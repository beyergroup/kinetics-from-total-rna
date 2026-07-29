from enum import StrEnum
from typing import Optional, NamedTuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
import torch
import torch.nn as nn

from rna_kinetics.data import GeneData, GlobalGeneData


def safe_exp(x: torch.Tensor, output_threshold: float = 1e20) -> torch.Tensor:
    output_threshold_tensor = torch.tensor(output_threshold, dtype=x.dtype, device=x.device)
    input_threshold = torch.log(output_threshold_tensor)
    return torch.where(
        x <= input_threshold,
        torch.exp(x),
        output_threshold_tensor + output_threshold_tensor * (x - input_threshold)
    )


class TestableParameters(StrEnum):
    ALPHA = 'alpha'
    BETA = 'beta'
    GAMMA = 'gamma'


class LRTSpecification(NamedTuple):
    num_features_reduced_matrix: int
    tested_parameter: TestableParameters
    tested_intron: Optional[str] = None


def _fit_poisson_glm(X: np.ndarray, y: np.ndarray, offset: np.ndarray) -> np.ndarray:
    """
    Initialise alpha via Poisson GLM (log link, statsmodels IRLS).
    X: (num_samples, num_features), y: (num_samples,) or (num_samples, num_genes), offset: same shape as y.
    Returns alpha of shape (num_features,) or (num_genes, num_features).
    Falls back to zeros for any gene where the GLM fails (e.g. all-zero counts).
    """
    if y.ndim == 1:
        try:
            return sm.GLM(y, X, family=sm.families.Poisson(), offset=offset).fit(disp=False).params
        except Exception:
            return np.zeros(X.shape[1])
    G = y.shape[1]
    alpha = np.zeros((G, X.shape[1]))
    for g in range(G):
        try:
            alpha[g] = sm.GLM(y[:, g], X, family=sm.families.Poisson(),
                              offset=offset[:, g]).fit(disp=False).params
        except Exception:
            pass
    return alpha


class RNAKineticsModel(nn.Module):

    def __init__(self,
                 feature_names: list[str],
                 intron_names: list[str],
                 intron_specific_lfc: bool,
                 lrt_specification: Optional[LRTSpecification] = None
                 ):
        super().__init__()
        self.feature_names = feature_names
        self.intron_names = intron_names

        num_features = len(feature_names)
        num_introns = len(intron_names)

        self.lfc_init_over_deg = nn.Parameter(torch.zeros(num_features))
        self.intercept_exon = nn.Parameter(torch.zeros(1))

        if intron_specific_lfc:
            self.lfc_elong_over_deg = nn.Parameter(torch.zeros(num_features, num_introns))
            self.lfc_splice_over_deg = nn.Parameter(torch.zeros(num_features, num_introns))
        else:
            self.lfc_elong_over_deg = nn.Parameter(torch.zeros(num_features, 1))
            self.lfc_splice_over_deg = nn.Parameter(torch.zeros(num_features, 1))
        self.intron_specific_lfc = intron_specific_lfc

        self.intercept_intron = nn.Parameter(torch.zeros(num_introns))
        self.intercept_pi_logit = nn.Parameter(torch.zeros(num_introns))

        self.lrt_specification = lrt_specification
        self.tested_intron_index: Optional[int] = None
        if self.lrt_specification is not None:
            self.reduced_lfc = nn.Parameter(torch.zeros(self.lrt_specification.num_features_reduced_matrix))
            self.tested_intron_index: Optional[
                int] = None if not intron_specific_lfc or self.lrt_specification.tested_parameter == TestableParameters.ALPHA else self.intron_names.index(
                self.lrt_specification.tested_intron)

    def initialize_parameters(self,
                              gene_data: GeneData,
                              library_sizes: torch.Tensor,
                              design_matrix: torch.Tensor,
                              pi_eps: float = 0.01,
                              num_pi_grid_points: int = 20) -> None:
        with torch.no_grad():
            intercept_exon_scalar = torch.log(gene_data.exon_reads.mean() / library_sizes.mean())
            self.intercept_exon.data[:] = intercept_exon_scalar  # preserve shape [1]

            offset = (intercept_exon_scalar + torch.log(library_sizes) + gene_data.isoform_length_offset).cpu().numpy()
            glm_lfc = _fit_poisson_glm(design_matrix.cpu().numpy(), gene_data.exon_reads.cpu().numpy(), offset)
            self.lfc_init_over_deg.data.copy_(torch.from_numpy(glm_lfc).to(self.lfc_init_over_deg.dtype))

            aggregated_coverage = gene_data.coverage.sum(dim=0)
            pi_grid = torch.linspace(
                pi_eps,
                1 - pi_eps,
                num_pi_grid_points,
                device=aggregated_coverage.device,
                dtype=aggregated_coverage.dtype,
            )
            coverage_loss = CoverageLoss(num_position_coverage=aggregated_coverage.shape[1])
            coverage_loss_grid = coverage_loss.loss_for_pi_grid(pi_grid, aggregated_coverage)
            best_pi = pi_grid[coverage_loss_grid.argmin(dim=0)]
            self.intercept_pi_logit.data.copy_(torch.logit(best_pi, eps=pi_eps))

            intercept_intron_vector = torch.log(gene_data.intron_reads.mean(dim=0) / library_sizes.mean() / 2)
            self.intercept_intron.data.copy_(intercept_intron_vector)

    def forward(self,
                design_matrix: torch.Tensor,
                log_library_sizes: torch.Tensor,
                isoform_length_offset: torch.Tensor,
                reduced_design_matrix: Optional[torch.Tensor] = None):

        if self.lrt_specification is None:
            gene_expression_term = design_matrix @ self.lfc_init_over_deg
            speed_term = design_matrix @ self.lfc_elong_over_deg
            splicing_term = design_matrix @ self.lfc_splice_over_deg
        else:
            if reduced_design_matrix is None:
                raise ValueError(
                    "reduced_design_matrix must be provided in the LRT mode (i.e. when lrt_specification is not None.")
            gene_expression_term = reduced_design_matrix @ self.reduced_lfc if self.lrt_specification.tested_parameter == TestableParameters.ALPHA else design_matrix @ self.lfc_init_over_deg

            if self.intron_specific_lfc:
                speed_term = design_matrix @ self.lfc_elong_over_deg
                splicing_term = design_matrix @ self.lfc_splice_over_deg
                if self.lrt_specification.tested_parameter == TestableParameters.BETA:
                    speed_term[:, self.tested_intron_index] = reduced_design_matrix @ self.reduced_lfc
                elif self.lrt_specification.tested_parameter == TestableParameters.GAMMA:
                    splicing_term[:, self.tested_intron_index] = reduced_design_matrix @ self.reduced_lfc

            else:
                speed_term = (reduced_design_matrix @ self.reduced_lfc).unsqueeze(
                    1) if self.lrt_specification.tested_parameter == TestableParameters.BETA else design_matrix @ self.lfc_elong_over_deg
                splicing_term = (reduced_design_matrix @ self.reduced_lfc).unsqueeze(
                    1) if self.lrt_specification.tested_parameter == TestableParameters.GAMMA else design_matrix @ self.lfc_splice_over_deg

        predicted_log_reads_exon = self.intercept_exon + log_library_sizes + isoform_length_offset + gene_expression_term

        pi = torch.sigmoid(self.intercept_pi_logit - speed_term - splicing_term)

        intron_gene_expression_term = self.intercept_intron + log_library_sizes.unsqueeze(
            1) + gene_expression_term.unsqueeze(1)
        reads_intronic_polymerases = safe_exp(intron_gene_expression_term + self.intercept_pi_logit - speed_term)
        reads_unspliced_transcripts = safe_exp(intron_gene_expression_term + splicing_term)
        predicted_reads_intron = reads_intronic_polymerases + reads_unspliced_transcripts

        return safe_exp(predicted_log_reads_exon), predicted_reads_intron, pi

    def get_param_df(self) -> pd.DataFrame:
        model_parameters = dict(self.named_parameters())
        parameter_data: list[dict] = []
        for param_name, param_value in model_parameters.items():
            if param_name == 'intercept_exon':
                parameter_data.append({'parameter_type': param_name,
                                       'intron_name': None,
                                       'feature_name': None,
                                       'value': param_value.item()})
            elif param_name == 'alpha':
                for feature_index, feature_name in enumerate(self.feature_names):
                    parameter_data.append({'parameter_type': param_name,
                                           'intron_name': None,
                                           'feature_name': feature_name,
                                           'value': param_value[feature_index].item()})
            elif param_name in ('beta', 'gamma'):
                for feature_index, feature_name in enumerate(self.feature_names):
                    if self.intron_specific_lfc:
                        for intron_index, intron_name in enumerate(self.intron_names):
                            parameter_data.append({'parameter_type': param_name,
                                                   'intron_name': intron_name,
                                                   'feature_name': feature_name,
                                                   'value': param_value[feature_index, intron_index].item()})
                    else:
                        parameter_data.append({'parameter_type': param_name,
                                               'intron_name': None,
                                               'feature_name': feature_name,
                                               'value': param_value[feature_index, 0].item()})

            elif param_name in ('intercept_intron', 'theta'):
                for intron_index, intron_name in enumerate(self.intron_names):
                    parameter_data.append({'parameter_type': param_name,
                                           'intron_name': intron_name,
                                           'feature_name': None,
                                           'value': param_value[intron_index].item()})
            else:
                raise RuntimeError(f"Unexpected parameter name: {param_name}")
        df_param = pd.DataFrame(data=parameter_data)
        return df_param


class IntronCoverageModel(nn.Module):

    def __init__(self,
                 feature_names: list[str],
                 intron_names: list[str]):
        super().__init__()
        self.feature_names = feature_names
        self.intron_names = intron_names

        num_features = len(feature_names)
        num_introns = len(intron_names)

        self.lfc = nn.Parameter(torch.zeros(num_features, 1))
        self.theta = nn.Parameter(torch.zeros(num_introns))

    def initialize_theta(self,
                         coverage: torch.Tensor,
                         pi_eps: float = 0.01,
                         num_pi_grid_points: int = 20) -> None:
        aggregated_coverage = coverage.sum(dim=0)
        pi_grid = torch.linspace(
            pi_eps,
            1 - pi_eps,
            num_pi_grid_points,
            device=aggregated_coverage.device,
            dtype=aggregated_coverage.dtype,
        )
        coverage_loss = CoverageLoss(num_position_coverage=aggregated_coverage.shape[1])
        coverage_loss_grid = coverage_loss.loss_for_pi_grid(pi_grid, aggregated_coverage)
        best_pi = pi_grid[coverage_loss_grid.argmin(dim=0)]
        self.theta.data.copy_(torch.logit(best_pi, eps=pi_eps))

    def forward(self, design_matrix: torch.Tensor):
        logit_term = design_matrix @ self.lfc
        pi = torch.sigmoid(self.theta + logit_term)
        return pi

    def get_param_df(self) -> pd.DataFrame:
        model_parameters = dict(self.named_parameters())
        parameter_data: list[dict] = []
        for param_name, param_value in model_parameters.items():
            if param_name == 'lfc':
                for feature_index, feature_name in enumerate(self.feature_names):
                    parameter_data.append({'parameter_type': param_name,
                                           'intron_name': None,
                                           'feature_name': feature_name,
                                           'value': param_value[feature_index, 0].item()})
            elif param_name == 'theta':
                for intron_index, intron_name in enumerate(self.intron_names):
                    parameter_data.append({'parameter_type': param_name,
                                           'intron_name': intron_name,
                                           'feature_name': None,
                                           'value': param_value[intron_index].item()})
            else:
                raise RuntimeError(f"Unexpected parameter name: {param_name}")
        df_param = pd.DataFrame(data=parameter_data)
        return df_param


class GlobalPol2Model(nn.Module):

    def __init__(self,
                 feature_names: list[str],
                 gene_names: list[str],
                 intron_names: list[str],
                 gene_idx: torch.Tensor,
                 lrt_specification: Optional[LRTSpecification] = None,
                 ):
        super().__init__()
        self.feature_names = feature_names
        self.gene_names = gene_names
        self.intron_names = intron_names

        num_features = len(feature_names)
        num_genes = len(gene_names)
        num_introns = len(intron_names)

        self.alpha = nn.Parameter(torch.zeros(num_genes, num_features))
        self.beta = nn.Parameter(torch.zeros(num_features))
        self.gamma = nn.Parameter(torch.zeros(num_features))
        self.intercept_exon = nn.Parameter(torch.zeros(num_genes))
        self.intercept_intron = nn.Parameter(torch.zeros(num_introns))
        self.theta = nn.Parameter(torch.zeros(num_introns))

        self.register_buffer('gene_idx', gene_idx)

        self.lrt_specification = lrt_specification
        if lrt_specification is not None:
            if lrt_specification.tested_parameter not in (TestableParameters.BETA, TestableParameters.GAMMA):
                raise ValueError(
                    f"GlobalPol2Model LRT only supports BETA and GAMMA, got {lrt_specification.tested_parameter!r}.")
            self.reduced_lfc = nn.Parameter(torch.zeros(lrt_specification.num_features_reduced_matrix))

    def initialize_parameters(self,
                              global_gene_data: GlobalGeneData,
                              library_sizes: torch.Tensor,
                              design_matrix: torch.Tensor,
                              pi_eps: float = 0.01,
                              num_pi_grid_points: int = 20) -> None:
        with torch.no_grad():
            self.intercept_exon.data.copy_(
                torch.log(global_gene_data.exon_reads.mean(dim=0) / library_sizes.mean())
            )

            log_lib = torch.log(library_sizes)
            offset = (self.intercept_exon.unsqueeze(0) + log_lib.unsqueeze(1)
                      + global_gene_data.isoform_length_offset).cpu().numpy()  # (num_samples, num_genes)
            alpha_init = _fit_poisson_glm(
                design_matrix.cpu().numpy(), global_gene_data.exon_reads.cpu().numpy(), offset
            )  # (num_genes, num_features)
            self.alpha.data.copy_(torch.from_numpy(alpha_init).to(self.alpha.dtype))

            self.intercept_intron.data.copy_(
                torch.log(global_gene_data.intron_reads.mean(dim=0) / library_sizes.mean() / 2)
            )
            aggregated_coverage = global_gene_data.coverage.sum(dim=0)  # (num_introns, num_bins)
            pi_grid = torch.linspace(pi_eps, 1 - pi_eps, num_pi_grid_points,
                                     device=aggregated_coverage.device,
                                     dtype=aggregated_coverage.dtype)
            coverage_loss = CoverageLoss(num_position_coverage=aggregated_coverage.shape[1])
            coverage_loss_grid = coverage_loss.loss_for_pi_grid(pi_grid, aggregated_coverage)
            best_pi = pi_grid[coverage_loss_grid.argmin(dim=0)]
            self.theta.data.copy_(torch.logit(best_pi, eps=pi_eps))

    def forward(self,
                design_matrix: torch.Tensor,
                log_library_sizes: torch.Tensor,
                isoform_length_offset: torch.Tensor,
                reduced_design_matrix: Optional[torch.Tensor] = None):

        gene_expression_term = design_matrix @ self.alpha.T  # (num_samples, num_genes)

        predicted_log_reads_exon = (
                self.intercept_exon
                + log_library_sizes.unsqueeze(1)
                + isoform_length_offset
                + gene_expression_term
        )

        if self.lrt_specification is None:
            speed_term = (design_matrix @ self.beta).unsqueeze(1)
            splicing_term = (design_matrix @ self.gamma).unsqueeze(1)
        else:
            if reduced_design_matrix is None:
                raise ValueError("reduced_design_matrix must be provided in LRT mode.")
            reduced_term = (reduced_design_matrix @ self.reduced_lfc).unsqueeze(1)
            speed_term = (
                reduced_term if self.lrt_specification.tested_parameter == TestableParameters.BETA
                else (design_matrix @ self.beta).unsqueeze(1)
            )
            splicing_term = (
                reduced_term if self.lrt_specification.tested_parameter == TestableParameters.GAMMA
                else (design_matrix @ self.gamma).unsqueeze(1)
            )

        gene_expr_per_intron = gene_expression_term[:, self.gene_idx]  # (num_samples, num_introns)

        intron_gene_expression_term = (
                self.intercept_intron
                + log_library_sizes.unsqueeze(1)
                + gene_expr_per_intron
        )

        pi = torch.sigmoid(self.theta - speed_term - splicing_term)
        reads_intronic_polymerases = safe_exp(intron_gene_expression_term + self.theta - speed_term)
        reads_unspliced_transcripts = safe_exp(intron_gene_expression_term + splicing_term)
        predicted_reads_intron = reads_intronic_polymerases + reads_unspliced_transcripts

        return safe_exp(predicted_log_reads_exon), predicted_reads_intron, pi

    def get_param_df(self) -> pd.DataFrame:
        parameter_data: list[dict] = []
        for gene_index, gene_name in enumerate(self.gene_names):
            parameter_data.append({
                'parameter_type': 'intercept_exon',
                'gene_name': gene_name,
                'feature_name': None,
                'intron_name': None,
                'value': self.intercept_exon[gene_index].item(),
            })
            for feature_index, feature_name in enumerate(self.feature_names):
                parameter_data.append({
                    'parameter_type': 'alpha',
                    'gene_name': gene_name,
                    'feature_name': feature_name,
                    'intron_name': None,
                    'value': self.alpha[gene_index, feature_index].item(),
                })
        for feature_index, feature_name in enumerate(self.feature_names):
            for param_name in ('beta', 'gamma'):
                parameter_data.append({
                    'parameter_type': param_name,
                    'gene_name': None,
                    'feature_name': feature_name,
                    'intron_name': None,
                    'value': getattr(self, param_name)[feature_index].item(),
                })
        for intron_index, intron_name in enumerate(self.intron_names):
            for param_name in ('intercept_intron', 'theta'):
                parameter_data.append({
                    'parameter_type': param_name,
                    'gene_name': None,
                    'feature_name': None,
                    'intron_name': intron_name,
                    'value': getattr(self, param_name)[intron_index].item(),
                })
        return pd.DataFrame(parameter_data)


class CoverageLoss(nn.Module):

    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        locations = torch.linspace(start=1 / (2 * num_position_coverage),
                                   end=1 - 1 / (2 * num_position_coverage),
                                   steps=num_position_coverage)
        location_term = 1 - 2 * locations
        self.register_buffer("location_term", location_term)

    def forward(self, pi, coverage):
        loss_per_location = -torch.log(1 + pi.unsqueeze(2) * self.location_term)
        return torch.sum(loss_per_location * coverage)

    def loss_for_pi_grid(self, candidate_pi: torch.Tensor, aggregated_coverage: torch.Tensor) -> torch.Tensor:
        """
        Compute coverage loss for a grid of candidate pi values.
        Args:
            candidate_pi: Tensor of shape (k,), candidate values of pi.
            aggregated_coverage: Tensor of shape (i, b), coverage summed over samples.
        Returns:
            Tensor of shape (k, i) containing the total loss for each candidate pi
            and each intron.
        """
        loss_terms = -torch.log(
            1 + candidate_pi[:, None] * self.location_term[None, :]
        )  # shape (k, b)
        return torch.einsum("kb,ib->ki", loss_terms, aggregated_coverage)


class RNAKineticsLoss(nn.Module):
    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        self.loss_function_exon_counts = nn.PoissonNLLLoss(log_input=False, full=True, reduction='sum')
        self.loss_function_intron_counts = nn.PoissonNLLLoss(log_input=False, full=True, reduction='sum')
        self.loss_function_intron_coverage = CoverageLoss(num_position_coverage=num_position_coverage)

    def forward(self,
                reads_exon: torch.Tensor,
                reads_intron: torch.Tensor,
                intron_coverage: torch.Tensor,
                predicted_reads_exon: torch.Tensor,
                predicted_reads_intron: torch.Tensor,
                predicted_pi: torch.Tensor):
        loss_exon_counts = self.loss_function_exon_counts(predicted_reads_exon,
                                                          reads_exon)
        loss_intron_counts = self.loss_function_intron_counts(predicted_reads_intron, reads_intron)
        loss_intron_coverage = self.loss_function_intron_coverage(predicted_pi, intron_coverage)

        total_loss = loss_exon_counts + loss_intron_counts + loss_intron_coverage
        return total_loss
