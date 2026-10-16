#!/usr/bin/env Rscript

# Estimates the prior covariance of the rate ratio LFCs, per design column, and converts it into
# the L2 penalty weights the regularized refit applies.
#
# The penalty the refit adds is `0.5 * t(theta) %*% Lambda %*% theta` per design column, where
# `theta` holds that column's three parametrized ratios. That is the same thing as a zero-centred
# Gaussian prior on `theta` with covariance `Lambda^-1`, so what has to be estimated is that
# covariance. The matching estimator is a multivariate random-effects model over genes, fitted by
# REML: each gene's three estimates are a draw from the prior, observed with that gene's own
# sampling covariance, as `y_g ~ N(mu, Psi + V_g)`.
#
# REML is restricted, or residual, maximum likelihood (Patterson & Thompson 1971, Biometrika
# 58:545-554). It maximises the likelihood of residual contrasts, combinations of the data whose
# distribution does not involve the fixed effects, so the estimated covariance is corrected for the
# degrees of freedom spent estimating the mean; in one dimension it is the difference between
# dividing by n and by n - 1. Plain maximum likelihood is biased low for a variance component and
# can push one onto the zero boundary, which is why "reml" is passed rather than "ml".
#
# A gene whose estimates are explained by noise contributes almost nothing, which is why there is
# no bound on the standard errors: the fit discounts an imprecise estimate through that weight
# without needing a threshold, and stays stable even for standard errors many orders of magnitude
# above the rest. A finite, strictly positive standard error is still required: a zero sampling
# variance would be an infinite weight, and an infinite one makes the fit fail outright.
#
# Fitting one joint covariance rather than six separate per-ratio priors is what keeps the result
# coherent. Six independently fitted marginals need not describe any joint distribution, and on
# intron-specific fits they do not: they can imply a correlation above one. Every ratio's prior is
# instead read off the single fitted matrix, so consistency holds by construction.
#
# Notation: `Psi`, `tau`, `mu` are mixmeta's, for the between-study covariance, its standard
# deviations, and the pooled mean of a random-effects meta-analysis (Sera, Armstrong, Blangiardo &
# Gasparrini 2019, Stat. Med. 38:5429-5444); genes play the role of studies. `lambda`, `prior_sd`
# and `prior_cov_*` are this pipeline's own names and are columns of the output. Nothing here uses
# the per-rate symbols the model itself is written with.

library(argparse)
library(tidyverse)
library(mixmeta)

parser <- ArgumentParser()
parser$add_argument("--ratio_lfcs",
                    required = TRUE,
                    help = "TSV with one row per rate ratio, intron and design column.")
parser$add_argument("--output_folder", default = ".",
                    help = "Where to write the output TSV with regularization coefficients.")
parser$add_argument("--max_abs_lfc", type = "double", default = 10 * log(2),
                    help = paste("Estimates with |lfc| above this (natural log scale) do not",
                                 "contribute to the prior fit. They are still regularized: the",
                                 "fitted prior is applied to every gene regardless. The default is",
                                 "a 1024-fold change in the rate ratio, the same bound DESeq2",
                                 "applies as abs(x) < 10 on its log2 scale before estimating its",
                                 "own coefficient prior. Beyond it the estimates come from fits",
                                 "whose likelihood has no interior maximum rather than from",
                                 "measured effects."))

args <- parser$parse_args()

output_folder <- args$output_folder
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# The basis the penalty acts on. The other three ratios are differences of these, so their priors
# follow from the same matrix instead of being fitted.
BASIS <- c("lfc_init_over_deg", "lfc_elong_over_deg", "lfc_splice_over_deg")

# Each remaining ratio as the pair of basis ratios it is the difference of. This is what turns the
# derived ratios' reported standard errors into the sampling covariances between the basis
# estimates, through cov(a, b) = (var(a) + var(b) - var(a - b)) / 2. All three variances are
# marginals of one gene's inverse Fisher information, so the reconstruction is exact and the
# resulting block is positive definite.
DIFFERENCE_OF <- list(lfc_init_over_elong   = c("lfc_init_over_deg", "lfc_elong_over_deg"),
                      lfc_init_over_splice  = c("lfc_init_over_deg", "lfc_splice_over_deg"),
                      lfc_elong_over_splice = c("lfc_elong_over_deg", "lfc_splice_over_deg"))

# Every ratio as a contrast of the basis, which is what lets one fitted matrix yield all six
# priors.
CONTRAST <- c(lapply(set_names(seq_along(BASIS), BASIS), function(i) as.numeric(seq_along(BASIS) == i)),
              lapply(DIFFERENCE_OF,
                     function(pair) as.numeric(BASIS == pair[1]) - as.numeric(BASIS == pair[2])))

# With no usable estimates there is nothing to learn a prior from, so fall back to a prior so wide
# that it imposes no shrinkage. DESeq2 does the same, returning a prior variance of 1e6 when its
# own finiteness filter empties the input.
UNINFORMATIVE_PRIOR_VARIANCE <- 1e6

# REML needs at least two genes to separate between-gene covariance from sampling covariance.
MIN_UNITS_FOR_FIT <- 2

ratio_lfcs <- read_tsv(args$ratio_lfcs, show_col_types = FALSE) |>
  mutate(usable = lfc_identifiable &
           !training_diverged_full_model &
           is.finite(lfc) &
           is.finite(lfc_se) &
           lfc_se > 0 &
           abs(lfc) <= args$max_abs_lfc)

# One row per unit of the fit, with every ratio's estimate and variance in its own column. A unit
# is a (gene, intron) pair when any ratio is intron-specific, and a gene otherwise. A ratio the
# model shares across a gene's introns is reported once with a null intron_name, so it joins on the
# gene alone and is broadcast onto each of that gene's introns.
build_fit_table <- function(df_feature) {
  intron_level <- df_feature |> filter(!is.na(intron_name))
  units <- if (nrow(intron_level) > 0) {
    distinct(intron_level, gene_name, intron_name)
  } else {
    distinct(df_feature, gene_name) |> mutate(intron_name = NA_character_)
  }
  for (ratio_name in names(CONTRAST)) {
    one_ratio <- df_feature |>
      filter(ratio == ratio_name) |>
      transmute(gene_name,
                intron_name,
                !!paste0("lfc_of_", ratio_name) := if_else(usable, lfc, NA_real_),
                !!paste0("var_of_", ratio_name) := if_else(usable, lfc_se^2, NA_real_))
    is_shared <- all(is.na(one_ratio$intron_name))
    units <- left_join(units,
                       if (is_shared) select(one_ratio, -intron_name) else one_ratio,
                       by = if (is_shared) "gene_name" else c("gene_name", "intron_name"))
  }
  units
}

# The estimates and their sampling covariances, in the form mixmeta takes.
#
# A basis ratio is dropped for a unit when its own estimate is unusable, and also when a covariance
# it needs is unavailable, which happens only where the bound excluded a difference while sparing
# both its terms. Both members of such a pair are dropped, which needs no tie-break. mixmeta uses
# whichever outcomes remain, so a partly dropped unit still contributes.
build_estimates_and_covariances <- function(fit_table) {
  num_units <- nrow(fit_table)
  estimates <- as.matrix(fit_table[, paste0("lfc_of_", BASIS)])
  colnames(estimates) <- BASIS
  variances <- as.matrix(fit_table[, paste0("var_of_", BASIS)])
  covariances <- matrix(NA_real_, num_units, length(DIFFERENCE_OF))
  num_pairs_dropped <- 0

  for (position in seq_along(DIFFERENCE_OF)) {
    pair <- DIFFERENCE_OF[[position]]
    first <- match(pair[1], BASIS)
    second <- match(pair[2], BASIS)
    difference_variance <- fit_table[[paste0("var_of_", names(DIFFERENCE_OF)[position])]]
    covariances[, position] <- (variances[, first] + variances[, second] - difference_variance) / 2
    unavailable <- is.na(covariances[, position]) &
      !is.na(variances[, first]) & !is.na(variances[, second])
    num_pairs_dropped <- num_pairs_dropped + sum(unavailable)
    estimates[unavailable, c(first, second)] <- NA_real_
    variances[unavailable, c(first, second)] <- NA_real_
  }

  # Entries belonging to a dropped outcome are never read, since mixmeta removes that outcome's row
  # and column, so filling them with zero is only to keep each block well formed.
  sampling_covariances <- lapply(seq_len(num_units), function(unit) {
    block <- diag(replace_na(variances[unit, ], 0), nrow = length(BASIS))
    for (position in seq_along(DIFFERENCE_OF)) {
      first <- match(DIFFERENCE_OF[[position]][1], BASIS)
      second <- match(DIFFERENCE_OF[[position]][2], BASIS)
      block[first, second] <- block[second, first] <- replace_na(covariances[unit, position], 0)
    }
    block
  })
  list(estimates = estimates,
       sampling_covariances = sampling_covariances,
       num_pairs_dropped = num_pairs_dropped)
}

# The fitted prior for one design column, as one row per rate ratio.
#
# A basis ratio with no usable estimates anywhere carries no information about its own direction,
# and mixmeta cannot estimate a variance for it. Those directions are dropped and the reduced model
# is fitted on the rest, so the ratios lying in the fitted subspace still get a prior and only the
# ones needing a dropped direction fall back.
estimate_prior_for_feature <- function(df_feature) {
  feature <- df_feature$feature_name[1]
  fit_table <- build_fit_table(df_feature)
  prepared <- build_estimates_and_covariances(fit_table)
  estimates <- prepared$estimates
  is_fittable <- colSums(!is.na(estimates)) >= MIN_UNITS_FOR_FIT

  counts <- tibble(num_units_used = sum(rowSums(!is.na(estimates)) > 0),
                   num_units_total = nrow(fit_table),
                   num_estimates_used = sum(!is.na(estimates)),
                   num_pairs_dropped = prepared$num_pairs_dropped)
  uninformative <- function(ratio_names, status) {
    tibble(feature_name = feature,
           parameter_type = ratio_names,
           prior_sd = sqrt(UNINFORMATIVE_PRIOR_VARIANCE),
           tau = NA_real_,
           mu = NA_real_,
           prior_cov_init_over_deg = NA_real_,
           prior_cov_elong_over_deg = NA_real_,
           prior_cov_splice_over_deg = NA_real_,
           status = status) |>
      bind_cols(counts)
  }
  if (!any(is_fittable)) {
    return(uninformative(names(CONTRAST), "no_usable_estimates"))
  }

  fitted_index <- which(is_fittable)
  fit <- tryCatch(mixmeta(estimates[, fitted_index, drop = FALSE] ~ 1,
                          S = lapply(prepared$sampling_covariances,
                                     function(block) block[fitted_index, fitted_index, drop = FALSE]),
                          method = "reml"),
                  error = function(e) NULL)
  if (is.null(fit)) {
    return(uninformative(names(CONTRAST), "fit_failed"))
  }

  between_gene_covariance <- as.matrix(fit$Psi)
  pooled_mean <- as.numeric(coef(fit))
  # The penalty shrinks each LFC toward zero rather than toward the fitted mean, so the prior
  # covariance is the second moment about zero, Psi + mu %*% t(mu). Shrinking toward the fitted
  # mean would break two things: the effect size is reported next to an LRT p-value for an LFC of
  # zero, so a gene with no effect would carry a non-zero regularized one, and a mean fitted per
  # ratio is not consistent with the differences between ratios. `mu` is reported because a large
  # value means either a genuine genome-wide shift or a normalisation problem.
  prior_covariance <- between_gene_covariance + outer(pooled_mean, pooled_mean)
  is_singular <- min(eigen(prior_covariance, only.values = TRUE)$values) <= 0

  rows <- map_dfr(names(CONTRAST), function(ratio_name) {
    contrast <- CONTRAST[[ratio_name]]
    # A ratio is estimable only if it puts no weight on a dropped direction.
    if (any(contrast[!is_fittable] != 0)) {
      return(uninformative(ratio_name, "not_estimable"))
    }
    fitted_contrast <- contrast[fitted_index]
    covariance_with_basis <- set_names(rep(NA_real_, length(BASIS)), BASIS)
    covariance_with_basis[fitted_index] <- as.numeric(fitted_contrast %*% prior_covariance)
    tibble(feature_name = feature,
           parameter_type = ratio_name,
           prior_sd = sqrt(as.numeric(t(fitted_contrast) %*% prior_covariance %*% fitted_contrast)),
           tau = sqrt(as.numeric(t(fitted_contrast) %*% between_gene_covariance %*% fitted_contrast)),
           mu = sum(fitted_contrast * pooled_mean),
           prior_cov_init_over_deg = covariance_with_basis[[1]],
           prior_cov_elong_over_deg = covariance_with_basis[[2]],
           prior_cov_splice_over_deg = covariance_with_basis[[3]],
           # A prior of zero width in a direction asserts those effects are exactly zero, which
           # gives an infinite penalty and a non-finite loss in the refit. It is a failed fit, not
           # a finding.
           status = if (is_singular || !is.finite(prior_sd) || prior_sd <= 0)
             "degenerate_prior" else "ok") |>
      bind_cols(counts)
  })
  rows
}

regularization_coefficients_df <- ratio_lfcs |>
  group_split(feature_name) |>
  map_dfr(estimate_prior_for_feature) |>
  mutate(lambda = 1 / (2 * prior_sd^2)) |>
  select(parameter_type, feature_name, prior_sd, lambda, tau, mu,
         starts_with("prior_cov_"), status, starts_with("num_")) |>
  arrange(parameter_type, feature_name)

degenerate <- regularization_coefficients_df |>
  filter(status == "degenerate_prior") |>
  distinct(feature_name)
if (nrow(degenerate) > 0) {
  stop("Prior covariance estimated as singular for design column(s): ",
       paste(degenerate$feature_name, collapse = ", "),
       ". This would put an infinite penalty on those coefficients.")
}

write_tsv(regularization_coefficients_df,
          file.path(output_folder, "regularization_coefficients.tsv"))

print(as.data.frame(regularization_coefficients_df), digits = 4)
