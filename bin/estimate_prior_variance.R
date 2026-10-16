#!/usr/bin/env Rscript

# Estimates the prior variance of each rate ratio's LFC, per design column, and converts it into
# the L2 penalty weight the regularized refit applies.
#
# The penalty the refit adds is `lambda * lfc^2` per design column, which is the same thing as a
# zero-centred Gaussian prior on that LFC with variance `1 / (2 * lambda)`. So what has to be
# estimated is the variance of that Gaussian, and the matching estimator is the random-effects
# variance of the per-gene estimates, fitted by REML: each gene's LFC is a draw from the prior,
# observed with its own sampling error, as `lfc_g ~ N(mu, tau^2 + lfc_se_g^2)`. A gene whose
# estimate is explained by noise therefore contributes almost nothing, which is why there is no
# bound on `lfc_se`: REML discounts an imprecise estimate through that weight without needing a
# threshold; it stays stable even for standard errors many orders of magnitude above the rest.
# A finite, strictly positive `lfc_se` is still required: a zero sampling variance would be an
# infinite weight, and an infinite one makes rma fail outright.
#
# Notation: `tau` and `mu` are metafor's, for the between-study standard deviation and the pooled
# mean of a random-effects meta-analysis (Viechtbauer 2010, J. Stat. Softw. 36(3)); genes play the
# role of studies. `lambda` and `prior_sd` are this pipeline's own names, and both are columns of
# the output. Nothing here uses the per-rate symbols the model itself is written with.

library(argparse)
library(tidyverse)
library(metafor)

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

# With no informative estimates there is nothing to learn a prior from, so fall back to a prior so
# wide that it imposes no shrinkage. DESeq2 does the same, returning a prior variance of 1e6 when
# its own finiteness filter empties the input.
UNINFORMATIVE_PRIOR_VARIANCE <- 1e6

# REML needs at least two estimates to separate between-gene variance from sampling variance.
MIN_ESTIMATES_FOR_FIT <- 2

ratio_lfcs <- read_tsv(args$ratio_lfcs, show_col_types = FALSE)

estimate_prior_variance <- function(df_group) {
  usable <- df_group |>
    filter(lfc_identifiable,
           !training_diverged_full_model,
           is.finite(lfc),
           is.finite(lfc_se),
           lfc_se > 0)
  num_before_bound <- nrow(usable)
  usable <- usable |>
    filter(abs(lfc) <= args$max_abs_lfc)

  # The fallback values go in first so that an early return still carries a status. `status` is
  # a column rather than a warning because it qualifies one specific lambda: anyone reading that
  # number needs to know, alongside it, whether a fit succeeded at all.
  result <- tibble(num_estimates_used = nrow(usable),
                   num_estimates_above_lfc_bound = num_before_bound - nrow(usable),
                   num_estimates_total = nrow(df_group),
                   tau = NA_real_,
                   mu = NA_real_,
                   prior_sd = sqrt(UNINFORMATIVE_PRIOR_VARIANCE),
                   status = "no_informative_estimates")
  if (nrow(usable) < MIN_ESTIMATES_FOR_FIT) {
    return(result)
  }

  fit <- tryCatch(rma(yi = usable$lfc, vi = usable$lfc_se^2, method = "REML"),
                  error = function(e) NULL)
  if (is.null(fit)) {
    result$status <- "fit_failed"
    return(result)
  }

  # rma fits a pooled mean as well as a between-gene variance, but the penalty shrinks each LFC
  # toward zero rather than toward that mean, so the variance the prior needs is the second moment
  # about zero, `tau^2 + mu^2`. Shrinking toward the fitted mean would break two things: the
  # effect size is reported next to an LRT p-value for an LFC of zero, so a gene with no effect
  # would carry a non-zero regularized one, and the means fitted per ratio are not mutually
  # consistent, so the derived ratios would stop being differences of the parametrized ones. `mu`
  # is reported because a large value means either a genuine genome-wide shift or a normalisation
  # problem.
  result$tau <- sqrt(fit$tau2)
  result$mu <- as.numeric(coef(fit))
  result$prior_sd <- sqrt(fit$tau2 + result$mu^2)
  result$status <- "ok"

  # A prior of exactly zero width asserts that every effect is exactly zero, which gives an
  # infinite penalty and a non-finite loss in the refit. It is a failed fit, not a finding.
  if (!is.finite(result$prior_sd) || result$prior_sd <= 0) {
    result$status <- "degenerate_prior"
  }
  result
}

regularization_coefficients_df <- ratio_lfcs |>
  group_by(ratio, feature_name) |>
  group_modify(~estimate_prior_variance(.x)) |>
  ungroup() |>
  mutate(lambda = 1 / (2 * prior_sd^2)) |>
  # The refit joins on parameter_type, the name this vocabulary carries everywhere except in
  # ratio_lfcs.tsv.
  rename(parameter_type = ratio) |>
  select(parameter_type, feature_name, prior_sd, lambda, tau, mu, status,
         num_estimates_used, num_estimates_above_lfc_bound, num_estimates_total) |>
  arrange(parameter_type, feature_name)

degenerate <- regularization_coefficients_df |> filter(status == "degenerate_prior")
if (nrow(degenerate) > 0) {
  stop("Prior variance estimated as zero for: ",
       paste(degenerate$parameter_type, degenerate$feature_name, sep = "/", collapse = ", "),
       ". This would put an infinite penalty on those coefficients.")
}

write_tsv(regularization_coefficients_df,
          file.path(output_folder, "regularization_coefficients.tsv"))

print(as.data.frame(regularization_coefficients_df), digits = 4)
