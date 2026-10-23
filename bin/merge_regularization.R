#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("--regularization_chunks_folder",
                    required = TRUE,
                    help = "Folder containing TSV files with regularization chunk results.")
parser$add_argument("--test_results",
                    required = TRUE,
                    help = "TSV file with test results before regularization.")
parser$add_argument("--output_folder",
                    default = '.',
                    help = 'Path to output folder')

args <- parser$parse_args()

output_folder <- args$output_folder
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

merge_chunks <- function(prefix, description) {
  chunk_files <- list.files(path = args$regularization_chunks_folder,
                            pattern = paste0("^", prefix, ".*\\.tsv$"),
                            full.names = TRUE)
  if (length(chunk_files) == 0) {
    stop("No matching TSV files with ", description, " found in the folder: ",
         args$regularization_chunks_folder)
  }
  sort(chunk_files) |>
    map(read_tsv) |>
    keep(~nrow(.x) > 0) |>
    list_rbind()
}

regularized_model_parameters_merged_df <- merge_chunks("regularized_model_parameters",
                                                       "regularized model parameters")
write_tsv(regularized_model_parameters_merged_df,
          file.path(output_folder, 'regularized_model_parameters.tsv'))

# Renamed before joining: test_results already carries an `lfc`, the unregularized one.
regularized_ratio_lfcs <- merge_chunks("regularized_ratio_lfcs",
                                       "regularized rate ratio LFCs") |>
  rename(regularized_lfc = lfc)

test_results <- read_tsv(args$test_results)

# The regularized effect size is a contrast, the difference between two design columns, while
# regularized_ratio_lfcs holds one column at a time. So each side of the contrast is joined
# separately and the two are subtracted, mirroring how Python assembles the unregularized `lfc`.
#
# A contrast names only one column when the other side is a factor's reference level, or when it
# is a continuous covariate. create_design_matrices.R writes an empty string there, which arrives
# as NA, and that side has to contribute zero. It does: feature_name is never NA in the ratio
# LFCs, so a NA key matches nothing and the replacement below takes over. The default NA matching
# has to stay, because intron_name is legitimately NA on both sides whenever a ratio is shared
# across a gene's introns.
#
# Joining against the ratio LFCs rather than the parameter table is what makes this correct for
# all six ratios: the parameter table has no row for the three that are not parameters, so the
# join used to yield NA for them and the replacement below turned it into a regularized LFC of
# exactly zero.
lfc_of_column <- function(column) {
  test_results |>
    left_join(regularized_ratio_lfcs,
              by = c("rate_ratio",
                     "gene_name",
                     "intron_name",
                     setNames("feature_name", column))) |>
    mutate(regularized_lfc = replace_na(regularized_lfc, 0.0)) |>
    pull(regularized_lfc)
}

test_results$lfc_regularized <- lfc_of_column("lfc_column_positive") -
  lfc_of_column("lfc_column_negative")

write_tsv(test_results, file.path(output_folder, 'raw_test_results.tsv'))