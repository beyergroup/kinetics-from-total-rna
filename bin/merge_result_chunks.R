#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)


parser <- ArgumentParser()
parser$add_argument("--input_folder",
                    required = TRUE,
                    help = "Folder containing TSV files to merge.")
parser$add_argument("--output_folder",
                    default = '.',
                    help = 'Path to output folder')

args <- parser$parse_args()

input_folder <- args$input_folder
output_folder <- args$output_folder

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# The fitting processes write one file per gene chunk, named <prefix>_<chunk_name>.tsv. Chunks
# holding no modelable genes produce an empty file, which read_tsv cannot type, so those are
# dropped rather than bound.
merge_chunks <- function(prefix, output_name, description) {
  chunk_files <- list.files(path = input_folder,
                            pattern = paste0("^", prefix, ".*\\.tsv$"),
                            full.names = TRUE)
  if (length(chunk_files) == 0) {
    stop("No matching TSV files with ", description, " found in the input folder: ", input_folder)
  }
  merged_df <- sort(chunk_files) |>
    map(read_tsv) |>
    keep(~nrow(.x) > 0) |>
    list_rbind()
  write_tsv(merged_df, file.path(output_folder, output_name))
}

merge_chunks("model_parameters", "model_parameters.tsv", "model parameters")
merge_chunks("ratio_lfcs", "ratio_lfcs.tsv", "rate ratio LFCs")
merge_chunks("test_results", "test_results_before_regularization.tsv", "test results")
