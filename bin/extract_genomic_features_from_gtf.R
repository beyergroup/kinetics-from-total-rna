#!/usr/bin/env Rscript

library(argparse)
library(rtracklayer)
library(GenomicRanges)
library(BiocParallel)
library(readr)

parser <- ArgumentParser()
parser$add_argument("--gtf",
                    type = "character",
                    required = TRUE,
                    help = "Path to the input GTF file.")
parser$add_argument("--threads",
                    type = "integer",
                    default = 4,
                    help = "Number of threads used by BiocParallel.")
parser$add_argument("--output_folder",
                    type = "character",
                    default = '.')
parser$add_argument("--min_intron_length",
                    type = "integer",
                    default = 50)
parser$add_argument("--min_constitutive_exon_length",
                    type = "integer",
                    default = 20)
parser$add_argument("--tss_blacklist_margin",
                    type = "integer",
                    default = 50,
                    help = paste("Blacklist an intron if its 3' end lies within this many bp",
                                 "upstream of an alternative transcription start site (first exon",
                                 "of another isoform of the same gene). Such introns carry a",
                                 "spurious 3' coverage spike from alt-promoter Pol II loading.",
                                 "Set to 0 to disable blacklisting."))

args <- parser$parse_args()

register(SnowParam(workers = args$threads, type = "SOCK", progressbar = FALSE))

output_folder <- args$output_folder
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

gtf <- rtracklayer::import(args$gtf)

genes <- gtf[gtf$type == "gene"]
exons <- gtf[gtf$type == "exon"]

genes <- genes[!is.na(genes$gene_id)]
exons <- exons[!is.na(exons$gene_id)]

if (anyDuplicated(genes$gene_id) > 0) {
  stop("Duplicated gene_id entries found in 'gene' features.")
}

# Extract gene IDs
gene_type_column <- if ("gene_biotype" %in% names(mcols(genes))) {
  "gene_biotype"
} else if ("gene_type" %in% names(mcols(genes))) {
  "gene_type"
} else {
  stop("Neither 'gene_biotype' nor 'gene_type' found in GTF attributes.")
}

all_gene_ids <- sort(unique(genes$gene_id))
protein_coding_gene_ids <- sort(unique(genes$gene_id[mcols(genes)[[gene_type_column]] == "protein_coding"]))

readr::write_csv(data.frame(gene_id = all_gene_ids),
                 file.path(output_folder, "all_genes.csv"))
readr::write_csv(data.frame(gene_id = protein_coding_gene_ids),
                 file.path(output_folder, "protein_coding_genes.csv"))

# Different annotations (Ensemble vs. Gencode) use different naming for UTRs
utr_types <- intersect(c("five_prime_utr", "three_prime_utr", "UTR"), unique(gtf$type))

exons_and_utr <- gtf[gtf$type %in% c("exon", utr_types)]
exons_and_utr <- exons_and_utr[!is.na(exons_and_utr$gene_id)]

exons_by_gene <- split(exons, exons$gene_id)
gene_ranges_by_gene <- split(genes, genes$gene_id)
exons_and_utr_by_gene <- split(exons_and_utr, exons_and_utr$gene_id)

# Transcription start sites (5' end of each transcript), used to blacklist introns whose
# 3' boundary abuts an alternative TSS. When a gene has an alternative promoter, the first
# exon of a shorter isoform starts inside the gene body; the setdiff-based intron of the
# longer isoform then ends at that TSS (a non-splice boundary), and alt-promoter Pol II
# loading / 5' heterogeneity piles reads up just upstream of it, producing a spurious 3' spike.
exons_with_transcript <- exons[!is.na(exons$transcript_id)]
transcript_spans <- unlist(range(split(exons_with_transcript, exons_with_transcript$transcript_id)))
tss_points <- resize(transcript_spans, width = 1L, fix = "start")
tss_points$gene_id <- exons_with_transcript$gene_id[
  match(names(tss_points), exons_with_transcript$transcript_id)]

# Extract introns
introns_by_gene <- bplapply(
  names(gene_ranges_by_gene),
  function(gene_id, gene_ranges_by_gene, exons_and_utr_by_gene) {

    gene_range <- range(gene_ranges_by_gene[[gene_id]])
    exons_and_utr_ranges <- exons_and_utr_by_gene[[gene_id]]

    if (is.null(exons_and_utr_ranges) || length(exons_and_utr_ranges) == 0) {
      warning("No exon/UTR for gene ", gene_id)
      return(GRanges())
    }

    exons_and_utr_ranges <- GenomicRanges::reduce(exons_and_utr_ranges, ignore.strand = FALSE)
    GenomicRanges::setdiff(gene_range, exons_and_utr_ranges, ignore.strand = FALSE)
  },
  gene_ranges_by_gene = gene_ranges_by_gene,
  exons_and_utr_by_gene = exons_and_utr_by_gene
)
names(introns_by_gene) <- names(gene_ranges_by_gene)


# Extract constitutive exons
constitutive_exons_by_gene <- bplapply(exons_by_gene, function(gene_exons) {
  if (length(gene_exons) == 0) return(GRanges())

  transcripts_ids <- unique(gene_exons$transcript_id)
  num_transcripts <- length(transcripts_ids)

  atomic_exon_ranges <- disjoin(gene_exons, with.revmap = TRUE, ignore.strand = FALSE)

  atomic_exon_to_transcripts <- mcols(atomic_exon_ranges)$revmap

  keep_mask <- vapply(atomic_exon_to_transcripts, function(exon_transcript_indices) {
    length(unique(gene_exons$transcript_id[exon_transcript_indices])) == num_transcripts
  }, logical(1))
  atomic_exon_ranges[keep_mask]
})


constitutive_exons <- unlist(GRangesList(constitutive_exons_by_gene), use.names = TRUE)
constitutive_exons$gene_id <- names(constitutive_exons)
names(constitutive_exons) <- NULL
constitutive_exons <- constitutive_exons[width(constitutive_exons) >= args$min_constitutive_exon_length]

introns <- unlist(GRangesList(introns_by_gene), use.names = TRUE)
introns$gene_id <- names(introns)
names(introns) <- NULL
introns <- introns[width(introns) >= args$min_intron_length]


features_list <- list(
  introns = introns,
  constitutive_exons = constitutive_exons
)

for (feature_name in names(features_list)) {
  features <- features_list[[feature_name]]
  num_overlapping_genes <- countOverlaps(features, genes, ignore.strand = FALSE)
  features <- features[num_overlapping_genes == 1]

  # Assign number in the 5' to 3' direction
  features <- features[order(features$gene_id, start(features), end(features))]
  feature_index_in_gene <- ave(seq_along(features),
                               features$gene_id,
                               FUN = seq_along)
  num_features_in_gene <- ave(seq_along(features),
                              features$gene_id,
                              FUN = length)
  mask_minus_strand <- as.character(strand(features)) == "-"

  features$number_in_gene_5_to_3_prime <- feature_index_in_gene
  features$number_in_gene_5_to_3_prime[mask_minus_strand] <- num_features_in_gene[mask_minus_strand] + 1 - feature_index_in_gene[mask_minus_strand]
  features$name <- paste0(features$gene_id, "_", features$number_in_gene_5_to_3_prime)

  feature_genes <- genes[match(features$gene_id, genes$gene_id)]

  features$mid <- (start(features) + end(features)) / 2
  features$mid_relative_5_to_3_prime <- (features$mid - start(feature_genes)) / width(feature_genes)
  features$mid_relative_5_to_3_prime[mask_minus_strand] <- 1 - features$mid_relative_5_to_3_prime[mask_minus_strand]


  features$score <- 0
  features <- sort(features)


  if (feature_name == "constitutive_exons") {
    mcols(features)$constitutive_exon_id <- features$name
    features$type <- "constitutive_exon"

    rtracklayer::export(features, file.path(output_folder, "constitutive_exons.bed"), format = "BED")
    rtracklayer::export(features, file.path(output_folder, "constitutive_exons.gtf"), format = "GTF")

  } else if (feature_name == "introns") {
    mcols(features)$intron_id <- features$name
    features$type <- "intron"

    # Flag introns whose 3' end lies within tss_blacklist_margin bp upstream of an alternative
    # TSS of the same gene (see TSS computation above). The intron numbering / IDs are assigned
    # over the full set before this split, so kept introns keep stable positional IDs (the
    # numbering may have gaps where blacklisted introns were removed).
    blacklist_reason <- rep(NA_character_, length(features))
    if (args$tss_blacklist_margin > 0) {
      intron_downstream <- flank(features, width = args$tss_blacklist_margin, start = FALSE)
      hits <- findOverlaps(intron_downstream, tss_points, ignore.strand = FALSE)
      hits <- hits[features$gene_id[queryHits(hits)] == tss_points$gene_id[subjectHits(hits)]]
      if (length(hits) > 0) {
        reason_by_intron <- tapply(names(tss_points)[subjectHits(hits)],
                                   queryHits(hits),
                                   function(tx) paste(sort(unique(tx)), collapse = ","))
        blacklist_reason[as.integer(names(reason_by_intron))] <- unname(reason_by_intron)
      }
    }

    blacklist_mask <- !is.na(blacklist_reason)
    kept_introns <- features[!blacklist_mask]
    blacklisted_introns <- features[blacklist_mask]
    blacklisted_introns$alt_tss_transcript <- blacklist_reason[blacklist_mask]

    message(sprintf(
      "Introns: %d total, %d blacklisted (alt-TSS within %d bp of 3' end), %d kept.",
      length(features), length(blacklisted_introns), args$tss_blacklist_margin, length(kept_introns)))

    rtracklayer::export(kept_introns, file.path(output_folder, "introns.bed"), format = "BED")
    rtracklayer::export(kept_introns, file.path(output_folder, "introns.gtf"), format = "GTF")
    rtracklayer::export(blacklisted_introns, file.path(output_folder, "blacklisted_introns.bed"), format = "BED")
    rtracklayer::export(blacklisted_introns, file.path(output_folder, "blacklisted_introns.gtf"), format = "GTF")
  }

}

