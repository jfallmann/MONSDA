#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(R.utils)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript CPM_TPM_from_counts.R counts.tsv.gz")
}

infile <- args[1]
cpmoutfile <- args[2]
tpmoutfile <- args[3]

# ---- Read data ----
message("Reading ", infile)
counts <- read.delim(infile, check.names = FALSE)

# Extract columns
gene_lengths <- counts$Length
gene_ids <- counts$ID
count_matrix <- counts[ , -(1:2)]

# ---- CPM calculation ----
# CPM = (counts / library_size) * 1e6
lib_sizes <- colSums(count_matrix)
cpm <- sweep(count_matrix, 2, lib_sizes, FUN = "/") * 1e6

# ---- TPM calculation ----
# 1. Normalize counts by gene length in kb
rate <- sweep(count_matrix, 1, gene_lengths / 1000, FUN = "/")

# 2. Scale so that sum of rates per sample = 1e6
scaling_factors <- colSums(rate)
tpm <- sweep(rate, 2, scaling_factors, FUN = "/") * 1e6

# ---- Add back gene IDs ----
cpm_out <- data.frame(ID = gene_ids, cpm, check.names = FALSE)
tpm_out <- data.frame(ID = gene_ids, tpm, check.names = FALSE)

# ---- Write output ----
write.table(cpm_out, gzfile(cpmoutfile), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_out, gzfile(tpmoutfile), sep = "\t", quote = FALSE, row.names = FALSE)
