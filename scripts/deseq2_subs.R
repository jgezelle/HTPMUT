#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  cat("Usage:\n  deseq2_substitutions.R <substitution_counts.tsv> <out_prefix>\n")
  quit(status=1)
}

in_tsv <- args[1]
out_prefix <- args[2]

dt <- fread(in_tsv)

# Expect columns from your parser:
# sample, code, virus, condition, replicate, contig, pos, ref, alt, alt_count, depth, feature
req <- c("sample","virus","condition","feature","alt_count","pos","ref","alt")
miss <- setdiff(req, names(dt))
if (length(miss) > 0) stop(paste("Missing columns:", paste(miss, collapse=", ")))

dt <- dt[virus %in% c("ZIKV","DENV") & condition %in% c("R","X")]

run_one <- function(vv) {
  d <- dt[virus == vv]

  # collapse any accidental duplicates by summing
  d <- d[, .(alt_count = sum(alt_count),
             pos = pos[1],
             ref = ref[1],
             alt = alt[1]),
         by=.(feature, sample, condition)]

  # build count matrix: rows = feature, cols = sample
  mat <- dcast(d, feature + pos + ref + alt ~ sample, value.var="alt_count", fill=0)

  feature_info <- mat[, .(feature, pos, ref, alt)]
  rownames_mat <- mat$feature
  count_mat <- as.matrix(mat[, -(1:4)])
  rownames(count_mat) <- rownames_mat

  # coldata
  samples <- colnames(count_mat)
  cond <- unique(d[, .(sample, condition)])[match(samples, sample), condition]
  coldata <- data.frame(row.names=samples, condition=factor(cond, levels=c("R","X")))

  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData=count_mat, colData=coldata, design=~condition)

  # Basic filtering: drop features with near-zero signal across all samples
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  dds <- DESeq(dds)

  res <- results(dds, contrast=c("condition","X","R"))
  resdf <- as.data.frame(res)
  resdf$feature <- rownames(resdf)

  # add pos/ref/alt back
  resdf <- merge(resdf, as.data.frame(feature_info), by="feature", all.x=TRUE)

  # write
  out_file <- paste0(out_prefix, "_", vv, "_DESeq2.tsv")
  fwrite(resdf, out_file, sep="\t")

  cat("Wrote:", out_file, " (n=", nrow(resdf), ")\n")
}

run_one("ZIKV")
run_one("DENV")