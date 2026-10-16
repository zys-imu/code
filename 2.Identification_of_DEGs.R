#!/usr/bin/env Rscript

# Differential analysis for HepG2 vs hepatocyte.
# Expected layout under project_dir:
#   input/hepatocyte_1.bed, input/hepatocyte_2.bed
#   input/HepG2_1.bed, input/HepG2_2.bed
#   result/DEGs/

options(stringsAsFactors = FALSE, warn = 1)
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("R/Bioconductor package 'DESeq2' is required.", call. = FALSE)
}
suppressPackageStartupMessages(library(DESeq2))

# Generic path; change this or use --project-dir.
project_dir <- "path/to/project"
input_dir <- file.path(project_dir, "input")
result_dir <- file.path(project_dir, "result", "DEGs")
sample_files <- c(
  hepatocyte_1 = file.path(input_dir, "hepatocyte_1.bed"),
  hepatocyte_2 = file.path(input_dir, "hepatocyte_2.bed"),
  HepG2_1 = file.path(input_dir, "HepG2_1.bed"),
  HepG2_2 = file.path(input_dir, "HepG2_2.bed")
)
condition_levels <- c("hepatocyte", "HepG2")
fdr_cutoff <- 0.01
log2fc_cutoff <- 1

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
project_arg <- get_arg("--project-dir")
if (!is.null(project_arg)) {
  project_dir <- project_arg
  input_dir <- file.path(project_dir, "input")
  result_dir <- file.path(project_dir, "result", "DEGs")
  sample_files <- c(
    hepatocyte_1 = file.path(input_dir, "hepatocyte_1.bed"),
    hepatocyte_2 = file.path(input_dir, "hepatocyte_2.bed"),
    HepG2_1 = file.path(input_dir, "HepG2_1.bed"),
    HepG2_2 = file.path(input_dir, "HepG2_2.bed")
  )
}

check_file <- function(path, label = "Input file") {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}
read_sample <- function(path, sample_name) {
  check_file(path, paste0("Sample ", sample_name))
  x <- read.table(path, header = TRUE, sep = "", quote = "",
                  comment.char = "", check.names = FALSE,
                  stringsAsFactors = FALSE)
  if (ncol(x) < 2L) stop("Sample file has fewer than 2 columns: ", path, call. = FALSE)
  gene_col <- match("geneName", names(x)); if (is.na(gene_col)) gene_col <- 1L
  value_cols <- setdiff(seq_len(ncol(x)), gene_col)
  if (length(value_cols) != 1L)
    stop("Each sample file must have exactly one value column besides geneName: ", path, call. = FALSE)
  gene <- as.character(x[[gene_col]])
  value <- suppressWarnings(as.numeric(x[[value_cols]]))
  if (anyNA(gene) || any(!nzchar(gene)) || anyDuplicated(gene))
    stop("Gene IDs are missing or duplicated in: ", path, call. = FALSE)
  if (anyNA(value) || any(!is.finite(value)) || any(value < 0))
    stop("Counts must be finite and non-negative in: ", path, call. = FALSE)
  data.frame(geneName = gene, value = value, check.names = FALSE)
}

sample_data <- Map(read_sample, unname(sample_files), names(sample_files))
common_genes <- Reduce(intersect, lapply(sample_data, function(x) x$geneName))
if (!length(common_genes)) stop("No common genes found across all four samples.", call. = FALSE)
mycounts <- do.call(cbind, lapply(sample_data, function(x) x$value[match(common_genes, x$geneName)]))
rownames(mycounts) <- common_genes; colnames(mycounts) <- names(sample_files)
if (any(abs(mycounts - round(mycounts)) > 1e-8))
  warning("Non-integer values detected; rounding for DESeq2. Use raw counts when possible.")
mycounts <- round(mycounts); storage.mode(mycounts) <- "integer"
condition <- factor(c("hepatocyte", "hepatocyte", "HepG2", "HepG2"), levels = condition_levels)
col_data <- data.frame(row.names = colnames(mycounts), condition = condition)
keep <- rowSums(mycounts) > 0L; mycounts <- mycounts[keep, , drop = FALSE]
if (!nrow(mycounts)) stop("No genes remain after removing all-zero rows.", call. = FALSE)

dds <- DESeqDataSetFromMatrix(countData = mycounts, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "HepG2", "hepatocyte"))
res <- res[order(res$padj, res$pvalue, na.last = TRUE), ]

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
write.table(as.data.frame(res), file.path(result_dir, "All_results.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(as.data.frame(subset(res, !is.na(padj) & padj < fdr_cutoff & abs(log2FoldChange) > log2fc_cutoff)), file.path(result_dir, "All_DEGs.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(as.data.frame(subset(res, !is.na(padj) & padj < fdr_cutoff & log2FoldChange > log2fc_cutoff)), file.path(result_dir, "Up_DEGs.txt"), sep = "\t", quote = FALSE, col.names = NA)
write.table(as.data.frame(subset(res, !is.na(padj) & padj < fdr_cutoff & log2FoldChange < -log2fc_cutoff)), file.path(result_dir, "Down_DEGs.txt"), sep = "\t", quote = FALSE, col.names = NA)
message("Analysis complete: HepG2 vs hepatocyte.")
message("Results: ", normalizePath(result_dir, winslash = "/", mustWork = FALSE))
