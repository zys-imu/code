#!/usr/bin/env Rscript

# Integrated promoter-window generation and DEGseq mapping pipeline.
# Directory layout expected under project_dir:
#   raw_data/RefSeq_gene/RefSeq-noNR_noRep.txt
#   raw_data/bed_format/<sample>.bed
# Outputs:
#   raw_data/RefSeq_gene/TSS_up_*.txt and TSS_down_*.txt
#   pretreated_data/HepG2/H3K4me1/<sample>_<direction>_<window>.bed

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  if (!("--merge-only" %in% commandArgs(trailingOnly = TRUE)) &&
      !requireNamespace("DEGseq", quietly = TRUE))
    stop("R package 'DEGseq' is required.", call. = FALSE)
  # Use DEGseq::getGeneExp explicitly; merge-only needs no DEGseq installation.
})

# ------------------------- user configuration -----------------------------
# Use a generic project path. Change this line or pass --project-dir.
project_dir <- "path/to/project"
samples <- c("HepG2_H3K4me1_1", "HepG2_H3K4me1_2")
windows <- 1:50
directions <- c("up", "down")
tss_offset <- 5000
bin_width <- 100
min_overlap_percent <- 0.5
skip_existing <- FALSE
overwrite_tss <- TRUE
merge_only <- FALSE
merge_value_col <- 2L # DEGseq column 2: raw count; column 1: gene ID

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
project_arg <- get_arg("--project-dir")
if (!is.null(project_arg)) project_dir <- project_arg
if ("--skip-existing" %in% args) skip_existing <- TRUE
if ("--keep-tss" %in% args) overwrite_tss <- FALSE
if ("--merge-only" %in% args) merge_only <- TRUE

ref_dir <- file.path(project_dir, "raw_data", "RefSeq_gene")
bed_dir <- file.path(project_dir, "raw_data", "bed_format")
output_dir <- file.path(project_dir, "pretreated_data", "HepG2", "H3K4me1")
refseq_file <- file.path(ref_dir, "RefSeq-noNR_noRep.txt")

check_file <- function(path, label) {
  if (!file.exists(path))
    stop(label, " not found: ", path, call. = FALSE)
}

make_window <- function(tss, bin) {
  start0 <- tss - tss_offset
  start <- if (bin == 1L) start0 else start0 + bin_width * (bin - 1L) + 1
  end <- start0 + bin_width * bin
  list(start = start, end = end)
}

read_refseq <- function(path) {
  check_file(path, "RefSeq input")
  x <- read.table(path, header = FALSE, sep = "", quote = "",
                  comment.char = "", fill = TRUE, check.names = FALSE)
  if (ncol(x) < 5L) stop("RefSeq input needs at least 5 columns.", call. = FALSE)
  tss <- suppressWarnings(as.numeric(x[[5L]]))
  if (anyNA(tss)) stop("The 5th RefSeq column must contain numeric TSS coordinates.",
                       call. = FALSE)
  x
}

write_tss_files <- function(refseq) {
  dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", 100L)
  for (bin in 1:100) {
    w <- make_window(as.numeric(refseq[[5L]]), bin)
    rows[[bin]] <- paste(
      refseq[[1L]], refseq[[2L]], refseq[[3L]], refseq[[4L]],
      w[["start"]], w[["end"]], w[["start"]], w[["end"]], 1,
      paste0(w[["start"]], ","), paste0(w[["end"]], ","), sep = "\t")
  }
  names(rows) <- c(paste0("TSS_up_", 50:1), paste0("TSS_down_", 1:50))
  for (i in seq_along(rows)) {
    path <- file.path(ref_dir, paste0(names(rows)[i], ".txt"))
    if (overwrite_tss || !file.exists(path))
      writeLines(rows[[i]], path, useBytes = TRUE)
    else message("KEEP  ", basename(path))
  }
}

run_mapping <- function(sample, direction, window) {
  bed_file <- file.path(bed_dir, paste0(sample, ".bed"))
  ref_file <- file.path(ref_dir, sprintf("TSS_%s_%d.txt", direction, window))
  output_file <- file.path(output_dir, sprintf("%s_%s_%d.bed",
                                               sample, direction, window))
  check_file(bed_file, "BED input")
  check_file(ref_file, "TSS reference")
  if (skip_existing && file.exists(output_file) &&
      isTRUE(file.info(output_file)$size > 0)) {
    message("SKIP  ", basename(output_file))
    return(invisible(NULL))
  }
  tryCatch(
    DEGseq::getGeneExp(list(bed_file), refFlat = ref_file,
                       output = output_file,
                       min.overlapPercent = min_overlap_percent),
    error = function(e) stop(sprintf("Mapping failed for %s/%s/%d: %s",
                                     sample, direction, window,
                                     conditionMessage(e)), call. = FALSE)
  )
  message("DONE  ", output_file)
}

merge_sample <- function(sample) {
  labels <- c(paste0("up_", 50:1), paste0("down_", 1:50))
  paths <- file.path(output_dir, paste0(sample, "_", labels, ".bed"))
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files))
    stop("Missing window files:\n", paste(missing_files, collapse = "\n"),
         call. = FALSE)

  ids <- NULL
  counts <- NULL
  for (i in seq_along(paths)) {
    tab <- read.table(paths[i], header = TRUE, sep = "\t", quote = "\"",
                      comment.char = "", check.names = FALSE,
                      colClasses = "character", row.names = NULL)
    if (!nrow(tab) || ncol(tab) < merge_value_col)
      stop("Empty or invalid mapping table: ", paths[i], call. = FALSE)
    gene <- tab[[1L]]
    if (anyNA(gene) || any(!nzchar(gene)) || anyDuplicated(gene))
      stop("Missing or duplicate gene IDs in: ", paths[i], call. = FALSE)
    values <- suppressWarnings(as.numeric(tab[[merge_value_col]]))
    if (any(!is.finite(values)) || any(values < 0))
      stop("Invalid count values in: ", paths[i], call. = FALSE)
    if (is.null(ids)) {
      ids <- gene
      counts <- matrix(NA_real_, nrow = length(ids), ncol = length(labels),
                       dimnames = list(NULL, labels))
    }
    if (length(gene) != length(ids) || !setequal(gene, ids))
      stop("Gene set differs from up_50 in: ", paths[i], call. = FALSE)
    counts[, i] <- values[match(ids, gene)]
  }
  merged <- data.frame(gene_id = ids, counts, check.names = FALSE)
  # Extensionless tab-separated filenames match the requested sample names.
  destination <- file.path(output_dir, sample)
  write.table(merged, file = destination, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  message("MERGED  ", destination, " (", nrow(merged), " genes x 100 windows)")
  invisible(destination)
}

if (!merge_only) {
check_file(refseq_file, "RefSeq input")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
refseq <- read_refseq(refseq_file)
write_tss_files(refseq)

jobs <- expand.grid(sample = samples, direction = directions, window = windows,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
message("Starting ", nrow(jobs), " mapping jobs...")
for (i in seq_len(nrow(jobs))) {
  j <- jobs[i, , drop = FALSE]
  run_mapping(j$sample, j$direction, as.integer(j$window))
}
message("Finished ", nrow(jobs), " mapping jobs.")
}

for (sample in samples) merge_sample(sample)

