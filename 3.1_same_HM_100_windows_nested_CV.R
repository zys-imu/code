#!/usr/bin/env Rscript
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1L], mustWork = FALSE))
source(file.path(script_dir, "hm_nested_cv_common.R"))
load_hm_packages()

project_dir <- "path/to/project"
signal_dir <- file.path(project_dir, "signals")
deg_dir <- file.path(project_dir, "DEGs")
output_dir <- file.path(project_dir, "results", "3.1_same_HM_100_windows")
hms <- c("H2AFZ","H3K4me1","H3K4me2","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","H3K79me2","H4K20me1")
methods <- c("RF","SVM","GLM")
outer_k <- 10L; inner_k <- 5L; seed <- 1L; max_features <- 10L
args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) { i <- match(flag,args); if (is.na(i)||i==length(args)) default else args[[i+1L]] }
if (!is.null(z <- get_arg("--project-dir"))) { project_dir <- z; signal_dir <- file.path(z,"signals"); deg_dir <- file.path(z,"DEGs"); output_dir <- file.path(z,"results","3.1_same_HM_100_windows") }
if (!is.null(z <- get_arg("--signal-dir"))) signal_dir <- z
if (!is.null(z <- get_arg("--deg-dir"))) deg_dir <- z
if (!is.null(z <- get_arg("--output-dir"))) output_dir <- z
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

labels <- make_labels(file.path(deg_dir,"Up_DEGs.txt"), file.path(deg_dir,"Down_DEGs.txt"))
changes <- make_hm_changes(signal_dir, hms)
all_results <- list()
for (hm in hms) {
  z <- changes[[hm]]
  dat <- merge(labels, z, by="gene_id")
  x <- as.matrix(dat[, setdiff(names(dat),c("gene_id","y")), drop=FALSE]); y <- dat$y
  for (method in methods) {
    fit <- nested_oof(x, y, method, selector="greedy", outer_k=outer_k, inner_k=inner_k, seed=seed, max_features=max_features)
    all_results[[length(all_results)+1L]] <- data.frame(HM=hm, method=method, fit$metrics)
    write.table(data.frame(gene_id=dat$gene_id, label=y, prediction=fit$prediction), file.path(output_dir,paste0(hm,"_",method,"_OOF_predictions.txt")), sep="\t", quote=FALSE, row.names=FALSE)
  }
}
write.table(do.call(rbind,all_results), file.path(output_dir,"nested_CV_metrics.txt"), sep="\t", quote=FALSE, row.names=FALSE)

