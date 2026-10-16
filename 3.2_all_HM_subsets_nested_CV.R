#!/usr/bin/env Rscript
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1L], mustWork = FALSE))
source(file.path(script_dir, "hm_nested_cv_common.R"))
load_hm_packages()

project_dir <- "path/to/project"
signal_dir <- file.path(project_dir, "signals"); deg_dir <- file.path(project_dir, "DEGs")
output_dir <- file.path(project_dir, "results", "3.2_all_HM_subsets")
hms <- c("H2AFZ","H3K4me1","H3K4me2","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","H3K79me2","H4K20me1")
methods <- c("RF","SVM","GLM"); outer_k <- 10L; inner_k <- 5L; seed <- 1L
max_hm_subset_size <- length(hms)
args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) { i <- match(flag,args); if (is.na(i)||i==length(args)) default else args[[i+1L]] }
if (!is.null(z <- get_arg("--project-dir"))) { project_dir<-z; signal_dir<-file.path(z,"signals"); deg_dir<-file.path(z,"DEGs"); output_dir<-file.path(z,"results","3.2_all_HM_subsets") }
if (!is.null(z <- get_arg("--max-hm-subset"))) max_hm_subset_size <- as.integer(z)
if (!is.null(z <- get_arg("--signal-dir"))) signal_dir <- z
if (!is.null(z <- get_arg("--deg-dir"))) deg_dir <- z
if (!is.null(z <- get_arg("--output-dir"))) output_dir <- z
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

labels <- make_labels(file.path(deg_dir,"Up_DEGs.txt"), file.path(deg_dir,"Down_DEGs.txt"))
changes <- make_hm_changes(signal_dir, hms)
merged <- labels
for (hm in hms) merged <- merge(merged, changes[[hm]], by="gene_id")
results <- list()
cols <- unlist(lapply(hms, function(h) grep(paste0("^",h,"_"), names(merged), value=TRUE)))
x <- as.matrix(merged[,cols,drop=FALSE]); y <- merged$y
for (method in methods) {
  # All HM subsets are candidates inside inner CV; the outer test fold is untouched.
  fit <- nested_oof(x,y,method,selector="subset",outer_k=outer_k,inner_k=inner_k,seed=seed,
                    max_subset_size=max_hm_subset_size)
  selected <- vapply(fit$selected, paste, collapse=";", FUN.VALUE=character(1))
  results[[length(results)+1L]] <- data.frame(HM_subset_selected=paste(selected,collapse=" | "), method=method, fit$metrics)
  write.table(data.frame(gene_id=merged$gene_id,label=y,prediction=fit$prediction), file.path(output_dir,paste0(method,"_OOF_predictions.txt")), sep="\t", quote=FALSE, row.names=FALSE)
}
write.table(do.call(rbind,results),file.path(output_dir,"nested_CV_metrics.txt"),sep="\t",quote=FALSE,row.names=FALSE)

