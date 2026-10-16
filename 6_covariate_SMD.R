#!/usr/bin/env Rscript
# Input: tab-separated table with headers:
# gene_id, group, gene_length, exon_length, expression, tss_gc,
# H3K27ac, H3K4me1, DNase
# group: differential / non_differential
# All numeric columns must be on the intended analysis scale BEFORE input.
# No log/standardization is silently applied. GC may be a fraction or percent,
# but must use one consistent unit. One row per gene.
project_dir <- "path/to/project"
args <- commandArgs(TRUE)
i <- match("--project-dir", args)
if (!is.na(i)) {
  if (i == length(args)) stop("--project-dir requires a value.")
  project_dir <- args[i + 1L]
}
input_file <- file.path(project_dir, "input", "chromatin_gene_table.tsv")
output_dir <- file.path(project_dir, "results", "chromatin_adjustment")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
dat <- read.delim(input_file, check.names=FALSE, colClasses="character")
covariates <- c("gene_length","exon_length","expression","tss_gc")
outcomes <- c("H3K27ac","H3K4me1","DNase")
required <- c("gene_id","group",covariates,outcomes)
if (!all(required %in% names(dat))) stop("Missing columns: ",paste(setdiff(required,names(dat)),collapse=", "))
if (anyNA(dat$gene_id) || any(!nzchar(dat$gene_id)) || anyDuplicated(dat$gene_id))
  stop("gene_id must be unique and nonempty.")
if (any(!is.na(dat$group) & !dat$group %in% c("non_differential","differential")))
  stop("group must be differential or non_differential.")
for (v in c(covariates,outcomes)) {
  original <- dat[[v]]
  dat[[v]] <- suppressWarnings(as.numeric(original))
  if (any(!is.na(original) & nzchar(original) & is.na(dat[[v]])))
    stop("Non-numeric values in ",v)
}
# A single complete-case cohort for S3, all outcomes and every model stage.
keep <- !is.na(dat$group) & Reduce(`&`,lapply(dat[c(covariates,outcomes)],is.finite))
write.table(data.frame(gene_id=dat$gene_id,included=keep),
            file.path(output_dir,"complete_case_audit.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
dat <- dat[keep,,drop=FALSE]
dat$group <- factor(dat$group,levels=c("non_differential","differential"))
if (any(table(dat$group)<2L)) stop("Each group needs at least two complete-case genes.")
write_tsv <- function(x,name) write.table(x,file.path(output_dir,name),
                                        sep="\t",quote=FALSE,row.names=FALSE,na="NA")

# Unweighted pooled within-group SD: sqrt((SD_diff^2 + SD_non_diff^2)/2).
# SMD is differential minus non-differential, computed without outcome fitting.
results <- do.call(rbind,lapply(covariates,function(v) {
  a <- dat[dat$group=="differential",v]
  b <- dat[dat$group=="non_differential",v]
  pooled <- sqrt((var(a)+var(b))/2)
  data.frame(covariate=v,n_differential=length(a),n_non_differential=length(b),
             mean_differential=mean(a),mean_non_differential=mean(b),
             sd_differential=sd(a),sd_non_differential=sd(b),
             SMD=if(pooled>0) (mean(a)-mean(b))/pooled else NA_real_)
}))
write_tsv(results,"Supplementary_Table_S3.tsv")
capture.output(sessionInfo(),file=file.path(output_dir,"S3_sessionInfo.txt"))
message("Saved S3: ",output_dir)
