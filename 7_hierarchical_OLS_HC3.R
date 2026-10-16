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

if(!requireNamespace("sandwich",quietly=TRUE)) stop("Install package sandwich.")
# Hierarchy is an explicit configurable choice:
# M0: group; M1: + gene/exon lengths; M2: + expression; M3: + GC.
# All stages and outcomes use exactly the complete-case cohort defined above.
stages <- list(M0=character(), M1=c("gene_length","exon_length"),
               M2=c("gene_length","exon_length","expression"),
               M3=covariates)
results <- list()
for(stage in names(stages)) {
  for(outcome in outcomes) {
    model <- lm(reformulate(c("group",stages[[stage]]),response=outcome),data=dat)
    if(model$rank<ncol(model.matrix(model))) stop("Rank-deficient model: ",stage,"/",outcome)
    if(df.residual(model)<=0 || any(hatvalues(model)>=1-1e-10))
      stop("Insufficient residual degrees of freedom or unit leverage.")
    vc <- sandwich::vcovHC(model,type="HC3")
    term <- "groupdifferential"
    beta <- unname(coef(model)[term]);se <- sqrt(vc[term,term])
    if(!is.finite(beta) || !is.finite(se) || se<=0) stop("Invalid HC3 estimate.")
    # Two-sided asymptotic Wald z inference, consistently used for CI and p.
    p <- 2*pnorm(-abs(beta/se))
    results[[length(results)+1L]] <- data.frame(model=stage,outcome=outcome,n=nrow(dat),
      beta=beta,HC3_SE=se,CI_lower=beta-qnorm(.975)*se,
      CI_upper=beta+qnorm(.975)*se,p_value=p,formula=paste(deparse(formula(model)),collapse=" "))
  }
}
results <- do.call(rbind,results)
# BH correction across the THREE outcomes separately within EACH model.
results$BH_FDR <- NA_real_
for(stage in names(stages)) {
  k <- results$model==stage
  results$BH_FDR[k] <- p.adjust(results$p_value[k],method="BH")
}
write_tsv(results,"Supplementary_Table_S4.tsv")
attenuation <- do.call(rbind,lapply(outcomes,function(outcome) {
  b0 <- results$beta[results$model=="M0" & results$outcome==outcome]
  b3 <- results$beta[results$model=="M3" & results$outcome==outcome]
  data.frame(outcome=outcome,beta_unadjusted=b0,beta_fully_adjusted=b3,
             attenuation_percent=if(abs(b0)>.Machine$double.eps) 100*(1-b3/b0) else NA_real_)
}))
# Signed attenuation: a sign reversal may give >100%; report with coefficients.
write_tsv(attenuation,"coefficient_attenuation.tsv")
capture.output(sessionInfo(),file=file.path(output_dir,"S4_sessionInfo.txt"))
message("Saved S4 and attenuation: ",output_dir)
