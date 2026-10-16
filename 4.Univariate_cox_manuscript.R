#!/usr/bin/env Rscript
# Manuscript section 1.5: zero-expression filter then continuous univariate Cox.
script_file<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(script_file)),"survival_pipeline_common.R"))
if(!requireNamespace("survival",quietly=TRUE)) stop("Install survival.")
project_dir<-project_arg()
input_dir<-file.path(project_dir,"input")
output_dir<-file.path(project_dir,"results","univariate_cox")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
data<-read_cohort(file.path(input_dir,"TCGA_expression.tsv"),file.path(input_dir,"TCGA_clinical.tsv"))
# Headerless candidate list, first column = differential H3K79me2 signal gene ID.
candidate<-unique(read.delim(file.path(input_dir,"H3K79me2_candidate_genes.txt"),
                           header=FALSE,colClasses="character")[[1]])
present<-intersect(candidate,rownames(data$expr))
write_tsv(data.frame(gene=setdiff(candidate,present)),
          file.path(output_dir,"missing_candidate_genes.tsv"))
x<-data$expr[present,,drop=FALSE]
zero_fraction<-rowMeans(x==0)
eligible<-zero_fraction<=0.5
write_tsv(data.frame(gene=present,zero_fraction=zero_fraction,retained=eligible),
          file.path(output_dir,"expression_filter.tsv"))
x<-x[eligible,,drop=FALSE]
if(!nrow(x)) stop("No eligible candidate genes.")
results<-lapply(rownames(x),function(g) {
  d<-data$clinical; d$expression<-as.numeric(x[g,])
  ans<-data.frame(gene=g,n=nrow(d),events=sum(d$OS),HR=NA_real_,
                  CI_lower=NA_real_,CI_upper=NA_real_,p_value=NA_real_,status="ok")
  if(stats::sd(d$expression)==0) { ans$status<-"constant_expression"; return(ans) }
  warnings<-character()
  fit<-tryCatch(withCallingHandlers(
    survival::coxph(survival::Surv(OS.time,OS)~expression,data=d),
    warning=function(w) {warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e)e)
  if(inherits(fit,"error")) {ans$status<-conditionMessage(fit);return(ans)}
  s<-summary(fit)
  ans$HR<-s$conf.int[1,"exp(coef)"];ans$CI_lower<-s$conf.int[1,"lower .95"]
  ans$CI_upper<-s$conf.int[1,"upper .95"];ans$p_value<-s$coefficients[1,"Pr(>|z|)"]
  if(length(warnings)) ans$status<-paste(warnings,collapse="; ")
  if(any(!is.finite(c(ans$HR,ans$CI_lower,ans$CI_upper,ans$p_value))))
    ans$status<-"nonfinite_estimate"
  ans
})
results<-do.call(rbind,results)
results$BH_FDR<-p.adjust(results$p_value,method="BH")
results<-results[order(results$p_value,na.last=TRUE),]
# Manuscript specifies raw p < 0.05, not FDR/log-rank/5-year curve filtering.
seeds<-results[results$status=="ok" & !is.na(results$p_value) & results$p_value<0.05,,drop=FALSE]
write_tsv(results,file.path(output_dir,"all_univariate_cox.tsv"))
write_tsv(seeds,file.path(output_dir,"seed_genes.tsv"))
saveRDS(list(expr=x,clinical=data$clinical,seeds=seeds$gene),
        file.path(output_dir,"training_data.rds"))
capture.output(sessionInfo(),file=file.path(output_dir,"sessionInfo.txt"))
message("Seed genes: ",nrow(seeds),"; output: ",output_dir)
