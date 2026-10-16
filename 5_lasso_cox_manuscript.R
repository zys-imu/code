#!/usr/bin/env Rscript
# Manuscript section 1.5: LASSO Cox -> multivariable Cox -> risk score.
script_file<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(script_file)),"survival_pipeline_common.R"))
for(pkg in c("survival","glmnet")) if(!requireNamespace(pkg,quietly=TRUE)) stop("Install ",pkg)
project_dir<-project_arg()
output_dir<-file.path(project_dir,"results","lasso_cox")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
data<-readRDS(file.path(project_dir,"results","univariate_cox","training_data.rds"))
genes<-data$seeds
if(length(genes)<2) stop("LASSO requires at least two eligible seed genes.")
x<-t(data$expr[genes,,drop=FALSE]);cl<-data$clinical
if(!identical(rownames(x),cl$Sample)) stop("Patient order mismatch.")
# These choices are configurable because the manuscript does not specify them.
seed<-1006L
cv_folds<-10L
lambda_choice<-"lambda.min"
set.seed(seed)
if(sum(cl$OS)<cv_folds || nrow(cl)<2*cv_folds)
  stop("Insufficient patients/events for configured CV folds; reduce cv_folds explicitly.")
foldid<-integer(nrow(cl))
for(event in c(0,1)) {
  idx<-which(cl$OS==event)
  foldid[idx]<-sample(rep(seq_len(cv_folds),length.out=length(idx)))
}
y<-survival::Surv(cl$OS.time,cl$OS)
cv_fit<-glmnet::cv.glmnet(x=x,y=y,family="cox",alpha=1,nlambda=100,
                         lambda.min.ratio=0.2,foldid=foldid,nfolds=cv_folds)
beta<-as.matrix(stats::coef(cv_fit,s=lambda_choice))
selected<-rownames(beta)[beta[,1]!=0]
write_tsv(data.frame(gene=rownames(beta),lasso_coefficient=beta[,1]),
          file.path(output_dir,"lasso_all_coefficients.tsv"))
write_tsv(data.frame(Sample=cl$Sample,fold=foldid),file.path(output_dir,"lambda_cv_folds.tsv"))
saveRDS(cv_fit,file.path(output_dir,"lasso_cv_model.rds"))
pdf(file.path(output_dir,"lasso_cross_validation.pdf"));plot(cv_fit);dev.off()
pdf(file.path(output_dir,"lasso_coefficient_paths.pdf"));plot(cv_fit$glmnet.fit,xvar="lambda");dev.off()
write_tsv(data.frame(gene=selected),file.path(output_dir,"lasso_selected_genes.tsv"))
if(!length(selected)) stop("Chosen lambda selects zero genes; no risk signature fitted.")
# Internal column names avoid formula failures for gene symbols with punctuation.
xx<-x[,selected,drop=FALSE]; colnames(xx)<-paste0("g",seq_along(selected))
d<-data.frame(OS.time=cl$OS.time,OS=cl$OS,xx,check.names=FALSE)
fit<-survival::coxph(survival::Surv(OS.time,OS)~.,data=d,x=TRUE,y=TRUE)
b<-stats::coef(fit)
if(any(!is.finite(b))) stop("Multivariable Cox coefficients are not identifiable.")
s<-summary(fit)
coef_table<-data.frame(gene=selected,coefficient=as.numeric(b),
  HR=s$conf.int[,"exp(coef)"],CI_lower=s$conf.int[,"lower .95"],
  CI_upper=s$conf.int[,"upper .95"],p_value=s$coefficients[,"Pr(>|z|)"])
write_tsv(coef_table,file.path(output_dir,"multivariable_cox_coefficients.tsv"))
# Uncentered linear predictor: sum(expression * multivariable Cox coefficient).
risk<-as.numeric(x[,selected,drop=FALSE] %*% b)
cutoff<-median(risk)
write_tsv(data.frame(Sample=cl$Sample,OS.time=cl$OS.time,OS=cl$OS,
          risk_score=risk,risk_group=ifelse(risk>cutoff,"High","Low")),
          file.path(output_dir,"TCGA_risk_scores.tsv"))
write_tsv(data.frame(lambda_choice=lambda_choice,lambda=cv_fit[[lambda_choice]],
          alpha=1,nlambda=100,lambda_min_ratio=0.2,cv_folds=cv_folds,
          seed=seed,training_median_cutoff=cutoff,n_selected=length(selected)),
          file.path(output_dir,"model_settings.tsv"))
saveRDS(list(cox=fit,genes=selected,coefficients=b,cutoff=cutoff,lambda_choice=lambda_choice),
        file.path(output_dir,"prognostic_signature.rds"))
capture.output(sessionInfo(),file=file.path(output_dir,"sessionInfo.txt"))
message("Selected ",length(selected)," genes. Risk scores are development-cohort scores, not external validation.")
