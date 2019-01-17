#!/bin/Rscript
library(stringr)
library("edgeR")
setwd(dir = "E:/zky/m6Apipe/bin/readscount/") #test
args<-commandArgs(T) # Rscript get_htseq_matrix.R aligner_tools eg. bwa
output_name <- str_c("treated_vs_control_","bwa","_edgeR.csv")

#合并表达矩阵
control_file = str_c("control_input_","bwa",".count") 
treated_file = str_c("treated_input_","bwa",".count")
control_database=read.table(control_file,sep="\t",header=T,row.names=1)
treated_database=read.table(treated_file,sep="\t",header=T,row.names=1)
combined_database <- cbind(control_database,treated_database)
group <- factor(c(rep("control",ncol(control_database)), rep("treated",ncol(treated_database)))) #设定因子

#edgeR
y <- DGEList(counts=combined_database,group=group)
#rownames(y) <- rownames(combined_database)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#To perform quasi-likelihood F-tests:
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)
#To perform likelihood ratio tests:
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
write.csv(combined_database, file = output_name ) #test
