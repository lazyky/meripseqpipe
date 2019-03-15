#!/bin/Rscript
## Rscript edgeR.R aligner_tools designfile gtf eg. Rscript edgeR.R tophat2 designfile_single.txt
### designfile: filename, control_or_treated, input_or_ip, situation(default 1 is CONTROL_SITUATION else are TREATED_SITUATION)
### TREATED_SITUATION_STARTPOINT is the default situation check point
library(stringr)
library("edgeR")
args<-commandArgs(T) 
TREATED_SITUATION_STARTPOINT <- args[1] 
aligner_tools_name <- args[2]
designfile <- args[3]
TREATED_SITUATION_STARTPOINT <- as.numeric(TREATED_SITUATION_STARTPOINT)

# setting CONTROL_SITUATION and TREATED_SITUATION 
## default 1 is CONTROL_SITUATION else are TREATED_SITUATION
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
CONTROL_SITUATION <- c()
TREATED_SITUATION <- c()
if (TREATED_SITUATION_STARTPOINT-1 >= 1){
  for (i in c(1:(TREATED_SITUATION_STARTPOINT-1))){
  CONTROL_SITUATION <- c(CONTROL_SITUATION,str_c("_",i,"_"))
  }
}
for (i in c(TREATED_SITUATION_STARTPOINT:max(as.numeric(designtable$situation)))){
  TREATED_SITUATION <- c(TREATED_SITUATION,str_c("_",i,"_"))
}
#合并表达矩阵
## Running edgeR and rename the output name
for(i in CONTROL_SITUATION){
  for(j in TREATED_SITUATION){
    control_database=read.table(str_c("htseq_situation", i, aligner_tools_name, "_input.count"), header = TRUE, row.names = 1)
    treated_database=read.table(str_c("htseq_situation", j, aligner_tools_name, "_input.count"), header = TRUE, row.names = 1)
    combined_database <- cbind(control_database,treated_database)
    group <- factor(c(rep("control",ncol(control_database)), rep("treated",ncol(treated_database)))) #设定因子
    y <- DGEList(counts=combined_database,group=group)
    rownames(y) <- rownames(combined_database)
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
    ### set output_name
    output_name <- str_c("edgeR_situation",i ,j ,aligner_tools_name,".csv")
    write.csv(combined_database, file = output_name)
  }
}


