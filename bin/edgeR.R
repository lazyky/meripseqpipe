#!/bin/Rscript
## Rscript edgeR.R aligner_tools designfile gtf eg. Rscript edgeR.R tophat2 designfile_single.txt
### designfile: filename, control_or_treated, input_or_ip, group(default 1 is CONTROL_SITUATION else are TREATED_SITUATION)
### TREATED_SITUATION_STARTPOINT is the default group check point
library("edgeR")
args<-commandArgs(T) 
designfile <- args[1]

designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
# Running edgeR by compare.file
## while there are only 2 groups, running edgeR without compare.file
if(length(unique(designtable$Group)) == 2){
  # Combine expression matrix
  ## Running edgeR and rename the output name
  group_id_1 <- unique(designtable$Group)[1]
  group_id_2 <- unique(designtable$Group)[2]
  control_database = read.table(paste0("htseq_group_", group_id_1, "_input.count"), header = TRUE, row.names = 1)
  treated_database = read.table(paste0("htseq_group_", group_id_2, "_input.count"), header = TRUE, row.names = 1)
  combined_database <- cbind(control_database,treated_database)
  group <- factor(c(rep(group_id_1,ncol(control_database)), rep(group_id_2,ncol(treated_database)))) #setting factors
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
  output_name <- paste0("edgeR_group_",group_id_1, "_",group_id_2)
  write.csv(combined_database, file = paste0(output_name,".matirx") )
  write.csv(qlf$table, file = paste0(output_name, "_qlf.csv"))
  write.csv(lrt$table, file = paste0(output_name, "_lrt.csv"))
}else if(length(unique(designtable$Group)) < 2){
  stop( "The count of Group is less than two, please check your designfile.")
}else{
  print("multi-groups compare will coming soon")
}