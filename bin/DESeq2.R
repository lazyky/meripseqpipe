#!/bin/Rscript
## Rscript DESeq2.R aligner_tools designfile gtf eg. Rscript DESeq2.R tophat2 designfile_single.txt
### designfile: filename, control_or_treated, input_or_ip, group(default 1 is CONTROL_SITUATION else are TREATED_SITUATION)
### TREATED_SITUATION_STARTPOINT is the default group check point
library(DESeq2)
args<-commandArgs(T) 
designfile <- args[1]
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
# Running DEseq2 by compare.file
## while there are only 2 groups, running DEseq2 without compare.file
if(length(unique(designtable$Group)) == 2){
  # Combine expression matrix
  ## Running DEseq2 and rename the output name
  group_id_1 <- unique(designtable$Group)[1]
  group_id_2 <- unique(designtable$Group)[2]
  control_database = read.table(paste0("htseq_group_", group_id_1, "_input.count"), header = TRUE, row.names = 1)
  treated_database = read.table(paste0("htseq_group_", group_id_2, "_input.count"), header = TRUE, row.names = 1)
  combined_database <- cbind(control_database,treated_database)
  condition <- factor(c(rep(group_id_1,ncol(control_database)), rep(group_id_2,ncol(treated_database)))) #setting factors
  ### assign gene names
  colData <- data.frame(row.names=colnames(combined_database), condition)
  dds <- DESeqDataSetFromMatrix(countData = combined_database,colData = colData,design = ~ condition)
  rownames(dds) <- rownames(combined_database)
  #dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  res <- results(dds)
  table(res$padj <0.05)
  res <- res[order(res$padj),]
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
  resdata2=resdata[resdata$log2FoldChange > 1|resdata$log2FoldChange < -1, ]
  ### set output_name
  output_name <- paste0("Deseq2_group_",group_id_1, "_",group_id_2)
  write.csv(resdata, file = paste0(output_name, ".csv"))
  write.csv(resdata2,file = paste0(output_name, "_log2.csv"),row.names =FALSE)
}else if(length(unique(designtable$Group)) < 2){
  stop( "The count of Group is less than two, please check your designfile.")
}else{
  print("multi-groups compare will coming soon")
}
