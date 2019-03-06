#!/bin/Rscript
## Rscript DESeq2.R aligner_tools designfile gtf eg. Rscript DESeq2.R tophat2 designfile_single.txt
### designfile: filename, control_or_treated, input_or_ip, situation(default 0 is CONTROL_SITUATION else are TREATED_SITUATION)
### TREATED_SITUATION_STARTPOINT is the default situation check point
library(stringr)
library(DESeq2)
args<-commandArgs(T) 
TREATED_SITUATION_STARTPOINT <- args[1] 
aligner_tools_name <- args[2]
designfile <- args[3]
TREATED_SITUATION_STARTPOINT <- as.numeric(TREATED_SITUATION_STARTPOINT)

# setting CONTROL_SITUATION and TREATED_SITUATION 
## default 0 is CONTROL_SITUATION else are TREATED_SITUATION
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
CONTROL_SITUATION <- c()
TREATED_SITUATION <- c()
for (i in c(0:(TREATED_SITUATION_STARTPOINT-1))){
  CONTROL_SITUATION <- c(CONTROL_SITUATION,str_c("_",i,"_"))
}
for (i in c(TREATED_SITUATION_STARTPOINT:max(designtable$situation))){
  TREATED_SITUATION <- c(TREATED_SITUATION,str_c("_",i,"_"))
}

## Running DEseq2 and rename the output name
for(i in CONTROL_SITUATION){
  for(j in TREATED_SITUATION){
    control_database=read.table(str_c("htseq_situation", i, aligner_tools_name, "_input.count"), header = TRUE, row.names = 1)
    treated_database=read.table(str_c("htseq_situation", j, aligner_tools_name, "_input.count"), header = TRUE, row.names = 1)
    combined_database <- cbind(control_database,treated_database)
    condition <- factor(c(rep("control",ncol(control_database)), rep("treated",ncol(treated_database)))) #设定因子
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
    resdata2=resdata[resdata$log2FoldChange > 1|resdata$log2FoldChange < -1,]
    ### set output_name
    output_name <- str_c("Deseq2_situation",i ,j ,aligner_tools_name,".csv")
    write.csv(resdata, file = output_name)
    write.csv(resdata2,file = str_c("Deseq2_situation",i ,j ,aligner_tools_name,"_log2.csv"),row.names =FALSE)
  }
}







