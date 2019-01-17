#!/bin/Rscript
library(stringr)
library("DESeq2")
#setwd(dir = "E:/zky/m6Apipe/bin/readscount/") #test
args<-commandArgs(T) # Rscript get_htseq_matrix.R aligner_tools eg. bwa
output_name <- str_c("treated_vs_control_",args[1],"_deseq2.csv")

control_file = str_c("control_input_",args[1],".count") 
treated_file = str_c("treated_input_",args[1],".count")
control_database=read.table(control_file,sep="\t",header=T,row.names=1)
treated_database=read.table(treated_file,sep="\t",header=T,row.names=1)
combined_database <- cbind(control_database,treated_database)
condition <- factor(c(rep("control",ncol(control_database)), rep("treated",ncol(treated_database)))) #设定因子

#assign gene names
colData <- data.frame(row.names=colnames(combined_database), condition)
dds <- DESeqDataSetFromMatrix(countData = combined_database,colData = colData,design = ~ condition)
rownames(dds) <- rownames(combined_database)
#dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)
table(res$padj <0.05)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata, file = output_name )
#
resdata2=resdata[resdata$log2FoldChange > 1|resdata$log2FoldChange < -1,]
write.csv(resdata2,file = str_c("treated_vs_control_",args[1],"_log2.csv"),row.names =FALSE)

