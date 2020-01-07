#!/bin/Rscript
## generate DESeq2 logfc matrix 
anno.exp.m6a <- subset(m6a.anno.matrix,select= c(PeakRegion,ID))
final.count.table <- NULL
htseq_list <- dir(pattern = "count",full.names = T)
combined_htseq_count <- read.table(htseq_list[1], header = TRUE, row.names = 1, check.names = FALSE)
for (file in htseq_list[-1]){
  combined_htseq_count <- cbind(combined_htseq_count,read.table(file, header = TRUE, row.names = 1, check.names = FALSE))
}
colnames(combined_htseq_count) <- unlist(lapply(strsplit(colnames(combined_htseq_count),split = "_"),FUN = function(x){x[1]}))
combined_htseq_count$ID <- rownames(combined_htseq_count)
final.count.table <- merge(anno.exp.m6a,combined_htseq_count,by= "ID")

sample_id <- designtable$Sample_ID[1]
ip_count_file <- grep(paste0("[.]",sample_id,"[.]ip"),count_filelist,value = TRUE)
ip.matrix <- read.table(file = ip_count_file, sep = "\t", row.names = NULL, header = T, check.names=F)
rownames.ip.matrix <- ip.matrix$PeakName
ip.matrix <- subset(ip.matrix , select= 5)
colnames(ip.matrix) <- paste0(sample_id,".ip")
for(sample_id in designtable$Sample_ID[-1]){
  ip_count_file <- grep(paste0("[.]",sample_id,"[.]ip"),count_filelist,value = TRUE)
  ip_count_table <- read.table(file = ip_count_file, sep = "\t", row.names = NULL, header = T, check.names=F)
  ip_count_table <- subset(ip_count_table , select= 5)
  colnames(ip_count_table) <- paste0(sample_id,".ip")
  ip.matrix <- cbind(ip.matrix,ip_count_table)
}
row.names(ip.matrix) <- rownames.ip.matrix
ip.matrix$PeakRegion <- rownames(ip.matrix)
final.count.table <- merge(final.count.table,ip.matrix,by= "PeakRegion")

load("deseq2.Rdata")
library("DESeq2")
deseq2.count.table <- subset(final.count.table,select= colnames(final.count.table)[c(-1,-2)])
rownames(deseq2.count.table) <- final.count.table$PeakRegion
coldata <- data.frame(row.names = colnames(deseq2.count.table),group = colnames(deseq2.count.table) ,sample = unlist(lapply(strsplit(colnames(deseq2.count.table),split = "[.]"),FUN = function(x){x[1]})))
final.deseq2.logfc.matrix <- subset(final.count.table, select = PeakRegion)
for (sample_id in unique(coldata$sample)) {
  coldata.sample <- subset(coldata,sample == sample_id,group)
  coldata.sample$group <- unlist(lapply(strsplit(as.character(coldata.sample$group),split = "[.]"),function(x){x[2]}))
  coldata.sample$group <- factor(coldata.sample$group)
  inf.dds <- DESeq2::DESeqDataSetFromMatrix(countData = deseq2.count.table[,rownames(coldata.sample)],colData = coldata.sample,design = ~group)
  inf.dds.LRT <- DESeq2::DESeq(inf.dds)
  head(deseq2.count.table[,rownames(coldata.sample)])
  results <- DESeq2::results(inf.dds.LRT,constract=c("group","input","ip"))
  results <- data.frame(PeakRegion = rownames(results),sample_id = results$log2FoldChange)
  colnames(results)[2] <- sample_id
  final.deseq2.logfc.matrix <- merge(final.deseq2.logfc.matrix,results,by= "PeakRegion")
}
2^head(final.deseq2.logfc.matrix)
exp()final.deseq2.logfc.matrix
rownames(final.deseq2.logfc.matrix) <- final.deseq2.logfc.matrix$PeakRegion
final.deseq2.logfc.matrix <- final.deseq2.logfc.matrix[,-1]
final.deseq2.logfc.matrix <- 2^(final.deseq2.logfc.matrix)
save(final.deseq2.logfc.matrix,file = "deseq.quantification.matrix.RData")