#!/bin/Rscript
## Rscript GLM_edgeR_DM.R <desginfile> <compare_str> <THREAD_NUM>
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)

#####GLM model###
#edgeR
library(edgeR)
library(DESeq2)
library(readr)
library(dplyr)
library(BiocParallel)
#load data
designfile <- args[1]
compare_str <- args[2]
THREAD_NUM <- args[3]
register(MulticoreParam(THREAD_NUM))

designtable <- read.csv(designfile, head = TRUE, stringsAsFactors=FALSE, colClasses = c("character"), check.names=F)
design.matrix <- as.data.frame(designtable$Group)
rownames(design.matrix) <- designtable$Sample_ID
colnames(design.matrix) <- "Condition"

# Get the information of groups from compare_str
if(length(unique(design.matrix$Condition)) < 2){
  stop( "The count of Group is less than two, please check your designfile.")
}else if( compare_str == "two_group" ){
  # Get the information without compare_str beacause of only two groups
  group_id_1 <- unique(design.matrix$Condition)[1]
  group_id_2 <- unique(design.matrix$Condition)[2]
}else{
  # Running MeTDiff quantification with compare_str
  group_id_1 <- strsplit(as.character(compare_str), "_vs_")[[1]][1]
  group_id_2 <- strsplit(as.character(compare_str), "_vs_")[[1]][2]
}
design.matrix <- subset(design.matrix, Condition == group_id_1 | Condition == group_id_2 )
design.matrix$Condition <- factor(design.matrix$Condition,labels = c("control","treatment"))
filelist = list.files(path = "count/peak_count",pattern = ".count",full.names = T)
## Generate the matrix of peaks count
rpkm_peaks_list <- NULL
for(sample_id in row.names(design.matrix)){
  input_count_file <- grep(paste0("[.]",sample_id,"[.]input"),filelist,value = TRUE)
  input_count_table <- read.table(file = input_count_file, sep = "\t", row.names = NULL,header = T)

  ip_count_file <- grep(paste0("[.]",sample_id,"[.]ip"),filelist,value = TRUE)
  ip_count_table <- read.table(file = ip_count_file, sep = "\t", row.names = NULL, header = T)
  rpkm <- cbind(input_count_table[,5],ip_count_table[,5])
  colnames(rpkm) <- c(paste0(sample_id,".input"),paste0(sample_id,".ip"))
  rpkm_peaks_list <- cbind(rpkm_peaks_list,rpkm)
}
rownames(rpkm_peaks_list) <- ip_count_table$PeakName

## generate design matrix
design.matrix$m6A <- "input"
design.matrix$sample_id <- paste0(rownames(design.matrix),".input")
design.matrix_ip <- design.matrix
design.matrix_ip$m6A <- "IP"
design.matrix_ip$sample_id <- paste0(rownames(design.matrix_ip),".ip")
design.matrix <- rbind(design.matrix,design.matrix_ip)
rownames(design.matrix) <- design.matrix$sample_id
design.matrix$m6A <- factor(design.matrix$m6A)
design.matrix <- design.matrix[colnames(rpkm_peaks_list),]

design1 <- model.matrix(~ Condition + m6A + Condition:m6A, data=design.matrix)
design2 <- model.matrix(~ Condition + m6A , data=design.matrix)

#construct deglst
deg_lst <- DGEList(counts = rpkm_peaks_list ,group = design.matrix$Condition)
deg_lst <- calcNormFactors(deg_lst)

#disper and test  
deg_lst1<-estimateDisp(deg_lst,design1)
deg_lst2<-estimateDisp(deg_lst,design2)
#GLM fit

fit1 <- glmFit(deg_lst1, design1)
colnames(fit1)
#fit1$deviance
fit2 <- glmFit(deg_lst2, design2)
#fit2$deviance
##chip test nested interaction
deviance.chip.test <- chisq.test(as.numeric(fit1$deviance), as.numeric(fit2$deviance))
print(deviance.chip.test)
## 
lrt <- glmLRT(fit1,coef = 4)
coeffients.df = as.data.frame(lrt$coefficients)
coeffients.df$id = rownames(coeffients.df)
colnames(coeffients.df) = c("base","condition","m6A","meth_in_condition","id")
results2 <- as.data.frame(cbind(lrt$table$PValue,p.adjust(lrt$table$PValue)))
head(results)
head(results2)
sum(results2$V1<0.05)
sum(results$edger.p<0.05)
topTags(lrt)
glm_DM.res <- as.data.frame(lrt$table)
glm_DM.res$id <- rownames(glm_DM.res)

glm_DM.res.df <- data.frame(glm_DM.res, FDR = p.adjust(glm_DM.res$PValue, method="BH"))%>%
  inner_join(coeffients.df,by="id")
colnames(glm_DM.res.df)
glm_DM.res.df = glm_DM.res.df%>%mutate(meth_status = case_when( meth_in_condition > 0 & PValue <= 0.05 ~ "hyper",
                                                                meth_in_condition < 0 & PValue <= 0.05 ~ "hypo",
                                                              TRUE ~ "non_diff"))

table(glm_DM.res.df$meth_status)


## DESeq2
dds <- DESeqDataSetFromMatrix(countData = rpkm_peaks_list,
                              colData = design.matrix,
                              design= ~ m6A + Condition + Condition:m6A )
dds <- DESeq(dds,parallel = T)
substr(names(mcols(dds)),1,10) 
res <- results(dds,name = resultsNames(dds)[4])
sum(res$padj < 0.1, na.rm=TRUE)

dds2 <- DESeqDataSetFromMatrix(countData = rpkm_peaks_list,
                              colData = design.matrix,
                              design=  ~ m6A + Condition )
dds2 <- DESeq(dds2,parallel = T)
resultsNames(dds) # lists the coefficients

deseq.deviance.chip.test <- chisq.test(mcols(dds)$deviance, mcols(dds2)$deviance)
print(deseq.deviance.chip.test)

deseq_results <- run.deseq2(rpkm_peaks_list,design.matrix)
edgeR_results <- run.edger(rpkm_peaks_list,design.matrix)
qnb_results <- run.qnb(rpkm_peaks_list,design.matrix)
metdiff.results <- run.metdiff(rpkm_peaks_list,design.matrix)
diffm6A_results <- cbind(deseq_results,edgeR_results,qnb_results,metdiff_results)
rownames(diffm6A_results) <- rownames(rpkm_peaks_list)

print(paste0("deseq_diffm6A: ",sum(diffm6A_results$deseq2.p<=0.05,na.rm = T)))
print(paste0("edgeR_diffm6A: ",sum(diffm6A_results$edger.p<=0.05,na.rm = T)))
print(paste0("qnb_diffm6A: ",sum(diffm6A_results$qnb.p<=0.05,na.rm = T)))
print(paste0("metdiff_diffm6A: ",sum(diffm6A_results$metdiff.p<=0.05,na.rm = T)))
write.table(diffm6A_results,file = "deqm_diffm6A_results",sep = "\t",quote = F)

