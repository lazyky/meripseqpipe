#!/bin/Rscript
## Rscript arranged_result.R <designfile> <comparefile> <diffm6A_mode>
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)
args <- commandArgs(T)
designfile <- args[1]#"formatted_designfile.txt"
comparefile <- args[2]#"compare_info.txt"
diffm6A_mode <- args[3]#"QNB"

## generate design matrix
designtable <- read.csv(designfile, head = TRUE, stringsAsFactors=FALSE, colClasses = c("character"))
design.matrix <- as.matrix(designtable$Group)
rownames(design.matrix) <- designtable$Sample_ID
colnames(design.matrix) <- "Type"

## generate expression matrix
htseq.filelist = grep("htseq",list.files(path = "./",pattern = "input.count"), value = T)
expression.matrix <- NULL
for( file in htseq.filelist ){
  tmp.expression.table <- as.matrix(read.table(file, header = TRUE, row.names = 1))
  expression.matrix <- cbind(expression.matrix, tmp.expression.table)
}
expression.matrix <- expression.matrix[c(-nrow(expression.matrix):-(nrow(expression.matrix)-4)),]
colnames(expression.matrix) <- as.matrix(lapply(strsplit(colnames(expression.matrix),".input"), function(x){ x[1]}))

## generate m6A matrix
m6a.matrix <- as.matrix(read.table(file = grep("quantification.matrix",x = list.files(),value = T), header = T, row.names = 1))
annotation.file <- list.files(pattern = "merged_peaks.anno.txt")
anno.merged.bed <- read.table(annotation.file, header = F, sep = "\t", stringsAsFactors = F, quote = "")[,c(4,15,11,12,13,14,17)]
colnames(anno.merged.bed) <- c("PeakRegion","ID","Gene_symbol","Coding","Location","Count","RNA_type")
annotation.info <- as.matrix(anno.merged.bed$ID)
rownames(annotation.info) <- anno.merged.bed$PeakRegion
colnames(annotation.info) <- "ID"
m6a.anno.matrix <- as.data.frame(m6a.matrix[which(row.names(m6a.matrix) %in% rownames(annotation.info)),])
m6a.anno.matrix <- cbind(annotation.info[rownames(m6a.anno.matrix),],m6a.anno.matrix)
colnames(m6a.anno.matrix)[1] <- "ID"
m6a.anno.matrix <- aggregate(m6a.anno.matrix, by=list(m6a.anno.matrix$ID), FUN=mean)
m6a.anno.matrix <- subset(m6a.anno.matrix, select = -ID)
colnames(m6a.anno.matrix)[1] <- "ID"

## generate diffm6A list
diffm6A.filelist <- grep("_diffm6A_",list.files(pattern = ".txt"), value = T)
compare.list <- read.csv(comparefile,header = F,stringsAsFactors = F)
diffm6A.list <- NULL
diffm6A.anno.list <-NULL
for( compare_str in compare.list ){
  diffm6A.list[[compare_str]] <- read.table(grep(sub("_vs_","_",compare_str), diffm6A.filelist, value = T),header = T,row.names = 1)
  if( diffm6A_mode == "QNB" ){
    colnames(diffm6A.list[[compare_str]]) <- c("p.treated","p.control","log2FC","log2.OR","pvalue","qvalue","padj")
  }else if( diffm6A_mode == "MeTDiff" ){
    
  }else if( diffm6A_mode == "MATK" ){
    diffm6A.list[[compare_str]]$padj = p.adjust(diffm6A.list[[compare_str]]$pvalue, method = "BH")
    diffm6A.list[[compare_str]] <- subset(diffm6A.list[[compare_str]], select = -ID)
    diffm6A.list[[compare_str]] <- subset(diffm6A.list[[compare_str]], select = -GeneName)
  }else if( diffm6A_mode == "bedtools" ){
    
  }else{
    stop("Please check your setting of quantification_mode")
  }
  ## Adding the ID for Peak Region
  diffm6A.anno.list[[compare_str]] <- diffm6A.list[[compare_str]][which(row.names(diffm6A.list[[compare_str]]) %in% rownames(annotation.info)),]
  diffm6A.anno.list[[compare_str]] <- cbind(annotation.info[rownames(diffm6A.anno.list[[compare_str]]),],diffm6A.anno.list[[compare_str]])
  colnames(diffm6A.anno.list[[compare_str]])[1] <- "ID"
  # diffm6A.anno.list[[compare_str]] <- aggregate(diffm6A.anno.list[[compare_str]], by=list(diffm6A.anno.list[[compare_str]]$ID), FUN=mean)
  # diffm6A.anno.list[[compare_str]] <- subset(diffm6A.anno.list[[compare_str]], select = -ID)
  # colnames(diffm6A.anno.list[[compare_str]])[1] <- "ID"
}
## generate QC list
QC.filelist <- grep("unanno.txt",grep(".anno.txt", list.files(pattern = "merged_group"), value = T), value = T, invert = T)
QC.list <- NULL
for( file in QC.filelist ){
  tmp.QC <-  read.table(file , header = F, sep = "\t", quote = "")[,c(1,2,3,15,11,12,13,14,17)]
  colnames(tmp.QC) <- c("Chr","ChrStart","ChrEnd","ID","Gene_symbol","Coding","Location","Count","RNA_type")
  name.QC <- strsplit(strsplit(file,"_group_")[[1]][2],".anno.txt")[[1]][1]
  QC.list[[name.QC]] <- tmp.QC
}

## generate the motif list of ATCG matrix
motif.filelist = dir(".",pattern = "motif3.motif",recursive = TRUE)
motif.list <- NULL
motif.pvalue <- NULL
for( file in motif.filelist ){
  motif_matrix <- read.delim(file,header = F,sep = "\t",stringsAsFactors = F)
  group_name <- strsplit(strsplit(file,split = c("_homer/homerResults"))[[1]][1],split = "_group_")[[1]][2]
  pvalue <- strsplit(motif_matrix[1,6], split = ":")[[1]][4]
  motif_matrix <- motif_matrix[-1,c(-5,-6)]
  colnames(motif_matrix) <- c("A","C","G","T")
  rownames(motif_matrix) <- c(1:nrow(motif_matrix))
  motif_matrix <- as.matrix(t(motif_matrix))
  motif.list[[group_name]] <- motif_matrix
  motif.pvalue[[group_name]] <- pvalue
}
motif.pvalue = as.data.frame(t(motif.pvalue))

## save variable for m6Aviewer
write.table(expression.matrix,file = "expression.matrix")
write.table(m6a.matrix, file = "m6a.matrix")
write.table(m6a.anno.matrix,file= "m6a.anno.matrix")
#write.table(diffm6A.list, file= "diffm6A.list")
#write.table(diffm6A.anno.list, file = "diffm6A.anno.list")

save(design.matrix, expression.matrix, m6a.matrix,m6a.anno.matrix, 
     diffm6A.list,diffm6A.anno.list, QC.list, compare.list, annotation.info,
     motif.list, motif.pvalue,
     file = paste(diffm6A_mode,"arranged_results.m6APipe", sep = "_"))