#!/bin/Rscript
## Rscript arranged_result.R <designfile> <comparefile> <diffm6A_mode> <expression_mode> <mergePeaks_mode>
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)
args <- commandArgs(T)
#args <- c("formatted_designfile.txt", "compare_info", "DESeq2", "edgeR", "macs2")
designfile <- args[1]#"formatted_designfile.txt"
comparefile <- args[2]#"compare_info"
diffm6A_mode <- args[3]#"QNB"
rnaseq_mode <- args[4]#"DESeq2"
peakMerged_mode <- args[5]
options(stringsAsFactors = F)

## generate design matrix
compare.list <- read.csv(comparefile,header = F, check.names=F)
designtable <- read.csv(designfile, head = TRUE, colClasses = c("character"), check.names=F)
design.matrix <- as.matrix(designtable$Group)
rownames(design.matrix) <- designtable$Sample_ID
colnames(design.matrix) <- "Type"

## generate QC.peaks list
QC.reads.stats.file <- dir(".",pattern = "multiqc_rseqc_bam_stat.txt",recursive = TRUE)
QC.reads.stats <- read.table(QC.reads.stats.file, header =F, check.names=F)
colnames(QC.reads.stats)  <- QC.reads.stats[1,]
QC.reads.stats <- QC.reads.stats[-1,]

QC.reads.distribution.file <- dir(".",pattern = "multiqc_rseqc_read_distribution.txt",recursive = TRUE)
QC.reads.distribution <- read.table(QC.reads.distribution.file, sep = "\t", header =F, check.names=F)
colnames(QC.reads.distribution)  <- QC.reads.distribution[1,]
QC.reads.distribution <- QC.reads.distribution[-1,]
rownames(QC.reads.distribution)=QC.reads.distribution[,1]
QC.reads.distribution <- QC.reads.distribution[, grep("tag_pct",colnames(QC.reads.distribution))]

## generate QC.peaks list
QC.peaks.filelist <- grep("unanno.txt",grep(".anno.txt", list.files(pattern = "merged_group"), value = T), value = T, invert = T)
QC.peaks.list <- NULL
for( file in QC.peaks.filelist ){
  tmp.QC <-  read.table(file, header = F, sep = "\t", quote = "")[,c(1,2,3,15,11,12,13,14,17)]
  colnames(tmp.QC) <- c("Chr","ChrStart","ChrEnd","ID","Gene_symbol","Coding","Location","Relative_distance","RNA_type")
  name.QC <- strsplit(strsplit(file,"_group_")[[1]][2],".anno.txt")[[1]][1]
  QC.peaks.list[[name.QC]] <- tmp.QC
}

## generate the motif list of ATCG matrix
QC.motif.filelist = dir(".",pattern = "motif[1,2,3].motif",recursive = TRUE)
QC.motif.list <- NULL
QC.motif.pvalue <- NULL
for( group_name in as.vector(unique(design.matrix)[,1]) ){
  group.motif.list <- NULL
  group.motif.pvalue <- NULL
  for (file in grep(paste0(group_name, "_homer/homerResults"), QC.motif.filelist, value = T) ){
    motif_matrix <- read.delim(file,header = F,sep = "\t", check.names=F)
    motif_pvalue <- strsplit(motif_matrix[1,6], split = ":")[[1]][4]
    motif_matrix <- motif_matrix[-1,c(-5,-6)]
    colnames(motif_matrix) <- c("A","C","G","T")
    rownames(motif_matrix) <- c(1:nrow(motif_matrix))
    motif_matrix <- as.matrix(t(motif_matrix))
    motif_name <- strsplit(strsplit(file,split = c("_homer/homerResults/"))[[1]][2],split = "[.]motif")[[1]][1]
    group.motif.list[[motif_name]] <- motif_matrix
    group.motif.pvalue[[motif_name]] <- motif_pvalue
  }
  QC.motif.list[[group_name]] <- group.motif.list
  QC.motif.pvalue[[group_name]] <- group.motif.pvalue
}

## generate peak Visualization
annotation.file <- list.files(pattern = "merged_allpeaks.anno.txt")
annotation.info <- read.table(annotation.file, header = F, sep = "\t", quote = "", check.names=F)[,c(4,15,11)]
colnames(annotation.info) <- c("PeakRegion","ID","Gene_symbol")
m6a.peaks.file <- list.files(pattern = "merged_allpeaks.bed$")
m6a.peaks.table <- read.table(m6a.peaks.file, header = F, sep = "\t", quote = "", check.names=F)
colnames(m6a.peaks.table) <- c("Chr","ChrStart","ChrEnd","PeakRegion","pvalue")
m6a.peaks.table = merge(x = m6a.peaks.table,y = annotation.info,by = "PeakRegion",all.x = TRUE)
m6a.sites.file <- list.files(pattern = "m6A_sites_merged.bed")
m6a.sites.table <- read.table(m6a.sites.file, header = F, sep = "\t", quote = "", check.names=F)
colnames(m6a.sites.table) <- c("Chr","ChrStart","ChrEnd","Gene_symbol&ID","Strand","Score","Group","Sequence")

expression.matrix <- NULL
diffexpression.list <- NULL
if (rnaseq_mode != "none"){
  ## generate expression matrix
  htseq.filelist = grep("htseq",list.files(path = "./",pattern = "input.count"), value = T)
  for( file in htseq.filelist ){
    tmp.expression.table <- as.matrix(read.table(file, sep = "\t", header = TRUE, row.names = 1, check.names=F))
    expression.matrix <- cbind(expression.matrix, tmp.expression.table)
  }
  colnames(expression.matrix) <- as.matrix(lapply(strsplit(colnames(expression.matrix),".input"), function(x){ x[1]}))
  
  ## generate diff_expression list
  diffexpression.filelist <- grep(rnaseq_mode,list.files(pattern = ".csv"), value = T)
  for( compare_str in compare.list ){
    diffexpression.list[[compare_str]] <- read.csv(grep(sub("_vs_","_",compare_str), diffexpression.filelist, value = T), header = T, check.names=F)
    colnames(diffexpression.list[[compare_str]])[1] <- "ID"
  }
}

## generate m6A matrix
m6a.anno.matrix <- read.delim(file = grep("quantification.matrix",x = list.files(),value = T), header = T, sep = "\t", row.names = 1, check.names=F)
m6a.anno.matrix <- cbind(PeakRegion = row.names(m6a.anno.matrix), m6a.anno.matrix)

## generate diffm6A list
diffm6A.filelist <- grep("_diffm6A_",list.files(pattern = ".txt"), value = T)
diffm6A.list <- NULL
for( compare_str in compare.list ){
  diffm6A.list[[compare_str]] <- read.table(grep(sub("_vs_","_",compare_str), diffm6A.filelist, value = T),header = T,row.names = 1, check.names=F)
  if( diffm6A_mode == "MATK" ){
    diffm6A.list[[compare_str]]$padj = p.adjust(diffm6A.list[[compare_str]]$pvalue, method = "BH")
    diffm6A.list[[compare_str]] <- diffm6A.list[[compare_str]][,-seq(1,2)]
  }
  diffm6A.list[[compare_str]]$PeakRegion <- rownames(diffm6A.list[[compare_str]])
  diffm6A.list[[compare_str]] <- merge(x = annotation.info,y = diffm6A.list[[compare_str]],by = "PeakRegion", all.y = TRUE)
}

## save variable for m6Aviewer
write.table(expression.matrix,file = "expression.matrix",quote=F)
write.table(m6a.anno.matrix,file= "m6a.anno.matrix",quote=F)
#write.table(diffm6A.list, file= "diffm6A.list")
#write.table(diffm6A.anno.list, file = "diffm6A.anno.list")

save(design.matrix, compare.list, 
     QC.reads.stats, QC.reads.distribution,
     QC.motif.list, QC.motif.pvalue,QC.peaks.list,
     m6a.peaks.table, m6a.sites.table,
     expression.matrix, m6a.anno.matrix, 
     diffexpression.list, diffm6A.list,
     file = paste0(peakMerged_mode,"_",diffm6A_mode,"_",rnaseq_mode,"_arranged_results_",Sys.Date(),".m6APipe"))

