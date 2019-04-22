# Rscript arranged_result.R designfile comparefile
args <- commandArgs(T)
designfile <- args[1]#"formatted_designfile.txt"
comparefile <- args[2]#"compare_info.txt"
merged_peak_mode <- args[3]#"bedtools"
diffm6A_mode <- args[4]#"QNB"

## generate design matrix
designtable <- read.csv(designfile, head = TRUE, stringsAsFactors=FALSE, colClasses = c("character"))
design.matrix <- as.matrix(designtable$Group)
rownames(design.matrix) <- designtable$Sample_ID
colnames(design.matrix) <- "Type"

## generate expression matrix
htseq.filelist = grep("htseq",list.files(path = "./",pattern = "input.count"), value = T)
expression.matrix <- NULL
class(expression.matrix)
for( file in htseq.filelist ){
  tmp.expression.table <- as.matrix(read.table(file, header = TRUE, row.names = 1))
  expression.matrix <- cbind(expression.matrix, tmp.expression.table)
}
expression.matrix <- expression.matrix[c(-nrow(expression.matrix):-(nrow(expression.matrix)-5)),]
colnames(expression.matrix) <- as.matrix(lapply(strsplit(colnames(expression.matrix),".input"), function(x){ x[1]}))

## generate m6A matrix
m6a.matrix <- as.matrix(read.table(file = grep("quantification.matrix",x = list.files(),value = T), header = T, row.names = 1))
anno_merged_bed <- read.table(list.files(pattern = "_merged_peaks.refSeq.all.bed"), header = F, sep = "\t", quote = "")
colnames(anno_merged_bed) <- c("Chr_site","SiteStart","SiteEnd","Peak_name","Chr","ChrStart","ChrEnd","ENST_ID","Count","Strand","ID")

## generate diffm6A list
diffm6A.filelist <- grep("_diffm6A_",list.files(pattern = ".txt"), value = T)
compare.list <- read.csv(comparefile,header = F)
diffm6A.list <- NULL
for( compare_str in compare.list ){
  compare_str = "control_vs_METLT14_KD"
  sub("_vs_","_",compare_str)
  diffm6A.list[[compare_str]] <- read.table(grep(sub("_vs_","_",compare_str), diffm6A.filelist, value = T))
  if( diffm6A_mode == "QNB"){
    colnames(diffm6A.list[[compare_str]]) <- c("p.treated","p.control","log2FC","log2.OR","pvalue","q","padj")
  }
}

## generate QC list
QC.filelist <- grep("unanno.txt",grep(".anno.txt", list.files(pattern = paste0(merged_peak_mode,"_group")), value = T), value = T, invert = T)
QC.list <- NULL
for( file in QC.filelist ){
  tmp.QC <-  read.table(file , header = F, sep = "\t", quote = "")[,c(1,2,3,15,11,12,13,14,17)]
  colnames(tmp.QC) <- c("Chr","ChrStart","ChrEnd","ID","Gene_symbol","Coding","Location","Count","RNA_type")
  name.QC <- strsplit(strsplit(file,"_group_")[[1]][2],".anno.txt")[[1]][1]
  QC.list[[name.QC]] <- tmp.QC
}

## save variable for m6Aviewer
save(design.matrix, expression.matrix, m6a.matrix, diffm6A.list, QC.list, compare.list, anno_merged_bed,
     file ="arranged_results.m6APipe")