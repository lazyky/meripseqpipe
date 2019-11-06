library(DESeq2)
library(BiocParallel)

#load data
args <- commandArgs(T)
#args <- c("formatted_designfile.txt","shGFPa549_vs_shMettl3a549", "10")
designfile <- args[1]
compare_str <- args[2]
THREAD_NUM <- as.numeric(args[3])
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
filelist = list.files(path = ".",pattern = ".count",full.names = T)
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

run.deseq2 <- function(cnts,meta){
  inf.dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnts,colData = meta,design = ~Condition+m6A+Condition:m6A)
  inf.dds.LRT <- DESeq2::DESeq(inf.dds,betaPrior=FALSE, test="LRT",
                               full=~Condition+m6A+Condition:m6A,reduced=~Condition+m6A)    
  inf.dds.res <- DESeq2::results(inf.dds.LRT)
  results <- inf.dds.res
  colnames(results) <- c("baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj")
  return(results)
}
results <- run.deseq2(rpkm_peaks_list,design.matrix)
write.table(results,file = paste0("DESeq2_diffm6A_",group_id_1, "_",group_id_2,".txt") ,sep = "\t",quote = F)
