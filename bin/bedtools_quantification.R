## Rscript arranged_results.R aligners_tools_name peak_calling_tools_name 
args <- commandArgs(T)
designfile <- args[1]
bam_stat_summary <- args[2]#"bam_stat_summary.txt"
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
bam_stat_table <- read.table(bam_stat_summary,row.names = 1)
filelist =list.files(path = "./",pattern = ".count")
rpkm_peaks_list <- NULL
for(sample_id in designtable$Sample_ID){
  input_count_file <- grep(paste0(sample_id,".input"),filelist,value = TRUE)
  input_count_table <- read.table(file = input_count_file, sep = "\t", row.names = NULL,header = T)
  bam_stat_index = grep(paste0(sample_id,".input"),rownames(bam_stat_table))
  input_rpkm =  apply(input_count_table,1,function(x) (as.numeric(x[5])/(as.numeric(x[3])-as.numeric(x[2]))*1000/bam_stat_table[bam_stat_index,]*1000000)+1)
  
  ip_count_file <- grep(paste0(sample_id,".ip"),filelist,value = TRUE)
  ip_count_table <- read.table(file = ip_count_file, sep = "\t", row.names = NULL, header = T)
  bam_stat_index = grep(paste0(sample_id,".ip"),rownames(bam_stat_table))
  ip_rpkm =  apply(ip_count_table,1,function(x) (as.numeric(x[5])/(as.numeric(x[3])-as.numeric(x[2]))*1000/bam_stat_table[bam_stat_index,]*1000000)+1)
  
  rpkm <- as.matrix(ip_rpkm/input_rpkm)
  colnames(rpkm) <- sample_id
  rpkm_peaks_list <- cbind(rpkm_peaks_list,rpkm)
}
rownames(rpkm_peaks_list) <- input_count_table[,4] 
write.table(rpkm_peaks_list,file = "bedtools_quantification.matrix")
