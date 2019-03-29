## Rscript arranged_results.R aligners_tools_name peak_calling_tools_name 
aligner_tools_name <- "aligners"#args[1]
peak_calling_tools_name <- "MATK" #args[2]
bam_stat_summary <- "bam_stat_summary.txt"
getwd()
setwd("/data2/yeying_by_zky/EBV_pipe/work/cf/8346833244f3b36ed0c12ba0a38809")
bam_stat_table <- read.table(bam_stat_summary,row.names = 1)
filelist = grep(paste0(aligner_tools_name,"_",peak_calling_tools_name),list.files(path = "./",pattern = "input.count"),value = TRUE)
rpkm_peaks_list <- NULL
rpkm_peaks_list_index <- 0
class(rpkm_peaks_list)
for (file in filelist){
  count_table <- read.table(file = file, sep = "\t", row.names = NULL)
  colnames(count_table) <- as.vector(as.matrix(count_table[1,]))
  count_table <- count_table[-1,]
  for (count_table_index in c(5:length(count_table[1,]))){
    bam_stat_index <- as.numeric(which( row.names(bam_stat_table) == colnames(count_table)[count_table_index] ))
    rpkm =  apply(count_table,1,function(x) (as.numeric(x[count_table_index])/(as.numeric(x[3])-as.numeric(x[2]))/bam_stat_table[bam_stat_index,]*1000))
    rpkm <- as.matrix(rpkm)
    colnames(rpkm) <- colnames(count_table)[count_table_index]
    rpkm_peaks_list_index = rpkm_peaks_list_index + 1
      rpkm_peaks_list <- cbind(rpkm_peaks_list,rpkm)
  }
}
write.table(rpkm_peaks_list,file = paste0(aligner_tools_name,"_",peak_calling_tools_name))