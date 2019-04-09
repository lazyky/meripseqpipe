
## get expression matrix
filelist = grep("htseq_group",list.files(path = "./",pattern = "aligners_input.count"),value = TRUE)
expression_matrix <- NULL
for (file in filelist){
  part_matrix <- as.matrix(read.table(file, header = TRUE, row.names = 1))
  expression_matrix <- cbind(expression_matrix,part_matrix)
}
expression_matrix_final <- as.data.frame(expression_matrix)

## get m6A matrix
m6A_matrix_final <- read.table("final_merged.matrix", header = TRUE, row.names = 1)
m6A_diff_result <- NULL

## get information for m6A QC
QC_for_m6A <- read.table("bedtools_merged_peaks.anno.txt", header = F, sep = "\t", quote = "")[,c(1,2,3,15,11,12,13,14,17)]
colnames(QC_for_m6A) <- c("Chr","ChrStart","ChrEnd","ENSG_ID","Gene_name","Coding","Location","Count","RNA_type")

## get peaks information
merged_bed <- read.table("bedtools_merged_peaks.bed", header = F, sep = "\t", quote = "")

## sample_group_info
expression_group_info <-  read.table("expression_group.file", header = TRUE, row.names = 1)
m6A_group_info <-  read.table("m6A_group.file", header = TRUE, row.names = 1)

## group_compare_info
group_compare_info <- read.table("group_compare_info.file", header = F, sep = ",")
length(rownames(group_compare_info))
matrix_list <- NULL
for (i in c(1:length(rownames(group_compare_info)))){
  colnames(m6A_matrix_final)
  tmp <- m6A_matrix_final[,c(1,3)]
}
## save variable for m6Aviewer
save(expression_matrix_final, m6A_matrix_final, QC_for_m6A, merged_bed, 
     m6A_diff_result, expression_group_info, m6A_group_info, 
     file ="arranged_results.m6APipe")
