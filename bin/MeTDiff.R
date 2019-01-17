library(stringr)
library(MeTDiff)
args<-commandArgs(T) # Rscript MeTDiff.R aligner_tools eg. bwa
output_bed_name <- str_c("metdiff_",args[1],".bed") #diff_peak.bed

#生成designfile的table
designfile <- read.csv("designfile.txt",head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
control_input <- subset(designfile, control_or_treated == "control" & ip_or_input == "input")
control_ip <- subset(designfile, control_or_treated == "control" & ip_or_input == "ip")
treated_input <- subset(designfile, control_or_treated == "treated" & ip_or_input == "input")
treated_ip <- subset(designfile, control_or_treated == "treated" & ip_or_input == "ip")
gtf <- "genes.gtf"

#通过designfile生成对应bam文件的name，并合并到向量里
control_input_bam = vector()
for (i in c(1:nrow(control_input))) {
    tmp_bam_name=str_c(str_c(as.character(control_input[i,]),collapse='_'),"_",args[1],"_sort.bam")
    control_input_bam = c(control_input_bam,tmp_bam_name)
}
control_ip_bam = vector()
for (i in c(1:nrow(control_ip))) {
    tmp_bam_name=str_c(str_c(as.character(control_ip[i,]),collapse='_'),"_",args[1],"_sort.bam")
    control_ip_bam = c(control_ip_bam,tmp_bam_name)
}
treated_input_bam = vector()
for (i in c(1:nrow(treated_input))) {
    tmp_bam_name=str_c(str_c(as.character(treated_input[i,]),collapse='_'),"_",args[1],"_sort.bam")
    treated_input_bam = c(treated_input_bam,tmp_bam_name)
}
treated_ip_bam = vector()
for (i in c(1:nrow(treated_ip))) {
    tmp_bam_name=str_c(str_c(as.character(treated_ip[i,]),collapse='_'),"_",args[1],"_sort.bam")
    treated_ip_bam = c(treated_ip_bam,tmp_bam_name)
}

metdiff(GENE_ANNO_GTF=gtf,
        IP_BAM = control_ip_bam,
        INPUT_BAM = control_input_bam,
        TREATED_IP_BAM = treated_ip_bam,
        TREATED_INPUT_BAM=treated_input_bam,
        EXPERIMENT_NAME=args[1])
bed_name <- str_c(args[1],"/diff_peak.bed")
file.rename( bed_name , output_bed_name )
