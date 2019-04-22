## Rscript MeTPeak.R aligner_tools designfile gtf eg. Rscript MeTPeak.R tophat2 designfile_single.txt genes.gtf
## designfile: filename, input_or_ip, situation(default 1 is CONTROL_SITUATION else are TREATED_SITUATION)
#!/bin/Rscript
library(MeTPeak)
library(parallel)
args <- commandArgs(T) 
designfile <- args[1]
gtf <- args[2]
THREAD_NUM <- as.numeric(args[3])
flag_peakCallingbygroup <- as.numeric(args[4])
##setting CONTROL_SITUATION and TREATED_SITUATION 
#default 1 is CONTROL_SITUATION else are TREATED_SITUATION
designtable <- read.csv(designfile,header = TRUE ,stringsAsFactors=FALSE, colClasses = c("character"))
##Traversing all situations
filelist = grep(".bai",list.files(path = "./",pattern = ".bam"),value = TRUE,invert = TRUE)
if(flag_peakCallingbygroup){
  bamlist <- NULL
  for(group_id in unique(designtable$Group)){
    input = grep(paste0("input_",group_id),filelist,value = TRUE)
    ip = grep(paste0("ip_",group_id),filelist,value = TRUE)
    bamlist[[group_id]] <- cbind(input,ip)
  }
  ##Running MeTPeak and rename the output name
  mclapply(unique(designtable$Group),function(x){
    group_id = x
    metpeak(GENE_ANNO_GTF = gtf,
            IP_BAM = bamlist[[group_id]][,2],
            INPUT_BAM = bamlist[[group_id]][,1],
            EXPERIMENT_NAME = paste0( "metpeak_",group_id )
    )
    bed_name <- paste0( "metpeak_",group_id ,"/peak.xls")
    output_bed_name <- paste0("metpeak_group_",group_id,"_normalized.bed") #peak.bed
    bed12.to.bed6 <- paste0("awk 'BEGIN{OFS=\"\t\"}NR>1{print $1,$2,$3,$4,-$13,$6,$7,$8,$9,$10,$11,$12}' ", bed_name," | bed12ToBed6 -i > ", output_bed_name)
    system(bed12.to.bed6)
    },
    mc.cores = THREAD_NUM)
}else{
  if(length(designtable$Sample_ID) > length(unique(designtable$Sample_ID))) stop("Sample_id is repeat")
  mclapply(unique(designtable$Sample_ID),function(x){
    sample_id=x
    sample_vector <- grep(sample_id,filelist,value = TRUE)
    metpeak(GENE_ANNO_GTF = gtf,
            IP_BAM = sample_vector[2],
            INPUT_BAM = sample_vector[1],
            EXPERIMENT_NAME = paste0( "metpeak_",sample_id )
    )
    bed_name <- paste0( "metpeak_",sample_id ,"/peak.xls")
    output_bed_name <- paste0("metpeak_",sample_id,"_normalized.bed") #peak.bed
    bed12.to.bed6 <- paste0("awk 'BEGIN{OFS=\"\t\"}NR>1{print $1,$2,$3,$4,-$13,$6,$7,$8,$9,$10,$11,$12}' ", bed_name," | bed12ToBed6 -i > ", output_bed_name)
    system(bed12.to.bed6)
  },
  mc.cores = THREAD_NUM)
}

