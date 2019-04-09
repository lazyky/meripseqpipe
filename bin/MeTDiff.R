## Rscript MeTDiff.R aligner_tools designfile gtf eg. Rscript MeTDiff.R tophat2 designfile_single.txt genes.gtf
## designfile: filename, control_or_treated, input_or_ip, situation(default 1 is CONTROL_SITUATION else are TREATED_SITUATION)
## TREATED_SITUATION_STARTPOINT is the default situation check point
#!/bin/Rscript
library(stringr)
library(MeTDiff)
args <- commandArgs(T)
designfile <- args[1]
gtf <- args[2]
comparefile <- args[3]

#setting CONTROL_SITUATION and TREATED_SITUATION 
#default 1 is CONTROL_SITUATION else are TREATED_SITUATION
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
filelist = grep(".bai",list.files(path = "./",pattern = ".bam"),value = TRUE,invert = TRUE)
  bamlist <- NULL
  for(group_id in designtable$Group){
    input = grep(str_c("input_",group_id),filelist,value = TRUE)
    ip = grep(str_c("ip_",group_id),filelist,value = TRUE)
    bamlist[[group_id]] <- cbind(input,ip)
  }
  ##Running MeTPeak and rename the output name
  for (group_id in designtable$Group){
    metpeak(GENE_ANNO_GTF = gtf,
            IP_BAM = bamlist[[group_id]][,2],
            INPUT_BAM = bamlist[[group_id]][,1],
            EXPERIMENT_NAME = str_c( "metpeak_",group_id )
    )
    control_bed_name <- str_c( "metpeak_",group_id ,"/peak.bed")
    output_control_bed_name <- str_c("metpeak_group_",group_id,".bed") #peak.bed
    file.rename( control_bed_name , output_control_bed_name )
  }
}
##Running MeTDiff and rename the output name
for(i in CONTROL_SITUATION){
  for(j in TREATED_SITUATION){
    output_pattern <- str_c("metdiff_situation",i,j,aligner_tools_name,"/")
    metdiff(GENE_ANNO_GTF=gtf,
            IP_BAM = bamlist[[i]][,2],
            INPUT_BAM = bamlist[[i]][,1],
            TREATED_IP_BAM = bamlist[[j]][,2],
            TREATED_INPUT_BAM = bamlist[[j]][,1],
            EXPERIMENT_NAME = output_pattern)
    #set output_name
    output_bed_name <- str_c("metdiff_situation",i,j,aligner_tools_name,".bed") #diff_peak.bed
    bed_name <- str_c(output_pattern,"/diff_peak.bed") #choose peak.bed diff_peak.bed
    file.rename( bed_name , output_bed_name )
  }
}
