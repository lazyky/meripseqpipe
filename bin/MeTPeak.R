## Rscript MeTPeak.R aligner_tools designfile gtf eg. Rscript MeTPeak.R tophat2 designfile_single.txt genes.gtf
## designfile: filename, input_or_ip, situation(default 1 is CONTROL_SITUATION else are TREATED_SITUATION)
#!/bin/Rscript
library(stringr)
library(MeTPeak)
args <- commandArgs(T) 
aligner_tools_name <- args[1]
designfile <- args[2]
gtf <- args[3]

##setting CONTROL_SITUATION and TREATED_SITUATION 
#default 1 is CONTROL_SITUATION else are TREATED_SITUATION
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
##Traversing all situations
filelist = grep(aligner_tools_name,grep(".bai",list.files(path = "./",pattern = ".bam"),value = TRUE,invert = TRUE),value = TRUE)
bamlist <- NULL
for(i in c(1:max(as.numeric(designtable$situation)))){
  a = grep(str_c("input_",i,"_"),filelist,value = TRUE)
  b = grep(str_c("ip_",i,"_"),filelist,value = TRUE)
  bamlist[[i]] <- cbind(a,b)
}

##Running MeTPeak and rename the output name
output_pattern <- str_c("metpeak_",aligner_tools_name,"_")
for (i in c(1:max(as.numeric(designtable$situation)))){
  metpeak(GENE_ANNO_GTF=gtf,
          IP_BAM = bamlist[[i]][,2],
          INPUT_BAM = bamlist[[i]][,1],
          EXPERIMENT_NAME = str_c( output_pattern,i )
  ) 
  control_bed_name <- str_c(output_pattern,i,"/peak.bed")
  output_control_bed_name <- str_c("metpeak_situation_",i,"_",aligner_tools_name,".bed") #peak.bed
  file.rename( control_bed_name , output_control_bed_name )
}
