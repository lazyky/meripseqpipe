## Rscript MeTPeak.R aligner_tools designfile gtf eg. Rscript MeTPeak.R tophat2 designfile_single.txt genes.gtf
## designfile: filename, control_or_treated, input_or_ip, situation(default 1 is CONTROL_SITUATION else are TREATED_SITUATION)
## TREATED_SITUATION_STARTPOINT is the default situation check point
#!/bin/Rscript
library(stringr)
library(MeTPeak)
args <- commandArgs(T) 
TREATED_SITUATION_STARTPOINT <- args[1]
aligner_tools_name <- args[2]
designfile <- args[3]
gtf <- args[4]
TREATED_SITUATION_STARTPOINT <- as.numeric(TREATED_SITUATION_STARTPOINT)

##setting CONTROL_SITUATION and TREATED_SITUATION 
#default 1 is CONTROL_SITUATION else are TREATED_SITUATION
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
CONTROL_SITUATION <- c()
TREATED_SITUATION <- c()
if (TREATED_SITUATION_STARTPOINT-1 >= 1){
  for (i in c(1:(TREATED_SITUATION_STARTPOINT-1))){
  CONTROL_SITUATION <- c(CONTROL_SITUATION,str_c("_",i,"_"))
  }
}
for (i in c(TREATED_SITUATION_STARTPOINT:max(as.numeric(designtable$situation)))){
  TREATED_SITUATION <- c(TREATED_SITUATION,str_c("_",i,"_"))
}

##Traversing all situations
filelist = grep(aligner_tools_name,grep(".bai",list.files(path = "./",pattern = ".bam"),value = TRUE,invert = TRUE),value = TRUE)
bamlist <- NULL
for(i in c(CONTROL_SITUATION,TREATED_SITUATION)){
  if (i %in% CONTROL_SITUATION){
    a = grep(str_c("input",i),filelist,value = TRUE)
    b = grep(str_c("ip",i),filelist,value = TRUE)
    bamlist[[i]] <- cbind(a,b)
  }
  if (i %in% TREATED_SITUATION){
    a = grep(str_c("input",i),filelist,value = TRUE)
    b = grep(str_c("ip",i),filelist,value = TRUE)
    bamlist[[i]] <- cbind(a,b)
  }
}

##Running MeTPeak and rename the output name
output_pattern <- str_c("metpeak_",aligner_tools_name,"_")
for (i in CONTROL_SITUATION){
  metpeak(GENE_ANNO_GTF=gtf,
          IP_BAM = bamlist[[i]][,2],
          INPUT_BAM = bamlist[[i]][,1],
          EXPERIMENT_NAME = str_c( output_pattern, "control",i )
  ) 
  control_bed_name <- str_c(output_pattern,"control",i,"/peak.bed")
  output_control_bed_name <- str_c("metpeak_situaion",i,aligner_tools_name,".bed") #peak.bed
  file.rename( control_bed_name , output_control_bed_name )
}
for (i in TREATED_SITUATION){
  metpeak(GENE_ANNO_GTF=gtf,
          IP_BAM = bamlist[[i]][,2],
          INPUT_BAM = bamlist[[i]][,1],
          EXPERIMENT_NAME = str_c( output_pattern, "treated",i )
  )
  treated_bed_name <- str_c(output_pattern,"treated",i,"/peak.bed")
  output_treated_bed_name <- str_c("metpeak_situaion",i,aligner_tools_name,".bed") #peak.bed
  file.rename( treated_bed_name , output_treated_bed_name )
}



