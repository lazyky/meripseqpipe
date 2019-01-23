## Rscript diffexomePeak.R aligner_tools eg. bwa
library(stringr)
library("exomePeak")
args <- commandArgs(T) 
aligner_tools_name <- args[1]
designfile <- args[2]
gtf <- args[3]
TREATED_SITUATION_STARTPOINT <- 1

##setting CONTROL_SITUATION and TREATED_SITUATION 
#default 0 is CONTROL_SITUATION else are TREATED_SITUATION
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
CONTROL_SITUATION <- c()
TREATED_SITUATION <- c()
for (i in c(0:(TREATED_SITUATION_STARTPOINT-1))){
  CONTROL_SITUATION <- c(CONTROL_SITUATION,str_c("_",i,"_"))
}
for (i in c(TREATED_SITUATION_STARTPOINT:max(designtable$situation))){
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

##Running MeTDiff and rename the output name
for(i in CONTROL_SITUATION){
  for(j in TREATED_SITUATION){
    output_pattern <- str_c("diffexomePeak_situation",i,j,aligner_tools_name)
    exomepeak(GENE_ANNO_GTF = gtf,
              OUTPUT_DIR = output_pattern,
              IP_BAM = bamlist[[i]][,2],
              INPUT_BAM = bamlist[[i]][,1],
              TREATED_IP_BAM = bamlist[[j]][,2],
              TREATED_INPUT_BAM = bamlist[[j]][,1]
    )
    #set output_name
    output_bed_name <- str_c("diffexomePeak_situation",i,j,aligner_tools_name,".bed") #con_sig_diff_peak.bed
    bed_name <- str_c(output_pattern,"/exomePeak_output/","con_sig_diff_peak.bed") #choose diff_peak.bed sig_diff_peak.bed con_sig_diff_peak
    file.rename( bed_name , output_bed_name )  
  }
}

