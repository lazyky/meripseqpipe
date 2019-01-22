## Rscript get_htseq_matrix.R aligner_tools designfile gtf eg. Rscript get_htseq_matrix.R tophat2 designfile_single.txt
## designfile: filename, control_or_treated, input_or_ip, situation(default 0 is CONTROL_SITUATION else are TREATED_SITUATION)
## TREATED_SITUATION_STARTPOINT is the default situation check point
#!/bin/Rscript
library(stringr)
library("QNB")
args<-commandArgs(T) 
aligner_tools_name <- args[1]
designfile <- args[2]
TREATED_SITUATION_STARTPOINT <- 1

#setting CONTROL_SITUATION and TREATED_SITUATION 
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

# When there are replicates under two conditions, we could select 
# mode="per-condition" or mode="pooled" to estimate the dispersion. 
# The default is mode="auto".
##Running MeTDiff and rename the output name
for(i in CONTROL_SITUATION){
  for(j in TREATED_SITUATION){
    control_ip <- read.table(str_c("situation", i, aligner_tools_name, "_ip.count"), header = TRUE, row.names = 1)
    control_input <- read.table(str_c("situation",i, aligner_tools_name, "_input.count"), header = TRUE, row.names = 1)
    treated_ip <- read.table(str_c("situation", j, aligner_tools_name, "_ip.count"), header = TRUE, row.names = 1)
    treated_input <- read.table(str_c("situation", j, aligner_tools_name, "_input.count"), header = TRUE, row.names = 1)
    #set output_name
    result <- qnbtest(control_ip, treated_ip,control_input,treated_input,mode="auto")
    output_name <- str_c("QNB_situation",i,j,"QNB.csv")
    write.csv(result, file = output_name)
  }
}

#dispersion.pdf
#dif_meth.xls
## End(Not run)