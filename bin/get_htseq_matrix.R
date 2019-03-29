## Rscript get_htseq_matrix.R aligner_tools designfile gtf eg. Rscript get_htseq_matrix.R tophat2 designfile_single.txt
## designfile: filename, control_or_treated, input_or_ip, situation(default 0 is CONTROL_SITUATION else are TREATED_SITUATION)
## TREATED_SITUATION_STARTPOINT is the default situation check point
#!/bin/Rscript
library(stringr)
args<-commandArgs(T)
aligner_tools_name <- args[1]
designfile <- args[2]

# setting CONTROL_SITUATION and TREATED_SITUATION 
## default 1 is CONTROL_SITUATION else are TREATED_SITUATION
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))

#Gene expression matrix for control_input_bwa(aligner)
htseq.files <- list.files("./",pattern = str_c(aligner_tools_name,".txt"))
for(i in c(1:max(as.numeric(designtable$situation)))){
  trans.htseq.input.count <- c()
  pc.names <- c()
  pc.samples <- c()
  for(pc in grep(str_c("input_",i,"_",aligner_tools_name),htseq.files,value = TRUE)){
    pc.exp <- read.table(pc,header=F,sep="\t",row.names=1,quote = "")
    trans.htseq.input.count <- cbind(trans.htseq.input.count,pc.exp[,1])
    pc.names <- rownames(pc.exp) #genes name
    pc.samples <- c(pc.samples,pc) #samples name
  }
  colnames(trans.htseq.input.count) <- pc.samples
  rownames(trans.htseq.input.count) <- pc.names  
  trans.htseq.ip.count <- c()
  pc.names <- c()
  pc.samples <- c()
  for(pc in grep(str_c("ip_",i,"_",aligner_tools_name),htseq.files,value = TRUE)){
    pc.exp <- read.table(pc,header=F,sep="\t",row.names=1,quote = "")
    trans.htseq.ip.count <- cbind(trans.htseq.ip.count,pc.exp[,1])
    pc.names <- rownames(pc.exp) #genes name
    pc.samples <- c(pc.samples,pc) #samples name
  }
  colnames(trans.htseq.ip.count) <- pc.samples
  rownames(trans.htseq.ip.count) <- pc.names  
  
  #parsing samplenames
  output_pattern = str_c("htseq_situation_",i,"_",aligner_tools_name)  #添加aligner
  write.table(trans.htseq.input.count, file = str_c(output_pattern,"_input.count") , sep ="\t", row.names =T,col.names =T)
  write.table(trans.htseq.ip.count, file = str_c(output_pattern,"_ip.count") , sep ="\t", row.names =T,col.names =T)
}

