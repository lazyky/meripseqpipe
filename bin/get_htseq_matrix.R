#!/bin/Rscript
## Rscript get_htseq_matrix.R aligner_tools designfile gtf eg. Rscript get_htseq_matrix.R tophat2 designfile_single.txt
## designfile: filename, control_or_treated, input_or_ip, group(default 0 is CONTROL_SITUATION else are TREATED_SITUATION)

library(parallel)
args<-commandArgs(T)
designfile <- args[1]
THREAD_NUM <- as.numeric(args[2])

designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
#Generate gene count matrix
htseq.files <- list.files("./",pattern = ".txt")
mclapply(unique(designtable$Group),function(x){
  group_id <- x
  trans.htseq.input.count <- c()
  pc.names <- c()
  pc.samples <- c()
  for(pc in grep(paste0(".input_",group_id),htseq.files,value = TRUE)){
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
  for(pc in grep(paste0(".ip_",group_id),htseq.files,value = TRUE)){
    pc.exp <- read.table(pc,header=F,sep="\t",row.names=1,quote = "")
    trans.htseq.ip.count <- cbind(trans.htseq.ip.count,pc.exp[,1])
    pc.names <- rownames(pc.exp) #genes name
    pc.samples <- c(pc.samples,pc) #samples name
  }
  colnames(trans.htseq.ip.count) <- pc.samples
  rownames(trans.htseq.ip.count) <- pc.names  
  
  #parsing samplenames
  output_pattern = paste0("htseq_group_",group_id)  #添加aligner
  trans.htseq.input.count <- trans.htseq.input.count[c(-nrow(trans.htseq.input.count):-(nrow(trans.htseq.input.count)-4)),]
  trans.htseq.ip.count <- trans.htseq.ip.count[c(-nrow(trans.htseq.ip.count):-(nrow(trans.htseq.ip.count)-4)),]
  write.table(trans.htseq.input.count, file = paste0(output_pattern,"_input.count") , sep ="\t", row.names =T,col.names =T)
  write.table(trans.htseq.ip.count, file = paste0(output_pattern,"_ip.count") , sep ="\t", row.names =T,col.names =T)
  },
  mc.cores = THREAD_NUM
)