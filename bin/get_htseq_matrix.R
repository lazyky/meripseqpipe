#!/bin/Rscript
## Rscript get_htseq_matrix.R designfile THREAD_NUM eg. Rscript get_htseq_matrix.R designfile_single.txt 10
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
  for(pc in grep(paste0(".input_",group_id,"[.]bam"),htseq.files,value = TRUE)){
    pc.exp <- read.table(pc,header=F,sep="\t",row.names=1,quote = "")
    trans.htseq.input.count <- cbind(trans.htseq.input.count,pc.exp[,1])
    pc.names <- rownames(pc.exp) #genes name
    pc.samples <- c(pc.samples,pc) #samples name
  }
  rownames(trans.htseq.input.count) <- pc.names  
  trans.htseq.input.count <- as.matrix(trans.htseq.input.count[c(-nrow(trans.htseq.input.count):-(nrow(trans.htseq.input.count)-4)),])
  colnames(trans.htseq.input.count) <- pc.samples
  #parsing samplenames
  output_pattern = paste0("htseq_group_",group_id)  #添加aligner
  write.table(trans.htseq.input.count, file = paste0(output_pattern,"_input.count") , sep ="\t", row.names =T,col.names =T)
},
mc.cores = THREAD_NUM
)
htseq.filelist = grep("htseq",list.files(path = "./",pattern = "input.count"), value = T)
expression.matrix <- NULL
for( file in htseq.filelist ){
  tmp.expression.table <- as.matrix(read.table(file, header = TRUE, row.names = 1, check.names=F))
  expression.matrix <- cbind(expression.matrix, tmp.expression.table)
}
expression.matrix <- expression.matrix[c(-nrow(expression.matrix):-(nrow(expression.matrix)-4)),]
colnames(expression.matrix) <- as.matrix(lapply(strsplit(colnames(expression.matrix),".input"), function(x){ x[1]}))
write.table(expression.matrix,file = "expression.matrix",quote=F)