#Zhao Qi 
library(stringr)
#setwd(dir = "E:/zky/m6Apipe/bin/htseq_count/") #test
args<-commandArgs(T) # Rscript get_htseq_matrix.R aligner_tools eg. bwa

#Gene expression matrix for control_input_bwa(aligner)

filename_pattern = str_c(str_c("*control_input","1",args[1],sep="_"),".txt")
htseq.files <- list.files("./",pattern = filename_pattern)
trans.htseq.count <- c()
pc.names <- c()
for(pc in htseq.files){
  pc.exp <- read.table(pc,header=F,sep="\t",row.names=1,quote = "")
  trans.htseq.count <- cbind(trans.htseq.count,pc.exp[,1])
  pc.names <- rownames(pc.exp) #基因名
}
#parsing samplenames
output_pattern = str_c("control_input",args[1],sep = "_")  #添加aligner
pc.samples <- str_c(unlist(lapply(htseq.files,function(x){unlist(strsplit(x,"_"))[1]})),"_",output_pattern) #添加列的samplename
colnames(trans.htseq.count) <- pc.samples
rownames(trans.htseq.count) <- pc.names
write.table(trans.htseq.count, file = str_c(output_pattern,".count") , sep ="\t", row.names =T,col.names =T)


#Gene expression matrix for control_ip_bwa(aligner)

filename_pattern = str_c(str_c("*control_ip","1",args[1],sep="_"),".txt")
htseq.files <- list.files("./",pattern = filename_pattern)
trans.htseq.count <- c()
pc.names <- c()
for(pc in htseq.files){
  pc.exp <- read.table(pc,header=F,sep="\t",row.names=1,quote = "")
  trans.htseq.count <- cbind(trans.htseq.count,pc.exp[,1])
  pc.names <- rownames(pc.exp)
}
#parsing samplenames
output_pattern = str_c("control_ip",args[1],sep = "_")
pc.samples <- str_c(unlist(lapply(htseq.files,function(x){unlist(strsplit(x,"_"))[1]})),"_",output_pattern)
colnames(trans.htseq.count) <- pc.samples
rownames(trans.htseq.count) <- pc.names
write.table(trans.htseq.count, file = str_c(output_pattern,".count") , sep ="\t", row.names =T,col.names =T)


#Gene expression matrix for treated_input_bwa(aligner)

filename_pattern = str_c(str_c("*treated_input","1",args[1],sep="_"),".txt")
htseq.files <- list.files("./",pattern = filename_pattern)
trans.htseq.count <- c()
pc.names <- c()
for(pc in htseq.files){
  pc.exp <- read.table(pc,header=F,sep="\t",row.names=1,quote = "")
  trans.htseq.count <- cbind(trans.htseq.count,pc.exp[,1])
  pc.names <- rownames(pc.exp)
}
#parsing samplenames
output_pattern = str_c("treated_input",args[1],sep = "_")
pc.samples <- str_c(unlist(lapply(htseq.files,function(x){unlist(strsplit(x,"_"))[1]})),"_",output_pattern)
colnames(trans.htseq.count) <- pc.samples
rownames(trans.htseq.count) <- pc.names
write.table(trans.htseq.count, file = str_c(output_pattern,".count") , sep ="\t", row.names =T,col.names =T)


#Gene expression matrix for treated_ip_bwa(aligner)

filename_pattern = str_c(str_c("*treated_ip","1",args[1],sep="_"),".txt")
htseq.files <- list.files("./",pattern = filename_pattern)
trans.htseq.count <- c()
pc.names <- c()
for(pc in htseq.files){
  pc.exp <- read.table(pc,header=F,sep="\t",row.names=1,quote = "")
  trans.htseq.count <- cbind(trans.htseq.count,pc.exp[,1])
  pc.names <- rownames(pc.exp)
}
#parsing samplenames
output_pattern = str_c("treated_ip",args[1],sep = "_")
pc.samples <- str_c(unlist(lapply(htseq.files,function(x){unlist(strsplit(x,"_"))[1]})),"_",output_pattern)
colnames(trans.htseq.count) <- pc.samples
rownames(trans.htseq.count) <- pc.names
write.table(trans.htseq.count, file = str_c(output_pattern,".count") , sep ="\t", row.names =T,col.names =T)

