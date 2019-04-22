## Rscript arranged_results.R aligners_tools_name peak_calling_tools_name 
library("QNB")
args <- commandArgs(T)
designfile <- args[1]
compare_str <- as.character(args[2])

designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
filelist =list.files(path = "./",pattern = ".count")
countlist <- NULL
#combine the matrix by groups
for(group_id in designtable$Group){
  input.count <- c()
  input.names <- c()
  input.samples <- c()
  for(input in grep(group_id, grep(".input",filelist,value = TRUE), value = T)){
    input.exp <- read.table(input,header=T,sep="\t",row.names= NULL,quote = "")
    input.count <- cbind(input.count,input.exp[,5])
    input.names <- input.exp[,4] #peaks name
    input.samples <- c(input.samples,input) #samples name
  }
  colnames(input.count) <- input.samples
  rownames(input.count) <- input.names
  countlist[[paste0(group_id,"_input")]] <- input.count 
  
  ip.count <- c()
  ip.names <- c()
  ip.samples <- c()
  for(ip in grep(group_id, grep(".ip",filelist,value = TRUE), value = T)){
    ip.exp <- read.table(ip,header=T,sep="\t",row.names= NULL,quote = "")
    ip.count <- cbind(ip.count,ip.exp[,5])
    ip.names <- ip.exp[,4] #peaks name
    ip.samples <- c(ip.samples,ip) #samples name
  }
  colnames(ip.count) <- ip.samples
  rownames(ip.count) <- ip.names
  countlist[[paste0(group_id,"_ip")]] <- ip.count
}

# Running QNB quantification
if(length(unique(designtable$Group)) < 2){
  stop( "The count of Group is less than two, please check your designfile.")
}else if( compare_str == "two_group" ){
  # Running QNB quantification without compare_str beacause of only two groups
  group_id_1 <- unique(designtable$Group)[1]
  group_id_2 <- unique(designtable$Group)[2]
}else{
  # Running QNB quantification with compare_str
  group_id_1 <- strsplit(as.character(compare_str), "_vs_")[[1]][1]
  group_id_2 <- strsplit(as.character(compare_str), "_vs_")[[1]][2]
}
meth1 = countlist[[paste0(group_id_1,"_ip")]]
meth2 = countlist[[paste0(group_id_2,"_ip")]]
unmeth1 = countlist[[paste0(group_id_1,"_input")]]
unmeth2 = countlist[[paste0(group_id_2,"_input")]]
output_name <- paste0("QNB_diffm6A_",group_id_1, "_",group_id_2)
dir.create(output_name)
result <- qnbtest(meth1, meth2, unmeth1, unmeth2, mode="auto", output.dir = output_name)
write.table(result,file = paste0(output_name,".txt"))