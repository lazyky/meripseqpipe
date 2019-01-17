## The function is currently defined as load library and specify the parameters

library("QNB")
library(stringr)
args<-commandArgs(T) # Rscript QNB.R aligner_tools eg. bwa
output_name <- str_c("treated_vs_control_",args[1],"_QNB.csv")

control_ip = read.table(str_c("control_ip_",args[1],".count"),header=TRUE,row.names=1)
treated_ip = read.table(str_c("treated_ip_",args[1],".count"),header=TRUE,row.names=1)
control_input = read.table(str_c("control_input_",args[1],".count"),header=TRUE,row.names=1)
treated_input = read.table(str_c("treated_input_",args[1],".count"),header=TRUE,row.names=1)

# When there are replicates under two conditions, we could select 
# mode="per-condition" or mode="pooled" to estimate the dispersion. 
# The default is mode="auto".

result = qnbtest(control_ip, treated_ip,control_input,treated_input,mode="auto")
write.csv(result, file = output_name )
#dispersion.pdf
#dif_meth.xls
## End(Not run)