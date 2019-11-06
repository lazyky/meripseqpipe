#!/bin/Rscript
## Rscript MeTDiff_diffm6A.R <designfile> <compare_str>
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)
library(MeTDiff)
args <- commandArgs(T)
designfile <- args[1]
compare_str <- args[2]

.help.digamma <- function(xx,alpha){
  Tm <- dim(xx)
  TT <- Tm[1]
  m <- Tm[2]
  
  res <- matrix(0,TT,1)
  for (ii in 1:m){
    res <- res + digamma(xx[,ii] + alpha)
  }
  return(res)
}

.help.trigamma <- function(xx,alpha){
  Tm <- dim(xx)
  if (is.null(Tm)){
    TT <- length(xx)
    m <- 1
  }else{
    TT <- Tm[1]
    m <- Tm[2]
  }
  res <- matrix(0,TT,1)
  for (ii in 1:m){
    res <- res + trigamma(xx[,ii] + alpha)
  }
  return(res)
}

.help.postprob <- function(dxx,dyy,dnn,xx,yy,alpha,beta){
  N <- length(alpha)
  Tm <- dim(xx)
  if (is.null(Tm)){
    TT <- length(xx)
    m <- 1
  }else{
    TT <- Tm[1]
    m <- Tm[2]
  }
  res <- matrix(0,TT,N)
  
  for (ii in 1:m){
    dnx <- as.matrix(dxx[,ii])
    dny <- as.matrix(dyy[,ii])
    dn <- as.matrix(dnn[,ii])
    x <- as.matrix(xx[,ii])
    y <- as.matrix(yy[,ii])
    res <- res + (dn-dnx-dny) %*% matrix(1,1,N) + lgamma(x %*% matrix(1,1,N) + matrix(1,TT) %*% alpha) - 
      lgamma(matrix(1,TT) %*% (alpha+beta) + (x+y) %*% matrix(1,1,N)) +
      lgamma(y %*% matrix(1,1,N) + matrix(1,TT) %*% beta) + lgamma(matrix(1,TT) %*% (alpha+beta)) - 
      lgamma(matrix(1,TT) %*% alpha) - lgamma(matrix(1,TT) %*% beta)
  }
  res <- exp(res)
}

.help.factorial <- function(count){
  #compute the log(count!)
  cm = max(count)
  if (is.null(ncol(count))){
    D <- 1
  }else{
    D <- ncol(count)
  }
  if(cm > 50000){
    dnorm <- as.matrix(lgamma(data.matrix(count+1)))
  }
  else{
    tmp  <- cumsum(rbind(0,log(as.matrix(1:max(count)))))
    dnorm <- matrix(tmp[data.matrix(count+1)],ncol=D)
  }
}

.betabinomial.lh <- function(x,y,Nit=40,Npara=1e-9){
  #   x <- as.matrix(x[peak,])
  #   y <- as.matrix(y[peak,])
  N <- 2 # number of states
  J <- matrix(0,N,1)
  H <- matrix(0,N,N)
  T <- nrow(x)
  IP_mean <- rowMeans(x)
  INPUT_mean <- rowMeans(y)
  nip = ncol(x)
  nin = ncol(y)
  # if the dimension for x and y does not match
  if (nip > nin) {
    avg_input <- round(matrix(rep(INPUT_mean,nip-nin),ncol=nip-nin))
    y <- cbind(y,avg_input) 
  }
  else if (nip < nin){
    avg_ip <- matrix(rep(IP_mean,nin-nip),ncol=nin-nip)
    x <- cbind(x,avg_ip) 
  }
  n <- x + y
  m <- ncol(x)
  rr <- x/n
  
  # use another method to initialize
  p1_e <- exp(sum( log(rr) )/(T*m))
  p2_e <- exp(sum( log(1-rr)/(T*m) ))
  alpha <- 1/2 *(1-p2_e)/(1-p1_e-p2_e ) # to avoid 0
  beta <-  1/2 *(1-p1_e)/(1-p1_e-p2_e )
  c = rbind(alpha,beta)
  # add break condition to avoid alpha is na   
  if ( !any(is.finite(beta)) | !is.finite(alpha) | any(beta <= 0) | any(alpha<= 0) ){
    return(list(logl=rnorm(1)*10000,alpha=c(1,1),beta=c(1,1))) 
  }
  for (nit in 1:Nit){
    J[1] <- T*digamma(sum(c))*m - sum( .help.digamma(as.matrix(n),sum(c)) ) + sum( .help.digamma(as.matrix(x),c[1]) ) - T*digamma(c[1])*m
    J[2] <- T*digamma(sum(c))*m - sum( .help.digamma(as.matrix(n),sum(c)) ) + sum( .help.digamma(as.matrix(y),c[2]) ) - T*digamma(c[2])*m 
    H[1,1] <- T*trigamma(sum(c))*m - sum(.help.trigamma(as.matrix(n),sum(c))) + sum(.help.trigamma(as.matrix(x),c[1])) - T*trigamma(c[1])*m
    H[2,2] <- T*trigamma(sum(c))*m - sum(.help.trigamma(as.matrix(n),sum(c))) + sum(.help.trigamma(as.matrix(y),c[2]))  - T*trigamma(c[2])*m    
    H[1,2] <- T*trigamma(sum(c))*m - sum(.help.trigamma(as.matrix(n),sum(c)))
    H[2,1] <- H[1,2]
    eigvalue <- eigen(H)$values
    
    if ( (any(beta < Npara)) | (any(alpha < Npara)) 
         | abs(eigvalue[1]/eigvalue[2]) > 1e12 | abs(eigvalue[1]/eigvalue[2]) < 1e-12
         | any(eigvalue==0) ){   break  }
    
    #     tmp_step <- -solve(H,tol=1e-20) %*% J
    tmp_step <- -solve(H, J) # using newton smoothing
    tmp <- c + tmp_step
    while(any(tmp <= 0)){
      #       warning(sprintf("Could not update the Newton step ...\n"))
      tmp_step <- tmp_step / 20
      tmp <- c + tmp_step
    }
    c <- tmp
    
  }
  #   caculate the likelihood
  alpha <- c[1]
  beta <- c[2]
  dnx <- .help.factorial(x)
  dny <- .help.factorial(y)
  dn <- .help.factorial(n)
  prob <- .help.postprob(dnx,dny,dn,x,y,alpha,beta)
  return(list(logl=sum(log(prob)),alpha=alpha,beta=beta))
  
}

# merge and compare two conditions
diff.call.module <- function(meth1,unmeth1,meth2,unmeth2){
  #x = untreated IP, y = untreated input, xx = treated IP, yy = treated input
  no_peak=length(meth1[,1]) #PEAK$loci2peak_merged[,1])
  pvalues <- rep(1,no_peak)
  log.fc <- rep(0,no_peak)
  for (ipeak in 1:no_peak) {
    if (ipeak%%1000 == 0){print(ipeak)}
    x = t(as.array(meth1[ipeak,]))
    y = t(as.matrix(unmeth1[ipeak,]))
    xx = t(as.matrix(meth2[ipeak,]))
    yy = t(as.matrix(unmeth2[ipeak,]))
    xxx = cbind(x,xx)
    yyy = cbind(y,yy)
    #BBtest
    logl1 <- .betabinomial.lh(x,y+1)
    logl2 <- .betabinomial.lh(xx,yy+1)
    logl3 <- .betabinomial.lh(xxx,yyy+1)
    tst <- (logl1$logl+logl2$logl-logl3$logl)*2
    pvalues[ipeak] <- 1 - pchisq(tst,2)
    log.fc[ipeak] <- log2( (sum(xx)+1)/(1+sum(yy)) * (1+sum(y))/(1+sum(x)) ) 
    
  }
  p <- pvalues
  fdr <- p.adjust(pvalues,method='fdr')

  DIFF <- list(fdr=fdr,pvalues=p,fc=log.fc)
  # result
  result =list()
  result$DIFF = DIFF
  return(result)
  
}

designtable <- read.csv(designfile, head = TRUE, stringsAsFactors=FALSE, colClasses = c("character"), check.names=F)
design.matrix <- as.data.frame(designtable$Group)
rownames(design.matrix) <- designtable$Sample_ID
colnames(design.matrix) <- "Condition"

# Get the information of groups from compare_str
if(length(unique(design.matrix$Condition)) < 2){
  stop( "The count of Group is less than two, please check your designfile.")
}else if( compare_str == "two_group" ){
  # Get the information without compare_str beacause of only two groups
  group_id_1 <- unique(design.matrix$Condition)[1]
  group_id_2 <- unique(design.matrix$Condition)[2]
}else{
  # Running MeTDiff quantification with compare_str
  group_id_1 <- strsplit(as.character(compare_str), "_vs_")[[1]][1]
  group_id_2 <- strsplit(as.character(compare_str), "_vs_")[[1]][2]
}
design.matrix <- subset(design.matrix, Condition == group_id_1 | Condition == group_id_2 )
design.matrix$Condition <- factor(design.matrix$Condition,labels = c("control","treatment"))
filelist = list.files(path = ".",pattern = ".count",full.names = T)
## Generate the matrix of peaks count
rpkm_peaks_list <- NULL
for(sample_id in row.names(design.matrix)){
  input_count_file <- grep(paste0("[.]",sample_id,"[.]input"),filelist,value = TRUE)
  input_count_table <- read.table(file = input_count_file, sep = "\t", row.names = NULL,header = T)

  ip_count_file <- grep(paste0("[.]",sample_id,"[.]ip"),filelist,value = TRUE)
  ip_count_table <- read.table(file = ip_count_file, sep = "\t", row.names = NULL, header = T)
  rpkm <- cbind(input_count_table[,5],ip_count_table[,5])
  colnames(rpkm) <- c(paste0(sample_id,".input"),paste0(sample_id,".ip"))
  rpkm_peaks_list <- cbind(rpkm_peaks_list,rpkm)
}
rownames(rpkm_peaks_list) <- ip_count_table$PeakName

## generate design matrix
design.matrix$m6A <- "input"
design.matrix$sample_id <- paste0(rownames(design.matrix),".input")
design.matrix_ip <- design.matrix
design.matrix_ip$m6A <- "IP"
design.matrix_ip$sample_id <- paste0(rownames(design.matrix_ip),".ip")
design.matrix <- rbind(design.matrix,design.matrix_ip)
rownames(design.matrix) <- design.matrix$sample_id
design.matrix$m6A <- factor(design.matrix$m6A)
design.matrix <- design.matrix[colnames(rpkm_peaks_list),]

cnts <- rpkm_peaks_list
meta <- design.matrix
run.metdiff <- function(cnts,meta){
  meth1 <- cnts[,which(meta$Condition == 'treatment' & meta$m6A == "IP")]
  meth2 <- cnts[,which(meta$Condition != 'treatment' & meta$m6A == "IP")]
  unmeth1 <- cnts[,which(meta$Condition == 'treatment' & meta$m6A == "input")]
  unmeth2 <- cnts[,which(meta$Condition != 'treatment' & meta$m6A == "input")]
  metdiff.result <- diff.call.module(meth1,unmeth1,meth2,unmeth2)
  results <- data.frame(log2FC= metdiff.result$DIFF$fc, pvalue = metdiff.result$DIFF$pvalues, padj = p.adjust(metdiff.result$DIFF$pvalues,"BH"))
  rownames(results) <- rownames(cnts)
  return(results)
}
results <- run.metdiff(rpkm_peaks_list,design.matrix)
write.table(results,file = paste0("MeTDiff_diffm6A_",group_id_1, "_",group_id_2,".txt") ,sep = "\t",quote = F)
