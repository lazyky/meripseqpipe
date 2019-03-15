args <- commandArgs(trailingOnly = TRUE)

m6A_anno=read.table(args[1],header=F,sep="\t",quote="")
m6A_unanno=read.table(args[2],header=F,sep="\t",quote="")

peak_freq=c(as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="5UTR"),14])),as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="CDS"),14]))+100,as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="3UTR"),14]))+200)
pdf(paste(args[3],".m6A.coding_topology.pdf",sep=""))
#plot(density(peak_freq),xlim=c(0,300),lwd=3,ylab="m6A peak density",xlab=NA,main=NA,col="darkblue",xaxt="n")
plot(table(peak_freq)/length(peak_freq),xlim=c(0,300),lwd=3,ylab="m6A coding peak density",xlab=NA,main=NA,col="darkblue",xaxt="n",type="l") #xingyang
abline(v=100,lty=2,col="black",lwd=1)
abline(v=200,lty=2,col="black",lwd=1)
text(x=c(50,150,250),y=-0.001,pos=2,labels=c("5' UTR","CDS","3' UTR"),xpd=TRUE,font=2)
dev.off()

#noncoding gene
exon_peak_freq=c(as.numeric(as.vector(m6A_anno[which(m6A_anno[,13]=="exon"),14])))
pdf(paste(args[3],".m6A.non-coding_topology.pdf",sep=""))
plot(table(exon_peak_freq)/length(exon_peak_freq),xlim=c(0,100),lwd=3,ylab="m6A non-coding peak density",xlab=NA,main=NA,col="darkblue",xaxt="n",type="l")
text(x=c(50),y=-0.001,pos=2,labels=c("exon"),xpd=TRUE,font=2)
dev.off()

color_index <- c('Protein-coding', 'Long non-coding RNA', 'Small non-coding RNA', 'Pseudogenes', 'Immunoglobulin/T-cell receptor gene segments', 'Unknown')
color_name  <- c('red', 'blue', 'green', 'yellow', 'pink', 'grey')

first.level.stat.temp <- as.data.frame(table(m6A_anno[,17]))
first.level.stat <- rbind (first.level.stat.temp, data.frame(Var1="Unknown", Freq=nrow(m6A_unanno)))
pdf(paste(args[3],".m6A.location.barplot.1st.pdf",sep=""))
type_number<-length(first.level.stat[,1])
def.par <- par()
layout(matrix(c(1,2),nr=1),widths=c(20,15))
op1 <- par(mar=c(4.1,5.1,4.1,0))
a<-barplot(first.level.stat[,2], ylab="peak number",col=color_name[match(first.level.stat[,1],color_index)] )
box()
axis(side=1,at=a,labels=1:type_number,padj=-2,tick=FALSE)
op1 <- par(mar=c(0,0,1.1,1.1))
plot(0,0,pch="",xlim=c(0,1),ylim=c(-3,length(first.level.stat[,1])+4),bty="n",axes=F,xlab="",ylab="")
y_offset <- 0;
if(type_number >= 25) { y_offset <- -1; }
if(type_number == 24) { y_offset <- 0.8; }
if(type_number <= 23) { y_offset = 2; }
for(i in type_number:1) {text(0.01, type_number-i+1+y_offset, paste(i,": ", first.level.stat[i,1],"  ","(", as.character(first.level.stat[i,2]),")", sep=""),pos=4,cex=0.6,col=color_name[which(color_index == first.level.stat[i,1])])}
par(op1)
#par(def.par)
title("Gene stat (first level)")
dev.off()
write.table(first.level.stat,"m6A.location.1nd.data.txt")

for (i in 1:length(m6A_anno[,1])) {m6A_anno[i,18] <- paste(m6A_anno[i,17],m6A_anno[i,16],sep="--")}
second.level.stat = as.data.frame(table(m6A_anno[,18]))
for (i in 1:length(second.level.stat[,1])) {gene_1st_type = strsplit(toString(second.level.stat[i,1]),split="--")[[1]][1]; second.level.stat[i,3]=color_name[which(color_index == gene_1st_type)]}
pdf(paste(args[3],".m6A.location.barplot.2nd.pdf",sep=""))
type_number<-length(second.level.stat[,1])
def.par <- par()
layout(matrix(c(1,2),nr=1),widths=c(20,25))
op1 <- par(mar=c(4.1,5.1,4.1,0))
a<-barplot(second.level.stat[,2], ylab="peak number",col=second.level.stat[,3])
box()
axis(side=1,at=a,labels=1:type_number,padj=-2,tick=FALSE)
op1 <- par(mar=c(0,0,1.1,1.1))
plot(0,0,pch="",xlim=c(0,1),ylim=c(-3,length(second.level.stat[,1])+4),bty="n",axes=F,xlab="",ylab="")
y_offset <- 0;
if(type_number >= 25) { y_offset <- -1; }
if(type_number == 24) { y_offset <- 0.8; }
if(type_number <= 23) { y_offset = 2; }
for(i in type_number:1) {text(0.01, type_number-i+1+y_offset, paste(i,": ", second.level.stat[i,1],"  ","(", as.character(second.level.stat[i,2]),")", sep=""),pos=4,cex=0.6,col=second.level.stat[i,3])}
par(op1)
#par(def.par)
title("Gene stat (second level)",adj=0)
dev.off()
write.table(second.level.stat,"m6A.location.2nd.data.txt")
