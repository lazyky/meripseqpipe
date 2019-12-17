#Find m6Am 5'UTR peaks
anno_5UTR=overlap.anno[which(overlap.anno$Gene.site=="5UTR"),]

gtf.temp=fread("/data/database/hg38/GENCODE/gencode.v25.annotation.gtf",sep="\t",skip = 5,data.table = F)
gtf.temp=cbind(gtf.temp,Transcript.id=strsplit2(strsplit2(gtf.temp$V9,split = "transcript_id ")[,2],split = ";")[,1])
gtf.temp$Transcript.id=gsub("\"","",gtf.temp$Transcript.id)

write.table(strsplit2(anno_5UTR$Peak.id,split=":|-"),"UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")
system("fastaFromBed -fi /data/database/hg38/genome.fa -bed UTR5.peak.bed -fo UTR5.peak.fa")
system("/data/software/homer/bin/homer2 find -i UTR5.peak.fa -m /data/xingyang/m6A_zhengjian/BCA.motif -p 5 > /data/xingyang/m6A_zhengjian/analysis/BCA_peak_offset.txt")
BCA_in_5UTR_offset=read.table("analysis/BCA_peak_offset.txt",header=F)

anno_5UTR=overlap.anno[unique(BCA_in_5UTR_offset$V1),]
anno_5UTR=merge(gtf.temp,anno_5UTR,by="Transcript.id",all.y=T)
anno_5UTR=anno_5UTR[which(anno_5UTR$V3=="UTR"),]
anno_5UTR$temp.start=anno_5UTR$V4-anno_5UTR$Start
anno_5UTR$temp.end=anno_5UTR$V5-anno_5UTR$Start
anno_5UTR[which(anno_5UTR$temp.start>0),"temp.start"]=1
anno_5UTR[which(anno_5UTR$temp.start<0),"temp.start"]=(-1)
anno_5UTR[which(anno_5UTR$temp.end>0),"temp.end"]=1
anno_5UTR[which(anno_5UTR$temp.end<0),"temp.end"]=(-1)
anno_5UTR=anno_5UTR[which((anno_5UTR$temp.start*anno_5UTR$temp.end)<=0),]
anno_5UTR.bed=cbind(anno_5UTR$V1,anno_5UTR$V4,anno_5UTR$V5,anno_5UTR$Peak.id,".",anno_5UTR$V7)
write.table(anno_5UTR.bed,"UTR5.peak.bed",row.names = F,col.names = F,quote = F,sep="\t")
system("fastaFromBed -fi /data/database/hg38/genome.fa -bed UTR5.peak.bed -s -name -fo UTR5.peak.fa")

temp.utr5=read.table("UTR5.peak.fa",sep="\n")
temp.utr5=cbind(temp.utr5,substr(temp.utr5[,1],1,1))

i=2
n=nrow(temp.utr5)
temp.utr5=cbind(temp.utr5,type=NA)
while(i<=n){
  if(temp.utr5[i,2]=="A"){
    temp.utr5[c(i-1,i),"type"]="m6Am"
  }
  i=i+2
}
m6Am=na.omit(temp.utr5)
m6Am=m6Am[grep(">",m6Am$V1),]
m6Am=gsub(">","",m6Am$V1)
m6Am=strsplit2(m6Am,split="[(]")[,1]