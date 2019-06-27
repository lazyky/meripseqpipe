## generate QC.peaks list
QC.peaks.filelist <- grep("unanno.txt",grep(".anno.txt", list.files(pattern = "group"), value = T), value = T, invert = T)
QC.peaks.list <- NULL
for( file in QC.peaks.filelist ){
  tmp.QC <-  read.table(file , header = F, sep = "\t", quote = "")[,c(1,2,3,15,11,12,13,14,17)]
  colnames(tmp.QC) <- c("Chr","ChrStart","ChrEnd","ID","Gene_symbol","Coding","Location","Count","RNA_type")
  name.QC <- strsplit(file,".anno.txt")[[1]][1]
  QC.peaks.list[[name.QC]] <- tmp.QC
}

qclist=QC.peaks.list
qclist$A=QC.peaks.list$A
qclist$B=QC.peaks.list$B

distribute_df <- NULL
for(i in 1:length(qclist)){
  peak_freq = c(as.numeric(as.vector(qclist[[i]][which(qclist[[i]][,7]=="5UTR"),8])),
                as.numeric(as.vector(qclist[[i]][which(qclist[[i]][,7]=="CDS"),8]))+100,
                as.numeric(as.vector(qclist[[i]][which(qclist[[i]][,7]=="3UTR"),8]))+200)
  freq = data.frame(Freq = peak_freq, group=factor(rep(paste(names(qclist[i])), length(peak_freq))))
  distribute_df = rbind(distribute_df, freq)
}
ggplot(distribute_df, aes(x=Freq, colour = group))+
  geom_line(stat = "density", size=1, adjust = 0.8)+
  scale_x_continuous(breaks = c(50,150,250), labels = c("5'UTR", "CDS", "3'UTR"))+ #axis labels
  labs(y="m6A coding peak density",x="Region of gene")+
  geom_vline(xintercept = c(100,200), linetype = "dashed")+
  theme_bw()+
  theme(panel.grid =element_blank(),#remove grid line
        axis.title.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.title.y = element_text(size = 15, angle = 90, face = "plain", colour = "black"),
        axis.ticks.x = element_blank()) #remove ticks

compare.list <- read.csv("compare.file",header = F,stringsAsFactors = F)
designtable <- read.csv("formatted_designfile.txt", head = TRUE, stringsAsFactors=FALSE, colClasses = c("character"))
design.matrix <- as.matrix(designtable$Group)
rownames(design.matrix) <- designtable$Sample_ID
colnames(design.matrix) <- "Type"
## generate peak Visualization
annotation.file <- list.files(pattern = "merged_peaks.anno.txt")
annotation.info <- read.table(annotation.file, header = F, sep = "\t", stringsAsFactors = F, quote = "")[,c(4,15,11)]
colnames(annotation.info) <- c("PeakRegion","ID","Gene_symbol")
m6a.peaks.file <- list.files(pattern = "merged_peaks.bed")
m6a.peaks.table <- read.table(m6a.peaks.file, header = F, sep = "\t", stringsAsFactors = F, quote = "")
colnames(m6a.peaks.table) <- c("Chr","ChrStart","ChrEnd","PeakRegion","pvalue")
m6a.peaks.table = merge(x = m6a.peaks.table,y = annotation.info,by = "PeakRegion",all.x = TRUE)
m6a.sites.file <- list.files(pattern = "m6A_sites_merged.bed")
m6a.sites.table <- read.table(m6a.sites.file, header = F, sep = "\t", stringsAsFactors = F, quote = "")
colnames(m6a.sites.table) <- c("Chr","ChrStart","ChrEnd","Gene_symbol","ID","Strand","Score","Group","Sequence")
## generate m6A matrix
m6a.matrix <- as.matrix(read.table(file = grep("quantification.matrix",x = list.files(),value = T), header = T, row.names = 1))
m6a.anno.matrix <- as.data.frame(m6a.matrix)
m6a.anno.matrix$PeakRegion <- rownames(m6a.matrix)
m6a.anno.matrix <- merge(x = annotation.info,y = m6a.anno.matrix,by = "PeakRegion",all.y = TRUE)

## generate diffm6A list
diffm6A.filelist <- grep("_diffm6A_",list.files(pattern = ".txt"), value = T)
diffm6A.list <- NULL
diffm6A.anno.list <-NULL
for( compare_str in compare.list ){
  diffm6A.list[[compare_str]] <- read.table(grep(sub("_vs_","_",compare_str), diffm6A.filelist, value = T),header = T,row.names = 1)
  if( diffm6A_mode == "QNB" ){
    colnames(diffm6A.list[[compare_str]]) <- c("p.treated","p.control","log2FC","log2.OR","pvalue","qvalue","padj")
    diffm6A.list[[compare_str]]$PeakRegion <- rownames(diffm6A.list[[compare_str]])
    diffm6A.list[[compare_str]] <- merge(x = annotation.info,y = diffm6A.list[[compare_str]],by = "PeakRegion", all.y = TRUE)
  }else if( diffm6A_mode == "MeTDiff" ){
    
  }else if( diffm6A_mode == "MATK" ){
    diffm6A.list[[compare_str]]$padj = p.adjust(diffm6A.list[[compare_str]]$pvalue, method = "BH")
  }else if( diffm6A_mode == "bedtools" ){
    
  }else{
    stop("Please check your setting of quantification_mode")
  }
  
}
pca_plot <- function(mat,coldt){
  sum = summary(prcomp(mat[,rownames(coldt)]))
  pcaData <- as.data.frame(sum$rotation)
  ggplot(pcaData, aes(PC1, PC2, color=coldt$Type)) +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(pcaData)))+
    xlab(paste("PC1","(",round(100*sum$importance[2,1],1),"%)",sep = "")) +
    ylab(paste("PC2","(",round(100*sum$importance[2,2],1),"%)",sep = "")) +
    scale_colour_hue("Type") +
    #  coord_fixed() +
    theme_bw()
  
}