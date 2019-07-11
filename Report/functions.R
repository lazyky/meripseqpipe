library(ggplot2)
library(ggseqlogo)
#pie plot
pie_plot <- function(dt){
  data1 <- table(dt[,9])
  data2 <- as.data.frame(data1)
  names(data2)[1] <- "type"
  data2 = data2[order(data2$Freq, decreasing = TRUE),]
  data2$type = factor(data2$type, levels = data2$type, ordered = T)
  myLabel = paste(data2$type, " ( ", round(data2$Freq/sum(data2$Freq)*100,2), "% )  ", sep = "")  
  p = ggplot(data2, aes(x = "", y = Freq, fill = type)) + 
    geom_bar(stat = "identity", width = 1) + 
    scale_fill_brewer(palette ="Set3",direction = 1,breaks = data2$type, labels = myLabel)+ 
    #scale_fill_discrete(breaks = data2$type, labels = myLabel)+
    coord_polar(theta = "y", direction = 1) + 
    annotate("text",x=-Inf,y=Inf,vjust=1.5,hjust=-.12,label="a")+ 
    labs(x = "", y = "", fill="Type") +
    theme_bw()+
    #geom_text(aes(x=1.2,y=sum(data2$Freq)-cumsum(data2$Freq)+data2$Freq/2,label=paste( round(data2$Freq/sum(data2$Freq)*100,2),"%",sep = "")))+
    theme(axis.ticks = element_blank()) + 
    theme(legend.title = element_blank(),legend.text = element_text(size = 15),
          axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank())#,legend.position = "bottom")
  return(p)
}

#multidistribution plot
distribution_plot <- function(qclist){
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
}

#motif
motif_plot <- function(motif, pval, rank){
  ggplot()+
    geom_logo(motif, method = "bits")+
    annotate("text", x=ncol(motif)-0.5, y=2.45, label=paste0("p = ",pval),size = 8)+
    ggtitle(rank)+
    theme(plot.title = element_text(hjust = 0, size = 8))+
    theme_logo()
}

#stack bar plot
stack_bar_plot <- function(dat){
  ggplot(dat, aes(x = sample,y = Num, fill = type))+
    ####position="stack"堆叠状
    geom_bar(stat ="identity",width = 0.2,position ="stack",show.legend = F)+
    geom_point(aes(x=sample, y=Num, color = type),alpha=0)+
    labs(x = "",y = "Number of Tags", title = "Read Distribution")+    
    guides(fill = guide_legend(reverse = T))+
    guides(color = guide_legend(reverse = T, override.aes = list(alpha=1, size=3)))+
    #scale_fill_brewer(palette ="Set2")+
    #scale_color_brewer(palette ="Set2")+
    theme_bw()+
    theme(plot.title = element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),  
          legend.title = element_blank(),              
          legend.text = element_text(face = "bold",margin = margin(r = 30, unit = "pt")),     
          legend.position = "bottom",
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line=element_line(colour="black"))+
    coord_flip()
}

fill_bar_plot <- function(dat){
  
  ggplot(dat, aes(x = sample,y = Pct, fill = type))+
    geom_bar(stat ="identity",width = 0.5,position = "fill")+
    #geom_point(aes(x = point2, y= point, color = type),alpha=0)+
    labs(x = "",y = "Number of Tags", title = "Read Distribution")+    
    guides(fill = guide_legend(reverse = T, nrow = 2, byrow = T,
                               override.aes = list(alpha=1, size=3)))+
    scale_fill_manual(values = c("#966279","#91E8E1","#F45B5B","#2B908F","#E4D354",
                                 "#F15C80","#8085E9","#F7A35C","#90ED7D","#434348",
                                 "#7CB5EC"))+
    #guides(color = guide_legend(reverse = T, override.aes = list(alpha=1, size=3)))+
    #scale_fill_brewer(palette ="Set2")+
    #scale_color_brewer(palette ="Set2")+
    theme_bw()+
    theme(plot.title = element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),  
          legend.title = element_blank(),              
          legend.text = element_text(face = "bold",margin = margin(r = 30, unit = "pt")),     
          legend.position = "bottom",
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line=element_line(colour="black"))+
    coord_flip()
}

#volcano plots
volcano_plot_dm = function(res, Sample_1 = "A", Sample_2 = "B", lfc = 0, pval = 0.05, groupname = ""){
  par(mar = c(5, 6, 5, 5))
  tab = data.frame(logFC = res$log2FC, negLogPval = -log10(res$pvalue)) 
  tab$gene_name = rownames(res)
  tab = na.omit(tab)
  tab<-tab%>%mutate(threshold = ifelse(logFC >= 0.58 & negLogPval > -log10(0.05) ,"B", ifelse(logFC<=-0.58 & negLogPval > -log10(0.05), "A", "C")))
  n_up = length(which(tab$threshold=="B"))
  n_down = length(which(tab$threshold=="A"))
  tab_order = tab[order(tab$negLogPval, decreasing = T),]
  ggplot(tab_order, aes(x=logFC, y=negLogPval)) +
    geom_point(aes(colour = threshold)) +
    scale_colour_manual(values = c("A"= "#619cff", "B"="#f8766d",  "C"= "#c8c8c8"),
                        labels=c(paste("Down: ", n_down, sep=""),paste("Up: ", n_up, sep = "") , "No sig"), name = NULL) +
    geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed") +
    geom_vline(aes(xintercept=-0.58), linetype="dashed") +
    geom_vline(aes(xintercept=0.58), linetype="dashed") +
    ggtitle(paste("Volcano Plot of Different Methylation in", groupname))+
    xlab(expression(paste(Log[2], " fold change", sep = ""))) +
    ylab(expression(paste(-Log[10], " adjusted P value", sep = ""))) +
    theme_bw() +
    theme(legend.position = 'top',
          plot.title = element_text(hjust = 0.5))
}

volcano_plot_de = function(res, Sample_1 = "A", Sample_2 = "B", lfc = 0, pval = 0.05, groupname = ""){
  par(mar = c(5, 6, 5, 5))
  tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj)) 
  tab$gene_name = rownames(res)
  tab = na.omit(tab)
  tab<-tab%>%mutate(threshold = ifelse(logFC >= 0.58 & negLogPval > -log10(0.05) ,"B", ifelse(logFC<=-0.58 & negLogPval > -log10(0.05), "A", "C")))
  n_up = length(which(tab$threshold=="B"))
  n_down = length(which(tab$threshold=="A"))
  tab_order = tab[order(tab$negLogPval, decreasing = T),]
  ggplot(tab_order, aes(x=logFC, y=negLogPval)) +
    geom_point(aes(colour = threshold)) +
    scale_colour_manual(values = c("A"= "#619cff", "B"="#f8766d",  "C"= "#c8c8c8"),
                        labels=c(paste("Down: ", n_down, sep=""),paste("Up: ", n_up, sep = "") , "No sig"), name = NULL) +
    geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed") +
    geom_vline(aes(xintercept=-0.58), linetype="dashed") +
    geom_vline(aes(xintercept=0.58), linetype="dashed") +
    ggtitle(paste("Volcano Plot of Different Expression in", groupname))+
    xlab(expression(paste(Log[2], " fold change", sep = ""))) +
    ylab(expression(paste(-Log[10], " adjusted P value", sep = ""))) +
    theme_bw() +
    theme(legend.position = 'top',
          plot.title = element_text(hjust = 0.5))
}

#heatmaps
heatmap_dm <- function(mat,coldt){
pheatmap(mat, cluster_rows=FALSE, show_rownames=F, cluster_cols=FALSE, annotation_col=coldt, 
         main = "Heatmap of Different Methylation", scale = "row")
}

heatmap_de <- function(mat, coldt){
  anno_color = c("#e34a33", "#2ca25f")
  names(anno_color) = c(levels(as.factor(coldt$Type)))
  anno_colors = list(Type = anno_color)
  pheatmap(mat, cluster_rows=FALSE, show_rownames=F, cluster_cols=FALSE, annotation_col=coldt, 
           color = colorRampPalette(c(rep('#1C2B6F',1),'black', rep('#E31E26',1)))(50), annotation_colors = anno_colors,
           main = "Heatmap of Different Expression", scale = "row")
}
#PCA
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

#multiplot
ggplot2.multiplot <- function(..., plotlist=NULL, cols=2) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    grid::viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}






