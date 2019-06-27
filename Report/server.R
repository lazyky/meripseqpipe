library(shiny)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(dplyr)
library(stringr)
library(knitr)
library(reshape2)
source("functions.R")
options(stringsAsFactors=F)


server <- function(input, output, session) {

#qc peaks-----------------------------------------  
  
  peak_dt <- eventReactive(input$peakbutton, {
    tabledt = QC.peaks.list[[which(names(QC.peaks.list)==input$peak)]]
    motif = t(apply(QC.motif.list[[which(names(QC.motif.list)==input$peak)]], 1, function(x)as.numeric(x)))
    pvalue = as.character(QC.motif.pvalue[1,which(names(QC.motif.list)==input$peak)])
    return(list(tabledt = tabledt, motif = motif, p = pvalue, list = QC.peaks.list))
  })
  
  output$peaktable <- renderDataTable({
    peak_dt()$tabledt
  })
  
  output$Pieplot <- renderPlot({
    pie_plot(peak_dt()$tabledt)
  })
  
  output$motif <- renderPlot({
    motif_plot(peak_dt()$motif, peak_dt()$p)
  })
  
  output$peakdistri <- renderPlot({
    distribution_plot(peak_dt()$list)
  })

#qc reads--------------------------------------------
  
  reads_dt <- eventReactive(input$readsbutton, {
    bardt0 = QC.reads.distribution
    bardt0 = cbind(sample = rownames(bardt0), bardt0)
    colnames(bardt0) = sapply(colnames(bardt0), function(x)gsub("_tag_pct","",x), USE.NAMES = F)
    bardt = melt(bardt0, id.vars = "sample", variable.name = "type", value.name = "Pct")
    bardt$type = factor(bardt$type, levels = c("other_intergenic","tes_down_10kb","tes_down_5kb","tes_down_1kb",
                                               "tss_up_10kb","tss_up_5kb","tss_up_1kb","introns","3_utr_exons",
                                               "5_utr_exons","cds_exons"))
    bardt$Pct = as.numeric(bardt$Pct)
    read_stat = QC.reads.stats
    return(list(bardt = bardt, read_stat = read_stat))
  })
  
  output$readsbarplot <- renderPlot({
    fill_bar_plot(reads_dt()$bardt)
  })
  
  output$readstable <- renderDataTable({
    reads_dt()$read_stat
  })
  
#different methylation-------------------------------
  
  dm_dt <- eventReactive(input$dmbutton,{
    matrix = m6a.anno.matrix
    dmres = diffm6A.list[[which(names(diffm6A.list)==input$dmgroup)]]
    group1 = strsplit(input$dmgroup, "_vs_")[[1]][1]
    group2 = strsplit(input$dmgroup, "_vs_")[[1]][2]
    coldata0 = as.data.frame(design.matrix)
    coldata = subset(coldata0, Type==group1|Type==group2)
    coldata$Type = as.factor(coldata$Type)
    dmg = subset(dmres, abs(log2FC)> 0.58 & pvalue < 0.05)
    matrix1 = matrix[,-c(2:3)]
    rownames(matrix1) = matrix1$PeakRegion
    matrix1 = matrix1[,-1]
    dm_mat = matrix1[row.names(dmg),rownames(coldata)]
    select <- dmg[order(dmg$log2FC, decreasing = TRUE), ] 
    dm_mat = log2(dm_mat+1)
    dm_mat = dm_mat[rownames(select),]
    dm_mat = na.omit(dm_mat)
    return(list(mat = matrix, res = dmres, dm_mat = dm_mat, coldata = coldata, name = group2))
  })
  
  output$dmmatrix <- renderDataTable({
    dm_dt()$mat
  })
  
  output$dmres <- renderDataTable({
    dm_dt()$res
  })
  
  output$dmvolcano <- renderPlot({
    volcano_plot_dm(res = dm_dt()$res, groupname = dm_dt()$name)
  })
 
  output$dmheatmap <- renderPlot({
    heatmap_dm(mat = dm_dt()$dm_mat, coldt = dm_dt()$coldata)
  })
  
  output$dmPCA <- renderPlot({
    pca_plot(mat = dm_dt()$dm_mat, coldt = dm_dt()$coldata)
  })
  
#different expression--------------------------------------------
  
  de_dt <- eventReactive(input$debutton,{
    matrix = expression.matrix
    deres = diffexpression.list[[which(names(diffexpression.list)==input$degroup)]]
    group1 = strsplit(input$degroup, "_vs_")[[1]][1]
    group2 = strsplit(input$degroup, "_vs_")[[1]][2]
    coldata0 = as.data.frame(design.matrix)
    coldata = subset(coldata0, Type==group1|Type==group2)
    coldata$Type = as.factor(coldata$Type)
    deg = subset(deres, abs(log2FoldChange)> 0.58 & padj < 0.05)
    rownames(deg) = deg$ID
    de_mat = matrix[row.names(deg),rownames(coldata)]
    select <- deg[order(deg$log2FoldChange, decreasing = TRUE), ] 
    de_mat = log2(de_mat+1)
    de_mat = de_mat[rownames(select),]
    de_mat = na.omit(de_mat)
    groupname = input$degroup
    return(list(mat = matrix, res = deres, de_mat = de_mat, coldata = coldata, name = group2))
  })
  
  output$dematrix <- renderDataTable({
    de_dt()$mat
  })
  
  output$deres <- renderDataTable({
    de_dt()$res
  })
  
  output$devolcano <- renderPlot({
    volcano_plot_de(res = de_dt()$res, groupname = de_dt()$name)
  })
  
  output$deheatmap <- renderPlot({
    heatmap_de(mat = de_dt()$de_mat, coldt = de_dt()$coldata)
  })
  
  output$dePCA <- renderPlot({
    pca_plot(mat = de_dt()$de_mat, coldt = de_dt()$coldata)
  })
  
  
  
  
  
  
  
  
}