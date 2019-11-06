library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
library(shinyBS)
library(shinycssloaders)
library(DT)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(pheatmap)
library(dplyr)
library(stringr)
library(knitr)
library(reshape2)
library(ggseqlogo)
options(stringsAsFactors=F)

server <- function(input, output, session) {

#qc peaks-----------------------------------------  
  
  peak_dt <- eventReactive(input$peakbutton, {
    tabledt = QC.peaks.list[[which(names(QC.peaks.list)==input$peak)]]
    motif1 = t(apply(QC.motif.list[[which(names(QC.motif.list)==input$peak)]][[1]], 1, function(x)as.numeric(x)))
    motif2 = t(apply(QC.motif.list[[which(names(QC.motif.list)==input$peak)]][[2]], 1, function(x)as.numeric(x)))
    motif3 = t(apply(QC.motif.list[[which(names(QC.motif.list)==input$peak)]][[3]], 1, function(x)as.numeric(x)))
    pvalue1 = QC.motif.pvalue[[which(names(QC.motif.list)==input$peak)]][[1]]
    pvalue2 = QC.motif.pvalue[[which(names(QC.motif.list)==input$peak)]][[2]]
    pvalue3 = QC.motif.pvalue[[which(names(QC.motif.list)==input$peak)]][[3]]
    motif_p1 = motif_plot(motif1,pvalue1, "1")
    motif_p2 = motif_plot(motif2,pvalue2, "2")
    motif_p3 = motif_plot(motif3,pvalue3, "3")
    motif_all = ggplot2.multiplot(motif_p1,motif_p2,motif_p3,cols = 1)
    return(list(tabledt = tabledt, motif1 = motif_p1, motif2 = motif_p2, motif3 = motif_p3,list = QC.peaks.list))
  })
  
  output$peaktable <- renderDataTable({
    peak_dt()$tabledt
  })
  
  output$Pieplot <- renderPlot({
    pie_plot(peak_dt()$tabledt)
  })
  
  output$motif <- renderPlot({
    ggplot2.multiplot(peak_dt()$motif1, peak_dt()$motif2, peak_dt()$motif3, cols = 1)
  })
  
  output$peakdistri <- renderPlot({
    distribution_plot(peak_dt()$list)
  })
  
  #download
  output$peak_download_table <- downloadHandler(
    
    filename = function(){
      paste("Peaks table", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file){
      
      write.csv(peak_dt()$tabledt, file, row.names = F)
      
    }
  )
  
  output$peak_download_MOTIF <- downloadHandler(
    
    filename = function(){
      paste("Motif", Sys.time(), ".", input$peak_download_MOTIFtype, sep = "")
    },
    
    content = function(file){
      if(input$peak_download_MOTIFtype == "png")
        png(file=file, height = 1200, width = 740, res = 90)
      else
        pdf(file=file, height = 11, width = 6)
      
      print(ggplot2.multiplot(peak_dt()$motif1, peak_dt()$motif2, peak_dt()$motif3, cols = 1))
      
      dev.off()
    }
  )
  
  output$peak_download_PIE <- downloadHandler(
    
    filename = function(){
      paste("Pie Plot", Sys.time(), ".", input$peak_download_PIEtype, sep = "")
    },
    
    content = function(file){
      if(input$peak_download_PIEtype == "png")
        png(file=file, height = 1000, width = 1000, res = 120)
      else
        pdf(file=file)
      
      print(pie_plot(peak_dt()$tabledt))
      
      dev.off()
    }
  )
  
  output$peak_download_PD <- downloadHandler(
    
    filename = function(){
      paste("Peak Distribution Plot", Sys.time(), ".", input$peak_download_PDtype, sep = "")
    },
    
    content = function(file){
      if(input$peak_download_PDtype == "png")
        png(file=file, height = 1200, width = 2500, res = 250)
      else
        pdf(file=file, height = 7, width = 15)
      
      print(distribution_plot(peak_dt()$list))
      
      dev.off()
    }
  )

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
  
  #download
  
  output$reads_download_table <- downloadHandler(
    
    filename = function(){
      paste("Reads table", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file){
      
      write.csv(reads_dt()$read_stat, file, row.names = F)
      
    }
  )
  
  output$reads_download_bar <- downloadHandler(
    
    filename = function(){
      paste("Reads Distribution", Sys.time(), ".", input$reads_download_bartype, sep = "")
    },
    
    content = function(file){
      if(input$reads_download_bartype == "png")
        png(file=file, height = 1500, width = 3700, res = 310)
      else
        pdf(file=file, height = 6, width = 15)
      
      print(fill_bar_plot(reads_dt()$bardt))
      
      dev.off()
    }
  )
  
  #igv-------------------------------------------------
  
  igv_dt <- eventReactive(input$IGVbutton,{
    m6apeak <- m6a.peaks.table
    m6asite <- m6a.sites.table
    return(list(m6apeak = m6apeak, m6asite = m6asite))
  })
  
  output$m6apeak <- renderDataTable(
    igv_dt()$m6apeak
  )
  
  output$m6asite <- renderDataTable(
    igv_dt()$m6asite
  )
  
  #download
  output$igv_download_peaktable <- downloadHandler(
    
    filename = function(){
      paste("m6A peak table", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file){
      
      write.csv(igv_dt()$m6apeak, file, row.names = F)
      
    }
  )
  
  output$igv_download_sitetable <- downloadHandler(
    
    filename = function(){
      paste("m6A sites table", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file){
      
      write.csv(igv_dt()$m6asite, file, row.names = F)
      
    }
  )
  
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
  
  #download
  output$dm_download_matrix <- downloadHandler(
    
    filename = function(){
      paste("m6A Quantification Matrix", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file){
      
      write.csv(dm_dt()$mat, file, row.names = F)
      
    }
  )
  
  output$dm_download_difftable <- downloadHandler(
    
    filename = function(){
      paste("Diffm6A Peaks Table", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file){
      
      write.csv(dm_dt()$res, file, row.names = F)
      
    }
  )
  
  output$DM_download_HM <- downloadHandler(
    
    filename = function(){
      paste("Heatmap of Differential Methylation", Sys.time(), ".", input$DM_download_HMtype, sep = "")
    },
    
    content = function(file){
      if(input$DM_download_HMtype == "png")
        png(file=file, height = 1000, width = 1100, res = 160)
      else
        pdf(file=file)
      
      print(heatmap_dm(mat = dm_dt()$dm_mat, coldt = dm_dt()$coldata))
      
      dev.off()
    }
  )
  
  output$DM_download_vol <- downloadHandler(
    
    filename = function(){
      paste("Volcano plot of Differential Methylation", Sys.time(), ".", input$DM_download_voltype, sep = "")
    },
    
    content = function(file){
      if(input$DM_download_voltype == "png")
        png(file=file, height = 1000, width = 1000, res = 190)
      else
        pdf(file=file)
      
      print(volcano_plot_dm(res = dm_dt()$res, groupname = dm_dt()$name))
      
      dev.off()
    }
  )
  
  output$DM_download_PCA <- downloadHandler(
    
    filename = function(){
      paste("PCA plot of Differential Methylation", Sys.time(), ".", input$DM_download_PCAtype, sep = "")
    },
    
    content = function(file){
      if(input$DM_download_PCAtype == "png")
        png(file=file, height = 700, width = 900, res = 130)
      else
        pdf(file=file, width = 8)
      
      print(pca_plot(mat = dm_dt()$dm_mat, coldt = dm_dt()$coldata))
      
      dev.off()
    }
  )
  
#different expression--------------------------------------------
  
  de_dt <- eventReactive(input$debutton,{
    matrix = expression.matrix
    deres = diffexpression.list[[which(names(diffexpression.list)==input$degroup)]]
    group1 = strsplit(input$degroup, "_vs_")[[1]][1]
    group2 = strsplit(input$degroup, "_vs_")[[1]][2]
    coldata0 = as.data.frame(design.matrix)
    coldata = subset(coldata0, Type==group1|Type==group2)
    coldata$Type = as.factor(coldata$Type)
    deg = subset(deres, abs(log2FoldChange)> 0.58 & pvalue < 0.05)
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
  
  #download
  output$de_download_matrix <- downloadHandler(
    
    filename = function(){
      paste("Expression Matrix", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file){
      
      write.csv(de_dt()$mat, file, row.names = F)
      
    }
  )
  
  output$de_download_difftable <- downloadHandler(
    
    filename = function(){
      paste("Differential Expression Table", Sys.time(), ".csv", sep = "")
    },
    
    content = function(file){
      
      write.csv(de_dt()$res, file, row.names = F)
      
    }
  )
  
  output$DE_download_HM <- downloadHandler(
    
    filename = function(){
      paste("Heatmap of Differential Expression", Sys.time(), ".", input$DE_download_HMtype, sep = "")
    },
    
    content = function(file){
      if(input$DE_download_HMtype == "png")
        png(file=file, height = 1000, width = 1100, res = 160)
      else
        pdf(file=file)
      
      print(heatmap_de(mat = de_dt()$de_mat, coldt = de_dt()$coldata))
      
      dev.off()
    }
  )
  
  output$DE_download_vol <- downloadHandler(
    
    filename = function(){
      paste("Volcano plot of Differential Expression", Sys.time(), ".", input$DE_download_voltype, sep = "")
    },
    
    content = function(file){
      if(input$DE_download_voltype == "png")
        png(file=file, height = 1000, width = 1000, res = 190)
      else
        pdf(file=file)
      
      print(volcano_plot_de(res = de_dt()$res, groupname = de_dt()$name))
      
      dev.off()
    }
  )
  
  output$DE_download_PCA <- downloadHandler(
    
    filename = function(){
      paste("PCA plot of Differential Expression", Sys.time(), ".", input$DE_download_PCAtype, sep = "")
    },
    
    content = function(file){
      if(input$DE_download_PCAtype == "png")
        png(file=file, height = 700, width = 900, res = 130)
      else
        pdf(file=file, width = 8)
      
      print(pca_plot(mat = de_dt()$de_mat, coldt = de_dt()$coldata))
      
      dev.off()
    }
  )
  
  
  
  
  
  
}