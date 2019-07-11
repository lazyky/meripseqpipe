library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
library(shinyBS)
library(shinycssloaders)
library(DT)

ui<-dashboardPage(skin = "black",
  dashboardHeader(title = "m6A Report"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName = "introduction",icon = icon("bullseye")),
      menuItem("Quality Control", tabName = "qc",icon = icon("th"),
               menuSubItem("Reads", tabName = "reads"),
               menuSubItem("m6A Peak", tabName = "peak")),
      menuItem("IGV", tabName = "IGV",icon = icon("connectdevelop")),
      menuItem("Different Methylation (DM)", tabName = "dm",icon = icon("bar-chart-o")),
      menuItem("Different Expression (DE)", tabName = "de",icon = icon("gg")),
      menuItem("Citation", tabName = "citation",icon = icon("bookmark")),
      menuItem("Contact", tabName = "contact us",icon = icon("podcast")),
      menuItem("Help", tabName = "help",icon = icon("question"))
    )
  ),
  #introduction------------------------------------------------------------
  
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "introduction",
        h1("Welcome to m6AReport"),
        fluidRow(
          box(width = 10, h2("Introduction"),hr(),
              p("m6A (N6-Methyladenosine) is a common methylated modification on RNA. It is mainly found in mRNA, and it can also be found in other type of RNA such as rRNA and lncRNA. 
                m6A is a significant content in epigenetic studies and plays an important role in basic life course. High-throughput sequencing, to be specific, MeRIP-seq, is the common method to study m6A. 
                At present, there are problems in the m6A analysis process, such as difficulty in quality control, difficulty in interpreting the analysis results, lack of visualization schemes, etc.. Moreover,  
                a relatively high level of programming skills is required when applying most of tools for analysing m6A data, but many biology researchers have little programming experience."),
              p("Here, we introduce m6A Analyzer, an interactive tool based on shiny-R for qualitative and quantitative analysis of m6A data, which has a graphical user interface (GUI), making data analysis easy and fast. m6A analyzer includes five analysis modules:", 
                strong("Quality Control, Differential Methylation Analysis, Differential Gene Expression Analysis, Correlation Analysis of m6A and Gene Expression"), "and",strong("Correlation Analysis of m6A and Translation.")))
        )),
  #quality control----------------------------------------------------------
      
      tabItem(
        tabName = "reads",h1("Quality Control of Your m6A Data: Reads"), hr(),
        fluidRow(
         column(width = 4,
                box(width = NULL,solidHeader = TRUE,title = "DATA",
                    p("Click the button to visualize the reads distribution and reads statistic of your data"),
                    actionButton("readsbutton", label = "GO PLOT!"))),
                  
         column(width = 8, 
             tabBox(width = NULL,title = "TYPES",id = "readsplot",side = "right",
              tabPanel("Reads Statistic", div(style = 'overflow-x: scroll', 
                                           conditionalPanel(
                                             condition = "input.readsbutton > 0",
                                             withSpinner(dataTableOutput("readstable"),type = 4),
                                             downloadButton("reads_download_table", label = "Download Table")))),
              tabPanel("Reads Distribution",
                       withSpinner(plotOutput("readsbarplot")),
                       conditionalPanel(
                         condition = "input.readsbutton >0",
                         radioButtons(
                           inputId = "reads_download_bartype",
                           label = "File Type", 
                           choices = c("png", "pdf"),
                           selected = "png",
                           inline = TRUE
                         ),
                         downloadButton("reads_download_bar", label = "DOWNLOAD PLOT")
                       ))
            )
          )
        )
      ),
  
  tabItem(
    tabName = "peak",h1("Quality Control of Your m6A Data: m6A peak"), hr(),
    fluidRow(
      column(width = 4,
             box(width = NULL,solidHeader = TRUE,title = "DATA",
                 selectInput("peak", label = "Group", 
                             choices = names(QC.peaks.list)),
                 actionButton("peakbutton", label = "GO PLOT!"))),
      
      column(width = 8, 
             tabBox(width = NULL,title = "TYPES",id = "peakplot",side = "right",
                    tabPanel("Data Display", div(style = 'overflow-x: scroll', 
                                                 conditionalPanel(
                                                   condition = "input.peakbutton > 0",
                                                   helpText("Data for genomic and position distribution"),
                                                   withSpinner(dataTableOutput("peaktable"),type = 4),
                                                  downloadButton("peak_download_table", label = "Download Table")))),
                    tabPanel("Motif",
                             withSpinner(plotOutput("motif", height = "800px")),
                             conditionalPanel(
                               condition = "input.peakbutton > 0 ",
                               radioButtons(
                                 inputId = "peak_download_MOTIFtype",
                                 label = "File Type", 
                                 choices = c("png", "pdf"),
                                 selected = "png",
                                 inline = TRUE
                               ),
                               downloadButton("peak_download_MOTIF", label = "DOWNLOAD PLOT")
                             )),
                    tabPanel("Genomic Distribution",
                             withSpinner(plotOutput("Pieplot")),
                             conditionalPanel(
                               condition = "input.peakbutton >0",
                               radioButtons(
                                 inputId = "peak_download_PIEtype",
                                 label = "File Type", 
                                 choices = c("png", "pdf"),
                                 selected = "png",
                                 inline = TRUE
                               ),
                               downloadButton("peak_download_PIE", label = "DOWNLOAD PLOT")
                             )),
                    tabPanel("Position Distribution",
                             withSpinner(plotOutput("peakdistri")),
                             conditionalPanel(
                               condition = "input.peakbutton >0",
                               radioButtons(
                                 inputId = "peak_download_PDtype",
                                 label = "File Type", 
                                 choices = c("png", "pdf"),
                                 selected = "png",
                                 inline = TRUE
                               ),
                               downloadButton("peak_download_PD", label = "DOWNLOAD PLOT")
                             ))
                    
             )
      )
    )
  ),
  
      #IGV------------------------------------------------------------------
  
  tabItem(
    tabName = "IGV", h1("IGV Browser"),hr(),
    fluidRow(
      tabBox(width = NULL,title = "TYPES", id = "tabset4", side = "right",
             tabPanel("Browser",
                      conditionalPanel(
                        condition="input.IGVbutton > 0",
                        tags$head(
                          includeScript("igv.min.js")
                        ),
                        div(id="igvDiv"),
                        tags$script(src="igv.js")
                      ),
                      p(strong("Click button to view IGV browser and m6A data tables")),
                      actionButton("IGVbutton", label = "GO!")),
             tabPanel("m6A peaks table",
                      div(style = 'overflow-x: scroll',
                          conditionalPanel(
                            condition="input.IGVbutton > 0",
                            withSpinner(dataTableOutput("m6apeak"),type = 4),
                            downloadButton("igv_download_peaktable", label = "Download Table")))),      
             tabPanel("m6A sites table",
                      div(style = 'overflow-x: scroll',
                          conditionalPanel(
                            condition="input.IGVbutton > 0",
                            withSpinner(dataTableOutput("m6asite"),type = 4),
                            downloadButton("igv_download_sitetable", label = "Download Table"))))
             
      )
      
    )),
      
      #different methylation------------------------------------------------------------------
      
      tabItem(
        tabName = "dm", 
        h1("Different Methylation"),hr(),
        fluidRow(
          column(width = 4,
                 box(width = NULL,solidHeader = TRUE,title = "DATA",
                     selectInput("dmgroup", label = "Contrast Group",
                                 choices = names(diffm6A.list)),
                     actionButton("dmbutton", "GO PLOT!"))),
          
          column(width = 8,
            tabBox(width = NULL,title = "TYPES", id = "dmplot", side = "right",
                   tabPanel("m6A Quantification Matrix",div(style = 'overflow-x: scroll',
                                                            conditionalPanel(
                                                            condition="input.dmbutton > 0",
                                                            withSpinner(dataTableOutput("dmmatrix"),type = 4),
                                                            downloadButton("dm_download_matrix", label = "Download Table")))),
                   tabPanel("Diffm6A Peaks Table",div(style = 'overflow-x: scroll',
                                                      conditionalPanel(
                                                        condition="input.dmbutton > 0",
                                                        withSpinner(dataTableOutput("dmres"),type = 4),
                                                        downloadButton("dm_download_difftable", label = "Download Table")))),
                   tabPanel("Heat Map",
                            withSpinner(plotOutput("dmheatmap", height = "650px")),
                            conditionalPanel(
                              condition = "input.dmbutton >0",
                              radioButtons(
                                inputId = "DM_download_HMtype",
                                label = "File Type", 
                                choices = c("png", "pdf"),
                                selected = "png",
                                inline = TRUE
                              ),
                              downloadButton("DM_download_HM", label = "DOWNLOAD PLOT")
                            )
                            ),
                   tabPanel("Volcano Plot",
                            withSpinner(plotOutput("dmvolcano", height = "650px")),
                            conditionalPanel(
                              condition = "input.dmbutton >0",
                              radioButtons(
                                inputId = "DM_download_voltype",
                                label = "File Type", 
                                choices = c("png", "pdf"),
                                selected = "png",
                                inline = TRUE
                              ),
                              downloadButton("DM_download_vol", label = "DOWNLOAD PLOT")
                            )),
                        
                   tabPanel("PCA Plot",
                            withSpinner(plotOutput("dmPCA", height = "650px")),
                            conditionalPanel(
                              condition = "input.dmbutton >0",
                              radioButtons(
                                inputId = "DM_download_PCAtype",
                                label = "File Type", 
                                choices = c("png", "pdf"),
                                selected = "png",
                                inline = TRUE
                              ),
                              downloadButton("DM_download_PCA", label = "DOWNLOAD PLOT")))
                                 
              )
            )
          )
        ),
      
      #different expression------------------------------------------------------------------
      
      tabItem(
        tabName = "de", h1("Different Expression"),hr(),
        fluidRow(
          column(width = 4,
                 box(width = NULL, solidHeader = TRUE, title = "DATA",
                         selectInput("degroup", label = "Contrast Group",
                                     choices = names(diffexpression.list)),
                         actionButton("debutton", "GO PLOT!")
                  )
                ),
          
          column(width = 8,
                 tabBox(width = NULL,title = "TYPES", id = "deplot", side = "right",
                        tabPanel("Expression Matrix",
                                 div(style = 'overflow-x: scroll',
                                     conditionalPanel(
                                     condition="input.debutton > 0",
                                     withSpinner(dataTableOutput("dematrix"),type = 4),
                                     downloadButton("de_download_matrix", label = "Download Table")))),
                        tabPanel("Differential Expression Table",
                                 div(style = 'overflow-x: scroll',
                                     conditionalPanel(
                                     condition="input.debutton > 0",
                                     withSpinner(dataTableOutput("deres"),type = 4),
                                     downloadButton("de_download_difftable", label = "Download Table")))),
                        tabPanel("Heat Map",
                                 withSpinner(plotOutput("deheatmap", height = "650px")),
                                 conditionalPanel(
                                   condition = "input.debutton >0",
                                   radioButtons(
                                     inputId = "DE_download_HMtype",
                                     label = "File Type", 
                                     choices = c("png", "pdf"),
                                     selected = "png",
                                     inline = TRUE
                                   ),
                                   downloadButton("DE_download_HM", label = "DOWNLOAD PLOT")
                                 )),
                        tabPanel("Volcano Plot",
                                 withSpinner(plotOutput("devolcano", height = "650px")),
                                 conditionalPanel(
                                   condition = "input.debutton >0",
                                   radioButtons(
                                     inputId = "DE_download_voltype",
                                     label = "File Type", 
                                     choices = c("png", "pdf"),
                                     selected = "png",
                                     inline = TRUE
                                   ),
                                   downloadButton("DE_download_vol", label = "DOWNLOAD PLOT")
                                 )),
                        tabPanel("PCA Plot",
                                 withSpinner(plotOutput("dePCA", height = "650px")),
                                 conditionalPanel(
                                   condition = "input.debutton >0",
                                   radioButtons(
                                     inputId = "DE_download_PCAtype",
                                     label = "File Type", 
                                     choices = c("png", "pdf"),
                                     selected = "png",
                                     inline = TRUE
                                   ),
                                   downloadButton("DE_download_PCA", label = "DOWNLOAD PLOT")))
          )
        )
      ) 
    ),
      
      
    #--------------------------------------------------------------------------
      tabItem(tabName = "citation",h1("CITATION"),
              fluidRow(
                column(width=12,
                       box(width = NULL),
                       box(width = NULL),
                       box(width = NULL))
              )),
      
      tabItem(tabName = "contact",h1("CONTACT US"),
              fluidRow(
                column(width=12,
                       box(width = 10),
                       box(width = 10))
              )),
      tabItem(tabName = "help")
    )
  )
)

