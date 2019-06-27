library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
library(shinyBS)
library(shinycssloaders)


ui<-dashboardPage(skin = "black",
  dashboardHeader(title = "m6A Viewer"),
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
    useShinyjs(),
    tabItems(
      tabItem(
        tabName = "introduction",
        h1("Welcome to m6A Viewer"),
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
                                             withSpinner(dataTableOutput("readstable"),type = 4)))),
              tabPanel("Reads Distribution",
                       withSpinner(plotOutput("readsbarplot")),
                       conditionalPanel(
                         condition = "input.readsbutton >0",
                         awesomeRadio(
                           inputId = "reads_download_bartype",
                           label = "File Type", 
                           choices = c("png", "pdf"),
                           selected = "png",
                           inline = TRUE
                         ),
                         downloadBttn("reads_download_bar", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")
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
                                                   helpText("Data for genomic and position distribution")),
                                                 withSpinner(dataTableOutput("peaktable"),type = 4))),
                    tabPanel("Motif",
                             withSpinner(plotOutput("motif")),
                             conditionalPanel(
                               condition = "input.peakbutton > 0 ",
                               awesomeRadio(
                                 inputId = "peak_download_MOTIFtype",
                                 label = "File Type", 
                                 choices = c("png", "pdf"),
                                 selected = "png",
                                 inline = TRUE
                               ),
                               downloadBttn("peak_download_MOTIF", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")
                             )),
                    tabPanel("Genomic Distribution",
                             withSpinner(plotOutput("Pieplot")),
                             conditionalPanel(
                               condition = "input.peakbutton >0",
                               awesomeRadio(
                                 inputId = "peak_download_PIEtype",
                                 label = "File Type", 
                                 choices = c("png", "pdf"),
                                 selected = "png",
                                 inline = TRUE
                               ),
                               downloadBttn("peak_download_PIE", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")
                             )),
                    tabPanel("Position Distribution",
                             withSpinner(plotOutput("peakdistri")),
                             conditionalPanel(
                               condition = "input.peakbutton >0",
                               awesomeRadio(
                                 inputId = "peak_download_PDtype",
                                 label = "File Type", 
                                 choices = c("png", "pdf"),
                                 selected = "png",
                                 inline = TRUE
                               ),
                               downloadBttn("peak_download_PD", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")
                             ))
                    
             )
      )
    )
  ),
  
      #IGV------------------------------------------------------------------
  
      tabItem(
          tabName = "IGV", h1("IGV"),hr(),
          fluidRow(
            column(width = 4,
                   box(width = NULL, solidHeader = TRUE, status = "primary", title = "DATA",
                       fileInput("file4", "Choose CSV File",
                                 multiple = TRUE,
                                 accept = c("text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".csv",
                                            ".matrix")),
                 
                 
                       checkboxInput("header4", "Header", TRUE),
                 
                 
                       radioButtons("sep4", "Separator",
                                     choices = c(Comma = ",",
                                                Semicolon = ";",
                                                Tab = "\t"),
                                     selected = "\t"),
                 
                 
                       radioButtons("quote4", "Quote",
                                    choices = c(None = "",
                                                "Double Quote" = '"',
                                                "Single Quote" = "'"),
                                    selected = '"'),
                 
                 
                       radioButtons("disp4", "Display",
                                    choices = c(Head = "head",
                                                All = "all"),
                                    selected = "head"),
                       actionButton("dedmbutton", "GO PLOT!"),
                       uiOutput("selectfile4"))),
             
            column(width = 8,
                   tabBox(width = NULL,title = "TYPES", id = "tabset4", side = "right",
                          tabPanel("Venn Diagram",
                                   plotOutput("venn_dedm")),
                          tabPanel("Scatter Plot",
                                   plotOutput("Scatter_dedm")),      
                          tabPanel("ECDF Plot",
                                   plotOutput("ecdf_dedm")),
                          tabPanel("Box Plot",
                                   plotOutput("box_dedm")),
                          tabPanel("Enrichment Analysis",
                                   plotOutput("go_dedm"))
                   )
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
                                                            withSpinner(dataTableOutput("dmmatrix"),type = 4))),
                   tabPanel("Diffm6A Peaks Table",div(style = 'overflow-x: scroll',
                                                      withSpinner(dataTableOutput("dmres"),type = 4))),
                   tabPanel("Heat Map",
                            withSpinner(plotOutput("dmheatmap")),
                            conditionalPanel(
                              condition = "input.dmbutton >0",
                              awesomeRadio(
                                inputId = "DM_download_HMtype",
                                label = "File Type", 
                                choices = c("png", "pdf"),
                                selected = "png",
                                inline = TRUE
                              ),
                              downloadBttn("DM_download_HM", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")
                            )
                            ),
                   tabPanel("Volcano Plot",
                            withSpinner(plotOutput("dmvolcano")),
                            conditionalPanel(
                              condition = "input.dmbutton >0",
                              awesomeRadio(
                                inputId = "DM_download_voltype",
                                label = "File Type", 
                                choices = c("png", "pdf"),
                                selected = "png",
                                inline = TRUE
                              ),
                              downloadBttn("DM_download_vol", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")
                            )),
                        
                   tabPanel("PCA Plot",
                            withSpinner(plotOutput("dmPCA")),
                            conditionalPanel(
                              condition = "input.dmbutton >0",
                              awesomeRadio(
                                inputId = "DM_download_PCAtype",
                                label = "File Type", 
                                choices = c("png", "pdf"),
                                selected = "png",
                                inline = TRUE
                              ),
                              downloadBttn("DM_download_PCA", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")))
                                 
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
                                     withSpinner(dataTableOutput("dematrix"),type = 4))),
                        tabPanel("Differential Expression Table",
                                 div(style = 'overflow-x: scroll',
                                     withSpinner(dataTableOutput("deres"),type = 4))),
                        tabPanel("Heat Map",
                                 withSpinner(plotOutput("deheatmap")),
                                 conditionalPanel(
                                   condition = "input.debutton >0",
                                   awesomeRadio(
                                     inputId = "DE_download_HMtype",
                                     label = "File Type", 
                                     choices = c("png", "pdf"),
                                     selected = "png",
                                     inline = TRUE
                                   ),
                                   downloadBttn("DE_download_HM", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")
                                 )),
                        tabPanel("Volcano Plot",
                                 withSpinner(plotOutput("devolcano")),
                                 conditionalPanel(
                                   condition = "input.debutton >0",
                                   awesomeRadio(
                                     inputId = "DE_download_voltype",
                                     label = "File Type", 
                                     choices = c("png", "pdf"),
                                     selected = "png",
                                     inline = TRUE
                                   ),
                                   downloadBttn("DE_download_vol", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")
                                 )),
                        tabPanel("PCA Plot",
                                 withSpinner(plotOutput("dePCA")),
                                 conditionalPanel(
                                   condition = "input.debutton >0",
                                   awesomeRadio(
                                     inputId = "DE_download_PCAtype",
                                     label = "File Type", 
                                     choices = c("png", "pdf"),
                                     selected = "png",
                                     inline = TRUE
                                   ),
                                   downloadBttn("DE_download_PCA", style = "material-flat", size = "sm", label = "DOWNLOAD PLOT")))
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

