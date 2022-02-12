# ###########################################################
# #           This file contain UI widgets.                ##
# ###########################################################



####################
#pain
###################
pain_panel <- tabPanel(
  titlePanel(h4("BPS")),
  div(id="sidePanelGroupCmprsn_pain",sidebarPanel(
    width = 2,
    actionButton(inputId = "completeDATA_pain", label = "Complete Data",icon = icon("play-circle")),
    h6("OR", align = "center"),
    div(id = "geneNames_pain_form",selectizeInput(
      'geneNames_pain',
      label = 'Subset',
      choices = "",
      selected = NULL , multiple = TRUE,
      options = list(maxItems = 10,'plugins' = list('remove_button'),
                     'create' = TRUE,
                     'persist' = FALSE)
    )),
    hr(style = "border-color: grey"),
    fluidRow(column(12,numericInput("width_Pain", "Plot Width", min = 8, max = 80, value = 8))),
    fluidRow(column(12,numericInput("height_Pain", "Plot Height", min = 8, max = 80, value = 8))),
    )
    ),
  
  
  mainPanel(
    tabsetPanel(type = "tabs",
        ####################
        #Tables
        ###################
        tabPanel("Tables",
                 
            tabsetPanel(type = "tabs",
                        ######Raw Read counts######
                        tabPanel("Read Counts",
                                 
                                 fluidRow(
                                   #column(1,offset = 0),
                                   column(12,align="center",DT::dataTableOutput("painTable")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; padding: 10px;')
                                 )),
                        ######edgeR results######
                        tabPanel("edgeR result (all genes)",
                                 fluidRow(
                                   column(1,offset = 0),
                                   column(11,align="center",DT::dataTableOutput("painTable_edger")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; padding: 10px;')
                                 )),
                        ######edgeR results(DEGs genes onyl)######
                        tabPanel("edgeR result (DEGs only)", 
                                 fluidRow(style = "font-size: 80%",
                                   column(1,offset = 0),
                                   column(11,align="center",DT::dataTableOutput("painTable_edger_degs")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; ,padding: 10px;')
                                 ))
            )
        ),
        ####################
        #Plots
        ###################
        tabPanel("Plots",        
                 tabsetPanel(type = "tabs",
                 ######Volcano Plots######
                 tabPanel("Volcano Plot",
                          br(),
                          fluidRow( 
                            column(2,downloadButton(outputId = "down_volcanoPain", label = "    .pdf"))
                          ),
                          fluidRow(
                            column(1,offset = 0),
                            column(11,plotOutput("volcanoPain")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;padding: 10px;'),
                          ),       
            ),
                 

                 # ######Box Plots######
                 tabPanel("Box Plot",
                          br(),
                          fluidRow( 
                            column(2,downloadButton(outputId = "down_boxplot_pain", label = "    .pdf"))
                          ),
                          fluidRow(
                            column(1,offset = 0),
                            column(11,plotOutput("boxplot_pain")%>% withSpinner(color="#0dc5c1")),
                          ), 
                 ),
                #######Line Plots######
                tabPanel("Line Plot",
                         fluidRow(
                           column(1,offset = 0),
                           column(11,plotlyOutput("lineplot_pain")%>% withSpinner(color="#0dc5c1")),
                         )
                )



        )
    ),
    ####################
    #Enrichment results
    ###################
    tabPanel("GO Enrichment Plots",        
             tabsetPanel(type = "tabs",
                         
                         # ######GSEA######
                         tabPanel("GSEA Ridge Plot",
                                  br(),
                                  fluidRow( 
                                    column(2,downloadButton(outputId = "down_gsea_ridgePlot_pain", label = "  .pdf"))
                                  ),
                                  fluidRow(
                                    column(1,offset = 0),
                                    column(11,plotOutput("gsea_ridgePlot_pain")%>% withSpinner(color="#0dc5c1")),
                                  )

                                  ),
                         
                         ######GO ORA Treemap Plot######
                         tabPanel("ORA Treemap",
                                  br(),
                                  fluidRow( 
                                    column(2,downloadButton(outputId = "down_goORATreemap_pain", label = "  .pdf"))
                                  ),
                                  br(),
                                  fluidRow(
                                    column(1,offset = 0),
                                    column(11,plotOutput("goORATreemap_pain")%>% withSpinner(color="#0dc5c1")),
                                  )),
                         
                         ######GO ORA Interactive Plot######
                         tabPanel("ORA COW Plot",
                                  br(),
                                  fluidRow( 
                                    column(2,downloadButton(outputId = "down_goORA_Cowplot_pain", label = "    .pdf")),
                                    column(4,offset = 0),
                                    column(2,disabled(actionButton(inputId = "makeMeInteractivebtn_pain", label = "Interactive Plot", icon = icon("play-circle"))),),
                                    column(1,offset = 0),
                                    column(2,actionButton(inputId = "showSideBarPanel_pain", label = "Show Sidebar Panel", icon = icon("play-circle")),)
                                    
                                  ),
                                  
                                  fluidRow(
                                    column(1,offset = 0),
                                    column(11,plotOutput("goORA_Cowplot_pain")%>% withSpinner(color="#0dc5c1")),
                                          ),
                                    br(), br(),br(),br(),br(),
                                      br(), br(),br(),br(),br(),
                                      br(), br(),br(),br(),br(),
                                      fluidRow(
                                        plotOutput("interactiveGoPlot_pain")
                                      ),
                                    )

             )
    )
    )
  )
  
)

####################
#oab
###################
oab_panel <- tabPanel(
  titlePanel(h4("DO")),
  div(id="sidePanelGroupCmprsn_OAB",sidebarPanel(
    width = 2,
    actionButton(inputId = "completeDATA_OAB", label = "Complete Data",icon = icon("play-circle")),
    h6("OR", align = "center"),
    div(id = "geneNames_OAB_form",selectizeInput(
      'geneNames_OAB',
      label = 'Subset',
      choices = "",
      selected = NULL , multiple = TRUE,
      options = list(maxItems = 10,'plugins' = list('remove_button'),
                     'create' = TRUE,
                     'persist' = FALSE)
    )),
    hr(style = "border-color: grey"),
    fluidRow(column(12,numericInput("width_OAB", "Plot Width", min = 8, max = 80, value = 8))),
    fluidRow(column(12,numericInput("height_OAB", "Plot Height", min = 8, max = 80, value = 8))),
  )
  ),
  
  
  mainPanel(
    tabsetPanel(type = "tabs",
                ####################
                #Tables
                ###################
                tabPanel("Tables",
                         
                         tabsetPanel(type = "tabs",
                                     ######Raw Read counts######
                                     tabPanel("Read Counts",
                                              
                                              fluidRow(style = "font-size: 85%",
                                                column(12,align="center",DT::dataTableOutput("oabTable")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; padding: 10px;')
                                              )),
                                     ######edgeR results######
                                     tabPanel("edgeR result (all genes)",
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,align="center",DT::dataTableOutput("oabTable_edger")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; padding: 10px;')
                                              )),
                                     ######edgeR results(DEGs genes onyl)######
                                     tabPanel("edgeR result (DEGs only)", 
                                              fluidRow(style = "font-size: 80%",
                                                       column(1,offset = 0),
                                                       column(11,align="center",DT::dataTableOutput("oabTable_edger_degs")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; ,padding: 10px;')
                                              ))
                         )
                ),
                ####################
                #Plots
                ###################
                tabPanel("Plots",        
                         tabsetPanel(type = "tabs",
                                     ######Volcano Plots######
                                     tabPanel("Volcano Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_volcano_oab", label = "    .pdf"))
                                              ),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("volcano_oab")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;padding: 10px;'),
                                              ),       
                                              
                                     ),
                                     
                                     
                                     # ######Box Plots######
                                     tabPanel("Box Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_boxplot_OAB", label = "    .pdf"))
                                              ),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("boxplot_OAB")%>% withSpinner(color="#0dc5c1")),
                                              ), 
                                     ),
                                     #######Line Plots######
                                     tabPanel("Line Plot",
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotlyOutput("lineplot_OAB")%>% withSpinner(color="#0dc5c1")),
                                              )
                                     )
                                     
                                     
                                     
                         )
                ),
                ####################
                #Enrichment results
                ###################
                tabPanel("GO Enrichment Plots",        
                         tabsetPanel(type = "tabs",
                                     
                                     # ######GSEA######
                                     tabPanel("GSEA Ridge Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_gsea_ridgePlot_OAB", label = "  .pdf"))
                                              ),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("gsea_ridgePlot_OAB")%>% withSpinner(color="#0dc5c1")),
                                              )
                                     ),
                                     
                                     ######GO ORA Treemap Plot######
                                     tabPanel("ORA Treemap",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_goORATreemap_OAB", label = "  .pdf"))
                                              ),
                                              br(),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("goORATreemap_OAB")%>% withSpinner(color="#0dc5c1")),
                                              )                                  

                                     )
                                     ,
                                     ######GO ORA Interactive Plot######
                                     tabPanel("ORA COW Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_goORA_Cowplot_OAB", label = "    .pdf")),
                                                column(4,offset = 0),
                                                column(2,disabled(actionButton(inputId = "makeMeInteractivebtn_OAB", label = "Interactive Plot", icon = icon("play-circle"))),),
                                                column(1,offset = 0),
                                                column(2,actionButton(inputId = "showSideBarPanel_OAB", label = "Show Sidebar Panel", icon = icon("play-circle")),)
                                                
                                              ),
                                              
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("goORA_Cowplot_OAB")%>% withSpinner(color="#0dc5c1")),
                                              ),
                                              br(), br(),br(),br(),br(),
                                              br(), br(),br(),br(),br(),
                                              br(), br(),br(),br(),br(),
                                              fluidRow(
                                                plotOutput("interactiveGoPlot_OAB")
                                              ),
                                     )
                                     
                         )
                )
    )
  )
  
)

####################
#AC
###################
AC_panel <- tabPanel(
  titlePanel(h4("UA")),
  div(id="sidePanelGroupCmprsn_AC",sidebarPanel(
    width = 2,
    actionButton(inputId = "completeDATA_AC", label = "Complete Data",icon = icon("play-circle")),
    h6("OR", align = "center"),
    div(id = "geneNames_AC_form",selectizeInput(
      'geneNames_AC',
      label = 'Subset',
      choices = "",
      selected = NULL , multiple = TRUE,
      options = list(maxItems = 10,'plugins' = list('remove_button'),
                     'create' = TRUE,
                     'persist' = FALSE)
    )),
    hr(style = "border-color: grey"),
    fluidRow(column(12,numericInput("width_AC", "Plot Width", min = 8, max = 80, value = 8))),
    fluidRow(column(12,numericInput("height_AC", "Plot Height", min = 8, max = 80, value = 8))),
  )
  ),
  
  
  mainPanel(
    tabsetPanel(type = "tabs",
                ####################
                #Tables
                ###################
                tabPanel("Tables",
                         
                         tabsetPanel(type = "tabs",
                                     ######Raw Read counts######
                                     tabPanel("Read Counts",
                                              
                                              fluidRow(style = "font-size: 85%",
                                                       #column(1,offset = 0),
                                                       column(12,align="center",DT::dataTableOutput("ACTable")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; padding: 10px;')
                                              )),
                                     ######edgeR results######
                                     tabPanel("edgeR result (all genes)",
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,align="center",DT::dataTableOutput("ACTable_edger")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; padding: 10px;')
                                              )),
                                     ######edgeR results(DEGs genes onyl)######
                                     tabPanel("edgeR result (DEGs only)", 
                                              fluidRow(style = "font-size: 80%",
                                                       column(1,offset = 0),
                                                       column(11,align="center",DT::dataTableOutput("ACTable_edger_degs")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; ,padding: 10px;')
                                              ))
                         )
                ),
                ####################
                #Plots
                ###################
                tabPanel("Plots",        
                         tabsetPanel(type = "tabs",
                                     ######Volcano Plots######
                                     tabPanel("Volcano Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_volcano_AC", label = "    .pdf"))
                                              ),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("volcano_AC")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;padding: 10px;'),
                                              ),       
                                     ),
                                     
                                     
                                     # ######Box Plots######
                                     tabPanel("Box Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_boxplot_AC", label = "    .pdf"))
                                              ),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("boxplot_AC")%>% withSpinner(color="#0dc5c1")),
                                              ), 
                                     ),
                                     #######Line Plots######
                                     tabPanel("Line Plot",
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotlyOutput("lineplot_AC")%>% withSpinner(color="#0dc5c1")),
                                              )
                                     )
                                     
                                     
                                     
                         )
                ),
                ####################
                #Enrichment results
                ###################
                tabPanel("GO Enrichment Plots",        
                         tabsetPanel(type = "tabs",
                                     
                                     # ######GSEA######
                                     tabPanel("GSEA Ridge Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_gsea_ridgePlot_AC", label = "  .pdf"))
                                              ),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("gsea_ridgePlot_AC")%>% withSpinner(color="#0dc5c1")),
                                              )

                                     ),
                                     
                                     ######GO ORA Treemap Plot######
                                     tabPanel("ORA Treemap",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_goORATreemap_AC", label = "  .pdf"))
                                              ),
                                              br(),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("goORATreemap_AC")%>% withSpinner(color="#0dc5c1")),
                                              )                                  
                                         ),
                                     ######GO ORA Interactive Plot######
                                     tabPanel("ORA COW Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_goORA_Cowplot_AC", label = "    .pdf")),
                                                column(4,offset = 0),
                                                column(2,disabled(actionButton(inputId = "makeMeInteractivebtn_AC", label = "Interactive Plot", icon = icon("play-circle"))),),
                                                column(1,offset = 0),
                                                column(2,actionButton(inputId = "showSideBarPanel_AC", label = "Show Sidebar Panel", icon = icon("play-circle")),)
                                                
                                              ),
                                              
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("goORA_Cowplot_AC")%>% withSpinner(color="#0dc5c1")),
                                              ),
                                              br(), br(),br(),br(),br(),
                                              br(), br(),br(),br(),br(),
                                              br(), br(),br(),br(),br(),
                                              fluidRow(
                                                plotOutput("interactiveGoPlot_AC")
                                              ),
                                     )
                                     
                         )
                )
    )
  )
  
)

####################
#BOO
###################
BOO_panel <- tabPanel(
  titlePanel(h4("BOO")),
  div(id="sidePanelGroupCmprsn_BOO",sidebarPanel(
    width = 2,
    actionButton(inputId = "completeDATA_BOO", label = "Complete Data",icon = icon("play-circle")),
    h6("OR", align = "center"),
    div(id = "geneNames_BOO_form",selectizeInput(
      'geneNames_BOO',
      label = 'Subset',
      choices = "",
      selected = NULL , multiple = TRUE,
      options = list(maxItems = 10,'plugins' = list('remove_button'),
                     'create' = TRUE,
                     'persist' = FALSE)
    )),
    hr(style = "border-color: grey"),
    fluidRow(column(12,numericInput("width_BOO", "Plot Width", min = 8, max = 80, value = 8))),
    fluidRow(column(12,numericInput("height_BOO", "Plot Height", min = 8, max = 80, value = 8))),
  )
  ),
  
  
  mainPanel(
    tabsetPanel(type = "tabs",
                ####################
                #Tables
                ###################
                tabPanel("Tables",
                         
                         tabsetPanel(type = "tabs",
                                     ######Raw Read counts######
                                     tabPanel("Read Counts",
                                              
                                              fluidRow(style = "font-size: 85%",
                                                       #column(1,offset = 0),
                                                       column(12,align="center",DT::dataTableOutput("BOOTable")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; padding: 10px;')
                                              )),
                                     ######edgeR results######
                                     tabPanel("edgeR result (all genes)",
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,align="center",DT::dataTableOutput("BOOTable_edger")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; padding: 10px;')
                                              )),
                                     ######edgeR results(DEGs genes onyl)######
                                     tabPanel("edgeR result (DEGs only)", 
                                              fluidRow(style = "font-size: 80%",
                                                       column(1,offset = 0),
                                                       column(11,align="center",DT::dataTableOutput("BOOTable_edger_degs")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;border:1px solid; ,padding: 10px;')
                                              ))
                         )
                ),
                ####################
                #Plots
                ###################
                tabPanel("Plots",        
                         tabsetPanel(type = "tabs",
                                     ######Volcano Plots######
                                     tabPanel("Volcano Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_volcano_BOO", label = "    .pdf"))
                                              ),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("volcano_BOO")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px;padding: 10px;'),
                                              ),       

                                     ),
                                     
                                     
                                     # ######Box Plots######
                                     tabPanel("Box Plot",
   
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_boxplot_BOO", label = "    .pdf"))
                                              ),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("boxplot_BOO")%>% withSpinner(color="#0dc5c1")),
                                              ), 
                                     ),
                                     #######Line Plots######
                                     tabPanel("Line Plot",
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotlyOutput("lineplot_BOO")%>% withSpinner(color="#0dc5c1")),
                                              )
                                     )
                                     
                                     
                                     
                         )
                ),
                ####################
                #Enrichment results
                ###################
                tabPanel("GO Enrichment Plots",        
                         tabsetPanel(type = "tabs",
                                     
                                     # ######GSEA######
                                     tabPanel("GSEA Ridge Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_gsea_ridgePlot_BOO", label = "  .pdf"))
                                              ),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("gsea_ridgePlot_BOO")%>% withSpinner(color="#0dc5c1")),
                                              )

                                              
                                     ),
                                     
                                     ######GO ORA Treemap Plot######
                                     tabPanel("ORA Treemap",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_goORATreemap_BOO", label = "  .pdf"))
                                              ),
                                              br(),
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("goORATreemap_BOO")%>% withSpinner(color="#0dc5c1")),
                                              )                                  
                                              
                                              
                                     )
                                     ,
                                     ######GO ORA Interactive Plot######
                                     tabPanel("ORA COW Plot",
                                              br(),
                                              fluidRow( 
                                                column(2,downloadButton(outputId = "down_goORA_Cowplot_BOO", label = "    .pdf")),
                                                column(4,offset = 0),
                                                column(2,disabled(actionButton(inputId = "makeMeInteractivebtn_BOO", label = "Interactive Plot", icon = icon("play-circle"))),),
                                                column(1,offset = 0),
                                                column(2,actionButton(inputId = "showSideBarPanel_BOO", label = "Show Sidebar Panel", icon = icon("play-circle")),)
                                                
                                              ),
                                              
                                              fluidRow(
                                                column(1,offset = 0),
                                                column(11,plotOutput("goORA_Cowplot_BOO")%>% withSpinner(color="#0dc5c1")),
                                              ),
                                              br(), br(),br(),br(),br(),
                                              br(), br(),br(),br(),br(),
                                              br(), br(),br(),br(),br(),
                                              fluidRow(
                                                plotOutput("interactiveGoPlot_BOO")
                                              ),
                                     )
                                     
                         )
                )
    )
  )
  
)

####################
#goEnrichment_panel
###################
goEnrichment_panel <- tabPanel(
  titlePanel(h4("Functional Enrichment Exploration")),
  div(id="sidePanelGO_enrich",sidebarPanel(
    width = 2,
    
    selectizeInput(
      'groupName',
      label = 'Select Group',
      choices = c("BPS vs Control","DO vs Control","UA vs Control","BOO vs Control"),
      selected = "BPS vs Control"
    ),
    selectizeInput(
      'enrichmentType',
      label = 'Enrichment Type',
      choices = c("ORA","GSEA"),
      selected = "ORA"
    ),
    selectizeInput(
      'goTermName',
      label = 'Select GO Term',
      choices = "",
      selected = NULL
    ),
    
    actionButton(inputId = "plotVolcanoFor_goTerm", label = "Plot",icon = icon("play-circle")),
    
    hr(style = "border-color: grey"),
    fluidRow(column(12,numericInput("width_goTerm", "Plot Width", min = 8, max = 80, value = 8))),
    fluidRow(column(12,numericInput("height_goTerm", "Plot Height", min = 8, max = 80, value = 8))),
    downloadButton(outputId = "down_volcanoGoTerm", label = " Download")
  )
  ),

  mainPanel(fluidRow(
    column(1,offset = 0),
    column(11,plotOutput("volcano_plot_ForGOTerm")%>% withSpinner(color="#0dc5c1")),
  ),
  br(),br(),br(),
  fluidRow(
    column(12,align="center",DT::dataTableOutput("table_ForGOTerm")%>% withSpinner(color="#0dc5c1"), style='margin-bottom:30px; padding: 10px;')
  )
  )  
)


##########################################################
#                        about                          ##
###########################################################

about <- tabPanel(
  
  titlePanel(h4("About")),

    # General info.
    tabPanel(
      "Overview",
      tags$h2("About LUTD-Base"),
      tags$p(HTML("LUTD-Base is an easy, interactive, and user-friendly shiny 
                  based web application that allows the other researcher to 
                  deeply explore the bulk RNA-Seq data, differential expression,
                  functional enrichment analysis results, etc. from multiple 
                  LUTD related studies carried out by our group.")),
      
      tags$p(HTML("LUTD-Base web application is developed at
                  <a href=\"http://www.urofun.ch/\" target=\"_blank\"> Functional Urology</a>
                  lab at the 
                  <a href=\"https://www.unibe.ch/index_eng.html/\" target=\"_blank\"> University of Bern</a>.
      Source code is available at <a href=\"https://github.com/FunctionalUrology/LUTD-Base/\" target=\"_blank\"> https://github.com/FunctionalUrology/LUTD-Base</a>.")),
      
      tags$h4("Application"),
      tags$p(HTML("LUTD-Base contains 5 main tabs, one for each LUTD subtype 
                  and a Functional Enrichment Exploration tab. Each tab consists
                  of 3 different panels named Tables, Plots, GO Enrichment plots. 
                  Users can explore the raw read counts and edgeR results using 
                  the table panel. It allows users to explore the read counts 
                  and differential expression analysis results for the 
                  subset of genes using the plot panel. GO Enrichment plots 
                  provide static as well as interactive plots to explore 
                  functional enrichment analysis results. All the plots and 
                  tables are available with a download option.")),
      
      tags$h4("Contact"),
      tags$p(HTML("Bug reports and feature requests can be communicated via:")),
      tags$ul(
        tags$li(HTML("Github: <a href=\"https://github.com/FunctionalUrology/LUTD-Base/issues\" target=\"_blank\">https://github.com/FunctionalUrology/LUTD-Base/issues</a>")),
        tags$li(HTML("Email: <a href=\"mailto:akshay.akshay@dbmr.unibe.ch\" target=\"_blank\">akshay.akshay@dbmr.unibe.ch</a> <br>
                     &emsp; &emsp; &ensp; 
                     <a href=\"mailto:Ali.HashemiGheinani@childrens.harvard.edu\" target=\"_blank\">Ali.HashemiGheinani@childrens.harvard.edu</a>"))
        ),
      
      tags$h4("Citation"),
      tags$p(HTML("If LUTD-Base helps you in any way, please cite the following articles:")),
      tags$ul(
        tags$li(HTML("Gheinani, Ali Hashemi et al. “Characterization of 
                      miRNA-regulated networks, hubs of signaling, and 
                      biomarkers in obstruction-induced bladder dysfunction.” 
                      JCI insight vol. 2,2 e89560. 26 Jan. 2017, doi:
                     <a href=\"https://doi.org/10.1172/jci.insight.89560\" target=\"_blank\">10.1172/jci.insight.89560</a>")),
        
        tags$li(HTML("Gheinani, Ali Hashemi et al. “Integrated mRNA-miRNA
                     transcriptome analysis of bladder biopsies from patients 
                     with bladder pain syndrome identifies signaling alterations
                     contributing to the disease pathogenesis.” 
                     BMC urology vol. 21,1 172. 7 Dec. 2021, doi:
                     <a href=\"https://doi.org/10.1186/s12894-021-00934-0\" target=\"_blank\">10.1186/s12894-021-00934-0</a>"))
      )
      
      ),
  
  
    
  
)


#combine everything
ui <- fluidPage(
  theme = shinytheme("journal"),
  navbarPage(
    pain_panel,
    oab_panel,
    AC_panel,
    BOO_panel,
    goEnrichment_panel,
    about,
    useShinyjs(),
   title= tags$h2(p(em("LUTD-Base")),style="color: #4d1a00;"),
   # Create Right Side Logo/Image with Link       
   tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><a href=\"https://www.unibe.ch/index_ger.html\" target=\"_blank\" ><img src=\"https://www.unibe.ch/media/logo_unibern@2x.png\" alt=\"alt\" style=\"float:right;width:142px;height:108.86 px;padding-top:2px;padding-bottom:2px;\"> </a>`</div>');
    console.log(header)")
   )
))

