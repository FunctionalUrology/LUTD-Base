
# ########################################################################
#                             load data file                             #
# ########################################################################
load(file="data1.RData")

interactive_heatmap_Input=readRDS(file = paste(getwd(),"output/interactive_heatmap_Input.rds",sep="/"))
ht=interactive_heatmap_Input$ht
ht_OAB=interactive_heatmap_Input$ht_OAB
ht_AC=interactive_heatmap_Input$ht_AC
ht_BOO=interactive_heatmap_Input$ht_BOO

# ########################################################################
#                               Define functions                         #
# ########################################################################

##############
#Gens Info
##############
getVolcano_plots_selectedGenes<-function(data,select)
{
  if(length(select)==0)
  {
    return("error")
  }
  options(ggrepel.max.overlaps = Inf)
  temp=data[select,]
  p=EnhancedVolcano(temp,
                    lab=rownames(temp),
                    x = 'logFC',
                    y = 'FDR',
                    legendLabels = c("NS", expression(Log[2]~FC),
                                     "FDR",
                                     expression(FDR ~and~ log[2]~FC)),
                    selectLab =select,
                    subtitle = NULL,
                    pCutoff = 0.1,
                    FCcutoff = 0.5,
                    axisLabSize = 12,
                    legendLabSize = 7,
                    pointSize = 2,
                    labSize = 4,
                    xlab = bquote(~Log[2]~ 'fold change(±0.5)'),
                    colAlpha = 1,
                    legendIconSize = 5.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.5,
                    ylab = bquote(~-Log[10] ~ italic(FDR(0.1))))
  return(p)
}



##########
#boxplot
##########
getBoxplot<- function(genes,readCountMatrix,groupName)
{
  if(length(genes)==0)
  {
    return("error")
  }
  
  plots=list()
  for(gene in genes)
  {
    # create a dataset
    data <- data.frame(
      group=c( rep("Control",6), rep(groupName,6))
    )
    data$readCount=as.numeric(as.vector(readCountMatrix[gene,]))
    
    plots[[gene]]=ggplot( data,aes(x=group, y=readCount,fill=group)) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      theme_ipsum(base_family = 'Helvetica') +
      theme(
        legend.position="none",
        plot.title = element_text(size=11,hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
      ) +
      ggtitle(gene) +
      labs(y="read count",x="")
    
  }
  return(plots)
  
}



##########
#lineplot
##########

getLineplot<- function(genes,readCountMatrix)
{
  if(length(genes)==0)
  {
    return("error")
  }
  
  readCount=c()
  for(gene in genes)
  {
    readCount=c(readCount,as.numeric(as.vector(readCountMatrix[gene,])))
  }
  
  data <- data.frame(matrix(nrow=length(genes)*12,ncol=3))
  colnames(data)=c("sample","readCount","gene")
  
  data$sample=rep(colnames(readCountMatrix) , length(genes))
  data$gene=rep(genes , each=12)
  data$readCount=readCount
  fig <- plot_ly(data, x = ~sample, y = ~readCount, color = ~gene)
  fig <- fig %>% add_lines()
  return(fig)
  
}

##########
#get Volcano plots for selected go Term
##########
getVolcanoFor_goTerm<- function(groupNAME,enrichmentType,term)
{
  
  if (groupNAME=="BPS vs Control")
  {
    edgerR_result_table=edgeR_res_all_pain
    
    #type of enrichment
    if (enrichmentType=="ORA")
    {goORA_table=as.data.frame(enrichmentResults_pain[["go_over"]])}
    else{
      goORA_table=as.data.frame(gsea_pain_symbol)
      colnames(goORA_table)[11] <- "geneID"
    }}
  
  else if (groupNAME=="DO vs Control"){
    edgerR_result_table=edgeR_res_all_OAB
    
    #type of enrichment
    if (enrichmentType=="ORA")
    {goORA_table=as.data.frame(enrichmentResults_OAB[["go_over"]])}
    else{
      goORA_table=as.data.frame(gsea_OAB_symbol)
      colnames(goORA_table)[11] <- "geneID"
    }}
  
  else if (groupNAME=="UA vs Control"){
    edgerR_result_table=edgeR_res_all_AC
    
    #type of enrichment
    if (enrichmentType=="ORA")
    {goORA_table=as.data.frame(enrichmentResults_AC[["go_over"]])}
    else{
      goORA_table=as.data.frame(gsea_AC_symbol)
      colnames(goORA_table)[11] <- "geneID"
    }}
  
  else {
    edgerR_result_table=edgeR_res_all_BOO
    
    #type of enrichment
    if (enrichmentType=="ORA")
    {goORA_table=as.data.frame(enrichmentResults_BOO[["go_over"]])}
    else{
      goORA_table=as.data.frame(gsea_BOO_symbol)
      colnames(goORA_table)[11] <- "geneID"
    }}
  
  geneInvolved=goORA_table[goORA_table['Description'] ==term ,"geneID"]
  geneInvolved=unlist(strsplit(geneInvolved, "/"),use.names=FALSE)
  
  return(getVolcano_plots_selectedGenes(edgerR_result_table,geneInvolved))
}

##########
#get table for selected go Term
##########

getTable_goTerm<- function(groupNAME,enrichmentType,term)
{
  
  if (groupNAME=="BPS vs Control")
  {
    edgerR_result_table=edgeR_res
    rawCount_table=countdata
    
    #type of enrichment
    if (enrichmentType=="ORA")
    {goORA_table=as.data.frame(enrichmentResults_pain[["go_over"]])}
    else{
      goORA_table=as.data.frame(gsea_pain_symbol)
      colnames(goORA_table)[11] <- "geneID"
    }}
  
  else if (groupNAME=="DO vs Control"){
    edgerR_result_table=edgeR_res_OAB
    rawCount_table=countdata.OAB
    
    #type of enrichment
    if (enrichmentType=="ORA")
    {goORA_table=as.data.frame(enrichmentResults_OAB[["go_over"]])}
    else{
      goORA_table=as.data.frame(gsea_OAB_symbol)
      colnames(goORA_table)[11] <- "geneID"
    }}
  
  else if (groupNAME=="UA vs Control"){
    edgerR_result_table=edgeR_res_AC
    rawCount_table=countdata.AC
    
    #type of enrichment
    if (enrichmentType=="ORA")
    {goORA_table=as.data.frame(enrichmentResults_AC[["go_over"]])}
    else{
      goORA_table=as.data.frame(gsea_AC_symbol)
      colnames(goORA_table)[11] <- "geneID"
    }}
  
  else {
    edgerR_result_table=edgeR_res_BOO
    rawCount_table=countdata.BOO
    
    #type of enrichment
    if (enrichmentType=="ORA")
    {goORA_table=as.data.frame(enrichmentResults_BOO[["go_over"]])}
    else{
      goORA_table=as.data.frame(gsea_BOO_symbol)
      colnames(goORA_table)[11] <- "geneID"
    }}
  
  geneInvolved=goORA_table[goORA_table['Description'] ==term ,"geneID"]
  geneInvolved=unlist(strsplit(geneInvolved, "/"),use.names=FALSE)
  
  table_ForGOTerm=merge(rawCount_table[geneInvolved,], edgerR_result_table[geneInvolved,c("log2FoldChange","FDR","description")],by=0, all=TRUE)
  rownames(table_ForGOTerm)=table_ForGOTerm$Row.names
  table_ForGOTerm=table_ForGOTerm[2:length(table_ForGOTerm)]
  return(table_ForGOTerm)
}





###########################################################
##This part contain server side functions for computation.##
###########################################################

 
  server <- function(input,output, session){
    
    #Provide gene names for drop down in pain
    observe({
      updateSelectizeInput(session, 'geneNames_pain', choices = sort(rownames(edgeR_res_all_pain))[27:length(rownames(edgeR_res_all_pain))], server = TRUE)
    })
    
    #Provide gene names for drop down in oab
    observe({
      updateSelectizeInput(session, 'geneNames_OAB', choices = sort(rownames(edgeR_res_all_OAB)), server = TRUE)
    })
    
    #Provide gene names for drop down inr AC
    observe({
      updateSelectizeInput(session, 'geneNames_AC', choices = sort(rownames(edgeR_res_all_AC)), server = TRUE)
    })
    
    #Provide gene names for drop down in BOO
    observe({
      updateSelectizeInput(session, 'geneNames_BOO', choices = sort(rownames(edgeR_res_all_BOO)), server = TRUE)
    })
    

    ####################
    #Pain data panel
    ###################
    #common error message for empty gene list
    errorPlot<<-ggplot() +xlim(c(-5,5)) + ylim(c(-5,5)) +
      geom_text(size=8, aes(x = 0, y = 0, label = "Please select atleast one gene.")) + theme_void()
    
    
    #1. Tables
    output$painTable <- DT::renderDT(NULL)
    output$painTable_edger <- DT::renderDT(NULL)
    output$painTable_edger_degs <- DT::renderDT(NULL)
    
    colnameForEdgerRes=c("log2FoldChange","logCPM","PValue","FDR","description")
    options = list( 
      dom = "Blfrtip"
      , buttons = 
        list(list(
          extend = "csv",
          filename= 'dataTable',
          text='.CSV  <span class="glyphicon glyphicon-download-alt"></span>'
        )))
    
    output$painTable = DT::renderDataTable({datatable( data = countdata,
                                                       extensions = 'Buttons',
                                                       options = options)})
    
    output$painTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_pain),
                                                             extensions = 'Buttons',
                                                             options = options)})
    
    output$painTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res[colnameForEdgerRes],
                                                                  extensions = 'Buttons',
                                                                  options = options)},width=1000)
    
    observeEvent(input$geneNames_pain, {
        output$painTable = DT::renderDataTable({datatable( data = countdata[input$geneNames_pain,],
                                                           extensions = 'Buttons',
                                                           options = options)})
        
        output$painTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_pain)[input$geneNames_pain,],
                                                                 extensions = 'Buttons',
                                                                 options = options)
          })
        output$painTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res[input$geneNames_pain,colnameForEdgerRes],
                                                                      extensions = 'Buttons',
                                                                      options = options)
          })

    })
    
    restartFN <- observeEvent(input$completeDATA_pain, {
      output$painTable = DT::renderDataTable({datatable( data = countdata,
                                                         extensions = 'Buttons',
                                                         options = options)
        })
      output$painTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_pain),
                                                               extensions = 'Buttons',
                                                               options = options)
        })
      output$painTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res[colnameForEdgerRes],
                                                                    extensions = 'Buttons',
                                                                    options = options)
        })
      shinyjs::reset("geneNames_pain_form")
    })
    
    #2. Plots
    output$volcanoPain <- DT::renderDT(NULL)
    output$boxplot_pain <- DT::renderDT(NULL)
    output$lineplot_pain <- DT::renderDT(NULL)
    
    #volcano
          # Download
    output$down_volcanoPain<-downloadHandler(
      filename <- function(){paste0("volcano-plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_Pain, height = input$height_Pain)
        p=getVolcano_plots_selectedGenes(edgeR_res_all_pain,input$geneNames_pain)
        if(p=="error"){grid.arrange(volcanoPain)}
        else{grid.arrange(p)}
        dev.off()
        }
    )
    
    output$volcanoPain<-renderPlot({volcanoPain})
    observeEvent(input$geneNames_pain, {
        output$volcanoPain<-renderPlot({
          p=getVolcano_plots_selectedGenes(edgeR_res_all_pain,input$geneNames_pain)
          if(p=="error")
          {
            output$volcanoPain<-renderPlot({volcanoPain})
            
          }
          else{p}
          })
      
    })
    restartFN <- observeEvent(input$completeDATA_pain, {
      output$volcanoPain<-renderPlot({volcanoPain})
    })
    
    #random genes for initial boxplot and line plot
    randomGenes<<-rownames(edgeR_res_all_pain)[30:31]
    
    #boxplot
        # Download
        output$down_boxplot_pain<-downloadHandler(
          filename <- function(){paste0("boxplot.pdf")},
          content = function(file){
            pdf(file = file,onefile=TRUE, width = input$width_Pain, height = input$height_Pain)
            plots=getBoxplot(input$geneNames_pain,countdata,"BPS")
            if(plots=="error"){
              p=getBoxplot(randomGenes,countdata,"BPS")
              do.call("grid.arrange", c(p, ncol=2))
              }
            else{do.call("grid.arrange", c(plots, ncol=2))}
            dev.off()
          }
        )
        
    output$boxplot_pain<-renderPlot({
      p=getBoxplot(randomGenes,countdata,"BPS")
      do.call("grid.arrange", c(p, ncol=2))
    },height = 300)
    
    observeEvent(input$geneNames_pain, {
        output$boxplot_pain<-renderPlot({
          plots=getBoxplot(input$geneNames_pain,countdata,"BPS")
          
          if(plots=="error") #raise error if selectinput is empty
          {
            p=getBoxplot(randomGenes,countdata,"BPS")
            do.call("grid.arrange", c(p, ncol=2))
          }
          else{do.call("grid.arrange", c(plots, ncol=2))}
          },height =round(length(input$geneNames_pain)/2+0.1)*300)
        })
      

    #lineplot

    output$lineplot_pain<-renderPlotly({getLineplot(randomGenes,countdata)})
    
    observeEvent(input$geneNames_pain, {
      plots=getLineplot(input$geneNames_pain,countdata)
      output$lineplot_pain<-renderPlotly({plots})})
    

    
    #2. GO Enrichment
    output$gsea_ridgePlot_pain <- DT::renderDT(NULL)
    output$goORATreemap_pain <- DT::renderDT(NULL)
    output$goORA_Cowplot_pain <- DT::renderDT(NULL)
    
    #gsea
    output$down_gsea_ridgePlot_pain<-downloadHandler(
      filename <- function(){paste0("GO-GSEA-Ridge-Plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_Pain, height = input$height_Pain)
        grid.arrange(ridgeplot(gsea_pain,showCategory = 50)+
          theme(text=element_text(size=8)))
        dev.off()
      }
    )
    
    output$gsea_ridgePlot_pain<-renderPlot({ridgeplot(gsea_pain,showCategory = 50)+
        theme(text=element_text(size=8))},
        width = 900,height = 900)
    
    #treemap
    output$down_goORATreemap_pain<-downloadHandler(
      filename <- function(){paste0("GO-ORA-Treemap.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_Pain, height = input$height_Pain)
        treemapPlot(reducedTerms_pain)
        dev.off()
      }
    )
    
    output$goORATreemap_pain<-renderPlot({treemapPlot(reducedTerms_pain)},width = 900,height = 550)
    
    #cowplot
    output$down_goORA_Cowplot_pain<-downloadHandler(
      filename <- function(){paste0("GO-ORA-COW-Plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_Pain, height = input$height_Pain)
        draw(ht)
        dev.off()
      }
    )
    
    output$goORA_Cowplot_pain<-renderPlot({print(ht)},height = 2000)
    enable("makeMeInteractivebtn_pain")
    #make interactive and hide side bar panel
    observeEvent(input$makeMeInteractivebtn_pain, {
      hide(id = "sidePanelGroupCmprsn_pain")
      output$interactiveGoPlot_pain<-renderPlot({
        InteractiveComplexHeatmapModal(input, output, session,ht,cancel_action = "hide")
      })
    })
    
    #Show side bar pane
    observeEvent(input$showSideBarPanel_pain, {
      show(id = "sidePanelGroupCmprsn_pain")})
    
    ####################
    #oab data panel
    ###################

    
    #1. Tables
    output$oabTable <- DT::renderDT(NULL)
    output$oabTable_edger <- DT::renderDT(NULL)
    output$oabTable_edger_degs <- DT::renderDT(NULL)
    
    colnameForEdgerRes=c("log2FoldChange","logCPM","PValue","FDR","description")
    options = list( 
      dom = "Blfrtip"
      , buttons = 
        list(list(
          extend = "csv",
          filename= 'dataTable',
          text='.CSV  <span class="glyphicon glyphicon-download-alt"></span>'
        )))
    
    output$oabTable = DT::renderDataTable({datatable( data = countdata.OAB,
                                                       extensions = 'Buttons',
                                                       options = options)})
    
    output$oabTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_OAB),
                                                             extensions = 'Buttons',
                                                             options = options)})
    
    output$oabTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res_OAB[colnameForEdgerRes],
                                                                  extensions = 'Buttons',
                                                                  options = options)},width=1000)
    
    observeEvent(input$geneNames_OAB, {
      output$oabTable = DT::renderDataTable({datatable( data = countdata.OAB[input$geneNames_OAB,],
                                                         extensions = 'Buttons',
                                                         options = options)})
      
      output$oabTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_OAB)[input$geneNames_OAB,],
                                                               extensions = 'Buttons',
                                                               options = options)
      })
      output$oabTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res_OAB[input$geneNames_OAB,colnameForEdgerRes],
                                                                    extensions = 'Buttons',
                                                                    options = options)
      })
      
    })
    
    restartFN <- observeEvent(input$completeDATA_OAB, {
      output$oabTable = DT::renderDataTable({datatable( data = countdata.OAB,
                                                         extensions = 'Buttons',
                                                         options = options)
      })
      output$oabTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_OAB),
                                                               extensions = 'Buttons',
                                                               options = options)
      })
      output$oabTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res_OAB[colnameForEdgerRes],
                                                                    extensions = 'Buttons',
                                                                    options = options)
      })
      shinyjs::reset("geneNames_OAB_form")
      
    })
    
    #2. Plots
    output$volcano_oab <- DT::renderDT(NULL)
    output$boxplot_OAB <- DT::renderDT(NULL)
    output$lineplot_OAB <- DT::renderDT(NULL)
    
    #volcano
    # Download
    output$down_volcano_oab<-downloadHandler(
      filename <- function(){paste0("volcano-plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_OAB, height = input$height_OAB)
        p=getVolcano_plots_selectedGenes(edgeR_res_all_OAB,input$geneNames_OAB)
        if(p=="error"){grid.arrange(volcanoOAB)}
        else{grid.arrange(p)}
        dev.off()
      }
    )
    
    output$volcano_oab<-renderPlot({volcanoOAB})
    observeEvent(input$geneNames_OAB, {
      output$volcano_oab<-renderPlot({
        p=getVolcano_plots_selectedGenes(edgeR_res_all_OAB,input$geneNames_OAB)
        if(p=="error")
        {
          output$volcano_oab<-renderPlot({volcanoOAB})
          
        }
        else{p}
      })
      
    })
    restartFN <- observeEvent(input$completeDATA_OAB, {
      output$volcano_oab<-renderPlot({volcanoOAB})
    })
    
    #random genes for intial boxplot and line plot
    randomGenes<<-rownames(edgeR_res_all_OAB)[30:31]
    
    #boxplot
    # Download
    output$down_boxplot_OAB<-downloadHandler(
      filename <- function(){paste0("boxplot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_OAB, height = input$height_OAB)
        plots=getBoxplot(input$geneNames_OAB,countdata.OAB,"DO")
        if(plots=="error"){
          p=getBoxplot(randomGenes,countdata.OAB,"DO")
          do.call("grid.arrange", c(p, ncol=2))
        }
        else{do.call("grid.arrange", c(plots, ncol=2))}
        dev.off()
      }
    )
    
    output$boxplot_OAB<-renderPlot({
      p=getBoxplot(randomGenes,countdata.OAB,"DO")
      do.call("grid.arrange", c(p, ncol=2))
    },height = 300)
    
    observeEvent(input$geneNames_OAB, {
      output$boxplot_OAB<-renderPlot({
        plots=getBoxplot(input$geneNames_OAB,countdata.OAB,"DO")
        
        if(plots=="error") #raise error if selectinput is empty
        {
          p=getBoxplot(randomGenes,countdata.OAB,"DO")
          do.call("grid.arrange", c(p, ncol=2))
        }
        else{do.call("grid.arrange", c(plots, ncol=2))}
      },height =round(length(input$geneNames_OAB)/2+0.1)*300)
    })
    
    
    #lineplot
    
    output$lineplot_OAB<-renderPlotly({getLineplot(randomGenes,countdata.OAB)})
    
    observeEvent(input$geneNames_OAB, {
      plots=getLineplot(input$geneNames_OAB,countdata.OAB)
      output$lineplot_OAB<-renderPlotly({plots})})
    
    
    
    #2. GO Enrichment
    output$gsea_ridgePlot_OAB <- DT::renderDT(NULL)
    output$goORATreemap_OAB <- DT::renderDT(NULL)
    output$goORA_Cowplot_OAB <- DT::renderDT(NULL)
    
    #gsea
    output$down_gsea_ridgePlot_OAB<-downloadHandler(
      filename <- function(){paste0("GO-GSEA-Ridge-Plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_OAB, height = input$height_OAB)
        grid.arrange(ridgeplot(gsea_OAB,showCategory = 50)+
                       theme(text=element_text(size=8)))
        dev.off()
      }
    )
    
    output$gsea_ridgePlot_OAB<-renderPlot({ridgeplot(gsea_OAB,showCategory = 50)+
        theme(text=element_text(size=8))},
        width = 900,height = 1800)
    
    #treemap
    output$down_goORATreemap_OAB<-downloadHandler(
      filename <- function(){paste0("GO-ORA-Treemap.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_OAB, height = input$height_OAB)
        treemapPlot(reducedTerms_OAB)
        dev.off()
      }
    )
    
    output$goORATreemap_OAB<-renderPlot({treemapPlot(reducedTerms_OAB)},width = 900,height = 550)
    
    #cowplot
    output$down_goORA_Cowplot_OAB<-downloadHandler(
      filename <- function(){paste0("GO-ORA-COW-Plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_OAB, height = input$height_OAB)
        draw(ht_OAB)
        dev.off()
      }
    )
    
    output$goORA_Cowplot_OAB<-renderPlot({print(ht_OAB)},height = 2000)
    enable("makeMeInteractivebtn_OAB")
    #make interactive and hide side bar panel
    observeEvent(input$makeMeInteractivebtn_OAB, {
      hide(id = "sidePanelGroupCmprsn_OAB")
      output$interactiveGoPlot_OAB<-renderPlot({
        InteractiveComplexHeatmapModal(input, output, session,ht_OAB,cancel_action = "hide")
      })
    })
    
    #Show side bar panel
    observeEvent(input$showSideBarPanel_OAB, {
      show(id = "sidePanelGroupCmprsn_OAB")})
    
    
    ####################
    #AC data
    ###################
    
    
    #1. Tables
    output$ACTable <- DT::renderDT(NULL)
    output$ACTable_edger <- DT::renderDT(NULL)
    output$ACTable_edger_degs <- DT::renderDT(NULL)
    
    colnameForEdgerRes=c("log2FoldChange","logCPM","PValue","FDR","description")
    options = list( 
      dom = "Blfrtip"
      , buttons = 
        list(list(
          extend = "csv",
          filename= 'dataTable',
          text='.CSV  <span class="glyphicon glyphicon-download-alt"></span>'
        )))
    
    output$ACTable = DT::renderDataTable({datatable( data = countdata.AC,
                                                     extensions = 'Buttons',
                                                     options = options)})
    
    output$ACTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_AC),
                                                           extensions = 'Buttons',
                                                           options = options)})
    
    output$ACTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res_AC[colnameForEdgerRes],
                                                                extensions = 'Buttons',
                                                                options = options)},width=1000)
    
    observeEvent(input$geneNames_AC, {
      output$ACTable = DT::renderDataTable({datatable( data = countdata.AC[input$geneNames_AC,],
                                                       extensions = 'Buttons',
                                                       options = options)})
      
      output$ACTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_AC)[input$geneNames_AC,],
                                                             extensions = 'Buttons',
                                                             options = options)
      })
      output$ACTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res_AC[input$geneNames_AC,colnameForEdgerRes],
                                                                  extensions = 'Buttons',
                                                                  options = options)
      })
      
    })
    
    restartFN <- observeEvent(input$completeDATA_AC, {
      output$ACTable = DT::renderDataTable({datatable( data = countdata.AC,
                                                       extensions = 'Buttons',
                                                       options = options)
      })
      output$ACTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_AC),
                                                             extensions = 'Buttons',
                                                             options = options)
      })
      output$ACTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res_AC[colnameForEdgerRes],
                                                                  extensions = 'Buttons',
                                                                  options = options)
      })
      shinyjs::reset("geneNames_AC_form")
      
    })
    
    #2. Plots
    output$volcano_AC <- DT::renderDT(NULL)
    output$boxplot_AC <- DT::renderDT(NULL)
    output$lineplot_AC <- DT::renderDT(NULL)
    
    #volcano
    # Download
    output$down_volcano_AC<-downloadHandler(
      filename <- function(){paste0("volcano-plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_AC, height = input$height_AC)
        p=getVolcano_plots_selectedGenes(edgeR_res_all_AC,input$geneNames_AC)
        if(p=="error"){grid.arrange(volcanoAC)}
        else{grid.arrange(p)}
        dev.off()
      }
    )
    
    output$volcano_AC<-renderPlot({volcanoAC})
    observeEvent(input$geneNames_AC, {
      output$volcano_AC<-renderPlot({
        p=getVolcano_plots_selectedGenes(edgeR_res_all_AC,input$geneNames_AC)
        if(p=="error")
        {
          output$volcano_AC<-renderPlot({volcanoAC})
          
        }
        else{p}
      })
      
    })
    restartFN <- observeEvent(input$completeDATA_AC, {
      output$volcano_AC<-renderPlot({volcanoAC})
    })
    
    #random genes for initial boxplot and line plot
    randomGenes<<-rownames(edgeR_res_all_AC)[30:31]
    
    #boxplot
    # Download
    output$down_boxplot_AC<-downloadHandler(
      filename <- function(){paste0("boxplot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_AC, height = input$height_AC)
        plots=getBoxplot(input$geneNames_AC,countdata.AC,"UA")
        if(plots=="error"){
          p=getBoxplot(randomGenes,countdata.AC,"UA")
          do.call("grid.arrange", c(p, ncol=2))
        }
        else{do.call("grid.arrange", c(plots, ncol=2))}
        dev.off()
      }
    )
    
    output$boxplot_AC<-renderPlot({
      p=getBoxplot(randomGenes,countdata.AC,"UA")
      do.call("grid.arrange", c(p, ncol=2))
    },height = 300)
    
    observeEvent(input$geneNames_AC, {
      output$boxplot_AC<-renderPlot({
        plots=getBoxplot(input$geneNames_AC,countdata.AC,"UA")
        
        if(plots=="error") #raise error if selectinput is empty
        {
          p=getBoxplot(randomGenes,countdata.AC,"UA")
          do.call("grid.arrange", c(p, ncol=2))
        }
        else{do.call("grid.arrange", c(plots, ncol=2))}
      },height =round(length(input$geneNames_AC)/2+0.1)*300)
    })
    
    
    #lineplot
    
    output$lineplot_AC<-renderPlotly({getLineplot(randomGenes,countdata.AC)})
    
    observeEvent(input$geneNames_AC, {
      plots=getLineplot(input$geneNames_AC,countdata.AC)
      output$lineplot_AC<-renderPlotly({plots})})
    
    
    
    #2. GO Enrichment
    output$gsea_ridgePlot_AC <- DT::renderDT(NULL)
    output$goORATreemap_AC <- DT::renderDT(NULL)
    output$goORA_Cowplot_AC <- DT::renderDT(NULL)
    
    #gsea
    output$down_gsea_ridgePlot_AC<-downloadHandler(
      filename <- function(){paste0("GO-GSEA-Ridge-Plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_AC, height = input$height_AC)
        grid.arrange(ridgeplot(gsea_AC,showCategory = 50)+
                       theme(text=element_text(size=8))+
                       ggtitle("Top 50: visit Go Enrichment Exploration
 panel for more"))
        dev.off()
      }
    )
    
    output$gsea_ridgePlot_AC<-renderPlot({ridgeplot(gsea_AC,showCategory = 50)+
        theme(text=element_text(size=8))+
        ggtitle("Top 50: visit Go Enrichment Exploration
 panel for more")},
        width = 900,height = 1800)
    
    #treemap
    output$down_goORATreemap_AC<-downloadHandler(
      filename <- function(){paste0("GO-ORA-Treemap.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_AC, height = input$height_AC)
        treemapPlot(reducedTerms_AC)
        dev.off()
      }
    )
    
    output$goORATreemap_AC<-renderPlot({treemapPlot(reducedTerms_AC)},width = 900,height = 550)
    
    #cowplot
    output$down_goORA_Cowplot_AC<-downloadHandler(
      filename <- function(){paste0("GO-ORA-COW-Plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_AC, height = input$height_AC)
        draw(ht_AC)
        dev.off()
      }
    )
    
    output$goORA_Cowplot_AC<-renderPlot({print(ht_AC)},height = 2000)
    enable("makeMeInteractivebtn_AC")
    #make interactive and hide side bar panel
    observeEvent(input$makeMeInteractivebtn_AC, {
      hide(id = "sidePanelGroupCmprsn_AC")
      output$interactiveGoPlot_AC<-renderPlot({
        InteractiveComplexHeatmapModal(input, output, session,ht_AC,cancel_action = "hide")
      })
    })
    
    #Show side bar panel
    observeEvent(input$showSideBarPanel_AC, {
      show(id = "sidePanelGroupCmprsn_AC")})
    
    
    ####################
    #BOO data
    ###################
    
    
    #1. Tables
    output$BOOTable <- DT::renderDT(NULL)
    output$BOOTable_edger <- DT::renderDT(NULL)
    output$BOOTable_edger_degs <- DT::renderDT(NULL)
    
    colnameForEdgerRes=c("log2FoldChange","logCPM","PValue","FDR","description")
    options = list( 
      dom = "Blfrtip"
      , buttons = 
        list(list(
          extend = "csv",
          filename= 'dataTable',
          text='.CSV  <span class="glyphicon glyphicon-download-alt"></span>'
        )))
    
    output$BOOTable = DT::renderDataTable({datatable( data = countdata.BOO,
                                                      extensions = 'Buttons',
                                                      options = options)})
    
    output$BOOTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_BOO),
                                                            extensions = 'Buttons',
                                                            options = options)})
    
    output$BOOTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res_BOO[colnameForEdgerRes],
                                                                 extensions = 'Buttons',
                                                                 options = options)},width=1000)
    
    observeEvent(input$geneNames_BOO, {
      output$BOOTable = DT::renderDataTable({datatable( data = countdata.BOO[input$geneNames_BOO,],
                                                        extensions = 'Buttons',
                                                        options = options)})
      
      output$BOOTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_BOO)[input$geneNames_BOO,],
                                                              extensions = 'Buttons',
                                                              options = options)
      })
      output$BOOTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res_BOO[input$geneNames_BOO,colnameForEdgerRes],
                                                                   extensions = 'Buttons',
                                                                   options = options)
      })
      
    })
    
    restartFN <- observeEvent(input$completeDATA_BOO, {
      output$BOOTable = DT::renderDataTable({datatable( data = countdata.BOO,
                                                        extensions = 'Buttons',
                                                        options = options)
      })
      output$BOOTable_edger = DT::renderDataTable({datatable( data = as.data.frame(edgeR_res_all_BOO),
                                                              extensions = 'Buttons',
                                                              options = options)
      })
      output$BOOTable_edger_degs = DT::renderDataTable({datatable( data = edgeR_res_BOO[colnameForEdgerRes],
                                                                   extensions = 'Buttons',
                                                                   options = options)
      })
      shinyjs::reset("geneNames_BOO_form")
      
    })
    
    #2. Plots
    output$volcano_BOO <- DT::renderDT(NULL)
    output$boxplot_BOO <- DT::renderDT(NULL)
    output$lineplot_BOO <- DT::renderDT(NULL)
    
    #volcano
    # Download
    output$down_volcano_BOO<-downloadHandler(
      filename <- function(){paste0("volcano-plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_BOO, height = input$height_BOO)
        p=getVolcano_plots_selectedGenes(edgeR_res_all_BOO,input$geneNames_BOO)
        if(p=="error"){grid.arrange(volcanoBOO)}
        else{grid.arrange(p)}
        dev.off()
      }
    )
    
    output$volcano_BOO<-renderPlot({volcanoBOO})
    observeEvent(input$geneNames_BOO, {
      output$volcano_BOO<-renderPlot({
        p=getVolcano_plots_selectedGenes(edgeR_res_all_BOO,input$geneNames_BOO)
        if(p=="error")
        {
          output$volcano_BOO<-renderPlot({volcanoBOO})
          
        }
        else{p}
      })
      
    })
    restartFN <- observeEvent(input$completeDATA_BOO, {
      output$volcano_BOO<-renderPlot({volcanoBOO})
    })
    
    #random genes for bocplot and line plot
    randomGenes<<-rownames(edgeR_res_all_BOO)[30:31]
    
    #boxplot
    # Download
    output$down_boxplot_BOO<-downloadHandler(
      filename <- function(){paste0("boxplot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_BOO, height = input$height_BOO)
        plots=getBoxplot(input$geneNames_BOO,countdata.BOO,"BOO")
        if(plots=="error"){
          p=getBoxplot(randomGenes,countdata.BOO,"BOO")
          do.call("grid.arrange", c(p, ncol=2))
        }
        else{do.call("grid.arrange", c(plots, ncol=2))}
        dev.off()
      }
    )
    
    output$boxplot_BOO<-renderPlot({
      p=getBoxplot(randomGenes,countdata.BOO,"BOO")
      do.call("grid.arrange", c(p, ncol=2))
    },height = 300)
    
    observeEvent(input$geneNames_BOO, {
      output$boxplot_BOO<-renderPlot({
        plots=getBoxplot(input$geneNames_BOO,countdata.BOO,"BOO")
        
        if(plots=="error") #raise error if selectinput is empty
        {
          p=getBoxplot(randomGenes,countdata.BOO,"BOO")
          do.call("grid.arrange", c(p, ncol=2))
        }
        else{do.call("grid.arrange", c(plots, ncol=2))}
      },height =round(length(input$geneNames_BOO)/2+0.1)*300)
    })
    
    
    #lineplot
    
    output$lineplot_BOO<-renderPlotly({getLineplot(randomGenes,countdata.BOO)})
    
    observeEvent(input$geneNames_BOO, {
      plots=getLineplot(input$geneNames_BOO,countdata.BOO)
      output$lineplot_BOO<-renderPlotly({plots})})
    
    
    
    #2. GO Enrichment
    output$gsea_ridgePlot_BOO <- DT::renderDT(NULL)
    output$goORATreemap_BOO <- DT::renderDT(NULL)
    output$goORA_Cowplot_BOO <- DT::renderDT(NULL)
    
    #gsea
    output$down_gsea_ridgePlot_BOO<-downloadHandler(
      filename <- function(){paste0("GO-GSEA-Ridge-Plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_BOO, height = input$height_BOO)
        grid.arrange(ridgeplot(gsea_BOO,showCategory = 50)+
                       theme(text=element_text(size=8))+
                       ggtitle("Top 50: visit Go Enrichment Exploration
 panel for more"))
        dev.off()
      }
    )
    
    output$gsea_ridgePlot_BOO<-renderPlot({ridgeplot(gsea_BOO,showCategory = 50)+
        theme(text=element_text(size=8))+
        ggtitle("Top 50: visit Go Enrichment Exploration
 panel for more")},
        width = 900,height = 1800)
    
    #treemap
    output$down_goORATreemap_BOO<-downloadHandler(
      filename <- function(){paste0("GO-ORA-Treemap.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_BOO, height = input$height_BOO)
        treemapPlot(reducedTerms_BOO)
        dev.off()
      }
    )
    
    output$goORATreemap_BOO<-renderPlot({treemapPlot(reducedTerms_BOO)},width = 900,height = 550)
    
    #cowplot
    output$down_goORA_Cowplot_BOO<-downloadHandler(
      filename <- function(){paste0("GO-ORA-COW-Plot.pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_BOO, height = input$height_BOO)
        draw(ht_BOO)
        dev.off()
      }
    )
    
    output$goORA_Cowplot_BOO<-renderPlot({print(ht_BOO)},height = 2000)
    enable("makeMeInteractivebtn_BOO")
    #make interactive and hide side bar panel
    observeEvent(input$makeMeInteractivebtn_BOO, {
      hide(id = "sidePanelGroupCmprsn_BOO")
      output$interactiveGoPlot_BOO<-renderPlot({
        InteractiveComplexHeatmapModal(input, output, session,ht_BOO,cancel_action = "hide")
      })
    })
    
    #Show side bar panel
    observeEvent(input$showSideBarPanel_BOO, {
      show(id = "sidePanelGroupCmprsn_BOO")})
    
    
    ####################
    #GO Exploratory analysis
    ###################
    output$volcano_plot_ForGOTerm <- DT::renderDT(NULL)
    output$table_ForGOTerm <- DT::renderDT(NULL)
    
    #update drop down list of GO terms based on group name
    observeEvent(input$groupName,{
      if(input$groupName=="BPS vs Control")
      {
        #check enrichment type
        if(input$enrichmentType=="ORA")
        {updateSelectizeInput(session, 'goTermName',choices = sort(enrichmentResults_pain[["go_over"]]$Description), server = TRUE) }
        else{updateSelectizeInput(session, 'goTermName',choices = sort(gsea_pain$Description), server = TRUE)}
      }
      
      else if(input$groupName=="DO vs Control"){
        #check enrichment type
        if(input$enrichmentType=="ORA")
        {updateSelectizeInput(session, 'goTermName',choices = sort(enrichmentResults_OAB[["go_over"]]$Description), server = TRUE) }
        else{updateSelectizeInput(session, 'goTermName',choices = sort(gsea_OAB$Description), server = TRUE)}
      }
      
      else if(input$groupName=="UA vs Control"){
        #check enrichment type
        if(input$enrichmentType=="ORA")
        {updateSelectizeInput(session, 'goTermName',choices = sort(enrichmentResults_AC[["go_over"]]$Description), server = TRUE) }
        else{updateSelectizeInput(session, 'goTermName',choices = sort(gsea_AC$Description), server = TRUE)}
      }
      
      else{
        #check enrichment type
        if(input$enrichmentType=="ORA")
        {updateSelectizeInput(session, 'goTermName',choices = sort(enrichmentResults_BOO[["go_over"]]$Description), server = TRUE) }
        else{updateSelectizeInput(session, 'goTermName',choices = sort(gsea_BOO$Description), server = TRUE)}
      }
      
    })
    
    #update dropdown list of GO terms based on enrichment type
    observeEvent(input$enrichmentType,{
      if(input$enrichmentType=="ORA")
      {
        #check group type
        if(input$groupName=="BPS vs Control")
          {updateSelectizeInput(session, 'goTermName',choices = sort(enrichmentResults_pain[["go_over"]]$Description), server = TRUE) }
        else if(input$groupName=="DO vs Control")
          {updateSelectizeInput(session, 'goTermName',choices = sort(enrichmentResults_OAB[["go_over"]]$Description), server = TRUE)}
        else if(input$groupName=="UA vs Control")
        {updateSelectizeInput(session, 'goTermName',choices = sort(enrichmentResults_AC[["go_over"]]$Description), server = TRUE)}
        
        else{updateSelectizeInput(session, 'goTermName',choices = sort(enrichmentResults_BOO[["go_over"]]$Description), server = TRUE)}
        
      }
      
      else{
        #check group type
        if(input$groupName=="BPS vs Control")
        {updateSelectizeInput(session, 'goTermName',choices = sort(gsea_pain$Description), server = TRUE) }
        
        else if(input$groupName=="DO vs Control")
          {updateSelectizeInput(session, 'goTermName',choices = sort(gsea_OAB$Description), server = TRUE)}
        
        else if(input$groupName=="UA vs Control")
        {updateSelectizeInput(session, 'goTermName',choices = sort(gsea_AC$Description), server = TRUE)}
        
        else 
        {updateSelectizeInput(session, 'goTermName',choices = sort(gsea_BOO$Description), server = TRUE)}
        
        
        }
    })
    
    #download
    output$down_volcanoGoTerm<-downloadHandler(
      filename <- function(){paste(input$goTermName,"in",input$groupName,".pdf")},
      content = function(file){
        pdf(file = file,onefile=TRUE, width = input$width_OAB, height = input$height_OAB)
        grid.arrange(getVolcanoFor_goTerm(input$groupName,input$enrichmentType,input$goTermName))
        dev.off()
      }
    )
    
    #plot volcano
    restartFN <- observeEvent(input$plotVolcanoFor_goTerm, {
      output$volcano_plot_ForGOTerm<-renderPlot({
        getVolcanoFor_goTerm(input$groupName,input$enrichmentType,input$goTermName)
      })
      
    output$table_ForGOTerm = DT::renderDataTable({datatable( data = getTable_goTerm(input$groupName,input$enrichmentType,input$goTermName),
                                                       extensions = 'Buttons',
                                                       caption = "Data", 
                                                       options = options)})
    })
    

    
    
}