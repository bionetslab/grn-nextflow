#!/usr/bin/env Rscript

################################################################################
### Imports
# Order of imports is important here.
# Seurat must be imported before networkD3 and htmlwidgets.
# Otherwise, the JS function of these libraries will be masked and this leads to errors
# renv::activate('/home/bionets-og86asub/Documents/netmap/')
library(Seurat)
library(shiny)
library(ggplot2)
library(optparse)
library(shinyWidgets)
library(DT)
library(data.table)
library(shinythemes)
library(shinyhelper)
library(dplyr)
library(stringr)
library(cowplot)
library(igraph)
library(cyjShiny)
library(htmlwidgets)
library(graph)
library(jsonlite)
require(dplyr)
require(tidyr)
require(rcartocolor)
require(colorspace)
require(presto)
################################################################################

source('style_generator.R')
source('auxiliary_functions.R')

option_list <- list(
  make_option(
    c("-r", "--results.path"),
    type = "character",
    default = "",
    help = "Path to results"
  ),
  make_option(
    c("-s", "--selection"),
    type = "character",
    default = "",
    help = "Configurated selection"
  ),
  make_option(
    c("-f", "--seurat.file"),
    type = "character",
    default = "",
    help = "Seurat data file"
  ),
  make_option(
    c("--dgrntools"),
    type = "character",
    default = "No tools were chosen for DGRN Inference",
    help = "used diffGRN tools in the pipeline"
  ),
  make_option(
    c("--grntools"),
    type = "character",
    default = "No tools were chosen for GRN Inference",
    help = "used GRN tools in the pipeline"
  ),
  make_option(
    c("-m", "--mode"),
    type = "character",
    default = "seurat",
    help = "Specifies which file type was used as an input for the pipeline"
  ),
  make_option(
    c("-n", "--n.samples"),
    type = "integer",
    default = 100,
    help = "Number of meta cells used in the pipeline (Needed for comparison plot)"
  ),
  make_option(
    c("--p.missing"),
    type = "integer",
    default = 10,
    help = "Percentage of 0 allowed per gene",
    metavar = "number"
  ),
  make_option(
    c("--metacells"),
    default = TRUE,
    action = 'store_true',
    help = "Metacell aggreation has been used"
  )
)
# changing opt$grntoools and opt$dgrntools to be suitable as input for the Shiny app
opt <- parse_args(OptionParser(option_list = option_list))
metacells <- opt$metacells
if (opt$grntools != "No tools were chosen for GRN Inference") {
  opt$grntools <- c("NA", opt$grntools)
}

metacells<-FALSE
opt$results.path <-
  "/home/bionets-og86asub/Documents/external_analyses/huiqin/"
opt$selection <-"NSCLC_BRCA,RNA,cluster_type,NSCLC_M0,NSCLC_M0-NSCLC_BRCA,RNA,cluster_type,BRCA_M0,BRCA_M0-NSCLC_PRAD,RNA,cluster_type,NSCLC_M0,NSCLC_M0-NSCLC_PRAD,RNA,cluster_type,PRAD_M0,BRCA_M0-BRCA_PRAD,RNA,cluster_type,BRCA_M0,BRCA_M0-BRCA_PRAD,RNA,cluster_type,PRAD_M0,PRAD_M0"
opt$seurat.file <-
  "/home/bionets-og86asub/Documents/external_analyses/huiqin/cancer.rds"
opt$dgrntools <- "boostdiff,zscores,diffcoex"
opt$dgrntools <- c(strsplit(opt$dgrntools, ",")[[1]])
opt$grntools <- "grnboost2"
# 

# opt$results.path <-
#   "/home/bionets-og86asub/Documents/netmap/data/misc/"
# opt$selection <-"C5_C8,SCT,cluster:genotype,C8_WT,8:wt-C5_C8,SCT,cluster:genotype,C8_KO,8:ko"
# opt$seurat.file <-
#   "/home/bionets-og86asub/Documents/netmap/data/misc/ko_vs_wt.rds"
# opt$dgrntools <- "boostdiff,zscores,diffcoex"
# opt$dgrntools <- c(strsplit(opt$dgrntools, ",")[[1]])
# opt$grntools <- "grnboost2"
# opt$assay<-'SCT'


# parse selection
selections <-sapply(strsplit(opt$selection, "-")[[1]], function(x) {strsplit(x, ",")[[1]]})
selections <- as.data.table(rbind(t(selections)))
colnames(selections) <- c("key", "assay", "group_var", "condition", "variables")

all_networks<-read_all_networks(selections, opt$results.path)
color.map<-make_color_map(all_networks, selections)


  

########### Loading Seurat object, filtering the correct cells and performing differential testing
adata <- readRDS(opt$seurat.file)
# saves a bit of computation time
gene_names <- rownames(adata)
capitalized_gene_names <- toupper(gene_names)
# adata <- adata$all

group.var <- selections$group_var[1]
merge_cols<-str_split(group.var, ':')[[1]]
df <- as.data.table(adata@meta.data) %>% 
  unite(x, merge_cols, sep = ":", remove = FALSE)
adata@meta.data[[group.var]]<-df$x
  
if (metacells){
metacell.seurat<-read_metacell_files(selections)
print(metacell.seurat)
}
  
Idents(adata) <- group.var

<<<<<<< HEAD
create_graph <- function(all_networks, de.genes = NULL) {
  nodes <-
    data.table(
	      id = c(all_networks$target, all_networks$source),
      tool = c(all_networks$tool, all_networks$tool),
      keys = c(all_networks$key, all_networks$key)
    )
=======
# get all factors and columns that were used in the selection but were not of class factor
factors <- unique(c(names(which(sapply(adata@meta.data, class)=='factor')), selections$group_var[1]))
factors<-intersect(factors, colnames(adata@meta.data))
>>>>>>> c903c443582958ff11f589676d4faa83c76c17b7

Idents(adata) <- group.var
de.genes <- FindAllMarkers(adata)
de.genes$name <- as.character(de.genes$cluster)
colnames(de.genes)[2] <- "value"
#de.genes<-NULL

<<<<<<< HEAD
read_files <- function(file_list, key) {
  all_data <- list()
  
  for (file in file_list) {
    tryCatch(
      {
        data <- fread(file)
        if (!is.null(data) && nrow(data) > 0) {
          data$tool <- basename(dirname(file))
          data$key <- basename(dirname(dirname(file)))
          all_data[[basename(dirname(file))]] <- data
        } else {
          print("Data empty")
        }
      },
      error = function(e) {
        cat("Error reading", file, ": ", conditionMessage(e), "\n")
      }
    )
  }
  all_data <- rbindlist(all_data)
  return(all_data)
}

make_color_map<-function(all_networks){
  print('Generating colors')
  n.colors<-nrow(unique(all_networks[, .(interaction)]))
  color.pal<-rcartocolor::carto_pal(n=max(3, 2*n.colors), 'Prism')
  color.pal<-color.pal[seq(1, 2*n.colors, 2)]
 
  color.df<-data.table(color.pal[1:n.colors], unique(all_networks[, .(interaction)]))
  colnames(color.df)<-c('color', 'interaction')
  color.df<-color.df %>% pivot_longer(-c('interaction'))
  print(color.df)
  combos<-unique(all_networks[, .(interaction, tool)])
  print(combos)
  #combos<-unique(cross_join(combos, combos)[,.(interaction.x,tool.y)])
  colnames(combos)<-c('interaction', 'tool')
  print(combos)
  combos<-merge(combos, color.df, by = 'interaction', allow.cartesian=T)
  combos<-combos %>% group_by(interaction) %>% mutate(row_count = row_number()) %>% as.data.table()
  combos$value<-sapply(1:nrow(combos), function(i) lighten(combos$value[i], (combos$row_count[i]-1)*(1/max(combos$row_count)))) 
  combos<-as.data.table(combos)
  combos<-combos[, .(tool, interaction, value)]
  
  cl<-list()
  possible_styles<-c('solid', 'dotted', 'dashed')
  dashes<- ceiling(sqrt(length(unique(all_networks$key))-1))
  #dash.styles<-c('NA', 'NA',as.character(sapply(1:dashes, function(x) sapply(1:3, function(y) paste0('[', 2*x, ', ', 2*y, ']')))))
  #dash.width<- c('NA', 'NA', sapply(1:dashes, function(x) sapply(1:3, function(y) 2*x)))
  #dash.gap<- c('NA', 'NA', sapply(1:dashes, function(x) sapply(1:3, function(y) 2*y)))
  possible_styles<-c('solid', 'dotted', rep('dashed', dashes*dashes))

  for (i in 1:length(unique(all_networks$key))){
    c<-combos
    c$arrow_style<-possible_styles[i]
    c$key<-unique(all_networks$key)[i]
    cl[[i]]<-c
  }
  combos<-rbindlist(cl)
  print(combos)
  return(combos)
}

get_color<-function(color.df, t){
  print(color.df)
  print(t)
  values<- unique(color.df[(tool==t)]$value)
  print(values)
  names(values)<-unique(color.df[(tool==t)]$interaction)
  return(values)
  }




# parsing the selection into a configuration data.table that has all necessary information to enable all features of the shiny app.
# every row corresponds to ONE comparison in the input configuration file (or the one comparison of tsv mode)
configuration <-
  data.table(
    key = NULL,
    assay = NULL,
    group_var = NULL,
    condition_names = NULL,
    tools = NULL,
    network_files = NULL
  )

if (opt$mode == "tsv") {
  # get tools
  tools <-
    list.dirs(
      path = file.path(opt$results.path, opt$selection),
      full.names = FALSE,
      recursive = FALSE
    )
  tool_names <- paste(tools, collapse = ",")
  
  # get network files
  all_files <-
    list.files(
      path = file.path(opt$results.path, opt$selection),
      full.names = TRUE,
      recursive = TRUE
    )
  
  network_files <-
    c(grep("aggregated_filtered_network", all_files, value = TRUE))
  network_files <- paste(network_files, collapse = ",")
  
  # get condition names
  files <-
    list.files(
      path = file.path(opt$results.path, opt$selection),
      full.names = FALSE,
      recursive = FALSE
    )
  cond_files <- c(grep("out_", files, value = TRUE))
  cond_name1 <- gsub("out_", "", cond_files[1])
  cond_name1 <- gsub(".tsv", "", cond_name1)
  cond_name2 <- gsub("out_", "", cond_files[2])
  cond_name2 <- gsub(".tsv", "", cond_name2)
  cond_names <- paste(c(cond_name1, cond_name2), collapse = ",")
  
  configuration <- rbind(
    configuration,
    list(
      key = opt$selection,
      assay = "RNA",
      group_var = NA,
      condition_names = cond_names,
      tools = tool_names,
      network_files = network_files
    )
  )
} else {
  # parse selection
  selections <-
    sapply(strsplit(opt$selection, "-")[[1]], function(x) {
      strsplit(x, ",")[[1]]
    })
  selections <- as.data.table(rbind(t(selections)))
  colnames(selections) <-
    c("key", "assay", "group_var", "condition", "condition2")
  
  an <- list()
  for (k in unique(selections$key)) {
    all_files <-
      list.files(
        path = file.path(opt$results.path, k),
        full.names = TRUE,
        recursive = TRUE
      )
    network_files <-
      c(grep("aggregated_filtered_network", all_files, value = TRUE))
    all_networks <- read_files(network_files, k)
    an[[k]] <- all_networks
  }
  all_networks <- rbindlist(an)
  
    colnames(all_networks) <-
    c(
      "target",
      "source",
      "weight",
      "interaction",
      "effect",
      "tool",
      "key"
    )
    
  color.map<-make_color_map(all_networks)

# get all factors and columns that were used in the selection but were not of class factor
factors <- unique(c(names(which(sapply(adata@meta.data, class)=='factor')), selections$group_var[1]))
factors<-intersect(factors, colnames(adata@meta.data))

Idents(adata) <- group.var
de.genes <- FindAllMarkers(adata)
de.genes$name <- as.character(de.genes$cluster)
colnames(de.genes)[2] <- "value"
#de.genes<-NULL

  # 
graph <- create_graph(all_networks, de.genes)
nodes <- graph[[1]]
edges <- graph[[2]]

# SET INPUT OPTIONS ----
# Generate the style JSON
style_json <- generate_style_json(color.map$tool, color.map$value, color.map$condition, color.map$key, color.map$arrow_style, color.map$dash.style)
# Write the style JSON to a file
writeLines(style_json, "style.js")
basicStyleFile <- "style.js"
styleList <- c("", "Basic" = "basicStyleFile")


# UI ----
ui <- shinyUI(fluidPage(
  tags$head(
    tags$link(
      rel = "stylesheet", type = "text/css",
      href = "http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css"
    ),
    tags$style("#cyjShiny{height:95vh !important;}")
  ),
  navbarPage(
    "Internet Xplorer",
    tabPanel(
      "net",
      sidebarLayout(
        sidebarPanel(
          fluidPage(fluidRow(
            column(
              12,
              fluidRow(
                selectInput(
                  "Key_pick",
                  "Choose key to display DiffGRN/GRN for:",
                  choices = unique(selections$key),
                  selected = unique(selections$key)[1],
                  multiple = TRUE,
                  selectize = TRUE,
                  width = NULL,
                  size = NULL
                ),
                selectInput(
                  "DiffGRN_pick",
                  "Choose DiffGRN to display:",
                  choices = unique(all_networks$tool),
                  selected = unique(all_networks$tool)[1],
                  multiple = TRUE,
                  selectize = TRUE,
                  width = NULL,
                  size = NULL
                ),
                selectInput('selectCondition', 
                            'Select Condition',
                            choices = unique(c(all_networks$interaction, selections$condition)),
                            selected = unique(selections$condition)[1],
                            multiple = TRUE,
                            selectize = TRUE,
                            width = NULL,
                            size = NULL)
              )
            )
          )),
          sliderInput('n.edges', 'Number of edges', 1, max(edges[, .N, by=c('tool', 'key')]$N),  max(edges[, .N, by=c('tool', 'key')]$N), step = 1, round = FALSE),
          
          selectInput(
            "doLayout",
            "Select Layout:",
            choices = c(
              "",
              "cose",
              "cola",
              "circle",
              "concentric",
              "breadthfirst",
              "grid",
              "random",
              "preset",
              "fcose"
            )
          ),
          selectInput(
            "setNodeAttributes",
            "Select Node Attribute:",
            choices = c('Select condition', unique(nodes$name)),
            selected = 'Select condition'
          ),
          actionButton('generate_network', 'Generate Network'),
          HTML("<br>"),
          plotOutput('legend'),
          uiOutput('arrow_legend'),
          
          HTML("<br>"),
          actionButton("sfn", "Select First Neighbor"),
          actionButton("fit", "Fit Graph"),
          actionButton("hideSelection", "Hide Selection"),
          actionButton("showOnlySelection", 'Show selection'),
          actionButton("fitSelected", "Fit Selected"),
          actionButton("clearSelection", "Unselect Nodes"),
          HTML("<br>"),
          # actionButton("loopConditions", "Loop Conditions"),
          # HTML("<br>"),
          actionButton("getSelectedNodes", "Get Selected Nodes"),
          HTML("<br>"),
          htmlOutput("selectedNodesDisplay"),
          HTML("<br>"),
          width = 2,
          selectInput(
            "select_genes",
            "Select Genes to display",
            character(0),
            multiple = TRUE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          ),
          actionButton("redrawGraph", "Draw Graph with selected nodes"),
        ),
        mainPanel(
          width = 10,
          fluidRow(
            column(8, cyjShinyOutput("cyjShiny")),
            column(
              4,
              fluidRow(
                plotOutput("violin_plot", height = "50vh"),
                downloadButton("downloadGRNViolinPlot", "Download Plot"),
                 column(6, uiOutput("geneLink")
                )
              ),
              fluidRow(
                plotOutput("linear_model_plot", height = "50vh"),
                downloadButton("downloadGRNLinearModelPlot", "Download Plot")
                )
              ),
              
              ),
              
            ),
          )
    ),
    tabPanel(
      "Gene Expression",
      h2("Gene Expression"),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "select_genes",
            "Select Genes to display",
            character(0),
            multiple = TRUE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          ),
          selectInput(
            "select_primary_grouping",
            "Select primary grouping variable",
            factors,
            selected = factors[1],
            multiple = FALSE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          ),
          selectInput(
            "select_secondary_grouping",
            "Select secondary grouping variable",
            factors,
            selected = factors[1],
            multiple = FALSE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          ),
          actionButton("plot_button", "Generate Plots")
        ),
        mainPanel(tabsetPanel(
          type = "tabs",
          tabPanel(
            "Umap",
            fluidRow(
              column(12, plotOutput("umap_plot")),
              column(12, plotOutput("standard_violin_plot"))
            )
          )
          # ,tabPanel(
          #   "Dotplot",
          #   fluidRow(column(
          #     12, plotOutput("dot_plot")
          #   ))
          # )
        ))
      )
    ) # sidebarLayout
  )
))

# SERVER ----
server <- function(input, output, session) {
  
  color.map<-reactiveVal(color.map)
  
  # Event observers
  observeEvent(input$fit, {
    fit(session, 80)
  })
  
  update_node_vals<-function(){
    attribute <- "value"
    if(input$setNodeAttributes == 'Select condition'){

      setNodeAttributes(
        session,
        attributeName = attribute,
        nodes = nodes[tool %in% input$DiffGRN_pick]$id,
        values = rep(0, length(nodes[tool %in% input$DiffGRN_pick]$id))
      )
    }
    else{
    if('value' %in% colnames(nodes)){
      updated.nodes <-
        nodes[tool %in% input$DiffGRN_pick &
                cluster %in% input$setNodeAttributes]$id
      updated.values <-
        nodes[tool %in% input$DiffGRN_pick &
                cluster %in% input$setNodeAttributes]$value
      difference <-
        setdiff(unique(nodes$id), nodes[tool %in% input$DiffGRN_pick &
                                          cluster %in% input$setNodeAttributes]$id)
      
      updated.nodes <- c(updated.nodes, difference)
      updated.values <- c(updated.values, rep(0, length(difference)))
      updated.values[is.infinite(updated.values)] <- 0
      updated.values[is.na(updated.values)]<-0
      setNodeAttributes(
        session,
        attributeName = attribute,
        nodes = updated.nodes,
        values = updated.values
      )
    }
    }
  }
  

  render_graph<-function(genes = NULL){
    
    k <- isolate(input$Key_pick)
    d <- isolate(input$DiffGRN_pick)
    c <- isolate(input$selectCondition)
    n.edges <- isolate(input$n.edges)
    print(c)
    edges <- edges %>% group_by(tool, key) %>% slice_max(n = n.edges, weight) %>% as.data.table()
    edges<-edges[(tool %in% d) & (key %in% k) & (interaction %in% c)]
    if (nrow(edges)<=0){
      return(NULL)
    }
    color.map(make_color_map(edges, selections))
    # Generate the style JSON
    style_json <- generate_style_json(color.map()$tool, color.map()$value, color.map()$interaction, color.map()$key, color.map()$arrow_style, color.map()$dash.style)
    # Write the style JSON to a file
    writeLines(style_json, "style.js")
    loadStyleFile(basicStyleFile)
    
    if (!is.null(genes)){
    edges<-edges[ (source %in% genes) & (target %in% genes)]
    }

    nodes<-nodes[id %in% c(edges$target, edges$source)]
    nodes<-nodes[(tool %in% d) & (key %in% k)]

    # Did not manage to init graph otherwise
    n_t<-nodes
    n_t$value<-0
    
    graph.json <-
      dataFramesToJSON(
        tbl.nodes = n_t,
        tbl.edges = edges
      )
    
    
    output$cyjShiny <- renderCyjShiny({
      cyjShiny(graph.json,
               layoutName = "cola",
               styleFile = basicStyleFile
      )
    })

    output$legend<-renderPlot({
      cm<-color.map()[tool %in% d]
      cm<-cm[tool %in% opt$dgrntools]
      gp<-ggplot(cm, aes(x = interaction, y = tool, fill = value)) +
        geom_tile() +
        scale_fill_identity() +
        labs(x = "", y = "Tool") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20), axis.text.y =element_text(size = 20, hjust = 0.5) ) +
        guides(fill = FALSE)  # Hide the legend + facet_wrap(~key)
      cm<-color.map()[tool %in% d]
      
      cm<-cm[tool %in% opt$grntools]
      rel_height <- length(intersect(d, opt$grntools))/length(intersect(c(opt$dgrntools, opt$grntools), d))
      print(rel_height)
      gp2<-ggplot(cm, aes(x = interaction, y = tool, fill = value)) +
        geom_tile() +
        xlab("")+ylab("")+
        scale_fill_identity() +
        theme_minimal() +
        theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.text.y =element_text(size = 20, hjust = 0.5)) +
        guides(fill = FALSE)  # Hide the legend + facet_wrap(~key)
      plot_grid(gp, gp2, nrow = 2, rel_heights = c(1-rel_height, rel_height))
        
    })
    
    output$arrow_legend<-renderUI({
      d <- isolate(input$DiffGRN_pick)
      styles<- unique(color.map()[, .(key, arrow_style)])
      tl_dot<-lapply(1:12, function(x) htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 2px; margin-right: 0px; margin-bottom: 6px", 'grey')))
      tl_dash<-lapply(1:4, function(x) htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", 'grey')))
      tl_sol<-list(htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 69px; margin-right: 5px; margin-bottom: 2px", 'grey')))

      
      
      tl<-list()
      for (s in 1:nrow(styles)){
        if(styles$arrow_style[s]=='solid'){
          tl[[s]]<-list(tl_sol, htmltools::span(styles$key[s]), htmltools::br())
        }
        if(styles$arrow_style[s]=='dotted'){
          tl[[s]]<-list(tl_dot, htmltools::span(styles$key[s]), htmltools::br())
        }
        if(styles$arrow_style[s]=='dashed'){
          tl[[s]]<-list(tl_dash, htmltools::span(styles$key[s]), htmltools::br())
        }
      }
      
    print(tagList(tl))
    })
    

  }
  
  
  observeEvent(input$generate_network, {
    render_graph()

  })
  
  

  
      
  observeEvent(input$redrawGraph, ignoreInit = TRUE, {
    genes<-isolate(input$select_genes)
    if (length(genes)>0){
    render_graph(genes)
    update_node_vals()
    }
    else{
      render_graph()
      update_node_vals()
      
    }
  })
  

  
  observeEvent(input$setNodeAttributes, {
    update_node_vals()
  })
  
  
  
  observeEvent(input$doLayout, ignoreInit = TRUE, {
    strategy <- input$doLayout
    doLayout(session, strategy)
    # session$sendCustomMessage(type="doLayout", message=list(input$doLayout))
  })
  
  
  observeEvent(input$sfn, ignoreInit = TRUE, {
    session$sendCustomMessage(type = "sfn", message = list())
  })
  
  observeEvent(input$hideSelection, ignoreInit = TRUE, {
    session$sendCustomMessage(type = "hideSelection", message = list())
  })
  
  observeEvent(input$showOnlySelection, ignoreInit = TRUE, {
    session$sendCustomMessage(type='invertSelection', message = list())
    session$sendCustomMessage(type = "hideSelection", message = list())
  })
  
  
  
  observeEvent(input$fitSelected, ignoreInit = TRUE, {
    fitSelected(session, 100)
  })
  
  observeEvent(input$getSelectedNodes, ignoreInit = TRUE, {
    output$selectedNodesDisplay <- renderText({
      " "
    })
    getSelectedNodes(session)
  })
  
  observeEvent(input$clearSelection, ignoreInit = TRUE, {
    session$sendCustomMessage(type = "clearSelection", message = list())
  })
  
  
  observeEvent(input$selectedNodes, {
    newNodes <- input$selectedNodes
    
    output$selectedNodesDisplay <- renderText({
      paste(newNodes)
    })
    
    input_gene_names <- sapply(toupper(input$selectedNodes), function(x) gene_names[grep(paste0("^", x, "$"), capitalized_gene_names)])
    input_gene_names <- as.character(input_gene_names)
    
    edit <- isolate(input$select_genes)
    if (length(edit) < 12) {
      selection <- c(edit, input_gene_names)
    } else {
      selection <- edit
    }
    
    updateSelectizeInput(
      session,
      "select_genes",
      choices = selection,
      selected = selection,
      server = TRUE
    )
  })
  

  

  
  # renders the geneCard link for the selected gene
  output$geneLink <- renderUI({
    url <-
      a(
        sprintf("Genecard %s link", input$selectedNodes[1]),
        href = sprintf(
          "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s",
          input$selectedNodes[1]
        ),
        target = "_blank"
      )
    tagList(sprintf("%s: ", input$selectedNodes[1]), url)
  })
  


  
  observeEvent(input$selectCondition, {
    k <- isolate(input$Key_pick)
    d <- isolate(input$DiffGRN_pick)
    c <- isolate(input$selectCondition)
    print(c)
    n.edges <- isolate(input$n.edges)
    edges<-edges[(tool %in% d) & (key %in% k) & (interaction %in% c)]
    if (nrow(edges)<=0){
      return(NULL)
    }
    color.map(make_color_map(edges, selections))
    print(color.map())
    print('Legend')
    output$legend<-renderPlot({
      ggplot(color.map()[tool %in% d], aes(x = interaction, y = tool, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_identity() +
        labs(x = "Condition", y = "Tool") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        guides(fill = FALSE)  # Hide the legend + facet_wrap(~key)
    })
  })

  # function for plotting a linear model of two conditions given the x,y expression data of both
  # plot_linear_model <- function(expr1_x, expr1_y, expr2_x, expr2_y, conditions) {
  plot_linear_model <- function(df, source, target) {
    cm<-merge(color.map(), selections, by.x = c('key', 'interaction'), by.y = c('key', 'condition'))
    cm<-cm[key %in%  isolate(input$Key_pick)]
    cm<-cm[!duplicated(cm[, .(key, variables)])]
    
    cols<-cm$value
    names(cols)<-cm$variables
    
    title <-'Linear model'
    plot <- ggplot(df, aes(x = x, y = y, colour = Condition)) +
      geom_point() +
      scale_colour_manual(values = cols) +
      geom_smooth(
        method = "lm",
        formula = "y ~ x",
        se = FALSE
      ) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
    
    plot <- plot +
      theme_classic() +
      labs(x = paste0(source, " expression"), y = paste0(target, " expression")) +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 15),
        legend.position = 'bottom'
      ) + expand_limits(x = 0, y = 0)
    return(list(plot = plot, title = title))
  }
  
  get_expression_data<-function(source, target, conditions, adata){
    elist<-list()
    
    conds<-isolate(input$selectCondition)
    conditions<-unique(selections[condition %in% conds]$variables)
    
    print(conditions)
    for (co in conditions){
      if (co %in% unique(adata@meta.data[, group.var])){
      e<-FetchData(object = adata, vars = group.var)
      expr1 <-as.data.table(t(as.matrix(adata[, which(e==co)][c(source, target)]@assays$RNA@layers$counts)))
      expr1$condition <- co
      elist[[co]]<-expr1
    }
    }
    
    expression <- rbindlist(elist)
    ag <-aggregate(. ~ condition, expression, function(x) {c(mean = mean(x), sd = sd(x))})
    stats <- ag %>% pivot_longer(-c("condition"))
    stats <- as.data.table(stats)[, 1:4]
    colnames(stats) <- c("condition", "gene", "mean", "std")
    
    expression$row <- 1:nrow(expression)
    expression <-
      expression %>%
      pivot_longer(-c("condition", "row")) %>%
      as.data.table()
    expression <- merge(expression,stats, by.x = c("condition", "name"), by.y = c("condition", "gene"))
    expression[, outlier := ifelse(value > (mean + 2 * std), "outlier", "normal")]

    key <- isolate(input$Key_pick)
    
    print('Generating colors')
    cm<-merge(color.map(), selections, by.x = c('key', 'interaction'), by.y = c('key', 'condition'))
    cm<-cm[key %in%  isolate(input$Key_pick)]
    cm<-cm[!duplicated(cm[, .(key, variables)])]
    
    cols<-cm$value
    names(cols)<-cm$variables
    
    expression$color<-sapply(expression$condition, function(x) cols[x])
    expression <- expression[outlier != "outlier"]
    expression <-
      expression %>%
      select(condition, name, value, color, row) %>%
      group_by(condition, color) %>%
      pivot_wider(names_from = name, values_from = c(value)) %>%
      select(-row) %>%
      as.data.table()
    colnames(expression) <- c("Condition", "cols", "x", "y")
    return(expression)
  }
  
  
  output$linear_model_plot <- renderPlot({
    if (is.null(input$selectedNodes[1]) || is.null(input$selectedNodes[2])) {
      return(NULL)
    }
    if (is.na(input$selectedNodes[1]) || is.na(input$selectedNodes[2])) {
      return(NULL)
    }
    
    s <- input$selectedNodes[1]
    t <- input$selectedNodes[2]
    k<-isolate(input$Key_pick)
    tl<-isolate(input$DiffGRN_pick)
    
    if(nrow(edges[source==s & target==t & key %in% k & tool %in% tl]) >= 1){
      source<-s
      target<-t
    }
    else{
      source<-t
      target<-s
    }
    
    conditions <- isolate(input$selectCondition)
    print(conditions)
    if (!metacells) {
      expression<-get_expression_data(source, target, conditions, adata)
    }
    else{


      expression<-get_expression_data(source, target, conditions, metacell.seurat)
      print(expression)
    }
    # expr1_x <- as.numeric(displayed_network$cond1_metacells_data[displayed_network$cond1_metacells_data$Gene == source, -1])
    # expr1_y <- as.numeric(displayed_network$cond1_metacells_data[displayed_network$cond1_metacells_data$Gene == target, -1])
    # expr2_x <- as.numeric(displayed_network$cond2_metacells_data[displayed_network$cond2_metacells_data$Gene == source, -1])
    # expr2_y <- as.numeric(displayed_network$cond2_metacells_data[displayed_network$cond2_metacells_data$Gene == target, -1])
    # res <- plot_linear_model(expr1_x, expr1_y, expr2_x, expr2_y, conditions)

    res <- plot_linear_model(expression, source, target)
    
    output$downloadGRNLinearModelPlot <-
      download_plot(res[[1]], res[[2]])
    
    res[[1]]
  })
  
  get_comp_selections <- function(comp_selec1, comp_selec2) {
    col_idx <- 1
    col_names <-
      unique(lapply(strsplit(comp_selec1, ":"), function(x) {
        x[1]
      }))
    comp1_select.cells <- c()
    cells <- c()
    for (selection in comp_selec1) {
      s <- strsplit(selection, ":")[[1]]
      col <- s[1]
      if (col != col_names[col_idx]) {
        if (length(comp1_select.cells) == 0) {
          comp1_select.cells <- cells
        } else {
          comp1_select.cells <- intersect(comp1_select.cells, cells)
        }
        cells <- c()
        col_idx <- col_idx + 1
      }
      
      val <- str_sub(s[2], 2) # removing the whitespace for good looks :)
      cells <- c(cells, which(adata@meta.data[, col] == val))
    }
    if (length(col_names) > 1) {
      comp1_select.cells <- intersect(comp1_select.cells, cells)
    } else {
      comp1_select.cells <- cells
    }
    
    col_idx <- 1
    col_names <-
      unique(lapply(strsplit(comp_selec2, ":"), function(x) {
        x[1]
      }))
    comp2_select.cells <- c()
    cells <- c()
    for (selection in comp_selec2) {
      s <- strsplit(selection, ":")[[1]]
      col <- s[1]
      if (col != col_names[col_idx]) {
        if (length(comp2_select.cells) == 0) {
          comp2_select.cells <- cells
        } else {
          comp2_select.cells <- intersect(comp2_select.cells, cells)
        }
        cells <- c()
        col_idx <- col_idx + 1
      }
      val <-
        str_sub(s[2], 2) # removing the whitespace that is there for good looks :)
      cells <- c(cells, which(adata@meta.data[, col] == val))
    }
    if (length(col_names) > 1) {
      comp2_select.cells <- intersect(comp2_select.cells, cells)
    } else {
      comp2_select.cells <- cells
    }
    
    return(
      list(
        comp1_select.cells = comp1_select.cells,
        comp2_select.cells = comp2_select.cells
      )
    )
  }
  

  
  
  
  # Violin and linear model plot for the shown DiffGRN
  output$violin_plot <- renderPlot({
    print('Violin plot')
    if (is.null(input$selectedNodes[1])) {
      print('No nodes selected')
      return(NULL)
    }
    print('null')
    key <- isolate(input$Key_pick)

    print('Generating colors')
    cm<-merge(color.map(), selections, by.x = c('key', 'interaction'), by.y = c('key', 'condition'))
    print(cm)
    cm<-cm[key %in%  isolate(input$Key_pick)]
    cm<-cm[!duplicated(cm[, .(key, variables)])]

    cols<-cm$value
    names(cols)<-cm$variables
    
    print('probelm')
    gn<-isolate(input$selectedNodes)[1]
    
    input_gene_name <- gene_names[grep(paste0("^", toupper(gn), "$"), capitalized_gene_names)]
    print(input_gene_name)
    plot <-
      VlnPlot(adata,
              features = input_gene_name,
              #group.by = group_var,
              idents = cm$variables,
              cols = cols,
              assay = opt$assay)+
      theme(legend.position = "bottom")

    output$downloadGRNViolinPlot <-
      download_plot(
        plot,
        sprintf(
          "%s || %s vs. %s",
          isolate(input$selectedNodes)[1],
          condition_names[1],
          condition_names[2]
        )
      )
    plot
  })
  
  
  # function for downloading a plot
  download_plot <- function(plot, filename) {
    downloadHandler(
      filename = function() {
        sprintf("%s.png", filename) # Filename for the downloaded plot
      },
      content = function(file) {
        # Save the plot as a file
        ggsave(file, plot = plot)
      }
    )
  }
  
  
  
  
  ########################## SECOND TAB (Seurat object information) ################################
  umap_plot <- eventReactive(input$plot_button, {
    umap_plot <-
      try(FeaturePlot(adata, features = input$select_genes, ncol = 6),
          silent = TRUE
      )
    if (class(umap_plot) == "try-error") {
      showNotification(
        "No dimensionality plots were found in the seurat object. If you want to show any, please add them to the seurat object.",
        type = "warning"
      )
    } else {
      umap_plot
    }
  })
  output$umap_plot <- renderPlot({
    umap_plot()
  })
  
  standard_violin_plot <- eventReactive(input$plot_button, {
    # Seurat's VlnPlot does not show the legend if more than one gene is plotted (see https://github.com/satijalab/seurat/issues/2598) -> using cowplot's plot_grid with list of violin plots
    mygenes<-intersect(input$select_genes, rownames(adata))
    plots <-
      VlnPlot(
        adata,
        features = mygenes,
        group.by = input$select_primary_grouping,
        split.by = input$select_secondary_grouping,
        ncol = 6,
        combine = FALSE
      )
    do.call(plot_grid, plots)
    # VlnPlot(adata, features = input$select_genes, group.by=input$select_primary_grouping,  split.by = input$select_secondary_grouping, ncol=6)
  })
  output$standard_violin_plot <-
    renderPlot({
      standard_violin_plot()
    })
  
  dot_plot <- eventReactive(input$plot_button, {
    # Seurat's VlnPlot does not show the legend if more than one gene is plotted (see https://github.com/satijalab/seurat/issues/2598) -> using cowplot's plot_grid with list of violin plots
    mygenes<-intersect(input$select_genes, rownames(adata))
    DotPlot(
      adata,
      features = mygenes,
      group.by = input$select_primary_grouping
    )
  })
  output$dot_plot <- renderPlot({
    dot_plot()
  })
  
  # updateSelectizeInput(session, "select_genes", choices = gene_names, server = TRUE)
  
  
  ##########################
}

# RUN SHINY APP ----
shinyApp(ui = ui, server = server)
