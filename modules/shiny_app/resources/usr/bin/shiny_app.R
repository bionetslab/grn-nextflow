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
################################################################################
print("Hi")

#source('style_generator.R')
# Function to generate the style JSON based on tools and colors
generate_style_json <- function(tools, colors, conditions, keys, arrow_styles, dash_styles) {
  edge_styles <- c()
  
  for (i in seq_along(tools)) {
    tool <- tools[i]
    color <- colors[i]
    condition <-conditions[i]
    key <- keys[i]
    arrow_style<-arrow_styles[i]
    
    if (length(unique(keys))>3){
      print(unique(keys))
      print('many arrow stke')
      arrow_style<-'solid'
      
      edge_style <- paste0('{
      "selector": "edge[tool = \'', tool, '\'][interaction = \'', condition, '\'][key = \'', key, '\'][effect < 0]",
      "css": {
        "line-color": "', color, '",
        "target-arrow-color": "', color, '",
        "target-arrow-shape": "tee",
        "line-style": "', arrow_style, '"
        }
    }')
      edge_styles <- c(edge_styles, edge_style)
      
      edge_style <- paste0('{
      "selector": "edge[tool = \'', tool, '\'][interaction = \'', condition, '\'][key = \'', key, '\'][effect >= 0]",
      "css": {
        "line-color": "', color, '",
        "target-arrow-color": "', color, '",
        "target-arrow-shape": "triangle",
        "line-style": "', arrow_style, '"

      }
    }')
      
      edge_styles <- c(edge_styles, edge_style)
    }
    else{
      edge_style <- paste0('{
      "selector": "edge[tool = \'', tool, '\'][interaction = \'', condition, '\'][key = \'', key, '\'][effect < 0]",
      "css": {
        "line-color": "', color, '",
        "target-arrow-color": "', color, '",
        "target-arrow-shape": "tee",
        "line-style": "', arrow_style, '"
        }
    }')
      edge_styles <- c(edge_styles, edge_style)
      
      edge_style <- paste0('{
      "selector": "edge[tool = \'', tool, '\'][interaction = \'', condition, '\'][key = \'', key, '\'][effect >= 0]",
      "css": {
        "line-color": "', color, '",
        "target-arrow-color": "', color, '",
        "target-arrow-shape": "triangle",
        "line-style": "', arrow_style, '"

      }
    }')
      e<-'{
        "selector": "edge[key]",
        "css": {
          "label": "data(key)",
          "width": 3
        }
      }'
      edge_styles <- c(edge_styles, edge_style)
    }
  }
  if (length(unique(keys))>3){
    
    e<-'{
        "selector": "edge[key]",
        "css": {
          "label": "data(key)",
          "width": 3
        }
      }'
    print(e)
    edge_styles<-c(edge_styles, e)
    
  }
  
  style_json <- paste0('[', 
                       '{
                         "selector": "node",
                         "css": {
                           "border-width": "2px",
                           "width": "data(size)",
                           "height": "data(size)",
                           "content": "data(id)"
                         }
                       },
                       {
                         "selector": "edge",
                         "css": {
                           "curve-style": "bezier",
                           "control-point-step-size": 40,
                           "target-arrow-shape": "triangle"
                         }
                       },',
                       paste(edge_styles, collapse = ","),
                       ',{
                         "selector": "node:selected",
                         "css": {
                           "overlay-opacity": 0.3,
                           "overlay-color": "gray",
                           "background-color": "blue"
                         }
                       },
                       {
                         "selector": "node[value <= 0]",
                         "css": {
                           "background-color": "mapData(value, -5, 0, green, white)"
                         }
                       },
                       {
                         "selector": "node[value > 0]",
                         "css": {
                           "background-color": "mapData(value, 0, 10, white, red)"
                         }
                       }
                     ]')
  
  return(style_json)
}


option_list <- list(
  make_option(
    c("-r", "--results.path"),
    type = "character",
    default = "",
    help = "Path to results"
  ),
  make_option(
    c("-p", "--project.path"),
    type = "character",
    default = "",
    help = "Project path"
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
  )
)
# changing opt$grntoools and opt$dgrntools to be suitable as input for the Shiny app
opt <- parse_args(OptionParser(option_list = option_list))




# loading the modified network3d lib
# library(networkD3, lib.loc=sprintf("%s/lib", opt$project.path))
metacells <- FALSE
#opt$results.path <-
#  "/home/bionets-og86asub/Documents/external_analyses/huiqin/"
#opt$project.path <-
#  "/home/bionets-og86asub/Documents/netmap/workflows"
#opt$selection <-
#  "M0,RNA,cluster_type,NSCLC_M0,NSCLC_M0-M0,RNA,cluster_type,BRCA_M0,BRCA_M0-M1,RNA,cluster_type,NSCLC_M0,NSCLC_M0-M1,RNA,cluster_type,PRAD_M0,PRAD_M0-M2,RNA,cluster_type,NSCLC_M0,NSCLC_M0-M2,RNA,cluster_type,BRCA_M0,BRCA_M0-M3,RNA,cluster_type,NSCLC_M0,NSCLC_M0-M3,RNA,cluster_type,PRAD_M0,PRAD_M0"
#opt$seurat.file <-
#  "/home/bionets-og86asub/Documents/external_analyses/huiqin/cancer.rds"
#opt$dgrntools <- "boostdiff,zscores,diffcoex"
#opt$dgrntools <- c(strsplit(opt$dgrntools, ",")[[1]])
#opt$grntools <- "grnboost2"
if (opt$grntools != "No tools were chosen for GRN Inference") {
  opt$grntools <- c("NA", opt$grntools)
}
# opt$m <- "seurat"

create_graph <- function(all_networks, de.genes = NULL) {
  nodes <-
    data.table(
      id = c(all_networks$target, all_networks$regulator),
      tool = c(all_networks$tool, all_networks$tool),
      keys = c(all_networks$key, all_networks$key)
    )
  # Necessary because key is reserved keywort in df.
  colnames(nodes)<-c('id', 'tool', 'key')
  nodes <- unique(nodes)
  if (!is.null(de.genes)) {
    nodes <-
      merge(nodes,
            de.genes,
            by.x = "id",
            by.y = "gene",
            all.x = TRUE
      )
  }
  edges <- all_networks
  colnames(edges) <-
    c(
      "target",
      "source",
      "weight",
      "interaction",
      "effect",
      "tool",
      "key"
    )
  return(list(nodes, edges))
}

read_files <- function(file_list, key) {
  all_data <- list()
  
  for (file in file_list) {
    tryCatch(
      {
        data <- fread(file)
        if (!is.null(data) && nrow(data) > 0) {
          print(basename(dirname(file)))
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
  print(all_data)
  all_data <- rbindlist(all_data)
  return(all_data)
}

make_color_map<-function(all_networks){
  n.colors<-nrow(unique(all_networks[, .(interaction)]))
  color.pal<-rcartocolor::carto_pal(n=max(3, 2*n.colors), 'Prism')
  color.pal<-color.pal[seq(1, 2*n.colors, 2)]
 
  color.df<-data.table(color.pal[1:n.colors], unique(all_networks[, .(interaction)]))
  colnames(color.df)<-c('color', 'interaction')
  color.df<-color.df %>% pivot_longer(-c('interaction'))
  
  combos<-unique(all_networks[, .(interaction, tool)])
  combos<-unique(cross_join(combos, combos)[,.(interaction.x,tool.y)])
  colnames(combos)<-c('interaction', 'tool')
  combos<-merge(combos, color.df, by = 'interaction', allow.cartesian=T)
  combos<-combos %>% group_by(interaction) %>% mutate(row_count = row_number())
  combos$value<-sapply(1:nrow(combos), function(i) lighten(combos$value[i], (combos$row_count[i]-1)*(1/max(combos$row_count)))) 
  combos<-as.data.table(combos)
  combos<-combos[, .(tool, interaction, value)]
  
  cl<-list()
  possible_styles<-c('solid', 'dotted', 'dashed')
  dashes<- ceiling(sqrt(length(unique(all_networks$key))-1))
  dash.styles<-c('NA', 'NA',as.character(sapply(1:dashes, function(x) sapply(1:3, function(y) paste0('[', 2*x, ', ', 2*y, ']')))))
  dash.width<- c('NA', 'NA', sapply(1:dashes, function(x) sapply(1:3, function(y) 2*x)))
  dash.gap<- c('NA', 'NA', sapply(1:dashes, function(x) sapply(1:3, function(y) 2*y)))
  possible_styles<-c('solid', 'dotted', rep('dashed', dashes*dashes))

  for (i in 1:length(unique(all_networks$key))){
    c<-combos
    c$arrow_style<-possible_styles[i]
    c$key<-unique(all_networks$key)[i]
    c$dash.style<-dash.styles[i]
    c$dash.width<-dash.width[i]
    c$dash.gap<-dash.gap[i]
    
    cl[[i]]<-c
  }
  combos<-rbindlist(cl)
  return(combos)
}

get_color<-function(color.df, t){
  print('Getting color')
  print(color.df[(tool==t) ])
  values<- unique(color.df[(tool==t)]$value)
  print(values)
  names(values)<-unique(color.df$interaction)
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
  print(color.map)

  # Generate the style JSON
  style_json <- generate_style_json(color.map$tool, color.map$value, color.map$condition, color.map$key, color.map$arrow_style, color.map$dash.style)
  # Write the style JSON to a file
  writeLines(style_json, "style.js")
  
  
  ########### Loading Seurat object, filtering the correct cells and performing differential testing
  adata <- readRDS(opt$seurat.file)
  # saves a bit of computation time
  gene_names <- rownames(adata)
  capitilized_gene_names <- toupper(gene_names)
  # adata <- adata$all
  
  group.var <- selections$group_var[1]
  Idents(adata) <- group.var

  # get all factors and columns that were used in the selection but were not of class factor
  factors <- unique(c(names(which(sapply(adata@meta.data, class)=='factor')), selections$group_var[1]))
  # get all unique values for the factors
  factor_vals <- sapply(factors, function(factor) { unique(unlist(adata@meta.data[, factor], recursive = FALSE)) })
  # remove the levels from the entries that are factors
  factor_vals <- sapply(factor_vals, function(x) { if(class(x) == 'factor') { droplevels(x) } else { x }})
  # add name of factor as key to the values (needed for identification inside of shiny)
  factor_vals <- sapply(names(factor_vals), function(name) { paste0(rep(name, length(factor_vals[[name]])), ': ', factor_vals[[name]])})
  
  
  # Idents(adata) <- "cluster_type"
  # de.genes <- FindAllMarkers(adata)
  # de.genes$gene <- rownames(de.genes)
  # de.genes$name <- de.genes$cluster
  # colnames(de.genes)[2] <- "value"
  de.genes<-NULL
  # 
  graph <- create_graph(all_networks, de.genes)
  nodes <- graph[[1]]
  edges <- graph[[2]]
}

basicStyleFile <- "style.js"

# SET INPUT OPTIONS ----
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
                            choices = unique(all_networks$interaction),
                            selected = unique(all_networks$interaction)[1],
                            multiple = TRUE,
                            selectize = TRUE,
                            width = NULL,
                            size = NULL)
              )
            )
          )),
          sliderInput('n.edges', 'Number of edges', 1, nrow(edges),  nrow(edges), step = 1, round = FALSE),
          
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
            choices = unique(nodes$name)
          ),
          actionButton('generate_network', 'Generate Network'),
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
          ),
          tabPanel(
            "Dotplot",
            fluidRow(column(
              12, plotOutput("dot_plot")
            ))
          )
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
  

  render_graph<-function(genes = NULL){
    
    print('Rendering graph')
    k <- isolate(input$Key_pick)
    d <- isolate(input$DiffGRN_pick)
    c <- isolate(input$selectCondition)
    n.edges <- isolate(input$n.edges)
    
    print(c)
    edges <- edges %>% group_by(tool, key) %>% slice_max(n = n.edges, weight) %>% as.data.table()
    edges<-edges[(tool %in% d) & (key %in% k) & (interaction %in% c)]
    color.map(make_color_map(edges))
    print(color.map)
    # Generate the style JSON
    style_json <- generate_style_json(color.map()$tool, color.map()$value, color.map()$interaction, color.map()$key, color.map()$arrow_style, color.map()$dash.style)
    # Write the style JSON to a file
    writeLines(style_json, "style.js")
    loadStyleFile(basicStyleFile)
    
    if (!is.null(genes)){
    edges<-edges[ (source %in% genes) & (target %in% genes)]
    }

    print(d)
    print(edges)
    nodes<-nodes[id %in% c(edges$target, edges$source)]
    print(nodes)
    nodes<-nodes[(tool %in% d) & (key %in% k)]
    
    print(k)
    print(edges)

    
    graph.json <-
      dataFramesToJSON(
        tbl.nodes = nodes,
        tbl.edges = edges
      )

    output$cyjShiny <- renderCyjShiny({
      cyjShiny(graph.json,
               layoutName = "cola",
               styleFile = basicStyleFile
      )
    })

    
    output$legend<-renderPlot({
      ggplot(color.map()[tool %in% d], aes(x = interaction, y = tool, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_identity() +
        labs(x = "Condition", y = "Tool") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        guides(fill = FALSE)  # Hide the legend + facet_wrap(~key)
    })

  }
  
  
  observeEvent(input$generate_network, {
    print('Trigger')
    render_graph()
  })
  
  

  
      
  observeEvent(input$redrawGraph, ignoreInit = TRUE, {
    genes<-isolate(input$select_genes)
    if (length(genes)>0){
    render_graph(genes)
    }
    else{
      render_graph()
    }
  })
  
  observeEvent(input$setNodeAttributes, {
    attribute <- "value"
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
    isolate(
      setNodeAttributes(
        session,
        attributeName = attribute,
        nodes = updated.nodes,
        values = updated.values
      )
    )
    }
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
    print("Get selecte node")
    output$selectedNodesDisplay <- renderText({
      " "
    })
    getSelectedNodes(session)
  })
  
  observeEvent(input$clearSelection, ignoreInit = TRUE, {
    session$sendCustomMessage(type = "clearSelection", message = list())
  })
  
  
  observeEvent(input$selectedNodes, {
    print("Selected Node")
    newNodes <- input$selectedNodes
    
    output$selectedNodesDisplay <- renderText({
      paste(newNodes)
    })
    
    input_gene_names <- sapply(toupper(input$selectedNodes), function(x) gene_names[grep(paste0("^", x, "$"), capitilized_gene_names)])
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
  
  # loads the data of the uploadable gene list
  observeEvent(input$comparison_grn_file, {
    comparison_network <-
      read.table(input$comparison_grn_file$datapath,
                 sep = ",",
                 header = TRUE
      )[, 2:3]
    comparison_nodes <-
      unique(c(comparison_network[, 1], comparison_network[, 2]))
    displayed_network$fn$x$links$existsInComparisonGRN <-
      do.call(paste0, displayed_network$links[, 1:2]) %in% do.call(paste0, comparison_network)
    displayed_network$fn$x$nodes$existsInComparisonGRN <-
      displayed_network$fn$x$nodes$name %in% comparison_nodes
  })
  
  # observer events for showing the overlap (node, edges) of the uploadable gene list and the displayed network
  observeEvent(input$comparison_grn_switch_edges, {
    if (is.null(input$comparison_grn_file)) {
      return(NULL)
    }
    displayed_network$fn$x$links$comparisonGRNSwitch <-
      input$comparison_grn_switch_edges
  })
  
  observeEvent(input$selectCondition, {
    k <- isolate(input$Key_pick)
    d <- isolate(input$DiffGRN_pick)
    c <- isolate(input$selectCondition)
    n.edges <- isolate(input$n.edges)
    edges<-edges[(tool %in% d) & (key %in% k) & (interaction %in% c)]
    color.map(make_color_map(edges))
    
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
  
  observeEvent(input$comparison_grn_switch_nodes, {
    if (is.null(input$comparison_grn_file)) {
      return(NULL)
    }
    displayed_network$fn$x$nodes$comparisonGRNSwitch <-
      input$comparison_grn_switch_nodes
  })
  
  # function for plotting a linear model of two conditions given the x,y expression data of both
  # plot_linear_model <- function(expr1_x, expr1_y, expr2_x, expr2_y, conditions) {
  plot_linear_model <- function(df, source, target) {
    print(df)
    mapped_cols <- unique(df$cols)
    names(mapped_cols) <- unique(df$Condition)
    print("create Plot")
    title <-'Linear model'
    plot <- ggplot(df, aes(x = x, y = y, colour = Condition)) +
      geom_point() +
      scale_colour_manual(values = mapped_cols) +
      geom_smooth(
        method = "lm",
        formula = "y ~ x",
        se = FALSE
      ) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
    
    print('add title')
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

    for (co in conditions){
      if (co %in% unique(adata$cluster_type)){
      print(co)
      expr1 <-
        as.data.table(t(as.matrix(
          subset(adata, subset = cluster_type == co)[c(source, target)]@assays$RNA@counts
        )))
      
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
    print(color.map())
    colors<-get_color(color.map(), t= isolate(input$DiffGRN_pick[1]))
    print(colors)
    expression$color<-sapply(expression$condition, function(x) colors[x])
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
    
    print("Linear Model plot")
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
    
    if (!metacells) {
      expression<-get_expression_data(source, target, conditions, adata)
    }
    else{
      print('FAAAIIIL')
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
  
  # Comparison violin and linear model plots if a filter was set
  
  meta_cell_creation <- function(subset) {
    n.samples <- 100
    p.missing <- 10
    n.cells <- nrow(subset@meta.data)
    # Compute number of cells to aggregate
    cells.p.metasample <- nrow(subset@meta.data) / opt$n.samples
    # randomly assign each of the cells to a group
    set.seed(1)
    subset@meta.data$meta.cell <-
      sample(
        nrow(subset@meta.data),
        size = nrow(subset@meta.data),
        replace = FALSE
      ) %% opt$n.samples
    # Set the ident to the newly created meta.cell variable
    Idents(subset) <- "meta.cell"
    # Aggregate the expression
    agg <-
      AggregateExpression(
        subset,
        slot = "counts",
        return.seurat = T,
        assays = configuration[configuration$key == input$Key_pick]$assay
      )
    agg@assays[[opt$assay]]@counts <-
      agg@assays[[opt$assay]]@counts / cells.p.metasample
    
    agg <- NormalizeData(agg)
    result.data.frame <- agg@assays[[opt$assay]]@data
    
    row.names <- rownames(result.data.frame)
    column.names <-
      paste0(
        paste0(input$Key_pick, collapse = "_"),
        "_",
        colnames(result.data.frame)
      )
    result.data.frame <- as.data.table(result.data.frame)
    result.data.frame <- cbind(row.names, result.data.frame)
    colnames(result.data.frame) <- c("Gene", column.names)
    select <-
      which(rowSums(result.data.frame == 0) / (ncol(result.data.frame) - 1) <
              (opt$p.missing / 100))
    result.data.frame <- result.data.frame[select, ]
    result.data.frame$Gene <- toupper(result.data.frame$Gene)
    return(result.data.frame)
  }
  
  comparison_linear_model_plot <-
    eventReactive(input$compare_button, {
      if (opt$mode == "tsv" ||
          is.null(input$id_source) || is.null(input$id_target)) {
        return(NULL)
      }
      comp_selec1 <- input$comp_selec_cond1
      comp_selec2 <- input$comp_selec_cond2
      cond_name1 <- input$comp_name_cond1
      cond_name2 <- input$comp_name_cond2
      conditition_names <- c(cond_name1, cond_name2)
      
      comp_select.cells <-
        get_comp_selections(comp_selec1, comp_selec2)
      comp1_select.cells <- comp_select.cells[[1]]
      comp2_select.cells <- comp_select.cells[[2]]
      
      comp1_subset <- subset(adata, cells = comp1_select.cells)
      comp2_subset <- subset(adata, cells = comp2_select.cells)
      
      comp1_meta_cell_df <- meta_cell_creation(comp1_subset)
      comp2_meta_cell_df <- meta_cell_creation(comp2_subset)
      
      df_gene_names <-
        intersect(comp1_meta_cell_df$Gene, comp2_meta_cell_df$Gene)
      comp1_meta_cell_df <-
        comp1_meta_cell_df[Gene %in% df_gene_names]
      comp2_meta_cell_df <-
        comp2_meta_cell_df[Gene %in% df_gene_names]
      
      comp1_meta_cell_df$Gene <- toupper(comp1_meta_cell_df$Gene)
      comp2_meta_cell_df$Gene <- toupper(comp2_meta_cell_df$Gene)
      
      expr1_x <-
        as.numeric(comp1_meta_cell_df[Gene == toupper(isolate(input$id_source)), -1])
      expr1_y <-
        as.numeric(comp1_meta_cell_df[Gene == toupper(isolate(input$id_target)), -1])
      expr2_x <-
        as.numeric(comp2_meta_cell_df[Gene == toupper(isolate(input$id_source)), -1])
      expr2_y <-
        as.numeric(comp2_meta_cell_df[Gene == toupper(isolate(input$id_target)), -1])
      
      # outlier detection
      res <- plot_linear_model(expr1_x, expr1_y, expr2_x, expr2_y)
      output$downloadGRNLinearModelPlot <-
        download_plot(res[[1]], res[[2]])
      
      res[[1]]
    })
  
  # Violin and linear model plot for the shown DiffGRN
  output$violin_plot <- renderPlot({
    if (is.null(input$selectedNodes[1]) || opt$mode == "tsv") {
      return(NULL)
    }
    key <- isolate(input$Key_pick)
    condition_names <-
      unique(selections[key %in% isolate(input$Key_pick)]$condition)[1:2]
    print(condition_names)
    
    cols <-get_color(color.map(), t= isolate(input$DiffGRN_pick[1]))
    #group_var <- selections[key %in% isolate(input$Key_pick)]$group_var
    #print(group_var)
    conds<-isolate(input$selectCondition)
    print(cols)
    input_gene_name <-
      gene_names[grep(paste0("^", toupper(input$selectedNodes[1]), "$"), capitilized_gene_names)]
    plot <-
      VlnPlot(adata,
              features = input_gene_name,
              #group.by = group_var,
              idents = conds,
              cols = cols
      ) + theme(legend.position = "bottom")
    output$downloadGRNViolinPlot <-
      download_plot(
        plot,
        sprintf(
          "%s || %s vs. %s",
          input$selectedNodes[1],
          condition_names[1],
          condition_names[2]
        )
      )
    plot
  })
  
  output$comparison_linear_model_plot <-
    renderPlot({
      comparison_linear_model_plot()
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
    print("VLN")
    # Seurat's VlnPlot does not show the legend if more than one gene is plotted (see https://github.com/satijalab/seurat/issues/2598) -> using cowplot's plot_grid with list of violin plots
    plots <-
      VlnPlot(
        adata,
        features = input$select_genes,
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
    DotPlot(
      adata,
      features = input$select_genes,
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
