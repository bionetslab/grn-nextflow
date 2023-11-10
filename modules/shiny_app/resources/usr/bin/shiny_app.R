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
#(htmlwidgets)
library(optparse)
library(shinyWidgets)
library(DT)
library(data.table)
library(shinythemes)
library(shinyhelper)
library(dplyr)
library(stringr)
library(cowplot)
################################################################################



option_list <- list(
  make_option(c("-r", "--results.path"), type='character',
              default="", help="Path to results"),
  make_option(c("-p", "--project.path"), type='character',
              default="", help="Project path"),
  make_option(c("-s", "--selection"), type='character',
              default="", help="Configurated selection"),
  make_option(c("-f", "--seurat.file"), type='character',
              default="", help="Seurat data file"),
  make_option(c("--dgrntools"), type='character',
              default="No tools were chosen for DGRN Inference", help="used diffGRN tools in the pipeline"),
  make_option(c("--grntools"), type='character',
              default="No tools were chosen for GRN Inference", help="used GRN tools in the pipeline"),
  make_option(c("-m", "--mode"), type='character',
              default="seurat", help="Specifies which file type was used as an input for the pipeline"),
  make_option(c("-n", "--n.samples"), type="integer",
              default=100, help="Number of meta cells used in the pipeline (Needed for comparison plot)")
)
# changing opt$grntoools and opt$dgrntools to be suitable as input for the Shiny app
opt <- parse_args(OptionParser(option_list=option_list))
if (opt$grntools != "No tools were chosen for GRN Inference") {
  opt$grntools <- c("NA", opt$grntools)
}
opt$dgrntools <- c(strsplit(opt$dgrntools, ",")[[1]])

# loading the modified network3d lib
library(networkD3, lib.loc=sprintf("%s/lib", opt$project.path))

# colors for the conditions
colors <- c("#cc79a7", "#009e73", "#0072b2")

# parsing the selection into a configuration data.table that has all necessary information to enable all features of the shiny app.
# every row corresponds to ONE comparison in the input configuration file (or the one comparison of tsv mode)
configuration <- data.table(key = NULL, assay = NULL, group_var = NULL, condition_names = NULL, tools = NULL, network_files = NULL)

if (opt$mode == 'tsv') {
  # get tools
  tools <- list.dirs(path = file.path(opt$results.path, opt$selection), full.names = FALSE, recursive = FALSE)
  tool_names <- paste(tools, collapse = ",")

  # get network files
  all_files <- list.files(path = file.path(opt$results.path, opt$selection), full.names = TRUE, recursive = TRUE)
  network_files <- c(grep("aggregated_filtered_network", all_files, value = TRUE))
  network_files <- paste(network_files, collapse = ",")

  # get condition names
  files <- list.files(path = file.path(opt$results.path, opt$selection), full.names = FALSE, recursive = FALSE)
  print(files)
  cond_files <- c(grep('out_', files, value = TRUE))
  cond_name1 <- gsub('out_', '', cond_files[1])
  cond_name1 <- gsub('.tsv', '', cond_name1)
  cond_name2 <- gsub('out_', '', cond_files[2])
  cond_name2 <- gsub('.tsv', '', cond_name2)
  cond_names <- paste(c(cond_name1, cond_name2), collapse = ",")

  configuration <- rbind(configuration,
    list(key = opt$selection, assay = 'RNA', group_var = NA, condition_names = cond_names, tools = tool_names, network_files = network_files))
  print(configuration)
} else {

  # parse selection
  selections <- strsplit(opt$selection, "-")[[1]]
  # get condition names and all used group variables (they can change between the two conditions) first 
  all_group_vars <- c()
  conditions <- list()
  for (s in selections) {
    split_string <- strsplit(s, ",")[[1]]
    key <- split_string[1]
    condition_name <- split_string[4]
    conditions[[key]] <- c(conditions[[key]], condition_name)
    all_group_vars <-c(all_group_vars, strsplit(split_string[3], ":")[[1]])
  }

  # parsing the selection
  for (s in selections) {
    print(s)
    split_string <- strsplit(s, ",")[[1]]
    key <- split_string[1]
    assay <- split_string[2]
    group_var <- split_string[3]
    split_group_var <- strsplit(group_var, ":")[[1]]
    # indicates which group_vars are used in each selection of the comparison
    group_var_idx <- paste(match(split_group_var, all_group_vars), collapse=":")
    group_vals <- split_string[5:(length(split_string))]
    
    # parsing all information. "," separates the different values (key, assay, group_var, ...). 
    group_vals <- paste(group_vals, collapse=",")
    all_files <- list.files(path = file.path(opt$results.path, key), full.names = TRUE, recursive = TRUE)
    network_files <- c(grep("aggregated_filtered_network", all_files, value = TRUE))
    network_files <- paste(network_files, collapse = ",")

    tools <- list.dirs(path = file.path(opt$results.path, key), full.names = FALSE, recursive = FALSE)
    tools <- paste(tools, collapse = ",")

    cond <- sort(paste(conditions[[key]], collapse = ","), decreasing = FALSE)

    config <- list(key = key, assay = assay, group_var = group_var, group_var_idx = group_var_idx, group_vals = group_vals, tools = tools, condition_names = cond, network_files = network_files)
    configuration <- rbind(configuration, config)
  }

  # Postprocessing configuration data.table so that only 1 row exists per given comparison in the configuration file
  # For one column of a row: "-" separates between values of selections. ":" separates list values of one value type (e.g. multiple group_vars are separated by :). 
  configuration <- unique(configuration)
  duplicates <- configuration[duplicated(configuration$key),]
  # Must be because of differing group_vars ->
  if (nrow(duplicates) > 0) {
    for (i in 1:nrow(duplicates)) {
      dup_row <- configuration[i*2,]
      orig_row <- configuration[(i*2) - 1, ]

      configuration[(i*2) - 1, ]$group_vals <- paste(orig_row$group_vals, dup_row$group_vals, sep = "-")
      # used to map the used group vars of the selection to the list of all used group vars.
      configuration[(i*2) - 1, ]$group_var_idx <- paste(orig_row$group_var_idx, dup_row$group_var_idx, sep = "-")
      configuration[(i*2) - 1, ]$group_var <- paste(orig_row$group_var, dup_row$group_var, sep=":")
    }
    # removing duplicate rows has to be done after merging the group_vars
    configuration <- configuration[-sapply(1:nrow(duplicates), function(x) {x * 2}),]
  }
  ########### Loading Seurat object, filtering the correct cells and performing differential testing
  adata <- readRDS(opt$seurat.file)
  # saves a bit of computation time
  gene_names <- rownames(adata)
  capitilized_gene_names <- toupper(gene_names)
  # adata <- adata$all
  group.var<-unique(strsplit(configuration[1]$group_var, ':')[[1]])
  Idents(adata)<-group.var
  
  # get all factors and columns that were used in the selection but were not of class factor
  factors <- unique(c(names(which(sapply(adata@meta.data, class)=='factor')), unlist(strsplit(configuration$group_var, ':')[[1]])))
  # get all unique values for the factors
  factor_vals <- sapply(factors, function(factor) { unique(unlist(adata@meta.data[, factor], recursive = FALSE)) })
  # remove the levels from the entries that are factors
  factor_vals <- sapply(factor_vals, function(x) { if(class(x) == 'factor') { droplevels(x) } else { x }})
  # add name of factor as key to the values (needed for identification inside of shiny)
  factor_vals <- sapply(names(factor_vals), function(name) { paste0(rep(name, length(factor_vals[[name]])), ': ', factor_vals[[name]])})
}

# Creating the UI of the Seurat object
ui <- 
  navbarPage(
    theme = shinytheme("cerulean"), title = "InterNet Xplorer", 
    # Moving the notifications into the top left corner for better visibility
    tags$head(
      tags$style(
        HTML(".shiny-notification {
             position:fixed;
             top: calc(25%);
             left: calc(25%);
             }
             "
            )
        )
    ),
    tabPanel("Differential Network Analysis",
      sidebarLayout(
        mainPanel(
          fluidRow(class="network",
            # Legend of the displayed network                
            htmltools::div(
              style = "padding: 10px; background-color: white;",
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 71px; margin-right: 5px; margin-bottom: 6px", colors[1])),
              textOutput("condition1_activator", inline = TRUE),
              htmltools::br(),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
              textOutput("condition1_repressor", inline = TRUE),
              htmltools::br(),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 71px; margin-right: 5px; margin-bottom: 6px", colors[2])),
              textOutput("condition2_activator", inline = TRUE),
              htmltools::br(),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
              textOutput("condition2_repressor", inline = TRUE),
              htmltools::br(),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 71px; margin-right: 5px; margin-bottom: 6px", colors[3])),
              sprintf("Base GRN Activator"),
              htmltools::br(),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[3])),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[3])),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[3])),
              htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[3])),
              sprintf("Base GRN Repressor"),
            ),
            # Displayed DGRN
            forceNetworkOutput("net"),
            hr(style = "border-top: 0.5px solid #000000; opacity: 0.2;")
          ),
          # Option switches
          fluidRow(class="selection switches",
            column(3,
              selectInput(
                'Key_pick',
                "Choose key to display DiffGRN/GRN for:",
                choices = c(configuration$key),
                selected = configuration$key[1],
                multiple = FALSE,
                selectize = TRUE,
                width = NULL,
                size = NULL
              ),
              selectInput(
                'DiffGRN_pick',
                "Choose DiffGRN to display:",
                choices = opt$dgrntools,
                selected = opt$dgrntools[1],
                multiple = FALSE,
                selectize = TRUE,
                width = NULL,
                size = NULL
              ),
              selectInput(
                'GRN_pick',
                'Choose GRN to display:',
                choices = c(opt$grntools),
                selected = opt$grntools[1],
                multiple = FALSE,
                selectize = TRUE,
                width = NULL,
                size = NULL
              ),
            ),
            column(3,
              helper(
                fileInput('comparison_grn_file', 'Upload edge list to show overlap',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    '.csv'
                  )
                ),
                type = "inline",
                content = c("You can upload a GRN (column 1: Source genes, column 2: Target genes) as a csv file for comparison to this DiffGRN.",
                            "By clicking on either of the checkboxes below the upload you can compare the nodes or the edges.",
                            "Opaque nodes/edges are NOT contained in the uploaded GRN.")
              ),
              checkboxInput(inputId = "comparison_grn_switch_nodes", label = "Compare nodes to uploaded file", value = FALSE),
              checkboxInput(inputId = "comparison_grn_switch_edges", label = "Compare edges to uploaded file", value = FALSE)
            ),
            column(3,
              textInput('comp_name_cond1', 'Name the first comparison condition', value = "Comparison condition 1", width = NULL, placeholder = NULL),
              uiOutput("comp_selec_cond1"),
              textInput('comp_name_cond2', 'Name the second comparison condition', value = "Comparison condition 2", width = NULL, placeholder = NULL),
              uiOutput("comp_selec_cond2"),
              uiOutput("compare_button"),
            ),
          ),
          width = 6
        ),
        # Violin and linear model plots
        sidebarPanel(
          fluidRow(
            column(6,
              plotOutput("violin_plot"),
              column(6,
                downloadButton("downloadGRNViolinPlot", "Download Plot"),
              ),
              column(6,
                uiOutput("geneLink")
              ),
            ),
            column(6, 
              plotOutput("linear_model_plot"),
              downloadButton("downloadGRNLinearModelPlot", "Download Plot")
              ),
          ),
          uiOutput("comparison_plots"),
          width = 6
        ), 
      ),
    ),
    ### Second Tab (View Seurat object information) ### 
    tabPanel(
      "Gene Expression", 
      h2("Gene Expression"),
        sidebarLayout(
          uiOutput('secondTab_sidepanel'),
          mainPanel(
            tabsetPanel(
              type = 'tabs',
              tabPanel(
                'Umap',
                fluidRow(12, 
                  column(12, plotOutput("umap_plot")), 
                  column(12, plotOutput('standard_violin_plot'))
                )
              ),
              tabPanel('Dotplot',
                fluidRow(12, 
                  column(12, plotOutput("dot_plot"))
                )
              )      
            )
          )
        )
    )
  )

# Backend of the shiny app
server <- function(input, output, session) {
  # Needed to show helper buttons in the shiny app
  observe_helpers(withMathJax = TRUE)

  output$secondTab_sidepanel <- renderUI({
    if (opt$mode == 'tsv') {
      return(NULL)
    }
    sidebarPanel(
      selectInput(
        'select_genes',
        'Select Genes to display',
        character(0),
        multiple = TRUE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      ),
      selectInput(
        'select_primary_grouping',
        'Select primary grouping variable',
        factors,
        selected = factors[1],
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      ),
      selectInput(
        'select_secondary_grouping',
        'Select secondary grouping variable',
        factors,
        selected = factors[1],
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      ),
      actionButton('plot_button', 'Generate Plots')
    )
  })

  # The forcenetwork needs to be a reactive value to change the node/link values
  output$mode <- renderText({
    paste0(opt$mode)
  })

  # UI of the comparison button/plots. This is only shown if there is something to compare to. For this reason, it has to be coded on the server side.
  output$compare_button <- renderUI({
    if(opt$mode == "tsv") {
      return(NULL)
    }
    helper(
      actionButton("compare_button", "Compare plots!"),
      type = "inline",
      content = c(
        "How to compare these results to a different case:",
        "1. Select a case for comparison. Availabe comparison cases are:",
        "* Armstrong vs. Docile, Spleen, day 10",
        "* Armstrong vs. Docile, Spleen, day 28",
        "* Armstrong vs. Docile, Liver, day 10",
        "* Armstrong vs. Docile, Liver, day 28",
        "2. Select the cell clusters of the dataset to use:",
        "* Standard: cluster 1, 2"
      )
    )
  })


  output$comparison_plots <- renderUI({
    if(opt$mode == "tsv"){
      return(NULL)
    }
    fluidRow(
      column(6, 
        plotOutput("comparison_violin_plot"),
        downloadButton("downloadComparisonViolinPlot", "Download Plot")
        ),      
      column(6, 
        plotOutput("comparison_linear_model_plot"),
        downloadButton("downloadComparisonLinearModelPlot", "Download Plot")
        ),
    )
  })

  # all needed reactive values of the displayed networks.
  # diffGRN_network_data: Network data of current DGRN network
  # GRN_network_data:     Network data of current GRN network (NULL if "NA" selected or if it was not calculated in the pipeline)
  # network_data:         Merged network data of diffGRN_network_data and GRN_network_data
  # fn:                   Rendered network
  # links:                Links of the rendered network (needed for comparison to uploadable gene list)
  # diffConditions:       Differential Conditions of the current selected key
  # grnConditions:        "Condition" of the current selected GRN. WIP: remove it because it is probably unneccesary
  # adata:                adata object. Stored separately to add "group", "comparison" column for the violin/linear model plots to get the correct cells.
  # cond1_metacells_data: Metacell data of differential condition 1 (preloaded when selecting a key to save computation time)
  # cond2_metacells_data: Metacell data of differential condition 2 (preloaded when selecting a key to save computation time)
  # metacells_maxVal:     Maximum value of the metacell data. Used to keep a consistent range of the axis for the violin plot. The linear model plots do NOT have the same range because of visibility.

  if (opt$mode == 'tsv') {
    displayed_network <- reactiveValues(
      diffGRN_network_data = NULL, GRN_network_data = NULL, network_data = NULL, 
      fn = NULL, links = NULL, diffConditions = "", grnConditions = "", adata = NULL, 
      cond1_metacells_data = NULL, cond2_metacells_data = NULL, metacells_maxVal = 0
    )
  } else {
    displayed_network <- reactiveValues(
      diffGRN_network_data = NULL, GRN_network_data = NULL, network_data = NULL, 
      fn = NULL, links = NULL, diffConditions = "", grnConditions = "", adata = adata, 
      cond1_metacells_data = NULL, cond2_metacells_data = NULL, metacells_maxVal = 0
    )
  }

  # creates the displayed network using network3D. For more information see https://christophergandrud.github.io/networkD3/
  create_network <- function() {
    links <- data.frame(source = displayed_network$network_data[, 2], target = displayed_network$network_data[, 1], width=displayed_network$network_data[, 3]*100, stringsAsFactors=F)
    nodes <- data.frame(name = unique(c(unique(links$source), unique(links$target))), group = 1)
    displayed_network$links <- links # this is required to show the accurate link overlap with the upload edge list functionality because the node names are transformed into id's to create the force network

    nodes$id = seq(0, nrow(nodes)-1)
    # replacing node names with numeric values in links for plotting function
    i <- 0
    for(node in nodes$name) {
      links$source[links$source == node] = i
      links$target[links$target == node] = i
      i <- i + 1
    }

    network <- list(links = links, nodes = nodes)
    if (input$GRN_pick == "NA") {
      conditions <- c(displayed_network$diffConditions) 
    } else {
      conditions <- c(displayed_network$diffConditions, displayed_network$grnConditions) 
    }
    for(i in 1:length(conditions)) {
      displayed_network$network_data$condition <- ifelse(displayed_network$network_data$condition == conditions[i], colors[i], displayed_network$network_data$condition)
    }

    fn <- forceNetwork(Links = network$links, Nodes = network$nodes, 
                        Source = 'source', Target = 'target', 
                        NodeID = 'name', Group = 'group', Value="width", opacity = 1, arrows=TRUE, 
                        linkColour=displayed_network$network_data$condition, opacityNoHover = 0.6, zoom=TRUE, colourScale = JS('d3.scaleOrdinal().domain(["1"]).range(["#000000"])'))

    # adding linear model information to force network
    fn$x$links$effect <- displayed_network$network_data[,5]

    # adding switch information as it is currently set
    fn$x$links$comparisonGRNSwitch <- input$comparison_grn_switch_edges
    fn$x$nodes$comparisonGRNSwitch <- input$comparison_grn_switch_nodes
    if (!is.null(input$comparison_grn_file)) {
      comparison_network <- read.table(input$comparison_grn_file$datapath, sep = ",", header = TRUE)[, 2:3]
      comparison_nodes <- unique(c(comparison_network[, 1], comparison_network[, 2]))
      fn$x$links$existsInComparisonGRN <- do.call(paste0, displayed_network$links[, 1:2]) %in% do.call(paste0, comparison_network)
      fn$x$nodes$existsInComparisonGRN <- fn$x$nodes$name %in% comparison_nodes
    }
    displayed_network$fn <- fn
  }
  
  ## Observer event for selecting the input key. Each row in the configuration data.frame (computed at the top) is identified by ONE UNIQUE key (if the config file does not use the same keys twice. If so, the behaviour is undefined)
  observeEvent(input$Key_pick, {
  
    # getting all network files of the correct row
    network.file <- strsplit(configuration[configuration$key == input$Key_pick,]$network_files, ",")[[1]][1]
    # trying to load the data of the first network file. If no edges were found using the selected approach (mainly happens with z_scores), nothing will be done.
    tmp <- try(read.table(file = network.file, header = TRUE))
    if (class(tmp) == "try-error") {
      showNotification(paste0("No edges were found for the condition ", input$DiffGRN_pick, "! Displaying network of previously selected condition."), type = "warning")
    } else {
      displayed_network$network_data <- tmp
      # This is a remnant of the idea of a GRN only mode. WIP to integrate this.
      if(input$DiffGRN_pick != "No tools were chosen for DGRN Inference") { # -> initial selected network will be a differential network
        displayed_network$diffGRN_network_data <- displayed_network$network_data
        displayed_network$diffConditions <- strsplit(configuration[configuration$key == input$Key_pick,]$condition_names, ",")[[1]]

        if (opt$mode != 'tsv') {
          # creating a column "group" in the adata object to match the conditions to the cells
          group_var.all <- strsplit(configuration[configuration$key == input$Key_pick,]$group_var, ":")[[1]]
          group_vals.all <- strsplit(configuration[configuration$key == input$Key_pick,]$group_vals, "-")[[1]]
          group_idx.all <- strsplit(configuration[configuration$key == input$Key_pick,]$group_var_idx, "-")[[1]]

          # calculate grouping column for the current selected key
          displayed_network$adata$ShinyGroup <- "NA"
          displayed_network$adata$ShinyComparison <- "NA"
          comp_selec1 <- list()
          comp_selec2 <- list()
          for (i in 1:length(group_idx.all)) {
            group_idx <- as.numeric(strsplit(group_idx.all[i], ":")[[1]])
            group_var <- group_var.all[c(group_idx)]
            group_val.selection <- strsplit(group_vals.all[i], ",")[[1]]
            filter_vector <- nrow(adata)
            idx <- 1
            for (col_name in group_var) {
              
              col_values <- displayed_network$adata@meta.data[, col_name]
              group_val <- unique(unlist(lapply(strsplit(group_val.selection, ':'), `[[`, idx)))
              if (i == 1) {
                comp_selec1[[col_name]] <- group_val
              } else {
                comp_selec2[[col_name]] <- group_val
              }
              group_val <- type.convert(group_val, as.is = TRUE)
              if (length(group_val) == 1) { 
                filter_vector <- filter_vector & (col_values == group_val)
              } else {
                filter_vector <- filter_vector & (col_values %in% group_val)
              }
              idx <- idx + 1
            }
            displayed_network$adata$ShinyGroup[filter_vector] <- gsub(pattern = '_', replacement = ' ', displayed_network$diffConditions[i])        
            # set the identifications for the comparison selections
            if (i == 1) {
              comp_selec1 <- unlist(sapply(names(comp_selec1), function(name) { paste0(rep(name, length(comp_selec1[[name]])), ': ', comp_selec1[[name]])}))
            } else {
              comp_selec2 <- unlist(sapply(names(comp_selec2), function(name) { paste0(rep(name, length(comp_selec2[[name]])), ': ', comp_selec2[[name]])}))
            }
            if (i == 1) {
              output$comp_selec_cond1 <- renderUI({
                virtualSelectInput(
                  'comp_selec_cond1',
                  label = "Choose different factors for the first condition for comparison:",
                  choices = factor_vals,
                  selected = comp_selec1,
                  multiple = TRUE,
                  width = NULL,
                  size = NULL
                )
              }) 
            } else {
              output$comp_selec_cond2 <- renderUI({
                virtualSelectInput(
                  'comp_selec_cond2',
                  label = "Choose different factors for the second condition for comparison:",
                  choices = factor_vals,
                  selected = comp_selec2,
                  multiple = TRUE,
                  width = NULL,
                  size = NULL
                )
              })
            }
          }
        }
      } else { # WIP: Fix grn only mode
        displayed_network$GRN_network_data <- displayed_network$network_data
      }

      # render text for the network legend according to the condition names of the selected network
      output$condition1_activator <- renderText({
        paste("Stronger in ", gsub(displayed_network$diffConditions[1], pattern="_", replacement=" "), " (Activator)")
      })
      output$condition1_repressor <- renderText({
        paste("Stronger in ", gsub(displayed_network$diffConditions[1], pattern="_", replacement=" "), " (Repressor)")
      })
      output$condition2_activator <- renderText({
        paste("Stronger in ", gsub(displayed_network$diffConditions[2], pattern="_", replacement=" "), " (Activator)")
      })
      output$condition2_repressor <- renderText({
        paste("Stronger in ", gsub(displayed_network$diffConditions[2], pattern="_", replacement=" "), " (Repressor)")
      })
      
      updateTextInput(session, "DiffGRN_pick", value=opt$dgrntools[1])
      updateTextInput(session, "GRN_pick", value=opt$grntools[1])

      # loading the metacell expression data for the two differential conditions
      # TODO: Debug this if no metacells were selected.
      gexp_path <- paste(opt$results.path, input$Key_pick, sep="/")
      all_files <- list.files(gexp_path)
      out_files <- grep("^out", all_files, value = TRUE)
      out_files1_condition_name <- gsub(pattern="out_", "", out_files[1])
      out_files1_condition_name <- gsub(pattern=".tsv", "", out_files1_condition_name)
      if (out_files1_condition_name[1] == displayed_network$diffConditions[1]) {
        displayed_network$cond1_metacells_data <- read.table(file = paste(gexp_path, out_files[1], sep="/"), header = TRUE, sep = "\t")
        displayed_network$cond2_metacells_data <- read.table(file = paste(gexp_path, out_files[2], sep="/"), header = TRUE, sep = "\t")
      } else {
        displayed_network$cond2_metacells_data <- read.table(file = paste(gexp_path, out_files[1], sep="/"), header = TRUE, sep = "\t")
        displayed_network$cond1_metacells_data <- read.table(file = paste(gexp_path, out_files[2], sep="/"), header = TRUE, sep = "\t")
      }

      # Computing the maximum value of the metacells
      displayed_network$metacells_maxVal <- ceiling(
        max(
          c(
            max(displayed_network$cond1_metacells_data[,2:ncol(displayed_network$cond1_metacells_data)]), 
            max(displayed_network$cond2_metacells_data[,2:ncol(displayed_network$cond2_metacells_data)])
          )
        )
      )
      create_network()
    }
  })
  
  # observer Event for the DGRN input button
  observeEvent(input$DiffGRN_pick, {
    # Nothing happens in GRN only mode (GRN only mode is WIP)
    if (input$DiffGRN_pick != "No tools were chosen for DGRN Inference") {
      # loading the corresponding network files
      all_network.files <- strsplit(configuration[configuration$key == input$Key_pick,]$network_files, ",")[[1]]
      # finding the correct network based on input$Key_pick and input$DiffGRN_pick
      network.file <- grep(paste(input$Key_pick, input$DiffGRN_pick, sep="/"), all_network.files, value = TRUE)
      # trying to load the data. Nothing happens if no edges were found using the selected approach
      tmp <- try(read.table(file = network.file, header = TRUE))
      if (class(tmp) == "try-error") {
        showNotification(paste0("No edges were found for the condition ", input$DiffGRN_pick, "! Displaying network of previously selected condition."), type = "warning")
      } else {
        # Changing the values of displayed network to the new network data.
        displayed_network$diffGRN_network_data <- tmp
        if (!is.null(displayed_network$GRN_network_data)) {
          displayed_network$network_data <- rbind(displayed_network$diffGRN_network_data, displayed_network$GRN_network_data)
          displayed_network$network_data[[1]] <-  str_to_title(displayed_network$network_data[[1]])
          displayed_network$network_data[[2]] <-  str_to_title(displayed_network$network_data[[2]])
        } else {
          displayed_network$network_data <- displayed_network$diffGRN_network_data
        }
        # creating the new network
        create_network()
      }
    }
  })
  
  # Observer event for the GRN input button
  observeEvent(input$GRN_pick, {
    # Only do something if a GRN was computed in the pipeline
    if (input$GRN_pick != "No tools were chosen for GRN Inference") {  
      if (input$GRN_pick != "NA") {
        all_network.files <- strsplit(configuration[configuration$key == input$Key_pick,]$network_files, ",")[[1]]
        # finding the correct network based on input$Key_pick and input$GRN_pick
        network.file <- grep(paste(input$Key_pick, input$GRN_pick, sep="/"), all_network.files, value = TRUE)
        # Changing the values of displayed network to the new network data.
        displayed_network$GRN_network_data <- read.table(file = network.file, header = TRUE)
        displayed_network$grnConditions <- displayed_network$GRN_network_data$condition
        if (!is.null(displayed_network$diffGRN_network_data)) {
          displayed_network$network_data <- rbind(displayed_network$diffGRN_network_data, displayed_network$GRN_network_data)
          displayed_network$network_data[[1]] <-  toupper(displayed_network$network_data[[1]])
          displayed_network$network_data[[2]] <-  toupper(displayed_network$network_data[[2]])
        } else {
          displayed_network$network_data <- displayed_network$GRN_network_data
        }
      } else {
        displayed_network$GRN_network_data <- NULL
        displayed_network$network_data <- displayed_network$diffGRN_network_data
      }
      # creating the new network
      create_network()
    }
  })

# renders the displayed forcenetwork and modifies some of the underlying JS code (see comments in the string for more details)
  output$net <- renderForceNetwork(
    fnrender <- htmlwidgets::onRender(
      displayed_network$fn,
      '
      function(el, x) {
        d3.selectAll(".link").style("stroke-dasharray", function(d) { 
          // Effect for repressor or inhibitor    
          if (d.effect <= 0){
            return "4px";
          } else {
            return;
          }
        });
        
        // Changing shiny input values to interact with them
        d3.selectAll(".link").on("click", function(d) {
          Shiny.onInputChange("id_target", d.target.name);
          Shiny.onInputChange("id_source", d.source.name);
        });

        // Changing node clickaction based on genecard switch
        d3.selectAll(".node").on("click", function(d) {
          Shiny.onInputChange("id_node", d.name);
          if (d.geneCard_switch) {
            window.open("https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + d.name);
          }
        });

        // Removing mouseover effect from nodes 
        d3.selectAll(".node").on("mouseover", null);

        // Changing edge opacity according to matching edges in boostdiff and uploaded comparison GRN 
        d3.selectAll(".link").style("stroke-opacity", function(d) {
          if (d.comparisonGRNSwitch) {
            if(d.existsInComparisonGRN) {
              return 1;
            } else {
              return .2;
            }
          } else {
            return 1;
          }
        });

        // Change node color to light gray if it is not included in the uploaded comparison GRN
        d3.selectAll(".node").selectAll("circle").style("fill", function(d) {
          if (d.comparisonGRNSwitch) {
            if(d.existsInComparisonGRN) {
              return "rgb(0, 0, 0)";
            } else {
              return "rgb(211,211,211)";
            }
          } else {
            return "rgb(0, 0, 0)";
          }
        });
      }
      '
    )
  )
 
  # renders the geneCard link for the selected gene
  output$geneLink <- renderUI({
    url <- a(sprintf("Genecard %s link", input$id_node), href=sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", input$id_node), target="_blank")
    tagList(sprintf("%s: ", input$id_node), url)
  })

  # loads the data of the uploadable gene list
  observeEvent(input$comparison_grn_file, {
    comparison_network <- read.table(input$comparison_grn_file$datapath, sep = ",", header = TRUE)[, 2:3]
    comparison_nodes <- unique(c(comparison_network[, 1], comparison_network[, 2]))
    displayed_network$fn$x$links$existsInComparisonGRN <- do.call(paste0, displayed_network$links[, 1:2]) %in% do.call(paste0, comparison_network)
    displayed_network$fn$x$nodes$existsInComparisonGRN <- displayed_network$fn$x$nodes$name %in% comparison_nodes 
  })

  # observer events for showing the overlap (node, edges) of the uploadable gene list and the displayed network
  observeEvent(input$comparison_grn_switch_edges, {
    if(is.null(input$comparison_grn_file)) {
      return(NULL)
    }
    displayed_network$fn$x$links$comparisonGRNSwitch <- input$comparison_grn_switch_edges
  })

  observeEvent(input$comparison_grn_switch_nodes, {
    if(is.null(input$comparison_grn_file)) {
      return(NULL)
    }
    displayed_network$fn$x$nodes$comparisonGRNSwitch <- input$comparison_grn_switch_nodes
  })
  
  # function for plotting a linear model of two conditions given the x,y expression data of both
  plot_linear_model <- function(expr1_x, expr1_y, expr2_x, expr2_y) {
    # remove outliers for linear model plotting (|Z_score| > 2 <-> outlier)
    expr1x_mean <- mean(expr1_x)
    expr1y_mean <- mean(expr1_y)
    expr2x_mean <- mean(expr2_x)
    expr2y_mean <- mean(expr2_y)

    expr1x_sd <- sd(expr1_x)
    expr1y_sd <- sd(expr1_y)
    expr2x_sd <- sd(expr2_x)
    expr2y_sd <- sd(expr2_y)

    expr1x_mask <- !(expr1_x > (expr1x_mean + 2*expr1x_sd) | expr1_x < (expr1x_mean - 2*expr1x_sd)) 
    expr1y_mask <- !(expr1_y > (expr1y_mean + 2*expr1y_sd) | expr1_y < (expr1y_mean - 2*expr1y_sd)) 
    expr2x_mask <- !(expr2_x > (expr2x_mean + 2*expr2x_sd) | expr2_x < (expr2x_mean - 2*expr2x_sd))
    expr2y_mask <- !(expr2_y > (expr2y_mean + 2*expr2y_sd) | expr2_y < (expr2y_mean - 2*expr2y_sd))

    expr1_mask <- expr1x_mask & expr1y_mask 
    expr2_mask <- expr2x_mask & expr2y_mask 

    expr1_x <- expr1_x[expr1_mask]
    expr1_y <- expr1_y[expr1_mask]
    expr2_x <- expr2_x[expr2_mask]
    expr2_y <- expr2_y[expr2_mask]

    # prepare plot and plot linear model for the two conditions
    condition_names<-sapply(displayed_network$diffConditions, function(x) gsub(pattern = '_', replacement = ' ', x))
    conditions <- c(rep(condition_names[1], length(expr1_x)), rep(condition_names[2], length(expr2_x)))
    cols <- c(rep(colors[1], length(expr1_x)), rep(colors[2], length(expr2_x)))
    title <- sprintf("%s -> %s || %s vs. %s", input$id_source, input$id_target, condition_names[1], condition_names[2])

    df <- data.frame(
      x = c(expr1_x, expr2_x),
      y = c(expr1_y, expr2_y),
      Condition = conditions,
      col = cols
    )
    x_max <- ceiling(max(df$x))
    x_min <- floor(min(df$x))
    y_max <- ceiling(max(df$y))
    y_min <- floor(min(df$y))

    mapped_cols <- colors[c(1:2)]
    names(mapped_cols) <- condition_names
    plot <- ggplot(df, aes(x = x, y = y, colour = Condition)) +
      geom_point() +
      scale_colour_manual(values = mapped_cols) +
      geom_smooth(method = "lm", formula = "y ~ x", se = FALSE) + 
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))

    plot <- plot + 
      theme_classic() +
      labs(x = "Source gene", y = "Target gene") + 
      theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            plot.title = element_text(size = 15)) +
      scale_x_continuous(limits = c(x_min, x_max)) + scale_y_continuous(limits = c(y_min, y_max))
    return(list(plot=plot, title=title))
  }
  

# Violin and linear model plot for the shown DiffGRN
  output$violin_plot <- renderPlot({
    if (is.null(input$id_node) || opt$mode == 'tsv') {
      return(NULL)
    }
    condition_names <- sapply(displayed_network$diffConditions, function(x) gsub(pattern = '_', replacement = ' ', x))
    cols <- colors[c(1,2)]
    names(cols) <- condition_names
    input_gene_name <- gene_names[grep(paste0('^',input$id_node,'$'), capitilized_gene_names)]
    asub <- subset(displayed_network$adata, subset = ShinyGroup != 'NA')
    plot <- VlnPlot(asub, features = input_gene_name, group.by = 'ShinyGroup', cols = cols, y.max=displayed_network$metacells_maxVal)
    output$downloadGRNViolinPlot <- download_plot(plot, sprintf("%s || %s vs. %s", input$id_node, condition_names[1], condition_names[2]))
    plot
  })

  output$linear_model_plot <- renderPlot({    
    if (is.null(input$id_source)) {
      return(NULL)
    }

    expr1_x <- as.numeric(displayed_network$cond1_metacells_data[displayed_network$cond1_metacells_data$Gene == input$id_source, -1])    
    expr1_y <- as.numeric(displayed_network$cond1_metacells_data[displayed_network$cond1_metacells_data$Gene == input$id_target, -1])    
    expr2_x <- as.numeric(displayed_network$cond2_metacells_data[displayed_network$cond2_metacells_data$Gene == input$id_source, -1])    
    expr2_y <- as.numeric(displayed_network$cond2_metacells_data[displayed_network$cond2_metacells_data$Gene == input$id_target, -1])
    
    res <- plot_linear_model(expr1_x, expr1_y, expr2_x, expr2_y)
    output$downloadGRNLinearModelPlot <- download_plot(res[[1]], res[[2]])

    res[[1]]
  })
  
  get_comp_selections <- function(comp_selec1, comp_selec2) {
    
    col_idx <- 1
    col_names <- unique(lapply(strsplit(comp_selec1, ':'), function(x) { x[1] }))
    comp1_select.cells <- c()
    cells <- c()
    for (selection in comp_selec1) {
      s <- strsplit(selection, ':')[[1]]
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
      val <- str_sub(s[2], 2) # removing the whitespace that is there for good looks :)
      cells <- c(cells, which(adata@meta.data[, col] == val))
    }
    comp1_select.cells <- intersect(comp1_select.cells, cells)

    col_idx <- 1
    col_names <- unique(lapply(strsplit(comp_selec2, ':'), function(x) { x[1] }))
    comp2_select.cells <- c()
    cells <- c()
    for (selection in comp_selec2) {
      s <- strsplit(selection, ':')[[1]]
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
      val <- str_sub(s[2], 2) # removing the whitespace that is there for good looks :)
      cells <- c(cells, which(adata@meta.data[, col] == val))
    }
    comp2_select.cells <- intersect(comp2_select.cells, cells)

    return(list(comp1_select.cells = comp1_select.cells, comp2_select.cells = comp2_select.cells))
  }

  # Comparison violin and linear model plots if a filter was set
  comparison_violin_plot <- eventReactive(input$compare_button, {
    # Can only make comparison plots if the input file type was a Seurat or Scanpy object
    if(opt$mode == "tsv" || is.null(input$id_node)) {
      return(NULL)
    }

    cond_name1 <- input$comp_name_cond1
    cond_name2 <- input$comp_name_cond2  
    condition_names <- c(cond_name1, cond_name2)
    comp_selec1 <- input$comp_selec_cond1
    comp_selec2 <- input$comp_selec_cond2
    if (length(comp_selec1) == 0 || length(comp_selec2) == 0) {
      return(NULL)
    }    
    
    comp_select.cells <- get_comp_selections(comp_selec1, comp_selec2)
    comp1_select.cells <- comp_select.cells$comp1_select.cells
    comp2_select.cells <- comp_select.cells$comp2_select.cells
    
    displayed_network$adata$ShinyComparison[comp1_select.cells] <- cond_name1
    displayed_network$adata$ShinyComparison[comp2_select.cells] <- cond_name2
    input_gene_name <- gene_names[grep(paste0('^',input$id_node,'$'), capitilized_gene_names)]
    cols <- colors[c(1,2)]
    names(cols) <- condition_names
    asub <- subset(displayed_network$adata, subset = ShinyComparison != 'NA')
    plot <- VlnPlot(asub, features = input_gene_name, group.by = 'ShinyComparison', cols = cols, y.max=displayed_network$metacells_maxVal)
    output$downloadComparisonViolinPlot <- download_plot(plot, sprintf("%s || %s vs. %s", input$id_node, condition_names[1], condition_names[2]))
    plot
  })

  output$comparison_violin_plot <- renderPlot({comparison_violin_plot()})
  
  meta_cell_creation <- function(subset) {
      n.samples <- 100 
      p.missing <- 10  
      print(subset)
      n.cells<-nrow(subset@meta.data)
      # Compute number of cells to aggregate
      cells.p.metasample<-nrow(subset@meta.data)/n.samples
      # randomly assign each of the cells to a group
      set.seed(1)
      subset@meta.data$meta.cell<- sample(nrow(subset@meta.data), size = nrow(subset@meta.data), replace = FALSE) %% n.samples
      # Set the ident to the newly created meta.cell variable
      Idents(subset)<-"meta.cell"
      # Aggregate the expression
      print(configuration[configuration$key == input$Key_pick]$assay)
      agg<-AggregateExpression(subset, slot = "counts", return.seurat = T, assays = configuration[configuration$key == input$Key_pick]$assay)
      agg@assays[[configuration[configuration$key == input$Key_pick]$assay]]@counts <-agg@assays[[configuration[configuration$key == input$Key_pick]$assay]]@counts / cells.p.metasample
      agg <- NormalizeData(object = agg, assay = configuration[configuration$key == input$Key_pick]$assay)
      # export the results
      result.data.frame <- agg@assays[[configuration[configuration$key == input$Key_pick]$assay]]@data
      row.names<-rownames(result.data.frame)
      column.names<-paste0(paste0(input$Key_pick, collapse='_'), '_', colnames(result.data.frame))
      result.data.frame<-as.data.table(result.data.frame)
      result.data.frame<-cbind(row.names, result.data.frame)
      colnames(result.data.frame)<-c('Gene', column.names)
      select<-which(rowSums(result.data.frame==0)/(ncol(result.data.frame)-1)<(p.missing/100))
      
      result.data.frame<-result.data.frame[select, ]
      result.data.frame$Gene <- toupper(result.data.frame$Gene)
      return(result.data.frame)
  }

  comparison_linear_model_plot <- eventReactive(input$compare_button, {
    if(opt$mode == "tsv") {
      return(NULL)
    }
    comp_selec1 <- input$comp_selec_cond1
    comp_selec2 <- input$comp_selec_cond2
    cond_name1 <- input$comp_name_cond1
    cond_name2 <- input$comp_name_cond2   
    conditition_names <- c(cond_name1, cond_name2)

    comp_select.cells <- get_comp_selections(comp_selec1, comp_selec2)
    comp1_select.cells <- comp_select.cells[[1]]
    comp2_select.cells <- comp_select.cells[[2]]

    comp1_subset <- subset(adata, cells = comp1_select.cells)
    comp2_subset <- subset(adata, cells = comp2_select.cells) 

    comp1_meta_cell_df <- meta_cell_creation(comp1_subset)
    comp2_meta_cell_df <- meta_cell_creation(comp2_subset)

    df_gene_names <- intersect(comp1_meta_cell_df$Gene, comp2_meta_cell_df$Gene)
    comp1_meta_cell_df <- comp1_meta_cell_df[Gene %in% df_gene_names]
    comp2_meta_cell_df <- comp2_meta_cell_df[Gene %in% df_gene_names]      
    
    comp1_meta_cell_df$Gene <- toupper(comp1_meta_cell_df$Gene)
    comp2_meta_cell_df$Gene <- toupper(comp2_meta_cell_df$Gene)

    expr1_x <- as.numeric(comp1_meta_cell_df[Gene == toupper(isolate(input$id_source)), -1])
    expr1_y <- as.numeric(comp1_meta_cell_df[Gene == toupper(isolate(input$id_target)), -1])    
    expr2_x <- as.numeric(comp2_meta_cell_df[Gene == toupper(isolate(input$id_source)), -1])    
    expr2_y <- as.numeric(comp2_meta_cell_df[Gene == toupper(isolate(input$id_target)), -1])
    
    # outlier detection
    res <- plot_linear_model(expr1_x, expr1_y, expr2_x, expr2_y)
    output$downloadGRNLinearModelPlot <- download_plot(res[[1]], res[[2]])

    res[[1]]
  })

  output$comparison_linear_model_plot <- renderPlot({comparison_linear_model_plot()})

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
    FeaturePlot(adata, features = input$select_genes, ncol=6)
  })
  output$umap_plot <- renderPlot({umap_plot()})

  standard_violin_plot <- eventReactive(input$plot_button, {
    # Seurat's VlnPlot does not show the legend if more than one gene is plotted (see https://github.com/satijalab/seurat/issues/2598) -> using cowplot's plot_grid with list of violin plots
    plots <- VlnPlot(adata, features = input$select_genes, group.by=input$select_primary_grouping,  split.by = input$select_secondary_grouping, ncol=6, combine=FALSE)
    do.call(plot_grid, plots)
    # VlnPlot(adata, features = input$select_genes, group.by=input$select_primary_grouping,  split.by = input$select_secondary_grouping, ncol=6)
  })  
  output$standard_violin_plot <- renderPlot({standard_violin_plot()})

  dot_plot <- eventReactive(input$plot_button, {
    DotPlot(adata, features = input$select_genes, group.by = input$select_primary_grouping)
  })
  output$dot_plot <- renderPlot({dot_plot()})
  
  # updateSelectizeInput(session, "select_genes", choices = gene_names, server = TRUE)

  observeEvent(input$id_node, {
    input_gene_name <- gene_names[grep(paste0('^',input$id_node,'$'), capitilized_gene_names)]
    edit <- isolate(input$select_genes)
    if (length(edit)<12){
      selection<-c(edit, input_gene_name)
    }else{
      selection<-edit
    }
    updateSelectizeInput(session,'select_genes', choices = selection, selected = selection, server = TRUE)
  })
  ##########################
}

shinyApp(ui = ui, server = server) 

