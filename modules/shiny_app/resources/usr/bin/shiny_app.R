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
  split_string <- strsplit(s, ",")[[1]]
  key <- split_string[1]
  assay <- split_string[2]
  group_var <- split_string[3]
  split_group_var <- strsplit(group_var, ":")[[1]]
  # indicates which group_vars are used in each selection of the comparison
  group_var_idx <- paste(match(split_group_var, all_group_vars), collapse=":")

  # if no filter was used, add NULL values.
  filter <- NA
  filter_val <- NA
  # if the filter is set, the selection string was modified (by the pipeline) to have filter as last keyword
  if (split_string[length(split_string)] == "filter") {
    group_vals <- split_string[5:(length(split_string) - 3)]

    filter <- split_string[length(split_string) - 2]
    filter_val <- split_string[length(split_string) - 1]
  } else {
    group_vals <- split_string[5:(length(split_string) - 2)]
  }
  # parsing all information
  group_vals <- paste(group_vals, collapse=",")
  all_files <- list.files(path = file.path(opt$results.path, key), full.names = TRUE, recursive = TRUE)
  network_files <- c(grep("aggregated_filtered_network", all_files, value = TRUE))
  network_files <- paste(network_files, collapse = ",")

  tools <- list.dirs(path = file.path(opt$results.path, key), full.names = FALSE, recursive = FALSE)
  tools <- paste(tools, collapse = ",")

  cond <- sort(paste(conditions[[key]], collapse = ","), decreasing = FALSE)

  config <- list(key = key, assay = assay, group_var = group_var, group_var_idx = group_var_idx, group_vals = group_vals, tools = tools, condition_names = cond, network_files = network_files, filter = filter, filter_val = filter_val)
  configuration <- rbind(configuration, config)
}

# Changing list of data.table into 
configuration <- unique(configuration)
duplicates <- configuration[duplicated(configuration$key),]
# Must be because of differing group_vars ->
if (nrow(duplicates) > 0) {
  for (i in 1:nrow(duplicates)) {
    dup_row <- configuration[i*2,]
    orig_row <- configuration[(i*2) - 1, ]

    configuration[(i*2) - 1, ]$group_vals <- paste(orig_row$group_vals, dup_row$group_vals, sep = "-")
    configuration[(i*2) - 1, ]$group_var_idx <- paste(orig_row$group_var_idx, dup_row$group_var_idx, sep = "-")
    configuration[(i*2) - 1, ]$group_var <- paste(orig_row$group_var, dup_row$group_var, sep=":")
  }
  # removing duplicate rows has to be done after merging the group_vars
  configuration <- configuration[-sapply(1:nrow(duplicates), function(x) {x * 2}),]
}
########### Loading Seurat object, filtering the correct cells and performing differential testing
if (opt$mode != "tsv") {
  adata <- readRDS(opt$seurat.file)
  # adata <- adata$all
  group.var<-unique(strsplit(configuration[1]$group_var, ':')[[1]])
  Idents(adata)<-group.var
}
ui <- 
  navbarPage(
    theme = shinytheme("cerulean"), title = "InterNet Xplorer", 
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
            # adding legend                
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
            # adding DiffGRN
            forceNetworkOutput("net"),
            hr(style = "border-top: 0.5px solid #000000; opacity: 0.2;")
          ),
          # adding Option switches
          fluidRow(class="selection switches",
            # h3("Options:"),
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
              uiOutput("select_conditions"),
              uiOutput("select_filter_vals"),
              uiOutput("compare_button"),
            ),
          ),
          width = 6
        ),
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
    tabPanel(
      "Gene Expression", 
      h2("Gene Expression"),
        sidebarLayout(
          sidebarPanel(
            selectInput(
              'select_genes',
              'Select Genes to display',
              rownames(adata),
              selected = NULL,
              multiple = TRUE,
              selectize = TRUE,
              width = NULL,
              size = NULL
            ),
            selectInput(
              'select_primary_grouping',
              'Select primary grouping variable',
              colnames(adata@meta.data),
              selected = adata@meta.data[1],
              multiple = FALSE,
              selectize = TRUE,
              width = NULL,
              size = NULL
            ),
            selectInput(
              'select_secondary_grouping',
              'Select secondary grouping variable',
              colnames(adata@meta.data),
              selected = adata@meta.data[1],
              multiple = FALSE,
              selectize = TRUE,
              width = NULL,
              size = NULL
            )
          ),
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

server <- function(input, output, session) {
  observe_helpers(withMathJax = TRUE)
#   # The forcenetwork needs to be a reactive value to change the node/link values
  output$mode <- renderText({
    paste0(opt$mode)
  })

  output$comparison_plots <- renderUI({
    if(opt$mode == "tsv" | is.na(configuration[configuration$key == input$Key_pick,]$filter)){
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

  output$compare_button <- renderUI({
    if(opt$mode == "tsv" | is.na(configuration[configuration$key == input$Key_pick,]$filter)) {
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

  displayed_network <- reactiveValues(
    diffGRN_network_data = NULL, GRN_network_data = NULL, network_data = NULL, 
    fn = NULL, links = NULL, diffConditions = "", grnConditions = "", adata = adata, 
    cond1_metacells_data = NULL, cond2_metacells_data = NULL, metacells_minVal = 0, metacells_maxVal = 0
  )

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
  
  ## Observer event for selecting the comparison key 
  observeEvent(input$Key_pick, {
    # setting base network if key is switched
    network.file <- strsplit(configuration[configuration$key == input$Key_pick,]$network_files, ",")[[1]][1]
    # finding the correct network based on input$Key_pick and input$DiffGRN_pick
    tmp <- try(read.table(file = network.file, header = TRUE))
    if (class(tmp) == "try-error") {
      showNotification(paste0("No edges were found for the condition ", input$DiffGRN_pick, "! Displaying network of previously selected condition."), type = "warning")
    } else {
      displayed_network$network_data <- tmp
      if(input$DiffGRN_pick != "No tools were chosen for DGRN Inference") { # -> initial selected network will be a differential network
        displayed_network$diffGRN_network_data <- displayed_network$network_data
        displayed_network$diffConditions <- strsplit(configuration[configuration$key == input$Key_pick,]$condition_names, ",")[[1]]

        # creating a column in the adata object to match the conditions to the cells
        group_var.all <- strsplit(configuration[configuration$key == input$Key_pick,]$group_var, ":")[[1]]
        group_vals.all <- strsplit(configuration[configuration$key == input$Key_pick,]$group_vals, "-")[[1]]
        group_idx.all <- strsplit(configuration[configuration$key == input$Key_pick,]$group_var_idx, "-")[[1]]

        displayed_network$adata$group <- "NA"
        for (i in 1:length(group_idx.all)) {
          group_idx <- as.numeric(strsplit(group_idx.all[i], ":")[[1]])
          group_var <- group_var.all[c(group_idx)]
          group_val.selection <- strsplit(group_vals.all[i], ",")[[1]]
          filter_vector <- nrow(adata)
          idx <- 1
          for (col_name in group_var) {
            col_values <- adata@meta.data[, col_name]
            group_val <- strsplit(group_val.selection, ":")[[1]][idx]
            filter_vector <- filter_vector & (col_values == group_val)
            idx <- idx + 1
          }
          displayed_network$adata$group[filter_vector] <- gsub(pattern = '_', replacement = ' ', displayed_network$diffConditions[i])        
        }
      } else { # TODO: Fix grn only mode
        displayed_network$GRN_network_data <- displayed_network$network_data
      }
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
      
      output$select_conditions <- renderUI({
        if(opt$mode == "tsv" || is.na(configuration[configuration$key == input$Key_pick,]$filter)) {
          return(NULL)
        } else {
          selectInput(
            "filter_key_pick",
            label = "Choose a key to compare to.",
            choices = configuration$key, 
            multiple = FALSE,
            width = NULL,
            size = NULL
          )
        }
      })
      output$select_filter_vals <- renderUI({
        if(opt$mode == "tsv" || is.na(configuration[configuration$key == input$Key_pick,]$filter)) {
          return(NULL)
        } else {
          selectInput(
            'filter_value',
            label = "Choose Clusters for comparison:",
            choices = unique(adata[[configuration[configuration$key == input$Key_pick,]$filter]]),
            selected = NULL,
            multiple = TRUE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          )
        }
      })

      updateTextInput(session, "DiffGRN_pick", value=opt$dgrntools[1])
      updateTextInput(session, "GRN_pick", value=opt$grntools[1])

      gexp_path <- paste(opt$results.path, input$Key_pick, sep="/")
      all_files <- list.files(gexp_path)
      out_files <- grep("^out", all_files, value = TRUE)
      out_files1_condition_name <- gsub(pattern="out", "", out_files[1])
      out_files1_condition_name <- gsub(pattern=".tsv", "", out_files1_condition_name)
      if (out_files1_condition_name[1] == displayed_network$diffConditions[1]) {
        displayed_network$cond1_metacells_data <- read.table(file = paste(gexp_path, out_files[1], sep="/"), header = TRUE, sep = "\t")
        displayed_network$cond2_metacells_data <- read.table(file = paste(gexp_path, out_files[2], sep="/"), header = TRUE, sep = "\t")
      } else {
        displayed_network$cond2_metacells_data <- read.table(file = paste(gexp_path, out_files[1], sep="/"), header = TRUE, sep = "\t")
        displayed_network$cond1_metacells_data <- read.table(file = paste(gexp_path, out_files[2], sep="/"), header = TRUE, sep = "\t")
      }
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
  
  observeEvent(input$DiffGRN_pick, {
    if (input$DiffGRN_pick != "No tools were chosen for DGRN Inference") {
      all_network.files <- strsplit(configuration[configuration$key == input$Key_pick,]$network_files, ",")[[1]]
      # finding the correct network based on input$Key_pick and input$DiffGRN_pick
      network.file <- grep(paste(input$Key_pick, input$DiffGRN_pick, sep="/"), all_network.files, value = TRUE)
      tmp <- try(read.table(file = network.file, header = TRUE))
      if (class(tmp) == "try-error") {
        showNotification(paste0("No edges were found for the condition ", input$DiffGRN_pick, "! Displaying network of previously selected condition."), type = "warning")
      } else {
        displayed_network$diffGRN_network_data <- tmp
        displayed_network$diffConditions <- strsplit(configuration[configuration$key == input$Key_pick,]$condition_names, ",")[[1]]
        group_var.all <- strsplit(configuration[configuration$key == input$Key_pick,]$group_var, ":")[[1]]
        group_vals.all <- strsplit(configuration[configuration$key == input$Key_pick,]$group_vals, "-")[[1]]
        group_idx.all <- strsplit(configuration[configuration$key == input$Key_pick,]$group_var_idx, "-")[[1]]

        displayed_network$adata$group <- "NA"
        for (i in 1:length(group_idx.all)) {
          group_idx <- as.numeric(strsplit(group_idx.all[i], ":")[[1]])
          group_var <- group_var.all[c(group_idx)]
          group_val.selection <- strsplit(group_vals.all[i], ",")[[1]]
          filter_vector <- nrow(adata)
          idx <- 1
          for (col_name in group_var) {
            col_values <- adata@meta.data[, col_name]
            group_val <- strsplit(group_val.selection, ":")[[1]][idx]
            filter_vector <- filter_vector & (col_values == group_val)
            idx <- idx + 1
          }
          displayed_network$adata$group[filter_vector] <- gsub(pattern = '_', replacement = ' ', displayed_network$diffConditions[i])        
        }
        if (!is.null(displayed_network$GRN_network_data)) {
          displayed_network$network_data <- rbind(displayed_network$diffGRN_network_data, displayed_network$GRN_network_data)
          displayed_network$network_data[[1]] <-  str_to_title(displayed_network$network_data[[1]])
          displayed_network$network_data[[2]] <-  str_to_title(displayed_network$network_data[[2]])
        } else {
          displayed_network$network_data <- displayed_network$diffGRN_network_data
        }
        create_network()
      }
    }
  })
  
  observeEvent(input$GRN_pick, {
    if (input$GRN_pick != "No tools were chosen for GRN Inference") {  
      if (input$GRN_pick != "NA") {
        all_network.files <- strsplit(configuration[configuration$key == input$Key_pick,]$network_files, ",")[[1]]
        # finding the correct network based on input$Key_pick and input$DiffGRN_pick
        network.file <- grep(paste(input$Key_pick, input$GRN_pick, sep="/"), all_network.files, value = TRUE)
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
      create_network()
    }
  })

#    ### render displayed force network ###
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
 
  output$geneLink <- renderUI({
    url <- a(sprintf("Genecard %s link", input$id_node), href=sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", input$id_node), target="_blank")
    tagList(sprintf("%s: ", input$id_node), url)
  })

  observeEvent(input$comparison_grn_file, {
    comparison_network <- read.table(input$comparison_grn_file$datapath, sep = ",", header = TRUE)[, 2:3]
    comparison_nodes <- unique(c(comparison_network[, 1], comparison_network[, 2]))
    displayed_network$fn$x$links$existsInComparisonGRN <- do.call(paste0, displayed_network$links[, 1:2]) %in% do.call(paste0, comparison_network)
    displayed_network$fn$x$nodes$existsInComparisonGRN <- displayed_network$fn$x$nodes$name %in% comparison_nodes 
  })

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
  
  plot_linear_model <- function(expr1_x, expr1_y, expr2_x, expr2_y) {
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

    if (sum(is.na(df)) > 0.9 * nrow(df)) { # dont plot anything if there are more than 90% NA values in the dataframe
      return(NULL)
    }
    cols <- df$col
    names(cols) <- df$Condition
    plot <- ggplot(df, aes(x = x, y = y, colour = Condition)) +
      geom_point() +
      scale_colour_manual(values = cols) +
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
    if (is.null(input$id_node)) {
      return(NULL)
    }
    Idents(displayed_network$adata) <- 'group'
    # asub <- 
    filter_name <- configuration[configuration$key == input$Key_pick,]$filter
    if (!is.na(filter_name)) {
      filter_val <- strsplit(configuration[configuration$key == input$Key_pick,]$filter_val, ":")[[1]]
      asub <- displayed_network$adata[, displayed_network$adata@meta.data[[filter_name]] %in% filter_val]
    } else {
      asub <- displayed_network$adata
    }
    condition_names <- sapply(displayed_network$diffConditions, function(x) gsub(pattern = '_', replacement = ' ', x))
    cols <- colors[c(1,2)]
    names(cols) <- condition_names
    adata_gene_name <- rownames(asub)[which(toupper(rownames(asub)) == input$id_node)]
    plot <- VlnPlot(asub, features = adata_gene_name, idents = condition_names, cols = cols, y.max=displayed_network$metacells_maxVal)
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
  

#   ### Comparison Plots ###
  comparison_violin_plot <- eventReactive(input$compare_button, {
    # Can only make comparison plots if the input file type was a Seurat or Scanpy object
    if(opt$mode == "tsv" || is.na(configuration[configuration$key == input$Key_pick,]$filter)) {
      return(NULL)
    }
    if (is.null(input$id_node) || 
        is.null(input$filter_key_pick) || 
        is.null(input$filter_value)
    ) { # only plot if all three values are set and atleast 2 cases are selected
      return(NULL)
    }
    
    filter_name <- configuration[configuration$key == input$Key_pick,]$filter
    picked_filter_key <- configuration[configuration$key == input$filter_key_pick,]$filter
    if (filter_name != picked_filter_key) {
      showNotification(paste0("Selected filter key has not the same set filter in the configuration file. Therefore, these conditions cannot be compared!"), type = "error")
      return(NULL)
    }
    filter_val <- input$filter_value
    group_var.all <- strsplit(configuration[configuration$key == input$filter_key_pick,]$group_var, ":")[[1]]
    group_vals.all <- strsplit(configuration[configuration$key == input$filter_key_pick,]$group_vals, "-")[[1]]
    group_idx.all <- strsplit(configuration[configuration$key == input$filter_key_pick,]$group_var_idx, "-")[[1]]
    condition_names <- strsplit(configuration[configuration$key == input$filter_key_pick,]$condition_names, ",")[[1]]
    condition_names <- sapply(condition_names, function(x) gsub(pattern = '_', replacement = ' ', x))
    adata$comparison <- "NA"
    for (i in 1:length(group_idx.all)) {
      group_idx <- as.numeric(strsplit(group_idx.all[i], ":")[[1]])
      group_var <- group_var.all[c(group_idx)]
      group_val.selection <- strsplit(group_vals.all[i], ",")[[1]]
      filter_vector <- nrow(adata)
      idx <- 1
      for (col_name in group_var) {
        col_values <- adata@meta.data[, col_name]
        group_val <- strsplit(group_val.selection, ":")[[1]][idx]
        filter_vector <- filter_vector & (col_values == group_val)
        idx <- idx + 1
      }
      adata$comparison[filter_vector] <- condition_names[i]        
    }

    Idents(adata) <- 'comparison'
    asub <- subset(adata, cells = which(adata@meta.data[, filter_name] %in% as.numeric(filter_val)))
    adata_gene_name <- rownames(asub)[which(toupper(rownames(asub)) == input$id_node)]
    plot <- try(VlnPlot(asub, features = adata_gene_name, idents = condition_names, cols = colors, y.max=displayed_network$metacells_maxVal), silent=TRUE)
    output$downloadComparisonViolinPlot <- download_plot(plot, sprintf("%s || %s vs. %s", input$id_node, condition_names[1], condition_names[2]))
    plot
  })

  output$comparison_violin_plot <- renderPlot({comparison_violin_plot()})
  
  comparison_linear_model_plot <- eventReactive(input$compare_button, {
    if(opt$mode == "tsv" || is.na(configuration[configuration$key == input$Key_pick,]$filter)) {
      return(NULL)
    }
    if (is.null(input$id_source) || 
        is.null(input$filter_key_pick) || 
        is.null(input$filter_value)
    ) { # only plot if all three values are set
      return(NULL)
    }


    filter_name <- configuration[configuration$key == input$Key_pick,]$filter
    picked_filter_key <- configuration[configuration$key == input$filter_key_pick,]$filter
    if (filter_name != picked_filter_key) {
      showNotification(paste0("Selected filter key has not the same set filter in the configuration file. Therefore, these conditions cannot be compared!"), type = "error")
      return(NULL)
    }
    filter_val <- input$filter_value
    group_var.all <- strsplit(configuration[configuration$key == input$filter_key_pick,]$group_var, ":")[[1]]
    selection <- strsplit(configuration[configuration$key == input$filter_key_pick,]$group_vals, "-")[[1]]
    group_idx.all <- strsplit(configuration[configuration$key == input$filter_key_pick,]$group_var_idx, "-")[[1]]
    condition_names <- strsplit(configuration[configuration$key == input$filter_key_pick,]$condition_names, ",")[[1]]
    condition_names <- sapply(condition_names, function(x) gsub(pattern = '_', replacement = ' ', x))
    assay <- configuration[configuration$key == input$filter_key_pick,]$assay
    n.samples <- 100 
    p.missing <- 50            

    cond1_data <- NULL
    cond2_data <- NULL
    idx <- 1
    for (i in seq_len(length(selection))) {
      s <- selection[i]
      s <- strsplit(s, ',')[[1]]
      s <- lapply(s, function(x) strsplit(x, ":")[[1]])

      group_idx <- as.numeric(strsplit(group_idx.all[idx], ":")[[1]])
      group_var <- group_var.all[c(group_idx)]
      idx <- idx + 1 

      meta.cell.dfs<-list()
      meta.cell.df<-NULL
      select.cells.all <- c()
      for (g in 1:length(s)){
        select.cells<-list()
        # first criterion
        select.cells <- which(adata@meta.data[, group_var[1]]==s[[g]][1])

        # select metadata categories
        if (length(group_var) > 1) {
          for (j in 2:length(group_var)){
            select.cells<-intersect(select.cells, which(adata@meta.data[, group_var[j]]==s[[g]][j]))
          }
        }
        select.cells<-intersect(select.cells, which(adata@meta.data[, filter_name] %in% as.numeric(filter_val)))
        select.cells.all <- c(select.cells.all, unlist(select.cells))
      }
      select.cells.all <- unique(select.cells.all)
      subset<-subset(adata, cells = select.cells.all)        
      # Find number of barcodes in object
      n.cells<-nrow(subset@meta.data)
      # Compute number of cells to aggregate
      cells.p.metasample<-nrow(subset@meta.data)/n.samples
      # randomly assign each of the cells to a group
      set.seed(1)
      subset@meta.data$meta.cell<- sample(nrow(subset@meta.data), size = nrow(subset@meta.data), replace = FALSE) %% n.samples
      # Set the ident to the newly created meta.cell variable
      Idents(subset)<-"meta.cell"
      # Aggregate the expression
      agg<-AggregateExpression(subset, slot = "counts", return.seurat = T, assays = opt$assay)
      agg@assays[[opt$assay]]@counts <-agg@assays[[opt$assay]]@counts / cells.p.metasample
      agg <- NormalizeData(object = agg, assay = opt$assay)
      # export the results

      result.data.frame <- agg@assays[[opt$assay]]@data
      row.names<-rownames(result.data.frame)
      column.names<-paste0(paste0(s[[g]], collapse='_'), '_', colnames(result.data.frame))
      result.data.frame<-as.data.table(result.data.frame)
      result.data.frame<-cbind(row.names, result.data.frame)
      colnames(result.data.frame)<-c('Gene', column.names)
      if(is.null(meta.cell.df)){
        meta.cell.df<-result.data.frame
      }
      else{
        meta.cell.df<-merge(meta.cell.df, result.data.frame, by = 'Gene')
      }
      select<-which(rowSums(meta.cell.df==0)/(ncol(meta.cell.df)-1)<(p.missing/100))
      meta.cell.df<-meta.cell.df[select, ]
      meta.cell.df$Gene <- toupper(meta.cell.df$Gene)
      if(is.null(cond1_data)) {
        cond1_data <- meta.cell.df     
      } else {
        cond2_data <- meta.cell.df
      }
    }
    
    gene_names <- intersect(cond1_data$Gene, cond2_data$Gene)
    cond1_data <- cond1_data[Gene %in% gene_names]
    cond2_data <- cond2_data[Gene %in% gene_names]      

    expr1_x <- as.numeric(cond1_data[Gene == toupper(isolate(input$id_source)), -1])
    expr1_y <- as.numeric(cond1_data[Gene == toupper(isolate(input$id_target)), -1])    
    expr2_x <- as.numeric(cond2_data[Gene == toupper(isolate(input$id_source)), -1])    
    expr2_y <- as.numeric(cond2_data[Gene == toupper(isolate(input$id_target)), -1])
    
    # outlier detection
     
    res <- plot_linear_model(expr1_x, expr1_y, expr2_x, expr2_y)
    output$downloadGRNLinearModelPlot <- download_plot(res[[1]], res[[2]])

    res[[1]]
  })

  output$comparison_linear_model_plot <- renderPlot({comparison_linear_model_plot()})

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

  ##########################
  output$umap_plot<-renderPlot({
    FeaturePlot(adata, features = input$select_genes, ncol=6)
  })
  
  output$standard_violin_plot<-renderPlot({
    # Seurat's VlnPlot does not show the legend if more than one gene is plotted (see https://github.com/satijalab/seurat/issues/2598) -> using cowplot's plot_grid with list of violin plots
    plots <- VlnPlot(adata, features = input$select_genes, group.by=input$select_primary_grouping,  split.by = input$select_secondary_grouping, ncol=6, combine=FALSE)
    do.call(plot_grid, plots)
    # VlnPlot(adata, features = input$select_genes, group.by=input$select_primary_grouping,  split.by = input$select_secondary_grouping, ncol=6)
  })
  
  output$dot_plot<-renderPlot({
    DotPlot(adata, features = input$select_genes, group.by = input$select_primary_grouping)
  })
  
  observe({
    x = input$id_node
    edit <- isolate(input$select_genes)
    if (length(edit)<12){
      selection<-c(edit, x)
    }else{
      selection<-edit
    }
    updateSelectInput(session,'select_genes',
                      choices = rownames(adata),
                      selected = selection)
  })
  ##########################
}

shinyApp(ui = ui, server = server) 

