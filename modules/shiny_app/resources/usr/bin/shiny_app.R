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

opt <- parse_args(OptionParser(option_list=option_list))
if (opt$grntools != "No tools were chosen for GRN Inference") {
  opt$grntools <- c("NA", opt$grntools)
}
library(networkD3, lib.loc=sprintf("%s/lib", opt$project.path))
colors <- c("#cc79a7", "#009e73", "#0072b2")
configuration <- data.table(key = NULL, assay = NULL, group_var = NULL, condition_names = NULL, tools = NULL, network_files = NULL)

# parse selection
selections <- strsplit(opt$selection, "-")[[1]]
# get condition names first

conditions <- list()
for (s in selections) {
  key <- strsplit(s, ",")[[1]][1]
  condition_name <- strsplit(s, ",")[[1]][4]
  conditions[[key]] <- c(conditions[[key]], condition_name)
}


for (s in selections) {
  split_string <- strsplit(s, ",")[[1]]
  key <- split_string[1]
  assay <- split_string[2]
  group_var <- split_string[3]

  filter <- NULL
  filter_val <- NULL
  
  if (split_string[length(split_string)] == "filter") {
    filter <- split_string[length(split_string) - 2]
    filter_val <- split_string[length(split_string) - 1]
  }

  all_files <- list.files(path = file.path(opt$results.path, key), full.names = TRUE, recursive = TRUE)
  network_files <- c(grep("aggregated_filtered_network", all_files, value = TRUE))
  network_files <- paste(network_files, collapse = ",")

  tools <- list.dirs(path = file.path(opt$results.path, key), full.names = FALSE, recursive = FALSE)
  tools <- paste(tools, collapse = ",")

  cond <- sort(paste(conditions[[key]], collapse = ","), decreasing = FALSE)

  config <- list(key = key, assay = assay, group_var = group_var, tools = tools, condition_names = cond, network_files = network_files, filter = filter, filter_val = filter_val)
  configuration <- rbind(configuration, config)
}

configuration <- unique(configuration)

########### Loading Seurat object, filtering the correct cells and performing differential testing
adata <- readRDS(opt$seurat.file)
# adata <- adata$all

# Set grouping var to Armstrong vs Docile
group.var<-strsplit(configuration[1]$group_var, ':')[[1]]

Idents(adata)<-group.var
#### ADD SOME CONVENIENCE VARIABLES
factor_groups<-names(which(sapply(adata@meta.data, class)=='factor'))
###########

ui <- 
  navbarPage(
    theme = shinytheme("cerulean"), title = "InterNet Xplorer", 
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
              sprintf("Base GRN edge"),
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
                choices = c(opt$dgrntools),
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
              factor_groups,
              selected = adata@meta.data[1],
              multiple = FALSE,
              selectize = TRUE,
              width = NULL,
              size = NULL
            ),
            selectInput(
              'select_secondary_grouping',
              'Select secondary grouping variable',
              factor_groups,
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


  output$comparison_plots <- renderUI({
    if(opt$mode != "seurat"){
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
    if(opt$mode != "seurat") {
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

  displayed_network <- reactiveValues(display_condition = NULL, data_cond1 = NULL, data_cond2 = NULL, diffGRN_network_data = NULL, GRN_network_data = NULL, network_data = NULL, fn = NULL, conditions = NULL, links = NULL)

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
    conditions <- sort(unique(displayed_network$network_data$condition), decreasing = FALSE) # sorting to get matching colours with matching conditions for all plots, needs to be changed later on if other networks are displayed to change legend
    # print(displayed_network$network_data$condition)
    for(i in 1:length(conditions)) {
      displayed_network$network_data$condition <- ifelse(displayed_network$network_data$condition == conditions[i], colors[i], displayed_network$network_data$condition)
      # print(displayed_network$network_data$condition)
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
    print("Input Key pick")
    network.files <- strsplit(configuration[configuration$key == input$Key_pick,]$network_files, ",")[[1]]
    displayed_network$network_data <- read.table(file = network.files[1], header = TRUE)
    if(input$DiffGRN_pick != "No tools were chosen for DGRN Inference") { # -> initial selected network will be a differential network
      displayed_network$diffGRN_network_data <- displayed_network$network_data
    } else {
      displayed_network$GRN_network_data <- displayed_network$network_data
    }
    displayed_network$conditions <- strsplit(configuration[configuration$key == input$Key_pick,]$condition_names, ",")[[1]]
    output$condition1_activator <- renderText({
      paste("Stronger in ", gsub(displayed_network$conditions[1], pattern="_", replacement=" "), " (Activator)")
    })
    output$condition1_repressor <- renderText({
      paste("Stronger in ", gsub(displayed_network$conditions[1], pattern="_", replacement=" "), " (Repressor)")
    })
    output$condition2_activator <- renderText({
      paste("Stronger in ", gsub(displayed_network$conditions[2], pattern="_", replacement=" "), " (Activator)")
    })
    output$condition2_repressor <- renderText({
      paste("Stronger in ", gsub(displayed_network$conditions[2], pattern="_", replacement=" "), " (Repressor)")
    })
    
    output$select_conditions <- renderUI({
      if(opt$mode != "seurat") {
        return(NULL)
      }
      selectizeInput(
        "selection_dropdown",
        label = "Choose 2 conditions to compare to.",
        choices = unique(adata$group), ## Hardcoded change?
        multiple = TRUE,
        options = list(maxItems = 2)
      )
    })
    output$select_filter_vals <- renderUI({
      if(opt$mode != "seurat") {
        return(NULL)
      }
      selectInput(
        'clusters_pick',
        "Choose Clusters for comparison:",
        choices = unique(adata[[configuration[configuration$key == input$Key_pick,]$filter]]),
        selected = NULL,
        multiple = TRUE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      )
    })

    updateTextInput(session, "DiffGRN_pick", value=opt$dgrntools[1])
    updateTextInput(session, "GRN_pick", value=opt$grntools[1])
    create_network()
  })
  
  observeEvent(input$DiffGRN_pick, {
    print("DiffGRN Pick")
    if (input$DiffGRN_pick != "No tools were chosen for DGRN Inference") {
      all_network.files <- strsplit(configuration[configuration$key == input$Key_pick,]$network_files, ",")[[1]]
      # finding the correct network based on input$Key_pick and input$DiffGRN_pick
      network.file <- grep(paste(input$Key_pick, input$DiffGRN_pick, sep="/"), all_network.files, value = TRUE)
      displayed_network$diffGRN_network_data <- read.table(file = network.file, header = TRUE)
      displayed_network$conditions <- strsplit(configuration[configuration$key == input$Key_pick,]$condition_names, ",")[[1]]
      if (!is.null(displayed_network$GRN_network_data)) {
        displayed_network$network_data <- rbind(displayed_network$diffGRN_network_data, displayed_network$GRN_network_data)
        displayed_network$network_data[[1]] <-  str_to_title(displayed_network$network_data[[1]])
        displayed_network$network_data[[2]] <-  str_to_title(displayed_network$network_data[[2]])
      } else {
        displayed_network$network_data <- displayed_network$diffGRN_network_data
      }
      create_network()
    }
  })
  
  observeEvent(input$GRN_pick, {
    print("GRN Pick")
    if (input$GRN_pick != "No tools were chosen for GRN Inference") {  
      if (input$GRN_pick != "NA") {
        all_network.files <- strsplit(configuration[configuration$key == input$Key_pick,]$network_files, ",")[[1]]
        # finding the correct network based on input$Key_pick and input$DiffGRN_pick
        network.file <- grep(paste(input$Key_pick, input$GRN_pick, sep="/"), all_network.files, value = TRUE)
        displayed_network$GRN_network_data <- read.table(file = network.file, header = TRUE)
        ## TODO!: Fix this in the pipeline so that the network matrices all have the same order/name of columns
        # adjusting column order and name to fit boostdiff column names and order 
        displayed_network$GRN_network_data["condition"] <- input$GRN_pick
        displayed_network$GRN_network_data <- displayed_network$GRN_network_data[, c(1,2,3,5,4)]
        colnames(displayed_network$GRN_network_data)[colnames(displayed_network$GRN_network_data) == "TF"] <- "regulator"
        ##
        if (!is.null(displayed_network$diffGRN_network_data)) {
          displayed_network$network_data <- rbind(displayed_network$diffGRN_network_data, displayed_network$GRN_network_data)
          displayed_network$network_data[[1]] <-  str_to_title(displayed_network$network_data[[1]])
          displayed_network$network_data[[2]] <-  str_to_title(displayed_network$network_data[[2]])
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
  
  plot_linear_model <- function(x, y, conditions, title) {
    df <- data.frame(
      x = x,
      y = y,
      Condition = conditions
    )
    if (sum(is.na(df)) > 0.9 * nrow(df)) { # dont plot anything if there are more than 90% NA values in the dataframe
      return(NULL)
    }
    # print(df)
    plot <- ggplot(df, aes(x = x, y = y, colour = Condition)) +
      geom_point() +
      scale_colour_manual(values = unname(colors)) +
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
      scale_x_continuous(limits = c(0, 5.5)) + scale_y_continuous(limits = c(0, 5.5))
    return(plot)
  }
  

# Violin and linear model plot for the shown DiffGRN
  output$violin_plot <- renderPlot({
    if (is.null(input$id_node)) {
      print("Select a node!")
      return(NULL)
    }
    Idents(adata)<-'group' # TODO this is hardcoded potentially problematic
    # asub <- 
    filter_name <- configuration[configuration$key == input$Key_pick,]$filter
    filter_val <- strsplit(configuration[configuration$key == input$Key_pick,]$filter_val, ":")[[1]]
    asub <- adata[, adata@meta.data[[filter_name]] %in% filter_val]
    condition_names<-sapply(displayed_network$conditions, function(x) gsub(pattern = '_', replacement = ' ', x))
    plot<-VlnPlot(asub, features = str_to_title(isolate(input$id_node)), idents = condition_names, cols = colors, y.max=4)
    output$downloadGRNViolinPlot <- download_plot(plot, sprintf("%s || %s vs. %s", input$id_node, condition_names[1], condition_names[2]))
    plot
  })

  output$linear_model_plot <- renderPlot({    
    if (is.null(input$id_source)) {
      print("Select an edge!")
      return(NULL)
    }

    gexp_path <- paste(opt$results.path, input$Key_pick, sep="/")
    all_files <- list.files(gexp_path)
    out_files <- sort(grep("^out", all_files, value = TRUE), decreasing = FALSE)

    gexp_condition1 <- read.table(file = paste(gexp_path, out_files[1], sep="/"), header = TRUE, sep = "\t")
    gexp_condition2 <- read.table(file = paste(gexp_path, out_files[2], sep="/"), header = TRUE, sep = "\t")

    expr1_x <- as.numeric(gexp_condition1[gexp_condition1$Gene == input$id_source, -1])    
    expr1_y <- as.numeric(gexp_condition1[gexp_condition1$Gene == input$id_target, -1])    
    expr2_x <- as.numeric(gexp_condition2[gexp_condition2$Gene == input$id_source, -1])    
    expr2_y <- as.numeric(gexp_condition2[gexp_condition2$Gene == input$id_target, -1])
    x <- c(expr1_x, expr2_x)
    y <- c(expr1_y, expr2_y)
    
    condition_name1 <- gsub(displayed_network$conditions[1], pattern="_", replacement=" ")
    condition_name2 <- gsub(displayed_network$conditions[2], pattern="_", replacement=" ")

    conditions <- c(rep(condition_name1, length(expr1_x)), rep(condition_name2, length(expr2_x)))
    title <- sprintf("%s -> %s || %s vs. %s", input$id_source, input$id_target, condition_name1, condition_name2)

    plot <- plot_linear_model(x, y, conditions, title)
    output$downloadGRNLinearModelPlot <- download_plot(plot, title)

    plot
  })
  

  ### Comparison Plots ###
  output$comparison_violin_plot <- renderPlot({
    # Can only make comparison plots if the input file type was a Seurat object
    if(opt$mode != "seurat") {
      return(NULL)
    }
    if (is.null(input$id_node) || 
        is.null(input$selection_dropdown) || 
        is.null(input$clusters_pick) ||
        length(input$selection_dropdown) < 2 
    ) { # only plot if all three values are set and atleast 2 cases are selected
      print("Select a node, two cases in the selection dropdown menu and atleast one cluster!")
      return(NULL)
    }
    
    filter_name <- configuration[configuration$key == input$Key_pick,]$filter
    filter_val <- input$clusters_pick

    Idents(adata) <- 'group' # TODO this is hardcoded potentially problematic
    asub <- adata[, adata@meta.data[[filter_name]] %in% filter_val]
    condition_names<-sapply(isolate(input$selection_dropdown), function(x) gsub(pattern = '_', replacement = ' ', x))
    condition_names <- sort(condition_names, decreasing=FALSE)
    plot <- VlnPlot(asub, features = str_to_title(isolate(input$id_node)), idents = condition_names, cols = colors, y.max=4)
    output$downloadComparisonViolinPlot <- download_plot(plot, sprintf("%s || %s vs. %s", input$id_node, condition_names[1], condition_names[2]))
    plot
  })
  

  output$comparison_linear_model_plot <- renderPlot({
    if(opt$mode != "seurat") {
      return(NULL)
    }
    if (is.null(input$id_source) || 
        is.null(input$selection_dropdown) || 
        is.null(input$clusters_pick) ||
        length(input$selection_dropdown) < 2 
    ) { # only plot if all three values are set
      print("Select an edge, two cases in the selection dropdown menu and atleast one cluster!")
      return(NULL)
    }
    input$compare_button

    cond1_data <- NULL
    cond2_data <- NULL
    
    condition_names <- isolate(input$selection_dropdown)

    selection<-c()
    for(i in isolate(input$selection_dropdown)){
      relevant<- as.character(unique(as.data.table(adata@meta.data)[group==i]$sample))
      relevant<- sapply(relevant, function(x) gsub(x, pattern=' ', replacement=':'))
      relevant<-paste0(relevant, collapse = ',')
      selection<-c(selection, relevant)
    }
    print(selection)

    print(input$clusters_pick)
    filter_val <- strsplit(input$clusters_pick, ",")[[1]]
    filter_name <- configuration[configuration$key == input$Key_pick,]$filter
    n.samples <- opt$n.samples 
    p.missing <- 50            
    
    group.var <- strsplit(configuration[configuration$key == input$Key_pick,]$group_var, ":")[[1]]
    assay <- configuration[configuration$key == input$Key_pick,]$assay
    # print(clusters)
    for (s in selection) {
      s <- strsplit(s, ',')
      s <- unname(sapply(s, function(x) strsplit(x ,":")))
      
      meta.cell.dfs<-list()
      meta.cell.df<-NULL
      for (g in 1:length(s)){
        # Create a subset of the required data
        select.cells<-list()
        # first criterion
        select.cells <- which(adata@meta.data[, group.var[1]]==s[[g]][1])
        
        # select metadata categories
        for (i in 2:length(group.var)){
          select.cells<-intersect(select.cells, which(adata@meta.data[, group.var[i]]==s[[g]][i]))
        }
        
        # select clusters (multiple)
        select.cells<-intersect(select.cells, which(adata@meta.data[, filter_name] %in% as.numeric(filter_val)))
        
        
        subset<-subset(adata, cells = select.cells)
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
        agg<-AggregateExpression(subset, slot = "data", return.seurat = T, assays = assay)

        # export the results
        result.data.frame<-agg@assays[[assay]]@data
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
      }
      select<-which(rowSums(meta.cell.df==0)/(ncol(meta.cell.df)-1)<(p.missing/100))
      meta.cell.df<-meta.cell.df[select, ]
      if(is.null(cond1_data)) {
        cond1_data <- meta.cell.df     
      } else {
        cond2_data <- meta.cell.df
      }
    }

    gene_names <- intersect(cond1_data$Gene, cond2_data$Gene)
    cond1_data <- cond1_data[Gene %in% gene_names]
    cond2_data <- cond2_data[Gene %in% gene_names]      

    expr1_x <- as.numeric(cond1_data[Gene == isolate(input$id_source), -1])
    expr1_y <- as.numeric(cond1_data[Gene == isolate(input$id_target), -1])    
    expr2_x <- as.numeric(cond2_data[Gene == isolate(input$id_source), -1])    
    expr2_y <- as.numeric(cond2_data[Gene == isolate(input$id_target), -1])
    
    x <- c(expr1_x, expr2_x)
    y <- c(expr1_y, expr2_y)
    conditions <- c(rep(condition_names[1], length(expr1_x)), rep(condition_names[2], length(expr2_x)))
    condition_names <- sort(condition_names, decreasing = FALSE) # just to get the same order of conditions names in the title
    title <- sprintf("%s -> %s || %s vs. %s", input$id_source, input$id_target, condition_names[1], condition_names[2])
      
    plot <- plot_linear_model(x, y, conditions, title)

    output$downloadComparisonLinearModelPlot <- download_plot(plot, title)
    plot
  })

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
    x = str_to_title(input$id_node)
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

