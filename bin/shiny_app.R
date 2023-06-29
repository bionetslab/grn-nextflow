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
library(networkD3)
#(htmlwidgets)
library(optparse)
library(shinyWidgets)
library(DT)
library(data.table)
library(shinythemes)
library(shinyhelper)
library(dplyr)
library(cowplot)
################################################################################

option_list <- list(
  make_option(c("-n", "--network.file"), type = 'character',
              default='', help="Network data file"),
  make_option(c("-s", "--selection"), type = 'character',
              default='', help="selection (hardcoded parameter in boostdiff.nf)"),
  make_option(c("-c", "--metaCells_condition1.file"), type = 'character',
              default='', help="Gene expression file of meta cells of condition 1"),
  make_option(c("-d", "--metaCells_condition2.file"), type = 'character',
              default='', help="Gene expression file of of meta cells of condition 2"),
  make_option(c("-f", "--seurat.file"), type = 'character',
              default='', help="Path to Seurat file"),
  make_option(c("-g", "--group.var"), type = 'character', 
              help="Grouping variable(s) in meta.data object separated by colon",
              default = 'sample'),
    make_option(c("-a", "--assay"), type = 'character',
              help = "Assay to use", default = 'SCT'))
opt <- parse_args(OptionParser(option_list=option_list))


# opt$network.file<- '/home/bionets-og86asub/Documents/netmap/workflows/results/Arm_vs_Doc_D28:Spleen/aggregated_filtered_network.txt' 
# opt$metaCells_condition1.file<- '/home/bionets-og86asub/Documents/netmap/workflows/results/Arm_vs_Doc_D28:Spleen/out_Doc_Spleen_d28.tsv' 
# opt$metaCells_condition2.file<- '/home/bionets-og86asub/Documents/netmap/workflows/results/Arm_vs_Doc_D28:Spleen/out_Arm_Spleen_d28.tsv' 
# opt$seurat.file<-'/home/bionets-og86asub/Documents/netmap/data/misc/seurat_cd4_micro.rds'
# opt$group.var<-'infection:tissue:subject:time' 
# opt$selection<- 'Doc:Spleen:1:d28,Doc:Spleen:2:d28,Doc:Spleen:4:d28,Doc_Spleen_d28,Arm_vs_Doc_D28:Spleen,cluster,1:2-Arm:Spleen:1:d28,Arm:Spleen:3:d28,Arm:Spleen:5:d28,Arm_Spleen_d28,Arm_vs_Doc_D28:Spleen,cluster,1:2' 
# opt$scaNet.folder<-'/home/bionets-og86asub/Documents/netmap/workflows/scaNet_data'


########### Creating force directed network
# Loading differential network
network_data <- read.table(file = opt$network.file, sep = "\t", header = TRUE)
# Creating force network
links <- data.frame(source = network_data[, 2], target = network_data[, 1], width=network_data[, 3]*100, stringsAsFactors=F)
nodes <- data.frame(name = unique(c(unique(links$source), unique(links$target))), group = 1)
nodes$id = seq(0, nrow(nodes)-1)
# replacing node names with numeric values in links for plotting function
i <- 0
for(node in nodes$name) {
  links$source[links$source == node] = i
  links$target[links$target == node] = i
  i <- i + 1
}

network <- list(links = links, nodes = nodes)
colors <- c("#cc79a7", "#009e73")

conditions <- sort(unique(network_data$condition), decreasing = FALSE) # sorting to get matching colours with matching conditions for all plots
network_data$condition <- ifelse(network_data$condition == conditions[1], colors[1], network_data$condition)
network_data$condition <- ifelse(network_data$condition == conditions[2], colors[2], network_data$condition)

fn <- forceNetwork(Links = network$links, Nodes = network$nodes, 
                   Source = 'source', Target = 'target', 
                   NodeID = 'name', Group = 'group', Value="width", opacity = 1, arrows=TRUE, 
                   linkColour=network_data$condition, opacityNoHover = 0.6, zoom=TRUE, colourScale = JS('d3.scaleOrdinal().domain(["1"]).range(["#000000"])'))

# adding switch information to access it later on
fn$x$nodes$geneCard_switch <- FALSE
fn$x$links$scanet_switch <- FALSE
# adding linear model information to force network
fn$x$links$effect <- network_data[,5]

# scaNet dummy test
# dummy_scaNet_data <- network_data[,1:2]
# num_rows_to_remove <- round(nrow(dummy_scaNet_data) * 0.8)
# rows_to_remove <- sample(nrow(dummy_scaNet_data), num_rows_to_remove, replace = FALSE)
# dummy_scaNet_data <- dummy_scaNet_data[-rows_to_remove, ]

# # TODO: Load correct ScaNet data
# scaNet_arm_expr <- read.table(paste(opt$scaNet.folder, "M1_Arm_GRN.csv", sep = "/"), sep = ",", header = TRUE)[, 2:3]
# # # scaNet_arm_expr$condition <- "Arm"
# scaNet_doc_expr <- read.table(paste(opt$scaNet.folder, "M1_Doc_GRN.csv", sep = "/"), sep = ",", header = TRUE)[, 2:3]
# # scaNet_doc_expr$condition <- "Doc"
# scaNet_data <- rbind(scaNet_arm_expr, scaNet_doc_expr)

# fn$x$links$existsInScanet <- do.call(paste0, network_data[,1:2]) %in% do.call(paste0, scaNet_data)

###########

########### Loading the necessary tsv files for the gene expression violin plot and the linear model plot of the forcenetwork (done separately from the seurat object because bulking is random)
# reading in gene expression
# gexp_condition1_name <- tail(unlist(strsplit(opt$gexp_condition1.file, "/")), n=1) # removing .tsv

gexp_condition1_name <- gsub(basename(opt$metaCells_condition1.file), pattern="\\.tsv$|\\.csv)", replacement="")
gexp_condition1_name <- gsub(gexp_condition1_name, pattern="^out_", replacement="")

gexp_condition2_name <- gsub(basename(opt$metaCells_condition2.file), pattern="\\.tsv$|\\.csv)", replacement="")
gexp_condition2_name <- gsub(gexp_condition2_name, pattern="^out_", replacement="")


names(colors)<-c(gsub(pattern = '_', replacement = ' ',  gexp_condition2_name), gsub(pattern = '_', replacement = ' ',  gexp_condition1_name))

gexp_condition_names <- sort(c(gexp_condition1_name, gexp_condition2_name), decreasing = FALSE) # just to get the matching ordering with matching colours :)
if (gexp_condition_names[1] == gexp_condition1_name) {  
  gexp_condition1 <- read.table(file = opt$metaCells_condition1.file, sep = "\t", header = TRUE)
  gexp_condition2 <- read.table(file = opt$metaCells_condition2.file, sep = "\t", header = TRUE)
} else {
  gexp_condition1 <- read.table(file = opt$metaCells_condition2.file, sep = "\t", header = TRUE)
  gexp_condition2 <- read.table(file = opt$metaCells_condition1.file, sep = "\t", header = TRUE)
}
###########

########### Loading Seurat object, filtering the correct cells and performing differential testing
adata <- readRDS(opt$seurat.file)
adata <- adata$all

# Set grouping var to Armstrong vs Docile
opt$group.var<-strsplit(opt$group.var, ':')[[1]]

Idents(adata)<-opt$group.var


#### ADD SOME CONVENIENCE VARIABLES
factor_groups<-names(which(sapply(adata@meta.data, class)=='factor'))
# print(factor_groups) 
###########

ui <- 
  navbarPage(theme = shinytheme("cerulean"), title = "Results", 
             tabPanel("Boostdiff Results",
                      sidebarLayout(
                        mainPanel(
                          fluidRow(class="network",
                                   # adding legend                
                                   htmltools::div(
                                     style = "padding: 10px; background-color: white;",
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 71px; margin-right: 5px; margin-bottom: 6px", colors[1])),
                                     sprintf("Stronger in %s (Activator)",conditions[1]),
                                     htmltools::br(),
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
                                     sprintf("Stronger in %s (Repressor)", conditions[1]),
                                     htmltools::br(),
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 71px; margin-right: 5px; margin-bottom: 6px", colors[2])),
                                     sprintf("Stronger in %s (Activator)", conditions[2]),
                                     htmltools::br(),
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
                                     htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
                                     sprintf("Stronger in %s (Repressor)", conditions[2]),
                                     htmltools::br(),
                                   ),
                                   # adding DiffGRN
                                   forceNetworkOutput("net"),
                          ),
                          # adding Option switches
                          fluidRow(class="selection switches",
                                    column(3,
                                          checkboxInput(inputId = "geneCard_switch", label = "Enable GeneCard redirection", value = FALSE),
                                          helper(
                                            fileInput('comparison_grn_file', 'Choose GRN to compare to',
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
                                          selectizeInput(
                                            "selection_dropdown",
                                            label = "Choose 2 conditions to compare to.",
                                            choices = unique(adata$group), ## Hardcoded change?
                                            multiple = TRUE,
                                            options = list(maxItems = 2)
                                          ),
                                          selectInput(
                                            'clusters_pick',
                                            "Choose Clusters for comparison:",
                                            choices = unique(adata$cluster),
                                            selected = NULL,
                                            multiple = TRUE,
                                            selectize = TRUE,
                                            width = NULL,
                                            size = NULL
                                          ),
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
                                      ),
                          ),
                          width = 6
                        ),
                        sidebarPanel(
                          fluidRow(
                            column(6, 
                              plotOutput("violin_plot"),
                              downloadButton("downloadGRNViolinPlot", "Download Plot")
                            ),
                            column(6, 
                              plotOutput("linear_model_plot"),
                              downloadButton("downloadGRNLinearModelPlot", "Download Plot")
                              ),
                          ),
                          fluidRow(
                            column(6, 
                              plotOutput("comparison_violin_plot"),
                              downloadButton("downloadComparisonViolinPlot", "Download Plot")
                              ),      
                            column(6, 
                              plotOutput("comparison_linear_model_plot"),
                              downloadButton("downloadComparisonLinearModelPlot", "Download Plot")
                              ),
                          ),
                          width = 6
                        ), 
                      ),
             ),
             tabPanel("Seurat Object Analysis", 
                      h2("Seurat Object Analysis"),
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
                          tabsetPanel(type = 'tabs',
                                      tabPanel('Umap',
                                               fluidRow(12, 
                                                        column(12, plotOutput("umap_plot")), 
                                                        column(12, plotOutput('standard_violin_plot'))
                                               )
                                      ) ,
                                      tabPanel('Dotplot',
                                               fluidRow(12, 
                                                        column(12, plotOutput("dot_plot"))
                                               )
                                      )      
                          )
                          
                          
                        )
                      )
             ),
  )

server <- function(input, output, session) {
  observe_helpers(withMathJax = TRUE)
  # The forcenetwork needs to be a reactive value to change the node/link values
  network <- reactiveValues(fn = fn)
  # Observer event for the geneCard switch
  observeEvent(input$geneCard_switch, {
    network$fn$x$nodes$geneCard_switch <- input$geneCard_switch
  })
  # Observer event for the Scanet switch

  observeEvent(input$comparison_grn_file, {
    comparison_network <- read.table(input$comparison_grn_file$datapath, sep = ",", header = TRUE)[, 2:3]
    comparison_nodes <- unique(c(comparison_network[, 1], comparison_network[, 2]))
    network$fn$x$links$existsInComparisonGRN <- do.call(paste0, network$fn$x$links[, 1:2]) %in% do.call(paste0, comparison_network)
    network$fn$x$nodes$existsInComparisonGRN <- network$fn$x$nodes$name %in% comparison_nodes 
  })

  observeEvent(input$comparison_grn_switch_edges, {
    if(is.null(input$comparison_grn_file)) {
      return(NULL)
    }
    network$fn$x$links$comparisonGRNSwitch <- input$comparison_grn_switch_edges
  })

  observeEvent(input$comparison_grn_switch_nodes, {
    if(is.null(input$comparison_grn_file)) {
      return(NULL)
    }
    network$fn$x$nodes$comparisonGRNSwitch <- input$comparison_grn_switch_nodes
  })
  # render force network
  output$net <- renderForceNetwork(
    fnrender <- htmlwidgets::onRender(
      network$fn,
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
  
  plot_linear_model <- function(x, y, conditions, title) {
    df <- data.frame(
      x = x,
      y = y,
      condition = conditions
    )
    if (sum(is.na(df)) > 0.9 * nrow(df)) { # dont plot anything if there are more than 90% NA values in the dataframe
      return(NULL)
    }
    # print(df)
    plot <- ggplot(df, aes(x = x, y = y, colour = condition)) +
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
    condition_names<-sapply(gexp_condition_names, function(x) gsub(pattern = '_', replacement = ' ', x))
    plot<-VlnPlot(adata, features = isolate(input$id_node), idents = condition_names, cols = colors, y.max=4)
    output$downloadGRNViolinPlot <- download_plot(plot, sprintf("%s || %s vs. %s", input$id_node, condition_names[1], condition_names[2]))
    plot
  })

  output$linear_model_plot <- renderPlot({    
    if (is.null(input$id_source)) {
      print("Select an edge!")
      return(NULL)
    }
    expr1_x <- as.numeric(gexp_condition1[gexp_condition1$Gene == input$id_source, -1])    
    expr1_y <- as.numeric(gexp_condition1[gexp_condition1$Gene == input$id_target, -1])    
    expr2_x <- as.numeric(gexp_condition2[gexp_condition2$Gene == input$id_source, -1])    
    expr2_y <- as.numeric(gexp_condition2[gexp_condition2$Gene == input$id_target, -1])
    x <- c(expr1_x, expr2_x)
    y <- c(expr1_y, expr2_y)
    conditions <- c(rep(gexp_condition_names[1], length(expr1_x)), rep(gexp_condition_names[2], length(expr2_x)))
    title <- sprintf("%s -> %s || %s vs. %s", input$id_source, input$id_target, gexp_condition_names[1], gexp_condition_names[2])

    plot <- plot_linear_model(x, y, conditions, title)
    output$downloadGRNLinearModelPlot <- download_plot(plot, title)

    plot
  })
  

  ### Comparison Plots ###
  output$comparison_violin_plot <- renderPlot({
    if (is.null(input$id_node) || 
        is.null(input$selection_dropdown) || 
        is.null(input$clusters_pick) ||
        length(input$selection_dropdown) < 2 
    ) { # only plot if all three values are set and atleast 2 cases are selected
      print("Select a node, two cases in the selection dropdown menu and atleast one cluster!")
      return(NULL)
    }
    input$compare_button

    cluster.name <- 'cluster' # !!! This only works if the cluster_id is also 'cluster' in boostdiff.nf -> TODO: make it flexible
    Idents(adata)<-'group' # TODO this is hardcoded potentially problematic
    asub<-subset(adata, subset = cluster %in% isolate(input$clusters_pick))
    condition_names<-sapply(isolate(input$selection_dropdown), function(x) gsub(pattern = '_', replacement = ' ', x))
    condition_names <- sort(condition_names, decreasing=FALSE)
    cols <- colors
    names(cols) <- c(gsub(pattern = '_', replacement = ' ',  condition_names[1]), gsub(pattern = '_', replacement = ' ',  condition_names[2])) 
    plot <- VlnPlot(asub, features = isolate(input$id_node), idents = condition_names, cols = cols, y.max=4)
    output$downloadComparisonViolinPlot <- download_plot(plot, sprintf("%s || %s vs. %s", input$id_node, condition_names[1], condition_names[2]))
    plot
  })
  

  output$comparison_linear_model_plot <- renderPlot({
    if (is.null(input$id_source) || 
        is.null(input$selection_dropdown) || 
        is.null(input$clusters_pick) ||
        length(input$selection_dropdown) < 2 
    ) { # only plot if all three values are set
      print("Select an edge, two cases in the selection dropdown menu and atleast one cluster!")
      return(NULL)
    }
    input$compare_button
    doc_data <- NULL
    arm_data <- NULL
    
    condition_names <- isolate(input$selection_dropdown)

    selection<-c()
    for(i in isolate(input$selection_dropdown)){
      relevant<- as.character(unique(as.data.table(adata@meta.data)[group==i]$sample))
      relevant<- sapply(relevant, function(x) gsub(x, pattern=' ', replacement=':'))
      relevant<-paste0(relevant, collapse = ',')
      selection<-c(selection, relevant)
    }
    # print(selection)
    
    clusters <- isolate(input$clusters_pick)
    cluster.name <- 'cluster' # !!! This only works if the cluster_id is also 'cluster' in boostdiff.nf -> TODO: make it flexible
    n.samples <- 30           # TODO: Same as line above, make it flexible (also hardcoded in boostdiff.nf), however does not break as line above if not flexible and changed
    p.missing <- 10           # Hardcoded in create_metacells.R (not set in boostdiff.nf) 
    
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
        select.cells <- which(adata@meta.data[, opt$group.var[1]]==s[[g]][1])
        # select metadata categories
        for (i in 2:length(opt$group.var)){
          select.cells<-intersect(select.cells, which(adata@meta.data[, opt$group.var[i]]==s[[g]][i]))
        }
        
        # select clusters (multiple)
        select.cells<-intersect(select.cells, which(adata@meta.data[, cluster.name] %in% as.numeric(clusters)))
        
        
        subset<-subset(adata, cells = select.cells)
        # Find number of barcodes in object
        n.cells<-nrow(subset@meta.data)
        # Compute number of cells to aggregate
        cells.p.metasample<-nrow(subset@meta.data)/n.samples
        # randomly assign each of the cells to a group
        subset@meta.data$meta.cell<- sample(nrow(subset@meta.data), size = nrow(subset@meta.data), replace = FALSE) %% n.samples
        # Set the ident to the newly created meta.cell variable
        Idents(subset)<-"meta.cell"
        # Aggregate the expression
        agg<-AggregateExpression(subset, slot = "data", return.seurat = T, assays = opt$assay)

        # export the results
        result.data.frame<-agg@assays[[opt$assay]]@data
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
      if(is.null(doc_data)) {
        doc_data <- meta.cell.df     
      } else {
        arm_data <- meta.cell.df
      }
    }
    
    gene_names <- intersect(doc_data$Gene, arm_data$Gene)
    doc_data <- doc_data[Gene %in% gene_names]
    arm_data <- arm_data[Gene %in% gene_names]      
    
    expr1_x <- as.numeric(doc_data[Gene == isolate(input$id_source), -1])    
    expr1_y <- as.numeric(doc_data[Gene == isolate(input$id_target), -1])    
    expr2_x <- as.numeric(arm_data[Gene == isolate(input$id_source), -1])    
    expr2_y <- as.numeric(arm_data[Gene == isolate(input$id_target), -1])
    
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

