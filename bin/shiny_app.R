#!/usr/bin/env Rscript

################################################################################
### Imports
# Order of imports is important here.
# Seurat must be imported before networkD3 and htmlwidgets.
# Otherwise, the JS function of these libraries will be masked and this leads to errors
library(Seurat)
library(shiny)
library(ggplot2)
library(networkD3)
library(htmlwidgets)
library(optparse)
library(shinyWidgets)
library(DT)
library(data.table)
library(shinythemes)
library(shinyhelper)
library(dplyr)
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
    make_option(c("-a", "--scaNet.folder"), type = 'character',
              help="Path to scaNet data folder", default='')
)
opt <- parse_args(OptionParser(option_list=option_list))

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
dummy_scaNet_data <- network_data[,1:2]
num_rows_to_remove <- round(nrow(dummy_scaNet_data) * 0.8)
rows_to_remove <- sample(nrow(dummy_scaNet_data), num_rows_to_remove, replace = FALSE)
dummy_scaNet_data <- dummy_scaNet_data[-rows_to_remove, ]

# TODO: Load correct ScaNet data
scaNet_arm_expr <- read.table(paste(opt$scaNet.folder, "M1_Arm_GRN.csv", sep = "/"), sep = ",", header = TRUE)[, 2:3]
# # scaNet_arm_expr$condition <- "Arm"
scaNet_doc_expr <- read.table(paste(opt$scaNet.folder, "M1_Doc_GRN.csv", sep = "/"), sep = ",", header = TRUE)[, 2:3]
# scaNet_doc_expr$condition <- "Doc"
scaNet_data <- rbind(scaNet_arm_expr, scaNet_doc_expr)

fn$x$links$existsInScanet <- do.call(paste0, network_data[,1:2]) %in% do.call(paste0, scaNet_data)

###########

########### Loading the necessary tsv files for the gene expression violin plot and the linear model plot of the forcenetwork (done separately from the seurat object because bulking is random)
# reading in gene expression
# gexp_condition1_name <- tail(unlist(strsplit(opt$gexp_condition1.file, "/")), n=1) # removing .tsv
gexp_condition1_name <- substr(opt$metaCells_condition1.file, 1, nchar(opt$metaCells_condition1.file) - 4) # removing .tsv
gexp_condition1_name <- substr(gexp_condition1_name, 5, nchar(gexp_condition1_name))

# gexp_condition2_name <- tail(unlist(strsplit(opt$gexp_condition2.file, "/")), n=1) 
gexp_condition2_name <- substr(opt$metaCells_condition2.file, 1, nchar(opt$metaCells_condition2.file) - 4) # removing .tsv
gexp_condition2_name <- substr(gexp_condition2_name, 5, nchar(gexp_condition2_name)) # removing _out


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
opt$group.var<-strsplit(opt$group.var, ':')[[1]]
Idents(adata)<- opt$group.var

# opt$selection <- strsplit(opt$selection, "-")[[1]]
# opt$group.var <- strsplit(opt$group.var, ':')[[1]]
# Idents(adata_tmp)<- opt$group.var

# cell_subsets <- c()
# subset_seurat_object <- NULL
# # select the correct cells
# for (selection in opt$selection) {
#   criteria <- strsplit(selection, ",")[[1]]
#   index <- length(criteria) - 4 # not the cluster.ids, cluster.name, outputfolder or key
#   selec <- criteria[1:index] 
#   selec <- unname(sapply(selec, function(x) strsplit(x ,":")))
#   # key <- criteria[length(criteria) - 3]
#   cluster.ids <- strsplit(criteria[length(criteria)], ":")[[1]]
#   cluster.name <- criteria[length(criteria)-1]

#   for (g in 1:length(selec)){
#     # Create a subset of the required data
#     select.cells <- list()
#     # first criterion
#     select.cells <- which(adata_tmp@meta.data[, opt$group.var[1]]==selec[[g]][1])
#     # select metadata categories
#     for (i in 2:length(opt$group.var)){
#       select.cells <- intersect(select.cells, which(adata_tmp@meta.data[, opt$group.var[i]]==selec[[g]][i]))
#     }

#     # select clusters (multiple)
#     select.cells<-intersect(select.cells, which(adata_tmp@meta.data[, cluster.name] %in% as.numeric(cluster.ids)))
#     cell_subsets <- c(cell_subsets, select.cells)

#   }
# }      

# subset_seurat_object <- subset(adata, cells = cell_subsets)

# # perform differential testing
# markers <- FindMarkers(subset_seurat_object, ident.1="1", ident.2="2") 
###########
# Shiny app UI and server logic:
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
            materialSwitch(inputId = "geneCard_switch", label = "Enable GeneCard redirection", status = "danger", value = TRUE, right=TRUE),
            materialSwitch(inputId = "scanet_switch", label = "Compare to Scanet results (Opaque edges are not found in Scanet)", status = "danger", right=TRUE),
          ),
          column(3,
            selectInput(
              "selection_dropdown", 
              "Choose case for comparison:", 
              choices = c("NA", "Arm vs Doc D10 Spleen", "Arm vs Doc D28 Spleen", "Arm vs Doc D10 Liver", "Arm vs Doc D28 Liver"),
            )
          ),
          column(3,
            pickerInput(
              inputId = "clusters_pick", 
              label = "Choose Clusters for comparison:", 
              choices = seq(1, 10), 
              options = pickerOptions(
                actionsBox = TRUE, 
                size = 10,
                selectedTextFormat = "count > 3"
              ), 
              multiple = TRUE
            )
          ),
          column(2,
              helper(
                actionButton("compare_button", "Compare plots!"),
                type = "markdown",
                content = "Clickhelp"
              )
          ),              
        ),
        width = 6
      ),
      sidebarPanel(
        fluidRow(
          column(6, plotOutput("violin_plot")),
          column(6, plotOutput("linear_model_plot")),
        ),
        fluidRow(
          column(6, plotOutput("comparison_violin_plot")),      
          column(6, plotOutput("comparison_linear_model_plot")),
        ),
        width = 6
      ), 
    ),
  ),
  tabPanel("Seurat Object Analysis", 
    h2("Seurat Object Analysis")
  ),
)




        # fluidRow(class="network",
        #   column(6, 
        #     # adding legend                
        #     htmltools::div(
        #       style = "padding: 10px; background-color: white;",
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 71px; margin-right: 5px; margin-bottom: 6px", colors[1])),
        #       sprintf("Stronger in %s (Activator)",conditions[1]),
        #       htmltools::br(),
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
        #       sprintf("Stronger in %s (Repressor)", conditions[1]),
        #       htmltools::br(),
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 71px; margin-right: 5px; margin-bottom: 6px", colors[2])),
        #       sprintf("Stronger in %s (Activator)", conditions[2]),
        #       htmltools::br(),
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
        #       htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[2])),
        #       sprintf("Stronger in %s (Repressor)", conditions[2]),
        #       htmltools::br(),
        #     ),
        #     # adding DiffGRN
        #     forceNetworkOutput("net"),
        #     # adding Option switches
        #     fluidRow(class="selection switches",
        #       column(3,
        #         materialSwitch(inputId = "geneCard_switch", label = "Enable GeneCard redirection", status = "danger", value = TRUE, right=TRUE),
        #         materialSwitch(inputId = "scanet_switch", label = "Compare to Scanet results (Opaque edges are not found in Scanet)", status = "danger", right=TRUE),
        #       ),
        #       column(3,
        #         selectInput(
        #           "selection_dropdown", 
        #           "Choose case for comparison:", 
        #           choices = c("NA", "Arm vs Doc D10 Spleen", "Arm vs Doc D28 Spleen", "Arm vs Doc D10 Liver", "Arm vs Doc D28 Liver"),
        #         ),
        #       ),
        #       column(3,
        #         pickerInput(
        #           inputId = "clusters_pick", 
        #           label = "Choose Clusters for comparison:", 
        #           choices = seq(1, 10), 
        #           options = pickerOptions(
        #             actionsBox = TRUE, 
        #             size = 10,
        #             selectedTextFormat = "count > 3"
        #           ), 
        #           multiple = TRUE
        #         )
        #       ),
        #       column(2,
        #           actionButton("compare_button", "Compare plots!")
        #       ),              
        #     ),
        #   ),               
    #       column(3,
    #         fluidRow(class = "table row",
    #           plotOutput("violin_plot")),
    #         fluidRow(class = "table row",
    #           plotOutput("linear_model_plot")),
    #       ),
    #       column(3,
    #         fluidRow(class = "table row",
    #           plotOutput("comparison_violin_plot")),
    #         fluidRow(class = "table row",
    #           plotOutput("comparison_linear_model_plot")),
    #       ),
    #     ),
    #   ),
#     ),
#     tabPanel("Seurat Object Analysis",       
#       h2("Seurat Object Analysis")
#     )
#   ))
# )

server <- function(input, output) {
  observe_helpers(withMathJax = TRUE)
  # The forcenetwork needs to be a reactive value to change the node/link values
  network <- reactiveValues(fn = fn)
  # Observer event for the geneCard switch
  observeEvent(input$geneCard_switch, {
    network$fn$x$nodes$geneCard_switch <- input$geneCard_switch
  })
  # Observer event for the Scanet switch
  observeEvent(input$scanet_switch, {
    network$fn$x$links$scanet_switch <- input$scanet_switch
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
        // Changing edge opacity according to matching edges in boostdiff and scanet iff scanet switch is on
        d3.selectAll(".link").style("stroke-opacity", function(d) {
          if (d.scanet_switch) {
            if(d.existsInScanet) {
              return 1;
            } else {
              return .2;
            }
          } else {
            return 1;
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

    plot <- ggplot(df, aes(x = x, y = y, colour = condition)) +
      geom_point() +
      scale_colour_manual(values = colors) +
      geom_smooth(method = "lm", formula = "y ~ x", se = FALSE) + 
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))

    plot <- plot + theme_classic()
    return(plot)
  }

  plot_violin_plot <- function(expr_data, conditions, title) {
    df <- data.frame(
      gexpr = expr_data,
      condition = conditions
    )

     plot <- ggplot(df, aes(x = condition, y = gexpr, fill = condition)) +
      geom_violin() +
      scale_fill_manual(values = colors) +
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
      xlab("Condition") +
      ylab("Gene expression") 
    # Idea: "save plot switch -> saves all plots that you select"
    # ggsave(filename = "dnjac15_gene_expression_plot.png", plot = plot, width = 6, height = 4, dpi = 300)

    plot <- plot + theme_classic()
    return(plot)
  }
  # render Gene expression plot of the force network
  output$violin_plot <- renderPlot({
    expr1 <- as.numeric(gexp_condition1[gexp_condition1$Gene == input$id_node, -1])
    expr2 <- as.numeric(gexp_condition2[gexp_condition2$Gene == input$id_node, -1])
    expr_data <- c(expr1, expr2)
    conditions <- c(rep(gexp_condition_names[1], length(expr1)), rep(gexp_condition_names[2], length(expr2)))
    title <- sprintf("%s || %s vs %s", input$id_node, gexp_condition_names[1], gexp_condition_names[2])
    
    plot_violin_plot(expr_data, conditions, title)
  }) 

  # render linear model plot of the force network
  output$linear_model_plot <- renderPlot({    
    expr1_x <- as.numeric(gexp_condition1[gexp_condition1$Gene == input$id_source, -1])    
    expr1_y <- as.numeric(gexp_condition1[gexp_condition1$Gene == input$id_target, -1])    
    expr2_x <- as.numeric(gexp_condition2[gexp_condition2$Gene == input$id_source, -1])    
    expr2_y <- as.numeric(gexp_condition2[gexp_condition2$Gene == input$id_target, -1])
    x <- c(expr1_x, expr2_x)
    y <- c(expr1_y, expr2_y)
    conditions <- c(rep(gexp_condition_names[1], length(expr1_x)), rep(gexp_condition_names[2], length(expr2_x)))
    title <- sprintf("%s -> %s || %s vs. %s", input$id_source, input$id_target, gexp_condition_names[1], gexp_condition_names[2])

    plot_linear_model(x, y, conditions, title)
  })

  comparison_data <- reactiveValues(arm_data = NULL, doc_data = NULL, condition_names = NULL)
  # comparing to other Results:
  observeEvent(input$compare_button, {
    comparison_data$doc_data <- NULL
    comparison_data$arm_data <- NULL
    selection <- NULL
    if(input$selection_dropdown == "NA" || is.null(input$clusters_pick)) {
      return()
    } else if (input$selection_dropdown == "Arm vs Doc D10 Spleen") {
      selection <- c("Doc:Spleen:1:d10,Doc:Spleen:3:d10,Doc:Spleen:5:d10", "Arm:Spleen:2:d10,Arm:Spleen:3:d10,Arm:Spleen:4:d10")
      comparison_data$condition_names <- c("Doc_Spleen_d10", "Arm_Spleen_d10")
    } else if (input$selection_dropdown == "Arm vs Doc D28 Spleen") {
      selection <- c("Doc:Spleen:1:d28,Doc:Spleen:2:d28,Doc:Spleen:4:d28", "Arm:Spleen:1:d28,Arm:Spleen:3:d28,Arm:Spleen:5:d28")
      comparison_data$condition_names <- c("Doc_Spleen_d28", "Arm_Spleen_d28")
    } else if (input$selection_dropdown == "Arm vs Doc D10 Liver") {
      selection <- c("Doc:Liver:1:d10,Doc:Liver:3:d10,Doc:Liver:5:d10", "Arm:Liver:2:d10,Arm:Liver:3:d10,Arm:Liver:4:d10")
      comparison_data$condition_names <- c("Doc_Liver_d10", "Arm_Liver_d10")
    } else if (input$selection_dropdown == "Arm vs Doc D28 Liver") {
      selection <- c("Doc:Liver:1:d28,Doc:Liver:2:d28,Doc:Liver:4:d28", "Arm:Liver:1:d28,Arm:Liver:3:d28,Arm:Liver:5:d28")
      comparison_data$condition_names <- c("Doc_Liver_d28", "Arm_Liver_d28")
    } 
    
    clusters <- strsplit(input$clusters_pick, ",")
    cluster.name <- 'cluster' # !!! This only works if the cluster_id is also 'cluster' in boostdiff.nf -> TODO: make it flexible
    n.samples <- 30           # TODO: Same as line above, make it flexible (also hardcoded in boostdiff.nf), however does not break as line above if not flexible and changed
    p.missing <- 10           # Hardcoded in create_metacells.R (not set in boostdiff.nf) 
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
        agg<-AggregateExpression(subset, slot = "data", return.seurat = T, assays = 'SCT')

        # export the results
        result.data.frame<-agg@assays$SCT@data
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
      if(is.null(comparison_data$doc_data)) {
        comparison_data$doc_data <- meta.cell.df     
      } else {
        comparison_data$arm_data <- meta.cell.df
      }
    }
    gene_names <- intersect(comparison_data$doc_data$Gene, comparison_data$arm_data$Gene)
    comparison_data$doc_data <- comparison_data$doc_data[Gene %in% gene_names]
    comparison_data$arm_data <- comparison_data$arm_data[Gene %in% gene_names]      
  })  

  output$comparison_violin_plot <- renderPlot({
    if(!is.null(comparison_data$arm_data)) {
      expr1 <- as.numeric(comparison_data$doc_data[Gene == input$id_node, -1])
      expr2 <- as.numeric(comparison_data$arm_data[Gene == input$id_node, -1])
      expr_data <- c(expr1, expr2)
      conditions <- c(rep(comparison_data$condition_names[1], length(expr1)), rep(comparison_data$condition_names[2], length(expr2)))
      condition_names <- sort(comparison_data$condition_names, decreasing = FALSE) # just to get the same order of conditions names in the title
      title <- sprintf("%s || %s vs. %s", input$id_node, condition_names[1], condition_names[2])

      plot_violin_plot(expr_data, conditions, title)
    }
  })

  output$comparison_linear_model_plot <- renderPlot({
    if(!is.null(comparison_data$arm_data)) {
      expr1_x <- as.numeric(comparison_data$doc_data[Gene == input$id_source, -1])    
      expr1_y <- as.numeric(comparison_data$doc_data[Gene == input$id_target, -1])    
      expr2_x <- as.numeric(comparison_data$arm_data[Gene == input$id_source, -1])    
      expr2_y <- as.numeric(comparison_data$arm_data[Gene == input$id_target, -1])
      x <- c(expr1_x, expr2_x)
      y <- c(expr1_y, expr2_y)
      conditions <- c(rep(comparison_data$condition_names[1], length(expr1_x)), rep(comparison_data$condition_names[2], length(expr2_x)))
      condition_names <- sort(comparison_data$condition_names, decreasing = FALSE) # just to get the same order of conditions names in the title
      title <- sprintf("%s -> %s || %s vs. %s", input$id_source, input$id_target, condition_names[1], condition_names[2])

      plot_linear_model(x, y, conditions, title)
    }
  })

  # rendering Seurat diffexp test plot
  # output$diff_exp_table <- renderDT(
  #   datatable(
  #     markers, rownames = TRUE,
  #     selection = list(selected = 1, target = 'row')                                
  #   )
  # )
}

shinyApp(ui = ui, server = server) 

