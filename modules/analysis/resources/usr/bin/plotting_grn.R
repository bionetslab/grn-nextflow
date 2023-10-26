#!/usr/bin/env Rscript

################################################################################
# This script constructs a force directed graph visualization of a differential GRN
# using the networkD3 package  
# 
# Params: 
# Input file (tsv file containg the differential GRN)
# Output file (file name of output)
# Output folder (path to output folder)
#
################################################################################

### Imports
# library(networkD3)
library(htmlwidgets)
library(optparse)
library(data.table)

################################################################################

option_list <- list(
  make_option(c("-o", "--output.folder"), type = 'character',
              help="Output folder"), 
  make_option(c("-i", "--input.file"), type = 'character',
              help="Input file"),
  make_option(c("-n", "--output.file"), type = 'character',
              help="Output file name"),
  make_option(c("-p", "--project.folder"), type = 'character',
              help="Path to project base directory")      
        )
opt <- parse_args(OptionParser(option_list=option_list))

library(networkD3, lib.loc=sprintf("%s/lib", opt$project.folder))

colors <- c("#cc79a7", "#009e73", "#0072b2")

tmp <- try(read.table(file = opt$input.file, header = TRUE, sep = '\t'), silent=TRUE)
if (class(tmp) == 'try-error') {
  dir.create(opt$output.folder)
  writeLines(
    paste("<html>\n<head>\n</head>\n<body>\n<h1>",
          "No edges wer found using this approach",
          "</h1>\n</body>\n</html>")
    , opt$output.file)
} else {

  data <- read.table(file = opt$input.file, sep = "\t", header=TRUE)
  links <- data.frame(source = data[, 2], target = data[, 1], width=1, stringsAsFactors=F)
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

  conditions <- unique(data$condition)
  data$condition <- ifelse(data$condition == conditions[1], colors[1], data$condition)

  if (length(conditions) > 1) {
    legend <- htmltools::div(style = "padding: 10px; background-color: white;",
      htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 71px; margin-right: 5px; margin-bottom: 6px", colors[1])),
      sprintf("Stronger in %s (Activator)", conditions[1]),
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
      sprintf("Stronger in %s (Repressor)", conditions[2])
    )
    data$condition <- ifelse(data$condition == conditions[2], colors[2], data$condition)
  } else {
    legend <- htmltools::div(
      style = "padding: 10px; background-color: white;",
      htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 71px; margin-right: 5px; margin-bottom: 6px", colors[1])),
      "Activator",
      htmltools::br(),
      htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
      htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
      htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
      htmltools::div(style = sprintf("display: inline-block; background-color: %s; height: 5px; width: 10px; margin-right: 5px; margin-bottom: 6px", colors[1])),
      "Repressor",
      htmltools::br())
  }
  # create network
  fn <- forceNetwork(Links = network$links, Nodes = network$nodes, 
               Source = 'source', Target = 'target', 
               NodeID = 'name', Group = 'group', Value="width", opacity = 1, arrows=TRUE, 
               linkColour=data$condition, opacityNoHover = 0.6, zoom=TRUE, colourScale = JS('d3.scaleOrdinal().domain(["1"]).range(["#000000"])'))

  # MyClickScript <- 'alert("You clicked " + d.name + " which is in row " +
  #        (d.index + 1) +  " of your original R data frame");'

  fn$x$links$effect <- data[,5]
  fn <- prependContent(fn, legend)

  # add information about inhibition and activation on edges
  fnrender <- onRender(
    fn,
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

    }
    '
  )

  dir.create(opt$output.folder)
  saveNetwork(fnrender, opt$output.file)  
}  



# fnrender <- htmlwidgets::onRender(
#   fn,
#   '

#   //function(el, x) {
#   //  d3.select(el).selectAll(".node")
#   //    .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })

#   //}

#   function(el, x) {
#         d3.select(el).selectAll(".node").attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; });

#     d3.select(el).selectAll(".node text")
#       .attr("dy", ".31em")
#       .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
#       .attr("transform", function(d) { return d.x < 180 ? "translate(8)" : "rotate(180)translate(-8)"; })
#       .style("font", x.options.fontSize + "px " + x.options.fontFamily)
#       .style("opacity", x.options.opacity)
#       .style("fill", x.options.textColour)
#       .text(function(d) { return d.x; });



#     // d3.selectAll(".link").style("marker-end", "url") <- TODO: code a HTML marker that is a flat line or a half circle 
#   }
#   '
# )

# saveNetwork(fn, "testing.html")  



# data <- read.table(file = './diff_network.tsv', sep="\t", header=TRUE)


# links <- data.frame(source = data[, 3], target = data[, 2], width=1, stringsAsFactors=F)
# nodes <- data.frame(name = c(unique(links$source), unique(links$target)), group = 1)
# nodes$id = seq(0, 21)

# # replacing node names with numeric values in links for plotting function
# i <- 0
# for(node in nodes$name) {
#     links$source[links$source == node] = i
#     links$target[links$target == node] = i
#     i <- i + 1
# }

# network <- list(links = links, nodes = nodes)
# data$color <- ifelse(data$color == "darkred", "#ff6db6", data$color)
# data$color <- ifelse(data$color == "darkgreen", "#24ff24", data$color)


# # MyClickScript <- 'alert("You clicked " + d.name + " which is in row " +
# #        (d.index + 1) +  " of your original R data frame");'

# fn <- forceNetwork(Links = network$links, Nodes = network$nodes, 
#              Source = 'source', Target = 'target', 
#              NodeID = 'name', Group = 'group', Value="width", opacity = 1, arrows=TRUE, linkColour=data[,5], opacityNoHover = 0.6, zoom=TRUE, colourScale = JS('d3.scaleOrdinal().domain(["1"]).range(["#000000"])'))

# legend <- htmltools::div(
#   style = "padding: 10px; background-color: white;",
#   htmltools::div(style = "display: inline-block; background-color: #24ff24; height: 2px; width: 10px; margin-right: 5px; margin-bottom: 6px"),
#   "Stronger in Armstrong",
#   htmltools::br(),
#   htmltools::div(style = "display: inline-block; background-color: #ff6db6; height: 2px; width: 10px; margin-right: 5px; margin-bottom: 6px"),
#   "Stronger in Docile"
# )

# fn <- htmlwidgets::prependContent(fn, legend)


# fnrender <- htmlwidgets::onRender(
#   fn,
#   '

#   //function(el, x) {
#   //  d3.select(el).selectAll(".node")
#   //    .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })

#   //}

#   function(el, x) {
#         d3.select(el).selectAll(".node").attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; });

#     d3.select(el).selectAll(".node text")
#       .attr("dy", ".31em")
#       .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
#       .attr("transform", function(d) { return d.x < 180 ? "translate(8)" : "rotate(180)translate(-8)"; })
#       .style("font", x.options.fontSize + "px " + x.options.fontFamily)
#       .style("opacity", x.options.opacity)
#       .style("fill", x.options.textColour)
#       .text(function(d) { return d.x; });



#     // d3.selectAll(".link").style("marker-end", "url") <- TODO: code a HTML marker that is a flat line or a half circle 
#   }
#   '
# )

# saveNetwork(fn, "testing.html")  


# target <- data[, 2]
# src <- data[, 3]

# print(data[data$target=="TOX",])
# print(data[data$src=="PDS5A",])
# # Armstrong:
# # data$color <- ifelse(data$color == "darkred", "#009292", data$color)
# # data$color <- ifelse(data$color == "darkgreen", "#920000", data$color)


# # Docile:
# data$color <- ifelse(data$color == "darkred", "#009292", data$color)
# data$color <- ifelse(data$color == "darkgreen", "#920000", data$color)



# net<-workData <- data.frame(source=src, target=target)
# head(networkData)
# sn <- simpleNetwork(networkData, linkColour=data[,5], nodeColour="#000000", charge=-30, fontSize=12, opacity=0.8, linkDistance = 100, zoom=TRUE)
# # networkData$value <- 1

# # nodes <- data.frame(ID=unique(networkData))
# # nodes$group <- 1

# # fn <- forceNetwork(networkData, nodes, Source = "source", Target = "target", Value = "value", NodeID = "ID", Group = "group", opacity=0.9)
# # saveNetwork(fn, "testing.html")  

# saveNetwork(sn, "net.html")
# # webshot("sn.html", "net.png", vheight=1024, vwidth=768)
