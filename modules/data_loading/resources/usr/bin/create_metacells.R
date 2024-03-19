#!/usr/bin/env Rscript


################################################################################
# This script extracts scaled data matrices from a Seurat object for different 
# subgroups. Users need to specify the relevant meta variabiable(s) in
# meta.data and the corresponding selection criterion.
# 
# The data is aggregated to form meta-cells sampled at random from a specified
# subgroup/cluster as specified in the input. The user can specify a target 
# number of cells (default: 100 per group). 
#
# The data is saved as a dense data frame of the format genes x metacells (tab 
# separated). The first column name will be called 'Gene'. Genes without valid
# measures will be 0-imputed.
# 
# Params: Input file name (rds file containing Seurat object)
# Output filename (optional)
# Selector and selection criteria
# Number of target meta cells.
#
################################################################################

### Imports
require(Seurat)
require(optparse)
require(tidyverse)
require(data.table)

################################################################################


option_list <- list( 
  make_option(c("-f", "--input.file"), type = 'character',
              help="Input file"),
  make_option(c("-n", "--n.samples"), type="integer", default=100, 
              help="Number of meta cells to generate",
              metavar="number"),
  make_option(c("-p", "--p.missing"), type="integer", default=10, 
              help="Percentage of 0 allowed per gene",
              metavar="number"),
  make_option(c("-s", "--selection"), type = 'character', default="", 
              help="Selection criteria separated by colon"),
  make_option(c("-g", "--group.var"), type = 'character',
              help="Grouping variable(s) in meta.data object separated by colon",
              default = 'sample'),
  make_option(c("-o", "--output.file"), type = 'character',
              help="Output file"),
  make_option(c("-a", "--assay"), type = 'character', default="",
              help="Assay used"),
  make_option(c("-m", "--mode"), type = 'character',
              help="data loading mode. Available modes: 'tsv', 'seurat', 'anndata'", default="seurat"),
  make_option(c("--key"), type = 'character',
              help="Key that specifies the selection"),
  make_option(c('--bulk_RNA'), type = 'logical',
              default = T, help = 'If bulk RNA is used (provide Normalized count matrix in tsv mode)'),
  make_option(c('--create_metacells'), type = 'logical',
              default = T, help = 'Determines if metacells are created')        
  )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
opt$input.file<-'/home/bionets-og86asub/Documents/netmap/data/misc/ko_vs_wt.rds'
n.samples<-100
opt$group.var<- "cluster:genotype"
opt$selection<- "8:wt,8:ko"
opt$cluster.name<-'cluster'
opt$cluster.ids<-'1:2'
opt$clusters<-c(1,2)
opt$output.file<-'/home/bionets-og86asub/Documents/netmap/data/Doc_Spleen_d10.tsv'
opt$assay<-'SCT'
n.samples<-opt$n.samples
opt$key<-'test'

if (opt$mode == "seurat" || opt$mode == "anndata") {

  print('Selected parameters')
  print(opt$input.file)
  print(opt$n.samples)
  print(opt$selection)
  print(opt$group.var)
  print(opt$output.file)
  print(opt$assay)

  opt$group.var<-strsplit(opt$group.var, ':')[[1]]
  opt$selection<-strsplit(opt$selection, ',')[[1]]
  opt$selection<-unname(sapply(opt$selection, function(x) strsplit(x ,":")))
  #load data and extract seurat object
  adata<-readRDS(opt$input.file)
  #adata<-adata$all

  # Set the ident to the required identity class
  Idents(adata)<- opt$group.var
  DefaultAssay(adata)<-opt$assay
  meta.cell.df<-NULL
  # For each group:
  select.cells.all <- c()
  for (g in 1:length(opt$selection)){
    # Create a subset of the required data
    select.cells<-list()
    # first criterion
    select.cells <- which(adata@meta.data[, opt$group.var[1]]==opt$selection[[g]][1])
    # select metadata categories
    if (length(opt$group.var) > 1) {
      for (i in 2:length(opt$group.var)){
        select.cells<-intersect(select.cells, which(adata@meta.data[, opt$group.var[i]]==opt$selection[[g]][i]))
      }
    }
    # print(select.cells)
    select.cells.all <- c(select.cells.all, unlist(select.cells))
  }  
  subset<-subset(adata, cells = select.cells.all)
  if (opt$create_metacells) {
    # Find number of barcodes in object
    n.cells<-nrow(subset@meta.data)
    # Compute number of cells to aggregate
    cells.p.metasample<-nrow(subset@meta.data)/n.samples
    # randomly assign each of the cells to a group
    set.seed(1)
    subset@meta.data$meta.cell<- sample(nrow(subset@meta.data), size = nrow(subset@meta.data), replace = FALSE) %% n.samples
    # Set the ident to the newly created meta.cell variable
    Idents(subset)<-"meta.cell"
    # Aggregate the count expression
    agg<-AggregateExpression(subset, return.seurat = T, assays = opt$assay)
    agg[[opt$assay]]@layers$counts <- agg[[opt$assay]]@layers$counts / cells.p.metasample
    agg <- NormalizeData(agg)
  } else {
    agg <- NormalizeData(subset)
  }
  result.data.frame <- agg[[opt$assay]]@layers$data
  
  row.names<-rownames(agg)
  column.names<-paste0(paste0(opt$key, collapse='_'), '_', colnames(agg))
  result.data.frame<-as.data.table(result.data.frame)
  result.data.frame<-cbind(row.names, result.data.frame)
  colnames(result.data.frame)<-c('Gene', column.names)
  # filter out genes that have 1 - opt$p.missing counts that are 0
  select <- which(rowSums(result.data.frame==0)/(ncol(result.data.frame)-1)<(opt$p.missing/100))
  result.data.frame <- result.data.frame[select, ]
  # save the aggregated data frame into a tsv file
  fwrite(result.data.frame, file = file.path(opt$output.file), sep='\t')
} else if (opt$mode == "tsv") {

  expression_matrix <- read.table(opt$input.file, header = TRUE, row.names = 1, sep = "\t")

  subset <- CreateSeuratObject(counts = expression_matrix)  
  if (opt$bulk_RNA) {
    result.data.frame <- subset[['RNA']]@counts
    row.names<-rownames(result.data.frame)
    column.names<-paste0(paste0(opt$key, collapse='_'), '_', colnames(result.data.frame))
    result.data.frame<-as.data.table(result.data.frame)
    result.data.frame<-cbind(row.names, result.data.frame)
    colnames(result.data.frame)<-c('Gene', column.names)

    select<-which(rowSums(result.data.frame==0)/(ncol(result.data.frame)-1)<(opt$p.missing/100))
    result.data.frame<-result.data.frame[select, ]
    # save the aggregated data frame into a tsv sheet
    fwrite(result.data.frame, file = file.path(opt$output.file), sep='\t')
  } else {
    # Find number of barcodes in object
    n.cells <- nrow(subset)
    # Compute number of cells to aggregate
    cells.p.metasample <- n.cells/n.samples
    # randomly assign each of the cells to a group
    set.seed(1)
    subset@meta.data$meta.cell <- sample(nrow(subset@meta.data), size = nrow(subset@meta.data), replace = FALSE) %% n.samples
    # Set the ident to the newly created meta.cell variable
    Idents(subset) <- "meta.cell"
    # Aggregate the expression
    agg<-AggregateExpression(subset, slot = "counts", return.seurat = T, assays = 'RNA') 
    agg[['RNA']]@counts <- agg[['RNA']]@counts / cells.p.metasample
    result.data.frame <- agg[['RNA']]@counts
    row.names<-rownames(result.data.frame)
    column.names<-paste0(paste0(opt$key, collapse='_'), '_', colnames(result.data.frame))
    result.data.frame<-as.data.table(result.data.frame)
    result.data.frame<-cbind(row.names, result.data.frame)
    colnames(result.data.frame)<-c('Gene', column.names)

    select<-which(rowSums(result.data.frame==0)/(ncol(result.data.frame)-1)<(opt$p.missing/100))
    result.data.frame<-result.data.frame[select, ]
    # save the aggregated data frame into a tsv sheet
    fwrite(result.data.frame, file = file.path(opt$output.file), sep='\t')
  }

} 