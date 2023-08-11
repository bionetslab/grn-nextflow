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
  make_option(c("-s", "--selection"), type = 'character', 
              help="Selection criteria separated by colon"),
  make_option(c("-l", "--cluster.name"), type = 'character', 
              help="Selection criteria separated by colon"),
  make_option(c("-k", "--cluster.ids"), type = 'character', 
              help="cluster ids"),
  make_option(c("-g", "--group.var"), type = 'character', 
              help="Grouping variable(s) in meta.data object separated by colon",
              default = 'sample'),
  make_option(c("-o", "--output.file"), type = 'character',
              help="Output file"),
  make_option(c("-a", "--assay"), type = 'character',
              help="Assay used"))

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
# opt$input.file<-'/home/bionets-og86asub/Documents/netmap/data/ga_an0228_10x_deepseq_filtered_smarta_merged_tissue_integrated_rep_timepoint_infection_filtered_seurat.rds'
# n.samples<-100
# opt$group.var<- "infection:tissue:subject:time"
# opt$selection<- "Doc:Spleen:1:d10,Doc:Spleen:3:d10"
# opt$cluster.name<-'cluster'
# opt$cluster.ids<-'1:2'
# opt$clusters<-c(1,2)
# opt$output.file<-'/home/bionets-og86asub/Documents/netmap/data/Doc_Spleen_d10.tsv'



print('Selected parameters')
print(opt$input.file)
print(opt$n.samples)
print(opt$selection)
print(opt$group.var)
print(opt$output.file)
print(opt$assay)
print(opt$cluster.name)

if (opt$cluster.name == "null" || opt$cluster.ids == "null") {
  opt$cluster.name <- NULL
  opt$cluster.ids <- NULL
}

n.samples<-opt$n.samples
opt$group.var<-strsplit(opt$group.var, ':')[[1]]
opt$selection<-strsplit(opt$selection, ',')[[1]]
opt$selection<-unname(sapply(opt$selection, function(x) strsplit(x ,":")))
if (!is.null(opt$cluster.name) && !is.null(opt$cluster.ids)) {
  opt$cluster.ids<-strsplit(opt$cluster.ids, ':')[[1]]
}
#load data and extract seurat object
adata<-readRDS(opt$input.file)
# adata<-adata$all


# Set the ident to the required identity class
Idents(adata)<- opt$group.var
meta.cell.df<-NULL
# For each group:
for (g in 1:length(opt$selection)){
# Create a subset of the required data
select.cells<-list()
# first criterion
select.cells <- which(adata@meta.data[, opt$group.var[1]]==opt$selection[[g]][1])
# select metadata categories
for (i in 2:length(opt$group.var)){
  select.cells<-intersect(select.cells, which(adata@meta.data[, opt$group.var[i]]==opt$selection[[g]][i]))
}

# select clusters (multiple)
if (!is.null(opt$cluster.name) && !is.null(opt$cluster.ids)) {
  select.cells<-intersect(select.cells, which(adata@meta.data[, opt$cluster.name] %in% as.numeric(opt$cluster.ids)))
}

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
column.names<-paste0(paste0(opt$selection[[g]], collapse='_'), '_', colnames(result.data.frame))
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

select<-which(rowSums(meta.cell.df==0)/(ncol(meta.cell.df)-1)<(opt$p.missing/100))
meta.cell.df<-meta.cell.df[select, ]
# save the aggregated data frame into a tsv sheet
fwrite(meta.cell.df, file = file.path(opt$output.file), sep='\t')