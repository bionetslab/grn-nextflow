#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(dcanr)
library(igraph)
library(EBcoexpress)

option_list <- list(
    make_option(c("-c", "--input.file.1"), type = 'character', 
                help='Input file'),
    make_option(c("-d", "--input.file.2"), type = 'character',
                help='Input file'),
    make_option(c("-o", "--output.path"), type = 'character',
                help='output path')
)

opt <- parse_args(OptionParser(option_list=option_list))

opt$input.file.1 <- 'work/ed/6890859d5432ca859982d05ce9da5e/out_Arm_Spleen_d10.tsv'
opt$input.file.2 <- 'work/ed/6890859d5432ca859982d05ce9da5e/out_Doc_Spleen_d10.tsv'

data.1 <- read.table(file = opt$input.file.1, sep = "\t", header = TRUE)
data.1 <- as.data.frame(sapply(data.1, as.numeric))
data.2 <- read.table(file = opt$input.file.2, sep = "\t", header = TRUE)
data.2 <- as.data.frame(sapply(data.2, as.numeric))


conditions <- c(rep(1, times=ncol(data.1)), rep(2, times=ncol(data.2)))
data <- cbind(data.1, data.2)

scores <- dcScore(data, conditions, dc.method = 'ebcoexpress', cor.method = 'spearman', ebcoexpress.useBWMC=FALSE)
raw_p <- dcTest(scores, data, conditions)
adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')
dcnet <- dcNetwork(scores, adj_p)
edgedf <- as_data_frame(dcnet, what = 'edges')
print(head(edgedf))
