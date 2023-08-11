#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(dcanr)
library(igraph)
library(networkD3)

option_list <- list(
    make_option(c("-c", "--input.file.1"), type = 'character', 
                help='Input file'),
    make_option(c("-d", "--input.file.2"), type = 'character',
                help='Input file'),
    make_option(c("-o", "--output.file"), type = 'character',
                help='output file')
)

opt <- parse_args(OptionParser(option_list=option_list))

opt$input.file.1 <- "/home/nicolai/Documents/Arbeit/NetMap/new/boostdiff-nextflow/results/Arm_vs_Doc_D28:Spleen/out_Arm_Spleen_d28.tsv"
opt$input.file.2 <- "/home/nicolai/Documents/Arbeit/NetMap/new/boostdiff-nextflow/results/Arm_vs_Doc_D28:Spleen/out_Doc_Spleen_d28.tsv"

data.1 <- read.table(file = opt$input.file.1, sep = "\t", header = TRUE)
data.2 <- read.table(file = opt$input.file.2, sep = "\t", header = TRUE)

conditions <- c(rep(1, times=ncol(data.1)-1), rep(2, times=ncol(data.2)-1)) # -1 because of GENE column that has to be ignored
# merging data
data <- cbind(data.1, data.2[, 2:ncol(data.2)])
rownames(data) <- data[, 1]
data <- data[, 2:ncol(data)]
names(conditions) <- (colnames(data))

# running dcanr pipeline
z_scores <- dcScore(data, conditions, dc.method = 'zscore', cor.method = 'spearman')
raw_p <- dcTest(z_scores, data, conditions)
adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')
dcnet <- dcNetwork(z_scores, adj_p)
edgedf <- as_data_frame(dcnet, what = 'edges')

sorted_indices <- sort(abs(edgedf$score), decreasing=TRUE, index.return=TRUE) # sort based on highest magnitude -> outliers in distribution 
sorted_edgedf_by_score <- edgedf[match(sorted_indices$ix, rownames(edgedf)),]
write.table(sorted_edgedf_by_score, file = opt$output.file, sep='\t')

