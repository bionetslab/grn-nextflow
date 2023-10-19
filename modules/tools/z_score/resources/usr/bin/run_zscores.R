#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(dcanr)
library(igraph)
library(stringr)

option_list <- list(
    make_option(c("-c", "--input.file.1"), type = 'character', 
                help='Input file'),
    make_option(c("-d", "--input.file.2"), type = 'character',
                help='Input file'),
    make_option(c("-o", "--output.file"), type = 'character',
                default="result_zscores.tsv", help='output file'),
    make_option(c("-n", "--top_n_edges"), type="integer", 
                default=100, help="Filter network for top_n_edges"),
    make_option(c("-a", "--alpha"), type='integer',
                default=95, help="remove (1-alpha)% values")
)

opt <- parse_args(OptionParser(option_list=option_list))

opt$input.file.1 <- "/home/nicolai/Documents/Arbeit/InterNet_Xplorer/src/grn-nextflow/results/Arm_vs_Doc_D10:Liver/out_Arm_Liver_d10.tsv"
opt$input.file.2 <- "/home/nicolai/Documents/Arbeit/InterNet_Xplorer/src/grn-nextflow/results/Arm_vs_Doc_D10:Liver/out_Doc_Liver_d10.tsv"

cond_name_1 <- sub('.*out_', "", opt$input.file.1)
cond_name_2 <- sub('.*out_', "", opt$input.file.2)

cond_name_1 <- str_replace(cond_name_1, ".tsv", "")
cond_name_1 <- str_replace(cond_name_1, ".tsv", "")

data.1 <- read.table(file = opt$input.file.1, sep = "\t", header = TRUE)
data.2 <- read.table(file = opt$input.file.2, sep = "\t", header = TRUE)
# print(rownames(data.1))
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

edgedf$abs_score <- abs(edgedf$score)
edgedf$condition <- ifelse(edgedf$score >= 0, cond_name_1, cond_name_2)
edgedf <- setorder(edgedf, -abs_score)
edgedf <- edgedf[1:opt$top_n_edges,]

edgedf <- edgedf[, !(names(edgedf) %in% c("abs_score", "color"))]
colnames(edgedf) <- c("target", "regulator", "weight", "condition")
edgedf$effect <- 0

for(row in 1:nrow(edgedf)) {
    target <- edgedf[row, 1]
    regulator <- edgedf[row, 2]
    condition <- edgedf[row, 4]

    if (condition == cond_name_1) {
        target_data <- as.numeric(data.1[data.1$Gene == target,2:ncol(data.1)]) 
        regulator_data <- as.numeric(data.1[data.1$Gene == regulator,2:ncol(data.1)])        
    } else {
        target_data <- as.numeric(data.2[data.2$Gene == target,2:ncol(data.2)]) 
        regulator_data <- as.numeric(data.2[data.2$Gene == regulator,2:ncol(data.2)])        
    }
    # alpha capping
    regulator_threshold <- quantile(regulator_data, probs=0.95)
    target_threshold <- quantile(target_data, probs=0.95)
    regulator_mask <- regulator_data < regulator_threshold
    target_mask <- target_data < target_threshold
    mask <- regulator_mask & target_mask

    regulator_data <- regulator_data[mask]
    target_data <- target_data[mask]  
    model <- lm(formula = regulator_data ~ target_data)
    edgedf[row,]$effect <- as.numeric(model$coefficients[2])
}

fwrite(edgedf, opt$output.file, sep="\t")
