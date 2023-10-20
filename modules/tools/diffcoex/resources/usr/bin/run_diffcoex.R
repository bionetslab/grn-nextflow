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
                help='output file'),
    make_option(c("-n", "--top_n_edges"), type="integer", 
                default=100, help="Filter network for top_n_edges")
)

opt <- parse_args(OptionParser(option_list=option_list))

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

scores <- dcScore(data, conditions, dc.method = 'diffcoex', cor.method = 'spearman')
raw_p <- dcTest(scores, data, conditions)
adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')
dcnet <- dcNetwork(scores, adj_p)
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

    # outlier detection (|Z-score| > 2 -> outlier)
    target_mean <- mean(target_data)
    regulator_mean <- mean(regulator_data)

    target_sd <- sd(target_data)
    regulator_sd <- sd(regulator_data)

    target_mask <- !(target_data > (target_mean + 2*target_sd) | target_data < (target_mean - 2*target_sd)) 
    regulator_mask <- !(regulator_data > (regulator_mean + 2*regulator_sd) | regulator_data < (regulator_mean - 2*regulator_sd))

    mask <- regulator_mask & target_mask 

    model <- lm(formula = regulator_data ~ target_data)
    edgedf[row,]$effect <- as.numeric(model$coefficients[2])
}


fwrite(edgedf, paste0(opt$output.file, "aggregated_filtered_network_diffcoex.txt"), sep="\t")
