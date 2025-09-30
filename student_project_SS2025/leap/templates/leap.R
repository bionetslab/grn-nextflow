#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(igraph)
library(dcanr)

# Since LEAP is an archived CRAN package it is not available via conda and has to be installes sepreately
if (!require("LEAP", quietly = TRUE)) {
  install.packages('https://cran.r-project.org/src/contrib/Archive/LEAP/LEAP_0.2.tar.gz', 
                   type = 'source', repos = NULL, quiet = TRUE)
}

library(LEAP)

inFile <- "${sample_rds}"
prefix <- "${prefix}"
maxLag <- 0.33 # default maxLag value form Beeline implementation
outFile <- paste0(prefix, "_leap.csv")


# input expression data
inputExpr <- readRDS(inFile)
rownames(inputExpr) <- sub("_", ".", rownames(inputExpr))

inputExpr <- as.Seurat(inputExpr, data = NULL)
expr_matrix <- GetAssayData(inputExpr, layer = "counts")
expr_dense <- as.matrix(expr_matrix)
expr_df <- as.data.frame(expr_dense)
gene_names <- rownames(expr_df) 
# Run LEAP's compute Max. Absolute Correlation
# MAC_cutoff is set to zero to get a score for all TFs
# max_lag_prop is set to the max. recommended value from the paper's supplementary file
# Link to paper: https://academic.oup.com/bioinformatics/article/33/5/764/2557687

MAC_results = MAC_counter(data = expr_df, max_lag_prop=maxLag, MAC_cutoff = 0, 
                          file_name = "temp", lag_matrix = FALSE, symmetric = FALSE)

# Write output to a file
TF <- gene_names[MAC_results[,'Row gene index']]
target <- gene_names[MAC_results[,'Column gene index']]
importance <- MAC_results[,'Correlation']
outDF <- data.frame(TF, target, importance)
write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)