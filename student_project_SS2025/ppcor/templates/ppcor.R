#!/usr/bin/env Rscript

library(ppcor)
library(tidyverse)
library(Seurat)
library(igraph)
library(dcanr)

inFile <- "${sample_rds}"
prefix <- "${prefix}"
outFile <- paste0(prefix, "_ppcor.csv")

# input expression data
inputExpr <- readRDS(inFile)
rownames(inputExpr) <- sub("_", ".", rownames(inputExpr))

inputExpr <- as.Seurat(inputExpr, data = NULL)
expr_matrix <- GetAssayData(inputExpr, slot = "counts")
expr_dense <- as.matrix(expr_matrix)
expr_df <- as.data.frame(expr_dense)
gene_names <- rownames(expr_df) 


# Run pcor using spearman's correlation as mentioned in the PNI paper 
# Link to paper: https://www.pnas.org/content/114/23/5822

pcorResults=  pcor(x= t(expr_df), method = "spearman")

# Write output to a file
# https://stackoverflow.com/questions/38664241/ranking-and-counting-matrix-elements-in-r
TF <- gene_names[c(row(pcorResults[["estimate"]]))]
target <- gene_names[c(col(pcorResults[["estimate"]]))]
importance <- c(pcorResults[["estimate"]])
pValue <-  c(pcorResults[["p.value"]])
DF = data.frame(TF, target, importance, pValue)
outDF <- DF[order(DF[["importance"]], decreasing=TRUE), ]
write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)