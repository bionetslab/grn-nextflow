#!/usr/bin/env Rscript

library(sceasy)
library(reticulate)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input.file"), type="character",
              help="Input anndata file in h5ad format")
)

opt <- parse_args(OptionParser(option_list=option_list))

sceasy::convertFormat(opt$input.file, from="anndata", to="seurat", outFile='seurat_file.rds')
