#!/usr/bin/env Rscript

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input.file"), type="character",
              help="Input anndata file in h5ad format")
)

opt <- parse_args(OptionParser(option_list=option_list))

Convert(opt$input.file, dest = 'h5seurat', overwrite = TRUE)
seurat_v4 <- LoadH5Seurat(paste0(gsub('.h5ad', '', opt$input.file), '.h5seurat'))

counts <- seurat_v4[['RNA']]$counts
meta.data <- seurat_v4@meta.data
seurat_v5 <- CreateSeuratObject(counts = counts, meta.data = meta.data)

saveRDS(seurat_v5, 'seurat_object.rds')
