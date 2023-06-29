#!/usr/bin/env Rscript
renv::activate('/home/bionets-og86asub/Documents/netmap/')
require(Seurat)
library(SeuratData)
library(SeuratDisk)


input.file<-'/home/bionets-og86asub/Documents/netmap/data/misc/ga_an0228_10x_deepseq_filtered_smarta_merged_tissue_integrated_rep_timepoint_infection_filtered_seurat.rds'
adata<-readRDS(input.file)

adata<-adata$all

# Remove intermediate analysis layers
DefaultAssay(adata)<-'SCT'
adata[['RNA']]<-NULL
adata[['integrated']]<-NULL
adata[['SCTsplit']]<-NULL
adata[['SCTmerged']]<-NULL

output.file<-'/home/bionets-og86asub/Documents/netmap/data/misc/seurat_cd4_small.rds'
saveRDS(adata, file = output.file)

adata<-readRDS(output.file)
ada_sub<-subset(adata, downsample=1000)

output.file<-'/home/bionets-og86asub/Documents/netmap/data/misc/seurat_cd4_micro.rds'
ada_sub<- readRDS( file = output.file)

adata[['SCT']]<-NULL
SaveH5Seurat(adata, filename = "/home/bionets-og86asub/Documents/netmap/data/mouse.h5Seurat")
Convert("/home/bionets-og86asub/Documents/netmap/data/mouse.h5Seurat", dest = "h5ad")


VlnPlot(adata, features = 'Tox')

DimPlot(adata,)

adata[['integrated']]@data



gc()
