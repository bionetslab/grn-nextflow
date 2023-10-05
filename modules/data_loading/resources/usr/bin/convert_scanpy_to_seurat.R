library(reticulate)

option_list <- list(
  make_option(c("-f", "--file"), type='character',
              default="", help="path to scanpy file")
)

opt <- parse_args(OptionParser(option_list=option_list))

ad <- import("anndata", convert = FALSE)
dataset_ad <- ad$read_h5ad(file)
dataset_seurat <- Convert(dataset_ad, to = "seurat")

fwrite()
