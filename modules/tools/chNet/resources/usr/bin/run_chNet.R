library(chNet)
library(data.table)
library(optparse)
library(stringr)
library(igraph)

option_list <- list(
    make_option(
      c("-c", "--input.file.1"), 
      type = 'character', 
      help='Input file'
    ),
    make_option(
      c("-d", "--input.file.2"), 
      type = 'character',
      help='Input file'
    ),
    make_option(
      c("-o", "--output.file"), 
      type = 'character',
      help='output file'
    ),
    make_option(
      c('-l', '--lambda'),
      type = 'double',
      help = 'hyperparamter for thresholding',
      default = 2.85
    ),
    make_option(
      c('--n_cpus'),
      type = 'integer',
      help = 'Number of cpus to use',
      default = 10
    ),
    make_option(
      c('--run_parallel'),
      type = 'bool',
      help = 'Run method in parallel',
      default = T
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt$n_cpus)
opt$n_cpus <- 10
s
cond_name_1 <- sub('.*out_', "", opt$input.file.1)
cond_name_2 <- sub('.*out_', "", opt$input.file.2)

cond_name_1 <- str_replace(cond_name_1, ".tsv", "")
cond_name_2 <- str_replace(cond_name_2, ".tsv", "")

data.1 <- as.matrix(fread(opt$input.file.1), rownames=1)
data.2 <- as.matrix(fread(opt$input.file.2), rownames=1)

groups <- c(rep(cond_name_1, ncol(data.1)), rep(cond_name_2, ncol(data.2)))
data <- t(cbind(data.1, data.2))

result <- chNet(data, groups, lambar = opt$lambda, parallel = opt$run_parallel, nCpus = opt$n_cpus)

net <- as_data_frame(result$Diff.net, 'edges')

fwrite(net, paste0(opt$output.file, "aggregated_filtered_network_zscores.txt"), sep="\t")
