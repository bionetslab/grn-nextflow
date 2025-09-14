args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
repnum <- as.numeric(args[2])
output_file <- args[3]

if (is.na(output_file) || output_file == "") {
    stop("Missing output filename for meanA.txt")
}

# Read first A matrix from subdirectory structure
first_path <- file.path(dir, "out_1", "A.txt")
if (!file.exists(first_path)) {
    stop(paste("File not found:", first_path))
}
meanA <- as.matrix(read.table(first_path))

# Accumulate the rest
for(i in 2:repnum) {
    this_path <- file.path(dir, paste0("out_", i), "A.txt")
    if (!file.exists(this_path)) {
        stop(paste("File not found:", this_path))
    }
    tmp <- as.matrix(read.table(this_path))
    meanA <- meanA + tmp
}

# Average
meanA <- meanA / repnum

# Write to file
write.table(meanA, file=output_file, sep="\t", col.names=FALSE, row.names=FALSE)


