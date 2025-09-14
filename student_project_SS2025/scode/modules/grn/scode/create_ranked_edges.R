args <- commandArgs(trailingOnly = TRUE)
mean_a_file <- args[1]
output_file <- args[2]

# Read the mean A matrix
A <- as.matrix(read.table(mean_a_file))

# Get dimensions
n_genes <- nrow(A)

# Create edge list
edges <- data.frame(
  source = integer(),
  target = integer(),
  weight = numeric(),
  stringsAsFactors = FALSE
)

# Extract edges (excluding diagonal)
for (i in 1:n_genes) {
  for (j in 1:n_genes) {
    if (i != j) {
      edges <- rbind(edges, data.frame(
        source = i,
        target = j,
        weight = A[i, j]
      ))
    }
  }
}

# Rank edges by absolute weight (strongest interactions first)
edges$abs_weight <- abs(edges$weight)
edges <- edges[order(edges$abs_weight, decreasing = TRUE), ]

# Remove the temporary abs_weight column
edges$abs_weight <- NULL

# Add rank column
edges$rank <- 1:nrow(edges)

# Reorder columns
edges <- edges[, c("rank", "source", "target", "weight")]

# Write to CSV
write.csv(edges, file = output_file, row.names = FALSE)
