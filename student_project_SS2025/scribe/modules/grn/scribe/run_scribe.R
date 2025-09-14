#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
arg <- list()
for (i in seq(1, length(args), by = 2)) {
  key <- gsub("^--", "", args[i])
  val <- args[i + 1]
  arg[[key]] <- val
}

expr_path <- arg[["expr"]]
genes_path <- arg[["genes"]]
pt_path <- arg[["pseudotime"]]
threads <- as.integer(arg[["threads"]])
if (is.na(threads)) threads <- 1
set.seed(as.integer(arg[["seed"]]))

edges_out <- arg[["edges"]]

# Load data
expr <- as.matrix(read.csv(expr_path, header = TRUE, check.names = FALSE))
genes <- readLines(genes_path, warn = FALSE)
pt <- read.csv(pt_path, header = TRUE)[, 1]
if (ncol(expr) != length(genes)) {
  stop("Number of columns in expr does not match genes length")
}

# Z-score normalize genes
zscore <- function(x) {
  s <- sd(x)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x)) / s
}
X <- apply(expr, 2, zscore)

# Order by pseudotime
ord <- order(pt)
X <- X[ord, , drop = FALSE]
pt <- pt[ord]

# Build Gaussian kernel weights along pseudotime
gaussian_weights <- function(t, bandwidth) {
  n <- length(t)
  W <- matrix(0, n, n)
  for (i in 1:n) {
    d <- (t - t[i]) / bandwidth
    W[i, ] <- exp(-0.5 * d * d)
  }
  W
}

n <- nrow(X)
p <- ncol(X)
bw <- sd(pt)
if (is.na(bw) || bw == 0) bw <- 1
W <- gaussian_weights(pt, bw)

lag <- 1
if (n <= lag + 5) lag <- 0

# Weighted correlation scoring (dependency-light)
score_gene_pair <- function(y, x, w) {
  wy <- y - sum(w * y) / sum(w)
  wx <- x - sum(w * x) / sum(w)
  num <- sum(w * wx * wy)
  den <- sqrt(sum(w * wx * wx) * sum(w * wy * wy))
  if (den > 0) abs(num / den) else 0
}

edges <- NULL
if (lag == 0) {
  # contemporaneous
  for (ti in 1:p) {
    y <- X[, ti]
    scores <- numeric(p)
    for (j in 1:p) {
      if (j == ti) {
        scores[j] <- 0
        next
      }
      scores[j] <- score_gene_pair(y, X[, j], diag(W))
    }
    o <- order(scores, decreasing = TRUE)
    keep <- which(scores[o] > 0)
    if (length(keep) > 0) {
      o <- o[keep]
      edges <- rbind(edges, data.frame(regulator = genes[o], target = genes[ti], score = scores[o], stringsAsFactors = FALSE))
    }
  }
} else {
  # lagged: target at t depends on regulator at t-lag
  Tn <- n - lag
  wdiag <- diag(W[(1 + lag):n, 1:(n - lag)])
  if (length(wdiag) != Tn) wdiag <- rep(1, Tn)
  for (ti in 1:p) {
    y <- X[(lag + 1):n, ti]
    scores <- numeric(p)
    for (j in 1:p) {
      if (j == ti) {
        scores[j] <- 0
        next
      }
      xj <- X[1:(n - lag), j]
      scores[j] <- score_gene_pair(y, xj, wdiag)
    }
    o <- order(scores, decreasing = TRUE)
    keep <- which(scores[o] > 0)
    if (length(keep) > 0) {
      o <- o[keep]
      edges <- rbind(edges, data.frame(regulator = genes[o], target = genes[ti], score = scores[o], stringsAsFactors = FALSE))
    }
  }
}

if (is.null(edges)) {
  edges <- data.frame(regulator = character(0), target = character(0), score = numeric(0))
}
edges <- edges[edges$regulator != edges$target, ]
edges <- edges[order(-edges$score), ]
write.table(edges, file = edges_out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("[OK] Wrote", edges_out, "\n")
