#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
arg <- list()
for (i in seq(1, length(args), by=2)) {
  key <- gsub("^--", "", args[i])
  val <- args[i+1]
  arg[[key]] <- val
}
edges_path <- arg[["edges"]]
out_path <- arg[["out"]]

edges <- tryCatch(read.table(edges_path, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) data.frame())
if (nrow(edges) == 0) {
  write.table(data.frame(), file=out_path, sep="\t", row.names=FALSE, quote=FALSE)
  quit(save="no")
}

genes <- sort(unique(c(edges$regulator, edges$target)))
idx <- function(v) match(v, genes)
n <- length(genes)
A <- matrix(0, n, n, dimnames=list(genes, genes))

for (i in 1:nrow(edges)) {
  r <- edges$regulator[i]; t <- edges$target[i]; s <- edges$score[i]
  A[idx(r), idx(t)] <- max(A[idx(r), idx(t)], s)
}

mx <- max(A)
if (mx > 0) A <- A / mx

out <- data.frame(regulator=character(0), target=character(0), weight=numeric(0), stringsAsFactors=FALSE)
for (i in 1:n) for (j in 1:n) {
  if (A[i,j] > 0) out <- rbind(out, data.frame(regulator=genes[i], target=genes[j], weight=A[i,j], stringsAsFactors=FALSE))
}
out <- out[order(-out$weight), ]
write.table(out, file=out_path, sep="\t", row.names=FALSE, quote=FALSE)
cat("[OK] Wrote", out_path, "\n")
