library(MASS)

args <- commandArgs(trailingOnly = T)
fdata <- args[1]
ftime <- args[2]
dir <- args[3]
tfnum <- as.numeric(args[4])
pnum <- as.numeric(args[5])
cnum <- as.numeric(args[6])
maxite <- as.numeric(args[7])

maxB <- 2.0
minB <- -10.0

# Ensure subdirectory exists
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

# Read and validate input data
X <- as.matrix(read.table(fdata, sep = "\t"))[1:tfnum, 1:cnum]
pseudotime <- read.table(ftime, sep = "\t")[1:cnum, 2]

# Validate data
if (any(is.na(X)) || any(is.infinite(X))) {
    stop("Expression data contains NA or infinite values")
}
if (any(is.na(pseudotime)) || any(is.infinite(pseudotime))) {
    stop("Pseudotime data contains NA or infinite values")
}

# Normalize pseudotime
pseudotime <- pseudotime / max(pseudotime)

# Initialize matrices
W <- matrix(0, tfnum, pnum)
Z <- matrix(0, pnum, cnum)
WZ <- matrix(0, tfnum, cnum)

# Initialize B values
new_B <- runif(pnum, min = minB, max = maxB)
old_B <- new_B
RSS <- sum(X^2) # Initialize with sum of squares of X

sample_Z <- function() {
    for (i in 1:pnum) {
        for (j in 1:cnum) {
            Z[i, j] <<- exp(new_B[i] * pseudotime[j]) + runif(1, min = -0.001, max = 0.001)
        }
    }
}

for (ite in 1:maxite) {
    target <- sample(1:pnum, 1)
    new_B[target] <- runif(1, min = minB, max = maxB)

    if (ite == maxite) {
        new_B <- old_B
    }

    sample_Z()

    # Fit linear models with error handling
    for (i in 1:tfnum) {
        tryCatch(
            {
                X.lm <- lm(X[i, ] ~ t(Z) - 1)
                coeffs <- X.lm$coefficients

                # Handle NA coefficients
                if (any(is.na(coeffs))) {
                    # Use previous W values or zeros if first iteration
                    if (ite == 1) {
                        W[i, ] <- rep(0, pnum)
                    }
                    # else keep previous W[i,] values
                } else {
                    W[i, ] <- coeffs
                }
            },
            error = function(e) {
                # On error, keep previous W values or use zeros
                if (ite == 1) {
                    W[i, ] <- rep(0, pnum)
                }
            }
        )
    }

    # Calculate WZ
    WZ <- W %*% Z

    # Calculate RSS with error handling
    tmp_RSS <- tryCatch(
        {
            sum((X - WZ)^2)
        },
        error = function(e) {
            return(Inf)
        }
    )

    # Handle NA/NaN values in RSS
    if (is.na(tmp_RSS) || is.nan(tmp_RSS) || is.infinite(tmp_RSS)) {
        tmp_RSS <- Inf
    }

    # Update if improvement
    if (tmp_RSS < RSS && is.finite(tmp_RSS)) {
        RSS <- tmp_RSS
        old_B[target] <- new_B[target]
    } else {
        new_B[target] <- old_B[target]
    }
}

# Final calculations
B <- diag(new_B)
invW <- tryCatch(
    {
        ginv(W)
    },
    error = function(e) {
        # If ginv fails, use identity matrix
        diag(ncol(W))
    }
)

A <- W %*% B %*% invW

# Ensure output directory exists and write results
write.table(RSS, file.path(dir, "RSS.txt"), row.names = F, col.names = F, sep = "\t")
write.table(W, file.path(dir, "W.txt"), row.names = F, col.names = F, sep = "\t")
write.table(B, file.path(dir, "B.txt"), row.names = F, col.names = F, sep = "\t")
write.table(A, file.path(dir, "A.txt"), row.names = F, col.names = F, sep = "\t")
