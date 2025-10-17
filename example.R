# Example usage of Meff_based_FDR

source("Meff_Based_FDR.R")

set.seed(42)
n <- 200; m <- 6
X <- matrix(rnorm(n*m), n, m)
colnames(X) <- paste0("domain_", 1:m)

# fake p-values (e.g., from m regressions)
pi <- runif(m)
names(pi) <- colnames(X)

res <- Meff_based_FDR(X, pi)
print(res)
cat("\nEstimated Meff:", attr(res, "Meff"), "\n")
