# Meff-based multiple testing (Li & Ji)

This repo provides a single R function `Meff_based_FDR()` to compute multiple testing thresholds that account for correlation among tests using the Li & Ji (2005) effective number of tests (Meff). It reports both standard (S) and Meff-adjusted (E) versions of Bonferroni and FDR procedures and returns significance calls.

- Reference: Li, J., & Ji, L. (2005). Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix. Heredity, 95(3), 221–227. https://doi.org/10.1038/sj.hdy.6800717

## Installation / usage

You can use the function by sourcing the script:

```r
source("Meff_Based_FDR.R")
```

Then call `Meff_based_FDR()`:

```r
# n subjects, m tests
set.seed(123)
n <- 100; m <- 5
X <- matrix(rnorm(n*m), n, m)
colnames(X) <- paste0("test", 1:m)

# p-values for the m tests (vector or data.frame)
pi <- runif(m)
names(pi) <- colnames(X)

res <- Meff_based_FDR(X, pi, alpha = 0.05,
                      method = c("Bonferroni_S", "Bonferroni_E",
                                 "FDR_S", "FDR_E"))
print(res)
attr(res, "Meff")  # the estimated effective number of tests
```

## Function documentation

`Meff_based_FDR(data, pi, alpha = 0.05, method = c("Bonferroni_S", "Bonferroni_E", "FDR_S", "FDR_E"))`

Arguments:
- `data`: An n x m numeric matrix or data.frame of the correlated variables used across the m tests (e.g., m cognitive domain scores across n subjects).
- `pi`: P-values for the m tests. Either a numeric vector of length m (optionally named), or a data.frame with p-values in the second column (compatible with the original code). Length/rows must match `ncol(data)`.
- `alpha`: Overall error rate (default 0.05).
- `method`: One or more of `c("Bonferroni_S", "Bonferroni_E", "FDR_S", "FDR_E")`. Selects which threshold columns are included in the result. Defaults to all.

Returns: A data.frame with m rows (sorted by p ascending) containing:
- The p-values (with any labels retained)
- One column per requested threshold method
- Logical columns `Sig_*` indicating significance for each method

Notes:
- `_S` indicates the standard procedure; `_E` indicates the Meff-adjusted version.
- Missing values in `data` are handled via pairwise complete observations when forming the correlation matrix.
- Columns in `data` must have non-zero variance; otherwise Meff cannot be computed reliably.
- P-values must be in [0, 1]. Non-finite or out-of-range values will error.

Helper:
- `Meff_from_data(data)` returns just the estimated Meff (clamped to [1, m]).

## Interpreting the output

- The table is sorted from smallest to largest p-value. For FDR-style procedures, significance calls will be a run of TRUEs up to a cutoff index for each method.
- The attributes `Meff`, `alpha`, and `methods` are attached to the returned data.frame for convenience.

## Minimal example with regression p-values

Suppose you fit m separate regressions (one per cognitive domain score), and collect a p-value per model. Put the m domain scores in an `n x m` matrix `data`, and the m p-values in a vector `pi` with the same order as the columns of `data`.

```r
source("Meff_Based_FDR.R")

# data: n x m matrix of domain scores
data <- as.matrix(your_domain_scores_matrix)  # n x m

# pi: vector of m p-values (order aligns with columns of data)
pi <- your_p_values_vector  # length m

out <- Meff_based_FDR(data, pi)
print(out)
```

## License

MIT
