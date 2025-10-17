#' Meff-based FDR and Bonferroni thresholds
#'
#' Estimate the effective number of independent tests (Meff; Li & Ji, 2005)
#' from the correlation among m tests, then compute multiple-comparison
#' thresholds using both standard (S) and effective (E) versions of
#' Bonferroni and FDR procedures. Returns sorted p-values, thresholds, and
#' boolean significance calls.
#'
#' @references
#' Li, J., & Ji, L. (2005). Adjusting multiple testing in multilocus analyses
#' using the eigenvalues of a correlation matrix. Heredity, 95(3), 221–227.
#' <https://doi.org/10.1038/sj.hdy.6800717>
#'
#' @param data An n x m numeric matrix or data.frame of the correlated
#'   variables used across the m tests (e.g., m cognitive scores across n
#'   subjects). Columns must have non-zero variance.
#' @param pi P-values for the m tests. Either a numeric vector of length m
#'   (optionally named), or a data.frame with p-values in the second column
#'   (compatible with legacy code). Length/rows must match ncol(data).
#' @param alpha Overall error rate (default 0.05).
#' @param method Character vector: one or more of
#'   c("Bonferroni_S", "Bonferroni_E", "FDR_S", "FDR_E"). Selects which
#'   threshold columns are included in the result. Defaults to all.
#'
#' @return A data.frame with m rows (sorted by p ascending) containing:
#'   - p-values and optional labels
#'   - one column per requested threshold method
#'   - logical columns `Sig_*` indicating significance for each method
#'   Attributes `Meff`, `alpha`, and `methods` are attached.
#'
#' @export
#' @examples
#' set.seed(123)
#' n <- 100; m <- 5
#' X <- matrix(rnorm(n*m), n, m)
#' colnames(X) <- paste0("test", 1:m)
#' pvals <- stats::runif(m); names(pvals) <- colnames(X)
#' Meff_based_FDR(X, pvals)

Meff_based_FDR <- function(data, pi, alpha = 0.05,
                           method = c("Bonferroni_S", "Bonferroni_E",
                                       "FDR_S", "FDR_E")) {
  # Validate inputs
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.numeric(data)) stop("'data' must be numeric matrix/data.frame")
  if (length(dim(data)) != 2) stop("'data' must be 2D (n x m)")
  M <- ncol(data)
  if (M < 1) stop("'data' must have at least one column (m >= 1)")

  # Check alpha
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number in (0, 1)")
  }

  # Check columns have non-zero variance (or stop with a helpful message)
  sds <- apply(data, 2, stats::sd, na.rm = TRUE)
  if (any(!is.finite(sds) | sds == 0)) {
    bad <- which(!is.finite(sds) | sds == 0)
    stop(sprintf("One or more columns in 'data' have zero or undefined variance (columns: %s). Remove or adjust these columns before computing Meff.",
                 paste(bad, collapse = ", ")))
  }

  # Normalize p input to a data.frame with at least two columns where 2nd is p
  # Keep optional labels if provided
  if (is.numeric(pi)) {
    if (length(pi) != M) stop("Length of 'pi' must equal ncol(data)")
    if (any(!is.finite(pi))) stop("'pi' contains non-finite values")
    if (any(pi < 0 | pi > 1)) stop("All p-values in 'pi' must be in [0, 1]")
    lab <- names(pi)
    p_df <- data.frame(label = if (!is.null(lab)) lab else seq_len(M), p = as.numeric(pi),
                       stringsAsFactors = FALSE)
  } else if (is.data.frame(pi)) {
    if (ncol(pi) < 2) stop("'pi' data.frame must have at least 2 columns with p-values in column 2")
    if (nrow(pi) != M) stop("nrow(pi) must equal ncol(data)")
    p_df <- pi
    names(p_df)[2] <- "p"
    if (any(!is.finite(p_df[["p"]]))) stop("p-values column in 'pi' contains non-finite values")
    if (any(p_df[["p"]] < 0 | p_df[["p"]] > 1)) stop("All p-values in 'pi' must be in [0, 1]")
  } else {
    stop("'pi' must be a numeric vector or a data.frame")
  }

  # Step 1: estimate effective number of independent tests (Li & Ji, Eq. 5)
  cormat <- stats::cor(data, use = "pairwise.complete.obs")
  if (any(!is.finite(cormat))) {
    stop("Correlation matrix contains non-finite values; check for columns with too few non-missing observations or zero variance.")
  }
  ei <- eigen(cormat, only.values = TRUE, symmetric = TRUE)
  v <- ei$values
  Meff_raw <- sum(as.integer(v >= 1) + (v - as.integer(v)))
  # Clamp to [1, M] to avoid numerical issues
  Meff <- max(1, min(M, Meff_raw))

  # Step 2: sort p-values ascending
  ord <- order(p_df[["p"]])
  p_df <- p_df[ord, , drop = FALSE]
  i <- seq_len(M)

  # Step 3: compute thresholds
  # Columns in fixed order; we'll subset by 'method' later
  # Guard for M == 1 in FDR_E formula to avoid division by zero
  FDR_E_vec <- if (M > 1) {
    alpha / Meff + (i - 1) / (M - 1) * (alpha - alpha / Meff)
  } else {
    rep(alpha / Meff, M)
  }
  thresh <- cbind(
    Bonferroni_S = rep(1 - (1 - alpha)^(1 / M), M),
    Bonferroni_E = rep(1 - (1 - alpha)^(1 / Meff), M),
    FDR_S        = alpha * i / M,
    FDR_E        = FDR_E_vec
  )
  # Ensure numeric precision readable
  thresh <- signif(thresh, digits = 6)

  # Step 4: significance calls for each method
  # For per-row varying thresholds (FDR_*), compare row-wise
  calls <- sapply(seq_len(ncol(thresh)), function(j) {
    tcol <- thresh[, j]
    pj <- p_df[["p"]]
    # number of p <= threshold for this method
    idx <- which(pj <= tcol)
    sig <- rep(FALSE, M)
    if (length(idx)) sig[seq_len(max(idx))] <- TRUE
    sig
  })
  colnames(calls) <- colnames(thresh)

  # Filter to requested methods (keep order as requested if valid)
  method <- match.arg(method, choices = colnames(thresh), several.ok = TRUE)
  keep_cols <- match(method, colnames(thresh))
  thresh <- thresh[, keep_cols, drop = FALSE]
  calls  <- calls[,  keep_cols, drop = FALSE]

  # Build output: p-values + thresholds + significance
  # For clarity, append "Sig_" prefix to boolean columns
  sig_df <- as.data.frame(calls)
  names(sig_df) <- paste0("Sig_", colnames(sig_df))

  out <- cbind(p_df, as.data.frame(thresh), sig_df)
  rownames(out) <- NULL
  attr(out, "Meff") <- Meff
  attr(out, "alpha") <- alpha
  attr(out, "methods") <- method
  return(out)
}

#' Compute the effective number of independent tests (Meff) from data
#'
#' Convenience wrapper to estimate Meff (Li & Ji, 2005) directly from a
#' numeric matrix/data.frame using the same validation rules as
#' `Meff_based_FDR`.
#'
#' @inheritParams Meff_based_FDR
#' @return A single numeric value: estimated Meff in [1, m].
#' @export
Meff_from_data <- function(data) {
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.numeric(data)) stop("'data' must be numeric matrix/data.frame")
  if (length(dim(data)) != 2) stop("'data' must be 2D (n x m)")
  M <- ncol(data)
  if (M < 1) stop("'data' must have at least one column (m >= 1)")
  sds <- apply(data, 2, stats::sd, na.rm = TRUE)
  if (any(!is.finite(sds) | sds == 0)) {
    bad <- which(!is.finite(sds) | sds == 0)
    stop(sprintf("One or more columns in 'data' have zero or undefined variance (columns: %s). Remove or adjust these columns before computing Meff.",
                 paste(bad, collapse = ", ")))
  }
  cormat <- stats::cor(data, use = "pairwise.complete.obs")
  if (any(!is.finite(cormat))) {
    stop("Correlation matrix contains non-finite values; check for columns with too few non-missing observations or zero variance.")
  }
  v <- eigen(cormat, only.values = TRUE, symmetric = TRUE)$values
  Meff_raw <- sum(as.integer(v >= 1) + (v - as.integer(v)))
  max(1, min(M, Meff_raw))
}
