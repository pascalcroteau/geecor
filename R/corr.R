

#' Convert Estimated Correlation Parameters to a Working Correlation Matrix
#'
#' Reconstructs the full working correlation matrix implied by a vector of
#' estimated `alpha` parameters, for any structure supported by
#' `geefit()`, including `"ar_m"`.
#'
#' @param corstr the working correlation structure.
#' @param alpha the estimated correlation parameter(s).
#' @param maxwave total number of possible waves. Required for every
#'   structure except `"block_exchangeable"`.
#' @param bandwidth truncation lag, for `"banded_toeplitz"`,
#'   `"banded_unstructured"`, and `"banded_exchangeable"`.
#' @param subgroup length-`maxwave` vector of subgroup membership, for
#'   `"nested_exchangeable"`.
#' @param block length-`maxwave` vector of group membership, for
#'   `"pairwise_grouped_exchangeable"`.
#' @param n_individual,n_period number of individuals and periods per
#'   cluster, for `"block_exchangeable"` (Li et al.).
#' @param Mv autoregressive order, for `"ar_m"`.
#'
#' @return a symmetric correlation matrix with 1's on the diagonal.
#'
#' @seealso [build_zcor()], [geeglm_arm()]
#' @noRd
build_corr_matrix <- function(corstr, alpha, maxwave = NULL, bandwidth = NULL,
                              subgroup = NULL, block = NULL,
                              n_individual = NULL, n_period = NULL,
                              Mv = NULL) {
  switch(corstr,
         "independence" = diag(maxwave),
         "exchangeable" = {
           stopifnot(length(alpha) == 1)
           R <- matrix(alpha, maxwave, maxwave); diag(R) <- 1; R
         },
         "ar1" = {
           stopifnot(length(alpha) == 1)
           alpha^abs(outer(1:maxwave, 1:maxwave, "-"))
         },
         "unstructured"        = build_corr_unstructured(alpha, maxwave),
         "toeplitz"            = build_corr_toeplitz(alpha, maxwave),
         "banded-toeplitz"     = build_corr_banded_toeplitz(alpha, maxwave, bandwidth),
         "banded-unstructured" = build_corr_banded_unstructured(alpha, maxwave, bandwidth),
         "banded-exchangeable" = build_corr_banded_exchangeable(alpha, maxwave, bandwidth),
         "nested-exchangeable" = build_corr_nested_exch(alpha, subgroup),
         "pairwise-grouped-exchangeable" = build_corr_pairwise_grouped(alpha, block),
         "block-exchangeable"  = build_corr_block_exch(alpha, n_individual, n_period),
         "ar-m"                = build_corr_ar_m(alpha, Mv, maxwave),
         stop("Unknown corstr: ", corstr)
  )
}

# ---- Helpers, one per structure --------------------------------------------

build_corr_unstructured <- function(alpha, maxwave) {
  pairs <- combn(maxwave, 2)
  stopifnot(length(alpha) == ncol(pairs))
  R <- diag(maxwave)
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]; j <- pairs[2, k]
    R[i, j] <- R[j, i] <- alpha[k]
  }
  R
}

build_corr_toeplitz <- function(alpha, maxwave) {
  stopifnot(length(alpha) == maxwave - 1)
  toeplitz(c(1, alpha))
}

build_corr_banded_toeplitz <- function(alpha, maxwave, bandwidth) {
  stopifnot(length(alpha) == bandwidth)
  full <- c(alpha, rep(0, maxwave - 1 - bandwidth))
  toeplitz(c(1, full))
}

build_corr_banded_unstructured <- function(alpha, maxwave, bandwidth) {
  pairs <- combn(maxwave, 2)
  lags <- abs(pairs[2, ] - pairs[1, ])
  cols <- which(lags <= bandwidth)
  stopifnot(length(alpha) == length(cols))
  R <- diag(maxwave)
  for (m in seq_along(cols)) {
    k <- cols[m]; i <- pairs[1, k]; j <- pairs[2, k]
    R[i, j] <- R[j, i] <- alpha[m]
  }
  R
}

build_corr_banded_exchangeable <- function(alpha, maxwave, bandwidth) {
  stopifnot(length(alpha) == 1)
  R <- diag(maxwave)
  for (i in 1:maxwave) for (j in 1:maxwave) {
    if (i != j && abs(i - j) <= bandwidth) R[i, j] <- alpha
  }
  R
}

build_corr_nested_exch <- function(alpha, subgroup) {
  stopifnot(length(alpha) == 2)
  maxwave <- length(subgroup)
  R <- diag(maxwave)
  for (i in 1:maxwave) for (j in 1:maxwave) if (i != j) {
    R[i, j] <- if (subgroup[i] == subgroup[j]) alpha[1] else alpha[2]
  }
  R
}

build_corr_pairwise_grouped <- function(alpha, block) {
  maxwave <- length(block)
  R <- diag(maxwave)
  for (i in 1:maxwave) for (j in 1:maxwave) if (i != j) {
    bi <- block[i]; bj <- block[j]
    nm <- if (bi == bj) paste0("intra_bloc", bi)
    else paste0("inter_", min(bi, bj), "_", max(bi, bj))
    if (!nm %in% names(alpha)) stop("Missing alpha for pair type: ", nm)
    R[i, j] <- alpha[[nm]]
  }
  R
}

build_corr_block_exch <- function(alpha, n_individual, n_period) {
  stopifnot(length(alpha) == 3)
  n <- n_individual * n_period
  period     <- rep(1:n_period, each = n_individual)
  individual <- rep(1:n_individual, times = n_period)
  R <- diag(n)
  for (i in 1:n) for (j in 1:n) if (i != j) {
    same_period <- period[i] == period[j]
    same_indiv  <- individual[i] == individual[j]
    R[i, j] <- if (same_period && !same_indiv) alpha[1]
    else if (!same_period && !same_indiv) alpha[2]
    else alpha[3]
  }
  R
}

build_corr_ar_m <- function(alpha, Mv, maxwave) {
  # stopifnot(length(alpha) == Mv)
  # phi <- phi_of_rho(alpha)                       # défini dans geeglm_arm.R
  # rho_all <- extend_ar_corr(phi, alpha, maxwave)  # défini dans geeglm_arm.R
  # toeplitz(c(1, rho_all))
  toeplitz(c(1, alpha))
}





