


# Comment les récupérer, sans perdre la vitesse
#
# Étape 1 — un seul fit geeglm supplémentaire, une fois à la convergence (pas à chaque itération) pour obtenir la variance sandwich de ρ1, ..., ρMv directement depuis geepack :
#
#   r
# fit_band_final <- geeglm(formula, id = id_v, waves = waves_v, data = data, family = family,
#                          corstr = "userdefined",
#                          zcor = build_banded_toeplitz_zcor(id_v, waves_v, k = Mv, maxwave = maxwave))
#
# rho_hat  <- fit_band_final$geese$alpha
# V_rho    <- fit_band_final$geese$valpha   # matrice de covariance sandwich de rho (nom de champ à vérifier)
# se_rho   <- sqrt(diag(V_rho))
# wald_rho <- rho_hat / se_rho
# p_rho    <- 2 * (1 - pnorm(abs(wald_rho)))
#
# Ça ne coûte qu'un seul appel geeglm de plus dans tout l'algorithme (pas un par itération), donc ça préserve l'essentiel du gain de vitesse du dernier correctif.
#
# Étape 2 — propager l'incertitude à φ par la méthode delta. φ = R(ρ)⁻¹ρ est une fonction non linéaire de ρ (puisque R(ρ) elle-même dépend de ρ) — il faut donc le Jacobien ∂φ/∂ρ, calculé numériquement, puis appliquer la méthode delta :
#
#   r
# library(numDeriv)
#
phi_of_alpha <- function(alpha) {
  Mv <- length(alpha)
  R_mat <- if (Mv == 1) matrix(1, 1, 1) else toeplitz(c(1, alpha[1:(Mv - 1)]))
  solve(R_mat, alpha)
}
#
# J <- jacobian(phi_of_rho, rho_hat)
# V_phi <- J %*% V_rho %*% t(J)
# se_phi <- sqrt(diag(V_phi))
# wald_phi <- phi_hat / se_phi
# p_phi <- 2 * (1 - pnorm(abs(wald_phi)))
#
# Étape 3 — stocker tout ça dans l'objet retourné, pour que summary.geeglm_arm() puisse l'afficher :
#
#   r
# fit$ar_rho_se    <- se_rho
# fit$ar_rho_wald  <- wald_rho
# fit$ar_rho_pval  <- p_rho
# fit$ar_phi_se    <- se_phi
# fit$ar_phi_wald  <- wald_phi
# fit$ar_phi_pval  <- p_phi
#
# Et dans summary.geeglm_arm(), présentez un vrai tableau, dans le même esprit que printCoefmat() :
#
#   r
# summary.geeglm_arm <- function(object, ...) {
#   s <- NextMethod()
#   s$ar_table <- data.frame(
#     Estimate  = object$ar_rho_by_lag_used,  # rho_1..rho_Mv seulement (pas les lags extrapolés)
#     Std.Error = object$ar_rho_se,
#     Wald      = object$ar_rho_wald,
#     `Pr(>|W|)` = object$ar_rho_pval,
#     row.names = paste0("lag", seq_along(object$ar_rho_se)),
#     check.names = FALSE
#   )
#   class(s) <- c("summary.geeglm_arm", class(s))
#   s
# }
# Une limite importante à documenter
#
# Ce tableau de significativité ne concerne que ρ1, ..., ρMv (et par extension φ) — pas les corrélations extrapolées au-delà de Mv (ar_correlation_by_lag[Mv+1:]), puisque celles-ci sont des fonctions déterministes de φ sans incertitude propre distincte (leur incertitude est déjà entièrement capturée par celle de φ, mais leur propre erreur-type demanderait une seconde application de la méthode delta si vous vouliez l'afficher explicitement — faisable, mais pas indispensable pour l'inférence sur φ lui-même).





# rho_full_from_rho_hat <- function(rho_hat, maxwave) {
#   phi <- phi_of_rho(rho_hat)
#   extend_ar_corr(phi, rho_hat, maxwave)
# }
alpha_all_from_alpha <- function(alpha, maxwave) {
  phi <- phi_of_alpha(alpha)
  extend_ar_corr(phi, alpha, maxwave)
}
#
# J_full <- jacobian(function(r) rho_full_from_rho_hat(r, maxwave), rho_hat)  # (maxwave-1) x Mv
# V_rho_full <- J_full %*% V_rho %*% t(J_full)   # V_rho vient de fit_band_final$geese$valpha
#
# se_rho_full   <- sqrt(diag(V_rho_full))
# wald_rho_full <- rho_all / se_rho_full          # rho_all = ar_correlation_by_lag déjà calculé
# p_rho_full    <- 2 * (1 - pnorm(abs(wald_rho_full)))





d_phi_d_alpha <- function(alpha, phi) {
  Mv <- length(alpha)
  R_mat <- if (Mv == 1) matrix(1, 1, 1) else toeplitz(c(1, alpha[1:(Mv - 1)]))
  R_inv <- solve(R_mat)
  J <- matrix(0, Mv, Mv)
  for (k in seq_len(Mv)) {
    e_k <- numeric(Mv); e_k[k] <- 1
    if (k < Mv) {
      Ek <- matrix(0, Mv, Mv)
      for (i in 1:Mv) for (j in 1:Mv) if (abs(i - j) == k) Ek[i, j] <- 1
      J[, k] <- R_inv %*% (e_k - Ek %*% phi)
    } else {
      J[, k] <- R_inv %*% e_k   # dR/drho_Mv = 0
    }
  }
  J
}




extend_ar_corr <- function(phi, alpha, maxwave) {
  Mv <- length(phi)
  alpha_ext <- c(1, alpha, rep(NA_real_, maxwave - 1 - Mv))  # alpha_ext[i] = lag (i-1)
  if (maxwave - 1 > Mv) {
    for (lag in (Mv + 1):(maxwave - 1)) {
      idx <- lag + 1
      alpha_ext[idx] <- sum(phi * alpha_ext[(idx - 1):(idx - Mv)])
    }
  }
  alpha_ext[-1]  # lag 1..(maxwave-1)
}





extend_ar_corr_jacobian <- function(alpha, maxwave) {
  Mv <- length(alpha)
  phi <- phi_of_alpha(alpha)
  Jphi <- d_phi_d_alpha(alpha, phi)     # Mv x Mv

  n_lags <- maxwave - 1
  rho_ext <- c(alpha, rep(NA_real_, n_lags - Mv))
  J <- matrix(0, n_lags, Mv)
  J[1:Mv, ] <- diag(Mv)              # rho_ext[m] = rho_m pour m <= Mv -> dérivée = identité

  if (n_lags > Mv) {
    for (m in (Mv + 1):n_lags) {
      idx_prev <- (m - 1):(m - Mv)
      rho_ext[m] <- sum(phi * rho_ext[idx_prev])

      term1 <- as.numeric(rho_ext[idx_prev] %*% Jphi)
      term2 <- as.numeric(phi %*% J[idx_prev, , drop = FALSE])
      J[m, ] <- term1 + term2
    }
  }
  list(rho_ext = rho_ext, J = J)
}





pearson_alpha_by_lag <- function(id, waves, resid_pearson, k, p) {
  id_f <- factor(id, levels = unique(id))
  resid_by_cluster <- split(resid_pearson, id_f)
  wave_by_cluster  <- split(waves, id_f)

  sums <- numeric(k); counts <- numeric(k)
  for (ci in seq_along(resid_by_cluster)) {
    e <- resid_by_cluster[[ci]]; w <- wave_by_cluster[[ci]]
    n <- length(e)
    if (n < 2) next
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      lag <- abs(w[i] - w[j])
      if (lag >= 1 && lag <= k) {
        sums[lag] <- sums[lag] + e[i] * e[j]
        counts[lag] <- counts[lag] + 1
      }
    }
  }
  phi_hat <- sum(resid_pearson^2) / (length(resid_pearson) - p)
  (sums / (counts - p)) / phi_hat
}







#' Fit GEE Models with an Autoregressive Working Correlation of Order Mv > 1
#'
#' `geeglm_arm()` fits a GEE model with an autoregressive working correlation
#' structure of order `Mv` (`"AR-M"` in the **gee** package terminology),
#' generalizing [geepack::geeglm()]'s native `corstr = "ar1"` (which only
#' covers `Mv = 1`) to any order.
#'
#' Unlike every `build_*_zcor()`-based structure in this package, the AR-M
#' correlation is a non-linear function of its underlying parameters and cannot
#' be expressed as a fixed linear combination of indicator columns, so it is fit
#' through an iterative procedure rather than a single call to
#' `geeglm(..., corstr = "userdefined")`.
#'
#' @details
#' # Algorithm
#'
#' The fitting procedure alternates, each iteration, between:
#'
#' 1. **Moment-based lag correlations.** Given the current mean-model fit,
#'    Pearson residuals are used to estimate `rho_1, ..., rho_Mv` — the
#'    correlation at each lag from 1 to `Mv` — exactly as in the toeplitz
#'    correlation structure.
#'
#' 2. **Yule-Walker.** The `Mv` autoregressive coefficients `phi` are
#'    recovered from `rho_1, ..., rho_Mv` by solving the Yule-Walker
#'    equations `R(rho) %*% phi = rho`, where `R(rho)` is the `Mv x Mv`
#'    Toeplitz matrix built from `rho_1, ..., rho_{Mv-1}`.
#'
#' 3. **Extrapolation.** Correlations beyond lag `Mv` are extended via the
#'    autoregressive recursion `rho_k = sum(phi * rho_{k-1}, ..., rho_{k-Mv})`
#'    — unlike the banded structures in this package, AR-M correlation does
#'    **not** vanish beyond lag `Mv`; it continues to decay according to
#'    this recursion.
#'
#' 4. **Mean model update.** A single Fisher-scoring step updates the
#'    regression coefficients under the resulting fixed correlation matrix
#'    (`corstr = "fixed"`), mirroring the alternating (not nested) update
#'    used by the classical `gee::gee(corstr = "AR-M")` implementation.
#'
#' # Inference
#'
#' After convergence, two additional model fits are used purely to obtain
#' correctly-aligned point estimates and sandwich standard errors (not
#' repeated at every iteration, to keep the cost negligible):
#'
#' - a `corstr = "userdefined"` fit with `zcor` from
#'   `build_banded_toeplitz_zcor(k = Mv)`, supplying `rho_1, ..., rho_Mv`
#'   and their sandwich covariance matrix together, from the same fit;
#'
#' - a final `corstr = "fixed"` fit, run to full convergence (not a single
#'   step) from the already-converged coefficients, supplying `eta` and
#'   its sandwich covariance matrix.
#'
#' Standard errors for `phi` and for the extrapolated correlations beyond
#' lag `Mv` are obtained from the sandwich covariance of `rho_1, ...,
#' rho_Mv` by the delta method, using the closed-form Jacobians of `phi =
#' R(rho)^-1 rho` and of the extrapolation recursion (no numerical
#' differentiation).
#'
#' # Limitations
#'
#' - No stationarity check is performed on the estimated `phi`; the
#'   implied correlations may in principle exceed 1 in absolute value for
#'   pathological data. Consider checking `polyroot(c(1, -phi))` after
#'   fitting.
#' - Standard errors for the extrapolated correlations tend to grow with
#'   the lag, reflecting accumulated uncertainty through the recursion —
#'   this is expected, not a symptom of a fitting problem.
#'
#' @param formula a two-sided formula, as in [geepack::geeglm()].
#' @param id cluster identifier variable.
#' @param waves variable giving the true time/position of each observation
#'   within its cluster.
#' @param data a data frame, sorted by `id` (contiguous per cluster) and,
#'   within each cluster, by `waves` in increasing order.
#' @param family as in [geepack::geeglm()].
#' @param Mv order of the autoregressive structure (`Mv >= 1`; `Mv = 1` is
#'   equivalent to `corstr = "ar1"`, though estimated by this function's own
#'   procedure rather than natively).
#' @param max_iter maximum number of outer iterations. Defaults to 50.
#' @param tol convergence tolerance on the implied correlation vector
#'   between successive iterations. Defaults to `1e-6`.
#'
#' @return
#' An object of class `c("geefit", "geeglm"`. In addition:
#'
#' * `ar_coefficients`: the estimated autoregressive coefficients `phi_1, ..., phi_Mv`.
#'
#' * `ar_correlation_by_lag`: the implied working correlation at every lag from
#'   1 to `max(waves) - 1`, including lags beyond `Mv`.
#'
#' * `ar_rho_se`, `ar_phi_se`: sandwich standard errors for
#'   `ar_correlation_by_lag` and `ar_coefficients`, respectively (delta method).
#'
#' @references
#' Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. *Biometrika*, 73(1), 13-22.
#'
#' @examples
#' \dontrun{
#' fit <- geeglm_arm(y ~ x1 + x2, id = id, waves = waves, data = mydata,
#'                    family = gaussian, Mv = 2)
#' summary(fit)
#' }
#'
#' @noRd
geeglm_arm <- function(formula, id, waves, data, family = gaussian, Mv,
                       max_iter = 50, tol = 1e-6) {

  cl <- match.call()

  # env <- new.env(parent = parent.frame())

  id_v <- eval(substitute(id), data)
  waves_v <- eval(substitute(waves), data)
  waves_v <- as.integer(factor(waves_v))
  maxwave <- max(waves_v)
  stopifnot(Mv >= 1, Mv <= maxwave - 1)

  zcor.unstr <- build_unstructured_zcor(id_v, waves_v, maxwave)
  pairs <- combn(maxwave, 2)
  lags  <- abs(pairs[2, ] - pairs[1, ])


  geeargs <- setdiff(names(formals(geefit)), "...")
  m <- match(geeargs, names(cl), 0L)
  cl_fix <- cl[c(1L, m)]
  cl_toep <- cl[c(1L, m)]

  cl_fix[[1]] <- quote(geefit)
  cl_toep[[1]] <- quote(geefit)


  cl_fix <- rlang::call_modify(cl_fix, corstr = "fixed", zcor = quote(zcor_fixed),
                               etastart = quote(eta_current),
                               control = geese.control(maxit = 1),
                               Mv = rlang::zap())
  cl_toep <- rlang::call_modify(cl_toep, corstr = "banded-toeplitz", bandwidth = Mv,
                                Mv = rlang::zap())

  alpha_all <- rep(0, maxwave - 1)  # initialisation : indépendance
  fit <- NULL
  eta_current <- NULL
  # env$eta_current <- eta_current

  for (iter in seq_len(max_iter)) {
    zcor_fixed <- as.numeric(zcor.unstr %*% alpha_all[lags])

    # fit <- geeglm(formula, id = id_v, waves = waves_v, data = data, family = family,
    #               corstr = "fixed", zcor = zcor_fixed)
    # fit <- geeglm(formula, id = id, waves = time, data = data, family = family,
    #               corstr = "fixed", zcor = zcor_fixed)
    fit <- eval(cl_fix,
                list(eta_current = eta_current,
                     zcor_fixed = zcor_fixed),
                parent.frame())
    # fit <- geefit(formula, id = id, waves = waves, data = data, family = family,
    #               corstr = "fixed", zcor = zcor_fixed, etastart = eta_current,
    #               control = geese.control(maxit = 1))

    eta_current <- predict(fit, type = "link")


    # zcor_band <- build_banded_toeplitz_zcor(id_v, waves_v,
    #                                         bandwidth = Mv, maxwave = maxwave)
    # fit_band <- geeglm(formula, id = id, waves = time, data = data, family = family,
    #                    corstr = "userdefined", zcor = zcor_band)
    # fit_band <- eval(cl_toep)
    # rho_1_to_Mv <- fit_band$geese$alpha

    resid_p <- residuals(fit, type = "pearson")   # résidus du fit "fixed" courant, MÊME beta
    alpha <- pearson_alpha_by_lag(id_v, waves_v, resid_p, k = Mv, p = length(coef(fit)))


    # R_mat <- if (Mv == 1) matrix(1, 1, 1) else toeplitz(c(1, alpha[1:(Mv - 1)]))
    # phi <- solve(R_mat, alpha)
    phi <- phi_of_alpha(alpha)

    alpha_new <- extend_ar_corr(phi, alpha, maxwave)
    if (max(abs(alpha_new - alpha_all)) < tol) { alpha_all <- alpha_new; break }
    alpha_all <- alpha_new
  }


  fit_band <- eval(cl_toep, parent.frame())
  alpha <- fit_band$geese$alpha
  V_alpha    <- fit_band$geese$valpha
  # alpha_all <- extend_ar_corr(phi, alpha, maxwave)
  res <- extend_ar_corr_jacobian(alpha, maxwave)
  alpha_all <- res$rho_ext          # remplace l'appel séparé à extend_ar_corr()
  J_full <- res$J
  names(alpha_all) <- paste("alpha", 1:(maxwave-1), sep = ":")

  # J_full <- numDeriv::jacobian(
  #   function(a) alpha_all_from_alpha(a, maxwave),
  #   alpha
  # )  # (maxwave-1) x Mv
  # V_alpha_all <- J_full %*% V_alpha %*% t(J_full)
  Vch <- chol(V_alpha, pivot = TRUE)
  Vrk <- attr(Vch, "rank")
  piv <- attr(Vch, "pivot")
  U <- Vch[1:Vrk, , drop = FALSE]
  Vh <- tcrossprod(J_full[, piv, drop = FALSE], U)
  V_alpha_all <- tcrossprod(Vh)


  fit$ar_coefficients <- phi
  fit$ar_correlation_by_lag <- alpha_all
  fit$ar_iterations <- iter
  fit$geese$alpha <- alpha_all
  fit$geese$valpha <- V_alpha_all
  fit$geese$zcor.names <- names(alpha_all)
  fit
}









geecorfit_arm <- function(x, y, id, waves, family = gaussian,
                          weights, subset, na.action, start = NULL,
                          etastart, mustart, offset, control = geese.control(...),
                          contrasts = NULL, Mv = NULL, max_iter = 50, tol = 1e-6,
                          ...)
{

  cl <- match.call()

  # env <- new.env(parent = parent.frame())

  # id_v <- eval(substitute(id), data)
  # waves_v <- eval(substitute(waves), data)
  # waves_v <- as.integer(factor(waves_v))
  waves_v <- as.integer(factor(waves))
  maxwave <- max(waves_v)
  stopifnot(Mv >= 1, Mv <= maxwave - 1)

  zcor.unstr <- build_unstructured_zcor(id, waves_v, maxwave)
  pairs <- combn(maxwave, 2)
  lags  <- abs(pairs[2, ] - pairs[1, ])


  geeargs <- setdiff(names(formals(geecor_fit)), "...")
  m <- match(geeargs, names(cl), 0L)
  cl_fix <- cl[c(1L, m)]
  cl_toep <- cl[c(1L, m)]

  cl_fix[[1]] <- quote(geecor_fit)
  cl_toep[[1]] <- quote(geecor_fit)


  cl_fix <- rlang::call_modify(cl_fix, corstr = "fixed",
                               waves = quote(waves_v),
                               zcor = quote(zcor_fixed),
                               etastart = quote(eta_current),
                               control = geese.control(maxit = 1),
                               Mv = rlang::zap())
  cl_toep <- rlang::call_modify(cl_toep, corstr = "banded-toeplitz",
                                waves = quote(waves_v),
                                bandwidth = Mv,
                                Mv = rlang::zap())

  alpha_all <- rep(0, maxwave - 1)  # initialisation : indépendance
  fit <- NULL
  eta_current <- NULL
  # env$eta_current <- eta_current
# return(list(waves_v, maxwave, cl_fix,
#             zcor.unstr,
#             zcor_fixed = as.numeric(zcor.unstr %*% alpha_all[lags])))
  for (iter in seq_len(max_iter)) {
    zcor_fixed <- as.numeric(zcor.unstr %*% alpha_all[lags])

    # fit <- geeglm(formula, id = id_v, waves = waves_v, data = data, family = family,
    #               corstr = "fixed", zcor = zcor_fixed)
    # fit <- geeglm(formula, id = id, waves = time, data = data, family = family,
    #               corstr = "fixed", zcor = zcor_fixed)
    fit <- eval(cl_fix,
                list(waves_v = waves_v,
                     eta_current = eta_current,
                     zcor_fixed = zcor_fixed),
                parent.frame())
    # fit <- geefit(formula, id = id, waves = waves, data = data, family = family,
    #               corstr = "fixed", zcor = zcor_fixed, etastart = eta_current,
    #               control = geese.control(maxit = 1))

    beta_current <- fit$beta
    eta_current <- as.numeric(x %*% beta_current) + offset
    # eta_current <- predict(fit, type = "link")


    # zcor_band <- build_banded_toeplitz_zcor(id_v, waves_v,
    #                                         bandwidth = Mv, maxwave = maxwave)
    # fit_band <- geeglm(formula, id = id, waves = time, data = data, family = family,
    #                    corstr = "userdefined", zcor = zcor_band)
    # fit_band <- eval(cl_toep)
    # rho_1_to_Mv <- fit_band$geese$alpha

    # resid_p <- residuals(fit, type = "pearson")   # résidus du fit "fixed" courant, MÊME beta
    mu <- family$linkinv(eta_current)
    resid_p <- (y - mu)*weights/sqrt(family$variance(mu))
    alpha <- pearson_alpha_by_lag(id, waves_v, resid_p, k = Mv, p = length(beta_current))


    # R_mat <- if (Mv == 1) matrix(1, 1, 1) else toeplitz(c(1, alpha[1:(Mv - 1)]))
    # phi <- solve(R_mat, alpha)
    phi <- phi_of_alpha(alpha)

    alpha_new <- extend_ar_corr(phi, alpha, maxwave)
    if (max(abs(alpha_new - alpha_all)) < tol) { alpha_all <- alpha_new; break }
    alpha_all <- alpha_new
  }


  fit_band <- eval(cl_toep, list(waves_v = waves_v), parent.frame())
  alpha <- fit_band$alpha
  V_alpha    <- fit_band$valpha
  # # alpha_all <- extend_ar_corr(phi, alpha, maxwave)
  res <- extend_ar_corr_jacobian(alpha, maxwave)
  alpha_all <- res$rho_ext          # remplace l'appel séparé à extend_ar_corr()
  J_full <- res$J
  names(alpha_all) <- paste("alpha", 1:(maxwave-1), sep = ":")
  #
  # # J_full <- numDeriv::jacobian(
  # #   function(a) alpha_all_from_alpha(a, maxwave),
  # #   alpha
  # # )  # (maxwave-1) x Mv
  # # V_alpha_all <- J_full %*% V_alpha %*% t(J_full)
  Vch <- chol(V_alpha, pivot = TRUE)
  Vrk <- attr(Vch, "rank")
  piv <- attr(Vch, "pivot")
  U <- Vch[1:Vrk, , drop = FALSE]
  Vh <- tcrossprod(J_full[, piv, drop = FALSE], U)
  V_alpha_all <- tcrossprod(Vh)
  #
  #
  # fit$ar_coefficients <- phi
  # fit$ar_correlation_by_lag <- alpha_all
  # fit$ar_iterations <- iter
  fit$alpha <- alpha_all
  fit$valpha <- V_alpha_all
  fit$zcor.names <- names(alpha_all)
  fit
}
