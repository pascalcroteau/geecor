gendat <- function() {

  id <- gl(5, 4, 20)

  visit <- rep(1:4, 5)

  y <- rnorm(20)

  data.frame(
    y,
    id,
    visit
  )[c(-2, -9), ]

}





six <- data.frame(
  case  = rep(1:16, each = 4),
  city  = rep(c("portage","kingston","kingston","portage","kingston","portage",
                "kingston","portage","portage","kingston","kingston","portage",
                "kingston","portage","kingston","portage"), each = 4),
  age   = rep(c(9, 10, 11, 12), times = 16),
  smoke = c(0,0,0,0, 1,2,2,2, 0,0,1,1, 0,0,0,1, 0,1,1,1, 0,1,1,1, 1,1,0,0,
            1,1,1,2, 2,2,1,1, 0,0,0,1, 1,0,0,0, 1,0,0,0, 1,0,1,1, 1,2,1,2,
            1,1,1,2, 1,1,2,1),
  wheeze = c(1,1,1,0, 1,1,0,0, 1,0,0,0, 0,1,1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 1,0,0,0, 0,0,0,0, 1,0,1,1, 0,0,0,0, 0,1,1,1, 0,0,0,1,
             0,0,0,1, 1,1,0,0)
)
six$city <- factor(six$city, levels = c("portage", "kingston"))





# -------- Generic simulator (known Sigma matrix) --------------------------
# Simulates n_clusters clusters of size `ncol(Sigma)`, with mean vector `mu` and
# exact correlation matrix `Sigma` (variance 1). Introduces imbalance by
# randomly removing a fraction of the observations (robustness test).

simulate_from_sigma <- function(n_clusters, Sigma, mu = 0, prop_missing = 0.15, seed = 1) {
  set.seed(seed)
  eig <- eigen(Sigma)
  if (any(eig$values <= 1e-8)) {
    stop("Sigma n'est pas définie positive - ajustez les valeurs d'alpha choisies.")
  }
  L <- chol(Sigma)  # Sigma = t(L) %*% L
  n <- ncol(Sigma)

  out <- vector("list", n_clusters)
  for (i in seq_len(n_clusters)) {
    z <- rnorm(n)
    y <- mu + as.numeric(z %*% L)
    waves <- 1:n
    if (prop_missing > 0) {
      keep <- rbinom(n, 1, 1 - prop_missing) == 1
      keep[1] <- TRUE  # garder au moins une observation
      y <- y[keep]; waves <- waves[keep]
    }
    out[[i]] <- data.frame(id = i, waves = waves, y = y)
  }
  do.call(rbind, out)
}



add_missing <- function(dat, seed = 2)
{
  dat_na <- dat
  set.seed(seed)
  na_idx <- sample(seq_len(nrow(dat_na)),
                   size = round(0.05 * nrow(dat_na)))
  dat_na$y[na_idx] <- NA

  dat_na
}




# ---------------------------------------------------------------------
# banded-toeplitz.
# ---------------------------------------------------------------------

true_alpha_band <- c(0.5, 0.3)
Sigma_band <- diag(6)
for (i in 1:6) for (j in 1:6) {
  lag <- abs(i - j)
  if (lag == 1) Sigma_band[i, j] <- 0.5
  if (lag == 2) Sigma_band[i, j] <- 0.3
}

dat_band <- simulate_from_sigma(n_clusters = 400, Sigma = Sigma_band)
dat_band_na <- add_missing(dat_band)




# ---------------------------------------------------------------------
# banded-unstructured
# ---------------------------------------------------------------------

maxwave_bu <- 5
k_bu <- 2

# true values, one per distinct pair with lag <= 2, in the order of combn()
# filtered by lag <= k -- this order must match the colnames(zcor) produced by
# build_banded_unstructured_zcor
true_pairs_bu <- list(
  "1_2" = 0.50,  # lag 1
  "1_3" = 0.30,  # lag 2
  "2_3" = 0.45,  # lag 1
  "2_4" = 0.25,  # lag 2
  "3_4" = 0.40,  # lag 1
  "3_5" = 0.20,  # lag 2
  "4_5" = 0.35   # lag 1
)

Sigma_bu <- diag(maxwave_bu)
pairs_bu <- combn(maxwave_bu, 2)
lags_bu  <- abs(pairs_bu[2, ] - pairs_bu[1, ])
for (c in which(lags_bu <= k_bu)) {
  i <- pairs_bu[1, c]; j <- pairs_bu[2, c]
  val <- true_pairs_bu[[paste0(i, "_", j)]]
  Sigma_bu[i, j] <- val
  Sigma_bu[j, i] <- val
}

dat_bu <- simulate_from_sigma(n_clusters = 2000, Sigma = Sigma_bu)
dat_bu_na <- add_missing(dat_bu)




# ---------------------------------------------------------------------
# banded-exchangeable
# ---------------------------------------------------------------------


true_alpha_bdex <- 0.4
Sigma_bdex <- diag(6)
for (i in 1:6) for (j in 1:6) if (abs(i - j) %in% c(1, 2)) Sigma_bdex[i, j] <- 0.4

dat_bdex <- simulate_from_sigma(n_clusters = 400, Sigma = Sigma_bdex)
dat_bdex_na <- add_missing(dat_bdex)



# ---------------------
# nested-exchangeable:
# ---------------------


true_alpha_nested <- c(alpha1 = 0.5, alpha2 = 0.2)
subgroup_nested <- c(1, 1, 1, 2, 2, 2)
Sigma_nested <- diag(6)
for (i in 1:6) for (j in 1:6) if (i != j) {
  Sigma_nested[i, j] <- if (subgroup_nested[i] == subgroup_nested[j]) 0.5 else 0.2
}

dat_nested <- simulate_from_sigma(n_clusters = 400, Sigma = Sigma_nested)
dat_nested <- dplyr::mutate(dat_nested,
                            sub = dplyr::if_else(waves <= 3, 1, 2))

dat_nested_na <- add_missing(dat_nested)



# ---------------------------------------------------------------------
# pairwise-grouped-exchangeable: The answer must be close to 0.5 0.4 0.3 0.2 0.15 0.25
# ---------------------------------------------------------------------


block_pge <- c(1, 1, 2, 2, 3, 3)
true_params_pge <- c(intra1 = 0.5, intra2 = 0.4, intra3 = 0.3,
                     inter12 = 0.2, inter13 = 0.15, inter23 = 0.25)
Sigma_pge <- diag(6)
for (i in 1:6) for (j in 1:6) if (i != j) {
  bi <- block_pge[i]; bj <- block_pge[j]
  Sigma_pge[i, j] <- if (bi == bj) {
    switch(bi, true_params_pge["intra1"], true_params_pge["intra2"], true_params_pge["intra3"])
  } else if (all(sort(c(bi, bj)) == c(1, 2))) true_params_pge["inter12"]
  else if (all(sort(c(bi, bj)) == c(1, 3))) true_params_pge["inter13"]
  else true_params_pge["inter23"]
}

dat_pge <- simulate_from_sigma(n_clusters = 500, Sigma = Sigma_pge)
dat_pge <- dplyr::mutate(dat_pge,
                         blk = dplyr::case_when(waves <= 2 ~ 1,
                                                waves <= 4 ~ 2,
                                                TRUE ~ 3))


dat_pge_na <- add_missing(dat_pge)

# ---------------------------------------------------------------------
# block-exchangeable: The answer must be close to 0.3 0.15 0.25
# ---------------------------------------------------------------------

sigma_a2 <- 0.15   # cluster -> contribue à alpha1, alpha2, alpha3 (= sigma_a2)
sigma_b2 <- 0.15   # période x cluster -> contribue en plus à alpha1 seulement
sigma_c2 <- 0.10   # individu (persistant) -> contribue en plus à alpha3 seulement
sigma_e2 <- 0.60   # résiduel
sigma_tot2 <- sigma_a2 + sigma_b2 + sigma_c2 + sigma_e2

true_alpha1 <- (sigma_a2 + sigma_b2) / sigma_tot2  # intra-période
true_alpha2 <- sigma_a2 / sigma_tot2                # inter-période, diff individu
true_alpha3 <- (sigma_a2 + sigma_c2) / sigma_tot2   # même individu, périodes diff

simulate_block_exch_li <- function(n_clusters, n_indiv = 3, n_periods = 3, seed = 1) {
  set.seed(seed)
  out <- vector("list", n_clusters)
  for (i in seq_len(n_clusters)) {
    a_i <- rnorm(1, 0, sqrt(sigma_a2))
    b_ik <- rnorm(n_periods, 0, sqrt(sigma_b2))
    c_ij <- rnorm(n_indiv, 0, sqrt(sigma_c2))
    rows <- list()
    for (k in 1:n_periods) for (j in 1:n_indiv) {
      e <- rnorm(1, 0, sqrt(sigma_e2))
      y <- a_i + b_ik[k] + c_ij[j] + e
      rows[[length(rows) + 1]] <- data.frame(id = i, period = k, individual = j, y = y)
    }
    out[[i]] <- do.call(rbind, rows)
  }
  do.call(rbind, out)
}

dat_li <- simulate_block_exch_li(n_clusters = 300)


dat_li_na <- add_missing(dat_li)


