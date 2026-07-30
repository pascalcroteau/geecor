
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


set.seed(88)
dat <- gendat()






# ---- Simulateur générique (matrice Sigma connue) --------------------------
# Simule n_clusters clusters de taille ncol(Sigma), avec vecteur moyen mu et
# matrice de corrélation Sigma exacte (variance 1). Introduit du déséquilibre
# en retirant aléatoirement une fraction des observations (test robustesse).
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





test_that(
  "testing structures",
  {

    # ---------------------------------------------------------------------
    # Toeplitz - example from geepack
    # ---------------------------------------------------------------------

    zcor <- genZcor(clusz = table(dat$id), waves = dat$visit, corstrv=4)

    # defining the Toeplitz structure
    zcor.toep     <- matrix(NA, nrow(zcor), 3)
    zcor.toep[,1] <- apply(zcor[,c(1, 4, 6)], 1, sum)
    zcor.toep[,2] <- apply(zcor[,c(2, 5)], 1, sum)
    zcor.toep[,3] <- zcor[,3]

    fit <- geefit(y ~ 1, data = dat, id = id, waves = visit,
                  corstr = "toeplitz")

    ref <- geeglm(y ~ 1,id = id, data = dat,
                    corstr = "userdefined", zcor = zcor.toep)

    expect_identical(coef(fit), coef(ref))
    expect_identical(fit$geese$alpha, ref$geese$alpha)
    expect_identical(fit$geese$gamma, ref$geese$gamma)





    # ---------------------------------------------------------------------
    # m-dependent - example from https://www.sfu.ca/sasdoc/sashtml/stat/chap29/sect7.htm
    # The original is for exchangeable, but in this case m-dependent with
    # mdep 3 is equivalent
    # ---------------------------------------------------------------------

    #
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



    ref <- geeglm(wheeze ~ city + age + smoke, id = case, data = six,
                  waves = age, family = binomial, corstr = "exchangeable")

    fite <- geefit(wheeze ~ city + age + smoke, id = case, data = six,
                   waves = age, family = binomial, corstr = "exchangeable")

    fite_smp <- geefit(wheeze ~ city + age + smoke, id = case,
                       data = six |> dplyr::slice_sample(prop = 1),
                       waves = age, family = binomial, corstr = "exchangeable")

    fitm <- geefit(wheeze ~ city + age + smoke, id = case, data = six,
                   waves = age, family = binomial, corstr = "m-dependent",
                   mdep = 3)

    fitm_smp <- geefit(wheeze ~ city + age + smoke, id = case,
                   data = six |> dplyr::slice_sample(prop = 1), waves = age,
                   family = binomial, corstr = "m-dependent", mdep = 3)



    expect_identical(coef(fite), coef(ref))
    expect_identical(fite$geese$alpha, ref$geese$alpha)
    expect_identical(fite$geese$gamma, ref$geese$gamma)


    expect_identical(coef(fite_smp), coef(ref))
    expect_identical(fite_smp$geese$alpha, ref$geese$alpha)
    expect_identical(fite_smp$geese$gamma, ref$geese$gamma)


    expect_identical(coef(fitm), coef(ref))
    expect_identical(fitm$geese$alpha, ref$geese$alpha, ignore_attr = TRUE)
    expect_identical(fitm$geese$gamma, ref$geese$gamma)


    expect_identical(coef(fitm_smp), coef(ref))
    expect_identical(fitm_smp$geese$alpha, ref$geese$alpha, ignore_attr = TRUE)
    expect_identical(fitm_smp$geese$gamma, ref$geese$gamma)





    # ---------------------------------------------------------------------
    # banded-toeplitz. The answer must be close to 0.5 0.3
    # ---------------------------------------------------------------------


    true_alpha_band <- c(0.5, 0.3)
    Sigma_band <- diag(6)
    for (i in 1:6) for (j in 1:6) {
      lag <- abs(i - j)
      if (lag == 1) Sigma_band[i, j] <- 0.5
      if (lag == 2) Sigma_band[i, j] <- 0.3
    }

    dat_band <- simulate_from_sigma(n_clusters = 400, Sigma = Sigma_band)
    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_band, family = gaussian,
             corstr = "banded-toeplitz", waves = waves, bandwidth = 2)
      )




    # ---------------------------------------------------------------------
    # banded-toeplitz. The answer must be close to true_pairs_bu
    # ---------------------------------------------------------------------
    maxwave_bu <- 5
    k_bu <- 2

    # valeurs vraies, une par paire distincte de lag <= 2, dans l'ordre de combn()
    # filtré sur lag <= k -- cet ordre doit correspondre à colnames(zcor) produit
    # par build_banded_unstructured_zcor
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
    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_bu, family = gaussian,
             corstr = "banded-unstructured", waves = waves, bandwidth = k_bu)
    )




    # ---------------------------------------------------------------------
    # m-dependent. The answer must be close to 0.4
    # ---------------------------------------------------------------------


    true_alpha_mdep <- 0.4
    Sigma_mdep <- diag(6)
    for (i in 1:6) for (j in 1:6) if (abs(i - j) %in% c(1, 2)) Sigma_mdep[i, j] <- 0.4

    dat_mdep <- simulate_from_sigma(n_clusters = 400, Sigma = Sigma_mdep)
    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_mdep, family = gaussian,
             corstr = "m-dependent", mdep = 2, waves = waves)
    )



    # ---------------------------------------------------------------------
    # m-dependent, example from geeM.
    # We verified that the answers are close to 4 digits before taking the
    # snapshot
    # ---------------------------------------------------------------------



    resplogit_ref <- geeM::geem(resp ~ age + smoke + age:smoke, id=id,
                                data = ohio, family = binomial,
                                corstr = "m-dep" , Mv=1)
    expect_snapshot(
      geefit(resp ~ age + smoke + age:smoke, id=id,
             data = ohio, family = binomial, waves = age,
             corstr = "m-dependent" , mdep = 1)
    )



    # ---------------------------------------------------------------------
    # nested-exchangeable: The answer must be close to 0.5 0.2
    # ---------------------------------------------------------------------


    true_alpha_nested <- c(alpha1 = 0.5, alpha2 = 0.2)
    subgroup_nested <- c(1, 1, 1, 2, 2, 2)
    Sigma_nested <- diag(6)
    for (i in 1:6) for (j in 1:6) if (i != j) {
      Sigma_nested[i, j] <- if (subgroup_nested[i] == subgroup_nested[j]) 0.5 else 0.2
    }

    dat_nested <- simulate_from_sigma(n_clusters = 400, Sigma = Sigma_nested)
    dat_nested <- dplyr::mutate(dat_nested,
                                sub = dplyr::if_else(waves <= 3, 1, 2))
    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_nested, family = gaussian,
             corstr = "nested-exchangeable", waves = waves, subgroup = sub)
    )



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

    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_pge, family = gaussian,
             corstr = "pairwise-grouped-exchangeable", waves = waves,
             block = blk)
    )



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

    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_li, family = gaussian,
             corstr = "block-exchangeable", waves = period,
             individual = individual)
    )



  })
