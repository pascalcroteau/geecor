# trié dat par id puis waves
# pour grouped_exchangeable: par id (cluster), puis period puis subj_id


build_toeplitz_zcor <- function(id, waves, maxwave = NULL)
{
  clusz <- as.integer(table(factor(id, levels = unique(id))))
  # wvs <- as.integer(factor(waves))
  wvs <- waves
  # if (is.null(maxwave)) maxwave <- max(waves)
  if (is.null(maxwave)) maxwave <- max(wvs)

  # matrice de design non structurée (gère nativement le déséquilibre)
  # zcor.unstr <- genZcor(clusz = clusz, waves = waves, corstrv = 4)
  zcor.unstr <- genZcor(clusz = clusz, waves = wvs, corstrv = 4)

  # correspondance colonne <-> paire (s,t), s<t, dans l'ordre de combn()
  pairs <- combn(maxwave, 2)
  lags  <- abs(pairs[2, ] - pairs[1, ])

  maxlag <- maxwave - 1
  zcor.toep <- sapply(seq_len(maxlag), function(l) {
    cols <- which(lags == l)
    if (length(cols) == 1) zcor.unstr[, cols] else rowSums(zcor.unstr[, cols, drop = FALSE])
  })

  zcor.toep
}




build_banded_toeplitz_zcor <- function(id, waves, bandwidth, maxwave = NULL)
{
  # id      : identifiant du cluster/sujet
  # waves   : temps de mesure réel (pas l'ordre des lignes)
  # bandwidth       : largeur de bande (nombre de lags distincts à estimer, bandwidth >= 1)
  # maxwave : nb total de vagues possibles (déduit des données si non fourni)

  clusz <- as.integer(table(factor(id, levels = unique(id))))
  # wvs <- as.integer(factor(waves))
  wvs <- waves
  # if (is.null(maxwave)) maxwave <- max(waves)
  if (is.null(maxwave)) maxwave <- max(wvs)

  stopifnot(bandwidth >= 1, bandwidth <= maxwave - 1)

  # matrice de design non structurée (gère nativement le déséquilibre)
  # zcor.unstr <- genZcor(clusz = clusz, waves = waves, corstrv = 4)
  zcor.unstr <- genZcor(clusz = clusz, waves = wvs, corstrv = 4)

  # correspondance colonne <-> paire (s,t), s<t, ordre de combn()
  pairs <- combn(maxwave, 2)
  lags  <- abs(pairs[2, ] - pairs[1, ])

  # une colonne par lag de 1 à bandwidth ; les lags > bandwidth ne sont assignés à aucune
  # colonne -> ligne nulle -> corrélation fixée à 0 pour ces paires
  zcor.band <- sapply(seq_len(bandwidth), function(l) {
    cols <- which(lags == l)
    if (length(cols) == 1) zcor.unstr[, cols] else rowSums(zcor.unstr[, cols, drop = FALSE])
  })

  zcor.band
}




build_banded_unstructured_zcor <- function(id, waves, bandwidth, maxwave = NULL) {
  # id      : identifiant du cluster/sujet
  # waves   : temps de mesure réel (pas l'ordre des lignes)
  # bandwidth       : largeur de bande (lag maximal ayant une corrélation non nulle)
  # maxwave : nb total de vagues possibles (déduit des données si non fourni)

  clusz <- as.integer(table(factor(id, levels = unique(id))))
  # wvs <- as.integer(factor(waves))
  wvs <- waves
  # if (is.null(maxwave)) maxwave <- max(waves)
  if (is.null(maxwave)) maxwave <- max(wvs)

  stopifnot(bandwidth >= 1, bandwidth <= maxwave - 1)

  # zcor.unstr <- genZcor(clusz = clusz, waves = waves, corstrv = 4)
  zcor.unstr <- genZcor(clusz = clusz, waves = wvs, corstrv = 4)

  pairs <- combn(maxwave, 2)
  lags  <- abs(pairs[2, ] - pairs[1, ])

  # on garde chaque paire de lag <= bandwidth comme sa PROPRE colonne (pas de regroupement)
  cols <- which(lags <= bandwidth)
  zcor.band.unstr <- zcor.unstr[, cols, drop = FALSE]

  colnames(zcor.band.unstr) <- paste0("pair_", pairs[1, cols], "_", pairs[2, cols])
  zcor.band.unstr
}





build_mdep_common_zcor <- function(id, waves, m, maxwave = NULL) {
  # id      : identifiant du cluster/sujet
  # waves   : temps de mesure réel (pas l'ordre des lignes)
  # m       : ordre de dépendance ; un seul paramètre commun pour les lags 1..m
  # maxwave : nb total de vagues possibles (déduit des données si non fourni)

  clusz <- as.integer(table(factor(id, levels = unique(id))))
  # wvs <- as.integer(factor(waves))
  wvs <- waves
  # if (is.null(maxwave)) maxwave <- max(waves)
  if (is.null(maxwave)) maxwave <- max(wvs)

  stopifnot(m >= 1, m <= maxwave - 1)

  # zcor.unstr <- genZcor(clusz = clusz, waves = waves, corstrv = 4)
  zcor.unstr <- genZcor(clusz = clusz, waves = wvs, corstrv = 4)

  pairs <- combn(maxwave, 2)
  lags  <- abs(pairs[2, ] - pairs[1, ])

  # toutes les paires avec lag <= m regroupées dans UNE seule colonne
  # (un seul paramètre) ; lags > m -> exclues -> corrélation fixée à 0
  cols <- which(lags <= m)
  zcor.mdep <- matrix(rowSums(zcor.unstr[, cols, drop = FALSE]), ncol = 1)

  zcor.mdep
}





build_nested_exch_zcor <- function(id, waves, subgroup, maxwave = NULL) {
  # id       : identifiant du patient/cluster (le niveau GEE le plus haut)
  # waves    : temps/position réel de chaque observation (pas l'ordre des lignes)
  # subgroup : vecteur de longueur maxwave donnant le sous-groupe de chaque
  #            vague (ex: c("G","G","D","D") pour 2 mesures oeil gauche,
  #            2 mesures oeil droit ; ou le numéro de visite si plusieurs
  #            mesures par visite)
  # maxwave  : nb total de vagues possibles (déduit de waves si non fourni)

  clusz <- as.integer(table(factor(id, levels = unique(id))))
  # wvs <- as.integer(factor(waves))
  wvs <- waves
  # if (is.null(maxwave)) maxwave <- max(waves)
  if (is.null(maxwave)) maxwave <- max(wvs)

  stopifnot(length(subgroup) == maxwave)

  # zcor.unstr <- genZcor(clusz = clusz, waves = waves, corstrv = 4)
  zcor.unstr <- genZcor(clusz = clusz, waves = wvs, corstrv = 4)

  pairs <- combn(maxwave, 2)
  same_subgroup <- subgroup[pairs[1, ]] == subgroup[pairs[2, ]]

  col_alpha1 <- rowSums(zcor.unstr[, same_subgroup, drop = FALSE])   # même sous-groupe
  col_alpha2 <- rowSums(zcor.unstr[, !same_subgroup, drop = FALSE])  # sous-groupe différent

  cbind(alpha1 = col_alpha1, alpha2 = col_alpha2)
}






build_pairwise_grouped_exch_zcor <- function(id, waves, block, maxwave = NULL) {
  # id      : identifiant du cluster/sujet
  # waves   : temps de mesure réel (pas l'ordre des lignes)
  # block   : vecteur de longueur maxwave donnant le numéro de bloc pour
  #           chaque vague (ex: c(1,1,2,2,2,3) pour 3 blocs)
  # maxwave : nb total de vagues possibles (déduit de waves si non fourni)

  clusz <- as.integer(table(factor(id, levels = unique(id))))
  # wvs <- as.integer(factor(waves))
  wvs <- waves
  # if (is.null(maxwave)) maxwave <- max(waves)
  if (is.null(maxwave)) maxwave <- max(wvs)

  stopifnot(length(block) == maxwave)

  # zcor.unstr <- genZcor(clusz = clusz, waves = waves, corstrv = 4)
  zcor.unstr <- genZcor(clusz = clusz, waves = wvs, corstrv = 4)

  pairs <- combn(maxwave, 2)
  block_s <- block[pairs[1, ]]
  block_t <- block[pairs[2, ]]

  blocks_uniq <- sort(unique(block))
  nblocks <- length(blocks_uniq)

  # toutes les paires de blocs (b1, b2) avec b1 <= b2 : b1==b2 -> intra-bloc,
  # b1 < b2 -> inter-blocs pour cette paire de blocs précise
  block_pairs <- combn(nblocks, 2)  # paires de blocs distincts (b1 < b2)

  # colonnes intra-bloc : une par bloc
  intra_cols <- lapply(blocks_uniq, function(b) {
    cols <- which(block_s == b & block_t == b)
    if (length(cols) == 1) zcor.unstr[, cols] else rowSums(zcor.unstr[, cols, drop = FALSE])
  })
  names(intra_cols) <- paste0("intra_bloc", blocks_uniq)

  # colonnes inter-blocs : une par paire de blocs distincte (b1, b2)
  inter_cols <- lapply(seq_len(ncol(block_pairs)), function(k) {
    b1 <- blocks_uniq[block_pairs[1, k]]
    b2 <- blocks_uniq[block_pairs[2, k]]
    cols <- which((block_s == b1 & block_t == b2) | (block_s == b2 & block_t == b1))
    if (length(cols) == 1) zcor.unstr[, cols] else rowSums(zcor.unstr[, cols, drop = FALSE])
  })
  names(inter_cols) <- apply(block_pairs, 2, function(k) {
    paste0("inter_", blocks_uniq[k[1]], "_", blocks_uniq[k[2]])
  })

  zcor.block <- do.call(cbind, c(intra_cols, inter_cols))
  zcor.block
}




build_block_exch_li_zcor <- function(id, period, individual) {
  # id         : identifiant du cluster (l'unité GEE, ex: le cluster/site dans
  #              un essai en grappes)
  # period     : période/vague de chaque observation (peut se répéter au sein
  #              d'un même cluster, puisque plusieurs individus partagent une
  #              même période)
  # individual : identifiant de l'individu, unique au sein d'un cluster, qui
  #              PERSISTE à travers les périodes pour les individus suivis en
  #              cohorte (même sujet revu à chaque période), ou qui est
  #              distinct à chaque observation si le design est en coupe
  #              transversale pour cet individu

  stopifnot(length(id) == length(period), length(id) == length(individual))

  id_f <- factor(id, levels = unique(id))  # préserve l'ordre d'apparition
  clusters <- split(seq_along(id), id_f)

  zcor_list <- lapply(clusters, function(idx) {
    n <- length(idx)
    if (n < 2) return(NULL)

    pr <- period[idx]
    iv <- individual[idx]

    pairs <- combn(n, 2)
    same_period <- pr[pairs[1, ]] == pr[pairs[2, ]]
    same_indiv  <- iv[pairs[1, ]] == iv[pairs[2, ]]

    if (any(same_period & same_indiv)) {
      stop("Paire avec même période ET même individu détectée ",
           "(observation dupliquée ?) au sein d'un cluster.")
    }

    alpha1 <- as.numeric(same_period & !same_indiv)   # intra-période, individus différents
    alpha2 <- as.numeric(!same_period & !same_indiv)  # inter-période, individus différents
    alpha3 <- as.numeric(same_indiv & !same_period)   # même individu, périodes différentes

    cbind(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3)
  })

  do.call(rbind, zcor_list)
}






build_zcor <- function(corstr, id, waves, maxwave = NULL, bandwidth = NULL,
                       m = NULL, subgroup = NULL, block = NULL,
                       individual = NULL)
{
  switch(corstr,
         "toeplitz" = build_toeplitz_zcor(id, waves),
         "banded-toeplitz" = build_banded_toeplitz_zcor(id, waves,
                                                        bandwidth = bandwidth),
         "banded-unstructured" = build_banded_unstructured_zcor(
           id, waves, bandwidth = bandwidth
           ),
         "m-dependent" = build_mdep_common_zcor(id, waves, m = m),
         "nested-exchangeable" = build_nested_exch_zcor(id, waves,
                                                        subgroup = subgroup),
         "pairwise-grouped-exchangeable" = build_pairwise_grouped_exch_zcor(
             id,
             waves,
             block = block),
         "block-exchangeable" = build_block_exch_li_zcor(
           id,
           waves,
           individual
           )
         )
}





























