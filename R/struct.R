# trié dat par id puis waves
# pour grouped_exchangeable: par id (cluster), puis period puis subj_id



build_toeplitz_zcor <- function(zcor.unstr, maxwave, pairs, lags)
{
  maxlag <- maxwave - 1
  zcor.toep <- sapply(seq_len(maxlag), function(l) {
    cols <- which(lags == l)
    if (length(cols) == 1) zcor.unstr[, cols] else rowSums(zcor.unstr[, cols, drop = FALSE])
  })

  zcor.toep
}




build_banded_toeplitz_zcor <- function(zcor.unstr, maxwave, pairs, lags,
                                       bandwidth)
{
  # id      : identifiant du cluster/sujet
  # waves   : temps de mesure réel (pas l'ordre des lignes)
  # bandwidth       : largeur de bande (nombre de lags distincts à estimer, bandwidth >= 1)
  # maxwave : nb total de vagues possibles (déduit des données si non fourni)

  stopifnot(bandwidth >= 1, bandwidth <= maxwave - 1)

  # une colonne par lag de 1 à bandwidth ; les lags > bandwidth ne sont assignés à aucune
  # colonne -> ligne nulle -> corrélation fixée à 0 pour ces paires
  zcor.band <- sapply(seq_len(bandwidth), function(l) {
    cols <- which(lags == l)
    if (length(cols) == 1) zcor.unstr[, cols] else rowSums(zcor.unstr[, cols, drop = FALSE])
  })

  zcor.band
}






build_banded_unstructured_zcor <- function(zcor.unstr, maxwave, pairs, lags,
                                           bandwidth) {
  # id      : identifiant du cluster/sujet
  # waves   : temps de mesure réel (pas l'ordre des lignes)
  # bandwidth       : largeur de bande (lag maximal ayant une corrélation non nulle)
  # maxwave : nb total de vagues possibles (déduit des données si non fourni)

  stopifnot(bandwidth >= 1, bandwidth <= maxwave - 1)

  # on garde chaque paire de lag <= bandwidth comme sa PROPRE colonne (pas de regroupement)
  cols <- which(lags <= bandwidth)
  zcor.band.unstr <- zcor.unstr[, cols, drop = FALSE]

  colnames(zcor.band.unstr) <- paste0("pair_", pairs[1, cols], "_", pairs[2, cols])
  zcor.band.unstr
}






build_banded_exchangeable_zcor <- function(zcor.unstr, maxwave, pairs, lags,
                                           bandwidth) {
  # id      : identifiant du cluster/sujet
  # waves   : temps de mesure réel (pas l'ordre des lignes)
  # bandwidth       : ordre de dépendance ; un seul paramètre commun pour les lags 1..m
  # maxwave : nb total de vagues possibles (déduit des données si non fourni)

  stopifnot(bandwidth >= 1, bandwidth <= maxwave - 1)

  # toutes les paires avec lag <= bandwidth regroupées dans UNE seule colonne
  # (un seul paramètre) ; lags > bandwidth -> exclues -> corrélation fixée à 0
  cols <- which(lags <= bandwidth)
  zcor.mdep <- matrix(rowSums(zcor.unstr[, cols, drop = FALSE]), ncol = 1)

  zcor.mdep
}





build_nested_exch_zcor <- function(zcor.unstr, maxwave, pairs, lags, subgroup) {
  # id       : identifiant du patient/cluster (le niveau GEE le plus haut)
  # waves    : temps/position réel de chaque observation (pas l'ordre des lignes)
  # subgroup : vecteur de longueur maxwave donnant le sous-groupe de chaque
  #            vague (ex: c("G","G","D","D") pour 2 mesures oeil gauche,
  #            2 mesures oeil droit ; ou le numéro de visite si plusieurs
  #            mesures par visite)
  # maxwave  : nb total de vagues possibles (déduit de waves si non fourni)

  stopifnot(length(subgroup) == maxwave)
  same_subgroup <- subgroup[pairs[1, ]] == subgroup[pairs[2, ]]

  col_alpha1 <- rowSums(zcor.unstr[, same_subgroup, drop = FALSE])   # même sous-groupe
  col_alpha2 <- rowSums(zcor.unstr[, !same_subgroup, drop = FALSE])  # sous-groupe différent

  cbind(alpha1 = col_alpha1, alpha2 = col_alpha2)
}






build_pairwise_grouped_exch_zcor <- function(zcor.unstr, maxwave, pairs, lags,
                                             block) {
  # id      : identifiant du cluster/sujet
  # waves   : temps de mesure réel (pas l'ordre des lignes)
  # block   : vecteur de longueur maxwave donnant le numéro de bloc pour
  #           chaque vague (ex: c(1,1,2,2,2,3) pour 3 blocs)
  # maxwave : nb total de vagues possibles (déduit de waves si non fourni)

  stopifnot(length(block) == maxwave)

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
  #              individual identifier, unique within a cluster, that persists across periods for individuals followed in a cohort (same subject revisited at each period), or that is distinct at each observation if the design is cross-sectional for that individual

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





build_unstructured_zcor <- function(id, waves, maxwave = NULL) {
  if (is.null(maxwave)) maxwave <- max(waves)

  id_f <- factor(id, levels = unique(id))
  clusters <- split(waves, id_f)

  pairs <- combn(maxwave, 2)
  ncol_z <- ncol(pairs)
  pair_index <- function(v1, v2) which(pairs[1, ] == v1 & pairs[2, ] == v2)

  rows <- list()
  for (w in clusters) {
    n <- length(w)
    if (n < 2) next
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      v1 <- min(w[i], w[j]); v2 <- max(w[i], w[j])
      col <- pair_index(v1, v2)
      row <- numeric(ncol_z)
      row[col] <- 1
      rows[[length(rows) + 1]] <- row
    }
  }
  zcor <- do.call(rbind, rows)
  colnames(zcor) <- paste0("alpha.", pairs[1, ], ":", pairs[2, ])
  zcor
}






#' Build a `zcor` Matrix for Extended Working Correlation Structures
#'
#' `build_zcor()` builds the `zcor` matrix based on `corstr`, ready to be passed
#' to [geepack::geeglm()] as
#' `geeglm(..., corstr = "userdefined", zcor = build_zcor(...))`.
#'
#' @param corstr character string naming the working correlation
#'   structure. One of `"independence"`, `"exchangeable"`, `"ar1"`,
#'   `"unstructured"`, `r nice_collapse(builders)`. See **Details** for the
#'   argument(s) required by each.
#' @param id cluster identifier (the GEE unit, e.g., the cluster/site in a
#'   cluster randomized trial)
#' @param waves actual time/position of each observation (not the row order)
#'   within each cluster. For `"block-exchangeable"`, can recur within the same
#'   cluster, since several individuals share the same period (wave).
#' @param maxwave total number of possible waves. Deduced from `max(waves)`
#'   when `NULL`. Used by every structure except `"block_exchangeable"`.
#' @param bandwidth truncation lag `k`: correlation is fixed at 0 beyond this
#'   lag. Required when `corstr` is `"banded-toeplitz"`, `"banded-exchangeable"`
#'   or `"banded-unstructured"`; ignored otherwise.
#' @param m integer; for `"m-dependent"`, same as `"bandwidth"` in
#'   `"banded-toeplitz"`. Required when `corstr = "m-dependent"`; ignored
#'   otherwise.
#' @param subgroup vector of length `maxwave` giving the subgroup of each
#'   wave. Required when `corstr = "nested_exchangeable"`; ignored otherwise.
#'
#'     e.g., c("L","L","R","R") for 2 left-eye measurements, 2 right-eye
#'     measurements for 4 waves; or the visit number if there are multiple
#'     measurements per visit.
#' @param block vector of length `maxwave` giving the group of each wave.
#'   Required when `corstr = "pairwise_grouped_exchangeable"`; ignored
#'   otherwise.
#'
#'     e.g.,: c(1,1,2,2,2,3) for 3 blocks and 6 waves.
#' @param individual identifier of the individual within each cluster,
#'   persisting across periods for individuals followed in a cohort design.
#'   Required when `corstr = "block_exchangeable"`; ignored otherwise.
#'
#'
#' @details
#' `corstr` determines both the dispatch target and which of the remaining
#' arguments are required:
#'
#' | `corstr`                          | requires                   |
#' | :-------------------------------- | :------------------------- |
#' | `"toeplitz"`                      | `id`, `waves`              |
#' | `"banded-toeplitz"`               | `id`, `waves`, `bandwidth` |
#' | `"banded-unstructured"`           | `id`, `waves`, `bandwidth` |
#' | `"banded-exchangeable"`           | `id`, `waves`, `bandwidth` |
#' | `"m-dependent"`                   | `id`, `waves`, `m`         |
#' | `"nested-exchangeable"`           | `id`, `waves`, `subgroup`  |
#' | `"pairwise-grouped-exchangeable"` | `id`, `waves`, `block`     |
#' | `"block-exchangeable"`            | `id`, `waves`, `individual`|
#'
#' `data` must be sorted by `id` (contiguous per cluster) and, within each
#' cluster, in the same row order used to build `waves` (and `individual`,
#' when applicable) — see the package vignette for details on required
#' sorting per structure.
#'
#' @return
#' A numeric matrix suitable for the `zcor` argument of
#' `geeglm(..., corstr = "userdefined")`: one row per observed pair within a
#' cluster, one column per correlation parameter to estimate.
#'
#' @seealso [geepack::geeglm()], [geepack::genZcor()]
#'
#' @export
#'
build_zcor <- function(corstr, id, waves, maxwave = NULL, bandwidth = NULL,
                       m = NULL, subgroup = NULL, block = NULL,
                       individual = NULL)
{
  if (is.null(maxwave)) maxwave <- max(waves)

  # unstructured design matrix (natively handles imbalance)
  zcor.unstr <- build_unstructured_zcor(id, waves, maxwave)

  # correspondence between column and pair (s, t), s < t, in the order of combn()
  pairs <- combn(maxwave, 2)
  lags  <- abs(pairs[2, ] - pairs[1, ])



  switch(corstr,
         "toeplitz" = build_toeplitz_zcor(zcor.unstr, maxwave, pairs, lags),
         "m-dependent" = build_banded_toeplitz_zcor(zcor.unstr, maxwave,
                                                    pairs, lags,
                                                    bandwidth = m),
         "banded-toeplitz" = build_banded_toeplitz_zcor(zcor.unstr, maxwave,
                                                        pairs, lags,
                                                        bandwidth = bandwidth),
         "banded-unstructured" = build_banded_unstructured_zcor(
           zcor.unstr, maxwave, pairs, lags, bandwidth = bandwidth
           ),
         "banded-exchangeable" = build_banded_exchangeable_zcor(
           zcor.unstr, maxwave, pairs, lags, bandwidth = bandwidth),
         "nested-exchangeable" = build_nested_exch_zcor(zcor.unstr, maxwave,
                                                        pairs, lags,
                                                        subgroup = subgroup),
         "pairwise-grouped-exchangeable" = build_pairwise_grouped_exch_zcor(
           zcor.unstr, maxwave, pairs, lags, block = block),
         "block-exchangeable" = build_block_exch_li_zcor(
           id,
           waves,
           individual
           ),
         "unrecognized correlation structure"
         )
}





























