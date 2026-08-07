# DOCUMENTS PROBLEMS WITH ANOVA (USE CAR::ANOVA), DROP1/ADD1/STEP

builders <- c("toeplitz", "banded-toeplitz", "banded-unstructured",
              "banded-exchangeable", "m-dependent", "nested-exchangeable",
              "pairwise-grouped-exchangeable", "block-exchangeable")
other_structs <- "ar-m"


#
# @section Usage differences with [geepack::geeglm()]:
# waves mandatory|  auto sort


#' Fit GEE Models with Extended Working Correlation Structures
#'
#' @description
#' `geefit()` is a wrapper around [geepack::geeglm()] that extends the set of
#' available working correlation structures beyond those natively supported
#' by **geepack**.
#'
#' @details
#' In addition to supporting more correlation structures, the data is
#' automatically sorted by `id` (contiguous per cluster) and, within
#' each cluster, by `waves` in increasing order.
#'
#' Filtering of the missing values, as well as other subsetting via `subset`,
#' is also performed by `geefit()`. See section **`"Missing Values"`**.
#'
#' Based on the `corstr` argument, `geefit()` dispatches to the appropriate
#' internal `build_*_zcor()` helper to construct the `zcor` matrix passed to
#' `geeglm(..., corstr = "userdefined")`, or falls back to a native
#' **geepack** structure when one already exists.
#'
#' Each working correlation structure defines how pairs of observations
#' within a cluster share (or don't share) correlation parameters. The
#' sections below document every structure supported by `geefit()`, along
#' with equivalent names used in other software, where such an equivalent
#' exists.
#'
#' ## Native geepack structures
#'
#' These are passed through unchanged to [geepack::geeglm()]:
#'
#' - **`"independence"`** — No correlation is assumed between observations
#'   within a cluster. Equivalent to `corstr = "independence"` in the
#'   **gee** package and `TYPE=IND` in SAS PROC GENMOD.
#'
#' - **`"exchangeable"`** — A single correlation parameter is shared by
#'   every pair of observations within a cluster, regardless of lag.
#'   Equivalent to `corstr = "exchangeable"` in **gee** and `TYPE=CS`
#'   (compound symmetry) in SAS.
#'
#' - **`"ar1"`** — First-order autoregressive structure: `corr(lag = k) =
#'   rho^k`. Equivalent to `corstr = "AR-M"` with `Mv = 1` in **gee**, and
#'   `TYPE=AR(1)` in SAS.
#'
#' - **`"unstructured"`** — Every pair of observations has its own, freely
#'   estimated correlation parameter. Equivalent to `corstr =
#'   "unstructured"` in **gee** and `TYPE=UN` in SAS.
#'
#' ## Extended structures
#'
#' These are dispatched to the corresponding `build_*_zcor()` helper:
#'
#' - **`"toeplitz"`** — Correlation depends only on the lag between two
#'   observations, with one freely estimated parameter per lag (no
#'   truncation). Equivalent to `TYPE=TOEP` in SAS.
#'
#' - **`"banded-toeplitz"`** — Toeplitz structure truncated at a maximum
#'   lag `k`: one parameter per lag from 1 to `k`, and correlation fixed
#'   at 0 beyond `k`. Equivalent to `corstr = "stat_M_dep"` (with `Mv =
#'   k`) in **gee**, and `TYPE=TOEP(k+1)` or MDEP(k+1)` in SAS. Also referred to
#'   as "stationary M-dependent" in the literature.
#'
#' - **`"banded-exchangeable"`** — A single shared correlation parameter for
#'   every pair with lag `<= m`, and 0 beyond. Distinct from both
#'   `"banded-toeplitz"` and `"banded_unstructured"` above; no directly
#'   equivalent named option is known in SAS or the **gee** package.
#'
#' - **`"banded-unstructured"`** — Same truncation as
#'   `"banded_toeplitz"` (0 beyond lag `k`), but without pooling
#'   observations within a band: every pair with lag `<= k` keeps its own
#'   separate parameter. Equivalent to `corstr = "non_stat_M_dep"` (with
#'   `Mv = k`) in **gee**. Also referred to as "nonstationary M-dependent"
#'   in the literature.
#'
#' - **`"m-dependent"`** — alias for **`"banded-toeplitz"`**.
#'
#' - **`"nested-exchangeable"`** — Observations are partitioned into
#'   subgroups (e.g., left/right eye, or visit-level clustering of
#'   repeated measurements). Two parameters: one shared correlation for
#'   pairs within the same subgroup, one shared correlation for pairs
#'   across different subgroups. This is the cross-sectional (non-cohort)
#'   case described by Li, Turner, and Preisser (2018) for cluster
#'   randomized trials.
#'
#' - **`"pairwise-grouped-exchangeable"`** — Generalization of
#'   `"nested_exchangeable"` to more than two groups, with a dedicated
#'   correlation parameter for each group (intra-group) and a dedicated
#'   parameter for each distinct pair of groups (inter-group). For `k`
#'   groups this yields `k(k+1)/2` parameters. Not a structure with a
#'   standard name in the literature; this naming is specific to this
#'   package.
#'
#' - **`"block-exchangeable"`** — Three-parameter structure for cohort
#'   designs (the same individual followed across multiple periods
#'   within a cluster): `alpha1` for pairs in the same period (different
#'   individuals), `alpha2` for pairs in different periods (different
#'   individuals), and `alpha3` for pairs of the same individual across
#'   different periods. As defined by Li, Turner, and Preisser (2018);
#'   requires an additional `individual` identifier beyond `id` and
#'   `waves`. Reduces to `"nested_exchangeable"` when no individual is
#'   followed across periods.
#'
#' ## Structures requiring iterative fitting
#'
#' Unlike every structure above, `"ar-m"` cannot be expressed as a fixed
#' linear combination of indicator columns, and is therefore not a simple
#' `zcor` build:
#'
#' - **`"ar-m"`** — Autoregressive structure of order `Mv > 1`. Unlike the
#'   banded structures above, correlation does not vanish beyond lag
#'   `Mv`; it continues to decay according to the autoregressive
#'   recursion (Yule-Walker equations) implied by the `Mv` autoregressive
#'   coefficients. Equivalent to `corstr = "AR-M"` with `Mv > 1` in the
#'   **gee** package. Because the correlation is a non-linear function of
#'   the underlying parameters, it is instead fit via `geeglm_arm()`, an
#'   iterative procedure alternating between banded-Toeplitz estimation
#'   and a Yule-Walker update.
#'
#'
#' # Missing Values
#'
#' **Very first step**
#'
#'
#' @param formula a two-sided formula, as in [geepack::geeglm()].
#' @param data a data frame.
#' @param id  cluster identifier column in `data`.
#'
#' `r lifecycle::badge("experimental")`
#'
#'     Can also be a vector of clusters, the same length as the number of rows
#'     of `data`, prior to removing missing values or subsetting. Must
#'     corresponds to the order of `data`.
#' @param waves column in `data`; variable giving the true time/position of each
#'  observation within its cluster. Required for all correlation structures to
#'  carry out sorting.
#'
#' `r lifecycle::badge("experimental
#'
#'     Can also be a vector, the same length as the number of rows of `data`,
#'     prior to removing missing values or subsetting. Must corresponds to the
#'     order of `data`.
#' @param family as in [geepack::geeglm()].
#' @param corstr character string naming the working correlation
#'   structure. One of `"independence"`, `"exchangeable"`, `"ar1"`,
#'   `"unstructured"`, `r nice_collapse(builders)`.
#' @param zcor user-defined correlation structure. Mandatory if `corstr` is
#'   `"userdefined"`, ignored for all other structures.
#' @param weights See corresponding documentation to glm
#' @param subset See corresponding documentation to glm
#' @param na.action See corresponding documentation to glm
#' @param start See corresponding documentation to glm
#' @param etastart See corresponding documentation to glm
#' @param mustart See corresponding documentation to glm
#' @param offset See corresponding documentation to glm
#' @param control See corresponding documentation to geeglm
#' @param method See corresponding documentation to glm
#' @param contrasts See corresponding documentation to glm
#' @param std.err See corresponding documentation to geeglm
#'
#' @param bandwidth truncation lag `k`: correlation is fixed at 0 beyond this
#'   lag. Required when `corstr` is `"banded-toeplitz"`, `"banded-exchangeable"`
#'   or `"banded-unstructured"`.
#' @param mdep integer; same as `"bandwidth"` in `"banded-toeplitz"`. Required
#'   when `corstr = "m-dependent"`.
#' @param subgroup Column in `data`; variable giving the subgroup of each
#'   observation within its cluster. Required when `corstr` is
#'   `"nested-exchangeable"`.
#'
#' `r lifecycle::badge("experimental")`
#'
#'     Can also be a vector, the same length as the number of rows of `data`,
#'     prior to removing missing values or subsetting. Must corresponds to the
#'     order of `data`.
#' @param block Column in `data`; variable giving the block id of each
#'    observation within its cluster. Required when `corstr` is
#'    `"pairwise-grouped-exchangeable"`.
#'
#' `r lifecycle::badge("experimental")`
#'
#'     Can also be a vector, the same length as the number of rows of `data`,
#'     prior to removing missing values or subsetting. Must corresponds to the
#'     order of `data`.
#' @param individual Column in `data`; variable giving the individual id of
#'    each observation within its cluster. Required when `corstr` is
#'    `"group-exchangeable"`.
#'
#'     This is not the same as the cluster ID provided with
#'     `id`, but the individual identifier, unique within a cluster, that
#'     persists across periods (waves) for individuals followed in a cohort
#'     (same subject revisited at each period), or that is distinct at each
#'     observation if the design is cross-sectional for that individual.
#'
#' `r lifecycle::badge("experimental")`
#'
#'     Can also be a vector, the same length as the number of rows of `data`,
#'     prior to removing missing values or subsetting. Must corresponds to the
#'     order of `data`.
#' @param ... further arguments passed to or from other methods.
#'
#' @references
#' Li, F., Turner, E. L., & Preisser, J. S. (2018). Sample size
#' determination for GEE analyses of stepped wedge cluster randomized
#' trials. *Biometrics*, 74(4), 1450-1458.
#'
#' Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. *Biometrika*, 73(1), 13-22.
#'
#' @export
geefit <- function(formula, data, id, waves = NULL, family = gaussian,
                   corstr = "independence", zcor = NULL,
                   weights, subset, na.action, start = NULL, etastart, mustart,
                   offset, control = geese.control(...), method = "glm.fit",
                   contrasts = NULL, std.err = "san.se",
                   bandwidth = NULL, Mv = NULL, mdep = NULL,
                   subgroup = NULL, block = NULL, individual = NULL,
                   ...)
{
  cl <- match.call()
  origcall <- cl

  corstr <- match.arg(corstr,
                      c(c("independence", "exchangeable", "ar1", "unstructured",
                          "fixed", "userdefined"),
                        c(builders, other_structs)))

  idCol <- rlang::ensym(id)
  wavesCol <- rlang::enquo(waves)
  subgroupCol <- rlang::enquo(subgroup)
  blockCol <- rlang::enquo(block)
  indCol <- rlang::enquo(individual)

  subExpr <- rlang::enquo(subset)
  naExpr <- rlang::enquo(na.action)


  if (rlang::quo_is_null(wavesCol)) {
    stop("'waves' must be provided.")
  }


  if (corstr == "ar-m") {

    a <- match(c("formula", "data", "id", "waves", "family", "Mv"),
               names(cl), 0L)

    cl <- cl[c(1L, a)]
    cl[[1]] <- quote(geeglm_arm)

    cl <- rlang::call_modify(cl, data = quote(data))


  } else {

    geeargs <- setdiff(names(formals(geepack::geeglm)), "...")
    m <- match(geeargs, names(cl), 0L)
    cl <- cl[c(1L, m)]

    cl[[1]] <- quote(geepack::geeglm)

  }


  dat_nms <- colnames(data)

  id_name <- rlang::as_name(idCol)
  if (!id_name %in% dat_nms) {
    data[[id_name]] <- id
  }


  waves_name <- rlang::as_name(wavesCol)
  if (!waves_name %in% dat_nms) {
    data[[waves_name]] <- waves
  }


  if (corstr == "nested-exchangeable") {
    if (rlang::quo_is_null(subgroupCol)) {
      stop("'subgroup' must be provided")
    }

    subgroup_name <- rlang::as_name(subgroupCol)
    if (!subgroup_name %in% dat_nms) {
      data[[subgroup_name]] <- subgroup
    }
  }


  if (corstr == "pairwise-grouped-exchangeable") {
    if (rlang::quo_is_null(blockCol)) {
      stop("'block' must be provided")
    }

    block_name <- rlang::as_name(blockCol)
    if (!block_name %in% dat_nms) {
      data[[block_name]] <- block
    }
  }


  if (corstr == "block-exchangeable") {
    if (rlang::quo_is_null(indCol)) {
      stop("'individual' must be provided")
    }

    ind_name <- rlang::as_name(indCol)
    if (!ind_name %in% dat_nms) {
      data[[ind_name]] <- individual
    }
  }




  if (!rlang::quo_is_missing(subExpr)) {
    keep <- eval(rlang::quo_get_expr(subExpr), data)
    data <- data[keep, ]
  }


  if (rlang::quo_is_missing(naExpr)) {
    na.action <- if (!is.null(naa <- attr(data, "na.action")) &&
                     mode(naa) != "numeric")
      naa
    else {
      naa <- getOption("na.action")
      if (!is.null(naa) && is.character(naa)) {
        naa <- get(naa, mode = "function",
                   envir = as.environment("package:stats"))
      }
      if (is.null(naa)) naa <- na.fail
      naa
    }
  }

  if (anyNA(data)) {
    added_vars <- tidyselect::eval_select(
      rlang::expr(c(!!idCol, !!wavesCol, !!subgroupCol, !!blockCol, !!indCol)),
      data,
      strict = FALSE
    )
    vars <- union(all.vars(formula), names(added_vars))
    data <- dplyr::select(data, tidyselect::all_of(vars))
    data <- naa(data)

#     upd_fm <- reformulate(c(".", encodeString(names(added_vars), quote = "`")),
#                           quote(.))
#     mf_formula <- update.formula(formula, upd_fm)
# print(naa)
#     print(rlang::new_formula(NULL, rlang::expr(. + !!idCol + !!wavesCol)))
#     print(update.formula(formula, . ~ . + rlang::quo_get_expr(idCol))   )
#     data <- na.action(data)
  }


  data <- dplyr::arrange(data, !!idCol, !!wavesCol)


  if (corstr %in% c("independence", "exchangeable", "ar1", "unstructured",
                    "fixed", "userdefined")) {

    return(
      eval(
        rlang::call_modify(cl, data = quote(data), waves = NULL),
        envir = list(data = data),
        enclos = parent.frame()
      )
    )
  }




  subgroup <- block <- individual <- NULL


  if (corstr == "nested-exchangeable") {

    dsub <- dplyr::select(data, !!wavesCol, !!subgroupCol)
    dssub <- dplyr::arrange(dplyr::distinct(dsub), !!wavesCol)
    subgroup <- dplyr::pull(dssub, !!subgroupCol)

  } else if (corstr == "pairwise-grouped-exchangeable") {

    dbl <- dplyr::select(data, !!wavesCol, !!blockCol)
    dbbl <- dplyr::arrange(dplyr::distinct(dbl), !!wavesCol)
    block <- dplyr::pull(dbbl, !!blockCol)

  } else if (corstr == "block-exchangeable") {

    individual <- data[[ind_name]]
  }

  waves_v <- as.integer(factor(data[[waves_name]]))

  if (corstr == "ar-m") {

    out <- eval(cl, envir = list(data = data), enclos = parent.frame())

  } else {

    zcor <- build_zcor(corstr,
                       id = data[[id_name]],
                       waves = waves_v,
                       bandwidth = bandwidth, m = mdep,
                       subgroup = subgroup, block = block,
                       individual = individual)


    out <- eval(rlang::call_modify(cl, data = quote(data), waves = NULL,
                                   corstr = "userdefined", zcor = quote(zcor),
                                   subset = rlang::zap()),
                envir = list(data = data,
                             zcor = zcor),
                enclos = parent.frame())

  }





  out$.corstruct <- corstr
  out$call <- origcall
  out$waves <- waves_v
  out$.corparams <- list(bandwidth = bandwidth, mdep = mdep,
                         subgroup = subgroup, block = block,
                         individual = individual)

  class(out) <- c("geecor", class(out))

  out
}




#' @noRd
geecor_fit <- function(x, y, id, waves, family = gaussian,
                       corstr = "independence", zcor = NULL,
                       weights, subset, na.action, start = NULL,
                       etastart, mustart, offset, control = geese.control(...),
                       contrasts = NULL, bandwidth = NULL, mdep = NULL,
                       subgroup = NULL, block = NULL, individual = NULL, ...)
{
  corstr <- match.arg(corstr,
                      c(c("independence", "exchangeable", "ar1", "unstructured",
                          "fixed", "userdefined"),
                        builders))


  N <- NROW(y)
  soffset <- rep(0, N)
  zsca <- matrix(1, N, 1)

  if (corstr == "userdefined" & is.null(zcor)) {
    stop("'zcor' must be provided for corstr == 'userdefined'")
  }


  if (missing(offset)) offset <- rep(0, N)
  if (missing(weights)) weights <- rep(1, N)
  if (family$family == "binomial") {
    if (is.matrix(y) && ncol(y) == 2) {
      weights <- apply(y, 1, sum)
      y <- y[, 1]/weights
    }
  }


  if (!corstr %in% c("independence", "exchangeable", "ar1", "unstructured",
                    "fixed")) {
    zcor <- build_zcor(corstr, id = id, waves = waves, bandwidth = bandwidth,
                       m = mdep, subgroup = subgroup, block = block,
                       individual = individual)
    corstr <- "userdefined"
  }



  ans <- geese.fit(x = x, y = y, id = id, offset = offset, soffset = soffset,
                   weights = weights, waves = waves, zsca = zsca, zcor = zcor,
                   corp = NULL, control = control, b = start, alpha = NULL,
                   gm = NULL, family = family, mean.link = NULL,
                   variance = NULL, cor.link = "identity", sca.link = "identity",
                   link.same = TRUE, scale.fix = FALSE,
                   scale.value = 1, corstr = corstr)

  ans
  # ans <- c(ans, list(call = call, formula = formula))
  # class(ans) <- "geese"
  # ans$X <- x
  # ans$id <- id
  # ans$weights <- weights
  #
  #
  # out <- glmFit
  # toDelete <- c("R", "deviance", "aic", "null.deviance", "iter",
  #               "df.null", "converged", "boundary")
  # out[match(toDelete, names(out))] <- NULL
  # out$method <- "geese.fit"
  # out$geese <- ans
  # out$weights <- ans$weights
  # out$coefficients <- ans$beta
  # out$offset <- offset
  # if (is.null(out$offset)) {
  #   out$linear.predictors <- ans$X %*% ans$beta
  # }
  # else {
  #   out$linear.predictors <- out$offset + ans$X %*% ans$beta
  # }
  # out$fitted.values <- family(out)$linkinv(out$linear.predictors)
  # out$modelInfo <- ans$model
  # out$id <- ans$id
  # out$call <- ans$call
  # out$corstr <- ans$model$corstr
  # out$cor.link <- ans$model$cor.link
  # out$control <- ans$control
  # out$std.err <- std.err
  # class(out) <- c("geeglm", "gee", "glm", "lm")
  # out
}


