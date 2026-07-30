# DOCUMENTS PROBLEMS WITH ANOVA (USE CAR::ANOVA), DROP1/ADD1/STEP

builders <- c("toeplitz", "banded-toeplitz", "banded-unstructured",
              "m-dependent", "nested-exchangeable",
              "pairwise-grouped-exchangeable", "block-exchangeable")



# subgroup: a column in data. it should correspond, for each wave, to the
#           subgroup it belongs to
#' @export
geefit <- function(formula, data, id, waves = NULL, family = gaussian,
                   corstr = "independence", zcor = NULL,
                   weights, subset, na.action, start = NULL, etastart, mustart,
                   offset, control = geese.control(...), method = "glm.fit",
                   contrasts = NULL, scale.fix = FALSE, scale.value = 1,
                   std.err = "san.se", bandwidth = NULL, mdep = NULL,
                   subgroup = NULL, block = NULL, individual = NULL,
                   ...)
{
  cl <- match.call()
  origcall <- cl

  corstr <- match.arg(corstr,
                      c(c("independence", "exchangeable", "ar1", "unstructured",
                          "fixed", "userdefined"),
                        builders))

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


  geeargs <- setdiff(names(formals(geepack::geeglm)), "...")
  m <- match(geeargs, names(cl), 0L)
  cl <- cl[c(1L, m)]

  cl[[1]] <- quote(geepack::geeglm)


  if (!rlang::quo_is_missing(subExpr)) {
    keep <- eval(rlang::quo_get_expr(subExpr), data)
    data <- data[keep, ]
  }


  if (!rlang::quo_is_missing(naExpr)) {
    keep <- eval(rlang::quo_get_expr(subExpr), data)
    data <- data[keep, ]
  }


  data <- dplyr::arrange(data, !!idCol, !!wavesCol)
  # data <- dplyr::mutate(data,
  #                       .wave_std = dplyr::row_number())


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
    if (rlang::quo_is_null(subgroupCol)) {
      stop("'subgroup' must be provided")
    }

    dsub <- dplyr::select(data, !!wavesCol, !!subgroupCol)
    dssub <- dplyr::arrange(dplyr::distinct(dsub), !!wavesCol)
    subgroup <- dplyr::pull(dssub, !!subgroupCol)

  } else if (corstr == "pairwise-grouped-exchangeable") {

    if (rlang::quo_is_null(blockCol)) {
      stop("'block' must be provided")
    }

    dbl <- dplyr::select(data, !!wavesCol, !!blockCol)
    dbbl <- dplyr::arrange(dplyr::distinct(dbl), !!wavesCol)
    block <- dplyr::pull(dbbl, !!blockCol)

  } else if (corstr == "block-exchangeable") {

    if (rlang::quo_is_null(indCol)) {
      stop("'individual' must be provided")
    }

    individual <- data[[rlang::as_name(indCol)]]
  }

  waves_v <- as.integer(factor(data[[rlang::as_name(wavesCol)]]))

  zcor <- build_zcor(corstr,
                     id = data[[rlang::as_name(idCol)]],
                     waves = waves_v,
                     # waves = data[[rlang::as_name(wavesCol)]],
                     bandwidth = bandwidth, m = mdep,
                     subgroup = subgroup, block = block,
                     individual = individual)


  out <- eval(rlang::call_modify(cl, data = quote(data), waves = NULL,
                                 corstr = "userdefined", zcor = quote(zcor),
                                 subset = rlang::zap()),
              envir = list(data = data,
                           zcor = zcor),
              enclos = parent.frame())


  out$.corstruct <- corstr
  out$call <- origcall
  # out$waves <- data[[rlang::as_name(wavesCol)]]
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


