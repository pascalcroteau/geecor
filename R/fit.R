

builders <- c("toeplitz", "banded-toeplitz", "m-dependent",
              "nested-exchangeable", "pairwise-grouped-exchangeable",
              "block-exchangeable")



# subgroup: a column in data. it should correspond, for each wave, to the
#           subgroup it belongs to
#' @export
geefit <- function(formula, data, family = gaussian,
                   weights, subset, na.action, start = NULL, etastart, mustart,
                   offset, control = geese.control(...), method = "glm.fit",
                   contrasts = NULL, id, waves = NULL, zcor = NULL,
                   corstr = "independence",
                   scale.fix = FALSE, scale.value = 1, std.err = "san.se",
                   bandwidth = NULL, mdep = NULL, subgroup = NULL, block = NULL,
                   individual = NULL,
                   ...)
{
  cl <- match.call()

  corstr <- match.arg(corstr,
                      c(c("independence", "exchangeable", "ar1", "unstructured",
                          "userdefined"),
                        builders))

  idCol <- rlang::ensym(id)
  wavesCol <- rlang::enquo(waves)
  subgroupCol <- rlang::enquo(subgroup)
  blockCol <- rlang::enquo(block)
  indCol <- rlang::enquo(individual)

  subExpr <- rlang::enquo(subset)


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


  data <- dplyr::arrange(data, !!idCol, !!wavesCol)
  # data <- dplyr::mutate(data,
  #                       .wave_std = dplyr::row_number())


  if (corstr %in% c("independence", "exchangeable", "ar1", "unstructured",
                    "userdefined")) {

    return(
      eval(
        rlang::call_modify(cl, data = quote(data), waves = NULL),
        envir = list(data = data),
        enclos = parent.frame()
      )
    )
  }



  zcor <- switch(corstr,
                 "toeplitz" = build_toeplitz_zcor(data[[rlang::as_name(idCol)]],
                                                  data[[rlang::as_name(wavesCol)]]),
                 "banded-toeplitz" = build_banded_toeplitz_zcor(
                   data[[rlang::as_name(idCol)]],
                   data[[rlang::as_name(wavesCol)]],
                   bandwidth = bandwidth),
                 "m-dependent" = build_mdep_common_zcor(
                   data[[rlang::as_name(idCol)]],
                   data[[rlang::as_name(wavesCol)]],
                   m = mdep),
                 "nested-exchangeable" = {
                   if (rlang::quo_is_null(subgroupCol)) {
                     stop("'subgroup' must be provided")
                   }

                   dsub <- dplyr::select(data, !!wavesCol, !!subgroupCol)
                   dssub <- dplyr::arrange(dplyr::distinct(dsub), !!wavesCol)
                   subg <- dplyr::pull(dssub, !!subgroupCol)

                   build_nested_exch_zcor(
                     data[[rlang::as_name(idCol)]],
                     data[[rlang::as_name(wavesCol)]],
                     subgroup = subg)
                 },
                 "pairwise-grouped-exchangeable" = {
                   if (rlang::quo_is_null(blockCol)) {
                     stop("'block' must be provided")
                   }

                   dbl <- dplyr::select(data, !!wavesCol, !!blockCol)
                   dbbl <- dplyr::arrange(dplyr::distinct(dbl), !!wavesCol)
                   bl <- dplyr::pull(dbbl, !!blockCol)

                   build_pairwise_grouped_exch_zcor(
                     data[[rlang::as_name(idCol)]],
                     data[[rlang::as_name(wavesCol)]],
                     block = bl)
                 },
                 "block-exchangeable" = build_block_exch_li_zcor(
                   data[[rlang::as_name(idCol)]],
                   data[[rlang::as_name(wavesCol)]],
                   data[[rlang::as_name(indCol)]]))



  out <- eval(rlang::call_modify(cl, data = quote(data), waves = NULL,
                                 corstr = "userdefined", zcor = quote(zcor),
                                 subset = rlang::zap()),
              envir = list(data = data,
                           zcor = zcor),
              enclos = parent.frame())


  out$corstr <- corstr

  class(out) <- c("geecor", class(out))

  out
}
