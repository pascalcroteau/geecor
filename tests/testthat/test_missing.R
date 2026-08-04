# library(geepack)
#
# # ============================================================================
# # NOTE : ce script n'a pas pu être exécuté ni vérifié (pas de R disponible
# # dans l'environnement où il a été écrit). Testez-le et signalez toute erreur.
# # ============================================================================
#
# # ============================================================================
# # PARTIE 1 : structures NATIVES de geepack (independence/exchangeable/ar1/
# #            unstructured) - jeu de données public "respiratory" (Davis 1991,
# #            inclus dans geepack), avec comparaison à des coefficients publiés
# #            (SAS) pour une variante de ce même jeu, ET un test de cohérence
# #            interne geeglm(na.omit) vs geefit() face aux valeurs manquantes.
# # ============================================================================
#
#
test_that(
  "testing missingness",
  {

    # ---------------------------------------------------------------------
    # Native structures - example from geepack
    # ---------------------------------------------------------------------


    respiratory <- geepack::respiratory
    respiratory$center <- factor(respiratory$center)

    gee.ind <- geeglm(outcome ~ center + treat + age + baseline, data = respiratory,
                      id = id, family = binomial(), corstr = "independence")
    gee.exc <- geeglm(outcome ~ center + treat + age + baseline, data = respiratory,
                      id = id, family = binomial(), corstr = "exchangeable")
    gee.uns <- geeglm(outcome ~ center + treat + age + baseline, data = respiratory,
                      id = id, family = binomial(), corstr = "unstructured")
    gee.ar1 <- geeglm(outcome ~ center + treat + age + baseline, data = respiratory,
                      id = id, family = binomial(), corstr = "ar1")




    geecor.ind <- geeglm(outcome ~ center + treat + age + baseline,
                         data = respiratory,
                         id = id, waves = visit,
                         family = binomial(), corstr = "independence")
    geecor.exc <- geeglm(outcome ~ center + treat + age + baseline,
                         data = respiratory,
                         id = id, waves = visit,
                         family = binomial(), corstr = "exchangeable")
    geecor.uns <- geeglm(outcome ~ center + treat + age + baseline,
                         data = respiratory,
                         id = id, waves = visit,
                         family = binomial(), corstr = "unstructured")
    geecor.ar1 <- geeglm(outcome ~ center + treat + age + baseline,
                         data = respiratory,
                         id = id, waves = visit,
                         family = binomial(), corstr = "ar1")



    expect_identical(coef(gee.ind), coef(geecor.ind))
    expect_identical(gee.ind$geese$alpha, geecor.ind$geese$alpha)
    expect_identical(coef(gee.exc), coef(geecor.exc))
    expect_identical(gee.exc$geese$alpha, geecor.exc$geese$alpha)
    expect_identical(coef(gee.uns), coef(geecor.uns))
    expect_identical(gee.uns$geese$alpha, geecor.uns$geese$alpha)
    expect_identical(coef(gee.ar1), coef(geecor.ar1))
    expect_identical(gee.ar1$geese$alpha, geecor.ar1$geese$alpha)




    set.seed(1)
    respiratory_na <- respiratory
    na_rows <- sample(seq_len(nrow(respiratory_na)),
                      size = round(0.1 * nrow(respiratory_na)))
    respiratory_na$age[na_rows] <- NA


    # geeglm
    respiratory_clean <- respiratory_na[complete.cases(
      respiratory_na[, c("id", "outcome", "center", "treat", "age", "baseline")]), ]
    respiratory_clean <- respiratory_clean[order(respiratory_clean$id), ]

    fit_manual <- geeglm(outcome ~ center + treat + age + baseline,
                         data = respiratory_clean,
                         id = id, family = binomial(), corstr = "exchangeable")

    # geecor
    fit_geefit <- geefit(outcome ~ center + treat + age + baseline,
                         id = id, waves = visit,
                         data = respiratory_na, family = binomial(),
                         corstr = "exchangeable")

    expect_equal(coef(fit_manual), coef(fit_geefit))
    expect_identical(fit_manual$geese$alpha, fit_geefit$geese$alpha)







    # ---------------------------------------------------------------------
    # New structures -
    # ---------------------------------------------------------------------
    #


    dat_band_clean <- dat_band_na[complete.cases(dat_band_na), ]
    dat_band_clean <- dat_band_clean[order(dat_band_clean$id,
                                           dat_band_clean$waves), ]

    zcor_manual <- build_zcor("banded-toeplitz", id = dat_band_clean$id,
                              waves = dat_band_clean$waves, bandwidth = 2)
    fit_manual_band <- geeglm(y ~ 1, id = id, data = dat_band_clean, family = gaussian,
                              corstr = "userdefined", zcor = zcor_manual)


    fit_geefit_band <- geefit(y ~ 1, id = id, waves = waves, data = dat_band_na,
                              family = gaussian, corstr = "banded-toeplitz",
                              bandwidth = 2)



    expect_identical(coef(fit_manual_band), coef(fit_geefit_band))
    expect_identical(fit_manual_band$geese$alpha, fit_geefit_band$geese$alpha)



    # ---------------------------------------------------------------------
    # banded-unstructured.
    # ---------------------------------------------------------------------

    dat_bu_clean <- dat_bu_na[complete.cases(dat_bu_na), ]
    dat_bu_clean <- dat_bu_clean[order(dat_bu_clean$id,
                                       dat_bu_clean$waves), ]

    zcor_manual <- build_zcor("banded-unstructured", id = dat_bu_clean$id,
                              waves = dat_bu_clean$waves, bandwidth = 2)
    fit_manual_bu <- geeglm(y ~ 1, id = id, data = dat_bu_clean, family = gaussian,
                              corstr = "userdefined", zcor = zcor_manual)


    fit_geefit_bu <- geefit(y ~ 1, id = id, waves = waves, data = dat_bu_na,
                              family = gaussian, corstr = "banded-unstructured",
                              bandwidth = 2)



    expect_identical(coef(fit_manual_bu), coef(fit_geefit_bu))
    expect_identical(fit_manual_bu$geese$alpha, fit_geefit_bu$geese$alpha)




    # ---------------------------------------------------------------------
    # banded-exchangeable. The answer must be close to 0.4
    # ---------------------------------------------------------------------


    dat_bdex_clean <- dat_bdex_na[complete.cases(dat_bdex_na), ]
    dat_bdex_clean <- dat_bdex_clean[order(dat_bdex_clean$id,
                                           dat_bdex_clean$waves), ]

    zcor_manual <- build_zcor("banded-exchangeable", id = dat_bdex_clean$id,
                              waves = dat_bdex_clean$waves, bandwidth = 2)
    fit_manual_bdex <- geeglm(y ~ 1, id = id, data = dat_bdex_clean, family = gaussian,
                              corstr = "userdefined", zcor = zcor_manual)


    fit_geefit_bdex <- geefit(y ~ 1, id = id, waves = waves, data = dat_bdex_na,
                              family = gaussian, corstr = "banded-exchangeable",
                              bandwidth = 2)



    expect_identical(coef(fit_manual_bdex), coef(fit_geefit_bdex))
    expect_identical(fit_manual_bdex$geese$alpha, fit_geefit_bdex$geese$alpha)




    # ---------------------------------------------------------------------
    # nested-exchangeable: The answer must be close to 0.5 0.2
    # ---------------------------------------------------------------------

    dat_nested_clean <- dat_nested_na[complete.cases(dat_nested_na), ]
    dat_nested_clean <- dat_nested_clean[order(dat_nested_clean$id,
                                               dat_nested_clean$waves), ]

    zcor_manual <- build_zcor("nested-exchangeable", id = dat_nested_clean$id,
                              waves = dat_nested_clean$waves,
                              subgroup = dplyr::distinct(dat_nested_clean, waves, sub)$sub)
    fit_manual_nested <- geeglm(y ~ 1, id = id, data = dat_nested_clean, family = gaussian,
                                corstr = "userdefined", zcor = zcor_manual)


    fit_geefit_nested <- geefit(y ~ 1, id = id, waves = waves, data = dat_nested_na,
                                family = gaussian, corstr = "nested-exchangeable",
                                subgroup = sub)



    expect_identical(coef(fit_manual_nested), coef(fit_geefit_nested))
    expect_identical(fit_manual_nested$geese$alpha, fit_geefit_nested$geese$alpha)




    # ---------------------------------------------------------------------
    # pairwise-grouped-exchangeable: The answer must be close to 0.5 0.4 0.3 0.2 0.15 0.25
    # ---------------------------------------------------------------------

    dat_pge_clean <- dat_pge_na[complete.cases(dat_pge_na), ]
    dat_pge_clean <- dat_pge_clean[order(dat_pge_clean$id,
                                         dat_pge_clean$waves), ]

    zcor_manual <- build_zcor("pairwise-grouped-exchangeable", id = dat_pge_clean$id,
                              waves = dat_pge_clean$waves,
                              block = dplyr::distinct(dat_pge_clean, waves, blk)$blk)
    fit_manual_pge <- geeglm(y ~ 1, id = id, data = dat_pge_clean, family = gaussian,
                             corstr = "userdefined", zcor = zcor_manual)


    fit_geefit_pge <- geefit(y ~ 1, id = id, waves = waves, data = dat_pge_na,
                             family = gaussian, corstr = "pairwise-grouped-exchangeable",
                             block = blk)



    expect_identical(coef(fit_manual_pge), coef(fit_geefit_pge))
    expect_identical(fit_manual_pge$geese$alpha, fit_geefit_pge$geese$alpha)


    # ---------------------------------------------------------------------
    # block-exchangeable: The answer must be close to 0.3 0.15 0.25
    # ---------------------------------------------------------------------

    dat_li_clean <- dat_li_na[complete.cases(dat_li_na), ]
    dat_li_clean <- dat_li_clean[order(dat_li_clean$id,
                                       dat_li_clean$period), ]

    zcor_manual <- build_zcor("block-exchangeable", id = dat_li_clean$id,
                              waves = dat_li_clean$period,
                              individual = dat_li_clean$individual)
    fit_manual_li<- geeglm(y ~ 1, id = id, data = dat_li_clean, family = gaussian,
                           corstr = "userdefined", zcor = zcor_manual)


    fit_geefit_li <- geefit(y ~ 1, id = id, waves = period, data = dat_li_na,
                            family = gaussian, corstr = "block-exchangeable",
                            individual = individual)



    expect_identical(coef(fit_manual_li), coef(fit_geefit_li))
    expect_identical(fit_manual_li$geese$alpha, fit_geefit_li$geese$alpha)
  })

