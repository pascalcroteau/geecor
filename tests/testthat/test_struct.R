









test_that(
  "testing structures",
  {

    # ---------------------------------------------------------------------
    # Toeplitz - example from geepack
    # ---------------------------------------------------------------------

    set.seed(88)
    dat <- gendat()

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






    ref <- geeglm(wheeze ~ city + age + smoke, id = case, data = six,
                  waves = age, family = binomial, corstr = "exchangeable")

    fite <- geefit(wheeze ~ city + age + smoke, id = case, data = six,
                   waves = age, family = binomial, corstr = "exchangeable")

    fite_smp <- geefit(wheeze ~ city + age + smoke, id = case,
                       data = six |> dplyr::slice_sample(prop = 1),
                       waves = age, family = binomial, corstr = "exchangeable")

    fitm <- geefit(wheeze ~ city + age + smoke, id = case, data = six,
                   waves = age, family = binomial, corstr = "banded-exchangeable",
                   bandwidth = 3)

    fitm_smp <- geefit(wheeze ~ city + age + smoke, id = case,
                   data = six |> dplyr::slice_sample(prop = 1), waves = age,
                   family = binomial, corstr = "banded-exchangeable", bandwidth = 3)



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

    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_band, family = gaussian,
             corstr = "banded-toeplitz", waves = waves, bandwidth = 2)
      )


    # ---------------------------------------------------------------------
    # banded-unstructured. The answer must be close to true_pairs_bu
    # ---------------------------------------------------------------------

    # true_pairs_bu <- list(
    #   "1_2" = 0.50,  # lag 1
    #   "1_3" = 0.30,  # lag 2
    #   "2_3" = 0.45,  # lag 1
    #   "2_4" = 0.25,  # lag 2
    #   "3_4" = 0.40,  # lag 1
    #   "3_5" = 0.20,  # lag 2
    #   "4_5" = 0.35   # lag 1
    # )
    #

    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_bu, family = gaussian,
             corstr = "banded-unstructured", waves = waves, bandwidth = k_bu)
    )



    # ---------------------------------------------------------------------
    # banded-exchangeable The answer must be close to 0.4
    # ---------------------------------------------------------------------

    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_bdex, family = gaussian,
             corstr = "banded-exchangeable", bandwidth = 2, waves = waves)
    )



    # ---------------------------------------------------------------------
    # banded-exchangeable, example from geeM.
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


    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_nested, family = gaussian,
             corstr = "nested-exchangeable", waves = waves, subgroup = sub)
    )



    # ---------------------------------------------------------------------
    # pairwise-grouped-exchangeable: The answer must be close to 0.5 0.4 0.3 0.2 0.15 0.25
    # ---------------------------------------------------------------------


    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_pge, family = gaussian,
             corstr = "pairwise-grouped-exchangeable", waves = waves,
             block = blk)
    )



    # ---------------------------------------------------------------------
    # block-exchangeable: The answer must be close to 0.3 0.15 0.25
    # ---------------------------------------------------------------------

    expect_snapshot(
      geefit(y ~ 1, id = id, data = dat_li, family = gaussian,
             corstr = "block-exchangeable", waves = period,
             individual = individual)
    )



  })
