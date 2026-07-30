test_that(
  "testing other functions",
  {


    # ---------------------------------------------------------------------
    # testing 'update()'
    # ---------------------------------------------------------------------


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



    fit_full <- geefit(wheeze ~ city + age + smoke, id = case, data = six,
                       waves = age, family = binomial, corstr = "m-dependent",
                       mdep = 3)

    fit_upd <- update(fit_full, . ~ . - age)

    fit_ref <- geefit(wheeze ~ city + smoke, id = case, data = six,
                      waves = age, family = binomial, corstr = "m-dependent",
                      mdep = 3)
    expect_identical(coef(fit_upd), coef(fit_ref))
    expect_identical(fit_upd$geese$alpha, fit_ref$geese$alpha)
    expect_identical(fit_upd$geese$gamma, fit_ref$geese$gamma)




    # ---------------------------------------------------------------------
    # testing 'getCall()'
    # ---------------------------------------------------------------------
    expect_identical(getCall(fit_full),
                     quote(geefit(formula = wheeze ~ city + age + smoke,
                                  data = six, id = case, waves = age,
                                  family = binomial, corstr = "m-dependent",
                                  mdep = 3)))


    # ---------------------------------------------------------------------
    # testing 'geecor_fit()'
    # ---------------------------------------------------------------------

    fit_call <- fit_full$call
    fit_call[[1]] <- quote(geecor_fit)
    X <- model.matrix(fit_full)
    fit_call <- rlang::call_modify(fit_call, x = quote(X),
                                   y = quote(fit_full$y),
                                   id = quote(fit_full$id),
                                   corstr = quote(fit_full$.corstruct),
                                   waves = quote(fit_full$waves),
                                   family = quote(fit_full$family),
                                   weights = quote(fit_full$prior.weights),
                                   start = quote(fit_full$start),
                                   offset = quote(fit_full$offset),
                                   control = quote(fit_full$control),
                                   !!!fit_full$.corparams,
                                   formula = rlang::zap())
    fit_full_recons <- eval(fit_call)
    expect_identical(coef(fit_full), fit_full_recons$beta)
    expect_identical(fit_full$geese$alpha, fit_full_recons$alpha)


    # ---------------------------------------------------------------------
    # testing 'anova()'
    # ---------------------------------------------------------------------
    aov1 <- anova(geefit(wheeze ~ city*smoke + age, id = case, data = six,
                         waves = age, family = binomial, corstr = "m-dependent",
                         mdep = 3))

    citystats <- geefit(wheeze ~ city, id = case, data = six,
                     waves = age, family = binomial, corstr = "m-dependent",
                     mdep = 3)
    smokestats <- geefit(wheeze ~ city + smoke, id = case, data = six,
                       waves = age, family = binomial, corstr = "m-dependent",
                       mdep = 3)
    agestats <- geefit(wheeze ~ city + smoke + age, id = case, data = six,
                       waves = age, family = binomial, corstr = "m-dependent",
                       mdep = 3)
    allstats <- geefit(wheeze ~ city*smoke + age, id = case, data = six,
                       waves = age, family = binomial, corstr = "m-dependent",
                       mdep = 3)


    expect_equal(aov1["city", ]$X2,
                 summary(citystats)$coefficients["citykingston", ]$Wald)
    expect_equal(aov1["city", ]$`P(>|Chi|)`,
                 summary(citystats)$coefficients["citykingston", ]$`Pr(>|W|)`)

    expect_equal(aov1["smoke", ]$X2,
                 summary(smokestats)$coefficients["smoke", ]$Wald)
    expect_equal(aov1["smoke", ]$`P(>|Chi|)`,
                 summary(smokestats)$coefficients["smoke", ]$`Pr(>|W|)`)

    expect_equal(aov1["age", ]$X2,
                 summary(agestats)$coefficients["age", ]$Wald)
    expect_equal(aov1["age", ]$`P(>|Chi|)`,
                 summary(agestats)$coefficients["age", ]$`Pr(>|W|)`)

    expect_equal(aov1["city:smoke", ]$X2,
                 summary(allstats)$coefficients["citykingston:smoke", ]$Wald)
    expect_equal(aov1["city:smoke", ]$`P(>|Chi|)`,
                 summary(allstats)$coefficients["citykingston:smoke", ]$`Pr(>|W|)`)
  })

