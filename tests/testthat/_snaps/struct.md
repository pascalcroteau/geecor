# testing structures

    Code
      geefit(y ~ 1, id = id, data = dat_band, family = gaussian, corstr = "banded-toeplitz",
      waves = waves, bandwidth = 2)
    Output
      
      Call:
      geepack::geeglm(formula = y ~ 1, family = gaussian, data = data, 
          id = id, waves = NULL, zcor = zcor, corstr = "userdefined")
      
      Coefficients:
      (Intercept) 
           0.0224 
      
      Degrees of Freedom: 2073 Total (i.e. Null);  2072 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.08
      
      Correlation:  Structure = banded-toeplitz    Link = identity 
      Estimated Correlation Parameters:
      alpha:1 alpha:2 
        0.512   0.308 
      
      Number of clusters:   400   Maximum cluster size: 6 
      

---

    Code
      geefit(y ~ 1, id = id, data = dat_mdep, family = gaussian, corstr = "m-dependent",
      mdep = 2, waves = waves)
    Output
      
      Call:
      geepack::geeglm(formula = y ~ 1, family = gaussian, data = data, 
          id = id, waves = NULL, zcor = zcor, corstr = "userdefined")
      
      Coefficients:
      (Intercept) 
           0.0268 
      
      Degrees of Freedom: 2073 Total (i.e. Null);  2072 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.08
      
      Correlation:  Structure = m-dependent    Link = identity 
      Estimated Correlation Parameters:
      alpha:1 
         0.41 
      
      Number of clusters:   400   Maximum cluster size: 6 
      

---

    Code
      geefit(resp ~ age + smoke + age:smoke, id = id, data = ohio, family = binomial,
      waves = age, corstr = "m-dependent", mdep = 1)
    Output
      
      Call:
      geepack::geeglm(formula = resp ~ age + smoke + age:smoke, family = binomial, 
          data = data, id = id, waves = NULL, zcor = zcor, corstr = "userdefined")
      
      Coefficients:
      (Intercept)         age       smoke   age:smoke 
          -1.9259     -0.1576      0.3054      0.0994 
      
      Degrees of Freedom: 2148 Total (i.e. Null);  2144 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.02
      
      Correlation:  Structure = m-dependent    Link = identity 
      Estimated Correlation Parameters:
      alpha:1 
          0.4 
      
      Number of clusters:   537   Maximum cluster size: 4 
      

---

    Code
      geefit(y ~ 1, id = id, data = dat_nested, family = gaussian, corstr = "nested-exchangeable",
      waves = waves, subgroup = sub)
    Output
      
      Call:
      geepack::geeglm(formula = y ~ 1, family = gaussian, data = data, 
          id = id, waves = NULL, zcor = zcor, corstr = "userdefined")
      
      Coefficients:
      (Intercept) 
           0.0314 
      
      Degrees of Freedom: 2073 Total (i.e. Null);  2072 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.09
      
      Correlation:  Structure = nested-exchangeable    Link = identity 
      Estimated Correlation Parameters:
      alpha1 alpha2 
       0.526  0.212 
      
      Number of clusters:   400   Maximum cluster size: 6 
      

---

    Code
      geefit(y ~ 1, id = id, data = dat_pge, family = gaussian, corstr = "pairwise-grouped-exchangeable",
      waves = waves, block = blk)
    Output
      
      Call:
      geepack::geeglm(formula = y ~ 1, family = gaussian, data = data, 
          id = id, waves = NULL, zcor = zcor, corstr = "userdefined")
      
      Coefficients:
      (Intercept) 
           0.0102 
      
      Degrees of Freedom: 2598 Total (i.e. Null);  2597 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.05
      
      Correlation:  Structure = pairwise-grouped-exchangeable    Link = identity 
      Estimated Correlation Parameters:
      intra_bloc1 intra_bloc2 intra_bloc3   inter_1_2   inter_1_3   inter_2_3 
            0.522       0.410       0.312       0.219       0.170       0.259 
      
      Number of clusters:   500   Maximum cluster size: 6 
      

---

    Code
      geefit(y ~ 1, id = id, data = dat_li, family = gaussian, corstr = "block-exchangeable",
      waves = period, individual = individual)
    Output
      
      Call:
      geepack::geeglm(formula = y ~ 1, family = gaussian, data = data, 
          id = id, waves = NULL, zcor = zcor, corstr = "userdefined")
      
      Coefficients:
      (Intercept) 
          0.00325 
      
      Degrees of Freedom: 2700 Total (i.e. Null);  2699 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.06
      
      Correlation:  Structure = block-exchangeable    Link = identity 
      Estimated Correlation Parameters:
      alpha1 alpha2 alpha3 
       0.318  0.170  0.251 
      
      Number of clusters:   300   Maximum cluster size: 9 
      

