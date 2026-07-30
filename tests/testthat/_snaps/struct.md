# testing structures

    Code
      geefit(y ~ 1, id = id, data = dat_band, family = gaussian, corstr = "banded-toeplitz",
      waves = waves, bandwidth = 2)
    Output
      
      Call:
      geefit(formula = y ~ 1, data = dat_band, id = id, waves = waves, 
          family = gaussian, corstr = "banded-toeplitz", bandwidth = 2)
      
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
      geefit(y ~ 1, id = id, data = dat_bu, family = gaussian, corstr = "banded-unstructured",
      waves = waves, bandwidth = k_bu)
    Output
      
      Call:
      geefit(formula = y ~ 1, data = dat_bu, id = id, waves = waves, 
          family = gaussian, corstr = "banded-unstructured", bandwidth = k_bu)
      
      Coefficients:
      (Intercept) 
          0.00368 
      
      Degrees of Freedom: 8734 Total (i.e. Null);  8733 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 0.99
      
      Correlation:  Structure = banded-unstructured    Link = identity 
      Estimated Correlation Parameters:
      pair_1_2 pair_1_3 pair_2_3 pair_2_4 pair_3_4 pair_3_5 pair_4_5 
         0.517    0.318    0.431    0.252    0.373    0.201    0.344 
      
      Number of clusters:   2000   Maximum cluster size: 5 
      

---

    Code
      geefit(y ~ 1, id = id, data = dat_mdep, family = gaussian, corstr = "m-dependent",
      mdep = 2, waves = waves)
    Output
      
      Call:
      geefit(formula = y ~ 1, data = dat_mdep, id = id, waves = waves, 
          family = gaussian, corstr = "m-dependent", mdep = 2)
      
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
      geefit(formula = resp ~ age + smoke + age:smoke, data = ohio, 
          id = id, waves = age, family = binomial, corstr = "m-dependent", 
          mdep = 1)
      
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
      geefit(formula = y ~ 1, data = dat_nested, id = id, waves = waves, 
          family = gaussian, corstr = "nested-exchangeable", subgroup = sub)
      
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
      geefit(formula = y ~ 1, data = dat_pge, id = id, waves = waves, 
          family = gaussian, corstr = "pairwise-grouped-exchangeable", 
          block = blk)
      
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
      geefit(formula = y ~ 1, data = dat_li, id = id, waves = period, 
          family = gaussian, corstr = "block-exchangeable", individual = individual)
      
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
      

