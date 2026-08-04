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
       0.02239333 
      
      Degrees of Freedom: 2073 Total (i.e. Null);  2072 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.082758
      
      Correlation:  Structure = banded-toeplitz    Link = identity 
      Estimated Correlation Parameters:
        alpha:1   alpha:2 
      0.5116232 0.3084207 
      
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
      0.003683427 
      
      Degrees of Freedom: 8734 Total (i.e. Null);  8733 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 0.9896692
      
      Correlation:  Structure = banded-unstructured    Link = identity 
      Estimated Correlation Parameters:
       pair_1_2  pair_1_3  pair_2_3  pair_2_4  pair_3_4  pair_3_5  pair_4_5 
      0.5168216 0.3184708 0.4306707 0.2524204 0.3733138 0.2006915 0.3441165 
      
      Number of clusters:   2000   Maximum cluster size: 5 
      

---

    Code
      geefit(y ~ 1, id = id, data = dat_bdex, family = gaussian, corstr = "banded-exchangeable",
      bandwidth = 2, waves = waves)
    Output
      
      Call:
      geefit(formula = y ~ 1, data = dat_bdex, id = id, waves = waves, 
          family = gaussian, corstr = "banded-exchangeable", bandwidth = 2)
      
      Coefficients:
      (Intercept) 
       0.02682528 
      
      Degrees of Freedom: 2073 Total (i.e. Null);  2072 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.079079
      
      Correlation:  Structure = banded-exchangeable    Link = identity 
      Estimated Correlation Parameters:
        alpha:1 
      0.4103338 
      
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
      -1.92586991 -0.15764935  0.30538744  0.09943373 
      
      Degrees of Freedom: 2148 Total (i.e. Null);  2144 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.017311
      
      Correlation:  Structure = m-dependent    Link = identity 
      Estimated Correlation Parameters:
        alpha:1 
      0.3995063 
      
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
       0.03138783 
      
      Degrees of Freedom: 2073 Total (i.e. Null);  2072 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.086539
      
      Correlation:  Structure = nested-exchangeable    Link = identity 
      Estimated Correlation Parameters:
         alpha1    alpha2 
      0.5257682 0.2116866 
      
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
       0.01019327 
      
      Degrees of Freedom: 2598 Total (i.e. Null);  2597 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.053845
      
      Correlation:  Structure = pairwise-grouped-exchangeable    Link = identity 
      Estimated Correlation Parameters:
      intra_bloc1 intra_bloc2 intra_bloc3   inter_1_2   inter_1_3   inter_2_3 
        0.5219813   0.4100338   0.3120269   0.2186978   0.1697476   0.2589218 
      
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
      0.003248133 
      
      Degrees of Freedom: 2700 Total (i.e. Null);  2699 Residual
      
      Scale Link:                   identity
      Estimated Scale Parameters:  [1] 1.061016
      
      Correlation:  Structure = block-exchangeable    Link = identity 
      Estimated Correlation Parameters:
         alpha1    alpha2    alpha3 
      0.3181522 0.1696427 0.2511806 
      
      Number of clusters:   300   Maximum cluster size: 9 
      

