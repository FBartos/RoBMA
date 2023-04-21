context("(3) Model setup")
skip_on_cran()

# test model preview
test_that("Model preview works", {

  expect_equal(
    capture_output_lines(check_setup(models = FALSE), print = TRUE, width = 150),
    c("Robust Bayesian meta-analysis (set-up)"                                                                        ,
      "Components summary:"                                                                                           ,
      "              Models Prior prob."                                                                              ,
      "Effect         18/36       0.500"                                                                              ,
      "Heterogeneity  18/36       0.500"                                                                              ,
      "Bias           32/36       0.500"                                                                              )
  )

  expect_equal(
    capture_output_lines(check_setup(models = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian meta-analysis (set-up)"                                                                        ,
      "Components summary:"                                                                                           ,
      "              Models Prior prob."                                                                              ,
      "Effect         18/36       0.500"                                                                              ,
      "Heterogeneity  18/36       0.500"                                                                              ,
      "Bias           32/36       0.500"                                                                              ,
      ""                                                                                                              ,
      "Models overview:"                                                                                              ,
      " Model Prior Effect Prior Heterogeneity                         Prior Bias                         Prior prob.",
      "     1     Spike(0)            Spike(0)                                                                  0.125",
      "     2     Spike(0)            Spike(0)           omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "     3     Spike(0)            Spike(0)       omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "     4     Spike(0)            Spike(0)           omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "     5     Spike(0)            Spike(0)     omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010",
      "     6     Spike(0)            Spike(0)       omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "     7     Spike(0)            Spike(0) omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010",
      "     8     Spike(0)            Spike(0)                             PET ~ Cauchy(0, 1)[0, Inf]           0.031",
      "     9     Spike(0)            Spike(0)                           PEESE ~ Cauchy(0, 5)[0, Inf]           0.031",
      "    10     Spike(0)   InvGamma(1, 0.15)                                                                  0.125",
      "    11     Spike(0)   InvGamma(1, 0.15)           omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    12     Spike(0)   InvGamma(1, 0.15)       omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    13     Spike(0)   InvGamma(1, 0.15)           omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    14     Spike(0)   InvGamma(1, 0.15)     omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010",
      "    15     Spike(0)   InvGamma(1, 0.15)       omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    16     Spike(0)   InvGamma(1, 0.15) omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010",
      "    17     Spike(0)   InvGamma(1, 0.15)                             PET ~ Cauchy(0, 1)[0, Inf]           0.031",
      "    18     Spike(0)   InvGamma(1, 0.15)                           PEESE ~ Cauchy(0, 5)[0, Inf]           0.031",
      "    19 Normal(0, 1)            Spike(0)                                                                  0.125",
      "    20 Normal(0, 1)            Spike(0)           omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    21 Normal(0, 1)            Spike(0)       omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    22 Normal(0, 1)            Spike(0)           omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    23 Normal(0, 1)            Spike(0)     omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010",
      "    24 Normal(0, 1)            Spike(0)       omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    25 Normal(0, 1)            Spike(0) omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010",
      "    26 Normal(0, 1)            Spike(0)                             PET ~ Cauchy(0, 1)[0, Inf]           0.031",
      "    27 Normal(0, 1)            Spike(0)                           PEESE ~ Cauchy(0, 5)[0, Inf]           0.031",
      "    28 Normal(0, 1)   InvGamma(1, 0.15)                                                                  0.125",
      "    29 Normal(0, 1)   InvGamma(1, 0.15)           omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    30 Normal(0, 1)   InvGamma(1, 0.15)       omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    31 Normal(0, 1)   InvGamma(1, 0.15)           omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    32 Normal(0, 1)   InvGamma(1, 0.15)     omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010",
      "    33 Normal(0, 1)   InvGamma(1, 0.15)       omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    34 Normal(0, 1)   InvGamma(1, 0.15) omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010",
      "    35 Normal(0, 1)   InvGamma(1, 0.15)                             PET ~ Cauchy(0, 1)[0, Inf]           0.031",
      "    36 Normal(0, 1)   InvGamma(1, 0.15)                           PEESE ~ Cauchy(0, 5)[0, Inf]           0.031")
  )

  expect_equal(
    capture_output_lines(check_setup(model_type = "PSMA", models = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian meta-analysis (set-up)"                                                                        ,
      "Components summary:"                                                                                           ,
      "              Models Prior prob."                                                                              ,
      "Effect         18/36       0.500"                                                                              ,
      "Heterogeneity  18/36       0.500"                                                                              ,
      "Bias           32/36       0.500"                                                                              ,
      ""                                                                                                              ,
      "Models overview:"                                                                                              ,
      " Model Prior Effect Prior Heterogeneity                         Prior Bias                         Prior prob.",
      "     1     Spike(0)            Spike(0)                                                                  0.125",
      "     2     Spike(0)            Spike(0)           omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "     3     Spike(0)            Spike(0)       omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "     4     Spike(0)            Spike(0)           omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "     5     Spike(0)            Spike(0)     omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010",
      "     6     Spike(0)            Spike(0)       omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "     7     Spike(0)            Spike(0) omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010",
      "     8     Spike(0)            Spike(0)                             PET ~ Cauchy(0, 1)[0, Inf]           0.031",
      "     9     Spike(0)            Spike(0)                           PEESE ~ Cauchy(0, 5)[0, Inf]           0.031",
      "    10     Spike(0)   InvGamma(1, 0.15)                                                                  0.125",
      "    11     Spike(0)   InvGamma(1, 0.15)           omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    12     Spike(0)   InvGamma(1, 0.15)       omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    13     Spike(0)   InvGamma(1, 0.15)           omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    14     Spike(0)   InvGamma(1, 0.15)     omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010",
      "    15     Spike(0)   InvGamma(1, 0.15)       omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    16     Spike(0)   InvGamma(1, 0.15) omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010",
      "    17     Spike(0)   InvGamma(1, 0.15)                             PET ~ Cauchy(0, 1)[0, Inf]           0.031",
      "    18     Spike(0)   InvGamma(1, 0.15)                           PEESE ~ Cauchy(0, 5)[0, Inf]           0.031",
      "    19 Normal(0, 1)            Spike(0)                                                                  0.125",
      "    20 Normal(0, 1)            Spike(0)           omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    21 Normal(0, 1)            Spike(0)       omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    22 Normal(0, 1)            Spike(0)           omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    23 Normal(0, 1)            Spike(0)     omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010",
      "    24 Normal(0, 1)            Spike(0)       omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    25 Normal(0, 1)            Spike(0) omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010",
      "    26 Normal(0, 1)            Spike(0)                             PET ~ Cauchy(0, 1)[0, Inf]           0.031",
      "    27 Normal(0, 1)            Spike(0)                           PEESE ~ Cauchy(0, 5)[0, Inf]           0.031",
      "    28 Normal(0, 1)   InvGamma(1, 0.15)                                                                  0.125",
      "    29 Normal(0, 1)   InvGamma(1, 0.15)           omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    30 Normal(0, 1)   InvGamma(1, 0.15)       omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    31 Normal(0, 1)   InvGamma(1, 0.15)           omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010",
      "    32 Normal(0, 1)   InvGamma(1, 0.15)     omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010",
      "    33 Normal(0, 1)   InvGamma(1, 0.15)       omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010",
      "    34 Normal(0, 1)   InvGamma(1, 0.15) omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010",
      "    35 Normal(0, 1)   InvGamma(1, 0.15)                             PET ~ Cauchy(0, 1)[0, Inf]           0.031",
      "    36 Normal(0, 1)   InvGamma(1, 0.15)                           PEESE ~ Cauchy(0, 5)[0, Inf]           0.031")
  )

  expect_equal(
    capture_output_lines(check_setup(model_type = "PP", models = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian meta-analysis (set-up)"                                         ,
      "Components summary:"                                                             ,
      "              Models Prior prob."                                                ,
      "Effect          6/12       0.500"                                                ,
      "Heterogeneity   6/12       0.500"                                                ,
      "Bias            8/12       0.500"                                                ,
      ""                                                                                ,
      "Models overview:"                                                                ,
      " Model Prior Effect Prior Heterogeneity          Prior Bias          Prior prob.",
      "     1     Spike(0)            Spike(0)                                    0.125",
      "     2     Spike(0)            Spike(0)   PET ~ Cauchy(0, 1)[0, Inf]       0.062",
      "     3     Spike(0)            Spike(0) PEESE ~ Cauchy(0, 5)[0, Inf]       0.062",
      "     4     Spike(0)   InvGamma(1, 0.15)                                    0.125",
      "     5     Spike(0)   InvGamma(1, 0.15)   PET ~ Cauchy(0, 1)[0, Inf]       0.062",
      "     6     Spike(0)   InvGamma(1, 0.15) PEESE ~ Cauchy(0, 5)[0, Inf]       0.062",
      "     7 Normal(0, 1)            Spike(0)                                    0.125",
      "     8 Normal(0, 1)            Spike(0)   PET ~ Cauchy(0, 1)[0, Inf]       0.062",
      "     9 Normal(0, 1)            Spike(0) PEESE ~ Cauchy(0, 5)[0, Inf]       0.062",
      "    10 Normal(0, 1)   InvGamma(1, 0.15)                                    0.125",
      "    11 Normal(0, 1)   InvGamma(1, 0.15)   PET ~ Cauchy(0, 1)[0, Inf]       0.062",
      "    12 Normal(0, 1)   InvGamma(1, 0.15) PEESE ~ Cauchy(0, 5)[0, Inf]       0.062")
  )

  expect_equal(
    capture_output_lines(check_setup(model_type = "2w", models = TRUE), print = TRUE, width = 150),
    c("Robust Bayesian meta-analysis (set-up)"                                                                ,
      "Components summary:"                                                                                   ,
      "              Models Prior prob."                                                                      ,
      "Effect          6/12       0.500"                                                                      ,
      "Heterogeneity   6/12       0.500"                                                                      ,
      "Bias            8/12       0.500"                                                                      ,
      ""                                                                                                      ,
      "Models overview:"                                                                                      ,
      " Model Prior Effect Prior Heterogeneity                     Prior Bias                     Prior prob.",
      "     1     Spike(0)            Spike(0)                                                          0.125",
      "     2     Spike(0)            Spike(0)      omega[two-sided: .05] ~ CumDirichlet(1, 1)          0.062",
      "     3     Spike(0)            Spike(0)  omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)       0.062",
      "     4     Spike(0)   InvGamma(1, 0.15)                                                          0.125",
      "     5     Spike(0)   InvGamma(1, 0.15)      omega[two-sided: .05] ~ CumDirichlet(1, 1)          0.062",
      "     6     Spike(0)   InvGamma(1, 0.15)  omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)       0.062",
      "     7 Normal(0, 1)            Spike(0)                                                          0.125",
      "     8 Normal(0, 1)            Spike(0)      omega[two-sided: .05] ~ CumDirichlet(1, 1)          0.062",
      "     9 Normal(0, 1)            Spike(0)  omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)       0.062",
      "    10 Normal(0, 1)   InvGamma(1, 0.15)                                                          0.125",
      "    11 Normal(0, 1)   InvGamma(1, 0.15)      omega[two-sided: .05] ~ CumDirichlet(1, 1)          0.062",
      "    12 Normal(0, 1)   InvGamma(1, 0.15)  omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)       0.062")
  )
})

test_that("RoBMA.reg model preview works", {

  df_reg <- data.frame(
    d       = c(rep(-1, 5), rep(0, 5), rep(1, 5)),
    se      = rep(0.1, 15),
    mod_cat = c(rep("A", 5), rep("B", 5), rep("C", 5)),
    mod_con = c((1:15)/15)
  )

  expect_equal(
    capture_output_lines(check_setup.reg(~ mod_cat + mod_con, data = df_reg, models = FALSE), print = TRUE, width = 150),
    c("Robust Bayesian meta-regression (set-up)",
      "Components summary:"                     ,
      "               Models Prior prob."       ,
      "Effect         72/144       0.500"       ,
      "Heterogeneity  72/144       0.500"       ,
      "Bias          128/144       0.500"       ,
      ""                                        ,
      "Meta-regression components summary:"     ,
      "        Models Prior prob."              ,
      "mod_cat 72/144       0.500"              ,
      "mod_con 72/144       0.500"
    )
  )

  expect_equal(
    capture_output_lines(check_setup.reg(~ 1, data = df_reg, models = FALSE), print = TRUE, width = 150),
    c("Robust Bayesian meta-regression (set-up)",
      "Components summary:"                     ,
      "              Models Prior prob."        ,
      "Effect         18/36       0.500"        ,
      "Heterogeneity  18/36       0.500"        ,
      "Bias           32/36       0.500"        ,
      ""                                        ,
      "Meta-regression components summary:"     ,
      "[1] Models      Prior prob."             ,
      "<0 rows> (or 0-length row.names)"
  ))

  expect_equal(
    capture_output_lines(check_setup.reg(~ mod_cat + mod_con, data = df_reg, models = TRUE, priors_bias = NULL, priors_heterogeneity = NULL), print = TRUE, width = 150),
    c("Robust Bayesian meta-regression (set-up)"                                                                                  ,
      "Components summary:"                                                                                                       ,
      "              Models Prior prob."                                                                                          ,
      "Effect           4/8       0.500"                                                                                          ,
      "Heterogeneity    0/8       0.000"                                                                                          ,
      "Bias             0/8       0.000"                                                                                          ,
      ""                                                                                                                          ,
      "Meta-regression components summary:"                                                                                       ,
      "        Models Prior prob."                                                                                                ,
      "mod_cat    4/8       0.500"                                                                                                ,
      "mod_con    4/8       0.500"                                                                                                ,
      ""                                                                                                                          ,
      "Models overview:"                                                                                                          ,
      " Model Prior mod_cat                Prior mod_cat                Prior mod_cat  Prior Heterogeneity Prior Bias Prior prob.",
      "     1      Spike(0)                                   Spike(0)        Spike(0)            Spike(0)                  0.125",
      "     2  Normal(0, 1)                                   Spike(0)        Spike(0)            Spike(0)                  0.125",
      "     3      Spike(0) mean difference contrast: mNormal(0, 0.25)        Spike(0)            Spike(0)                  0.125",
      "     4  Normal(0, 1) mean difference contrast: mNormal(0, 0.25)        Spike(0)            Spike(0)                  0.125",
      "     5      Spike(0)                                   Spike(0) Normal(0, 0.25)            Spike(0)                  0.125",
      "     6  Normal(0, 1)                                   Spike(0) Normal(0, 0.25)            Spike(0)                  0.125",
      "     7      Spike(0) mean difference contrast: mNormal(0, 0.25) Normal(0, 0.25)            Spike(0)                  0.125",
      "     8  Normal(0, 1) mean difference contrast: mNormal(0, 0.25) Normal(0, 0.25)            Spike(0)                  0.125"
    )
  )

  expect_equal(
    capture_output_lines(check_setup.reg(~ mod_cat + mod_con, data = df_reg, models = TRUE,
                                         priors_bias = NULL, priors_heterogeneity = NULL,
                                         priors = list(
                                           "mod_cat" = prior_factor("beta", list(1, 1), contrast = "treatment", prior_weights = 1/2),
                                           "mod_con" = list(
                                             "null" = prior("normal", list(0,    0.01)),
                                             "alt"  = prior("normal", list(0.10, 0.30))
                                           )
                                         )), print = TRUE, width = 150),
    c("Robust Bayesian meta-regression (set-up)"                                                                        ,
      "Components summary:"                                                                                             ,
      "              Models Prior prob."                                                                                ,
      "Effect           4/8       0.500"                                                                                ,
      "Heterogeneity    0/8       0.000"                                                                                ,
      "Bias             0/8       0.000"                                                                                ,
      ""                                                                                                                ,
      "Meta-regression components summary:"                                                                             ,
      "        Models Prior prob."                                                                                      ,
      "mod_cat    4/8       0.667"                                                                                      ,
      "mod_con    4/8       0.500"                                                                                      ,
      ""                                                                                                                ,
      "Models overview:"                                                                                                ,
      " Model Prior mod_cat          Prior mod_cat           Prior mod_cat   Prior Heterogeneity Prior Bias Prior prob.",
      "     1      Spike(0)                       Spike(0)   Normal(0, 0.01)            Spike(0)                  0.167",
      "     2  Normal(0, 1)                       Spike(0)   Normal(0, 0.01)            Spike(0)                  0.167",
      "     3      Spike(0) treatment contrast: Beta(1, 1)   Normal(0, 0.01)            Spike(0)                  0.083",
      "     4  Normal(0, 1) treatment contrast: Beta(1, 1)   Normal(0, 0.01)            Spike(0)                  0.083",
      "     5      Spike(0)                       Spike(0)  Normal(0.1, 0.3)            Spike(0)                  0.167",
      "     6  Normal(0, 1)                       Spike(0)  Normal(0.1, 0.3)            Spike(0)                  0.167",
      "     7      Spike(0) treatment contrast: Beta(1, 1)  Normal(0.1, 0.3)            Spike(0)                  0.083",
      "     8  Normal(0, 1) treatment contrast: Beta(1, 1)  Normal(0.1, 0.3)            Spike(0)                  0.083"
    )
  )
})

test_that("Set autofit control works", {

  expect_error(set_autofit_control(max_Rhat = .99), "Checking 'autofit_control':\n\tThe 'max_Rhat' must be equal or higher than 1.")
  expect_error(set_autofit_control(min_ESS  =  -1), "Checking 'autofit_control':\n\tThe 'min_ESS' must be equal or higher than 0.")
  expect_error(set_autofit_control(max_error=  -1), "Checking 'autofit_control':\n\tThe 'max_error' must be equal or higher than 0.")
  expect_error(set_autofit_control(max_SD_error=  -.1), "Checking 'autofit_control':\n\tThe 'max_SD_error' must be equal or higher than 0.")
  expect_error(set_autofit_control(max_SD_error=  1.1), "Checking 'autofit_control':\n\tThe 'max_SD_error' must be equal or lower than 1.")
  expect_error(set_autofit_control(max_time = list(time = -1, unit = "secs")), "Checking 'autofit_control':\n\tThe 'max_time:time' must be equal or higher than 0.")
  expect_error(set_autofit_control(max_time = list(time = 10, unit = "maps")), "Checking 'autofit_control':\n\tThe 'maps' values are not recognized by the 'max_time:unit' argument.")

  expect_equal(set_autofit_control(), list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(max_Rhat = 1.01),  list(
    max_Rhat      = 1.01,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(min_ESS = 200),  list(
    max_Rhat      = 1.05,
    min_ESS       = 200,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(max_error = 0.01),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = 0.01,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(max_SD_error = 0.01),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = 0.01,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(max_time = list(time = 30, unit = "secs")),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 30, unit = "secs"),
    sample_extend = 1000
  ))

  expect_equal(set_autofit_control(sample_extend = 200),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    max_time      = list(time = 60, unit = "mins"),
    sample_extend = 200
  ))

})

test_that("Set convergence checks works", {


  expect_error(set_convergence_checks(max_Rhat = .99), "Checking 'convergence_checks':\n\tThe 'max_Rhat' must be equal or higher than 1.")
  expect_error(set_convergence_checks(min_ESS  =  -1), "Checking 'convergence_checks':\n\tThe 'min_ESS' must be equal or higher than 0.")
  expect_error(set_convergence_checks(max_error=  -1), "Checking 'convergence_checks':\n\tThe 'max_error' must be equal or higher than 0.")
  expect_error(set_convergence_checks(max_SD_error=  -.1), "Checking 'convergence_checks':\n\tThe 'max_SD_error' must be equal or higher than 0.")
  expect_error(set_convergence_checks(max_SD_error=  1.1), "Checking 'convergence_checks':\n\tThe 'max_SD_error' must be equal or lower than 1.")


  expect_equal(set_convergence_checks(), list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(max_Rhat = 1.01),  list(
    max_Rhat      = 1.01,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(min_ESS = 200),  list(
    max_Rhat      = 1.05,
    min_ESS       = 200,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(max_error = 0.01),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = 0.01,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(max_SD_error = 0.01),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = 0.01,
    remove_failed       = FALSE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(remove_failed = TRUE),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = TRUE,
    balance_probability = TRUE
  ))

  expect_equal(set_convergence_checks(balance_probability = FALSE),  list(
    max_Rhat      = 1.05,
    min_ESS       = 500,
    max_error     = NULL,
    max_SD_error  = NULL,
    remove_failed       = FALSE,
    balance_probability = FALSE
  ))

})

