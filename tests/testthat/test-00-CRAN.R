context("(0) Basic tests for CRAN")
# These are just a very rudimentary tests that don't require time or saved files.
# The full range of tests is run locally.


test_that("Basic functionality works", {

  RoBMA:::.load_RoBMA_module()

  # fit a default model
  d <- c(0.1, 0.2, 0.3)
  n <- c(20, 15, 10)
  fit_default <- suppressWarnings(RoBMA(d = d, n = n, chains = 1, burnin = 100, sample = 200, autofit = FALSE, seed = 1, algorithm = "ss"))

  expect_equal(TRUE, is.RoBMA(fit_default))

  expect_equal(
    capture_output_lines(fit_default, print = TRUE, width = 150),
    c("Call:"                                                                                                                                           ,
      "RoBMA(d = d, n = n, algorithm = \"ss\", chains = 1, sample = 200, "                                                                              ,
      "    burnin = 100, autofit = FALSE, seed = 1)"                                                                                                    ,
      ""                                                                                                                                                ,
      "Estimates:"                                                                                                                                      ,
      "               mu               tau    omega[0,0.025] omega[0.025,0.05]   omega[0.05,0.5]   omega[0.5,0.95] omega[0.95,0.975]    omega[0.975,1] ",
      "     -0.005022785       0.093402625       1.000000000       0.947454681       0.853829443       0.708481230       0.708800846       0.720853861 ",
      "              PET             PEESE "                                                                                                            ,
      "      0.025472293       0.027265737 "                                                                                                            )
  )

  expect_equal(
    capture_output_lines(summary(fit_default), print = TRUE, width = 150),
    c("Call:"                                                                                              ,
      "RoBMA(d = d, n = n, algorithm = \"ss\", chains = 1, sample = 200, "                                 ,
      "    burnin = 100, autofit = FALSE, seed = 1)"                                                       ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Components summary:"                                                                                ,
      "              Prior prob. Post. prob. Inclusion BF"                                                 ,
      "Effect              0.500       0.285        0.399"                                                 ,
      "Heterogeneity       0.500       0.490        0.961"                                                 ,
      "Bias                0.500       0.465        0.869"                                                 ,
      ""                                                                                                   ,
      "Model-averaged estimates:"                                                                          ,
      "                    Mean Median  0.025 0.975"                                                       ,
      "mu                -0.005  0.000 -0.586 0.474"                                                       ,
      "tau                0.093  0.000  0.000 0.451"                                                       ,
      "omega[0,0.025]     1.000  1.000  1.000 1.000"                                                       ,
      "omega[0.025,0.05]  0.947  1.000  0.534 1.000"                                                       ,
      "omega[0.05,0.5]    0.854  1.000  0.269 1.000"                                                       ,
      "omega[0.5,0.95]    0.708  1.000  0.017 1.000"                                                       ,
      "omega[0.95,0.975]  0.709  1.000  0.017 1.000"                                                       ,
      "omega[0.975,1]     0.721  1.000  0.017 1.000"                                                       ,
      "PET                0.025  0.000  0.000 0.379"                                                       ,
      "PEESE              0.027  0.000  0.000 0.016"                                                       ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
      "(Estimated publication weights omega correspond to one-sided p-values.)"                            ,
      "\033[0;31mESS 122 is lower than the set target (500).\033[0m"                                       )
  )


  df_reg <- data.frame(
    d       = d,
    se      = se_d(d, n),
    mod_con = scale(1:3),
    mod_cat = c("A", "A", "B")
  )

  fit_reg <- suppressWarnings(RoBMA.reg(~ mod_cat + mod_con, data = df_reg, priors_bias = NULL, chains = 1, burnin = 100, sample = 500, autofit = FALSE, seed = 1, algorithm = "ss"))

  expect_equal(
    capture_output_lines(summary(fit_reg), print = TRUE, width = 150),
    c("Call:"                                                                                              ,
      "RoBMA.reg(formula = ~mod_cat + mod_con, data = df_reg, priors_bias = NULL, "                        ,
      "    algorithm = \"ss\", chains = 1, sample = 500, burnin = 100, "                                   ,
      "    autofit = FALSE, seed = 1)"                                                                     ,
      ""                                                                                                   ,
      "Robust Bayesian meta-regression"                                                                    ,
      "Components summary:"                                                                                ,
      "              Prior prob. Post. prob. Inclusion BF"                                                 ,
      "Effect              0.500       0.282        0.393"                                                 ,
      "Heterogeneity       0.500       0.402        0.672"                                                 ,
      ""                                                                                                   ,
      "Meta-regression components summary:"                                                                ,
      "        Prior prob. Post. prob. Inclusion BF"                                                       ,
      "mod_cat       0.500       0.472        0.894"                                                       ,
      "mod_con       0.500       0.510        1.041"                                                       ,
      ""                                                                                                   ,
      "Model-averaged estimates:"                                                                          ,
      "     Mean Median  0.025 0.975"                                                                      ,
      "mu  0.052  0.000 -0.268 0.677"                                                                      ,
      "tau 0.090  0.000  0.000 0.476"                                                                      ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
      "(Estimated publication weights omega correspond to one-sided p-values.)"                            ,
      "\033[0;31mESS 193 is lower than the set target (500).\033[0m"                                       ,
      ""                                                                                                   ,
      "Model-averaged meta-regression estimates:"                                                          ,
      "                   Mean Median  0.025 0.975"                                                        ,
      "intercept         0.052  0.000 -0.268 0.677"                                                        ,
      "mod_cat [dif: A] -0.024  0.000 -0.356 0.284"                                                        ,
      "mod_cat [dif: B]  0.024  0.000 -0.284 0.356"                                                        ,
      "mod_con           0.009  0.000 -0.327 0.364"                                                        ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
      "\033[0;31mESS 193 is lower than the set target (500).\033[0m"
  ))

})
