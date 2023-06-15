context("(0) Basic tests for CRAN")
# These are just a very rudimentary tests that don't require time or saved files.
# The full range of tests is run locally.


test_that("Basic functionality works", {

  RoBMA:::.load_RoBMA_module()

  # fit a default model
  d <- c(0.1, 0.2, 0.3)
  n <- c(20, 15, 10)
  fit_default <- suppressWarnings(RoBMA(d = d, n = n, chains = 1, burnin = 50, sample = 100, autofit = FALSE, seed = 1))

  expect_equal(TRUE, is.RoBMA(fit_default))

  expect_equal(
    capture_output_lines(fit_default, print = TRUE, width = 150),
    c("Call:"                                                                                                                                           ,
      "RoBMA(d = d, n = n, chains = 1, sample = 100, burnin = 50, autofit = FALSE, "                                                                    ,
      "    seed = 1)"                                                                                                                                   ,
      ""                                                                                                                                                ,
      "Estimates:"                                                                                                                                      ,
      "               mu               tau    omega[0,0.025] omega[0.025,0.05]   omega[0.05,0.5]   omega[0.5,0.95] omega[0.95,0.975]    omega[0.975,1] ",
      "     -0.003244705       0.088661462       1.000000000       0.946723035       0.866019085       0.760510379       0.765383672       0.781386929 ",
      "              PET             PEESE "                                                                                                            ,
      "      0.044580252       0.045223594 " )
  )

  expect_equal(
    capture_output_lines(summary(fit_default), print = TRUE, width = 150),
    c("Call:"                                                                                        ,
      "RoBMA(d = d, n = n, chains = 1, sample = 100, burnin = 50, autofit = FALSE, "                 ,
      "    seed = 1)"                                                                                ,
      ""                                                                                             ,
      "Robust Bayesian meta-analysis"                                                                ,
      "Components summary:"                                                                          ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                    ,
      "Effect         18/36       0.500       0.288        0.405"                                    ,
      "Heterogeneity  18/36       0.500       0.412        0.702"                                    ,
      "Bias           32/36       0.500       0.460        0.852"                                    ,
      ""                                                                                             ,
      "Model-averaged estimates:"                                                                    ,
      "                    Mean Median  0.025 0.975"                                                 ,
      "mu                -0.003  0.000 -0.613 0.606"                                                 ,
      "tau                0.089  0.000  0.000 0.514"                                                 ,
      "omega[0,0.025]     1.000  1.000  1.000 1.000"                                                 ,
      "omega[0.025,0.05]  0.947  1.000  0.470 1.000"                                                 ,
      "omega[0.05,0.5]    0.866  1.000  0.238 1.000"                                                 ,
      "omega[0.5,0.95]    0.761  1.000  0.017 1.000"                                                 ,
      "omega[0.95,0.975]  0.765  1.000  0.017 1.000"                                                 ,
      "omega[0.975,1]     0.781  1.000  0.017 1.000"                                                 ,
      "PET                0.045  0.000  0.000 0.641"                                                 ,
      "PEESE              0.045  0.000  0.000 0.614"                                                 ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
      "(Estimated publication weights omega correspond to one-sided p-values.)"                      ,
      "\033[0;31mModel (2): ESS 45 is lower than the set target (500).\033[0m"                       ,
      "\033[0;31mModel (3): ESS 39 is lower than the set target (500).\033[0m"                       ,
      "\033[0;31mModel (4): ESS 29 is lower than the set target (500).\033[0m"                       ,
      "\033[0;31mModel (5): ESS 36 is lower than the set target (500).\033[0m"                       ,
      "\033[0;31mModel (6): ESS 28 is lower than the set target (500).\033[0m"                       ,
      "\033[0;31mThere were another 29 warnings. To see all warnings call 'check_RoBMA(fit)'.\033[0m")
  )


  df_reg <- data.frame(
    d       = d,
    se      = se_d(d, n),
    mod_con = scale(1:3),
    mod_cat = c("A", "A", "B")
  )

  fit_reg <- suppressWarnings(RoBMA.reg(~ mod_cat + mod_con, data = df_reg, priors_bias = NULL, chains = 1, burnin = 50, sample = 100, autofit = FALSE, seed = 1))

  expect_equal(
    capture_output_lines(summary(fit_reg), print = TRUE, width = 150),
    c("Call:"                                                                                              ,
      "RoBMA.reg(formula = ~mod_cat + mod_con, data = df_reg, priors_bias = NULL, "                        ,
      "    chains = 1, sample = 100, burnin = 50, autofit = FALSE, seed = 1)"                              ,
      ""                                                                                                   ,
      "Robust Bayesian meta-regression"                                                                    ,
      "Components summary:"                                                                                ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                          ,
      "Effect          8/16       0.500       0.267        0.364"                                          ,
      "Heterogeneity   8/16       0.500       0.401        0.671"                                          ,
      ""                                                                                                   ,
      "Meta-regression components summary:"                                                                ,
      "        Models Prior prob. Post. prob. Inclusion BF"                                                ,
      "mod_cat   8/16       0.500       0.470        0.886"                                                ,
      "mod_con   8/16       0.500       0.446        0.805"                                                ,
      ""                                                                                                   ,
      "Model-averaged estimates:"                                                                          ,
      "     Mean Median  0.025 0.975"                                                                      ,
      "mu  0.048  0.000 -0.300 0.667"                                                                      ,
      "tau 0.094  0.000  0.000 0.581"                                                                      ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
      "\033[0;31mModel (2): ESS 44 is lower than the set target (500).\033[0m"                             ,
      "\033[0;31mModel (3): ESS 45 is lower than the set target (500).\033[0m"                             ,
      "\033[0;31mModel (4): ESS 58 is lower than the set target (500).\033[0m"                             ,
      "\033[0;31mModel (5): ESS 4 is lower than the set target (500).\033[0m"                              ,
      "\033[0;31mModel (6): ESS 6 is lower than the set target (500).\033[0m"                              ,
      "\033[0;31mThere were another 9 warnings. To see all warnings call 'check_RoBMA(fit)'.\033[0m"       ,
      ""                                                                                                   ,
      "Model-averaged meta-regression estimates:"                                                          ,
      "                  Mean Median  0.025 0.975"                                                         ,
      "intercept        0.048  0.000 -0.300 0.667"                                                         ,
      "mod_cat [dif: A] 0.000  0.000 -0.357 0.459"                                                         ,
      "mod_cat [dif: B] 0.000  0.000 -0.459 0.357"                                                         ,
      "mod_con          0.005  0.000 -0.341 0.343"                                                         ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
      "\033[0;31mModel (2): ESS 44 is lower than the set target (500).\033[0m"                             ,
      "\033[0;31mModel (3): ESS 45 is lower than the set target (500).\033[0m"                             ,
      "\033[0;31mModel (4): ESS 58 is lower than the set target (500).\033[0m"                             ,
      "\033[0;31mModel (5): ESS 4 is lower than the set target (500).\033[0m"                              ,
      "\033[0;31mModel (6): ESS 6 is lower than the set target (500).\033[0m"                              ,
      "\033[0;31mThere were another 9 warnings. To see all warnings call 'check_RoBMA(fit)'.\033[0m"
  ))

})
