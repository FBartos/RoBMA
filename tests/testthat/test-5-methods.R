context("(5) Print and summary functions")
skip_on_cran()

# the summary tables and print functions are imported from BayesTools and tested henceforth
# test objects - assuming that the fit function worked properly
saved_files <- paste0("fit_", 1:13, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}

test_that("Print functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(saved_fits[[i]], print = TRUE, width = 150),
      read.table(file = file.path("../results/print", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }
})

test_that("Summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]]), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # all options
  expect_equal(
    capture_output_lines(summary(saved_fits[[1]], conditional = TRUE, logBF = TRUE, BF01 = TRUE, probs = c(0.10, 0.50, .90)), print = TRUE, width = 150),
    c("Call:"                                                                  ,
      "RoBMA(d = d, se = d_se, parallel = TRUE, seed = 1)"                     ,
      ""                                                                       ,
      "Robust Bayesian meta-analysis"                                          ,
      "Components summary:"                                                    ,
      "              Models Prior prob. Post. prob. Inclusion log(1/BF)"       ,
      "Effect         18/36       0.500       0.493               0.026"       ,
      "Heterogeneity  18/36       0.500       0.462               0.153"       ,
      "Bias           32/36       0.500       0.540              -0.159"       ,
      ""                                                                       ,
      "Model-averaged estimates:"                                              ,
      "                   Mean Median   0.1   0.5   0.9"                       ,
      "mu                0.196  0.000 0.000 0.000 0.618"                       ,
      "tau               0.109  0.000 0.000 0.000 0.327"                       ,
      "omega[0,0.025]    1.000  1.000 1.000 1.000 1.000"                       ,
      "omega[0.025,0.05] 0.924  1.000 0.671 1.000 1.000"                       ,
      "omega[0.05,0.5]   0.824  1.000 0.308 1.000 1.000"                       ,
      "omega[0.5,0.95]   0.759  1.000 0.126 1.000 1.000"                       ,
      "omega[0.95,0.975] 0.770  1.000 0.131 1.000 1.000"                       ,
      "omega[0.975,1]    0.801  1.000 0.136 1.000 1.000"                       ,
      "PET               0.110  0.000 0.000 0.000 0.401"                       ,
      "PEESE             0.087  0.000 0.000 0.000 0.000"                       ,
      "(Estimated publication weights omega correspond to one-sided p-values.)",
      ""                                                                       ,
      "Conditional estimates:"                                                 ,
      "                   Mean Median   0.1   0.5   0.9"                       ,
      "mu                0.398  0.414 0.054 0.414 0.732"                       ,
      "tau               0.242  0.170 0.061 0.170 0.489"                       ,
      "omega[0,0.025]    1.000  1.000 1.000 1.000 1.000"                       ,
      "omega[0.025,0.05] 0.784  0.888 0.371 0.888 1.000"                       ,
      "omega[0.05,0.5]   0.499  0.487 0.156 0.487 0.864"                       ,
      "omega[0.5,0.95]   0.313  0.247 0.039 0.247 0.710"                       ,
      "omega[0.95,0.975] 0.341  0.268 0.039 0.268 0.777"                       ,
      "omega[0.975,1]    0.435  0.304 0.040 0.304 1.000"                       ,
      "PET               0.825  0.753 0.167 0.753 1.564"                       ,
      "PEESE             1.690  1.540 0.347 1.540 3.201"                       ,
      "(Estimated publication weights omega correspond to one-sided p-values.)")
  )

  # different effect size measure
  expect_equal(
    capture_output_lines(summary(saved_fits[[1]], output_scale = "fishers_z"), print = TRUE, width = 150),
    c("Call:"                                                                                               ,
      "RoBMA(d = d, se = d_se, parallel = TRUE, seed = 1)"                                                  ,
      ""                                                                                                    ,
      "Robust Bayesian meta-analysis"                                                                       ,
      "Components summary:"                                                                                 ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                           ,
      "Effect         18/36       0.500       0.493        0.974"                                           ,
      "Heterogeneity  18/36       0.500       0.462        0.858"                                           ,
      "Bias           32/36       0.500       0.540        1.173"                                           ,
      ""                                                                                                    ,
      "Model-averaged estimates:"                                                                           ,
      "                   Mean Median  0.025 0.975"                                                         ,
      "mu                0.097  0.000 -0.032 0.399"                                                         ,
      "tau               0.055  0.000  0.000 0.312"                                                         ,
      "omega[0,0.025]    1.000  1.000  1.000 1.000"                                                         ,
      "omega[0.025,0.05] 0.924  1.000  0.309 1.000"                                                         ,
      "omega[0.05,0.5]   0.824  1.000  0.126 1.000"                                                         ,
      "omega[0.5,0.95]   0.759  1.000  0.027 1.000"                                                         ,
      "omega[0.95,0.975] 0.770  1.000  0.027 1.000"                                                         ,
      "omega[0.975,1]    0.801  1.000  0.027 1.000"                                                         ,
      "PET               0.110  0.000  0.000 1.310"                                                         ,
      "PEESE             0.175  0.000  0.000 3.142"                                                         ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale).",
      "(Estimated publication weights omega correspond to one-sided p-values.)" )
  )

  # no conditional, yet requested
  expect_equal(
    capture_output_lines(summary(saved_fits[[8]], conditional = TRUE), print = TRUE, width = 150),
    c( "Call:"                                                                      ,
       "RoBMA(d = d, se = d_se, priors_effect = NULL, priors_heterogeneity = NULL, ",
       "    priors_bias = NULL, parallel = TRUE, seed = 1)"                         ,
       ""                                                                           ,
       "Robust Bayesian meta-analysis"                                              ,
       "Components summary:"                                                        ,
       "              Models Prior prob. Post. prob. Inclusion BF"                  ,
       "Effect           0/1       0.000       0.000        0.000"                  ,
       "Heterogeneity    0/1       0.000       0.000        0.000"                  ,
       "Bias             0/1       0.000       0.000        0.000"                  ,
       ""                                                                           ,
       "Model-averaged estimates:"                                                  ,
       "     Mean Median 0.025 0.975"                                               ,
       "mu  0.000  0.000 0.000 0.000"                                               ,
       "tau 0.000  0.000 0.000 0.000"                                               ,
       ""                                                                           ,
       "Conditional estimates:"                                                     ,
       "[1] Mean   Median 0.025  0.975 "                                            ,
       "<0 rows> (or 0-length row.names)")
  )
})

test_that("Models summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]], type = "models"), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary.models", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # test short names
  expect_equal(
    capture_output_lines(summary(saved_fits[[11]], type = "models", short_name = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                                                    ,
      "RoBMA(y = d, se = d_se, priors_bias = list(prior_weightfunction(\"two-sided\", "                                          ,
      "    list(c(0.1), c(1, 1))), prior_PET(\"normal\", list(0, 1))), "                                                         ,
      "    parallel = TRUE, seed = 1)"                                                                                           ,
      ""                                                                                                                         ,
      "Robust Bayesian meta-analysis"                                                                                            ,
      "Models overview:"                                                                                                         ,
      " Model Prior Effect Prior Heterogeneity            Prior Bias           Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1         S(0)                S(0)                                       0.083        -2.90       0.043        0.499",
      "     2         S(0)                S(0) omega[2s: .1] ~ CumD(1, 1)            0.083        -2.60       0.059        0.686",
      "     3         S(0)                S(0)           PET ~ N(0, 1)[0, Inf]       0.083        -1.46       0.185        2.492",
      "     4         S(0)         Ig(1, 0.15)                                       0.083        -2.66       0.056        0.647",
      "     5         S(0)         Ig(1, 0.15) omega[2s: .1] ~ CumD(1, 1)            0.083        -2.61       0.058        0.683",
      "     6         S(0)         Ig(1, 0.15)           PET ~ N(0, 1)[0, Inf]       0.083        -1.73       0.141        1.799",
      "     7      N(0, 1)                S(0)                                       0.083        -2.01       0.107        1.313",
      "     8      N(0, 1)                S(0) omega[2s: .1] ~ CumD(1, 1)            0.083        -2.28       0.081        0.969",
      "     9      N(0, 1)                S(0)           PET ~ N(0, 1)[0, Inf]       0.083        -2.25       0.083        0.999",
      "    10      N(0, 1)         Ig(1, 0.15)                                       0.083        -2.38       0.073        0.870",
      "    11      N(0, 1)         Ig(1, 0.15) omega[2s: .1] ~ CumD(1, 1)            0.083        -2.70       0.053        0.621",
      "    12      N(0, 1)         Ig(1, 0.15)           PET ~ N(0, 1)[0, Inf]       0.083        -2.56       0.061        0.716"
    ))

  # test no spikes
  expect_equal(
    capture_output_lines(summary(saved_fits[[11]], type = "models", remove_spike_0 = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                                                                ,
      "RoBMA(y = d, se = d_se, priors_bias = list(prior_weightfunction(\"two-sided\", "                                                      ,
      "    list(c(0.1), c(1, 1))), prior_PET(\"normal\", list(0, 1))), "                                                                     ,
      "    parallel = TRUE, seed = 1)"                                                                                                       ,
      ""                                                                                                                                     ,
      "Robust Bayesian meta-analysis"                                                                                                        ,
      "Models overview:"                                                                                                                     ,
      " Model Prior Effect Prior Heterogeneity                  Prior Bias                 Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1                                                                                    0.083        -2.90       0.043        0.499",
      "     2                                  omega[two-sided: .1] ~ CumDirichlet(1, 1)         0.083        -2.60       0.059        0.686",
      "     3                                                   PET ~ Normal(0, 1)[0, Inf]       0.083        -1.46       0.185        2.492",
      "     4                InvGamma(1, 0.15)                                                   0.083        -2.66       0.056        0.647",
      "     5                InvGamma(1, 0.15) omega[two-sided: .1] ~ CumDirichlet(1, 1)         0.083        -2.61       0.058        0.683",
      "     6                InvGamma(1, 0.15)                  PET ~ Normal(0, 1)[0, Inf]       0.083        -1.73       0.141        1.799",
      "     7 Normal(0, 1)                                                                       0.083        -2.01       0.107        1.313",
      "     8 Normal(0, 1)                     omega[two-sided: .1] ~ CumDirichlet(1, 1)         0.083        -2.28       0.081        0.969",
      "     9 Normal(0, 1)                                      PET ~ Normal(0, 1)[0, Inf]       0.083        -2.25       0.083        0.999",
      "    10 Normal(0, 1)   InvGamma(1, 0.15)                                                   0.083        -2.38       0.073        0.870",
      "    11 Normal(0, 1)   InvGamma(1, 0.15) omega[two-sided: .1] ~ CumDirichlet(1, 1)         0.083        -2.70       0.053        0.621",
      "    12 Normal(0, 1)   InvGamma(1, 0.15)                  PET ~ Normal(0, 1)[0, Inf]       0.083        -2.56       0.061        0.716"
    ))
})

test_that("Diagnostics summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]], type = "diagnostics"), print = TRUE, width = 200),
      read.table(file = file.path("../results/summary.diagnostics", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }
})

test_that("Individual summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      capture_output_lines(summary(saved_fits[[i]], type = "individual"), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary.individual", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # different effect size measure
  expect_equal(
    capture_output_lines(summary(saved_fits[[10]], type = "individual", output_scale = "fishers_z"), print = TRUE, width = 150),
    c("Call:"                                                                                                  ,
      "RoBMA(d = d, se = d_se, priors_bias = prior_weightfunction(\"one-sided.fixed\", "                       ,
      "    list(c(0.1), c(1, 0.5))), priors_effect_null = NULL, priors_heterogeneity_null = NULL, "            ,
      "    priors_bias_null = NULL, parallel = TRUE, seed = 1)"                                                ,
      ""                                                                                                       ,
      "Robust Bayesian meta-analysis                                                                          ",
      " Model              1                        Parameter prior distributions"                             ,
      " Prior prob.    1.000                               mu ~ Normal(0, 1)     "                             ,
      " log(marglik)   -1.04                              tau ~ InvGamma(1, 0.15)"                             ,
      " Post. prob.    1.000             omega[one-sided: .1] = (0.5, 1)         "                             ,
      " Inclusion BF     Inf                                                     "                             ,
      ""                                                                                                       ,
      "Parameter estimates:"                                                                                   ,
      "              Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD  ESS R-hat"                     ,
      "mu           0.176 0.134 -0.097  0.177 0.428     0.00153          0.000 8076 1.000"                     ,
      "tau          0.106 0.102  0.019  0.075 0.369     0.00124          0.000 6744 1.000"                     ,
      "omega[0,0.1] 1.000 0.000  1.000  1.000 1.000          NA             NA   NA    NA"                     ,
      "omega[0.1,1] 0.500 0.000  0.500  0.500 0.500          NA             NA   NA    NA"                     ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))

  # different effect size measure (with PET and PEESE)
  expect_equal(
    capture_output_lines(summary(saved_fits[[3]], type = "individual", output_scale = "fishers_z"), print = TRUE, width = 150)[130:155],
    c(  " Model             11             Parameter prior distributions"                                     ,
        " Prior prob.    0.062                 mu ~ Normal(0, 1)        "                                     ,
        " log(marglik)   -0.87                tau ~ InvGamma(1, 0.15)   "                                     ,
        " Post. prob.    0.044                PET ~ Cauchy(0, 1)[0, Inf]"                                     ,
        " Inclusion BF   0.695                                          "                                     ,
        ""                                                                                                    ,
        "Parameter estimates:"                                                                                ,
        "     Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD  ESS R-hat"                           ,
        "mu  0.090 0.184 -0.335  0.108 0.403     0.00357          0.000 2783 1.003"                           ,
        "tau 0.115 0.121  0.018  0.078 0.429     0.00146          0.000 6854 1.000"                           ,
        "PET 0.730 0.646  0.023  0.556 2.439     0.01383          0.021 2184 1.002"                           ,
        "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale).",
        ""                                                                                                    ,
        "                                                               "                                     ,
        " Model             12             Parameter prior distributions"                                     ,
        " Prior prob.    0.062                 mu ~ Normal(0, 1)        "                                     ,
        " log(marglik)   -1.71                tau ~ InvGamma(1, 0.15)   "                                     ,
        " Post. prob.    0.019              PEESE ~ Cauchy(0, 5)[0, Inf]"                                     ,
        " Inclusion BF   0.291                                          "                                     ,
        ""                                                                                                    ,
        "Parameter estimates:"                                                                                ,
        "       Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD  ESS R-hat"                         ,
        "mu    0.102 0.179 -0.304  0.117 0.411     0.00311          0.000 3482 1.001"                         ,
        "tau   0.130 0.151  0.020  0.085 0.490     0.00228          0.000 4430 1.001"                         ,
        "PEESE 3.107 2.647  0.100  2.467 9.645     0.05316          0.000 2480 1.000"                         ,
        "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))

  # different effect size measure
  expect_equal(
    capture_output_lines(summary(saved_fits[[10]], type = "individual", output_scale = "fishers_z"), print = TRUE, width = 150),
    c("Call:"                                                                                                  ,
      "RoBMA(d = d, se = d_se, priors_bias = prior_weightfunction(\"one-sided.fixed\", "                       ,
      "    list(c(0.1), c(1, 0.5))), priors_effect_null = NULL, priors_heterogeneity_null = NULL, "            ,
      "    priors_bias_null = NULL, parallel = TRUE, seed = 1)"                                                ,
      ""                                                                                                       ,
      "Robust Bayesian meta-analysis                                                                          ",
      " Model              1                        Parameter prior distributions"                             ,
      " Prior prob.    1.000                               mu ~ Normal(0, 1)     "                             ,
      " log(marglik)   -1.04                              tau ~ InvGamma(1, 0.15)"                             ,
      " Post. prob.    1.000             omega[one-sided: .1] = (0.5, 1)         "                             ,
      " Inclusion BF     Inf                                                     "                             ,
      ""                                                                                                       ,
      "Parameter estimates:"                                                                                   ,
      "              Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD  ESS R-hat"                     ,
      "mu           0.176 0.134 -0.097  0.177 0.428     0.00153          0.000 8076 1.000"                     ,
      "tau          0.106 0.102  0.019  0.075 0.369     0.00124          0.000 6744 1.000"                     ,
      "omega[0,0.1] 1.000 0.000  1.000  1.000 1.000          NA             NA   NA    NA"                     ,
      "omega[0.1,1] 0.500 0.000  0.500  0.500 0.500          NA             NA   NA    NA"                     ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))
  # test short names
  expect_equal(
    capture_output_lines(summary(saved_fits[[10]], type = "individual", short_name = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                       ,
      "RoBMA(d = d, se = d_se, priors_bias = prior_weightfunction(\"one-sided.fixed\", "            ,
      "    list(c(0.1), c(1, 0.5))), priors_effect_null = NULL, priors_heterogeneity_null = NULL, " ,
      "    priors_bias_null = NULL, parallel = TRUE, seed = 1)"                                     ,
      ""                                                                                            ,
      "Robust Bayesian meta-analysis                                                               ",
      " Model              1             Parameter prior distributions"                             ,
      " Prior prob.    1.000                          mu ~ N(0, 1)    "                             ,
      " log(marglik)   -1.04                         tau ~ Ig(1, 0.15)"                             ,
      " Post. prob.    1.000               omega[1s: .1] = (0.5, 1)   "                             ,
      " Inclusion BF     Inf                                          "                             ,
      ""                                                                                            ,
      "Parameter estimates:"                                                                        ,
      "              Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD  ESS R-hat"          ,
      "mu           0.353 0.275 -0.195  0.357 0.883     0.00306          0.011 8076 1.000"          ,
      "tau          0.211 0.204  0.037  0.149 0.738     0.00249          0.012 6744 1.000"          ,
      "omega[0,0.1] 1.000 0.000  1.000  1.000 1.000          NA             NA   NA    NA"          ,
      "omega[0.1,1] 0.500 0.000  0.500  0.500 0.500          NA             NA   NA    NA"
    ))

  # test no spikes
  expect_equal(
    capture_output_lines(summary(saved_fits[[8]], type = "individual", remove_spike_0 = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                                                                     ,
      "RoBMA(d = d, se = d_se, priors_effect = NULL, priors_heterogeneity = NULL, "                                                               ,
      "    priors_bias = NULL, parallel = TRUE, seed = 1)"                                                                                        ,
      ""                                                                                                                                          ,
      "Robust Bayesian meta-analysis                                                               "                                              ,
      " Model              1             Parameter prior distributions"                                                                           ,
      " Prior prob.    1.000                                          "                                                                           ,
      " log(marglik)   -0.83                                          "                                                                           ,
      " Post. prob.    1.000                                          "                                                                           ,
      " Inclusion BF     Inf                                          "                                                                           ,
      ""                                                                                                                                          ,
      "Parameter estimates:"                                                                                                                      ,
      "[1] Mean           SD             lCI            Median         uCI            error(MCMC)    error(MCMC)/SD ESS            R-hat         ",
      "<0 rows> (or 0-length row.names)"
    ))


})

test_that("Interpret functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_fits)){
    expect_equal(
      gsub("(.{80})", "\\1\\\n", interpret(saved_fits[[i]])),
      read.table(file = file.path("../results/interpret", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # with transformation
  expect_equal(
    interpret(saved_fits[[1]], output_scale = "r"),
    "Robust Bayesian meta-analysis found weak evidence against the effect, BF_10 = 0.974, with mean model-averaged estimate correlation = 0.094, 95% CI [-0.032,  0.379]. Robust Bayesian meta-analysis found weak evidence against the heterogeneity, BF^rf = 0.858, with mean model-averaged estimate tau = 0.055, 95% CI [0.000, 0.312]. Robust Bayesian meta-analysis found weak evidence in favor of the publication bias, BF_pb = 1.17."
  )

})

#### creating / updating the test settings ####
if(FALSE){

  saved_fits       <- readRDS(file = "tests/results/saved_fits.RDS")

  # generate print files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(saved_fits[[i]], print = TRUE, width = 150), file = file.path("tests/results/print", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]]), print = TRUE, width = 150), file = file.path("tests/results/summary", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.models files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]], type = "models"), print = TRUE, width = 150), file = file.path("tests/results/summary.models", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.diagnostics files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]], type = "diagnostics"), print = TRUE, width = 200), file = file.path("tests/results/summary.diagnostics", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.individual files
  for(i in seq_along(saved_fits)){
    write.table(capture_output_lines(summary(saved_fits[[i]], type = "individual"), print = TRUE, width = 150), file = file.path("tests/results/summary.individual", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.individual files
  for(i in seq_along(saved_fits)){
    write.table(gsub("(.{80})", "\\1\\\n", interpret(saved_fits[[i]])), file = file.path("tests/results/interpret", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

}
