context("(5) Print and summary functions")
skip_on_cran()

# the summary tables and print functions are imported from BayesTools and tested henceforth
# test objects - assuming that the fit function worked properly
saved_files <- paste0("fit_", 1:16, ".RDS")
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
      "              Models Prior prob. Post. prob. log(Exclusion BF)"         ,
      "Effect         18/36       0.500       0.493             0.026"         ,
      "Heterogeneity  18/36       0.500       0.462             0.153"         ,
      "Bias           32/36       0.500       0.540            -0.159"         ,
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
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
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
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
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

  expect_equal(
    capture_output_lines(summary(saved_fits[[1]], output_scale = "r"), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                                    ,
      "RoBMA(d = d, se = d_se, parallel = TRUE, seed = 1)"                                                                                                                       ,
      ""                                                                                                                                                                         ,
      "Robust Bayesian meta-analysis"                                                                                                                                            ,
      "Components summary:"                                                                                                                                                      ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                                                                                                ,
      "Effect         18/36       0.500       0.493        0.974"                                                                                                                ,
      "Heterogeneity  18/36       0.500       0.462        0.858"                                                                                                                ,
      "Bias           32/36       0.500       0.540        1.173"                                                                                                                ,
      ""                                                                                                                                                                         ,
      "Model-averaged estimates:"                                                                                                                                                ,
      "                   Mean Median  0.025 0.975"                                                                                                                              ,
      "mu                0.094  0.000 -0.032 0.379"                                                                                                                              ,
      "tau               0.055  0.000  0.000 0.312"                                                                                                                              ,
      "omega[0,0.025]    1.000  1.000  1.000 1.000"                                                                                                                              ,
      "omega[0.025,0.05] 0.924  1.000  0.309 1.000"                                                                                                                              ,
      "omega[0.05,0.5]   0.824  1.000  0.126 1.000"                                                                                                                              ,
      "omega[0.5,0.95]   0.759  1.000  0.027 1.000"                                                                                                                              ,
      "omega[0.95,0.975] 0.770  1.000  0.027 1.000"                                                                                                                              ,
      "omega[0.975,1]    0.801  1.000  0.027 1.000"                                                                                                                              ,
      "PET               0.110  0.000  0.000 1.310"                                                                                                                              ,
      "PEESE             0.175  0.000  0.000 3.142"                                                                                                                              ,
      "The effect size estimates are summarized on the correlation scale and heterogeneity is summarized on the Fisher's z scale (priors were specified on the Cohen's d scale).",
      "(Estimated publication weights omega correspond to one-sided p-values.)"
      )
  )

  expect_equal(
    capture_output_lines(summary(saved_fits[[1]], output_scale = "OR"), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                      ,
      "RoBMA(d = d, se = d_se, parallel = TRUE, seed = 1)"                                                                                                         ,
      ""                                                                                                                                                           ,
      "Robust Bayesian meta-analysis"                                                                                                                              ,
      "Components summary:"                                                                                                                                        ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                                                                                  ,
      "Effect         18/36       0.500       0.493        0.974"                                                                                                  ,
      "Heterogeneity  18/36       0.500       0.462        0.858"                                                                                                  ,
      "Bias           32/36       0.500       0.540        1.173"                                                                                                  ,
      ""                                                                                                                                                           ,
      "Model-averaged estimates:"                                                                                                                                  ,
      "                   Mean Median 0.025 0.975"                                                                                                                 ,
      "mu                1.650  1.000 0.892 4.423"                                                                                                                 ,
      "tau               0.198  0.000 0.000 1.133"                                                                                                                 ,
      "omega[0,0.025]    1.000  1.000 1.000 1.000"                                                                                                                 ,
      "omega[0.025,0.05] 0.924  1.000 0.309 1.000"                                                                                                                 ,
      "omega[0.05,0.5]   0.824  1.000 0.126 1.000"                                                                                                                 ,
      "omega[0.5,0.95]   0.759  1.000 0.027 1.000"                                                                                                                 ,
      "omega[0.95,0.975] 0.770  1.000 0.027 1.000"                                                                                                                 ,
      "omega[0.975,1]    0.801  1.000 0.027 1.000"                                                                                                                 ,
      "PET               0.110  0.000 0.000 1.310"                                                                                                                 ,
      "PEESE             0.888  0.551 0.551 2.653"                                                                                                                 ,
      "The effect size estimates are summarized on the OR scale and heterogeneity is summarized on the logOR scale (priors were specified on the Cohen's d scale).",
      "(Estimated publication weights omega correspond to one-sided p-values.)"
    )
  )

  expect_equal(
    capture_output_lines(summary(saved_fits[[1]], output_scale = "logOR"), print = TRUE, width = 150),
    c("Call:"                                                                                            ,
      "RoBMA(d = d, se = d_se, parallel = TRUE, seed = 1)"                                               ,
      ""                                                                                                 ,
      "Robust Bayesian meta-analysis"                                                                    ,
      "Components summary:"                                                                              ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                        ,
      "Effect         18/36       0.500       0.493        0.974"                                        ,
      "Heterogeneity  18/36       0.500       0.462        0.858"                                        ,
      "Bias           32/36       0.500       0.540        1.173"                                        ,
      ""                                                                                                 ,
      "Model-averaged estimates:"                                                                        ,
      "                   Mean Median  0.025 0.975"                                                      ,
      "mu                0.355  0.000 -0.114 1.487"                                                      ,
      "tau               0.198  0.000  0.000 1.133"                                                      ,
      "omega[0,0.025]    1.000  1.000  1.000 1.000"                                                      ,
      "omega[0.025,0.05] 0.924  1.000  0.309 1.000"                                                      ,
      "omega[0.05,0.5]   0.824  1.000  0.126 1.000"                                                      ,
      "omega[0.5,0.95]   0.759  1.000  0.027 1.000"                                                      ,
      "omega[0.95,0.975] 0.770  1.000  0.027 1.000"                                                      ,
      "omega[0.975,1]    0.801  1.000  0.027 1.000"                                                      ,
      "PET               0.110  0.000  0.000 1.310"                                                      ,
      "PEESE             0.048  0.000  0.000 0.866"                                                      ,
      "The estimates are summarized on the log(OR) scale (priors were specified on the Cohen's d scale).",
      "(Estimated publication weights omega correspond to one-sided p-values.)"
    )
  )

  expect_equal( # try with BiBMA
    capture_output_lines(summary(saved_fits[[16]]), print = TRUE, width = 150),
    c("Call:"                                                                                          ,
      "BiBMA(x1 = 0:4, x2 = 2:6, n1 = rep(20, 5), n2 = rep(20, 5), chains = 2, "                       ,
      "    sample = 500, burnin = 250, adapt = 100, parallel = TRUE, "                                 ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                              ,
      ""                                                                                               ,
      "Bayesian model-averaged meta-analysis (binomial-normal model)"                                  ,
      "Components summary:"                                                                            ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                      ,
      "Effect           2/4       0.500       0.655        1.896"                                      ,
      "Heterogeneity    2/4       0.500       0.415        0.709"                                      ,
      ""                                                                                               ,
      "Model-averaged estimates:"                                                                      ,
      "      Mean Median  0.025 0.975"                                                                 ,
      "mu  -0.376 -0.295 -1.239 0.055"                                                                 ,
      "tau  0.150  0.000  0.000 0.792"                                                                 ,
      "The estimates are summarized on the log(OR) scale (priors were specified on the log(OR) scale)."
    )
  )

  expect_equal( # try with BiBMA
    capture_output_lines(summary(saved_fits[[16]], output_scale = "cohens_d"), print = TRUE, width = 150),
    c("Call:"                                                                                            ,
      "BiBMA(x1 = 0:4, x2 = 2:6, n1 = rep(20, 5), n2 = rep(20, 5), chains = 2, "                         ,
      "    sample = 500, burnin = 250, adapt = 100, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                  ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                ,
      ""                                                                                                 ,
      "Bayesian model-averaged meta-analysis (binomial-normal model)"                                    ,
      "Components summary:"                                                                              ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                        ,
      "Effect           2/4       0.500       0.655        1.896"                                        ,
      "Heterogeneity    2/4       0.500       0.415        0.709"                                        ,
      ""                                                                                                 ,
      "Model-averaged estimates:"                                                                        ,
      "      Mean Median  0.025 0.975"                                                                   ,
      "mu  -0.207 -0.162 -0.683 0.030"                                                                   ,
      "tau  0.083  0.000  0.000 0.437"                                                                   ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the log(OR) scale)."
    )
  )

  # no conditional, yet requested
  expect_equal(
    capture_output_lines(summary(saved_fits[[8]], conditional = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                                        ,
      "RoBMA(d = d, se = d_se, priors_effect = NULL, priors_heterogeneity = NULL, "                                  ,
      "    priors_bias = NULL, chains = 2, sample = 500, burnin = 250, "                                             ,
      "    adapt = 100, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, ",
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                            ,
      ""                                                                                                             ,
      "Robust Bayesian meta-analysis"                                                                                ,
      "Components summary:"                                                                                          ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                                    ,
      "Effect           0/1       0.000       0.000        0.000"                                                    ,
      "Heterogeneity    0/1       0.000       0.000        0.000"                                                    ,
      ""                                                                                                             ,
      "Model-averaged estimates:"                                                                                    ,
      "     Mean Median 0.025 0.975"                                                                                 ,
      "mu  0.000  0.000 0.000 0.000"                                                                                 ,
      "tau 0.000  0.000 0.000 0.000"                                                                                 ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."          ,
      ""                                                                                                             ,
      "Conditional estimates:"                                                                                       ,
      "[1] Mean   Median 0.025  0.975 "                                                                              ,
      "<0 rows> (or 0-length row.names)"                                                                             ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
  ))
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
      "    chains = 2, sample = 500, burnin = 250, adapt = 100, parallel = TRUE, "                                               ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                                          ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                                        ,
      ""                                                                                                                         ,
      "Robust Bayesian meta-analysis"                                                                                            ,
      "Models overview:"                                                                                                         ,
      " Model Prior Effect Prior Heterogeneity              Prior Bias         Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1         S(0)                S(0)                                       0.083        -2.90       0.043        0.499",
      "     2         S(0)                S(0) omega[2s: .1] ~ CumD(1, 1)            0.083        -2.60       0.059        0.687",
      "     3         S(0)                S(0)           PET ~ N(0, 1)[0, Inf]       0.083        -1.45       0.185        2.501",
      "     4         S(0)         Ig(1, 0.15)                                       0.083        -2.66       0.055        0.646",
      "     5         S(0)         Ig(1, 0.15) omega[2s: .1] ~ CumD(1, 1)            0.083        -2.61       0.058        0.683",
      "     6         S(0)         Ig(1, 0.15)           PET ~ N(0, 1)[0, Inf]       0.083        -1.73       0.140        1.796",
      "     7      N(0, 1)                S(0)                                       0.083        -2.00       0.107        1.316",
      "     8      N(0, 1)                S(0) omega[2s: .1] ~ CumD(1, 1)            0.083        -2.29       0.080        0.958",
      "     9      N(0, 1)                S(0)           PET ~ N(0, 1)[0, Inf]       0.083        -2.24       0.084        1.015",
      "    10      N(0, 1)         Ig(1, 0.15)                                       0.083        -2.40       0.072        0.857",
      "    11      N(0, 1)         Ig(1, 0.15) omega[2s: .1] ~ CumD(1, 1)            0.083        -2.67       0.055        0.638",
      "    12      N(0, 1)         Ig(1, 0.15)           PET ~ N(0, 1)[0, Inf]       0.083        -2.59       0.060        0.698"
    ))

  # test no spikes
  expect_equal(
    capture_output_lines(summary(saved_fits[[11]], type = "models", remove_spike_0 = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                                                                ,
      "RoBMA(y = d, se = d_se, priors_bias = list(prior_weightfunction(\"two-sided\", "                                                      ,
      "    list(c(0.1), c(1, 1))), prior_PET(\"normal\", list(0, 1))), "                                                                     ,
      "    chains = 2, sample = 500, burnin = 250, adapt = 100, parallel = TRUE, "                                                           ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                                                      ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                                                    ,
      ""                                                                                                                                     ,
      "Robust Bayesian meta-analysis"                                                                                                        ,
      "Models overview:"                                                                                                                     ,
      " Model Prior Effect Prior Heterogeneity                  Prior Bias                 Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1                                                                                    0.083        -2.90       0.043        0.499",
      "     2                                  omega[two-sided: .1] ~ CumDirichlet(1, 1)         0.083        -2.60       0.059        0.687",
      "     3                                                   PET ~ Normal(0, 1)[0, Inf]       0.083        -1.45       0.185        2.501",
      "     4                InvGamma(1, 0.15)                                                   0.083        -2.66       0.055        0.646",
      "     5                InvGamma(1, 0.15) omega[two-sided: .1] ~ CumDirichlet(1, 1)         0.083        -2.61       0.058        0.683",
      "     6                InvGamma(1, 0.15)                  PET ~ Normal(0, 1)[0, Inf]       0.083        -1.73       0.140        1.796",
      "     7 Normal(0, 1)                                                                       0.083        -2.00       0.107        1.316",
      "     8 Normal(0, 1)                     omega[two-sided: .1] ~ CumDirichlet(1, 1)         0.083        -2.29       0.080        0.958",
      "     9 Normal(0, 1)                                      PET ~ Normal(0, 1)[0, Inf]       0.083        -2.24       0.084        1.015",
      "    10 Normal(0, 1)   InvGamma(1, 0.15)                                                   0.083        -2.40       0.072        0.857",
      "    11 Normal(0, 1)   InvGamma(1, 0.15) omega[two-sided: .1] ~ CumDirichlet(1, 1)         0.083        -2.67       0.055        0.638",
      "    12 Normal(0, 1)   InvGamma(1, 0.15)                  PET ~ Normal(0, 1)[0, Inf]       0.083        -2.59       0.060        0.698"
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
    c("Call:"                                                                                                        ,
      "RoBMA(d = d, se = d_se, priors_bias = prior_weightfunction(\"one-sided.fixed\", "                             ,
      "    list(c(0.1), c(1, 0.5))), priors_effect_null = NULL, priors_heterogeneity_null = NULL, "                  ,
      "    priors_bias_null = NULL, chains = 2, sample = 500, burnin = 250, "                                        ,
      "    adapt = 100, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, ",
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                            ,
      ""                                                                                                             ,
      "Robust Bayesian meta-analysis                                                                          "      ,
      " Model              1                        Parameter prior distributions"                                   ,
      " Prior prob.    1.000                               mu ~ Normal(0, 1)     "                                   ,
      " log(marglik)   -0.92                              tau ~ InvGamma(1, 0.15)"                                   ,
      " Post. prob.    1.000             omega[one-sided: .1] = (0.5, 1)         "                                   ,
      " Inclusion BF     Inf                                                     "                                   ,
      ""                                                                                                             ,
      "Parameter estimates:"                                                                                         ,
      "              Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD ESS R-hat"                            ,
      "mu           0.158 0.136 -0.136  0.160 0.402     0.00632          0.046 484 1.001"                            ,
      "tau          0.111 0.111  0.021  0.078 0.407     0.00496          0.045 500 1.025"                            ,
      "omega[0,0.1] 1.000 0.000  1.000  1.000 1.000          NA             NA  NA    NA"                            ,
      "omega[0.1,1] 0.500 0.000  0.500  0.500 0.500          NA             NA  NA    NA"                            ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))

  # different effect size measure (with PET and PEESE)
  expect_equal(
    capture_output_lines(summary(saved_fits[[3]], type = "individual", output_scale = "fishers_z"), print = TRUE, width = 150)[132:158],
    c("                                                               "                                     ,
      " Model             11             Parameter prior distributions"                                     ,
      " Prior prob.    0.062                 mu ~ Normal(0, 1)        "                                     ,
      " log(marglik)   -0.88                tau ~ InvGamma(1, 0.15)   "                                     ,
      " Post. prob.    0.044                PET ~ Cauchy(0, 1)[0, Inf]"                                     ,
      " Inclusion BF   0.686                                          "                                     ,
      ""                                                                                                    ,
      "Parameter estimates:"                                                                                ,
      "     Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD ESS R-hat"                            ,
      "mu  0.098 0.171 -0.280  0.107 0.388     0.01005          0.059 300 1.004"                            ,
      "tau 0.112 0.109  0.019  0.075 0.387     0.00495          0.045 487 1.004"                            ,
      "PET 0.689 0.584  0.029  0.544 2.226     0.04007          0.069 212 1.001"                            ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale).",
      ""                                                                                                    ,
      "                                                               "                                     ,
      " Model             12             Parameter prior distributions"                                     ,
      " Prior prob.    0.062                 mu ~ Normal(0, 1)        "                                     ,
      " log(marglik)   -1.81                tau ~ InvGamma(1, 0.15)   "                                     ,
      " Post. prob.    0.017              PEESE ~ Cauchy(0, 5)[0, Inf]"                                     ,
      " Inclusion BF   0.265                                          "                                     ,
      ""                                                                                                    ,
      "Parameter estimates:"                                                                                ,
      "       Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD ESS R-hat"                          ,
      "mu    0.100 0.190 -0.333  0.113 0.436     0.01353          0.071 209 1.018"                          ,
      "tau   0.148 0.157  0.026  0.097 0.560     0.00923          0.059 288 1.068"                          ,
      "PEESE 2.962 2.675  0.095  2.222 9.717     0.17426          0.065 236 1.026"                          ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))

  # different effect size measure
  expect_equal(
    capture_output_lines(summary(saved_fits[[10]], type = "individual", output_scale = "fishers_z"), print = TRUE, width = 150),
    c("Call:"                                                                                                        ,
      "RoBMA(d = d, se = d_se, priors_bias = prior_weightfunction(\"one-sided.fixed\", "                             ,
      "    list(c(0.1), c(1, 0.5))), priors_effect_null = NULL, priors_heterogeneity_null = NULL, "                  ,
      "    priors_bias_null = NULL, chains = 2, sample = 500, burnin = 250, "                                        ,
      "    adapt = 100, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, ",
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                            ,
      ""                                                                                                             ,
      "Robust Bayesian meta-analysis                                                                          "      ,
      " Model              1                        Parameter prior distributions"                                   ,
      " Prior prob.    1.000                               mu ~ Normal(0, 1)     "                                   ,
      " log(marglik)   -0.92                              tau ~ InvGamma(1, 0.15)"                                   ,
      " Post. prob.    1.000             omega[one-sided: .1] = (0.5, 1)         "                                   ,
      " Inclusion BF     Inf                                                     "                                   ,
      ""                                                                                                             ,
      "Parameter estimates:"                                                                                         ,
      "              Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD ESS R-hat"                            ,
      "mu           0.158 0.136 -0.136  0.160 0.402     0.00632          0.046 484 1.001"                            ,
      "tau          0.111 0.111  0.021  0.078 0.407     0.00496          0.045 500 1.025"                            ,
      "omega[0,0.1] 1.000 0.000  1.000  1.000 1.000          NA             NA  NA    NA"                            ,
      "omega[0.1,1] 0.500 0.000  0.500  0.500 0.500          NA             NA  NA    NA"                            ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))
  # test short names
  expect_equal(
    capture_output_lines(summary(saved_fits[[10]], type = "individual", short_name = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                                        ,
      "RoBMA(d = d, se = d_se, priors_bias = prior_weightfunction(\"one-sided.fixed\", "                             ,
      "    list(c(0.1), c(1, 0.5))), priors_effect_null = NULL, priors_heterogeneity_null = NULL, "                  ,
      "    priors_bias_null = NULL, chains = 2, sample = 500, burnin = 250, "                                        ,
      "    adapt = 100, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, ",
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                            ,
      ""                                                                                                             ,
      "Robust Bayesian meta-analysis                                                               "                 ,
      " Model              1             Parameter prior distributions"                                              ,
      " Prior prob.    1.000                          mu ~ N(0, 1)    "                                              ,
      " log(marglik)   -0.92                         tau ~ Ig(1, 0.15)"                                              ,
      " Post. prob.    1.000               omega[1s: .1] = (0.5, 1)   "                                              ,
      " Inclusion BF     Inf                                          "                                              ,
      ""                                                                                                             ,
      "Parameter estimates:"                                                                                         ,
      "              Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD ESS R-hat"                            ,
      "mu           0.321 0.278 -0.272  0.322 0.827     0.01265          0.045 484 1.001"                            ,
      "tau          0.222 0.222  0.042  0.157 0.814     0.00993          0.045 500 1.025"                            ,
      "omega[0,0.1] 1.000 0.000  1.000  1.000 1.000          NA             NA  NA    NA"                            ,
      "omega[0.1,1] 0.500 0.000  0.500  0.500 0.500          NA             NA  NA    NA"                            ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    ))

  # test no spikes
  expect_equal(
    capture_output_lines(summary(saved_fits[[8]], type = "individual", remove_spike_0 = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                                                                     ,
      "RoBMA(d = d, se = d_se, priors_effect = NULL, priors_heterogeneity = NULL, "                                                               ,
      "    priors_bias = NULL, chains = 2, sample = 500, burnin = 250, "                                                                          ,
      "    adapt = 100, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                             ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                                                         ,
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
      "<0 rows> (or 0-length row.names)"                                                                                                          ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
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

test_that("Marginal summary functions work", {

  expect_error(marginal_summary(saved_fits[[1]]), "'marginal_summary' function is available only for RoBMA regression models")

  expect_equal(
    capture_output_lines(marginal_summary(saved_fits[[14]]), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                                               ,
      "RoBMA.reg(formula = ~mod_cat + mod_con, data = df_reg, priors_bias = NULL, "                                                                                                         ,
      "    chains = 2, sample = 500, burnin = 250, adapt = 100, parallel = TRUE, "                                                                                                          ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                                                                                                     ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                                                                                                   ,
      ""                                                                                                                                                                                    ,
      "Robust Bayesian meta-analysis"                                                                                                                                                       ,
      "Model-averaged marginal estimates:"                                                                                                                                                  ,
      "                Mean Median  0.025  0.975 Inclusion BF"                                                                                                                              ,
      "intercept      0.000  0.000  0.000  0.000        0.024"                                                                                                                              ,
      "mod_cat[A]    -0.847 -0.853 -0.996 -0.571          Inf"                                                                                                                              ,
      "mod_cat[B]     0.002  0.001 -0.064  0.080        0.102"                                                                                                                              ,
      "mod_cat[C]     0.845  0.870  0.569  1.021          Inf"                                                                                                                              ,
      "mod_con[-1SD] -0.087 -0.010 -0.299  0.000        7.962"                                                                                                                              ,
      "mod_con[0SD]   0.000  0.000  0.000  0.000        0.419"                                                                                                                              ,
      "mod_con[1SD]   0.087  0.012  0.000  0.299       10.197"                                                                                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."                                                                                 ,
      "\033[0;31mModel (13): ESS 6 is lower than the set target (10).\033[0m"                                                                                                               ,
      "\033[0;31mModel (15): ESS 5 is lower than the set target (10).\033[0m"                                                                                                               ,
      "\033[0;31mModel (16): R-hat 2.51 is larger than the set target (2).\033[0m"                                                                                                          ,
      "\033[0;31mmu_mod_cat[A]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m"                     ,
      "\033[0;31mmu_mod_cat[C]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m"                     ,
      "\033[0;31mmu_mod_con[0SD]: There is a considerable cluster of posterior samples at the exact null hypothesis values. The Savage-Dickey density ratio is likely to be invalid.\033[0m",
      "\033[0;31mmu_mod_con[0SD]: There is a considerable cluster of prior samples at the exact null hypothesis values. The Savage-Dickey density ratio is likely to be invalid.\033[0m"
    )
  )

  expect_equal(
    capture_output_lines(marginal_summary(saved_fits[[15]], conditional = TRUE, output_scale = "r"), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                             ,
      "RoBMA.reg(formula = ~mod_con, data = df_reg, priors = list(mod_con = list(null = prior(\"normal\", "                                                               ,
      "    list(0, 0.05)), alt = prior(\"normal\", list(0.3, 0.15)))), "                                                                                                  ,
      "    priors_heterogeneity = NULL, priors_bias = list(prior_weightfunction(distribution = \"two.sided\", "                                                           ,
      "        parameters = list(alpha = c(1, 1), steps = c(0.05)), "                                                                                                     ,
      "        prior_weights = 1/2), prior_PET(distribution = \"Cauchy\", "                                                                                               ,
      "        parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2)), "                                                                               ,
      "    priors_effect_null = NULL, chains = 2, sample = 500, burnin = 250, "                                                                                           ,
      "    adapt = 100, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                                                     ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                                                                                 ,
      ""                                                                                                                                                                  ,
      "Robust Bayesian meta-analysis"                                                                                                                                     ,
      "Model-averaged marginal estimates:"                                                                                                                                ,
      "                Mean Median  0.025  0.975 Inclusion BF"                                                                                                            ,
      "intercept     -0.020 -0.003 -0.216  0.022        0.029"                                                                                                            ,
      "mod_con[-1SD] -0.432 -0.422 -0.564 -0.395          Inf"                                                                                                            ,
      "mod_con[0SD]  -0.020 -0.003 -0.216  0.022        0.029"                                                                                                            ,
      "mod_con[1SD]   0.403  0.416  0.229  0.442          Inf"                                                                                                            ,
      "The estimates are summarized on the correlation scale (priors were specified on the Cohen's d scale)."                                                             ,
      "\033[0;31mModel (3): R-hat 2.133 is larger than the set target (2).\033[0m"                                                                                        ,
      "\033[0;31mModel (3): ESS 5 is lower than the set target (10).\033[0m"                                                                                              ,
      "\033[0;31mModel (3): MCMC error 1.33108 is larger than the set target (1).\033[0m"                                                                                 ,
      "\033[0;31mModel (6): ESS 7 is lower than the set target (10).\033[0m"                                                                                              ,
      "\033[0;31mmu_mod_con[-1SD]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m",
      "\033[0;31mmu_mod_con[1SD]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m" ,
      ""                                                                                                                                                                  ,
      "Conditional marginal estimates:"                                                                                                                                   ,
      "                Mean Median  0.025  0.975 Inclusion BF"                                                                                                            ,
      "intercept     -0.020 -0.003 -0.216  0.022        0.029"                                                                                                            ,
      "mod_con[-1SD] -0.432 -0.422 -0.564 -0.395          Inf"                                                                                                            ,
      "mod_con[0SD]  -0.020 -0.003 -0.216  0.022        0.029"                                                                                                            ,
      "mod_con[1SD]   0.403  0.416  0.229  0.442          Inf"                                                                                                            ,
      "The estimates are summarized on the correlation scale (priors were specified on the Cohen's d scale)."                                                             ,
      "\033[0;31mModel (3): R-hat 2.133 is larger than the set target (2).\033[0m"                                                                                        ,
      "\033[0;31mModel (3): ESS 5 is lower than the set target (10).\033[0m"                                                                                              ,
      "\033[0;31mModel (3): MCMC error 1.33108 is larger than the set target (1).\033[0m"                                                                                 ,
      "\033[0;31mModel (6): ESS 7 is lower than the set target (10).\033[0m"                                                                                              ,
      "\033[0;31mmu_mod_con[-1SD]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m",
      "\033[0;31mmu_mod_con[1SD]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m"
    )
  )


})

#### creating / updating the test settings ####
if(FALSE){

  saved_files <- paste0("fit_", 1:16, ".RDS")
  saved_fits  <- list()
  for(i in seq_along(saved_files)){
    saved_fits[[i]] <- readRDS(file = file.path("tests/results/fits", saved_files[i]))
  }

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
