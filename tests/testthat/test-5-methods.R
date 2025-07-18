context("(5) Print and summary functions")
skip_on_cran()

### Read all prefitted objects
# test objects - assuming that the fit function worked properly
temp_fits_dir <- Sys.getenv("ROBMA_TEST_FITS_DIR")
if (temp_fits_dir == "" || !dir.exists(temp_fits_dir)) {
  stop("Temporary fits directory not found. Run test-4-fit.R first.")
}

saved_files <- paste0("fit_", 1:16, ".RDS")
fits <- list()
for (i in 1:16) {
  fits[[i]] <- readRDS(file = file.path(temp_fits_dir, saved_files[i]))
}
names(fits) <- paste0("fit_", 1:16)


test_that("Print functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_files)){
    expect_equal(
      capture_output_lines(fits[[i]], print = TRUE, width = 150),
      read.table(file = file.path("../results/print", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }
})

test_that("Summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_files)){
    expect_equal(
      capture_output_lines(summary(fits[[i]]), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # all options
  expect_equal(
    capture_output_lines(summary(fits[["fit_1"]], conditional = TRUE, logBF = TRUE, BF01 = TRUE, probs = c(0.10, 0.50, .90)), print = TRUE, width = 150),
  c( "Call:"                                                                                                           ,
     "RoBMA(d = d, se = d_se, chains = 2, sample = 1000, burnin = 250, "                                               ,
     "    adapt = 250, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 1.25, ",
     "        min_ESS = 100, max_error = 1, max_SD_error = 1), seed = 1)"                                              ,
     ""                                                                                                                ,
     "Robust Bayesian meta-analysis"                                                                                   ,
     "Components summary:"                                                                                             ,
     "              Models Prior prob. Post. prob. log(Exclusion BF)"                                                  ,
     "Effect         18/36       0.500       0.494             0.025"                                                  ,
     "Heterogeneity  18/36       0.500       0.462             0.154"                                                 ,
     "Bias           32/36       0.500       0.537            -0.150"                                                 ,
     ""                                                                                                               ,
     "Model-averaged estimates:"                                                                                      ,
     "                   Mean Median   0.1   0.5   0.9"                                                               ,
     "mu                0.193  0.000 0.000 0.000 0.605"                                                               ,
     "tau               0.115  0.000 0.000 0.000 0.357"                                                               ,
     "omega[0,0.025]    1.000  1.000 1.000 1.000 1.000"                                                               ,
     "omega[0.025,0.05] 0.924  1.000 0.673 1.000 1.000"                                                               ,
     "omega[0.05,0.5]   0.823  1.000 0.309 1.000 1.000"                                                               ,
     "omega[0.5,0.95]   0.757  1.000 0.118 1.000 1.000"                                                               ,
     "omega[0.95,0.975] 0.768  1.000 0.123 1.000 1.000"                                                               ,
     "omega[0.975,1]    0.800  1.000 0.129 1.000 1.000"                                                               ,
     "PET               0.105  0.000 0.000 0.000 0.365"                                                               ,
     "PEESE             0.095  0.000 0.000 0.000 0.000"                                                               ,
     "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."            ,
     "(Estimated publication weights omega correspond to one-sided p-values.)"                                        ,
     ""                                                                                                               ,
     "Conditional estimates:"                                                                                         ,
     "                   Mean Median   0.1   0.5   0.9"                                                               ,
     "mu                0.399  0.412 0.059 0.412 0.724"                                                               ,
     "tau               0.251  0.181 0.064 0.181 0.508"                                                               ,
     "omega[0,0.025]    1.000  1.000 1.000 1.000 1.000"                                                               ,
     "omega[0.025,0.05] 0.783  0.889 0.366 0.889 1.000"                                                               ,
     "omega[0.05,0.5]   0.497  0.484 0.152 0.484 0.873"                                                               ,
     "omega[0.5,0.95]   0.310  0.243 0.039 0.243 0.706"                                                               ,
     "omega[0.95,0.975] 0.339  0.267 0.039 0.267 0.781"                                                               ,
     "omega[0.975,1]    0.435  0.305 0.039 0.305 1.000"                                                               ,
     "PET               0.789  0.706 0.160 0.706 1.506"                                                               ,
     "PEESE             1.678  1.522 0.341 1.522 3.175"                                                               ,
     "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."            ,
     "(Estimated publication weights omega correspond to one-sided p-values.)"
      )
  )

  # different effect size measure
  expect_equal(
    capture_output_lines(summary(fits[["fit_1"]], output_scale = "fishers_z"), print = TRUE, width = 150),
    c("Call:"                                                                                                           ,
      "RoBMA(d = d, se = d_se, chains = 2, sample = 1000, burnin = 250, "                                               ,
      "    adapt = 250, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 1.25, ",
      "        min_ESS = 100, max_error = 1, max_SD_error = 1), seed = 1)"                                              ,
      ""                                                                                                                ,
      "Robust Bayesian meta-analysis"                                                                                   ,
      "Components summary:"                                                                                             ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                                       ,
      "Effect         18/36       0.500       0.494        0.975"                                                       ,
      "Heterogeneity  18/36       0.500       0.462        0.858"                                                       ,
      "Bias           32/36       0.500       0.537        1.162"                                                       ,
      ""                                                                                                                ,
      "Model-averaged estimates:"                                                                                       ,
      "                   Mean Median  0.025 0.975"                                                                     ,
      "mu                0.095  0.000 -0.027 0.396"                                                                     ,
      "tau               0.058  0.000  0.000 0.330"                                                                     ,
      "omega[0,0.025]    1.000  1.000  1.000 1.000"                                                                     ,
      "omega[0.025,0.05] 0.924  1.000  0.311 1.000"                                                                     ,
      "omega[0.05,0.5]   0.823  1.000  0.130 1.000"                                                                     ,
      "omega[0.5,0.95]   0.757  1.000  0.027 1.000"                                                                     ,
      "omega[0.95,0.975] 0.768  1.000  0.028 1.000"                                                                     ,
      "omega[0.975,1]    0.800  1.000  0.028 1.000"                                                                     ,
      "PET               0.105  0.000  0.000 1.259"                                                                     ,
      "PEESE             0.190  0.000  0.000 3.359"                                                                     ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."            ,
      "(Estimated publication weights omega correspond to one-sided p-values.)"
      )
  )

  expect_equal(
    capture_output_lines(summary(fits[["fit_1"]], output_scale = "r"), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                                    ,
      "RoBMA(d = d, se = d_se, chains = 2, sample = 1000, burnin = 250, "                                                                                                        ,
      "    adapt = 250, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 1.25, "                                                         ,
      "        min_ESS = 100, max_error = 1, max_SD_error = 1), seed = 1)"                                                                                                       ,
      ""                                                                                                                                                                         ,
      "Robust Bayesian meta-analysis"                                                                                                                                            ,
      "Components summary:"                                                                                                                                                      ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                                                                                                ,
      "Effect         18/36       0.500       0.494        0.975"                                                                                                                ,
      "Heterogeneity  18/36       0.500       0.462        0.858"                                                                                                                ,
      "Bias           32/36       0.500       0.537        1.162"                                                                                                                ,
      ""                                                                                                                                                                         ,
      "Model-averaged estimates:"                                                                                                                                                ,
      "                   Mean Median  0.025 0.975"                                                                                                                              ,
      "mu                0.093  0.000 -0.027 0.376"                                                                                                                              ,
      "tau               0.058  0.000  0.000 0.330"                                                                                                                              ,
      "omega[0,0.025]    1.000  1.000  1.000 1.000"                                                                                                                              ,
      "omega[0.025,0.05] 0.924  1.000  0.311 1.000"                                                                                                                              ,
      "omega[0.05,0.5]   0.823  1.000  0.130 1.000"                                                                                                                              ,
      "omega[0.5,0.95]   0.757  1.000  0.027 1.000"                                                                                                                              ,
      "omega[0.95,0.975] 0.768  1.000  0.028 1.000"                                                                                                                              ,
      "omega[0.975,1]    0.800  1.000  0.028 1.000"                                                                                                                              ,
      "PET               0.105  0.000  0.000 1.259"                                                                                                                              ,
      "PEESE             0.190  0.000  0.000 3.359"                                                                                                                              ,
      "The effect size estimates are summarized on the correlation scale and heterogeneity is summarized on the Fisher's z scale (priors were specified on the Cohen's d scale).",
      "(Estimated publication weights omega correspond to one-sided p-values.)"
      )
  )

  expect_equal(
    capture_output_lines(summary(fits[["fit_1"]], output_scale = "OR"), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                      ,
      "RoBMA(d = d, se = d_se, chains = 2, sample = 1000, burnin = 250, "                                                                                          ,
      "    adapt = 250, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 1.25, "                                           ,
      "        min_ESS = 100, max_error = 1, max_SD_error = 1), seed = 1)"                                                                                         ,
      ""                                                                                                                                                           ,
      "Robust Bayesian meta-analysis"                                                                                                                              ,
      "Components summary:"                                                                                                                                        ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                                                                                  ,
      "Effect         18/36       0.500       0.494        0.975"                                                                                                  ,
      "Heterogeneity  18/36       0.500       0.462        0.858"                                                                                                  ,
      "Bias           32/36       0.500       0.537        1.162"                                                                                                  ,
      ""                                                                                                                                                           ,
      "Model-averaged estimates:"                                                                                                                                  ,
      "                   Mean Median 0.025 0.975"                                                                                                                 ,
      "mu                1.636  1.000 0.906 4.366"                                                                                                                 ,
      "tau               0.209  0.000 0.000 1.195"                                                                                                                 ,
      "omega[0,0.025]    1.000  1.000 1.000 1.000"                                                                                                                 ,
      "omega[0.025,0.05] 0.924  1.000 0.311 1.000"                                                                                                                 ,
      "omega[0.05,0.5]   0.823  1.000 0.130 1.000"                                                                                                                 ,
      "omega[0.5,0.95]   0.757  1.000 0.027 1.000"                                                                                                                 ,
      "omega[0.95,0.975] 0.768  1.000 0.028 1.000"                                                                                                                 ,
      "omega[0.975,1]    0.800  1.000 0.028 1.000"                                                                                                                 ,
      "PET               0.105  0.000 0.000 1.259"                                                                                                                 ,
      "PEESE             1.372  0.551 0.551 2.956"                                                                                                                 ,
      "The effect size estimates are summarized on the OR scale and heterogeneity is summarized on the logOR scale (priors were specified on the Cohen's d scale).",
      "(Estimated publication weights omega correspond to one-sided p-values.)"
    )
  )

  expect_equal(
    capture_output_lines(summary(fits[["fit_1"]], output_scale = "logOR"), print = TRUE, width = 150),
    c("Call:"                                                                                                           ,
      "RoBMA(d = d, se = d_se, chains = 2, sample = 1000, burnin = 250, "                                               ,
      "    adapt = 250, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 1.25, ",
      "        min_ESS = 100, max_error = 1, max_SD_error = 1), seed = 1)"                                              ,
      ""                                                                                                                ,
      "Robust Bayesian meta-analysis"                                                                                   ,
      "Components summary:"                                                                                             ,
      "              Models Prior prob. Post. prob. Inclusion BF"                                                       ,
      "Effect         18/36       0.500       0.494        0.975"                                                       ,
      "Heterogeneity  18/36       0.500       0.462        0.858"                                                       ,
      "Bias           32/36       0.500       0.537        1.162"                                                       ,
      ""                                                                                                                ,
      "Model-averaged estimates:"                                                                                       ,
      "                   Mean Median  0.025 0.975"                                                                     ,
      "mu                0.351  0.000 -0.098 1.474"                                                                     ,
      "tau               0.209  0.000  0.000 1.195"                                                                     ,
      "omega[0,0.025]    1.000  1.000  1.000 1.000"                                                                     ,
      "omega[0.025,0.05] 0.924  1.000  0.311 1.000"                                                                     ,
      "omega[0.05,0.5]   0.823  1.000  0.130 1.000"                                                                     ,
      "omega[0.5,0.95]   0.757  1.000  0.027 1.000"                                                                     ,
      "omega[0.95,0.975] 0.768  1.000  0.028 1.000"                                                                     ,
      "omega[0.975,1]    0.800  1.000  0.028 1.000"                                                                     ,
      "PET               0.105  0.000  0.000 1.259"                                                                     ,
      "PEESE             0.052  0.000  0.000 0.926"                                                                     ,
      "The estimates are summarized on the log(OR) scale (priors were specified on the Cohen's d scale)."               ,
      "(Estimated publication weights omega correspond to one-sided p-values.)"
    )
  )

  expect_equal( # try with BiBMA
    capture_output_lines(summary(fits[["fit_16"]]), print = TRUE, width = 150),
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
    capture_output_lines(summary(fits[["fit_16"]], output_scale = "cohens_d"), print = TRUE, width = 150),
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
    capture_output_lines(summary(fits[["fit_8"]], conditional = TRUE), print = TRUE, width = 150),
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
  for(i in 1:length(saved_files)){
    expect_equal(
      capture_output_lines(summary(fits[[i]], type = "models"), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary.models", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # test short names
  expect_equal(
    capture_output_lines(summary(fits[["fit_1"]], type = "models", short_name = TRUE), print = TRUE, width = 150)[1:15],
    c("Call:"                                                                                                                                ,
      "RoBMA(d = d, se = d_se, chains = 2, sample = 1000, burnin = 250, "                                                                    ,
      "    adapt = 250, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 1.25, "                     ,
      "        min_ESS = 100, max_error = 1, max_SD_error = 1), seed = 1)"                                                                   ,
      ""                                                                                                                                     ,
      "Robust Bayesian meta-analysis"                                                                                                        ,
      "Models overview:"                                                                                                                     ,
      " Model Prior Effect Prior Heterogeneity                  Prior Bias                 Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1         S(0)                S(0)                                                   0.125        -0.83       0.076        0.580",
      "     2         S(0)                S(0)           omega[2s: .05] ~ CumD(1, 1)             0.010        -0.30       0.011        1.049",
      "     3         S(0)                S(0)       omega[2s: .1, .05] ~ CumD(1, 1, 1)          0.010        -0.40       0.010        0.948",
      "     4         S(0)                S(0)           omega[1s: .05] ~ CumD(1, 1)             0.010        -0.28       0.011        1.066",
      "     5         S(0)                S(0)     omega[1s: .05, .025] ~ CumD(1, 1, 1)          0.010         0.02       0.015        1.440",
      "     6         S(0)                S(0)       omega[1s: .5, .05] ~ CumD(1, 1, 1)          0.010         0.55       0.025        2.472",
      "     7         S(0)                S(0) omega[1s: .5, .05, .025] ~ CumD(1, 1, 1, 1)       0.010         0.79       0.032        3.170"
    ))

  # test no spikes
  expect_equal(
    capture_output_lines(summary(fits[["fit_1"]], type = "models", remove_spike_0 = TRUE), print = TRUE, width = 150)[1:15],
    c("Call:"                                                                                                                                               ,
      "RoBMA(d = d, se = d_se, chains = 2, sample = 1000, burnin = 250, "                                                                                   ,
      "    adapt = 250, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 1.25, "                                    ,
      "        min_ESS = 100, max_error = 1, max_SD_error = 1), seed = 1)"                                                                                  ,
      ""                                                                                                                                                    ,
      "Robust Bayesian meta-analysis"                                                                                                                       ,
      "Models overview:"                                                                                                                                    ,
      " Model Prior Effect Prior Heterogeneity                         Prior Bias                         Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1                                                                                                   0.125        -0.83       0.076        0.580",
      "     2                                            omega[two-sided: .05] ~ CumDirichlet(1, 1)             0.010        -0.30       0.011        1.049",
      "     3                                        omega[two-sided: .1, .05] ~ CumDirichlet(1, 1, 1)          0.010        -0.40       0.010        0.948",
      "     4                                            omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.010        -0.28       0.011        1.066",
      "     5                                      omega[one-sided: .05, .025] ~ CumDirichlet(1, 1, 1)          0.010         0.02       0.015        1.440",
      "     6                                        omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.010         0.55       0.025        2.472",
      "     7                                  omega[one-sided: .5, .05, .025] ~ CumDirichlet(1, 1, 1, 1)       0.010         0.79       0.032        3.170"
    ))
})

test_that("Diagnostics summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_files)){
    expect_equal(
      capture_output_lines(summary(fits[[i]], type = "diagnostics"), print = TRUE, width = 200),
      read.table(file = file.path("../results/summary.diagnostics", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }
})

test_that("Individual summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_files)){
    expect_equal(
      capture_output_lines(summary(fits[[i]], type = "individual"), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary.individual", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # different effect size measure
  expect_equal(
    capture_output_lines(summary(fits[["fit_10"]], type = "individual", output_scale = "fishers_z"), print = TRUE, width = 150),
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
      "mu           0.158 0.136 -0.136  0.160 0.402     0.00633          0.046 485 1.001"                            ,
      "tau          0.111 0.111  0.021  0.078 0.407     0.00498          0.045 500 1.027"                            ,
      "omega[0,0.1] 1.000 0.000  1.000  1.000 1.000          NA             NA  NA    NA"                            ,
      "omega[0.1,1] 0.500 0.000  0.500  0.500 0.500     0.00000             NA   0    NA"                            ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))

  # different effect size measure (with PET and PEESE)
  expect_equal(
    capture_output_lines(summary(fits[["fit_3"]], type = "individual", output_scale = "fishers_z"), print = TRUE, width = 150)[132:158],
    c("                                                               "                                     ,
      " Model             11             Parameter prior distributions"                                     ,
      " Prior prob.    0.062                 mu ~ Normal(0, 1)        "                                     ,
      " log(marglik)   -0.88                tau ~ InvGamma(1, 0.15)   "                                     ,
      " Post. prob.    0.044                PET ~ Cauchy(0, 1)[0, Inf]"                                     ,
      " Inclusion BF   0.686                                          "                                     ,
      ""                                                                                                    ,
      "Parameter estimates:"                                                                                ,
      "     Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD ESS R-hat"                            ,
      "mu  0.098 0.171 -0.280  0.107 0.388     0.01048          0.061 277 1.003"                            ,
      "tau 0.112 0.109  0.019  0.075 0.387     0.00515          0.047 487 1.002"                            ,
      "PET 0.689 0.584  0.029  0.544 2.226     0.04006          0.069 212 1.001"                            ,
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
      "mu    0.100 0.190 -0.333  0.113 0.436     0.01309          0.069 210 1.018"                          ,
      "tau   0.148 0.157  0.026  0.097 0.560     0.01629          0.104 288 1.089"                          ,
      "PEESE 2.962 2.675  0.095  2.222 9.717     0.18092          0.068 236 1.020"                          ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))

  # different effect size measure
  expect_equal(
    capture_output_lines(summary(fits[["fit_10"]], type = "individual", output_scale = "fishers_z"), print = TRUE, width = 150),
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
      "mu           0.158 0.136 -0.136  0.160 0.402     0.00633          0.046 485 1.001"                            ,
      "tau          0.111 0.111  0.021  0.078 0.407     0.00498          0.045 500 1.027"                            ,
      "omega[0,0.1] 1.000 0.000  1.000  1.000 1.000          NA             NA  NA    NA"                            ,
      "omega[0.1,1] 0.500 0.000  0.500  0.500 0.500     0.00000             NA   0    NA"                            ,
      "The estimates are summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))
  # test short names
  expect_equal(
    capture_output_lines(summary(fits[["fit_10"]], type = "individual", short_name = TRUE), print = TRUE, width = 150),
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
      "mu           0.321 0.278 -0.272  0.322 0.827     0.01295          0.047 484 1.001"                            ,
      "tau          0.222 0.222  0.042  0.157 0.814     0.00997          0.045 500 1.027"                            ,
      "omega[0,0.1] 1.000 0.000  1.000  1.000 1.000          NA             NA  NA    NA"                            ,
      "omega[0.1,1] 0.500 0.000  0.500  0.500 0.500     0.00000             NA   0    NA"                            ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    ))

  # test no spikes
  expect_equal(
    capture_output_lines(summary(fits[["fit_8"]], type = "individual", remove_spike_0 = TRUE), print = TRUE, width = 150),
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
  for(i in 1:length(saved_files)){
    expect_equal(
      gsub("(.{80})", "\\1\\\n", interpret(fits[[i]])),
      read.table(file = file.path("../results/interpret", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # with transformation
  expect_equal(
    interpret(fits[["fit_1"]], output_scale = "r"),
    "Robust Bayesian meta-analysis found weak evidence against the effect, BF_10 = 0.975, with mean model-averaged estimate correlation = 0.093, 95% CI [-0.027,  0.376]. Robust Bayesian meta-analysis found weak evidence against the heterogeneity, BF^rf = 0.858, with mean model-averaged estimate tau = 0.058, 95% CI [0.000, 0.330]. Robust Bayesian meta-analysis found weak evidence in favor of the publication bias, BF_pb = 1.16."
  )

})

test_that("Marginal summary functions work", {

  expect_error(marginal_summary(fits[["fit_1"]]), "'marginal_summary' function is available only for RoBMA regression models")

  expect_equal(
    capture_output_lines(marginal_summary(fits[["fit_14"]]), print = TRUE, width = 150),
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
      "mod_cat[B]     0.002  0.001 -0.064  0.080        0.079"                                                                                                                              ,
      "mod_cat[C]     0.845  0.870  0.569  1.021          Inf"                                                                                                                              ,
      "mod_con[-1SD] -0.087 -0.010 -0.299  0.000        1.627"                                                                                                                              ,
      "mod_con[0SD]   0.000  0.000  0.000  0.000        0.004"                                                                                                                              ,
      "mod_con[1SD]   0.087  0.012  0.000  0.299        1.775"                                                                                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."                                                                                 ,
      "\033[0;31mModel (6): R-hat 2.409 is larger than the set target (2).\033[0m"                                                                                                          ,
      "\033[0;31mModel (7): R-hat 2.124 is larger than the set target (2).\033[0m"                                                                                                          ,
      "\033[0;31mModel (8): R-hat 2.04 is larger than the set target (2).\033[0m"                                                                                                           ,
      "\033[0;31mModel (13): ESS 6 is lower than the set target (10).\033[0m"                                                                                                               ,
      "\033[0;31mModel (14): R-hat 3.9 is larger than the set target (2).\033[0m"                                                                                                           ,
      "\033[0;31mThere were another 1 warnings. To see all warnings call 'check_RoBMA(fit)'.\033[0m"                                                                                        ,
      "\033[0;31mmu_mod_cat[A]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m"                     ,
      "\033[0;31mmu_mod_cat[C]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m"                     ,
      "\033[0;31mmu_mod_con[0SD]: There is a considerable cluster of posterior samples at the exact null hypothesis values. The Savage-Dickey density ratio is likely to be invalid.\033[0m",
      "\033[0;31mmu_mod_con[0SD]: There is a considerable cluster of prior samples at the exact null hypothesis values. The Savage-Dickey density ratio is likely to be invalid.\033[0m"
    )
  )

  expect_equal(
    capture_output_lines(marginal_summary(fits[["fit_15"]], conditional = TRUE, output_scale = "r"), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                             ,
      "RoBMA.reg(formula = ~mod_con, data = df_reg, priors = list(mod_con = list(null = prior(\"normal\", "                                                               ,
      "    list(0, 0.05)), alt = prior(\"normal\", list(0.3, 0.15)))), "                                                                                                  ,
      "    priors_heterogeneity = NULL, priors_bias = list(prior_weightfunction(distribution = \"two.sided\", "                                                           ,
      "        parameters = list(alpha = c(1, 1), steps = c(0.05)), "                                                                                                     ,
      "        prior_weights = 1/2), prior_PET(distribution = \"Cauchy\", "                                                                                               ,
      "        parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2)), "                                                                               ,
      "    priors_effect_null = NULL, algorithm = \"ss\", chains = 2, "                                                                                                   ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                                                                                  ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                                                                                   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                                                                                 ,
      ""                                                                                                                                                                  ,
      "Robust Bayesian meta-analysis"                                                                                                                                     ,
      "Model-averaged marginal estimates:"                                                                                                                                ,
      "                Mean Median  0.025  0.975 Inclusion BF"                                                                                                            ,
      "intercept     -0.004 -0.002 -0.048  0.023           NA"                                                                                                            ,
      "mod_con[-1SD] -0.422 -0.421 -0.457 -0.395          Inf"                                                                                                            ,
      "mod_con[0SD]  -0.004 -0.002 -0.048  0.023        0.025"                                                                                                            ,
      "mod_con[1SD]   0.416  0.418  0.380  0.443          Inf"                                                                                                            ,
      "The estimates are summarized on the correlation scale (priors were specified on the Cohen's d scale)."                                                             ,
      "\033[0;31mmu_mod_con[-1SD]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m",
      "\033[0;31mmu_mod_con[1SD]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m" ,
      ""                                                                                                                                                                  ,
      "Conditional marginal estimates:"                                                                                                                                   ,
      "                Mean Median  0.025  0.975 Inclusion BF"                                                                                                            ,
      "intercept     -0.004 -0.002 -0.048  0.023           NA"                                                                                                            ,
      "mod_con[-1SD] -0.422 -0.421 -0.457 -0.395          Inf"                                                                                                            ,
      "mod_con[0SD]  -0.004 -0.002 -0.048  0.023        0.025"                                                                                                            ,
      "mod_con[1SD]   0.416  0.418  0.380  0.443          Inf"                                                                                                            ,
      "The estimates are summarized on the correlation scale (priors were specified on the Cohen's d scale)."                                                             ,
      "\033[0;31mmu_mod_con[-1SD]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m",
      "\033[0;31mmu_mod_con[1SD]: Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.\033[0m"
    )
  )
})

# test heterogeneity summary function
test_that("Heterogeneity summary functions work", {

  # testing consistency across all model specifications
  for(i in 1:length(saved_files)){
    expect_equal(
      capture_output_lines(summary_heterogeneity(fits[[i]]), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary_heterogeneity", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  # testing consistency across all model specifications
  for(i in 1:length(saved_files)){
    expect_equal(
      capture_output_lines(summary_heterogeneity(fits[[i]], type = "i"), print = TRUE, width = 150),
      read.table(file = file.path("../results/summary_heterogeneity.individual", paste0(i, ".txt")), header = FALSE, blank.lines.skip = FALSE)[,1])
  }

  expect_equal(
    capture_output_lines(summary_heterogeneity(fits[["fit_1"]], conditional = TRUE, output_scale = "logOR", probs = c(0.025, 0.5)), print = TRUE, width = 150),
    c("Call:"                                                                                                           ,
      "RoBMA(d = d, se = d_se, chains = 2, sample = 1000, burnin = 250, "                                               ,
      "    adapt = 250, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 1.25, ",
      "        min_ESS = 100, max_error = 1, max_SD_error = 1), seed = 1)"                                              ,
      ""                                                                                                                ,
      "Robust Bayesian meta-analysis"                                                                                   ,
      "Model-averaged heterogeneity estimates:"                                                                         ,
      "       Mean Median  0.025   0.5"                                                                                 ,
      "PI    0.354  0.163 -0.768 0.163"                                                                                 ,
      "tau   0.209  0.000  0.000 0.000"                                                                                ,
      "tau2  0.051  0.000  0.000 0.000"                                                                                ,
      "I2   10.424  0.000  0.000 0.000"                                                                                ,
      "H2    1.289  1.000  1.000 1.000"                                                                                ,
      "The prediction interval (PI) is summarized on the log(OR) scale."                                               ,
      "The absolute heterogeneity (tau, tau^2) is summarized on the log(OR) scale."                                    ,
      "The relative heterogeneity indicies (I^2 and H^2) were computed on the Fisher's z scale."                       ,
      ""                                                                                                               ,
      "Conditional heterogeneity estimates:"                                                                           ,
      "       Mean Median  0.025    0.5"                                                                               ,
      "PI    0.715  0.747 -0.570  0.747"                                                                               ,
      "tau   0.454  0.328  0.072  0.328"                                                                               ,
      "tau2  0.108  0.030  0.001  0.030"                                                                               ,
      "I2   22.526 14.392  0.812 14.392"                                                                               ,
      "H2    1.612  1.168  1.008  1.168"                                                                               ,
      "The prediction interval (PI) is summarized on the log(OR) scale."                                               ,
      "The absolute heterogeneity (tau, tau^2) is summarized on the log(OR) scale."                                    ,
      "The relative heterogeneity indicies (I^2 and H^2) were computed on the Fisher's z scale."
    ))

})

# test effect size summary function
test_that("Effect size summary functions work", {

  # testing for consistency among pooled vs adjusted for standard models
  expect_equivalent(
    as.data.frame(pooled_effect(fits[["fit_15"]])[["estimates"]]),
    as.data.frame(adjusted_effect(fits[["fit_15"]])[["estimates"]])
  )
  expect_equivalent(
    as.data.frame(pooled_effect(fits[["fit_15"]], conditional = TRUE)[["estimates_conditional"]]),
    as.data.frame(adjusted_effect(fits[["fit_15"]], conditional = TRUE)[["estimates_conditional"]])
  )

  expect_equal(
    capture_output_lines(pooled_effect(fits[["fit_15"]], conditional = TRUE, output_scale = "logOR", probs = c(0.025, 0.5)), print = TRUE, width = 150),
    c("Call:"                                                                                                  ,
      "RoBMA.reg(formula = ~mod_con, data = df_reg, priors = list(mod_con = list(null = prior(\"normal\", "    ,
      "    list(0, 0.05)), alt = prior(\"normal\", list(0.3, 0.15)))), "                                       ,
      "    priors_heterogeneity = NULL, priors_bias = list(prior_weightfunction(distribution = \"two.sided\", ",
      "        parameters = list(alpha = c(1, 1), steps = c(0.05)), "                                          ,
      "        prior_weights = 1/2), prior_PET(distribution = \"Cauchy\", "                                    ,
      "        parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2)), "                    ,
      "    priors_effect_null = NULL, algorithm = \"ss\", chains = 2, "                                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                       ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                        ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                      ,
      ""                                                                                                       ,
      "Robust Bayesian meta-regression"                                                                        ,
      "Model-averaged pooled effect estimate:"                                                                 ,
      "           Mean Median  0.025    0.5"                                                                   ,
      "estimate -0.013 -0.007 -0.175 -0.007"                                                                   ,
      "PI       -0.013 -0.007 -0.175 -0.007"                                                                   ,
      ""                                                                                                       ,
      "Conditional pooled effect estimate:"                                                                    ,
      "           Mean Median  0.025    0.5"                                                                   ,
      "estimate -0.013 -0.007 -0.175 -0.007"                                                                   ,
      "PI       -0.013 -0.007 -0.175 -0.007"                                                                   ,
      "The estimates are summarized on the log(OR) scale (priors were specified on the Cohen's d scale)."
    ))

  expect_equal(
    capture_output_lines(adjusted_effect(fits[["fit_15"]], conditional = TRUE, output_scale = "r", probs = c(0.025, 0.5)), print = TRUE, width = 150),
    c("Call:"                                                                                                                                                                    ,
      "RoBMA.reg(formula = ~mod_con, data = df_reg, priors = list(mod_con = list(null = prior(\"normal\", "                                                                      ,
      "    list(0, 0.05)), alt = prior(\"normal\", list(0.3, 0.15)))), "                                                                                                         ,
      "    priors_heterogeneity = NULL, priors_bias = list(prior_weightfunction(distribution = \"two.sided\", "                                                                  ,
      "        parameters = list(alpha = c(1, 1), steps = c(0.05)), "                                                                                                            ,
      "        prior_weights = 1/2), prior_PET(distribution = \"Cauchy\", "                                                                                                      ,
      "        parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2)), "                                                                                      ,
      "    priors_effect_null = NULL, algorithm = \"ss\", chains = 2, "                                                                                                          ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                                                                                         ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                                                                                          ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                                                                                        ,
      ""                                                                                                                                                                         ,
      "Robust Bayesian meta-regression"                                                                                                                                          ,
      "Model-averaged adjusted effect estimate:"                                                                                                                                 ,
      "           Mean Median  0.025    0.5"                                                                                                                                     ,
      "estimate -0.004 -0.002 -0.048 -0.002"                                                                                                                                     ,
      "PI       -0.004 -0.002 -0.048 -0.002"                                                                                                                                     ,
      ""                                                                                                                                                                         ,
      "Conditional adjusted effect estimate:"                                                                                                                                    ,
      "           Mean Median  0.025    0.5"                                                                                                                                     ,
      "estimate -0.004 -0.002 -0.048 -0.002"                                                                                                                                     ,
      "PI       -0.004 -0.002 -0.048 -0.002"                                                                                                                                     ,
      "The effect size estimates are summarized on the correlation scale and heterogeneity is summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))

})

# test posterior extraction
test_that("Posterior extraction works", {

  temp_summary <- summary(fits[["fit_15"]], conditional = TRUE)

  temp_samples <- extract_posterior(fits[["fit_15"]], parameter = "mu")
  expect_equal(mean(temp_samples), temp_summary$estimates["mu", "Mean"])

  temp_samples <- extract_posterior(fits[["fit_15"]], parameter = "mu", conditional = TRUE)
  expect_equal(mean(temp_samples), temp_summary$estimates_conditional["mu", "Mean"])

  temp_samples <- extract_posterior(fits[["fit_15"]], parameter = "tau")
  expect_equal(mean(temp_samples), temp_summary$estimates["tau", "Mean"])

  temp_samples <- extract_posterior(fits[["fit_15"]], parameter = "omega")
  expect_equal(unname(apply(temp_samples, 2, mean)), temp_summary$estimates[colnames(temp_samples), "Mean"])

  temp_samples <- extract_posterior(fits[["fit_15"]], parameter = "PET")
  expect_equal(mean(temp_samples), temp_summary$estimates["PET", "Mean"])

  temp_summary <- summary(fits[["fit_15"]], conditional = TRUE, output_scale = "OR")

  temp_samples <- extract_posterior(fits[["fit_15"]], parameter = "mu", output_scale = "OR")
  expect_equal(mean(temp_samples), temp_summary$estimates["mu", "Mean"])

  temp_samples <- extract_posterior(fits[["fit_15"]], parameter = "mu", conditional = TRUE, output_scale = "OR")
  expect_equal(mean(temp_samples), temp_summary$estimates_conditional["mu", "Mean"])
})

# test posterior extraction
test_that("True effects summary", {

  expect_equal(
    capture_output_lines(true_effects(fits[["fit_4"]], conditional = TRUE), print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "True effect estimates:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.192  0.045 -0.051 0.795"                                                              ,
      "estimate[2] 0.213  0.105 -0.016 0.797"                                                              ,
      "estimate[3] 0.266  0.227  0.000 0.795"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
      ""                                                                                                   ,
      "Conditional true effect estimates:"                                                                 ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.391  0.402 -0.146 0.889"                                                              ,
      "estimate[2] 0.405  0.413 -0.134 0.889"                                                              ,
      "estimate[3] 0.438  0.453 -0.071 0.883"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    ))

  expect_equal(
    capture_output_lines(true_effects(fits[["fit_15"]], conditional = TRUE), print = TRUE, width = 150),
    c(
      "Call:"                                                                                                  ,
      "RoBMA.reg(formula = ~mod_con, data = df_reg, priors = list(mod_con = list(null = prior(\"normal\", "    ,
      "    list(0, 0.05)), alt = prior(\"normal\", list(0.3, 0.15)))), "                                       ,
      "    priors_heterogeneity = NULL, priors_bias = list(prior_weightfunction(distribution = \"two.sided\", ",
      "        parameters = list(alpha = c(1, 1), steps = c(0.05)), "                                          ,
      "        prior_weights = 1/2), prior_PET(distribution = \"Cauchy\", "                                    ,
      "        parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2)), "                    ,
      "    priors_effect_null = NULL, algorithm = \"ss\", chains = 2, "                                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                       ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                        ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                      ,
      ""                                                                                                       ,
      "Robust Bayesian meta-regression"                                                                        ,
      "True effect estimates:"                                                                                 ,
      "               Mean Median  0.025  0.975"                                                               ,
      "estimate[1]  -1.585 -1.581 -1.731 -1.472"                                                               ,
      "estimate[2]  -1.329 -1.325 -1.457 -1.232"                                                               ,
      "estimate[3]  -1.088 -1.084 -1.202 -1.004"                                                               ,
      "estimate[4]  -0.858 -0.854 -0.961 -0.786"                                                               ,
      "estimate[5]  -0.637 -0.634 -0.737 -0.574"                                                               ,
      "estimate[6]  -0.423 -0.420 -0.517 -0.365"                                                               ,
      "estimate[7]  -0.214 -0.211 -0.305 -0.160"                                                               ,
      "estimate[8]  -0.007 -0.004 -0.097  0.047"                                                               ,
      "estimate[9]   0.200  0.203  0.112  0.255"                                                               ,
      "estimate[10]  0.409  0.412  0.317  0.468"                                                               ,
      "estimate[11]  0.622  0.626  0.529  0.689"                                                               ,
      "estimate[12]  0.842  0.846  0.742  0.918"                                                               ,
      "estimate[13]  1.071  1.075  0.964  1.160"                                                               ,
      "estimate[14]  1.312  1.315  1.193  1.413"                                                               ,
      "estimate[15]  1.567  1.569  1.432  1.686"                                                               ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."    ,
      ""                                                                                                       ,
      "Conditional true effect estimates:"                                                                     ,
      "               Mean Median  0.025  0.975"                                                               ,
      "estimate[1]  -1.585 -1.581 -1.731 -1.472"                                                               ,
      "estimate[2]  -1.329 -1.325 -1.457 -1.232"                                                               ,
      "estimate[3]  -1.088 -1.084 -1.202 -1.004"                                                               ,
      "estimate[4]  -0.858 -0.854 -0.961 -0.786"                                                               ,
      "estimate[5]  -0.637 -0.634 -0.737 -0.574"                                                               ,
      "estimate[6]  -0.423 -0.420 -0.517 -0.365"                                                               ,
      "estimate[7]  -0.214 -0.211 -0.305 -0.160"                                                               ,
      "estimate[8]  -0.007 -0.004 -0.097  0.047"                                                               ,
      "estimate[9]   0.200  0.203  0.112  0.255"                                                               ,
      "estimate[10]  0.409  0.412  0.317  0.468"                                                               ,
      "estimate[11]  0.622  0.626  0.529  0.689"                                                               ,
      "estimate[12]  0.842  0.846  0.742  0.918"                                                               ,
      "estimate[13]  1.071  1.075  0.964  1.160"                                                               ,
      "estimate[14]  1.312  1.315  1.193  1.413"                                                               ,
      "estimate[15]  1.567  1.569  1.432  1.686"                                                               ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    ))

  expect_equal(
    capture_output_lines(true_effects(fits[["fit_15"]], conditional = TRUE, output_scale = "logOR", probs = c(0.025, 0.5)), print = TRUE, width = 150),
    c(
      "Call:"                                                                                                  ,
      "RoBMA.reg(formula = ~mod_con, data = df_reg, priors = list(mod_con = list(null = prior(\"normal\", "    ,
      "    list(0, 0.05)), alt = prior(\"normal\", list(0.3, 0.15)))), "                                       ,
      "    priors_heterogeneity = NULL, priors_bias = list(prior_weightfunction(distribution = \"two.sided\", ",
      "        parameters = list(alpha = c(1, 1), steps = c(0.05)), "                                          ,
      "        prior_weights = 1/2), prior_PET(distribution = \"Cauchy\", "                                    ,
      "        parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2)), "                    ,
      "    priors_effect_null = NULL, algorithm = \"ss\", chains = 2, "                                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                       ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                        ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                      ,
      ""                                                                                                       ,
      "Robust Bayesian meta-regression"                                                                        ,
      "True effect estimates:"                                                                                 ,
      "               Mean Median  0.025    0.5"                                                               ,
      "estimate[1]  -2.875 -2.868 -3.139 -2.868"                                                               ,
      "estimate[2]  -2.411 -2.404 -2.642 -2.404"                                                               ,
      "estimate[3]  -1.973 -1.966 -2.180 -1.966"                                                               ,
      "estimate[4]  -1.556 -1.549 -1.744 -1.549"                                                               ,
      "estimate[5]  -1.156 -1.149 -1.337 -1.149"                                                               ,
      "estimate[6]  -0.768 -0.761 -0.938 -0.761"                                                               ,
      "estimate[7]  -0.389 -0.382 -0.554 -0.382"                                                               ,
      "estimate[8]  -0.013 -0.007 -0.175 -0.007"                                                               ,
      "estimate[9]   0.362  0.368  0.203  0.368"                                                               ,
      "estimate[10]  0.741  0.747  0.575  0.747"                                                               ,
      "estimate[11]  1.129  1.135  0.959  1.135"                                                               ,
      "estimate[12]  1.528  1.534  1.347  1.534"                                                               ,
      "estimate[13]  1.943  1.949  1.748  1.949"                                                               ,
      "estimate[14]  2.380  2.385  2.164  2.385"                                                               ,
      "estimate[15]  2.841  2.845  2.596  2.845"                                                               ,
      "The estimates are summarized on the log(OR) scale (priors were specified on the Cohen's d scale)."      ,
      ""                                                                                                       ,
      "Conditional true effect estimates:"                                                                     ,
      "               Mean Median  0.025    0.5"                                                               ,
      "estimate[1]  -2.875 -2.868 -3.139 -2.868"                                                               ,
      "estimate[2]  -2.411 -2.404 -2.642 -2.404"                                                               ,
      "estimate[3]  -1.973 -1.966 -2.180 -1.966"                                                               ,
      "estimate[4]  -1.556 -1.549 -1.744 -1.549"                                                               ,
      "estimate[5]  -1.156 -1.149 -1.337 -1.149"                                                               ,
      "estimate[6]  -0.768 -0.761 -0.938 -0.761"                                                               ,
      "estimate[7]  -0.389 -0.382 -0.554 -0.382"                                                               ,
      "estimate[8]  -0.013 -0.007 -0.175 -0.007"                                                               ,
      "estimate[9]   0.362  0.368  0.203  0.368"                                                               ,
      "estimate[10]  0.741  0.747  0.575  0.747"                                                               ,
      "estimate[11]  1.129  1.135  0.959  1.135"                                                               ,
      "estimate[12]  1.528  1.534  1.347  1.534"                                                               ,
      "estimate[13]  1.943  1.949  1.748  1.949"                                                               ,
      "estimate[14]  2.380  2.385  2.164  2.385"                                                               ,
      "estimate[15]  2.841  2.845  2.596  2.845"                                                               ,
      "The estimates are summarized on the log(OR) scale (priors were specified on the Cohen's d scale)."
    ))


})


#### creating / updating the test settings ####
if(FALSE){

  temp_fits_dir <- Sys.getenv("ROBMA_TEST_FITS_DIR")

  saved_files <- paste0("fit_", 1:16, ".RDS")
  fits <- list()
  for (i in 1:16) {
    fits[[i]] <- readRDS(file = file.path(temp_fits_dir, saved_files[i]))
  }
  names(fits) <- paste0("fit_", 1:16)

  # generate print files
  for(i in seq_along(fits)){
    write.table(capture_output_lines(fits[[i]], print = TRUE, width = 150), file = file.path("tests/results/print", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary files
  for(i in seq_along(fits)){
    write.table(capture_output_lines(summary(fits[[i]]), print = TRUE, width = 150), file = file.path("tests/results/summary", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.models files
  for(i in seq_along(fits)){
    write.table(capture_output_lines(summary(fits[[i]], type = "models"), print = TRUE, width = 150), file = file.path("tests/results/summary.models", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.diagnostics files
  for(i in seq_along(fits)){
    write.table(capture_output_lines(summary(fits[[i]], type = "diagnostics"), print = TRUE, width = 200), file = file.path("tests/results/summary.diagnostics", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.individual files
  for(i in seq_along(fits)){
    write.table(capture_output_lines(summary(fits[[i]], type = "individual"), print = TRUE, width = 150), file = file.path("tests/results/summary.individual", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary.individual files
  for(i in seq_along(fits)){
    write.table(gsub("(.{80})", "\\1\\\n", interpret(fits[[i]])), file = file.path("tests/results/interpret", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary_heterogeneity files
  for(i in seq_along(fits)){
    write.table(capture_output_lines(summary_heterogeneity(fits[[i]]), print = TRUE, width = 150), file = file.path("tests/results/summary_heterogeneity", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }

  # generate summary_heterogeneity.individual files
  for(i in seq_along(fits)){
    write.table(capture_output_lines(summary_heterogeneity(fits[[i]], type = "i"), print = TRUE, width = 150), file = file.path("tests/results/summary_heterogeneity.individual", paste0(i, ".txt")), row.names = FALSE, col.names = FALSE)
  }
}
