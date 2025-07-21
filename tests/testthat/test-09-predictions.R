context("(9) Prediction functions")
skip_on_cran()

### Read all prefitted objects
# test objects - assuming that the fit function worked properly
temp_fits_dir <- Sys.getenv("ROBMA_TEST_FITS_DIR")
if (temp_fits_dir == "" || !dir.exists(temp_fits_dir)) {
  stop("Temporary fits directory not found. Run test-4-fit.R first.")
}

saved_files <- paste0("fit_", 1:18, ".RDS")
fits <- list()
for (i in 1:18) {
  fits[[i]] <- readRDS(file = file.path(temp_fits_dir, saved_files[i]))
}
names(fits) <- paste0("fit_", 1:18)


test_that("Meta-analysis prediction", {

  set.seed(1)
  pred1 <- predict(fits[["fit_4"]], incorporate_publication_bias = FALSE)
  pred2 <- predict(fits[["fit_4"]], incorporate_publication_bias = TRUE)
  pred3 <- predict(fits[["fit_4"]], incorporate_publication_bias = FALSE, output_scale = "r")
  pred4 <- predict(fits[["fit_4"]], incorporate_publication_bias = FALSE, type = "terms")
  pred5 <- predict(fits[["fit_4"]], incorporate_publication_bias = FALSE, type = "effect")
  pred6 <- predict(fits[["fit_4"]], incorporate_publication_bias = FALSE, newdata = list(r = c(0, 0.1, 0.2, 0.5), se = c(0.25, 0.25, 0.25, 0.25)))
  pred7 <- predict(fits[["fit_4"]], incorporate_publication_bias = TRUE , newdata = list(r = c(0, 0.1, 0.2, 0.5), se = c(0.25, 0.25, 0.25, 0.25)))
  pred8 <- predict(fits[["fit_4"]], incorporate_publication_bias = FALSE, newdata = list(r = c(0, 0.1, 0.2, 0.5), se = c(0.25, 0.25, 0.25, 0.25)), type = "terms")
  pred9 <- predict(fits[["fit_4"]], incorporate_publication_bias = TRUE, type = "terms")
  pred10<- predict(fits[["fit_4"]], incorporate_publication_bias = TRUE, type = "effect")

  expect_equal(
    capture_output_lines(pred1, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.186  0.171 -1.244 1.648"                                                              ,
      "estimate[2] 0.184  0.165 -1.058 1.427"                                                              ,
      "estimate[3] 0.175  0.153 -0.690 1.125"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )


  expect_equal(
    capture_output_lines(pred2, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.404  0.378 -1.135 1.982"                                                              ,
      "estimate[2] 0.368  0.340 -0.922 1.672"                                                              ,
      "estimate[3] 0.267  0.256 -0.666 1.188"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred3, print = TRUE, width = 150),
    c(
      "Call:"                                                                                                                                                                    ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                                                                                              ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                                                                                         ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                                                                                          ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                                                                                        ,
      ""                                                                                                                                                                         ,
      "Robust Bayesian meta-analysis"                                                                                                                                            ,
      "Posterior predictions:"                                                                                                                                                   ,
      "             Mean Median  0.025 0.975"                                                                                                                                    ,
      "estimate[1] 0.079  0.087 -0.514 0.621"                                                                                                                                    ,
      "estimate[2] 0.085  0.089 -0.447 0.592"                                                                                                                                    ,
      "estimate[3] 0.083  0.076 -0.341 0.491"                                                                                                                                    ,
      "The effect size estimates are summarized on the correlation scale and heterogeneity is summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred4, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "       Mean Median  0.025 0.975"                                                              ,
      "mu[1] 0.182  0.000 -0.103 0.837"                                                              ,
      "mu[2] 0.182  0.000 -0.103 0.837"                                                              ,
      "mu[3] 0.182  0.000 -0.103 0.837"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred5, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "          Mean Median  0.025 0.975"                                                              ,
      "theta[1] 0.183  0.060 -0.456 0.906"                                                              ,
      "theta[2] 0.183  0.062 -0.457 0.915"                                                              ,
      "theta[3] 0.178  0.058 -0.449 0.917"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred6, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.189  0.179 -1.138 1.626"                                                              ,
      "estimate[2] 0.176  0.164 -1.194 1.571"                                                              ,
      "estimate[3] 0.192  0.191 -1.234 1.703"                                                              ,
      "estimate[4] 0.211  0.214 -1.966 2.399"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred7, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.411  0.386 -1.055 1.881"                                                              ,
      "estimate[2] 0.406  0.377 -1.114 1.894"                                                              ,
      "estimate[3] 0.384  0.357 -1.135 1.928"                                                              ,
      "estimate[4] 0.587  0.461 -1.852 3.400"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred8, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "       Mean Median  0.025 0.975"                                                              ,
      "mu[1] 0.182  0.000 -0.103 0.837"                                                              ,
      "mu[2] 0.182  0.000 -0.103 0.837"                                                              ,
      "mu[3] 0.182  0.000 -0.103 0.837"                                                              ,
      "mu[4] 0.182  0.000 -0.103 0.837"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred9, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "       Mean Median  0.025 0.975"                                                              ,
      "mu[1] 0.285  0.196 -0.019 1.024"                                                              ,
      "mu[2] 0.264  0.180 -0.020 0.914"                                                              ,
      "mu[3] 0.226  0.125 -0.044 0.844"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred10, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, "                        ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                   ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                    ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "          Mean Median  0.025 0.975"                                                              ,
      "theta[1] 0.298  0.247 -0.311 1.119"                                                              ,
      "theta[2] 0.271  0.227 -0.385 1.032"                                                              ,
      "theta[3] 0.240  0.184 -0.357 0.954"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

})

test_that("Meta-regression prediction", {

  set.seed(1)
  pred1 <- predict(fits[["fit_15"]], incorporate_publication_bias = FALSE)
  pred2 <- predict(fits[["fit_15"]], incorporate_publication_bias = TRUE)
  pred3 <- predict(fits[["fit_15"]], incorporate_publication_bias = FALSE, output_scale = "r")
  pred4 <- predict(fits[["fit_15"]], incorporate_publication_bias = FALSE, type = "terms")  # no heterogeneity = should be the same
  pred5 <- predict(fits[["fit_15"]], incorporate_publication_bias = FALSE, type = "effect") # no heterogeneity = should be the same
  newdf <-  data.frame(
    d       = 0,
    se      = 1,
    mod_con = 2
  )
  pred6 <- predict(fits[["fit_15"]], incorporate_publication_bias = FALSE, newdata = newdf)
  pred7 <- predict(fits[["fit_15"]], incorporate_publication_bias = TRUE , newdata = newdf)
  pred8 <- predict(fits[["fit_15"]], incorporate_publication_bias = FALSE, newdata = newdf, type = "terms")
  pred9 <- predict(fits[["fit_15"]], incorporate_publication_bias = FALSE, newdata = newdf, type = "effect")
  pred10<- predict(fits[["fit_15"]], conditional = TRUE , newdata = newdf)
  pred11<- predict(fits[["fit_15"]], incorporate_publication_bias = TRUE, type = "terms")  # no heterogeneity = should be the same
  pred12<- predict(fits[["fit_15"]], incorporate_publication_bias = TRUE, type = "effect") # no heterogeneity = should be the same

  expect_equal(
    capture_output_lines(pred1, print = TRUE, width = 150),
    c(
      "Call:",
      "RoBMA.reg(formula = ~mod_con, data = df_reg, priors = list(mod_con = list(null = prior(\"normal\", ",
      "    list(0, 0.05)), alt = prior(\"normal\", list(0.3, 0.15)))), ",
      "    priors_heterogeneity = NULL, priors_bias = list(prior_weightfunction(distribution = \"two.sided\", ",
      "        parameters = list(alpha = c(1, 1), steps = c(0.05)), ",
      "        prior_weights = 1/2), prior_PET(distribution = \"Cauchy\", ",
      "        parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2)), ",
      "    priors_effect_null = NULL, algorithm = \"ss\", chains = 2, ",
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, ",
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, ",
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)",
      "",
      "Robust Bayesian meta-regression",
      "Posterior predictions:",
      "               Mean Median  0.025  0.975",
      "estimate[1]  -1.587 -1.585 -1.863 -1.340",
      "estimate[2]  -1.332 -1.328 -1.584 -1.101",
      "estimate[3]  -1.091 -1.090 -1.323 -0.868",
      "estimate[4]  -0.857 -0.854 -1.085 -0.636",
      "estimate[5]  -0.635 -0.630 -0.855 -0.420",
      "estimate[6]  -0.425 -0.424 -0.640 -0.212",
      "estimate[7]  -0.213 -0.213 -0.425 -0.007",
      "estimate[8]  -0.010 -0.008 -0.219  0.192",
      "estimate[9]   0.201  0.201 -0.013  0.413",
      "estimate[10]  0.408  0.410  0.191  0.617",
      "estimate[11]  0.621  0.622  0.404  0.843",
      "estimate[12]  0.843  0.845  0.618  1.058",
      "estimate[13]  1.072  1.073  0.842  1.303",
      "estimate[14]  1.311  1.312  1.077  1.545",
      "estimate[15]  1.570  1.569  1.317  1.824",
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred2, print = TRUE, width = 150),
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
      "Posterior predictions:"                                                                                 ,
      "               Mean Median  0.025  0.975"                                                               ,
      "estimate[1]  -1.576 -1.576 -1.820 -1.336"                                                               ,
      "estimate[2]  -1.325 -1.325 -1.558 -1.099"                                                               ,
      "estimate[3]  -1.081 -1.078 -1.307 -0.866"                                                               ,
      "estimate[4]  -0.847 -0.845 -1.064 -0.639"                                                               ,
      "estimate[5]  -0.630 -0.627 -0.844 -0.427"                                                               ,
      "estimate[6]  -0.417 -0.416 -0.626 -0.213"                                                               ,
      "estimate[7]  -0.217 -0.218 -0.419 -0.017"                                                               ,
      "estimate[8]   0.000 -0.001 -0.224  0.226"                                                               ,
      "estimate[9]   0.220  0.223  0.014  0.418"                                                               ,
      "estimate[10]  0.418  0.417  0.218  0.630"                                                               ,
      "estimate[11]  0.631  0.632  0.420  0.838"                                                               ,
      "estimate[12]  0.851  0.849  0.635  1.074"                                                               ,
      "estimate[13]  1.079  1.077  0.861  1.305"                                                               ,
      "estimate[14]  1.321  1.319  1.094  1.552"                                                               ,
      "estimate[15]  1.577  1.573  1.342  1.826"                                                               ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred3, print = TRUE, width = 150),
    c(
      "Call:"                                                                                                                                                                    ,
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
      "Posterior predictions:"                                                                                                                                                   ,
      "               Mean Median  0.025  0.975"                                                                                                                                 ,
      "estimate[1]  -0.620 -0.621 -0.681 -0.558"                                                                                                                                 ,
      "estimate[2]  -0.553 -0.553 -0.623 -0.482"                                                                                                                                 ,
      "estimate[3]  -0.476 -0.476 -0.553 -0.396"                                                                                                                                 ,
      "estimate[4]  -0.393 -0.394 -0.477 -0.308"                                                                                                                                 ,
      "estimate[5]  -0.303 -0.303 -0.397 -0.208"                                                                                                                                 ,
      "estimate[6]  -0.207 -0.206 -0.307 -0.109"                                                                                                                                 ,
      "estimate[7]  -0.106 -0.106 -0.208 -0.006"                                                                                                                                 ,
      "estimate[8]  -0.003 -0.002 -0.104  0.098"                                                                                                                                 ,
      "estimate[9]   0.099  0.099 -0.008  0.200"                                                                                                                                 ,
      "estimate[10]  0.199  0.199  0.100  0.297"                                                                                                                                 ,
      "estimate[11]  0.296  0.297  0.201  0.388"                                                                                                                                 ,
      "estimate[12]  0.386  0.388  0.297  0.469"                                                                                                                                 ,
      "estimate[13]  0.471  0.472  0.392  0.546"                                                                                                                                 ,
      "estimate[14]  0.547  0.548  0.473  0.614"                                                                                                                                 ,
      "estimate[15]  0.616  0.618  0.551  0.674"                                                                                                                                 ,
      "The effect size estimates are summarized on the correlation scale and heterogeneity is summarized on the Fisher's z scale (priors were specified on the Cohen's d scale)."
    ))

  expect_equal(
    capture_output_lines(pred4, print = TRUE, width = 150),
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
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                       ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                     ,
      ""                                                                                                      ,
      "Robust Bayesian meta-regression"                                                                       ,
      "Posterior predictions:"                                                                                ,
      "         Mean Median  0.025  0.975"                                                              ,
      "mu[1]  -1.585 -1.581 -1.731 -1.472"                                                              ,
      "mu[2]  -1.329 -1.325 -1.457 -1.232"                                                              ,
      "mu[3]  -1.088 -1.084 -1.202 -1.004"                                                              ,
      "mu[4]  -0.858 -0.854 -0.961 -0.786"                                                              ,
      "mu[5]  -0.637 -0.634 -0.737 -0.574"                                                              ,
      "mu[6]  -0.423 -0.420 -0.517 -0.365"                                                              ,
      "mu[7]  -0.214 -0.211 -0.305 -0.160"                                                              ,
      "mu[8]  -0.007 -0.004 -0.097  0.047"                                                              ,
      "mu[9]   0.200  0.203  0.112  0.255"                                                              ,
      "mu[10]  0.409  0.412  0.317  0.468"                                                              ,
      "mu[11]  0.622  0.626  0.529  0.689"                                                              ,
      "mu[12]  0.842  0.846  0.742  0.918"                                                              ,
      "mu[13]  1.071  1.075  0.964  1.160"                                                              ,
      "mu[14]  1.312  1.315  1.193  1.413"                                                              ,
      "mu[15]  1.567  1.569  1.432  1.686"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred5, print = TRUE, width = 150),
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
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                       ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                     ,
      ""                                                                                                      ,
      "Robust Bayesian meta-regression"                                                                       ,
      "Posterior predictions:"                                                                                ,
      "            Mean Median  0.025  0.975"                                                              ,
      "theta[1]  -1.585 -1.581 -1.731 -1.472"                                                              ,
      "theta[2]  -1.329 -1.325 -1.457 -1.232"                                                              ,
      "theta[3]  -1.088 -1.084 -1.202 -1.004"                                                              ,
      "theta[4]  -0.858 -0.854 -0.961 -0.786"                                                              ,
      "theta[5]  -0.637 -0.634 -0.737 -0.574"                                                              ,
      "theta[6]  -0.423 -0.420 -0.517 -0.365"                                                              ,
      "theta[7]  -0.214 -0.211 -0.305 -0.160"                                                              ,
      "theta[8]  -0.007 -0.004 -0.097  0.047"                                                              ,
      "theta[9]   0.200  0.203  0.112  0.255"                                                              ,
      "theta[10]  0.409  0.412  0.317  0.468"                                                              ,
      "theta[11]  0.622  0.626  0.529  0.689"                                                              ,
      "theta[12]  0.842  0.846  0.742  0.918"                                                              ,
      "theta[13]  1.071  1.075  0.964  1.160"                                                              ,
      "theta[14]  1.312  1.315  1.193  1.413"                                                              ,
      "theta[15]  1.567  1.569  1.432  1.686"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred6, print = TRUE, width = 150),
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
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                       ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                     ,
      ""                                                                                                      ,
      "Robust Bayesian meta-regression"                                                                       ,
      "Posterior predictions:"                                                                                ,
      "             Mean Median  0.025  0.975"                                                                ,
      "estimate[1] 3.467  2.093 -2.463 18.218"                                                                ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred7, print = TRUE, width = 150),
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
      "Posterior predictions:"                                                                                 ,
      "             Mean Median  0.025  0.975"                                                                 ,
      "estimate[1] 4.617  2.815 -2.413 22.111"                                                                 ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred8, print = TRUE, width = 150),
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
      "Posterior predictions:"                                                                                 ,
      "       Mean Median 0.025 0.975"                                                                   ,
      "mu[1] 2.112  2.113 1.938 2.274"                                                                   ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred9, print = TRUE, width = 150),
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
      "Posterior predictions:"                                                                                 ,
      "          Mean Median 0.025 0.975"                                                                   ,
      "theta[1] 2.112  2.113 1.938 2.274"                                                                   ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred10, print = TRUE, width = 150),
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
      "Posterior predictions:"                                                                                 ,
      "             Mean Median  0.025  0.975"                                                                 ,
      "estimate[1] 4.802  2.953 -2.453 23.402"                                                                 ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."    ,
      ""                                                                                                       ,
      "Conditional posterior predictions:"                                                                     ,
      "             Mean Median  0.025  0.975"                                                                 ,
      "estimate[1] 4.802  2.953 -2.453 23.402"                                                                 ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred11, print = TRUE, width = 150),
    c(
      "Call:",
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
      "Posterior predictions:"                                                                                 ,
      "         Mean Median  0.025  0.975"                                                               ,
      "mu[1]  -1.576 -1.576 -1.687 -1.469"                                                               ,
      "mu[2]  -1.321 -1.321 -1.414 -1.230"                                                               ,
      "mu[3]  -1.079 -1.079 -1.160 -1.001"                                                               ,
      "mu[4]  -0.850 -0.849 -0.918 -0.783"                                                               ,
      "mu[5]  -0.629 -0.629 -0.689 -0.571"                                                               ,
      "mu[6]  -0.415 -0.416 -0.468 -0.362"                                                               ,
      "mu[7]  -0.206 -0.207 -0.255 -0.156"                                                               ,
      "mu[8]   0.001  0.000 -0.048  0.049"                                                               ,
      "mu[9]   0.208  0.207  0.157  0.257"                                                               ,
      "mu[10]  0.417  0.417  0.363  0.470"                                                               ,
      "mu[11]  0.630  0.630  0.570  0.690"                                                               ,
      "mu[12]  0.850  0.850  0.782  0.920"                                                               ,
      "mu[13]  1.080  1.080  1.001  1.161"                                                               ,
      "mu[14]  1.320  1.320  1.228  1.416"                                                               ,
      "mu[15]  1.575  1.574  1.466  1.689"                                                               ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred12, print = TRUE, width = 150),
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
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                       ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                     ,
      ""                                                                                                      ,
      "Robust Bayesian meta-regression"                                                                       ,
      "Posterior predictions:"                                                                                ,
      "            Mean Median  0.025  0.975"                                                              ,
      "theta[1]  -1.576 -1.576 -1.687 -1.469"                                                              ,
      "theta[2]  -1.321 -1.321 -1.414 -1.230"                                                              ,
      "theta[3]  -1.079 -1.079 -1.160 -1.001"                                                              ,
      "theta[4]  -0.850 -0.849 -0.918 -0.783"                                                              ,
      "theta[5]  -0.629 -0.629 -0.689 -0.571"                                                              ,
      "theta[6]  -0.415 -0.416 -0.468 -0.362"                                                              ,
      "theta[7]  -0.206 -0.207 -0.255 -0.156"                                                              ,
      "theta[8]   0.001  0.000 -0.048  0.049"                                                              ,
      "theta[9]   0.208  0.207  0.157  0.257"                                                              ,
      "theta[10]  0.417  0.417  0.363  0.470"                                                              ,
      "theta[11]  0.630  0.630  0.570  0.690"                                                              ,
      "theta[12]  0.850  0.850  0.782  0.920"                                                              ,
      "theta[13]  1.080  1.080  1.001  1.161"                                                              ,
      "theta[14]  1.320  1.320  1.228  1.416"                                                              ,
      "theta[15]  1.575  1.574  1.466  1.689"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )
})

test_that("Multilevel Meta-analysis prediction", {

  set.seed(1)
  pred1 <- predict(fits[["fit_18"]], incorporate_publication_bias = FALSE)
  pred2 <- predict(fits[["fit_18"]], incorporate_publication_bias = TRUE)
  pred3 <- predict(fits[["fit_18"]], incorporate_publication_bias = FALSE, conditional = TRUE)
  pred4 <- predict(fits[["fit_18"]], incorporate_publication_bias = FALSE, type = "terms")
  pred5 <- predict(fits[["fit_18"]], incorporate_publication_bias = FALSE, type = "effect")
  pred6 <- predict(fits[["fit_18"]], incorporate_publication_bias = FALSE, newdata = list(z = c(0, 0.1), se = c(0.25, 0.25)))
  pred7 <- predict(fits[["fit_18"]], incorporate_publication_bias = TRUE , newdata = list(z = c(0, 0.1), se = c(0.25, 0.25)))
  pred8 <- predict(fits[["fit_18"]], incorporate_publication_bias = FALSE, newdata = list(z = c(0, 0.1), se = c(0.25, 0.25)), type = "terms")
  pred9 <- predict(fits[["fit_18"]], incorporate_publication_bias = TRUE, type = "terms")
  pred10<- predict(fits[["fit_18"]], incorporate_publication_bias = TRUE, type = "effect")

  expect_equal(
    capture_output_lines(pred1, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.203  0.190 -1.085 1.670"                                                              ,
      "estimate[2] 0.179  0.155 -0.995 1.434"                                                              ,
      "estimate[3] 0.183  0.172 -0.741 1.159"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )


  expect_equal(
    capture_output_lines(pred2, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.357  0.380 -1.211 1.853"                                                              ,
      "estimate[2] 0.377  0.383 -0.927 1.694"                                                              ,
      "estimate[3] 0.289  0.290 -0.508 1.139"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred3, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.158  0.187 -1.103 1.488"                                                              ,
      "estimate[2] 0.182  0.115 -1.122 1.446"                                                              ,
      "estimate[3] 0.189  0.154 -0.848 1.277"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale).",
      ""                                                                                                   ,
      "Conditional posterior predictions:"                                                                 ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.353  0.323 -0.744 1.752"                                                              ,
      "estimate[2] 0.400  0.378 -0.845 1.476"                                                              ,
      "estimate[3] 0.445  0.412 -0.511 1.398"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    ))

  expect_equal(
    capture_output_lines(pred4, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "       Mean Median  0.025 0.975"                                                              ,
      "mu[1] 0.179  0.000 -0.097 0.779"                                                              ,
      "mu[2] 0.179  0.000 -0.097 0.779"                                                              ,
      "mu[3] 0.179  0.000 -0.097 0.779"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred5, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "          Mean Median  0.025 0.975"                                                              ,
      "theta[1] 0.171  0.067 -0.514 0.916"                                                              ,
      "theta[2] 0.183  0.080 -0.436 0.926"                                                              ,
      "theta[3] 0.169  0.084 -0.612 0.908"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred6, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.152  0.147 -1.056 1.344"                                                              ,
      "estimate[2] 0.182  0.136 -1.185 1.401"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred7, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.343  0.290 -0.872 1.763"                                                              ,
      "estimate[2] 0.279  0.274 -0.846 1.459"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred8, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "       Mean Median  0.025 0.975"                                                              ,
      "mu[1] 0.179  0.000 -0.097 0.779"                                                              ,
      "mu[2] 0.179  0.000 -0.097 0.779"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred9, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "       Mean Median  0.025 0.975"                                                              ,
      "mu[1] 0.279  0.215 -0.014 0.959"                                                              ,
      "mu[2] 0.258  0.193 -0.031 0.854"                                                              ,
      "mu[3] 0.221  0.151 -0.040 0.788"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred10, print = TRUE, width = 150),
    c(
      "Call:"                                                                                              ,
      "RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                               ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                                ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "   ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                  ,
      ""                                                                                                   ,
      "Robust Bayesian meta-analysis"                                                                      ,
      "Posterior predictions:"                                                                             ,
      "          Mean Median  0.025 0.975"                                                              ,
      "theta[1] 0.290  0.264 -0.285 0.975"                                                              ,
      "theta[2] 0.270  0.252 -0.268 0.882"                                                              ,
      "theta[3] 0.272  0.240 -0.184 0.924"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

})
