context("(10) Prediction functions")
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
      "             Mean Median  0.025 0.975"                                                              ,
      "estimate[1] 0.182  0.000 -0.103 0.837"                                                              ,
      "estimate[2] 0.182  0.000 -0.103 0.837"                                                              ,
      "estimate[3] 0.182  0.000 -0.103 0.837"                                                              ,
      "The estimates are summarized on the Cohen's d scale (priors were specified on the Cohen's d scale)."
    )
  )

  expect_equal(
    capture_output_lines(pred5, print = TRUE, width = 150),
    c(
    )
  )

  expect_equal(
    capture_output_lines(pred6, print = TRUE, width = 150),
    c(
    )
  )

  expect_equal(
    capture_output_lines(pred7, print = TRUE, width = 150),
    c(
    )
  )

  expect_equal(
    capture_output_lines(pred8, print = TRUE, width = 150),
    c(
    )
  )

})
