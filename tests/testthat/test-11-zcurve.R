context("(11) Z-curve Diagnostics")
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


test_that("zcurve functions function works", {

  # meta-analysis
  zfit4 <- as_zcurve(fits[[4]])
  zfit6 <- as_zcurve(fits[[6]]) # flips to negative effects

  # meta-regression
  zfit15 <- as_zcurve(fits[[15]])

  # hierarchical
  zfit18 <- as_zcurve(fits[[18]])

  # check summary functions
  expect_equal(
    capture_output_lines(suppressWarnings(summary(zfit4)), print = TRUE, width = 150),
    c(
      "Call:"                                                                                 ,
      "as_zcurve: RoBMA(r = r, n = n, model_type = \"PSMA\", algorithm = \"ss\", chains = 2, ",
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                      ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "       ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                     ,
      ""                                                                                      ,
      "Z-curve"                                                                               ,
      "Model-averaged estimates:"                                                             ,
      "           Mean Median 0.025  0.975"                                                   ,
      "EDR       0.158  0.098 0.050  0.536"                                                   ,
      "Soric FDR 0.545  0.484 0.046  1.000"                                                   ,
      "Missing N 2.313  0.000 0.000 17.343"                                                   ,
      "Estimated using 3 estimates, 1 significant (ODR = 0.33, 95 CI [0.02, 0.87])."
    ))

  expect_equal(
    capture_output_lines(suppressWarnings(summary(zfit6, conditional = TRUE)), print = TRUE, width = 150),
    c(
      "Call:"                                                                                                        ,
      "as_zcurve: RoBMA(d = -d, se = d_se, effect_direction = \"negative\", model_type = \"PSMA\", "                 ,
      "    algorithm = \"ss\", chains = 2, sample = 2500, burnin = 1000, "                                           ,
      "    adapt = 500, parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, ",
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                            ,
      ""                                                                                                             ,
      "Z-curve"                                                                                                      ,
      "Model-averaged estimates:"                                                                                    ,
      "           Mean Median 0.025  0.975"                                                                          ,
      "EDR       0.156  0.094 0.050  0.511"                                                                          ,
      "Soric FDR 0.552  0.507 0.050  1.000"                                                                          ,
      "Missing N 2.262  0.000 0.000 18.736"                                                                          ,
      "Estimated using 3 estimates, 1 significant (ODR = 0.33, 95 CI [0.02, 0.87])."                                 ,
      ""                                                                                                             ,
      "Conditional estimates:"                                                                                       ,
      "           Mean Median 0.025 0.975"                                                                           ,
      "EDR       0.226  0.196 0.052 0.538"                                                                           ,
      "Soric FDR 0.301  0.216 0.045 0.968"                                                                           ,
      "Missing N 1.189  0.000 0.000 9.469"                                                                           ,
      "Estimated using 3 estimates, 1 significant (ODR = 0.33, 95 CI [0.02, 0.87])."
    ))

  expect_equal(
    capture_output_lines(suppressWarnings(summary(zfit15)), print = TRUE, width = 150),
    c(
      "Call:"                                                                                                         ,
      "as_zcurve: RoBMA.reg(formula = ~mod_con, data = df_reg, priors = list(mod_con = list(null = prior(\"normal\", ",
      "    list(0, 0.05)), alt = prior(\"normal\", list(0.3, 0.15)))), "                                              ,
      "    priors_heterogeneity = NULL, priors_bias = list(prior_weightfunction(distribution = \"two.sided\", "       ,
      "        parameters = list(alpha = c(1, 1), steps = c(0.05)), "                                                 ,
      "        prior_weights = 1/2), prior_PET(distribution = \"Cauchy\", "                                           ,
      "        parameters = list(0, 1), truncation = list(0, Inf), prior_weights = 1/2)), "                           ,
      "    priors_effect_null = NULL, algorithm = \"ss\", chains = 2, "                                               ,
      "    sample = 2500, burnin = 1000, adapt = 500, parallel = TRUE, "                                              ,
      "    autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, "                               ,
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                                             ,
      ""                                                                                                              ,
      "Z-curve"                                                                                                       ,
      "Model-averaged estimates:"                                                                                     ,
      "           Mean Median 0.025 0.975"                                                                            ,
      "EDR       0.834  0.871 0.616 0.879"                                                                            ,
      "Soric FDR 0.011  0.008 0.007 0.033"                                                                            ,
      "Missing N 1.181  0.000 0.000 8.326"                                                                            ,
      "Estimated using 15 estimates, 14 significant (ODR = 0.93, 95 CI [0.66, 1.00])."
    )
  )

  expect_equal(
    capture_output_lines(suppressWarnings(summary(zfit18)), print = TRUE, width = 150),
    c(
      "Call:"                                                                                           ,
      "as_zcurve: RoBMA(d = d, se = d_se, study_ids = c(1, 1, 2), algorithm = \"ss\", "                 ,
      "    chains = 1, sample = 500, burnin = 250, adapt = 100, thin = 2, "                             ,
      "    parallel = TRUE, autofit = FALSE, convergence_checks = set_convergence_checks(max_Rhat = 2, ",
      "        min_ESS = 10, max_error = 1, max_SD_error = 1), seed = 1)"                               ,
      ""                                                                                                ,
      "Z-curve"                                                                                         ,
      "Model-averaged estimates:"                                                                       ,
      "           Mean Median 0.025  0.975"                                                             ,
      "EDR       0.146  0.089 0.050  0.472"                                                             ,
      "Soric FDR 0.558  0.538 0.059  1.000"                                                             ,
      "Missing N 2.023  0.000 0.000 14.278"                                                             ,
      "Estimated using 3 estimates, 1 significant (ODR = 0.33, 95 CI [0.02, 0.87])."
    )
  )

  # check plotting functions
  vdiffr::expect_doppelganger(paste0("plot_zcurve_",1), function() plot(zfit4, by.hist = 0.5))
  vdiffr::expect_doppelganger(paste0("plot_zcurve_",2), function() plot(zfit4, by.hist = 1))
  vdiffr::expect_doppelganger(paste0("plot_zcurve_",3), function() plot(zfit6, by.hist = 1))
  vdiffr::expect_doppelganger(paste0("plot_zcurve_",4), function() plot(zfit6, plot_fit = FALSE))
  vdiffr::expect_doppelganger(paste0("plot_zcurve_",5), function() plot(zfit6, by.hist = 0.5, conditional = TRUE, plot_CI = FALSE, plot_extrapolation = FALSE, plot_thresholds = FALSE))
  vdiffr::expect_doppelganger(paste0("plot_zcurve_",6), suppressMessages(plot(zfit15, by.hist = 0.5, plot_type = "ggplot")))

  vdiffr::expect_doppelganger(paste0("plot_zcurve_",7), function() hist(zfit4))
  vdiffr::expect_doppelganger(paste0("plot_zcurve_",8), function() {
    hist(zfit4, by = 1)
    lines(zfit4, col = "green")
  })
  vdiffr::expect_doppelganger(paste0("plot_zcurve_",9), function() {
    out <- hist(zfit4, plot_type = "ggplot")
  })
  vdiffr::expect_doppelganger(paste0("plot_zcurve_",10), function() {
    out <- hist(zfit4, plot_type = "ggplot")
    out <- out + lines(zfit4, plot_type = "ggplot", extrapolate = TRUE, col = "green")
    out
  })

})
