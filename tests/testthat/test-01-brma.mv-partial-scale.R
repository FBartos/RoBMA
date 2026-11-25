context("Live fit of a known-V brma.mv with a scale formula for one random component")

source(testthat::test_path("common-functions.R"))
skip_on_cran()

# A live, uncached fit: the model compiles only when every formula syntax
# that BayesTools::JAGS_fit() joins ends with a complete line (BayesTools
# 0.3.1.128). The location formula ends with the standard deviation of the
# marginalized 'effect' component, followed by the scale formula of 'study';
# with BayesTools 0.3.1.127 JAGS stopped with a syntax error. The syntax
# itself is checked without fitting in
# test-00-input-data-mv-02-random-formulas.R.
test_that("known-V brma.mv with a scale formula for one of two random components fits", {

  dat <- data.frame(
    yi     = c(0.08, 0.13, 0.18, 0.20, 0.01, 0.05),
    study  = rep(c("s1", "s2", "s3"), each = 2L),
    effect = rep(c("a", "b"), 3L),
    x      = c(0, 1, 0, 1, 0, 1)
  )
  V <- kronecker(diag(3L), matrix(c(0.04, 0.018, 0.018, 0.05), nrow = 2L))

  fit <- suppressWarnings(brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = list(study = ~ 1 | study, effect = ~ 1 | study:effect),
    scale                     = list(study = ~ x),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    chains                    = 1,
    sample                    = 500,
    burnin                    = 200,
    adapt                     = 200,
    seed                      = 1,
    silent                    = TRUE,
    convergence_checks        = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
  ))

  expect_s3_class(fit, "brma.mv")
  expect_true(isTRUE(fit[["fit"]][["has_posterior"]]))
  expect_identical(
    .data_scale_formula_parameters(fit[["data"]]),
    c(study = "log_tau_study")
  )
  # The compiled model has the standard deviation of the marginalized
  # component on its own line, followed by the scale formula.
  model <- strsplit(
    paste(as.character(fit[["fit"]][["model"]]), collapse = "\n"),
    "\n", fixed = TRUE
  )[[1L]]
  model   <- trimws(model[nzchar(trimws(model))])
  sd_line <- which(model == "mu__xREx__effect_xRE_STDx[1] = mu__xREx__effect_intercept")
  expect_length(sd_line, 1L)
  expect_identical(model[sd_line + 1L], "for(i in 1:N_log_tau_study){")
  draws <- as.matrix(.get_posterior_samples(fit[["fit"]]))
  expect_true(all(c("log_tau_study_intercept", "log_tau_study_x",
                    "mu__xREx__effect_intercept") %in% colnames(draws)))
  expect_true(all(is.finite(draws[, c("log_tau_study_intercept", "log_tau_study_x")])))

  # The summary reports the scale formula of 'study' and the heterogeneity
  # of the marginalized 'effect' component.
  fit_summary <- summary(fit)
  expect_s3_class(fit_summary, "summary.brma")
  summary_table <- as.data.frame(fit_summary)
  expect_setequal(
    summary_table[["parameter"]][summary_table[["component"]] == "scale"],
    c("(study: tau) exp(intercept)", "(study: tau) x")
  )
  expect_true("effect: tau" %in% summary_table[["parameter"]])
  expect_true(all(is.finite(summary_table[["Mean"]])))
  expect_output(print(fit_summary))

  # hypothesis(): qCMDE/IWMDE are unavailable for component-specific scale
  # formulas, so the scale-formula quantities are tested with KDE.
  quantities <- hypothesis_quantities(fit)
  scale_rows <- quantities[quantities[["component"]] == "scale", , drop = FALSE]
  expect_setequal(
    unique(scale_rows[["parameter"]]),
    c("log_tau_study_intercept", "log_tau_study_x")
  )
  expect_true(all(scale_rows[["point_test_methods"]] == "KDE"))
  expect_error(
    hypothesis(fit, "log_tau_study_x = 0", component = "scale"),
    class = "RoBMA_density_method_scale_components"
  )
  point <- suppressWarnings(hypothesis(
    fit, "log_tau_study_x = 0", component = "scale", density_method = "KDE"
  ))
  expect_s3_class(point, "BayesTools_hypothesis_BF")
  expect_true(is.finite(attr(point, "raw_BF")))
  direction <- hypothesis(
    fit, "log_tau_study_x > 0", component = "scale", density_method = "KDE"
  )
  expect_s3_class(direction, "BayesTools_hypothesis_BF")
  expect_true(is.finite(attr(direction, "raw_BF")))
  expect_output(print(point))
})
