context("Residuals")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "residuals")

# list & load all fits
skip_if_no_fits()
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


# ============================================================================ #
# Test: Simple Meta-Analysis Residuals
# ============================================================================ #

test_that("Residuals for simple meta-analysis match metafor", {

  skip_if_not_installed("metafor")

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Residuals: compare brma vs metafor
  # --------------------------------------------------

  # metafor: residuals returns the difference between yi and fitted values
  metafor_resid <- residuals(fit_metafor)

  # brma: using residuals method
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  # comparison (allow MCMC tolerance)
  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma residuals should match metafor residuals")

  # --------------------------------------------------
  # Verify residuals sum to approximately zero (for simple model)
  # --------------------------------------------------

  expect_equal(mean(brma_resid_mean), 0, tolerance = 0.1,
               info = "residuals should sum approximately to zero for simple model")

  # --------------------------------------------------
  # Verify residuals match yi - fitted
  # --------------------------------------------------

  yi <- fit_brma$data$outcome$yi
  pooled_mu <- pooled_effect(fit_brma)$summary["mu", "Mean"]
  expected_resid <- yi - pooled_mu

  expect_equal(brma_resid_mean, expected_resid, tolerance = 0.05,
               info = "residuals should equal yi - pooled effect")
})


# ============================================================================ #
# Test: Meta-Regression Residuals
# ============================================================================ #

test_that("Residuals for meta-regression match metafor", {

  skip_if_not_installed("metafor")

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Residuals: compare brma vs metafor
  # --------------------------------------------------

  # metafor: residuals for meta-regression
  metafor_resid <- residuals(fit_metafor)

  # brma: using residuals method
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  # comparison (allow MCMC tolerance)
  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.10,
               info = "brma residuals should match metafor residuals for meta-regression")
})


# ============================================================================ #
# Test: Location-Scale Model Residuals
# ============================================================================ #

test_that("Residuals for location-scale model match metafor", {

  skip_if_not_installed("metafor")

  name        <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Residuals: compare brma vs metafor
  # --------------------------------------------------

  # metafor: residuals for location-scale model
  metafor_resid <- residuals(fit_metafor)

  # brma: using residuals method
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  # comparison (allow MCMC tolerance)
  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma residuals should match metafor residuals for location-scale model")
})


# ============================================================================ #
# Test: 3-Level Model Residuals
# ============================================================================ #

test_that("Residuals for 3-level model match expected values", {

  skip_if_not_installed("metafor")

  name        <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Residuals computation
  # --------------------------------------------------

  # metafor: residuals for multilevel model
  metafor_resid <- residuals(fit_metafor)

  # brma: using residuals method
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  # comparison (allow MCMC tolerance)
  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma residuals should match metafor residuals for 3-level model")

})


# ============================================================================ #
# Test: Selection Model Residuals
# ============================================================================ #

test_that("Residuals for selection model are computed correctly", {

  skip_if_not_installed("metafor")

  name        <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Residuals computation
  # --------------------------------------------------

  # metafor: residuals for multilevel model
  metafor_resid <- residuals(fit_metafor)

  # brma: using residuals method
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  # comparison (allow MCMC tolerance)
  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma residuals should match metafor residuals for selection model")
})


# ============================================================================ #
# Test: PET Model Residuals
# ============================================================================ #

test_that("Residuals for PET model are computed correctly", {

  skip_if_not_installed("metafor")

  name        <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Residuals computation
  # --------------------------------------------------

  # metafor: residuals for multilevel model
  metafor_resid <- residuals(fit_metafor)

  # brma: using residuals method
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  # comparison (allow MCMC tolerance)
  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma residuals should match metafor residuals for selection model")
})


# ============================================================================ #
# Test: GLMM Model Residuals
# ============================================================================ #

test_that("Residuals for GLMM model are computed correctly", {

  skip_if_not_installed("metafor")

  name        <- "nielweise2008_glmm"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Residuals computation
  # --------------------------------------------------

  # metafor: residuals for multilevel model
  metafor_resid <- residuals(fit_metafor)

  # brma: using residuals method
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  # comparison (allow MCMC tolerance)
  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma residuals should match metafor residuals for GLMM model")
})


# ============================================================================ #
# Test: Residuals Function Interface
# ============================================================================ #

test_that("Residuals function has correct interface", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_samples = TRUE returns matrix
  # --------------------------------------------------

  samples_resid <- residuals(fit_brma, as_samples = TRUE)
  expect_true(is.matrix(samples_resid),
              info = "residuals with as_samples = TRUE should return matrix")

  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(ncol(samples_resid), n_studies,
               info = "residuals should return one column per study")

  # --------------------------------------------------
  # Test probs argument
  # --------------------------------------------------

  resid_90ci <- residuals(fit_brma, probs = c(0.05, 0.95))
  ci_names <- colnames(resid_90ci$summary)
  # BayesTools may format as "5%" or "0.05" - check for either format
  has_lower <- any(grepl("5%|0\\.05", ci_names))
  has_upper <- any(grepl("95%|0\\.95", ci_names))
  expect_true(has_lower && has_upper,
              info = "probs argument should control CI quantiles")

  # --------------------------------------------------
  # Test return class
  # --------------------------------------------------

  resid_result <- residuals(fit_brma)
  expect_s3_class(resid_result, "brma.residuals")
  expect_true("summary" %in% names(resid_result),
              info = "result should contain summary element")
  expect_true("data" %in% names(resid_result),
              info = "result should contain data element")

})
