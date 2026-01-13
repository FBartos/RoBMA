context("Residuals")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list & load all fits
skip_if_no_fits()
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


# ============================================================================ #
# Test: Residuals Match Metafor
# ============================================================================ #

test_that("Residuals for simple meta-analysis match metafor", {

  skip_if_not_installed("metafor")

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Residuals: compare brma vs metafor
  # --------------------------------------------------

  # metafor: residuals
  metafor_resid <- stats::residuals(fit_metafor)

  # brma: residuals
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  # comparison
  expect_equal(brma_resid_mean, metafor_resid, tolerance = 0.05,
               info = "brma residuals should match metafor residuals")

})


test_that("Residuals for meta-regression match metafor", {

  skip_if_not_installed("metafor")

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Residuals: compare brma vs metafor
  # --------------------------------------------------

  # metafor: residuals
  metafor_resid <- stats::residuals(fit_metafor)

  # brma: residuals
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  # comparison (higher tolerance for regression models)
  expect_equal(brma_resid_mean, metafor_resid, tolerance = 0.1,
               info = "brma residuals should match metafor residuals for meta-regression")

})


test_that("Residuals with bias_adjusted argument work correctly", {

  skip_if_not_installed("metafor")

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test bias_adjusted parameter
  # --------------------------------------------------

  # bias-adjusted residuals (default)
  resid_bias_adj <- residuals(fit_brma, bias_adjusted = TRUE)

  # non-bias-adjusted residuals
  resid_no_bias <- residuals(fit_brma, bias_adjusted = FALSE)

  # both should return brma.predict objects
  expect_s3_class(resid_bias_adj, "brma.predict")
  expect_s3_class(resid_no_bias, "brma.predict")

  # both should have summary tables
  expect_true("summary" %in% names(resid_bias_adj))
  expect_true("summary" %in% names(resid_no_bias))

  # both should have same number of residuals as observations
  K <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(resid_bias_adj$summary), K)
  expect_equal(nrow(resid_no_bias$summary), K)

})


test_that("Residuals as_samples returns matrix", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_samples parameter
  # --------------------------------------------------

  # get residual samples
  resid_samples <- residuals(fit_brma, as_samples = TRUE)

  # should be a matrix
  expect_true(is.matrix(resid_samples))

  # should have correct dimensions (S samples x K observations)
  K <- nrow(fit_brma$data$outcome)
  expect_equal(ncol(resid_samples), K)
  expect_true(nrow(resid_samples) > 1000)  # should have many posterior samples

  # column names should be residual[1], residual[2], etc.
  expected_names <- paste0("residual[", seq_len(K), "]")
  expect_equal(colnames(resid_samples), expected_names)

})
