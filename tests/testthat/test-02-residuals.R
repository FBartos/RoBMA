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
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson <- residuals(fit_brma, type = "pearson")
  brma_pearson_mean <- brma_pearson$summary[, "Mean"]

  expect_equal(brma_pearson_mean, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard <- residuals(fit_brma, type = "rstandard")
  brma_rstandard_mean <- brma_rstandard$summary[, "Mean"]

  expect_equal(brma_rstandard_mean, as.vector(metafor_rstandard$z), tolerance = 0.05,
               info = "brma rstandard residuals should match metafor")
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
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.10,
               info = "brma outcome residuals should match metafor for meta-regression")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson <- residuals(fit_brma, type = "pearson")
  brma_pearson_mean <- brma_pearson$summary[, "Mean"]

  expect_equal(brma_pearson_mean, as.vector(metafor_pearson), tolerance = 0.10,
               info = "brma pearson residuals should match metafor for meta-regression")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard <- residuals(fit_brma, type = "rstandard")
  brma_rstandard_mean <- brma_rstandard$summary[, "Mean"]

  expect_equal(brma_rstandard_mean, as.vector(metafor_rstandard$z), tolerance = 0.10,
               info = "brma rstandard residuals should match metafor for meta-regression")
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
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for location-scale model")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson <- residuals(fit_brma, type = "pearson")
  brma_pearson_mean <- brma_pearson$summary[, "Mean"]

  expect_equal(brma_pearson_mean, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for location-scale model")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard <- residuals(fit_brma, type = "rstandard")
  brma_rstandard_mean <- brma_rstandard$summary[, "Mean"]

  expect_equal(brma_rstandard_mean, as.vector(metafor_rstandard$z), tolerance = 0.05,
               info = "brma rstandard residuals should match metafor for location-scale model")
})


# ============================================================================ #
# Test: 3-Level Model Residuals
# ============================================================================ #

test_that("Residuals for 3-level model match metafor", {

  skip_if_not_installed("metafor")

  name        <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for 3-level model")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson <- residuals(fit_brma, type = "pearson")
  brma_pearson_mean <- brma_pearson$summary[, "Mean"]

  expect_equal(brma_pearson_mean, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for 3-level model")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard <- residuals(fit_brma, type = "rstandard")
  brma_rstandard_mean <- brma_rstandard$summary[, "Mean"]

  expect_equal(brma_rstandard_mean, as.vector(metafor_rstandard$z), tolerance = 0.05,
               info = "brma rstandard residuals should match metafor for 3-level model")
})


# ============================================================================ #
# Test: Selection Model Residuals (Positive Effects)
# ============================================================================ #

test_that("Residuals for selection model (positive) match metafor", {

  skip_if_not_installed("metafor")

  name        <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for selection model (positive)")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson <- residuals(fit_brma, type = "pearson")
  brma_pearson_mean <- brma_pearson$summary[, "Mean"]

  expect_equal(brma_pearson_mean, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for selection model (positive)")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------
# TODO: not implemented in metafor
#   metafor_rstandard <- rstandard(fit_metafor)
#   brma_rstandard <- residuals(fit_brma, type = "rstandard")
#   brma_rstandard_mean <- brma_rstandard$summary[, "Mean"]
#
#   expect_equal(brma_rstandard_mean, as.vector(metafor_rstandard$z), tolerance = 0.05,
#                info = "brma rstandard residuals should match metafor for selection model (positive)")
})


# ============================================================================ #
# Test: Selection Model Residuals (Negative Effects)
# ============================================================================ #

test_that("Residuals for selection model (negative) match metafor", {

  skip_if_not_installed("metafor")

  name        <- "dat.lehmann2018-3PSM_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for selection model (negative)")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson <- residuals(fit_brma, type = "pearson")
  brma_pearson_mean <- brma_pearson$summary[, "Mean"]

  expect_equal(brma_pearson_mean, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for selection model (negative)")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------
  # TODO: not implemented in metafor
  #   metafor_rstandard <- rstandard(fit_metafor)
  #   brma_rstandard <- residuals(fit_brma, type = "rstandard")
  #   brma_rstandard_mean <- brma_rstandard$summary[, "Mean"]
  #
  #   expect_equal(brma_rstandard_mean, as.vector(metafor_rstandard$z), tolerance = 0.05,
  #                info = "brma rstandard residuals should match metafor for selection model (negative)")
})


# ============================================================================ #
# Test: PET Model Residuals (Positive Effects)
# ============================================================================ #

test_that("Residuals for PET model (positive) match metafor", {

  skip_if_not_installed("metafor")

  name        <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for PET model (positive)")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson <- residuals(fit_brma, type = "pearson")
  brma_pearson_mean <- brma_pearson$summary[, "Mean"]

  expect_equal(brma_pearson_mean, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for PET model (positive)")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard <- residuals(fit_brma, type = "rstandard")
  brma_rstandard_mean <- brma_rstandard$summary[, "Mean"]

  expect_equal(brma_rstandard_mean, as.vector(metafor_rstandard$z), tolerance = 0.05,
               info = "brma rstandard residuals should match metafor for PET model (positive)")
})


# ============================================================================ #
# Test: PET Model Residuals (Negative Effects)
# ============================================================================ #

test_that("Residuals for PET model (negative) match metafor", {

  skip_if_not_installed("metafor")

  name        <- "dat.lehmann2018-PET_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for PET model (negative)")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson <- residuals(fit_brma, type = "pearson")
  brma_pearson_mean <- brma_pearson$summary[, "Mean"]

  expect_equal(brma_pearson_mean, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for PET model (negative)")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard <- residuals(fit_brma, type = "rstandard")
  brma_rstandard_mean <- brma_rstandard$summary[, "Mean"]

  expect_equal(brma_rstandard_mean, as.vector(metafor_rstandard$z), tolerance = 0.05,
               info = "brma rstandard residuals should match metafor for PET model (negative)")
})


# ============================================================================ #
# Test: GLMM Model Residuals
# ============================================================================ #

test_that("Residuals for GLMM model are computed correctly", {

  skip_if_not_installed("metafor")

  name        <- "nielweise2008_glmm"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  n_studies <- nrow(fit_brma$data$outcome)

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid <- residuals(fit_brma)
  brma_resid_mean <- brma_resid$summary[, "Mean"]

  expect_equal(brma_resid_mean, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for GLMM model")

  # --------------------------------------------------
  # Pearson residuals: check finite and correct length
  # (metafor GLMM may not support pearson residuals directly)
  # --------------------------------------------------

  # TODO: these do not exist for metafor
  brma_pearson <- residuals(fit_brma, type = "pearson")
  brma_pearson_mean <- brma_pearson$summary[, "Mean"]

  expect_true(all(is.finite(brma_pearson_mean)),
              info = "pearson residuals should be finite for GLMM model")
  expect_equal(length(brma_pearson_mean), n_studies,
               info = "pearson residuals should have one value per study")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: check finite and correct length
  # --------------------------------------------------

  # TODO: these do not exist for metafor
  brma_rstandard <- residuals(fit_brma, type = "rstandard")
  brma_rstandard_mean <- brma_rstandard$summary[, "Mean"]

  expect_true(all(is.finite(brma_rstandard_mean)),
              info = "rstandard residuals should be finite for GLMM model")
  expect_equal(length(brma_rstandard_mean), n_studies,
               info = "rstandard residuals should have one value per study")
})


# ============================================================================ #
# Test: Residuals Function Interface
# ============================================================================ #

test_that("Residuals function has correct interface", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  n_studies <- nrow(fit_brma$data$outcome)

  # --------------------------------------------------
  # Test type argument accepts all valid values
  # --------------------------------------------------

  expect_no_error(residuals(fit_brma))
  expect_no_error(residuals(fit_brma, type = "pearson"))
  expect_no_error(residuals(fit_brma, type = "rstandard"))

  # --------------------------------------------------
  # Test as_samples = TRUE returns matrix for all types
  # --------------------------------------------------

  samples_outcome   <- residuals(fit_brma,   as_samples = TRUE)
  samples_pearson   <- residuals(fit_brma, type = "pearson",   as_samples = TRUE)
  samples_rstandard <- residuals(fit_brma, type = "rstandard", as_samples = TRUE)

  expect_true(is.matrix(samples_outcome),
              info = "outcome residuals with as_samples = TRUE should return matrix")
  expect_true(is.matrix(samples_pearson),
              info = "pearson residuals with as_samples = TRUE should return matrix")
  expect_true(is.matrix(samples_rstandard),
              info = "rstandard residuals with as_samples = TRUE should return matrix")

  # all should have same dimensions
  expect_equal(dim(samples_outcome), dim(samples_pearson))
  expect_equal(dim(samples_outcome), dim(samples_rstandard))
  expect_equal(ncol(samples_outcome), n_studies,
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

  # --------------------------------------------------
  # Test that all residual types return finite values
  # --------------------------------------------------

  expect_true(all(is.finite(samples_outcome)),
              info = "outcome residuals should all be finite")
  expect_true(all(is.finite(samples_pearson)),
              info = "pearson residuals should all be finite")
  expect_true(all(is.finite(samples_rstandard)),
              info = "rstandard residuals should all be finite")

  # --------------------------------------------------
  # Test rstandard() wrapper function
  # --------------------------------------------------

  rstandard_result <- rstandard(fit_brma)
  residuals_rstandard <- residuals(fit_brma, type = "rstandard")

  expect_equal(rstandard_result$summary, residuals_rstandard$summary,
               info = "rstandard() should return same result as residuals(type='rstandard')")})
