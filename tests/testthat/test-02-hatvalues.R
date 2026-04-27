context("Hatvalues")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "hatvalues")

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

expect_hatvalues_vector <- function(x, n) {

  expect_type(x, "double")
  expect_equal(length(x), n)
  expect_true(all(is.finite(x)))
  expect_true(all(x >= 0 & x <= 1 + sqrt(.Machine$double.eps)))
}


# ============================================================================ #
# Test: Simple Meta-Analysis Hatvalues
# ============================================================================ #

test_that("Hatvalues for simple meta-analysis match metafor", {

  name <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_hat <- hatvalues(fit_metafor)
  brma_hat    <- hatvalues(fit_brma)

  expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor")
})


# ============================================================================ #
# Test: Meta-Regression Hatvalues
# ============================================================================ #

test_that("Hatvalues for meta-regression match metafor", {

  name <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_hat <- hatvalues(fit_metafor)
  brma_hat    <- hatvalues(fit_brma)

  expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for meta-regression")
})

test_that("Hatvalues for meta-regression with interaction match metafor", {

  name <- "bcg_meta-regression4"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_hat <- hatvalues(fit_metafor)
  brma_hat    <- hatvalues(fit_brma)

  expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for meta-regression")
})

# ============================================================================ #
# Test: 3-Level Model Hatvalues
# ============================================================================ #

test_that("Hatvalues for 3-level model match metafor", {

  name <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_hat <- hatvalues(fit_metafor)
  brma_hat    <- hatvalues(fit_brma)

  expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for 3-level model")
})

test_that("Hatvalues for 3-level meta-regression match metafor", {

  name <- "konstantopoulos2011_3lvl2"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_hat <- hatvalues(fit_metafor)
  brma_hat    <- hatvalues(fit_brma)

  expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for 3-level model")
})

# ============================================================================ #
# Test: PET Models Hatvalues
# ============================================================================ #

test_that("Hatvalues for PET models work", {

  name <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_hat <- hatvalues(fit_metafor)
  brma_hat    <- hatvalues(fit_brma)

  expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for PET model")
})

test_that("Hatvalues for PET meta-regression work", {

  name <- "dat.lehmann2018-PETreg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_hat <- hatvalues(fit_metafor)
  brma_hat    <- hatvalues(fit_brma)

  expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for PET model")
})

test_that("Hatvalues for PET models work  (negative)", {

  name <- "dat.lehmann2018-PET_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_hat <- hatvalues(fit_metafor)
  brma_hat    <- hatvalues(fit_brma)

  expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for PET model")
})

# ============================================================================ #
# Test: glmm/bselmod Hatvalues are not available
# ============================================================================ #

test_that("Hatvalues should error for other models", {

  name <- "nielweise2008_glmm"
  fit_brma <- fits[[name]]

  expect_error(hatvalues(fit_brma), "available")

  name <- "dat.lehmann2018-3PSM"
  fit_brma <- fits[[name]]

  expect_error(hatvalues(fit_brma), "available")
})

# ============================================================================ #
# Test: Model-Averaging Hatvalues
# ============================================================================ #

test_that("Hatvalues for BMA.norm model-averaging fits are internally consistent", {

  fit_simple <- fits[["dat.lehmann2018_BMA.norm"]]
  n_simple   <- nrow(fit_simple[["data"]][["outcome"]])
  hat_simple <- hatvalues(fit_simple)
  expect_hatvalues_vector(hat_simple, n_simple)

  fit_mods <- fits[["dat.lehmann2018_BMA.norm_mods"]]
  n_mods   <- nrow(fit_mods[["data"]][["outcome"]])
  hat_mods <- hatvalues(fit_mods)
  expect_hatvalues_vector(hat_mods, n_mods)
})

test_that("Hatvalues error for BMA.glmm and RoBMA model-averaging fits", {

  names <- c("bcg_BMA.glmm", "dat.lehmann2018_RoBMA")
  skip_if_missing_fits(names)

  expect_error(
    hatvalues(fits[["bcg_BMA.glmm"]]),
    "normal outcome models"
  )
  expect_error(
    hatvalues(fits[["dat.lehmann2018_RoBMA"]]),
    "selection models"
  )
})
