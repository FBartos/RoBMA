context("Influence Diagnostics")
# these test implicitly test
# - cooks.distance
# - dffits
# - covratio
# - dfbetas
# - rstudent
# - hatvalues

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list & load all fits
skip_if_no_fits()
skip_if_not_installed("metafor")

fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


# ============================================================================ #
# Test: Simple Meta-Analysis Influence
# ============================================================================ #

test_that("Influence stats for simple meta-analysis match metafor", {

  name <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # Compute Influence
  inf_metafor <- influence(fit_metafor)
  inf_brma    <- influence(fit_brma)

  # Check the individual components
  expect_equal(inf_metafor$inf$rstudent, inf_brma$inf$rstudent,  tolerance = 0.05, info = "rstudent matches")
  expect_equal(inf_metafor$inf$dffits,   inf_brma$inf$dffits,    tolerance = 0.05, info = "dffits matches")
  expect_equal(inf_metafor$inf$cook.d,   inf_brma$inf$cook.d,    tolerance = 0.05, info = "cook.d matches")
  expect_equal(inf_metafor$inf$cov.r,    inf_brma$inf$cov.r,     tolerance = 0.05, info = "cov.r matches")
  expect_equal(inf_metafor$inf$tau2.del, inf_brma$inf$tau.del^2, tolerance = 0.05, info = "tau2.del matches")
  expect_equal(inf_metafor$inf$hat,      inf_brma$inf$hat,       tolerance = 0.05, info = "hat matches")
  expect_equal(inf_metafor$dfbs$intrcpt, inf_brma$dfbs$mu,       tolerance = 0.10, info = "dfbetas matches")

  # Check standalone functions
  expect_equal(inf_metafor$inf$dffits, dffits(fit_brma),         tolerance = 0.05)
  expect_equal(inf_metafor$inf$cook.d, cooks.distance(fit_brma), tolerance = 0.05)
  expect_equal(inf_metafor$inf$cov.r,  covratio(fit_brma),       tolerance = 0.05)
})


# ============================================================================ #
# Test: Meta-Regression Influence
# ============================================================================ #

test_that("Influence stats for meta-regression match metafor", {

  name <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # Compute Influence
  inf_metafor <- influence(fit_metafor)
  inf_brma    <- suppressWarnings(influence(fit_brma))

  # Skip outliers not well approximated by loo
  skip_case <- c(4, 6, 13)

  # Check the individual components
  expect_equal(inf_metafor$inf$rstudent[-skip_case], inf_brma$inf$rstudent[-skip_case],  tolerance = 0.15, info = "rstudent matches")
  expect_equal(inf_metafor$inf$dffits[-skip_case],   inf_brma$inf$dffits[-skip_case],    tolerance = 0.10, info = "dffits matches")
  expect_equal(inf_metafor$inf$cook.d[-skip_case],   inf_brma$inf$cook.d[-skip_case],    tolerance = 0.10, info = "cook.d matches")
  # skip covariance.ratio check - the metafors approach does not seen to account for uncertrainty in tau
  # expect_equal(inf_metafor$inf$cov.r[-skip_case],    inf_brma$inf$cov.r[-skip_case],     tolerance = 0.10, info = "cov.r matches")
  expect_equal(inf_metafor$inf$tau2.del[-skip_case], inf_brma$inf$tau.del[-skip_case]^2, tolerance = 0.10, info = "tau2.del matches")
  expect_equal(inf_metafor$inf$hat[-skip_case],      inf_brma$inf$hat[-skip_case],       tolerance = 0.10, info = "hat matches")
})


# ============================================================================ #
# Test: 3-Level Model Influence
# ============================================================================ #

test_that("Influence stats for 3-level model match metafor", {

  name <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # There is no influence for multivariate models (check the separate functions)
  expect_equal(as.vector(cooks.distance(fit_metafor)), cooks.distance(fit_brma), tolerance = 0.05)
  # no dffits
  # no cov.r
})

# ============================================================================ #
# Test: Selection Model Influence
# ============================================================================ #

test_that("Influence stats for 3-level model match metafor", {

  name <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # There is no influence for selection models nor seperate functions
  inf_brma <- influence(fit_brma)
  expect_true(!is.null(inf_brma))
})

