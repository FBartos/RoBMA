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

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

expect_influence_object <- function(x, n, inf_cols, min_dfbs_cols = 1) {

  expect_s3_class(x, "infl.brma")
  expect_true(all(inf_cols %in% names(x[["inf"]])))
  expect_equal(nrow(x[["inf"]]), n)
  expect_true(all(is.finite(as.matrix(x[["inf"]][, inf_cols, drop = FALSE]))))
  expect_s3_class(x[["dfbs"]], "data.frame")
  expect_equal(nrow(x[["dfbs"]]), n)
  expect_gte(ncol(x[["dfbs"]]), min_dfbs_cols)
  expect_true(all(is.finite(as.matrix(x[["dfbs"]]))))
}


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

test_that("Influence stats for meta-regression with interaction match metafor", {

  name <- "bcg_meta-regression4"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # Compute Influence
  inf_metafor <- influence(fit_metafor)
  inf_brma    <- suppressWarnings(influence(fit_brma))

  # Skip outliers not well approximated by loo
  skip_case <- c(4, 5, 6, 8, 10)

  # Check the individual components
  expect_equal(inf_metafor$inf$rstudent[-skip_case], inf_brma$inf$rstudent[-skip_case],  tolerance = 0.10, info = "rstudent matches")
  expect_equal(inf_metafor$inf$dffits[-skip_case],   inf_brma$inf$dffits[-skip_case],    tolerance = 0.10, info = "dffits matches")
  # expect_equal(inf_metafor$inf$cook.d[-skip_case],   inf_brma$inf$cook.d[-skip_case],    tolerance = 0.10, info = "cook.d matches") # TODO: check
  # skip covariance.ratio check - the metafors approach does not seen to account for uncertrainty in tau
  # expect_equal(inf_metafor$inf$cov.r[-skip_case],    inf_brma$inf$cov.r[-skip_case],     tolerance = 0.10, info = "cov.r matches")
  # expect_equal(inf_metafor$inf$tau2.del[-skip_case], inf_brma$inf$tau.del[-skip_case]^2, tolerance = 0.10, info = "tau2.del matches")
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
  expect_equal(as.vector(cooks.distance(fit_metafor)), suppressWarnings(cooks.distance(fit_brma)), tolerance = 0.05)
  # no dffits
  # no cov.r
})

test_that("Influence stats for 3-level meta-regression match metafor", {

  name <- "konstantopoulos2011_3lvl2"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # There is no influence for multivariate models (check the separate functions)
  expect_equal(as.vector(cooks.distance(fit_metafor)), suppressWarnings(cooks.distance(fit_brma)), tolerance = 0.05)
  # no dffits
  # no cov.r
})

# ============================================================================ #
# Test: Selection Model Influence
# ============================================================================ #

test_that("Influence stats for selection model work", {

  name <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # There is no influence for selection models nor seperate functions
  inf_brma <- suppressWarnings(influence(fit_brma))
  expect_true(!is.null(inf_brma))
})

test_that("DFFITS errors for unsupported model families", {

  skip_if_not("bcg_glmm" %in% names(fits), "GLMM cached fit not available.")
  skip_if_not("bcg_BMA.glmm" %in% names(fits), "BMA.glmm cached fit not available.")
  skip_if_not("dat.lehmann2018-3PSM" %in% names(fits), "Selection cached fit not available.")
  skip_if_not("dat.lehmann2018_RoBMA" %in% names(fits), "RoBMA cached fit not available.")

  expect_error(
    dffits(fits[["bcg_glmm"]]),
    "only available for normal outcome models"
  )
  expect_error(
    dffits(fits[["bcg_BMA.glmm"]]),
    "only available for normal outcome models"
  )
  expect_error(
    dffits(fits[["dat.lehmann2018-3PSM"]]),
    "not available for selection models"
  )
  expect_error(
    dffits(fits[["dat.lehmann2018_RoBMA"]]),
    "not available for selection models"
  )
})

# ============================================================================ #
# Test: Model-Averaging Influence Diagnostics
# ============================================================================ #

test_that("Influence stats for BMA.norm model-averaging fit are internally consistent", {

  name <- "dat.lehmann2018_BMA.norm"
  fit_brma <- fits[[name]]
  n        <- nrow(fit_brma[["data"]][["outcome"]])
  inf_brma <- suppressWarnings(influence(fit_brma))

  expect_influence_object(
    inf_brma,
    n,
    inf_cols = c("rstudent", "dffits", "cook.d", "cov.r", "tau.del", "hat")
  )
})

test_that("Influence stats for BMA.glmm model-averaging fit are internally consistent", {

  name <- "bcg_BMA.glmm"
  fit_brma <- fits[[name]]
  n        <- nrow(fit_brma[["data"]][["outcome"]])
  inf_brma <- suppressWarnings(influence(fit_brma))

  expect_influence_object(
    inf_brma,
    n,
    inf_cols = c("rstudent", "cov.r", "tau.del")
  )
  expect_false(any(c("dffits", "cook.d", "hat") %in% names(inf_brma[["inf"]])))
  expect_error(cooks.distance(fit_brma), "normal outcome models")
})

test_that("Influence stats for RoBMA model-averaging fit are internally consistent", {

  name <- "dat.lehmann2018_RoBMA"
  fit_brma <- fits[[name]]
  n        <- nrow(fit_brma[["data"]][["outcome"]])
  inf_brma <- suppressWarnings(influence(fit_brma))

  expect_influence_object(
    inf_brma,
    n,
    inf_cols = c("rstudent", "cov.r", "tau.del")
  )
  expect_false(any(c("dffits", "cook.d", "hat") %in% names(inf_brma[["inf"]])))
  expect_error(cooks.distance(fit_brma), "selection models")
})
