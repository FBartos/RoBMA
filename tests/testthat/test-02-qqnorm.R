context("QQ Normal plot")

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
# Test: Simple Meta-Analysis QQ Plot
# ============================================================================ #

test_that("QQ plot for simple meta-analysis matches metafor structure", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Visual comparison: side-by-side (rstandard for comparability)
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_simple_comparison_rstandard", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    qqnorm(fit_metafor, main = "metafor", type = "rstandard", ylim = c(-3, 3))
    qqnorm(fit_brma, plot_type = "base", main = "brma", type = "rstandard", ylim = c(-3, 3))
  })

  # --------------------------------------------------
  # Default (rstudent) - brma only
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_simple_rstudent_base", function() {
    qqnorm(fit_brma, plot_type = "base")
  })

  vdiffr::expect_doppelganger(
    "qqnorm_simple_rstudent_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: Meta-Regression QQ Plot
# ============================================================================ #

test_that("QQ plot for meta-regression works correctly", {

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Visual comparison: side-by-side (rstandard for comparability)
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_regression_comparison_rstandard", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    qqnorm(fit_metafor, main = "metafor", type = "rstandard", ylim = c(-3, 3))
    qqnorm(fit_brma, plot_type = "base", main = "brma", type = "rstandard", ylim = c(-3, 3))
  })

  vdiffr::expect_doppelganger(
    "qqnorm_regression_rstudent_ggplot",
    suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  )
})


# ============================================================================ #
# Test: Location-Scale Model QQ Plot
# ============================================================================ #

test_that("QQ plot for location-scale model works correctly", {

  name     <- "bangertdrowns2004_location-scale"
  fit_brma <- fits[[name]]
  set.seed(1)

  vdiffr::expect_doppelganger("qqnorm_scale_rstandard_base", function() {
    qqnorm(fit_brma, plot_type = "base", type = "rstandard")
  })

  vdiffr::expect_doppelganger(
    "qqnorm_scale_rstudent_ggplot",
    suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  )
})


# ============================================================================ #
# Test: 3-Level Model QQ Plot
# ============================================================================ #

test_that("QQ plot for 3-level model works correctly", {

  name     <- "konstantopoulos2011_3lvl"
  fit_brma <- fits[[name]]
  set.seed(1)

  vdiffr::expect_doppelganger("qqnorm_3lvl_base", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })

  vdiffr::expect_doppelganger(
    "qqnorm_3lvl_ggplot",
    suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  )
})


# ============================================================================ #
# Test: GLMM Model QQ Plot
# ============================================================================ #

test_that("QQ plot for GLMM model works correctly", {

  name     <- "nielweise2008_glmm"
  fit_brma <- fits[[name]]
  set.seed(1)

  # rstudent only (rstandard not available for GLMM)
  vdiffr::expect_doppelganger("qqnorm_glmm_base", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })

  vdiffr::expect_doppelganger(
    "qqnorm_glmm_ggplot",
    suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  )

  # rstandard should error for GLMM models
  expect_error(
    qqnorm(fit_brma, type = "rstandard"),
    info = "rstandard should error for GLMM models"
  )
})


# ============================================================================ #
# Test: PET Model QQ Plot
# ============================================================================ #

test_that("QQ plot for PET works correctly", {

  name        <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Visual comparison: side-by-side (rstandard for comparability)
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_PET_base", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    qqnorm(fit_metafor, main = "metafor", type = "rstudent", ylim = c(-3, 3))
    suppressWarnings(qqnorm(fit_brma, plot_type = "base", main = "brma", type = "rstudent", ylim = c(-3, 3)))
  })
})

# ============================================================================ #
# Test: Selection Model QQ Plot
# ============================================================================ #

test_that("QQ plot for selection model works correctly", {

  name     <- "dat.lehmann2018-3PSM"
  fit_brma <- fits[[name]]
  set.seed(1)

  vdiffr::expect_doppelganger("qqnorm_selmodel_base", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })

  # skip the ggplot version for to save test time
  # vdiffr::expect_doppelganger(
  #   "qqnorm_selmodel_ggplot",
  #   suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  # )

  # rstandard should error for selection models
  expect_error(
    qqnorm(fit_brma, type = "rstandard"),
    info = "rstandard should error for selection models"
  )
})


# ============================================================================ #
# Test: QQ Plot Interface
# ============================================================================ #

test_that("QQ plot has correct interface", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Test as_data = TRUE returns list with expected components
  # --------------------------------------------------

  qq_data <- qqnorm(fit_brma, as_data = TRUE, type = "rstandard")

  expect_true(is.list(qq_data),
    info = "as_data = TRUE should return a list"
  )

  expected_components <- c(
    "x", "y", "points", "refline", "envelope",
    "xlim", "ylim", "xlab", "ylab"
  )
  expect_true(all(expected_components %in% names(qq_data)),
    info = "QQ data should contain all expected components"
  )

  # Check points data.frame structure
  expect_true(is.data.frame(qq_data$points),
    info = "points should be a data.frame"
  )
  expect_true(all(c("x", "y") %in% names(qq_data$points)),
    info = "points should have x and y columns"
  )

  # Check number of points matches number of observations
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(qq_data$points), n_studies,
    info = "number of points should match number of observations"
  )

  # Check x values are sorted (theoretical quantiles)
  expect_equal(qq_data$x, sort(qq_data$x),
    info = "theoretical quantiles should be sorted"
  )

  # Check y values are sorted (sorted residuals)
  expect_equal(qq_data$y, sort(qq_data$y),
    info = "sample quantiles should be sorted"
  )

  # --------------------------------------------------
  # Test envelope = FALSE suppresses envelope
  # --------------------------------------------------

  qq_data_no_env <- qqnorm(fit_brma, as_data = TRUE, envelope = FALSE, type = "rstandard")
  expect_null(qq_data_no_env$envelope,
    info = "envelope should be NULL when envelope = FALSE"
  )

  # --------------------------------------------------
  # Test envelope = TRUE includes envelope
  # --------------------------------------------------

  expect_true(is.data.frame(qq_data$envelope),
    info = "envelope should be a data.frame when envelope = TRUE"
  )
  expect_true(all(c("x", "lower", "upper") %in% names(qq_data$envelope)),
    info = "envelope should have x, lower, upper columns"
  )

  # --------------------------------------------------
  # Test error on invalid plot_type
  # --------------------------------------------------

  expect_error(qqnorm(fit_brma, plot_type = "invalid"),
    info = "should error on invalid plot_type"
  )

  # --------------------------------------------------
  # Test error on invalid type
  # --------------------------------------------------

  expect_error(qqnorm(fit_brma, type = "invalid"),
    info = "should error on invalid type"
  )
})


# ============================================================================ #
# Test: QQ Plot Customization
# ============================================================================ #

test_that("QQ plot customization works", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Test custom point aesthetics
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_custom_points_base", function() {
    qqnorm(fit_brma, plot_type = "base", pch = 21, col = "blue",
           bg = "lightblue", cex = 1.5, type = "rstandard")
  })

  vdiffr::expect_doppelganger(
    "qqnorm_custom_points_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", pch = 21, col = "blue",
           bg = "lightblue", size = 3, type = "rstandard")
  )

  # --------------------------------------------------
  # Test custom axis labels and title
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_custom_labels_base", function() {
    qqnorm(fit_brma, plot_type = "base", xlab = "Expected", ylab = "Observed",
           main = "QQ Plot", type = "rstandard")
  })

  vdiffr::expect_doppelganger(
    "qqnorm_custom_labels_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", xlab = "Expected", ylab = "Observed",
           main = "QQ Plot", type = "rstandard")
  )

  # --------------------------------------------------
  # Test without envelope
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_no_envelope_base", function() {
    qqnorm(fit_brma, plot_type = "base", envelope = FALSE, type = "rstandard")
  })

  vdiffr::expect_doppelganger(
    "qqnorm_no_envelope_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", envelope = FALSE, type = "rstandard")
  )

  # --------------------------------------------------
  # Test custom envelope level
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_custom_level_base", function() {
    qqnorm(fit_brma, plot_type = "base", level = 99, type = "rstandard")
  })

  vdiffr::expect_doppelganger(
    "qqnorm_custom_level_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", level = 99, type = "rstandard")
  )

  # --------------------------------------------------
  # Test Bonferroni correction
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_bonferroni_base", function() {
    qqnorm(fit_brma, plot_type = "base", bonferroni = TRUE, type = "rstandard")
  })

  vdiffr::expect_doppelganger(
    "qqnorm_bonferroni_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", bonferroni = TRUE, type = "rstandard")
  )

  # --------------------------------------------------
  # Test suppressing envelope shading
  # --------------------------------------------------

  vdiffr::expect_doppelganger("qqnorm_no_shade_base", function() {
    qqnorm(fit_brma, plot_type = "base", shade = NA, type = "rstandard")
  })

  vdiffr::expect_doppelganger(
    "qqnorm_no_shade_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", shade = NA, type = "rstandard")
  )
})
