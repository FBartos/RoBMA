context("Regression plot (bubble plot)")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list & load all fits
skip_if_no_fits()
skip_if_not_installed("metafor")
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()

# TODO:
# - missing border of the shade (no shade deletes everything)
# - the predicticted line and bands should extend all the way to the end of plotting range
# - remove these warning messages:
# - should use the quantile function (as the funnel plot) instead of doing random simulation via predict
# ```
# Warning messages:
#  1: In data.frame(x = c(at_pred, rev(at_pred)), y = c(pred_lower, rev(pred_upper)),  :
#  row names were found from a short variable and have been discarded
#  2: In data.frame(x = c(at_pred, rev(at_pred)), y = c(pi_lower, rev(pi_upper)),  :
#  row names were found from a short variable and have been discarded
# ```

# ============================================================================ #
# Test: Continuous Moderator Regression Plot
# ============================================================================ #

test_that("Regression plot for continuous moderator matches metafor structure", {

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_continuous_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::regplot(fit_metafor, mod = "year", main = "metafor", ylim = c(-2, 0.5))
    regplot(fit_brma, mod = "year", plot_type = "base", main = "brma", ylim = c(-2, 0.5))
  })

  vdiffr::expect_doppelganger(
    "regplot_continuous_brma_ggplot",
    regplot(fit_brma, mod = "year", plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: Continuous Moderator with Prediction Intervals
# ============================================================================ #

test_that("Regression plot with prediction intervals works", {

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Visual tests with PI bands
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_continuous_pi_base", function() {
    regplot(fit_brma, mod = "year", plot_type = "base", pi = TRUE, main = "with PI")
  })

  vdiffr::expect_doppelganger(
    "regplot_continuous_pi_ggplot",
    regplot(fit_brma, mod = "year", plot_type = "ggplot", pi = TRUE)
  )
})


# ============================================================================ #
# Test: Categorical Moderator Regression Plot
# ============================================================================ #

test_that("Regression plot for categorical moderator works", {

  name        <- "bcg_meta-regression-categorical"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # TODO: implement
  # skip if this fit doesn't exist
  skip("Categorical meta-regression fit not added yet")

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_categorical_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::regplot(fit_metafor, mod = "year", main = "metafor")
    regplot(fit_brma, mod = "year", plot_type = "base", main = "brma")
  })

  vdiffr::expect_doppelganger(
    "regplot_categorical_brma_ggplot",
    regplot(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: Location-Scale Model Regression Plot
# ============================================================================ #

test_that("Regression plot for location-scale model works", {

  name        <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual tests
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_scale_brma_base", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::regplot(fit_metafor, mod = "ni100", main = "metafor", ylim = c(-1, 1.5))
    regplot(fit_brma, mod = "ni100", plot_type = "base", main = "brma", ylim = c(-1, 1.5))
  })

  vdiffr::expect_doppelganger(
    "regplot_scale_brma_ggplot",
    regplot(fit_brma, plot_type = "ggplot", mod = "ni100")
  )
})


# ============================================================================ #
# Test: Regression Plot Options
# ============================================================================ #

test_that("Regression plot has correct interface", {

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_data = TRUE returns list with expected components
  # --------------------------------------------------

  regplot_data <- regplot(fit_brma, mod = "year", as_data = TRUE)

  expect_true(is.list(regplot_data),
    info = "as_data = TRUE should return a list"
  )

  expected_components <- c(
    "points", "pred", "ci", "xlim", "ylim", "xlab", "ylab", "mod_type"
  )
  expect_true(all(expected_components %in% names(regplot_data)),
    info = "regplot data should contain all expected components"
  )

  # Check points data.frame structure
  expect_true(is.data.frame(regplot_data$points),
    info = "points should be a data.frame"
  )
  expect_true(all(c("x", "y", "size") %in% names(regplot_data$points)),
    info = "points should have x, y, and size columns"
  )

  # Check continuous interval-band structure
  band_columns <- c("x", "y", "lower", "upper", "xpred")
  n_pred       <- nrow(regplot_data$pred)

  for (band_name in c("ci", "pi", "si")) {
    band_data <- regplot_data[[band_name]]

    expect_true(all(band_columns %in% names(band_data)),
      info = paste0(band_name, " should have x, y, lower, upper, and xpred columns")
    )
    expect_equal(nrow(band_data), 2 * n_pred,
      info = paste0(band_name, " should contain one row per polygon vertex")
    )
    expect_true(all(vapply(band_data[band_columns], length, integer(1)) == nrow(band_data)),
      info = paste0(band_name, " columns should all have matching lengths")
    )
    expect_equal(band_data$x, band_data$xpred,
      info = paste0(band_name, " x coordinates should match the stored prediction x values")
    )
  }

  # Check number of points matches number of studies
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(regplot_data$points), n_studies,
    info = "number of points should match number of studies"
  )

  # --------------------------------------------------
  # Test error on intercept-only model
  # --------------------------------------------------

  name_simple <- "bcg_meta-analysis"
  fit_simple  <- fits[[name_simple]]

  expect_error(regplot(fit_simple),
    regexp = "regplot requires a model with moderators",
    info = "should error on intercept-only model"
  )

  # --------------------------------------------------
  # Test error on invalid plot_type
  # --------------------------------------------------

  expect_error(regplot(fit_brma, plot_type = "invalid"),
    info = "should error on invalid plot_type"
  )
})


# ============================================================================ #
# Test: Regression Plot Customization
# ============================================================================ #

test_that("Regression plot customization works", {

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test custom point aesthetics
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_custom_points_base", function() {
    regplot(fit_brma, mod = "year", plot_type = "base", pch = 21, col = "blue", bg = "lightblue")
  })

  vdiffr::expect_doppelganger(
    "regplot_custom_points_ggplot",
    regplot(fit_brma, mod = "year", plot_type = "ggplot", pch = 21, col = "blue", bg = "lightblue")
  )

  # --------------------------------------------------
  # Test custom line and band styling
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_custom_bands_base", function() {
    regplot(fit_brma, mod = "year", plot_type = "base", lcol = "darkblue", col.ci = "lightblue",
            pi = TRUE, col.pi = "lightyellow")
  })

  vdiffr::expect_doppelganger(
    "regplot_custom_bands_ggplot",
    regplot(fit_brma, mod = "year", plot_type = "ggplot", lcol = "darkblue", col.ci = "lightblue",
            pi = TRUE, col.pi = "lightyellow")
  )

  # --------------------------------------------------
  # Test reference line
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_refline_base", function() {
    regplot(fit_brma, mod = "year", plot_type = "base", refline = 0)
  })

  vdiffr::expect_doppelganger(
    "regplot_refline_ggplot",
    regplot(fit_brma, mod = "year", plot_type = "ggplot", refline = 0)
  )

  # --------------------------------------------------
  # Test custom axis labels and title
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_custom_labels_base", function() {
    regplot(fit_brma, mod = "ablat", plot_type = "base",
            xlab = "Absolute Latitude", ylab = "Log Risk Ratio",
            main = "BCG Meta-Regression")
  })

  vdiffr::expect_doppelganger(
    "regplot_custom_labels_ggplot",
    regplot(fit_brma, mod = "ablat", plot_type = "ggplot",
            xlab = "Absolute Latitude", ylab = "Log Risk Ratio",
            main = "BCG Meta-Regression")
  )

  # --------------------------------------------------
  # Test no shading
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_no_shade_base", function() {
    regplot(fit_brma, mod = "ablat", plot_type = "base", shade = FALSE)
  })

  vdiffr::expect_doppelganger(
    "regplot_no_shade_ggplot",
    regplot(fit_brma, mod = "ablat", plot_type = "ggplot", shade = FALSE)
  )

  # --------------------------------------------------
  # Test custom axis limits
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_custom_limits_base", function() {
    regplot(fit_brma, mod = "ablat", plot_type = "base", xlim = c(0, 60), ylim = c(-3, 1))
  })

  vdiffr::expect_doppelganger(
    "regplot_custom_limits_ggplot",
    regplot(fit_brma, mod = "ablat", plot_type = "ggplot", xlim = c(0, 60), ylim = c(-3, 1))
  )

  # --------------------------------------------------
  # Test CI only (no PI)
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_ci_only_base", function() {
    regplot(fit_brma, mod = "ablat", plot_type = "base", ci = TRUE, pi = FALSE)
  })

  # --------------------------------------------------
  # Test neither CI nor PI
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_no_bands_base", function() {
    regplot(fit_brma, mod = "ablat", plot_type = "base", ci = FALSE, pi = FALSE)
  })
})


# ============================================================================ #
# Test: Regression Plot Confidence Level
# ============================================================================ #

test_that("Regression plot confidence level works", {

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test different confidence levels
  # --------------------------------------------------

  vdiffr::expect_doppelganger("regplot_level_90_base", function() {
    regplot(fit_brma, mod = "ablat", plot_type = "base", level = 90, main = "90% CI")
  })

  vdiffr::expect_doppelganger("regplot_level_99_base", function() {
    regplot(fit_brma, mod = "ablat", plot_type = "base", level = 99, main = "99% CI")
  })
})
