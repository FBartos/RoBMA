context("Z-curve plot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list & load all fits
skip_if_no_fits()
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


# ============================================================================ #
# Test: Simple Meta-Analysis Z-Curve
# ============================================================================ #

test_that("Z-curve plot for simple meta-analysis works", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]
  zc       <- as_zcurve(fit_brma)

  # --------------------------------------------------
  # Visual tests
  # --------------------------------------------------

  vdiffr::expect_doppelganger("zcurve_simple_base", function() {
    suppressMessages(plot(zc, plot_type = "base"))
  })

  vdiffr::expect_doppelganger(
    "zcurve_simple_ggplot",
    suppressMessages(plot(zc, plot_type = "ggplot"))
  )
})


# ============================================================================ #
# Test: Z-Curve Customization
# ============================================================================ #

test_that("Z-curve plot customization works", {

  name     <- "bcg_meta-analysis"
  zc       <- as_zcurve(fits[[name]]$models[[1]])

  # --------------------------------------------------
  # Test custom styles
  # --------------------------------------------------

  vdiffr::expect_doppelganger("zcurve_custom_base", function() {
    plot(zc, plot_type = "base",
         plot_fit = TRUE, plot_CI = TRUE,
         col = "blue", lwd = 2, lty = 2,           # line args
         dots_hist = list(fill = "lightblue", color = "blue"),
         main = "Custom Z-Curve"
    )
  })

  # For ggplot, we pass specific list args for hist/lines as per .get_dots_* in implementation
  # or global args that get mapped
  vdiffr::expect_doppelganger(
    "zcurve_custom_ggplot",
    suppressMessages(plot(
      zc, plot_type = "ggplot",
      dots_hist = list(fill = "lightblue", color = "blue"),
      dots_thresholds = list(color = "red", linetype = "dashed"),
      main = "Custom Z-Curve GGplot"
    ))

  )

  # --------------------------------------------------
  # Test components only (hist / lines)
  # --------------------------------------------------

  vdiffr::expect_doppelganger("zcurve_hist_only_base", function() {
    suppressMessages(hist(zc, plot_type = "base", main = "Hist Only"))
  })

  vdiffr::expect_doppelganger("zcurve_lines_only_base", function() {
    # lines() adds to existing plot usually, but here we test the function
    # so we create an empty plot and add lines
    plot(0, 0, type = "n", xlim = c(-6, 6), ylim = c(0, 0.5), main = "Lines Only")
    lines(zc, plot_type = "base", col = "purple")
  })

})


# ============================================================================ #
# Test: Meta-Regression Z-Curve
# ============================================================================ #

test_that("Z-curve plot for meta-regression works", {

  name <- "bcg_meta-regression"
  zc   <- as_zcurve(fits[[name]])

  vdiffr::expect_doppelganger("zcurve_regression_base", function() {
    suppressMessages(plot(zc, plot_type = "base", main = "Meta-Regression Z-Curve"))
  })

})

# ============================================================================ #
# Test: Selection Models Z-Curve
# ============================================================================ #

test_that("Z-curve plot for selection model (positive) works", {

  name <- "dat.lehmann2018-3PSM"
  zc   <- as_zcurve(fits[[name]])

  vdiffr::expect_doppelganger("zcurve_selection_pos_base", function() {
    suppressMessages(plot(zc, plot_type = "base", main = "Selection Model (Pos) Z-Curve"))
  })

})

test_that("Z-curve plot for selection model (negative) works", {

  name <- "dat.lehmann2018-3PSM_neg"
  zc   <- as_zcurve(fits[[name]])

  vdiffr::expect_doppelganger("zcurve_selection_neg_base", function() {
    suppressMessages(plot(zc, plot_type = "base", main = "Selection Model (Neg) Z-Curve"))
  })

})

# ============================================================================ #
# Test: PET Models Z-Curve
# ============================================================================ #

test_that("Z-curve plot for PET model (positive) works", {

  name <- "dat.lehmann2018-PET"
  zc   <- as_zcurve(fits[[name]])

  vdiffr::expect_doppelganger("zcurve_PET_pos_base", function() {
    suppressMessages(plot(zc, plot_type = "base", main = "PET Model (Pos) Z-Curve"))
  })

})

test_that("Z-curve plot for PET model (negative) works", {

  name <- "dat.lehmann2018-PET_neg"
  zc   <- as_zcurve(fits[[name]])

  vdiffr::expect_doppelganger("zcurve_PET_neg_base", function() {
    suppressMessages(plot(zc, plot_type = "base", main = "PET Model (Neg) Z-Curve"))
  })

})

# ============================================================================ #
# Test: Multilevel Models Z-Curve
# ============================================================================ #

test_that("Z-curve plot for multilevel model works", {

  name <- "konstantopoulos2011_3lvl"
  zc   <- as_zcurve(fits[[name]])

  vdiffr::expect_doppelganger("zcurve_multilevel_base", function() {
    suppressMessages(plot(zc, plot_type = "base", main = "Multilevel Model Z-Curve"))
  })

})

