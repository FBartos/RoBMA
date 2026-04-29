context("Z-curve plot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)


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
  zc       <- as_zcurve(fits[[name]])

  # --------------------------------------------------
  # Test custom styles
  # --------------------------------------------------

  vdiffr::expect_doppelganger("zcurve_custom_base", function() {
    suppressMessages(plot(
      zc, plot_type = "base",
      plot_fit = TRUE, plot_ci = TRUE,
      lwd = 2, lty = 2,           # line args
      dots_hist = list(col = "lightblue"),
      main = "Custom Z-Curve"
    ))
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

test_that("Z-curve handles RoBMA bias-mixture branches", {

  name <- "dat.lehmann2018_RoBMA"
  skip_if_not(name %in% names(fits), "RoBMA cached fit not available.")

  fit               <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit[["fit"]])
  if (nrow(posterior_samples) > 50) {
    selected_ind      <- round(seq(from = 1, to = nrow(posterior_samples), length.out = 50))
    posterior_samples <- posterior_samples[selected_ind, , drop = FALSE]
  }
  selection         <- .zcurve_selection_context(
    object            = fit,
    data              = fit[["data"]],
    priors            = fit[["priors"]],
    posterior_samples = posterior_samples,
    is_weightfunction = .is_weightfunction(fit)
  )
  weighted_rows <- which(!selection[["use_normal"]])
  skip_if_not(length(weighted_rows) > 0, "No weightfunction posterior rows in cached fit.")

  selection_args <- .zcurve_selection_args(
    selection = selection,
    row       = weighted_rows[1],
    estimate  = 1,
    n         = 3
  )
  branch      <- selection[["bias_indicator"]][weighted_rows[1]]
  active_cuts <- selection[["fit_data"]][["crit_yi_mapping_max"]][branch]

  expect_equal(nrow(selection_args[["omega"]]), 3)
  expect_equal(ncol(selection_args[["omega"]]), active_cuts + 1)
  expect_length(selection_args[["crit_yi"]], active_cuts)

  zc <- as_zcurve(fit, max_samples = 50)
  expect_true(all(is.finite(zc[["zcurve"]][["estimates"]][["EDR"]])))
  expect_true(all(zc[["zcurve"]][["estimates"]][["EDR"]] >= 0))
  expect_true(all(zc[["zcurve"]][["estimates"]][["EDR"]] <= 1))
  expect_true(all(is.finite(zc[["zcurve"]][["estimates"]][["weights"]])))

  fitted_density <- lines(
    zc, as_data = TRUE, max_samples = 50, plot_ci = FALSE,
    extrapolate = FALSE, length.out = 25
  )
  extrapolated_density <- lines(
    zc, as_data = TRUE, max_samples = 50, plot_ci = FALSE,
    extrapolate = TRUE, length.out = 25
  )

  expect_true(all(is.finite(unlist(fitted_density[c("y", "y_lCI", "y_uCI")]))))
  expect_true(all(is.finite(unlist(extrapolated_density[c("y", "y_lCI", "y_uCI")]))))
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

