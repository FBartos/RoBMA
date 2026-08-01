context("zplot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))

# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

zplot_smoke_samples <- test_profile_value(100L, 1000L)

.test_as_zplot <- function(..., max_samples = 1000) {

  as_zplot(..., max_samples = max_samples)
}

.test_zplot <- function(..., summary_max_samples = 1000, max_samples = 1000) {

  zplot(..., summary_max_samples = summary_max_samples, max_samples = max_samples)
}

.test_plot_zplot <- function(..., max_samples = 1000) {

  plot(..., max_samples = max_samples)
}

.test_lines_zplot <- function(..., max_samples = 1000) {

  lines(..., max_samples = max_samples)
}

.zplot_test_area <- function(df) {

  x <- df[["x"]]
  y <- df[["y"]]

  return(sum(diff(x) * (y[-length(y)] + y[-1L]) / 2))
}


test_that("zplot rejects GLMM fits early", {

  name <- "bcg_BMA.glmm"
  skip_if_missing_fits(name)

  expect_error(
    .test_as_zplot(fits[[name]]),
    "normal outcome models"
  )
})


test_that("zplot creates reusable objects and plots directly", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit <- fits[[name]]
  zp  <- .test_as_zplot(fit, max_samples = zplot_smoke_samples)

  expect_s3_class(zp, "zplot_brma")
  expect_named(zp[["zplot"]], c("estimates", "data"))
  expect_true(.is_ggplot(.test_zplot(
    fit,
    plot_type           = "ggplot",
    summary_max_samples = zplot_smoke_samples,
    max_samples         = zplot_smoke_samples
  )))
  expect_true(.is_ggplot(.test_zplot(
    zp,
    plot_type   = "ggplot",
    max_samples = zplot_smoke_samples
  )))
})

test_that("zplot allows descriptive known-V brma.mv displays", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  zc <- .test_as_zplot(
    fits[[name]],
    max_samples = zplot_smoke_samples
  )

  expect_s3_class(zc, "zplot_brma")
  expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["EDR"]])))
  expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["weights"]])))
  expect_equal(
    length(zc[["zplot"]][["data"]][["z"]]),
    nobs(fits[[name]])
  )

  density <- .test_lines_zplot(
    zc,
    as_data     = TRUE,
    max_samples = zplot_smoke_samples,
    plot_ci     = FALSE,
    length.out  = 25
  )
  expect_true(all(is.finite(density[["y"]])))
})


test_that("zplot known-V brma.mv uses estimate-target predictive components", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma          <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
  selected_rows     <- seq_len(min(25L, nrow(posterior_samples)))
  posterior_samples <- posterior_samples[selected_rows, , drop = FALSE]

  setup <- .estimate_likelihood_setup_from_parts(
    fit               = fit_brma[["fit"]],
    data              = fit_brma[["data"]],
    priors            = fit_brma[["priors"]],
    posterior_samples = posterior_samples,
    unit              = "estimate",
    data_hash         = .get_outcome_hash(fit_brma),
    bias_adjusted     = TRUE,
    caller            = "test"
  )
  components <- .zplot_predictive_components(
    object            = fit_brma,
    posterior_samples = posterior_samples,
    extrapolate       = TRUE
  )

  expect_equal(components[["mu"]], setup[["mu"]])
  expect_equal(components[["tau_within"]], .known_v_extra_sd_from_setup(setup))
  expect_equal(components[["sei"]], setup[["selection_sei"]])
})


# ============================================================================ #
# Test: Simple Meta-Analysis Zplot
# ============================================================================ #

test_that("zplot for simple meta-analysis renders base and ggplot output", {

  name     <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  zc       <- .test_as_zplot(fit_brma)

  # --------------------------------------------------
  # Visual tests
  # --------------------------------------------------

  expect_vdiffr_snapshot("zplot_simple_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base"))
  })

  expect_vdiffr_snapshot(
    "zplot_simple_ggplot",
    suppressMessages(.test_plot_zplot(zc, plot_type = "ggplot"))
  )
})


# ============================================================================ #
# Test: Zplot Customization
# ============================================================================ #

test_that("zplot customization snapshots are stable", {

  skip_if_not_full_visuals("Customization snapshots are visual-gallery coverage.")

  name     <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  zc       <- .test_as_zplot(fits[[name]])

  # --------------------------------------------------
  # Test custom styles
  # --------------------------------------------------

  expect_vdiffr_snapshot("zplot_custom_base", function() {
    suppressMessages(.test_plot_zplot(
      zc, plot_type = "base",
      plot_fit = TRUE, plot_ci = TRUE,
      lwd = 2, lty = 2,           # line args
      dots_hist = list(col = "lightblue"),
      main = "Custom Zplot"
    ))
  })

  # For ggplot, we pass specific list args for hist/lines as per .get_dots_* in implementation
  # or global args that get mapped
  expect_vdiffr_snapshot(
    "zplot_custom_ggplot",
    suppressMessages(.test_plot_zplot(
      zc, plot_type = "ggplot",
      dots_hist = list(fill = "lightblue", color = "blue"),
      dots_thresholds = list(color = "red", linetype = "dashed"),
      main = "Custom Zplot GGplot"
    ))

  )

  # --------------------------------------------------
  # Test components only (hist / lines)
  # --------------------------------------------------

  expect_vdiffr_snapshot("zplot_hist_only_base", function() {
    suppressMessages(hist(zc, plot_type = "base", main = "Hist Only"))
  })

  expect_vdiffr_snapshot("zplot_lines_only_base", function() {
    # lines() adds to existing plot usually, but here we test the function
    # so we create an empty plot and add lines
    plot(0, 0, type = "n", xlim = c(-6, 6), ylim = c(0, 0.5), main = "Lines Only")
    .test_lines_zplot(zc, plot_type = "base", col = "purple")
  })

})


# ============================================================================ #
# Test: Meta-Regression Zplot
# ============================================================================ #

test_that("zplot for meta-regression renders base output", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_regression_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "Meta-Regression Zplot"))
  })

})

# ============================================================================ #
# Test: Selection Models Zplot
# ============================================================================ #

test_that("zplot for positive-direction selection model renders base output", {

  name <- "dat.lehmann2018-3PSM"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_selection_pos_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "Selection Model (Pos) Zplot"))
  })

})

test_that("zplot for negative-direction selection model renders base output", {

  skip_if_not_full_visuals("Negative-direction selection zplot is gallery coverage.")

  name <- "dat.lehmann2018-3PSM_neg"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_selection_neg_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "Selection Model (Neg) Zplot"))
  })

})

test_that("zplot handles RoBMA bias-mixture branches", {

  name <- "dat.lehmann2018_RoBMA"
  skip_if_not(name %in% names(fits), "RoBMA cached fit not available.")

  max_samples        <- test_profile_value(100L, 1000L)
  integration_points <- test_profile_value(401L, 2001L)
  fit               <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit[["fit"]])
  if (nrow(posterior_samples) > max_samples) {
    selected_ind      <- round(seq(from = 1, to = nrow(posterior_samples), length.out = max_samples))
    posterior_samples <- posterior_samples[selected_ind, , drop = FALSE]
  }
  selection         <- .zplot_selection_context(
    object            = fit,
    posterior_samples = posterior_samples,
    is_weightfunction = .is_weightfunction(fit)
  )
  weighted_rows <- which(!selection[["use_normal"]])
  skip_if_not(length(weighted_rows) > 0, "No weightfunction posterior rows in cached fit.")

  selection_args <- .zplot_selection_args(
    selection = selection,
    row       = weighted_rows[1],
    estimate  = 1,
    n         = 3
  )
  active_cuts <- selection[["n_bins"]] - 1L

  expect_equal(nrow(selection_args[["omega"]]), 3)
  expect_equal(ncol(selection_args[["omega"]]), selection[["n_bins"]])
  expect_length(selection_args[["crit_yi"]], active_cuts)

  zc <- .test_as_zplot(fit, max_samples = max_samples)
  expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["EDR"]])))
  expect_true(all(zc[["zplot"]][["estimates"]][["EDR"]] >= 0))
  expect_true(all(zc[["zplot"]][["estimates"]][["EDR"]] <= 1))
  expect_true(all(is.finite(zc[["zplot"]][["estimates"]][["weights"]])))

  fitted_density <- .test_lines_zplot(
    zc, as_data = TRUE, max_samples = max_samples, plot_ci = FALSE,
    extrapolate = FALSE, length.out = 25
  )
  extrapolated_density <- .test_lines_zplot(
    zc, as_data = TRUE, max_samples = max_samples, plot_ci = FALSE,
    extrapolate = TRUE, length.out = 25
  )

  expect_true(all(is.finite(unlist(fitted_density[c("y", "y_lCI", "y_uCI")]))))
  expect_true(all(is.finite(unlist(extrapolated_density[c("y", "y_lCI", "y_uCI")]))))

  fitted_area       <- .zplot_test_area(.test_lines_zplot(
    zc,
    as_data     = TRUE,
    max_samples = max_samples,
    plot_ci     = FALSE,
    extrapolate = FALSE,
    from        = -20,
    to          = 20,
    length.out  = integration_points
  ))
  extrapolated_area <- .zplot_test_area(.test_lines_zplot(
    zc,
    as_data     = TRUE,
    max_samples = max_samples,
    plot_ci     = FALSE,
    extrapolate = TRUE,
    from        = -20,
    to          = 20,
    length.out  = integration_points
  ))
  expected_area     <- mean(zc[["zplot"]][["estimates"]][["weights"]])

  expect_equal(fitted_area, 1, tolerance = 0.01)
  expect_equal(extrapolated_area, expected_area, tolerance = 0.01)
  expect_gt(extrapolated_area, fitted_area)
})

# ============================================================================ #
# Test: PET Models Zplot
# ============================================================================ #

test_that("zplot for positive-direction PET model renders base output", {

  name <- "dat.lehmann2018-PET"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_PET_pos_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "PET Model (Pos) Zplot"))
  })

})

test_that("zplot for negative-direction PET model renders base output", {

  skip_if_not_full_visuals("Negative-direction PET zplot is gallery coverage.")

  name <- "dat.lehmann2018-PET_neg"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_PET_neg_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "PET Model (Neg) Zplot"))
  })

})

# ============================================================================ #
# Test: Multilevel Models Zplot
# ============================================================================ #

test_that("zplot for multilevel model renders base output", {

  skip_if_not_full_visuals("Multilevel zplot is gallery coverage.")

  name <- "konstantopoulos2011_3lvl"
  skip_if_missing_fits(name)

  zc   <- .test_as_zplot(fits[[name]])

  expect_vdiffr_snapshot("zplot_multilevel_base", function() {
    suppressMessages(.test_plot_zplot(zc, plot_type = "base", main = "Multilevel Model Zplot"))
  })

})



