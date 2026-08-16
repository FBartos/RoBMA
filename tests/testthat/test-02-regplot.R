context("Regression plot (bubble plot)")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))

test_that("random-slope regplot variance replays its prediction grid", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    x     = c(0, 1, 0, 1),
    study = c("s1", "s1", "s2", "s2")
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.04, nrow(dat)),
    mods                      = ~ x,
    scale                     = ~ x,
    random                    = ~ diag(1 + x | study),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- matrix(
    c(
      0.20, 0,      0.25, 0.75,
      0.30, log(2), 0.40, 0.60
    ),
    nrow     = 2L,
    byrow    = TRUE,
    dimnames = list(NULL, c(
      "log_tau_intercept",
      "log_tau_x",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]"
    ))
  )
  grid_data <- data.frame(x = c(0, 1))

  tau <- rbind(c(0.20, 0.20), c(0.30, 0.60))
  z   <- (grid_data[["x"]] - mean(dat[["x"]])) / stats::sd(dat[["x"]])
  allocation <- rbind(c(0.50, 1.50), c(0.80, 1.20))
  expected <- tau * sqrt(
    allocation[, 1L] + tcrossprod(allocation[, 2L], z^2)
  )

  out <- .regplot_mv_random_sd_samples(
    x                 = object,
    newdata           = grid_data,
    posterior_samples = posterior_samples
  )

  expect_equal(unname(out), unname(expected), tolerance = 1e-12)
})


test_that("regplot prediction precision follows or overrides precision moderators", {

  vi <- c(0.01, 0.04)
  vi_grid <- .regplot_prediction_precision(
    newdata       = data.frame(vi = vi),
    mod_name      = "vi",
    n_pred        = length(vi),
    default_sei   = 1,
    reference_sei = NULL
  )
  expect_equal(vi_grid[["sei"]], sqrt(vi))
  expect_equal(vi_grid[["newdata"]][["vi"]], vi)

  sei <- c(0.1, 0.2)
  sei_grid <- .regplot_prediction_precision(
    newdata       = data.frame(sei = sei),
    mod_name      = "sei",
    n_pred        = length(sei),
    default_sei   = 1,
    reference_sei = NULL
  )
  expect_equal(sei_grid[["sei"]], sei)
  expect_equal(sei_grid[["newdata"]][["sei"]], sei)

  override <- .regplot_prediction_precision(
    newdata       = data.frame(year = 1:2, vi = vi),
    mod_name      = "year",
    n_pred        = length(vi),
    default_sei   = 1,
    reference_sei = 0.3
  )
  expect_equal(override[["sei"]], rep(0.3, length(vi)))
  expect_equal(override[["newdata"]][["vi"]], rep(0.3^2, length(vi)))
})

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

.test_regplot <- function(..., max_samples = 1000) {

  regplot(..., max_samples = max_samples)
}

# ============================================================================ #
# Test: Continuous Moderator Regression Plot
# ============================================================================ #

test_that("Regression plot for continuous moderator matches metafor", {

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_continuous_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::regplot(fit_metafor, mod = "year", main = "metafor", ylim = c(-2, 0.5))
    .test_regplot(fit_brma, mod = "year", plot_type = "base", main = "brma", ylim = c(-2, 0.5))
  })

  expect_vdiffr_snapshot(
    "regplot_continuous_brma_ggplot",
    .test_regplot(fit_brma, mod = "year", plot_type = "ggplot")
  )
})

test_that("Regression plot renders prediction intervals", {

  skip_if_not_full_visuals("Prediction-interval regplot variants are gallery coverage.")

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Visual tests with PI bands
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_continuous_pi_base", function() {
    .test_regplot(fit_brma, mod = "year", plot_type = "base", pi = TRUE, main = "with PI")
  })

  expect_vdiffr_snapshot(
    "regplot_continuous_pi_ggplot",
    .test_regplot(fit_brma, mod = "year", plot_type = "ggplot", pi = TRUE)
  )
})

test_that("Regression plot renders categorical dummy moderator", {

  name        <- "bcg_meta-regression2"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  expect_vdiffr_snapshot("regplot_categorical_brma", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    .test_regplot(fit_brma, mod = "alloc", plot_type = "base", main = "brma")
  })

  expect_vdiffr_snapshot(
    "regplot_categorical_brma_ggplot",
    .test_regplot(fit_brma, plot_type = "ggplot")
  )
})

test_that("Regression plot renders alternative categorical coding", {

  skip_if_not_full_visuals("Alternative coding snapshots duplicate categorical regplot coverage.")

  fit_brma1 <- fits[["bcg_meta-regression2"]]
  fit_brma2 <- fits[["bcg_meta-regression2b"]]

  expect_vdiffr_snapshot("regplot_categorical_coding_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    .test_regplot(fit_brma1, mod = "alloc", plot_type = "base", main = "brma", si = TRUE)
    .test_regplot(fit_brma2, mod = "alloc", plot_type = "base", main = "brma", si = TRUE)
  })
})

test_that("Regression plot renders interaction moderators", {

  skip_if_not_full_visuals("Interaction regplot variants are gallery coverage.")

  fit_brma3 <- fits[["bcg_meta-regression3"]]
  fit_brma4 <- fits[["bcg_meta-regression4"]]

  expect_vdiffr_snapshot("regplot_interaction_coding_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    .test_regplot(fit_brma3, mod = "year",            plot_type = "base", by = "alloc")
    .test_regplot(fit_brma4, mod = "year_before1969", plot_type = "base", by = "alloc")
  })

  expect_vdiffr_snapshot(
    "regplot_interaction_coding_ggplot",
    .test_regplot(fit_brma4, mod = "year_before1969", plot_type = "ggplot", by = "alloc")
  )
})

# ============================================================================ #
# Test: Location-Scale Model Regression Plot
# ============================================================================ #

test_that("Regression plot for location-scale model matches metafor view", {

  skip_if_not_full_visuals("Location-scale regplot variants are gallery coverage.")

  name        <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual tests
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_scale_brma_base", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::regplot(fit_metafor, mod = "ni100", main = "metafor", ylim = c(-1, 1.5))
    .test_regplot(fit_brma, mod = "ni100", plot_type = "base", main = "brma", ylim = c(-1, 1.5))
  })

  expect_vdiffr_snapshot(
    "regplot_scale_brma_ggplot",
    .test_regplot(fit_brma, plot_type = "ggplot", mod = "ni100")
  )
})

test_that("Regression plot thins location-scale prediction intervals", {

  name <- "bangertdrowns2004_location-scale"
  skip_if_missing_fits(name)

  fit_brma     <- fits[[name]]
  regplot_data <- NULL

  expect_no_error({
    regplot_data <- .test_regplot(
      fit_brma,
      mod         = "ni100",
      pi          = TRUE,
      as_data     = TRUE,
      at          = c(0, 3, 6),
      max_samples = 200
    )
  })

  expect_equal(nrow(regplot_data[["pi"]]), 2 * nrow(regplot_data[["pred"]]))
  expect_true(all(is.finite(regplot_data[["pi"]][["lower"]])))
  expect_true(all(is.finite(regplot_data[["pi"]][["upper"]])))
})

test_that("Regression plot supports scale-only moderators", {

  name <- "dat.lehmann2018_BMA.norm_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  regplot_data <- .test_regplot(
    fit_brma,
    mod         = "Preregistered",
    si          = TRUE,
    as_data     = TRUE,
    max_samples = 200
  )

  expect_identical(regplot_data[["mod_name"]], "Preregistered")
  expect_equal(
    regplot_data[["pred"]][["y"]],
    rep(regplot_data[["pred"]][["y"]][1L], nrow(regplot_data[["pred"]]))
  )
  expect_equal(
    regplot_data[["ci"]][["lower"]],
    rep(regplot_data[["ci"]][["lower"]][1L], nrow(regplot_data[["ci"]]))
  )
  expect_equal(
    regplot_data[["ci"]][["upper"]],
    rep(regplot_data[["ci"]][["upper"]][1L], nrow(regplot_data[["ci"]]))
  )
  si_width <- regplot_data[["si"]][["upper"]] -
    regplot_data[["si"]][["lower"]]
  expect_gt(diff(range(si_width)), 0)
})

test_that("Regression plot supports pointwise brma.mv random-formula PI and SI", {

  name <- "brma.mv_block_mvn_random_mods_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  regplot_data <- .test_regplot(
    fit_brma,
    mod         = "x",
    pi          = TRUE,
    si          = TRUE,
    as_data     = TRUE,
    at          = c(0, 1),
    max_samples = 80
  )

  expect_equal(nrow(regplot_data[["pi"]]), 2 * nrow(regplot_data[["pred"]]))
  expect_equal(nrow(regplot_data[["si"]]), 2 * nrow(regplot_data[["pred"]]))
  expect_true(all(is.finite(regplot_data[["pi"]][["lower"]])))
  expect_true(all(is.finite(regplot_data[["pi"]][["upper"]])))
  expect_true(all(is.finite(regplot_data[["si"]][["lower"]])))
  expect_true(all(is.finite(regplot_data[["si"]][["upper"]])))
  expect_equal(
    regplot_data[["sei"]],
    stats::median(.outcome_data_sei(fit_brma))
  )
})

# ============================================================================ #
# Test: 3-Level Model Regression Plot
# ============================================================================ #

test_that("Regression plot for multilevel model renders fitted moderator", {

  name        <- "konstantopoulos2011_3lvl2"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  mod         <- "vi"

  expect_no_error({
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::regplot(fit_metafor, mod = mod, main = "metafor", ylim = c(-1, 1.5))
    .test_regplot(fit_brma, mod = mod, plot_type = "base", main = "brma", ylim = c(-1, 1.5))
  })

  expect_s3_class(.test_regplot(fit_brma, plot_type = "ggplot", mod = mod), "ggplot")

  plot_data <- .test_regplot(
    fit_brma,
    mod      = mod,
    si       = TRUE,
    as_data  = TRUE
  )
  expect_equal(plot_data[["sei"]], sqrt(plot_data[["pred"]][["x"]]))
  expect_true(all(is.finite(plot_data[["si"]][["lower"]])))
  expect_true(all(is.finite(plot_data[["si"]][["upper"]])))
  expect_error(
    .test_regplot(fit_brma, mod = mod, sei = 0.1),
    "cannot be supplied"
  )

  sei_grid <- c(0.1, 0.2)
  sei_newdata <- .regplot_add_dummy_outcome(
    x       = fit_brma,
    newdata = data.frame(sei = sei_grid),
    n_pred  = length(sei_grid),
    sei     = 1
  )
  expect_equal(sei_newdata[["sei"]], sei_grid)
  expect_equal(sei_newdata[["vi"]], sei_grid^2)
})

# ============================================================================ #
# Test: Selection Regression Plot
# ============================================================================ #

test_that("Regression plot for selection model renders fitted moderator", {

  name        <- "dat.lehmann2018-3PSMreg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  mod         <- "Preregistered"

  expect_no_error({
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::regplot(fit_metafor, mod = mod, main = "metafor", ylim = c(-1, 1.5))
    .test_regplot(fit_brma, mod = mod, plot_type = "base", main = "brma", ylim = c(-1, 1.5))
  })

  expect_s3_class(.test_regplot(fit_brma, plot_type = "ggplot", mod = mod), "ggplot")
})

# ============================================================================ #
# Test: PET Model Regression Plot
# ============================================================================ #

test_that("Regression plot for PET model renders fitted moderator", {

  name        <- "dat.lehmann2018-PETreg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  mod         <- "Preregistered"

  expect_no_error({
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::regplot(fit_metafor, mod = mod, main = "metafor", ylim = c(-1, 1.5))
    .test_regplot(fit_brma, mod = mod, plot_type = "base", main = "brma", ylim = c(-1, 1.5))
  })

  expect_s3_class(.test_regplot(fit_brma, plot_type = "ggplot", mod = mod), "ggplot")
})

# ============================================================================ #
# Test: BMA.norm Model Regression Plot
# ============================================================================ #

test_that("Regression plot for BMA.norm meta-regression renders base output", {

  name     <- "dat.lehmann2018_BMA.norm_mods"
  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("regplot_BMAreg", function() {
    suppressWarnings(.test_regplot(fit_brma, plot_type = "base"))
  })
})

# ============================================================================ #
# Test: BMA.glmm Model Regression Plot
# ============================================================================ #

test_that("Regression plot for BMA.glmm model renders base output", {

  skip_if_not_full_visuals("BMA.glmm regplot smoke is gallery coverage.")

  name     <- "bcg_BMA.glmm_3lvl_location_scale"
  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("regplot_BMA.glmm", function() {
    suppressWarnings(.test_regplot(fit_brma, mod = "year", plot_type = "base"))
  })
})

# ============================================================================ #
# Test: RoBMA Model Regression Plot
# ============================================================================ #

test_that("Regression plot for RoBMA meta-regression renders base output", {

  name     <- "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("regplot_RoBMA_complex", function() {
    suppressWarnings(.test_regplot(fit_brma, mod = "Preregistered", plot_type = "base"))
  })
})

# ============================================================================ #
# Test: Regression Plot Options
# ============================================================================ #

test_that("Regression plot data and argument validation are stable", {

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_data = TRUE returns list with expected components
  # --------------------------------------------------

  regplot_data <- .test_regplot(fit_brma, mod = "year", pi = TRUE, si = TRUE, as_data = TRUE)
  expect_equal(regplot_data[["sei"]], stats::median(RoBMA:::.outcome_data_sei(fit_brma)))

  regplot_data_custom_sei <- .test_regplot(
    fit_brma,
    mod     = "year",
    si      = TRUE,
    sei     = .123,
    as_data = TRUE
  )
  expect_equal(regplot_data_custom_sei[["sei"]], .123)

  expect_true(is.list(regplot_data),
    info = "as_data = TRUE returns a list"
  )

  expected_components <- c(
    "points", "pred", "ci", "xlim", "ylim", "xlab", "ylab", "mod_type"
  )
  expect_true(all(expected_components %in% names(regplot_data)),
    info = "regplot data contains all expected components"
  )

  # Check points data.frame structure
  expect_true(is.data.frame(regplot_data$points),
    info = "points are returned as a data.frame"
  )
  expect_true(all(c("x", "y", "size") %in% names(regplot_data$points)),
    info = "points contain x, y, and size columns"
  )

  # Check continuous interval-band structure
  band_columns <- c("x", "y", "lower", "upper", "xpred")
  n_pred       <- nrow(regplot_data$pred)

  for (band_name in c("ci", "pi", "si")) {
    band_data <- regplot_data[[band_name]]

    expect_true(all(band_columns %in% names(band_data)),
      info = paste0(band_name, " contains x, y, lower, upper, and xpred columns")
    )
    expect_equal(nrow(band_data), 2 * n_pred,
      info = paste0(band_name, " contains one row per polygon vertex")
    )
    expect_true(all(vapply(band_data[band_columns], length, integer(1)) == nrow(band_data)),
      info = paste0(band_name, " columns all have matching lengths")
    )
    expect_equal(band_data$x, band_data$xpred,
      info = paste0(band_name, " x coordinates match the stored prediction x values")
    )
  }

  # Check number of points matches number of studies
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(regplot_data$points), n_studies,
    info = "number of points matches number of studies"
  )

  # --------------------------------------------------
  # Test error on intercept-only model
  # --------------------------------------------------

  name_simple <- "bcg_meta-analysis"
  fit_simple  <- fits[[name_simple]]

  expect_error(.test_regplot(fit_simple),
    regexp = "regplot requires a model with moderators",
    info = "intercept-only model is rejected"
  )

  # --------------------------------------------------
  # Test error on invalid plot_type
  # --------------------------------------------------

  expect_error(.test_regplot(fit_brma, plot_type = "invalid"),
    info = "invalid plot_type is rejected"
  )

  expect_error(.test_regplot(fit_brma, sei = 0),
    info = "non-positive reference sei is rejected"
  )

  expect_error(.test_regplot(fit_brma, mod = "year", at = "bad"),
    regexp = "at",
    info = "non-numeric prediction grid is rejected"
  )
  expect_error(.test_regplot(fit_brma, mod = "year", plim = c(3, 1)),
    regexp = "plim",
    info = "non-increasing point-size limits are rejected"
  )
  expect_error(.test_regplot(fit_brma, mod = "year", psize = 1:2),
    regexp = "psize",
    info = "wrong-length point sizes are rejected"
  )
  expect_error(.test_regplot(fit_brma, mod = "year", lwd = 0),
    regexp = "lwd",
    info = "non-positive line width is rejected"
  )

  fit_categorical <- fits[["bcg_meta-regression2"]]
  expect_error(.test_regplot(fit_categorical, mod = "alloc", at = 1),
    regexp = "continuous moderators",
    info = "at is rejected for categorical moderators"
  )
})

test_that("Regression plot transforms effect-size data", {

  skip_if_not_full_visuals("Transformation regplot snapshots are gallery coverage.")

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("regplot_transform_exp_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    .test_regplot(fit_brma, mod = "year", plot_type = "base", main = "log RR")
    .test_regplot(fit_brma, mod = "year", transf = exp, plot_type = "base", main = "RR")
  })
})

# ============================================================================ #
# Test: Regression Plot Customization
# ============================================================================ #

test_that("Regression plot customization snapshots are stable", {

  skip_if_not_full_visuals("Customization snapshots are visual-gallery coverage.")

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test custom point aesthetics
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_custom_points_base", function() {
    .test_regplot(fit_brma, mod = "year", plot_type = "base", pch = 21, col = "blue", bg = "lightblue")
  })

  expect_vdiffr_snapshot(
    "regplot_custom_points_ggplot",
    .test_regplot(fit_brma, mod = "year", plot_type = "ggplot", pch = 21, col = "blue", bg = "lightblue")
  )

  # --------------------------------------------------
  # Test custom line and band styling
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_custom_bands_base", function() {
    .test_regplot(fit_brma, mod = "year", plot_type = "base", lcol = "darkblue", col.ci = "lightblue",
            pi = TRUE, col.pi = "lightyellow")
  })

  expect_vdiffr_snapshot(
    "regplot_custom_bands_ggplot",
    .test_regplot(fit_brma, mod = "year", plot_type = "ggplot", lcol = "darkblue", col.ci = "lightblue",
            pi = TRUE, col.pi = "lightyellow")
  )

  # --------------------------------------------------
  # Test reference line
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_refline_base", function() {
    .test_regplot(fit_brma, mod = "year", plot_type = "base", refline = 0)
  })

  expect_vdiffr_snapshot(
    "regplot_refline_ggplot",
    .test_regplot(fit_brma, mod = "year", plot_type = "ggplot", refline = 0)
  )

  # --------------------------------------------------
  # Test custom axis labels and title
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_custom_labels_base", function() {
    .test_regplot(fit_brma, mod = "ablat", plot_type = "base",
            xlab = "Absolute Latitude", ylab = "Log Risk Ratio",
            main = "BCG Meta-Regression")
  })

  expect_vdiffr_snapshot(
    "regplot_custom_labels_ggplot",
    .test_regplot(fit_brma, mod = "ablat", plot_type = "ggplot",
            xlab = "Absolute Latitude", ylab = "Log Risk Ratio",
            main = "BCG Meta-Regression")
  )

  # --------------------------------------------------
  # Test no shading
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_no_shade_base", function() {
    .test_regplot(fit_brma, mod = "ablat", plot_type = "base", shade = FALSE)
  })

  expect_vdiffr_snapshot(
    "regplot_no_shade_ggplot",
    .test_regplot(fit_brma, mod = "ablat", plot_type = "ggplot", shade = FALSE)
  )

  # --------------------------------------------------
  # Test custom axis limits
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_custom_limits_base", function() {
    .test_regplot(fit_brma, mod = "ablat", plot_type = "base", xlim = c(0, 60), ylim = c(-3, 1))
  })

  expect_vdiffr_snapshot(
    "regplot_custom_limits_ggplot",
    .test_regplot(fit_brma, mod = "ablat", plot_type = "ggplot", xlim = c(0, 60), ylim = c(-3, 1))
  )

  # --------------------------------------------------
  # Test CI only (no PI)
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_ci_only_base", function() {
    .test_regplot(fit_brma, mod = "ablat", plot_type = "base", ci = TRUE, pi = FALSE)
  })

  # --------------------------------------------------
  # Test neither CI nor PI
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_no_bands_base", function() {
    .test_regplot(fit_brma, mod = "ablat", plot_type = "base", ci = FALSE, pi = FALSE)
  })
})

# ============================================================================ #
# Test: Regression Plot Confidence Level
# ============================================================================ #

test_that("Regression plot confidence-level snapshots are stable", {

  skip_if_not_full_visuals("Confidence-level snapshots are visual-gallery coverage.")

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test different confidence levels
  # --------------------------------------------------

  expect_vdiffr_snapshot("regplot_level_90_base", function() {
    .test_regplot(fit_brma, mod = "ablat", plot_type = "base", level = 90, main = "90% CI")
  })

  expect_vdiffr_snapshot("regplot_level_99_base", function() {
    .test_regplot(fit_brma, mod = "ablat", plot_type = "base", level = 99, main = "99% CI")
  })
})

test_that("native regplot mixture intervals match R fallback", {

  skip_if_not(.has_native_regplot_mixture())

  set.seed(1024)

  S     <- 200
  K     <- 4
  probs <- c(.025, .975)

  mean_samples <- matrix(stats::rnorm(S * K), nrow = S, ncol = K)
  sd_samples   <- matrix(stats::runif(S * K, .05, 1.25), nrow = S, ncol = K)

  expected <- .regplot_mixture_interval_quantiles_r(
    mean_samples = mean_samples,
    sd_samples   = sd_samples,
    probs        = probs
  )
  actual <- RoBMA:::.regplot_mixture_interval_quantiles(
    mean_samples = mean_samples,
    sd_samples   = sd_samples,
    probs        = probs
  )

  expect_equal(actual, expected, tolerance = 1e-4)
})
