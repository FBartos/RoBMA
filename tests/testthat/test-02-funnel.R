context("Funnel plot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)


# ============================================================================ #
# Test: Simple Meta-Analysis Funnel Plot
# ============================================================================ #

test_that("Funnel plot for simple meta-analysis matches metafor structure", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_simple_comparison_no_tau", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-3, 3), ylim = c(0, 0.8))
    funnel(fit_brma, plot_type = "base", xlim = c(-3, 3), ylim = c(0, 0.8), main = "brma", sampling_heterogeneity = FALSE)
  })

  vdiffr::expect_doppelganger("funnel_simple_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", addtau2 = TRUE, xlim = c(-3, 3), ylim = c(0, 0.8))
    funnel(fit_brma, plot_type = "base", xlim = c(-3, 3), ylim = c(0, 0.8), main = "brma")
  })

  vdiffr::expect_doppelganger(
    "funnel_simple_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: Meta-Regression Funnel Plot
# ============================================================================ #

test_that("Funnel plot for meta-regression works correctly", {

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_regression_comparison-standard", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstandard")
    funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstandard")
  })

  vdiffr::expect_doppelganger("funnel_regression_comparison-rstudent", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstudent")
    suppressWarnings(funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstudent"))
  })

  vdiffr::expect_doppelganger(
    "funnel_regression_brma_ggplot",
    suppressWarnings(funnel(fit_brma, plot_type = "ggplot"))
  )
})

test_that("Funnel plot for meta-regression with interaction works correctly", {

  name        <- "bcg_meta-regression4"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_regression4_comparison-standard", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstandard")
    suppressWarnings(funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstandard"))
  })

  vdiffr::expect_doppelganger("funnel_regression4_comparison-rstudent", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 1.2), xlim = c(-2, 2), type = "rstudent")
    suppressWarnings(funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 1.2), xlim = c(-2, 2), type = "rstudent"))
  })
})

# ============================================================================ #
# Test: Location-Scale Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for location-scale model works correctly", {

  name <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_scale_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 0.6), type = "rstandard")
    funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 0.6), sampling_heterogeneity = FALSE, type = "rstandard")
  })

  vdiffr::expect_doppelganger(
    "funnel_scale_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: 3-Level Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for 3-level model works correctly", {

  name <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_3lvl_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-1, 1.5), ylim = c(0.4, 0))
    funnel(fit_brma, plot_type = "base", main = "brma", sampling_heterogeneity = FALSE, xlim = c(-1, 1.5), ylim = c(0.4, 0))
  })

  # TODO:
  # add cluster-level residuals once implemented

  vdiffr::expect_doppelganger(
    "funnel_3lvl_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})

test_that("Funnel plot for 3-level meta-regression works correctly", {

  name <- "konstantopoulos2011_3lvl2"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_3lvl2_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-1, 1.5), ylim = c(0.5, 0), type = "rstandard")
    suppressWarnings(suppressWarnings(funnel(fit_brma, plot_type = "base", type = "rstandard", conditioning_depth = "marginal", main = "brma", sampling_heterogeneity = FALSE, xlim = c(-1, 1.5), ylim = c(0.5, 0))))
  })

  # TODO:
  # add cluster-level residuals once implemented

  vdiffr::expect_doppelganger(
    "funnel_3lvl2_brma_ggplot",
    suppressWarnings(funnel(fit_brma, plot_type = "ggplot"))
  )
})

# ============================================================================ #
# Test: GLMM Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for GLMM model works correctly", {

  name <- "nielweise2008_glmm"
  fit_brma <- fits[[name]]

  # there is no funnel plot for metafor
  vdiffr::expect_doppelganger(
    "funnel_glmm_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})

test_that("Funnel plot for GLMM meta-regression works correctly", {

  name <- "bcg_glmm_reg"
  fit_brma <- fits[[name]]

  # there is no funnel plot for metafor
  vdiffr::expect_doppelganger(
    "funnel_glmm_reg_ggplot",
    suppressWarnings(funnel(fit_brma, plot_type = "ggplot"))
  )
})

# ============================================================================ #
# Test: Selection Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for selection model works correctly", {

  name <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_selmodel_pos_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0))
    funnel(fit_brma, plot_type = "base", main = "brma", xlim = c(-2, 2), ylim = c(0.8, 0), sampling_bias = FALSE, sampling_heterogeneity = FALSE)
  })

  vdiffr::expect_doppelganger(
    "funnel_selmodel_pos_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot", sampling_bias = TRUE, sampling_heterogeneity = TRUE, xlim = c(-2, 2))
  )
})

test_that("Funnel plot for selection meta-regression works correctly", {

  name <- "dat.lehmann2018-3PSMreg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_selmodelreg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    # not available for metafor
    # par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    # metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0))
    suppressWarnings(funnel(fit_brma, plot_type = "base", main = "brma", xlim = c(-2, 2), ylim = c(0.8, 0), sampling_bias = FALSE, sampling_heterogeneity = FALSE))
  })
})

test_that("Funnel plot for selection model (negative) works correctly", {

  name <- "dat.lehmann2018-3PSM_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_selmodel_neg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0))
    funnel(fit_brma, plot_type = "base", main = "brma", xlim = c(-2, 2), ylim = c(0.8, 0), sampling_bias = FALSE, sampling_heterogeneity = FALSE)
  })
})

# ============================================================================ #
# Test: PET Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for PET model works correctly", {

  name <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_PET_pos_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    # the residuals need to be selected specifically because bPET is not treated as a regression
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0), type = "rstudent")
    suppressWarnings(funnel(fit_brma, plot_type = "base", sampling_bias = FALSE, sampling_heterogeneity = FALSE, residual = TRUE, type = "rstudent", xlim = c(-2, 2), ylim = c(0.8, 0)))
  })

  vdiffr::expect_doppelganger(
    "funnel_PET_pos_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})

test_that("Funnel plot for PET meta-regression works correctly", {

  name <- "dat.lehmann2018-PETreg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_PETreg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    # the residuals need to be selected specifically because bPET is not treated as a regression
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0), type = "rstudent")
    suppressWarnings(funnel(fit_brma, plot_type = "base", sampling_bias = FALSE, sampling_heterogeneity = FALSE, residual = TRUE, type = "rstudent", xlim = c(-2, 2), ylim = c(0.8, 0)))
  })
})

test_that("Funnel plot for PET model (negative) works correctly", {

  name <- "dat.lehmann2018-PET_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_PET_neg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    # the residuals need to be selected specifically because bPET is not treated as a regression
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0), type = "rstudent")
    suppressWarnings(funnel(fit_brma, plot_type = "base", sampling_bias = FALSE, sampling_heterogeneity = FALSE, residual = TRUE, type = "rstudent", xlim = c(-2, 2), ylim = c(0.8, 0)))
  })

  vdiffr::expect_doppelganger(
    "funnel_PET_neg_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: BMA.norm Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for BMA.norm model works correctly", {

  name     <- "dat.lehmann2018_BMA.norm"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("funnel_BMA", function() {
    suppressWarnings(funnel(fit_brma, plot_type = "base", sampling_heterogeneity = TRUE))
  })
})

test_that("Funnel plot for BMA.norm meta-regression works correctly", {

  name     <- "dat.lehmann2018_BMA.norm_mods"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("funnel_BMAreg", function() {
    suppressWarnings(funnel(fit_brma, plot_type = "base", sampling_heterogeneity = TRUE))
  })
})

# ============================================================================ #
# Test: BMA.glmm Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for BMA.glmm model works correctly", {

  name     <- "bcg_BMA.glmm_3lvl_location_scale"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("funnel_BMA.glmm", function() {
    suppressWarnings(funnel(fit_brma, plot_type = "base"))
  })
})

# ============================================================================ #
# Test: RoBMA Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for RoBMA model works correctly", {

  name     <- "dat.lehmann2018_RoBMA"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("funnel_RoBMA", function() {
    suppressWarnings(funnel(fit_brma, plot_type = "base", sampling_heterogeneity = TRUE, sampling_bias = TRUE))
  })
})

test_that("Funnel plot for RoBMA meta-regression works correctly", {

  name     <- "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("funnel_RoBMA_complex", function() {
    suppressWarnings(funnel(fit_brma, plot_type = "base", type = "LOO-PIT"))
  })
})

# ============================================================================ #
# Test: Funnel Plot Options
# ============================================================================ #

test_that("Funnel plot has correct interface", {

  name <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_data = TRUE returns list with expected components
  # --------------------------------------------------

  funnel_data <- funnel(fit_brma, as_data = TRUE)

  expect_true(is.list(funnel_data),
    info = "as_data = TRUE should return a list"
  )

  expected_components <- c(
    "points", "funnel", "funnel_edge1", "funnel_edge2",
    "background", "x_range", "y_range"
  )
  expect_true(all(expected_components %in% names(funnel_data)),
    info = "funnel data should contain all expected components"
  )

  # Check points data.frame structure
  expect_true(is.data.frame(funnel_data$points),
    info = "points should be a data.frame"
  )
  expect_true(all(c("x", "y") %in% names(funnel_data$points)),
    info = "points should have x and y columns"
  )

  # Check number of points matches number of studies
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(funnel_data$points), n_studies,
    info = "number of points should match number of studies"
  )

  # --------------------------------------------------
  # Test error on invalid plot_type
  # --------------------------------------------------

  expect_error(funnel(fit_brma, plot_type = "invalid"),
    info = "should error on invalid plot_type"
  )

  # --------------------------------------------------
  # Test error on invalid sampling_heterogeneity
  # --------------------------------------------------

  expect_error(funnel(fit_brma, sampling_heterogeneity = "yes"),
    info = "should error on invalid sampling_heterogeneity"
  )
})

# ============================================================================ #
# Test: Funnel Plot Customization
# ============================================================================ #

test_that("Funnel plot customization works", {

  name <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test custom point aesthetics
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_custom_points_base", function() {
    funnel(fit_brma, plot_type = "base", pch = 21, col = "blue", bg = "lightblue", cex = 1.5)
  })

  vdiffr::expect_doppelganger(
    "funnel_custom_points_ggplot",
    funnel(fit_brma, plot_type = "ggplot", pch = 19, col = "blue", bg = "lightblue", size = 3)
  )

  # --------------------------------------------------
  # Test custom funnel region styling
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_custom_regions_base", function() {
    funnel(fit_brma, plot_type = "base", back = "lightgrey", shade = "lightyellow", lty = "dashed")
  })

  vdiffr::expect_doppelganger(
    "funnel_custom_regions_ggplot",
    funnel(fit_brma, plot_type = "ggplot", back = "lightblue", shade = "lightyellow", lty = "dashed")
  )

  # --------------------------------------------------
  # Test suppressing background/shade
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_no_background_base", function() {
    funnel(fit_brma, plot_type = "base", back = NA, shade = "white")
  })

  vdiffr::expect_doppelganger(
    "funnel_no_shade_ggplot",
    funnel(fit_brma, plot_type = "ggplot", back = "grey", shade = NA)
  )

  # --------------------------------------------------
  # Test custom axis labels and title
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_custom_labels_base", function() {
    funnel(fit_brma, plot_type = "base", xlab = "Effect Size Residual", ylab = "SE", main = "Funnel Plot")
  })

  vdiffr::expect_doppelganger(
    "funnel_custom_labels_ggplot",
    funnel(fit_brma, plot_type = "ggplot", xlab = "Effect Size Residual", ylab = "SE", main = "Funnel Plot")
  )

  # --------------------------------------------------
  # Test line color customization
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_custom_lines_base", function() {
    funnel(fit_brma, plot_type = "base", col.line = "darkgrey", col.refline = "red", lty = "solid")
  })

  vdiffr::expect_doppelganger(
    "funnel_custom_lines_ggplot",
    funnel(fit_brma, plot_type = "ggplot", col.line = "darkgrey", col.refline = "red", lty = "solid")
  )
})

