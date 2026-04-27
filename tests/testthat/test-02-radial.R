context("Radial (Galbraith) plot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)


# ============================================================================ #
# Test: Simple Meta-Analysis Radial Plot
# ============================================================================ #

test_that("Radial plot for simple meta-analysis matches metafor structure", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("radial_simple_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 5))
    metafor::radial(fit_metafor, main = "metafor")
    radial(fit_brma, plot_type = "base", main = "brma")
  })

  vdiffr::expect_doppelganger("radial_simple_centered_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 5))
    metafor::radial(fit_metafor, center = TRUE, main = "metafor")
    radial(fit_brma, center = TRUE, plot_type = "base", main = "brma")
  })

  vdiffr::expect_doppelganger(
    "radial_simple_brma_ggplot",
    radial(fit_brma, plot_type = "ggplot")
  )

  vdiffr::expect_doppelganger(
    "radial_simple_centered_brma_ggplot",
    radial(fit_brma, center = TRUE, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: 3-Level Model Radial Plot
# ============================================================================ #

test_that("Radial plot for 3-level model works correctly", {

  name     <- "konstantopoulos2011_3lvl"
  fit_brma <- fits[[name]]

  # metafor does not support radial for rma.mv, brma-only tests

  vdiffr::expect_doppelganger("radial_3lvl_brma_base", function() {
    radial(fit_brma, plot_type = "base")
  })

  vdiffr::expect_doppelganger(
    "radial_3lvl_brma_ggplot",
    radial(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: GLMM Model Radial Plot
# ============================================================================ #

test_that("Radial plot for GLMM model works correctly", {

  name     <- "nielweise2008_glmm"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("radial_glmm_brma_base", function() {
    radial(fit_brma, plot_type = "base")
  })

  vdiffr::expect_doppelganger(
    "radial_glmm_brma_ggplot",
    radial(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: Selection Model Radial Plot
# ============================================================================ #

test_that("Radial plot for selection model works correctly", {

  name     <- "dat.lehmann2018-3PSM"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("radial_selmodel_brma_base", function() {
    radial(fit_brma, plot_type = "base")
  })

  vdiffr::expect_doppelganger(
    "radial_selmodel_brma_ggplot",
    radial(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: PET Model Radial Plot
# ============================================================================ #

test_that("Radial plot for PET model works correctly", {

  name     <- "dat.lehmann2018-PET"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("radial_PET_brma_base", function() {
    radial(fit_brma, plot_type = "base")
  })

  vdiffr::expect_doppelganger(
    "radial_PET_brma_ggplot",
    radial(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: BMA.norm Model Radial Plot
# ============================================================================ #

test_that("Radial plot for BMA.norm model works correctly", {

  name     <- "dat.lehmann2018_BMA.norm"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("radial_BMA", function() {
    suppressWarnings(radial(fit_brma, plot_type = "base"))
  })
})

# ============================================================================ #
# Test: BMA.glmm Model Radial Plot
# ============================================================================ #

test_that("Radial plot for BMA.glmm model works correctly", {

  name     <- "bcg_BMA.glmm"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("radial_BMA.glmm", function() {
    suppressWarnings(radial(fit_brma, plot_type = "base"))
  })
})

# ============================================================================ #
# Test: RoBMA Model Radial Plot
# ============================================================================ #

test_that("Radial plot for RoBMA model works correctly", {

  name     <- "dat.lehmann2018_RoBMA"
  fit_brma <- fits[[name]]

  vdiffr::expect_doppelganger("radial_RoBMA", function() {
    suppressWarnings(radial(fit_brma, plot_type = "base"))
  })
})

# ============================================================================ #
# Test: Error on Unsupported Models
# ============================================================================ #

test_that("Radial plot errors on unsupported model types", {

  # meta-regression (has moderators)
  expect_error(
    radial(fits[["bcg_meta-regression"]]),
    "moderators",
    info = "should error on meta-regression model"
  )

  # location-scale model (has both mods and scale, moderator check fires first)
  expect_error(
    radial(fits[["bangertdrowns2004_location-scale"]]),
    "moderators",
    info = "should error on location-scale model"
  )
})


# ============================================================================ #
# Test: Radial Plot Interface
# ============================================================================ #

test_that("Radial plot has correct interface", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_data = TRUE returns list with expected components
  # --------------------------------------------------

  radial_data <- radial(fit_brma, as_data = TRUE)

  expect_true(is.list(radial_data),
    info = "as_data = TRUE should return a list"
  )

  expected_components <- c(
    "points", "refline", "band", "arc", "arc_ticks",
    "ci_arc", "ci_ticks", "xlim", "zlim"
  )
  expect_true(all(expected_components %in% names(radial_data)),
    info = "radial data should contain all expected components"
  )

  # Check points data.frame structure
  expect_true(is.data.frame(radial_data$points),
    info = "points should be a data.frame"
  )
  expect_true(all(c("x", "z") %in% names(radial_data$points)),
    info = "points should have x and z columns"
  )

  # Check number of points matches number of studies
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(radial_data$points), n_studies,
    info = "number of points should match number of studies"
  )

  # --------------------------------------------------
  # Test center changes data values
  # --------------------------------------------------

  radial_data_centered <- radial(fit_brma, center = TRUE, as_data = TRUE)

  expect_true(radial_data_centered$center,
    info = "center flag should be TRUE when centering"
  )

  # centered z-values should differ from non-centered
  expect_false(
    all(abs(radial_data$points$z - radial_data_centered$points$z) < 1e-10),
    info = "centering should change z-axis values"
  )

  # --------------------------------------------------
  # Test error on invalid plot_type
  # --------------------------------------------------

  expect_error(radial(fit_brma, plot_type = "invalid"),
    info = "should error on invalid plot_type"
  )

  # --------------------------------------------------
  # Test galbraith alias works
  # --------------------------------------------------

  galbraith_data <- galbraith(fit_brma, as_data = TRUE)
  expect_equal(radial_data, galbraith_data,
    info = "galbraith() should produce same output as radial()"
  )
})


# ============================================================================ #
# Test: Radial Plot Customization
# ============================================================================ #

test_that("Radial plot customization works", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test custom point aesthetics
  # --------------------------------------------------

  vdiffr::expect_doppelganger("radial_custom_points_base", function() {
    radial(fit_brma, plot_type = "base", pch = 21, col = "blue", bg = "lightblue", cex = 1.5)
  })

  vdiffr::expect_doppelganger(
    "radial_custom_points_ggplot",
    radial(fit_brma, plot_type = "ggplot", pch = 21, col = "blue", bg = "lightblue", size = 3)
  )

  # --------------------------------------------------
  # Test custom axis labels and title
  # --------------------------------------------------

  vdiffr::expect_doppelganger("radial_custom_labels_base", function() {
    radial(fit_brma, plot_type = "base", xlab = "Precision", zlab = "z-score", main = "Radial Plot")
  })

  vdiffr::expect_doppelganger(
    "radial_custom_labels_ggplot",
    radial(fit_brma, plot_type = "ggplot", xlab = "Precision", zlab = "z-score", main = "Radial Plot")
  )

  # --------------------------------------------------
  # Test suppressing background
  # --------------------------------------------------

  vdiffr::expect_doppelganger("radial_no_background_base", function() {
    radial(fit_brma, plot_type = "base", back = NA)
  })

  vdiffr::expect_doppelganger(
    "radial_no_background_ggplot",
    radial(fit_brma, plot_type = "ggplot", back = NA)
  )

  # --------------------------------------------------
  # Test custom level
  # --------------------------------------------------

  vdiffr::expect_doppelganger("radial_custom_level_base", function() {
    radial(fit_brma, plot_type = "base", level = 99)
  })

  vdiffr::expect_doppelganger(
    "radial_custom_level_ggplot",
    radial(fit_brma, plot_type = "ggplot", level = 99)
  )
})
