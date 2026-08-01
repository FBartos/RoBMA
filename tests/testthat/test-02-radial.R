context("Radial (Galbraith) plot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))

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

  skip_if_not_full_visuals(
    "Default mathematical-label geometry needs maintainer review on the current graphics toolchain."
  )

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("radial_simple_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 5))
    metafor::radial(fit_metafor, main = "metafor")
    radial(fit_brma, plot_type = "base", main = "brma")
  })

  expect_vdiffr_snapshot("radial_simple_centered_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 5))
    metafor::radial(fit_metafor, center = TRUE, main = "metafor")
    radial(fit_brma, center = TRUE, plot_type = "base", main = "brma")
  })

  expect_vdiffr_snapshot(
    "radial_simple_brma_ggplot",
    radial(fit_brma, plot_type = "ggplot")
  )

  expect_vdiffr_snapshot(
    "radial_simple_centered_brma_ggplot",
    radial(fit_brma, center = TRUE, plot_type = "ggplot")
  )
})

test_that("Radial plot retains representative visual coverage", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("radial_custom_labels_base", function() {
    radial(
      fit_brma,
      plot_type = "base",
      xlab      = "Precision",
      zlab      = "z-score",
      main      = "Radial Plot"
    )
  })

  expect_vdiffr_snapshot(
    "radial_custom_labels_ggplot",
    radial(
      fit_brma,
      plot_type = "ggplot",
      xlab      = "Precision",
      zlab      = "z-score",
      main      = "Radial Plot"
    )
  )
})

# ============================================================================ #
# Test: 3-Level Model Radial Plot
# ============================================================================ #

test_that("Radial plot rejects 3-level models", {

  name     <- "konstantopoulos2011_3lvl"
  fit_brma <- fits[[name]]

  expect_error(
    radial(fit_brma, plot_type = "base"),
    "multilevel",
    info = "3-level models are rejected"
  )
})

# ============================================================================ #
# Test: GLMM Model Radial Plot
# ============================================================================ #

test_that("Radial plot rejects GLMM models", {

  name     <- "nielweise2008_glmm"
  fit_brma <- fits[[name]]

  expect_error(
    radial(fit_brma, plot_type = "base"),
    "normal outcome models",
    info = "GLMM models are rejected"
  )
})

# ============================================================================ #
# Test: Selection Model Radial Plot
# ============================================================================ #

test_that("Radial plot rejects selection models", {

  name     <- "dat.lehmann2018-3PSM"
  fit_brma <- fits[[name]]

  expect_error(
    radial(fit_brma, plot_type = "base"),
    "publication-bias",
    info = "selection models are rejected"
  )
})

# ============================================================================ #
# Test: PET Model Radial Plot
# ============================================================================ #

test_that("Radial plot rejects PET models", {

  name     <- "dat.lehmann2018-PET"
  fit_brma <- fits[[name]]

  expect_error(
    radial(fit_brma, plot_type = "base"),
    "publication-bias",
    info = "PET models are rejected"
  )
})

# ============================================================================ #
# Test: BMA.norm Model Radial Plot
# ============================================================================ #

test_that("Radial plot rejects BMA.norm model-averaging objects", {

  name     <- "dat.lehmann2018_BMA.norm"
  fit_brma <- fits[[name]]

  expect_error(
    radial(fit_brma, plot_type = "base"),
    "model-averaging",
    info = "BMA.norm objects are rejected"
  )
})

# ============================================================================ #
# Test: BMA.glmm Model Radial Plot
# ============================================================================ #

test_that("Radial plot rejects BMA.glmm model-averaging objects", {

  name     <- "bcg_BMA.glmm"
  fit_brma <- fits[[name]]

  expect_error(
    radial(fit_brma, plot_type = "base"),
    "normal outcome models",
    info = "BMA.glmm objects are rejected"
  )
})

# ============================================================================ #
# Test: RoBMA Model Radial Plot
# ============================================================================ #

test_that("Radial plot rejects RoBMA model-averaging objects", {

  name     <- "dat.lehmann2018_RoBMA"
  fit_brma <- fits[[name]]

  expect_error(
    radial(fit_brma, plot_type = "base"),
    "model-averaging",
    info = "RoBMA objects are rejected"
  )
})

# ============================================================================ #
# Test: Error on Unsupported Models
# ============================================================================ #

test_that("Radial plot errors on unsupported model types", {

  # meta-regression (has moderators)
  expect_error(
    radial(fits[["bcg_meta-regression"]]),
    "moderators",
    info = "meta-regression model is rejected"
  )

  # location-scale model (has both mods and scale, moderator check fires first)
  expect_error(
    radial(fits[["bangertdrowns2004_location-scale"]]),
    "moderators",
    info = "location-scale model is rejected"
  )

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)
  expect_error(
    radial(fits[[name]]),
    "known-V",
    info = "known-V brma.mv model is rejected"
  )
})


# ============================================================================ #
# Test: Radial Plot Interface
# ============================================================================ #

test_that("Radial plot data and alias interface are stable", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_data = TRUE returns list with expected components
  # --------------------------------------------------

  radial_data <- radial(fit_brma, as_data = TRUE)

  expect_true(is.list(radial_data),
    info = "as_data = TRUE returns a list"
  )

  expected_components <- c(
    "points", "refline", "band", "arc", "arc_ticks",
    "ci_arc", "ci_ticks", "xlim", "zlim"
  )
  expect_true(all(expected_components %in% names(radial_data)),
    info = "radial data contains all expected components"
  )

  # Check points data.frame structure
  expect_true(is.data.frame(radial_data$points),
    info = "points are returned as a data.frame"
  )
  expect_true(all(c("x", "z") %in% names(radial_data$points)),
    info = "points contain x and z columns"
  )

  # Check number of points matches number of studies
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(radial_data$points), n_studies,
    info = "number of points matches number of studies"
  )

  # --------------------------------------------------
  # Test center changes data values
  # --------------------------------------------------

  radial_data_centered <- radial(fit_brma, center = TRUE, as_data = TRUE)

  expect_true(radial_data_centered$center,
    info = "center flag is TRUE when centering"
  )

  # centered z-values should differ from non-centered
  expect_false(
    all(abs(radial_data$points$z - radial_data_centered$points$z) < 1e-10),
    info = "centering changes z-axis values"
  )

  # --------------------------------------------------
  # Test error on invalid plot_type
  # --------------------------------------------------

  expect_error(radial(fit_brma, plot_type = "invalid"),
    info = "invalid plot_type is rejected"
  )

  expect_error(radial(fit_brma, level = 100),
    regexp = "level",
    info = "invalid confidence level is rejected"
  )
  expect_error(radial(fit_brma, xlim = c(1, 1)),
    regexp = "xlim",
    info = "non-increasing xlim is rejected"
  )
  expect_error(radial(fit_brma, arc.res = 1),
    regexp = "arc.res",
    info = "invalid arc resolution is rejected"
  )
  expect_error(radial(fit_brma, transf = 1),
    regexp = "transf",
    info = "invalid transformation is rejected"
  )

  # --------------------------------------------------
  # Test galbraith alias works
  # --------------------------------------------------

  galbraith_data <- galbraith(fit_brma, as_data = TRUE)
  expect_equal(radial_data, galbraith_data,
    info = "galbraith() produces same output as radial()"
  )
})


# ============================================================================ #
# Test: Radial Plot Customization
# ============================================================================ #

test_that("Radial plot customization snapshots are stable", {

  skip_if_not_full_visuals("Customization snapshots are visual-gallery coverage.")

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test custom point aesthetics
  # --------------------------------------------------

  expect_vdiffr_snapshot("radial_custom_points_base", function() {
    radial(fit_brma, plot_type = "base", pch = 21, col = "blue", bg = "lightblue", cex = 1.5)
  })

  expect_vdiffr_snapshot(
    "radial_custom_points_ggplot",
    radial(fit_brma, plot_type = "ggplot", pch = 21, col = "blue", bg = "lightblue", size = 3)
  )

  # --------------------------------------------------
  # Test suppressing background
  # --------------------------------------------------

  expect_vdiffr_snapshot("radial_no_background_base", function() {
    radial(fit_brma, plot_type = "base", back = NA)
  })

  expect_vdiffr_snapshot(
    "radial_no_background_ggplot",
    radial(fit_brma, plot_type = "ggplot", back = NA)
  )

  # --------------------------------------------------
  # Test custom level
  # --------------------------------------------------

  expect_vdiffr_snapshot("radial_custom_level_base", function() {
    radial(fit_brma, plot_type = "base", level = 99)
  })

  expect_vdiffr_snapshot(
    "radial_custom_level_ggplot",
    radial(fit_brma, plot_type = "ggplot", level = 99)
  )
})
