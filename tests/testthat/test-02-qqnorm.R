context("Q-Q normal plot")

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
# Test: Simple Meta-Analysis Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot for simple meta-analysis matches metafor structure", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Visual comparison: side-by-side (rstandard for comparability)
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_simple_comparison_rstandard", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    qqnorm(fit_metafor, main = "metafor", type = "rstandard", ylim = c(-3, 3))
    qqnorm(fit_brma, plot_type = "base", main = "brma", type = "rstandard", ylim = c(-3, 3))
  })

  # --------------------------------------------------
  # Default (rstudent) - brma only
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_simple_rstudent_base", function() {
    qqnorm(fit_brma, plot_type = "base")
  })

  expect_vdiffr_snapshot(
    "qqnorm_simple_rstudent_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot")
  )
})

test_that("Q-Q plot supports known-V brma.mv residual targets", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  qq_student <- suppressWarnings(qqnorm(fit_brma, as_data = TRUE))
  qq_standard <- qqnorm(fit_brma, as_data = TRUE, type = "rstandard")

  expect_true(is.list(qq_student))
  expect_true(is.list(qq_standard))
  expect_equal(nrow(qq_student[["points"]]), nobs(fit_brma))
  expect_equal(nrow(qq_standard[["points"]]), nobs(fit_brma))
})

# ============================================================================ #
# Test: Meta-Regression Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot for meta-regression matches metafor residual quantiles", {

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Visual comparison: side-by-side (rstandard for comparability)
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_regression_comparison_rstandard", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    qqnorm(fit_metafor, main = "metafor", type = "rstandard", ylim = c(-3, 3))
    qqnorm(fit_brma, plot_type = "base", main = "brma", type = "rstandard", ylim = c(-3, 3))
  })

  expect_vdiffr_snapshot(
    "qqnorm_regression_rstudent_ggplot",
    suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  )
})

test_that("Q-Q plot for interaction meta-regression renders residual quantiles", {

  skip_if_not_full_visuals("Interaction Q-Q variants duplicate the core meta-regression visual.")

  name        <- "bcg_meta-regression4"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Visual comparison: side-by-side (rstandard for comparability)
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_regression4_comparison_rstandard", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    qqnorm(fit_metafor, main = "metafor", type = "rstandard", ylim = c(-3, 3))
    qqnorm(fit_brma, plot_type = "base", main = "brma", type = "rstandard", ylim = c(-3, 3))
  })

  expect_vdiffr_snapshot(
    "qqnorm_regression4_rstudent_ggplot",
    suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  )
})

# ============================================================================ #
# Test: Location-Scale Model Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot for location-scale model renders residual quantiles", {

  skip_if_not_full_visuals("Location-scale Q-Q variants are gallery coverage.")

  name     <- "bangertdrowns2004_location-scale"
  fit_brma <- fits[[name]]
  set.seed(1)

  expect_vdiffr_snapshot("qqnorm_scale_rstandard_base", function() {
    qqnorm(fit_brma, plot_type = "base", type = "rstandard")
  })

  expect_vdiffr_snapshot(
    "qqnorm_scale_rstudent_ggplot",
    suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  )
})

# ============================================================================ #
# Test: 3-Level Model Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot for 3-level model renders residual quantiles", {

  name     <- "konstantopoulos2011_3lvl"
  fit_brma <- fits[[name]]
  set.seed(1)

  expect_vdiffr_snapshot("qqnorm_3lvl_base", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })

  expect_vdiffr_snapshot(
    "qqnorm_3lvl_ggplot",
    suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  )
})

test_that("Q-Q plot for 3-level meta-regression renders residual quantiles", {

  skip_if_not_full_visuals("3-level meta-regression duplicates the default multilevel Q-Q visual.")

  name     <- "konstantopoulos2011_3lvl2"
  fit_brma <- fits[[name]]
  set.seed(1)

  expect_vdiffr_snapshot("qqnorm_3lvl2_base", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })

  expect_vdiffr_snapshot(
    "qqnorm_3lvl2_ggplot",
    suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  )
})

# ============================================================================ #
# Test: GLMM Model Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot rejects GLMMs without a discrete PIT convention", {

  name     <- "nielweise2008_glmm"
  fit_brma <- fits[[name]]
  expect_error(
    rstudent(fit_brma),
    "discrete PIT convention"
  )
  expect_error(
    qqnorm(fit_brma, envelope = FALSE, as_data = TRUE),
    "discrete PIT convention"
  )

  # rstandard should error for GLMM models
  expect_error(
    qqnorm(fit_brma, type = "rstandard"),
    info = "rstandard residuals are rejected for GLMM models"
  )
})

test_that("Q-Q plot rejects GLMM meta-regression without a PIT convention", {

  skip_if_not_full_visuals("GLMM meta-regression duplicates the default GLMM Q-Q visual.")

  name     <- "bcg_glmm_reg"
  fit_brma <- fits[[name]]
  expect_error(
    qqnorm(fit_brma),
    "discrete PIT convention"
  )

  # rstandard should error for GLMM models
  expect_error(
    qqnorm(fit_brma, type = "rstandard"),
    info = "rstandard residuals are rejected for GLMM models"
  )
})

# ============================================================================ #
# Test: PET Model Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot for PET model matches metafor residual quantiles", {

  name        <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Visual comparison: side-by-side (rstandard for comparability)
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_PET_base", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    qqnorm(fit_metafor, main = "metafor", type = "rstudent", ylim = c(-3, 3))
    suppressWarnings(qqnorm(fit_brma, plot_type = "base", main = "brma", type = "rstudent", ylim = c(-3, 3)))
  })
})

test_that("Q-Q plot for PET meta-regression matches metafor residual quantiles", {

  skip_if_not_full_visuals("PET meta-regression duplicates the default PET Q-Q visual.")

  name        <- "dat.lehmann2018-PETreg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Visual comparison: side-by-side (rstandard for comparability)
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_PETreg_base", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    qqnorm(fit_metafor, main = "metafor", type = "rstudent", ylim = c(-3, 3))
    suppressWarnings(qqnorm(fit_brma, plot_type = "base", main = "brma", type = "rstudent", ylim = c(-3, 3)))
  })
})

# ============================================================================ #
# Test: Selection Model Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot for selection model renders residual quantiles", {

  name     <- "dat.lehmann2018-3PSM"
  fit_brma <- fits[[name]]

  # rstandard should error for selection models
  expect_error(
    qqnorm(fit_brma, type = "rstandard"),
    info = "rstandard residuals are rejected for selection models"
  )
})

test_that("Selection Q-Q endpoint geometry is retained for certification", {

  skip_if_not_full_visuals(
    "Exact PIT-tail evaluation shifted one plotted point and needs maintainer review."
  )

  name     <- "dat.lehmann2018-3PSM"
  fit_brma <- fits[[name]]
  set.seed(1)

  expect_vdiffr_snapshot("qqnorm_selmodel_base", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })
})

test_that("Q-Q plot for selection meta-regression renders residual quantiles", {

  skip_if_not_full_visuals("Selection meta-regression duplicates the default selection Q-Q visual.")

  name     <- "dat.lehmann2018-3PSMreg"
  fit_brma <- fits[[name]]
  set.seed(1)

  expect_vdiffr_snapshot("qqnorm_selmodel_reg_base", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })

  # skip the ggplot version for to save test time
  # expect_vdiffr_snapshot(
  #   "qqnorm_selmodel_ggplot",
  #   suppressWarnings(qqnorm(fit_brma, plot_type = "ggplot"))
  # )

  # rstandard should error for selection models
  expect_error(
    qqnorm(fit_brma, type = "rstandard"),
    info = "rstandard residuals are rejected for selection models"
  )
})

# ============================================================================ #
# Test: BMA.norm Model Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot for BMA.norm model renders base output", {

  name     <- "dat.lehmann2018_BMA.norm"
  fit_brma <- fits[[name]]
  set.seed(1)

  expect_vdiffr_snapshot("qqnorm_BMA", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })
})

test_that("Q-Q plot for BMA.norm meta-regression renders base output", {

  skip_if_not_full_visuals("BMA meta-regression duplicates the default BMA Q-Q smoke test.")

  name     <- "dat.lehmann2018_BMA.norm_mods"
  fit_brma <- fits[[name]]
  set.seed(1)

  expect_vdiffr_snapshot("qqnorm_BMAreg", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })
})

# ============================================================================ #
# Test: BMA.glmm Model Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot rejects BMA.glmm without a discrete PIT convention", {

  name     <- "bcg_BMA.glmm_3lvl_location_scale"
  fit_brma <- fits[[name]]
  expect_error(qqnorm(fit_brma), "discrete PIT convention")
})

# ============================================================================ #
# Test: RoBMA Model Q-Q Plot
# ============================================================================ #

test_that("Q-Q plot for RoBMA model renders base output", {

  name     <- "dat.lehmann2018_RoBMA"
  fit_brma <- fits[[name]]
  set.seed(1)

  expect_vdiffr_snapshot("qqnorm_RoBMA", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base"))
  })
})

test_that("Q-Q plot for RoBMA meta-regression renders LOO-PIT output", {

  skip_if_not_full_visuals("RoBMA meta-regression duplicates the default RoBMA Q-Q smoke test.")

  name     <- "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  fit_brma <- fits[[name]]
  set.seed(1)

  expect_vdiffr_snapshot("qqnorm_RoBMA_complex", function() {
    suppressWarnings(qqnorm(fit_brma, plot_type = "base", type = "LOO-PIT"))
  })
})


# ============================================================================ #
# Test: Q-Q Plot Interface
# ============================================================================ #

test_that("Q-Q plot data and argument validation are stable", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Test as_data = TRUE returns list with expected components
  # --------------------------------------------------

  qq_data <- qqnorm(fit_brma, as_data = TRUE, type = "rstandard")

  expect_true(is.list(qq_data),
    info = "as_data = TRUE returns a list"
  )

  expected_components <- c(
    "x", "y", "points", "refline", "envelope",
    "xlim", "ylim", "xlab", "ylab"
  )
  expect_true(all(expected_components %in% names(qq_data)),
    info = "Q-Q data contains all expected components"
  )

  # Check points data.frame structure
  expect_true(is.data.frame(qq_data$points),
    info = "points are returned as a data.frame"
  )
  expect_true(all(c("x", "y") %in% names(qq_data$points)),
    info = "points contain x and y columns"
  )

  # Check number of points matches number of observations
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(qq_data$points), n_studies,
    info = "number of points matches number of observations"
  )

  # Check x values are sorted (theoretical quantiles)
  expect_equal(qq_data$x, sort(qq_data$x),
    info = "theoretical quantiles are sorted"
  )

  # Check y values are sorted (sorted residuals)
  expect_equal(qq_data$y, sort(qq_data$y),
    info = "sample quantiles are sorted"
  )

  # --------------------------------------------------
  # Test envelope = FALSE suppresses envelope
  # --------------------------------------------------

  qq_data_no_env <- qqnorm(fit_brma, as_data = TRUE, envelope = FALSE, type = "rstandard")
  expect_null(qq_data_no_env$envelope,
    info = "envelope is NULL when envelope = FALSE"
  )

  # --------------------------------------------------
  # Test envelope = TRUE includes envelope
  # --------------------------------------------------

  expect_true(is.data.frame(qq_data$envelope),
    info = "envelope is a data.frame when envelope = TRUE"
  )
  expect_true(all(c("x", "lower", "upper") %in% names(qq_data$envelope)),
    info = "envelope contains x, lower, upper columns"
  )

  set.seed(1)
  qq_env_1 <- qqnorm(
    fit_brma,
    as_data  = TRUE,
    type     = "rstandard",
    reps     = 10,
    envelope = TRUE
  )[["envelope"]]
  set.seed(999)
  qq_env_2 <- qqnorm(
    fit_brma,
    as_data  = TRUE,
    type     = "rstandard",
    reps     = 100,
    envelope = TRUE
  )[["envelope"]]
  expect_equal(qq_env_1, qq_env_2,
    info = "closed-form QQ envelope is deterministic and ignores reps"
  )

  # --------------------------------------------------
  # Test error on invalid plot_type
  # --------------------------------------------------

  expect_error(qqnorm(fit_brma, plot_type = "invalid"),
    info = "invalid plot_type is rejected"
  )

  # --------------------------------------------------
  # Test error on invalid type
  # --------------------------------------------------

  expect_error(qqnorm(fit_brma, type = "invalid"),
    info = "invalid type is rejected"
  )

  expect_error(
    qqnorm(fit_brma, type = "rstandard", max_samples = 9),
    "'max_samples' must be at least 10",
    fixed = TRUE,
    info  = "invalid maximum posterior sample count is rejected"
  )
})

# ============================================================================ #
# Test: Q-Q Plot Customization
# ============================================================================ #

test_that("Q-Q plot customization snapshots are stable", {

  skip_if_not_full_visuals("Customization snapshots are visual-gallery coverage.")

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]
  set.seed(1)

  # --------------------------------------------------
  # Test custom point aesthetics
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_custom_points_base", function() {
    qqnorm(fit_brma, plot_type = "base", pch = 21, col = "blue",
           bg = "lightblue", cex = 1.5, type = "rstandard")
  })

  expect_vdiffr_snapshot(
    "qqnorm_custom_points_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", pch = 21, col = "blue",
           bg = "lightblue", size = 3, type = "rstandard")
  )

  # --------------------------------------------------
  # Test custom axis labels and title
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_custom_labels_base", function() {
    qqnorm(fit_brma, plot_type = "base", xlab = "Expected", ylab = "Observed",
           main = "QQ Plot", type = "rstandard")
  })

  expect_vdiffr_snapshot(
    "qqnorm_custom_labels_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", xlab = "Expected", ylab = "Observed",
           main = "QQ Plot", type = "rstandard")
  )

  # --------------------------------------------------
  # Test without envelope
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_no_envelope_base", function() {
    qqnorm(fit_brma, plot_type = "base", envelope = FALSE, type = "rstandard")
  })

  expect_vdiffr_snapshot(
    "qqnorm_no_envelope_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", envelope = FALSE, type = "rstandard")
  )

  # --------------------------------------------------
  # Test custom envelope level
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_custom_level_base", function() {
    qqnorm(fit_brma, plot_type = "base", level = 99, type = "rstandard")
  })

  expect_vdiffr_snapshot(
    "qqnorm_custom_level_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", level = 99, type = "rstandard")
  )

  # --------------------------------------------------
  # Test Bonferroni correction
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_bonferroni_base", function() {
    qqnorm(fit_brma, plot_type = "base", bonferroni = TRUE, type = "rstandard")
  })

  expect_vdiffr_snapshot(
    "qqnorm_bonferroni_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", bonferroni = TRUE, type = "rstandard")
  )

  # --------------------------------------------------
  # Test suppressing envelope shading
  # --------------------------------------------------

  expect_vdiffr_snapshot("qqnorm_no_shade_base", function() {
    qqnorm(fit_brma, plot_type = "base", shade = NA, type = "rstandard")
  })

  expect_vdiffr_snapshot(
    "qqnorm_no_shade_ggplot",
    qqnorm(fit_brma, plot_type = "ggplot", shade = NA, type = "rstandard")
  )
})
