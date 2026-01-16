context("Funnel plot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip("UNDER UPDATE")
# list & load all fits
skip_if_no_fits()
skip_if_not_installed("metafor")
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


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

  vdiffr::expect_doppelganger("funnel_simple_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", addtau2 = TRUE, xlim = c(-3, 3), ylim = c(0, 0.8))
    funnel(fit_brma, plot_type = "base", xlim = c(-3, 3), ylim = c(0, 0.8), main = "brma")
  })

  vdiffr::expect_doppelganger("funnel_simple_brma_ggplot",
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

  vdiffr::expect_doppelganger("funnel_regression_comparison", function() {
    # the main difference is that metafor visualizes standardized residuals
    # while brma visualizes in how metafor incorporates heterogeneity into residuals
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor")
    funnel(fit_brma, plot_type = "base", sampling_heterogeneity = FALSE, main = "brma")
  })

  vdiffr::expect_doppelganger("funnel_regression_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: Location-Scale Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for location-scale model works correctly", {

  name        <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_scale_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 0.6))
    funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 0.6), sampling_heterogeneity = FALSE)
  })

  vdiffr::expect_doppelganger("funnel_scale_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: 3-Level Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for 3-level model works correctly", {

  name        <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_3lvl_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor")
    funnel(fit_brma, plot_type = "base")
    title(main = "brma")
  })

  vdiffr::expect_doppelganger("funnel_3lvl_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: GLMM Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for GLMM model works correctly", {

  name     <- "nielweise2008_glmm"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Residuals computation for GLMM uses same approach as residuals.brma
  # --------------------------------------------------

  # brma: get funnel plot data
  brma_funnel_data <- funnel(fit_brma, as_data = TRUE)
  brma_resid <- brma_funnel_data$points$x

  # brma: verify consistency with residuals function
  brma_resid_direct <- residuals(fit_brma)
  brma_resid_mean <- brma_resid_direct$summary[, "Mean"]

  expect_equal(brma_resid, brma_resid_mean, tolerance = 0.001,
               info = "funnel residuals should match residuals() output for GLMM")

  # --------------------------------------------------
  # Verify standard errors are computed from raw data
  # --------------------------------------------------

  # For GLMM, SE should be approximated from counts
  brma_se <- brma_funnel_data$points$y
  expect_true(all(brma_se > 0),
              info = "SE should be positive for all observations")
  expect_true(all(is.finite(brma_se)),
              info = "SE should be finite for all observations")

  # --------------------------------------------------
  # Visual comparison (GLMM only - no metafor funnel equivalent)
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_glmm_brma_base", function() {
    funnel(fit_brma, plot_type = "base")
  })

  vdiffr::expect_doppelganger("funnel_glmm_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: Selection Model Funnel Plot (Positive Effects)
# ============================================================================ #

test_that("Funnel plot for selection model (positive) works correctly", {

  name        <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_selmodel_pos_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor")
    funnel(fit_brma, plot_type = "base")
    title(main = "brma")
  })

  vdiffr::expect_doppelganger("funnel_selmodel_pos_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: Selection Model Funnel Plot (Negative Effects)
# ============================================================================ #

test_that("Funnel plot for selection model (negative) works correctly", {

  name        <- "dat.lehmann2018-3PSM_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_selmodel_neg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor")
    funnel(fit_brma, plot_type = "base")
    title(main = "brma")
  })

  vdiffr::expect_doppelganger("funnel_selmodel_neg_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: PET Model Funnel Plot (Positive Effects)
# ============================================================================ #

test_that("Funnel plot for PET model (positive) works correctly", {

  name        <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_PET_pos_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor")
    funnel(fit_brma, plot_type = "base")
    title(main = "brma")
  })

  vdiffr::expect_doppelganger("funnel_PET_pos_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: PET Model Funnel Plot (Negative Effects)
# ============================================================================ #

test_that("Funnel plot for PET model (negative) works correctly", {

  name        <- "dat.lehmann2018-PET_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_PET_neg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor")
    funnel(fit_brma, plot_type = "base")
    title(main = "brma")
  })

  vdiffr::expect_doppelganger("funnel_PET_neg_brma_ggplot",
    funnel(fit_brma, plot_type = "ggplot")
  )
})


# ============================================================================ #
# Test: Funnel Plot Options
# ============================================================================ #

test_that("Funnel plot options work correctly", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test sampling_heterogeneity = FALSE
  # --------------------------------------------------

  # with heterogeneity
  data_with_het <- funnel(fit_brma, sampling_heterogeneity = TRUE, as_data = TRUE)

  # without heterogeneity
  data_no_het <- funnel(fit_brma, sampling_heterogeneity = FALSE, as_data = TRUE)

  # funnel should be narrower without heterogeneity
  width_with_het <- max(data_with_het$funnel$x) - min(data_with_het$funnel$x)
  width_no_het   <- max(data_no_het$funnel$x)   - min(data_no_het$funnel$x)

  expect_true(width_no_het < width_with_het,
              info = "funnel should be narrower without heterogeneity")

  # --------------------------------------------------
  # Visual comparison: with vs without heterogeneity
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_het_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    funnel(fit_brma, sampling_heterogeneity = TRUE, plot_type = "base")
    title(main = "with heterogeneity")
    funnel(fit_brma, sampling_heterogeneity = FALSE, plot_type = "base")
    title(main = "without heterogeneity")
  })

  # --------------------------------------------------
  # Test plot_type works
  # --------------------------------------------------

  # base plot should return NULL invisibly
  result_base <- funnel(fit_brma, plot_type = "base")
  expect_null(result_base)

  # ggplot should return ggplot object
  result_gg <- funnel(fit_brma, plot_type = "ggplot")
  expect_s3_class(result_gg, "ggplot")
})


# ============================================================================ #
# Test: Funnel Plot Interface
# ============================================================================ #

test_that("Funnel plot has correct interface", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_data = TRUE returns list with expected components
  # --------------------------------------------------

  funnel_data <- funnel(fit_brma, as_data = TRUE)

  expect_true(is.list(funnel_data),
              info = "as_data = TRUE should return a list")

  expected_components <- c("points", "funnel", "funnel_edge1", "funnel_edge2",
                          "background", "x_range", "y_range")
  expect_true(all(expected_components %in% names(funnel_data)),
              info = "funnel data should contain all expected components")

  # Check points data.frame structure
  expect_true(is.data.frame(funnel_data$points),
              info = "points should be a data.frame")
  expect_true(all(c("x", "y") %in% names(funnel_data$points)),
              info = "points should have x and y columns")

  # Check number of points matches number of studies
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(funnel_data$points), n_studies,
               info = "number of points should match number of studies")

  # --------------------------------------------------
  # Test error on invalid plot_type
  # --------------------------------------------------

  expect_error(funnel(fit_brma, plot_type = "invalid"),
               info = "should error on invalid plot_type")

  # --------------------------------------------------
  # Test error on invalid sampling_heterogeneity
  # --------------------------------------------------

  expect_error(funnel(fit_brma, sampling_heterogeneity = "yes"),
               info = "should error on invalid sampling_heterogeneity")
})


# ============================================================================ #
# Test: Funnel Plot Customization
# ============================================================================ #

test_that("Funnel plot customization works", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test custom point aesthetics
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_custom_points_base", function() {
    funnel(fit_brma, plot_type = "base", pch = 19, col = "blue", bg = "lightblue", cex = 1.5)
  })

  vdiffr::expect_doppelganger("funnel_custom_points_ggplot",
    funnel(fit_brma, plot_type = "ggplot", pch = 19, col = "blue", bg = "lightblue", size = 3)
  )

  # --------------------------------------------------
  # Test custom funnel region styling
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_custom_regions_base", function() {
    funnel(fit_brma, plot_type = "base", back = "lightgrey", shade = "lightyellow", lty = "dashed")
  })

  vdiffr::expect_doppelganger("funnel_custom_regions_ggplot",
    funnel(fit_brma, plot_type = "ggplot", back = "lightgrey", shade = "lightyellow", lty = "dashed")
  )

  # --------------------------------------------------
  # Test suppressing background/shade
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_no_background_base", function() {
    funnel(fit_brma, plot_type = "base", back = NA, shade = "white")
  })

  vdiffr::expect_doppelganger("funnel_no_shade_ggplot",
    funnel(fit_brma, plot_type = "ggplot", back = "grey", shade = NA)
  )

  # --------------------------------------------------
  # Test custom axis labels and title
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_custom_labels_base", function() {
    funnel(fit_brma, plot_type = "base", xlab = "Effect Size Residual", ylab = "SE", main = "Funnel Plot")
  })

  vdiffr::expect_doppelganger("funnel_custom_labels_ggplot",
    funnel(fit_brma, plot_type = "ggplot", xlab = "Effect Size Residual", ylab = "SE", main = "Funnel Plot")
  )

  # --------------------------------------------------
  # Test custom axis ranges
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_custom_range_base", function() {
    funnel(fit_brma, plot_type = "base", xlim = c(-2, 2), ylim = c(0, 1))
  })

  vdiffr::expect_doppelganger("funnel_custom_range_ggplot",
    funnel(fit_brma, plot_type = "ggplot", xlim = c(-2, 2), ylim = c(0, 1))
  )

  # --------------------------------------------------
  # Test line color customization
  # --------------------------------------------------

  vdiffr::expect_doppelganger("funnel_custom_lines_base", function() {
    funnel(fit_brma, plot_type = "base", col.line = "darkgrey", col.refline = "red", lty = "solid")
  })

  vdiffr::expect_doppelganger("funnel_custom_lines_ggplot",
    funnel(fit_brma, plot_type = "ggplot", col.line = "darkgrey", col.refline = "red", lty = "solid")
  )
})


# ============================================================================ #
# Test: Standard Error Computation for Different Outcome Types
# ============================================================================ #

test_that("Standard errors are computed correctly for normal models", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # Get SE from funnel
  funnel_data <- funnel(fit_brma, as_data = TRUE)
  brma_se <- funnel_data$points$y

  # Get SE directly from data
  expected_se <- fit_brma$data$outcome$sei

  expect_equal(brma_se, expected_se, tolerance = 0.001,
               info = "SE for normal models should match sei from data")
})


test_that("Standard errors are computed correctly for binomial GLMM models", {

  name     <- "bcg_glmm"
  fit_brma <- fits[[name]]

  # Get SE from funnel
  funnel_data <- funnel(fit_brma, as_data = TRUE)
  brma_se <- funnel_data$points$y

  # Compute expected SE using metafor::escalc formula
  # SE(logOR) = sqrt(1/ai + 1/bi + 1/ci + 1/di) with zero-cell adjustment
  # For binomial GLMM, bi = n1i - ai, di = n2i - ci
  outcome_data <- fit_brma$data$outcome
  ai <- outcome_data$ai
  bi <- outcome_data$n1i - outcome_data$ai
  ci <- outcome_data$ci
  di <- outcome_data$n2i - outcome_data$ci

  # Apply zero-cell adjustment where needed
  needs_adj <- (ai == 0 | bi == 0 | ci == 0 | di == 0)
  ai_adj <- ai + 0.5 * needs_adj
  bi_adj <- bi + 0.5 * needs_adj
  ci_adj <- ci + 0.5 * needs_adj
  di_adj <- di + 0.5 * needs_adj

  expected_se <- sqrt(1/ai_adj + 1/bi_adj + 1/ci_adj + 1/di_adj)

  expect_equal(brma_se, expected_se, tolerance = 0.001,
               info = "SE for binomial GLMM should match escalc formula")

  # All SE should be positive and finite
  expect_true(all(brma_se > 0), info = "All SE should be positive")
  expect_true(all(is.finite(brma_se)), info = "All SE should be finite")
})


test_that("Standard errors are computed correctly for Poisson GLMM models", {

  name     <- "nielweise2008_glmm"
  fit_brma <- fits[[name]]

  # Get SE from funnel
  funnel_data <- funnel(fit_brma, as_data = TRUE)
  brma_se <- funnel_data$points$y

  # Compute expected SE using metafor::escalc formula
  # SE(logIRR) = sqrt(1/x1i + 1/x2i) with zero-cell adjustment
  outcome_data <- fit_brma$data$outcome
  x1i <- outcome_data$x1i
  x2i <- outcome_data$x2i

  # Apply zero-cell adjustment where needed
  needs_adj <- (x1i == 0 | x2i == 0)
  x1i_adj <- x1i + 0.5 * needs_adj
  x2i_adj <- x2i + 0.5 * needs_adj

  expected_se <- sqrt(1/x1i_adj + 1/x2i_adj)

  expect_equal(brma_se, expected_se, tolerance = 0.001,
               info = "SE for Poisson GLMM should match escalc formula")

  # All SE should be positive and finite
  expect_true(all(brma_se > 0), info = "All SE should be positive")
  expect_true(all(is.finite(brma_se)), info = "All SE should be finite")
})


# ============================================================================ #
# Test: Funnel Residuals Match metafor
# ============================================================================ #

test_that("Funnel residuals match metafor across model types", {

  # Define models and their tolerances
  models <- list(
    list(name = "bcg_meta-analysis",           tol = 0.05),
    list(name = "bcg_meta-regression",         tol = 0.10),
    list(name = "bangertdrowns2004_location-scale", tol = 0.05),
    list(name = "konstantopoulos2011_3lvl",    tol = 0.05),
    list(name = "dat.lehmann2018-3PSM",        tol = 0.05),
    list(name = "dat.lehmann2018-3PSM_neg",    tol = 0.05),
    list(name = "dat.lehmann2018-PET",         tol = 0.05),
    list(name = "dat.lehmann2018-PET_neg",     tol = 0.10)
  )

  for (m in models) {
    fit_metafor <- info[[m$name]][["metafor"]]
    fit_brma    <- fits[[m$name]]

    metafor_resid    <- residuals(fit_metafor)
    brma_funnel_data <- funnel(fit_brma, as_data = TRUE)
    brma_resid       <- brma_funnel_data$points$x

    expect_equal(brma_resid, as.vector(metafor_resid), tolerance = m$tol,
                 info = paste("brma funnel residuals should match metafor for", m$name))
  }
})


# ============================================================================ #
# Test: Funnel Standard Errors Match metafor
# ============================================================================ #

test_that("Funnel standard errors match metafor for normal models", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_se       <- sqrt(fit_metafor$vi)
  brma_funnel_data <- funnel(fit_brma, as_data = TRUE)
  brma_se          <- brma_funnel_data$points$y

  expect_equal(brma_se, metafor_se, tolerance = 0.001,
               info = "brma funnel SE should match metafor SE")
})
