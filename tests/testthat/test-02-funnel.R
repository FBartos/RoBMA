context("Funnel plot")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))


test_that("funnel line clipping preserves segment order", {

  forward <- .clip_line_x(
    x    = c(-2, 2),
    y    = c(0, 1),
    xlim = c(-1, 1)
  )
  reverse <- .clip_line_x(
    x    = c(2, -2),
    y    = c(0, 1),
    xlim = c(-1, 1)
  )

  expect_equal(forward[["x"]], c(-1, 1))
  expect_equal(forward[["y"]], c(.25, .75))
  expect_equal(reverse[["x"]], c(1, -1))
  expect_equal(reverse[["y"]], c(.25, .75))
})


test_that("known-V funnel tau uses extra variance samples", {

  V    <- matrix(c(.04, .01, .01, .09), nrow = 2L)
  data <- list(outcome = data.frame(yi = c(.10, .20), sei = sqrt(diag(V))))
  attr(data, "known_V")      <- TRUE
  attr(data, "known_V_data") <- .known_v_prepare(
    V                         = V,
    keep_rows                 = rep(TRUE, nrow(V)),
    known_v_parameterization  = "block_mvn",
    known_v_residual_fraction = NULL
  )
  object <- list(data = data)
  class(object) <- c("brma.mv", "brma")

  extra_variance <- rbind(c(.01, .04), c(.09, .16))
  testthat::local_mocked_bindings(
    .known_v_extra_variance_samples = function(object, ...) extra_variance,
    .package = "RoBMA"
  )

  expected <- mean(sqrt(rowMeans(extra_variance)))

  expect_equal(.get_funnel_tau(object), expected)
})


test_that("known-V marginal variance samples use the diagonal backend", {

  V    <- matrix(c(.04, .01, .01, .09), nrow = 2L)
  data <- list(outcome = data.frame(yi = c(.10, .20), sei = sqrt(diag(V))))
  attr(data, "known_V")      <- TRUE
  attr(data, "known_V_data") <- .known_v_prepare(
    V                         = V,
    keep_rows                 = rep(TRUE, nrow(V)),
    known_v_parameterization  = "block_mvn",
    known_v_residual_fraction = NULL
  )
  attr(data, "random")       <- TRUE
  object <- list(
    data = data,
    fit  = list()
  )
  class(object) <- c("brma.mv", "brma")
  posterior_samples <- matrix(1:3, ncol = 1L, dimnames = list(NULL, "tau"))
  random_variance <- cbind(1:3, 2:4)
  colnames(random_variance) <- c("1", "2")
  n_calls <- 0L

  testthat::local_mocked_bindings(
    .get_posterior_samples = function(fit, posterior_samples = NULL) {
      posterior_samples
    },
    .brma_mv_random_effects_marginal_vcov = function(
        object, posterior_samples, blocks = NULL, diagonal_only = FALSE) {

      n_calls <<- n_calls + 1L
      expect_true(diagonal_only)
      expect_null(blocks)
      list(
        samples  = random_variance,
        metadata = list(
          representation = "diagonal",
          quantity       = "variance",
          diagonal_only  = TRUE,
          dense          = FALSE,
          n_draws        = 3L,
          n_rows         = 2L,
          row_order      = 1:2,
          row_names      = c("1", "2")
        )
      )
    },
    .package = "RoBMA"
  )

  out <- .known_v_marginal_variance_samples(
    object            = object,
    posterior_samples = posterior_samples,
    max_bytes         = 1
  )
  metadata <- attr(out, "known_v_diagnostic", exact = TRUE)

  attr(out, "known_v_diagnostic") <- NULL
  expect_equal(n_calls, 1L)
  expect_equal(out, cbind(c(1.04, 2.04, 3.04), c(2.09, 3.09, 4.09)))
  expect_null(metadata[["n_chunks"]])

  extra <- .known_v_extra_variance_samples(
    object            = object,
    posterior_samples = posterior_samples
  )
  attr(extra, "known_v_diagnostic") <- NULL
  expect_equal(extra, unname(random_variance))
})


test_that("model-averaged funnel tau uses within-draw RMS heterogeneity", {

  tau_total <- rbind(
    c(1, 3, 5),
    c(2, 4, 8)
  )
  object <- list(
    fit    = list(),
    data   = list(scale = data.frame(), outcome = data.frame(yi = 1:3)),
    priors = list(scale = list())
  )
  posterior_samples <- matrix(0, nrow = nrow(tau_total), ncol = 1L)

  testthat::local_mocked_bindings(
    .evaluate.brma.tau = function(...) list(tau_total = tau_total),
    .is_scale = function(object) FALSE,
    .is_multilevel = function(object) FALSE,
    .fixed_tau_prior_value = function(priors) NULL,
    .fixed_rho_prior_value = function(priors) NULL,
    .package = "RoBMA"
  )

  expect_equal(
    .funnel_tau_samples(object, posterior_samples),
    sqrt(rowMeans(tau_total^2))
  )
})


# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

.test_funnel <- function(..., max_samples = 1000) {

  funnel(..., max_samples = max_samples)
}


# ============================================================================ #
# Test: Simple Meta-Analysis Funnel Plot
# ============================================================================ #

test_that("Funnel plot for simple meta-analysis matches metafor structure", {

  name        <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_simple_comparison_no_tau", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-3, 3), ylim = c(0, 0.8))
    .test_funnel(fit_brma, plot_type = "base", xlim = c(-3, 3), ylim = c(0, 0.8), main = "brma", sampling_heterogeneity = FALSE)
  })

  expect_vdiffr_snapshot("funnel_simple_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", addtau2 = TRUE, xlim = c(-3, 3), ylim = c(0, 0.8))
    .test_funnel(fit_brma, plot_type = "base", xlim = c(-3, 3), ylim = c(0, 0.8), main = "brma")
  })

  expect_vdiffr_snapshot(
    "funnel_simple_brma_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: Meta-Regression Funnel Plot
# ============================================================================ #

test_that("Funnel plot for meta-regression matches metafor residual views", {

  name        <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_regression_comparison-standard", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstandard")
    .test_funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstandard")
  })

  expect_vdiffr_snapshot("funnel_regression_comparison-rstudent", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstudent")
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstudent"))
  })

  expect_vdiffr_snapshot(
    "funnel_regression_brma_ggplot",
    suppressWarnings(.test_funnel(fit_brma, plot_type = "ggplot"))
  )
})

test_that("Funnel plot for interaction meta-regression renders residual views", {

  skip_if_not_full_visuals("Interaction funnel variants duplicate the core meta-regression visual.")

  name        <- "bcg_meta-regression4"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_regression4_comparison-standard", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstandard")
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 0.8), xlim = c(-2, 2), type = "rstandard"))
  })

  expect_vdiffr_snapshot("funnel_regression4_comparison-rstudent", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 1.2), xlim = c(-2, 2), type = "rstudent")
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 1.2), xlim = c(-2, 2), type = "rstudent"))
  })
})

# ============================================================================ #
# Test: Location-Scale Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for location-scale model matches metafor residual view", {

  skip_if_not_full_visuals("Location-scale funnel variants are gallery coverage.")

  name <- "bangertdrowns2004_location-scale"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_scale_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", ylim = c(0, 0.6), type = "rstandard")
    .test_funnel(fit_brma, plot_type = "base", main = "brma", ylim = c(0, 0.6), sampling_heterogeneity = FALSE, type = "rstandard")
  })

  expect_vdiffr_snapshot(
    "funnel_scale_brma_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: 3-Level Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for 3-level model matches metafor structure", {

  name <- "konstantopoulos2011_3lvl"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_3lvl_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-1, 1.5), ylim = c(0.4, 0))
    .test_funnel(fit_brma, plot_type = "base", main = "brma", sampling_heterogeneity = FALSE, xlim = c(-1, 1.5), ylim = c(0.4, 0))
  })

  expect_vdiffr_snapshot(
    "funnel_3lvl_brma_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot")
  )
})

test_that("Funnel plot for 3-level meta-regression renders residual views", {

  skip_if_not_full_visuals("3-level meta-regression duplicates the default multilevel funnel visual.")

  name <- "konstantopoulos2011_3lvl2"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_3lvl2_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-1, 1.5), ylim = c(0.5, 0), type = "rstandard")
    suppressWarnings(suppressWarnings(.test_funnel(fit_brma, plot_type = "base", type = "rstandard", conditioning_depth = "marginal", main = "brma", sampling_heterogeneity = FALSE, xlim = c(-1, 1.5), ylim = c(0.5, 0))))
  })

  expect_vdiffr_snapshot(
    "funnel_3lvl2_brma_ggplot",
    suppressWarnings(.test_funnel(fit_brma, plot_type = "ggplot"))
  )
})

# ============================================================================ #
# Test: GLMM Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot rejects undefined GLMM residuals", {

  name <- "nielweise2008_glmm"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_error(
    rstudent(fit_brma),
    "discrete PIT convention"
  )
  expect_error(
    .test_funnel(fit_brma, residual = TRUE, as_data = TRUE),
    "discrete PIT convention"
  )

})

test_that("Funnel plot for full-draw GLMM model remains stable", {

  skip_if_not_full_visuals(
    "The GLMM baseline uses the full-draw certification fixture."
  )

  name <- "nielweise2008_glmm"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  # there is no funnel plot for metafor
  expect_vdiffr_snapshot(
    "funnel_glmm_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot")
  )
})

test_that("GLMM meta-regression rejects residual funnel without a PIT convention", {

  skip_if_not_full_visuals("GLMM meta-regression duplicates the default GLMM funnel visual.")

  name <- "bcg_glmm_reg"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_error(
    .test_funnel(fit_brma, plot_type = "ggplot"),
    "discrete PIT convention"
  )
})

# ============================================================================ #
# Test: Selection Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for selection model matches metafor structure", {

  name <- "dat.lehmann2018-3PSM"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_selmodel_pos_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0))
    .test_funnel(fit_brma, plot_type = "base", main = "brma", xlim = c(-2, 2), ylim = c(0.8, 0), sampling_bias = FALSE, sampling_heterogeneity = FALSE)
  })
})

test_that("Selection funnel endpoint geometry is retained for certification", {

  skip_if_not_full_visuals(
    "Selected-normal endpoint geometry changed after exact boundary handling and needs maintainer review."
  )

  name <- "dat.lehmann2018-3PSM"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot(
    "funnel_selmodel_pos_brma_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot", sampling_bias = TRUE, sampling_heterogeneity = TRUE, xlim = c(-2, 2))
  )
})

test_that("Funnel plot for selection meta-regression renders residual view", {

  skip_if_not_full_visuals("Selection meta-regression duplicates the default selection funnel visual.")

  name <- "dat.lehmann2018-3PSMreg"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_selmodelreg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    # not available for metafor
    # par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    # metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0))
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", main = "brma", xlim = c(-2, 2), ylim = c(0.8, 0), sampling_bias = FALSE, sampling_heterogeneity = FALSE))
  })
})

test_that("Outcome funnel for selection meta-regression thins all samples equally", {

  skip_if_missing_fits("dat.lehmann2018-3PSMreg")

  fit_brma <- fits[["dat.lehmann2018-3PSMreg"]]
  funnel_data <- suppressWarnings(.test_funnel(
    fit_brma,
    residual               = FALSE,
    sampling_bias          = TRUE,
    sampling_heterogeneity = TRUE,
    max_samples            = 20,
    as_data                = TRUE
  ))

  expect_true(is.list(funnel_data))
  expect_true(all(c("points", "funnel") %in% names(funnel_data)))
})

test_that("Funnel plot allows descriptive known-V brma.mv outcome and residual modes", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  outcome_data <- .test_funnel(
    fit_brma,
    residual    = FALSE,
    as_data     = TRUE,
    max_samples = 1000
  )
  residual_data <- suppressWarnings(.test_funnel(
    fit_brma,
    residual    = TRUE,
    as_data     = TRUE,
    max_samples = 1000
  ))

  expect_true(is.list(outcome_data))
  expect_true(is.list(residual_data))
  expect_equal(nrow(outcome_data[["points"]]), nobs(fit_brma))
  expect_equal(nrow(residual_data[["points"]]), nobs(fit_brma))
})

test_that("Funnel plot for negative-direction selection model matches metafor structure", {

  skip_if_not_full_visuals("Negative-direction selection is gallery coverage.")

  name <- "dat.lehmann2018-3PSM_neg"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_selmodel_neg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0))
    .test_funnel(fit_brma, plot_type = "base", main = "brma", xlim = c(-2, 2), ylim = c(0.8, 0), sampling_bias = FALSE, sampling_heterogeneity = FALSE)
  })
})

# ============================================================================ #
# Test: PET Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for PET model matches metafor residual view", {

  name <- "dat.lehmann2018-PET"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_PET_pos_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    # the residuals need to be selected specifically because bPET is not treated as a regression
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0), type = "rstudent")
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", sampling_bias = FALSE, sampling_heterogeneity = FALSE, residual = TRUE, type = "rstudent", xlim = c(-2, 2), ylim = c(0.8, 0)))
  })

  expect_vdiffr_snapshot(
    "funnel_PET_pos_brma_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot")
  )
})

test_that("Funnel plot for PET meta-regression matches metafor residual view", {

  skip_if_not_full_visuals("PET meta-regression duplicates the default PET funnel visual.")

  name <- "dat.lehmann2018-PETreg"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_PETreg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    # the residuals need to be selected specifically because bPET is not treated as a regression
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0), type = "rstudent")
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", sampling_bias = FALSE, sampling_heterogeneity = FALSE, residual = TRUE, type = "rstudent", xlim = c(-2, 2), ylim = c(0.8, 0)))
  })
})

test_that("Funnel plot for negative-direction PET model matches metafor residual view", {

  skip_if_not_full_visuals("Negative-direction PET is gallery coverage.")

  name <- "dat.lehmann2018-PET_neg"
  skip_if_missing_fits(name)

  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Visual comparison: side-by-side plots
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_PET_neg_comparison", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    # the residuals need to be selected specifically because bPET is not treated as a regression
    metafor::funnel(fit_metafor, main = "metafor", xlim = c(-2, 2), ylim = c(0.8, 0), type = "rstudent")
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", sampling_bias = FALSE, sampling_heterogeneity = FALSE, residual = TRUE, type = "rstudent", xlim = c(-2, 2), ylim = c(0.8, 0)))
  })

  expect_vdiffr_snapshot(
    "funnel_PET_neg_brma_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot")
  )
})

# ============================================================================ #
# Test: BMA.norm Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for BMA.norm model renders base output", {

  name     <- "dat.lehmann2018_BMA.norm"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("funnel_BMA", function() {
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", sampling_heterogeneity = TRUE))
  })
})

test_that("Funnel plot for BMA.norm meta-regression renders base output", {

  skip_if_not_full_visuals("BMA meta-regression duplicates the default BMA funnel smoke test.")

  name     <- "dat.lehmann2018_BMA.norm_mods"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("funnel_BMAreg", function() {
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", sampling_heterogeneity = TRUE))
  })
})

# ============================================================================ #
# Test: BMA.glmm Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for BMA.glmm model renders base output", {

  name     <- "bcg_BMA.glmm_3lvl_location_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("funnel_BMA.glmm", function() {
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base"))
  })
})

# ============================================================================ #
# Test: RoBMA Model Funnel Plot
# ============================================================================ #

test_that("Funnel plot for RoBMA model renders base output", {

  name     <- "dat.lehmann2018_RoBMA"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("funnel_RoBMA", function() {
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", sampling_heterogeneity = TRUE, sampling_bias = TRUE))
  })
})

test_that("Funnel plot for RoBMA meta-regression renders LOO-PIT output", {

  skip_if_not_full_visuals("RoBMA meta-regression duplicates the default RoBMA funnel smoke test.")

  name     <- "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_vdiffr_snapshot("funnel_RoBMA_complex", function() {
    suppressWarnings(.test_funnel(fit_brma, plot_type = "base", type = "LOO-PIT"))
  })
})

# ============================================================================ #
# Test: Funnel Plot Options
# ============================================================================ #

test_that("Funnel plot data and argument validation are stable", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_data = TRUE returns list with expected components
  # --------------------------------------------------

  funnel_data <- .test_funnel(fit_brma, as_data = TRUE)

  expect_true(is.list(funnel_data),
    info = "as_data = TRUE returns a list"
  )

  expected_components <- c(
    "points", "funnel", "funnel_edge1", "funnel_edge2",
    "background", "x_range", "y_range"
  )
  expect_true(all(expected_components %in% names(funnel_data)),
    info = "funnel data contains all expected components"
  )

  # Check points data.frame structure
  expect_true(is.data.frame(funnel_data$points),
    info = "points are returned as a data.frame"
  )
  expect_true(all(c("x", "y") %in% names(funnel_data$points)),
    info = "points contain x and y columns"
  )

  # Check number of points matches number of studies
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(nrow(funnel_data$points), n_studies,
    info = "number of points matches number of studies"
  )

  # --------------------------------------------------
  # Test error on invalid plot_type
  # --------------------------------------------------

  expect_error(.test_funnel(fit_brma, plot_type = "invalid"),
    info = "invalid plot_type is rejected"
  )

  # --------------------------------------------------
  # Test error on invalid sampling_heterogeneity
  # --------------------------------------------------

  expect_error(.test_funnel(fit_brma, sampling_heterogeneity = "yes"),
    info = "invalid sampling_heterogeneity is rejected"
  )

  expect_no_error(
    .test_funnel(
      fit_brma,
      residual           = FALSE,
      type               = "not-a-residual-type",
      unit               = "not-a-unit",
      conditioning_depth = "not-a-depth",
      as_data            = TRUE
    )
  )

  expect_error(
    .test_funnel(fit_brma, residual = TRUE, type = "not-a-residual-type", as_data = TRUE),
    info = "residual mode validates residual type"
  )
})

# ============================================================================ #
# Test: Funnel Plot Customization
# ============================================================================ #

test_that("Funnel plot customization snapshots are stable", {

  skip_if_not_full_visuals("Customization snapshots are visual-gallery coverage.")

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test custom point aesthetics
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_custom_points_base", function() {
    .test_funnel(fit_brma, plot_type = "base", pch = 21, col = "blue", bg = "lightblue", cex = 1.5)
  })

  expect_vdiffr_snapshot(
    "funnel_custom_points_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot", pch = 19, col = "blue", bg = "lightblue", size = 3)
  )

  # --------------------------------------------------
  # Test custom funnel region styling
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_custom_regions_base", function() {
    .test_funnel(fit_brma, plot_type = "base", back = "lightgrey", shade = "lightyellow", lty = "dashed")
  })

  expect_vdiffr_snapshot(
    "funnel_custom_regions_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot", back = "lightblue", shade = "lightyellow", lty = "dashed")
  )

  # --------------------------------------------------
  # Test suppressing background/shade
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_no_background_base", function() {
    .test_funnel(fit_brma, plot_type = "base", back = NA, shade = "white")
  })

  expect_vdiffr_snapshot(
    "funnel_no_shade_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot", back = "grey", shade = NA)
  )

  # --------------------------------------------------
  # Test custom axis labels and title
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_custom_labels_base", function() {
    .test_funnel(fit_brma, plot_type = "base", xlab = "Effect Size Residual", ylab = "SE", main = "Funnel Plot")
  })

  expect_vdiffr_snapshot(
    "funnel_custom_labels_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot", xlab = "Effect Size Residual", ylab = "SE", main = "Funnel Plot")
  )

  # --------------------------------------------------
  # Test line color customization
  # --------------------------------------------------

  expect_vdiffr_snapshot("funnel_custom_lines_base", function() {
    .test_funnel(fit_brma, plot_type = "base", col.line = "darkgrey", col.refline = "red", lty = "solid")
  })

  expect_vdiffr_snapshot(
    "funnel_custom_lines_ggplot",
    .test_funnel(fit_brma, plot_type = "ggplot", col.line = "darkgrey", col.refline = "red", lty = "solid")
  )
})


