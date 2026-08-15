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


test_that("LOO-PIT funnels project normalized residuals by predictive scale", {

  expected_z  <- c(1.25, -2)
  expected_se <- c(.4, .8)
  object      <- list()

  testthat::local_mocked_bindings(
    rstudent.brma = function(model, unit = "estimate", ...) {
      data.frame(
        resid = c(10, -10),
        se    = expected_se,
        z     = expected_z
      )
    },
    .package = "RoBMA"
  )

  loo_pit <- .funnel_data_residual(
    x                  = object,
    type               = "LOO-PIT",
    unit               = "estimate",
    conditioning_depth = "marginal",
    dots               = .set_dots_funnel(list())
  )
  rstudent_alias <- .funnel_data_residual(
    x                  = object,
    type               = "rstudent",
    unit               = "estimate",
    conditioning_depth = "marginal",
    dots               = .set_dots_funnel(list())
  )

  expect_equal(loo_pit[["points"]][["x"]], expected_z * expected_se)
  expect_equal(loo_pit[["points"]][["y"]], expected_se)
  expect_equal(
    loo_pit[["points"]][["x"]] / loo_pit[["points"]][["y"]],
    expected_z
  )
  expect_equal(rstudent_alias[["points"]], loo_pit[["points"]])
  expect_identical(loo_pit[["xlab"]], "LOO-PIT Residual")
})


test_that("known-V radial tau uses extra variance samples", {

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

  expect_equal(.get_radial_tau(object), expected)
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


test_that("radial tau uses within-draw RMS heterogeneity", {

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
    .radial_tau_samples(object, posterior_samples),
    sqrt(rowMeans(tau_total^2))
  )
})


test_that("outcome funnel eligibility requires row-invariant heterogeneity", {

  posterior_samples <- matrix(
    c(.1, .2, .3),
    ncol = 1L,
    dimnames = list(NULL, "mu")
  )
  object <- list(fit = list())

  testthat::local_mocked_bindings(
    .get_posterior_samples = function(fit, posterior_samples = NULL) {
      posterior_samples
    },
    .funnel_row_heterogeneity_samples = function(object, posterior_samples) {
      cbind(c(.2, .3, .4), c(.2, .3, .4))
    },
    .package = "RoBMA"
  )

  common <- .funnel_common_heterogeneity(object)
  expect_true(common[["common"]])
  expect_equal(common[["tau"]], c(.2, .3, .4))

  testthat::local_mocked_bindings(
    .funnel_row_heterogeneity_samples = function(object, posterior_samples) {
      cbind(c(.2, .3, .4), c(.2, .35, .4))
    },
    .package = "RoBMA"
  )

  varying <- .funnel_common_heterogeneity(object)
  expect_false(varying[["common"]])
  expect_null(varying[["tau"]])
})


test_that("random-formula funnel eligibility uses row-marginal ZGZ variance", {

  dat <- data.frame(
    yi    = c(.1, .2, .3, .4),
    study = c("s1", "s1", "s2", "s2"),
    x     = c(0, 1, 2, 4)
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(.04, 4)),
    data                      = dat,
    random                    = ~ diag(1 + x | study),
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- rbind(
    c(1, .5, .5),
    c(2, .8, .2)
  )
  colnames(posterior_samples) <- c(
    "mu__xRE_ALLOCx_allocation__total_sd",
    "mu__xRE_ALLOCx_allocation__weight[1]",
    "mu__xRE_ALLOCx_allocation__weight[2]"
  )

  heterogeneity <- .funnel_row_heterogeneity_samples(
    object            = object,
    posterior_samples = posterior_samples
  )

  expect_false(all(heterogeneity == heterogeneity[, 1L]))
  expect_false(.funnel_common_heterogeneity(
    object            = object,
    posterior_samples = posterior_samples
  )[["common"]])
})


test_that("plug-in funnels average complete joint-model CDFs", {

  posterior_samples <- cbind(
    mu    = c(0, 2, 10, 14),
    tau   = c(1, 3, 2, 4),
    model = c(1, 1, 2, 2)
  )
  common_heterogeneity <- list(
    posterior_samples = posterior_samples,
    tau                = posterior_samples[, "tau"]
  )

  testthat::local_mocked_bindings(
    .funnel_joint_model_groups = function(x, posterior_samples) {
      posterior_samples[, "model"]
    },
    .funnel_setup_from_samples = function(
        x, posterior_samples, tau_samples, sampling_heterogeneity,
        sampling_bias, weights) {
      list(
        posterior_samples = posterior_samples,
        tau                = tau_samples,
        weights            = weights
      )
    },
    .package = "RoBMA"
  )

  setup <- .funnel_sampling_setup(
    x                      = list(),
    sampling_heterogeneity = TRUE,
    sampling_bias          = TRUE,
    max_samples            = 10,
    estimand               = "plugin",
    common_heterogeneity   = common_heterogeneity
  )

  expect_equal(unname(setup[["posterior_samples"]][, "mu"]), c(1, 12))
  expect_equal(unname(setup[["tau"]]), c(2, 3))
  expect_equal(unname(setup[["weights"]]), c(.5, .5))
})


test_that("weighted funnel CDF mixes models before inversion", {

  setup <- list(
    mu                = c(0, 2),
    tau               = c(1, 1),
    PET               = c(0, 0),
    PEESE             = c(0, 0),
    is_weightfunction = c(FALSE, FALSE),
    selection         = NULL,
    weights           = c(.75, .25)
  )

  actual <- .funnel_model_averaged_cdf(
    q                = .5,
    se               = .3,
    setup            = setup,
    effect_direction = "positive"
  )
  total_sd <- sqrt(1 + .3^2)
  expected <- .75 * stats::pnorm(.5, 0, total_sd) +
    .25 * stats::pnorm(.5, 2, total_sd)

  expect_equal(actual, expected)
})


test_that("posterior-predictive contours integrate continuous uncertainty", {

  plugin <- list(
    mu                = 0,
    tau               = .5,
    PET               = 0,
    PEESE             = 0,
    is_weightfunction = FALSE,
    selection         = NULL,
    weights           = 1
  )
  posterior <- list(
    mu                = c(-1, 1),
    tau               = c(.5, .5),
    PET               = c(0, 0),
    PEESE             = c(0, 0),
    is_weightfunction = c(FALSE, FALSE),
    selection         = NULL,
    weights           = NULL
  )

  plugin_quantiles <- .get_funnel_quantiles_from_setup(
    se_sequence      = .2,
    setup            = plugin,
    effect_direction = "positive"
  )
  posterior_quantiles <- .get_funnel_quantiles_from_setup(
    se_sequence      = .2,
    setup            = posterior,
    effect_direction = "positive"
  )

  expect_lt(posterior_quantiles[["lower"]], plugin_quantiles[["lower"]])
  expect_gt(posterior_quantiles[["upper"]], plugin_quantiles[["upper"]])
})


test_that("funnel routes study-specific heterogeneity to residual mode", {

  object <- structure(list(), class = "brma")
  testthat::local_mocked_bindings(
    .is_mods  = function(object) FALSE,
    .is_scale = function(object) FALSE,
    .funnel_common_heterogeneity = function(object, posterior_samples = NULL) {
      list(common = FALSE)
    },
    .check_unit_conditioning_depth = function(...) invisible(TRUE),
    .funnel_data_residual = function(...) list(mode = "residual"),
    .package = "RoBMA"
  )

  expect_identical(funnel.brma(object, as_data = TRUE)[["mode"]], "residual")
  expect_error(
    funnel.brma(object, residual = FALSE, as_data = TRUE),
    "common marginal heterogeneity"
  )
})


test_that("bfunnel rejects targets without a common normal outcome scale", {

  object <- structure(list(), class = "brma")
  testthat::local_mocked_bindings(
    .outcome_type = function(object) "norm",
    .is_mods      = function(object) FALSE,
    .is_scale     = function(object) FALSE,
    .funnel_common_heterogeneity = function(object, posterior_samples = NULL) {
      list(common = FALSE)
    },
    .package = "RoBMA"
  )

  expect_error(
    bfunnel.brma(object, as_data = TRUE),
    "common marginal heterogeneity"
  )

  testthat::local_mocked_bindings(
    .outcome_type = function(object) "bin",
    .package = "RoBMA"
  )
  expect_error(
    bfunnel.brma(object, as_data = TRUE),
    "only for normal outcome models"
  )
})


test_that("GLMM outcome funnels disclose their descriptive approximation", {

  data <- list(outcome = data.frame(ai = c(1, 2), ci = c(2, 1)))
  attr(data, "outcome_type") <- "bin"
  object <- structure(list(data = data), class = c("brma.glmm", "brma"))

  testthat::local_mocked_bindings(
    .effect_direction     = function(object) "positive",
    .outcome_data_yi      = function(object) c(-.1, .1),
    .outcome_data_sei     = function(object) c(.2, .3),
    .get_funnel_quantiles = function(x, se_sequence, ...) {
      list(
        lower = rep(-1, length(se_sequence)),
        upper = rep(1, length(se_sequence)),
        mid   = rep(.25, length(se_sequence))
      )
    },
    .package = "RoBMA"
  )

  expect_warning(
    funnel_data <- .funnel_data_outcome(
      x                      = object,
      sampling_heterogeneity = TRUE,
      sampling_bias          = TRUE,
      max_samples            = 1000,
      estimand               = "plugin",
      common_heterogeneity   = list(),
      dots                   = .set_dots_funnel(list())
    ),
    "descriptive normal effect-size approximation"
  )
  expect_identical(
    funnel_data[["xlab"]],
    "Observed Effect Size (Descriptive Normal Approximation)"
  )
  expect_identical(
    funnel_data[["approximation"]],
    list(
      method                              = "normal_effect_size_approximation",
      descriptive                         = TRUE,
      fitted_discrete_likelihood_coverage = FALSE
    )
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
    suppressWarnings(.test_funnel(fit_brma, plot_type = "ggplot"))
  )
})

test_that("GLMM meta-regression exposes only a descriptive raw residual funnel", {

  name <- "bcg_glmm_reg"
  skip_if_missing_fits(name)

  fit_brma          <- fits[[name]]
  expected_residual <- residuals(fit_brma, type = "outcome")
  expected_se       <- .outcome_data_sei(fit_brma)

  expect_error(
    .test_funnel(fit_brma, plot_type = "ggplot"),
    "discrete PIT convention"
  )

  funnel_data <- .test_funnel(
    fit_brma,
    plot_type = "ggplot",
    type      = "outcome",
    as_data   = TRUE,
    xlim      = c(-100, 100),
    ylim      = c(0, max(expected_se))
  )

  expect_equal(
    unname(funnel_data[["points"]][["x"]]),
    unname(expected_residual)
  )
  expect_equal(funnel_data[["points"]][["y"]], expected_se)
  expect_equal(
    funnel_data[["funnel_edge1"]][["x"]],
    stats::qnorm(0.025) * funnel_data[["funnel_edge1"]][["y"]]
  )
  expect_equal(
    funnel_data[["funnel_edge2"]][["x"]],
    stats::qnorm(0.975) * funnel_data[["funnel_edge2"]][["y"]]
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

test_that("Outcome funnel rejects selection meta-regression", {

  skip_if_missing_fits("dat.lehmann2018-3PSMreg")

  fit_brma <- fits[["dat.lehmann2018-3PSMreg"]]
  expect_error(
    .test_funnel(
      fit_brma,
      residual               = FALSE,
      sampling_bias          = TRUE,
      sampling_heterogeneity = TRUE,
      max_samples            = 20,
      as_data                = TRUE
    ),
    "not supported.*location or scale predictors"
  )
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

test_that("BMA.glmm funnel respects the discrete PIT policy", {

  name     <- "bcg_BMA.glmm_3lvl_location_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_error(
    .test_funnel(fit_brma, plot_type = "base"),
    "discrete PIT convention"
  )

  expect_vdiffr_snapshot("funnel_BMA.glmm", function() {
    suppressWarnings(.test_funnel(
      fit_brma,
      plot_type = "base",
      type      = "outcome"
    ))
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
  bfunnel_data <- bfunnel(fit_brma, as_data = TRUE, max_samples = 100)

  expect_true(is.list(funnel_data),
    info = "as_data = TRUE returns a list"
  )
  expect_identical(funnel_data[["estimand"]], "plugin")
  expect_identical(bfunnel_data[["estimand"]], "posterior_predictive")

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


