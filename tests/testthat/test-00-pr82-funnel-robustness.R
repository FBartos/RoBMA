test_that("common funnel heterogeneity tolerates only double roundoff", {

  tau <- matrix(c(1, 2, 1 + .Machine$double.eps, 2), nrow = 2L)
  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) matrix(0, 2L, 1L),
    .funnel_row_heterogeneity_samples = function(...) tau,
    .funnel_selection_source_variances = function(...) NULL,
    .package = "RoBMA"
  )
  expect_true(.funnel_common_heterogeneity(list())$common)
  tau[1L, 2L] <- 1 + 1e-8
  expect_false(.funnel_common_heterogeneity(list())$common)
  tau[1L, 2L] <- NA_real_
  expect_error(.funnel_common_heterogeneity(list()),
               "finite row-marginal heterogeneity", fixed = TRUE)
})

test_that("conditioned funnel quadrature does not silently cap unresolved widths", {

  condition <- tryCatch(.funnel_conditioned_order(c(1, 0), c(0.1, 1)), error = identity)
  expect_identical(conditionMessage(condition), paste0(
    "Selected funnel contours were rejected by diagnostics: the conditioned-to-integrated ",
    "standard deviation ratio is not finite. Use 'sampling_bias = FALSE' for an ",
    "unselected reference."
  ))
  expect_identical(condition$maximum, 401L)
  expect_lt(.funnel_conditioned_order(c(1, 0), c(0.1, 0)), 401L)
})

test_that("funnel clipping handles singleton edges and rejects invalid coordinates", {

  expect_equal(.clip_line_x(2, 1, c(0, 1)), data.frame(x = numeric(), y = numeric()))
  expect_equal(.clip_line_x(0.5, 1, c(0, 1)), data.frame(x = 0.5, y = 1))
  expect_equal(.clip_line_x(c(-1, 2), c(0, 3), c(0, 1)),
               data.frame(x = c(0, 1), y = c(1, 2)))
  expect_error(.clip_line_x(c(0, NA), c(0, 1), c(0, 1)),
               "Funnel contour coordinates must be finite.", fixed = TRUE)
})

test_that("funnel data rejects missing standard errors before constructing axes", {

  testthat::local_mocked_bindings(
    .effect_direction = function(...) "positive",
    .outcome_type = function(...) "norm",
    .outcome_data_yi = function(...) c(0, 1),
    .outcome_data_sei = function(...) c(0.1, NA_real_),
    .package = "RoBMA"
  )
  expect_error(
    .funnel_data_outcome(list(), TRUE, TRUE, 10, "plugin", list(), list()),
    "Funnel-plot standard errors must be finite and non-negative.", fixed = TRUE
  )
})

test_that("funnel paths retain disjoint visible runs in SE order", {

  witness <- .clip_line_x(c(0, 2, 0), c(0, 1, 2), c(-1, 1))
  expect_equal(witness, data.frame(x = c(0, 1, NA, 1, 0),
    y = c(0, .5, NA, 1.5, 2)))
  disjoint <- .clip_line_x(c(-2, 2, -2), 0:2, c(-1, 1))
  expect_equal(disjoint, data.frame(x = c(-1, 1, NA, 1, -1),
    y = c(.25, .75, NA, 1.25, 1.75)))
  expect_equal(.clip_line_x(c(-2, 0, 2), 0:2, c(-1, 1)),
    data.frame(x = c(-1, 0, 1), y = c(.5, 1, 1.5)))
  expect_equal(.clip_line_x(rep(0, 3), 0:2, c(-1, 1)),
    data.frame(x = rep(0, 3), y = 0:2))
  expect_equal(.clip_line_x(c(2, 3, 2), 0:2, c(-1, 1)),
    data.frame(x = numeric(), y = numeric()))
  expect_identical(.clip_line_x(c(0, 2, 0), 0:2, c(1, -1)), witness)

  testthat::skip_if_not_installed("ggplot2")
  polygon <- data.frame(x = c(-1, -1, 1, 1), y = c(0, 2, 2, 0))
  data <- list(points = data.frame(x = 0, y = 1), funnel = polygon,
    background = polygon, funnel_edge1 = witness, funnel_edge2 = disjoint,
    refline = data.frame(x = c(.1, .1), y = c(0, 2)),
    xlim = c(-1, 1), ylim = c(0, 2), xlab = "Effect", ylab = "SE")
  dots <- .set_dots_funnel(list(shade = NA, back = NA))
  plot <- .funnel_plot_ggplot(data, dots)
  expect_s3_class(plot$layers[[2L]]$geom, "GeomPath")
  expect_identical(plot$layers[[2L]]$data, witness)
  expect_identical(plot$layers[[3L]]$data, disjoint)
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[2L]]$x, witness$x)
  expect_equal(built$data[[2L]]$y, -witness$y)
  expect_identical(is.na(built$data[[2L]]$x), is.na(witness$x))
  .with_temp_plot_device(expect_null(.funnel_plot_base(data, dots)))
})

test_that("funnel styling validates only the documented numeric controls", {

  for (argument in c("xlim", "ylim")) {
    for (value in list(NA_real_, c(0, Inf), c(0, 1, 2))) {
      expect_error(.set_dots_funnel(stats::setNames(list(value), argument)),
        paste0("'", argument, "'"), fixed = TRUE)
    }
  }
  expect_error(.set_dots_funnel(list(ylim = c(-.1, .2))), "'ylim'", fixed = TRUE)
  expect_error(.set_dots_funnel(list(refline = Inf)), "'refline'", fixed = TRUE)
  expect_error(.set_dots_funnel(list(refline = 1:2)), "'refline'", fixed = TRUE)
  expect_error(.set_dots_funnel(list(as_data = NA)), "'as_data'", fixed = TRUE)
  expect_error(.set_dots_funnel(list(as_data = NULL)), "as_data", fixed = TRUE)
  expect_identical(.set_dots_funnel(list(ylim = c(.3, 0)))$ylim, c(.3, 0))
  expect_identical(.set_dots_funnel(list(shade = "skyblue"))$shade, "skyblue")
  expect_true(is.na(.set_dots_funnel(list(shade = NA))$shade))
})

test_that("only automatic ineligible funnel commonness skips source splitting", {

  tau <- cbind(c(1, 1), c(2, 2))
  error <- structure(list(message = "Declared source failure.", call = NULL),
    class = c("review_source_error", "error", "condition"))
  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) matrix(0, 2L, 1L),
    .funnel_row_heterogeneity_samples = function(...) tau,
    .funnel_selection_source_variances = function(...) stop(error),
    .package = "RoBMA")
  automatic <- .funnel_common_heterogeneity(list(), commonness_only = TRUE)
  expect_false(automatic$common)
  expect_null(automatic$sources)
  expect_identical(tryCatch(.funnel_common_heterogeneity(list()), error = identity), error)
  tau[,] <- 1
  expect_identical(tryCatch(.funnel_common_heterogeneity(list(), commonness_only = TRUE),
    error = identity), error)
})

test_that("funnel source matrices reject malformed shapes before recycling", {

  samples <- matrix(seq_len(6L), 2L, 3L)
  testthat::local_mocked_bindings(
    .is_data_random = function(...) TRUE,
    .brma_mv_random_effects_marginal_vcov = function(...) list(samples = samples),
    .package = "RoBMA")
  evaluate <- function() .funnel_selection_source_row_variance(
    object = list(data = list()), posterior_samples = matrix(0, 2L, 1L),
    source = list(name = "declared"), S = 2L, K = 3L)
  expect_identical(evaluate(), samples)
  samples <- array(0, c(2L, 3L, 3L))
  samples[, 1L, 1L] <- 1
  samples[, 2L, 2L] <- 2
  samples[, 3L, 3L] <- 3
  expect_equal(evaluate(), matrix(rep(1:3, each = 2L), 2L))
  for (invalid in list(1:6, matrix(0, 3L, 2L), array(0, c(2L, 3L, 2L)))) {
    samples <- invalid
    expect_error(evaluate(), "Funnel random-source", fixed = TRUE)
  }
})

test_that("automatic public funnel routing avoids ineligible source splitting", {

  replay <- review_heterogeneity_har_replay(scaled = FALSE, direct = TRUE,
    n = 24L, weights = c(.2, .8), rho = 0, time = c(1, 2, 1, 2))
  BayesTools::JAGS_validate_fit_contract(replay$object$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  tau <- .funnel_row_heterogeneity_samples(replay$object, replay$samples)
  expect_identical(dim(tau), c(24L, 4L))
  expect_equal(as.numeric(tau), rep(2 * sqrt(c(.4, 1.6)), each = 24L, times = 2L),
    tolerance = 64 * .Machine$double.eps)
  error <- structure(list(message = "Declared source failure.", call = NULL),
    class = c("review_source_error", "error", "condition"))
  sentinel <- list(review_public_residual_route = TRUE)
  # Routing only: the covariance/commonness owner is real; residual density
  # production is replaced after admission and is not asserted numerically.
  testthat::local_mocked_bindings(
    .funnel_selection_source_variances = function(...) stop(error),
    .funnel_data_residual = function(...) sentinel, .package = "RoBMA")
  expect_identical(tryCatch(funnel.brma(replay$object, as_data = TRUE),
    error = identity), sentinel)
  expect_identical(tryCatch(funnel.brma(replay$object, residual = FALSE,
    as_data = TRUE), error = identity), error)
})
