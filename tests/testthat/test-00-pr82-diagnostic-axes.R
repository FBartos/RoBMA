test_that("QQ limits contain one-sided residual ranges and reject missing residuals", {

  z <- c(0.5, 1, 2)
  testthat::local_mocked_bindings(
    residuals.brma = function(...) z,
    .package = "RoBMA"
  )
  get_data <- function() .qqnorm_data(
    x = list(), type = "LOO-PIT", unit = "estimate", conditioning_depth = "marginal",
    envelope = FALSE, conf_level = 95, bonferroni = FALSE, reps = 10,
    smooth = FALSE, max_samples = 10, xlim = NULL, ylim = NULL,
    xlab = NULL, ylab = NULL, dots = list()
  )
  out <- get_data()
  expect_lt(out$ylim[1L], min(z))
  expect_gt(out$ylim[2L], max(z))
  z <- -z
  out <- get_data()
  expect_lt(out$ylim[1L], min(z))
  expect_gt(out$ylim[2L], max(z))
  z[1L] <- NA_real_
  expect_error(get_data(), paste0("Normal QQ plots are unavailable because the requested ",
    "standardized residuals are non-finite. Inspect 'residuals()' for the requested target."), fixed = TRUE)
})

test_that("radial confidence limits survive scientific-notation probabilities", {

  values <- matrix(seq(-1, 1, length.out = 100), ncol = 1L,
                   dimnames = list(NULL, "mu"))
  testthat::local_mocked_bindings(
    pooled_effect = function(x, probs) .new_effect_brma_samples(
      samples = values, n_chains = 1L, n_iter = 100L, title = "test", probs = probs
    ),
    .outcome_data_yi = function(...) c(-0.2, 0.1, 0.4),
    .outcome_data_vi = function(...) c(0.04, 0.09, 0.16),
    .get_radial_tau_rows = function(...) rep(0.1, 3),
    .package = "RoBMA"
  )
  out <- .radial_data(
    x = list(), center = FALSE, xlim = NULL, zlim = NULL, xlab = NULL,
    zlab = NULL, atz = NULL, aty = NULL, steps = 5, level = 99.99,
    digits = 2, transf = NULL, targs = NULL, dots = .set_dots_radial(arc.res = 2)
  )
  expect_equal(out$ci_values[c(1L, 3L)],
               unname(stats::quantile(values, c(0.00005, 0.99995))))
  expect_gte(nrow(out$ci_arc), 2L)
})

test_that("QQ confidence levels and limits are checked before residual evaluation", {

  object <- brma(yi = c(.1, .2, .3), sei = rep(.2, 3), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  for (level in c(0, 100)) {
    expect_error(qqnorm.brma(object, conf_level = level),
      "'conf_level' must be greater than 0 and less than 100.", fixed = TRUE)
  }
  for (level in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(qqnorm.brma(object, conf_level = level),
      "'conf_level' must be a finite numeric vector.", fixed = TRUE)
  }
  for (argument in c("xlim", "ylim")) {
    args <- stats::setNames(list(c(0, Inf)), argument)
    expect_error(do.call(qqnorm.brma, c(list(y = object), args)),
      paste0("'", argument, "'"), fixed = TRUE)
    args <- stats::setNames(list(c(1, 0)), argument)
    expect_error(do.call(qqnorm.brma, c(list(y = object), args)),
      paste0("'", argument, "' must be an increasing numeric vector of length 2."), fixed = TRUE)
  }
  expect_error(qqnorm.brma(object, xlab = 1), "'xlab' must be a single character string")
  expect_error(qqnorm.brma(object, as_data = NA), "as_data")
  expect_error(qqnorm.brma(object, las = 5), "'las' must be one of 0, 1, 2, or 3.", fixed = TRUE)
  expect_error(qqnorm.brma(object, cex = -1), "'cex'")
  expect_type(.set_dots_qqnorm(list(cex = c(.5, 1, 1.5), size = 0)), "list")
})

test_that("row-varying radial coordinates use a declared random source", {

  replay <- review_heterogeneity_har_replay(scaled = FALSE, direct = TRUE,
    n = 24L, weights = c(.2, .8), rho = 0, time = c(1, 2, 1, 2))
  BayesTools::JAGS_validate_fit_contract(replay$object$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  tau <- rep(2 * sqrt(c(.4, 1.6)), 2L)
  actual <- radial.brma(replay$object, as_data = TRUE)
  expected_x <- 1 / sqrt(.04 + tau^2)
  expect_equal(actual$points$x, expected_x, tolerance = 1e-14)
  expect_equal(actual$points$z, seq(.1, .4, length.out = 4L) * expected_x,
    tolerance = 1e-14)
  expect_identical(actual$xlab, expression(x[i] == 1 / sqrt(v[i] + tau[i]^2)))
  scaled <- review_heterogeneity_har_replay(scaled = TRUE)
  expect_error(radial.brma(scaled$object, as_data = TRUE),
    "Radial plots cannot be drawn for models with scale regression.", fixed = TRUE)
})

test_that("radial labels follow declared designs even when row means coincide", {

  replay <- review_heterogeneity_har_replay(scaled = FALSE, direct = TRUE,
    n = 24L, weights = c(.5, .5), rho = 0, time = c(1, 2, 1, 2))
  object <- replay$object
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  tau <- .get_radial_tau_rows(object)
  expect_equal(tau, rep(tau[1L], length(tau)), tolerance = 0)
  actual <- radial.brma(object, as_data = TRUE)
  expect_identical(actual$xlab, expression(x[i] == 1 / sqrt(v[i] + tau[i]^2)))
  expect_equal(actual$points$x, 1 / sqrt(.04 + tau^2), tolerance = 0)

  setup <- review_compile_object(brma(yi = c(.1, .2, .3, .4), sei = rep(.2, 4L),
    measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("point", list(location = 0)),
    prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
    only_priors = TRUE, silent = TRUE))
  priors <- .create_fit_priors(setup$object$data, setup$object$priors)
  columns <- review_prior_columns(priors)
  expect_setequal(columns, c("mu", "tau"))
  samples <- matrix(0, 24L, length(columns), dimnames = list(NULL, columns))
  for (name in names(priors)) {
    stopifnot(BayesTools::is.prior.point(priors[[name]]))
    samples[, name] <- priors[[name]]$parameters$location
  }
  ordinary <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(ordinary$fit, requires = c(
    "name_encoding", "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  actual <- radial.brma(ordinary, as_data = TRUE)
  expect_identical(actual$xlab, expression(x[i] == 1 / sqrt(v[i] + tau^2)))
})

test_that("QQ plots omit exactly fitted estimates with undefined standardized residuals", {

  object <- brma(yi = c(.1, .2, .3), sei = rep(.2, 3), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  testthat::local_mocked_bindings(
    .check_unit_conditioning_depth = function(...) invisible(TRUE),
    rstandard.brma = function(...) data.frame(resid = c(.1, 0, -.2), se = c(.1, 0, .1),
                                              z = c(1, NA, -2)),
    .package = "RoBMA"
  )
  expect_message(
    data <- qqnorm.brma(object, type = "rstandard", envelope = FALSE, as_data = TRUE),
    "1 exactly fitted estimate(s) with an undefined standardized residual are omitted.",
    fixed = TRUE
  )
  expect_identical(data$K, 2L)
  expect_identical(data$y, c(-2, 1))

  testthat::local_mocked_bindings(
    rstandard.brma = function(...) data.frame(resid = 0, se = 0, z = NaN),
    .package = "RoBMA"
  )
  expect_error(qqnorm.brma(object, type = "rstandard", envelope = FALSE, as_data = TRUE),
               "non-finite")
})

test_that("QQ point styles follow their estimates through omission and sorting", {

  testthat::skip_if_not_installed("ggplot2")
  object <- brma(yi = c(.1, .2, .3, .4), sei = rep(.2, 4), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  testthat::local_mocked_bindings(
    .check_unit_conditioning_depth = function(...) invisible(TRUE),
    rstandard.brma = function(...) data.frame(resid = c(.1, 0, -.2, .05),
                                              se = c(.1, 0, .1, .1), z = c(1, NA, -2, .5)),
    .package = "RoBMA"
  )
  expect_message(
    plot <- qqnorm.brma(object, type = "rstandard", envelope = FALSE, plot_type = "ggplot",
                        size = c(1, 2, 3, 4), col = c("red", "green", "blue", "black")),
    "omitted"
  )
  built  <- ggplot2::ggplot_build(plot)
  points <- built$data[[which(vapply(plot$layers, function(layer) {
    inherits(layer$geom, "GeomPoint")
  }, logical(1)))[1L]]]
  # rows 3, 4, 1 in increasing z order; row 2 is omitted
  expect_identical(points$size, c(3, 4, 1))
  expect_identical(points$colour, c("blue", "black", "red"))
})
