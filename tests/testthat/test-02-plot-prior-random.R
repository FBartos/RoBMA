source(testthat::test_path("common-functions.R"))

# plot_prior() and print_prior() of the random-effect selectors of two fitted
# brma.mv models prepared by test-01-plot-prior-random.R: a total SD split by a
# Dirichlet(2, 3) allocation between two named blocks (the declared priors the
# expectations use), and a block whose SD is modelled on the scale. See
# test-00-plot-prior-random.R for the selectors of the other fits.

test_that("random SD, total and share selectors of an allocation show their priors", {

  skip_on_cran()
  fit <- load_fit("plot_prior_random_allocated")
  prior_list <- attr(fit[["fit"]], "prior_list")
  total   <- prior_list[["mu__xRE_ALLOCx_heterogeneity__allocation_sd"]]
  weights <- prior_list[["mu__xRE_ALLOCx_heterogeneity__weight"]]
  alpha   <- weights[["parameters"]][["alpha"]]
  expect_equal(alpha, c(2, 3))

  # a block-less SD is ambiguous between the blocks
  expect_error(plot_prior(fit, "tau"), class = "BayesTools_parameter_ambiguous")
  expect_error(print_prior(fit, "tau"), class = "BayesTools_parameter_ambiguous")

  # the total SD is the declared SD prior
  curve <- .random_prior_curve(plot_prior(fit, "tau_total", plot_type = "ggplot"))
  .expect_random_prior_curve(
    curve, exp(BayesTools::lpdf(total, curve[["x"]])), 1e-8
  )
  printed <- capture.output(returned <- print_prior(fit, "tau_total"))
  expect_identical(returned, total)
  expect_true(all(paste0("  ", print(total, silent = TRUE)) %in% printed))

  # the variance of the total is the square of the declared total SD
  curve <- .random_prior_curve(plot_prior(fit, "tau2_total", plot_type = "ggplot"))
  .expect_random_prior_curve(
    curve,
    exp(BayesTools::lpdf(total, sqrt(curve[["x"]]))) / (2 * sqrt(curve[["x"]])),
    1e-6
  )
  printed <- capture.output(print_prior(fit, "tau2_total"))
  expect_match(printed[[1L]], "tau2_total = tau_total^2", fixed = TRUE)
  expect_true(all(paste0("    ", print(total, silent = TRUE)) %in% printed))

  # the share of the total variance is the Beta marginal of the Dirichlet
  # weights, whose derivation names the declared weights prior
  for (component in 1:2) {
    selector <- c("tau2_prop(Study)", "tau2_prop(Crit)")[[component]]
    curve <- .random_prior_curve(plot_prior(fit, selector, plot_type = "ggplot"))
    expect_gt(nrow(curve), 500L)
    .expect_random_prior_curve(
      curve,
      stats::dbeta(curve[["x"]], alpha[[component]], sum(alpha) - alpha[[component]]),
      1e-6, info = selector
    )
    printed <- capture.output(print_prior(fit, selector))
    expect_match(
      printed[[1L]],
      paste0("w[", component, "], the share of the total variance of allocation 'heterogeneity'"),
      fixed = TRUE
    )
    expect_true(all(paste0("    ", print(weights, silent = TRUE)) %in% printed))
  }

  # the SD of a component is the total SD times sqrt(w[k]): the plotted law
  # reproduces the distribution function of draws simulated from the declared
  # priors, within their Monte Carlo error
  set.seed(1)
  n_draws <- 2e5
  draws   <- abs(stats::rnorm(n_draws, 0, 0.4)) *
    sqrt(stats::rbeta(n_draws, alpha[[1L]], alpha[[2L]]))
  curve   <- ggplot2::ggplot_build(
    plot_prior(fit, "Study: tau", plot_type = "ggplot")
  )[["data"]][[1L]]
  curve   <- curve[is.finite(curve[["y"]]), c("x", "y")]
  curve   <- curve[order(curve[["x"]]), ]
  expect_gt(nrow(curve), 100L)
  cumulative <- c(0, cumsum(
    diff(curve[["x"]]) * (utils::head(curve[["y"]], -1L) + utils::tail(curve[["y"]], -1L)) / 2
  ))
  expect_equal(max(cumulative), 1, tolerance = 5e-3)
  points <- unname(stats::quantile(draws, seq(0.1, 0.9, by = 0.2)))
  from_curve <- stats::approx(curve[["x"]], cumulative, xout = points, ties = "ordered")[["y"]]
  from_draws <- vapply(points, function(point) mean(draws <= point), numeric(1))
  expect_lt(
    max(abs(from_curve - from_draws) / sqrt(from_draws * (1 - from_draws) / n_draws)), 4
  )

  printed <- capture.output(returned <- print_prior(fit, "Study: tau"))
  expect_match(printed[[1L]], "Study: tau = sigma_total * sqrt(w[1])", fixed = TRUE)
  expect_true(all(paste0("    ", print(total, silent = TRUE)) %in% printed))
  expect_true(all(paste0("    ", print(weights, silent = TRUE)) %in% printed))
  printed <- capture.output(print_prior(fit, "Crit: tau2"))
  expect_match(printed[[1L]], "Crit: tau2 = (sigma_total * sqrt(w[2]))^2", fixed = TRUE)
  .random_prior_device(expect_null(plot_prior(fit, "Study: tau")))
})

test_that("the SD of a scale-modelled block keeps resolving to the scale intercept", {

  skip_on_cran()
  fit <- load_fit("plot_prior_random_scale")

  selected <- .select_plot_prior_parameter(
    fit, "tau", NULL, NULL, component = "auto"
  )
  expect_identical(selected[["source"]], "scale")
  expect_identical(selected[["label"]], "log_tau_intercept")
  expect_true(BayesTools::is.prior(selected[["prior"]]))
  printed <- capture.output(print_prior(fit, "tau"))
  expect_identical(printed[[1L]], "log_tau_intercept:")
  .random_prior_device(expect_null(plot_prior(fit, "tau")))
})
