test_that("paired scalar z integration preserves fitted, extrapolated and inactive laws", {

  cutoff <- stats::qnorm(.025, lower.tail = FALSE)
  sei <- c(.2, .3)
  location <- matrix(c(.1, -.1, .2, -.15, .3, -.2), 3L)
  sd <- matrix(c(.4, .5, .3, .7, .9, .6), 3L)
  retained <- matrix(c(.8, 0, .2, .3, 1e-12, .5), 3L)
  extrapolated_location <- location
  extrapolated_location[3L, ] <- location[3L, ] + .1
  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .2)))
  context <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), sei,
    effect_direction = "positive", signed_data = FALSE)
  context$omega <- matrix(c(1, .2, 1, 1, 1, .2), 3L, byrow = TRUE)
  context$alpha <- numeric(3L)
  context$phack_kind <- integer(3L)
  context$kernel_mode <- c(SELKERNEL_STEP, SELKERNEL_STEP, SELKERNEL_NORMAL)
  context$vector_rule <- integer(3L)
  context$use_normal <- c(FALSE, FALSE, TRUE)
  z <- sort(unique(c(seq(-6, 10, length.out = 241L), cutoff)))
  selected <- match(c(z[1L], z[51L], cutoff, z[181L], z[length(z)]), z)
  control <- set_selection_likelihood_control(relative_tolerance = 1e-7)
  result <- .zplot_latent_mixture(z, location, sd, retained, sei, context, FALSE,
    control, extrapolated_location)

  # Gaussian conditioning gives an independent adaptive integral for H.
  reference <- function(value, column) {
    total <- sqrt(sd[1L, column]^2 + retained[1L, column]^2)
    conditional_mean <- location[1L, column] + (retained[1L, column] / total)^2 *
      (value * sei[column] - location[1L, column])
    conditional_sd <- sd[1L, column] * retained[1L, column] / total
    H <- stats::integrate(function(u) stats::dnorm(u) / (.2 + .8 * stats::pnorm(
      cutoff * sei[column], conditional_mean + conditional_sd * u,
      sd[1L, column], lower.tail = FALSE)), -Inf, Inf, rel.tol = 1e-11)$value
    sei[column] * stats::dnorm(value * sei[column], location[1L, column], total) * H
  }
  expected <- vapply(z[selected], function(value) mean(vapply(seq_along(sei),
    function(column) reference(value, column), numeric(1L))), numeric(1L))
  expect_equal(result$extrapolated[1L, selected], expected, tolerance = 2e-7)
  expect_equal(result$fitted[1L, selected], expected * ifelse(z[selected] >= cutoff, 1, .2),
    tolerance = 2e-7)
  for (row in 2:3) {
    normal <- function(means) vapply(z, function(value) mean(sei * stats::dnorm(
      value * sei, means, sqrt(sd[row, ]^2 + retained[row, ]^2))), numeric(1L))
    expect_equal(result$fitted[row, ], normal(location[row, ]), tolerance = 1e-12)
    expect_equal(result$extrapolated[row, ], normal(extrapolated_location[row, ]), tolerance = 1e-12)
  }
  inverse_area <- vapply(seq_along(sei), function(column) {
    stats::integrate(function(u) stats::dnorm(u) / (.2 + .8 * stats::pnorm(
      cutoff * sei[column], location[1L, column] + retained[1L, column] * u,
      sd[1L, column], lower.tail = FALSE)), -Inf, Inf, rel.tol = 1e-11)$value
  }, numeric(1L))
  expect_equal(result$weights, c(mean(inverse_area), 1, 1), tolerance = 2e-7)
  tail <- .zplot_latent_mixture(2, location, sd, retained, sei, context, TRUE,
    control, extrapolated_location)
  tail_area <- vapply(seq_along(sei), function(column) {
    stats::integrate(function(u) {
      mu <- location[1L, column] + retained[1L, column] * u
      stats::dnorm(u) * (stats::pnorm(-2 * sei[column], mu, sd[1L, column]) +
        stats::pnorm(2 * sei[column], mu, sd[1L, column], lower.tail = FALSE)) /
        (.2 + .8 * stats::pnorm(cutoff * sei[column], mu, sd[1L, column], lower.tail = FALSE))
    }, -Inf, Inf, rel.tol = 1e-11)$value
  }, numeric(1L))
  expect_equal(tail$EDR[1L], mean(tail_area) / mean(inverse_area), tolerance = 2e-7)

  # Only the actual zero-weight interval justifies a structural zero display.
  zero_context <- BayesTools::selection_context_subset_observations(
    BayesTools::selection_context_subset_rows(context, 1L), 1L)
  zero_context$omega <- matrix(c(1, 0), 1L)
  zero <- .zplot_latent_mixture(seq(-6, -2, length.out = 241L), matrix(.1),
    matrix(.4), matrix(.1), .2, zero_context, FALSE, control, fitted_only = TRUE)
  expect_true(all(zero$fitted == 0))

  # The native integration diagnostic must remain active even if consecutive
  # displayed curves happen to coincide exactly.
  testthat::local_mocked_bindings(
    .gauss_hermite_nodes = function(...) list(nodes = 0, weights = 1, log_weights = 0),
    .zplot_selnorm_density_matrix = function(z_sequence, mean, ...) {
      structure(matrix(.5, nrow(mean), length(z_sequence)), relative_integration_error = rep(.02, nrow(mean)))
    },
    .package = "RoBMA"
  )
  expect_error(.zplot_latent_mixture(z, matrix(.1), matrix(.4), matrix(.8), .2,
    zero_context, FALSE, control, fitted_only = TRUE),
    paste0("Zplot marginal integration was rejected by diagnostics: relative integration error was 0.02. ",
      "Inspect the fitted selection and heterogeneity parameters."), fixed = TRUE)
})
