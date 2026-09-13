.zplot_context_test_selection <- function(sei, cutoff = 1.96, weights = c(1, .2),
                                          direction = "positive") {

  prior <- BayesTools::prior_weightfunction("one-sided",
    stats::pnorm(cutoff, lower.tail = FALSE), BayesTools::wf_fixed(weights),
    model = BayesTools::selection_model(weight_rule = "product", group = "paper"))
  context <- .selection_spec(list(outcome = list(bias = prior)), rep(0, length(sei)), sei, direction)
  context$omega <- matrix(context$fixed_omega[1L, ], 1L)
  context$alpha <- 0
  context$phack_kind <- 0L
  context$kernel_mode <- SELKERNEL_STEP
  context$vector_rule <- 0L
  context$use_normal <- FALSE
  BayesTools::selection_context_validate(context, n_samples = 1L)
}


test_that("context projection keeps the full bivariate normalizer inside the population mixture", {

  sei <- c(.7, 1.3)
  covariance <- matrix(c(.61, .364, .364, 1.81), 2L)
  mean <- c(.2, .2)
  context_sd <- sqrt(.3)
  weight <- .2
  cutoff <- 1.96 * sei
  sd <- sqrt(diag(covariance))
  schur <- diag(covariance) - covariance[1L, 2L]^2 / rev(diag(covariance))
  z <- c(-.5, 1.95, 1.97)
  inner_error <- outer_error <- 0
  # Independent QUADPACK reference: expand A(b) into two univariate tails
  # and a bivariate tail, then condition the focal density on its companion.
  # The original full covariance is used throughout; no native factor rule.
  acceptance <- function(location) {

    marginal <- stats::pnorm((location - cutoff) / sd)
    joint <- stats::integrate(function(u) stats::dnorm(u) * stats::pnorm(
      (location[2L] + covariance[2L, 1L] / sd[1L] * u - cutoff[2L]) / sqrt(schur[2L])),
      (cutoff[1L] - location[1L]) / sd[1L], Inf,
      rel.tol = 1e-10, abs.tol = 1e-12, subdivisions = 200L)
    stopifnot(identical(joint$message, "OK"))
    inner_error <<- max(inner_error, (1 - weight)^2 * joint$abs.error)
    weight^2 + weight * (1 - weight) * sum(marginal) + (1 - weight)^2 * joint$value
  }
  expected <- vapply(z, function(value) {
    integral <- stats::integrate(function(u) vapply(u, function(node) {
      population <- stats::dnorm(node)
      if (population == 0) return(0)
      location <- mean + context_sd * node
      x <- value * sei
      conditional_mean <- rev(location) + covariance[1L, 2L] / diag(covariance) * (x - location)
      companion <- weight + (1 - weight) * stats::pnorm(
        (conditional_mean - rev(cutoff)) / sqrt(rev(schur)))
      observed <- if (value >= 1.96) 1 else weight
      population * mean(sei * stats::dnorm(x, location, sd) * observed * companion) / acceptance(location)
    }, numeric(1L)), -Inf, Inf, rel.tol = 1e-9, abs.tol = 1e-11, subdivisions = 200L)
    stopifnot(identical(integral$message, "OK"))
    outer_error <<- max(outer_error, integral$abs.error)
    integral$value
  }, numeric(1L))
  reference_error <- outer_error + max(sei / sd) / sqrt(2 * pi) * inner_error / weight^4
  expect_lt(reference_error, 1e-7)
  control <- set_selection_likelihood_control(relative_tolerance = .001)
  for (sign in c(1, -1)) {
    context <- .zplot_context_test_selection(sei, direction = if (sign > 0) "positive" else "negative")
    actual <- .zplot_context_projection(sign * z, sign * mean, covariance,
      matrix(rep(sign * context_sd, 2L), 1L), sei, context, control)
    expect_type(actual, "list")
    expect_identical(dim(actual$density), c(1L, length(z)))
    expect_lte(actual$relative_error, control$relative_tolerance)
    expect_lte(actual$mass_error, control$relative_tolerance)
    # Reflection gives the same independent reference for negative direction.
    expect_lt(max(abs(as.numeric(actual$density) - expected)), 1e-3)
  }
})


test_that("point contexts agree with an analytic Gaussian orthant identity", {

  sei <- c(.7, 1.3)
  covariance <- matrix(c(.61, .364, .364, 1.81), 2L)
  sd <- sqrt(diag(covariance))
  weight <- .2
  correlation <- covariance[1L, 2L] / prod(sd)
  acceptance <- ((1 + weight) / 2)^2 + (1 - weight)^2 * asin(correlation) / (2 * pi)
  schur <- diag(covariance) - covariance[1L, 2L]^2 / rev(diag(covariance))
  z <- c(-.7, .3, 1.6)
  expected <- vapply(z, function(value) {
    x <- value * sei
    companion <- weight + (1 - weight) * stats::pnorm(
      covariance[1L, 2L] / diag(covariance) * x / sqrt(rev(schur)))
    observed <- if (value > 0) 1 else weight
    mean(sei * stats::dnorm(x, 0, sd) * observed * companion) / acceptance
  }, numeric(1L))
  control <- set_selection_likelihood_control(relative_tolerance = .001)
  actual <- .zplot_context_projection(z, c(0, 0), covariance, matrix(0, 1L, 2L),
    sei, .zplot_context_test_selection(sei, cutoff = 0), control)
  expect_type(actual, "list")
  expect_lte(actual$relative_error, control$relative_tolerance)
  expect_lte(actual$mass_error, control$relative_tolerance)
  expect_lt(max(abs(as.numeric(actual$density) - expected)), 1e-3)
})


test_that("constant selection weights retain the full Gaussian marginal", {

  sei <- c(.7, 1.3)
  mean <- c(.15, -.1)
  covariance <- matrix(c(.61, .364, .364, 1.81), 2L)
  loading <- matrix(c(.4, -.1, .2, .35), 2L)
  z <- c(-2, -.1, .3, 1.96, 3)
  context <- .zplot_context_test_selection(sei, weights = c(1, 1))
  actual <- .zplot_context_projection(z, mean, covariance, loading, sei,
    context, set_selection_likelihood_control())
  expected <- vapply(z, function(value) mean(sei * stats::dnorm(value * sei,
    mean, sqrt(diag(covariance) + colSums(loading^2)))), numeric(1L))
  expect_lt(max(abs(as.numeric(actual$density) - expected)), 1e-12)
  expect_identical(actual$relative_error, 0)
  expect_identical(actual$mass_error, 0)
  expect_false(context$use_normal)
})


test_that("input mixture pruning preserves raw anchors and reports bounded omission", {

  sei <- c(.8, 1.1, 1.4)
  correlation <- matrix(.4, 3L, 3L)
  diag(correlation) <- 1
  covariance <- correlation * tcrossprod(sei)
  means <- matrix(rep(seq(-.6, .6, length.out = 9L), 3L), ncol = 3L)
  rule <- .gauss_hermite_nodes(31L)
  context_weights <- rep(-log(nrow(means)), nrow(means))
  evaluate <- function(allowance, weights = context_weights) {

    .Call("RoBMA_selnorm_context_star_zplot", means,
      as.double(covariance[lower.tri(covariance, diag = TRUE)]), sei,
      c(5, .2), c(0, -Inf), c(Inf, 0), 1L, TRUE, NULL, weights,
      rule$nodes, rule$log_weights, seq(-4, 4, length.out = 65L), allowance,
      PACKAGE = "RoBMA")
  }
  set.seed(419L)
  before <- .Random.seed
  exact <- evaluate(0)
  pruned <- evaluate(1e-4)
  mass_only <- evaluate(1e-4, rep(-Inf, nrow(means)))
  expect_true(exact$available)
  expect_true(pruned$available)
  omitted <- attr(pruned$compression_error, "omitted_mass_error", exact = TRUE)
  expect_type(omitted, "double")
  expect_length(omitted, 1L)
  expect_true(is.finite(omitted))
  expect_gt(omitted, 0)
  expect_identical(attr(exact$compression_error, "omitted_mass_error", exact = TRUE), 0)
  expect_identical(pruned$mass, exact$mass)
  expect_identical(pruned$compact_log_normalizers, exact$compact_log_normalizers)
  expect_lte(max(abs(pruned$density - exact$density)), as.double(pruned$compression_error))
  expect_false(mass_only$available)
  expect_identical(mass_only$mass, 0)
  expect_identical(as.double(mass_only$compression_error), 0)
  expect_identical(attr(mass_only$compression_error, "omitted_mass_error", exact = TRUE), 0)
  expect_identical(mass_only$compact_log_normalizers, exact$compact_log_normalizers)
  expect_identical(.Random.seed, before)
})
