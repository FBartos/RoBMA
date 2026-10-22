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


.zplot_context_test_bivariate_density <- function(locations, covariance, sei, z, weight = .2) {

  sd <- sqrt(diag(covariance))
  schur <- diag(covariance) - covariance[1L, 2L]^2 / rev(diag(covariance))
  # Independent bivariate Gaussian probabilities give the complete event
  # normalizer for every retained realization, before averaging the curves.
  normalizers <- apply(locations, 1L, function(location) {
    joint <- as.numeric(mvtnorm::pmvnorm(lower = c(0, 0), upper = c(Inf, Inf),
      mean = location, sigma = covariance, algorithm = mvtnorm::TVPACK()))
    weight^2 + weight * (1 - weight) * sum(stats::pnorm(location / sd)) +
      (1 - weight)^2 * joint
  })
  matrix(vapply(z, function(value) {
    observed <- if (value >= 0) 1 else weight
    vapply(seq_len(nrow(locations)), function(row) {
      location <- locations[row, ]
      x <- value * sei
      conditional_mean <- rev(location) + covariance[1L, 2L] / diag(covariance) * (x - location)
      companion <- weight + (1 - weight) * stats::pnorm(conditional_mean / sqrt(rev(schur)))
      mean(sei * stats::dnorm(x, location, sd) * observed * companion) / normalizers[row]
    }, numeric(1L))
  }, numeric(nrow(locations))), nrow(locations), length(z))
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


test_that("a batch of point contexts is the average of their full-event densities", {

  # Each retained realization keeps its own full-event normalizer, so the batch
  # target is the average of ratios, checked here against the independent
  # bivariate reference. A batch of one is the single-context projection.
  sei <- c(.7, 1.3)
  covariance <- matrix(c(.61, .364, .364, 1.81), 2L)
  z <- c(-.7, .3, 1.6, 2.4)
  context <- .zplot_context_test_selection(sei, cutoff = 0)
  control <- set_selection_likelihood_control(relative_tolerance = .001)
  locations <- rbind(c(-.4, .1), c(0, 0), c(.3, .5), c(.8, -.2), c(-.1, .9))

  single <- .zplot_context_projection(z, locations[3L, ], covariance,
    matrix(0, 1L, 2L), sei, context, control)
  expect_identical(
    .zplot_context_projection(z, locations[3L, ], covariance, matrix(0, 1L, 2L),
      sei, context, control, context_means = locations[3L, , drop = FALSE]),
    single
  )

  batch <- .zplot_context_projection(z, colMeans(locations), covariance,
    matrix(0, 1L, 2L), sei, context, control, context_means = locations)
  expect_type(batch, "list")
  expect_lte(batch$relative_error, control$relative_tolerance)
  expect_lte(batch$mass_error, control$relative_tolerance)
  expected <- colMeans(.zplot_context_test_bivariate_density(locations, covariance, sei, z))
  expect_lt(max(abs(as.numeric(batch$density) - expected)),
            batch$integration_error[["absolute"]] + 1e-6)
  # The same nodes projected one at a time give the same curve.
  separate <- colMeans(do.call(rbind, lapply(seq_len(nrow(locations)), function(row) {
    .zplot_context_projection(z, locations[row, ], covariance, matrix(0, 1L, 2L),
      sei, context, control)$density[1L, ]
  })))
  expect_lt(max(abs(as.numeric(batch$density) - separate)), 1e-4 * max(separate))

  # A batch must stay a set of point contexts.
  expect_error(.zplot_context_projection(z, colMeans(locations), covariance,
    matrix(.1, 1L, 2L), sei, context, control, context_means = locations),
    "batched context projection inputs are invalid", fixed = TRUE)
})


test_that("equivalent selection partitions give identical checked context projections", {

  sei <- c(.7, 1.3)
  covariance <- matrix(c(.61, .364, .364, 1.81), 2L)
  cases <- list(
    list(cuts = .5, weights = c(1, .2),
      refined_cuts = c(.1, .5, .9), refined_weights = c(1, 1, .2, .2)),
    list(cuts = c(.25, .75), weights = c(1, .2, 1),
      refined_cuts = c(.05, .25, .5, .75, .95), refined_weights = c(1, 1, .2, .2, 1, 1))
  )
  for (case in cases) {
    for (direction in c("positive", "negative")) {
      context <- .zplot_context_test_selection(sei,
        cutoff = stats::qnorm(case$cuts, lower.tail = FALSE),
        weights = case$weights, direction = direction)
      refined <- .zplot_context_test_selection(sei,
        cutoff = stats::qnorm(case$refined_cuts, lower.tail = FALSE),
        weights = case$refined_weights, direction = direction)
      before <- BayesTools::selection_native_static_args(refined)
      evaluate <- function(selection) .zplot_context_projection(c(-2, -.6, 0, .7, 2),
        c(.15, -.1), covariance, matrix(.2, 1L, 2L), sei, selection,
        set_selection_likelihood_control(relative_tolerance = .001))
      expected <- evaluate(context)
      actual <- evaluate(refined)
      expect_type(actual, "list")
      expect_identical(actual, expected)
      expect_identical(BayesTools::selection_native_static_args(refined), before)
    }
  }
})


test_that("multidimensional retained contexts keep each full-event normalizer", {

  sei <- c(.7, 1.3)
  covariance <- matrix(c(.61, .364, .364, 1.81), 2L)
  latent <- diag(c(.02, .04))
  mean <- c(.15, -.1)
  weight <- .2
  z <- c(-.7, .3, 1.6)
  control <- set_selection_likelihood_control(points_per_scramble = 64L,
    max_points_per_scramble = 64L, scrambles = 4L, relative_tolerance = .01)
  uniforms <- BayesTools::selection_qmc_design(dimensions = 2L, points = 64L,
    scrambles = 4L, seed = control$seed)
  locations <- sweep(matrix(stats::qnorm(uniforms), 256L, 2L) %*% chol(latent), 2L, mean, "+")
  expected <- colMeans(.zplot_context_test_bivariate_density(locations, covariance, sei, z, weight))
  actual <- .zplot_full_event_context_draw(z, mean, covariance, latent, sei,
    .zplot_context_test_selection(sei, cutoff = 0), FALSE, control, list())
  # The reference uses exactly the same outer QMC nodes. This isolates the
  # checked inner approximation from uncertainty in the outer population mean.
  expect_equal(as.numeric(actual), expected, tolerance = 1e-3)
})


test_that("negligible conditional tails use the aggregate density error scale", {

  covariance <- matrix(c(.0115, .006544158, .006544158, .0076), 2L)
  sei <- sqrt(diag(covariance))
  z <- seq(-6, 10, .4)
  context <- .zplot_context_test_selection(sei, cutoff = 0)
  control <- set_selection_likelihood_control(points_per_scramble = 8L,
    max_points_per_scramble = 8L, scrambles = 2L)
  locations <- rbind(matrix(.1, 4L, 2L),
    matrix(c(2.359389, 2.051164), 4L, 2L, byrow = TRUE))
  uniforms <- array(NA_real_, c(2L, 8L, 2L))
  for (column in seq_len(2L)) uniforms[, , column] <-
    rep(stats::pnorm(locations[, column]), each = 2L)
  testthat::local_mocked_bindings(selection_qmc_design = function(dimensions, points, scrambles, seed) {
    stopifnot(dimensions == 2L, points == 8L, scrambles == 2L, seed == control$seed)
    uniforms
  }, .package = "BayesTools")
  locations <- matrix(stats::qnorm(uniforms), 16L, 2L)
  expected <- colMeans(.zplot_context_test_bivariate_density(locations, covariance, sei, z))
  # Standalone point-context evaluation retains its relative requirement.
  expect_null(.zplot_context_projection(z, locations[16L, ], covariance,
    matrix(0, 1L, 2L), sei, context, control))
  actual <- .zplot_full_event_context_draw(z, c(0, 0), covariance, diag(2L), sei,
    context, FALSE, control, list())
  expect_lt(max(abs(as.numeric(actual) - expected)), control$relative_tolerance * max(expected))
})


test_that("the aggregate curve rejects non-negligible conditional errors", {

  sei <- c(.7, 1.3)
  context <- .zplot_context_test_selection(sei, cutoff = 0)
  control <- set_selection_likelihood_control(points_per_scramble = 8L,
    max_points_per_scramble = 8L, scrambles = 2L)
  testthat::local_mocked_bindings(.zplot_context_projection = function(
      z, mean, covariance, context_factor, sei, selection, control, absolute_tolerance = NULL,
      context_means = NULL) {
    list(density = matrix(1e-4, 1L, length(z)), relative_error = Inf, mass_error = 0,
      integration_error = c(absolute = absolute_tolerance / 2))
  })
  expect_error(.zplot_full_event_context_draw(c(-1, 0, 1), c(0, 0),
    matrix(c(.61, .364, .364, 1.81), 2L), diag(2L), sei,
    context, FALSE, control, list()), "retained-context integration was rejected by diagnostics")
})


test_that("retained-context acceleration preserves diagnostic rejection", {

  sei <- c(.7, 1.3)
  context <- .zplot_context_test_selection(sei)
  control <- set_selection_likelihood_control()
  testthat::local_mocked_bindings(.zplot_context_projection = function(...) {
    list(density = matrix(.1, 1L, 3L), relative_error = 2 * control$relative_tolerance,
      mass_error = 0)
  })
  expect_error(.zplot_full_event_context_draw(c(-1, 0, 1), c(.1, .2),
    matrix(c(.61, .364, .364, 1.81), 2L), matrix(0, 2L, 2L), sei,
    context, FALSE, control, list()), "rejected by diagnostics")
})


test_that("unavailable context geometry retains the full-covariance fallback", {

  sei <- c(.7, 1.3)
  covariance <- matrix(c(.61, -.364, -.364, 1.81), 2L)
  sd <- sqrt(diag(covariance))
  weight <- .2
  acceptance <- ((1 + weight) / 2)^2 +
    (1 - weight)^2 * asin(covariance[1L, 2L] / prod(sd)) / (2 * pi)
  schur <- diag(covariance) - covariance[1L, 2L]^2 / rev(diag(covariance))
  z <- c(-.7, .3, 1.6)
  expected <- vapply(z, function(value) {
    x <- value * sei
    companion <- weight + (1 - weight) * stats::pnorm(
      covariance[1L, 2L] / diag(covariance) * x / sqrt(rev(schur)))
    observed <- if (value > 0) 1 else weight
    mean(sei * stats::dnorm(x, 0, sd) * observed * companion) / acceptance
  }, numeric(1L))
  context <- .zplot_context_test_selection(sei, cutoff = 0)
  control <- set_selection_likelihood_control(relative_tolerance = .001)
  expect_null(.zplot_context_projection(z, c(0, 0), covariance,
    matrix(0, 1L, 2L), sei, context, control))
  actual <- .zplot_full_event_context_draw(z, c(0, 0), covariance,
    matrix(0, 2L, 2L), sei, context, FALSE, control, list())
  expect_equal(as.numeric(actual), expected, tolerance = 1e-3)
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
  withr::local_seed(419L)
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


test_that("the context projection is invariant to the native thread count", {

  # Contexts and mixture rows are evaluated in parallel but reduced in their
  # original order, so a projection must not depend on the thread budget.
  sei <- c(.7, 1.3, .9, 1.1, .8, 1.4)
  # Compound symmetry with a positive residual on every row, which is what the
  # covariance envelope the projection builds requires.
  covariance <- diag(.5 * sei^2) + .3 * tcrossprod(sei)
  z <- seq(-3, 3, length.out = 61L)
  control <- set_selection_likelihood_control(relative_tolerance = .005)
  selection <- .zplot_context_test_selection(sei, cutoff = 1.96,
                                             weights = c(1, .35))
  previous <- RoBMA.get_option("native_threads")
  on.exit(RoBMA.options(native_threads = previous), add = TRUE)

  reference <- NULL
  for (threads in c(1L, 2L, 3L, 8L)) {
    RoBMA.options(native_threads = threads)
    actual <- .zplot_context_projection(z, rep(.2, 6L), covariance,
      matrix(sqrt(.3), 1L, 6L), sei, selection, control)
    expect_type(actual, "list")
    values <- c(as.numeric(actual$density), actual$relative_error,
                actual$mass_error, actual$mass_omission_error)
    expect_true(all(is.finite(values)))
    if (is.null(reference)) {
      reference <- values
      expect_gt(max(actual$density), 0)
    } else {
      expect_identical(values, reference)
    }
  }
})


# The transition rule applied at each row of `means`, the definition that
# .selection_transition_rule_resolves_rows() must reproduce in one pass.
.zplot_context_test_rowwise_rule <- function(nodes, means, residual_sd, loading, sei,
                                             selection, axes = seq_len(ncol(loading))) {

  all(vapply(seq_len(nrow(means)), function(row) {
    .selection_transition_rule_resolves(nodes, means[row, ], residual_sd, loading, sei,
      selection, axes)
  }, logical(1L)))
}


test_that("the vectorized transition rule decides every branch like the rule at each row", {

  # Nodes are dyadic, so a center is exactly a node, and the gaps are
  # 1, .5, .5, .5, .5, 1. With coefficient 1, selection SE .5 and cut 2 the
  # center is 1 - mean.
  grid <- c(-2, -1, -.5, 0, .5, 1, 2)
  step <- function(weights = c(1, .2), z_lower = c(2, 9), sign = 1, kernel = SELKERNEL_STEP) {
    list(kernel_mode = kernel, omega = matrix(weights, 1L), z_lower = z_lower, sign = sign)
  }
  decide <- function(expected, centers = 0, nodes = grid, residual_sd = .75, loading = matrix(1),
                     sei = .5, selection = step(), axes = seq_len(ncol(loading)), means = NULL) {

    if (is.null(means)) means <- matrix(1 - centers)
    actual <- .selection_transition_rule_resolves_rows(nodes, means, residual_sd, loading, sei,
      selection, axes)
    reference <- .zplot_context_test_rowwise_rule(nodes, means, residual_sd, loading, sei,
      selection, axes)
    expect_identical(actual, expected)
    expect_identical(reference, expected)
  }

  # Interior centers within one gap resolve; exterior centers are not owned here.
  decide(TRUE, c(0, .25, .75))
  decide(TRUE, c(2, -2, 3, -3, 100, -100))
  # One failing row fails the set wherever it is, so no row is skipped.
  for (position in 1:4) {
    centers <- c(0, .25, .75, .1)
    centers[position] <- 1.5
    decide(FALSE, centers)
  }
  # A center on a node takes the larger adjacent gap (1 below -1, .5 above).
  decide(FALSE, -1)
  decide(TRUE, -1, residual_sd = 1)
  decide(TRUE, 0)
  decide(TRUE, 0, residual_sd = .5)
  decide(FALSE, 0, residual_sd = .49)
  # Nodes are sorted first, which drops NA nodes; duplicated, infinite and too
  # few nodes fail.
  decide(TRUE, .25, nodes = rev(grid))
  decide(TRUE, .25, nodes = c(grid, NA))
  decide(FALSE, .25, nodes = c(grid, grid[1L]))
  decide(FALSE, .25, nodes = c(grid, Inf))
  decide(FALSE, .25, nodes = c(-Inf, grid))
  decide(FALSE, .25, nodes = 0)
  decide(FALSE, .25, nodes = numeric())
  # Nothing to resolve: no rows, a normal kernel, unchanged weights, no axes,
  # an unselected cut, an inactive element.
  decide(TRUE, nodes = c(grid, grid[1L]), means = matrix(0, 0L, 1L))
  decide(TRUE, 1.5, nodes = c(grid, grid[1L]), selection = step(kernel = SELKERNEL_NORMAL))
  decide(TRUE, 1.5, selection = step(weights = c(1, 1)))
  decide(TRUE, 1.5, axes = integer())
  decide(TRUE, 1.5, selection = step(z_lower = c(-Inf, 9)))
  decide(TRUE, 1.5, loading = matrix(0))
  # Invalid scales, coefficients and means on active elements fail; the same
  # values on an inactive element are ignored.
  for (residual_sd in list(0, -1, NA_real_, Inf)) decide(FALSE, 0, residual_sd = residual_sd)
  decide(FALSE, 0, loading = matrix(Inf))
  decide(TRUE, 0, loading = matrix(NA_real_))
  decide(TRUE, means = matrix(c(1, 1.25, .75, .9)))
  for (bad in c(NA, NaN, Inf, -Inf)) decide(FALSE, means = matrix(c(1, bad, .75, .9)))
  decide(TRUE, loading = rbind(1, 0), sei = c(.5, .5), residual_sd = c(10, NA),
         means = cbind(1 - c(0, .25), c(NA, Inf)))
  # The sign of the effect direction and each element's own SE place the center.
  decide(TRUE, selection = step(sign = -1), means = matrix(-1 - c(0, .25)))
  decide(FALSE, selection = step(sign = -1), means = matrix(-1 - c(0, 1.5)))
  decide(TRUE, loading = rbind(1, 2), sei = c(.5, 1), residual_sd = c(.75, 1.5),
         means = cbind(c(1, .75), c(-2, -2)))
  # A missing center fails (missing SE, or an infinite SE at a zero cut); an
  # infinite one is exterior.
  decide(FALSE, 0, sei = NA_real_)
  decide(FALSE, 0, selection = step(z_lower = c(0, 9)), sei = Inf)
  decide(TRUE, 0, sei = Inf)
  # An element coupled to another axis is decided by the node gap alone, not by
  # the means: the coupled gap here is 1 and the scale is compared to it.
  coupled <- cbind(c(1, 1), c(.5, 0))
  exterior <- matrix(1000, 3L, 2L)
  decide(FALSE, loading = coupled, residual_sd = c(.3, 5), sei = c(.5, .5), means = exterior)
  decide(TRUE, loading = coupled, residual_sd = c(1.2, 5), sei = c(.5, .5), means = exterior)
  # Scale .8 fails the axis width but resolves by its lattice width .471
  # (exponent 2 pi^2 .471 = 9.3 at gap 1); scale .3 has width .254 (5.0).
  decide(TRUE, loading = coupled, residual_sd = c(.8, 5), sei = c(.5, .5), means = exterior)
  # Only the intervals that reach [-3, 3] count: the widest gap here is 2.
  decide(TRUE, loading = coupled, residual_sd = c(2.5, 5), sei = c(.5, .5), means = exterior,
         nodes = c(-6, -4, grid, 4, 7))
  decide(TRUE, loading = coupled, residual_sd = c(.8, 5), sei = c(.5, .5), means = exterior,
         axes = integer())
})


test_that("the vectorized transition rule equals the rule at each row over a randomized grid", {

  withr::local_seed(2057L)
  # Dyadic node sets, regular and irregular, put centers exactly on nodes with
  # unequal gaps on either side; the Gauss-Hermite sets give generic centers.
  grid <- c(-2, -1, -.5, 0, .5, 1, 2)
  irregular <- c(-3, -2.5, -1, -.5, 0, 2, 2.25, 3.5)
  node_sets <- c(
    list(.gauss_hermite_nodes(15L)$nodes, .gauss_hermite_nodes(127L)$nodes, grid, irregular,
         rev(irregular)),
    list(c(grid, grid[1L]), c(grid, Inf), c(grid, NA), c(-Inf, grid), grid[4L], numeric()))
  node_prob <- c(.15, .05, .2, .35, .1, .03, .03, .03, .03, .02, .01)
  pool <- c(1, 2, .5, -1, -2, -.5)
  one <- function(values) values[sample.int(length(values), 1L)]
  n_cases <- 2400L
  actual <- reference <- logical(n_cases)
  for (case in seq_len(n_cases)) {
    k <- sample.int(4L, 1L)
    m <- sample.int(3L, 1L)
    loading <- matrix(0, k, m)
    for (element in seq_len(k)) {
      columns <- switch(one(c("zero", "single", "coupled")),
        zero = integer(), single = sample.int(m, 1L), coupled = sample.int(m, min(m, 2L)))
      for (column in columns) {
        loading[element, column] <- if (stats::runif(1L) < .6) one(pool) else 2 * stats::rnorm(1L)
      }
    }
    if (stats::runif(1L) < .04) loading[sample.int(k, 1L), sample.int(m, 1L)] <- Inf
    if (m == 1L && stats::runif(1L) < .04) loading[sample.int(k, 1L), 1L] <- NA_real_
    n_cuts <- sample.int(3L, 1L)
    omega <- sample(c(1, .5, .2), n_cuts + 1L, replace = TRUE)
    if (stats::runif(1L) < .2) omega[] <- 1
    if (n_cuts > 1L && stats::runif(1L) < .25) omega[2L] <- omega[1L]
    if (stats::runif(1L) < .02) omega[sample.int(n_cuts + 1L, 1L)] <- NA_real_
    z_lower <- sample(c(2, 1, .5, 0, -1, 3, 1.96, -Inf, Inf), n_cuts, replace = TRUE)
    if (stats::runif(1L) < .5) z_lower <- c(z_lower, 9)
    if (stats::runif(1L) < .03) z_lower <- z_lower[-length(z_lower)]
    kernel_mode <- one(list(SELKERNEL_STEP, SELKERNEL_STEP, SELKERNEL_STEP, SELKERNEL_NORMAL,
                            NULL))[[1L]]
    selection <- list(kernel_mode = kernel_mode, omega = matrix(omega, 1L), z_lower = z_lower,
      sign = one(c(1, -1)))
    sei <- if (stats::runif(1L) < .6) sample(c(.5, 1, 2, .25), k, replace = TRUE) else
      exp(stats::runif(k, -1, 1))
    if (stats::runif(1L) < .04) sei[sample.int(k, 1L)] <- one(c(NA_real_, Inf))
    if (stats::runif(1L) < .03) sei <- sei[1L]
    residual_sd <- if (stats::runif(1L) < .6) {
      sample(c(.25, .5, 1, 1.5, 2, 4, 10), k, replace = TRUE)
    } else exp(stats::runif(k, -3, 2))
    if (stats::runif(1L) < .05) residual_sd[sample.int(k, 1L)] <- one(c(0, NA_real_, Inf, -1))
    nodes <- node_sets[[sample.int(length(node_sets), 1L, prob = node_prob)]]
    centers <- if (any(is.finite(nodes))) nodes[is.finite(nodes)] else grid
    n_rows <- sample(0:4, 1L, prob = c(.04, .36, .3, .2, .1))
    means <- matrix(2 * stats::rnorm(n_rows * k), n_rows, k)
    for (row in seq_len(n_rows)) {
      # Put a transition of one element exactly on a node, far outside, or off.
      element <- sample.int(k, 1L)
      column <- which(loading[element, ] != 0)
      kind <- one(c("node", "node", "exterior", "none"))
      if (length(column) && kind != "none" && is.finite(sei[min(element, length(sei))]) &&
          any(is.finite(z_lower))) {
        column <- column[sample.int(length(column), 1L)]
        shift <- selection$sign * one(z_lower[is.finite(z_lower)]) * sei[element]
        center <- if (kind == "node") one(centers) else one(c(-1, 1)) * 1e3
        means[row, element] <- shift - loading[element, column] * center
      }
    }
    if (n_rows > 0L && stats::runif(1L) < .06) {
      means[sample.int(n_rows, 1L), sample.int(k, 1L)] <- one(c(NA_real_, Inf, -Inf, NaN))
    }
    axes <- switch(sample.int(5L, 1L, prob = c(.5, .15, .2, .1, .05)), seq_len(m), m, 1L,
      sort(sample.int(m, sample.int(m, 1L))), integer())
    actual[case] <- .selection_transition_rule_resolves_rows(nodes, means, residual_sd, loading,
      sei, selection, axes)
    reference[case] <- .zplot_context_test_rowwise_rule(nodes, means, residual_sd, loading,
      sei, selection, axes)
  }
  expect_identical(actual, reference)
  # Both decisions occur often, so agreement is not vacuous.
  expect_gt(mean(reference), .3)
  expect_lt(mean(reference), .85)
})


test_that("the context projection decides the inner transitions once per rule, not per node", {

  sei <- c(.7, 1.3)
  covariance <- matrix(c(.61, .364, .364, 1.81), 2L)
  control <- set_selection_likelihood_control(relative_tolerance = .001)
  context <- .zplot_context_test_selection(sei)
  row_wise <- .selection_transition_rule_resolves
  vectorized <- .selection_transition_rule_resolves_rows
  tally <- new.env(parent = emptyenv())
  tally$row_wise <- tally$rows <- list()
  tally$inner_axes <- logical()
  # Counts the calls of both rules while `code` runs.
  counting <- function(code) {

    tally$row_wise <- tally$rows <- list()
    tally$inner_axes <- logical()
    testthat::with_mocked_bindings(code,
      .selection_transition_rule_resolves = function(nodes, mean, residual_sd, loading, sei,
                                                     selection, axes = seq_len(ncol(loading))) {
        tally$row_wise[[length(tally$row_wise) + 1L]] <-
          c(nodes = length(nodes), axes = length(axes), axis = axes[[1L]], columns = ncol(loading))
        row_wise(nodes, mean, residual_sd, loading, sei, selection, axes)
      },
      .selection_transition_rule_resolves_rows = function(nodes, means, residual_sd, loading, sei,
                                                          selection, axes) {
        tally$rows[[length(tally$rows) + 1L]] <- c(nodes = length(nodes), rows = nrow(means))
        # The sweeps decide every axis but the last, which is the context.
        tally$inner_axes <- c(tally$inner_axes, identical(axes, seq_len(ncol(loading) - 1L)))
        vectorized(nodes, means, residual_sd, loading, sei, selection, axes)
      })
  }

  # A continuous context: the numerator and denominator sweeps over the outer
  # nodes are vectorized, once per (rule, outer rule); only the outer axis
  # itself is decided by the row-wise rule.
  projected <- counting(.zplot_context_projection(c(-.5, 1.95, 1.97), c(.2, .2), covariance,
    matrix(sqrt(.3), 1L, 2L), sei, context, control))
  expect_type(projected, "list")
  outer_calls <- do.call(rbind, tally$row_wise)
  sweeps <- do.call(rbind, tally$rows)
  expect_gt(nrow(outer_calls), 0L)
  expect_true(all(outer_calls[, "axes"] == 1L & outer_calls[, "axis"] == outer_calls[, "columns"]))
  expect_false(anyDuplicated(sweeps) > 0L)
  # The same inner rule is decided again for another outer rule, so the cache
  # key includes the outer rule.
  expect_true(anyDuplicated(sweeps[, "nodes"]) > 0L)
  # The first decision is the finest rule of every sweep, in a single call.
  expect_identical(unname(sweeps[1L, ]), c(1023L, 1023L))
  expect_true(all(sweeps[, "rows"] > 1L))
  expect_true(all(tally$inner_axes))
  # The result is that of the row-wise definition.
  reference <- testthat::with_mocked_bindings(
    .zplot_context_projection(c(-.5, 1.95, 1.97), c(.2, .2), covariance,
      matrix(sqrt(.3), 1L, 2L), sei, context, control),
    .selection_transition_rule_resolves_rows = .zplot_context_test_rowwise_rule)
  expect_identical(projected, reference)

  # A batch of point contexts has no outer axis: the rule is never called row by row.
  locations <- rbind(c(-.4, .1), c(0, 0), c(.3, .5), c(.8, -.2), c(-.1, .9))
  point_context <- .zplot_context_test_selection(sei, cutoff = 0)
  projected <- counting(.zplot_context_projection(c(-.7, .3, 1.6), colMeans(locations),
    covariance, matrix(0, 1L, 2L), sei, point_context, control, context_means = locations))
  expect_type(projected, "list")
  expect_length(tally$row_wise, 0L)
  sweeps <- do.call(rbind, tally$rows)
  expect_false(anyDuplicated(sweeps) > 0L)
  expect_true(all(sweeps[, "rows"] == nrow(locations)))
  expect_true(all(tally$inner_axes))
  reference <- testthat::with_mocked_bindings(
    .zplot_context_projection(c(-.7, .3, 1.6), colMeans(locations), covariance,
      matrix(0, 1L, 2L), sei, point_context, control, context_means = locations),
    .selection_transition_rule_resolves_rows = .zplot_context_test_rowwise_rule)
  expect_identical(projected, reference)
})
