.zplot_threshold_batch_context <- function(S, sei, direction = "positive", normal = FALSE) {

  prior <- BayesTools::prior_weightfunction("one-sided", c(.025, .05),
    BayesTools::wf_fixed(c(1, .6, .2)))
  selection <- .selection_spec(list(outcome = list(bias = prior)),
    rep(0, length(sei)), sei, effect_direction = direction, signed_data = FALSE)
  selection$omega <- matrix(rep(c(1, .6, .2), S), S, byrow = TRUE)
  selection$alpha <- numeric(S)
  selection$phack_kind <- integer(S)
  selection$use_normal <- rep_len(normal, S)
  selection$kernel_mode <- ifelse(selection$use_normal, SELKERNEL_NORMAL, SELKERNEL_STEP)
  selection$vector_rule <- integer(S)
  selection
}


# Keep the original R node reductions as an exact arithmetic reference.
.zplot_threshold_batch_reference <- function(z_threshold, mean, sd, sei,
                                              selection_context, mean_extrapolated,
                                              latent_sd, quadrature, fitted_only) {

  fitted <- extrapolated <- weights <- numeric(nrow(mean))
  for (j in seq_along(quadrature$nodes)) {
    fit <- .zplot_selnorm_threshold_summary(z_threshold,
      mean + latent_sd * quadrature$nodes[j], sd, sei, selection_context, FALSE)
    fitted <- fitted + quadrature$weights[j] * fit$EDR
    if (!fitted_only) {
      ext <- .zplot_selnorm_threshold_summary(z_threshold,
        mean_extrapolated + latent_sd * quadrature$nodes[j], sd, sei,
        selection_context, TRUE)
      extrapolated <- extrapolated + quadrature$weights[j] * ext$EDR * ext$weights
      weights <- weights + quadrature$weights[j] * ext$weights
    }
  }
  list(fitted = fitted, weights = weights, extrapolated = extrapolated)
}


test_that("threshold quadrature preserves scalar calls and R reduction order", {

  sei <- c(.3, .5, .3, .3, .8, .5)
  S <- 6L
  mean <- matrix(seq(-.4, .7, length.out = S * length(sei)), S)
  mean[, 3L] <- mean[, 1L]
  sd <- matrix(rep(sei + .2, each = S), S)
  latent_sd <- matrix(seq(0, .3, length.out = length(mean)), S)
  latent_sd[, 3L] <- latent_sd[, 1L]
  rule <- list(nodes = c(-3, -.7, 0, .2, 2.5), weights = c(0, .1, .2, .3, .4))
  for (direction in c("positive", "negative")) {
    for (normal in list(FALSE, TRUE, c(FALSE, TRUE))) {
      selection <- .zplot_threshold_batch_context(S, sei, direction, normal)
      selection$omega[2L, ] <- selection$omega[2L, ] * 1e280
      selection$omega[3L, ] <- selection$omega[3L, ] * 1e-280
      selection$omega[4L, ] <- c(1, 0, 0)
      selection$omega[5L, ] <- c(.1, .8, .2)
      for (fitted_only in c(FALSE, TRUE)) {
        args <- list(1.96, mean, sd, sei, selection, mean + .125,
          latent_sd, rule, fitted_only)
        expect_identical(do.call(.zplot_selnorm_threshold_quadrature, args),
          do.call(.zplot_threshold_batch_reference, args))
      }
    }
  }
})


test_that("threshold quadrature agrees with independent conditional Gaussian integrals", {

  for (direction in c("positive", "negative")) {
    selection <- .zplot_threshold_batch_context(1L, .5, direction)
    threshold <- 1.96 * .5
    residual_sd <- .7
    latent_sd <- .3
    fitted_mean <- -.2
    extrapolated_mean <- .35
    sign <- if (direction == "positive") 1 else -1
    lower <- selection$z_lower * .5
    upper <- selection$z_upper * .5
    omega <- as.numeric(selection$omega)
    interval <- function(lo, hi, location) {
      if (lo >= hi) return(0)
      stats::pnorm(hi, location, residual_sd) - stats::pnorm(lo, location, residual_sd)
    }
    components <- function(u) {
      location <- sign * (fitted_mean + latent_sd * u)
      location_ext <- sign * (extrapolated_mean + latent_sd * u)
      mass <- vapply(seq_along(omega), function(b) interval(lower[b], upper[b], location), numeric(1))
      mass_ext <- vapply(seq_along(omega), function(b) interval(lower[b], upper[b], location_ext), numeric(1))
      tail <- vapply(seq_along(omega), function(b) {
        interval(lower[b], min(upper[b], -threshold), location) +
          interval(max(lower[b], threshold), upper[b], location)
      }, numeric(1))
      norm <- sum(omega * mass)
      norm_ext <- sum(omega * mass_ext)
      c(fitted = sum(omega * tail) / norm,
        weights = 1 / norm_ext,
        extrapolated = (stats::pnorm(-threshold, location_ext, residual_sd) +
          stats::pnorm(threshold, location_ext, residual_sd, lower.tail = FALSE)) / norm_ext)
    }
    reference <- vapply(seq_len(3L), function(component) {
      stats::integrate(function(u) vapply(u, function(value) {
        components(value)[component] * stats::dnorm(value)
      }, numeric(1)), -9, 9, rel.tol = 1e-10)$value
    }, numeric(1))
    actual <- .zplot_selnorm_threshold_quadrature(1.96, matrix(fitted_mean),
      matrix(residual_sd), .5, selection, matrix(extrapolated_mean),
      matrix(latent_sd), .gauss_hermite_nodes(127L), FALSE)
    expect_equal(unlist(actual, use.names = FALSE), reference, tolerance = 2e-10)
  }
})


test_that("threshold quadrature preserves nonfinite values and validates its new shapes", {

  selection <- .zplot_threshold_batch_context(3L, c(.3, .3))
  args <- list(z_threshold = 1.96, mean = matrix(c(NA, Inf, 0), 3L, 2L),
    sd = matrix(c(.5, .5, 0), 3L, 2L), sei = c(.3, .3),
    selection_context = selection, mean_extrapolated = matrix(c(0, -Inf, NA), 3L, 2L),
    latent_sd = matrix(c(.1, 0, Inf), 3L, 2L),
    quadrature = list(nodes = c(-1, 0, 1), weights = c(.25, .5, .25)), fitted_only = FALSE)
  expect_identical(do.call(.zplot_selnorm_threshold_quadrature, args),
    do.call(.zplot_threshold_batch_reference, args))
  invalid <- args
  invalid$latent_sd <- matrix(0, 2L, 2L)
  expect_error(do.call(.zplot_selnorm_threshold_quadrature, invalid), "latent_sd", fixed = TRUE)
  invalid <- args
  invalid$mean_extrapolated <- matrix(0, 2L, 2L)
  expect_error(do.call(.zplot_selnorm_threshold_quadrature, invalid), "mean_extrapolated", fixed = TRUE)
  invalid <- args
  invalid$quadrature$weights <- numeric()
  expect_error(do.call(.zplot_selnorm_threshold_quadrature, invalid), "same positive length", fixed = TRUE)
  invalid <- args
  invalid$fitted_only <- NA
  expect_error(do.call(.zplot_selnorm_threshold_quadrature, invalid), "must not be missing", fixed = TRUE)
  invalid <- args
  invalid$selection_context$omega[1L, 1L] <- -1
  expect_error(do.call(.zplot_selnorm_threshold_quadrature, invalid), "omega", fixed = TRUE)
})


test_that("batched threshold rules preserve the adaptive R controller", {

  S <- 5L
  sei <- c(.3, .5, .3)
  selection <- .zplot_threshold_batch_context(S, sei, normal = c(FALSE, TRUE))
  mean <- matrix(seq(-.2, .3, length.out = S * length(sei)), S)
  sd <- matrix(rep(sei + .2, each = S), S)
  latent_sd <- matrix(rep(seq(.05, .5, length.out = S), length(sei)), S)
  for (fitted_only in c(FALSE, TRUE)) {
    evaluate <- function() .zplot_latent_mixture(1.96, mean, sd, latent_sd, sei,
      selection, TRUE, set_selection_likelihood_control(), mean + .125, fitted_only)
    reference <- testthat::with_mocked_bindings(evaluate(),
      .zplot_selnorm_threshold_quadrature = .zplot_threshold_batch_reference,
      .package = "RoBMA")
    expect_identical(evaluate(), reference)
  }
})


test_that("repeated transition columns preserve valid and invalid row gates", {

  sei <- rep(c(.3, .5), 4L)
  names(sei) <- paste0("estimate", seq_along(sei))
  S <- 6L
  selection <- .zplot_threshold_batch_context(S, sei, normal = c(FALSE, TRUE))
  mean <- matrix(rep(c(0, -.1, .2, NA, Inf, -0), length(sei)), S)
  sd <- matrix(.5, S, length(sei))
  latent <- matrix(rep(c(.1, .2, 0, .3, .5, -0), length(sei)), S)
  # Matching SE alone must never identify a shared transition predicate.
  mean[1L, 3L] <- .3
  sd[3L, 5L] <- .05
  latent[5L, 7L] <- .01
  original <- .selection_scalar_transition_rule_resolves
  calls <- 0L
  testthat::with_mocked_bindings(
    .zplot_latent_transition_resolves(c(-2, -1, 0, 1, 2), mean, sd, latent, sei, selection),
    .selection_scalar_transition_rule_resolves = function(...) {
      calls <<- calls + 1L
      original(...)
    }, .package = "RoBMA")
  expect_identical(calls, 5L)
  for (sign in c(1L, -1L)) for (nodes in list(c(-2, -1, 0, 1, 2),
      .gauss_hermite_nodes(31L)$nodes, c(NA_real_, 0), numeric())) {
    selection$sign <- sign
    reference <- rep(TRUE, S)
    for (column in seq_along(sei)) {
      reference <- reference & .selection_scalar_transition_rule_resolves(nodes,
        mean[, column], sd[, column], latent[, column], sei[column], selection, seq_len(S))
    }
    expect_identical(.zplot_latent_transition_resolves(nodes, mean, sd, latent, sei, selection),
      reference)
  }
})


test_that("threshold quadrature uses the same row reductions on actual native threads", {

  source(testthat::test_path("helper-selection-kernel.R"), local = TRUE)
  .with_kernel_thread_worker({
    inputs <- .kernel_thread_inputs(S = 256L)
    selection <- .kernel_thread_selection_context(inputs)
    rule <- .gauss_hermite_nodes(31L)
    .expect_thread_invariant("zcurve threshold quadrature", function() {
      .zplot_selnorm_threshold_quadrature(1.96, inputs$mean, inputs$sd, inputs$sei,
        selection, inputs$mean + .125, matrix(.1, inputs$S, inputs$K), rule, FALSE)
    }, threads = c(1L, 2L), rows = inputs$S,
      work_per_row = inputs$K * inputs$spec$n_bins * length(rule$nodes) * 2L)
  })
})
