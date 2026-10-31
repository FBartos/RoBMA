test_that("fitted Gaussian uncertainty is equivariant to effect-size units", {

  S <- 20L
  draw_at_scale <- function(scale) {
    set.seed(826)
    .evaluate.brma.true_effects_posterior.norm(
      mu_samples = matrix(.1 * scale, S, 2L),
      tau_within = matrix(c(.3, .5) * scale, S, 2L, byrow = TRUE),
      yi = c(.8, -.1) * scale,
      sei = c(.2, .4) * scale
    ) / scale
  }
  # Gaussian conjugacy is unchanged by a common change of effect-size units.
  # The SD formula must not square or multiply representable variances again.
  reference <- draw_at_scale(1)
  expect_equal(draw_at_scale(1e100), reference, tolerance = 1e-14)
  expect_equal(draw_at_scale(1e-100), reference, tolerance = 1e-14)
  expect_equal(draw_at_scale(1e154), reference, tolerance = 1e-14)
})

test_that("scalar Gaussian conditionals retain represented means and normal assignments", {

  S <- 20000L
  for (scale in c(1, 1e100, 1e-100, 1e154)) {
    set.seed(2302)
    normals <- matrix(stats::rnorm(S * 3L), S, 3L)
    expected_seed <- .Random.seed
    set.seed(2302)
    value <- .evaluate.brma.true_effects_posterior.norm(
      matrix(c(0, -.2, -.2) * scale, S, 3L, byrow = TRUE),
      matrix(c(1, 0, .7) * scale, S, 3L, byrow = TRUE),
      yi = c(8, 2, 2) * scale, sei = c(1, .7, 0) * scale,
      bias_offset = matrix(c(0, -.3, .3) * scale, S, 3L, byrow = TRUE))
    expect_equal(value[, 1L] / scale, 4 + normals[, 1L] / sqrt(2), tolerance = 1e-14)
    expect_equal(value[, 2L], rep(-.2 * scale, S), tolerance = 0)
    expect_equal(value[, 3L], rep(2 * scale - .3 * scale, S), tolerance = 0)
    expect_identical(.Random.seed, expected_seed)
  }
  expect_error(.evaluate.brma.true_effects.norm(matrix(0, 1L, 1L),
    matrix(1e-200, 1L, 1L), yi = 1e200, sei = 1, same_data = TRUE),
    "Gaussian conditional arithmetic is numerically unavailable.", fixed = TRUE)
  conditional_mean <- function(mu, yi, tau = 1, sei = 1) {
    .evaluate.brma.true_effects.norm(matrix(mu, 1L, 1L), matrix(tau, 1L, 1L),
      yi = yi, sei = sei, same_data = TRUE)
  }
  expect_error(conditional_mean(2^-1074, 2^-1074),
    "Gaussian conditional arithmetic is numerically unavailable.", fixed = TRUE)
  expect_equal(conditional_mean(0, 0), matrix(0, 1L, 1L), tolerance = 0)
  expect_equal(conditional_mean(2^-1022, -2^-1022), matrix(0, 1L, 1L), tolerance = 0)
  expect_equal(conditional_mean(2^-1074, 0, tau = 0), matrix(2^-1074, 1L, 1L), tolerance = 0)
  expect_equal(conditional_mean(0, 2^-1074, sei = 0), matrix(2^-1074, 1L, 1L), tolerance = 0)
  expect_equal(conditional_mean(2^-1074, 0, tau = 0, sei = 0),
    matrix(2^-1074, 1L, 1L), tolerance = 0)
})

test_that("one-estimate GLMM replication accepts scalar JAGS nuisance names", {

  binomial <- brma.glmm(ai = 2, ci = 1, n1i = 10, n2i = 12, measure = "OR",
    only_priors = TRUE)
  poisson <- brma.glmm(x1i = 2, x2i = 1, t1i = 10, t2i = 12, measure = "IRR",
    only_priors = TRUE)
  samples <- cbind(mu = c(-.2, .1, .3, .5), tau = .4,
    theta = c(-.5, 0, .5, 1), pi = c(.1, .2, .3, .4), phi = log(c(.1, .2, .3, .4)))
  effect <- samples[, "mu"] + samples[, "tau"] * samples[, "theta"]
  set.seed(827)
  expected_binomial <- cbind(
    stats::rbinom(4L, 10, stats::plogis(stats::qlogis(samples[, "pi"]) + effect / 2)),
    stats::rbinom(4L, 12, stats::plogis(stats::qlogis(samples[, "pi"]) - effect / 2))
  )
  set.seed(827)
  actual_binomial <- predict(binomial, type = "response", conditioning_depth = "estimate",
    as_measure = FALSE, quiet = TRUE, .posterior_samples = samples)
  expect_equal(unname(as.matrix(actual_binomial)), expected_binomial, tolerance = 0)

  set.seed(828)
  expected_poisson <- cbind(
    stats::rpois(4L, exp(samples[, "phi"] + effect / 2 + log(10))),
    stats::rpois(4L, exp(samples[, "phi"] - effect / 2 + log(12)))
  )
  set.seed(828)
  actual_poisson <- predict(poisson, type = "response", conditioning_depth = "estimate",
    as_measure = FALSE, quiet = TRUE, .posterior_samples = samples)
  expect_equal(unname(as.matrix(actual_poisson)), expected_poisson, tolerance = 0)
})
