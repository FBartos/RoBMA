test_that("retained sampling does not enter the scalar marginal selection CDF", {

  # With tau = 0, conditioned sampling, and positive selection weights, the
  # weights cancel and the marginal law is ordinary normal. Folding the
  # sampling variance into the scalar selected CDF incorrectly tilts that law.
  object <- bselmodel(
    yi = c(.2, .5), sei = c(.2, .3), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE,
    selection = BayesTools::selection_model(known_sampling_variance = "condition")
  )
  expected <- matrix(c(.3, .7), 1L)
  testthat::local_mocked_bindings(
    .cdf_lik_estimate.brma = function(object) expected,
    .package = "RoBMA"
  )
  expect_equal(unname(.cdf.brma(object, conditioning_depth = "estimate")), expected)
  expect_error(
    .cdf.brma(object, conditioning_depth = "marginal"),
    "Joint-selection CDF evaluation is unavailable at this conditioning depth.",
    fixed = TRUE
  )
})

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
