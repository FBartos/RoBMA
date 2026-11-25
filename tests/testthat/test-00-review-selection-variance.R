test_that("selected scalar and matrix moments return centered variance directly", {

  prior <- prior_weightfunction("one-sided", pnorm(-1), wf_fixed(c(1, .5)))
  context <- .selection_spec(list(outcome = list(bias = prior)), yi = 1e8,
    sei = 1e8, effect_direction = "positive", signed_data = FALSE)
  context$omega <- matrix(rep(c(1, .5), each = 2), 2, 2)
  context$alpha <- c(0, 0)
  context$phack_kind <- c(0L, 0L)
  context$kernel_mode <- rep(SELKERNEL_STEP, 2)
  means <- c(1e8 - 1, 1e8 + 2)
  a <- c(1, -2)
  mass <- pnorm(a, lower.tail = FALSE) + .5 * pnorm(a)
  shift <- .5 * dnorm(a) / mass
  expected <- 1 + .5 * a * dnorm(a) / mass - shift^2
  moments <- .selection_step_moments_matrix(matrix(means, 2), matrix(1, 2), 1e8, context)
  expect_equal(as.numeric(moments$variance), expected, tolerance = 1e-12)
  scalar <- .selection_step_moments_matrix(matrix(means[1], 1), matrix(1, 1), 1e8,
    BayesTools::selection_context_subset_rows(context, 1L))
  expect_equal(as.numeric(scalar$variance), expected[1], tolerance = 1e-12)
  summary <- .selection_step_weighted_summary(1e8, matrix(means, 2), matrix(1, 2),
    1e8, matrix(c(.25, .75), 2), context)
  offsets <- c(-1, 2) + shift
  center <- sum(c(.25, .75) * offsets)
  reference <- sum(c(.25, .75) * (expected + (offsets - center)^2))
  expect_lt(abs(summary$variance - reference), 1e-7)
})

test_that("outer sampling quadrature accumulates moments about its Gaussian location", {

  object <- bselmodel(yi = 1e8, sei = 1, measure = "GEN", effect_direction = "positive",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE,
    prior_effect = prior("normal", list(1e8, 1)), prior_heterogeneity = prior("point", list(1)),
    prior_bias = prior_weightfunction("one-sided", .5, wf_fixed(c(1, .2)),
      model = selection_model(known_sampling_variance = "condition")))
  samples <- cbind(mu = 1e8, tau = 1, theta = 0, sampling_z = 0)
  setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors,
    "estimate", NULL)
  # The unequal-weight STEP branch is retained. Its threshold is 0, over 1e8
  # Gaussian SDs from this location, so all omitted selection-tail effects are
  # below binary64 range. This is not active-selection translation invariance.
  result <- .selection_joint_estimate_targets(setup, c("mean", "variance"))
  expect_equal(as.numeric(result$mean), 1e8, tolerance = 0)
  expect_lt(abs(as.numeric(result$variance) - 2), 1e-7)
})
