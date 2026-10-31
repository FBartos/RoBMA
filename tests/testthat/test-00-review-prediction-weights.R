.prediction_weight_fixture <- function(weights = c(.25, 1.5, 4),
                                       selection = BayesTools::selection_model()) {

  bselmodel(yi = c(.8, -.1, .3), vi = c(.04, .09, .05), weights = weights,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .05,
      BayesTools::wf_fixed(c(1, .4)), model = selection))
}


.prediction_weight_state <- function(object, samples, type = "estimate", adjusted = FALSE,
                                     depth = "estimate") {

  context <- .predict_brma_context(object, newdata = NULL, V_new = NULL,
    type = type, conditioning_depth = depth, conditioning_depth_specified = TRUE,
    as_measure = TRUE, output_measure = NULL, transform = NULL, probs = c(.025, .975),
    bias_adjusted = adjusted, quiet = TRUE, conditional = FALSE,
    dots = list(.posterior_samples = samples))
  scale <- .predict_brma_scale_state(context)
  location <- .predict_brma_location_state(context, scale)
  list(context = context, scale = scale, location = location)
}


test_that("weighted fitted augmentation agrees with an independently normalized kernel", {

  object <- .prediction_weight_fixture()
  S <- 20000L
  samples <- cbind(mu = rep(.1, S), tau = .3)
  state <- .prediction_weight_state(object, samples)
  parts <- state$location$selection_parts
  data <- object$data$outcome
  means <- variances <- numeric(3L)
  for (i in 1:3) {
    publication_weight <- if (stats::pnorm(data$yi[i] / data$sei[i], lower.tail = FALSE) <= .05) 1 else .4
    kernel <- function(u) stats::dnorm(u, 0, .3) *
      stats::dnorm(data$yi[i], .1 + u, data$sei[i])^data$weights[i] * publication_weight
    mass <- stats::integrate(kernel, -Inf, Inf, rel.tol = 1e-11)$value
    mean <- stats::integrate(function(u) u * kernel(u), -Inf, Inf, rel.tol = 1e-11)$value / mass
    variances[i] <- stats::integrate(function(u) (u - mean)^2 * kernel(u),
      -Inf, Inf, rel.tol = 1e-11)$value / mass
    means[i] <- .1 + mean
  }
  actual_mean <- .predict_joint_selection_source_posterior(parts, data$yi, draw = FALSE)
  expect_lt(max(abs(actual_mean[1L, ] - means)), 1e-10)
  set.seed(13021)
  draws <- predict(object, type = "estimate", conditioning_depth = "estimate", quiet = TRUE,
    .posterior_samples = samples)
  expect_true(all(abs(colMeans(draws) - means) < 6 * sqrt(variances / S)))
  expect_true(all(abs(apply(draws, 2L, stats::var) - variances) <
    6 * variances * sqrt(2 / (S - 1))))
  expect_identical(attr(draws, "RoBMA_target")$latent_augmentation,
    list(policy = "normalized_powered_gaussian", likelihood_weights = data$weights))
  expect_equal(diag(parts$sampling_covariance_matrix), data$sei^2 / data$weights, tolerance = 0)
  response <- .predict_joint_selection_response_setup(state$context, state$location, state$scale)
  expect_equal(diag(response$sampling_covariance_matrix), data$sei^2, tolerance = 0)
})


test_that("weighted augmentation preserves branch offsets and raw response sampling", {

  bias <- BayesTools::prior_weightfunction("one-sided", .05, BayesTools::wf_fixed(c(1, .4)))
  for (direction in c("positive", "negative")) {
    sign <- if (direction == "positive") 1 else -1
    object <- RoBMA(yi = sign * c(.8, -.1, .3), vi = c(.04, .09, .05),
      weights = c(.25, 1.5, 4), measure = "GEN", effect_direction = direction,
      prior_unit_information_sd = 1, only_priors = TRUE,
      prior_effect = BayesTools::prior("normal", list(0, 1)), prior_effect_null = NULL,
      prior_heterogeneity = BayesTools::prior("normal", list(0, 1), list(0, Inf)),
      prior_heterogeneity_null = NULL,
      prior_bias = list(bias, BayesTools::prior_PET("normal", list(0, 1)),
        BayesTools::prior_PEESE("normal", list(0, 1))), prior_bias_null = BayesTools::prior_none())
    priors <- .create_fit_priors(object$data, object$priors)$bias
    per_branch <- 4000L
    indicator <- rep(seq_along(priors), each = per_branch)
    samples <- cbind(mu = sign * .1, tau = .3, bias_indicator = indicator,
      PET = .6 * vapply(priors, BayesTools::is.prior.PET, logical(1L))[indicator],
      PEESE = .8 * vapply(priors, BayesTools::is.prior.PEESE, logical(1L))[indicator],
      "omega[1]" = 1, "omega[2]" = .4)
    for (adjusted in c(FALSE, TRUE)) {
      state <- .prediction_weight_state(object, samples, adjusted = adjusted)
      parts <- state$location$selection_parts
      a <- object$data$outcome$sei^2 / object$data$outcome$weights
      expected_mean <- parts$latent_means +
        .3^2 / matrix(.3^2 + a, nrow(samples), 3L, byrow = TRUE) *
        (matrix(object$data$outcome$yi, nrow(samples), 3L, byrow = TRUE) - parts$means)
      variance <- .3^2 * a / (.3^2 + a)
      set.seed(13022)
      value <- predict(object, type = "estimate", conditioning_depth = "estimate",
        bias_adjusted = adjusted, quiet = TRUE, .posterior_samples = samples)
      for (branch in seq_along(priors)) {
        rows <- which(indicator == branch)
        expect_true(all(abs(colMeans(value[rows, , drop = FALSE]) - expected_mean[rows[1L], ]) <
          6 * sqrt(variance / per_branch)))
        expect_true(all(abs(apply(value[rows, , drop = FALSE], 2L, stats::var) - variance) <
          6 * variance * sqrt(2 / (per_branch - 1))))
      }
      response <- .predict_joint_selection_response_setup(state$context, state$location, state$scale)
      expect_equal(diag(response$sampling_covariance_matrix), object$data$outcome$sei^2, tolerance = 0)
    }
  }
  object <- .prediction_weight_fixture()
  samples <- cbind(mu = c(.1, .2, -.1), tau = .3)
  state <- .prediction_weight_state(object, samples, type = "response", adjusted = TRUE)
  set.seed(13023)
  truth <- .predict_brma_estimate_draws(state$context, state$location, state$scale)
  expected <- truth + matrix(stats::rnorm(9L), 3L, 3L) *
    matrix(object$data$outcome$sei, 3L, 3L, byrow = TRUE)
  expected_seed <- .Random.seed
  set.seed(13023)
  response <- predict(object, type = "response", conditioning_depth = "estimate",
    bias_adjusted = TRUE, quiet = TRUE, .posterior_samples = samples)
  expect_equal(unname(as.matrix(response)), expected, tolerance = 0)
  expect_identical(.Random.seed, expected_seed)
  expect_identical(attr(response, "RoBMA_target")$latent_augmentation$likelihood_weights,
    object$data$outcome$weights)
})


test_that("augmentation provenance follows reconstruction and represented variance", {

  omitted <- .prediction_weight_fixture(NULL)
  unit <- .prediction_weight_fixture(rep(1, 3L))
  samples <- cbind(mu = c(.1, .2), tau = .3)
  set.seed(13024)
  a <- predict(omitted, type = "estimate", conditioning_depth = "estimate", quiet = TRUE,
    .posterior_samples = samples)
  set.seed(13024)
  b <- predict(unit, type = "estimate", conditioning_depth = "estimate", quiet = TRUE,
    .posterior_samples = samples)
  expect_equal(as.numeric(a), as.numeric(b), tolerance = 0)
  expect_identical(attr(a, "RoBMA_target")$latent_augmentation,
    attr(b, "RoBMA_target")$latent_augmentation)
  expect_identical(attr(a, "RoBMA_target")$latent_augmentation$likelihood_weights, rep(1, 3L))
  zero <- predict(omitted, type = "estimate", conditioning_depth = "estimate", quiet = TRUE,
    .posterior_samples = cbind(mu = c(.1, .2), tau = 0))
  expect_equal(as.numeric(zero), as.numeric(matrix(c(.1, .2), 2L, 3L)), tolerance = 0)
  for (type in c("terms", "blup")) {
    value <- predict(omitted, type = type, quiet = TRUE, .posterior_samples = samples)
    expect_null(attr(value, "RoBMA_target")$latent_augmentation)
  }
  value <- predict(omitted, type = "estimate", conditioning_depth = "marginal", quiet = TRUE,
    .posterior_samples = samples)
  expect_null(attr(value, "RoBMA_target")$latent_augmentation)
  retained <- .prediction_weight_fixture(selection = BayesTools::selection_model(
    estimate_random_effects = "condition"))
  retained_samples <- cbind(samples, "theta[1]" = .2, "theta[2]" = -.1, "theta[3]" = .4)
  value <- predict(retained, type = "estimate", conditioning_depth = "estimate", quiet = TRUE,
    .posterior_samples = retained_samples)
  expect_null(attr(value, "RoBMA_target")$latent_augmentation)
  expect_null(.selection_postfit_target_metadata(omitted$data)$latent_augmentation)
  expect_null(.brma_mv_marglik_target_metadata(omitted)$latent_augmentation)
  tiny <- bselmodel(yi = c(.8, -.1, .3), vi = c(1e-100, .09, .05),
    weights = c(1e308, 1.5, 4), measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE)
  expect_error(.prediction_weight_state(tiny, samples),
    "Fitted Gaussian conditioning variance is numerically unavailable", fixed = TRUE)
})


test_that("retained sampling truth has no powered Gaussian augmentation marker", {

  samples <- cbind(mu = c(.1, .2), tau = .3,
    "theta[1]" = .4, "theta[2]" = -.2, "theta[3]" = .5,
    "sampling_z[1]" = .2, "sampling_z[2]" = -.3, "sampling_z[3]" = .1)
  omitted <- .prediction_weight_fixture(NULL, BayesTools::selection_model(
    known_sampling_variance = "condition"))
  unit <- .prediction_weight_fixture(rep(1, 3L), BayesTools::selection_model(
    known_sampling_variance = "condition"))
  sei <- omitted$data$outcome$sei
  y <- matrix(omitted$data$outcome$yi, 2L, 3L, byrow = TRUE)
  e0 <- matrix(c(.2, -.3, .1) * sei, 2L, 3L, byrow = TRUE)
  u0 <- matrix(c(.4, -.2, .5) * .3, 2L, 3L, byrow = TRUE)
  noise_weight <- matrix(sei^2 / (sei^2 + .3^2), 2L, 3L, byrow = TRUE)
  truth <- y - (e0 + noise_weight * (y - samples[, "mu"] - e0 - u0))
  for (type in c("estimate", "response")) for (adjusted in c(FALSE, TRUE)) {
    set.seed(13071)
    a <- predict(omitted, type = type, conditioning_depth = "estimate",
      bias_adjusted = adjusted, quiet = TRUE, .posterior_samples = samples)
    expected_seed <- .Random.seed
    set.seed(13071)
    b <- predict(unit, type = type, conditioning_depth = "estimate",
      bias_adjusted = adjusted, quiet = TRUE, .posterior_samples = samples)
    expect_equal(as.numeric(a), as.numeric(b), tolerance = 0)
    expect_identical(.Random.seed, expected_seed)
    expect_null(attr(a, "RoBMA_target")$latent_augmentation)
    expect_null(attr(b, "RoBMA_target")$latent_augmentation)
    if (type == "estimate") expect_equal(as.numeric(a), as.numeric(truth), tolerance = 1e-14)
    if (type == "response" && adjusted) {
      set.seed(13071)
      expected <- truth + matrix(stats::rnorm(6L), 2L, 3L) *
        matrix(sei, 2L, 3L, byrow = TRUE)
      expect_equal(as.numeric(a), as.numeric(expected), tolerance = 1e-14)
      expect_identical(.Random.seed, expected_seed)
    }
    integrated <- predict(.prediction_weight_fixture(NULL), type = type,
      conditioning_depth = "estimate", bias_adjusted = adjusted, quiet = TRUE,
      .posterior_samples = samples[, c("mu", "tau"), drop = FALSE])
    expect_identical(attr(integrated, "RoBMA_target")$latent_augmentation$policy,
      "normalized_powered_gaussian")
  }
})
