.estimate_selection_evidence_fixture <- function(mode, rule = "product",
                                                 direction = "positive") {

  dat <- data.frame(yi = c(.4, -.2, .3), sei = c(.8, 1.1, .6),
                    paper = rep("a", 3L))
  object <- bselmodel(
    yi = yi, sei = sei, data = dat, measure = "GEN",
    effect_direction = direction,
    prior_unit_information_sd = 1, only_priors = TRUE,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1),
                                  truncation = list(lower = 0)),
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(
        estimate_random_effects = mode, other_random_effects = "integrate",
        known_sampling_variance = "integrate", group = "paper", weight_rule = rule
      )
    )
  )
  samples <- rbind(c(.2, .4, .5, -1, 1.5), c(-.1, .7, -1, .2, .6))
  colnames(samples) <- c("mu", "tau", paste0("theta[", 1:3, "]"))
  if (mode == "integrate") samples <- samples[, c("mu", "tau"), drop = FALSE]
  context <- .iwmde_context_ensure_caches(structure(list(
    object = object, data = object$data, priors = object$priors,
    posterior_samples = samples,
    flat_prior_list = .create_fit_priors(object$data, object$priors),
    formula_fit = object$fit,
    formula_inputs = .iwmde_formula_inputs(object$data, object$priors),
    selection_spec = .iwmde_selection_spec(object$data, object$priors)
  ), class = "iwmde_context"))
  list(object = object, dat = dat, samples = samples, context = context)
}


test_that("estimate selection evidence uses the requested Gaussian kernel", {

  specifications <- expand.grid(mode = c("integrate", "condition"),
                                direction = c("positive", "negative"),
                                stringsAsFactors = FALSE)
  for (cell in seq_len(nrow(specifications))) {
    mode <- specifications$mode[cell]
    direction <- specifications$direction[cell]
    sign <- if (direction == "positive") 1 else -1
    fixture <- .estimate_selection_evidence_fixture(mode, direction = direction)
    object <- fixture$object
    samples <- fixture$samples
    dat <- fixture$dat
    setup <- .log_lik_posterior_setup(
      object$fit, samples, object$data, object$priors, "estimate", NULL
    )
    expected_mean <- matrix(samples[, "mu"], 2L, 3L)
    expected_sd <- matrix(dat$sei, 2L, 3L, byrow = TRUE)
    if (mode == "condition") {
      expected_mean <- expected_mean + samples[, 3:5] * samples[, "tau"]
    } else {
      expected_sd <- sqrt(expected_sd^2 + samples[, "tau"]^2)
    }
    dimnames(expected_mean) <- NULL
    # One-sided selection at p=.5 partitions outcomes by their declared direction.
    # This exact normal-CDF expression is independent of selection kernels.
    mass <- .2 + .8 * stats::pnorm(sign * expected_mean / expected_sd)
    expected <- stats::dnorm(matrix(dat$yi, 2L, 3L, byrow = TRUE),
                              expected_mean, expected_sd, log = TRUE) +
      matrix(log(ifelse(sign * dat$yi >= 0, 1, .2)), 2L, 3L, byrow = TRUE) - log(mass)
    expect_equal(unname(setup$mu), expected_mean, tolerance = 1e-14)
    expect_equal(setup$tau_within, matrix(samples[, "tau"], 2L, 3L),
                 tolerance = 0)
    expect_equal(.log_lik_estimate_from_setup(setup), expected, tolerance = 1e-12)
    expect_equal(.log_lik_estimate_sum_from_setup(setup), rowSums(expected),
                 tolerance = 1e-12)
    evaluated <- .log_lik_evaluated_setup(
      object$fit, object$data, object$priors, "estimate", NULL,
      mu_samples = matrix(samples[, "mu"], 2L, 3L),
      tau_within_samples = matrix(samples[, "tau"], 2L, 3L),
      tau_between_samples = NULL, posterior_samples = samples
    )
    expect_equal(.log_lik_estimate_sum_from_setup(evaluated), rowSums(expected),
                 tolerance = 1e-12)
    for (row in 1:2) {
      active <- .iwmde_active_setup(fixture$context, samples[row, ])
      parameters <- .iwmde_row_parameters(
        fixture$context, samples[row, ], active, state_scope = "local"
      )
      bridge <- .log_posterior(
        parameters = parameters, data = active$fit_data,
        is_scale = FALSE, is_multilevel = FALSE, is_weights = FALSE,
        is_known_v = FALSE, is_PET = FALSE, is_PEESE = FALSE,
        is_weightfunction = TRUE, effect_direction = direction, outcome_type = "norm",
        model_data = object$data
      )
      expect_equal(bridge, sum(expected[row, ]), tolerance = 1e-12)
    }
  }
})


test_that("conditional estimate q grids preserve latent priors and update their means", {

  fixture <- .estimate_selection_evidence_fixture("condition")
  context <- fixture$context
  samples <- fixture$samples
  values <- c(.15, .8)
  states <- .iwmde_row_states(context, 1:2, "tau", estimator = "q_grid_cmde")
  expect_identical(vapply(states, `[[`, character(1L), "state_scope"), rep("local", 2L))
  expect_true(all(vapply(states, function(state) "theta" %in% names(state$prior_list),
                        logical(1L))))
  expect_equal(lapply(states, function(state) state$parameters$theta),
               lapply(1:2, function(row) as.numeric(samples[row, 3:5])), tolerance = 0)
  actual <- .iwmde_log_q_grid(
    context, "tau", values, states, replacement = list(type = "scalar")
  )
  expected <- outer(values, 1:2, Vectorize(function(tau, row) {
    theta <- samples[row, 3:5]
    mean <- samples[row, "mu"] + tau * theta
    se <- fixture$dat$sei
    likelihood <- sum(stats::dnorm(fixture$dat$yi, mean, se, log = TRUE) +
      log(ifelse(fixture$dat$yi >= 0, 1, .2)) -
      log(.2 + .8 * stats::pnorm(mean / se)))
    likelihood + stats::dnorm(samples[row, "mu"], log = TRUE) +
      stats::dnorm(tau, log = TRUE) + log(2) + sum(stats::dnorm(theta, log = TRUE))
  }))
  expect_equal(actual, expected, tolerance = 1e-12)
})


test_that("conditional estimate deletion retains the original best-weight event", {

  fixture <- .estimate_selection_evidence_fixture("condition", "best")
  object <- fixture$object
  setup <- .log_lik_posterior_setup(
    object$fit, fixture$samples, object$data, object$priors, "estimate", NULL
  )
  # Delete estimates 1 and 3, retaining negative estimate 2. The event weight
  # is .2 exactly when both missing outcomes are also negative, otherwise 1.
  index <- c(1L, 3L)
  actual <- .selection_joint_deletion_loglik_from_setup(setup, list(index))[, 1L]
  mean <- setup$mu[, index, drop = FALSE]
  se <- matrix(fixture$dat$sei[index], 2L, 2L, byrow = TRUE)
  negative_mass <- stats::pnorm(-mean / se)
  normalizer <- 1 - .8 * negative_mass[, 1L] * negative_mass[, 2L]
  expected <- rowSums(stats::dnorm(
    matrix(fixture$dat$yi[index], 2L, 2L, byrow = TRUE), mean, se, log = TRUE
  )) - log(normalizer)
  expect_equal(actual, expected, tolerance = 1e-12)
})


test_that("Gaussian ensemble branches retain the fitted selection source representation", {

  diagonal <- c(.04, .09)
  loading <- c(.2, -.3)
  yi <- c(.1, -.2)
  object <- bselmodel.mv(
    yi = yi, V = known_v_factor(diagonal, matrix(loading, ncol = 1L)),
    data = data.frame(paper = c("a", "a")),
    random = NULL, measure = "GEN", prior_unit_information_sd = 1,
    effect_direction = "positive",
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(known_sampling_variance = "condition",
                                         group = "paper")
    ), only_priors = TRUE
  )
  samples <- matrix(c(.1, .4, -.2, .5), 1L,
                    dimnames = list(NULL, c("mu", paste0("sampling_z[", 1:3, "]"))))
  priors <- object$priors
  priors$outcome$bias <- BayesTools::prior_none()
  context <- list(object = object, data = object$data,
                  predictor_cache = new.env(parent = emptyenv()))
  # With no integrated random source, full sampling conditioning cancels
  # positive weights. Augmentation coordinates remain independent auxiliaries.
  actual <- .iwmde_log_lik_known_v_joint_sum_from_samples(
    context, samples, active_setup = list(priors = priors), unit = "estimate"
  )
  covariance <- diag(diagonal) + tcrossprod(loading)
  expected <- mvtnorm::dmvnorm(yi, rep(.1, 2L), covariance, log = TRUE)
  expect_equal(actual, expected, tolerance = 1e-12)
})


test_that("evaluated selected predictors retain estimate and cluster effects once", {

  object <- bselmodel(
    yi = c(.4, -.2), sei = c(.8, 1.1), cluster = c("a", "a"),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(estimate_random_effects = "condition")
    )
  )
  samples <- matrix(c(.2, .4, .3, .5, -1, 1.5), 1L,
    dimnames = list(NULL, c("mu", "tau", "rho", "gamma[1]", "theta[1]", "theta[2]")))
  within <- matrix(.4 * sqrt(.7), 1L, 2L)
  between <- matrix(.4 * sqrt(.3), 1L, 2L)
  expected_mean <- .2 + .5 * between + matrix(c(-1, 1.5), 1L) * within
  posterior <- .log_lik_posterior_setup(
    object$fit, samples, object$data, object$priors, "estimate", NULL
  )
  evaluated <- .log_lik_evaluated_setup(
    object$fit, object$data, object$priors, "estimate", NULL,
    mu_samples = matrix(.2, 1L, 2L), tau_within_samples = within,
    tau_between_samples = between, posterior_samples = samples
  )
  expected <- sum(stats::dnorm(c(.4, -.2), expected_mean, c(.8, 1.1), log = TRUE) +
    log(c(1, .2)) - log(.2 + .8 * stats::pnorm(expected_mean / c(.8, 1.1))))
  expect_equal(unname(posterior$mu), expected_mean, tolerance = 1e-14)
  expect_equal(unname(evaluated$mu), expected_mean, tolerance = 1e-14)
  expect_equal(.log_lik_estimate_sum_from_setup(evaluated), expected, tolerance = 1e-12)
})
