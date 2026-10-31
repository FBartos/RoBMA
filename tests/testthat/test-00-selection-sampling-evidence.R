.sampling_selection_evidence_fixture <- function(estimate = "integrate") {

  dat <- data.frame(yi = c(.4, -.2), sei = c(.8, 1.1), paper = c("a", "a"))
  object <- bselmodel(
    yi = yi, sei = sei, data = dat, measure = "GEN", effect_direction = "positive",
    prior_unit_information_sd = 1, only_priors = TRUE,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1), truncation = list(lower = 0)),
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(estimate_random_effects = estimate,
        known_sampling_variance = "condition", group = "paper")
    ),
    selection_control = set_selection_likelihood_control(relative_tolerance = 1e-6)
  )
  samples <- rbind(c(.2, .4, .5, -1, .3, -.4), c(-.1, .7, -.2, .6, -.5, .2))
  colnames(samples) <- c("mu", "tau", paste0("theta[", 1:2, "]"), paste0("sampling_z[", 1:2, "]"))
  context <- .iwmde_context_ensure_caches(structure(list(
    object = object, data = object$data, priors = object$priors, posterior_samples = samples,
    flat_prior_list = .create_fit_priors(object$data, object$priors), formula_fit = object$fit,
    formula_inputs = .iwmde_formula_inputs(object$data, object$priors),
    selection_spec = .iwmde_selection_spec(object$data, object$priors)
  ), class = "iwmde_context"))
  list(object = object, dat = dat, samples = samples, context = context)
}


test_that("sampling-conditioned joint evidence replays Gaussian augmentation at each candidate", {

  fixture <- .sampling_selection_evidence_fixture()
  object <- fixture$object
  samples <- fixture$samples
  yi <- fixture$dat$yi
  vi <- fixture$dat$sei^2
  oracle <- function(row, tau = row[["tau"]]) {

    e0 <- sqrt(vi) * row[paste0("sampling_z[", 1:2, "]")]
    u0 <- tau * row[paste0("theta[", 1:2, "]")]
    e <- e0 + vi / (vi + tau^2) * (yi - row[["mu"]] - u0 - e0)
    mass <- .2 + .8 * stats::pnorm((row[["mu"]] + e) / tau)
    sum(stats::dnorm(yi, row[["mu"]], sqrt(vi + tau^2), log = TRUE) +
          log(ifelse(yi >= 0, 1, .2)) - log(mass))
  }
  expected <- apply(samples, 1L, oracle)
  setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
  expect_equal(unname(setup$mu), matrix(samples[, "mu"], 2L, 2L), tolerance = 0)
  expect_equal(.log_lik_estimate_sum_from_setup(setup), expected, tolerance = 1e-12)
  for (row in 1:2) {
    active <- .iwmde_active_setup(fixture$context, samples[row, ])
    bridge <- .log_posterior(
      parameters = list(), data = active$fit_data, is_scale = FALSE,
      is_multilevel = FALSE, is_weights = FALSE, is_known_v = FALSE,
      is_PET = FALSE, is_PEESE = FALSE, is_weightfunction = TRUE,
      effect_direction = "positive", outcome_type = "norm", model_data = object$data,
      bridge_context = list(nodes = samples[row, ]), selection_fit = object$fit,
      selection_priors = object$priors
    )
    expect_equal(bridge, expected[row], tolerance = 1e-12)
  }
  states <- .iwmde_row_states(fixture$context, 1:2, "tau", estimator = "q_grid_cmde")
  expect_true(all(vapply(states, function(state) all(c("theta", "sampling_z") %in% names(state$prior_list)), logical(1L))))
  values <- c(.15, .8)
  grid <- .iwmde_log_q_grid(fixture$context, "tau", values, states, replacement = list(type = "scalar"))
  reference <- outer(values, 1:2, Vectorize(function(tau, row) {
    oracle(samples[row, ], tau) + stats::dnorm(samples[row, "mu"], log = TRUE) +
      stats::dnorm(tau, log = TRUE) + log(2) + sum(stats::dnorm(samples[row, 3:6], log = TRUE))
  }))
  expect_equal(grid, reference, tolerance = 1e-12)
})


test_that("sampling-conditioned row deletion integrates the deleted sampling error", {

  fixture <- .sampling_selection_evidence_fixture()
  object <- fixture$object
  setup <- .log_lik_posterior_setup(object$fit, fixture$samples, object$data, object$priors, "estimate", NULL)
  actual <- .log_lik_estimate_from_setup(setup)
  # Independent product sources reduce to scalar integrals over original
  # sampling errors, regardless of the fitted auxiliary coordinates.
  expected <- matrix(NA_real_, 2L, 2L)
  for (draw in 1:2) for (row in 1:2) {
    mu <- fixture$samples[draw, "mu"]
    tau <- fixture$samples[draw, "tau"]
    y <- fixture$dat$yi[row]
    se <- fixture$dat$sei[row]
    density <- stats::integrate(function(e) {
      stats::dnorm(e, 0, se) * stats::dnorm(y - mu - e, 0, tau) /
        (.2 + .8 * stats::pnorm((mu + e) / tau))
    }, -Inf, Inf, rel.tol = 1e-10)$value * ifelse(y >= 0, 1, .2)
    expected[draw, row] <- log(density)
  }
  expect_equal(actual, expected, tolerance = 2e-6)
})

test_that("independent candidate effects retain correlated sampling-error conditioning", {

  fixture <- .sampling_selection_evidence_fixture()
  setup <- .log_lik_posterior_setup(fixture$object$fit, fixture$samples,
    fixture$object$data, fixture$object$priors, "estimate", NULL)
  plan <- .data_selection_execution_plan(setup$data)
  plan$row_blocks <- list(1:2)
  plan$factor_ranks <- 0L
  attr(setup$data, "selection_execution_plan") <- plan
  V <- matrix(c(.64, .3, .3, 1.21), 2L)
  errors <- rbind(c(.2, -.3), c(-.4, .1))
  variances <- array(0, c(2L, 2L, 2L))
  for (draw in 1:2) diag(variances[draw, , ]) <- fixture$samples[draw, "tau"]^2
  state <- list(sampling_covariance = V, e = errors,
    integrated_covariance = variances, total_covariance = variances,
    baseline_mu = matrix(fixture$samples[, "mu"], 2L, 2L))
  selection <- .selection_joint_signed_context(setup, setup$yi)
  actual <- .selection_conditioned_sampling_independent_targets(setup, state,
    selection, c("log_density", "cdf"))
  expected <- lapply(actual, function(x) matrix(0, 2L, 2L))
  for (draw in 1:2) for (row in 1:2) {
    other <- 3L - row
    mu <- fixture$samples[draw, "mu"]
    tau <- fixture$samples[draw, "tau"]
    y <- setup$yi[row]
    conditional_mean <- V[row, other] / V[other, other] * errors[draw, other]
    conditional_sd <- sqrt(V[row, row] - V[row, other]^2 / V[other, other])
    expected$log_density[draw, row] <- log(stats::integrate(function(e) {
      stats::dnorm(e, conditional_mean, conditional_sd) *
        stats::dnorm(y, mu + e, tau) * ifelse(y >= 0, 1, .2) /
        (.2 + .8 * stats::pnorm((mu + e) / tau))
    }, -Inf, Inf, rel.tol = 1e-10)$value)
    expected$cdf[draw, row] <- stats::integrate(function(e) {
      numerator <- if (y < 0) .2 * stats::pnorm((y - mu - e) / tau) else
        stats::pnorm((y - mu - e) / tau) - .8 * stats::pnorm((-mu - e) / tau)
      stats::dnorm(e, conditional_mean, conditional_sd) * numerator /
        (.2 + .8 * stats::pnorm((mu + e) / tau))
    }, -Inf, Inf, rel.tol = 1e-10)$value
  }
  expect_equal(actual, expected, tolerance = 2e-6)
  setup$effect_direction <- "negative"
  setup$yi <- -setup$yi
  state$baseline_mu <- -state$baseline_mu
  state$e <- -state$e
  mirrored <- .selection_conditioned_sampling_independent_targets(setup, state,
    selection, c("log_density", "cdf"))
  expect_equal(mirrored$log_density, expected$log_density, tolerance = 2e-6)
  expect_equal(mirrored$cdf, 1 - expected$cdf, tolerance = 2e-6)
})

test_that("sampling deletion quadrature preserves narrow upper-tail intervals", {

  fixture <- .sampling_selection_evidence_fixture()
  fixture$samples[, "mu"] <- -6.4
  setup <- .log_lik_posterior_setup(fixture$object$fit, fixture$samples,
    fixture$object$data, fixture$object$priors, "estimate", NULL)
  plan <- .data_selection_execution_plan(setup$data)
  # Coarser blocking changes no probability law for these independent sources.
  plan$row_blocks <- list(1:2)
  plan$factor_ranks <- 0L
  plan$design_keys <- plan$design_keys[1L]
  attr(setup$data, "selection_execution_plan") <- plan
  testthat::local_mocked_bindings(
    .selection_conditioned_sampling_independent_targets = function(...) NULL,
    .package = "RoBMA"
  )
  actual <- .selection_conditioned_sampling_estimate_targets(setup, "log_density")$log_density
  reference <- matrix(0, 2L, 2L)
  for (draw in 1:2) for (row in 1:2) {
    mu <- fixture$samples[draw, "mu"]
    tau <- fixture$samples[draw, "tau"]
    se <- fixture$dat$sei[row]
    y <- fixture$dat$yi[row]
    total <- se^2 + tau^2
    conditional_mean <- mu + se^2 / total * (y - mu)
    conditional_sd <- sqrt(se^2 * tau^2 / total)
    reciprocal <- stats::integrate(function(z) {
      stats::dnorm(z) / (.2 + .8 * stats::pnorm((conditional_mean + conditional_sd * z) / tau))
    }, -Inf, Inf, rel.tol = 1e-10)$value
    reference[draw, row] <- stats::dnorm(y, mu, sqrt(total), log = TRUE) +
      log(ifelse(y >= 0, 1, .2)) + log(reciprocal)
  }
  expect_equal(unname(actual), reference, tolerance = 2e-6)
})


test_that("independent conditioned-sampling CDFs and moments match scalar Gaussian identities", {

  fixture <- .sampling_selection_evidence_fixture()
  object <- fixture$object
  setup <- .log_lik_posterior_setup(object$fit, fixture$samples, object$data, object$priors, "estimate", NULL)
  actual <- .selection_joint_estimate_targets(setup, c("cdf", "log_lower", "log_upper", "mean", "variance"))
  for (draw in 1:2) for (row in 1:2) {
    mu <- fixture$samples[draw, "mu"]
    tau <- fixture$samples[draw, "tau"]
    y <- fixture$dat$yi[row]
    se <- fixture$dat$sei[row]
    reference <- vapply(c("cdf", "mean", "second"), function(component) {
      stats::integrate(function(e) {
        m <- mu + e
        mass <- .2 + .8 * stats::pnorm(m / tau)
        value <- switch(component,
          cdf = if (y < 0) .2 * stats::pnorm((y - m) / tau) / mass else
            (stats::pnorm((y - m) / tau) - .8 * stats::pnorm(-m / tau)) / mass,
          mean = m + .8 * tau * stats::dnorm(m / tau) / mass,
          second = m^2 + tau^2 + .8 * m * tau * stats::dnorm(m / tau) / mass)
        value * stats::dnorm(e, 0, se)
      }, -Inf, Inf, rel.tol = 1e-10)$value
    }, numeric(1L))
    expect_equal(actual$cdf[draw, row], unname(reference["cdf"]), tolerance = 2e-6)
    expect_equal(actual$log_lower[draw, row], log(unname(reference["cdf"])), tolerance = 2e-6)
    expect_equal(actual$log_upper[draw, row], log1p(-unname(reference["cdf"])), tolerance = 2e-6)
    expect_equal(actual$mean[draw, row], unname(reference["mean"]), tolerance = 2e-6)
    expect_equal(actual$variance[draw, row], unname(reference["second"] - reference["mean"]^2), tolerance = 2e-6)
  }
})


test_that("symmetric conditioned-sampling selection preserves a zero mean", {

  object <- bselmodel(yi = c(0, 0), sei = c(.8, 1.1), measure = "GEN",
    effect_direction = "negative", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction("two-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(known_sampling_variance = "condition")))
  samples <- matrix(c(0, .5, .2, -.3, .1, -.1), 1L,
    dimnames = list(NULL, c("mu", "tau", "theta[1]", "theta[2]", "sampling_z[1]", "sampling_z[2]")))
  setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
  actual <- .selection_joint_estimate_targets(setup, c("cdf", "mean"))
  expect_equal(actual$cdf, matrix(.5, 1L, 2L), tolerance = 1e-12)
  expect_equal(actual$mean, matrix(0, 1L, 2L), tolerance = 1e-12)
})


test_that("conditioned-sampling log tails stay finite after probability underflow", {

  object <- bselmodel(yi = c(-100, .2), sei = c(.8, 1.1), measure = "GEN",
    effect_direction = "positive", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(known_sampling_variance = "condition")))
  samples <- matrix(c(0, .5, .2, -.3, .1, -.1), 1L,
    dimnames = list(NULL, c("mu", "tau", "theta[1]", "theta[2]", "sampling_z[1]", "sampling_z[2]")))
  setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
  actual <- .selection_joint_estimate_targets(setup, c("log_density", "log_lower"))
  # At this negative tail the acceptance mass is .2 up to a Gaussian tail
  # below binary64 range, cancelling the observed .2 weight. Base R's
  # log-normal tail is an independent finite reference (ordinary CDF is zero).
  sd <- sqrt(.8^2 + .5^2)
  expect_equal(actual$log_density[1L, 1L], stats::dnorm(-100, 0, sd, log = TRUE), tolerance = 1e-9)
  expect_equal(actual$log_lower[1L, 1L], stats::pnorm(-100, 0, sd, log.p = TRUE), tolerance = 1e-9)
})


test_that("all-conditioned evidence cancels weights without using sampling auxiliaries", {

  fixture <- .sampling_selection_evidence_fixture("condition")
  object <- fixture$object
  samples <- fixture$samples
  setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
  means <- matrix(samples[, "mu"], 2L, 2L)
  sd <- sqrt(matrix(fixture$dat$sei^2, 2L, 2L, byrow = TRUE) + samples[, "tau"]^2)
  expected <- stats::dnorm(matrix(fixture$dat$yi, 2L, 2L, byrow = TRUE), means, sd, log = TRUE)
  expect_equal(.log_lik_estimate_sum_from_setup(setup), rowSums(expected), tolerance = 1e-12)
  y <- matrix(fixture$dat$yi, 2L, 2L, byrow = TRUE)
  vi <- matrix(fixture$dat$sei^2, 2L, 2L, byrow = TRUE)
  e0 <- sweep(samples[, paste0("sampling_z[", 1:2, "]")], 2L, fixture$dat$sei, "*")
  u0 <- samples[, paste0("theta[", 1:2, "]")] * samples[, "tau"]
  retained <- means + u0 + samples[, "tau"]^2 / (vi + samples[, "tau"]^2) *
    (y - means - u0 - e0)
  deletion <- stats::dnorm(y, retained, sqrt(vi), log = TRUE)
  expect_equal(.log_lik_estimate_from_setup(setup), unname(deletion), tolerance = 1e-12)
})


test_that("prepared scalar bridge evaluation preserves augmentation and direction", {

  for (estimate in c("integrate", "condition")) for (direction in c("positive", "negative")) {
    object <- bselmodel(yi = c(.4, -.2), sei = c(.8, 1.1), measure = "GEN",
      effect_direction = direction, prior_unit_information_sd = 1, only_priors = TRUE,
      prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
        BayesTools::wf_fixed(c(1, .2)), model = BayesTools::selection_model(
          estimate_random_effects = estimate, known_sampling_variance = "condition")))
    samples <- matrix(c(.2, .4, .5, -1, .3, -.4), 1L,
      dimnames = list(NULL, c("mu", "tau", "theta[1]", "theta[2]", "sampling_z[1]", "sampling_z[2]")))
    setup <- .log_lik_posterior_setup(object$fit, samples, object$data,
      object$priors, "estimate", NULL)
    fit_data <- .marglik_add_selection_bridge_data(.create_fit_data(object$data, object$priors),
      object$priors, direction, object$data)
    parameters <- list(mu = .2, tau = .4, theta = c(.5, -1), sampling_z = c(.3, -.4))
    plan <- .marglik_conditioned_sampling_plan(object$data)
    expect_false(is.null(plan))
    actual <- .marglik_conditioned_sampling_independent(parameters, fit_data, plan,
      FALSE, FALSE, FALSE, direction, NULL)
    sign <- if (direction == "negative") -1 else 1
    y <- c(.4, -.2)
    vi <- c(.8, 1.1)^2
    expected <- sum(stats::dnorm(y, .2, sqrt(vi + .4^2), log = TRUE))
    if (estimate == "integrate") {
      e0 <- sign * sqrt(vi) * parameters$sampling_z
      e <- e0 + vi / (vi + .4^2) * (y - .2 - e0 - .4 * parameters$theta)
      acceptance <- .2 + .8 * stats::pnorm(sign * (.2 + e) / .4)
      expected <- expected + sum(log(ifelse(sign * y >= 0, 1, .2)) - log(acceptance))
    }
    expect_equal(actual, expected, tolerance = 1e-12)
    expect_equal(actual, .selection_joint_loglik_from_setup(setup), tolerance = 1e-12)
    # Prepared formula outputs may have a separate mean and SD per row.
    parameters$mu <- c(.1, -.3)
    parameters$tau <- NULL
    parameters$log_tau <- log(c(.25, .7))
    actual_scale <- .marglik_conditioned_sampling_independent(parameters, fit_data, plan,
      TRUE, FALSE, FALSE, direction, NULL)
    expected_scale <- stats::dnorm(y, parameters$mu, sqrt(vi + c(.25, .7)^2), log = TRUE)
    if (estimate == "integrate") {
      e0 <- sign * sqrt(vi) * parameters$sampling_z
      e <- e0 + vi / (vi + c(.25, .7)^2) *
        (y - parameters$mu - e0 - c(.25, .7) * parameters$theta)
      acceptance <- .2 + .8 * stats::pnorm(sign * (parameters$mu + e) / c(.25, .7))
      expected_scale <- expected_scale + log(ifelse(sign * y >= 0, 1, .2)) - log(acceptance)
    }
    expect_equal(actual_scale, sum(expected_scale), tolerance = 1e-12)
  }
})


test_that("singleton sampling deletion distinguishes declared retained random context", {

  keys <- lapply(c("integrate", "condition"), function(estimate) {
    object <- .sampling_selection_evidence_fixture(estimate)$object
    key <- .current_predictive_target_key(object, "estimate")
    expect_identical(key$retained_context, "remaining_data")
    key
  })
  expect_false(identical(keys[[1L]]$retained_random_context, keys[[2L]]$retained_random_context))
  objects <- lapply(keys, function(key) {
    structure(list(), class = "psis_loo", RoBMA_target = key)
  })
  expect_error(.check_loo_compare_targets(objects), "different declared retained random contexts", fixed = TRUE)
  object <- bselmodel.mv(yi = c(.4, -.2), V = matrix(c(.64, .2, .2, 1.21), 2L),
    random = ~ 1 | estimate, data = data.frame(estimate = c("a", "b")),
    measure = "GEN", effect_direction = "positive", prior_unit_information_sd = 1,
    only_priors = TRUE, selection = selection_model(known_sampling_variance = "condition", group = "estimate"))
  expect_identical(.current_predictive_target_key(object, "estimate")$retained_context,
    "remaining_data_and_remaining_sampling_errors")
  attr(objects[[2L]], "RoBMA_target")$retained_context <- "remaining_data_and_remaining_sampling_errors"
  expect_error(.check_loo_compare_targets(objects),
    "LOO/WAIC objects with different data, unit, or retained-context targets cannot be compared.",
    fixed = TRUE)
})


test_that("independent deletion broadcasts a fixed heterogeneity scale", {

  fixture <- .sampling_selection_evidence_fixture("integrate")
  object <- fixture$object
  setup <- .log_lik_posterior_setup(object$fit, fixture$samples, object$data,
    object$priors, "estimate", NULL)
  setup$tau_within <- matrix(.4, 1L, setup$K)
  selection <- .selection_joint_signed_context(setup, setup$yi)
  actual <- .selection_conditioned_sampling_independent_targets(setup, NULL,
    selection, "log_density")$log_density
  expected <- matrix(NA_real_, setup$S, setup$K)
  for (draw in seq_len(setup$S)) for (row in seq_len(setup$K)) {
    mu <- setup$mu[draw, row]
    y <- setup$yi[row]
    se <- fixture$dat$sei[row]
    expected[draw, row] <- log(stats::integrate(function(e) {
      dnorm(e, 0, se) * dnorm(y, mu + e, .4) * ifelse(y >= 0, 1, .2) /
        (.2 + .8 * pnorm((mu + e) / .4))
    }, -Inf, Inf, rel.tol = 1e-10)$value)
  }
  expect_equal(actual, unname(expected), tolerance = 2e-6)
})

test_that("the deleted-row conditional mean keeps both variance shares", {

  # mean + variance / total * (observed - mean) loses other_variance / total *
  # mean once the share of X rounds to one: 2^80 with a variance ratio 2^-60
  # gave 0 instead of 2^20.
  expect_identical(.gaussian_summand_conditional_mean(2^80, 1, 2^-60, 0), 2^20)
  expect_identical(.gaussian_summand_conditional_mean(0, 2^-60, 1, 2^80), 2^20)
  expect_identical(.gaussian_summand_conditional_mean(.3, .5, 0, -1.2), -1.2)
  expect_identical(.gaussian_summand_conditional_mean(.3, 0, .5, -1.2), .3)

  # Mixed scales (locations 2^-95..2^95, variances 2^-190..2^190) against the
  # product form (other * mean + variance * observed) / total, relative to the
  # scale of the two terms.
  grid <- withr::with_seed(20261021, {
    n <- 2000L
    magnitude <- function(range) stats::runif(n, 1, 2) * 2^stats::runif(n, -range, range)
    sign <- function() sample(c(-1, 1), n, replace = TRUE)
    list(mean = sign() * magnitude(95), observed = sign() * magnitude(95),
         variance = magnitude(190), other = magnitude(190))
  }, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  total <- grid$variance + grid$other
  reference <- (grid$other * grid$mean + grid$variance * grid$observed) / total
  scale <- abs(grid$other / total * grid$mean) + abs(grid$variance / total * grid$observed)
  actual <- .gaussian_summand_conditional_mean(grid$mean, grid$variance, grid$other, grid$observed)
  expect_true(all(abs(actual - reference) <= 1e-14 * scale))
  # The grid contains rounded-to-one shares whose other term is not negligible.
  expect_true(any(grid$variance / total == 1 &
                    abs(grid$other / total * grid$mean) > 1e-14 * scale))

  # The scalar deletion route centres its quadrature there: the transition rule
  # receives 2^20 + y for mean 2^80 and integrated variance 2^-60 sei^2.
  fixture <- .sampling_selection_evidence_fixture("integrate")
  object <- fixture$object
  setup <- .log_lik_posterior_setup(object$fit, fixture$samples, object$data,
    object$priors, "estimate", NULL)
  setup$mu[] <- 2^80
  setup$tau_within <- matrix(fixture$dat$sei * 2^-30, 1L, setup$K)
  selection <- .selection_joint_signed_context(setup, setup$yi)
  original <- .selection_scalar_transition_rule_resolves
  centres <- list()
  testthat::local_mocked_bindings(
    .selection_scalar_transition_rule_resolves = function(nodes, mean, residual_sd, coefficient,
                                                          sei, ...) {
      if (all(mean < 2^79)) centres[[length(centres) + 1L]] <<- list(mean = mean, sei = sei)
      original(nodes, mean, residual_sd, coefficient, sei, ...)
    }, .package = "RoBMA")
  .selection_conditioned_sampling_independent_targets(setup, NULL, selection,
    c("log_density", "mean"))
  rows <- vapply(centres, function(centre) match(centre$sei, fixture$dat$sei), integer(1L))
  expect_setequal(rows, 1:2)
  for (i in seq_along(centres)) {
    expect_equal(centres[[i]]$mean, rep(2^20 + setup$yi[rows[i]], length(centres[[i]]$mean)),
                 tolerance = 1e-15)
  }
})

test_that("multi-row sampling deletion keeps both variance shares of its conditional mean", {

  # Deleting rows 1:2 conditions their sampling errors on the retained error
  # 2^80 (conditional mean .3 * 2^80 each) under an integrated random covariance
  # of 2^-60. Its share times that mean (about 2^18) is what
  # m_e + gain (y - m_e) loses once the gain rounds to the identity.
  dat <- data.frame(yi = c(.4, -.2, .9), sei = c(.8, 1.1, 1), paper = "a")
  object <- bselmodel(
    yi = yi, sei = sei, data = dat, measure = "GEN", effect_direction = "positive",
    prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(estimate_random_effects = "integrate",
        known_sampling_variance = "condition", group = "paper")))
  samples <- matrix(c(.1, .5, rep(0, 6L)), 1L, dimnames = list(NULL,
    c("mu", "tau", paste0("theta[", 1:3, "]"), paste0("sampling_z[", 1:3, "]"))))
  setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
  V <- matrix(.3, 3L, 3L) + diag(c(.34, .91, .7))
  state <- list(e = matrix(c(0, 0, 2^80), 1L), baseline_mu = matrix(0, 1L, 3L),
    sampling_covariance = V, integrated_covariance = array(diag(2^-60, 3L), c(1L, 3L, 3L)),
    log_normalizer = matrix(0, 1L, length(.data_selection_execution_plan(setup$data)$row_blocks)))
  centres <- list()
  testthat::local_mocked_bindings(
    .selection_conditioned_sampling_state = function(...) state,
    .selection_conditioned_sampling_normalizer = function(setup, means) {
      # Only the draws given the deleted outcomes lie below the 2^80 scale.
      if (all(abs(means[, 1:2]) < 2^60)) centres[[length(centres) + 1L]] <<- means[, 1:2, drop = FALSE]
      list(log_mass = rep(0, nrow(means)))
    }, .package = "RoBMA")
  .selection_conditioned_sampling_deletion_loglik(setup, list(1:2))
  sampling <- .selection_deleted_gaussian_block(V, state$e[1L, ], 1:2)
  total <- sampling$covariance + diag(2^-60, 2L)
  expected <- as.vector(diag(2^-60, 2L) %*% solve(total, sampling$mean) +
    sampling$covariance %*% solve(total, dat$yi[1:2]))
  expect_gt(length(centres), 0L)
  for (centre in centres) expect_equal(colMeans(centre), expected, tolerance = 1e-12)
})

test_that("native scalar deletion guards its own returned conditional rule", {

  evaluate <- function(orders, variance = .16, mode = 1L, omega = c(1, .2)) {
    rule <- .selection_joint_cluster_quadrature_rules(orders)
    .Call("RoBMA_selnorm_sampling_deletion_loglik_batch", 0, matrix(.2), .64,
      matrix(variance), matrix(.64 + variance), .7, matrix(omega, 1L),
      c(0, -Inf), c(Inf, 0), 1L, mode, FALSE,
      as.double(rule$nodes), as.double(rule$log_weights), as.integer(rule$orders), 1,
      PACKAGE = "RoBMA")
  }
  expect_true(is.na(evaluate(c(3L, 7L))[1L]))
  actual <- evaluate(c(3L, 7L, 15L))[1L]
  conditional_mean <- .2 + .64 / .8 * (0 - .2)
  conditional_sd <- sqrt(.64 * .16 / .8)
  expected <- stats::dnorm(0, .2, sqrt(.8), log = TRUE) +
    log(stats::integrate(function(u) stats::dnorm(u) /
      (.2 + .8 * stats::pnorm((conditional_mean + conditional_sd * u) / .4)),
      -Inf, Inf, rel.tol = 1e-11)$value)
  expect_equal(actual, expected, tolerance = 1e-4)
  expect_equal(evaluate(c(3L, 7L), variance = 0)[1L], stats::dnorm(0, .2, .8, log = TRUE))
  expect_equal(evaluate(c(3L, 7L), mode = 0L)[1L], stats::dnorm(0, .2, sqrt(.8), log = TRUE))
  expect_equal(evaluate(c(3L, 7L), omega = c(1, 1))[1L], stats::dnorm(0, .2, sqrt(.8), log = TRUE))
})

test_that("native and dependent-block deletion centre their integrals on both variance shares", {

  # The native batch: mean 2^79, sampling variance 1, integrated variance 2^-60 and
  # selection SE 2^19 put the share-form centre 2^-60 * 2^79 = 2^19 on the cutoff
  # z = 1, which order 7 cannot resolve; the lost-share form centred at y = 0, far
  # outside the rule. The Gaussian term dominates the value itself.
  native <- function(orders) {
    rule <- .selection_joint_cluster_quadrature_rules(orders)
    .Call("RoBMA_selnorm_sampling_deletion_loglik_batch", 0, matrix(2^79), 1,
      matrix(2^-60), matrix(1 + 2^-60), 2^19, matrix(c(1, .2), 1L),
      c(1, -Inf), c(Inf, 1), 1L, 1L, FALSE,
      as.double(rule$nodes), as.double(rule$log_weights), as.integer(rule$orders), 1,
      PACKAGE = "RoBMA")[1L]
  }
  expect_true(is.na(native(c(3L, 7L))))
  expect_true(is.finite(native(c(3L, 7L, 15L))))

  # The dependent-block route deletes row 1 of a three-row block whose retained
  # sampling error is 2^80 under an integrated covariance of 2^-60: its reciprocal
  # integral is centred on .gaussian_summand_conditional_mean(), about 5.0e5, not
  # on the lost-share value 0.
  dat <- data.frame(yi = c(.4, -.2, .9), sei = c(.8, 1.1, 1), paper = "a")
  object <- bselmodel(
    yi = yi, sei = sei, data = dat, measure = "GEN", effect_direction = "positive",
    prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(estimate_random_effects = "integrate",
        known_sampling_variance = "condition", group = "paper")))
  samples <- matrix(c(.1, .5, rep(0, 6L)), 1L, dimnames = list(NULL,
    c("mu", "tau", paste0("theta[", 1:3, "]"), paste0("sampling_z[", 1:3, "]"))))
  setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
  V <- matrix(.3, 3L, 3L) + diag(c(.34, .91, .7))
  state <- list(e = matrix(c(0, 0, 2^80), 1L), baseline_mu = matrix(0, 1L, 3L),
    sampling_covariance = V, integrated_covariance = array(diag(2^-60, 3L), c(1L, 3L, 3L)),
    total_covariance = array(V + diag(2^-60, 3L), c(1L, 3L, 3L)),
    log_normalizer = matrix(0, 1L, length(.data_selection_execution_plan(setup$data)$row_blocks)))
  centres <- numeric()
  testthat::local_mocked_bindings(
    .selection_conditioned_sampling_state = function(...) state,
    .selection_conditioned_sampling_independent_targets = function(...) NULL,
    .selection_conditioned_sampling_normalizer = function(setup, means) {
      if (all(abs(means[, 1L]) < 2^60)) centres <<- c(centres, means[, 1L])
      list(log_mass = rep(0, nrow(means)))
    }, .package = "RoBMA")
  .selection_conditioned_sampling_estimate_targets_chunk(setup, "log_density")
  sampling <- .selection_deleted_gaussian_coordinates(V, state$e[1L, ])
  centre <- .gaussian_summand_conditional_mean(sampling$mean[1L], sampling$variance[1L], 2^-60, dat$yi[1L])
  expect_gt(centre, 4e5)
  expect_true(any(abs(centres - centre) <= 1e-9 * centre))

  # Its lower tail at y = .4 underflows, so the tail-factored ratio integrates the
  # reciprocal at proposed deleted values. At this scale those round to 0 or
  # -2^25 from the location; the lost-share form would centre on them, the share
  # form adds the random-effect share of the retained error, about 5.0e5.
  conditional_mean <- .gaussian_summand_conditional_mean
  shares <- list()
  centres <- numeric()
  testthat::local_mocked_bindings(
    .gaussian_summand_conditional_mean = function(mean, variance, other, observed) {
      value <- conditional_mean(mean, variance, other, observed)
      shares[[length(shares) + 1L]] <<- c(mean = mean, variance = variance, other = other, value = value)
      value
    }, .package = "RoBMA")
  lower <- .selection_conditioned_sampling_estimate_targets_chunk(setup, "log_lower")
  expect_true(is.finite(lower$log_lower[1L, 1L]))
  # Row 1's calls; the normalizer mock records the first column only.
  shares <- matrix(as.numeric(unlist(shares)), ncol = 4L, byrow = TRUE,
    dimnames = list(NULL, c("mean", "variance", "other", "value")))
  shares <- shares[shares[, "mean"] == sampling$mean[1L], , drop = FALSE]
  expect_gt(nrow(shares), 0L)
  expect_true(all(shares[, "variance"] == sampling$variance[1L] & shares[, "other"] == 2^-60))
  expect_true(all(vapply(shares[, "value"], function(value) any(abs(centres - value) <= 1e-9 * abs(value)),
    logical(1L))))
  expect_gt(max(shares[, "value"]), 4e5)
})


test_that("fixed-zero random bridge replay preserves its fitted source map", {

  point <- function(value) BayesTools::prior("point", list(value))
  yi <- c(.4, -.2)
  se <- c(.8, 1.1)
  object <- bselmodel.mv(yi = yi, V = diag(se^2), random = ~ 1 | study,
    data = data.frame(study = c("a", "b")), measure = "GEN", effect_direction = "positive",
    prior_unit_information_sd = 1, prior_effect = point(.1), prior_heterogeneity = point(0),
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
      model = BayesTools::selection_model(known_sampling_variance = "condition", group = "study")),
    only_priors = TRUE)
  fit_priors <- .create_fit_priors(object$data, object$priors)
  fit <- structure(list(), formula_design = object$formula_design,
    prior_list = c(object$formula_design$mu$prior_list, fit_priors))
  object$fit <- fit
  bridge <- .marglik_fixed_zero_random_setup(object, fit, fit_priors)
  nodes <- c(mu_intercept = .1, mu__xREx__study_intercept = 0,
    "mu__xREx__study_xRE_Zx[1,1]" = .2, "mu__xREx__study_xRE_Zx[2,1]" = -.3,
    "sampling_z[1]" = .4, "sampling_z[2]" = -.5)
  fit_data <- .create_fit_data(object$data, object$priors)
  actual <- .log_posterior(parameters = list(), data = fit_data,
    is_scale = FALSE, is_multilevel = FALSE, is_weights = FALSE, is_known_v = TRUE,
    is_random = TRUE, is_PET = FALSE, is_PEESE = FALSE, is_weightfunction = TRUE,
    effect_direction = "positive", outcome_type = "norm", model_data = object$data,
    bridge_context = list(nodes = nodes), selection_fit = bridge$fit, selection_priors = object$priors)
  expect_equal(actual, sum(stats::dnorm(yi, .1, se, log = TRUE)), tolerance = 1e-12)
  expect_identical(bridge$fit_priors, fit_priors)
})


test_that("sampling-conditioned deletion retains a density with singular integrated random covariance", {

  yi <- c(.4, -.2)
  se <- c(.8, 1.1)
  object <- bselmodel(
    yi = yi, sei = se, cluster = c("a", "a"), measure = "GEN",
    effect_direction = "positive", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
      BayesTools::wf_fixed(c(1, .2)), model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "condition")),
    selection_control = set_selection_likelihood_control(relative_tolerance = 1e-6)
  )
  samples <- matrix(c(.1, .5, 1, .3, -.4, .2), 1L,
    dimnames = list(NULL, c("mu", "tau", "rho", "gamma[1]", "sampling_z[1]", "sampling_z[2]")))
  setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
  state <- .selection_conditioned_sampling_state(setup)
  expect_equal(.selection_covariance_draw(state$integrated_covariance, 1L), matrix(.25, 2L, 2L), tolerance = 0)
  actual <- .log_lik_estimate_from_setup(setup)
  # A shared random intercept gives a rank-one integrated covariance. Once
  # the other row and its sampling error are retained, that intercept is known;
  # the deleted row still varies through its own sampling error.
  mass <- function(e1, e2) {
    cuts <- cbind(-.1 - e1, -.1 - e2)
    .04 + .16 * stats::pnorm(-apply(cuts, 1L, min) / .5) +
      .8 * stats::pnorm(-apply(cuts, 1L, max) / .5)
  }
  expected <- numeric(2L)
  for (row in 1:2) {
    other <- 3L - row
    other_error <- state$e[1, other]
    intercept <- yi[other] - .1 - other_error
    error <- yi[row] - .1 - intercept
    denominator <- stats::integrate(function(e) {
      stats::dnorm(e, 0, se[row]) * ifelse(.1 + intercept + e >= 0, 1, .2) /
        mass(e, rep(other_error, length(e)))
    }, -Inf, Inf, rel.tol = 1e-9, subdivisions = 1000L)$value
    expected[row] <- stats::dnorm(error, 0, se[row], log = TRUE) +
      log(ifelse(yi[row] >= 0, 1, .2)) - log(mass(error, other_error)) - log(denominator)
  }
  expect_equal(as.numeric(actual), expected, tolerance = 2e-6)

  # Deleting the complete publication integrates its sampling vector. Because
  # its random covariance has rank one, an independent scalar conditional
  # random-intercept integral supplies a reference for the vector QMC target.
  attr(setup$data, "selection_execution_plan") <- .selection_joint_execution_plan_with_control(
    .data_selection_execution_plan(setup$data),
    set_selection_likelihood_control(relative_tolerance = .002, max_points_per_scramble = 65536L)
  )
  joint <- .selection_conditioned_sampling_deletion_loglik(setup, list(1:2))[1L, 1L]
  variance <- 1 / (1 / .25 + sum(1 / se^2))
  mean <- variance * sum((yi - .1) / se^2)
  reciprocal <- stats::integrate(function(u) {
    stats::dnorm(u, mean, sqrt(variance)) / mass(yi[1] - .1 - u, yi[2] - .1 - u)
  }, -Inf, Inf, rel.tol = 1e-10)$value
  reference <- mvtnorm::dmvnorm(yi, rep(.1, 2L), diag(se^2) + matrix(.25, 2L, 2L), log = TRUE) +
    log(.2) + log(reciprocal)
  expect_equal(joint, reference, tolerance = .006)
})


test_that("conditioned-sampling deletion distinguishes deterministic and failed numerical targets", {

  fixture <- .sampling_selection_evidence_fixture()
  object <- fixture$object
  setup <- .log_lik_posterior_setup(object$fit, fixture$samples[1L, , drop = FALSE],
    object$data, object$priors, "estimate", NULL)
  plan <- .data_selection_execution_plan(setup$data)
  plan$row_blocks <- list(1:2)
  attr(setup$data, "selection_execution_plan") <- plan
  candidate <- tcrossprod(c(.5, -.5))
  state <- list(e = matrix(-.1, 1L, 2L), baseline_mu = setup$mu,
    sampling_covariance = matrix(1, 2L, 2L), integrated_covariance = array(candidate, c(1L, 2L, 2L)),
    total_covariance = array(matrix(1, 2L, 2L) + candidate, c(1L, 2L, 2L)))
  testthat::local_mocked_bindings(.selection_conditioned_sampling_state = function(setup) state,
    .selection_conditioned_sampling_normalizer = function(...) list(log_mass = 0), .package = "RoBMA")
  error <- tryCatch(.selection_joint_estimate_targets(setup, "log_density"), error = identity)
  expect_identical(conditionMessage(error), paste0(
    "Selection estimate deletion is unavailable because the deleted outcome is determined by retained sampling errors and outcomes. ",
    "Use 'known_sampling_variance = \"integrate\"' when fitting to obtain sampling-marginal deletion scores."
  ))
  expect_null(conditionCall(error))
  state$sampling_covariance <- diag(1, 2L)
  state$total_covariance <- array(diag(1, 2L) + candidate, c(1L, 2L, 2L))
  testthat::local_mocked_bindings(integrate = function(...) list(value = 1, abs.error = 1, message = "OK"), .package = "stats")
  error <- tryCatch(.selection_joint_estimate_targets(setup, "log_density"), error = identity)
  expect_identical(conditionMessage(error), paste0(
    "Selection estimate deletion was rejected by diagnostics: relative integration error was 1. ",
    "Increase 'max_points_per_scramble' in 'selection_control = set_selection_likelihood_control()' when fitting."
  ))
  expect_null(conditionCall(error))
})


test_that("singleton deletion evaluates the existing scalar law without auxiliary state", {

  components <- c("log_density", "cdf", "log_lower", "log_upper", "mean", "variance")
  cases <- list()
  for (estimate in c("integrate", "condition")) {
    for (direction in c("positive", "negative")) {
      fixture <- .sampling_selection_evidence_fixture(estimate)
      object <- fixture$object
      samples <- rbind(fixture$samples, fixture$samples[1L, , drop = FALSE])
      samples[3L, "tau"] <- 0
      if (direction == "negative") {
        attr(object$data, "effect_direction") <- direction
        object$data$outcome$yi <- -object$data$outcome$yi
        samples[, c("mu", "theta[1]", "theta[2]")] <-
          -samples[, c("mu", "theta[1]", "theta[2]")]
      }
      setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors,
        "estimate", NULL)
      signed_y <- if (direction == "negative") -setup$yi else setup$yi
      selection <- .selection_joint_signed_context(setup, signed_y)
      state <- .selection_conditioned_sampling_state(setup)
      expected <- .selection_conditioned_sampling_independent_targets(
        setup, state, selection, components)
      attr(expected, "dependency_blocks") <- .data_selection_execution_plan(setup$data)$row_blocks
      if (estimate == "condition") {
        expect_null(.selection_conditioned_sampling_independent_targets(
          setup, NULL, selection, components))
        expect_equal(.selection_conditioned_sampling_estimate_targets(setup, components),
          expected, tolerance = 1e-12)
      } else cases[[paste(estimate, direction)]] <- list(setup = setup, expected = expected)
      # Hard weights retain the general path's support/boundary checks.
      selection$omega[, 2L] <- 0
      expect_null(.selection_conditioned_sampling_independent_targets(
        setup, NULL, selection, components))
    }
  }
  testthat::local_mocked_bindings(
    .selection_conditioned_sampling_state = function(...) {
      stop("Singleton deletion must not construct auxiliary state.")
    }, .package = "RoBMA")
  for (case in cases) {
    # Dropping mathematically irrelevant auxiliary columns also checks that the
    # optimized route is independent of their names and monitored realizations.
    case$setup$posterior_samples <- case$setup$posterior_samples[, c("mu", "tau"), drop = FALSE]
    actual <- .selection_conditioned_sampling_estimate_targets(case$setup, components)
    expect_equal(actual, case$expected, tolerance = 1e-12)
    # The native log-density-only route may stop before the joint six-target
    # rule; both retain the fixture's 1e-6 relative integration criterion.
    density <- .selection_conditioned_sampling_estimate_targets(case$setup, "log_density")
    expect_equal(density$log_density, case$expected$log_density, tolerance = 2e-6)
  }
})


test_that("native singleton deletion preserves tails and per-cell quadrature rejection", {

  native <- function(y, mu, tau, omega, modes, orders = SELNORM_CLUSTER_QUADRATURE_ORDERS,
                     tolerance = 1e-6) {

    S <- length(mu)
    K <- length(y)
    quadrature <- .selection_joint_cluster_quadrature_rules(orders)
    variance <- matrix(tau^2, S, K)
    .Call("RoBMA_selnorm_sampling_deletion_loglik_batch",
      as.numeric(y), matrix(mu, S, K), rep(.8^2, K), variance, variance + .8^2,
      rep(.8, K), cbind(1, omega), c(0, -Inf), c(Inf, 0),
      as.integer(ifelse(y >= 0, 1L, 2L)), as.integer(modes), TRUE,
      as.numeric(quadrature$nodes), as.numeric(quadrature$log_weights),
      as.integer(quadrature$orders), as.numeric(tolerance), PACKAGE = "RoBMA")
  }
  y <- c(-100, 100)
  mu <- c(0, .1, -.2)
  tau <- c(.5, 0, .4)
  actual <- native(y, mu, tau, c(1e-320, .2, .02),
    c(SELKERNEL_STEP, SELKERNEL_STEP, SELKERNEL_NORMAL))
  reference <- vapply(y, function(value) stats::dnorm(value, mu,
    sqrt(.8^2 + tau^2), log = TRUE), numeric(3L))
  # At these tails, inverse acceptance cancels the observed weight up to
  # Gaussian tail terms far below binary64 range, even for the subnormal
  # positive weight. Zero heterogeneity and normal branches cancel exactly.
  expect_true(all(is.finite(actual)))
  expect_equal(actual, reference, tolerance = 1e-9)

  # This deliberately short internal rule sequence cannot resolve the first
  # selected cell at the requested tolerance; the normal cell remains valid.
  rejected <- native(0, c(0, 0), c(.05, .05), c(1e-4, 1e-4),
    c(SELKERNEL_STEP, SELKERNEL_NORMAL), orders = c(15L, 31L), tolerance = 1e-12)
  expect_true(is.na(rejected[1L, 1L]))
  expect_equal(rejected[2L, 1L], stats::dnorm(0, 0, sqrt(.8^2 + .05^2), log = TRUE),
    tolerance = 0)
  resolved <- native(0, 0, .05, 1e-4, SELKERNEL_STEP)
  ratio <- .8 / sqrt(.8^2 + .05^2)
  reciprocal <- stats::integrate(function(z) {
    stats::dnorm(z) / (1e-4 + (1 - 1e-4) * stats::pnorm(ratio * z))
  }, -Inf, Inf, rel.tol = 1e-11)$value
  reference <- stats::dnorm(0, 0, sqrt(.8^2 + .05^2), log = TRUE) + log(reciprocal)
  expect_equal(as.numeric(resolved), reference, tolerance = 2e-6)
})


test_that("native singleton deletion retains nonmonotone step weights", {

  yi <- c(-.2, .1, .8)
  sei <- c(.2, .25, .3)
  mu <- .1
  tau <- .35
  object <- bselmodel(yi = yi, sei = sei, measure = "GEN", effect_direction = "positive",
    prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction("one-sided", c(.025, .5),
      BayesTools::wf_fixed(c(1, 3, .2)),
      model = BayesTools::selection_model(known_sampling_variance = "condition")),
    selection_control = set_selection_likelihood_control(relative_tolerance = 1e-6))
  samples <- matrix(c(mu, tau), 1L, dimnames = list(NULL, c("mu", "tau")))
  setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors,
    "estimate", NULL)
  actual <- .selection_conditioned_sampling_estimate_targets(setup, "log_density")$log_density
  reference <- vapply(seq_along(yi), function(i) {
    cutoff <- sei[i] * stats::qnorm(.975)
    observed_weight <- if (yi[i] >= cutoff) 1 else if (yi[i] >= 0) 3 else .2
    log(stats::integrate(function(e) {
      upper <- stats::pnorm((mu + e - cutoff) / tau)
      lower <- stats::pnorm(-(mu + e) / tau)
      mass <- upper + 3 * (stats::pnorm((mu + e) / tau) - upper) + .2 * lower
      stats::dnorm(e, 0, sei[i]) * stats::dnorm(yi[i], mu + e, tau) * observed_weight / mass
    }, -Inf, Inf, rel.tol = 1e-10)$value)
  }, numeric(1L))
  expect_equal(as.numeric(actual), reference, tolerance = 2e-6)
})
