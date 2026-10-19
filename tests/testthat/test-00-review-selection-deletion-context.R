.review_sampling_deletion_setup <- function(retained = TRUE, cluster = FALSE, weight = 1) {

  args <- list(yi = c(.4, -.25), sei = c(1, 1), measure = "GEN",
    effect_direction = "positive", prior_unit_information_sd = 1,
    prior_effect = prior("point", list(0)), prior_heterogeneity = prior("point", list(1)),
    prior_bias = prior_weightfunction("one-sided", .5, wf_fixed(c(1, weight)),
      model = selection_model(estimate_random_effects = if (retained) "condition" else "integrate",
        known_sampling_variance = "condition")), only_priors = TRUE, silent = TRUE)
  if (cluster) {
    args$cluster <- factor(c("a", "b"))
    args$prior_heterogeneity_allocation <- prior("point", list(.36))
  }
  object <- do.call(bselmodel, args)
  samples <- cbind(mu = 0, tau = 1, `theta[1]` = .25, `theta[2]` = -.4,
    `sampling_z[1]` = .3, `sampling_z[2]` = -.2)
  if (cluster) samples <- cbind(samples, rho = .36, `gamma[1]` = .5, `gamma[2]` = -.4)
  list(object = object, setup = .log_lik_posterior_setup(object$fit, samples,
    object$data, object$priors, "estimate", NULL))
}

test_that("internal deletion covariance guards distinguish ordinary and retained-sampling remedies", {

  # Public admission rejects singular total covariance. This is an internal
  # failed-factorization diagnostic control using otherwise valid live states.
  for (weight in c(1, .2)) {
    fixture <- .review_sampling_deletion_setup(retained = FALSE, weight = weight)
    state <- .selection_conditioned_sampling_state(fixture$setup)
    testthat::with_mocked_bindings({
      result <- tryCatch(.selection_conditioned_sampling_deletion_loglik_chunk(
        fixture$setup, list(1L)), error = identity)
      expect_s3_class(result, "error")
      expected <- if (weight == 1) {
        "Selection deletion is unavailable because the deleted outcome vector has no density conditional on retained outcomes. Inspect the supplied sampling covariance and random-effect structure."
      } else {
        "Selection deletion is unavailable because the deleted outcome vector has no density conditional on retained sampling errors and outcomes. Use 'known_sampling_variance = \"integrate\"' when fitting to obtain sampling-marginal deletion scores."
      }
      expect_identical(conditionMessage(result), expected)
    }, .selection_conditioned_sampling_state = function(setup) state,
    .covariance_cholesky = function(...) NULL, .package = "RoBMA")
  }
})

test_that("estimate deletion retains declared random context when selection cancels", {

  for (retained in c(FALSE, TRUE)) {
    fixture <- .review_sampling_deletion_setup(retained)
    setup <- fixture$setup
    # Independent Gaussian augmentation, separately from the production state:
    # E0=(.3,-.2), U0=(.25,-.4), conditional correction (y-E0-U0)/2.
    baseline <- c(.25, -.4) + (c(.4, -.25) - c(.3, -.2) - c(.25, -.4)) / 2
    mean <- if (retained) baseline else c(0, 0)
    variance <- if (retained) c(1, 1) else c(2, 2)
    actual <- .selection_joint_estimate_targets(setup,
      c("mean", "variance", "log_density", "cdf", "log_lower", "log_upper"))
    expected <- list(mean = mean, variance = variance,
      log_density = dnorm(setup$yi, mean, sqrt(variance), log = TRUE),
      cdf = pnorm(setup$yi, mean, sqrt(variance)),
      log_lower = pnorm(setup$yi, mean, sqrt(variance), log.p = TRUE),
      log_upper = pnorm(setup$yi, mean, sqrt(variance), lower.tail = FALSE, log.p = TRUE))
    for (component in names(expected)) expect_equal(as.numeric(actual[[component]]),
      expected[[component]], tolerance = 1e-12)
  }
})

test_that("singleton-cluster deletion keeps its distinct joint cancellation law", {

  fixture <- .review_sampling_deletion_setup(retained = FALSE, cluster = TRUE)
  setup <- fixture$setup
  old_cluster <- dnorm(setup$yi, 0, sqrt(2), log = TRUE)
  cluster <- .selection_conditioned_sampling_deletion_loglik(setup, setup$cluster)
  estimate <- .selection_joint_estimate_targets(setup, c("mean", "variance", "log_density"))
  e0 <- c(.3, -.2)
  u0 <- .8 * c(.25, -.4)
  b0 <- .6 * c(.5, -.4)
  baseline <- b0 + .36 / 2 * (setup$yi - b0 - u0 - e0)
  expect_equal(as.numeric(cluster), old_cluster, tolerance = 1e-12)
  expect_equal(as.numeric(estimate$mean), baseline, tolerance = 1e-12)
  expect_equal(as.numeric(estimate$variance), c(1.64, 1.64), tolerance = 1e-12)
  expect_equal(as.numeric(estimate$log_density), dnorm(setup$yi, baseline, sqrt(1.64), log = TRUE), tolerance = 1e-12)
})

test_that("estimate context identity and legacy normalization share the score owner", {

  for (retained in c(FALSE, TRUE)) {
    fixture <- .review_sampling_deletion_setup(retained)
    object <- fixture$object
    setup <- fixture$setup
    score <- .log_lik_from_posterior_samples(object$fit, setup$posterior_samples,
      object$data, object$priors, unit = "estimate", add_metadata = TRUE,
      data_hash = .get_outcome_hash(object))
    metadata <- attr(score, "RoBMA_target")
    cdf <- .cdf_lik_estimate.brma(object, setup)
    key <- .current_predictive_target_key(object, "estimate")
    expect_identical(metadata$sampling_deletion_policy, "declared_retained_context")
    expect_identical(attr(cdf, "RoBMA_target")$retained_random_context,
      metadata$retained_random_context)
    expect_identical(.predictive_target_fingerprint(metadata), .predictive_target_fingerprint(key))
    expect_identical(metadata$retained_random_context == "none", !retained)
    legacy <- metadata
    legacy[c("sampling_deletion_policy", "retained_random_context")] <- NULL
    legacy$ordinary_branch_handling <- "all_local_Gaussian_sources_integrated_when_selection_cancels"
    before <- serialize(legacy, NULL)
    if (retained) {
      expect_null(.predictive_target_fingerprint(legacy))
      expect_error(.check_cached_predictive_target(object, legacy, "estimate", "LOO"), "Recompute")
    } else {
      expect_identical(.predictive_target_fingerprint(legacy), .predictive_target_fingerprint(key))
      expect_silent(.check_cached_predictive_target(object, legacy, "estimate", "LOO"))
    }
    expect_identical(serialize(legacy, NULL), before)
  }
  cluster <- .review_sampling_deletion_setup(retained = FALSE, cluster = TRUE)$object
  retained <- .review_sampling_deletion_setup(retained = TRUE)$object
  expect_false(identical(.current_predictive_target_key(cluster, "estimate")$retained_random_context,
    .current_predictive_target_key(retained, "estimate")$retained_random_context))
  expect_null(.current_predictive_target_key(cluster, "cluster")$sampling_deletion_policy)
})


test_that("correlated sampling deletion is continuous at constant selection weights", {

  # Complete compiler-declared supplied states, not a fitted posterior chain.
  # With retained e=0 and integrated G=I, the remaining sampling coordinate
  # mixes h~N(0,.75) with the selected N(h,1) observation kernel.
  reference <- function(weight) {
    if (weight == 1) return(c(log_density = stats::dnorm(0, sd = sqrt(1.75), log = TRUE),
      cdf = .5, log_lower = log(.5), log_upper = log(.5), mean = 0, variance = 1.75))
    integrate_h <- function(f) stats::integrate(function(h) f(h) * stats::dnorm(h, sd = sqrt(.75)),
      -Inf, Inf, rel.tol = 1e-10, abs.tol = 1e-12, subdivisions = 1000L)$value
    mass <- function(h) weight * stats::pnorm(-h) + stats::pnorm(h)
    shift <- function(h) (1 - weight) * stats::dnorm(h) / mass(h)
    mean_h <- function(h) h + shift(h)
    variance_h <- function(h) 1 - (1 - weight) * h * stats::dnorm(h) / mass(h) - shift(h)^2
    mean <- integrate_h(mean_h)
    lower <- integrate_h(function(h) weight * stats::pnorm(-h) / mass(h))
    upper <- integrate_h(function(h) stats::pnorm(h) / mass(h))
    c(log_density = log(integrate_h(function(h) stats::dnorm(h) / mass(h))),
      cdf = lower, log_lower = log(lower), log_upper = log(upper), mean = mean,
      variance = integrate_h(function(h) variance_h(h) + (mean_h(h) - mean)^2))
  }
  for (weight in c(1, 1 - 1e-6)) {
    setup <- review_compile_object(bselmodel.mv(yi = c(0, 0),
      V = matrix(c(1, .5, .5, 1), 2L), random = ~ 1 | estimate,
      data = data.frame(estimate = factor(1:2)), measure = "GEN", effect_direction = "positive",
      prior_effect = BayesTools::prior("point", list(location = 0)),
      prior_heterogeneity = BayesTools::prior_random(sd = BayesTools::prior("point", list(location = 1))),
      prior_bias = BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, weight)),
        model = BayesTools::selection_model(known_sampling_variance = "condition", weight_rule = "product")),
      prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE))
    priors <- unlist(lapply(setup$compiled, `[[`, "prior_list"), recursive = FALSE)
    names(priors) <- unlist(lapply(setup$compiled, function(part) names(part$prior_list)), use.names = FALSE)
    priors <- c(priors, .create_fit_priors(setup$object$data, setup$object$priors))
    columns <- setdiff(review_prior_columns(priors), "omega")
    samples <- matrix(0, 24L, length(columns), dimnames = list(NULL, columns))
    for (name in names(priors)) if (BayesTools::is.prior.point(priors[[name]])) {
      samples[, name] <- priors[[name]]$parameters$location
    }
    terms <- .formula_design_random_effects_by_mode(setup$object$formula_design$mu, "sampled")
    extras <- unlist(lapply(terms, function(term) {
      base <- paste0(term$parameter_stem, "_xRE_Zx")
      expect_true(base %in% setup$compiled$mu$add_parameters)
      unlist(lapply(seq_len(term$n_columns), function(column) {
        paste0(base, "[", seq_len(term$n_groups), ",", column, "]")
      }))
    }), use.names = FALSE)
    samples <- cbind(samples, matrix(0, 24L, length(extras), dimnames = list(NULL, extras)))
    generated <- BayesTools::JAGS_evaluate_deterministic(review_attach_replay(setup, samples)$fit, draws = samples)
    for (name in colnames(generated)) {
      if (name %in% colnames(samples)) samples[, name] <- generated[, name]
      else samples <- cbind(samples, generated[, name, drop = FALSE])
    }
    required <- c(unlist(lapply(BayesTools::JAGS_to_monitor(priors), function(name) {
      if (name == "omega") paste0("omega[", 1:2, "]")
      else if (name == "sampling_z") paste0("sampling_z[", 1:2, "]") else name
    })), extras)
    expect_true(all(required %in% colnames(samples)))
    object <- review_attach_replay(setup, samples)
    expect_identical(class(object), class(setup$object))
    BayesTools::JAGS_validate_fit_contract(object$fit, requires = c("name_encoding", "formula_name_map",
      "formula_design", "parameter_map", "draw_geometry"))
    likelihood <- .estimate_likelihood_setup.brma(object)
    state <- .selection_conditioned_sampling_state(likelihood)
    expect_true(all(state$e == 0))
    expect_true(all(state$baseline_mu == 0))
    expect_equal(unname(.selection_covariance_draw(state$integrated_covariance, 1L)), diag(2L), tolerance = 0)
    expected <- reference(weight)
    actual <- .selection_joint_estimate_targets(likelihood, names(expected))
    tolerance <- if (weight == 1) 1e-12 else 2e-6
    for (component in names(expected)) {
      expect_true(max(abs(actual[[component]] - expected[[component]])) <= tolerance, info = component)
    }
    expect_equal(as.numeric(log_lik(object)), as.numeric(actual$log_density), tolerance = 0)
  }
})
