.review_sampling_deletion_setup <- function(retained = TRUE, cluster = FALSE) {

  args <- list(yi = c(.4, -.25), sei = c(1, 1), measure = "GEN",
    effect_direction = "positive", prior_unit_information_sd = 1,
    prior_effect = prior("point", list(0)), prior_heterogeneity = prior("point", list(1)),
    prior_bias = prior_weightfunction("one-sided", .5, wf_fixed(c(1, 1)),
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
