test_that("single-cluster pseudo-BMA retains its target dimension with unequal draws", {

  # Original C2 probes A/C: compatible synthetic stored LOO scores, not MCMC fits.
  object <- brma(yi = c(.3, -.2, .5), sei = c(.1, .2, .15),
    cluster = rep("a", 3), measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE)
  make_stored <- function(S, shift) {
    out <- object
    scores <- matrix(-1 + .01 * sin(seq_len(S)) + shift, S, 1L)
    result <- loo::loo(scores, r_eff = 1, save_psis = TRUE, cores = 1)
    attr(result, "RoBMA_target") <- .current_predictive_target_key(out, "cluster")
    out[["loo"]] <- list(cluster = result)
    out
  }
  a <- make_stored(200L, 0)
  b <- make_stored(300L, .2)
  c <- make_stored(200L, .2)
  pointwise <- c(a[["loo"]][["cluster"]][["pointwise"]][, "elpd_loo"],
                 b[["loo"]][["cluster"]][["pointwise"]][, "elpd_loo"])
  expected <- exp(pointwise - max(pointwise))
  expected <- expected / sum(expected)
  actual <- loo_model_weights(a, b, unit = "cluster", method = "pseudobma", BB = FALSE)
  expect_identical(names(actual), c("a", "b"))
  expect_s3_class(actual, "pseudobma_weights")
  expect_equal(as.numeric(actual), unname(expected), tolerance = 1e-14)
  # Equal-count delegation remains a separate control; upstream single-target
  # stacking has its own recorded dependency limitation and is not changed here.
  equal <- loo_model_weights(a, c, unit = "cluster", method = "pseudobma", BB = FALSE)
  equal_pointwise <- c(a[["loo"]][["cluster"]][["pointwise"]][, "elpd_loo"],
                       c[["loo"]][["cluster"]][["pointwise"]][, "elpd_loo"])
  equal_expected <- exp(equal_pointwise - max(equal_pointwise))
  expect_s3_class(equal, "pseudobma_weights")
  expect_equal(as.numeric(equal), unname(equal_expected / sum(equal_expected)), tolerance = 1e-14)
})


test_that("cluster log scores expose user labels while retaining row partitions", {

  object <- brma(yi = c(.1, .2, .3), sei = c(.2, .3, .4),
    cluster = c("zeta", "zeta", "alpha"), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE)
  samples <- cbind(mu = c(.05, .1), tau = .3, rho = .4,
    "gamma[1]" = c(-.2, .3), "gamma[2]" = c(.5, -.4))
  scores <- .log_lik_from_posterior_samples(object[["fit"]], samples,
    object[["data"]], object[["priors"]], unit = "cluster", add_metadata = TRUE,
    data_hash = .get_outcome_hash(object))
  target <- attr(scores, "RoBMA_target", exact = TRUE)
  expect_identical(colnames(scores), c("log_lik_cluster[zeta]", "log_lik_cluster[alpha]"))
  expect_identical(target[["targets"]], c("zeta", "alpha"))
  expect_identical(target[["cluster_partition"]], list(1:2, 3L))
  setup <- .log_lik_posterior_setup(object[["fit"]], samples,
    object[["data"]], object[["priors"]], "estimate", NULL)
  expect_error(.log_lik_estimate_sum_from_setup(setup), paste0(
    "Internal error: a multilevel normal joint log-likelihood requires ",
    "the cluster likelihood evaluator, not an estimate-unit sum."
  ), fixed = TRUE)
})


test_that("check_loo uses stored draw counts and reports observed Pareto values", {

  object <- brma(yi = c(.1, -.2), sei = c(.3, .4), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE)
  scores <- cbind(-1 + .01 * sin(seq_len(100L)), -2 + .02 * cos(seq_len(100L)))
  result <- loo::loo(scores, r_eff = 1, save_psis = TRUE, cores = 1)
  attr(result, "RoBMA_target") <- .current_predictive_target_key(object, "estimate")
  # Controlled diagnostics on a real loo object isolate the warning boundary.
  result[["diagnostics"]][["pareto_k"]] <- c(.6, .5)
  object[["loo"]] <- list(estimate = result)
  expect_warning(check_loo(object), paste0(
    "Some Pareto k values are high, indicating potentially unreliable LOO ",
    "diagnostics for estimates: 1 (k = 0.6). ",
    "Inspect the LOO result with 'loo(object, unit = \"estimate\")'."
  ), fixed = TRUE)
  result[["diagnostics"]][["pareto_k"]] <- c(.5, .5)
  object[["loo"]][["estimate"]] <- result
  expect_warning(check_loo(object), NA)
})


test_that("non-finite log scores and invalid relative ESS have distinct diagnostics", {

  object <- brma(yi = c(.1, -.2), sei = c(.3, .4), measure = "GEN",
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("point", list(0)),
    prior_unit_information_sd = 1, only_priors = TRUE)
  setup <- review_compile_object(object)
  parameter <- if (is.null(setup[["object"]][["formula_design"]][["mu"]])) {
    "mu"
  } else review_coefficient_name(setup, "intercept")
  draws <- matrix(c(1e200, -1e200, 1e200, -1e200), ncol = 1L,
                  dimnames = list(NULL, parameter))
  # Real public constructor metadata and explicit synthetic rows, not MCMC fits.
  nonfinite <- review_attach_replay(setup, draws)
  expect_true(all(!is.finite(log_lik(nonfinite))))
  expect_error(add_loo(nonfinite), paste0(
    "LOO diagnostics are unavailable because log-likelihood draws contain ",
    "non-finite values. Inspect 'log_lik()' before retrying."
  ), fixed = TRUE)
  draws[, 1L] <- c(.05, .1, .2, .15)
  finite <- review_attach_replay(setup, draws)
  expect_true(all(is.finite(log_lik(finite))))
  expect_error(add_loo(finite, r_eff = NA_real_),
    "'r_eff' must be finite and positive for every varying observation.", fixed = TRUE)
})

test_that("check_loo retains small-sample and capped Pareto boundaries", {

  object <- brma(yi = c(.1, -.2), sei = c(.3, .4), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE)
  for (S in c(50L, 4000L)) {
    scores <- cbind(-1 + .01 * sin(seq_len(S)), -2 + .02 * cos(seq_len(S)))
    result <- loo::loo(scores, r_eff = 1, save_psis = TRUE, cores = 1)
    attr(result, "RoBMA_target") <- .current_predictive_target_key(object, "estimate")
    boundary <- if (S == 50L) 1 - 1 / log10(50) else .7
    result$diagnostics$pareto_k <- rep(boundary, 2L)
    object$loo <- list(estimate = result)
    expect_warning(check_loo(object), NA)
    result$diagnostics$pareto_k <- c(boundary + .01, boundary)
    object$loo$estimate <- result
    expected <- paste0("Some Pareto k values are high, indicating potentially unreliable LOO ",
      "diagnostics for estimates: 1 (k = ", format(boundary + .01, digits = 4), "). ",
      "Inspect the LOO result with 'loo(object, unit = \"estimate\")'.")
    expect_warning(check_loo(object), expected, fixed = TRUE)
  }
})
