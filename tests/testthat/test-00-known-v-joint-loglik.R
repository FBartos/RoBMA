# test-00-known-v-joint-loglik.R

test_that("known-V joint log-likelihood uses block MVN density", {

  V    <- matrix(c(.04, .015, .015, .09), nrow = 2L)
  data <- list(outcome = data.frame(yi = c(.10, -.20), sei = c(.20, .30)))
  attr(data, "known_V")      <- TRUE
  attr(data, "known_V_data") <- .known_v_prepare(
    V                         = V,
    keep_rows                 = rep(TRUE, nrow(V)),
    known_v_parameterization  = "block_mvn",
    known_v_residual_fraction = NULL
  )
  attr(data, "random")       <- FALSE

  setup <- list(
    outcome_type      = "norm",
    is_weightfunction = FALSE,
    weights           = NULL,
    data              = data,
    K                 = 2L,
    S                 = 2L,
    yi                = c(.10, -.20),
    mu                = matrix(c(.02, -.10, .06, -.16), nrow = 2L, byrow = TRUE),
    tau_within        = matrix(c(.05, .08, .04, .06), nrow = 2L, byrow = TRUE),
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 2L, ncol = 0L),
    marginalized_random_source_samples = NULL
  )

  out <- .log_lik_known_v_joint_sum_from_setup(setup)
  ref <- vapply(seq_len(setup[["S"]]), function(s) {
    covariance <- V + diag(setup[["tau_within"]][s, ]^2, nrow = 2L)
    .marglik_mvn_log_density(
      y          = setup[["yi"]],
      mean       = setup[["mu"]][s, ],
      covariance = covariance
    )
  }, numeric(1))

  expect_equal(out, ref)
  expect_equal(.log_lik_estimate_sum_from_setup(setup), ref)
  expect_false(isTRUE(all.equal(
    out,
    rowSums(.log_lik_known_v_estimate_target_from_setup(setup))
  )))

  singular_setup <- setup
  singular_setup[["data"]] <- data
  singular_setup[["data"]][["outcome"]][["sei"]] <- c(1, 1)
  attr(singular_setup[["data"]], "known_V_data") <-
    .known_v_prepare(
      V                         = matrix(1, nrow = 2L, ncol = 2L),
      keep_rows                 = rep(TRUE, 2L),
      known_v_parameterization  = "block_mvn",
      known_v_residual_fraction = NULL,
      warn_singular             = FALSE
    )
  singular_setup[["tau_within"]] <- matrix(0, nrow = 2L, ncol = 2L)
  expect_error(
    .log_lik_known_v_joint_sum_from_setup(singular_setup),
    "Known-V joint likelihood covariance"
  )
})


test_that("diagonal known-V estimate log-likelihood is exactly vectorized", {

  sampling_variance <- c(.04, .09, .16, .25)
  V <- diag(sampling_variance)
  data <- list(outcome = data.frame(
    yi  = c(.10, -.20, .30, -.40),
    sei = sqrt(sampling_variance)
  ))
  attr(data, "known_V") <- TRUE
  attr(data, "known_V_data") <- .known_v_prepare(
    V                         = V,
    keep_rows                 = rep(TRUE, nrow(V)),
    known_v_parameterization  = "block_mvn",
    known_v_residual_fraction = NULL
  )
  attr(data, "random") <- FALSE

  mu <- matrix(
    c(.02, -.10, .06, -.16, .04, -.12, .08, -.18, .06, -.14, .10, -.20),
    nrow = 3L,
    byrow = TRUE
  )
  tau <- matrix(
    c(.05, .08, .04, .06, .06, .09, .05, .07, .07, .10, .06, .08),
    nrow = 3L,
    byrow = TRUE
  )
  setup <- list(
    outcome_type      = "norm",
    is_weightfunction = FALSE,
    weights           = NULL,
    data              = data,
    K                 = 4L,
    S                 = 3L,
    yi                = data[["outcome"]][["yi"]],
    mu                = mu,
    tau_within        = tau,
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 3L, ncol = 0L),
    marginalized_random_source_samples = NULL
  )

  expected <- matrix(NA_real_, nrow = 3L, ncol = 4L)
  for (s in seq_len(nrow(expected))) {
    for (k in seq_len(ncol(expected))) {
      variance <- sampling_variance[[k]] + tau[s, k]^2
      residual <- data[["outcome"]][["yi"]][[k]] - mu[s, k]
      expected[s, k] <- -0.5 * (
        log(2 * pi * variance) + residual^2 / variance
      )
    }
  }

  expect_identical(
    .log_lik_known_v_estimate_target_from_setup(setup),
    expected
  )
  expect_equal(
    .log_lik_known_v_joint_sum_from_setup(setup),
    rowSums(expected),
    tolerance = 1e-13
  )

  block_data <- .known_v_dependency_block_data(data, setup[["K"]])
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(setup[["yi"]]),
    .marglik_known_v_covariance_matrix(.data_known_v_data(data)),
    list(),
    lapply(block_data, `[[`, "index"),
    PACKAGE = "RoBMA"
  )
  expect_identical(attr(plan, "low_rank_blocks"), 0L)
  expect_identical(attr(plan, "root_dense_blocks"), 4L)
})


test_that("native covariance plan returns exact Schur conditional densities", {

  y <- c(0.1, -0.2, 0.3, -0.1, 0.05)
  mean <- c(0.02, -0.03, 0.08, -0.04, 0.01)
  sampling_covariance <- matrix(0, nrow = 5L, ncol = 5L)
  sampling_covariance[1:3, 1:3] <- matrix(
    c(0.08, 0.02, 0.01, 0.02, 0.09, 0.015, 0.01, 0.015, 0.07),
    nrow = 3L,
    byrow = TRUE
  )
  sampling_covariance[4:5, 4:5] <- matrix(
    c(0.06, 0.012, 0.012, 0.05),
    nrow = 2L
  )
  extra_variance <- c(0.01, 0.02, 0.015, 0.025, 0.018)
  blocks <- list(1:3, 4:5)
  plan <- .Call(
    "RoBMA_known_v_covariance_plan_create",
    as.double(y),
    sampling_covariance,
    list(),
    blocks,
    PACKAGE = "RoBMA"
  )
  actual <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    plan,
    as.double(mean),
    list(),
    as.double(extra_variance),
    PACKAGE = "RoBMA"
  )
  expected <- numeric(length(y))
  for (idx in blocks) {
    covariance <- sampling_covariance[idx, idx, drop = FALSE] +
      diag(extra_variance[idx], nrow = length(idx))
    expected[idx] <- .log_lik_known_v_component_conditional(
      yi         = y[idx],
      mu         = mean[idx],
      covariance = covariance
    )
  }

  expect_equal(actual, expected, tolerance = 1e-12)

  means <- rbind(mean, mean + c(0.01, -0.02, 0.03, -0.01, 0.02))
  extra_variances <- rbind(extra_variance, extra_variance * 1.25)
  states <- rep(list(list()), nrow(means))
  summary <- .marglik_covariance_plan_conditional_summary_batch(
    cache                    = NULL,
    y                        = y,
    means                    = means,
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = list(),
    random_covariance_states = states,
    block_indices            = blocks,
    extra_variances          = extra_variances
  )
  actual_batch <- .marglik_covariance_plan_conditional_loglik_batch(
    cache                    = NULL,
    y                        = y,
    means                    = means,
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = list(),
    random_covariance_states = states,
    block_indices            = blocks,
    extra_variances          = extra_variances
  )
  expected_residual <- expected_variance <- matrix(
    NA_real_,
    nrow = nrow(means),
    ncol = length(y)
  )
  for (draw in seq_len(nrow(means))) {
    for (idx in blocks) {
      distribution <- .known_v_component_conditional_distribution(
        yi         = y[idx],
        mu         = means[draw, idx],
        covariance = sampling_covariance[idx, idx, drop = FALSE] +
          diag(extra_variances[draw, idx], nrow = length(idx))
      )
      expected_residual[draw, idx] <- distribution[["residual"]]
      expected_variance[draw, idx] <- distribution[["variance"]]
    }
  }
  expected_batch <- -0.5 * (
    log(2 * pi * expected_variance) +
      expected_residual^2 / expected_variance
  )

  expect_equal(summary[["residual"]], expected_residual, tolerance = 1e-12)
  expect_equal(summary[["variance"]], expected_variance, tolerance = 1e-12)
  expect_equal(actual_batch, expected_batch, tolerance = 1e-12)
})


test_that("known-V factor plans recover batched precision right-hand sides", {

  sampling_covariance <- matrix(c(
    1.4, 0.2, 0.1, 0.0,
    0.2, 1.1, 0.3, 0.1,
    0.1, 0.3, 1.6, 0.2,
    0.0, 0.1, 0.2, 1.3
  ), nrow = 4L, byrow = TRUE)
  extra_variances <- rbind(
    c(0.10, 0.20, 0.15, 0.05),
    c(0.30, 0.05, 0.25, 0.10)
  )
  rhs <- cbind(c(0.3, -0.5, 0.8, 0.1), 1, c(-1, -0.2, 0.6, 1.4))
  plan_data <- list(
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = list(),
    random_covariance_states = rep(list(list()), nrow(extra_variances)),
    block_indices            = list(seq_len(nrow(rhs))),
    extra_variances          = extra_variances
  )

  observed <- .known_v_covariance_plan_precision_rhs_batch(
    plan_data = plan_data,
    rhs       = rhs
  )
  expected <- array(NA_real_, dim = dim(observed))
  for (draw in seq_len(nrow(extra_variances))) {
    covariance <- sampling_covariance + diag(extra_variances[draw, ])
    expected[draw, , ] <- solve(covariance, rhs)
  }

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("latent known-V metadata use the exact low-rank covariance plan", {

  loading <- c(0.2, 0.3, -0.1, 0.4)
  sampling_covariance <- tcrossprod(loading)
  known_V <- .known_v_prepare(
    V                         = sampling_covariance,
    keep_rows                 = rep(TRUE, length(loading)),
    known_v_parameterization  = "auto",
    known_v_residual_fraction = NULL,
    warn_singular             = FALSE
  )
  factor <- .known_v_latent_sampling_factor_plan(known_V)
  extra_variances <- rbind(
    c(0.10, 0.20, 0.15, 0.05),
    c(0.30, 0.05, 0.25, 0.10)
  )
  rhs <- cbind(c(0.3, -0.5, 0.8, 0.1), 1)
  plan_data <- list(
    sampling_covariance      = factor[["sampling_covariance"]],
    random_covariance_plans  = list(factor[["factor_plan"]]),
    random_covariance_states = rep(
      list(list(factor[["factor_state"]])),
      nrow(extra_variances)
    ),
    block_indices   = list(seq_along(loading)),
    extra_variances = extra_variances
  )

  observed <- .known_v_covariance_plan_precision_rhs_batch(
    plan_data = plan_data,
    rhs       = rhs
  )
  expected <- array(NA_real_, dim = dim(observed))
  for (draw in seq_len(nrow(extra_variances))) {
    covariance <- sampling_covariance + diag(extra_variances[draw, ])
    expected[draw, , ] <- solve(covariance, rhs)
  }

  expect_equal(
    tcrossprod(factor[["factor_plan"]][["model_matrix"]]),
    sampling_covariance,
    tolerance = 0
  )
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_null(.known_v_latent_sampling_factor_plan(.known_v_prepare(
    V                         = sampling_covariance + diag(0.01, 4L),
    keep_rows                 = rep(TRUE, 4L),
    known_v_parameterization  = "block_mvn",
    known_v_residual_fraction = NULL,
    warn_singular             = FALSE
  )))
})


test_that("known-V factor-plan GLS matches exact dense projections", {

  sampling_covariance <- matrix(c(
    1.4, 0.2, 0.1, 0.0,
    0.2, 1.1, 0.3, 0.1,
    0.1, 0.3, 1.6, 0.2,
    0.0, 0.1, 0.2, 1.3
  ), nrow = 4L, byrow = TRUE)
  extra_variances <- rbind(
    c(0.10, 0.20, 0.15, 0.05),
    c(0.30, 0.05, 0.25, 0.10)
  )
  y <- c(0.3, -0.5, 0.8, 0.1)
  X <- cbind(1, c(-1, -0.2, 0.6, 1.4))
  plan_data <- list(
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = list(),
    random_covariance_states = rep(list(list()), nrow(extra_variances)),
    block_indices            = list(seq_len(length(y))),
    extra_variances          = extra_variances,
    covariance_diagonal      = sweep(
      extra_variances,
      2L,
      diag(sampling_covariance),
      "+"
    )
  )
  testthat::local_mocked_bindings(
    .known_v_marginal_factor_plan = function(...) plan_data,
    .package = "RoBMA"
  )

  observed <- .known_v_factor_gls_projection_batch(
    object            = list(),
    posterior_samples = matrix(0, nrow = 2L, ncol = 1L),
    known_V           = list(),
    X                 = X,
    y                 = y,
    return_full_H     = TRUE,
    return_se         = TRUE,
    return_resid      = TRUE
  )

  for (draw in seq_len(nrow(extra_variances))) {
    covariance <- sampling_covariance + diag(extra_variances[draw, ])
    expected <- .known_v_gls_projection(X, y, covariance)
    A <- diag(length(y)) - expected[["H"]]

    expect_equal(observed[["H"]][draw, , ], expected[["H"]],
                 tolerance = 1e-12)
    expect_equal(observed[["H_diag"]][draw, ], diag(expected[["H"]]),
                 tolerance = 1e-12)
    expect_equal(observed[["M_diag"]][draw, ], diag(covariance),
                 tolerance = 1e-12)
    expect_equal(observed[["residual"]][draw, ], expected[["residual"]],
                 tolerance = 1e-12)
    expect_equal(
      observed[["residual_variance"]][draw, ],
      rowSums((A %*% expected[["covariance_factor"]])^2),
      tolerance = 1e-12
    )
  }
})


test_that("rank-one known V retains sub-ULP diagonal variance", {

  V    <- matrix(1, nrow = 2L, ncol = 2L)
  data <- list(outcome = data.frame(yi = c(0.10, -0.20), sei = c(1, 1)))
  attr(data, "known_V")      <- TRUE
  attr(data, "known_V_data") <- .known_v_prepare(
    V                         = V,
    keep_rows                 = rep(TRUE, 2L),
    known_v_parameterization  = "block_mvn",
    known_v_residual_fraction = NULL,
    warn_singular             = FALSE
  )
  attr(data, "random") <- FALSE

  diagonal <- rep(1e-18, 2L)
  setup <- list(
    outcome_type      = "norm",
    is_weightfunction = FALSE,
    weights           = NULL,
    data              = data,
    K                 = 2L,
    S                 = 1L,
    yi                = c(0.10, -0.20),
    mu                = matrix(0, nrow = 1L, ncol = 2L),
    tau_within        = matrix(sqrt(diagonal), nrow = 1L),
    effect_direction  = "positive",
    posterior_samples = matrix(numeric(0), nrow = 1L, ncol = 0L),
    marginalized_random_source_samples = NULL
  )

  determinant <- diagonal[[1L]] * (2 + diagonal[[1L]])
  residual    <- setup[["yi"]]
  quadratic   <- (
    (1 + diagonal[[1L]]) * sum(residual^2) -
      2 * prod(residual)
  ) / determinant
  expected_joint <- -0.5 * (
    2 * log(2 * pi) + log(determinant) + quadratic
  )
  expected_variance <- diagonal[[1L]] *
    (2 + diagonal[[1L]]) / (1 + diagonal[[1L]])

  expect_equal(
    .log_lik_known_v_joint_sum_from_setup(setup),
    expected_joint,
    tolerance = 1e-12
  )
  distribution <- .known_v_estimate_target_summary_from_setup(setup)
  expect_equal(
    distribution[["variance"]],
    matrix(expected_variance, nrow = 1L, ncol = 2L),
    tolerance = 1e-12
  )
  expect_true(all(is.finite(distribution[["log_lower"]])))
  expect_true(all(is.finite(distribution[["log_upper"]])))

  unequal_diagonal <- c(1e-18, 2e-18)
  setup[["tau_within"]] <- matrix(sqrt(unequal_diagonal), nrow = 1L)
  expected <- .known_v_diagonal_rank_one_conditional(
    yi       = setup[["yi"]],
    mu       = setup[["mu"]][1L, ],
    diagonal = unequal_diagonal,
    rank_one = c(1, 1)
  )
  distribution <- .known_v_estimate_target_summary_from_setup(
    setup      = setup,
    components = c("mean", "variance")
  )
  expect_equal(
    distribution[["variance"]] / expected[["variance"]],
    matrix(1, nrow = 1L, ncol = 2L),
    tolerance = 1e-12
  )
  expect_equal(
    matrix(setup[["yi"]], nrow = 1L) - distribution[["mean"]],
    matrix(expected[["residual"]], nrow = 1L),
    tolerance = 1e-12
  )

  projection <- .known_v_gls_projection_blocks(
    X              = matrix(1, nrow = 2L, ncol = 1L),
    y              = setup[["yi"]],
    known_V        = .data_known_v_data(data),
    extra_variance = diagonal
  )
  expect_true(all(is.finite(unlist(projection))))
})


test_that("evaluated known-V random log-likelihood requires conditioned mu", {

  dat <- data.frame(
    yi    = c(.10, -.20),
    study = c("s1", "s2")
  )
  object <- brma.mv(
    yi                         = yi,
    V                          = diag(c(.04, .09)),
    random                     = ~ 1 | study,
    data                       = dat,
    measure                    = "GEN",
    marginalize_estimate_level = FALSE,
    prior_unit_information_sd  = 1,
    only_priors                = TRUE
  )

  mu                <- matrix(0, nrow = 1L, ncol = 2L)
  tau_within        <- matrix(0, nrow = 1L, ncol = 2L)
  posterior_samples <- matrix(numeric(0), nrow = 1L, ncol = 0L)

  expect_error(
    .log_lik_known_v_joint_sum_from_evaluated_predictors(
      fit                = object[["fit"]],
      data               = object[["data"]],
      priors             = object[["priors"]],
      mu_samples         = mu,
      tau_within_samples = tau_within,
      posterior_samples  = posterior_samples
    ),
    "sampled random effects included"
  )
  expect_silent(
    .log_lik_known_v_joint_sum_from_evaluated_predictors(
      fit                         = object[["fit"]],
      data                        = object[["data"]],
      priors                      = object[["priors"]],
      mu_samples                  = mu,
      tau_within_samples          = tau_within,
      posterior_samples           = posterior_samples,
      random_effects_conditioning = "included_in_mu"
    )
  )
})


test_that("evaluated known-V marginalized scale maps component row source", {

  dat <- data.frame(
    yi   = c(0.10, 0.20, -0.05, 0.15),
    type = factor(
      c("RCT", "RCT", "cohort", "cohort"),
      levels = c("RCT", "cohort")
    )
  )
  dat$study   <- factor(seq_len(nrow(dat)))
  dat$cohort  <- as.numeric(dat$type == "cohort")
  dat$bias_id <- factor("cohort_bias")
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    random                    = list(
      coh_bias    = ~ diag(0 + cohort | bias_id),
      ran_effects = ~ 1 | study
    ),
    scale                     = list(ran_effects = ~ type),
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  expect_equal(.data_sampled_random_effect_blocks(object[["data"]]), "coh_bias")

  mu_samples         <- matrix(0, nrow = 1L, ncol = 4L)
  tau_within_samples <- matrix(c(0.10, 0.20, 0.30, 0.40), nrow = 1L)
  posterior_samples  <- matrix(numeric(0), nrow = 1L, ncol = 0L)

  source_samples <- .known_v_marginalized_random_source_samples_from_tau(
    data               = object[["data"]],
    tau_within_samples = tau_within_samples
  )
  expect_named(source_samples, "tau_ran_effects")
  expect_equal(source_samples[["tau_ran_effects"]], tau_within_samples)

  formula_args <- .create_jags_formula_args(
    data   = object[["data"]],
    priors = object[["priors"]]
  )
  expect_true("tau_ran_effects" %in% formula_args[["add_parameters"]])

  log_lik <- .log_lik_known_v_joint_sum_from_evaluated_predictors(
    fit                         = object[["fit"]],
    data                        = object[["data"]],
    priors                      = object[["priors"]],
    mu_samples                  = mu_samples,
    tau_within_samples          = tau_within_samples,
    posterior_samples           = posterior_samples,
    random_effects_conditioning = "included_in_mu"
  )
  expected <- .marglik_mvn_log_density(
    y          = dat[["yi"]],
    mean       = rep(0, 4),
    covariance = diag(rep(0.04, 4) + as.numeric(tau_within_samples)^2)
  )

  expect_equal(log_lik, expected, tolerance = 1e-12)

  bridge_context <- structure(
    list(nodes = stats::setNames(
      as.numeric(tau_within_samples),
      paste0("tau_ran_effects[", seq_len(ncol(tau_within_samples)), "]")
    )),
    class = c("BayesTools_bridge_context", "list")
  )
  expect_equal(
    .log_posterior(
      parameters                 = list(mu = 0),
      data                       = .create_fit_data(object[["data"]], object[["priors"]]),
      is_mods                    = FALSE,
      is_scale                   = TRUE,
      is_random                  = TRUE,
      is_multilevel              = FALSE,
      is_weights                 = FALSE,
      is_known_v                 = TRUE,
      is_PET                     = FALSE,
      is_PEESE                   = FALSE,
      is_weightfunction          = FALSE,
      effect_direction           = "positive",
      outcome_type               = "norm",
      model_data                 = object[["data"]],
      bridge_context             = bridge_context
    ),
    expected,
    tolerance = 1e-12
  )
})


test_that("IWMDE evaluated known-V likelihood matches joint MVN oracle", {

  V <- matrix(c(0.04, 0.015, 0.015, 0.09), nrow = 2L)
  object <- brma.mv(
    yi                         = c(0.10, -0.20),
    V                          = V,
    random                     = ~ 1 | study,
    data                       = data.frame(
      yi    = c(0.10, -0.20),
      study = c("s1", "s2")
    ),
    known_v_parameterization   = "block_mvn",
    measure                    = "GEN",
    marginalize_estimate_level = FALSE,
    prior_unit_information_sd  = 1,
    only_priors                = TRUE
  )
  mu_samples  <- matrix(c(0.02, -0.10, 0.06, -0.16), nrow = 2L, byrow = TRUE)
  tau_samples <- matrix(c(0.05, 0.08, 0.04, 0.06), nrow = 2L, byrow = TRUE)
  context <- list(
    object = object,
    data   = object[["data"]]
  )
  active_setup <- list(
    priors            = object[["priors"]],
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    is_weightfunction = FALSE
  )

  posterior_samples <- matrix(numeric(0), nrow = 2L, ncol = 0L)
  setup <- .log_lik_evaluated_setup(
    fit                         = object[["fit"]],
    data                        = object[["data"]],
    priors                      = object[["priors"]],
    unit                        = "estimate",
    data_hash                   = NULL,
    mu_samples                  = mu_samples,
    tau_within_samples          = tau_samples,
    tau_between_samples         = NULL,
    posterior_samples           = posterior_samples,
    random_effects_conditioning = "included_in_mu"
  )
  expected <- vapply(seq_len(nrow(mu_samples)), function(s) {
    .marglik_mvn_log_density(
      y          = object[["data"]][["outcome"]][["yi"]],
      mean       = mu_samples[s, ],
      covariance = V
    )
  }, numeric(1))

  expect_true(.iwmde_uses_known_v_joint_likelihood(context))
  expect_equal(
    .log_lik_known_v_joint_sum_from_setup(setup),
    expected,
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    rowSums(.log_lik_known_v_estimate_target_from_setup(setup)),
    expected,
    tolerance = 1e-8
  )))
  expect_equal(
    .iwmde_log_lik_from_evaluated_predictors_sum_active_branch(
      context            = context,
      active_setup       = active_setup,
      mu_samples         = mu_samples,
      tau_within_samples = tau_samples,
      posterior_samples  = posterior_samples,
      unit               = "estimate"
    ),
    expected,
    tolerance = 1e-12
  )
})
