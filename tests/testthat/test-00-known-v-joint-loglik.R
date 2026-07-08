# test-00-known-v-joint-loglik.R

test_that("known-V joint log-likelihood uses block MVN density", {

  V    <- matrix(c(.04, .015, .015, .09), nrow = 2L)
  data <- list(outcome = data.frame(yi = c(.10, -.20), sei = c(.20, .30)))
  attr(data, "known_V")      <- TRUE
  attr(data, "known_V_data") <- list(
    V             = V,
    block_indices = list(1:2)
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
  expect_false(isTRUE(all.equal(
    out,
    rowSums(.log_lik_known_v_estimate_target_from_setup(setup))
  )))
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
    yi       = c(0.10, 0.20, -0.05, 0.15),
    estimate = paste0("e", 1:4),
    x        = c(0, 1, 0, 1)
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, 4)),
    random                    = list(effect = ~ 1 | estimate),
    scale                     = list(effect = ~ x),
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  mu_samples         <- matrix(0, nrow = 1L, ncol = 4L)
  tau_within_samples <- matrix(c(0.10, 0.20, 0.30, 0.40), nrow = 1L)
  posterior_samples  <- matrix(numeric(0), nrow = 1L, ncol = 0L)

  source_samples <- .known_v_marginalized_random_source_samples_from_tau(
    data               = object[["data"]],
    tau_within_samples = tau_within_samples
  )
  expect_named(source_samples, "tau_effect")
  expect_equal(source_samples[["tau_effect"]], tau_within_samples)

  formula_args <- .create_jags_formula_args(
    data   = object[["data"]],
    priors = object[["priors"]]
  )
  expect_true("tau_effect" %in% formula_args[["add_parameters"]])

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
      paste0("tau_effect[", seq_len(ncol(tau_within_samples)), "]")
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
      known_v_parameterization   = "block_mvn",
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
