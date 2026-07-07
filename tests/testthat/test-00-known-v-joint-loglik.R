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
