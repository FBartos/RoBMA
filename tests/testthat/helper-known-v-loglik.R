.test_log_lik_known_v_joint_sum_from_evaluated_predictors <- function(
    fit, data, priors, mu_samples, tau_within_samples,
    tau_between_samples = NULL, posterior_samples = NULL, unit = "estimate",
    data_hash = NULL, random_effects_conditioning = "none") {

  unit  <- .normalize_unit(unit)
  setup <- .log_lik_evaluated_setup(
    fit                         = fit,
    data                        = data,
    priors                      = priors,
    unit                        = unit,
    data_hash                   = data_hash,
    mu_samples                  = mu_samples,
    tau_within_samples          = tau_within_samples,
    tau_between_samples         = tau_between_samples,
    posterior_samples           = posterior_samples,
    random_effects_conditioning = random_effects_conditioning
  )

  return(.log_lik_known_v_joint_sum_from_setup(setup))
}
