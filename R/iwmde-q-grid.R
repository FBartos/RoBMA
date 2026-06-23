# ============================================================================ #
# IWMDE Q-Grid Dispatch and GLMM Batches
# ============================================================================ #

.iwmde_log_q_grid <- function(context, parameter, values, row_states,
                              replacement) {

  out <- .iwmde_log_q_grid_predictor_batch(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )

  if (is.matrix(out) &&
      nrow(out) == length(values) &&
      ncol(out) == length(row_states)) {
    return(out)
  }

  out <- .iwmde_log_q_grid_glmm_conditional_batch(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )

  if (is.matrix(out) &&
      nrow(out) == length(values) &&
      ncol(out) == length(row_states)) {
    return(out)
  }

  out <- .iwmde_log_q_grid_batch(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )

  if (is.matrix(out) &&
      nrow(out) == length(values) &&
      ncol(out) == length(row_states)) {
    return(out)
  }

  return(.iwmde_log_q_grid_scalar(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  ))
}


.iwmde_log_q_grid_glmm_conditional_batch <- function(context, parameter, values,
                                                     row_states, replacement) {

  if (!.data_outcome_type(context[["data"]]) %in% c("bin", "pois")) {
    return(NULL)
  }

  return(.iwmde_log_q_grid_from_samples(
    context         = context,
    parameter       = parameter,
    values          = values,
    row_states      = row_states,
    replacement     = replacement,
    likelihood_mode = "conditional",
    log_lik_fun     = function(samples, active_setup) {
      .iwmde_glmm_conditional_log_likelihood_samples(
        context      = context,
        samples      = samples,
        active_setup = active_setup
      )
    }
  ))
}


.iwmde_glmm_conditional_log_likelihood_samples <- function(context, samples,
                                                           active_setup) {

  data         <- context[["data"]]
  outcome_type <- .data_outcome_type(data)
  if (!outcome_type %in% c("bin", "pois")) {
    return(NULL)
  }

  K             <- nrow(data[["outcome"]])
  is_mods       <- .is_data_mods(data)
  is_scale      <- .is_data_scale(data)
  is_multilevel <- .is_data_multilevel(data)

  mu_samples <- .evaluate.brma.mu(
    fit               = context[["object"]][["fit"]],
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (is_mods) .create_fit_formula_list(data = data, "mods") else NULL,
    mods_priors       = active_setup[["priors"]][["mods"]],
    is_mods           = is_mods,
    is_PET            = active_setup[["is_PET"]],
    is_PEESE          = active_setup[["is_PEESE"]],
    effect_direction  = .data_effect_direction(data),
    bias_adjusted     = FALSE,
    K                 = K,
    posterior_samples = samples,
    priors            = active_setup[["priors"]]
  )
  tau_result <- .evaluate.brma.tau(
    fit               = context[["object"]][["fit"]],
    scale_data        = data[["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = data, "scale") else NULL,
    scale_priors      = active_setup[["priors"]][["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = samples,
    allow_missing_tau = .fixed_tau_prior_value(active_setup[["priors"]])
  )

  return(.log_lik_glmm_conditional_sum_from_evaluated_predictors(
    fit                 = context[["object"]][["fit"]],
    data                = data,
    posterior_samples   = samples,
    mu_samples          = mu_samples,
    tau_within_samples  = tau_result[["tau_within"]],
    tau_between_samples = tau_result[["tau_between"]]
  ))
}


.iwmde_log_q_grid_predictor_batch <- function(context, parameter, values,
                                              row_states, replacement) {

  if (length(values) == 0L || length(row_states) == 0L) {
    return(matrix(-Inf, nrow = length(values), ncol = length(row_states)))
  }

  likelihood_modes <- vapply(row_states, function(state) {
    state[["likelihood_mode"]]
  }, character(1))
  if (!all(likelihood_modes == "marginal")) {
    return(NULL)
  }

  out  <- matrix(-Inf, nrow = length(values), ncol = length(row_states))
  keys <- vapply(row_states, function(state) {
    .iwmde_state_active_key(context, state)
  }, character(1))

  for (key in unique(keys)) {
    state_cols   <- which(keys == key)
    group_states <- row_states[state_cols]
    group_out    <- .iwmde_log_q_grid_predictor_group(
      context     = context,
      parameter   = parameter,
      values      = values,
      row_states  = group_states,
      replacement = replacement
    )

    if (!is.matrix(group_out)) {
      return(NULL)
    }

    out[, state_cols] <- group_out
  }

  return(out)
}


.iwmde_log_q_grid_predictor_group <- function(context, parameter, values,
                                              row_states, replacement) {

  active_setup <- row_states[[1L]][["active_setup"]]
  unit         <- if (.is_data_multilevel(context[["data"]])) {
    "cluster"
  } else {
    "estimate"
  }
  setup <- .iwmde_predictor_setup(
    context      = context,
    row_states   = row_states,
    active_setup = active_setup,
    unit         = unit
  )
  basis <- .iwmde_predictor_update_basis(
    context     = context,
    parameter   = parameter,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup
  )
  if (is.null(basis)) {
    return(NULL)
  }
  basis <- .iwmde_predictor_resolve_formula_basis(
    context      = context,
    active_setup = active_setup,
    setup        = setup,
    basis        = basis,
    row_states   = row_states
  )

  normal_out <- .iwmde_log_q_grid_normal_location_group(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  if (is.matrix(normal_out)) {
    return(normal_out)
  }

  candidates <- .iwmde_predictor_candidates(
    context     = context,
    active_setup = active_setup,
    setup       = setup,
    basis       = basis,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )
  posterior_samples <- candidates[["posterior_samples"]]
  log_prior <- .iwmde_predictor_log_prior(
    context             = context,
    parameter           = parameter,
    values              = values,
    row_states          = row_states,
    replacement         = replacement,
    replacement_samples = candidates[["replacement_samples"]]
  )
  if (is.null(log_prior) || length(log_prior) != nrow(candidates[["mu"]])) {
    return(NULL)
  }

  log_lik <- .iwmde_log_lik_from_evaluated_predictors_sum_active_branch(
    context             = context,
    active_setup        = active_setup,
    mu_samples          = candidates[["mu"]],
    tau_within_samples  = candidates[["tau_within"]],
    tau_between_samples = candidates[["tau_between"]],
    posterior_samples   = posterior_samples,
    unit                = unit
  )
  log_q <- log_lik + log_prior
  log_q[!candidates[["valid"]]] <- -Inf

  return(matrix(log_q, nrow = length(values), ncol = length(row_states)))
}
