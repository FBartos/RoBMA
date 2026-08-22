# ============================================================================ #
# Exact coefficient-state random-covariance q-grid
# ============================================================================ #

.iwmde_log_q_grid_known_v_random_factor <- function(
    context, parameter, values, row_states, replacement, active_setup) {

  update <- replacement[["covariance_update"]]
  if (!.iwmde_uses_known_v_random_marginal_likelihood(context) ||
      !inherits(
        update,
        "BayesTools_random_effects_marginal_update_plan"
      ) || !update[["family"]] %in% c("factor", "markov") ||
      !replacement[["type"]] %in% c("primitive", "scalar") ||
      .is_data_weights(context[["data"]])) {
    return(NULL)
  }

  finite <- is.finite(values)
  if (sum(finite) < 2L) {
    return(NULL)
  }
  setup <- .iwmde_known_v_random_marginal_setup(context)
  rows  <- vapply(row_states, `[[`, integer(1), "row_index")
  samples <- context[["posterior_samples"]][rows, , drop = FALSE]
  mu_samples <- .iwmde_predictor_evaluate_fixed_mu(
    context      = context,
    active_setup = active_setup,
    samples      = samples
  )
  data <- context[["data"]]
  yi   <- data[["outcome"]][["yi"]]
  if (.data_effect_direction(data) == "negative") {
    yi         <- -yi
    mu_samples <- -mu_samples
  }

  inputs <- .brma_mv_random_effects_marginal_inputs(
    object            = context[["object"]],
    posterior_samples = samples,
    data              = data
  )
  random_factors <- .brma_mv_random_effects_marginal_factor_states(
    object            = context[["object"]],
    posterior_samples = samples,
    blocks            = setup[["blocks"]],
    row_blocks        = setup[["dependency_blocks"]],
    data              = data,
    inputs            = inputs
  )
  update_grid <- BayesTools::random_effects_marginal_update_grid(
    fit               = inputs[["formula_fit"]],
    update            = update,
    values            = values[finite],
    posterior_samples = samples,
    prior_list        = inputs[["location_priors"]]
  )
  factor_names <- names(random_factors[["factor_plans"]])
  factor_index <- match(update_grid[["block"]], factor_names)
  if (is.na(factor_index)) {
    stop(
      "Random-effect update grid does not match the compiled factor plan.",
      call. = FALSE
    )
  }
  update_grid[["factor_index"]] <- as.integer(factor_index)

  S <- length(row_states)
  K <- length(yi)
  evaluated <- .marglik_covariance_plan_factor_grid_loglik(
    cache                    = setup[["covariance_plan_cache"]],
    y                        = as.double(yi),
    means                    = mu_samples,
    sampling_covariance      = setup[["sampling_covariance"]],
    random_covariance_plans  = random_factors[["factor_plans"]],
    random_covariance_states = random_factors[["factor_states"]],
    block_indices            = setup[["dependency_blocks"]],
    extra_variances          = matrix(0, nrow = S, ncol = K),
    update_grid              = update_grid
  )
  if (!is.matrix(evaluated) ||
      !identical(dim(evaluated), c(sum(finite), S))) {
    return(NULL)
  }
  out <- matrix(-Inf, nrow = length(values), ncol = S)
  out[finite, ] <- evaluated

  log_prior <- .iwmde_predictor_log_prior(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )
  if (is.null(log_prior) || length(log_prior) != length(values) * S) {
    return(NULL)
  }

  out + matrix(log_prior, nrow = length(values), ncol = S)
}
