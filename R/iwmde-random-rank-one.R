# ============================================================================ #
# Exact coefficient-state rank-one variance q-grid
# ============================================================================ #
# Joint-selection counterpart of the affine mean-shift fast path: when the
# replacement varies a single random component's SD with every other
# coordinate fixed, the dense-block selection normalizers are interpolated
# along the certified rank-one covariance family instead of being re-integrated
# at every grid point. See `R/selection-covariance-grid.R` for the derivation
# and the error certificates.

.iwmde_joint_rank_one_variance_log_likelihood <- function(
    context, parameter, values, samples, active_setup, batch, replacement) {

  shared <- context[["covariance_grid"]]
  if (!is.environment(shared) || isTRUE(shared$untracked)) {
    return(NULL)
  }
  update <- shared[["update"]]
  if (!inherits(update, "BayesTools_random_effects_marginal_update_plan")) {
    return(NULL)
  }
  candidates <- batch[["candidates"]]
  positions <- batch[["valid_positions"]]
  state_index <- candidates[["state_index"]][positions]
  grid_index <- candidates[["grid_index"]][positions]
  row_values <- values[grid_index]
  if (any(!is.finite(row_values)) || any(row_values < 0)) {
    return(NULL)
  }
  # The certified family is parameterized by the component SD; the affine
  # coefficient must be its exact square for the transport bound.
  coefficients <- tryCatch(
    BayesTools::parameter_transform_forward(
      row_values, update[["coefficient_transform"]]),
    error = function(e) rep(NA_real_, length(row_values))
  )
  if (any(!is.finite(coefficients)) ||
      any(coefficients != row_values * row_values)) {
    return(NULL)
  }

  conditioned_random_effects <- .iwmde_conditioned_random_effects_from_latent(
    context           = context,
    posterior_samples = samples,
    unit              = "estimate"
  )
  setup <- .log_lik_posterior_setup(
    fit                        = context[["object"]][["fit"]],
    posterior_samples          = samples,
    data                       = context[["data"]],
    priors                     = active_setup[["priors"]],
    unit                       = "estimate",
    conditioned_random_effects = conditioned_random_effects
  )
  setup[["covariance_grid"]] <- list(
    shared       = shared,
    state_index  = state_index,
    values       = row_values,
    row_index    = vapply(batch[["row_states"]][state_index], `[[`, integer(1), "row_index"),
    context      = context,
    active_setup = active_setup
  )

  result <- tryCatch(
    .log_lik_known_v_joint_sum_from_setup(setup),
    error = function(e) {
      shared$untracked <- TRUE
      NULL
    }
  )
  if (is.null(result) || length(result) != nrow(samples) ||
      any(!is.finite(result))) {
    return(NULL)
  }
  result
}
