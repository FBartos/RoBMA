# ============================================================================ #
# IWMDE Known-V Predictor Kernels
# ============================================================================ #

# Evaluate known-V tau candidates while reusing each block factorization across
# posterior rows.
.iwmde_log_q_grid_normal_known_v_tau_group <- function(
    context, parameter, values, row_states, replacement, setup, basis) {

  data <- context[["data"]]
  if (!identical(parameter, "tau") ||
      !.iwmde_uses_known_v_joint_likelihood(context) ||
      .is_data_random(data) ||
      .is_data_multilevel(data) ||
      isTRUE(setup[["is_weightfunction"]]) ||
      !is.null(setup[["weights"]]) ||
      !identical(basis[["scale_update"]], "tau") ||
      !is.null(basis[["mu_basis"]]) ||
      !is.null(basis[["log_tau_basis"]]) ||
      isTRUE(basis[["formula_mu"]]) ||
      isTRUE(basis[["formula_logtau"]])) {
    return(NULL)
  }

  known_V <- .data_known_v_data(data)
  G       <- length(values)
  S       <- length(row_states)
  K       <- setup[["K"]]
  if (.known_v_nrow(known_V) != K ||
      !identical(dim(setup[["mu"]]), c(S, K))) {
    return(NULL)
  }

  log_prior <- .iwmde_predictor_log_prior(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )
  if (is.null(log_prior) || length(log_prior) != G * S) {
    return(NULL)
  }

  block_data        <- .known_v_dependency_block_data(data, K)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")
  if (any(vapply(
    block_covariances,
    function(covariance) {
      !is.null(.covariance_exact_rank_one_factor(covariance)) &&
        nrow(covariance) > 1L
    },
    logical(1)
  ))) {
    return(NULL)
  }

  residual <- matrix(
    setup[["yi"]],
    nrow  = S,
    ncol  = K,
    byrow = TRUE
  ) - setup[["mu"]]
  log_lik <- matrix(-Inf, nrow = G, ncol = S)
  valid   <- is.finite(values) & values >= 0

  for (g in which(valid)) {
    candidate_log_lik <- numeric(S)
    for (block in seq_along(block_indices)) {
      idx        <- block_indices[[block]]
      covariance <- block_covariances[[block]] +
        diag(values[[g]]^2, nrow = length(idx))
      chol_covariance <- .known_v_chol_covariance(
        covariance = covariance,
        context    = "joint likelihood"
      )
      residual_whitened <- backsolve(
        chol_covariance,
        t(residual[, idx, drop = FALSE]),
        transpose = TRUE
      )
      candidate_log_lik <- candidate_log_lik - .5 * (
        length(idx) * log(2 * pi) +
          2 * sum(log(diag(chol_covariance))) +
          colSums(residual_whitened^2)
      )
    }
    log_lik[g, ] <- candidate_log_lik
  }

  return(log_lik + matrix(log_prior, nrow = G, ncol = S))
}


.iwmde_normal_location_likelihood_change_known_v <- function(context, setup,
                                                              yi, mu,
                                                              mu_basis) {

  if (isTRUE(setup[["is_weightfunction"]]) ||
      !is.null(setup[["weights"]])) {
    return(NULL)
  }

  known_V <- .data_known_v_data(context[["data"]])
  K       <- setup[["K"]]
  S       <- setup[["S"]]
  if (.known_v_nrow(known_V) != K ||
      !identical(dim(mu), c(S, K)) ||
      !identical(dim(mu_basis), c(S, K))) {
    return(NULL)
  }

  block_data        <- .known_v_dependency_block_data(context[["data"]], K)
  block_indices     <- lapply(block_data, `[[`, "index")
  block_covariances <- lapply(block_data, `[[`, "covariance")
  if (any(vapply(
    block_covariances,
    function(covariance) {
      !is.null(.covariance_exact_rank_one_factor(covariance)) &&
        nrow(covariance) > 1L
    },
    logical(1)
  ))) {
    return(NULL)
  }

  extra_variance <- .known_v_extra_variance_from_setup(setup)
  linear         <- numeric(S)
  quadratic      <- numeric(S)

  for (block in seq_along(block_indices)) {
    idx        <- block_indices[[block]]
    covariance <- block_covariances[[block]]
    residual   <- matrix(yi[idx], nrow = S, ncol = length(idx), byrow = TRUE) -
      mu[, idx, drop = FALSE]
    basis      <- mu_basis[, idx, drop = FALSE]
    extra      <- extra_variance[, idx, drop = FALSE]

    invariant_extra <- S == 1L || all(extra == matrix(
      extra[1L, ],
      nrow  = S,
      ncol  = length(idx),
      byrow = TRUE
    ))
    if (invariant_extra) {
      covariance <- covariance + diag(extra[1L, ], nrow = length(idx))
      chol_covariance <- .known_v_chol_covariance(
        covariance = covariance,
        context    = "joint likelihood"
      )
      residual_whitened <- backsolve(
        chol_covariance,
        t(residual),
        transpose = TRUE
      )
      basis_whitened <- backsolve(
        chol_covariance,
        t(basis),
        transpose = TRUE
      )
      linear    <- linear + colSums(residual_whitened * basis_whitened)
      quadratic <- quadratic + colSums(basis_whitened^2)
      next
    }

    for (s in seq_len(S)) {
      row_covariance <- covariance +
        diag(extra[s, ], nrow = length(idx))
      chol_covariance <- .known_v_chol_covariance(
        covariance = row_covariance,
        context    = "joint likelihood"
      )
      residual_whitened <- backsolve(
        chol_covariance,
        residual[s, ],
        transpose = TRUE
      )
      basis_whitened <- backsolve(
        chol_covariance,
        basis[s, ],
        transpose = TRUE
      )
      linear[[s]]    <- linear[[s]] +
        sum(residual_whitened * basis_whitened)
      quadratic[[s]] <- quadratic[[s]] + sum(basis_whitened^2)
    }
  }

  if (any(!is.finite(linear)) ||
      any(!is.finite(quadratic)) ||
      any(quadratic < 0)) {
    return(NULL)
  }

  return(list(
    linear    = linear,
    quadratic = quadratic
  ))
}
