# ============================================================================ #
# Exact affine random-covariance q-grid
# ============================================================================ #

.iwmde_log_q_grid_known_v_random_affine <- function(
    context, parameter, values, row_states, replacement, active_setup) {

  update <- replacement[["covariance_update"]]
  coefficient_input <- if (identical(
    replacement[["type"]],
    "random_component_sd"
  )) {
    "quantity"
  } else {
    "source"
  }
  if (!.iwmde_uses_known_v_random_marginal_likelihood(context) ||
      !inherits(
        update,
        "BayesTools_random_effects_marginal_update_plan"
      ) || !identical(update[["family"]], "affine") ||
      !identical(update[["coefficient_input"]], coefficient_input) ||
      .is_data_weights(context[["data"]])) {
    return(NULL)
  }

  setup <- .iwmde_known_v_random_marginal_setup(context)

  coefficients <- .iwmde_random_affine_coefficients(values, update)
  finite <- is.finite(values) & is.finite(coefficients)
  if (sum(finite) < 2L ||
      length(unique(coefficients[finite])) < 2L) {
    return(NULL)
  }

  rows <- vapply(row_states, `[[`, integer(1), "row_index")
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

  S     <- length(row_states)
  K     <- length(yi)
  out   <- matrix(-Inf, nrow = length(values), ncol = S)
  invariant <- .iwmde_random_affine_invariant_covariance(
    update = update,
    setup  = setup,
    K      = K
  )
  if (!is.null(invariant)) {
    evaluated <- .iwmde_random_affine_log_likelihood_chunk(
      base_covariance       = invariant[["base_covariance"]],
      update_covariance     = invariant[["update_covariance"]],
      reference_coefficient = invariant[["reference_coefficient"]],
      coefficients          = coefficients[finite],
      means                 = mu_samples,
      outcome               = yi,
      blocks                = setup[["dependency_blocks"]],
      plan_cache            = setup[["affine_plan_cache"]],
      sampling_covariance   = setup[["sampling_covariance"]]
    )
    if (!is.matrix(evaluated) ||
        !identical(dim(evaluated), c(sum(finite), S))) {
      return(NULL)
    }
    out[finite, ] <- evaluated
  } else {
    finite_rows <- which(finite)
    anchor_rows <- finite_rows[c(
      which.min(coefficients[finite]),
      which.max(coefficients[finite])
    )]
    anchor_values <- values[anchor_rows]
    anchor_coefficients <- coefficients[anchor_rows]
    if (anchor_coefficients[[1L]] == anchor_coefficients[[2L]]) {
      return(NULL)
    }
    max_bytes <- .known_v_covariance_max_bytes()
    chunk_bytes <- if (is.infinite(max_bytes)) Inf else max_bytes / 2
    chunks <- .known_v_covariance_chunk_indices(
      S         = S,
      K         = K,
      max_bytes = chunk_bytes
    )

    for (chunk in chunks) {
      anchor_samples <- lapply(anchor_values, function(value) {
        .iwmde_random_affine_replacement_samples(
          context     = context,
          parameter   = parameter,
          value       = value,
          row_states  = row_states[chunk],
          replacement = replacement
        )
      })
      if (any(vapply(anchor_samples, is.null, logical(1)))) {
        return(NULL)
      }
      covariance <- lapply(anchor_samples, function(anchor_sample) {
        result <- .brma_mv_random_effects_marginal_vcov(
          object            = context[["object"]],
          posterior_samples = anchor_sample
        )
        result[["samples"]]
      })
      expected_dim <- c(length(chunk), K, K)
      valid_covariance <- vapply(covariance, function(value) {
        is.array(value) && identical(dim(value), expected_dim) &&
          all(is.finite(value))
      }, logical(1))
      if (!all(valid_covariance)) {
        return(NULL)
      }
      coefficient_difference <-
        anchor_coefficients[[2L]] - anchor_coefficients[[1L]]
      update_covariance <-
        (covariance[[2L]] - covariance[[1L]]) /
        coefficient_difference
      base_covariance <- covariance[[1L]]
      for (j in seq_along(chunk)) {
        base_covariance[j, , ] <-
          base_covariance[j, , ] + setup[["sampling_covariance"]]
      }
      evaluated <- .iwmde_random_affine_log_likelihood_chunk(
        base_covariance       = base_covariance,
        update_covariance     = update_covariance,
        reference_coefficient = anchor_coefficients[[1L]],
        coefficients          = coefficients[finite],
        means                 = mu_samples[chunk, , drop = FALSE],
        outcome               = yi,
        blocks                = setup[["dependency_blocks"]],
        plan_cache            = setup[["affine_plan_cache"]],
        sampling_covariance   = setup[["sampling_covariance"]]
      )
      if (!is.matrix(evaluated) ||
          !identical(dim(evaluated), c(sum(finite), length(chunk)))) {
        return(NULL)
      }
      out[finite, chunk] <- evaluated
    }
  }

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


.iwmde_random_affine_invariant_covariance <- function(update, setup, K) {

  invariant <- update[["invariant_covariance"]]
  if (!is.list(invariant) ||
      !all(c(
        "reference_coefficient", "base_covariance", "update_covariance"
      ) %in% names(invariant))) {
    return(NULL)
  }
  reference <- invariant[["reference_coefficient"]]
  base      <- invariant[["base_covariance"]]
  basis     <- invariant[["update_covariance"]]
  sampling  <- setup[["sampling_covariance"]]
  expected  <- c(K, K)
  if (!is.numeric(reference) || length(reference) != 1L ||
      !is.finite(reference) || !is.matrix(base) || !is.matrix(basis) ||
      !is.matrix(sampling) || !identical(dim(base), expected) ||
      !identical(dim(basis), expected) || !identical(dim(sampling), expected) ||
      any(!is.finite(base)) || any(!is.finite(basis)) ||
      any(!is.finite(sampling)) || !identical(base, t(base)) ||
      !identical(basis, t(basis)) || !identical(sampling, t(sampling))) {
    return(NULL)
  }

  list(
    reference_coefficient = reference,
    base_covariance       = base + sampling,
    update_covariance     = basis
  )
}


.iwmde_random_affine_coefficients <- function(values, update) {

  transform <- update[["coefficient_transform"]]
  if (is.null(transform)) {
    return(rep(NA_real_, length(values)))
  }
  tryCatch(
    BayesTools::parameter_transform_forward(values, transform),
    error = function(e) rep(NA_real_, length(values))
  )
}


.iwmde_random_affine_replacement_samples <- function(
    context, parameter, value, row_states, replacement) {

  rows <- lapply(row_states, function(state) {
    replaced <- .iwmde_replace_row_for_value(
      context     = context,
      state       = state,
      parameter   = parameter,
      value       = value,
      replacement = replacement
    )
    if (!isTRUE(replaced[["valid"]])) {
      return(NULL)
    }
    replaced[["row"]]
  })
  if (any(vapply(rows, is.null, logical(1)))) {
    return(NULL)
  }
  out <- do.call(rbind, rows)
  storage.mode(out) <- "double"
  out
}


.iwmde_random_affine_log_likelihood_chunk <- function(
    base_covariance, update_covariance, reference_coefficient, coefficients,
    means, outcome, blocks, plan_cache, sampling_covariance) {

  draw_specific <- length(dim(base_covariance)) == 3L
  if (draw_specific != (length(dim(update_covariance)) == 3L)) {
    return(NULL)
  }

  .marglik_covariance_plan_affine_grid_loglik(
    cache                 = plan_cache,
    y                     = as.double(outcome),
    means                 = means,
    sampling_covariance   = sampling_covariance,
    block_indices         = blocks,
    base_covariances      = base_covariance,
    update_covariances    = update_covariance,
    reference_coefficient = reference_coefficient,
    coefficients          = coefficients
  )
}
