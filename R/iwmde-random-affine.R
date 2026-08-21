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
  if (!.iwmde_random_affine_has_sampling_dependence(setup)) {
    return(NULL)
  }

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
      blocks                = setup[["dependency_blocks"]]
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
        blocks                = setup[["dependency_blocks"]]
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


.iwmde_random_affine_has_sampling_dependence <- function(setup) {

  covariance <- setup[["sampling_covariance"]]
  if (!is.matrix(covariance) || nrow(covariance) != ncol(covariance)) {
    return(FALSE)
  }
  off_diagonal <- row(covariance) != col(covariance)
  any(covariance[off_diagonal] != 0)
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
    means, outcome, blocks) {

  draw_specific <- length(dim(base_covariance)) == 3L
  if (draw_specific != (length(dim(update_covariance)) == 3L)) {
    return(NULL)
  }
  S <- nrow(means)
  G <- length(coefficients)
  out <- matrix(0, nrow = G, ncol = S)
  invariant <- !draw_specific ||
    .iwmde_random_affine_covariance_invariant(
      base_covariance,
      update_covariance
    )

  for (block in blocks) {
    shared_plan <- NULL
    if (invariant) {
      shared_plan <- .iwmde_random_affine_spectral_plan(
        base_covariance = .iwmde_random_affine_covariance_block(
          base_covariance,
          1L,
          block
        ),
        update_covariance = .iwmde_random_affine_covariance_block(
          update_covariance,
          1L,
          block
        ),
        reference_coefficient = reference_coefficient,
        coefficients          = coefficients
      )
      if (is.null(shared_plan)) {
        return(NULL)
      }
    }
    for (s in seq_len(S)) {
      plan <- shared_plan
      if (is.null(plan)) {
        plan <- .iwmde_random_affine_spectral_plan(
          base_covariance = .iwmde_random_affine_covariance_block(
            base_covariance,
            s,
            block
          ),
          update_covariance = .iwmde_random_affine_covariance_block(
            update_covariance,
            s,
            block
          ),
          reference_coefficient = reference_coefficient,
          coefficients          = coefficients
        )
      }
      if (is.null(plan)) {
        return(NULL)
      }
      value <- .known_v_affine_log_likelihood(
        plan         = plan,
        residual     = outcome[block] - means[s, block],
        coefficients = coefficients
      )
      if (is.null(value)) {
        return(NULL)
      }
      out[, s] <- out[, s] + value
    }
  }
  out
}


.iwmde_random_affine_spectral_plan <- function(
    base_covariance, update_covariance, reference_coefficient,
    coefficients) {

  plan <- .known_v_affine_spectral_plan(
    base_covariance       = base_covariance,
    update_covariance     = update_covariance,
    reference_coefficient = reference_coefficient
  )
  if (!is.null(plan)) {
    return(plan)
  }

  candidates <- unique(coefficients[
    is.finite(coefficients) & coefficients != reference_coefficient
  ])
  if (length(candidates) == 0L) {
    return(NULL)
  }
  candidates <- candidates[order(abs(
    candidates - stats::median(candidates)
  ))]
  for (candidate in candidates) {
    candidate_covariance <- base_covariance +
      (candidate - reference_coefficient) * update_covariance
    plan <- .known_v_affine_spectral_plan(
      base_covariance       = candidate_covariance,
      update_covariance     = update_covariance,
      reference_coefficient = candidate
    )
    if (!is.null(plan)) {
      return(plan)
    }
  }

  return(NULL)
}


.iwmde_random_affine_covariance_invariant <- function(base_covariance,
                                                       update_covariance) {

  S <- dim(base_covariance)[1L]
  if (S < 2L) {
    return(TRUE)
  }
  all(vapply(2:S, function(s) {
    identical(
      base_covariance[s, , , drop = FALSE],
      base_covariance[1L, , , drop = FALSE]
    ) && identical(
      update_covariance[s, , , drop = FALSE],
      update_covariance[1L, , , drop = FALSE]
    )
  }, logical(1)))
}


.iwmde_random_affine_covariance_block <- function(x, draw, block) {

  if (is.matrix(x)) {
    return(x[block, block, drop = FALSE])
  }
  matrix(
    x[draw, block, block],
    nrow = length(block),
    ncol = length(block)
  )
}
