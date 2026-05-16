# ============================================================================ #
# IWMDE Replacement State
# ============================================================================ #

.iwmde_log_q_grid_scalar <- function(context, parameter, values, row_states,
                                     replacement) {

  out <- matrix(
    -Inf,
    nrow = length(values),
    ncol = length(row_states)
  )

  for (g in seq_along(values)) {
    for (i in seq_along(row_states)) {
      out[g, i] <- .iwmde_log_q_replacement(
        context     = context,
        state       = row_states[[i]],
        parameter   = parameter,
        value       = values[g],
        replacement = replacement
      )
    }
  }

  return(out)
}


.iwmde_log_q_grid_from_samples <- function(context, parameter, values,
                                           row_states, replacement,
                                           likelihood_mode, log_lik_fun) {

  if (length(values) == 0L || length(row_states) == 0L) {
    return(matrix(-Inf, nrow = length(values), ncol = length(row_states)))
  }

  likelihood_modes <- vapply(row_states, function(state) {
    state[["likelihood_mode"]]
  }, character(1))
  if (!all(likelihood_modes == likelihood_mode)) {
    return(NULL)
  }

  out  <- matrix(-Inf, nrow = length(values), ncol = length(row_states))
  keys <- vapply(row_states, function(state) {
    .iwmde_state_active_key(context, state)
  }, character(1))

  for (key in unique(keys)) {
    state_cols   <- which(keys == key)
    group_states <- row_states[state_cols]
    candidates   <- .iwmde_build_replacement_samples(
      context     = context,
      parameter   = parameter,
      values      = values,
      row_states  = group_states,
      replacement = replacement
    )

    valid <- candidates[["valid"]]
    if (!any(valid)) {
      next
    }

    valid_samples <- candidates[["samples"]][valid, , drop = FALSE]
    active_setup  <- group_states[[1L]][["active_setup"]]
    log_lik <- log_lik_fun(
      samples      = valid_samples,
      active_setup = active_setup
    )
    if (is.null(log_lik)) {
      return(NULL)
    }

    valid_positions <- which(valid)
    log_prior <- vapply(seq_along(valid_positions), function(i) {
      position <- valid_positions[i]
      state    <- group_states[[candidates[["state_index"]][position]]]

      .iwmde_log_prior_row(valid_samples[i, ], state[["prior_list"]])
    }, numeric(1))

    log_q <- log_lik + log_prior
    for (i in seq_along(valid_positions)) {
      position <- valid_positions[i]
      out[
        candidates[["grid_index"]][position],
        state_cols[candidates[["state_index"]][position]]
      ] <- log_q[i]
    }
  }

  return(out)
}


.iwmde_log_q_grid_batch <- function(context, parameter, values, row_states,
                                    replacement) {

  unit <- if (.is_data_multilevel(context[["data"]])) {
    "cluster"
  } else {
    "estimate"
  }

  return(.iwmde_log_q_grid_from_samples(
    context         = context,
    parameter       = parameter,
    values          = values,
    row_states      = row_states,
    replacement     = replacement,
    likelihood_mode = "marginal",
    log_lik_fun     = function(samples, active_setup) {
      likelihood_samples <- .iwmde_likelihood_posterior_samples(
        context      = context,
        samples      = samples,
        active_setup = active_setup
      )

      .iwmde_log_lik_from_posterior_samples_sum_active_branch(
        context           = context,
        posterior_samples = likelihood_samples,
        active_setup      = active_setup,
        unit              = unit
      )
    }
  ))
}


.iwmde_likelihood_posterior_samples <- function(context, samples,
                                                active_setup) {

  if (!isTRUE(active_setup[["is_weightfunction"]])) {
    return(samples)
  }

  global_spec <- context[["selection_spec"]]
  active_spec <- active_setup[["selection_spec"]]
  if (is.null(global_spec) || is.null(active_spec) ||
      .iwmde_same_p_cuts(global_spec[["p_cuts"]], active_spec[["p_cuts"]])) {
    return(samples)
  }

  omega <- .iwmde_indexed_parameter_matrix(samples, "omega")
  if (is.null(omega)) {
    log_omega <- .iwmde_indexed_parameter_matrix(samples, "log_omega")
    if (!is.null(log_omega)) {
      omega <- exp(log_omega)
    }
  }
  if (is.null(omega) || ncol(omega) != global_spec[["n_bins"]]) {
    return(samples)
  }

  active_omega <- .iwmde_collapse_omega_matrix(
    omega       = omega,
    global_cuts = global_spec[["p_cuts"]],
    active_cuts = active_spec[["p_cuts"]]
  )
  colnames(active_omega) <- paste0("omega[", seq_len(ncol(active_omega)), "]")

  keep <- !grepl("^(omega|log_omega)\\[[0-9]+\\]$", colnames(samples))
  out  <- cbind(samples[, keep, drop = FALSE], active_omega)

  return(out)
}


.iwmde_indexed_parameter_matrix <- function(samples, parameter) {

  cols <- grep(paste0("^", parameter, "\\[[0-9]+\\]$"), colnames(samples), value = TRUE)
  if (length(cols) == 0L) {
    return(NULL)
  }

  index <- as.integer(sub(
    paste0("^", parameter, "\\[([0-9]+)\\]$"),
    "\\1",
    cols
  ))
  cols <- cols[order(index)]

  return(samples[, cols, drop = FALSE])
}


.iwmde_build_replacement_samples <- function(context, parameter, values,
                                             row_states, replacement) {

  sample_names <- colnames(context[["posterior_samples"]])
  n_candidates <- length(values) * length(row_states)

  if (n_candidates == 0L) {
    samples <- matrix(NA_real_, nrow = n_candidates, ncol = length(sample_names))
    colnames(samples) <- sample_names
    return(list(
      samples     = samples,
      valid       = logical(n_candidates),
      grid_index  = integer(n_candidates),
      state_index = integer(n_candidates)
    ))
  }

  grid_index  <- rep(seq_along(values), times = length(row_states))
  state_index <- rep(seq_along(row_states), each = length(values))
  row_index   <- vapply(row_states, function(state) {
    state[["row_index"]]
  }, integer(1))
  samples <- context[["posterior_samples"]][
    row_index[state_index],
    sample_names,
    drop = FALSE
  ]
  valid <- rep(TRUE, n_candidates)

  if (identical(replacement[["type"]], "linear")) {
    valid <- rep(FALSE, n_candidates)

    for (i in seq_along(row_states)) {
      positions <- which(state_index == i)
      linear    <- .iwmde_linear_replacement_state(
        context     = context,
        state       = row_states[[i]],
        replacement = replacement
      )

      if (!isTRUE(linear[["valid"]])) {
        next
      }
      if (length(linear[["active_columns"]]) == 0L) {
        valid[positions] <- vapply(values, function(value) {
          isTRUE(all.equal(linear[["current"]], value, tolerance = 1e-10))
        }, logical(1))
        next
      }

      value_delta <- values - linear[["current"]]
      finite      <- is.finite(value_delta)
      valid[positions] <- finite
      if (any(finite)) {
        samples[positions[finite], linear[["active_columns"]]] <-
          samples[positions[finite], linear[["active_columns"]], drop = FALSE] +
          outer(
            value_delta[finite],
            linear[["coefficients"]]
          )
      }
    }
  } else if (parameter %in% sample_names) {
    samples[, parameter] <- values[grid_index]
  }

  return(list(
    samples     = samples,
    valid       = valid,
    grid_index  = grid_index,
    state_index = state_index
  ))
}


.iwmde_replace_row_for_value <- function(context, state, parameter, value,
                                         replacement) {

  out <- function(row, valid = TRUE) {
    return(list(row = row, valid = valid))
  }

  if (identical(replacement[["type"]], "linear")) {
    return(.iwmde_replace_linear_row(
      context     = context,
      state       = state,
      value       = value,
      replacement = replacement
    ))
  }

  row <- state[["row"]]
  row[[parameter]] <- value

  return(out(row))
}


.iwmde_replace_linear_row <- function(context, state, value, replacement) {

  out <- function(row, valid = TRUE) {
    return(list(row = row, valid = valid))
  }

  row    <- state[["row"]]
  linear <- .iwmde_linear_replacement_state(
    context     = context,
    state       = state,
    replacement = replacement
  )
  if (!isTRUE(linear[["valid"]]) || !is.finite(value)) {
    return(out(row, valid = FALSE))
  }

  if (length(linear[["active_columns"]]) == 0L) {
    valid <- isTRUE(all.equal(linear[["current"]], value, tolerance = 1e-10))
    return(out(row, valid = valid))
  }

  row[linear[["active_columns"]]] <- row[linear[["active_columns"]]] +
    (value - linear[["current"]]) * linear[["coefficients"]]

  return(out(row))
}


.iwmde_linear_replacement_state <- function(context, state, replacement) {

  weights <- replacement[["weights"]]
  key <- paste(
    "linear_replacement",
    state[["row_index"]],
    paste(names(weights), .iwmde_key_number(weights), sep = "=", collapse = ","),
    sep = "|"
  )
  if (exists(key, envir = context[["row_cache"]], inherits = FALSE)) {
    return(get(key, envir = context[["row_cache"]], inherits = FALSE))
  }

  row     <- state[["row"]]
  current <- .iwmde_linear_value_row(row, weights)
  if (!is.finite(current)) {
    out <- list(valid = FALSE, current = current)
    assign(key, out, envir = context[["row_cache"]])
    return(out)
  }

  active_columns <- .iwmde_linear_active_columns(context, row, weights)
  if (length(active_columns) == 0L) {
    out <- list(
      valid          = TRUE,
      current        = current,
      active_columns = character(),
      active_weights = numeric(),
      denominator    = 0,
      coefficients   = numeric()
    )
    assign(key, out, envir = context[["row_cache"]])
    return(out)
  }

  active_weights <- weights[active_columns]
  denominator    <- sum(active_weights^2)
  if (!is.finite(denominator) || denominator <= 0) {
    out <- list(valid = FALSE, current = current)
    assign(key, out, envir = context[["row_cache"]])
    return(out)
  }

  out <- list(
    valid          = TRUE,
    current        = current,
    active_columns = active_columns,
    active_weights = active_weights,
    denominator    = denominator,
    coefficients   = active_weights / denominator
  )
  assign(key, out, envir = context[["row_cache"]])

  return(out)
}


.iwmde_log_q_replacement <- function(context, state, parameter, value,
                                     replacement) {

  row <- state[["row"]]
  if (!identical(replacement[["type"]], "linear")) {
    row[[parameter]] <- value
  }
  replaced <- .iwmde_replace_parameters(
    context     = context,
    state       = state,
    parameter   = parameter,
    value       = value,
    row         = row,
    replacement = replacement
  )
  if (!isTRUE(replaced[["valid"]])) {
    return(-Inf)
  }
  row        <- replaced[["row"]]
  parameters <- replaced[["parameters"]]

  if (isTRUE(state[["use_focal_prior_delta"]]) &&
      isTRUE(replaced[["use_focal_prior_delta"]])) {
    log_lik <- .iwmde_log_likelihood_parameters(
      context         = context,
      parameters      = parameters,
      active_setup    = state[["active_setup"]],
      likelihood_mode = state[["likelihood_mode"]],
      row             = row
    )
    if (!is.finite(log_lik)) {
      return(-Inf)
    }

    focal_log_prior <- .iwmde_focal_log_prior(
      prior     = state[["focal_prior"]],
      value     = value,
      parameter = parameter
    )
    if (!is.finite(focal_log_prior)) {
      return(-Inf)
    }

    log_prior <- state[["baseline_log_prior"]] +
      focal_log_prior - state[["baseline_focal_log_prior"]]
    out <- log_lik + log_prior
    if (!is.finite(out)) {
      return(-Inf)
    }

    return(out)
  }

  return(.iwmde_log_q_state(
    context         = context,
    row             = row,
    active_setup    = state[["active_setup"]],
    parameters      = parameters,
    prior_list      = state[["prior_list"]],
    likelihood_mode = state[["likelihood_mode"]]
  ))
}


.iwmde_replace_parameters <- function(context, state, parameter, value, row,
                                      replacement) {

  out <- function(row, parameters, use_focal_prior_delta = TRUE,
                  valid = TRUE) {

    return(list(
      row                   = row,
      parameters            = parameters,
      use_focal_prior_delta = use_focal_prior_delta,
      valid                 = valid
    ))
  }

  if (identical(replacement[["type"]], "linear")) {
    return(.iwmde_replace_linear_parameters(
      context     = context,
      state       = state,
      value       = value,
      replacement = replacement
    ))
  }

  if (identical(replacement[["type"]], "scalar")) {
    parameters <- state[["parameters"]]
    name       <- replacement[["name"]]
    if (!is.null(parameters[[name]])) {
      parameters[[name]] <- value
      return(out(row, parameters))
    }
  }

  if (identical(replacement[["type"]], "indexed")) {
    parameters <- state[["parameters"]]
    name       <- replacement[["name"]]
    index      <- replacement[["index"]]
    if (!is.null(parameters[[name]]) && length(parameters[[name]]) >= index) {
      parameters[[name]][index] <- value
      return(out(row, parameters))
    }
  }

  return(out(
    row        = row,
    parameters = .iwmde_row_parameters(context, row, state[["active_setup"]])
  ))
}


.iwmde_replace_linear_parameters <- function(context, state, value,
                                             replacement) {

  out <- function(row, parameters, valid = TRUE) {

    return(list(
      row                   = row,
      parameters            = parameters,
      use_focal_prior_delta = FALSE,
      valid                 = valid
    ))
  }

  row    <- state[["row"]]
  linear <- .iwmde_linear_replacement_state(
    context     = context,
    state       = state,
    replacement = replacement
  )
  if (!isTRUE(linear[["valid"]]) || !is.finite(value)) {
    return(out(row, state[["parameters"]], valid = FALSE))
  }

  if (length(linear[["active_columns"]]) == 0L) {
    valid <- isTRUE(all.equal(linear[["current"]], value, tolerance = 1e-10))
    return(out(row, state[["parameters"]], valid = valid))
  }

  # Move along the minimum-norm direction that changes a' beta to the target
  # value while keeping the orthogonal complement fixed.
  row[linear[["active_columns"]]] <- row[linear[["active_columns"]]] +
    (value - linear[["current"]]) * linear[["coefficients"]]
  parameters <- .iwmde_row_parameters(context, row, state[["active_setup"]])

  return(out(row, parameters))
}


.iwmde_linear_value_row <- function(row, weights) {

  return(sum(as.numeric(row[names(weights)]) * weights))
}


.iwmde_linear_active_columns <- function(context, row, weights) {

  active <- vapply(names(weights), function(parameter) {
    state <- .iwmde_focal_prior_state(context, parameter, row)

    return(identical(state[["status"]], "continuous"))
  }, logical(1))

  return(names(weights)[active])
}


.iwmde_replacement_spec <- function(context, parameter,
                                    parameter_spec = NULL) {

  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "linear")) {
    return(list(type = "linear", weights = parameter_spec[["weights"]]))
  }

  data <- context[["data"]]
  if (.iwmde_parameter_is_eta(parameter)) {
    return(list(type = "fallback"))
  }
  if (.is_data_mods(data) && startsWith(parameter, "mu_")) {
    return(list(type = "fallback"))
  }
  if (.is_data_scale(data) && startsWith(parameter, "log_tau_")) {
    return(list(type = "fallback"))
  }

  if (parameter %in% c("mu", "tau", "log_tau", "rho", "PET", "PEESE")) {
    return(list(type = "scalar", name = parameter))
  }

  indexed <- regexec("^(.+)\\[([0-9]+)\\]$", parameter)
  match   <- regmatches(parameter, indexed)[[1]]
  if (length(match) == 3L && match[2] %in% c("gamma", "pi", "theta", "phi")) {
    return(list(
      type  = "indexed",
      name  = match[2],
      index = as.integer(match[3])
    ))
  }

  return(list(type = "fallback"))
}


