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

    valid_positions <- which(valid)
    valid_samples   <- candidates[["samples"]][valid, , drop = FALSE]
    active_setup    <- group_states[[1L]][["active_setup"]]
    log_lik <- log_lik_fun(
      samples      = valid_samples,
      active_setup = active_setup,
      batch        = list(
        candidates      = candidates,
        valid_positions = valid_positions,
        row_states      = group_states
      )
    )
    if (is.null(log_lik)) {
      return(NULL)
    }

    if (!is.numeric(log_lik) || length(log_lik) != length(valid_positions)) {
      warning(
        "Batched IWMDE likelihood returned an invalid length; falling back to scalar evaluation.",
        call. = FALSE
      )
      return(NULL)
    }

    log_prior <- .iwmde_replacement_log_prior(
      parameter       = parameter,
      values          = values,
      valid_samples   = valid_samples,
      valid_positions = valid_positions,
      candidates      = candidates,
      row_states      = group_states,
      replacement     = replacement
    )

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


.iwmde_replacement_log_prior <- function(parameter, values, valid_samples,
                                         valid_positions, candidates,
                                         row_states, replacement) {

  state_index <- candidates[["state_index"]][valid_positions]
  grid_index  <- candidates[["grid_index"]][valid_positions]
  log_prior   <- rep(NA_real_, length(valid_positions))

  use_delta <- vapply(row_states, function(state) {
    identical(replacement[["type"]], "scalar") &&
      isTRUE(state[["use_focal_prior_delta"]])
  }, logical(1))
  prior_list      <- row_states[[1L]][["prior_list"]]
  same_prior_list <- all(vapply(row_states, function(state) {
    identical(state[["prior_list"]], prior_list)
  }, logical(1)))
  if (all(use_delta[state_index])) {
    focal_prior      <- row_states[[1L]][["focal_prior"]]
    same_focal_prior <- all(vapply(row_states, function(state) {
      identical(state[["focal_prior"]], focal_prior)
    }, logical(1)))
    if (same_focal_prior) {
      baseline_log_prior <- vapply(
        row_states,
        `[[`,
        numeric(1),
        "baseline_log_prior"
      )
      baseline_focal_log_prior <- vapply(
        row_states,
        `[[`,
        numeric(1),
        "baseline_focal_log_prior"
      )
      focal_log_prior <- .iwmde_focal_log_prior_values(
        prior     = focal_prior,
        values    = values,
        parameter = parameter
      )
      return(
        baseline_log_prior[state_index] +
          focal_log_prior[grid_index] -
          baseline_focal_log_prior[state_index]
      )
    }
  }
  if (!any(use_delta[state_index]) && same_prior_list) {
    return(.iwmde_replacement_log_prior_rows(
      samples     = valid_samples,
      prior_list  = prior_list,
      replacement = replacement,
      evaluator   = row_states[[1L]][["prior_evaluator"]]
    ))
  }

  for (i in unique(state_index)) {
    positions <- which(state_index == i)
    state     <- row_states[[i]]

    if (use_delta[i]) {
      focal_log_prior <- .iwmde_focal_log_prior_values(
        prior     = state[["focal_prior"]],
        values    = values[grid_index[positions]],
        parameter = parameter
      )
      log_prior[positions] <- state[["baseline_log_prior"]] +
        focal_log_prior - state[["baseline_focal_log_prior"]]
      next
    }

    # Vector and Dirichlet priors must be evaluated one posterior row at a
    # time; their auxiliary coordinates are joint within a row.
    log_prior[positions] <- .iwmde_replacement_log_prior_rows(
      samples     = valid_samples[positions, , drop = FALSE],
      prior_list  = state[["prior_list"]],
      replacement = replacement,
      evaluator   = state[["prior_evaluator"]]
    )
  }

  return(log_prior)
}


.iwmde_replacement_log_prior_rows <- function(samples, prior_list,
                                               replacement,
                                               evaluator = NULL) {

  if (!identical(replacement[["type"]], "simplex_pair")) {
    return(.iwmde_log_prior_rows(samples, prior_list, evaluator = evaluator))
  }

  parameter <- replacement[["parameter"]]
  if (is.null(parameter)) {
    candidates <- names(prior_list)[vapply(prior_list, function(prior) {
      inherits(prior, "prior.simplex") &&
        identical(prior[["distribution"]], "dirichlet")
    }, logical(1))]
    if (length(candidates) == 1L) {
      parameter <- candidates[[1L]]
    }
  }
  prior     <- if (is.null(parameter)) NULL else prior_list[[parameter]]
  eta_names <- replacement[["auxiliary_columns"]]
  if (is.null(eta_names) && !is.null(prior)) {
    eta_names <- .iwmde_simplex_auxiliary_columns(
      parameter,
      length(prior[["parameters"]][["alpha"]])
    )
  }
  if (is.null(prior) || !inherits(prior, "prior.simplex") ||
      !identical(prior[["distribution"]], "dirichlet") ||
      length(prior[["parameters"]][["alpha"]]) != length(eta_names) ||
      !all(eta_names %in% colnames(samples))) {
    stop(
      "The simplex density target does not have its matching Dirichlet prior ",
      "and auxiliary coordinates.",
      call. = FALSE
    )
  }

  log_prior <- .iwmde_log_prior_rows(
    samples    = samples,
    prior_list = prior_list[setdiff(names(prior_list), parameter)]
  )
  alpha <- prior[["parameters"]][["alpha"]]
  for (i in seq_along(alpha)) {
    eta     <- samples[, eta_names[[i]]]
    invalid <- !is.finite(eta) | eta < 0
    term    <- stats::dgamma(eta, shape = alpha[[i]], rate = 1, log = TRUE)
    term[invalid] <- -Inf
    log_prior <- log_prior + term
  }

  return(log_prior)
}


.iwmde_log_prior_rows <- function(samples, prior_list, evaluator = NULL) {

  if (length(prior_list) == 0L) {
    return(numeric(nrow(samples)))
  }

  samples   <- .resolve_fixed_prior_sample_columns(samples, prior_list)
  if (is.null(evaluator)) {
    evaluator <- BayesTools::JAGS_marglik_priors_rows_evaluator(prior_list)
  }
  log_prior <- evaluator(samples)
  if (!is.numeric(log_prior) || length(log_prior) != nrow(samples) ||
      !is.null(dim(log_prior))) {
    stop("Batched IWMDE prior evaluation returned an invalid result.", call. = FALSE)
  }

  return(as.numeric(log_prior))
}


.iwmde_prior_rows_evaluator <- function(context, row, parameter, state_scope,
                                        prior_list) {

  context <- .iwmde_context_ensure_caches(context)
  key <- paste(
    .iwmde_active_key(context, row),
    state_scope,
    if (.iwmde_parameter_is_eta(parameter)) "eta" else "ordinary",
    sep = "|"
  )
  cache <- context[["prior_cache"]]
  if (exists(key, envir = cache, inherits = FALSE)) {
    cached <- get(key, envir = cache, inherits = FALSE)
    if (!identical(cached[["prior_list"]], prior_list)) {
      stop("Cached IWMDE prior structure does not match its active branch.", call. = FALSE)
    }
    return(cached[["evaluator"]])
  }

  evaluator <- BayesTools::JAGS_marglik_priors_rows_evaluator(prior_list)
  assign(
    key,
    list(prior_list = prior_list, evaluator = evaluator),
    envir = cache
  )
  return(evaluator)
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
    log_lik_fun     = function(samples, active_setup, batch) {
      fixed_mu_samples <- .iwmde_q_grid_known_v_random_fixed_mu(
        context      = context,
        active_setup = active_setup,
        unit         = unit,
        batch        = batch
      )
      .iwmde_log_lik_from_posterior_samples_sum_active_branch(
        context           = context,
        posterior_samples = samples,
        active_setup      = active_setup,
        unit              = unit,
        fixed_mu_samples  = fixed_mu_samples
      )
    }
  ))
}


.iwmde_q_grid_known_v_random_fixed_mu <- function(
    context, active_setup, unit, batch) {

  if (!.iwmde_uses_known_v_random_marginal_likelihood(context) ||
      !.iwmde_q_grid_fixed_mu_is_invariant(context, batch[["candidates"]])) {
    return(NULL)
  }

  valid_positions <- batch[["valid_positions"]]
  state_index <- batch[["candidates"]][["state_index"]][valid_positions]
  row_states  <- batch[["row_states"]]
  if (length(state_index) == 0L ||
      anyNA(state_index) || any(state_index < 1L) ||
      any(state_index > length(row_states))) {
    return(NULL)
  }

  setup <- .iwmde_predictor_setup(
    context      = context,
    row_states   = row_states,
    active_setup = active_setup,
    unit         = unit
  )
  fixed_mu <- setup[["mu"]]
  if (!is.matrix(fixed_mu) || nrow(fixed_mu) != length(row_states) ||
      any(!is.finite(fixed_mu))) {
    return(NULL)
  }

  return(fixed_mu[state_index, , drop = FALSE])
}


.iwmde_q_grid_fixed_mu_is_invariant <- function(context, candidates) {

  changed <- candidates[["changed_parameters"]]
  if (!is.character(changed) || length(changed) == 0L || anyNA(changed)) {
    return(FALSE)
  }

  if (any(changed %in% c("mu", "PET", "PEESE"))) {
    return(FALSE)
  }

  fit <- context[["formula_fit"]]
  if (is.null(fit) && !is.null(context[["object"]])) {
    fit <- context[["object"]][["fit"]]
  }
  if (is.null(fit) || !inherits(fit, "BayesTools_fit")) {
    return(FALSE)
  }

  dependencies <- tryCatch(
    BayesTools::JAGS_formula_coordinate_dependencies(
      fit,
      coordinates = unique(changed)
    ),
    error = function(e) NULL
  )
  if (is.null(dependencies)) {
    return(FALSE)
  }

  return(!any(dependencies[["formula_parameter"]] == "mu"))
}


.iwmde_likelihood_posterior_samples <- function(context, samples,
                                                active_setup) {

  samples <- .iwmde_localize_active_branch_samples(
    context      = context,
    samples      = samples,
    active_setup = active_setup
  )
  samples <- .iwmde_resolve_fixed_prior_sample_columns(
    context      = context,
    samples      = samples,
    active_setup = active_setup
  )

  if (!isTRUE(active_setup[["is_weightfunction"]])) {
    return(samples)
  }

  global_spec <- context[["selection_spec"]]
  active_spec <- active_setup[["selection_spec"]]
  if (is.null(global_spec) || is.null(active_spec)) {
    return(samples)
  }

  active_omega <- .iwmde_active_branch_omega_matrix(
    samples     = samples,
    global_spec = global_spec,
    active_spec = active_spec
  )
  if (is.null(active_omega)) {
    return(samples)
  }

  keep <- !.iwmde_indexed_any_parameter_columns(
    columns    = colnames(samples),
    parameters = unique(c(
      global_spec[["jags_omega"]],
      active_spec[["jags_omega"]]
    ))
  )
  out  <- cbind(samples[, keep, drop = FALSE], active_omega)

  return(out)
}


.iwmde_resolve_fixed_prior_sample_columns <- function(context, samples,
                                                      active_setup) {

  prior_list <- active_setup[["fit_priors"]]
  if (is.null(prior_list) &&
      !is.null(context[["data"]]) &&
      !is.null(active_setup[["priors"]]) &&
      .iwmde_can_create_fit_priors(context[["data"]])) {
    prior_list <- .create_fit_priors(
      data   = context[["data"]],
      priors = active_setup[["priors"]]
    )
  }

  return(.resolve_fixed_prior_sample_columns(samples, prior_list))
}


.iwmde_can_create_fit_priors <- function(data) {

  outcome_type <- .data_outcome_type(data)

  return(
    length(outcome_type) == 1L &&
      !is.na(outcome_type) &&
      outcome_type %in% c("norm", "bin", "pois")
  )
}


.iwmde_active_branch_omega_matrix <- function(samples, global_spec,
                                              active_spec) {

  if (is.null(global_spec[["jags_omega"]]) ||
      is.null(active_spec[["jags_omega"]])) {
    return(NULL)
  }

  active_omega <- BayesTools::JAGS_indexed_parameter_matrix(
    samples   = samples,
    parameter = active_spec[["jags_omega"]]
  )
  if (!is.null(active_omega) &&
      ncol(active_omega) == active_spec[["n_bins"]]) {
    colnames(active_omega) <- .iwmde_indexed_parameter_names(
      parameter = active_spec[["jags_omega"]],
      n         = active_spec[["n_bins"]]
    )
    return(active_omega)
  }

  omega <- BayesTools::JAGS_indexed_parameter_matrix(
    samples   = samples,
    parameter = global_spec[["jags_omega"]]
  )
  if (is.null(omega) || ncol(omega) != global_spec[["n_bins"]]) {
    return(NULL)
  }

  if (.iwmde_same_p_cuts(global_spec[["p_cuts"]], active_spec[["p_cuts"]])) {
    active_omega <- omega
  } else {
    active_omega <- .iwmde_collapse_omega_matrix(
      omega       = omega,
      global_cuts = global_spec[["p_cuts"]],
      active_cuts = active_spec[["p_cuts"]]
    )
  }
  colnames(active_omega) <- .iwmde_indexed_parameter_names(
    parameter = active_spec[["jags_omega"]],
    n         = ncol(active_omega)
  )

  return(active_omega)
}


.iwmde_localize_active_branch_samples <- function(context, samples,
                                                  active_setup) {

  # Active priors have already been selected from the product-space mixture.
  # Likelihood dispatch should therefore see local branch indicators.
  .iwmde_check_active_branch_priors(active_setup)

  columns <- colnames(samples)
  if (is.null(columns)) {
    return(samples)
  }

  indicator_names <- context[["indicator_names"]]
  if (is.null(indicator_names)) {
    indicator_names <- grep("(^|_)indicator$", columns, value = TRUE)
    indicator_names <- c(
      indicator_names,
      intersect("bias_indicator", columns)
    )
  }
  indicator_names <- unique(indicator_names)
  indicator_names <- intersect(indicator_names, columns)
  if (length(indicator_names) == 0L) {
    return(samples)
  }

  samples[, indicator_names] <- 1
  return(samples)
}


.iwmde_check_active_branch_priors <- function(active_setup) {

  if (.iwmde_prior_tree_has_mixture(active_setup[["priors"]])) {
    stop(
      "IWMDE active-branch likelihood requires priors localized to a single mixture component.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}


.iwmde_prior_tree_has_mixture <- function(x) {

  if (is.null(x)) {
    return(FALSE)
  }
  if (BayesTools::is.prior.mixture(x)) {
    return(TRUE)
  }
  if (BayesTools::is.prior(x)) {
    return(FALSE)
  }
  if (!is.list(x)) {
    return(FALSE)
  }

  return(any(vapply(x, .iwmde_prior_tree_has_mixture, logical(1))))
}


.iwmde_indexed_any_parameter_columns <- function(columns, parameters) {

  parameters <- parameters[!is.na(parameters) & nzchar(parameters)]
  if (length(parameters) == 0L) {
    return(rep(FALSE, length(columns)))
  }

  out <- rep(FALSE, length(columns))
  for (parameter in parameters) {
    out <- out | BayesTools::JAGS_indexed_parameter_columns(
      columns   = columns,
      parameter = parameter
    )
  }

  return(out)
}


.iwmde_indexed_parameter_names <- function(parameter, n) {

  return(paste0(parameter, "[", seq_len(n), "]"))
}


.iwmde_sync_invgamma_auxiliary_matrix <- function(context, samples,
                                                  parameters) {

  if (length(parameters) == 0L || is.null(colnames(samples))) {
    return(list(samples = samples, valid = rep(TRUE, nrow(samples))))
  }

  valid <- rep(TRUE, nrow(samples))
  for (parameter in unique(parameters)) {
    auxiliary <- .iwmde_invgamma_auxiliary_name(context, parameter)
    if (is.null(auxiliary) ||
        !parameter %in% colnames(samples) ||
        !auxiliary %in% colnames(samples)) {
      next
    }

    sync_rows <- .iwmde_invgamma_auxiliary_rows(
      context   = context,
      samples   = samples,
      parameter = parameter
    )
    if (length(sync_rows) == 0L) {
      next
    }

    values          <- samples[sync_rows, parameter]
    finite_positive <- is.finite(values) & values > 0
    samples[sync_rows[finite_positive], auxiliary] <-
      1 / values[finite_positive]
    samples[sync_rows[!finite_positive], auxiliary] <- NA_real_
    valid[sync_rows] <- valid[sync_rows] & finite_positive
  }

  return(list(samples = samples, valid = valid))
}


.iwmde_sync_invgamma_auxiliary_row <- function(context, row, parameters) {

  if (length(parameters) == 0L) {
    return(list(row = row, valid = TRUE))
  }

  valid <- TRUE
  for (parameter in unique(parameters)) {
    auxiliary <- .iwmde_invgamma_auxiliary_name(context, parameter)
    if (is.null(auxiliary) ||
        !parameter %in% names(row) ||
        !auxiliary %in% names(row) ||
        !.iwmde_row_uses_invgamma_prior(context, parameter, row)) {
      next
    }

    value <- as.numeric(row[[parameter]])
    if (length(value) == 1L && is.finite(value) && value > 0) {
      row[[auxiliary]] <- 1 / value
    } else {
      row[[auxiliary]] <- NA_real_
      valid <- FALSE
    }
  }

  return(list(row = row, valid = valid))
}


.iwmde_invgamma_auxiliary_name <- function(context, parameter) {

  if (grepl("\\[[0-9]+\\]$", parameter)) {
    return(NULL)
  }
  prior_name <- .iwmde_parameter_prior_name(context, parameter)
  if (is.null(prior_name) ||
      !.iwmde_prior_contains_invgamma(
        context[["flat_prior_list"]][[prior_name]]
      )) {
    return(NULL)
  }

  return(paste0("inv_", parameter))
}


.iwmde_invgamma_auxiliary_rows <- function(context, samples, parameter) {

  rows <- seq_len(nrow(samples))
  if (length(rows) == 0L) {
    return(rows)
  }

  active_keys <- .iwmde_active_keys_matrix(context, samples)
  unique_keys <- unique(active_keys)
  first_rows  <- match(unique_keys, active_keys)
  uses_by_key <- vapply(first_rows, function(i) {
    .iwmde_row_uses_invgamma_prior(context, parameter, samples[i, ])
  }, logical(1))

  return(rows[uses_by_key[match(active_keys, unique_keys)]])
}


.iwmde_active_keys_matrix <- function(context, samples) {

  indicator_names <- context[["indicator_names"]]
  if (is.null(indicator_names)) {
    indicator_names <- grep("(^|_)indicator$", colnames(samples), value = TRUE)
    indicator_names <- c(
      indicator_names,
      intersect("bias_indicator", colnames(samples))
    )
    indicator_names <- unique(indicator_names)
  }

  if (length(indicator_names) == 0L) {
    return(rep("all", nrow(samples)))
  }

  parts <- lapply(indicator_names, function(name) {
    if (!name %in% colnames(samples)) {
      return(rep(paste(name, "NA", sep = "="), nrow(samples)))
    }

    values <- vapply(samples[, name], function(value) {
      as.character(.iwmde_indicator_index(value, name))
    }, character(1))

    paste(name, values, sep = "=")
  })

  return(do.call(paste, c(parts, sep = "|")))
}


.iwmde_row_uses_invgamma_prior <- function(context, parameter, row) {

  prior <- .iwmde_focal_prior(context, parameter, row)

  return(.iwmde_prior_is_invgamma(prior))
}


.iwmde_prior_contains_invgamma <- function(prior) {

  if (is.null(prior) || BayesTools::is.prior.none(prior)) {
    return(FALSE)
  }
  if (BayesTools::is.prior.mixture(prior)) {
    return(any(vapply(seq_along(prior), function(i) {
      .iwmde_prior_contains_invgamma(prior[[i]])
    }, logical(1))))
  }

  return(.iwmde_prior_is_invgamma(prior))
}


.iwmde_prior_is_invgamma <- function(prior) {

  inherits(prior, "prior.simple") &&
    identical(prior[["distribution"]], "invgamma")
}


.iwmde_build_replacement_samples <- function(context, parameter, values,
                                             row_states, replacement) {

  state_scope <- unique(vapply(
    row_states,
    .iwmde_state_scope_value,
    character(1)
  ))
  source_samples <- context[["posterior_samples"]]
  if (length(state_scope) == 1L && identical(state_scope, "global")) {
    source_samples <- .iwmde_drop_local_latent_sample_columns(
      source_samples,
      context
    )
  }
  sample_names <- colnames(source_samples)
  n_candidates <- length(values) * length(row_states)

  if (n_candidates == 0L) {
    samples <- matrix(NA_real_, nrow = n_candidates, ncol = length(sample_names))
    colnames(samples) <- sample_names
    return(list(
      samples            = samples,
      valid              = logical(n_candidates),
      grid_index         = integer(n_candidates),
      state_index        = integer(n_candidates),
      changed_parameters = character()
    ))
  }

  grid_index  <- rep(seq_along(values), times = length(row_states))
  state_index <- rep(seq_along(row_states), each = length(values))
  row_index   <- vapply(row_states, function(state) {
    state[["row_index"]]
  }, integer(1))
  samples <- source_samples[
    row_index[state_index],
    sample_names,
    drop = FALSE
  ]
  valid <- rep(TRUE, n_candidates)
  changed_parameters <- character()

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
        valid[positions] <- values == linear[["current"]]
        next
      }

      value_delta <- values - linear[["current"]]
      finite      <- is.finite(value_delta)
      valid[positions] <- finite
      if (any(finite)) {
        changed_parameters <- unique(c(
          changed_parameters,
          linear[["active_columns"]]
        ))
        samples[positions[finite], linear[["active_columns"]]] <-
          samples[positions[finite], linear[["active_columns"]], drop = FALSE] +
          outer(
            value_delta[finite],
            linear[["coefficients"]]
          )
      }
    }
  } else if (identical(replacement[["type"]], "simplex_pair")) {
    n_targets         <- replacement[["n_targets"]]
    columns           <- paste0(
      replacement[["parameter"]],
      "[", seq_len(n_targets), "]"
    )
    auxiliary_columns <- replacement[["auxiliary_columns"]]
    index             <- replacement[["index"]]
    target            <- values[grid_index]
    eta_sum           <- rowSums(samples[, auxiliary_columns, drop = FALSE])
    other             <- setdiff(seq_len(n_targets), index)
    other_eta_sum     <- rowSums(
      samples[, auxiliary_columns[other], drop = FALSE]
    )
    valid             <- is.finite(target) & target >= 0 & target <= 1 &
      is.finite(eta_sum) & eta_sum > 0 &
      is.finite(other_eta_sum) & other_eta_sum > 0
    samples[valid, columns[[index]]] <- target[valid]
    samples[valid, auxiliary_columns[[index]]] <-
      eta_sum[valid] * target[valid]
    other_proportions <- samples[valid, auxiliary_columns[other], drop = FALSE] /
      other_eta_sum[valid]
    samples[valid, columns[other]] <-
      other_proportions * (1 - target[valid])
    samples[valid, auxiliary_columns[other]] <-
      other_proportions * eta_sum[valid] * (1 - target[valid])
    changed_parameters <- c(columns, auxiliary_columns)
  } else if (identical(replacement[["type"]], "random_component_sd")) {
    target     <- values[grid_index]
    multiplier <- .iwmde_random_component_sd_multiplier(
      samples,
      replacement[["factors"]]
    )
    source     <- replacement[["source_parameter"]]
    valid      <- is.finite(target) & target >= 0 &
      is.finite(multiplier) & multiplier > 0
    samples[valid, source] <- target[valid] / multiplier[valid]
    changed_parameters <- source
  } else if (parameter %in% sample_names) {
    samples[, parameter] <- values[grid_index]
    changed_parameters <- parameter
  }
  synced <- .iwmde_sync_invgamma_auxiliary_matrix(
    context    = context,
    samples    = samples,
    parameters = changed_parameters
  )
  samples <- synced[["samples"]]
  valid   <- valid & synced[["valid"]]
  synced <- .iwmde_sync_random_allocation_sd_matrix(
    context    = context,
    samples    = samples,
    parameters = changed_parameters
  )
  samples <- synced[["samples"]]
  valid   <- valid & synced[["valid"]]
  return(list(
    samples            = samples,
    valid              = valid,
    grid_index         = grid_index,
    state_index        = state_index,
    changed_parameters = unique(changed_parameters)
  ))
}


.iwmde_replace_row_for_value <- function(context, state, parameter, value,
                                         replacement) {

  out <- function(row, valid = TRUE, changed_parameters = character()) {

    if (isTRUE(valid) && length(changed_parameters) > 0L) {
      synced <- .iwmde_sync_random_allocation_sd_row(
        context    = context,
        row        = row,
        parameters = changed_parameters
      )
      row   <- synced[["row"]]
      valid <- synced[["valid"]]
    }

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

  if (identical(replacement[["type"]], "simplex_pair")) {
    n_targets         <- replacement[["n_targets"]]
    columns           <- paste0(
      replacement[["parameter"]],
      "[", seq_len(n_targets), "]"
    )
    auxiliary_columns <- replacement[["auxiliary_columns"]]
    index             <- replacement[["index"]]
    other             <- setdiff(seq_len(n_targets), index)
    if (!is.finite(value) || value < 0 || value > 1 ||
        !all(c(columns, auxiliary_columns) %in% names(state[["row"]]))) {
      return(out(state[["row"]], valid = FALSE))
    }
    row     <- state[["row"]]
    eta_sum <- sum(row[auxiliary_columns])
    other_eta_sum <- sum(row[auxiliary_columns[other]])
    if (!is.finite(eta_sum) || eta_sum <= 0 ||
        !is.finite(other_eta_sum) || other_eta_sum <= 0) {
      return(out(row, valid = FALSE))
    }
    other_proportions <- row[auxiliary_columns[other]] / other_eta_sum
    row[[columns[[index]]]] <- value
    row[[auxiliary_columns[[index]]]] <- eta_sum * value
    row[columns[other]] <- other_proportions * (1 - value)
    row[auxiliary_columns[other]] <-
      other_proportions * eta_sum * (1 - value)

    return(out(
      row,
      changed_parameters = c(columns, auxiliary_columns)
    ))
  }

  if (identical(replacement[["type"]], "random_component_sd")) {
    row <- state[["row"]]
    multiplier <- .iwmde_random_component_sd_multiplier(
      matrix(row, nrow = 1L, dimnames = list(NULL, names(row))),
      replacement[["factors"]]
    )
    source <- replacement[["source_parameter"]]
    if (!is.finite(value) || value < 0 || !is.finite(multiplier) ||
        multiplier <= 0 || !source %in% names(row)) {
      return(out(row, valid = FALSE))
    }
    row[[source]] <- value / multiplier
    synced <- .iwmde_sync_replacement_row(
      context    = context,
      row        = row,
      parameters = source
    )

    return(out(
      synced[["row"]],
      valid = synced[["valid"]]
    ))
  }

  row <- state[["row"]]
  row[[parameter]] <- value
  synced <- .iwmde_sync_replacement_row(
    context    = context,
    row        = row,
    parameters = parameter
  )

  return(out(
    synced[["row"]],
    valid = synced[["valid"]]
  ))
}


.iwmde_sync_random_allocation_sd_matrix <- function(context, samples,
                                                     parameters) {

  samples <- as.matrix(samples)
  valid   <- rep(TRUE, nrow(samples))
  object  <- context[["object"]]
  if (length(parameters) == 0L || is.null(object) ||
      is.null(context[["data"]]) || !.is_data_random(context[["data"]])) {
    return(list(samples = samples, valid = valid))
  }

  design <- .fitted_formula_design(object, "mu", required = FALSE)
  terms  <- design[["random_effects"]]
  if (length(terms) == 0L) {
    return(list(samples = samples, valid = valid))
  }

  K <- nrow(context[["data"]][["outcome"]])
  for (term in terms) {
    if (!.marginalized_random_effect_has_allocation(term)) {
      next
    }

    dependencies <- .iwmde_random_allocation_sd_dependencies(term)
    if (length(intersect(parameters, dependencies)) == 0L) {
      next
    }

    derived <- .marginalized_random_effect_allocated_sd_samples(
      term              = term,
      posterior_samples = samples,
      K                 = K
    )
    columns <- .iwmde_random_allocation_sd_columns(
      term      = term,
      samples   = samples,
      n_columns = ncol(derived)
    )
    if (length(columns) == 0L) {
      next
    }
    if (length(columns) != ncol(derived) || anyNA(columns)) {
      stop(
        "Cannot synchronize allocation-derived random-effect SD columns for ",
        "block '", term[["block_name"]], "'.",
        call. = FALSE
      )
    }

    finite <- rowSums(!is.finite(derived) | derived < 0) == 0L
    valid  <- valid & finite
    if (any(finite)) {
      samples[finite, columns] <- derived[finite, , drop = FALSE]
    }
  }

  return(list(samples = samples, valid = valid))
}


.iwmde_sync_random_allocation_sd_row <- function(context, row, parameters) {

  synced <- .iwmde_sync_random_allocation_sd_matrix(
    context    = context,
    samples    = matrix(
      as.numeric(row),
      nrow     = 1L,
      dimnames = list(NULL, names(row))
    ),
    parameters = parameters
  )

  return(list(
    row   = synced[["samples"]][1L, ],
    valid = synced[["valid"]][[1L]]
  ))
}


.iwmde_sync_replacement_row <- function(context, row, parameters) {

  synced <- .iwmde_sync_invgamma_auxiliary_row(
    context    = context,
    row        = row,
    parameters = parameters
  )
  if (!isTRUE(synced[["valid"]])) {
    return(synced)
  }

  .iwmde_sync_random_allocation_sd_row(
    context    = context,
    row        = synced[["row"]],
    parameters = parameters
  )
}


.iwmde_random_allocation_sd_dependencies <- function(term) {

  allocation <- term[["sd_binding"]][["allocations"]][[1L]]
  source     <- allocation[["source"]][["name"]]
  factors    <- .marginalized_random_effect_allocation_factors(
    term,
    all = TRUE
  )
  factor_columns <- vapply(factors, function(factor) {
    paste0(factor[["weight_name"]], "[", factor[["index"]], "]")
  }, character(1))

  dependencies <- unique(c(source, factor_columns))

  dependencies[!is.na(dependencies) & nzchar(dependencies)]
}


.iwmde_random_allocation_sd_columns <- function(term, samples, n_columns) {

  parameters <- unique(term[["sd_parameter_names"]])
  parameters <- parameters[!is.na(parameters) & nzchar(parameters)]
  if (length(parameters) == n_columns &&
      all(parameters %in% colnames(samples))) {
    return(parameters)
  }
  if (length(parameters) == 1L && n_columns > 1L) {
    indexed <- paste0(parameters, "[", seq_len(n_columns), "]")
    if (all(indexed %in% colnames(samples))) {
      return(indexed)
    }
  }

  present <- parameters[parameters %in% colnames(samples)]
  indexed_present <- if (length(parameters) == 1L) {
    colnames(samples)[BayesTools::JAGS_indexed_parameter_columns(
      columns   = colnames(samples),
      parameter = parameters
    )]
  } else {
    character()
  }
  if (length(c(present, indexed_present)) > 0L) {
    return(NA_character_)
  }

  character()
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
    valid <- linear[["current"]] == value
    return(out(row, valid = valid))
  }

  row[linear[["active_columns"]]] <- row[linear[["active_columns"]]] +
    (value - linear[["current"]]) * linear[["coefficients"]]
  synced <- .iwmde_sync_replacement_row(
    context    = context,
    row        = row,
    parameters = linear[["active_columns"]]
  )

  return(out(synced[["row"]], valid = synced[["valid"]]))
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
  current <- .iwmde_linear_value_row(context, row, weights)
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
    likelihood_row <- if (identical(.iwmde_state_scope_value(state), "global") &&
                          identical(state[["likelihood_mode"]], "marginal") &&
                          !.iwmde_marginal_likelihood_requires_row(context)) {
      NULL
    } else {
      row
    }
    log_lik <- .iwmde_log_likelihood_parameters(
      context         = context,
      parameters      = parameters,
      active_setup    = state[["active_setup"]],
      likelihood_mode = state[["likelihood_mode"]],
      row             = likelihood_row
    )
    log_lik <- .iwmde_scalar_log_density(log_lik)
    focal_log_prior <- .iwmde_focal_log_prior(
      prior     = state[["focal_prior"]],
      value     = value,
      parameter = parameter
    )
    focal_log_prior <- .iwmde_scalar_log_density(focal_log_prior)
    baseline_log_prior <- .iwmde_scalar_log_density(
      state[["baseline_log_prior"]]
    )
    baseline_focal_log_prior <- .iwmde_scalar_log_density(
      state[["baseline_focal_log_prior"]]
    )
    log_prior <- baseline_log_prior +
      focal_log_prior - baseline_focal_log_prior
    out <- log_lik + log_prior

    return(.iwmde_scalar_log_density(out))
  }

  return(.iwmde_log_q_state(
    context         = context,
    row             = row,
    active_setup    = state[["active_setup"]],
    parameters      = parameters,
    prior_list      = state[["prior_list"]],
    likelihood_mode = state[["likelihood_mode"]],
    replacement     = replacement
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

  if (is.null(state[["parameters"]])) {
    state[["parameters"]] <- .iwmde_row_parameters(
      context      = context,
      row          = state[["row"]],
      active_setup = state[["active_setup"]],
      state_scope  = .iwmde_state_scope_value(state)
    )
  }

  if (identical(replacement[["type"]], "linear")) {
    return(.iwmde_replace_linear_parameters(
      context     = context,
      state       = state,
      value       = value,
      replacement = replacement
    ))
  }

  if (replacement[["type"]] %in% c(
    "simplex_pair",
    "random_component_sd"
  )) {
    replaced <- .iwmde_replace_row_for_value(
      context     = context,
      state       = state,
      parameter   = parameter,
      value       = value,
      replacement = replacement
    )
    if (!isTRUE(replaced[["valid"]])) {
      return(out(replaced[["row"]], state[["parameters"]], valid = FALSE))
    }
    parameters <- .iwmde_row_parameters(
      context      = context,
      row          = replaced[["row"]],
      active_setup = state[["active_setup"]],
      state_scope  = .iwmde_state_scope_value(state)
    )

    return(out(
      row                   = replaced[["row"]],
      parameters            = parameters,
      use_focal_prior_delta = FALSE
    ))
  }

  if (identical(replacement[["type"]], "scalar")) {
    parameters <- state[["parameters"]]
    name       <- replacement[["name"]]
    if (!is.null(parameters[[name]])) {
      parameters[[name]] <- value
      synced <- .iwmde_sync_replacement_row(
        context    = context,
        row        = row,
        parameters = name
      )
      if (!isTRUE(synced[["valid"]])) {
        return(out(synced[["row"]], parameters, valid = FALSE))
      }
      parameters <- .iwmde_row_parameters(
        context      = context,
        row          = synced[["row"]],
        active_setup = state[["active_setup"]],
        state_scope  = .iwmde_state_scope_value(state)
      )
      return(out(synced[["row"]], parameters))
    }
  }

  if (identical(replacement[["type"]], "indexed")) {
    parameters <- state[["parameters"]]
    name       <- replacement[["name"]]
    index      <- replacement[["index"]]
    if (!is.null(parameters[[name]]) && length(parameters[[name]]) >= index) {
      parameters[[name]][index] <- value
      synced <- .iwmde_sync_replacement_row(
        context    = context,
        row        = row,
        parameters = parameter
      )
      if (!isTRUE(synced[["valid"]])) {
        return(out(synced[["row"]], parameters, valid = FALSE))
      }
      parameters <- .iwmde_row_parameters(
        context      = context,
        row          = synced[["row"]],
        active_setup = state[["active_setup"]],
        state_scope  = .iwmde_state_scope_value(state)
      )
      return(out(synced[["row"]], parameters))
    }
  }

  synced <- .iwmde_sync_replacement_row(
    context    = context,
    row        = row,
    parameters = parameter
  )
  if (!isTRUE(synced[["valid"]])) {
    return(out(synced[["row"]], state[["parameters"]], valid = FALSE))
  }

  return(out(
    row        = synced[["row"]],
    parameters = .iwmde_row_parameters(
      context      = context,
      row          = synced[["row"]],
      active_setup = state[["active_setup"]],
      state_scope  = .iwmde_state_scope_value(state)
    )
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
    valid <- linear[["current"]] == value
    return(out(row, state[["parameters"]], valid = valid))
  }

  # Move along the minimum-norm direction that changes a' beta to the target
  # value while keeping the orthogonal complement fixed.
  row[linear[["active_columns"]]] <- row[linear[["active_columns"]]] +
    (value - linear[["current"]]) * linear[["coefficients"]]
  synced <- .iwmde_sync_replacement_row(
    context    = context,
    row        = row,
    parameters = linear[["active_columns"]]
  )
  if (!isTRUE(synced[["valid"]])) {
    return(out(synced[["row"]], state[["parameters"]], valid = FALSE))
  }
  parameters <- .iwmde_row_parameters(
    context      = context,
    row          = synced[["row"]],
    active_setup = state[["active_setup"]],
    state_scope  = .iwmde_state_scope_value(state)
  )

  return(out(synced[["row"]], parameters))
}


.iwmde_linear_value_row <- function(context, row, weights) {

  values <- vapply(names(weights), function(parameter) {
    .iwmde_parameter_value_row(context, row, parameter)
  }, numeric(1))

  return(sum(values * weights))
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
  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "simplex_pair")) {
    return(list(
      type              = "simplex_pair",
      parameter         = parameter_spec[["parameter"]],
      index             = parameter_spec[["index"]],
      n_targets         = parameter_spec[["n_targets"]],
      target_columns    = parameter_spec[["target_columns"]],
      auxiliary_columns = parameter_spec[["auxiliary_columns"]]
    ))
  }
  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "random_component_sd")) {
    return(list(
      type              = "random_component_sd",
      source_parameter  = parameter_spec[["source_parameter"]],
      factors           = parameter_spec[["factors"]],
      target_columns    = parameter_spec[["target_columns"]],
      factor_columns    = parameter_spec[["factor_columns"]],
      auxiliary_columns = parameter_spec[["auxiliary_columns"]]
    ))
  }

  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "primitive") &&
      identical(parameter_spec[["target_columns"]], parameter) &&
      parameter %in% colnames(context[["posterior_samples"]])) {
    return(list(type = "scalar", name = parameter))
  }

  if (.iwmde_parameter_is_eta(parameter)) {
    return(list(type = "fallback"))
  }
  formula_parameter <- .iwmde_predictor_formula_parameter(
    context = context,
    column  = parameter
  )
  if (isTRUE(formula_parameter %in% c("mu", "log_tau"))) {
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
