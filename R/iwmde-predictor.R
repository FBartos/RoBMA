# ============================================================================ #
# IWMDE Predictor and Location Fast Paths
# ============================================================================ #

.iwmde_predictor_setup <- function(context, row_states, active_setup, unit) {

  rows <- vapply(row_states, function(state) {
    state[["row_index"]]
  }, integer(1))
  state_scope <- unique(vapply(
    row_states,
    .iwmde_state_scope_value,
    character(1)
  ))
  if (length(state_scope) != 1L) {
    state_scope <- "mixed"
  }
  key <- .iwmde_predictor_cache_key(
    prefix      = "predictor_setup",
    active_key  = .iwmde_state_active_key(context, row_states[[1L]]),
    unit        = unit,
    state_scope = state_scope,
    rows        = rows
  )

  if (exists(key, envir = context[["predictor_cache"]], inherits = FALSE)) {
    return(get(key, envir = context[["predictor_cache"]]))
  }

  samples <- context[["posterior_samples"]][rows, , drop = FALSE]
  if (identical(state_scope, "global")) {
    samples <- .iwmde_drop_local_latent_sample_columns(samples, context)
  }
  setup <- if (.iwmde_uses_known_v_random_marginal_likelihood(context) &&
               identical(unit, "estimate")) {
    .iwmde_predictor_setup_known_v_random(
      context           = context,
      posterior_samples = samples,
      active_setup      = active_setup,
      unit              = unit
    )
  } else {
    .iwmde_log_lik_posterior_setup_active_branch(
      context           = context,
      posterior_samples = samples,
      active_setup      = active_setup,
      unit              = unit,
      data_hash         = NULL
    )
  }

  assign(key, setup, envir = context[["predictor_cache"]])
  return(setup)
}


.iwmde_predictor_setup_known_v_random <- function(
    context, posterior_samples, active_setup, unit) {

  posterior_samples <- .iwmde_likelihood_posterior_samples(
    context      = context,
    samples      = posterior_samples,
    active_setup = active_setup
  )
  K <- nrow(context[["data"]][["outcome"]])

  return(.log_lik_posterior_setup(
    fit                        = context[["object"]][["fit"]],
    posterior_samples          = posterior_samples,
    data                       = context[["data"]],
    priors                     = active_setup[["priors"]],
    unit                       = unit,
    data_hash                  = NULL,
    conditioned_random_effects = matrix(
      0,
      nrow = nrow(posterior_samples),
      ncol = K
    )
  ))
}


.iwmde_predictor_cache_key <- function(prefix, ...) {

  bytes <- as.integer(serialize(list(...), NULL, version = 3))
  hash1 <- 5381
  hash2 <- 0

  for (byte in bytes) {
    hash1 <- (hash1 * 33 + byte) %% 2147483647
    hash2 <- (hash2 * 65599 + byte) %% 2147483629
  }

  return(paste(
    prefix,
    sprintf("%08x", as.integer(hash1)),
    sprintf("%08x", as.integer(hash2)),
    sep = "|"
  ))
}


.iwmde_predictor_update_basis <- function(context, parameter, row_states,
                                          replacement, setup) {

  if (identical(replacement[["type"]], "linear")) {
    return(.iwmde_predictor_linear_basis(
      context    = context,
      row_states = row_states,
      weights    = replacement[["weights"]],
      setup      = setup
    ))
  }

  if (!identical(replacement[["type"]], "scalar") &&
      !identical(replacement[["type"]], "fallback")) {
    return(NULL)
  }

  if (!parameter %in% colnames(setup[["posterior_samples"]])) {
    return(NULL)
  }

  basis <- .iwmde_predictor_column_basis(
    context = context,
    column  = parameter,
    state   = row_states[[1L]]
  )
  if (is.null(basis)) {
    return(NULL)
  }

  basis <- .iwmde_predictor_materialize_formula_basis(
    context = context,
    setup   = setup,
    basis   = basis,
    formula_mu_directions = if (isTRUE(basis[["formula_mu"]])) {
      matrix(
        1,
        nrow     = nrow(setup[["mu"]]),
        ncol     = 1L,
        dimnames = list(NULL, parameter)
      )
    } else {
      NULL
    },
    formula_logtau_directions = if (isTRUE(basis[["formula_logtau"]])) {
      matrix(
        1,
        nrow     = nrow(setup[["mu"]]),
        ncol     = 1L,
        dimnames = list(NULL, parameter)
      )
    } else {
      NULL
    }
  )
  basis <- .iwmde_predictor_expand_basis(
    basis = basis,
    S     = nrow(setup[["mu"]]),
    K     = ncol(setup[["mu"]])
  )
  basis[["current"]] <- as.numeric(setup[["posterior_samples"]][, parameter])
  return(basis)
}


.iwmde_predictor_linear_basis <- function(context, row_states, weights, setup) {

  S                            <- length(row_states)
  K                            <- ncol(setup[["mu"]])
  mu_basis                     <- matrix(0, nrow = S, ncol = K)
  log_basis                    <- matrix(0, nrow = S, ncol = K)
  has_mu                       <- FALSE
  has_logtau                   <- FALSE
  has_formula_mu               <- FALSE
  has_formula_logtau           <- FALSE
  formula_mu_directions        <- matrix(0, nrow = S, ncol = 0L)
  formula_logtau_directions    <- matrix(0, nrow = S, ncol = 0L)
  current                      <- numeric(S)

  for (i in seq_along(row_states)) {
    linear <- .iwmde_linear_replacement_state(
      context     = context,
      state       = row_states[[i]],
      replacement = list(type = "linear", weights = weights)
    )
    if (!isTRUE(linear[["valid"]])) {
      return(NULL)
    }

    current[i]     <- linear[["current"]]
    active_columns <- linear[["active_columns"]]
    if (length(active_columns) == 0L) {
      return(NULL)
    }

    active_weights <- linear[["active_weights"]]
    denominator    <- linear[["denominator"]]

    for (column in active_columns) {
      column_basis <- .iwmde_predictor_column_basis(
        context = context,
        column  = column,
        state   = row_states[[i]]
      )
      if (is.null(column_basis) ||
          !identical(column_basis[["scale_update"]], "none")) {
        return(NULL)
      }

      coefficient <- active_weights[[column]] / denominator
      if (isTRUE(column_basis[["formula_mu"]])) {
        has_formula_mu <- TRUE
        if (!column %in% colnames(formula_mu_directions)) {
          formula_mu_directions <- cbind(
            formula_mu_directions,
            matrix(
              0,
              nrow     = S,
              ncol     = 1L,
              dimnames = list(NULL, column)
            )
          )
        }
        formula_mu_directions[i, column] <- coefficient
        next
      }
      if (isTRUE(column_basis[["formula_logtau"]])) {
        has_formula_logtau <- TRUE
        if (!column %in% colnames(formula_logtau_directions)) {
          formula_logtau_directions <- cbind(
            formula_logtau_directions,
            matrix(
              0,
              nrow     = S,
              ncol     = 1L,
              dimnames = list(NULL, column)
            )
          )
        }
        formula_logtau_directions[i, column] <- coefficient
        next
      }

      if (!is.null(column_basis[["mu_basis"]])) {
        mu_basis[i, ] <- mu_basis[i, ] + coefficient * column_basis[["mu_basis"]]
        has_mu        <- TRUE
      }
      if (!is.null(column_basis[["log_tau_basis"]])) {
        log_basis[i, ] <- log_basis[i, ] + coefficient * column_basis[["log_tau_basis"]]
        has_logtau     <- TRUE
      }
    }
  }

  if ((has_formula_mu || has_formula_logtau) && (has_mu || has_logtau)) {
    return(NULL)
  }
  if (!has_mu && !has_logtau) {
    if (!has_formula_mu && !has_formula_logtau) {
      return(NULL)
    }
  }

  basis <- list(
    current                   = current,
    mu_basis                  = if (has_mu) mu_basis else NULL,
    log_tau_basis             = if (has_logtau) log_basis else NULL,
    scale_update              = "none",
    formula_mu                = has_formula_mu,
    formula_logtau            = has_formula_logtau,
    formula_mu_columns        = if (has_formula_mu) {
      colnames(formula_mu_directions)
    } else {
      NULL
    },
    formula_logtau_columns    = if (has_formula_logtau) {
      colnames(formula_logtau_directions)
    } else {
      NULL
    }
  )
  basis <- .iwmde_predictor_materialize_formula_basis(
    context                    = context,
    setup                      = setup,
    basis                      = basis,
    formula_mu_directions      = if (has_formula_mu) {
      formula_mu_directions
    } else {
      NULL
    },
    formula_logtau_directions = if (has_formula_logtau) {
      formula_logtau_directions
    } else {
      NULL
    }
  )

  return(basis)
}


.iwmde_predictor_column_basis <- function(context, column, state) {

  data         <- context[["data"]]
  active_setup <- state[["active_setup"]]
  K            <- nrow(data[["outcome"]])

  formula_parameter <- .iwmde_predictor_formula_parameter(
    context = context,
    column  = column
  )
  if (identical(formula_parameter, "mu")) {
    return(list(
      mu_basis           = NULL,
      log_tau_basis      = NULL,
      scale_update       = "none",
      formula_mu         = TRUE,
      formula_logtau     = FALSE,
      formula_mu_columns = column
    ))
  }
  if (identical(formula_parameter, "log_tau")) {
    return(list(
      mu_basis              = NULL,
      log_tau_basis         = NULL,
      scale_update          = "none",
      formula_mu            = FALSE,
      formula_logtau        = TRUE,
      formula_logtau_columns = column
    ))
  }

  if (column %in% c("mu", "mu_intercept")) {
    return(list(
      mu_basis      = rep(1, K),
      log_tau_basis = NULL,
      scale_update  = "none"
    ))
  }

  if (column == "PET" && isTRUE(active_setup[["is_PET"]])) {
    direction <- ifelse(.data_effect_direction(data) == "negative", -1, 1)
    return(list(
      mu_basis      = direction * data[["outcome"]][["sei"]],
      log_tau_basis = NULL,
      scale_update  = "none"
    ))
  }

  if (column == "PEESE" && isTRUE(active_setup[["is_PEESE"]])) {
    direction <- ifelse(.data_effect_direction(data) == "negative", -1, 1)
    return(list(
      mu_basis      = direction * data[["outcome"]][["sei"]]^2,
      log_tau_basis = NULL,
      scale_update  = "none"
    ))
  }

  if (column == "tau" && !.is_data_scale(data)) {
    if (.is_data_random(data) && .is_data_known_v(data)) {
      return(NULL)
    }
    return(list(
      mu_basis      = NULL,
      log_tau_basis = NULL,
      scale_update  = "tau"
    ))
  }

  if (column == "log_tau" && !.is_data_scale(data)) {
    return(list(
      mu_basis      = NULL,
      log_tau_basis = rep(1, K),
      scale_update  = "none"
    ))
  }

  if (column == "rho" && .is_data_multilevel(data)) {
    return(list(
      mu_basis      = NULL,
      log_tau_basis = NULL,
      scale_update  = "rho"
    ))
  }

  return(NULL)
}


.iwmde_predictor_formula_parameter <- function(context, column) {

  cache <- context[["predictor_cache"]]
  key   <- paste0("formula_parameter|", column)
  if (is.environment(cache) &&
      exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE)[["value"]])
  }

  fit <- context[["formula_fit"]]
  if (is.null(fit) && !is.null(context[["object"]])) {
    fit <- context[["object"]][["fit"]]
  }
  if (is.null(fit) || !inherits(fit, "BayesTools_fit")) {
    return(NULL)
  }
  coordinates <- BayesTools::parameter_coordinates(fit)
  matches <- coordinates[["coordinate_name"]] == column &
    coordinates[["role"]] == "fixed_coefficient" &
    !coordinates[["internal"]]
  if (sum(matches) != 1L) {
    parameter <- NULL
  } else {
    parameter <- coordinates[["formula_parameter"]][matches]
    if (!parameter %in% c("mu", "log_tau")) {
      parameter <- NULL
    }
  }
  if (is.environment(cache)) {
    assign(key, list(value = parameter), envir = cache)
  }

  return(parameter)
}


.iwmde_predictor_materialize_formula_basis <- function(
    context, setup, basis, formula_mu_directions = NULL,
    formula_logtau_directions = NULL) {

  directions <- list(
    mu      = formula_mu_directions,
    log_tau = formula_logtau_directions
  )
  fit <- context[["formula_fit"]]
  if (is.null(fit) && !is.null(context[["object"]])) {
    fit <- context[["object"]][["fit"]]
  }
  for (formula_parameter in names(directions)) {
    parameter_directions <- directions[[formula_parameter]]
    if (is.null(parameter_directions)) {
      next
    }
    result <- BayesTools::JAGS_formula_predictor_basis(
      fit               = fit,
      directions        = parameter_directions,
      posterior_samples = setup[["posterior_samples"]]
    )
    if (!identical(result[["status"]], "affine")) {
      basis[[paste0("formula_", formula_parameter, "_status")]] <-
        result[["status"]]
      basis[[paste0("formula_", formula_parameter, "_reason")]] <-
        result[["reason"]]
      next
    }
    if (!identical(result[["parameter"]], formula_parameter)) {
      stop(
        "Formula predictor basis owner does not match the selected model component.",
        call. = FALSE
      )
    }
    if (identical(formula_parameter, "mu")) {
      basis[["mu_basis"]]           <- result[["basis"]]
      basis[["formula_mu"]]         <- FALSE
      basis[["formula_mu_columns"]] <- NULL
      basis[["formula_mu_affine"]]  <- TRUE
    } else {
      basis[["log_tau_basis"]]          <- result[["basis"]]
      basis[["formula_logtau"]]         <- FALSE
      basis[["formula_logtau_columns"]] <- NULL
    }
  }

  return(basis)
}


.iwmde_predictor_expand_basis <- function(basis, S, K) {

  if (!is.null(basis[["mu_basis"]]) && is.null(dim(basis[["mu_basis"]]))) {
    basis[["mu_basis"]] <- matrix(
      basis[["mu_basis"]],
      nrow  = S,
      ncol  = K,
      byrow = TRUE
    )
  }
  if (!is.null(basis[["log_tau_basis"]]) &&
      is.null(dim(basis[["log_tau_basis"]]))) {
    basis[["log_tau_basis"]] <- matrix(
      basis[["log_tau_basis"]],
      nrow  = S,
      ncol  = K,
      byrow = TRUE
    )
  }

  return(basis)
}


.iwmde_log_q_grid_normal_location_group <- function(context, parameter, values,
                                                    row_states, replacement,
                                                    setup, basis) {

  active_setup <- row_states[[1L]][["active_setup"]]
  is_weightfunction <- isTRUE(active_setup[["is_weightfunction"]])
  if (.data_outcome_type(context[["data"]]) != "norm") {
    return(NULL)
  }
  if (identical(basis[["scale_update"]], "rho")) {
    return(.iwmde_log_q_grid_normal_cluster_rho_group(
      context     = context,
      parameter   = parameter,
      values      = values,
      row_states  = row_states,
      replacement = replacement,
      setup       = setup,
      basis       = basis
    ))
  }
  if (identical(basis[["scale_update"]], "tau")) {
    return(.iwmde_log_q_grid_normal_known_v_tau_group(
      context     = context,
      parameter   = parameter,
      values      = values,
      row_states  = row_states,
      replacement = replacement,
      setup       = setup,
      basis       = basis
    ))
  }
  if (isTRUE(basis[["formula_mu"]]) ||
      isTRUE(basis[["formula_logtau"]]) ||
      !identical(basis[["scale_update"]], "none") ||
      !is.null(basis[["log_tau_basis"]]) ||
      is.null(basis[["mu_basis"]])) {
    return(NULL)
  }
  if (is_weightfunction && .is_data_multilevel(context[["data"]])) {
    return(.iwmde_log_q_grid_selnorm_multilevel_location_group(
      context      = context,
      parameter    = parameter,
      values       = values,
      row_states   = row_states,
      replacement  = replacement,
      setup        = setup,
      basis        = basis,
      active_setup = active_setup
    ))
  }

  likelihood_change <- .iwmde_normal_location_likelihood_change_cached(
    context     = context,
    parameter   = parameter,
    row_states  = row_states,
    replacement = replacement,
    setup       = setup,
    basis       = basis
  )
  if (is.null(likelihood_change)) {
    return(NULL)
  }

  normalizer_change <- NULL
  if (is_weightfunction) {
    normalizer_change <- .iwmde_selected_normal_location_normalizer_change(
      context      = context,
      setup        = setup,
      basis        = basis,
      values       = values,
      row_states   = row_states,
      active_setup = active_setup
    )
    if (is.null(normalizer_change)) {
      return(NULL)
    }
  }

  log_prior <- .iwmde_predictor_log_prior(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )
  if (is.null(log_prior) ||
      length(log_prior) != length(values) * length(row_states)) {
    return(NULL)
  }

  G          <- length(values)
  S          <- length(row_states)
  row_index  <- rep(seq_len(S), each = G)
  grid_index <- rep(seq_len(G), times = S)
  delta      <- values[grid_index] - basis[["current"]][row_index]
  baseline <- vapply(row_states, function(state) {
    state[["baseline_log_lik"]]
  }, numeric(1))
  log_lik <- baseline[row_index] +
    likelihood_change[["linear"]][row_index] * delta -
    .5 * likelihood_change[["quadratic"]][row_index] * delta^2
  if (!is.null(normalizer_change)) {
    log_lik <- log_lik - normalizer_change
  }
  log_q <- log_lik + log_prior
  log_q[!is.finite(delta)] <- -Inf

  return(matrix(log_q, nrow = G, ncol = S))
}


.iwmde_log_q_grid_normal_cluster_rho_group <- function(
    context, parameter, values, row_states, replacement, setup, basis) {

  if (!identical(parameter, "rho") ||
      .is_data_known_v(context[["data"]]) ||
      !.is_data_multilevel(context[["data"]]) ||
      isTRUE(setup[["is_weightfunction"]]) ||
      !is.null(setup[["weights"]]) ||
      !identical(basis[["scale_update"]], "rho") ||
      isTRUE(basis[["formula_mu"]]) ||
      isTRUE(basis[["formula_logtau"]]) ||
      !is.null(basis[["mu_basis"]]) ||
      !is.null(basis[["log_tau_basis"]])) {
    return(NULL)
  }

  G <- length(values)
  S <- length(row_states)
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

  valid   <- is.finite(values) & values >= 0 & values <= 1
  log_lik <- matrix(-Inf, nrow = G, ncol = S)
  if (any(valid)) {
    evaluated <- .log_lik_cluster_norm_analytic_rho_grid_sum(
      setup = setup,
      yi    = setup[["yi"]],
      vi    = setup[["sei"]]^2,
      rho   = values[valid]
    )
    if (!is.matrix(evaluated) ||
        !identical(dim(evaluated), c(sum(valid), S)) ||
        any(!is.finite(evaluated))) {
      return(NULL)
    }
    log_lik[valid, ] <- evaluated
  }

  return(log_lik + matrix(log_prior, nrow = G, ncol = S))
}


.iwmde_normal_location_likelihood_change_cached <- function(
    context, parameter, row_states, replacement, setup, basis) {

  cache <- context[["predictor_cache"]]
  if (!.iwmde_uses_known_v_joint_likelihood(context) ||
      !is.environment(cache)) {
    return(.iwmde_normal_location_likelihood_change(
      context = context,
      setup   = setup,
      basis   = basis
    ))
  }

  rows <- vapply(row_states, function(state) {
    state[["row_index"]]
  }, integer(1))
  state_scope <- unique(vapply(
    row_states,
    .iwmde_state_scope_value,
    character(1)
  ))
  if (length(state_scope) != 1L) {
    state_scope <- "mixed"
  }
  key <- .iwmde_predictor_cache_key(
    prefix      = "known_v_location_likelihood_change",
    active_key  = .iwmde_state_active_key(context, row_states[[1L]]),
    parameter   = parameter,
    replacement = replacement,
    state_scope = state_scope,
    rows        = rows
  )
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }

  likelihood_change <- .iwmde_normal_location_likelihood_change(
    context = context,
    setup   = setup,
    basis   = basis
  )
  assign(key, likelihood_change, envir = cache)

  return(likelihood_change)
}


.iwmde_log_q_grid_selnorm_multilevel_location_group <- function(context,
                                                                 parameter,
                                                                 values,
                                                                 row_states,
                                                                 replacement,
                                                                 setup, basis,
                                                                 active_setup) {

  selection_context <- .iwmde_selection_context_active_branch(
    context           = context,
    active_setup      = active_setup,
    posterior_samples = setup[["posterior_samples"]]
  )
  if (is.null(selection_context)) {
    return(NULL)
  }

  log_prior <- .iwmde_predictor_log_prior(
    context     = context,
    parameter   = parameter,
    values      = values,
    row_states  = row_states,
    replacement = replacement
  )
  if (is.null(log_prior) ||
      length(log_prior) != length(values) * length(row_states)) {
    return(NULL)
  }

  log_lik <- .log_lik_cluster_selnorm_location_grid(
    setup             = setup,
    yi                = setup[["yi"]],
    sei               = setup[["sei"]],
    basis             = basis[["mu_basis"]],
    current           = basis[["current"]],
    values            = values,
    selection_context = selection_context
  )
  if (!is.matrix(log_lik) ||
      nrow(log_lik) != length(values) ||
      ncol(log_lik) != length(row_states)) {
    return(NULL)
  }
  quadrature_change <- attr(
    log_lik,
    "quadrature_relative_change",
    exact = TRUE
  )

  G          <- length(values)
  S          <- length(row_states)
  row_index  <- rep(seq_len(S), each = G)
  grid_index <- rep(seq_len(G), times = S)
  delta      <- values[grid_index] - basis[["current"]][row_index]
  log_q      <- as.vector(log_lik) + log_prior
  log_q[!is.finite(delta)] <- -Inf
  out <- matrix(log_q, nrow = G, ncol = S)
  if (!is.null(quadrature_change)) {
    quadrature_change <- suppressWarnings(as.numeric(quadrature_change))
    quadrature_change <- quadrature_change[is.finite(quadrature_change)]
    if (length(quadrature_change) > 0L) {
      attr(out, "max_quadrature_relative_change") <-
        max(quadrature_change)
    }
  }

  return(out)
}


.iwmde_selected_normal_location_normalizer_change <- function(context, setup,
                                                              basis, values,
                                                              row_states,
                                                              active_setup) {

  if (!is.null(basis[["log_tau_basis"]]) ||
      !identical(basis[["scale_update"]], "none") ||
      isTRUE(basis[["formula_logtau"]])) {
    return(NULL)
  }

  selection_context <- .iwmde_selection_context_active_branch(
    context           = context,
    active_setup      = active_setup,
    posterior_samples = setup[["posterior_samples"]]
  )
  if (is.null(selection_context) ||
      any(!selection_context[["kernel_mode"]] %in% c(SELKERNEL_NORMAL, SELKERNEL_STEP))) {
    return(NULL)
  }

  G          <- length(values)
  S          <- length(row_states)
  K          <- ncol(setup[["mu"]])
  row_index  <- rep(seq_len(S), each = G)
  grid_index <- rep(seq_len(G), times = S)
  delta      <- values[grid_index] - basis[["current"]][row_index]
  sei        <- matrix(setup[["sei"]], nrow = S, ncol = K, byrow = TRUE)
  sd         <- .root_sum_squares(setup[["tau_within"]], sei)

  current_log_norm <- .iwmde_selected_normal_current_log_norm(
    context           = context,
    setup             = setup,
    sd                = sd,
    selection_context = selection_context,
    row_states        = row_states
  )

  if (.has_native_selnorm_log_norm_delta()) {
    out <- .selnorm_kernel_log_norm_delta_grid(
      mean             = setup[["mu"]],
      sd               = sd,
      basis            = basis[["mu_basis"]],
      current_log_norm = current_log_norm,
      current          = basis[["current"]],
      values           = values,
      sei              = setup[["sei"]],
      weights          = setup[["weights"]],
      omega            = selection_context[["omega"]],
      selection_spec   = selection_context,
      alpha            = selection_context[["alpha"]],
      phack_kind       = selection_context[["phack_kind"]],
      kernel_mode      = selection_context[["kernel_mode"]]
    )
    if (!is.matrix(out) ||
        !identical(dim(out), c(G, S))) {
      return(NULL)
    }

    out[!is.finite(delta)] <- NA_real_
    if (any(is.na(out))) {
      return(NULL)
    }

    return(as.numeric(out))
  }

  mean <- setup[["mu"]][row_index, , drop = FALSE] +
    basis[["mu_basis"]][row_index, , drop = FALSE] * delta
  candidate_context <- BayesTools::selection_context_subset_rows(
    context = selection_context,
    rows    = row_index
  )
  candidate_log_norm <- .selection_step_log_norm_matrix(
    mean              = mean,
    sd                = sd[row_index, , drop = FALSE],
    sei               = setup[["sei"]],
    selection_context = candidate_context
  )

  log_norm_change <- candidate_log_norm -
    current_log_norm[row_index, , drop = FALSE]
  if (!is.null(setup[["weights"]])) {
    log_norm_change <- .apply_log_lik_weights(log_norm_change, setup[["weights"]])
  }

  out <- matrix(rowSums(log_norm_change), nrow = G, ncol = S)
  out[!is.finite(delta)] <- NA_real_

  if (any(is.na(out))) {
    return(NULL)
  }

  return(as.numeric(out))
}


.iwmde_selected_normal_current_log_norm <- function(context, setup, sd,
                                                    selection_context,
                                                    row_states) {

  rows <- vapply(row_states, function(state) {
    state[["row_index"]]
  }, integer(1))
  key <- .iwmde_predictor_cache_key(
    prefix                  = "selnorm_current_log_norm",
    active_key              = .iwmde_state_active_key(context, row_states[[1L]]),
    rows                    = rows,
    n_mu                    = ncol(setup[["mu"]]),
    n_bins                  = selection_context[["n_bins"]],
    telescope_probabilities = selection_context[["telescope_probabilities"]]
  )

  if (exists(key, envir = context[["predictor_cache"]], inherits = FALSE)) {
    return(get(key, envir = context[["predictor_cache"]], inherits = FALSE))
  }

  current_log_norm <- .selection_step_log_norm_matrix(
    mean              = setup[["mu"]],
    sd                = sd,
    sei               = setup[["sei"]],
    selection_context = selection_context
  )

  assign(key, current_log_norm, envir = context[["predictor_cache"]])
  return(current_log_norm)
}


.iwmde_normal_location_likelihood_change <- function(context, setup, basis) {

  data      <- context[["data"]]
  yi        <- setup[["yi"]]
  mu        <- setup[["mu"]]
  mu_basis  <- basis[["mu_basis"]]
  sei       <- setup[["sei"]]
  data_weights <- setup[["weights"]]

  if (setup[["effect_direction"]] == "negative") {
    yi       <- -yi
    mu       <- -mu
    mu_basis <- -mu_basis
  }

  if (.iwmde_uses_known_v_joint_likelihood(context)) {
    return(.iwmde_normal_location_likelihood_change_known_v(
      context  = context,
      setup    = setup,
      yi       = yi,
      mu       = mu,
      mu_basis = mu_basis
    ))
  }

  if (!.is_data_multilevel(data)) {
    return(.iwmde_normal_location_likelihood_change_estimate(
      yi           = yi,
      mu           = mu,
      mu_basis     = mu_basis,
      tau_within   = setup[["tau_within"]],
      sei          = sei,
      data_weights = data_weights
    ))
  }

  if (!is.null(data_weights)) {
    return(NULL)
  }

  return(.iwmde_normal_location_likelihood_change_cluster(
    yi          = yi,
    mu          = mu,
    mu_basis    = mu_basis,
    tau_within  = setup[["tau_within"]],
    tau_between = setup[["tau_between"]],
    sei         = sei,
    cluster     = setup[["cluster"]]
  ))
}


.iwmde_normal_location_likelihood_change_estimate <- function(yi, mu,
                                                              mu_basis,
                                                              tau_within,
                                                              sei,
                                                              data_weights) {

  S        <- nrow(mu)
  K        <- ncol(mu)
  variance <- tau_within^2 + matrix(sei^2, nrow = S, ncol = K, byrow = TRUE)
  if (!all(is.finite(variance)) || any(variance <= 0)) {
    return(NULL)
  }

  residual <- matrix(yi, nrow = S, ncol = K, byrow = TRUE) - mu
  if (is.null(data_weights)) {
    weight <- matrix(1, nrow = S, ncol = K)
  } else {
    weight <- matrix(data_weights, nrow = S, ncol = K, byrow = TRUE)
  }

  return(list(
    linear    = rowSums(weight * residual * mu_basis / variance),
    quadratic = rowSums(weight * mu_basis^2 / variance)
  ))
}


.iwmde_normal_location_likelihood_change_cluster <- function(yi, mu,
                                                             mu_basis,
                                                             tau_within,
                                                             tau_between,
                                                             sei, cluster) {

  S         <- nrow(mu)
  linear    <- numeric(S)
  quadratic <- numeric(S)

  for (g in seq_along(cluster)) {
    idx       <- cluster[[g]]
    diagonal  <- sweep(tau_within[, idx, drop = FALSE]^2, 2L, sei[idx]^2, "+")
    rank_one  <- tau_between[, idx, drop = FALSE]
    residual  <- matrix(yi[idx], nrow = S, ncol = length(idx), byrow = TRUE) -
      mu[, idx, drop = FALSE]
    basis     <- mu_basis[, idx, drop = FALSE]
    valid     <- rowSums(!is.finite(diagonal) | diagonal <= 0) == 0L
    if (!all(valid)) {
      return(NULL)
    }

    inv_diag <- 1 / diagonal
    denom    <- 1 + rowSums(rank_one^2 * inv_diag)
    if (any(!is.finite(denom) | denom <= 0)) {
      return(NULL)
    }

    basis_inv    <- basis * inv_diag
    residual_inv <- residual * inv_diag
    cross_residual <- rowSums(basis * residual_inv) -
      rowSums(rank_one * basis_inv) *
      rowSums(rank_one * residual_inv) / denom
    cross_basis <- rowSums(basis * basis_inv) -
      rowSums(rank_one * basis_inv)^2 / denom

    if (any(!is.finite(cross_residual)) || any(!is.finite(cross_basis))) {
      return(NULL)
    }

    linear    <- linear + cross_residual
    quadratic <- quadratic + cross_basis
  }

  return(list(
    linear    = linear,
    quadratic = quadratic
  ))
}


.iwmde_predictor_candidates <- function(context, active_setup, setup, basis,
                                        parameter, values, row_states,
                                        replacement) {

  S         <- nrow(setup[["mu"]])
  K         <- ncol(setup[["mu"]])
  G         <- length(values)
  row_index <- rep(seq_len(S), each = G)
  grid_index <- rep(seq_len(G), times = S)
  delta     <- values[grid_index] - basis[["current"]][row_index]

  replacement_samples <- NULL
  if (isTRUE(basis[["formula_mu"]]) ||
      isTRUE(basis[["formula_logtau"]])) {
    replacement_samples <- .iwmde_build_replacement_samples(
      context     = context,
      parameter   = parameter,
      values      = values,
      row_states  = row_states,
      replacement = replacement
    )
  }

  mu <- setup[["mu"]][row_index, , drop = FALSE]
  if (!is.null(basis[["mu_basis"]])) {
    mu <- mu + basis[["mu_basis"]][row_index, , drop = FALSE] * delta
  }
  if (isTRUE(basis[["formula_mu"]])) {
    mu <- .iwmde_predictor_evaluate_mu(
      context      = context,
      active_setup = active_setup,
      samples      = replacement_samples[["samples"]]
    )
  }

  tau_total     <- setup[["tau_total"]]
  is_multilevel <- isTRUE(setup[["is_multilevel"]])
  rho           <- .iwmde_predictor_rho(setup)
  valid         <- is.finite(values[grid_index]) & is.finite(delta)

  if (identical(basis[["scale_update"]], "tau")) {
    tau_value <- values[grid_index]
    valid     <- valid & tau_value >= 0
    tau_total <- matrix(tau_value, nrow = length(row_index), ncol = K)
  } else {
    tau_total <- tau_total[row_index, , drop = FALSE]
  }

  if (identical(basis[["scale_update"]], "rho")) {
    rho_value <- values[grid_index]
    valid     <- valid & is.finite(rho_value) &
      rho_value >= 0 & rho_value <= 1
    rho       <- rho_value
  } else {
    rho <- rho[row_index]
  }

  if (!is.null(basis[["log_tau_basis"]])) {
    log_tau <- log(tau_total)
    log_tau <- log_tau + basis[["log_tau_basis"]][row_index, , drop = FALSE] * delta
    valid   <- valid & is.finite(rowSums(log_tau))
    tau_total <- exp(log_tau)
  }
  if (isTRUE(basis[["formula_logtau"]])) {
    tau_result <- .iwmde_predictor_evaluate_tau(
      context      = context,
      active_setup = active_setup,
      samples      = replacement_samples[["samples"]]
    )
    tau_within  <- tau_result[["tau_within"]]
    tau_between <- tau_result[["tau_between"]]
  } else {
    valid <- valid &
      rowSums(!is.finite(tau_total) | tau_total < 0) == 0L
    tau_within  <- matrix(NA_real_, nrow = nrow(tau_total), ncol = K)
    tau_between <- matrix(NA_real_, nrow = nrow(tau_total), ncol = K)
    valid_rows  <- which(valid)
    if (length(valid_rows) > 0L) {
      tau_result <- .heterogeneity_components(
        tau_total     = tau_total[valid_rows, , drop = FALSE],
        rho           = if (is_multilevel) rho[valid_rows] else NULL,
        is_multilevel = is_multilevel,
        context       = "IWMDE predictor candidates"
      )
      tau_within[valid_rows, ]  <- tau_result[["tau_within"]]
      tau_between[valid_rows, ] <- tau_result[["tau_between"]]
    }
  }

  if (!is.null(replacement_samples)) {
    valid <- valid & replacement_samples[["valid"]]
    posterior_samples <- replacement_samples[["samples"]]
  } else {
    posterior_samples <- setup[["posterior_samples"]][row_index, , drop = FALSE]
  }

  valid <- valid & is.finite(rowSums(mu)) &
    is.finite(rowSums(tau_within)) &
    is.finite(rowSums(tau_between))

  return(list(
    mu                  = mu,
    tau_within          = tau_within,
    tau_between         = tau_between,
    row_index           = row_index,
    grid_index          = grid_index,
    valid               = valid,
    posterior_samples   = posterior_samples,
    replacement_samples = replacement_samples
  ))
}


.iwmde_predictor_evaluate_mu <- function(context, active_setup, samples) {

  mu_samples <- .iwmde_predictor_evaluate_fixed_mu(
    context      = context,
    active_setup = active_setup,
    samples      = samples
  )

  data <- context[["data"]]
  if (.is_data_known_v(data) &&
      length(.data_effective_sampled_random_effect_terms(data)) > 0L) {
    mu_samples <- mu_samples + .evaluate.brma.random_effects(
      fit               = context[["object"]][["fit"]],
      data              = data,
      priors            = active_setup[["priors"]],
      posterior_samples = samples,
      same_data         = TRUE,
      required          = TRUE
    )
  }

  return(mu_samples)
}


.iwmde_predictor_evaluate_fixed_mu <- function(context, active_setup, samples) {

  data         <- context[["data"]]
  mods_formula <- context[["formula_inputs"]][["mods"]][["formula"]]

  mu_samples <- .evaluate.brma.mu(
    fit               = context[["object"]][["fit"]],
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = mods_formula,
    mods_priors       = if (.is_data_random(data)) {
      active_setup[["priors"]][["location"]]
    } else {
      active_setup[["priors"]][["mods"]]
    },
    is_mods           = .is_data_mods(data),
    is_PET            = active_setup[["is_PET"]],
    is_PEESE          = active_setup[["is_PEESE"]],
    effect_direction  = .data_effect_direction(data),
    bias_adjusted     = FALSE,
    K                 = nrow(data[["outcome"]]),
    posterior_samples = samples,
    priors            = active_setup[["priors"]]
  )

  return(mu_samples)
}


.iwmde_predictor_evaluate_tau <- function(context, active_setup, samples) {

  data          <- context[["data"]]
  scale_formula <- context[["formula_inputs"]][["scale"]][["formula"]]

  return(.evaluate.brma.tau(
    fit               = context[["object"]][["fit"]],
    scale_data        = data[["scale"]],
    scale_formula     = scale_formula,
    scale_priors      = active_setup[["priors"]][["scale"]],
    is_scale          = .is_data_scale(data),
    is_multilevel     = .is_data_multilevel(data),
    K                 = nrow(data[["outcome"]]),
    posterior_samples = samples,
    fixed_tau         = .fixed_tau_prior_value(active_setup[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(active_setup[["priors"]])
  ))
}


.iwmde_predictor_rho <- function(setup) {

  if (!isTRUE(setup[["is_multilevel"]])) {
    return(rep(0, nrow(setup[["mu"]])))
  }

  return(.resolve_heterogeneity_allocation(
    rho       = setup[["rho"]],
    n_samples = nrow(setup[["mu"]]),
    context   = "IWMDE predictor setup"
  ))
}


.iwmde_predictor_log_prior <- function(context, parameter, values,
                                       row_states, replacement,
                                       replacement_samples = NULL) {

  G <- length(values)
  S <- length(row_states)

  use_delta <- !identical(replacement[["type"]], "linear") &&
    all(vapply(row_states, function(state) {
      isTRUE(state[["use_focal_prior_delta"]])
    }, logical(1)))

  if (use_delta) {
    focal_log_prior <- .iwmde_focal_log_prior_values(
      prior     = row_states[[1L]][["focal_prior"]],
      values    = values,
      parameter = parameter
    )
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
    out <- sweep(
      matrix(focal_log_prior, nrow = G, ncol = S),
      2L,
      baseline_log_prior,
      "+"
    )
    out <- sweep(out, 2L, baseline_focal_log_prior, "-")

    return(as.numeric(out))
  }

  if (identical(replacement[["type"]], "linear")) {
    out <- .iwmde_predictor_linear_log_prior_delta(
      context     = context,
      values      = values,
      row_states  = row_states,
      replacement = replacement
    )
    if (!is.null(out)) {
      return(out)
    }
  }

  candidates <- replacement_samples
  if (is.null(candidates)) {
    candidates <- .iwmde_build_replacement_samples(
      context     = context,
      parameter   = parameter,
      values      = values,
      row_states  = row_states,
      replacement = replacement
    )
  }
  out <- rep(-Inf, G * S)
  valid_positions <- which(candidates[["valid"]])
  if (length(valid_positions) == 0L) {
    return(out)
  }

  out[valid_positions] <- .iwmde_replacement_log_prior(
    parameter       = parameter,
    values          = values,
    valid_samples   = candidates[["samples"]][valid_positions, , drop = FALSE],
    valid_positions = valid_positions,
    candidates      = candidates,
    row_states      = row_states,
    replacement     = replacement
  )

  return(out)
}


.iwmde_predictor_linear_log_prior_delta <- function(context, values,
                                                    row_states,
                                                    replacement) {

  G       <- length(values)
  out     <- rep(-Inf, G * length(row_states))
  weights <- replacement[["weights"]]

  for (i in seq_along(row_states)) {
    state   <- row_states[[i]]
    row     <- state[["row"]]
    idx     <- (i - 1L) * G + seq_len(G)
    current <- .iwmde_linear_value_row(context, row, weights)
    if (!is.finite(current)) {
      next
    }

    active_columns <- .iwmde_linear_active_columns(context, row, weights)
    if (length(active_columns) == 0L) {
      valid <- values == current
      out[idx[valid]] <- state[["baseline_log_prior"]]
      next
    }

    active_weights <- weights[active_columns]
    denominator    <- sum(active_weights^2)
    if (!is.finite(denominator) || denominator <= 0) {
      return(NULL)
    }

    log_prior <- rep(state[["baseline_log_prior"]], G)
    for (column in active_columns) {
      prior <- .iwmde_focal_prior(context, column, row)
      if (!.iwmde_can_use_focal_prior_delta(prior)) {
        return(NULL)
      }
      if (is.null(prior)) {
        next
      }

      baseline <- .iwmde_focal_log_prior(
        prior     = prior,
        value     = row[[column]],
        parameter = column
      )
      if (!is.finite(baseline)) {
        return(NULL)
      }

      candidate_values <- row[[column]] +
        (values - current) * active_weights[[column]] / denominator
      candidate_log_prior <- .iwmde_focal_log_prior_values(
        prior     = prior,
        values    = candidate_values,
        parameter = column
      )
      log_prior <- log_prior + candidate_log_prior - baseline
    }

    out[idx] <- log_prior
  }

  return(out)
}
