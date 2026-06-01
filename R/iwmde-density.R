# ============================================================================ #
# IWMDE Row States, Density Aggregation, and Chen Weights
# ============================================================================ #

.iwmde_row_states <- function(context, rows, parameter = NULL,
                              parameter_spec = NULL) {

  lapply(rows, function(row) {
    .iwmde_row_state(
      context        = context,
      row_index      = row,
      parameter      = parameter,
      parameter_spec = parameter_spec
    )
  })
}


.iwmde_likelihood_mode <- function(parameter, parameter_spec = NULL) {

  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "linear")) {
    if (.iwmde_linear_uses_local_latent(parameter_spec[["weights"]])) {
      return("conditional")
    }

    return("marginal")
  }

  if (.iwmde_parameter_is_local_latent(parameter)) {
    return("conditional")
  }

  return("marginal")
}


.iwmde_linear_uses_local_latent <- function(weights) {

  return(any(vapply(
    names(weights),
    .iwmde_parameter_is_local_latent,
    logical(1)
  )))
}


.iwmde_parameter_is_local_latent <- function(parameter) {

  if (is.null(parameter) || length(parameter) != 1L || is.na(parameter)) {
    return(FALSE)
  }

  return(grepl("^(gamma|theta|pi|phi)\\[", parameter))
}


.iwmde_row_state <- function(context, row_index, parameter = NULL,
                             parameter_spec = NULL) {

  base_state      <- .iwmde_base_row_state(context, row_index)
  likelihood_mode <- .iwmde_likelihood_mode(parameter, parameter_spec)
  is_primitive    <- is.null(parameter_spec) ||
    identical(parameter_spec[["type"]], "primitive")
  prior_list <- .iwmde_active_flat_prior_list(
    context   = context,
    row       = base_state[["row"]],
    parameter = parameter
  )
  log_prior <- .iwmde_log_prior_row(base_state[["row"]], prior_list)
  if (is_primitive) {
    focal_prior <- .iwmde_focal_prior(context, parameter, base_state[["row"]])
    baseline_focal_log_prior <- .iwmde_focal_log_prior(
      prior     = focal_prior,
      value     = base_state[["row"]][[parameter]],
      parameter = parameter
    )
  } else {
    focal_prior              <- NULL
    baseline_focal_log_prior <- NA_real_
  }
  baseline_log_lik <- .iwmde_baseline_log_likelihood(
    context         = context,
    base_state      = base_state,
    likelihood_mode = likelihood_mode
  )
  baseline_log_q <- if (is.finite(baseline_log_lik) &&
                        is.finite(log_prior)) {
    baseline_log_lik + log_prior
  } else {
    -Inf
  }

  base_state[["prior_list"]]               <- prior_list
  base_state[["baseline_log_lik"]]         <- baseline_log_lik
  base_state[["baseline_log_prior"]]       <- log_prior
  base_state[["focal_prior"]]              <- focal_prior
  base_state[["baseline_focal_log_prior"]] <- baseline_focal_log_prior
  base_state[["use_focal_prior_delta"]]        <- is_primitive &&
    !.iwmde_parameter_is_eta(parameter) &&
    is.finite(baseline_focal_log_prior) &&
    .iwmde_can_use_focal_prior_delta(focal_prior)
  base_state[["baseline_log_q"]]               <- baseline_log_q
  base_state[["likelihood_mode"]]              <- likelihood_mode

  return(base_state)
}


.iwmde_baseline_log_likelihood <- function(context, base_state,
                                           likelihood_mode) {

  key <- paste(
    base_state[["row_index"]],
    .iwmde_active_key(context, base_state[["row"]]),
    likelihood_mode,
    sep = "|"
  )
  if (exists(key, envir = context[["likelihood_cache"]], inherits = FALSE)) {
    return(get(key, envir = context[["likelihood_cache"]]))
  }

  log_lik <- .iwmde_log_likelihood_parameters(
    context         = context,
    parameters      = base_state[["parameters"]],
    active_setup    = base_state[["active_setup"]],
    likelihood_mode = likelihood_mode,
    row             = base_state[["row"]]
  )

  assign(key, log_lik, envir = context[["likelihood_cache"]])
  return(log_lik)
}


.iwmde_base_row_state <- function(context, row_index) {

  key <- as.character(row_index)
  if (exists(key, envir = context[["row_cache"]], inherits = FALSE)) {
    return(get(key, envir = context[["row_cache"]]))
  }

  row          <- context[["posterior_samples"]][row_index, ]
  active_key   <- .iwmde_active_key(context, row)
  active_setup <- .iwmde_active_setup(context, row, active_key)
  parameters   <- .iwmde_row_parameters(context, row, active_setup)

  state <- list(
    row_index        = row_index,
    row              = row,
    active_key       = active_key,
    active_setup     = active_setup,
    parameters       = parameters
  )

  assign(key, state, envir = context[["row_cache"]])
  return(state)
}


.iwmde_state_active_key <- function(context, state) {

  if (!is.null(state[["active_key"]])) {
    return(state[["active_key"]])
  }

  return(.iwmde_active_key(context, state[["row"]]))
}


.iwmde_density_grid <- function(context, parameter, display_grid,
                                normalization_grid, transform, row_states,
                                active_mass, replacement) {

  n_display        <- length(display_grid)
  y                <- numeric(n_display)
  finite_terms     <- integer(n_display)
  max_log_ratio    <- numeric(n_display)
  ess              <- numeric(n_display)
  max_weight_share <- numeric(n_display)
  q_grid           <- c(display_grid, normalization_grid[["x"]])
  log_q_grid       <- .iwmde_log_q_grid(
    context     = context,
    parameter   = parameter,
    values      = q_grid,
    row_states  = row_states,
    replacement = replacement
  )
  norm_rows        <- length(display_grid) + seq_along(normalization_grid[["x"]])
  log_q_display    <- log_q_grid[seq_along(display_grid), , drop = FALSE]
  log_q_norm       <- log_q_grid[norm_rows, , drop = FALSE]

  log_normalizer <- .iwmde_log_trapz_columns(
    x     = normalization_grid[["z"]],
    log_y = log_q_norm + normalization_grid[["log_jacobian"]]
  )
  keep_rows      <- is.finite(log_normalizer)
  log_q_display  <- log_q_display[, keep_rows, drop = FALSE]
  log_q_norm     <- log_q_norm[, keep_rows, drop = FALSE]
  log_normalizer <- log_normalizer[keep_rows]

  if (length(log_normalizer) == 0L) {
    return(list(
      x                      = display_grid,
      y                      = y,
      finite_terms           = rep(0L, n_display),
      max_log_ratio          = rep(Inf, n_display),
      ess                    = rep(0, n_display),
      max_weight_share       = rep(1, n_display),
      mcse                   = rep(NA_real_, n_display),
      relative_mcse          = rep(NA_real_, n_display),
      log_normalizer         = log_normalizer,
      n_normalized_rows      = 0L,
      normalization_points   = length(normalization_grid[["x"]]),
      normalization_range    = range(normalization_grid[["x"]]),
      normalization_integral = NA_real_,
      normalization_scale    = transform[["type"]],
      integral_mcse          = NA_real_,
      integral_relative_mcse = NA_real_,
      batch_size             = NA_integer_,
      n_batches              = 0L,
      estimator              = "q_grid_cmde",
      weight_method          = "conditional_grid"
    ))
  }

  density_terms <- .iwmde_density_aggregate(
    log_terms   = sweep(log_q_display, 2L, log_normalizer, "-"),
    active_mass = active_mass,
    denominator = length(log_normalizer)
  )
  y                <- density_terms[["y"]]
  finite_terms     <- density_terms[["finite_terms"]]
  max_log_ratio    <- density_terms[["max_log_ratio"]]
  ess              <- density_terms[["ess"]]
  max_weight_share <- density_terms[["max_weight_share"]]
  contributions    <- density_terms[["contributions"]]

  mcse_data      <- .iwmde_batch_mcse(contributions)
  integral_mcse  <- .iwmde_integral_mcse(contributions, display_grid)
  norm_y         <- .iwmde_normalization_density(
    log_q_norm         = log_q_norm,
    log_normalizer     = log_normalizer,
    log_jacobian       = normalization_grid[["log_jacobian"]],
    normalization_grid = normalization_grid[["z"]],
    active_mass        = active_mass
  )

  return(list(
    x                      = display_grid,
    y                      = y,
    finite_terms           = finite_terms,
    max_log_ratio          = max_log_ratio,
    ess                    = ess,
    max_weight_share       = max_weight_share,
    mcse                   = mcse_data[["mcse"]],
    relative_mcse          = mcse_data[["relative_mcse"]],
    log_normalizer         = log_normalizer,
    n_normalized_rows      = length(log_normalizer),
    normalization_points   = length(normalization_grid[["x"]]),
    normalization_range    = range(normalization_grid[["x"]]),
    normalization_integral = .iwmde_trapz(normalization_grid[["z"]], norm_y),
    normalization_scale    = transform[["type"]],
    integral_mcse          = integral_mcse[["mcse"]],
    integral_relative_mcse = integral_mcse[["relative_mcse"]],
    batch_size             = mcse_data[["batch_size"]],
    n_batches              = mcse_data[["n_batches"]],
    estimator              = "q_grid_cmde",
    weight_method          = "conditional_grid"
  ))
}


.iwmde_density_iwmde <- function(context, parameter, parameter_spec,
                                 display_grid, row_states, active_rows,
                                 active_values, weight_rows, weight_values,
                                 support, active_mass, replacement) {

  n_display        <- length(display_grid)
  y                <- numeric(n_display)
  finite_terms     <- integer(n_display)
  max_log_ratio    <- numeric(n_display)
  ess              <- numeric(n_display)
  max_weight_share <- numeric(n_display)
  weight           <- .iwmde_chen_log_weight(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    active_rows    = active_rows,
    active_values  = active_values,
    weight_rows    = weight_rows,
    weight_values  = weight_values,
    support        = support
  )
  log_weight       <- weight[["log_weight"]]
  keep_rows <- is.finite(log_weight)
  row_states <- row_states[keep_rows]
  log_weight <- log_weight[keep_rows]

  if (length(log_weight) == 0L) {
    return(list(
      x                      = display_grid,
      y                      = y,
      finite_terms           = rep(0L, n_display),
      max_log_ratio          = rep(Inf, n_display),
      ess                    = rep(0, n_display),
      max_weight_share       = rep(1, n_display),
      mcse                   = rep(NA_real_, n_display),
      relative_mcse          = rep(NA_real_, n_display),
      log_normalizer         = numeric(),
      n_normalized_rows      = 0L,
      normalization_points   = 0L,
      normalization_range    = c(NA_real_, NA_real_),
      normalization_integral = NA_real_,
      normalization_scale    = NA_character_,
      integral_mcse          = NA_real_,
      integral_relative_mcse = NA_real_,
      batch_size             = NA_integer_,
      n_batches              = 0L,
      estimator              = "iwmde",
      weight_method          = weight[["method"]]
    ))
  }

  log_q_display <- .iwmde_log_q_grid(
    context     = context,
    parameter   = parameter,
    values      = display_grid,
    row_states  = row_states,
    replacement = replacement
  )
  baseline_log_q <- vapply(row_states, function(state) {
    state[["baseline_log_q"]]
  }, numeric(1))

  density_terms <- .iwmde_density_aggregate(
    log_terms = sweep(
      sweep(log_q_display, 2L, baseline_log_q, "-"),
      2L,
      log_weight,
      "+"
    ),
    active_mass = active_mass,
    denominator = length(log_weight)
  )
  y                <- density_terms[["y"]]
  finite_terms     <- density_terms[["finite_terms"]]
  max_log_ratio    <- density_terms[["max_log_ratio"]]
  ess              <- density_terms[["ess"]]
  max_weight_share <- density_terms[["max_weight_share"]]
  contributions    <- density_terms[["contributions"]]

  mcse_data     <- .iwmde_batch_mcse(contributions)
  integral_mcse <- .iwmde_integral_mcse(contributions, display_grid)

  return(list(
    x                      = display_grid,
    y                      = y,
    finite_terms           = finite_terms,
    max_log_ratio          = max_log_ratio,
    ess                    = ess,
    max_weight_share       = max_weight_share,
    mcse                   = mcse_data[["mcse"]],
    relative_mcse          = mcse_data[["relative_mcse"]],
    log_normalizer         = numeric(),
    n_normalized_rows      = length(log_weight),
    normalization_points   = 0L,
    normalization_range    = c(NA_real_, NA_real_),
    normalization_integral = NA_real_,
    normalization_scale    = NA_character_,
    integral_mcse          = integral_mcse[["mcse"]],
    integral_relative_mcse = integral_mcse[["relative_mcse"]],
    batch_size             = mcse_data[["batch_size"]],
    n_batches              = mcse_data[["n_batches"]],
    estimator              = "iwmde",
    weight_method          = weight[["method"]]
  ))
}


.iwmde_density_aggregate <- function(log_terms, active_mass, denominator) {

  if (!is.matrix(log_terms)) {
    log_terms <- as.matrix(log_terms)
  }

  n_grid           <- nrow(log_terms)
  y                <- numeric(n_grid)
  max_log_ratio    <- rep(Inf, n_grid)
  ess              <- numeric(n_grid)
  max_weight_share <- rep(1, n_grid)
  finite           <- is.finite(log_terms)
  finite_terms     <- rowSums(finite)
  contributions    <- matrix(0, nrow = n_grid, ncol = ncol(log_terms))

  active_rows <- finite_terms > 0L
  if (!any(active_rows)) {
    return(list(
      y                = y,
      finite_terms     = finite_terms,
      max_log_ratio    = max_log_ratio,
      ess              = ess,
      max_weight_share = max_weight_share,
      contributions    = contributions
    ))
  }

  safe_terms <- log_terms
  safe_terms[!finite] <- -Inf
  max_term <- safe_terms[cbind(
    seq_len(n_grid),
    max.col(safe_terms, ties.method = "first")
  )]

  scaled_terms <- exp(safe_terms - max_term)
  scaled_terms[!finite] <- 0
  scaled_terms[!active_rows, ] <- 0

  sum_scaled_terms <- rowSums(scaled_terms)
  sum_scaled_sq    <- rowSums(scaled_terms^2)
  max_scaled       <- scaled_terms[cbind(
    seq_len(n_grid),
    max.col(scaled_terms, ties.method = "first")
  )]

  y[active_rows] <- active_mass * exp(max_term[active_rows]) *
    sum_scaled_terms[active_rows] / denominator
  ess[active_rows] <- sum_scaled_terms[active_rows]^2 /
    sum_scaled_sq[active_rows]
  max_weight_share[active_rows] <- max_scaled[active_rows] /
    sum_scaled_terms[active_rows]
  contributions[finite] <- active_mass * exp(log_terms[finite])

  median_terms <- vapply(which(active_rows), function(row) {
    stats::median(log_terms[row, finite[row, ]])
  }, numeric(1))
  max_log_ratio[active_rows] <- max_term[active_rows] - median_terms

  return(list(
    y                = y,
    finite_terms     = finite_terms,
    max_log_ratio    = max_log_ratio,
    ess              = ess,
    max_weight_share = max_weight_share,
    contributions    = contributions
  ))
}


.iwmde_chen_log_weight <- function(context, parameter, parameter_spec,
                                   active_rows, active_values, weight_rows,
                                   weight_values, support) {

  if (all(is.finite(support))) {
    return(.iwmde_chen_power_log_weight(
      active_values = active_values,
      weight_values = weight_values,
      support       = support
    ))
  }
  if (is.finite(support[1])) {
    return(.iwmde_chen_exponential_log_weight(
      active_values = active_values,
      weight_values = weight_values,
      support       = support,
      side          = "lower"
    ))
  }
  if (is.finite(support[2])) {
    return(.iwmde_chen_exponential_log_weight(
      active_values = active_values,
      weight_values = weight_values,
      support       = support,
      side          = "upper"
    ))
  }

  return(.iwmde_chen_conditional_normal_log_weight(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    active_rows    = active_rows,
    active_values  = active_values,
    weight_rows    = weight_rows,
    weight_values  = weight_values
  ))
}


.iwmde_chen_power_log_weight <- function(active_values, weight_values,
                                         support) {

  out    <- rep(-Inf, length(active_values))
  lower  <- support[1]
  upper  <- support[2]
  width  <- upper - lower
  values <- weight_values[
    is.finite(weight_values) &
      weight_values > lower &
      weight_values < upper
  ]
  if (!is.finite(width) || width <= 0 || length(values) < 2L) {
    return(list(log_weight = out, method = "chen_power"))
  }

  location <- mean(values)
  if (!is.finite(location) || location <= lower || location >= upper) {
    return(list(log_weight = out, method = "chen_power"))
  }

  midpoint <- (lower + upper) / 2
  finite   <- is.finite(active_values) &
    active_values > lower &
    active_values < upper

  if (location >= midpoint) {
    alpha <- (location - lower) / (upper - location)
    if (!is.finite(alpha) || alpha <= 0) {
      return(list(log_weight = out, method = "chen_power"))
    }
    out[finite] <- log(alpha) +
      (alpha - 1) * log(active_values[finite] - lower) -
      alpha * log(width)
  } else {
    alpha <- (upper - location) / (location - lower)
    if (!is.finite(alpha) || alpha <= 0) {
      return(list(log_weight = out, method = "chen_power"))
    }
    out[finite] <- log(alpha) +
      (alpha - 1) * log(upper - active_values[finite]) -
      alpha * log(width)
  }

  return(list(log_weight = out, method = "chen_power"))
}


.iwmde_chen_exponential_log_weight <- function(active_values, weight_values,
                                               support, side) {

  out    <- rep(-Inf, length(active_values))
  values <- weight_values[is.finite(weight_values)]

  if (identical(side, "lower")) {
    lower  <- support[1]
    values <- values[values > lower]
    if (length(values) < 2L) {
      return(list(log_weight = out, method = "chen_exponential"))
    }
    scale <- mean(values) - lower
    keep  <- is.finite(active_values) & active_values > lower
    if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) {
      return(list(log_weight = out, method = "chen_exponential"))
    }
    rate      <- 1 / scale
    out[keep] <- log(rate) - rate * (active_values[keep] - lower)
  } else {
    upper  <- support[2]
    values <- values[values < upper]
    if (length(values) < 2L) {
      return(list(log_weight = out, method = "chen_exponential"))
    }
    scale <- upper - mean(values)
    keep  <- is.finite(active_values) & active_values < upper
    if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) {
      return(list(log_weight = out, method = "chen_exponential"))
    }
    rate      <- 1 / scale
    out[keep] <- log(rate) - rate * (upper - active_values[keep])
  }

  return(list(log_weight = out, method = "chen_exponential"))
}


.iwmde_chen_conditional_normal_log_weight <- function(context, parameter,
                                                      parameter_spec,
                                                      active_rows,
                                                      active_values,
                                                      weight_rows,
                                                      weight_values) {

  samples <- context[["posterior_samples"]]
  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    active_rows    = active_rows,
    weight_rows    = weight_rows
  )
  if (ncol(conditioning[["fit"]]) == 0L) {
    return(.iwmde_chen_marginal_normal_log_weight(
      active_values = active_values,
      weight_values = weight_values
    ))
  }

  x_fit  <- conditioning[["fit"]]
  x_eval <- conditioning[["eval"]]
  z_fit  <- weight_values

  fit_keep <- is.finite(z_fit) &
    stats::complete.cases(x_fit)
  if (sum(fit_keep) < 3L) {
    return(.iwmde_chen_marginal_normal_log_weight(
      active_values = active_values,
      weight_values = weight_values
    ))
  }

  x_fit <- x_fit[fit_keep, , drop = FALSE]
  z_fit <- z_fit[fit_keep]
  center <- colMeans(x_fit)
  scale  <- apply(x_fit, 2L, stats::sd)
  keep   <- is.finite(scale) & scale > sqrt(.Machine$double.eps)
  if (!any(keep)) {
    return(.iwmde_chen_marginal_normal_log_weight(
      active_values = active_values,
      weight_values = weight_values
    ))
  }

  x_fit  <- sweep(x_fit[, keep, drop = FALSE], 2L, center[keep], "-")
  x_fit  <- sweep(x_fit, 2L, scale[keep], "/")
  x_eval <- sweep(x_eval[, keep, drop = FALSE], 2L, center[keep], "-")
  x_eval <- sweep(x_eval, 2L, scale[keep], "/")

  cov_mat <- stats::cov(cbind(z_fit, x_fit))
  if (!all(is.finite(cov_mat))) {
    return(.iwmde_chen_marginal_normal_log_weight(
      active_values = active_values,
      weight_values = weight_values
    ))
  }

  p         <- ncol(x_fit)
  sigma_zz  <- cov_mat[1L, 1L]
  sigma_zx  <- matrix(cov_mat[1L, -1L], nrow = 1L)
  sigma_xx  <- cov_mat[-1L, -1L, drop = FALSE]
  ridge     <- max(1e-10, 1e-8 * mean(diag(sigma_xx)))
  inverse   <- try(
    chol2inv(chol(sigma_xx + diag(ridge, p))),
    silent = TRUE
  )
  if (inherits(inverse, "try-error")) {
    return(.iwmde_chen_marginal_normal_log_weight(
      active_values = active_values,
      weight_values = weight_values
    ))
  }

  beta     <- sigma_zx %*% inverse
  variance <- as.numeric(sigma_zz - beta %*% t(sigma_zx))
  if (!is.finite(variance) || variance <= sqrt(.Machine$double.eps)) {
    return(.iwmde_chen_marginal_normal_log_weight(
      active_values = active_values,
      weight_values = weight_values
    ))
  }

  out       <- rep(-Inf, length(active_values))
  eval_keep <- is.finite(active_values) & stats::complete.cases(x_eval)
  means     <- as.numeric(mean(z_fit) + x_eval[eval_keep, , drop = FALSE] %*%
    t(beta))
  out[eval_keep] <- stats::dnorm(
    active_values[eval_keep],
    mean = means,
    sd   = sqrt(variance),
    log  = TRUE
  )

  return(list(log_weight = out, method = "chen_conditional_normal"))
}


.iwmde_chen_marginal_normal_log_weight <- function(active_values,
                                                   weight_values) {

  out    <- rep(-Inf, length(active_values))
  values <- weight_values[is.finite(weight_values)]
  if (length(values) < 2L) {
    return(list(log_weight = out, method = "chen_marginal_normal"))
  }

  location <- mean(values)
  scale    <- stats::sd(values)
  finite   <- is.finite(active_values)
  if (!is.finite(location) || !is.finite(scale) ||
      scale <= sqrt(.Machine$double.eps)) {
    return(list(log_weight = out, method = "chen_marginal_normal"))
  }

  out[finite] <- stats::dnorm(
    active_values[finite],
    mean = location,
    sd   = scale,
    log  = TRUE
  )

  return(list(log_weight = out, method = "chen_marginal_normal"))
}


.iwmde_chen_conditioning_matrix <- function(context, parameter,
                                            parameter_spec, active_rows,
                                            weight_rows) {

  samples <- context[["posterior_samples"]]
  columns <- .iwmde_chen_conditioning_columns(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec
  )
  if (length(columns) == 0L) {
    return(list(
      fit     = matrix(numeric(), nrow = length(weight_rows), ncol = 0L),
      eval    = matrix(numeric(), nrow = length(active_rows), ncol = 0L),
      columns = character()
    ))
  }

  fit_values  <- matrix(NA_real_, nrow = length(weight_rows), ncol = length(columns))
  eval_values <- matrix(NA_real_, nrow = length(active_rows), ncol = length(columns))
  colnames(fit_values)  <- columns
  colnames(eval_values) <- columns

  for (column in columns) {
    transformed <- .iwmde_chen_transform_conditioning_column(
      context     = context,
      fit_values  = samples[weight_rows, column],
      eval_values = samples[active_rows, column],
      column      = column
    )
    fit_values[, column]  <- transformed[["fit"]]
    eval_values[, column] <- transformed[["eval"]]
  }

  keep <- vapply(seq_along(columns), function(i) {
    values <- fit_values[, i]
    finite <- is.finite(values)
    sum(finite) >= 3L &&
      stats::sd(values[finite]) > sqrt(.Machine$double.eps) &&
      length(unique(signif(values[finite], 12))) > 1L
  }, logical(1))

  if (!any(keep)) {
    return(list(
      fit     = matrix(numeric(), nrow = length(weight_rows), ncol = 0L),
      eval    = matrix(numeric(), nrow = length(active_rows), ncol = 0L),
      columns = character()
    ))
  }

  columns     <- columns[keep]
  fit_values  <- fit_values[, keep, drop = FALSE]
  eval_values <- eval_values[, keep, drop = FALSE]

  return(list(
    fit     = fit_values,
    eval    = eval_values,
    columns = columns
  ))
}


.iwmde_chen_conditioning_columns <- function(context, parameter,
                                             parameter_spec) {

  samples <- context[["posterior_samples"]]
  columns <- colnames(samples)
  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "linear")) {
    target_columns <- names(parameter_spec[["weights"]])
  } else {
    target_columns <- parameter
  }

  candidates <- columns[
    !columns %in% target_columns &
      !columns %in% context[["indicator_names"]] &
      !vapply(columns, .iwmde_parameter_is_indicator, logical(1)) &
      !vapply(columns, .iwmde_parameter_is_local_latent, logical(1)) &
      vapply(columns, function(column) {
        .iwmde_chen_is_global_conditioning_column(
          context = context,
          column  = column
        )
      }, logical(1))
  ]

  return(candidates)
}


.iwmde_chen_is_global_conditioning_column <- function(context, column) {

  return(
    column == "mu" ||
      startsWith(column, "mu_") ||
      grepl("^mu\\[[0-9]+\\]$", column) ||
      column %in% c("tau", "log_tau", "rho", "PET", "PEESE") ||
      startsWith(column, "tau_") ||
      startsWith(column, "log_tau_") ||
      .iwmde_parameter_is_weightfunction_coordinate(
        parameter       = column,
        context         = context,
        include_private = FALSE
      )
  )
}


.iwmde_chen_transform_conditioning_column <- function(context,
                                                      fit_values,
                                                      eval_values,
                                                      column) {

  if (column == "tau" || startsWith(column, "tau_")) {
    return(.iwmde_chen_transform_nonnegative(fit_values, eval_values))
  }
  if (column == "rho") {
    return(.iwmde_chen_transform_unit_interval(fit_values, eval_values))
  }
  omega_coordinates <- if (!is.null(context[["selection_spec"]])) {
    context[["selection_spec"]][["jags_omega"]]
  } else {
    character()
  }
  if (.iwmde_parameter_matches_coordinate(column, omega_coordinates)) {
    return(.iwmde_chen_transform_omega(fit_values, eval_values))
  }
  if (column == "log_tau" ||
      startsWith(column, "log_tau_")) {
    return(list(
      fit  = as.numeric(fit_values),
      eval = as.numeric(eval_values)
    ))
  }

  return(list(
    fit  = as.numeric(fit_values),
    eval = as.numeric(eval_values)
  ))
}


.iwmde_chen_transform_nonnegative <- function(fit_values, eval_values) {

  return(list(
    fit  = log1p(pmax(as.numeric(fit_values), 0)),
    eval = log1p(pmax(as.numeric(eval_values), 0))
  ))
}


.iwmde_chen_transform_unit_interval <- function(fit_values, eval_values) {

  eps <- sqrt(.Machine$double.eps)

  return(list(
    fit  = stats::qlogis(pmin(pmax(as.numeric(fit_values), eps), 1 - eps)),
    eval = stats::qlogis(pmin(pmax(as.numeric(eval_values), eps), 1 - eps))
  ))
}


.iwmde_chen_transform_omega <- function(fit_values, eval_values) {

  fit  <- as.numeric(fit_values)
  eval <- as.numeric(eval_values)
  all_values <- c(fit, eval)
  finite     <- all_values[is.finite(all_values)]
  if (length(finite) > 0L && min(finite) >= 0 &&
      max(finite) <= 1 + sqrt(.Machine$double.eps)) {
    return(.iwmde_chen_transform_unit_interval(fit, eval))
  }

  return(.iwmde_chen_transform_nonnegative(fit, eval))
}


.iwmde_batch_mcse <- function(contributions) {

  n <- ncol(contributions)
  if (n < 4L) {
    return(list(
      mcse          = rep(NA_real_, nrow(contributions)),
      relative_mcse = rep(NA_real_, nrow(contributions)),
      batch_size    = NA_integer_,
      n_batches     = 0L
    ))
  }

  batch_size <- max(2L, floor(sqrt(n)))
  n_batches  <- floor(n / batch_size)
  if (n_batches < 2L) {
    return(list(
      mcse          = rep(NA_real_, nrow(contributions)),
      relative_mcse = rep(NA_real_, nrow(contributions)),
      batch_size    = batch_size,
      n_batches     = n_batches
    ))
  }

  keep          <- seq_len(n_batches * batch_size)
  contributions <- contributions[, keep, drop = FALSE]
  batch_index   <- rep(seq_len(n_batches), each = batch_size)
  batch_means   <- matrix(NA_real_, nrow = nrow(contributions), ncol = n_batches)

  for (batch in seq_len(n_batches)) {
    batch_means[, batch] <- rowMeans(
      contributions[, batch_index == batch, drop = FALSE]
    )
  }

  mcse <- apply(batch_means, 1, stats::sd) / sqrt(n_batches)
  mean_contribution <- rowMeans(contributions)
  relative_mcse <- rep(NA_real_, length(mcse))
  positive <- mean_contribution > 0
  relative_mcse[positive] <- mcse[positive] / mean_contribution[positive]

  return(list(
    mcse          = mcse,
    relative_mcse = relative_mcse,
    batch_size    = batch_size,
    n_batches     = n_batches
  ))
}


.iwmde_integral_mcse <- function(contributions, x) {

  n <- ncol(contributions)
  if (n < 4L) {
    return(list(mcse = NA_real_, relative_mcse = NA_real_))
  }

  integrals <- .iwmde_trapz_columns(
    x = x,
    y = contributions
  )
  finite <- is.finite(integrals)
  if (sum(finite) < 4L) {
    return(list(mcse = NA_real_, relative_mcse = NA_real_))
  }

  integrals  <- integrals[finite]
  batch_size <- max(2L, floor(sqrt(length(integrals))))
  n_batches  <- floor(length(integrals) / batch_size)
  if (n_batches < 2L) {
    return(list(mcse = NA_real_, relative_mcse = NA_real_))
  }

  keep        <- seq_len(n_batches * batch_size)
  integrals   <- integrals[keep]
  batch_index <- rep(seq_len(n_batches), each = batch_size)
  batch_means <- vapply(seq_len(n_batches), function(batch) {
    mean(integrals[batch_index == batch])
  }, numeric(1))

  mcse          <- stats::sd(batch_means) / sqrt(n_batches)
  mean_integral <- mean(integrals)
  relative_mcse <- if (mean_integral > 0) mcse / mean_integral else NA_real_

  return(list(mcse = mcse, relative_mcse = relative_mcse))
}


.iwmde_normalization_density <- function(log_q_norm, log_normalizer,
                                         log_jacobian, normalization_grid,
                                         active_mass) {

  y <- numeric(length(normalization_grid))
  for (g in seq_along(normalization_grid)) {
    log_terms <- log_q_norm[g, ] + log_jacobian[g] - log_normalizer
    finite    <- is.finite(log_terms)
    if (any(finite)) {
      max_term         <- max(log_terms[finite])
      scaled_terms     <- exp(log_terms[finite] - max_term)
      y[g]             <- active_mass * exp(max_term) *
        sum(scaled_terms) / length(log_normalizer)
    }
  }

  return(y)
}


.iwmde_log_trapz <- function(x, log_y) {

  finite <- is.finite(log_y)
  if (sum(finite) < 2L) {
    return(-Inf)
  }

  max_log <- max(log_y[finite])
  y       <- exp(log_y - max_log)
  area    <- .iwmde_trapz(x, y)
  if (!is.finite(area) || area <= 0) {
    return(-Inf)
  }

  return(log(area) + max_log)
}


.iwmde_log_trapz_columns <- function(x, log_y) {

  log_y <- as.matrix(log_y)
  out   <- rep(-Inf, ncol(log_y))
  keep  <- colSums(is.finite(log_y)) >= 2L
  if (!any(keep)) {
    return(out)
  }

  log_y_keep <- log_y[, keep, drop = FALSE]
  max_log    <- apply(log_y_keep, 2L, max, na.rm = TRUE)
  y          <- exp(sweep(log_y_keep, 2L, max_log, "-"))
  y[!is.finite(y)] <- 0
  area       <- .iwmde_trapz_columns(x = x, y = y)
  valid      <- is.finite(area) & area > 0
  out[which(keep)[valid]] <- log(area[valid]) + max_log[valid]

  return(out)
}
