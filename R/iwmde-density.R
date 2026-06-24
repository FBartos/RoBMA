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


.iwmde_state_scope <- function(parameter, parameter_spec = NULL) {

  if (identical(.iwmde_likelihood_mode(parameter, parameter_spec), "conditional")) {
    return("local")
  }

  return("global")
}


.iwmde_state_scope_value <- function(state) {

  state_scope <- state[["state_scope"]]
  if (length(state_scope) != 1L ||
      !state_scope %in% c("local", "global")) {
    return("local")
  }

  return(state_scope)
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


.iwmde_prior_name_is_local_latent <- function(parameter) {

  if (is.null(parameter) || length(parameter) != 1L || is.na(parameter)) {
    return(FALSE)
  }

  return(parameter %in% c("gamma", "theta", "pi", "phi"))
}


.iwmde_drop_local_latent_sample_columns <- function(samples) {

  keep <- !vapply(
    colnames(samples),
    .iwmde_parameter_is_local_latent,
    logical(1)
  )

  return(samples[, keep, drop = FALSE])
}


.iwmde_row_state <- function(context, row_index, parameter = NULL,
                             parameter_spec = NULL) {

  likelihood_mode <- .iwmde_likelihood_mode(parameter, parameter_spec)
  state_scope     <- .iwmde_state_scope(parameter, parameter_spec)
  base_state      <- .iwmde_base_row_state(
    context     = context,
    row_index   = row_index,
    state_scope = state_scope
  )
  is_primitive    <- is.null(parameter_spec) ||
    identical(parameter_spec[["type"]], "primitive")
  prior_list <- .iwmde_active_flat_prior_list(
    context     = context,
    row         = base_state[["row"]],
    parameter   = parameter,
    state_scope = state_scope
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
  base_state[["state_scope"]]                  <- state_scope

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

  row <- if (identical(base_state[["state_scope"]], "global") &&
             identical(likelihood_mode, "marginal")) {
    NULL
  } else {
    base_state[["row"]]
  }

  log_lik <- .iwmde_log_likelihood_parameters(
    context         = context,
    parameters      = base_state[["parameters"]],
    active_setup    = base_state[["active_setup"]],
    likelihood_mode = likelihood_mode,
    row             = row
  )

  assign(key, log_lik, envir = context[["likelihood_cache"]])
  return(log_lik)
}


.iwmde_base_row_state <- function(context, row_index,
                                  state_scope = c("local", "global")) {

  state_scope <- match.arg(state_scope)
  key <- paste(row_index, state_scope, sep = "|")
  if (exists(key, envir = context[["row_cache"]], inherits = FALSE)) {
    return(get(key, envir = context[["row_cache"]]))
  }

  row          <- context[["posterior_samples"]][row_index, ]
  active_key   <- .iwmde_active_key(context, row)
  active_setup <- .iwmde_active_setup(context, row, active_key)
  parameters   <- .iwmde_row_parameters(
    context      = context,
    row          = row,
    active_setup = active_setup,
    state_scope  = state_scope
  )

  state <- list(
    row_index        = row_index,
    row              = row,
    active_key       = active_key,
    active_setup     = active_setup,
    parameters       = parameters,
    state_scope      = state_scope
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
                                active_mass, replacement,
                                n_candidate_rows = length(row_states)) {

  n_display        <- length(display_grid)
  y                <- numeric(n_display)
  finite_terms     <- integer(n_display)
  max_log_ratio    <- numeric(n_display)
  ess              <- numeric(n_display)
  max_weight_share <- numeric(n_display)
  n_input_rows     <- length(row_states)
  n_candidate_rows <- as.integer(n_candidate_rows[[1L]])
  if (!is.finite(n_candidate_rows) || n_candidate_rows < n_input_rows) {
    n_candidate_rows <- n_input_rows
  }
  normalizer_plan  <- .iwmde_qcmde_normalizer_plan(
    normalization_grid = normalization_grid,
    transform          = transform
  )
  final_grid       <- normalizer_plan[["final_grid"]]
  pilot_grid       <- normalizer_plan[["pilot_grid"]]
  all_grid         <- normalizer_plan[["all_grid"]]
  q_grid           <- c(display_grid, all_grid[["x"]])
  log_q_grid       <- .iwmde_log_q_grid(
    context     = context,
    parameter   = parameter,
    values      = q_grid,
    row_states  = row_states,
    replacement = replacement
  )
  quadrature_change <- attr(
    log_q_grid,
    "max_quadrature_relative_change",
    exact = TRUE
  )
  if (is.null(quadrature_change)) {
    quadrature_change <- NA_real_
  }
  all_grid_rows    <- length(display_grid) + seq_along(all_grid[["x"]])
  log_q_display    <- log_q_grid[seq_along(display_grid), , drop = FALSE]
  log_q_all        <- log_q_grid[all_grid_rows, , drop = FALSE]
  log_q_sequence   <- lapply(normalizer_plan[["grid_sequence"]], function(grid) {
    log_q_all[grid[["all_index"]], , drop = FALSE]
  })
  log_normalizer_sequence <- lapply(seq_along(log_q_sequence), function(i) {
    grid <- normalizer_plan[["grid_sequence"]][[i]]
    .iwmde_log_trapz_columns(
      x     = grid[["z"]],
      log_y = log_q_sequence[[i]] + grid[["log_jacobian"]]
    )
  })
  refinement <- .iwmde_qcmde_select_refinement(
    log_q_display           = log_q_display,
    log_normalizer_sequence = log_normalizer_sequence,
    active_mass             = active_mass,
    denominator             = n_candidate_rows
  )
  normalizer_plan[["pilot_index"]] <- refinement[["pilot_index"]]
  normalizer_plan[["final_index"]] <- refinement[["final_index"]]
  normalizer_plan[["validation_index"]] <- refinement[["validation_index"]]
  normalizer_plan[["n_refinement_steps"]] <- refinement[["n_refinement_steps"]]
  pilot_grid      <- normalizer_plan[["grid_sequence"]][[refinement[["pilot_index"]]]]
  final_grid      <- normalizer_plan[["grid_sequence"]][[refinement[["final_index"]]]]
  validation_grid <- normalizer_plan[["grid_sequence"]][[refinement[["validation_index"]]]]
  log_q_initial   <- log_q_sequence[[refinement[["pilot_index"]]]]
  log_q_final     <- log_q_sequence[[refinement[["final_index"]]]]
  log_q_validation <- log_q_sequence[[refinement[["validation_index"]]]]

  initial_log_normalizer <- log_normalizer_sequence[[refinement[["pilot_index"]]]]
  final_log_normalizer   <- log_normalizer_sequence[[refinement[["final_index"]]]]
  validation_log_normalizer <- log_normalizer_sequence[[refinement[["validation_index"]]]]
  initial_finite <- is.finite(initial_log_normalizer)
  final_finite   <- is.finite(final_log_normalizer)
  validation_finite <- is.finite(validation_log_normalizer)
  pilot_y        <- .iwmde_qcmde_pilot_density(
    log_q_display  = log_q_display,
    log_normalizer = initial_log_normalizer,
    keep_rows      = initial_finite,
    active_mass    = active_mass,
    denominator    = n_candidate_rows
  )
  validation_y <- .iwmde_qcmde_pilot_density(
    log_q_display  = log_q_display,
    log_normalizer = validation_log_normalizer,
    keep_rows      = validation_finite,
    active_mass    = active_mass,
    denominator    = n_candidate_rows
  )
  keep_rows      <- final_finite
  log_q_display  <- log_q_display[, keep_rows, drop = FALSE]
  log_q_initial  <- log_q_initial[, keep_rows, drop = FALSE]
  log_q_final    <- log_q_final[, keep_rows, drop = FALSE]
  log_normalizer <- final_log_normalizer[keep_rows]
  row_states     <- row_states[keep_rows]
  normalizer_change <- .iwmde_qcmde_normalizer_change(
    initial_log_normalizer = final_log_normalizer,
    final_log_normalizer   = validation_log_normalizer
  )
  n_normalized_rows <- length(log_normalizer)
  n_dropped_rows    <- max(0L, n_candidate_rows - n_normalized_rows)

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
      pilot_log_normalizer   = initial_log_normalizer[keep_rows],
      n_normalized_rows      = 0L,
      n_candidate_rows       = n_candidate_rows,
      n_evaluated_rows       = n_input_rows,
      n_dropped_rows         = n_dropped_rows,
      row_drop_fraction      = .iwmde_row_drop_fraction(
        n_candidate_rows  = n_candidate_rows,
        n_normalized_rows = 0L
      ),
      normalization_points              = length(final_grid[["x"]]),
      normalization_range               = range(final_grid[["x"]]),
      normalization_initial_points      = length(pilot_grid[["x"]]),
      normalization_initial_range       = range(pilot_grid[["x"]]),
      normalization_integral            = NA_real_,
      normalization_final_integral      = NA_real_,
      normalization_scale               = transform[["type"]],
      normalization_mass_ratio          = NA_real_,
      pilot_y                           = pilot_y,
      validation_y                      = validation_y,
      ordinate_relative_change          = rep(NA_real_, n_display),
      ordinate_log_change               = rep(NA_real_, n_display),
      max_normalizer_relative_change    = normalizer_change[["max"]],
      max_quadrature_relative_change    = quadrature_change,
      p95_normalizer_relative_change    = normalizer_change[["p95"]],
      median_normalizer_relative_change = normalizer_change[["median"]],
      normalization_refined_points      = length(final_grid[["x"]]),
      normalization_refined_range       = range(final_grid[["x"]]),
      n_rescued_normalizer              = sum(!initial_finite & final_finite),
      n_dropped_normalizer              = sum(!final_finite),
      n_initial_dropped_normalizer      = sum(!initial_finite),
      n_validation_dropped_normalizer   = sum(final_finite & !validation_finite),
      n_refinement_steps                = normalizer_plan[["n_refinement_steps"]],
      integral_mcse                     = NA_real_,
      integral_relative_mcse            = NA_real_,
      batch_size                        = NA_integer_,
      n_batches                         = 0L,
      estimator                         = "q_grid_cmde",
      weight_method                     = "conditional_grid"
    ))
  }

  density_terms <- .iwmde_density_aggregate(
    log_terms   = sweep(log_q_display, 2L, log_normalizer, "-"),
    active_mass = active_mass,
    denominator = n_candidate_rows
  )
  y                <- density_terms[["y"]]
  finite_terms     <- density_terms[["finite_terms"]]
  max_log_ratio    <- density_terms[["max_log_ratio"]]
  ess              <- density_terms[["ess"]]
  max_weight_share <- density_terms[["max_weight_share"]]
  contributions    <- density_terms[["contributions"]]

  mcse_data      <- .iwmde_batch_mcse(contributions)
  integral_mcse  <- .iwmde_integral_mcse(contributions, display_grid)
  norm_y_initial <- .iwmde_normalization_density(
    log_q_norm         = log_q_initial,
    log_normalizer     = log_normalizer,
    log_jacobian       = pilot_grid[["log_jacobian"]],
    normalization_grid = pilot_grid[["z"]],
    active_mass        = active_mass,
    denominator        = n_candidate_rows
  )
  norm_y_final <- .iwmde_normalization_density(
    log_q_norm         = log_q_final,
    log_normalizer     = log_normalizer,
    log_jacobian       = final_grid[["log_jacobian"]],
    normalization_grid = final_grid[["z"]],
    active_mass        = active_mass,
    denominator        = n_candidate_rows
  )
  ordinate_change <- .iwmde_qcmde_ordinate_change(
    pilot_y = y,
    final_y = validation_y
  )
  pilot_change <- .iwmde_qcmde_ordinate_change(
    pilot_y = pilot_y,
    final_y = y
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
    pilot_log_normalizer   = initial_log_normalizer[keep_rows],
    n_normalized_rows      = n_normalized_rows,
    n_candidate_rows       = n_candidate_rows,
    n_evaluated_rows       = n_input_rows,
    n_dropped_rows         = n_dropped_rows,
    row_drop_fraction      = .iwmde_row_drop_fraction(
      n_candidate_rows  = n_candidate_rows,
      n_normalized_rows = n_normalized_rows
    ),
    normalization_points              = length(final_grid[["x"]]),
    normalization_range               = range(final_grid[["x"]]),
    normalization_initial_points      = length(pilot_grid[["x"]]),
    normalization_initial_range       = range(pilot_grid[["x"]]),
    normalization_integral            = .iwmde_trapz(pilot_grid[["z"]], norm_y_initial),
    normalization_final_integral      = .iwmde_trapz(final_grid[["z"]], norm_y_final),
    normalization_scale               = transform[["type"]],
    normalization_mass_ratio          = 1,
    pilot_y                           = pilot_y,
    validation_y                      = validation_y,
    ordinate_relative_change          = ordinate_change[["relative"]],
    ordinate_log_change               = ordinate_change[["log"]],
    pilot_ordinate_relative_change    = pilot_change[["relative"]],
    pilot_ordinate_log_change         = pilot_change[["log"]],
    max_normalizer_relative_change    = normalizer_change[["max"]],
    max_quadrature_relative_change    = quadrature_change,
    p95_normalizer_relative_change    = normalizer_change[["p95"]],
    median_normalizer_relative_change = normalizer_change[["median"]],
    normalization_refined_points      = length(final_grid[["x"]]),
    normalization_refined_range       = range(final_grid[["x"]]),
    n_rescued_normalizer              = sum(!initial_finite & final_finite),
    n_dropped_normalizer              = sum(!final_finite),
    n_initial_dropped_normalizer      = sum(!initial_finite),
    n_validation_dropped_normalizer   = sum(final_finite & !validation_finite),
    n_refinement_steps                = normalizer_plan[["n_refinement_steps"]],
    integral_mcse                     = integral_mcse[["mcse"]],
    integral_relative_mcse            = integral_mcse[["relative_mcse"]],
    batch_size                        = mcse_data[["batch_size"]],
    n_batches                         = mcse_data[["n_batches"]],
    estimator                         = "q_grid_cmde",
    weight_method                     = "conditional_grid"
  ))
}


.iwmde_density_iwmde <- function(context, parameter, parameter_spec,
                                 display_grid, row_states, active_rows,
                                 active_values, weight_rows, weight_values,
                                 support, active_mass, replacement,
                                 normalization_grid = NULL,
                                 n_candidate_rows = length(row_states)) {

  n_display        <- length(display_grid)
  y                <- numeric(n_display)
  finite_terms     <- integer(n_display)
  max_log_ratio    <- numeric(n_display)
  ess              <- numeric(n_display)
  max_weight_share <- numeric(n_display)
  n_input_rows     <- length(row_states)
  n_candidate_rows <- as.integer(n_candidate_rows[[1L]])
  if (!is.finite(n_candidate_rows) || n_candidate_rows < n_input_rows) {
    n_candidate_rows <- n_input_rows
  }
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
  raw_log_weight    <- weight[["log_weight"]]
  keep_rows         <- is.finite(raw_log_weight)
  n_dropped_weight  <- sum(!keep_rows)
  row_states        <- row_states[keep_rows]
  log_weight        <- raw_log_weight[keep_rows]
  n_normalized_rows <- length(log_weight)
  n_dropped_rows    <- max(0L, n_candidate_rows - n_normalized_rows)

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
      n_candidate_rows       = n_candidate_rows,
      n_evaluated_rows       = n_input_rows,
      n_dropped_rows         = n_dropped_rows,
      row_drop_fraction      = .iwmde_row_drop_fraction(
        n_candidate_rows  = n_candidate_rows,
        n_normalized_rows = 0L
      ),
      normalization_points              = 0L,
      normalization_range               = c(NA_real_, NA_real_),
      normalization_integral            = NA_real_,
      normalization_scale               = NA_character_,
      normalization_mass_ratio          = NA_real_,
      max_normalizer_relative_change    = NA_real_,
      max_quadrature_relative_change    = NA_real_,
      median_normalizer_relative_change = NA_real_,
      normalization_refined_points      = 0L,
      normalization_refined_range       = c(NA_real_, NA_real_),
      n_dropped_weight                  = n_dropped_weight,
      weight_partitions                 = weight[["partitions"]],
      integral_mcse                     = NA_real_,
      integral_relative_mcse            = NA_real_,
      batch_size                        = NA_integer_,
      n_batches                         = 0L,
      estimator                         = "iwmde",
      weight_method                     = weight[["method"]]
    ))
  }

  has_normalization_grid <- !is.null(normalization_grid) &&
    length(normalization_grid[["x"]]) >= 2L
  q_grid <- if (has_normalization_grid) {
    c(display_grid, normalization_grid[["x"]])
  } else {
    display_grid
  }
  log_q_grid <- .iwmde_log_q_grid(
    context     = context,
    parameter   = parameter,
    values      = q_grid,
    row_states  = row_states,
    replacement = replacement
  )
  quadrature_change <- attr(
    log_q_grid,
    "max_quadrature_relative_change",
    exact = TRUE
  )
  if (is.null(quadrature_change)) {
    quadrature_change <- NA_real_
  }
  log_q_display <- log_q_grid[seq_len(n_display), , drop = FALSE]
  log_q_norm <- if (has_normalization_grid) {
    log_q_grid[n_display + seq_along(normalization_grid[["x"]]), , drop = FALSE]
  } else {
    NULL
  }
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
    denominator = n_candidate_rows
  )
  y                <- density_terms[["y"]]
  finite_terms     <- density_terms[["finite_terms"]]
  max_log_ratio    <- density_terms[["max_log_ratio"]]
  ess              <- density_terms[["ess"]]
  max_weight_share <- density_terms[["max_weight_share"]]
  contributions    <- density_terms[["contributions"]]

  normalization <- .iwmde_iwmde_normalization(
    normalization_grid = normalization_grid,
    log_q_norm         = log_q_norm,
    baseline_log_q     = baseline_log_q,
    log_weight         = log_weight,
    active_mass        = active_mass,
    denominator        = n_candidate_rows
  )

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
    n_normalized_rows      = n_normalized_rows,
    n_candidate_rows       = n_candidate_rows,
    n_evaluated_rows       = n_input_rows,
    n_dropped_rows         = n_dropped_rows,
    row_drop_fraction      = .iwmde_row_drop_fraction(
      n_candidate_rows  = n_candidate_rows,
      n_normalized_rows = n_normalized_rows
    ),
    normalization_points              = normalization[["points"]],
    normalization_range               = normalization[["range"]],
    normalization_integral            = normalization[["integral"]],
    normalization_scale               = normalization[["scale_type"]],
    normalization_mass_ratio          = normalization[["mass_ratio"]],
    max_normalizer_relative_change    = NA_real_,
    max_quadrature_relative_change    = quadrature_change,
    median_normalizer_relative_change = NA_real_,
    normalization_refined_points      = 0L,
    normalization_refined_range       = c(NA_real_, NA_real_),
    n_dropped_weight                  = n_dropped_weight,
    weight_partitions                 = weight[["partitions"]],
    integral_mcse                     = integral_mcse[["mcse"]],
    integral_relative_mcse            = integral_mcse[["relative_mcse"]],
    batch_size                        = mcse_data[["batch_size"]],
    n_batches                         = mcse_data[["n_batches"]],
    estimator                         = "iwmde",
    weight_method                     = weight[["method"]]
  ))
}


.iwmde_density_aggregate <- function(log_terms, active_mass, denominator) {

  if (!is.matrix(log_terms)) {
    log_terms <- as.matrix(log_terms)
  }
  denominator <- as.integer(denominator[[1L]])
  if (!is.finite(denominator) || denominator < ncol(log_terms)) {
    denominator <- ncol(log_terms)
  }

  n_grid           <- nrow(log_terms)
  y                <- numeric(n_grid)
  max_log_ratio    <- rep(Inf, n_grid)
  ess              <- numeric(n_grid)
  max_weight_share <- rep(1, n_grid)
  finite           <- is.finite(log_terms)
  finite_terms     <- rowSums(finite)
  contributions    <- matrix(0, nrow = n_grid, ncol = denominator)

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
  if (ncol(log_terms) > 0L) {
    contribution_terms          <- active_mass * exp(log_terms)
    contribution_terms[!finite] <- 0
    contributions[, seq_len(ncol(log_terms))] <- contribution_terms
  }

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


.iwmde_iwmde_normalization <- function(normalization_grid,
                                       log_q_norm,
                                       baseline_log_q,
                                       log_weight, active_mass,
                                       denominator) {

  empty <- list(
    mass_ratio = NA_real_,
    integral   = NA_real_,
    points     = 0L,
    range      = c(NA_real_, NA_real_),
    scale_type = NA_character_
  )
  if (is.null(normalization_grid) ||
      is.null(log_q_norm) ||
      length(normalization_grid[["x"]]) < 2L ||
      nrow(log_q_norm) != length(normalization_grid[["x"]]) ||
      ncol(log_q_norm) == 0L ||
      length(baseline_log_q) != ncol(log_q_norm) ||
      length(log_weight) != ncol(log_q_norm) ||
      length(log_weight) == 0L) {
    return(empty)
  }

  log_terms <- sweep(
    sweep(log_q_norm, 2L, baseline_log_q, "-"),
    2L,
    log_weight,
    "+"
  )
  log_terms <- sweep(
    log_terms,
    1L,
    normalization_grid[["log_jacobian"]],
    "+"
  )
  density_terms <- .iwmde_density_aggregate(
    log_terms   = log_terms,
    active_mass = active_mass,
    denominator = denominator
  )
  integral <- .iwmde_trapz(
    normalization_grid[["z"]],
    density_terms[["y"]]
  )
  mass_ratio <- if (is.finite(integral) && integral > 0) {
    active_mass / integral
  } else {
    NA_real_
  }

  return(list(
    mass_ratio = mass_ratio,
    integral   = integral,
    points     = length(normalization_grid[["x"]]),
    range      = range(normalization_grid[["x"]]),
    scale_type = "support_grid"
  ))
}


.iwmde_chen_log_weight <- function(context, parameter, parameter_spec,
                                   active_rows, active_values, weight_rows,
                                   weight_values, support) {

  active_supports <- .iwmde_parameter_row_supports(
    context        = context,
    parameter      = parameter,
    rows           = active_rows,
    parameter_spec = parameter_spec
  )
  weight_supports <- .iwmde_parameter_row_supports(
    context        = context,
    parameter      = parameter,
    rows           = weight_rows,
    parameter_spec = parameter_spec
  )
  active_keys  <- paste(
    .iwmde_chen_row_active_keys(context, active_rows),
    .iwmde_chen_support_keys(active_supports),
    sep = "||"
  )
  weight_keys  <- paste(
    .iwmde_chen_row_active_keys(context, weight_rows),
    .iwmde_chen_support_keys(weight_supports),
    sep = "||"
  )
  out          <- rep(-Inf, length(active_values))
  methods      <- character()
  partitions   <- list()

  for (key in unique(active_keys)) {
    active_index <- which(active_keys == key)
    weight_index <- which(weight_keys == key)
    if (length(active_index) == 0L || length(weight_index) == 0L) {
      next
    }

    row_support <- active_supports[active_index[1L], ]
    if (!.iwmde_chen_valid_support(row_support)) {
      next
    }
    if (!is.finite(row_support[1]) && !is.finite(row_support[2])) {
      weight <- .iwmde_chen_try_weight(
        expr = .iwmde_chen_conditional_normal_log_weight(
          context        = context,
          parameter      = parameter,
          parameter_spec = parameter_spec,
          active_rows    = active_rows[active_index],
          active_values  = active_values[active_index],
          weight_rows    = weight_rows[weight_index],
          weight_values  = weight_values[weight_index]
        ),
        fallback = .iwmde_chen_marginal_normal_log_weight(
          active_values = active_values[active_index],
          weight_values = weight_values[weight_index]
        )
      )
    } else if (all(is.finite(row_support))) {
      weight <- .iwmde_chen_try_weight(
        expr = .iwmde_chen_logit_conditional_normal_log_weight(
          context        = context,
          parameter      = parameter,
          parameter_spec = parameter_spec,
          active_rows    = active_rows[active_index],
          active_values  = active_values[active_index],
          weight_rows    = weight_rows[weight_index],
          weight_values  = weight_values[weight_index],
          support        = row_support
        ),
        fallback = .iwmde_chen_beta_log_weight(
          active_values = active_values[active_index],
          weight_values = weight_values[weight_index],
          support       = row_support
        )
      )
    } else {
      weight <- .iwmde_chen_gamma_log_weight_single(
        active_values = active_values[active_index],
        weight_values = weight_values[weight_index],
        support       = row_support
      )
    }

    out[active_index] <- weight[["log_weight"]]
    methods <- c(methods, weight[["method"]])
    partitions[[length(partitions) + 1L]] <- list(
      key            = key,
      support        = row_support,
      method         = weight[["method"]],
      n_eval_rows    = length(active_index),
      n_fit_rows     = length(weight_index),
      n_finite_terms = sum(is.finite(weight[["log_weight"]]))
    )
  }

  return(list(
    log_weight = out,
    method     = .iwmde_chen_method_label(methods, "chen"),
    partitions = partitions
  ))
}


.iwmde_chen_try_weight <- function(expr, fallback) {

  tryCatch(
    expr,
    error = function(e) fallback
  )
}


.iwmde_chen_row_active_keys <- function(context, rows) {

  if (is.null(context[["posterior_samples"]]) || length(rows) == 0L) {
    return(rep("all", length(rows)))
  }
  samples <- context[["posterior_samples"]]

  vapply(rows, function(row) {
    if (!is.finite(row) || row < 1L || row > nrow(samples)) {
      return("all")
    }

    row_values <- samples[row, , drop = FALSE]
    if (is.data.frame(row_values)) {
      row_values <- as.list(row_values[1L, , drop = FALSE])
    } else {
      row_values <- as.list(row_values[1L, ])
      names(row_values) <- colnames(samples)
    }

    .iwmde_active_key(context, row_values)
  }, character(1))
}


.iwmde_chen_valid_support <- function(support) {

  return(
    length(support) == 2L &&
      !any(is.na(support)) &&
      support[1] < support[2]
  )
}


.iwmde_chen_support_keys <- function(supports) {

  apply(supports, 1L, function(support) {
    paste(
      .iwmde_key_number(support[1]),
      .iwmde_key_number(support[2]),
      sep = ","
    )
  })
}


.iwmde_chen_method_label <- function(methods, fallback) {

  methods <- unique(methods[nzchar(methods)])
  if (length(methods) == 0L) {
    return(fallback)
  }
  if (length(methods) == 1L) {
    return(methods)
  }

  return(paste0("chen_mixed(", paste(methods, collapse = ","), ")"))
}


.iwmde_chen_gamma_log_weight_single <- function(active_values, weight_values,
                                                support) {

  out <- rep(-Inf, length(active_values))
  if (is.finite(support[1])) {
    distance_fit  <- weight_values - support[1]
    distance_eval <- active_values - support[1]
  } else if (is.finite(support[2])) {
    distance_fit  <- support[2] - weight_values
    distance_eval <- support[2] - active_values
  } else {
    return(list(log_weight = out, method = "chen_gamma"))
  }

  distance_fit <- distance_fit[
    is.finite(distance_fit) & distance_fit > 0
  ]
  keep <- is.finite(distance_eval) & distance_eval > 0
  if (length(distance_fit) < 3L || !any(keep)) {
    return(list(log_weight = out, method = "chen_gamma"))
  }

  location <- mean(distance_fit)
  variance <- stats::var(distance_fit)
  if (!is.finite(location) || !is.finite(variance) ||
      location <= 0 || variance <= sqrt(.Machine$double.eps)) {
    return(list(log_weight = out, method = "chen_gamma"))
  }

  shape <- location^2 / variance
  rate  <- location / variance
  if (!is.finite(shape) || !is.finite(rate) ||
      shape <= 0 || rate <= 0) {
    return(list(log_weight = out, method = "chen_gamma"))
  }

  out[keep] <- stats::dgamma(
    distance_eval[keep],
    shape = shape,
    rate  = rate,
    log   = TRUE
  )

  return(list(log_weight = out, method = "chen_gamma"))
}


.iwmde_chen_logit_conditional_normal_log_weight <- function(context,
                                                            parameter,
                                                            parameter_spec,
                                                            active_rows,
                                                            active_values,
                                                            weight_rows,
                                                            weight_values,
                                                            support) {

  active <- .iwmde_chen_logit_transform(active_values, support)
  weight <- .iwmde_chen_logit_transform(weight_values, support)
  out    <- rep(-Inf, length(active_values))

  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    active_rows    = active_rows,
    weight_rows    = weight_rows
  )
  if (ncol(conditioning[["fit"]]) == 0L) {
    return(.iwmde_chen_beta_log_weight(
      active_values = active_values,
      weight_values = weight_values,
      support       = support
    ))
  }

  x_fit  <- conditioning[["fit"]]
  x_eval <- conditioning[["eval"]]
  z_fit  <- weight[["z"]]

  fit_keep <- is.finite(z_fit) & stats::complete.cases(x_fit)
  if (sum(fit_keep) < 3L) {
    .iwmde_chen_conditional_stop("fewer than three finite conditioning rows")
  }

  x_fit <- x_fit[fit_keep, , drop = FALSE]
  z_fit <- z_fit[fit_keep]
  center <- colMeans(x_fit)
  scale  <- apply(x_fit, 2L, stats::sd)
  keep   <- is.finite(scale) & scale > sqrt(.Machine$double.eps)
  if (!any(keep)) {
    .iwmde_chen_conditional_stop("all conditioning columns have zero variance")
  }

  x_fit  <- sweep(x_fit[, keep, drop = FALSE], 2L, center[keep], "-")
  x_fit  <- sweep(x_fit, 2L, scale[keep], "/")
  x_eval <- sweep(x_eval[, keep, drop = FALSE], 2L, center[keep], "-")
  x_eval <- sweep(x_eval, 2L, scale[keep], "/")

  cov_mat <- stats::cov(cbind(z_fit, x_fit))
  if (!all(is.finite(cov_mat))) {
    .iwmde_chen_conditional_stop("the conditional covariance matrix is not finite")
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
    .iwmde_chen_conditional_stop("the conditional covariance matrix is not positive definite")
  }

  beta     <- sigma_zx %*% inverse
  variance <- as.numeric(sigma_zz - beta %*% t(sigma_zx))
  if (!is.finite(variance) || variance <= sqrt(.Machine$double.eps)) {
    .iwmde_chen_conditional_stop("the conditional variance is not positive")
  }

  eval_keep <- is.finite(active[["z"]]) &
    is.finite(active[["log_jacobian"]]) &
    stats::complete.cases(x_eval)
  means <- as.numeric(mean(z_fit) + x_eval[eval_keep, , drop = FALSE] %*%
    t(beta))
  out[eval_keep] <- stats::dnorm(
    active[["z"]][eval_keep],
    mean = means,
    sd   = sqrt(variance),
    log  = TRUE
  ) + active[["log_jacobian"]][eval_keep]

  return(list(log_weight = out, method = "chen_logit_conditional_normal"))
}


.iwmde_chen_logit_transform <- function(values, support) {

  lower <- support[1]
  upper <- support[2]
  width <- upper - lower
  prob  <- (values - lower) / width

  return(list(
    z            = stats::qlogis(prob),
    log_jacobian = -log(width) - log(prob) - log1p(-prob)
  ))
}


.iwmde_chen_beta_log_weight <- function(active_values, weight_values,
                                        support) {

  out   <- rep(-Inf, length(active_values))
  lower <- support[1]
  upper <- support[2]
  width <- upper - lower

  prob_fit <- (weight_values - lower) / width
  prob_fit <- prob_fit[
    is.finite(prob_fit) & prob_fit > 0 & prob_fit < 1
  ]
  prob_eval <- (active_values - lower) / width
  keep      <- is.finite(prob_eval) & prob_eval > 0 & prob_eval < 1
  if (length(prob_fit) < 3L || !any(keep)) {
    return(list(log_weight = out, method = "chen_beta"))
  }

  location <- mean(prob_fit)
  variance <- stats::var(prob_fit)
  maximum  <- location * (1 - location)
  if (!is.finite(location) || !is.finite(variance) ||
      location <= 0 || location >= 1 ||
      variance <= sqrt(.Machine$double.eps) ||
      variance >= maximum) {
    return(list(log_weight = out, method = "chen_beta"))
  }

  common <- maximum / variance - 1
  alpha  <- location * common
  beta   <- (1 - location) * common
  if (!is.finite(alpha) || !is.finite(beta) ||
      alpha <= 0 || beta <= 0) {
    return(list(log_weight = out, method = "chen_beta"))
  }

  out[keep] <- stats::dbeta(
    prob_eval[keep],
    shape1 = alpha,
    shape2 = beta,
    log    = TRUE
  ) - log(width)

  return(list(log_weight = out, method = "chen_beta"))
}


.iwmde_chen_conditional_stop <- function(reason) {

  stop(
    "IWMDE Chen conditional-normal weights are unavailable: ",
    reason,
    ".",
    call. = FALSE
  )
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
    .iwmde_chen_conditional_stop("fewer than three finite conditioning rows")
  }

  x_fit <- x_fit[fit_keep, , drop = FALSE]
  z_fit <- z_fit[fit_keep]
  center <- colMeans(x_fit)
  scale  <- apply(x_fit, 2L, stats::sd)
  keep   <- is.finite(scale) & scale > sqrt(.Machine$double.eps)
  if (!any(keep)) {
    .iwmde_chen_conditional_stop("all conditioning columns have zero variance")
  }

  x_fit  <- sweep(x_fit[, keep, drop = FALSE], 2L, center[keep], "-")
  x_fit  <- sweep(x_fit, 2L, scale[keep], "/")
  x_eval <- sweep(x_eval[, keep, drop = FALSE], 2L, center[keep], "-")
  x_eval <- sweep(x_eval, 2L, scale[keep], "/")

  cov_mat <- stats::cov(cbind(z_fit, x_fit))
  if (!all(is.finite(cov_mat))) {
    .iwmde_chen_conditional_stop("the conditional covariance matrix is not finite")
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
    .iwmde_chen_conditional_stop("the conditional covariance matrix is not positive definite")
  }

  beta     <- sigma_zx %*% inverse
  variance <- as.numeric(sigma_zz - beta %*% t(sigma_zx))
  if (!is.finite(variance) || variance <= sqrt(.Machine$double.eps)) {
    .iwmde_chen_conditional_stop("the conditional variance is not positive")
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

  context <- .iwmde_context_ensure_caches(context)
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
      fit_values  = .iwmde_parameter_column_values(
        context   = context,
        samples   = samples[weight_rows, , drop = FALSE],
        parameter = column
      ),
      eval_values = .iwmde_parameter_column_values(
        context   = context,
        samples   = samples[active_rows, , drop = FALSE],
        parameter = column
      ),
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


.iwmde_qcmde_normalizer_plan <- function(normalization_grid, transform) {

  grid_sequence <- .iwmde_qcmde_grid_sequence(
    normalization_grid = normalization_grid,
    transform          = transform
  )
  if (length(grid_sequence) == 0L) {
    grid_sequence <- list(normalization_grid)
  }

  all_z <- sort(unique(unlist(lapply(grid_sequence, `[[`, "z"),
                              use.names = FALSE)))
  all_x <- .iwmde_from_internal(all_z, transform)
  keep  <- is.finite(all_z) & is.finite(all_x)
  all_z <- all_z[keep]
  all_x <- all_x[keep]

  if (length(all_z) < 2L || any(diff(all_z) <= 0)) {
    all_z <- as.numeric(normalization_grid[["z"]])
    all_x <- as.numeric(normalization_grid[["x"]])
  }

  all_grid <- list(
    x            = all_x,
    z            = all_z,
    log_jacobian = .iwmde_log_jacobian(all_z, transform)
  )

  grid_sequence <- lapply(grid_sequence, function(grid) {
    index <- match(grid[["z"]], all_grid[["z"]])
    keep  <- !is.na(index)
    grid[["x"]]            <- grid[["x"]][keep]
    grid[["z"]]            <- grid[["z"]][keep]
    grid[["log_jacobian"]] <- grid[["log_jacobian"]][keep]
    grid[["all_index"]]    <- index[keep]
    grid
  })
  grid_sequence <- grid_sequence[
    vapply(grid_sequence, function(grid) length(grid[["z"]]) >= 2L, logical(1))
  ]
  if (length(grid_sequence) == 0L) {
    grid_sequence <- list(c(all_grid, list(all_index = seq_along(all_z))))
  }

  return(list(
    grid_sequence      = grid_sequence,
    all_grid           = all_grid,
    pilot_grid         = grid_sequence[[1L]],
    final_grid         = grid_sequence[[length(grid_sequence)]],
    pilot_index        = 1L,
    final_index        = length(grid_sequence),
    validation_index   = length(grid_sequence),
    n_refinement_steps = length(grid_sequence) - 1L
  ))
}


.iwmde_qcmde_grid_sequence <- function(normalization_grid, transform,
                                       max_refinement_steps = 3L) {

  z <- sort(unique(as.numeric(normalization_grid[["z"]])))
  z <- z[is.finite(z)]
  if (length(z) < 2L || any(diff(z) <= 0)) {
    return(list())
  }

  width <- diff(range(z))
  if (!is.finite(width) || width <= 0) {
    return(list())
  }

  n_base <- length(z)
  out    <- vector("list", max_refinement_steps + 1L)
  out[[1L]] <- .iwmde_qcmde_grid_from_z(z, transform)

  for (step in seq_len(max_refinement_steps)) {
    extension <- width * .iwmde_qcmde_extension_fraction(step)
    n_points  <- .iwmde_qcmde_refinement_points(
      n_base = n_base,
      step   = step
    )
    z_step <- seq(
      min(z) - extension,
      max(z) + extension,
      length.out = n_points
    )
    z_step <- sort(unique(c(z, z_step)))
    out[[step + 1L]] <- .iwmde_qcmde_grid_from_z(z_step, transform)
  }

  out <- out[!vapply(out, is.null, logical(1))]

  return(out)
}


.iwmde_qcmde_grid_from_z <- function(z, transform) {

  x    <- .iwmde_from_internal(z, transform)
  keep <- is.finite(z) & is.finite(x)
  z    <- z[keep]
  x    <- x[keep]
  if (length(z) < 2L || any(diff(z) <= 0)) {
    return(NULL)
  }

  return(list(
    x            = x,
    z            = z,
    log_jacobian = .iwmde_log_jacobian(z, transform)
  ))
}


.iwmde_qcmde_extension_fraction <- function(step) {

  fractions <- c(.25, .50, 1)
  step      <- min(length(fractions), max(1L, as.integer(step[[1L]])))

  return(fractions[[step]])
}


.iwmde_qcmde_refinement_points <- function(n_base, step) {

  multiplier <- c(1.5, 2, 2.5)
  step       <- min(length(multiplier), max(1L, as.integer(step[[1L]])))

  return(max(2L, as.integer(ceiling(n_base * multiplier[[step]]))))
}


.iwmde_qcmde_select_refinement <- function(log_q_display,
                                           log_normalizer_sequence,
                                           active_mass, denominator) {

  n_sequence <- length(log_normalizer_sequence)
  if (n_sequence <= 1L) {
    return(list(
      pilot_index        = 1L,
      final_index        = 1L,
      validation_index   = 1L,
      n_refinement_steps = 0L
    ))
  }

  if (n_sequence == 2L) {
    return(list(
      pilot_index        = 1L,
      final_index        = 2L,
      validation_index   = 2L,
      n_refinement_steps = 1L
    ))
  }

  for (index in seq.int(2L, n_sequence - 1L)) {
    candidate_y <- .iwmde_qcmde_density_from_normalizer(
      log_q_display  = log_q_display,
      log_normalizer = log_normalizer_sequence[[index]],
      active_mass    = active_mass,
      denominator    = denominator
    )
    validation_y <- .iwmde_qcmde_density_from_normalizer(
      log_q_display  = log_q_display,
      log_normalizer = log_normalizer_sequence[[index + 1L]],
      active_mass    = active_mass,
      denominator    = denominator
    )
    change <- .iwmde_qcmde_ordinate_change(
      pilot_y = candidate_y,
      final_y = validation_y
    )
    max_change <- .iwmde_max_or_na(change[["relative"]])
    if (is.finite(max_change) &&
        max_change <= .iwmde_qcmde_refinement_target()) {
      return(list(
        pilot_index        = index - 1L,
        final_index        = index,
        validation_index   = index + 1L,
        n_refinement_steps = index - 1L
      ))
    }
  }

  return(list(
    pilot_index        = max(1L, n_sequence - 2L),
    final_index        = n_sequence - 1L,
    validation_index   = n_sequence,
    n_refinement_steps = n_sequence - 2L
  ))
}


.iwmde_qcmde_density_from_normalizer <- function(log_q_display,
                                                 log_normalizer,
                                                 active_mass,
                                                 denominator) {

  keep_rows <- is.finite(log_normalizer)
  return(.iwmde_qcmde_pilot_density(
    log_q_display  = log_q_display,
    log_normalizer = log_normalizer,
    keep_rows      = keep_rows,
    active_mass    = active_mass,
    denominator    = denominator
  ))
}


.iwmde_qcmde_refinement_target <- function() {

  return(.025)
}


.iwmde_qcmde_pilot_density <- function(log_q_display, log_normalizer,
                                       keep_rows, active_mass,
                                       denominator) {

  if (!any(keep_rows)) {
    return(rep(NA_real_, nrow(log_q_display)))
  }

  density_terms <- .iwmde_density_aggregate(
    log_terms = sweep(
      log_q_display[, keep_rows, drop = FALSE],
      2L,
      log_normalizer[keep_rows],
      "-"
    ),
    active_mass = active_mass,
    denominator = denominator
  )

  return(density_terms[["y"]])
}


.iwmde_qcmde_normalizer_change <- function(initial_log_normalizer,
                                           final_log_normalizer) {

  finite <- is.finite(initial_log_normalizer) &
    is.finite(final_log_normalizer)
  if (!any(finite)) {
    return(list(max = NA_real_, p95 = NA_real_, median = NA_real_))
  }

  relative_change <- abs(expm1(
    final_log_normalizer[finite] - initial_log_normalizer[finite]
  ))

  return(list(
    max    = max(relative_change),
    p95    = stats::quantile(relative_change, .95, names = FALSE, type = 8),
    median = stats::median(relative_change)
  ))
}


.iwmde_qcmde_ordinate_change <- function(pilot_y, final_y) {

  relative_change <- rep(NA_real_, length(final_y))
  log_change      <- rep(NA_real_, length(final_y))

  positive_pilot <- is.finite(pilot_y) & pilot_y > 0
  finite_final   <- is.finite(final_y)
  relative_rows  <- positive_pilot & finite_final
  if (any(relative_rows)) {
    relative_change[relative_rows] <- abs(
      final_y[relative_rows] / pilot_y[relative_rows] - 1
    )
  }

  zero_to_positive <- is.finite(pilot_y) & pilot_y == 0 &
    finite_final & final_y > 0
  relative_change[zero_to_positive] <- Inf

  positive_final <- finite_final & final_y > 0
  log_rows       <- positive_pilot & positive_final
  if (any(log_rows)) {
    log_change[log_rows] <- abs(log(final_y[log_rows]) - log(pilot_y[log_rows]))
  }

  return(list(relative = relative_change, log = log_change))
}


.iwmde_normalization_density <- function(log_q_norm, log_normalizer,
                                         log_jacobian, normalization_grid,
                                         active_mass, denominator) {

  y <- numeric(length(normalization_grid))
  for (g in seq_along(normalization_grid)) {
    log_terms <- log_q_norm[g, ] + log_jacobian[g] - log_normalizer
    finite    <- is.finite(log_terms)
    if (any(finite)) {
      max_term         <- max(log_terms[finite])
      scaled_terms     <- exp(log_terms[finite] - max_term)
      y[g]             <- active_mass * exp(max_term) *
        sum(scaled_terms) / denominator
    }
  }

  return(y)
}


.iwmde_row_drop_fraction <- function(n_candidate_rows, n_normalized_rows) {

  n_candidate_rows  <- as.numeric(n_candidate_rows)[1]
  n_normalized_rows <- as.numeric(n_normalized_rows)[1]
  if (!is.finite(n_candidate_rows) || n_candidate_rows <= 0 ||
      !is.finite(n_normalized_rows) || n_normalized_rows < 0) {
    return(NA_real_)
  }

  return(min(1, max(0, 1 - n_normalized_rows / n_candidate_rows)))
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
