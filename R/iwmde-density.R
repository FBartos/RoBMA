# ============================================================================ #
# IWMDE Row States, Density Aggregation, and Chen Weights
# ============================================================================ #

.iwmde_row_states <- function(context, rows, parameter = NULL,
                              parameter_spec = NULL, estimator = NULL) {

  lapply(rows, function(row) {
    tryCatch(
      .iwmde_row_state(
        context        = context,
        row_index      = row,
        parameter      = parameter,
        parameter_spec = parameter_spec
      ),
      error = function(e) {
        if (inherits(e, "iwmde_construction_error") || is.null(estimator)) {
          stop(e)
        }
        .iwmde_stop_construction_failure(
          estimator = estimator,
          parameter = parameter,
          rows      = row,
          stage     = "baseline joint-density evaluation",
          detail    = conditionMessage(e)
        )
      }
    )
  })
}


.iwmde_stop_construction_failure <- function(estimator, parameter, rows,
                                             stage, detail = NULL,
                                             chain_coverage = NULL) {

  estimator_label <- if (identical(estimator, "q_grid_cmde") ||
                         identical(estimator, "qCMDE")) {
    "qCMDE"
  } else {
    "IWMDE"
  }
  rows <- unique(as.integer(rows))
  row_text <- if (length(rows) == 0L) {
    "an unknown posterior row"
  } else {
    shown <- utils::head(rows, 8L)
    suffix <- if (length(rows) > length(shown)) {
      paste0(" (and ", length(rows) - length(shown), " more)")
    } else {
      ""
    }
    paste0(
      "posterior row", if (length(rows) == 1L) " " else "s ",
      paste(shown, collapse = ", "), suffix
    )
  }
  detail_text <- if (!is.null(detail) && nzchar(as.character(detail)[[1L]])) {
    paste0(
      ": ",
      sub("[.]+$", "", gsub("[\r\n\t]+", " ", as.character(detail)[[1L]]))
    )
  } else {
    ""
  }

  message <- paste0(
    estimator_label, " construction failed for target '", parameter,
    "' at ", row_text, " during ", stage, detail_text, "."
  )
  condition <- structure(
    list(
      message        = message,
      call           = NULL,
      estimator      = estimator,
      target         = parameter,
      posterior_rows = rows,
      stage          = stage,
      detail         = detail,
      chain_coverage = chain_coverage
    ),
    class = c("iwmde_construction_error", "error", "condition")
  )
  stop(condition)
}


.iwmde_validate_log_grid <- function(log_q_grid, estimator, parameter, rows,
                                     n_values, stage) {

  if (!is.numeric(log_q_grid) || !is.matrix(log_q_grid) ||
      nrow(log_q_grid) != n_values || ncol(log_q_grid) != length(rows)) {
    .iwmde_stop_construction_failure(
      estimator = estimator,
      parameter = parameter,
      rows      = rows,
      stage     = stage,
      detail    = "joint log-density evaluation returned an invalid matrix"
    )
  }

  invalid <- is.na(log_q_grid) | log_q_grid == Inf
  if (any(invalid)) {
    bad_columns <- unique(which(invalid, arr.ind = TRUE)[, "col"])
    .iwmde_stop_construction_failure(
      estimator = estimator,
      parameter = parameter,
      rows      = rows[bad_columns],
      stage     = stage,
      detail    = "joint log density was undefined or positive-infinite"
    )
  }

  invisible(log_q_grid)
}
.iwmde_likelihood_mode <- function(parameter, parameter_spec = NULL,
                                   context = NULL) {

  if (.iwmde_context_uses_local_likelihood(context)) {
    return("conditional")
  }

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


.iwmde_context_uses_local_likelihood <- function(context) {

  if (is.null(context) || is.null(context[["data"]])) {
    return(FALSE)
  }

  return(.data_outcome_type(context[["data"]]) %in% c("bin", "pois"))
}


.iwmde_state_scope <- function(parameter, parameter_spec = NULL,
                               context = NULL) {

  if (identical(
    .iwmde_likelihood_mode(parameter, parameter_spec, context),
    "conditional"
  )) {
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

  likelihood_mode <- .iwmde_likelihood_mode(
    parameter      = parameter,
    parameter_spec = parameter_spec,
    context        = context
  )
  state_scope <- .iwmde_state_scope(
    parameter      = parameter,
    parameter_spec = parameter_spec,
    context        = context
  )
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
    !.iwmde_parameter_controls_sampled_random_sd(context, parameter) &&
    is.finite(baseline_focal_log_prior) &&
    .iwmde_can_use_focal_prior_delta(focal_prior)
  base_state[["baseline_log_q"]]               <- baseline_log_q
  base_state[["likelihood_mode"]]              <- likelihood_mode
  base_state[["state_scope"]]                  <- state_scope

  return(.iwmde_new_row_state(base_state))
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
             identical(likelihood_mode, "marginal") &&
             !.iwmde_marginal_likelihood_requires_row(context)) {
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


.iwmde_marginal_likelihood_requires_row <- function(context) {

  .is_data_known_v(context[["data"]]) && .is_data_random(context[["data"]])
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
                                estimator_rows = seq_along(row_states),
                                population_rows = estimator_rows,
                                chain_id = rep(1L, length(estimator_rows)),
                                expected_chain_ids = unique(chain_id),
                                n_candidate_rows = length(row_states)) {

  n_display        <- length(display_grid)
  y                <- numeric(n_display)
  finite_terms     <- integer(n_display)
  max_log_ratio    <- numeric(n_display)
  ess              <- numeric(n_display)
  max_weight_share <- numeric(n_display)
  n_input_rows     <- length(row_states)
  selected_rows    <- estimator_rows
  n_candidate_rows <- as.integer(n_candidate_rows[[1L]])
  if (!is.finite(n_candidate_rows) || n_candidate_rows != n_input_rows ||
      length(estimator_rows) != n_input_rows) {
    .iwmde_stop_construction_failure(
      estimator = "q_grid_cmde",
      parameter = parameter,
      rows      = estimator_rows,
      stage     = "row-contract validation",
      detail    = "candidate rows, estimator rows, and row states are inconsistent"
    )
  }
  normalizer_plan  <- .iwmde_qcmde_normalizer_plan(
    normalization_grid = normalization_grid,
    transform          = transform
  )
  final_grid       <- normalizer_plan[["final_grid"]]
  pilot_grid       <- normalizer_plan[["pilot_grid"]]
  all_grid         <- normalizer_plan[["all_grid"]]
  q_grid           <- c(display_grid, all_grid[["x"]])
  log_q_grid <- tryCatch(
    .iwmde_log_q_grid(
      context     = context,
      parameter   = parameter,
      values      = q_grid,
      row_states  = row_states,
      replacement = replacement
    ),
    error = function(e) {
      if (inherits(e, "iwmde_construction_error")) {
        stop(e)
      }
      .iwmde_stop_construction_failure(
        estimator = "q_grid_cmde",
        parameter = parameter,
        rows      = estimator_rows,
        stage     = "joint-density grid evaluation",
        detail    = conditionMessage(e)
      )
    }
  )
  .iwmde_validate_log_grid(
    log_q_grid = log_q_grid,
    estimator  = "q_grid_cmde",
    parameter  = parameter,
    rows       = estimator_rows,
    n_values   = length(q_grid),
    stage      = "joint-density grid evaluation"
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
  if (any(!final_finite)) {
    .iwmde_stop_construction_failure(
      estimator = "q_grid_cmde",
      parameter = parameter,
      rows      = estimator_rows[!final_finite],
      stage     = "conditional-density normalization",
      detail    = paste0(
        "no finite positive normalizer was obtained after ",
        normalizer_plan[["n_refinement_steps"]], " refinement step(s)"
      )
    )
  }
  if (any(!validation_finite)) {
    .iwmde_stop_construction_failure(
      estimator = "q_grid_cmde",
      parameter = parameter,
      rows      = estimator_rows[!validation_finite],
      stage     = "conditional-density normalization validation",
      detail    = "the validation grid did not produce a finite positive normalizer"
    )
  }
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
  estimator_rows <- estimator_rows[keep_rows]
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
      sampling_mcse          = rep(NA_real_, n_display),
      sampling_relative_mcse = rep(NA_real_, n_display),
      sampling_fraction      = if (length(population_rows) == 0L) {
        1
      } else {
        length(estimator_rows) / length(population_rows)
      },
      sampling_uncertainty_type = "finite_population_srswor",
      mcmc_uncertainty_scope = "selected_continuous_rows_only",
      mcmc_uncertainty_status = "unavailable",
      mcmc_uncertainty_reason = "no normalized continuous contributions",
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
      pilot_normalization_integral      = NA_real_,
      final_normalization_integral      = NA_real_,
      normalization_relative_error      = NA_real_,
      normalization_scale               = transform[["type"]],
      normalization_mass_ratio          = NA_real_,
      pilot_y                           = pilot_y,
      validation_y                      = validation_y,
      ordinate_relative_change          = rep(NA_real_, n_display),
      ordinate_log_change               = rep(NA_real_, n_display),
      pilot_ordinate_relative_change    = rep(NA_real_, n_display),
      pilot_ordinate_log_change         = rep(NA_real_, n_display),
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
    log_terms         = sweep(log_q_display, 2L, log_normalizer, "-"),
    active_mass       = active_mass,
    denominator       = n_candidate_rows,
    contribution_rows = estimator_rows,
    sampling_population_rows = population_rows,
    chain_id = chain_id[match(estimator_rows, selected_rows)],
    expected_chain_ids = expected_chain_ids
  )
  y                <- density_terms[["y"]]
  finite_terms     <- density_terms[["finite_terms"]]
  max_log_ratio    <- density_terms[["max_log_ratio"]]
  ess              <- density_terms[["ess"]]
  max_weight_share <- density_terms[["max_weight_share"]]
  contributions    <- density_terms[["contributions"]]

  mcse_data      <- .iwmde_batch_mcse(contributions)
  ess            <- mcse_data[["ess"]]
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
  pilot_normalization_integral <- .iwmde_trapz(
    pilot_grid[["z"]],
    norm_y_initial
  )
  final_normalization_integral <- .iwmde_trapz(
    final_grid[["z"]],
    norm_y_final
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
    sampling_mcse            = density_terms[["sampling_mcse"]],
    sampling_relative_mcse   = density_terms[["sampling_relative_mcse"]],
    sampling_fraction        = density_terms[["sampling_fraction"]],
    sampling_uncertainty_type =
      density_terms[["sampling_uncertainty_type"]],
    mcmc_uncertainty_scope   = mcse_data[["uncertainty_scope"]],
    mcmc_uncertainty_status  = mcse_data[["uncertainty_status"]],
    mcmc_uncertainty_reason  = mcse_data[["uncertainty_reason"]],
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
    pilot_normalization_integral      = pilot_normalization_integral,
    final_normalization_integral      = final_normalization_integral,
    normalization_relative_error      = abs(
      final_normalization_integral / active_mass - 1
    ),
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
                                 population_rows = active_rows,
                                 chain_id = rep(1L, length(active_rows)),
                                 expected_chain_ids = unique(chain_id),
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
  if (!is.finite(n_candidate_rows) || n_candidate_rows != n_input_rows ||
      length(active_rows) != n_input_rows ||
      length(active_values) != n_input_rows) {
    .iwmde_stop_construction_failure(
      estimator = "iwmde",
      parameter = parameter,
      rows      = active_rows,
      stage     = "row-contract validation",
      detail    = "candidate rows, active values, and row states are inconsistent"
    )
  }
  weight <- tryCatch(
    .iwmde_chen_log_weight(
      context        = context,
      parameter      = parameter,
      parameter_spec = parameter_spec,
      active_rows    = active_rows,
      active_values  = active_values,
      weight_rows    = weight_rows,
      weight_values  = weight_values,
      support        = support
    ),
    error = function(e) {
      if (inherits(e, "iwmde_construction_error")) {
        stop(e)
      }
      .iwmde_stop_construction_failure(
        estimator = "iwmde",
        parameter = parameter,
        rows      = active_rows,
        stage     = "proposal-density construction",
        detail    = conditionMessage(e)
      )
    }
  )
  weight_fallbacks <- list(
    count   = if (is.null(weight[["fallback_count"]])) {
      0L
    } else {
      weight[["fallback_count"]]
    },
    rows    = if (is.null(weight[["fallback_rows"]])) {
      0L
    } else {
      weight[["fallback_rows"]]
    },
    from    = if (is.null(weight[["fallback_from"]])) {
      character()
    } else {
      weight[["fallback_from"]]
    },
    reasons = if (is.null(weight[["fallback_reasons"]])) {
      structure(integer(), names = character())
    } else {
      weight[["fallback_reasons"]]
    }
  )
  raw_log_weight    <- weight[["log_weight"]]
  if (!is.numeric(raw_log_weight) ||
      length(raw_log_weight) != length(active_rows)) {
    .iwmde_stop_construction_failure(
      estimator = "iwmde",
      parameter = parameter,
      rows      = active_rows,
      stage     = "proposal-density construction",
      detail    = "the proposal log-density vector has an invalid length or type"
    )
  }
  keep_rows <- is.finite(raw_log_weight)
  if (any(!keep_rows)) {
    .iwmde_stop_construction_failure(
      estimator = "iwmde",
      parameter = parameter,
      rows      = active_rows[!keep_rows],
      stage     = "proposal-density construction",
      detail    = "the normalized proposal density was zero or non-finite at an evaluation row"
    )
  }
  contribution_rows <- active_rows[keep_rows]
  n_dropped_weight  <- 0L
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
      sampling_mcse          = rep(NA_real_, n_display),
      sampling_relative_mcse = rep(NA_real_, n_display),
      sampling_fraction      = if (length(population_rows) == 0L) {
        1
      } else {
        length(active_rows) / length(population_rows)
      },
      sampling_uncertainty_type = "finite_population_srswor",
      mcmc_uncertainty_scope = "selected_continuous_rows_only",
      mcmc_uncertainty_status = "unavailable",
      mcmc_uncertainty_reason = "no normalized continuous contributions",
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
      support_grid_normalization_integral = NA_real_,
      normalization_relative_error      = NA_real_,
      normalization_scale               = NA_character_,
      normalization_mass_ratio          = NA_real_,
      max_normalizer_relative_change    = NA_real_,
      max_quadrature_relative_change    = NA_real_,
      median_normalizer_relative_change = NA_real_,
      normalization_refined_points      = 0L,
      normalization_refined_range       = c(NA_real_, NA_real_),
      n_dropped_weight                  = n_dropped_weight,
      weight_partitions                 = weight[["partitions"]],
      n_weight_fallbacks                 = weight_fallbacks[["count"]],
      n_weight_fallback_rows             = weight_fallbacks[["rows"]],
      weight_fallback_from               = weight_fallbacks[["from"]],
      weight_fallback_reasons            = weight_fallbacks[["reasons"]],
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
  log_q_grid <- tryCatch(
    .iwmde_log_q_grid(
      context     = context,
      parameter   = parameter,
      values      = q_grid,
      row_states  = row_states,
      replacement = replacement
    ),
    error = function(e) {
      if (inherits(e, "iwmde_construction_error")) {
        stop(e)
      }
      .iwmde_stop_construction_failure(
        estimator = "iwmde",
        parameter = parameter,
        rows      = contribution_rows,
        stage     = "joint-density ordinate evaluation",
        detail    = conditionMessage(e)
      )
    }
  )
  .iwmde_validate_log_grid(
    log_q_grid = log_q_grid,
    estimator  = "iwmde",
    parameter  = parameter,
    rows       = contribution_rows,
    n_values   = length(q_grid),
    stage      = "joint-density ordinate evaluation"
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
    active_mass       = active_mass,
    denominator       = n_candidate_rows,
    contribution_rows = contribution_rows,
    sampling_population_rows = population_rows,
    chain_id = chain_id[match(contribution_rows, active_rows)],
    expected_chain_ids = expected_chain_ids
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
  ess           <- mcse_data[["ess"]]
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
    sampling_mcse            = density_terms[["sampling_mcse"]],
    sampling_relative_mcse   = density_terms[["sampling_relative_mcse"]],
    sampling_fraction        = density_terms[["sampling_fraction"]],
    sampling_uncertainty_type =
      density_terms[["sampling_uncertainty_type"]],
    mcmc_uncertainty_scope   = mcse_data[["uncertainty_scope"]],
    mcmc_uncertainty_status  = mcse_data[["uncertainty_status"]],
    mcmc_uncertainty_reason  = mcse_data[["uncertainty_reason"]],
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
    support_grid_normalization_integral = normalization[["integral"]],
    normalization_relative_error      = abs(
      normalization[["integral"]] / active_mass - 1
    ),
    normalization_scale               = normalization[["scale_type"]],
    normalization_mass_ratio          = normalization[["mass_ratio"]],
    max_normalizer_relative_change    = NA_real_,
    max_quadrature_relative_change    = quadrature_change,
    median_normalizer_relative_change = NA_real_,
    normalization_refined_points      = 0L,
    normalization_refined_range       = c(NA_real_, NA_real_),
    n_dropped_weight                  = n_dropped_weight,
    weight_partitions                 = weight[["partitions"]],
    n_weight_fallbacks                 = weight_fallbacks[["count"]],
    n_weight_fallback_rows             = weight_fallbacks[["rows"]],
    weight_fallback_from               = weight_fallbacks[["from"]],
    weight_fallback_reasons            = weight_fallbacks[["reasons"]],
    integral_mcse                     = integral_mcse[["mcse"]],
    integral_relative_mcse            = integral_mcse[["relative_mcse"]],
    batch_size                        = mcse_data[["batch_size"]],
    n_batches                         = mcse_data[["n_batches"]],
    estimator                         = "iwmde",
    weight_method                     = weight[["method"]]
  ))
}
