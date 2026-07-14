# ============================================================================ #
# IWMDE Estimate Facade
# ============================================================================ #

.iwmde_estimate <- function(context, parameter, density_method, density_control,
                            outputs = c("density", "ordinate"), values = NULL,
                            parameter_spec = NULL, metadata = NULL,
                            cache = NULL) {

  context <- .iwmde_context_ensure_caches(context)
  .iwmde_check_context_density_method_supported(context, density_method)
  outputs <- unique(match.arg(
    outputs,
    c("density", "ordinate"),
    several.ok = TRUE
  ))
  if (identical(outputs, "ordinate")) {
    return(.iwmde_estimate_adaptive_ordinate(
      context         = context,
      parameter       = parameter,
      density_method  = density_method,
      density_control = density_control,
      values          = values,
      parameter_spec  = parameter_spec,
      metadata        = metadata,
      cache           = cache
    ))
  }

  plan <- .iwmde_plan(
    context         = context,
    parameter       = parameter,
    density_method  = density_method,
    density_control = density_control,
    outputs         = outputs,
    values          = values,
    parameter_spec  = parameter_spec,
    metadata        = metadata
  )

  return(.iwmde_estimate_from_plan(
    context = context,
    plan    = plan,
    cache   = cache
  ))
}


.iwmde_estimate_adaptive_ordinate <- function(
    context, parameter, density_method, density_control, values,
    parameter_spec = NULL, metadata = NULL, cache = NULL) {

  control <- .iwmde_density_control_resolve(
    density_method  = density_method,
    density_control = density_control,
    purpose         = "ordinate"
  )
  hard_cap <- control[["max_samples"]]
  budget   <- min(control[["initial_samples"]], hard_cap)
  history  <- list()

  repeat {
    plan <- .iwmde_plan(
      context         = context,
      parameter       = parameter,
      density_method  = density_method,
      density_control = control,
      outputs         = "ordinate",
      values          = values,
      parameter_spec  = parameter_spec,
      metadata        = metadata,
      row_budget      = as.integer(budget)
    )
    out <- .iwmde_estimate_from_plan(
      context = context,
      plan    = plan,
      cache   = cache
    )

    diagnostic <- out[["diagnostics"]][["ordinate"]]
    metrics    <- .iwmde_adaptation_metrics(diagnostic)
    eligible   <- length(plan[["rows"]][["continuous_rows"]])
    achieved   <- plan[["rows"]][["n_candidate_rows"]]
    if (!is.finite(eligible) || eligible < 1L) {
      eligible <- achieved
    }
    effective_cap <- min(hard_cap, eligible)
    precision_target_met <- is.finite(metrics[["relative_mcse"]]) &&
      metrics[["relative_mcse"]] <= control[["target_relative_mcse"]]
    bf_grade_met <- is.null(.iwmde_diagnostics_bf_failure_reason(
      diagnostic[["diagnostics"]]
    ))
    target_met <- precision_target_met && bf_grade_met
    all_rows_used <- is.finite(eligible) && achieved >= eligible
    hard_cap_reached <- is.finite(hard_cap) && hard_cap < eligible &&
      achieved >= hard_cap

    history[[length(history) + 1L]] <- data.frame(
      requested_row_budget = as.integer(budget),
      evaluated_rows       = as.integer(achieved),
      relative_mcse        = metrics[["relative_mcse"]],
      ess                  = metrics[["ess"]],
      max_weight_share     = metrics[["max_weight_share"]],
      precision_target_met = precision_target_met,
      bf_grade_met         = bf_grade_met,
      stringsAsFactors     = FALSE
    )

    exhausted <- all_rows_used || hard_cap_reached ||
      !is.finite(effective_cap) || achieved >= effective_cap
    if (target_met || exhausted || !identical(diagnostic[["status"]], "ok")) {
      break
    }

    next_budget <- .iwmde_next_row_budget(
      current              = budget,
      relative_mcse        = metrics[["relative_mcse"]],
      target_relative_mcse = control[["target_relative_mcse"]],
      cap                  = effective_cap
    )
    if (!is.finite(next_budget) || next_budget <= budget) {
      break
    }
    budget <- next_budget
  }

  adaptation <- list(
    adaptive             = TRUE,
    initial_row_budget   = control[["initial_samples"]],
    achieved_row_budget  = as.integer(achieved),
    eligible_rows        = as.integer(eligible),
    hard_cap             = hard_cap,
    hard_cap_reached     = hard_cap_reached,
    all_rows_used        = all_rows_used,
    target_relative_mcse = control[["target_relative_mcse"]],
    precision_target_met = precision_target_met,
    bf_grade_met         = bf_grade_met,
    target_met           = target_met,
    n_steps              = length(history),
    history              = do.call(rbind, history)
  )

  return(.iwmde_estimate_attach_adaptation(out, adaptation))
}


.iwmde_adaptation_metrics <- function(diagnostic) {

  diagnostics <- diagnostic[["diagnostics"]]
  if (!is.list(diagnostics)) {
    diagnostics <- list()
  }

  return(list(
    relative_mcse = .iwmde_diagnostic_scalar_any(
      diagnostics,
      c("bf_relative_mcse", "relative_mcse")
    ),
    ess = .iwmde_diagnostic_scalar_any(
      diagnostics,
      c("bf_ess", "ess")
    ),
    max_weight_share = .iwmde_diagnostic_scalar_any(
      diagnostics,
      c("bf_max_weight_share", "max_weight_share")
    )
  ))
}


.iwmde_next_row_budget <- function(current, relative_mcse,
                                    target_relative_mcse, cap) {

  growth <- 2
  if (is.finite(relative_mcse) && relative_mcse > 0 &&
      is.finite(target_relative_mcse) && target_relative_mcse > 0) {
    projected <- 1.2 * (relative_mcse / target_relative_mcse)^2
    growth    <- min(4, max(1.5, projected))
  }
  next_budget <- max(current + 20, ceiling(current * growth))

  return(min(as.integer(next_budget), cap))
}


.iwmde_estimate_attach_adaptation <- function(estimate, adaptation) {

  diagnostic <- estimate[["diagnostics"]][["ordinate"]]
  if (is.list(diagnostic)) {
    diagnostic[["adaptation"]] <- adaptation
    if (is.list(diagnostic[["diagnostics"]])) {
      for (name in setdiff(names(adaptation), "history")) {
        diagnostic[["diagnostics"]][[name]] <- adaptation[[name]]
      }
    }
    estimate[["diagnostics"]][["ordinate"]] <- diagnostic
  }

  for (ordinate_name in c(
    "posterior_ordinate",
    "rejected_posterior_ordinate"
  )) {
    ordinate <- estimate[[ordinate_name]]
    if (is.list(ordinate) && is.list(ordinate[["diagnostics"]])) {
      for (name in setdiff(names(adaptation), "history")) {
        ordinate[["diagnostics"]][[name]] <- adaptation[[name]]
      }
      ordinate[["diagnostics"]][["adaptation_history"]] <-
        adaptation[["history"]]
      estimate[[ordinate_name]] <- ordinate
    }
  }

  estimate[["adaptation"]] <- adaptation

  return(estimate)
}


.iwmde_estimate_from_plan <- function(context, plan, cache = NULL) {

  key <- plan[["plan_key"]]
  if (.iwmde_estimate_cache_has(cache, key)) {
    return(.iwmde_estimate_cache_get(cache, key))
  }

  diagnostic_cache <- .iwmde_estimate_diagnostic_cache(cache)
  execution_cache  <- new.env(parent = emptyenv())
  density_diagnostic <- NULL
  ordinate_diagnostic <- NULL

  if (isTRUE(plan[["outputs"]][["need_density"]])) {
    density_diagnostic <- .iwmde_execute_plan_diagnostic(
      context          = context,
      plan             = plan,
      output           = "density",
      execution_cache  = execution_cache,
      diagnostic_cache = diagnostic_cache
    )
    density_diagnostic <- .iwmde_attach_plan_to_diagnostic(
      diagnostic = density_diagnostic,
      plan       = plan
    )
  }

  if (isTRUE(plan[["outputs"]][["need_ordinate"]])) {
    ordinate_diagnostic <- .iwmde_execute_plan_diagnostic(
      context          = context,
      plan             = plan,
      output           = "ordinate",
      execution_cache  = execution_cache,
      diagnostic_cache = diagnostic_cache
    )
    ordinate_diagnostic <- .iwmde_attach_plan_to_diagnostic(
      diagnostic = ordinate_diagnostic,
      plan       = plan
    )
  }

  out <- .iwmde_estimate_result(
    plan                = plan,
    density_diagnostic  = density_diagnostic,
    ordinate_diagnostic = ordinate_diagnostic
  )
  .iwmde_estimate_cache_set(cache, key, out)

  return(out)
}


.iwmde_attach_plan_to_diagnostic <- function(diagnostic, plan) {

  if (!is.list(diagnostic)) {
    return(diagnostic)
  }

  diagnostic[["plan"]] <- plan
  diagnostic[["iwmde_provenance"]] <- .iwmde_plan_provenance(plan)

  return(diagnostic)
}


.iwmde_estimate_result <- function(plan, density_diagnostic = NULL,
                                   ordinate_diagnostic = NULL) {

  density_attribute <- NULL
  ordinate_attribute <- NULL
  rejected_ordinate_attribute <- NULL

  if (!is.null(density_diagnostic)) {
    density_attribute <- .iwmde_posterior_density_attribute(
      diagnostic      = density_diagnostic,
      density_method  = plan[["density_method"]],
      density_control = plan[["control"]],
      metadata        = plan[["target"]][["metadata"]]
    )
  }
  if (!is.null(ordinate_diagnostic)) {
    candidate_ordinate <- .iwmde_posterior_ordinate_attribute(
      diagnostic      = ordinate_diagnostic,
      density_method  = plan[["density_method"]],
      density_control = plan[["control"]],
      metadata        = plan[["target"]][["metadata"]],
      allow_rejected  = TRUE
    )
    if (!is.null(candidate_ordinate) &&
        .iwmde_posterior_ordinate_supports_bf(candidate_ordinate)) {
      ordinate_attribute <- candidate_ordinate
    } else {
      rejected_ordinate_attribute <- candidate_ordinate
    }
  }

  out <- list(
    target = plan[["target"]],
    plan   = plan,
    density = if (is.null(density_diagnostic)) {
      NULL
    } else {
      density_diagnostic[["iwmde"]]
    },
    ordinates = if (is.null(ordinate_diagnostic)) {
      NULL
    } else {
      ordinate_diagnostic[["iwmde"]]
    },
    diagnostics = list(
      density  = density_diagnostic,
      ordinate = ordinate_diagnostic
    ),
    posterior_density            = density_attribute,
    posterior_ordinate           = ordinate_attribute,
    rejected_posterior_ordinate  = rejected_ordinate_attribute,
    provenance                   = .iwmde_plan_provenance(plan)
  )
  class(out) <- c("iwmde_estimate", "list")

  return(out)
}


.iwmde_execute_plan_diagnostic <- function(context, plan,
                                           output = c("density", "ordinate"),
                                           execution_cache = NULL,
                                           diagnostic_cache = NULL) {

  context <- .iwmde_context_ensure_caches(context)
  output    <- match.arg(output)
  parameter <- plan[["target"]][["parameter"]]
  cache_key <- .iwmde_plan_output_cache_key(plan, output)
  if (.iwmde_cache_has(diagnostic_cache, cache_key)) {
    return(.iwmde_relabel_diagnostic(
      diagnostic = .iwmde_cache_get(diagnostic_cache, cache_key),
      parameter  = parameter
    ))
  }

  if (identical(plan[["status"]], "unsupported")) {
    out <- .iwmde_unsupported(parameter, plan[["reason"]])
    .iwmde_cache_set(diagnostic_cache, cache_key, out)
    return(out)
  }
  if (identical(plan[["status"]], "point_only")) {
    if (identical(output, "density")) {
      out <- .iwmde_point_only_diagnostic(
        parameter = parameter,
        samples   = plan[["rows"]][["samples"]],
        component = plan[["rows"]][["component"]]
      )
    } else {
      out <- .iwmde_unsupported(
        parameter,
        "parameter has no continuous active posterior component"
      )
    }
    .iwmde_cache_set(diagnostic_cache, cache_key, out)
    return(out)
  }

  if (identical(output, "ordinate") &&
      length(plan[["outputs"]][["requested_values"]]) == 0L) {
    out <- .iwmde_unsupported(parameter, "no finite ordinate values")
    .iwmde_cache_set(diagnostic_cache, cache_key, out)
    return(out)
  }
  if (identical(output, "ordinate") &&
      length(plan[["grids"]][["requested_values"]]) == 0L) {
    out <- .iwmde_unsupported(parameter, "ordinate values are outside support")
    .iwmde_cache_set(diagnostic_cache, cache_key, out)
    return(out)
  }

  execution <- .iwmde_plan_row_execution(
    context         = context,
    plan            = plan,
    execution_cache = execution_cache
  )
  if (length(execution[["active_rows"]]) < 20L) {
    out <- .iwmde_unsupported(
      parameter,
      "fewer than 20 finite baseline log-q values"
    )
    .iwmde_cache_set(diagnostic_cache, cache_key, out)
    return(out)
  }

  density <- .iwmde_execute_plan_estimator(
    context   = context,
    plan      = plan,
    output    = output,
    execution = execution
  )
  if (identical(density[["status"]], "unsupported")) {
    .iwmde_cache_set(diagnostic_cache, cache_key, density)
    return(density)
  }
  if (density[["n_normalized_rows"]] < 20L) {
    out <- .iwmde_unsupported(
      parameter,
      "fewer than 20 rows had finite conditional normalizers"
    )
    .iwmde_cache_set(diagnostic_cache, cache_key, out)
    return(out)
  }
  if (min(density[["finite_terms"]]) == 0L) {
    out <- .iwmde_unsupported(
      parameter,
      paste0("at least one IWMDE ", .iwmde_plan_output_label(output),
             " had no finite importance terms")
    )
    .iwmde_cache_set(diagnostic_cache, cache_key, out)
    return(out)
  }

  out <- .iwmde_plan_diagnostic_result(
    plan      = plan,
    output    = output,
    execution = execution,
    density   = density
  )
  .iwmde_cache_set(diagnostic_cache, cache_key, out)

  return(out)
}


.iwmde_plan_active_key_counts <- function(row_states) {

  if (length(row_states) == 0L) {
    return(integer())
  }

  keys <- vapply(row_states, function(state) {
    key <- state[["active_key"]]
    if (length(key) != 1L || is.na(key) || !nzchar(key)) {
      return("all")
    }

    return(as.character(key))
  }, character(1))

  return(table(keys))
}


.iwmde_plan_output_cache_key <- function(plan, output) {

  return(.iwmde_hash("iwmde_plan_output", list(
    plan_key = plan[["plan_key"]],
    output   = output
  )))
}


.iwmde_plan_output_label <- function(output) {

  if (identical(output, "ordinate")) {
    return("ordinate")
  }

  return("grid point")
}


.iwmde_plan_row_execution <- function(context, plan, execution_cache = NULL) {

  key <- "row_execution"
  if (!is.null(execution_cache) &&
      exists(key, envir = execution_cache, inherits = FALSE)) {
    return(get(key, envir = execution_cache, inherits = FALSE))
  }

  rows <- plan[["rows"]]

  out <- list(
    active_rows      = rows[["estimator_rows"]],
    active_values    = rows[["estimator_values"]],
    row_states       = rows[["row_states"]],
    baseline_log_q   = rows[["baseline_log_q"]],
    n_dropped_log_q  = rows[["n_dropped_log_q"]],
    n_candidate_rows = rows[["n_denominator_rows"]]
  )
  if (!is.null(execution_cache)) {
    assign(key, out, envir = execution_cache)
  }

  return(out)
}


.iwmde_execute_plan_estimator <- function(context, plan, output, execution) {

  display_grid <- if (identical(output, "density")) {
    plan[["grids"]][["display_grid"]]
  } else {
    plan[["grids"]][["evaluation_values"]]
  }
  if (length(display_grid) == 0L) {
    return(.iwmde_unsupported(
      plan[["target"]][["parameter"]],
      "estimator display grid is empty"
    ))
  }

  if (identical(plan[["method"]], "q_grid_cmde")) {
    if (is.null(plan[["grids"]][["normalization_grid"]])) {
      return(.iwmde_unsupported(
        plan[["target"]][["parameter"]],
        "could not construct a normalization grid"
      ))
    }
    density <- .iwmde_density_grid(
      context            = context,
      parameter          = plan[["target"]][["parameter"]],
      display_grid       = display_grid,
      normalization_grid = plan[["grids"]][["normalization_grid"]],
      transform          = plan[["support"]][["transform"]],
      row_states         = execution[["row_states"]],
      active_mass        = plan[["rows"]][["active_mass"]],
      replacement        = plan[["replacement"]],
      n_candidate_rows   = execution[["n_candidate_rows"]]
    )
  } else {
    density <- .iwmde_density_iwmde(
      context            = context,
      parameter          = plan[["target"]][["parameter"]],
      parameter_spec     = plan[["execution_spec"]],
      display_grid       = display_grid,
      row_states         = execution[["row_states"]],
      active_rows        = execution[["active_rows"]],
      active_values      = execution[["active_values"]],
      weight_rows        = execution[["active_rows"]],
      weight_values      = execution[["active_values"]],
      support            = plan[["support"]][["support"]],
      active_mass        = plan[["rows"]][["active_mass"]],
      replacement        = plan[["replacement"]],
      normalization_grid = plan[["grids"]][["normalization_grid"]],
      n_candidate_rows   = execution[["n_candidate_rows"]]
    )
  }

  if (identical(output, "ordinate") && is.list(density)) {
    density[["evaluation_x"]] <- plan[["grids"]][["evaluation_values"]]
    density[["x"]]            <- plan[["grids"]][["requested_values"]]
  }

  return(.iwmde_new_density_result(
    fields    = density,
    estimator = plan[["method"]]
  ))
}


.iwmde_plan_diagnostic_result <- function(plan, output, execution, density) {

  parameter <- plan[["target"]][["parameter"]]
  rows      <- plan[["rows"]]
  is_density <- identical(output, "density")
  .iwmde_validate_density_result(
    density   = density,
    estimator = plan[["method"]]
  )

  plot_integral <- if (is_density) {
    .iwmde_trapz(density[["x"]], density[["y"]])
  } else {
    NA_real_
  }
  kde <- if (is_density) {
    .iwmde_kde(
      values   = rows[["continuous_values"]],
      xlim     = plan[["support"]][["density_xlim"]],
      n_points = plan[["control"]][["n_points"]],
      mass     = rows[["active_mass"]]
    )
  } else {
    NULL
  }
  hist_data <- if (is_density) {
    .iwmde_histogram(
      values = rows[["continuous_values"]],
      xlim   = plan[["support"]][["density_xlim"]],
      mass   = rows[["active_mass"]]
    )
  } else {
    NULL
  }

  bf_values <- if (is_density) {
    plan[["outputs"]][["requested_values"]]
  } else {
    plan[["grids"]][["requested_values"]]
  }
  bf_diagnostics <- .iwmde_density_bf_diagnostics(density, bf_values)
  active_key_counts <- .iwmde_plan_active_key_counts(execution[["row_states"]])

  diagnostics <- list(
    integral                    = plot_integral,
    plot_integral               = plot_integral,
    point_mass_total            = rows[["point_mass_total"]],
    plot_total_mass             = if (is_density) {
      plot_integral + rows[["point_mass_total"]]
    } else {
      NA_real_
    },
    display_mass_fraction       = if (is_density) {
      plot_integral / rows[["active_mass"]]
    } else {
      NA_real_
    },
    pilot_normalization_integral =
      density[["pilot_normalization_integral"]],
    final_normalization_integral =
      density[["final_normalization_integral"]],
    support_grid_normalization_integral =
      density[["support_grid_normalization_integral"]],
    normalization_relative_error =
      density[["normalization_relative_error"]],
    normalization_points        = density[["normalization_points"]],
    normalization_range         = density[["normalization_range"]],
    normalization_initial_points = density[["normalization_initial_points"]],
    normalization_initial_range = density[["normalization_initial_range"]],
    normalization_scale         = density[["normalization_scale"]],
    normalization_mass_ratio    = density[["normalization_mass_ratio"]],
    max_ordinate_relative_change =
      .iwmde_max_or_na(density[["ordinate_relative_change"]]),
    max_normalizer_relative_change =
      density[["max_normalizer_relative_change"]],
    max_quadrature_relative_change =
      density[["max_quadrature_relative_change"]],
    p95_normalizer_relative_change =
      density[["p95_normalizer_relative_change"]],
    median_normalizer_relative_change =
      density[["median_normalizer_relative_change"]],
    normalization_refined_points =
      density[["normalization_refined_points"]],
    normalization_refined_range =
      density[["normalization_refined_range"]],
    n_rescued_normalizer        = density[["n_rescued_normalizer"]],
    n_dropped_normalizer        = density[["n_dropped_normalizer"]],
    n_initial_dropped_normalizer =
      density[["n_initial_dropped_normalizer"]],
    n_validation_dropped_normalizer =
      density[["n_validation_dropped_normalizer"]],
    n_refinement_steps          = density[["n_refinement_steps"]],
    active_mass                 = rows[["active_mass"]],
    n_candidate_rows            = rows[["n_candidate_rows"]],
    n_denominator_rows          = rows[["n_denominator_rows"]],
    n_estimator_rows            = rows[["n_estimator_rows"]],
    n_evaluated_rows            = density[["n_evaluated_rows"]],
    n_dropped_rows              = density[["n_dropped_rows"]],
    row_drop_fraction           = density[["row_drop_fraction"]],
    n_active                    = length(execution[["active_rows"]]),
    n_active_state_keys         = length(active_key_counts),
    min_active_state_rows       = if (length(active_key_counts) == 0L) {
      NA_integer_
    } else {
      min(active_key_counts)
    },
    n_total                     = rows[["n_total"]],
    n_dropped_log_q             = execution[["n_dropped_log_q"]],
    n_dropped_weight            = density[["n_dropped_weight"]],
    weight_partitions           = density[["weight_partitions"]],
    n_weight_fallbacks          = if (is.null(density[["n_weight_fallbacks"]])) {
      0L
    } else {
      density[["n_weight_fallbacks"]]
    },
    n_weight_fallback_rows      = if (is.null(density[["n_weight_fallback_rows"]])) {
      0L
    } else {
      density[["n_weight_fallback_rows"]]
    },
    weight_fallback_from        = if (is.null(density[["weight_fallback_from"]])) {
      character()
    } else {
      density[["weight_fallback_from"]]
    },
    weight_fallback_reasons     = if (is.null(density[["weight_fallback_reasons"]])) {
      structure(integer(), names = character())
    } else {
      density[["weight_fallback_reasons"]]
    },
    max_log_ratio               = density[["max_log_ratio"]],
    min_finite_terms            = min(density[["finite_terms"]]),
    n_normalized_rows           = density[["n_normalized_rows"]],
    min_ess                     = min(density[["ess"]]),
    max_weight_share            = max(density[["max_weight_share"]]),
    max_mcse                    = .iwmde_max_or_na(density[["mcse"]]),
    max_relative_mcse           =
      .iwmde_max_or_na(density[["relative_mcse"]]),
    plot_integral_mcse          = density[["integral_mcse"]],
    plot_integral_relative_mcse = density[["integral_relative_mcse"]],
    batch_size                  = density[["batch_size"]],
    n_batches                   = density[["n_batches"]],
    estimator                   = density[["estimator"]],
    weight_method               = density[["weight_method"]],
    display_grid                = if (is_density) {
      plan[["control"]][["display_grid"]]
    } else {
      "ordinate"
    },
    bf_value                    = bf_diagnostics[["bf_value"]],
    bf_evaluation_value         = bf_diagnostics[["bf_evaluation_value"]],
    bf_included                 = bf_diagnostics[["bf_included"]],
    bf_grid_index               = bf_diagnostics[["bf_grid_index"]],
    bf_ordinate                 = bf_diagnostics[["bf_ordinate"]],
    bf_pilot_ordinate           = bf_diagnostics[["bf_pilot_ordinate"]],
    bf_validation_ordinate      =
      bf_diagnostics[["bf_validation_ordinate"]],
    bf_ordinate_relative_change =
      bf_diagnostics[["bf_ordinate_relative_change"]],
    bf_ordinate_log_change      =
      bf_diagnostics[["bf_ordinate_log_change"]],
    bf_pilot_ordinate_relative_change =
      bf_diagnostics[["bf_pilot_ordinate_relative_change"]],
    bf_pilot_ordinate_log_change =
      bf_diagnostics[["bf_pilot_ordinate_log_change"]],
    bf_mcse                     = bf_diagnostics[["bf_mcse"]],
    bf_relative_mcse            = bf_diagnostics[["bf_relative_mcse"]],
    bf_error_percent            = bf_diagnostics[["bf_error_percent"]],
    bf_finite_terms             = bf_diagnostics[["bf_finite_terms"]],
    bf_ess                      = bf_diagnostics[["bf_ess"]],
    bf_max_weight_share         = bf_diagnostics[["bf_max_weight_share"]],
    bf_max_log_ratio            = bf_diagnostics[["bf_max_log_ratio"]]
  )

  out <- .iwmde_new_diagnostic(list(
    parameter    = parameter,
    status       = "ok",
    samples      = rows[["samples"]],
    target_key   = plan[["target"]][["target_key"]],
    active_rows  = execution[["active_rows"]],
    active_mass  = rows[["active_mass"]],
    point_masses = rows[["point_masses"]],
    support      = plan[["support"]][["support"]],
    xlim         = if (is_density) {
      plan[["support"]][["density_xlim"]]
    } else {
      plan[["support"]][["xlim"]]
    },
    iwmde        = density,
    diagnostics  = diagnostics
  ))
  if (is_density) {
    out[["histogram"]] <- hist_data
    out[["kde"]]       <- kde
  }
  return(out)
}
