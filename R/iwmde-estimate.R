# ============================================================================ #
# IWMDE Estimate Facade
# ============================================================================ #

.iwmde_estimate <- function(context, parameter, density_method, density_control,
                            outputs = c("density", "ordinate"), values = NULL,
                            parameter_spec = NULL, metadata = NULL,
                            cache = NULL) {

  context <- .iwmde_context_ensure_caches(context)
  outputs <- unique(match.arg(
    outputs,
    c("density", "ordinate"),
    several.ok = TRUE
  ))
  ordinate_values <- NULL
  if ("ordinate" %in% outputs) {
    ordinate_values <- .iwmde_sorted_ordinate_values(values)
  }
  parameter_spec <- .iwmde_prepare_prior_ordinates(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    values         = ordinate_values
  )
  ordinate_warnings <- .iwmde_ordinate_prior_warnings(
    parameter       = parameter,
    prior_ordinates = parameter_spec[["prior_ordinates"]]
  )
  if (length(ordinate_warnings) > 0L) {
    warning(paste(ordinate_warnings, collapse = " "), call. = FALSE)
  }
  .iwmde_check_context_density_method_supported(context, density_method)
  if ("ordinate" %in% outputs) {
    out <- .iwmde_estimate_fixed_ordinate(
      context         = context,
      parameter       = parameter,
      density_method  = density_method,
      density_control = density_control,
      outputs         = outputs,
      values          = values,
      parameter_spec  = parameter_spec,
      metadata        = metadata,
      cache           = cache
    )
  } else {
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
    out <- .iwmde_estimate_from_plan(
      context = context,
      plan    = plan,
      cache   = cache
    )
  }

  return(out)
}


.iwmde_estimate_fixed_ordinate <- function(
    context, parameter, density_method, density_control, values,
    parameter_spec = NULL, metadata = NULL, cache = NULL,
    outputs = "ordinate") {

  control <- .iwmde_density_control_resolve(
    density_method  = density_method,
    density_control = density_control,
    purpose         = "ordinate"
  )
  plan <- .iwmde_plan(
    context         = context,
    parameter       = parameter,
    density_method  = density_method,
    density_control = control,
    outputs         = outputs,
    values          = values,
    parameter_spec  = parameter_spec,
    metadata        = metadata
  )
  out <- .iwmde_estimate_from_plan(
    context = context,
    plan    = plan,
    cache   = cache
  )

  diagnostic <- out[["diagnostics"]][["ordinate"]]
  eligible   <- length(plan[["rows"]][["continuous_rows"]])
  achieved   <- plan[["rows"]][["n_candidate_rows"]]
  if (!is.finite(eligible) || eligible < 1L) {
    eligible <- achieved
  }
  ordinate_entries <- c(
    .iwmde_posterior_ordinate_entries(out[["posterior_ordinate"]]),
    .iwmde_posterior_ordinate_entries(out[["rejected_posterior_ordinate"]])
  )
  if (length(ordinate_entries) == 0L) {
    ordinate_entries <- list(diagnostic)
  }
  ordinate_designs <- lapply(ordinate_entries, function(entry) {
    .iwmde_ordinate_sampling_design(
      diagnostic = entry,
      control    = control,
      eligible   = eligible,
      achieved   = achieved
    )
  })
  names(ordinate_designs) <- vapply(ordinate_entries, function(entry) {
    value <- .iwmde_ordinate_scalar(entry, "value")
    if (is.finite(value)) .iwmde_key_number(value) else "diagnostic"
  }, character(1))
  sampling_design <- .iwmde_ordinate_sampling_design_aggregate(
    ordinate_designs
  )
  out <- .iwmde_estimate_attach_sampling_design(
    estimate         = out,
    sampling_design  = sampling_design,
    ordinate_designs = ordinate_designs
  )

  return(out)
}


.iwmde_ordinate_sampling_design <- function(diagnostic, control,
                                             eligible, achieved) {

  metrics <- .iwmde_ordinate_precision_metrics(diagnostic)
  precision_target_met <- is.finite(metrics[["relative_mcse"]]) &&
    metrics[["relative_mcse"]] <= control[["target_relative_mcse"]]
  sampling_target_met <- is.finite(metrics[["sampling_relative_mcse"]]) &&
    metrics[["sampling_relative_mcse"]] <=
      control[["target_relative_mcse"]]
  numerical_grade_met <- if (!is.null(diagnostic[["ordinate"]])) {
    is.null(.iwmde_posterior_ordinate_bf_failure_reason(diagnostic))
  } else {
    is.null(.iwmde_diagnostics_bf_failure_reason(
      diagnostic[["diagnostics"]]
    ))
  }
  bf_grade_met <- numerical_grade_met && !identical(
    diagnostic[["diagnostics"]][["mixture_mcse_type"]],
    "worst_correlation_delta_upper_bound"
  )

  return(list(
    fixed_budget         = TRUE,
    requested_samples    = control[["samples"]],
    achieved_row_budget  = as.integer(achieved),
    eligible_rows        = as.integer(eligible),
    all_rows_used        = is.finite(eligible) && achieved >= eligible,
    target_relative_mcse = control[["target_relative_mcse"]],
    precision_target_met = precision_target_met,
    sampling_target_met  = sampling_target_met,
    bf_grade_met         = bf_grade_met,
    target_met           = precision_target_met && sampling_target_met &&
      bf_grade_met
  ))
}


.iwmde_ordinate_sampling_design_aggregate <- function(designs) {

  out <- designs[[1L]]
  fields <- c(
    "all_rows_used", "precision_target_met", "sampling_target_met",
    "bf_grade_met", "target_met"
  )
  for (field in fields) {
    out[[field]] <- all(vapply(designs, `[[`, logical(1), field))
  }

  return(out)
}


.iwmde_ordinate_precision_metrics <- function(diagnostic) {

  diagnostics <- diagnostic[["diagnostics"]]
  if (!is.list(diagnostics)) {
    diagnostics <- list()
  }

  return(list(
    relative_mcse = .iwmde_diagnostic_scalar_any(
      diagnostics,
      c("bf_relative_mcse", "relative_mcse")
    ),
    sampling_relative_mcse = .iwmde_diagnostic_scalar_any(
      diagnostics,
      c("bf_sampling_relative_mcse", "max_sampling_relative_mcse",
        "sampling_relative_mcse")
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


.iwmde_estimate_attach_sampling_design <- function(
    estimate, sampling_design, ordinate_designs = list()) {

  diagnostic <- estimate[["diagnostics"]][["ordinate"]]
  if (is.list(diagnostic)) {
    diagnostic[["sampling_design"]] <- sampling_design
    if (is.list(diagnostic[["diagnostics"]])) {
      for (name in names(sampling_design)) {
        diagnostic[["diagnostics"]][[name]] <- sampling_design[[name]]
      }
    }
    estimate[["diagnostics"]][["ordinate"]] <- diagnostic
  }

  for (ordinate_name in c(
    "posterior_ordinate",
    "rejected_posterior_ordinate"
  )) {
    ordinate <- estimate[[ordinate_name]]
    entries  <- .iwmde_posterior_ordinate_entries(ordinate)
    for (i in seq_along(entries)) {
      value <- .iwmde_ordinate_scalar(entries[[i]], "value")
      key   <- if (is.finite(value)) .iwmde_key_number(value) else NULL
      design <- if (!is.null(key) && key %in% names(ordinate_designs)) {
        ordinate_designs[[key]]
      } else {
        sampling_design
      }
      for (name in names(design)) {
        entries[[i]][["diagnostics"]][[name]] <- design[[name]]
      }
      warnings <- .iwmde_diagnostics_bf_warning(
        entries[[i]][["diagnostics"]]
      )
      if (length(warnings) == 0L) {
        entries[[i]][["diagnostics"]][["warning"]] <- NULL
      } else {
        entries[[i]][["diagnostics"]][["warning"]] <- warnings
      }
    }
    estimate[[ordinate_name]] <- .iwmde_posterior_ordinate_combine(entries)
  }

  estimate[["sampling_design"]] <- sampling_design

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
    if (length(plan[["grids"]][["requested_values"]]) <= 1L) {
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
    } else {
      candidate_ordinates <- .iwmde_posterior_ordinate_attributes(
        diagnostic      = ordinate_diagnostic,
        density_method  = plan[["density_method"]],
        density_control = plan[["control"]],
        metadata        = plan[["target"]][["metadata"]]
      )
      ordinate_attribute          <- candidate_ordinates[["accepted"]]
      rejected_ordinate_attribute <- candidate_ordinates[["rejected"]]
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
  context_chain_id <- context[["chain_id"]]
  if (is.null(context_chain_id)) {
    context_chain_id <- rep(1L, nrow(context[["posterior_samples"]]))
  }

  out <- list(
    active_rows      = rows[["estimator_rows"]],
    active_values    = rows[["estimator_values"]],
    population_rows  = rows[["population_rows"]],
    sampling_population_rows = rows[["continuous_rows"]],
    conditioned_chain_id = context_chain_id[rows[["population_rows"]]],
    chain_id         = context_chain_id[rows[["estimator_rows"]]],
    expected_chain_ids = rows[["chain_coverage"]][["expected_chain_ids"]],
    row_states       = rows[["row_states"]],
    baseline_log_q   = rows[["baseline_log_q"]],
    n_candidate_rows = rows[["n_denominator_rows"]]
  )
  if (identical(plan[["method"]], "iwmde")) {
    out[["proposal_weight"]] <- .iwmde_plan_proposal_weight(
      context         = context,
      plan            = plan,
      active_rows     = out[["active_rows"]],
      active_values   = out[["active_values"]],
      proposal_rows   = rows[["continuous_rows"]],
      proposal_values = rows[["continuous_values"]]
    )
  }
  if (!is.null(execution_cache)) {
    assign(key, out, envir = execution_cache)
  }

  return(out)
}


.iwmde_plan_proposal_weight <- function(context, plan, active_rows,
                                         active_values, proposal_rows,
                                         proposal_values) {

  tryCatch(
    .iwmde_chen_log_weight(
      context        = context,
      parameter      = plan[["target"]][["parameter"]],
      parameter_spec = plan[["execution_spec"]],
      active_rows    = active_rows,
      active_values  = active_values,
      weight_rows    = proposal_rows,
      weight_values  = proposal_values,
      support        = plan[["support"]][["support"]]
    ),
    error = function(e) {
      if (inherits(e, "iwmde_construction_error")) {
        stop(e)
      }
      .iwmde_stop_construction_failure(
        estimator = "iwmde",
        parameter = plan[["target"]][["parameter"]],
        rows      = proposal_rows,
        stage     = "proposal-density construction",
        detail    = conditionMessage(e)
      )
    }
  )
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

  mixed_output <- plan[["rows"]][["point_mass_total"]] > 0
  conditioned_rows <- if (mixed_output) {
    execution[["population_rows"]]
  } else {
    NULL
  }
  conditioned_chain_id <- if (mixed_output) {
    execution[["conditioned_chain_id"]]
  } else {
    NULL
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
      estimator_rows     = execution[["active_rows"]],
      population_rows    = execution[["sampling_population_rows"]],
      chain_id           = execution[["chain_id"]],
      expected_chain_ids = execution[["expected_chain_ids"]],
      conditioned_rows   = conditioned_rows,
      conditioned_chain_id = conditioned_chain_id,
      active_mass        = plan[["rows"]][["active_mass"]],
      replacement        = plan[["replacement"]],
      n_candidate_rows   = execution[["n_candidate_rows"]]
    )
  } else {
    density <- .iwmde_density_iwmde(
      context            = context,
      parameter          = plan[["target"]][["parameter"]],
      display_grid       = display_grid,
      row_states         = execution[["row_states"]],
      active_rows        = execution[["active_rows"]],
      active_values      = execution[["active_values"]],
      population_rows    = execution[["sampling_population_rows"]],
      chain_id           = execution[["chain_id"]],
      expected_chain_ids = execution[["expected_chain_ids"]],
      conditioned_rows   = conditioned_rows,
      conditioned_chain_id = conditioned_chain_id,
      proposal_weight    = execution[["proposal_weight"]],
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
  curve_diagnostics <- .iwmde_density_curve_diagnostics(
    density            = density,
    tail_probabilities = plan[["grids"]][["tail_probabilities"]],
    tail_values        = plan[["grids"]][["tail_values"]]
  )

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
    n_initial_dropped_normalizer =
      density[["n_initial_dropped_normalizer"]],
    n_refinement_steps          = density[["n_refinement_steps"]],
    active_mass                 = rows[["active_mass"]],
    n_candidate_rows            = rows[["n_candidate_rows"]],
    n_denominator_rows          = rows[["n_denominator_rows"]],
    n_estimator_rows            = rows[["n_estimator_rows"]],
    n_evaluated_rows            = density[["n_evaluated_rows"]],
    n_active                    = length(execution[["active_rows"]]),
    n_active_state_keys         = length(active_key_counts),
    min_active_state_rows       = if (length(active_key_counts) == 0L) {
      NA_integer_
    } else {
      min(active_key_counts)
    },
    n_total                     = rows[["n_total"]],
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
    min_ess                     = min(density[["ess"]]),
    max_weight_share            = max(density[["max_weight_share"]]),
    max_mcse                    = .iwmde_max_or_na(density[["mcse"]]),
    max_relative_mcse           =
      .iwmde_max_or_na(density[["relative_mcse"]]),
    max_active_branch_mcse      =
      .iwmde_max_or_na(density[["active_branch_mcse"]]),
    max_active_branch_relative_mcse =
      .iwmde_max_or_na(density[["active_branch_relative_mcse"]]),
    active_mass_mcse            = density[["active_mass_mcse"]],
    active_mass_relative_mcse   = density[["active_mass_relative_mcse"]],
    max_active_mass_component_mcse =
      .iwmde_max_or_na(density[["active_mass_component_mcse"]]),
    mixture_mcse_type           = density[["mixture_mcse_type"]],
    max_sampling_mcse           =
      .iwmde_max_or_na(density[["sampling_mcse"]]),
    max_sampling_relative_mcse  =
      .iwmde_max_or_na(density[["sampling_relative_mcse"]]),
    sampling_fraction           = density[["sampling_fraction"]],
    sampling_uncertainty_type   = density[["sampling_uncertainty_type"]],
    mcmc_uncertainty_scope      = density[["mcmc_uncertainty_scope"]],
    mcmc_uncertainty_status     = density[["mcmc_uncertainty_status"]],
    mcmc_uncertainty_reason     = density[["mcmc_uncertainty_reason"]],
    plot_scale_relative_mcse    =
      curve_diagnostics[["plot_scale_relative_mcse"]],
    plot_scale_sampling_relative_mcse =
      curve_diagnostics[["plot_scale_sampling_relative_mcse"]],
    bulk_probability_range      =
      curve_diagnostics[["bulk_probability_range"]],
    bulk_x_range                = curve_diagnostics[["bulk_x_range"]],
    bulk_max_relative_mcse      =
      curve_diagnostics[["bulk_max_relative_mcse"]],
    bulk_max_sampling_relative_mcse =
      curve_diagnostics[["bulk_max_sampling_relative_mcse"]],
    bulk_min_ess                = curve_diagnostics[["bulk_min_ess"]],
    bulk_max_weight_share       =
      curve_diagnostics[["bulk_max_weight_share"]],
    tail_probabilities          =
      curve_diagnostics[["tail_probabilities"]],
    tail_target_x               = curve_diagnostics[["tail_target_x"]],
    tail_evaluation_x           =
      curve_diagnostics[["tail_evaluation_x"]],
    tail_relative_mcse          =
      curve_diagnostics[["tail_relative_mcse"]],
    tail_ess                    = curve_diagnostics[["tail_ess"]],
    tail_max_weight_share       =
      curve_diagnostics[["tail_max_weight_share"]],
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
    prior_ordinates             = plan[["prior_ordinates"]],
    ordinate_warnings           = plan[["ordinate_warnings"]],
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
    bf_active_branch_mcse       =
      bf_diagnostics[["bf_active_branch_mcse"]],
    bf_active_branch_relative_mcse =
      bf_diagnostics[["bf_active_branch_relative_mcse"]],
    bf_active_mass_component_mcse =
      bf_diagnostics[["bf_active_mass_component_mcse"]],
    bf_sampling_mcse            = bf_diagnostics[["bf_sampling_mcse"]],
    bf_sampling_relative_mcse   =
      bf_diagnostics[["bf_sampling_relative_mcse"]],
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


.iwmde_density_curve_diagnostics <- function(density, tail_probabilities,
                                              tail_values) {

  empty <- list(
    plot_scale_relative_mcse = NA_real_,
    plot_scale_sampling_relative_mcse = NA_real_,
    bulk_probability_range   = c(NA_real_, NA_real_),
    bulk_x_range             = c(NA_real_, NA_real_),
    bulk_max_relative_mcse   = NA_real_,
    bulk_max_sampling_relative_mcse = NA_real_,
    bulk_min_ess             = NA_real_,
    bulk_max_weight_share    = NA_real_,
    tail_probabilities       = c(NA_real_, NA_real_),
    tail_target_x            = c(NA_real_, NA_real_),
    tail_evaluation_x        = c(NA_real_, NA_real_),
    tail_relative_mcse       = c(NA_real_, NA_real_),
    tail_ess                 = c(NA_real_, NA_real_),
    tail_max_weight_share    = c(NA_real_, NA_real_)
  )

  tail_probabilities <- as.numeric(tail_probabilities)
  tail_values        <- as.numeric(tail_values)
  if (length(tail_probabilities) != 2L ||
      length(tail_values) != 2L ||
      any(!is.finite(tail_probabilities)) ||
      any(!is.finite(tail_values)) ||
      tail_probabilities[1L] >= tail_probabilities[2L] ||
      tail_values[1L] > tail_values[2L]) {
    return(empty)
  }

  x                <- density[["x"]]
  y                <- density[["y"]]
  mcse             <- density[["mcse"]]
  relative_mcse    <- density[["relative_mcse"]]
  sampling_mcse    <- density[["sampling_mcse"]]
  sampling_relative_mcse <- density[["sampling_relative_mcse"]]
  ess              <- density[["ess"]]
  max_weight_share <- density[["max_weight_share"]]
  peak_density     <- .iwmde_max_or_na(y)
  if (length(x) == 0L || !is.finite(peak_density) || peak_density <= 0) {
    return(empty)
  }

  bulk <- x >= tail_values[1L] & x <= tail_values[2L]
  if (!any(bulk)) {
    return(empty)
  }
  tail_index <- vapply(tail_values, function(value) {
    which.min(abs(x - value))
  }, integer(1))

  bulk_ess <- ess[bulk]
  bulk_ess <- bulk_ess[is.finite(bulk_ess)]

  return(list(
    plot_scale_relative_mcse = .iwmde_max_or_na(mcse) / peak_density,
    plot_scale_sampling_relative_mcse =
      .iwmde_max_or_na(sampling_mcse) / peak_density,
    bulk_probability_range   = range(tail_probabilities),
    bulk_x_range             = range(tail_values),
    bulk_max_relative_mcse   = .iwmde_max_or_na(relative_mcse[bulk]),
    bulk_max_sampling_relative_mcse =
      .iwmde_max_or_na(sampling_relative_mcse[bulk]),
    bulk_min_ess             = if (length(bulk_ess) == 0L) {
      NA_real_
    } else {
      min(bulk_ess)
    },
    bulk_max_weight_share    = .iwmde_max_or_na(max_weight_share[bulk]),
    tail_probabilities       = tail_probabilities,
    tail_target_x            = tail_values,
    tail_evaluation_x        = x[tail_index],
    tail_relative_mcse       = relative_mcse[tail_index],
    tail_ess                 = ess[tail_index],
    tail_max_weight_share    = max_weight_share[tail_index]
  ))
}
