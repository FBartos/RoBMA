# ============================================================================ #
# IWMDE Density Methods, Metadata, and Attribute Adapters
# ============================================================================ #

.density_method_normalize <- function(density_method, allow_normal = FALSE) {

  allowed <- c("KDE", "qCMDE", "IWMDE")
  if (isTRUE(allow_normal)) {
    allowed <- c("KDE", "normal", "qCMDE", "IWMDE")
  }

  return(.density_method_match(density_method, allowed))
}


.density_method_normalize_precomputed <- function(density_method) {

  density_method <- .density_method_match(density_method, c("qCMDE", "IWMDE"))

  return(density_method)
}


.density_method_match <- function(density_method, allowed) {

  if (is.null(density_method)) {
    density_method <- allowed[[1L]]
  }

  if (is.character(density_method) &&
      length(density_method) > 0L &&
      !anyNA(density_method)) {
    density_method_lower <- tolower(density_method)
    allowed_lower        <- tolower(allowed)

    if (length(density_method) > 1L) {
      if (identical(density_method_lower, allowed_lower)) {
        density_method <- allowed
      }
    } else {
      index <- pmatch(density_method_lower, allowed_lower, duplicates.ok = TRUE)
      if (!is.na(index)) {
        density_method <- allowed[[index]]
      }
    }
  }

  return(BayesTools::posterior_density_method_match(
    method  = density_method,
    allowed = allowed,
    name    = "density_method"
  ))
}


.iwmde_posterior_metadata <- function(samples, parameter, level = NULL) {

  metadata <- c(
    list(
      parameter = parameter,
      level     = level
    ),
    .iwmde_sample_condition_metadata(samples)
  )
  metadata <- metadata[!vapply(metadata, is.null, logical(1))]

  return(metadata)
}


.iwmde_sample_condition_metadata <- function(samples,
                                            include_prior_density_context = FALSE,
                                            require_child_condition = FALSE) {

  conditional <- .iwmde_first_nonempty_attr(
    samples,
    c("effective_conditional", "conditional")
  )
  conditional_rule <- .iwmde_first_nonempty_attr(
    samples,
    c("effective_conditional_rule", "conditional_rule")
  )
  condition_key <- attr(samples, "condition_key", exact = TRUE)

  if (isTRUE(require_child_condition) &&
      (is.null(conditional) || is.null(conditional_rule) || is.null(condition_key))) {
    stop(
      "Conditional marginal-means qCMDE/IWMDE metadata is incomplete; recompute the marginal_means object with current BayesTools.",
      call. = FALSE
    )
  }

  metadata <- list(
    conditional              = conditional,
    conditional_rule         = conditional_rule,
    condition_key            = condition_key,
    condition_event          = attr(samples, "condition_event", exact = TRUE),
    resolved_condition_event = attr(samples, "resolved_condition_event", exact = TRUE)
  )
  if (isTRUE(include_prior_density_context)) {
    metadata[["prior_density_context"]] <- attr(
      samples,
      "prior_density_context",
      exact = TRUE
    )
  }
  metadata <- metadata[!vapply(metadata, is.null, logical(1))]

  return(metadata)
}


.iwmde_first_nonempty_attr <- function(x, names) {

  for (name in names) {
    value <- attr(x, name, exact = TRUE)
    if (!is.null(value) && length(value) > 0L) {
      return(value)
    }
  }

  return(NULL)
}


.density_method_uses_precomputed <- function(density_method, allow_normal = FALSE) {

  density_method <- .density_method_normalize(
    density_method = density_method,
    allow_normal   = allow_normal
  )

  return(BayesTools::posterior_density_method_uses_precomputed(density_method))
}


.density_method_iwmde_estimator <- function(density_method) {

  density_method <- .density_method_normalize(density_method)
  if (identical(density_method, "KDE")) {
    return(NULL)
  }
  if (identical(density_method, "qCMDE")) {
    return("q_grid_cmde")
  }

  return("iwmde")
}


.iwmde_posterior_density_attribute <- function(diagnostic, density_method,
                                               metadata = NULL,
                                               density_control = NULL) {

  if (!identical(diagnostic[["status"]], "ok") ||
      is.null(diagnostic[["iwmde"]])) {
    return(NULL)
  }

  density <- diagnostic[["iwmde"]]
  diagnostics <- diagnostic[["diagnostics"]]
  density_failure <- .iwmde_diagnostics_density_failure_reason(diagnostics)
  if (!is.null(density_failure)) {
    return(NULL)
  }
  density_warning <- .iwmde_diagnostics_density_warning(diagnostics)
  if (length(density_warning) > 0L) {
    diagnostics[["warning"]] <- density_warning
  }
  provenance <- .iwmde_result_provenance(
    diagnostic      = diagnostic,
    density_method  = density_method,
    metadata        = metadata,
    density_control = density_control,
    attribute       = "density"
  )
  args <- c(
    list(
      x              = density[["x"]],
      y              = density[["y"]],
      method         = density[["estimator"]],
      density_method = .density_method_normalize(density_method),
      diagnostics    = diagnostics,
      point_masses   = diagnostic[["point_masses"]],
      iwmde_provenance = provenance
    ),
    metadata
  )
  out <- do.call(
    BayesTools::posterior_density_attribute,
    args
  )

  return(out)
}


.iwmde_posterior_ordinate_attribute <- function(diagnostic, density_method,
                                                metadata = NULL,
                                                density_control = NULL,
                                                allow_rejected = FALSE) {

  if (!identical(diagnostic[["status"]], "ok") ||
      is.null(diagnostic[["diagnostics"]])) {
    return(NULL)
  }

  diagnostics <- diagnostic[["diagnostics"]]
  if (!isTRUE(diagnostics[["bf_included"]]) ||
      !is.finite(diagnostics[["bf_ordinate"]]) ||
      diagnostics[["bf_ordinate"]] <= 0) {
    return(NULL)
  }

  ordinate_diagnostics <- list(
    evaluation_value                  = diagnostics[["bf_evaluation_value"]],
    mcse                              = diagnostics[["bf_mcse"]],
    relative_mcse                     = diagnostics[["bf_relative_mcse"]],
    sampling_mcse                     = diagnostics[["bf_sampling_mcse"]],
    sampling_relative_mcse            =
      diagnostics[["bf_sampling_relative_mcse"]],
    sampling_fraction                 = diagnostics[["sampling_fraction"]],
    sampling_uncertainty_type         =
      diagnostics[["sampling_uncertainty_type"]],
    mcmc_uncertainty_scope            =
      diagnostics[["mcmc_uncertainty_scope"]],
    BF_error_percent                  = diagnostics[["bf_error_percent"]],
    finite_terms                      = diagnostics[["bf_finite_terms"]],
    ess                               = diagnostics[["bf_ess"]],
    max_weight_share                  = diagnostics[["bf_max_weight_share"]],
    max_log_ratio                     = diagnostics[["bf_max_log_ratio"]],
    active_mass                       = diagnostics[["active_mass"]],
    n_candidate_rows                  = diagnostics[["n_candidate_rows"]],
    n_denominator_rows                = diagnostics[["n_denominator_rows"]],
    n_estimator_rows                  = diagnostics[["n_estimator_rows"]],
    n_evaluated_rows                  = diagnostics[["n_evaluated_rows"]],
    n_normalized_rows                 = diagnostics[["n_normalized_rows"]],
    n_dropped_rows                    = diagnostics[["n_dropped_rows"]],
    row_drop_fraction                 = diagnostics[["row_drop_fraction"]],
    n_dropped_log_q                   = diagnostics[["n_dropped_log_q"]],
    n_dropped_weight                  = diagnostics[["n_dropped_weight"]],
    weight_partitions                 = diagnostics[["weight_partitions"]],
    n_weight_fallbacks                = diagnostics[["n_weight_fallbacks"]],
    n_weight_fallback_rows            = diagnostics[["n_weight_fallback_rows"]],
    weight_fallback_from              = diagnostics[["weight_fallback_from"]],
    weight_fallback_reasons           = diagnostics[["weight_fallback_reasons"]],
    n_dropped_normalizer              = diagnostics[["n_dropped_normalizer"]],
    n_validation_dropped_normalizer   = diagnostics[["n_validation_dropped_normalizer"]],
    pilot_ordinate                    = diagnostics[["bf_pilot_ordinate"]],
    validation_ordinate               = diagnostics[["bf_validation_ordinate"]],
    ordinate_relative_change          = diagnostics[["bf_ordinate_relative_change"]],
    ordinate_log_change               = diagnostics[["bf_ordinate_log_change"]],
    pilot_ordinate_relative_change    =
      diagnostics[["bf_pilot_ordinate_relative_change"]],
    pilot_ordinate_log_change         =
      diagnostics[["bf_pilot_ordinate_log_change"]],
    pilot_normalization_integral      =
      diagnostics[["pilot_normalization_integral"]],
    final_normalization_integral      =
      diagnostics[["final_normalization_integral"]],
    support_grid_normalization_integral =
      diagnostics[["support_grid_normalization_integral"]],
    normalization_relative_error      =
      diagnostics[["normalization_relative_error"]],
    normalization_mass_ratio          = diagnostics[["normalization_mass_ratio"]],
    normalization_range               = diagnostics[["normalization_range"]],
    normalization_initial_points      = diagnostics[["normalization_initial_points"]],
    normalization_initial_range       = diagnostics[["normalization_initial_range"]],
    max_normalizer_relative_change    = diagnostics[["max_normalizer_relative_change"]],
    p95_normalizer_relative_change    = diagnostics[["p95_normalizer_relative_change"]],
    median_normalizer_relative_change = diagnostics[["median_normalizer_relative_change"]],
    max_quadrature_relative_change    =
      diagnostics[["max_quadrature_relative_change"]],
    normalization_refined_points      = diagnostics[["normalization_refined_points"]],
    normalization_refined_range       = diagnostics[["normalization_refined_range"]],
    n_refinement_steps                = diagnostics[["n_refinement_steps"]],
    estimator                         = diagnostics[["estimator"]],
    weight_method                     = diagnostics[["weight_method"]],
    prior_ordinates                   = diagnostics[["prior_ordinates"]],
    ordinate_warnings                 = diagnostics[["ordinate_warnings"]]
  )
  ordinate_diagnostics <- .iwmde_compact_nulls(ordinate_diagnostics)
  ordinate_warning <- .iwmde_diagnostics_bf_warning(ordinate_diagnostics)
  if (length(ordinate_warning) > 0L) {
    ordinate_diagnostics[["warning"]] <- ordinate_warning
  }
  provenance <- .iwmde_result_provenance(
    diagnostic       = diagnostic,
    density_method   = density_method,
    metadata         = metadata,
    density_control  = density_control,
    value            = diagnostics[["bf_value"]],
    evaluation_value = diagnostics[["bf_evaluation_value"]],
    attribute        = "ordinate"
  )

  args <- c(
    list(
      value            = diagnostics[["bf_value"]],
      ordinate         = diagnostics[["bf_ordinate"]],
      method           = diagnostics[["estimator"]],
      density_method   = .density_method_normalize(density_method),
      diagnostics      = ordinate_diagnostics,
      evaluation_value = diagnostics[["bf_evaluation_value"]],
      iwmde_provenance = provenance
    ),
    metadata
  )
  out <- do.call(
    BayesTools::posterior_ordinate_attribute,
    args
  )

  if (!allow_rejected && !.iwmde_posterior_ordinate_supports_bf(out)) {
    return(NULL)
  }

  return(out)
}
