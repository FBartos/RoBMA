# ============================================================================ #
# IWMDE Diagnostic Policy and Presentation Adapters
# ============================================================================ #

.iwmde_posterior_ordinate_supports_bf <- function(posterior_ordinate) {

  BayesTools::posterior_ordinate_supports_bf(
    ordinate  = posterior_ordinate,
    validator = .iwmde_posterior_ordinate_bf_validator
  )
}


.iwmde_posterior_ordinate_bf_validator <- function(posterior_ordinate) {

  return(is.null(
    .iwmde_posterior_ordinate_bf_failure_reason(posterior_ordinate)
  ))
}


.iwmde_posterior_ordinate_bf_failure_reason <- function(posterior_ordinate) {

  diagnostics <- posterior_ordinate[["diagnostics"]]
  if (is.null(diagnostics) || !is.list(diagnostics)) {
    return("missing qCMDE/IWMDE BF diagnostics")
  }

  ordinate <- .iwmde_ordinate_scalar(posterior_ordinate, "ordinate")
  if (!is.finite(ordinate)) {
    ordinate <- .iwmde_ordinate_scalar(posterior_ordinate, "y")
  }
  if (!is.finite(ordinate) || ordinate <= 0) {
    return("posterior ordinate is zero or non-finite")
  }

  value <- .iwmde_ordinate_scalar(posterior_ordinate, "evaluation_value")
  if (!is.finite(value)) {
    value <- .iwmde_ordinate_scalar(posterior_ordinate, "value")
  }
  if (is.finite(value)) {
    diagnostics[["evaluation_value"]] <- value
  }

  return(.iwmde_diagnostics_bf_failure_reason(diagnostics))
}


.iwmde_diagnostics_bf_failure_reason <- function(diagnostics) {

  estimator <- diagnostics[["estimator"]]
  if (length(estimator) != 1L) {
    return("missing qCMDE/IWMDE estimator diagnostics")
  }
  estimator <- as.character(estimator)
  if (!estimator %in% c("q_grid_cmde", "iwmde")) {
    return("unknown qCMDE/IWMDE estimator")
  }

  ordinate <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("ordinate", "bf_ordinate")
  )
  if (!is.na(ordinate) && (!is.finite(ordinate) || ordinate <= 0)) {
    return("posterior ordinate is zero or non-finite")
  }

  relative_mcse <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("relative_mcse", "bf_relative_mcse")
  )
  if (!is.finite(relative_mcse) || relative_mcse < 0 ||
      relative_mcse >= .iwmde_bf_max_relative_mcse()) {
    return(paste0(
      "relative MCSE is ",
      .iwmde_percent(relative_mcse),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_bf_max_relative_mcse()),
      ")"
    ))
  }

  finite_terms <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("finite_terms", "bf_finite_terms")
  )
  if (!is.finite(finite_terms) ||
      finite_terms < .iwmde_bf_min_finite_terms()) {
    return(paste0(
      "only ", .iwmde_count(finite_terms),
      " finite importance terms are available (minimum ",
      .iwmde_count(.iwmde_bf_min_finite_terms()), ")"
    ))
  }

  ess <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("ess", "bf_ess")
  )
  if (!is.finite(ess) || ess < .iwmde_bf_min_ess()) {
    return(paste0(
      "effective sample size is ", .iwmde_count(ess),
      " (minimum ", .iwmde_count(.iwmde_bf_min_ess()), ")"
    ))
  }

  max_weight_share <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_weight_share", "bf_max_weight_share")
  )
  if (!is.finite(max_weight_share) ||
      max_weight_share >= .iwmde_bf_max_weight_share()) {
    return(paste0(
      "largest importance weight contributes ",
      .iwmde_percent(max_weight_share),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_bf_max_weight_share()),
      ")"
    ))
  }

  row_loss_reason <- .iwmde_diagnostics_row_loss_failure_reason(
    diagnostics = diagnostics,
    estimator   = estimator
  )
  if (!is.null(row_loss_reason)) {
    return(row_loss_reason)
  }

  if (identical(estimator, "q_grid_cmde")) {
    mass_reason <- .iwmde_diagnostics_mass_failure_reason(
      diagnostics = diagnostics,
      estimator   = estimator
    )
    if (!is.null(mass_reason)) {
      return(mass_reason)
    }
  } else if (identical(estimator, "iwmde")) {
    mass_reason <- .iwmde_diagnostics_mass_failure_reason(
      diagnostics = diagnostics,
      estimator   = estimator
    )
    if (!is.null(mass_reason)) {
      return(mass_reason)
    }
  }

  quadrature_reason <- .iwmde_diagnostics_quadrature_failure_reason(
    diagnostics
  )
  if (!is.null(quadrature_reason)) {
    return(quadrature_reason)
  }

  sampling_relative_mcse <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("sampling_relative_mcse", "bf_sampling_relative_mcse")
  )
  has_sampling_diagnostic <- any(c(
    "sampling_relative_mcse",
    "bf_sampling_relative_mcse"
  ) %in% names(diagnostics))
  if (has_sampling_diagnostic &&
      (!is.finite(sampling_relative_mcse) || sampling_relative_mcse < 0 ||
       sampling_relative_mcse >= .iwmde_bf_max_relative_mcse())) {
    return(paste0(
      "finite-population row-sampling relative error is ",
      .iwmde_percent(sampling_relative_mcse),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_bf_max_relative_mcse()),
      ")"
    ))
  }

  return(NULL)
}


.iwmde_diagnostics_density_failure_reason <- function(diagnostics) {

  if (is.null(diagnostics) || !is.list(diagnostics)) {
    return("missing qCMDE/IWMDE density diagnostics")
  }
  density_diagnostics <- .iwmde_density_diagnostics_as_ordinate(diagnostics)

  estimator <- diagnostics[["estimator"]]
  if (length(estimator) != 1L) {
    return("missing qCMDE/IWMDE estimator diagnostics")
  }
  estimator <- as.character(estimator)
  if (!estimator %in% c("q_grid_cmde", "iwmde")) {
    return("unknown qCMDE/IWMDE estimator")
  }

  reason <- .iwmde_diagnostics_density_sample_failure_reason(density_diagnostics)
  if (!is.null(reason)) {
    return(reason)
  }

  row_loss_reason <- .iwmde_diagnostics_row_loss_failure_reason(
    diagnostics = density_diagnostics,
    estimator   = estimator
  )
  if (!is.null(row_loss_reason)) {
    return(row_loss_reason)
  }

  mass_reason <- .iwmde_diagnostics_mass_failure_reason(
    diagnostics = density_diagnostics,
    estimator   = estimator
  )
  if (!is.null(mass_reason)) {
    return(mass_reason)
  }

  return(.iwmde_diagnostics_quadrature_failure_reason(density_diagnostics))
}


.iwmde_diagnostics_mass_failure_reason <- function(diagnostics, estimator) {

  normalization_error <- .iwmde_diagnostics_stability_relative_error(
    diagnostics = diagnostics,
    estimator   = estimator
  )
  if (!is.finite(normalization_error)) {
    return("normalization diagnostics are unavailable")
  }

  fail_tolerance <- .iwmde_bf_mass_fail_tolerance(estimator)
  if (normalization_error > fail_tolerance) {
    return(paste0(
      .iwmde_estimator_label(estimator),
      .iwmde_diagnostics_normalization_error_phrase(estimator),
      .iwmde_percent(normalization_error),
      " (maximum allowed ",
      .iwmde_percent(fail_tolerance),
      ")"
    ))
  }

  return(NULL)
}


.iwmde_diagnostics_density_warning <- function(diagnostics) {

  if (is.null(diagnostics) || !is.list(diagnostics)) {
    return(character())
  }

  density_diagnostics <- .iwmde_density_diagnostics_as_ordinate(diagnostics)
  estimator <- diagnostics[["estimator"]]
  if (length(estimator) != 1L) {
    return(character())
  }
  estimator <- as.character(estimator)
  if (!estimator %in% c("q_grid_cmde", "iwmde")) {
    return(character())
  }

  warnings <- .iwmde_diagnostics_density_sample_warning(density_diagnostics)
  fallback_warning <- .iwmde_diagnostics_weight_fallback_warning(
    density_diagnostics
  )
  normalization_warning <- .iwmde_diagnostics_density_mass_warning(
    diagnostics = density_diagnostics,
    estimator   = estimator
  )
  quadrature_warning <- .iwmde_diagnostics_quadrature_warning(density_diagnostics)
  warnings <- c(
    warnings,
    fallback_warning,
    normalization_warning,
    quadrature_warning
  )

  return(unique(warnings[nzchar(warnings)]))
}


.iwmde_diagnostics_weight_fallback_warning <- function(diagnostics) {

  fallback_count <- .iwmde_diagnostic_scalar(
    diagnostics,
    "n_weight_fallbacks"
  )
  if (!is.finite(fallback_count) || fallback_count <= 0) {
    return(character())
  }

  fallback_rows <- .iwmde_diagnostic_scalar(
    diagnostics,
    "n_weight_fallback_rows"
  )
  reasons <- diagnostics[["weight_fallback_reasons"]]
  reason_text <- if (is.numeric(reasons) && length(reasons) > 0L &&
                     !is.null(names(reasons))) {
    paste(names(reasons), collapse = "; ")
  } else {
    "conditional weight estimation was unavailable"
  }
  row_text <- if (is.finite(fallback_rows)) {
    paste0(" affecting ", .iwmde_count(fallback_rows), " estimator rows")
  } else {
    ""
  }

  return(paste0(
    "IWMDE used normalized fallback weights in ",
    .iwmde_count(fallback_count), " conditional partitions", row_text,
    ": ", reason_text, "."
  ))
}


.iwmde_diagnostics_density_sample_failure_reason <- function(diagnostics) {

  estimator_rows <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("n_estimator_rows", "n_normalized_rows", "n_evaluated_rows")
  )
  n_active_state_keys <- .iwmde_diagnostic_scalar(
    diagnostics,
    "n_active_state_keys"
  )
  if (is.finite(estimator_rows) &&
      is.finite(n_active_state_keys) &&
      n_active_state_keys > 1 &&
      estimator_rows < .iwmde_density_min_estimator_rows()) {
    return(paste0(
      "model-averaged density uses only ",
      .iwmde_count(estimator_rows),
      " estimator rows across ",
      .iwmde_count(n_active_state_keys),
      " active states (minimum ",
      .iwmde_count(.iwmde_density_min_estimator_rows()),
      ")"
    ))
  }

  metrics <- .iwmde_density_curve_metrics(diagnostics)
  plot_scale_sampling_mcse <-
    metrics[["plot_scale_sampling_relative_mcse"]]
  has_plot_scale_sampling <-
    "plot_scale_sampling_relative_mcse" %in% names(diagnostics)
  if (has_plot_scale_sampling &&
      (!is.finite(plot_scale_sampling_mcse) ||
      plot_scale_sampling_mcse < 0 ||
      plot_scale_sampling_mcse >= .iwmde_density_max_relative_mcse())) {
    return(paste0(
      "density plot-scale finite-population sampling MCSE is ",
      .iwmde_percent(plot_scale_sampling_mcse),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")"
    ))
  }
  bulk_sampling_relative_mcse <-
    metrics[["bulk_sampling_relative_mcse"]]
  has_bulk_sampling <-
    "bulk_max_sampling_relative_mcse" %in% names(diagnostics)
  if (has_bulk_sampling &&
      (!is.finite(bulk_sampling_relative_mcse) ||
      bulk_sampling_relative_mcse < 0 ||
      bulk_sampling_relative_mcse >= .iwmde_density_max_relative_mcse())) {
    return(paste0(
      "density bulk finite-population sampling relative MCSE is ",
      .iwmde_percent(bulk_sampling_relative_mcse),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")"
    ))
  }
  plot_scale_relative_mcse <- metrics[["plot_scale_relative_mcse"]]
  if (!is.finite(plot_scale_relative_mcse) ||
      plot_scale_relative_mcse < 0 ||
      plot_scale_relative_mcse >= .iwmde_density_max_relative_mcse()) {
    return(paste0(
      "density plot-scale MCSE is ",
      .iwmde_percent(plot_scale_relative_mcse),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")"
    ))
  }

  bulk_relative_mcse <- metrics[["bulk_relative_mcse"]]
  if (!is.finite(bulk_relative_mcse) || bulk_relative_mcse < 0 ||
      bulk_relative_mcse >= .iwmde_density_max_relative_mcse()) {
    return(paste0(
      "density bulk relative MCSE is ",
      .iwmde_percent(bulk_relative_mcse),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")"
    ))
  }

  ess <- metrics[["ess"]]
  min_ess <- .iwmde_density_min_ess(estimator_rows)
  if (!is.finite(ess) || ess < min_ess) {
    return(paste0(
      "density bulk effective sample size is ",
      .iwmde_count(ess),
      " (minimum ", .iwmde_count(min_ess), ")"
    ))
  }

  max_weight_share <- metrics[["max_weight_share"]]
  if (!is.finite(max_weight_share) ||
      max_weight_share >= .iwmde_density_max_weight_share()) {
    return(paste0(
      "largest bulk density importance weight contributes ",
      .iwmde_percent(max_weight_share),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_density_max_weight_share()),
      ")"
    ))
  }

  return(NULL)
}


.iwmde_diagnostics_density_sample_warning <- function(diagnostics) {

  warnings <- character()
  estimator_rows <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("n_estimator_rows", "n_normalized_rows", "n_evaluated_rows")
  )
  if (is.finite(estimator_rows) &&
      estimator_rows < .iwmde_density_warning_min_estimator_rows()) {
    warnings <- c(warnings, paste0(
      "Density uses only ",
      .iwmde_count(estimator_rows),
      " estimator rows",
      " (warning threshold ",
      .iwmde_count(.iwmde_density_warning_min_estimator_rows()),
      "; recommended minimum ",
      .iwmde_count(.iwmde_density_min_estimator_rows()),
      ")."
    ))
  }

  metrics <- .iwmde_density_curve_metrics(diagnostics)
  plot_scale_sampling_mcse <-
    metrics[["plot_scale_sampling_relative_mcse"]]
  if (is.finite(plot_scale_sampling_mcse) &&
      plot_scale_sampling_mcse >=
        .iwmde_density_warning_relative_mcse() &&
      plot_scale_sampling_mcse < .iwmde_density_max_relative_mcse()) {
    warnings <- c(warnings, paste0(
      "Density plot-scale finite-population sampling MCSE is ",
      .iwmde_percent(plot_scale_sampling_mcse),
      " (warning threshold ",
      .iwmde_percent(.iwmde_density_warning_relative_mcse()),
      "; rejection threshold ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")."
    ))
  }
  bulk_sampling_relative_mcse <-
    metrics[["bulk_sampling_relative_mcse"]]
  if (is.finite(bulk_sampling_relative_mcse) &&
      bulk_sampling_relative_mcse >=
        .iwmde_density_warning_relative_mcse() &&
      bulk_sampling_relative_mcse < .iwmde_density_max_relative_mcse()) {
    warnings <- c(warnings, paste0(
      "Density bulk finite-population sampling relative MCSE is ",
      .iwmde_percent(bulk_sampling_relative_mcse),
      " (warning threshold ",
      .iwmde_percent(.iwmde_density_warning_relative_mcse()),
      "; rejection threshold ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")."
    ))
  }
  plot_scale_relative_mcse <- metrics[["plot_scale_relative_mcse"]]
  if (is.finite(plot_scale_relative_mcse) &&
      plot_scale_relative_mcse >=
        .iwmde_density_warning_relative_mcse() &&
      plot_scale_relative_mcse < .iwmde_density_max_relative_mcse()) {
    warnings <- c(warnings, paste0(
      "Density plot-scale MCSE is ",
      .iwmde_percent(plot_scale_relative_mcse),
      " (warning threshold ",
      .iwmde_percent(.iwmde_density_warning_relative_mcse()),
      "; rejection threshold ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")."
    ))
  }

  bulk_relative_mcse <- metrics[["bulk_relative_mcse"]]
  if (is.finite(bulk_relative_mcse) &&
      bulk_relative_mcse >= .iwmde_density_warning_relative_mcse() &&
      bulk_relative_mcse < .iwmde_density_max_relative_mcse()) {
    warnings <- c(warnings, paste0(
      "Density bulk relative MCSE is ",
      .iwmde_percent(bulk_relative_mcse),
      " (warning threshold ",
      .iwmde_percent(.iwmde_density_warning_relative_mcse()),
      "; rejection threshold ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")."
    ))
  }

  ess <- metrics[["ess"]]
  min_ess <- .iwmde_density_min_ess(estimator_rows)
  warning_min_ess <- .iwmde_density_warning_min_ess(estimator_rows)
  if (is.finite(ess) &&
      ess >= min_ess &&
      ess < warning_min_ess) {
    warnings <- c(warnings, paste0(
      "Density bulk effective sample size is ",
      .iwmde_count(ess),
      " (warning threshold ",
      .iwmde_count(warning_min_ess),
      "; rejection threshold ",
      .iwmde_count(min_ess),
      ")."
    ))
  }

  max_weight_share <- metrics[["max_weight_share"]]
  if (is.finite(max_weight_share) &&
      max_weight_share >= .iwmde_density_warning_weight_share() &&
      max_weight_share < .iwmde_density_max_weight_share()) {
    warnings <- c(warnings, paste0(
      "Largest bulk density importance weight contributes ",
      .iwmde_percent(max_weight_share),
      " (warning threshold ",
      .iwmde_percent(.iwmde_density_warning_weight_share()),
      "; rejection threshold ",
      .iwmde_percent(.iwmde_density_max_weight_share()),
      ")."
    ))
  }

  return(warnings)
}


.iwmde_density_curve_metrics <- function(diagnostics) {

  return(list(
    plot_scale_relative_mcse = .iwmde_diagnostic_scalar(
      diagnostics,
      "plot_scale_relative_mcse"
    ),
    plot_scale_sampling_relative_mcse = .iwmde_diagnostic_scalar(
      diagnostics,
      "plot_scale_sampling_relative_mcse"
    ),
    bulk_relative_mcse = .iwmde_diagnostic_scalar(
      diagnostics,
      "bulk_max_relative_mcse"
    ),
    bulk_sampling_relative_mcse = .iwmde_diagnostic_scalar(
      diagnostics,
      "bulk_max_sampling_relative_mcse"
    ),
    ess = .iwmde_diagnostic_scalar(
      diagnostics,
      "bulk_min_ess"
    ),
    max_weight_share = .iwmde_diagnostic_scalar(
      diagnostics,
      "bulk_max_weight_share"
    )
  ))
}


.iwmde_diagnostics_density_mass_warning <- function(diagnostics, estimator) {

  warnings <- character()
  row_loss <- .iwmde_diagnostics_row_loss_fraction(diagnostics)
  if (is.finite(row_loss)) {
    warning_tolerance <- .iwmde_row_loss_warning_tolerance(estimator)
    fail_tolerance    <- .iwmde_row_loss_fail_tolerance(estimator)
    if (row_loss > warning_tolerance &&
        row_loss <= fail_tolerance) {
      warnings <- c(warnings, paste0(
        .iwmde_estimator_label(estimator),
        " dropped ",
        .iwmde_percent(row_loss),
        " of target rows during density estimation",
        " (warning threshold ",
        .iwmde_percent(warning_tolerance),
        "; rejection threshold ",
        .iwmde_percent(fail_tolerance),
        ")."
      ))
    }
  }

  normalization_error <- .iwmde_diagnostics_stability_relative_error(
    diagnostics = diagnostics,
    estimator   = estimator
  )
  warning_tolerance <- .iwmde_bf_mass_warning_tolerance(estimator)
  fail_tolerance    <- .iwmde_bf_mass_fail_tolerance(estimator)
  if (is.finite(normalization_error) &&
      normalization_error > warning_tolerance &&
      normalization_error <= fail_tolerance) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      .iwmde_diagnostics_normalization_error_phrase(estimator),
      .iwmde_percent(normalization_error),
      " (warning threshold ",
      .iwmde_percent(warning_tolerance),
      "; rejection threshold ",
      .iwmde_percent(fail_tolerance),
      ")."
    ))
  }

  return(warnings)
}


.iwmde_diagnostics_bf_warning <- function(diagnostics) {

  estimator <- diagnostics[["estimator"]]
  if (length(estimator) != 1L) {
    return(character())
  }
  estimator <- as.character(estimator)
  if (!estimator %in% c("q_grid_cmde", "iwmde")) {
    return(character())
  }

  warnings <- character()

  ordinate_warnings <- diagnostics[["ordinate_warnings"]]
  if (!is.null(ordinate_warnings)) {
    warnings <- c(warnings, as.character(ordinate_warnings))
  }
  warnings <- c(
    warnings,
    .iwmde_diagnostics_weight_fallback_warning(diagnostics)
  )

  relative_mcse <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("relative_mcse", "bf_relative_mcse")
  )
  if (is.finite(relative_mcse) &&
      relative_mcse >= .iwmde_bf_warning_relative_mcse() &&
      relative_mcse < .iwmde_bf_max_relative_mcse()) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      " relative MCSE is ",
      .iwmde_percent(relative_mcse)
    ))
  }

  sampling_relative_mcse <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("sampling_relative_mcse", "bf_sampling_relative_mcse")
  )
  if (is.finite(sampling_relative_mcse) &&
      sampling_relative_mcse >= .iwmde_bf_warning_relative_mcse() &&
      sampling_relative_mcse < .iwmde_bf_max_relative_mcse()) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      " finite-population row-sampling relative error is ",
      .iwmde_percent(sampling_relative_mcse)
    ))
  }

  finite_terms <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("finite_terms", "bf_finite_terms")
  )
  if (is.finite(finite_terms) &&
      finite_terms >= .iwmde_bf_min_finite_terms() &&
      finite_terms < .iwmde_bf_warning_min_finite_terms()) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      " uses only ",
      .iwmde_count(finite_terms),
      " finite importance terms",
      " (warning threshold ",
      .iwmde_count(.iwmde_bf_warning_min_finite_terms()),
      "; BF rejection threshold ",
      .iwmde_count(.iwmde_bf_min_finite_terms()), ")."
    ))
  }

  ess <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("ess", "bf_ess")
  )
  if (is.finite(ess) &&
      ess >= .iwmde_bf_min_ess() &&
      ess < .iwmde_bf_warning_min_ess()) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      " effective sample size is ",
      .iwmde_count(ess)
    ))
  }

  max_weight_share <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_weight_share", "bf_max_weight_share")
  )
  if (is.finite(max_weight_share) &&
      max_weight_share >= .iwmde_bf_warning_weight_share() &&
      max_weight_share < .iwmde_bf_max_weight_share()) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      " largest importance weight contributes ",
      .iwmde_percent(max_weight_share),
      " (warning threshold ",
      .iwmde_percent(.iwmde_bf_warning_weight_share()),
      "; BF rejection threshold ",
      .iwmde_percent(.iwmde_bf_max_weight_share()),
      ")."
    ))
  }

  row_loss <- .iwmde_diagnostics_row_loss_fraction(diagnostics)
  if (is.finite(row_loss)) {
    warning_tolerance <- .iwmde_row_loss_warning_tolerance(estimator)
    fail_tolerance    <- .iwmde_row_loss_fail_tolerance(estimator)
    if (row_loss > warning_tolerance &&
        row_loss <= fail_tolerance) {
      warnings <- c(warnings, paste0(
        .iwmde_estimator_label(estimator),
        " dropped ",
        .iwmde_percent(row_loss),
        " of target rows during density estimation",
        " (warning threshold ",
        .iwmde_percent(warning_tolerance),
        "; BF rejection threshold ",
        .iwmde_percent(fail_tolerance),
        ")."
      ))
    }
  }

  normalization_error <- .iwmde_diagnostics_stability_relative_error(
    diagnostics = diagnostics,
    estimator   = estimator
  )
  warning_tolerance <- .iwmde_bf_mass_warning_tolerance(estimator)
  fail_tolerance    <- .iwmde_bf_mass_fail_tolerance(estimator)
  if (is.finite(normalization_error) &&
      normalization_error > warning_tolerance &&
      normalization_error <= fail_tolerance) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      .iwmde_diagnostics_normalization_error_phrase(estimator),
      .iwmde_percent(normalization_error),
      " (warning threshold ",
      .iwmde_percent(warning_tolerance),
      "; BF rejection threshold ",
      .iwmde_percent(fail_tolerance),
      ")."
    ))
  }

  quadrature_warning <- .iwmde_diagnostics_quadrature_warning(diagnostics)
  if (length(quadrature_warning) > 0L) {
    warnings <- c(warnings, quadrature_warning)
  }

  return(unique(warnings[nzchar(warnings)]))
}


.iwmde_density_diagnostics_as_ordinate <- function(diagnostics) {

  out <- diagnostics
  out[["relative_mcse"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_relative_mcse", "relative_mcse")
  )
  out[["finite_terms"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("min_finite_terms", "finite_terms")
  )
  out[["ess"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("min_ess", "ess")
  )
  out[["max_weight_share"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_weight_share", "bf_max_weight_share")
  )
  out[["ordinate_relative_change"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_ordinate_relative_change", "ordinate_relative_change")
  )

  return(out)
}


.iwmde_diagnostics_quadrature_failure_reason <- function(diagnostics) {

  quadrature_change <- .iwmde_diagnostics_quadrature_relative_change(diagnostics)
  if (!is.finite(quadrature_change)) {
    return(NULL)
  }
  if (quadrature_change > .iwmde_quadrature_fail_tolerance()) {
    return(paste0(
      "adaptive quadrature sensitivity is ",
      .iwmde_percent(quadrature_change),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_quadrature_fail_tolerance()),
      ")"
    ))
  }

  return(NULL)
}


.iwmde_diagnostics_quadrature_warning <- function(diagnostics) {

  quadrature_change <- .iwmde_diagnostics_quadrature_relative_change(diagnostics)
  if (!is.finite(quadrature_change) ||
      quadrature_change <= .iwmde_quadrature_warning_tolerance() ||
      quadrature_change > .iwmde_quadrature_fail_tolerance()) {
    return(character())
  }

  return(paste0(
    "Adaptive quadrature sensitivity is ",
    .iwmde_percent(quadrature_change),
    " (warning threshold ",
    .iwmde_percent(.iwmde_quadrature_warning_tolerance()),
    "; rejection threshold ",
    .iwmde_percent(.iwmde_quadrature_fail_tolerance()),
    ")."
  ))
}


.iwmde_diagnostics_quadrature_relative_change <- function(diagnostics) {

  quadrature_change <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_quadrature_relative_change", "quadrature_relative_change")
  )
  if (!is.finite(quadrature_change) || quadrature_change < 0) {
    return(NA_real_)
  }

  return(quadrature_change)
}


.iwmde_diagnostics_row_loss_failure_reason <- function(diagnostics, estimator) {

  row_loss <- .iwmde_diagnostics_row_loss_fraction(diagnostics)
  if (!is.finite(row_loss)) {
    return(NULL)
  }

  fail_tolerance <- .iwmde_row_loss_fail_tolerance(estimator)
  if (row_loss > fail_tolerance) {
    return(paste0(
      .iwmde_estimator_label(estimator),
      " dropped ",
      .iwmde_percent(row_loss),
      " of target rows during density estimation",
      " (maximum allowed ",
      .iwmde_percent(fail_tolerance),
      ")"
    ))
  }

  return(NULL)
}


.iwmde_diagnostics_row_loss_fraction <- function(diagnostics) {

  row_loss <- .iwmde_diagnostic_scalar(diagnostics, "row_drop_fraction")
  if (is.finite(row_loss) && row_loss >= 0) {
    return(min(1, row_loss))
  }

  n_candidate_rows <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("n_candidate_rows", "n_active")
  )
  n_normalized_rows <- .iwmde_diagnostic_scalar(
    diagnostics,
    "n_normalized_rows"
  )
  if (!is.finite(n_candidate_rows) || n_candidate_rows <= 0 ||
      !is.finite(n_normalized_rows) || n_normalized_rows < 0) {
    return(NA_real_)
  }

  return(.iwmde_row_drop_fraction(
    n_candidate_rows  = n_candidate_rows,
    n_normalized_rows = n_normalized_rows
  ))
}


.iwmde_diagnostics_stability_relative_error <- function(diagnostics,
                                                        estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return(.iwmde_diagnostics_qcmde_ordinate_relative_change(diagnostics))
  }

  return(.iwmde_diagnostics_normalization_relative_error(diagnostics))
}


.iwmde_diagnostics_qcmde_ordinate_relative_change <- function(diagnostics) {

  ordinate_change <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("ordinate_relative_change", "bf_ordinate_relative_change")
  )
  if (!is.finite(ordinate_change) || ordinate_change < 0) {
    return(NA_real_)
  }

  return(ordinate_change)
}


.iwmde_diagnostics_normalization_error_phrase <- function(estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return(" posterior ordinate changes by ")
  }

  return(" normalization mass differs from active posterior mass by ")
}


.iwmde_diagnostics_normalization_relative_error <- function(diagnostics) {

  relative_error <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    "normalization_relative_error"
  )
  if (is.finite(relative_error) && relative_error >= 0) {
    return(relative_error)
  }

  normalization_integral <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c(
      "support_grid_normalization_integral",
      "final_normalization_integral"
    )
  )
  active_mass <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    "active_mass"
  )
  if (!is.finite(normalization_integral) || !is.finite(active_mass) ||
      active_mass <= 0) {
    return(NA_real_)
  }

  return(abs(normalization_integral / active_mass - 1))
}


.iwmde_ordinate_scalar <- function(posterior_ordinate, name) {

  value <- posterior_ordinate[[name]]
  if (is.null(value)) {
    return(NA_real_)
  }

  value <- suppressWarnings(as.numeric(value))[1]
  if (length(value) == 0L) {
    return(NA_real_)
  }

  return(value)
}

.iwmde_diagnostic_scalar <- function(diagnostics, name) {

  value <- diagnostics[[name]]
  if (is.null(value)) {
    return(NA_real_)
  }
  if (is.list(value)) {
    value <- unlist(value, use.names = FALSE)
  }

  value <- suppressWarnings(as.numeric(value))[1]
  if (length(value) == 0L) {
    return(NA_real_)
  }

  return(value)
}


.iwmde_diagnostic_scalar_any <- function(diagnostics, names) {

  for (name in names) {
    value <- .iwmde_diagnostic_scalar(diagnostics, name)
    if (!is.na(value)) {
      return(value)
    }
  }

  return(NA_real_)
}


.iwmde_estimator_label <- function(estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return("qCMDE")
  }
  if (identical(estimator, "iwmde")) {
    return("IWMDE")
  }

  return(as.character(estimator)[1])
}


.iwmde_percent <- function(value) {

  if (!is.finite(value)) {
    return("NA")
  }

  paste0(formatC(100 * value, digits = 3, format = "fg"), "%")
}


.iwmde_count <- function(value) {

  if (!is.finite(value)) {
    return("NA")
  }

  formatC(value, digits = 3, format = "fg")
}


.iwmde_posterior_ordinate_warnings <- function(posterior_ordinate) {

  entries <- .iwmde_posterior_ordinate_entries(posterior_ordinate)
  warnings <- unlist(lapply(entries, function(entry) {
    diagnostics <- entry[["diagnostics"]]
    if (is.null(diagnostics) || !is.list(diagnostics)) {
      return(character())
    }
    explicit <- diagnostics[["warning"]]
    if (!is.null(explicit)) {
      explicit <- as.character(explicit)
      explicit <- explicit[nzchar(explicit)]
      if (length(explicit) > 0L) {
        return(explicit)
      }
    }
    .iwmde_diagnostics_bf_warning(diagnostics)
  }), use.names = FALSE)

  unique(warnings[nzchar(warnings)])
}


.iwmde_posterior_ordinate_failure_reasons <- function(posterior_ordinate) {

  entries <- .iwmde_posterior_ordinate_entries(posterior_ordinate)
  reasons <- unlist(lapply(entries, function(entry) {
    reason <- .iwmde_posterior_ordinate_bf_failure_reason(entry)
    if (is.null(reason)) {
      return(character())
    }
    reason
  }), use.names = FALSE)

  unique(reasons[nzchar(reasons)])
}


.iwmde_posterior_ordinate_entries <- function(posterior_ordinate) {

  if (is.null(posterior_ordinate)) {
    return(list())
  }
  if (is.list(posterior_ordinate) &&
      !is.data.frame(posterior_ordinate[["ordinates"]]) &&
      is.list(posterior_ordinate[["ordinates"]]) &&
      is.null(posterior_ordinate[["ordinates"]][["x"]]) &&
      is.null(posterior_ordinate[["ordinates"]][["value"]]) &&
      is.null(posterior_ordinate[["ordinates"]][["null_hypothesis"]])) {
    return(posterior_ordinate[["ordinates"]])
  }

  return(list(posterior_ordinate))
}


.iwmde_bf_append_warning <- function(bf, posterior_ordinate) {

  warnings <- .iwmde_posterior_ordinate_warnings(posterior_ordinate)
  if (length(warnings) == 0L) {
    return(bf)
  }

  attr(bf, "warnings") <- unique(c(
    attr(bf, "warnings", exact = TRUE),
    warnings
  ))

  return(bf)
}


.iwmde_density_bf_diagnostics <- function(density, values) {

  out <- list(
    bf_value            = NA_real_,
    bf_evaluation_value = NA_real_,
    bf_included         = FALSE,
    bf_grid_index       = NA_integer_,
    bf_ordinate         = NA_real_,
    bf_pilot_ordinate   = NA_real_,
    bf_validation_ordinate = NA_real_,
    bf_ordinate_relative_change = NA_real_,
    bf_ordinate_log_change      = NA_real_,
    bf_pilot_ordinate_relative_change = NA_real_,
    bf_pilot_ordinate_log_change      = NA_real_,
    bf_mcse             = NA_real_,
    bf_relative_mcse    = NA_real_,
    bf_error_percent    = NA_real_,
    bf_finite_terms     = NA_integer_,
    bf_ess              = NA_real_,
    bf_max_weight_share = NA_real_,
    bf_max_log_ratio    = NA_real_
  )

  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(out)
  }
  out[["bf_value"]] <- values[1]

  x <- as.numeric(density[["x"]])
  if (length(x) == 0L) {
    return(out)
  }
  index <- which(x == out[["bf_value"]])
  if (length(index) == 0L) {
    return(out)
  }
  index <- index[1]

  out[["bf_included"]]         <- TRUE
  out[["bf_grid_index"]]       <- index
  out[["bf_evaluation_value"]] <- .iwmde_density_evaluation_value(density, index)
  out[["bf_ordinate"]]         <- .iwmde_density_index_value(density, "y", index)
  out[["bf_pilot_ordinate"]]   <- .iwmde_density_index_value(density, "pilot_y", index)
  out[["bf_validation_ordinate"]] <-
    .iwmde_density_index_value(density, "validation_y", index)
  out[["bf_ordinate_relative_change"]] <-
    .iwmde_density_index_value(density, "ordinate_relative_change", index)
  out[["bf_ordinate_log_change"]] <-
    .iwmde_density_index_value(density, "ordinate_log_change", index)
  out[["bf_pilot_ordinate_relative_change"]] <-
    .iwmde_density_index_value(density, "pilot_ordinate_relative_change", index)
  out[["bf_pilot_ordinate_log_change"]] <-
    .iwmde_density_index_value(density, "pilot_ordinate_log_change", index)
  out[["bf_mcse"]]             <- .iwmde_density_index_value(density, "mcse", index)
  out[["bf_relative_mcse"]]    <- .iwmde_density_index_value(density, "relative_mcse", index)
  if (!is.null(density[["sampling_mcse"]])) {
    out[["bf_sampling_mcse"]] <-
      .iwmde_density_index_value(density, "sampling_mcse", index)
  }
  if (!is.null(density[["sampling_relative_mcse"]])) {
    out[["bf_sampling_relative_mcse"]] <-
      .iwmde_density_index_value(density, "sampling_relative_mcse", index)
  }
  out[["bf_finite_terms"]]     <- as.integer(.iwmde_density_index_value(density, "finite_terms", index))
  out[["bf_ess"]]              <- .iwmde_density_index_value(density, "ess", index)
  out[["bf_max_weight_share"]] <- .iwmde_density_index_value(density, "max_weight_share", index)
  out[["bf_max_log_ratio"]]    <- .iwmde_density_index_value(density, "max_log_ratio", index)

  if (is.finite(out[["bf_relative_mcse"]])) {
    out[["bf_error_percent"]] <- 100 * out[["bf_relative_mcse"]]
  }

  return(out)
}


.iwmde_density_evaluation_value <- function(density, index) {

  value <- .iwmde_density_index_value(density, "evaluation_x", index)
  if (is.finite(value)) {
    return(value)
  }

  return(.iwmde_density_index_value(density, "x", index))
}


.iwmde_density_index_value <- function(density, name, index) {

  value <- density[[name]]
  if (is.null(value) || length(value) < index) {
    return(NA_real_)
  }

  value <- suppressWarnings(as.numeric(value[index]))[1]
  if (length(value) == 0L) {
    return(NA_real_)
  }

  return(value)
}
