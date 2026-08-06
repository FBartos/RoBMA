#' @title Extract qCMDE/IWMDE Density Diagnostics
#'
#' @description
#' Extracts the compact reliability diagnostics attached to a
#' qCMDE/IWMDE hypothesis result. The diagnostics describe the posterior
#' ordinate computation. They do not include uncertainty in the prior ordinate
#' or, for IWMDE, uncertainty from estimating its conditional weight function.
#'
#' @details
#' Point ordinates evaluate deterministic increasing posterior-row budgets,
#' beginning with `initial_samples`, until `target_relative_mcse` and every
#' BF-grade reliability gate pass, all eligible rows are used, or the
#' `max_samples` hard cap is reached. The table records the achieved budget,
#' `hard_cap_reached`, `all_rows_used`, adaptation steps, whether the precision
#' target was met, and whether the final ordinate passed every BF-grade gate.
#' Finite row budgets use a reproducible nested simple random sample without
#' replacement. `relative_mcse` and `ess` describe only the selected continuous
#' posterior-row sequence. They exclude uncertainty in the active/product-space
#' indicator mass. `sampling_relative_mcse` separately estimates
#' finite-population row-sampling uncertainty with its sampling fraction. Both
#' precision checks must pass before adaptive evaluation can stop.
#'
#' The reliability policy warns when relative MCSE is at least 5 percent, ESS is
#' below 100, the largest contribution share is at least 20 percent, fewer than
#' 100 finite contributions remain, or row loss exceeds the estimator-specific
#' warning threshold. It rejects an ordinate when relative MCSE is at least 25
#' percent, ESS is below 20, the largest contribution share is at least 50
#' percent, fewer than 20 finite contributions remain, or row loss exceeds its
#' rejection threshold. The corresponding thresholds are returned as columns
#' rather than being implicit.
#'
#' qCMDE's full method diagnostics distinguish `pilot_normalization_integral`
#' from `final_normalization_integral`; IWMDE uses
#' `support_grid_normalization_integral`. This compact accessor reports their
#' common normalization check as `normalization_relative_error`. The active
#' policy check is identified by `stability_metric`: qCMDE checks posterior
#' ordinate movement, whereas IWMDE checks normalization mass. The table also
#' reports the stability and adaptive-quadrature warning/rejection thresholds.
#'
#' @param object an object returned by [hypothesis()] with
#'   `density_method = "qCMDE"` or `density_method = "IWMDE"`, or an
#'   `RoBMA_density_ordinate_error` caught from a rejected point-ordinate
#'   computation.
#' @param ... unused.
#'
#' @return A data frame of class `RoBMA_density_diagnostics`, with one row per
#' computed point ordinate. Columns identify the estimator, parameter, requested
#' and evaluated values, schema/source provenance, row counts, active mass,
#' relative MCSE, ESS, largest contribution share, finite terms, normalization
#' and quadrature checks, adaptive-budget state, policy thresholds, weight
#' fallbacks, status, and warnings.
#'
#' The exact columns, in order, are `schema_version`, `algorithm_version`,
#' `source_fingerprint`, `estimator`, `density_method`, `parameter`, `level`,
#' `requested_value`, `evaluation_value`, `achieved_row_budget`, `eligible_rows`,
#' `evaluated_rows`, `retained_rows`, `finite_terms`, `row_drop_fraction`,
#' `active_mass`, `relative_mcse`, `sampling_relative_mcse`,
#' `sampling_fraction`, `mcmc_uncertainty_scope`,
#' `sampling_uncertainty_type`, `ess`, `max_weight_share`,
#' `normalization_relative_error`, `stability_metric`,
#' `stability_relative_error`, `ordinate_relative_change`,
#' `quadrature_relative_change`, `target_relative_mcse`,
#' `stability_warning_threshold`, `stability_rejection_threshold`,
#' `quadrature_warning_threshold`, `quadrature_rejection_threshold`,
#' `warning_relative_mcse`, `rejection_relative_mcse`,
#' `warning_min_finite_terms`, `rejection_min_finite_terms`, `warning_min_ess`,
#' `rejection_min_ess`, `warning_max_weight_share`,
#' `rejection_max_weight_share`, `warning_row_drop_fraction`,
#' `rejection_row_drop_fraction`, `hard_cap`, `hard_cap_reached`,
#' `all_rows_used`, `adaptation_steps`, `target_met`, `precision_target_met`,
#' `sampling_target_met`, `bf_grade_met`, `n_weight_fallbacks`,
#' `weight_fallback_reasons`, `status`, and `warnings`.
#'
#' @examples \dontrun{
#' result <- hypothesis(
#'   fit,
#'   hypothesis     = "mu != 0 vs mu = 0",
#'   density_method = "qCMDE"
#' )
#' density_diagnostics(result)
#'
#' failed <- tryCatch(
#'   hypothesis(
#'     fit,
#'     hypothesis     = "mu != 0 vs mu = 0",
#'     density_method = "IWMDE"
#'   ),
#'   RoBMA_density_ordinate_error = identity
#' )
#' density_diagnostics(failed)
#' }
#'
#' @seealso [hypothesis()]
#'
#' @export
density_diagnostics <- function(object, ...) {

  UseMethod("density_diagnostics")
}


#' @export
density_diagnostics.default <- function(object, ...) {

  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "density_diagnostics()"
  )
  stop(
    "No density diagnostics method is available for objects of class ",
    paste0("'", class(object), "'", collapse = ", "), ".",
    call. = FALSE
  )
}


.iwmde_stop_ordinate_unavailable <- function(message, estimate) {

  ordinate    <- estimate[["rejected_posterior_ordinate"]]
  diagnostics <- .iwmde_collect_ordinate_density_diagnostics(ordinate)
  condition <- structure(
    list(
      message             = message,
      call                = NULL,
      density_diagnostics = diagnostics
    ),
    class = c("RoBMA_density_ordinate_error", "error", "condition")
  )

  stop(condition)
}


.iwmde_collect_ordinate_density_diagnostics <- function(ordinate) {

  posterior <- numeric()
  attr(posterior, "posterior_ordinate") <- ordinate

  return(.iwmde_collect_public_density_diagnostics(posterior))
}


#' @export
density_diagnostics.BayesTools_hypothesis_BF <- function(object, ...) {

  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "density_diagnostics()"
  )
  diagnostics <- attr(object, "density_diagnostics", exact = TRUE)

  return(.density_diagnostics_validate(diagnostics))
}


#' @export
density_diagnostics.RoBMA_density_ordinate_error <- function(object, ...) {

  .warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "density_diagnostics()"
  )

  return(.density_diagnostics_validate(object[["density_diagnostics"]]))
}


.iwmde_collect_public_density_diagnostics <- function(posterior) {

  entries <- .iwmde_posterior_ordinate_entries(
    attr(posterior, "posterior_ordinate", exact = TRUE)
  )
  rows <- lapply(entries, .iwmde_public_density_diagnostic_row)

  if (is.list(posterior)) {
    child_rows <- lapply(posterior, .iwmde_collect_public_density_diagnostics)
    child_rows <- child_rows[vapply(child_rows, nrow, integer(1)) > 0L]
    rows       <- c(rows, child_rows)
  }
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (length(rows) == 0L) {
    return(.iwmde_empty_public_density_diagnostics())
  }

  out <- do.call(rbind, rows)
  out <- out[!duplicated(out), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("RoBMA_density_diagnostics", "data.frame")

  return(out)
}


.iwmde_public_density_diagnostic_row <- function(entry) {

  diagnostics <- entry[["diagnostics"]]
  if (!is.list(diagnostics)) {
    return(.iwmde_empty_public_density_diagnostics())
  }
  provenance <- .iwmde_attribute_provenance(entry)
  warnings   <- .iwmde_posterior_ordinate_warnings(entry)
  failure    <- .iwmde_posterior_ordinate_bf_failure_reason(entry)
  policy     <- .iwmde_diagnostic_policy()
  estimator  <- .iwmde_public_character(diagnostics[["estimator"]])
  stability_metric <- if (identical(estimator, "q_grid_cmde")) {
    "ordinate_relative_change"
  } else {
    "normalization_relative_error"
  }
  stability_error <- .iwmde_diagnostics_stability_relative_error(
    diagnostics = diagnostics,
    estimator   = estimator
  )
  stability_warning <- .iwmde_bf_mass_warning_tolerance(estimator)
  stability_rejection <- .iwmde_bf_mass_fail_tolerance(estimator)
  row_loss_warning <- .iwmde_row_loss_warning_tolerance(estimator)
  row_loss_rejection <- .iwmde_row_loss_fail_tolerance(estimator)

  data.frame(
    schema_version = .iwmde_public_character(
      provenance[["schema_version"]]
    ),
    algorithm_version = .iwmde_public_character(
      provenance[["algorithm_version"]]
    ),
    source_fingerprint = .iwmde_public_source_fingerprint(provenance),
    estimator = estimator,
    density_method = .iwmde_public_character(entry[["density_method"]]),
    parameter = .iwmde_public_character(entry[["parameter"]]),
    level = .iwmde_public_character(entry[["level"]]),
    requested_value = .iwmde_public_numeric(entry[["value"]]),
    evaluation_value = .iwmde_public_numeric(
      diagnostics[["evaluation_value"]]
    ),
    achieved_row_budget = .iwmde_public_integer_any(
      diagnostics,
      c("achieved_row_budget", "n_candidate_rows")
    ),
    eligible_rows = .iwmde_public_integer_any(
      diagnostics,
      c("eligible_rows", "n_total")
    ),
    evaluated_rows = .iwmde_public_integer_any(
      diagnostics,
      c("n_evaluated_rows", "n_candidate_rows")
    ),
    retained_rows = .iwmde_public_integer_any(
      diagnostics,
      c("n_normalized_rows", "n_estimator_rows")
    ),
    finite_terms = .iwmde_public_integer_any(
      diagnostics,
      "finite_terms"
    ),
    row_drop_fraction = .iwmde_diagnostics_row_loss_fraction(diagnostics),
    active_mass = .iwmde_public_numeric(diagnostics[["active_mass"]]),
    relative_mcse = .iwmde_public_numeric(diagnostics[["relative_mcse"]]),
    sampling_relative_mcse = .iwmde_public_numeric(
      diagnostics[["sampling_relative_mcse"]]
    ),
    sampling_fraction = .iwmde_public_numeric(
      diagnostics[["sampling_fraction"]]
    ),
    mcmc_uncertainty_scope = .iwmde_public_character(
      diagnostics[["mcmc_uncertainty_scope"]]
    ),
    sampling_uncertainty_type = .iwmde_public_character(
      diagnostics[["sampling_uncertainty_type"]]
    ),
    ess = .iwmde_public_numeric(diagnostics[["ess"]]),
    max_weight_share = .iwmde_public_numeric(
      diagnostics[["max_weight_share"]]
    ),
    normalization_relative_error = .iwmde_public_numeric(
      diagnostics[["normalization_relative_error"]]
    ),
    stability_metric = stability_metric,
    stability_relative_error = stability_error,
    ordinate_relative_change = .iwmde_public_numeric(
      diagnostics[["ordinate_relative_change"]]
    ),
    quadrature_relative_change = .iwmde_public_numeric(
      diagnostics[["max_quadrature_relative_change"]]
    ),
    target_relative_mcse = .iwmde_public_numeric_any(
      diagnostics,
      "target_relative_mcse",
      default = .05
    ),
    stability_warning_threshold = stability_warning,
    stability_rejection_threshold = stability_rejection,
    quadrature_warning_threshold = policy[["quadrature_warn"]],
    quadrature_rejection_threshold = policy[["quadrature_fail"]],
    warning_relative_mcse = policy[["warn_relative_mcse"]],
    rejection_relative_mcse = policy[["max_relative_mcse"]],
    warning_min_finite_terms = policy[["warn_min_finite_terms"]],
    rejection_min_finite_terms = policy[["min_finite_terms"]],
    warning_min_ess = policy[["warn_min_ess"]],
    rejection_min_ess = policy[["min_ess"]],
    warning_max_weight_share = policy[["warn_weight_share"]],
    rejection_max_weight_share = policy[["max_weight_share"]],
    warning_row_drop_fraction = row_loss_warning,
    rejection_row_drop_fraction = row_loss_rejection,
    hard_cap = .iwmde_public_numeric(diagnostics[["hard_cap"]]),
    hard_cap_reached = .iwmde_public_logical(
      diagnostics[["hard_cap_reached"]]
    ),
    all_rows_used = .iwmde_public_logical(diagnostics[["all_rows_used"]]),
    adaptation_steps = .iwmde_public_integer_any(
      diagnostics,
      "n_steps"
    ),
    target_met = .iwmde_public_logical(diagnostics[["target_met"]]),
    precision_target_met = .iwmde_public_logical(
      diagnostics[["precision_target_met"]]
    ),
    sampling_target_met = .iwmde_public_logical(
      diagnostics[["sampling_target_met"]]
    ),
    bf_grade_met = .iwmde_public_logical(diagnostics[["bf_grade_met"]]),
    n_weight_fallbacks = .iwmde_public_integer_any(
      diagnostics,
      "n_weight_fallbacks",
      default = 0L
    ),
    weight_fallback_reasons = .iwmde_public_named_counts(
      diagnostics[["weight_fallback_reasons"]]
    ),
    status = if (is.null(failure)) "ok" else "rejected",
    warnings = paste(warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
}


.iwmde_empty_public_density_diagnostics <- function() {

  out <- data.frame(
    schema_version = character(), algorithm_version = character(),
    source_fingerprint = character(), estimator = character(),
    density_method = character(), parameter = character(), level = character(),
    requested_value = numeric(), evaluation_value = numeric(),
    achieved_row_budget = integer(), eligible_rows = integer(),
    evaluated_rows = integer(), retained_rows = integer(), finite_terms = integer(),
    row_drop_fraction = numeric(), active_mass = numeric(),
    relative_mcse = numeric(), sampling_relative_mcse = numeric(),
    sampling_fraction = numeric(), mcmc_uncertainty_scope = character(),
    sampling_uncertainty_type = character(), ess = numeric(),
    max_weight_share = numeric(),
    normalization_relative_error = numeric(),
    stability_metric = character(), stability_relative_error = numeric(),
    ordinate_relative_change = numeric(), quadrature_relative_change = numeric(),
    target_relative_mcse = numeric(), stability_warning_threshold = numeric(),
    stability_rejection_threshold = numeric(),
    quadrature_warning_threshold = numeric(),
    quadrature_rejection_threshold = numeric(),
    warning_relative_mcse = numeric(),
    rejection_relative_mcse = numeric(), warning_min_finite_terms = numeric(),
    rejection_min_finite_terms = numeric(), warning_min_ess = numeric(),
    rejection_min_ess = numeric(), warning_max_weight_share = numeric(),
    rejection_max_weight_share = numeric(),
    warning_row_drop_fraction = numeric(),
    rejection_row_drop_fraction = numeric(), hard_cap = numeric(),
    hard_cap_reached = logical(), all_rows_used = logical(),
    adaptation_steps = integer(), target_met = logical(),
    precision_target_met = logical(), sampling_target_met = logical(),
    bf_grade_met = logical(),
    n_weight_fallbacks = integer(), weight_fallback_reasons = character(),
    status = character(), warnings = character(), stringsAsFactors = FALSE
  )
  class(out) <- c("RoBMA_density_diagnostics", "data.frame")

  return(out)
}


# Validate the stable public density-diagnostics schema.
.density_diagnostics_validate <- function(diagnostics) {

  template <- .iwmde_empty_public_density_diagnostics()
  if (!inherits(diagnostics, "RoBMA_density_diagnostics") ||
      !is.data.frame(diagnostics) ||
      !identical(names(diagnostics), names(template))) {
    stop(
      "Attached density diagnostics do not satisfy the RoBMA public schema.",
      call. = FALSE
    )
  }

  valid_types <- vapply(seq_along(template), function(i) {
    identical(typeof(diagnostics[[i]]), typeof(template[[i]]))
  }, logical(1))
  if (!all(valid_types)) {
    stop(
      "Attached density diagnostics do not satisfy the RoBMA public schema.",
      call. = FALSE
    )
  }

  return(diagnostics)
}


.iwmde_public_numeric <- function(value) {

  value <- suppressWarnings(as.numeric(value))[1L]
  if (length(value) == 0L) {
    return(NA_real_)
  }

  return(value)
}


.iwmde_public_numeric_any <- function(x, names, default = NA_real_) {

  value <- .iwmde_diagnostic_scalar_any(x, names)
  if (!is.finite(value)) {
    return(default)
  }

  return(value)
}


.iwmde_public_integer_any <- function(x, names, default = NA_integer_) {

  value <- .iwmde_public_numeric_any(x, names, default = default)
  if (!is.finite(value)) {
    return(as.integer(default))
  }

  return(as.integer(value))
}


.iwmde_public_character <- function(value) {

  if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) {
    return(NA_character_)
  }

  return(as.character(value[[1L]]))
}


.iwmde_public_logical <- function(value) {

  if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) {
    return(NA)
  }

  return(isTRUE(value[[1L]]))
}


.iwmde_public_named_counts <- function(x) {

  if (is.null(x) || length(x) == 0L) {
    return("")
  }
  labels <- names(x)
  if (is.null(labels)) {
    labels <- rep("fallback", length(x))
  }

  paste0(labels, "=", as.integer(x), collapse = "; ")
}


.iwmde_public_source_fingerprint <- function(provenance) {

  if (!is.list(provenance) || is.null(provenance[["source_fingerprint"]])) {
    return(NA_character_)
  }

  return(.iwmde_hash("iwmde_public_source", provenance[["source_fingerprint"]]))
}
