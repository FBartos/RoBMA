# ============================================================================ #
# IWMDE Internal Schemas
# ============================================================================ #

.iwmde_validate_required_fields <- function(x, required, schema) {

  if (!is.list(x)) {
    stop("Internal ", schema, " must be a list.", call. = FALSE)
  }
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(
      "Internal ", schema, " is missing required field(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }

  invisible(x)
}


.iwmde_validate_plan <- function(plan) {

  .iwmde_validate_required_fields(
    plan,
    c(
      "schema_version", "target", "outputs", "control", "row_budget", "method",
      "density_method", "source_fingerprint", "status", "rows"
    ),
    "iwmde_plan"
  )
  if (!identical(plan[["schema_version"]], .iwmde_schema_version())) {
    stop("Internal iwmde_plan has an incompatible schema version.",
         call. = FALSE)
  }
  scalar_character <- function(value) {

    is.character(value) && length(value) == 1L && !is.na(value)
  }
  if (!scalar_character(plan[["status"]]) ||
      !plan[["status"]] %in% c("ok", "point_only", "unsupported")) {
    stop("Internal iwmde_plan has an unknown status.", call. = FALSE)
  }
  if (!scalar_character(plan[["method"]]) ||
      !plan[["method"]] %in% c("q_grid_cmde", "iwmde")) {
    stop("Internal iwmde_plan has an unknown estimator.", call. = FALSE)
  }
  if (!scalar_character(plan[["density_method"]]) ||
      !plan[["density_method"]] %in% c("qCMDE", "IWMDE")) {
    stop("Internal iwmde_plan has an unknown density method.", call. = FALSE)
  }
  if (!is.numeric(plan[["row_budget"]]) ||
      length(plan[["row_budget"]]) != 1L || is.na(plan[["row_budget"]]) ||
      plan[["row_budget"]] < 20) {
    stop("Internal iwmde_plan has an invalid row budget.", call. = FALSE)
  }
  if (!is.list(plan[["target"]]) || !is.list(plan[["outputs"]]) ||
      !is.list(plan[["control"]]) || !is.list(plan[["rows"]]) ||
      !is.list(plan[["source_fingerprint"]])) {
    stop("Internal iwmde_plan has an invalid structural field.", call. = FALSE)
  }
  if (identical(plan[["status"]], "ok")) {
    .iwmde_validate_required_fields(
      plan,
      c("support", "replacement", "grids"),
      "executable iwmde_plan"
    )
  }

  invisible(plan)
}


.iwmde_new_plan <- function(fields) {

  fields[["schema_version"]] <- .iwmde_schema_version()
  class(fields) <- c("iwmde_plan", "list")
  .iwmde_validate_plan(fields)

  return(fields)
}


.iwmde_new_row_state <- function(fields) {

  if (!is.list(fields)) {
    stop("Internal IWMDE row state must be a list.", call. = FALSE)
  }
  fields[["schema_version"]] <- .iwmde_schema_version()
  class(fields) <- c("iwmde_row_state", "list")
  .iwmde_validate_row_state(fields)

  return(fields)
}


.iwmde_validate_row_state <- function(state) {

  .iwmde_validate_required_fields(
    state,
    c("schema_version", "baseline_log_q"),
    "IWMDE row state"
  )
  if (!identical(state[["schema_version"]], .iwmde_schema_version())) {
    stop("Internal IWMDE row state has an incompatible schema version.",
         call. = FALSE)
  }
  baseline_log_q <- state[["baseline_log_q"]]
  if (!is.numeric(baseline_log_q) || length(baseline_log_q) != 1L ||
      is.na(baseline_log_q)) {
    stop("Internal IWMDE row state has an invalid baseline log density.",
         call. = FALSE)
  }

  invisible(state)
}


.iwmde_validate_row_states <- function(row_states) {

  if (!is.list(row_states)) {
    stop("Internal IWMDE row states must be a list.", call. = FALSE)
  }
  if (length(row_states) == 0L) {
    return(invisible(row_states))
  }
  valid <- vapply(row_states, function(state) {
    tryCatch({
      .iwmde_validate_row_state(state)
      is.finite(state[["baseline_log_q"]])
    }, error = function(e) FALSE)
  }, logical(1))
  if (!all(valid)) {
    stop(
      "Internal IWMDE row states contain an invalid baseline log density.",
      call. = FALSE
    )
  }

  invisible(row_states)
}


.iwmde_validate_density_result <- function(density, estimator = NULL) {

  common_fields <- c(
    "schema_version", "x", "y", "finite_terms", "max_log_ratio", "ess",
    "max_weight_share", "mcse", "relative_mcse", "n_candidate_rows",
    "n_evaluated_rows", "n_normalized_rows", "n_dropped_rows",
    "row_drop_fraction", "normalization_points", "normalization_range",
    "normalization_relative_error", "normalization_scale",
    "normalization_mass_ratio", "max_normalizer_relative_change",
    "max_quadrature_relative_change", "median_normalizer_relative_change",
    "normalization_refined_points", "normalization_refined_range",
    "integral_mcse", "integral_relative_mcse", "batch_size", "n_batches",
    "estimator", "weight_method"
  )
  .iwmde_validate_required_fields(
    density,
    common_fields,
    "IWMDE density result"
  )
  if (!identical(density[["schema_version"]], .iwmde_schema_version())) {
    stop("Internal IWMDE density result has an incompatible schema version.",
         call. = FALSE)
  }
  if (!is.numeric(density[["x"]]) || !is.numeric(density[["y"]]) ||
      length(density[["x"]]) < 1L ||
      any(!is.finite(density[["x"]])) ||
      is.unsorted(density[["x"]], strictly = TRUE) ||
      any(!is.finite(density[["y"]])) || any(density[["y"]] < 0)) {
    stop("Internal IWMDE density result has a non-numeric grid.",
         call. = FALSE)
  }
  n_values <- length(density[["x"]])
  vector_fields <- c(
    "y", "finite_terms", "max_log_ratio", "ess", "max_weight_share",
    "mcse", "relative_mcse"
  )
  valid_lengths <- vapply(
    vector_fields,
    function(name) length(density[[name]]) == n_values,
    logical(1)
  )
  if (!all(valid_lengths)) {
    stop(
      "Internal IWMDE density result has inconsistent grid-field lengths.",
      call. = FALSE
    )
  }
  if (!all(vapply(
    density[vector_fields],
    is.numeric,
    logical(1)
  ))) {
    stop("Internal IWMDE density result has a non-numeric diagnostic field.",
         call. = FALSE)
  }
  finite_terms <- density[["finite_terms"]]
  if (any(!is.finite(finite_terms)) || any(finite_terms < 0) ||
      any(finite_terms != as.integer(finite_terms))) {
    stop("Internal IWMDE density result has invalid finite-term counts.",
         call. = FALSE)
  }
  if (any(!is.finite(density[["ess"]])) || any(density[["ess"]] < 0) ||
      any(density[["ess"]] > finite_terms + sqrt(.Machine$double.eps))) {
    stop("Internal IWMDE density result has invalid effective sample sizes.",
         call. = FALSE)
  }
  if (any(!is.finite(density[["max_weight_share"]])) ||
      any(density[["max_weight_share"]] < 0) ||
      any(density[["max_weight_share"]] > 1)) {
    stop("Internal IWMDE density result has invalid maximum weight shares.",
         call. = FALSE)
  }
  nonnegative_or_na <- function(value) {

    all(is.na(value) | (is.finite(value) & value >= 0))
  }
  if (!nonnegative_or_na(density[["mcse"]]) ||
      !nonnegative_or_na(density[["relative_mcse"]])) {
    stop("Internal IWMDE density result has invalid Monte Carlo errors.",
         call. = FALSE)
  }
  scalar_counts <- c(
    "n_candidate_rows", "n_evaluated_rows", "n_normalized_rows",
    "n_dropped_rows", "normalization_points", "normalization_refined_points",
    "n_batches"
  )
  if (!all(vapply(density[scalar_counts], function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value >= 0 && value == as.integer(value)
  }, logical(1)))) {
    stop("Internal IWMDE density result has an invalid row count.",
         call. = FALSE)
  }
  if (density[["n_normalized_rows"]] > density[["n_evaluated_rows"]] ||
      density[["n_evaluated_rows"]] > density[["n_candidate_rows"]] ||
      density[["n_dropped_rows"]] !=
        density[["n_candidate_rows"]] - density[["n_normalized_rows"]]) {
    stop("Internal IWMDE density result has inconsistent row counts.",
         call. = FALSE)
  }
  if (!is.numeric(density[["row_drop_fraction"]]) ||
      length(density[["row_drop_fraction"]]) != 1L ||
      !is.finite(density[["row_drop_fraction"]]) ||
      density[["row_drop_fraction"]] < 0 ||
      density[["row_drop_fraction"]] > 1) {
    stop("Internal IWMDE density result has an invalid row-drop fraction.",
         call. = FALSE)
  }
  if (!is.character(density[["estimator"]]) ||
      length(density[["estimator"]]) != 1L ||
      !density[["estimator"]] %in% c("q_grid_cmde", "iwmde")) {
    stop("Internal IWMDE density result has an unknown estimator.",
         call. = FALSE)
  }
  if (!is.character(density[["weight_method"]]) ||
      length(density[["weight_method"]]) != 1L ||
      is.na(density[["weight_method"]]) ||
      !nzchar(density[["weight_method"]])) {
    stop("Internal IWMDE density result has an invalid weight method.",
         call. = FALSE)
  }
  range_fields <- c("normalization_range", "normalization_refined_range")
  if (!all(vapply(density[range_fields], .iwmde_schema_numeric_range,
                  logical(1)))) {
    stop("Internal IWMDE density result has an invalid normalization range.",
         call. = FALSE)
  }
  common_nonnegative_fields <- c(
    "normalization_relative_error", "normalization_mass_ratio",
    "max_normalizer_relative_change", "max_quadrature_relative_change",
    "median_normalizer_relative_change", "integral_mcse",
    "integral_relative_mcse"
  )
  if (!all(vapply(
    density[common_nonnegative_fields],
    .iwmde_schema_nonnegative_scalar,
    logical(1)
  ))) {
    stop("Internal IWMDE density result has an invalid numerical diagnostic.",
         call. = FALSE)
  }
  normalization_scale <- density[["normalization_scale"]]
  if (!is.character(normalization_scale) ||
      length(normalization_scale) != 1L ||
      (!is.na(normalization_scale) && !nzchar(normalization_scale)) ||
      (density[["normalization_points"]] > 0L && is.na(normalization_scale))) {
    stop("Internal IWMDE density result has an invalid normalization scale.",
         call. = FALSE)
  }
  batch_size <- density[["batch_size"]]
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      (!is.na(batch_size) &&
        (!is.finite(batch_size) || batch_size < 1 ||
          batch_size != as.integer(batch_size)))) {
    stop("Internal IWMDE density result has an invalid batch size.",
         call. = FALSE)
  }
  if (!is.null(estimator) &&
      !identical(as.character(density[["estimator"]]), estimator)) {
    stop("Internal IWMDE density result changed estimator.", call. = FALSE)
  }
  estimator_fields <- if (identical(density[["estimator"]], "q_grid_cmde")) {
    c(
      "log_normalizer", "pilot_log_normalizer",
      "normalization_initial_points", "normalization_initial_range",
      "pilot_normalization_integral", "final_normalization_integral",
      "pilot_y", "validation_y", "ordinate_relative_change",
      "ordinate_log_change", "pilot_ordinate_relative_change",
      "pilot_ordinate_log_change", "p95_normalizer_relative_change",
      "n_rescued_normalizer", "n_dropped_normalizer",
      "n_initial_dropped_normalizer", "n_validation_dropped_normalizer",
      "n_refinement_steps"
    )
  } else {
    c(
      "log_normalizer", "support_grid_normalization_integral",
      "n_dropped_weight", "weight_partitions", "n_weight_fallbacks",
      "n_weight_fallback_rows", "weight_fallback_from",
      "weight_fallback_reasons"
    )
  }
  .iwmde_validate_required_fields(
    density,
    estimator_fields,
    paste(density[["estimator"]], "density result")
  )
  if (identical(density[["estimator"]], "q_grid_cmde")) {
    qcmde_vector_fields <- c(
      "pilot_y", "validation_y", "ordinate_relative_change",
      "ordinate_log_change", "pilot_ordinate_relative_change",
      "pilot_ordinate_log_change"
    )
    if (!all(vapply(density[qcmde_vector_fields], function(value) {
      is.numeric(value) && length(value) == n_values
    }, logical(1)))) {
      stop(
        "Internal qCMDE density result has inconsistent diagnostic vectors.",
        call. = FALSE
      )
    }
    if (any(!is.finite(density[["pilot_y"]])) ||
        any(density[["pilot_y"]] < 0) ||
        any(!is.finite(density[["validation_y"]])) ||
        any(density[["validation_y"]] < 0)) {
      stop("Internal qCMDE density result has invalid comparison densities.",
           call. = FALSE)
    }
    qcmde_change_fields <- c(
      "ordinate_relative_change", "ordinate_log_change",
      "pilot_ordinate_relative_change", "pilot_ordinate_log_change"
    )
    if (!all(vapply(density[qcmde_change_fields], function(value) {
      all(is.na(value) | (is.finite(value) & value >= 0))
    }, logical(1)))) {
      stop("Internal qCMDE density result has invalid comparison changes.",
           call. = FALSE)
    }
    if (!is.numeric(density[["log_normalizer"]]) ||
        !is.numeric(density[["pilot_log_normalizer"]]) ||
        length(density[["log_normalizer"]]) !=
          density[["n_normalized_rows"]] ||
        length(density[["pilot_log_normalizer"]]) !=
          density[["n_normalized_rows"]]) {
      stop(
        "Internal qCMDE density result has invalid normalizer vectors.",
        call. = FALSE
      )
    }
    if (any(!is.finite(density[["log_normalizer"]]))) {
      stop("Internal qCMDE density result has a non-finite final normalizer.",
           call. = FALSE)
    }
    qcmde_count_fields <- c(
      "normalization_initial_points", "n_rescued_normalizer",
      "n_dropped_normalizer", "n_initial_dropped_normalizer",
      "n_validation_dropped_normalizer", "n_refinement_steps"
    )
    if (!all(vapply(density[qcmde_count_fields], function(value) {
      is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value >= 0 && value == as.integer(value)
    }, logical(1)))) {
      stop("Internal qCMDE density result has an invalid count.", call. = FALSE)
    }
    qcmde_scalar_fields <- c(
      "pilot_normalization_integral", "final_normalization_integral",
      "p95_normalizer_relative_change"
    )
    if (!all(vapply(
      density[qcmde_scalar_fields],
      .iwmde_schema_nonnegative_scalar,
      logical(1)
    )) ||
        !.iwmde_schema_numeric_range(
          density[["normalization_initial_range"]],
          allow_missing = FALSE
        )) {
      stop("Internal qCMDE density result has invalid normalization metadata.",
           call. = FALSE)
    }
  } else {
    iwmde_count_fields <- c(
      "n_dropped_weight", "n_weight_fallbacks", "n_weight_fallback_rows"
    )
    if (!all(vapply(density[iwmde_count_fields], function(value) {
      is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value >= 0 && value == as.integer(value)
    }, logical(1)))) {
      stop("Internal IWMDE density result has an invalid fallback count.",
           call. = FALSE)
    }
    if (!is.numeric(density[["log_normalizer"]]) ||
        length(density[["log_normalizer"]]) != 0L ||
        !is.list(density[["weight_partitions"]]) ||
        !is.character(density[["weight_fallback_from"]]) ||
        anyNA(density[["weight_fallback_from"]]) ||
        any(!nzchar(density[["weight_fallback_from"]])) ||
        !.iwmde_schema_nonnegative_scalar(
          density[["support_grid_normalization_integral"]]
        )) {
      stop("Internal IWMDE density result has invalid fallback metadata.",
           call. = FALSE)
    }
    fallback_reasons <- density[["weight_fallback_reasons"]]
    valid_reasons <- is.numeric(fallback_reasons) &&
      all(is.finite(fallback_reasons)) &&
      all(fallback_reasons >= 0) &&
      all(fallback_reasons == as.integer(fallback_reasons)) &&
      (length(fallback_reasons) == 0L ||
        (!is.null(names(fallback_reasons)) &&
          all(!is.na(names(fallback_reasons))) &&
          all(nzchar(names(fallback_reasons)))))
    if (!valid_reasons ||
        sum(fallback_reasons) != density[["n_weight_fallbacks"]] ||
        density[["n_weight_fallback_rows"]] >
          density[["n_candidate_rows"]]) {
      stop("Internal IWMDE density result has inconsistent fallback metadata.",
           call. = FALSE)
    }
  }

  invisible(density)
}


.iwmde_schema_nonnegative_scalar <- function(value) {

  return(
    is.numeric(value) && length(value) == 1L &&
      (is.na(value) || (is.finite(value) && value >= 0))
  )
}


.iwmde_schema_numeric_range <- function(value, allow_missing = TRUE) {

  if (!is.numeric(value) || length(value) != 2L) {
    return(FALSE)
  }
  if (all(is.na(value))) {
    return(isTRUE(allow_missing))
  }

  return(all(is.finite(value)) && value[[1L]] <= value[[2L]])
}


.iwmde_new_density_result <- function(fields, estimator = NULL) {

  fields[["schema_version"]] <- .iwmde_schema_version()
  class(fields) <- c("iwmde_density_result", "list")
  .iwmde_validate_density_result(fields, estimator = estimator)

  return(fields)
}


.iwmde_new_diagnostic <- function(fields) {

  fields[["schema_version"]] <- .iwmde_schema_version()
  .iwmde_validate_required_fields(
    fields,
    c("schema_version", "parameter", "status"),
    "IWMDE diagnostic"
  )
  if (!identical(fields[["schema_version"]], .iwmde_schema_version())) {
    stop("Internal IWMDE diagnostic has an incompatible schema version.",
         call. = FALSE)
  }
  if (!is.character(fields[["status"]]) || length(fields[["status"]]) != 1L ||
      !fields[["status"]] %in% c("ok", "point_only", "unsupported")) {
    stop("Internal IWMDE diagnostic has an unknown status.", call. = FALSE)
  }
  if (!is.character(fields[["parameter"]]) ||
      length(fields[["parameter"]]) != 1L || is.na(fields[["parameter"]])) {
    stop("Internal IWMDE diagnostic has an invalid structural field.",
         call. = FALSE)
  }
  if (identical(fields[["status"]], "ok")) {
    .iwmde_validate_required_fields(
      fields,
      c("target_key", "iwmde", "diagnostics"),
      "successful IWMDE diagnostic"
    )
    if (!is.character(fields[["target_key"]]) ||
        length(fields[["target_key"]]) != 1L ||
        is.na(fields[["target_key"]]) || !is.list(fields[["iwmde"]]) ||
        !is.list(fields[["diagnostics"]])) {
      stop("Internal IWMDE diagnostic has an invalid structural field.",
           call. = FALSE)
    }
  } else if (identical(fields[["status"]], "point_only")) {
    .iwmde_validate_required_fields(
      fields,
      c("samples", "point_masses", "reason"),
      "point-only IWMDE diagnostic"
    )
  } else {
    .iwmde_validate_required_fields(
      fields,
      "reason",
      "unsupported IWMDE diagnostic"
    )
  }
  class(fields) <- c("iwmde_parameter_diagnostic", "list")

  return(fields)
}
