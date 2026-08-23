.hypothesis_brma_attach_iwmde <- function(object, posterior, parameter,
                                          parameter_label, hypothesis,
                                          conditional, n_points, samples,
                                          target_relative_mcse,
                                          normalization_points,
                                          normalization_prob, density_method,
                                          n_samples,
                                          parameter_spec = NULL,
                                          level_parameter_specs = list(),
                                          workspace = NULL) {

  point_refs <- .hypothesis_brma_point_refs(hypothesis, parameter)
  if (nrow(point_refs) == 0L) {
    return(posterior)
  }
  posterior <- .hypothesis_brma_keep_requested_ordinates(
    posterior  = posterior,
    point_refs = point_refs
  )

  .iwmde_check_point_ordinate_supported(object, density_method)

  sample_parameter <- .as_mixed_posteriors_parameters(object, parameter)
  raw_samples <- .brma_as_mixed_posteriors(
    object           = object,
    parameters       = sample_parameter,
    conditional      = conditional,
    conditional_rule = "AND",
    transform_scaled = FALSE,
    n_prior_samples  = n_samples
  )
  raw_sample_parameter <- .plot_brma_density_sample_parameter(
    samples          = raw_samples,
    parameter        = parameter,
    sample_parameter = sample_parameter
  )
  raw_posterior <- BayesTools::marginal_posterior(
    samples       = raw_samples,
    parameter     = raw_sample_parameter,
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = n_samples
  )

  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }

  resources      <- .iwmde_workspace_resources(object, workspace)
  context        <- resources[["context"]]
  estimate_cache <- resources[["estimate_cache"]]

  scalar_rows <- is.na(point_refs[["level"]])
  if (any(scalar_rows)) {
    posterior <- .hypothesis_brma_attach_iwmde_scalar(
      posterior                = posterior,
      raw_posterior            = raw_posterior,
      context                  = context,
      estimate_cache           = estimate_cache,
      parameter                = parameter,
      parameter_label          = parameter_label,
      value                    = unique(point_refs[["value"]][scalar_rows]),
      conditional              = conditional,
      n_points                 = n_points,
      samples                  = samples,
      target_relative_mcse     = target_relative_mcse,
      normalization_points     = normalization_points,
      normalization_prob       = normalization_prob,
      density_method           = density_method,
      parameter_spec           = parameter_spec
    )
  }

  for (i in which(!scalar_rows)) {
    ref <- point_refs[i, , drop = FALSE]
    posterior <- .hypothesis_brma_attach_iwmde_level(
      posterior            = posterior,
      raw_posterior        = raw_posterior,
      context              = context,
      estimate_cache       = estimate_cache,
      parameter            = parameter,
      level                = ref[["level"]],
      value                = ref[["value"]],
      conditional          = conditional,
      n_points             = n_points,
      samples              = samples,
      target_relative_mcse = target_relative_mcse,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      density_method       = density_method,
      parameter_spec       = level_parameter_specs[[ref[["level"]]]]
    )
  }

  return(posterior)
}


.hypothesis_brma_keep_requested_ordinates <- function(posterior, point_refs) {

  if (!is.list(posterior)) {
    values <- point_refs[["value"]][is.na(point_refs[["level"]])]
    attr(posterior, "posterior_ordinate") <-
      .iwmde_posterior_ordinate_keep_values(
        attr(posterior, "posterior_ordinate", exact = TRUE),
        values
      )
    return(posterior)
  }

  attr(posterior, "posterior_ordinate") <- NULL
  for (level in names(posterior)) {
    values <- point_refs[["value"]][
      !is.na(point_refs[["level"]]) & point_refs[["level"]] == level
    ]
    attr(posterior[[level]], "posterior_ordinate") <-
      .iwmde_posterior_ordinate_keep_values(
        attr(posterior[[level]], "posterior_ordinate", exact = TRUE),
        values
      )
  }

  return(posterior)
}


.hypothesis_brma_attach_iwmde_scalar <- function(
    posterior, raw_posterior, context, estimate_cache, parameter,
    parameter_label, value, conditional, n_points, samples,
    target_relative_mcse, normalization_points,
    normalization_prob, density_method, parameter_spec = NULL,
    display_transform = NULL) {

  if (is.list(raw_posterior) || is.list(posterior)) {
    stop(
      "qCMDE/IWMDE point hypotheses for factor parameters must specify ",
      "a level, e.g. '", parameter_label, "[level] = ", value, "'.",
      call. = FALSE
    )
  }

  raw_values <- as.numeric(raw_posterior)
  display_values <- as.numeric(posterior)
  if (is.null(parameter_spec) &&
      !.plot_brma_same_sample_scale(raw_values, display_values)) {
    stop("Could not align qCMDE/IWMDE ordinate to the displayed coefficient scale.",
         call. = FALSE)
  }
  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "unsupported_formula_transform")) {
    stop(parameter_spec[["reason"]], call. = FALSE)
  }
  if (is.null(parameter_spec)) {
    parameter_spec <- list(
      type          = "primitive",
      prior_density = attr(raw_posterior, "prior_density", exact = TRUE)
    )
  } else if (is.null(parameter_spec[["prior_density"]])) {
    parameter_spec[["prior_density"]] <- attr(
      raw_posterior,
      "prior_density",
      exact = TRUE
    )
  }
  parameter_spec[["conditional"]]      <- conditional
  parameter_spec[["conditional_rule"]] <- "AND"
  source_value <- if (is.null(display_transform)) {
    value
  } else {
    BayesTools::parameter_transform_inverse(value, display_transform)
  }
  source_keys <- vapply(source_value, .iwmde_key_number, character(1))

  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = parameter,
    density_method  = density_method,
    density_control = list(
      n_points             = n_points,
      samples              = samples,
      target_relative_mcse = target_relative_mcse,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      display_grid         = "ordinate"
    ),
    outputs        = "ordinate",
    values         = source_value,
    parameter_spec = parameter_spec,
    metadata       = .iwmde_posterior_metadata(
      samples   = posterior,
      parameter = parameter
    ),
    cache          = estimate_cache
  )
  diagnostic <- estimate[["diagnostics"]][["ordinate"]]
  marginal_parameter <- attr(posterior, "parameter", exact = TRUE)
  values <- .iwmde_sorted_ordinate_values(source_value)
  for (requested_source_value in values) {
    requested_index <- match(
      .iwmde_key_number(requested_source_value),
      source_keys
    )
    requested_value <- value[[requested_index]]
    ordinate <- .iwmde_posterior_ordinate_keep_values(
      posterior_ordinate = estimate[["posterior_ordinate"]],
      values             = requested_source_value
    )
    if (is.null(ordinate)) {
      .iwmde_stop_ordinate_unavailable(
        message = .hypothesis_brma_iwmde_ordinate_failure_message(
          density_method = density_method,
          target         = paste0(parameter, " = ", requested_value),
          diagnostic     = diagnostic,
          reason         = .hypothesis_brma_estimate_ordinate_reason(
            estimate = estimate,
            value    = requested_source_value
          )
        ),
        estimate = estimate
      )
    }
    if (!is.null(display_transform)) {
      ordinate <- .hypothesis_brma_transform_iwmde_ordinate(
        ordinate,
        display_transform
      )
      ordinate[["value"]] <- requested_value
    }
    if (is.character(marginal_parameter) &&
        length(marginal_parameter) == 1L &&
        !is.na(marginal_parameter) && nzchar(marginal_parameter)) {
      ordinate[["parameter"]] <- marginal_parameter
    }

    existing <- .iwmde_posterior_ordinate_drop_value(
      posterior_ordinate = attr(posterior, "posterior_ordinate", exact = TRUE),
      value              = requested_value
    )
    attr(posterior, "posterior_ordinate") <-
      BayesTools::posterior_ordinate_append(
        existing = existing,
        ordinate = ordinate
      )
  }

  return(posterior)
}


.hypothesis_brma_transform_iwmde_ordinate <- function(
    ordinate, display_transform) {

  source_evaluation_value <- ordinate[["evaluation_value"]]
  jacobian <- BayesTools::parameter_transform_jacobian(
    source_evaluation_value,
    display_transform
  )
  if (length(jacobian) != 1L || !is.finite(jacobian) || jacobian <= 0) {
    stop("Unsupported qCMDE/IWMDE display transform.", call. = FALSE)
  }
  ordinate[["value"]] <- BayesTools::parameter_transform_forward(
    ordinate[["value"]],
    display_transform
  )
  ordinate[["evaluation_value"]] <- BayesTools::parameter_transform_forward(
    ordinate[["evaluation_value"]],
    display_transform
  )
  ordinate[["ordinate"]] <- ordinate[["ordinate"]] / jacobian
  diagnostics <- ordinate[["diagnostics"]]
  value_fields <- intersect(
    c("evaluation_value"),
    names(diagnostics)
  )
  for (field in value_fields) {
    diagnostics[[field]] <- BayesTools::parameter_transform_forward(
      diagnostics[[field]],
      display_transform
    )
  }
  density_fields <- intersect(
    c(
      "mcse", "active_branch_mcse", "sampling_mcse",
      "pilot_ordinate", "validation_ordinate"
    ),
    names(diagnostics)
  )
  for (field in density_fields) {
    diagnostics[[field]] <- diagnostics[[field]] / jacobian
  }
  ordinate[["diagnostics"]] <- diagnostics
  provenance <- ordinate[["iwmde_provenance"]]
  if (is.list(provenance)) {
    for (field in intersect(c("value", "evaluation_value"), names(provenance))) {
      provenance[[field]] <- BayesTools::parameter_transform_forward(
        provenance[[field]],
        display_transform
      )
    }
    ordinate[["iwmde_provenance"]] <- provenance
  }

  ordinate
}


.hypothesis_brma_estimate_ordinate_reason <- function(estimate, value) {

  rejected <- .iwmde_posterior_ordinate_keep_values(
    posterior_ordinate = estimate[["rejected_posterior_ordinate"]],
    values             = value
  )
  if (!is.null(rejected)) {
    reason <- .iwmde_posterior_ordinate_bf_failure_reason(rejected)
    if (!is.null(reason) && nzchar(reason)) {
      return(reason)
    }
  }

  return(.hypothesis_brma_diagnostic_reason(
    estimate[["diagnostics"]][["ordinate"]]
  ))
}


.hypothesis_brma_attach_iwmde_level <- function(
    posterior, raw_posterior, context, estimate_cache, parameter, level,
    value, conditional, n_points, samples,
    target_relative_mcse, normalization_points, normalization_prob,
    density_method, parameter_spec = NULL) {

  if (!is.list(posterior) || !level %in% names(posterior)) {
    stop("Hypothesis references unknown level '", level,
         "' for parameter '", parameter, "'.", call. = FALSE)
  }
  if (!is.list(raw_posterior) || !level %in% names(raw_posterior)) {
    stop("Raw posterior for level '", level, "' is unavailable.",
         call. = FALSE)
  }
  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "unsupported_formula_transform")) {
    stop(parameter_spec[["reason"]], call. = FALSE)
  }

  weights <- if (is.null(parameter_spec)) {
    attr(raw_posterior[[level]], "linear_weights", exact = TRUE)
  } else {
    parameter_spec[["weights"]]
  }
  linear_weights <- .iwmde_linear_weights(weights)
  if (is.null(linear_weights) || length(linear_weights) == 0L) {
    reason <- if (is.null(linear_weights)) {
      "linear weights are unavailable"
    } else {
      "linear weights are all zero"
    }
    stop(
      "Precomputed ", density_method, " posterior ordinate is unavailable for '",
      parameter, "[", level, "] = ", value, "': ", reason,
      call. = FALSE
    )
  }

  raw_values     <- as.numeric(raw_posterior[[level]])
  display_values <- as.numeric(posterior[[level]])
  if (is.null(parameter_spec) &&
      !.plot_brma_same_sample_scale(raw_values, display_values)) {
    stop(
      "Could not align qCMDE/IWMDE ordinate for '", parameter, "[", level,
      "]' to the displayed scale.",
      call. = FALSE
    )
  }
  if (is.null(parameter_spec)) {
    parameter_spec <- list(
      type          = "linear",
      weights       = weights,
      prior_density = attr(
        raw_posterior[[level]],
        "prior_density",
        exact = TRUE
      )
    )
  }
  parameter_spec[["conditional"]] <- if (is.null(conditional)) {
    attr(raw_posterior[[level]], "effective_conditional", exact = TRUE)
  } else {
    conditional
  }
  parameter_spec[["conditional_rule"]] <- "AND"

  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = paste0(parameter, "[", level, "]"),
    density_method  = density_method,
    density_control = list(
      n_points             = n_points,
      samples              = samples,
      target_relative_mcse = target_relative_mcse,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      display_grid         = "ordinate"
    ),
    outputs        = "ordinate",
    values         = value,
    parameter_spec = parameter_spec,
    metadata       = .iwmde_posterior_metadata(
      samples   = posterior[[level]],
      parameter = parameter,
      level     = level
    ),
    cache          = estimate_cache
  )
  diagnostic <- estimate[["diagnostics"]][["ordinate"]]

  ordinate <- estimate[["posterior_ordinate"]]
  if (is.null(ordinate)) {
    .iwmde_stop_ordinate_unavailable(
      message = .hypothesis_brma_iwmde_ordinate_failure_message(
        density_method = density_method,
        target         = paste0(parameter, "[", level, "] = ", value),
        diagnostic     = diagnostic,
        reason         = .hypothesis_brma_diagnostic_reason(diagnostic)
      ),
      estimate = estimate
    )
  }

  existing <- .iwmde_posterior_ordinate_drop_value(
    posterior_ordinate = attr(posterior[[level]], "posterior_ordinate",
                              exact = TRUE),
    value              = value
  )
  attr(posterior[[level]], "posterior_ordinate") <- BayesTools::posterior_ordinate_append(
    existing = existing,
    ordinate = ordinate
  )

  return(posterior)
}


.hypothesis_brma_iwmde_ordinate_failure_message <- function(
    density_method, target, diagnostic, reason) {

  if (identical(diagnostic[["status"]], "ok")) {
    return(paste0(
      density_method, " posterior ordinate for '", target,
      "' was rejected by diagnostics: ", reason
    ))
  }

  return(paste0(
    density_method, " posterior ordinate is unavailable for '", target,
    "': ", reason
  ))
}


.hypothesis_brma_diagnostic_reason <- function(diagnostic) {

  reason <- diagnostic[["reason"]]
  if (is.null(reason) || !nzchar(reason)) {
    diagnostic_reason <- .iwmde_diagnostics_bf_failure_reason(
      diagnostic[["diagnostics"]]
    )
    if (!is.null(diagnostic_reason) && nzchar(diagnostic_reason)) {
      reason <- diagnostic_reason
    } else {
      reason <- "failed qCMDE/IWMDE numerical diagnostics"
    }
  }

  return(reason)
}


.hypothesis_brma_append_iwmde_warnings <- function(
    table, posterior, parameter = NULL) {

  diagnostics <- .iwmde_collect_public_density_diagnostics(posterior)
  if (nrow(diagnostics) > 0L) {
    attr(table, "density_diagnostics") <- diagnostics
  }

  warning_records <- .iwmde_collect_posterior_ordinate_warning_records(posterior)
  if (nrow(warning_records) == 0L) {
    return(table)
  }
  if (!is.null(parameter)) {
    warning_records[["parameter"]] <- parameter
  }

  existing <- attr(table, "warnings", exact = TRUE)
  warnings <- unlist(lapply(seq_len(nrow(warning_records)), function(i) {

    record <- warning_records[i, , drop = FALSE]
    rows   <- .hypothesis_brma_warning_rows(table, record)
    stats::setNames(rep(record[["warning"]], length(rows)), rows)
  }), use.names = TRUE)
  attr(table, "warnings") <- .hypothesis_brma_unique_named_warnings(
    c(existing, warnings)
  )

  return(table)
}


.hypothesis_brma_warning_rows <- function(table, record) {

  row_names <- rownames(table)
  if (is.null(row_names) || length(row_names) == 0L) {
    return(.hypothesis_brma_warning_row(table))
  }

  candidates <- seq_len(nrow(table))

  point_rows <- .hypothesis_brma_warning_point_rows(table, record)
  if (length(point_rows) > 0L) {
    candidates <- intersect(candidates, point_rows)
  }

  label <- record[["label"]]
  if (!is.na(label) && nzchar(label)) {
    label_rows <- unique(c(
      which(row_names == label),
      .hypothesis_brma_table_rows_with_labels(table, label)
    ))
    narrowed <- intersect(candidates, label_rows)
    if (length(narrowed) > 0L) {
      candidates <- narrowed
    }
  }

  if (length(candidates) == 0L) {
    return(.hypothesis_brma_warning_row(table))
  }

  row_names[candidates]
}


.hypothesis_brma_warning_point_rows <- function(table, record) {

  statements <- .hypothesis_brma_result_statements(table)
  if (is.null(statements) || length(statements) == 0L ||
      nrow(table) %% length(statements) != 0L) {
    point_labels <- .hypothesis_brma_warning_point_labels(table, record)
    return(.hypothesis_brma_table_rows_with_labels(table, point_labels))
  }

  statement_matches <- vapply(statements, function(statement) {

    any(vapply(
      list(statement[["left"]], statement[["right"]]),
      .hypothesis_brma_warning_side_matches,
      logical(1),
      record = record
    ))
  }, logical(1))
  if (!any(statement_matches)) {
    return(integer())
  }

  n_quantities <- nrow(table) %/% length(statements)
  rows <- unlist(lapply(which(statement_matches), function(i) {

    ((i - 1L) * n_quantities + 1L):(i * n_quantities)
  }), use.names = FALSE)

  rows
}


.hypothesis_brma_warning_point_labels <- function(table, record) {

  statements <- .hypothesis_brma_result_statements(table)
  if (is.null(statements) || length(statements) == 0L) {
    return(character())
  }

  labels <- unlist(lapply(statements, function(statement) {

    sides <- list(statement[["left"]], statement[["right"]])
    vapply(sides, function(side) {

      if (!.hypothesis_brma_warning_side_matches(side, record)) {
        return(NA_character_)
      }
      side[["label"]]
    }, character(1))
  }), use.names = FALSE)
  labels <- unique(labels[!is.na(labels) & nzchar(labels)])

  return(labels)
}


.hypothesis_brma_warning_side_matches <- function(side, record) {

  if (!is.list(side) || !identical(side[["type"]], "point")) {
    return(FALSE)
  }

  expression <- side[["expression"]]
  while (is.list(expression) &&
         identical(expression[["type"]], "parentheses")) {
    expression <- expression[["expression"]]
  }
  expr <- if (is.list(expression) &&
              identical(expression[["type"]], "symbol")) {
    expression[["name"]]
  } else if (is.list(expression) &&
             identical(expression[["type"]], "level_reference")) {
    paste0(expression[["parameter"]], "[", expression[["level"]], "]")
  } else {
    NULL
  }
  if (is.null(expr) || length(expr) != 1L || is.na(expr) || !nzchar(expr)) {
    return(FALSE)
  }

  label     <- record[["label"]]
  parameter <- record[["parameter"]]
  terms     <- unique(c(label, parameter))
  terms     <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) > 0L && !expr %in% terms) {
    return(FALSE)
  }

  value <- suppressWarnings(as.numeric(record[["value"]]))
  if (!is.finite(value)) {
    return(TRUE)
  }
  side_value <- suppressWarnings(as.numeric(side[["value"]]))
  if (!is.finite(side_value)) {
    return(FALSE)
  }

  identical(value, side_value)
}


.hypothesis_brma_table_rows_with_labels <- function(table, labels) {

  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) == 0L) {
    return(integer())
  }

  columns <- intersect(c("Alternative", "Null"), colnames(table))
  if (length(columns) == 0L) {
    return(integer())
  }

  matches <- logical(nrow(table))
  for (column in columns) {
    values <- as.character(table[[column]])
    matches <- matches | values %in% labels
    for (label in labels) {
      matches <- matches | grepl(label, values, fixed = TRUE)
    }
  }

  which(matches)
}


.hypothesis_brma_warning_row <- function(table) {

  row_names <- rownames(table)
  if (!is.null(row_names) && length(row_names) > 0L) {
    return(row_names[[1L]])
  }

  return("qCMDE/IWMDE")
}


.hypothesis_brma_unique_named_warnings <- function(warnings) {

  if (length(warnings) == 0L) {
    return(NULL)
  }

  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  if (length(warnings) == 0L) {
    return(NULL)
  }

  warning_names <- names(warnings)
  if (is.null(warning_names)) {
    warning_names <- rep("", length(warnings))
  }
  key <- paste(warning_names, warnings, sep = "\r")

  warnings[!duplicated(key)]
}


.iwmde_collect_posterior_ordinate_warning_records <- function(posterior) {

  records <- .iwmde_posterior_ordinate_warning_records(
    attr(posterior, "posterior_ordinate", exact = TRUE)
  )
  if (is.list(posterior)) {
    child_records <- lapply(posterior, function(sample) {
      .iwmde_collect_posterior_ordinate_warning_records(sample)
    })
    child_records <- child_records[vapply(child_records, nrow, integer(1)) > 0L]
    if (length(child_records) > 0L) {
      records <- rbind(records, do.call(rbind, child_records))
    }
  }

  records <- records[!duplicated(records), , drop = FALSE]
  rownames(records) <- NULL

  return(records)
}


.iwmde_posterior_ordinate_warning_records <- function(posterior_ordinate) {

  entries <- .iwmde_posterior_ordinate_entries(posterior_ordinate)
  if (length(entries) == 0L) {
    return(.iwmde_empty_posterior_ordinate_warning_records())
  }

  records <- lapply(entries, function(entry) {

    warnings <- .iwmde_posterior_ordinate_warnings(entry)
    if (length(warnings) == 0L) {
      return(.iwmde_empty_posterior_ordinate_warning_records())
    }

    data.frame(
      warning   = warnings,
      parameter = .hypothesis_brma_ordinate_character(entry, "parameter"),
      level     = .hypothesis_brma_ordinate_character(entry, "level"),
      label     = .hypothesis_brma_ordinate_label(entry),
      value     = .hypothesis_brma_ordinate_value(entry),
      stringsAsFactors = FALSE
    )
  })
  records <- records[vapply(records, nrow, integer(1)) > 0L]
  if (length(records) == 0L) {
    return(.iwmde_empty_posterior_ordinate_warning_records())
  }

  do.call(rbind, records)
}


.iwmde_empty_posterior_ordinate_warning_records <- function() {

  data.frame(
    warning   = character(),
    parameter = character(),
    level     = character(),
    label     = character(),
    value     = numeric(),
    stringsAsFactors = FALSE
  )
}


.hypothesis_brma_ordinate_character <- function(entry, name) {

  value <- entry[[name]]
  if (is.null(value) || length(value) == 0L) {
    return(NA_character_)
  }

  value <- as.character(value[[1L]])
  if (is.na(value) || !nzchar(value)) {
    return(NA_character_)
  }

  return(value)
}


.hypothesis_brma_ordinate_label <- function(entry) {

  parameter <- .hypothesis_brma_ordinate_character(entry, "parameter")
  level     <- .hypothesis_brma_ordinate_character(entry, "level")
  if (!is.na(parameter) && nzchar(parameter) &&
      !is.na(level) && nzchar(level)) {
    return(paste0(parameter, "[", level, "]"))
  }
  if (!is.na(parameter) && nzchar(parameter)) {
    return(parameter)
  }

  return(NA_character_)
}


.hypothesis_brma_ordinate_value <- function(entry) {

  for (name in c("value", "x", "null_hypothesis")) {
    value <- suppressWarnings(as.numeric(entry[[name]]))
    value <- value[is.finite(value)]
    if (length(value) > 0L) {
      return(value[[1L]])
    }
  }

  return(NA_real_)
}
