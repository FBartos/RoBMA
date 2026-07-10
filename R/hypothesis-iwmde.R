.hypothesis_brma_attach_iwmde <- function(object, posterior, parameter,
                                          parameter_label, hypothesis,
                                          conditional, n_points, max_samples,
                                          normalization_points,
                                          normalization_prob, density_method,
                                          n_samples) {

  point_refs <- .hypothesis_brma_point_refs(hypothesis, parameter)
  if (nrow(point_refs) == 0L) {
    return(posterior)
  }

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

  context        <- .iwmde_context(object)
  estimate_cache <- .iwmde_estimate_cache()

  for (i in seq_len(nrow(point_refs))) {
    ref <- point_refs[i, , drop = FALSE]
    if (is.na(ref[["level"]])) {
      posterior <- .hypothesis_brma_attach_iwmde_scalar(
        posterior                = posterior,
        raw_posterior            = raw_posterior,
        context                  = context,
        estimate_cache           = estimate_cache,
        parameter                = parameter,
        parameter_label          = parameter_label,
        value                    = ref[["value"]],
        conditional              = conditional,
        n_points                 = n_points,
        max_samples              = max_samples,
        normalization_points     = normalization_points,
        normalization_prob       = normalization_prob,
        density_method           = density_method
      )
    } else {
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
        max_samples          = max_samples,
        normalization_points = normalization_points,
        normalization_prob   = normalization_prob,
        density_method       = density_method
      )
    }
  }

  return(posterior)
}


.hypothesis_brma_attach_iwmde_scalar <- function(
    posterior, raw_posterior, context, estimate_cache, parameter,
    parameter_label, value, conditional, n_points, max_samples, normalization_points,
    normalization_prob, density_method) {

  if (is.list(raw_posterior) || is.list(posterior)) {
    stop(
      "qCMDE/IWMDE point hypotheses for factor parameters must specify ",
      "a level, e.g. '", parameter_label, "[level] = ", value, "'.",
      call. = FALSE
    )
  }

  raw_values <- as.numeric(raw_posterior)
  display_values <- as.numeric(posterior)
  transform <- .plot_brma_affine_sample_transform(raw_values, display_values)
  if (is.null(transform)) {
    stop("Could not align qCMDE/IWMDE ordinate to the displayed coefficient scale.",
         call. = FALSE)
  }

  raw_value <- (value - transform[["intercept"]]) / transform[["slope"]]
  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = parameter,
    density_method  = density_method,
    density_control = list(
      n_points             = n_points,
      max_samples          = max_samples,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      display_grid         = "ordinate"
    ),
    outputs        = "ordinate",
    values         = raw_value,
    parameter_spec = list(
      type             = "primitive",
      conditional      = conditional,
      conditional_rule = "AND"
    ),
    metadata       = .iwmde_posterior_metadata(
      samples   = posterior,
      parameter = parameter
    ),
    cache          = estimate_cache
  )
  diagnostic <- estimate[["diagnostics"]][["ordinate"]]

  ordinate <- estimate[["posterior_ordinate"]]
  if (is.null(ordinate)) {
    stop(
      "Precomputed ", density_method, " posterior ordinate is unavailable for '",
      parameter, " = ", value, "': ",
      .hypothesis_brma_diagnostic_reason(diagnostic),
      call. = FALSE
    )
  }

  ordinate <- .hypothesis_brma_transform_ordinate(
    ordinate  = ordinate,
    transform = transform,
    value     = value
  )
  existing <- .iwmde_posterior_ordinate_drop_value(
    posterior_ordinate = attr(posterior, "posterior_ordinate", exact = TRUE),
    value              = value
  )
  attr(posterior, "posterior_ordinate") <- BayesTools::posterior_ordinate_append(
    existing = existing,
    ordinate = ordinate
  )

  return(posterior)
}


.hypothesis_brma_attach_iwmde_level <- function(
    posterior, raw_posterior, context, estimate_cache, parameter, level,
    value, conditional, n_points, max_samples, normalization_points, normalization_prob,
    density_method) {

  if (!is.list(posterior) || !level %in% names(posterior)) {
    stop("Hypothesis references unknown level '", level,
         "' for parameter '", parameter, "'.", call. = FALSE)
  }
  if (!is.list(raw_posterior) || !level %in% names(raw_posterior)) {
    stop("Raw posterior for level '", level, "' is unavailable.",
         call. = FALSE)
  }

  weights <- attr(raw_posterior[[level]], "linear_weights", exact = TRUE)
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
  transform <- .plot_brma_affine_sample_transform(raw_values, display_values)
  if (is.null(transform)) {
    stop(
      "Could not align qCMDE/IWMDE ordinate for '", parameter, "[", level,
      "]' to the displayed scale.",
      call. = FALSE
    )
  }
  raw_value <- (value - transform[["intercept"]]) / transform[["slope"]]
  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = paste0(parameter, "[", level, "]"),
    density_method  = density_method,
    density_control = list(
      n_points             = n_points,
      max_samples          = max_samples,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      display_grid         = "ordinate"
    ),
    outputs        = "ordinate",
    values         = raw_value,
    parameter_spec = list(
      type             = "linear",
      weights          = weights,
      conditional      = if (is.null(conditional)) {
        attr(raw_posterior[[level]], "effective_conditional", exact = TRUE)
      } else {
        conditional
      },
      conditional_rule = "AND"
    ),
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
    stop(
      "Precomputed ", density_method, " posterior ordinate is unavailable for '",
      parameter, "[", level, "] = ", value, "': ",
      .hypothesis_brma_diagnostic_reason(diagnostic),
      call. = FALSE
    )
  }

  ordinate <- .hypothesis_brma_transform_ordinate(
    ordinate  = ordinate,
    transform = transform,
    value     = value
  )
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


.hypothesis_brma_attach_iwmde_log_intercept <- function(
    object, posterior, parameter, hypothesis, n_points, max_samples,
    normalization_points, normalization_prob, density_method) {

  point_refs <- .hypothesis_brma_point_refs(hypothesis, parameter)
  if (nrow(point_refs) == 0L) {
    return(posterior)
  }
  if (any(point_refs[["value"]] <= 0)) {
    stop(
      "qCMDE/IWMDE point hypotheses for a log-scale intercept require a ",
      "strictly positive value.",
      call. = FALSE
    )
  }
  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }

  context        <- .iwmde_context(object)
  estimate_cache <- .iwmde_estimate_cache()
  for (i in seq_len(nrow(point_refs))) {
    value <- point_refs[["value"]][i]
    estimate <- .iwmde_estimate(
      context         = context,
      parameter       = parameter,
      density_method  = density_method,
      density_control = list(
        n_points             = n_points,
        max_samples          = max_samples,
        normalization_points = normalization_points,
        normalization_prob   = normalization_prob,
        display_grid         = "ordinate"
      ),
      outputs        = "ordinate",
      values         = value,
      parameter_spec = list(
        type             = "scale_log_intercept",
        conditional      = NULL,
        conditional_rule = "AND"
      ),
      metadata = .iwmde_posterior_metadata(
        samples   = posterior,
        parameter = parameter
      ),
      cache = estimate_cache
    )
    ordinate <- estimate[["posterior_ordinate"]]
    if (is.null(ordinate)) {
      stop(
        "Precomputed ", density_method, " posterior ordinate is unavailable for '",
        parameter, " = ", value, "': ",
        .hypothesis_brma_diagnostic_reason(
          estimate[["diagnostics"]][["ordinate"]]
        ),
        call. = FALSE
      )
    }
    existing <- .iwmde_posterior_ordinate_drop_value(
      posterior_ordinate = attr(posterior, "posterior_ordinate", exact = TRUE),
      value              = value
    )
    attr(posterior, "posterior_ordinate") <- BayesTools::posterior_ordinate_append(
      existing = existing,
      ordinate = ordinate
    )
  }

  return(posterior)
}


.hypothesis_brma_transform_ordinate <- function(ordinate, transform, value) {

  slope <- abs(transform[["slope"]])
  ordinate[["value"]]    <- value
  ordinate[["ordinate"]] <- ordinate[["ordinate"]] / slope
  if (!is.null(ordinate[["evaluation_value"]])) {
    ordinate[["evaluation_value"]] <- transform[["intercept"]] +
      transform[["slope"]] * ordinate[["evaluation_value"]]
  }

  diagnostics <- ordinate[["diagnostics"]]
  if (!is.null(diagnostics[["evaluation_value"]])) {
    diagnostics[["evaluation_value"]] <- transform[["intercept"]] +
      transform[["slope"]] * diagnostics[["evaluation_value"]]
  }
  if (!is.null(diagnostics[["mcse"]])) {
    diagnostics[["mcse"]] <- diagnostics[["mcse"]] / slope
  }
  if (!is.null(diagnostics[["normalization_range"]])) {
    diagnostics[["normalization_range"]] <- transform[["intercept"]] +
      transform[["slope"]] * diagnostics[["normalization_range"]]
    diagnostics[["normalization_range"]] <- sort(diagnostics[["normalization_range"]])
  }
  diagnostics[["plot_scale_transform"]] <- transform
  ordinate[["diagnostics"]] <- diagnostics
  provenance <- ordinate[["iwmde_provenance"]]
  if (is.list(provenance)) {
    provenance[["result"]] <- c(
      provenance[["result"]],
      list(plot_scale_transform = transform)
    )
    ordinate[["iwmde_provenance"]] <- provenance
  }

  return(ordinate)
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
      reason <- "failed BF-grade diagnostics"
    }
  }

  return(reason)
}


.hypothesis_brma_append_iwmde_warnings <- function(table, posterior) {

  warning_records <- .iwmde_collect_posterior_ordinate_warning_records(posterior)
  if (nrow(warning_records) == 0L) {
    return(table)
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

  parsed <- attr(table, "parsed", exact = TRUE)
  if (is.null(parsed) || length(parsed) == 0L ||
      nrow(table) %% length(parsed) != 0L) {
    point_labels <- .hypothesis_brma_warning_point_labels(table, record)
    return(.hypothesis_brma_table_rows_with_labels(table, point_labels))
  }

  parsed_matches <- vapply(parsed, function(item) {

    any(vapply(
      list(item[["left"]], item[["right"]]),
      .hypothesis_brma_warning_side_matches,
      logical(1),
      record = record
    ))
  }, logical(1))
  if (!any(parsed_matches)) {
    return(integer())
  }

  n_quantities <- nrow(table) %/% length(parsed)
  rows <- unlist(lapply(which(parsed_matches), function(i) {

    ((i - 1L) * n_quantities + 1L):(i * n_quantities)
  }), use.names = FALSE)

  rows
}


.hypothesis_brma_warning_point_labels <- function(table, record) {

  parsed <- attr(table, "parsed", exact = TRUE)
  if (is.null(parsed) || length(parsed) == 0L) {
    return(character())
  }

  labels <- unlist(lapply(parsed, function(item) {

    sides <- list(item[["left"]], item[["right"]])
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

  expr <- side[["expr"]]
  if (is.null(expr) || length(expr) == 0L || is.na(expr[[1L]])) {
    return(FALSE)
  }
  expr <- as.character(expr[[1L]])

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

  isTRUE(all.equal(value, side_value, tolerance = sqrt(.Machine$double.eps)))
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
