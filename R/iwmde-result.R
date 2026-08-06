# ============================================================================ #
# IWMDE Result Schema and Provenance Helpers
# ============================================================================ #

.iwmde_schema_version <- function() {

  return("4")
}


.iwmde_algorithm_version <- function() {

  return("12")
}


.iwmde_package_version <- function(package) {

  tryCatch(
    as.character(utils::packageVersion(package)),
    error = function(e) NA_character_
  )
}


.iwmde_hash <- function(prefix, payload) {

  return(paste(prefix, rlang::hash(payload), sep = "|"))
}


.iwmde_compact_nulls <- function(x) {

  if (!is.list(x)) {
    return(x)
  }

  x <- lapply(x, .iwmde_compact_nulls)
  x <- x[!vapply(x, is.null, logical(1))]

  return(x)
}


.iwmde_density_control_provenance <- function(density_control) {

  if (is.null(density_control)) {
    return(NULL)
  }

  keep <- c(
    "n_points",
    "max_samples",
    "samples",
    "target_relative_mcse",
    "normalization_points",
    "normalization_prob",
    "display_grid"
  )
  density_control <- density_control[intersect(keep, names(density_control))]
  density_control <- .iwmde_compact_nulls(density_control)

  return(density_control)
}


.iwmde_target_provenance <- function(metadata) {

  if (is.null(metadata)) {
    return(NULL)
  }

  keep <- c(
    "parameter",
    "level",
    "conditional",
    "conditional_rule",
    "condition_key"
  )
  metadata <- metadata[intersect(keep, names(metadata))]
  metadata <- .iwmde_compact_nulls(metadata)

  return(metadata)
}


.iwmde_diagnostic_policy <- function() {

  return(list(
    warn_relative_mcse = .iwmde_bf_warning_relative_mcse(),
    warn_min_ess      = .iwmde_bf_warning_min_ess(),
    warn_weight_share = .iwmde_bf_warning_weight_share(),
    warn_min_finite_terms = .iwmde_bf_warning_min_finite_terms(),
    qcmde_warn        = .iwmde_bf_mass_warning_tolerance("q_grid_cmde"),
    qcmde_fail        = .iwmde_bf_mass_fail_tolerance("q_grid_cmde"),
    iwmde_warn        = .iwmde_bf_mass_warning_tolerance("iwmde"),
    iwmde_fail        = .iwmde_bf_mass_fail_tolerance("iwmde"),
    quadrature_warn   = .iwmde_quadrature_warning_tolerance(),
    quadrature_fail   = .iwmde_quadrature_fail_tolerance()
  ))
}


.iwmde_result_method <- function(diagnostic, density_method) {

  diagnostics <- diagnostic[["diagnostics"]]
  if (is.list(diagnostics) && !is.null(diagnostics[["estimator"]])) {
    return(as.character(diagnostics[["estimator"]])[[1L]])
  }

  density <- diagnostic[["iwmde"]]
  if (is.list(density) && !is.null(density[["estimator"]])) {
    return(as.character(density[["estimator"]])[[1L]])
  }

  return(.density_method_iwmde_estimator(density_method))
}


.iwmde_provenance_request <- function(density_method, method, metadata = NULL,
                                      density_control = NULL,
                                      value = NULL,
                                      attribute = c("density", "ordinate"),
                                      target_key = NULL,
                                      source_fingerprint = NULL,
                                      plan_key = NULL,
                                      prior_ordinates = NULL) {

  attribute <- match.arg(attribute)
  density_method <- .density_method_normalize(density_method)

  request <- list(
    schema_version    = .iwmde_schema_version(),
    algorithm_version = .iwmde_algorithm_version(),
    provenance_level  = "diagnostic_adapter",
    attribute         = attribute,
    density_method    = density_method,
    method            = method,
    internal_method   = method,
    density_control   = .iwmde_density_control_provenance(density_control),
    target            = .iwmde_target_provenance(metadata),
    target_key        = target_key,
    source_fingerprint = source_fingerprint,
    prior_ordinates   = prior_ordinates,
    requested_value   = if (is.null(value)) {
      NULL
    } else {
      .iwmde_key_number(value)
    }
  )
  request <- .iwmde_compact_nulls(request)

  request[["request_key"]] <- .iwmde_hash("iwmde_request", request)
  if (!is.null(plan_key)) {
    request[["plan_key"]] <- plan_key
  }

  return(request)
}


.iwmde_result_provenance <- function(diagnostic, density_method,
                                     metadata = NULL,
                                     density_control = NULL,
                                     value = NULL,
                                     evaluation_value = NULL,
                                     attribute = c("density", "ordinate")) {

  attribute <- match.arg(attribute)
  method    <- .iwmde_result_method(diagnostic, density_method)
  plan      <- diagnostic[["plan"]]
  source_fingerprint <- if (is.list(plan)) {
    plan[["source_fingerprint"]]
  } else {
    NULL
  }
  target_key <- if (is.list(plan)) {
    plan[["target"]][["target_key"]]
  } else {
    diagnostic[["target_key"]]
  }
  prior_ordinates <- if (identical(attribute, "ordinate") && is.list(plan)) {
    plan[["prior_ordinates"]]
  } else {
    NULL
  }
  provenance <- .iwmde_provenance_request(
    density_method   = density_method,
    method           = method,
    metadata         = metadata,
    density_control  = density_control,
    value            = value,
    attribute        = attribute,
    target_key       = target_key,
    source_fingerprint = source_fingerprint,
    plan_key           = if (is.list(plan)) plan[["plan_key"]] else NULL,
    prior_ordinates    = prior_ordinates
  )
  if (is.list(plan)) {
    provenance[["provenance_level"]] <- "iwmde_plan"
    provenance[["plan_key"]]         <- plan[["plan_key"]]
    if (identical(attribute, "ordinate")) {
      provenance[["prior_ordinates"]] <- prior_ordinates
    }
  }
  diagnostics <- diagnostic[["diagnostics"]]
  bf_grade <- if (identical(attribute[[1L]], "ordinate") &&
                  is.list(diagnostics)) {
    is.null(.iwmde_diagnostics_bf_failure_reason(diagnostics))
  } else {
    NA
  }

  provenance[["RoBMA_version"]]      <- .iwmde_package_version("RoBMA")
  provenance[["BayesTools_version"]] <- .iwmde_package_version("BayesTools")
  provenance[["diagnostic_policy"]]  <- .iwmde_diagnostic_policy()
  provenance[["result"]] <- .iwmde_compact_nulls(list(
    target_key       = diagnostic[["target_key"]],
    evaluation_value = if (is.null(evaluation_value)) {
      NULL
    } else {
      as.numeric(evaluation_value)
    },
    evaluation_value_key = if (is.null(evaluation_value)) {
      NULL
    } else {
      .iwmde_key_number(evaluation_value)
    },
    bf_grade         = bf_grade
  ))

  return(provenance)
}


.iwmde_attribute_provenance <- function(x) {

  if (is.null(x) || !is.list(x)) {
    return(NULL)
  }

  provenance <- x[["iwmde_provenance"]]
  if (!is.list(provenance) || is.null(provenance[["request_key"]])) {
    return(NULL)
  }

  return(provenance)
}


.iwmde_provenance_matches <- function(existing, requested) {

  existing <- .iwmde_attribute_provenance(existing)
  if (is.null(existing) || is.null(requested)) {
    return(FALSE)
  }

  identical(existing[["request_key"]], requested[["request_key"]])
}


.iwmde_attribute_method_matches_request <- function(attribute, provenance) {

  if (is.null(attribute) || !is.list(attribute) || is.null(provenance)) {
    return(FALSE)
  }
  if (!is.null(attribute[["density_method"]])) {
    density_method <- tryCatch(
      .density_method_normalize(attribute[["density_method"]]),
      error = function(e) NULL
    )
    if (!identical(density_method, provenance[["density_method"]])) {
      return(FALSE)
    }
  }
  if (!is.null(attribute[["method"]]) &&
      !identical(
        as.character(attribute[["method"]])[[1L]],
        provenance[["method"]]
      )) {
    return(FALSE)
  }

  return(TRUE)
}


.iwmde_posterior_density_matches_request <- function(posterior_density,
                                                     provenance) {

  if (is.null(posterior_density) || !is.list(posterior_density)) {
    return(FALSE)
  }
  if (!.iwmde_provenance_matches(posterior_density, provenance)) {
    return(FALSE)
  }
  if (!.iwmde_attribute_method_matches_request(posterior_density, provenance)) {
    return(FALSE)
  }
  diagnostics <- posterior_density[["diagnostics"]]
  if (!is.null(.iwmde_diagnostics_density_failure_reason(diagnostics))) {
    return(FALSE)
  }
  if (is.null(posterior_density[["x"]]) ||
      is.null(posterior_density[["y"]])) {
    return(FALSE)
  }

  x <- suppressWarnings(as.numeric(posterior_density[["x"]]))
  y <- suppressWarnings(as.numeric(posterior_density[["y"]]))

  return(
    length(x) == length(y) &&
      length(x) > 1L &&
      all(is.finite(x)) &&
      all(is.finite(y))
  )
}


.iwmde_ordinate_value_matches <- function(entry, value) {

  entry_value <- .iwmde_ordinate_scalar(entry, "value")
  value       <- suppressWarnings(as.numeric(value))[1]
  if (!is.finite(entry_value) || !is.finite(value)) {
    return(FALSE)
  }

  return(entry_value == value)
}


.iwmde_posterior_ordinate_matches_request <- function(posterior_ordinate,
                                                       value,
                                                       provenance) {

  entries <- .iwmde_posterior_ordinate_entries(posterior_ordinate)
  if (length(entries) == 0L) {
    return(FALSE)
  }

  any(vapply(entries, function(entry) {
    .iwmde_ordinate_value_matches(entry, value) &&
      .iwmde_provenance_matches(entry, provenance) &&
      .iwmde_attribute_method_matches_request(entry, provenance) &&
      .iwmde_posterior_ordinate_supports_bf(entry)
  }, logical(1)))
}


.iwmde_posterior_ordinate_drop_value <- function(posterior_ordinate, value) {

  entries <- .iwmde_posterior_ordinate_entries(posterior_ordinate)
  if (length(entries) == 0L) {
    return(NULL)
  }

  keep <- !vapply(entries, .iwmde_ordinate_value_matches, logical(1),
                  value = value)
  entries <- entries[keep]
  if (length(entries) == 0L) {
    return(NULL)
  }
  if (length(entries) == 1L) {
    return(entries[[1L]])
  }

  return(structure(
    list(status = "ok", ordinates = entries),
    class = c("BayesTools_posterior_ordinates", "list")
  ))
}


.iwmde_posterior_ordinate_keep_values <- function(posterior_ordinate, values) {

  entries <- .iwmde_posterior_ordinate_entries(posterior_ordinate)
  if (length(entries) == 0L || length(values) == 0L) {
    return(NULL)
  }

  keep <- vapply(entries, function(entry) {

    any(vapply(values, function(value) {

      .iwmde_ordinate_value_matches(entry, value)
    }, logical(1)))
  }, logical(1))
  entries <- entries[keep]
  if (length(entries) == 0L) {
    return(NULL)
  }
  if (length(entries) == 1L) {
    return(entries[[1L]])
  }

  return(structure(
    list(status = "ok", ordinates = entries),
    class = c("BayesTools_posterior_ordinates", "list")
  ))
}
