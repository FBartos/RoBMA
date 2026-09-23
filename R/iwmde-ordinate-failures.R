# Failed point computations are diagnostics, never BayesTools ordinate objects.
.iwmde_ordinate_numerical_status <- function(ordinate, log_ordinate = NULL,
                                             computed = NULL) {

  if (identical(computed, FALSE)) return("not_computed")
  value <- .iwmde_public_numeric(ordinate)
  log_value <- .iwmde_public_numeric(log_ordinate)
  if (is.finite(value) && value > 0) return("finite")
  if (is.finite(value) && value < 0) return("invalid_density")
  if (is.finite(log_value) && !is.na(value)) {
    represented <- exp(log_value)
    if (value == 0) return(if (represented == 0) "underflow" else "arithmetic_failure")
    if (value == Inf) return(if (represented == Inf) "overflow" else "arithmetic_failure")
  }
  if (!is.na(value) && value == 0) return("zero_without_log_evidence")
  "nonfinite"
}


.iwmde_ordinate_numerical_reason <- function(status, log_ordinate = NULL,
                                            reason = NULL) {

  log_value <- .iwmde_public_numeric(log_ordinate)
  if (!is.null(reason)) reason <- sub("[[:space:].]+$", "", reason)
  log_text <- if (is.finite(log_value)) {
    paste0(" (computed log ordinate ", format(log_value, digits = 8L), ")")
  } else ""
  switch(status,
    finite = NULL,
    underflow = paste0("posterior density underflowed to zero on the ordinary scale",
      log_text, ". This does not establish a mathematical zero or an unreliable estimate. ",
      "Ordinary-scale uncertainty is unavailable after this loss. ",
      "The point-null Bayes factor is unavailable on this scale; inspect 'density_diagnostics()'"),
    overflow = paste0("posterior density overflowed on the ordinary scale", log_text,
      ". The point-null Bayes factor is unavailable on this scale; inspect 'density_diagnostics()'"),
    arithmetic_failure = paste0("ordinary-scale posterior density arithmetic failed despite a representable log estimate",
      log_text, ". Inspect 'density_diagnostics()'"),
    invalid_density = paste0("posterior density computation returned a negative value. ",
      "Inspect 'density_diagnostics()' for the available numerical diagnostics"),
    zero_without_log_evidence = paste0("posterior density was computed as zero without a finite log estimate; ",
      "underflow and a mathematical zero cannot be distinguished. Inspect 'density_diagnostics()' and the target support"),
    not_computed = paste0(if (!is.null(reason) && nzchar(reason)) reason else
      "the requested posterior ordinate was not computed", ". Inspect 'density_diagnostics()' for the computation stage"),
    paste0("posterior density computation was non-finite without sufficient evidence to identify underflow or overflow. ",
      "Inspect 'density_diagnostics()' for the available numerical diagnostics")
  )
}


.iwmde_estimate_ordinate_failures <- function(plan, diagnostic,
                                             accepted = NULL) {

  values <- plan[["outputs"]][["requested_values"]]
  if (is.null(values)) values <- plan[["grids"]][["requested_values"]]
  values <- .iwmde_sorted_ordinate_values(values)
  rows <- lapply(values, function(value) {
    if (!is.null(.iwmde_posterior_ordinate_keep_values(accepted, value))) return(NULL)
    selected <- diagnostic
    if (identical(diagnostic[["status"]], "ok") && !is.null(diagnostic[["iwmde"]])) {
      selected <- .iwmde_select_ordinate_diagnostic(diagnostic, value)
    }
    raw <- selected[["diagnostics"]]
    computed <- if (!identical(selected[["status"]], "ok")) FALSE else raw[["bf_included"]]
    if (!is.null(raw[["bf_value"]]) && !isTRUE(raw[["bf_value"]] == value)) computed <- FALSE
    if (identical(computed, FALSE)) {
      raw[["bf_ordinate"]] <- NA_real_
      raw[["bf_log_ordinate"]] <- NA_real_
      raw[["bf_evaluation_value"]] <- NA_real_
    }
    normalized <- raw
    fields <- startsWith(names(raw), "bf_")
    names(normalized)[fields] <- substring(names(raw)[fields], 4L)
    normalized[["estimator"]] <- plan[["method"]]
    normalized[["prior_ordinates"]] <- .iwmde_prior_ordinates_select(plan[["prior_ordinates"]], value)
    normalized[["ordinate_warnings"]] <- .iwmde_ordinate_prior_warnings(
      plan[["target"]][["parameter"]], normalized[["prior_ordinates"]]
    )
    selected[["plan"]] <- plan
    metadata <- plan[["target"]][["metadata"]]
    entry <- c(list(value = value, ordinate = raw[["bf_ordinate"]],
      method = plan[["method"]], density_method = plan[["density_method"]],
      computed = computed, diagnostics = normalized,
      iwmde_provenance = .iwmde_result_provenance(selected, plan[["density_method"]],
        metadata = metadata, density_control = plan[["control"]], value = value,
        evaluation_value = raw[["bf_evaluation_value"]], attribute = "ordinate")), metadata)
    if (is.null(entry[["parameter"]])) entry[["parameter"]] <- plan[["target"]][["parameter"]]
    .iwmde_ordinate_failure_row(entry, selected[["reason"]])
  })
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) return(.iwmde_empty_public_density_diagnostics())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("RoBMA_density_diagnostics", "data.frame")
  out
}


.iwmde_ordinate_failure_row <- function(entry, reason = NULL, preserve_accuracy = FALSE) {

  status <- .iwmde_ordinate_numerical_status(entry[["ordinate"]],
    entry[["diagnostics"]][["log_ordinate"]], entry[["computed"]])
  entry[["diagnostics"]][["bf_grade_met"]] <- FALSE
  entry[["diagnostics"]][["target_met"]] <- FALSE
  if (!identical(status, "finite")) {
    # Retain raw fields on the estimate. Ordinary-scale arithmetic artifacts
    # are not evidence of exact MCSE or convergence. A display Jacobian can
    # lose range while source-scale relative diagnostics remain meaningful.
    fields <- c("mcse", "active_branch_mcse", "sampling_mcse")
    if (!preserve_accuracy) {
      fields <- c(fields, "relative_mcse", "active_branch_relative_mcse",
        "sampling_relative_mcse", "ess", "BF_error_percent",
        "ordinate_relative_change", "ordinate_log_change",
        "pilot_ordinate_relative_change", "pilot_ordinate_log_change")
      entry[["diagnostics"]][["precision_target_met"]] <- FALSE
      entry[["diagnostics"]][["sampling_target_met"]] <- FALSE
    }
    for (field in fields) entry[["diagnostics"]][[field]] <- NA_real_
  }
  entry[["failure_reason"]] <- if (identical(status, "finite")) {
    .iwmde_posterior_ordinate_bf_failure_reason(entry)
  } else .iwmde_ordinate_numerical_reason(status, entry[["diagnostics"]][["log_ordinate"]], reason)
  if (is.null(entry[["failure_reason"]])) {
    entry[["failure_reason"]] <- "posterior ordinate did not satisfy the Bayes-factor requirements"
  }
  out <- .iwmde_public_density_diagnostic_row(entry)
  class(out) <- c("RoBMA_density_diagnostics", "data.frame")
  out
}


.iwmde_ordinate_failure_record <- function(estimate, value) {

  records <- estimate[["ordinate_failures"]]
  if (is.null(records) || !nrow(records)) return(NULL)
  index <- which(records[["requested_value"]] == value)
  if (!length(index)) return(NULL)
  records[index[[1L]], , drop = FALSE]
}


.iwmde_stop_ordinate_numerical_failure <- function(ordinate, log_ordinate,
                                                  values = NULL,
                                                  finite_terms = NULL,
                                                  max_weight_share = NULL) {

  stop(structure(list(
    message = paste0("Posterior density aggregation failed on the ordinary scale; ",
      "computed log estimates were retained."),
    call = NULL, stage = "posterior-density aggregation",
    ordinate = ordinate, log_ordinate = log_ordinate, values = values,
    finite_terms = finite_terms, max_weight_share = max_weight_share
  ), class = c("iwmde_ordinate_numerical_error", "error", "condition")))
}


.iwmde_enrich_ordinate_error <- function(error, plan, context, parameter,
                                         density_method, control, values,
                                         metadata = NULL, prior_ordinates = NULL) {

  if (is.null(plan)) {
    plan <- list(target = list(parameter = parameter, metadata = metadata),
      outputs = list(requested_values = values), control = control,
      method = .density_method_iwmde_estimator(density_method), density_method = density_method,
      source_fingerprint = context[["source_fingerprint"]], prior_ordinates = prior_ordinates)
  }
  frames <- lapply(.iwmde_sorted_ordinate_values(values), function(value) {
    current_plan <- plan
    current_plan[["outputs"]][["requested_values"]] <- value
    diagnostic <- list(status = "unsupported", reason = conditionMessage(error))
    if (inherits(error, "iwmde_ordinate_numerical_error")) {
      index <- match(value, error[["values"]])
      if (!is.na(index)) {
        diagnostic <- list(status = "ok", diagnostics = list(
          bf_included = TRUE, bf_value = value, bf_evaluation_value = value,
          bf_ordinate = error[["ordinate"]][[index]],
          bf_log_ordinate = error[["log_ordinate"]][[index]],
          bf_finite_terms = error[["finite_terms"]][[index]],
          bf_max_weight_share = error[["max_weight_share"]][[index]]
        ))
      }
    }
    .iwmde_estimate_ordinate_failures(current_plan, diagnostic)
  })
  error[["density_diagnostics"]] <- if (length(frames)) do.call(rbind, frames) else
    .iwmde_empty_public_density_diagnostics()
  class(error) <- unique(c("RoBMA_density_ordinate_error", class(error)))
  stop(error)
}
