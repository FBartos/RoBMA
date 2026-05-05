# ============================================================================ #
# diagnostics-helpers.R
# ============================================================================ #
#
# Shared diagnostic helpers used by residual, influence, and deletion methods.
#
# ============================================================================ #


.diagnostic_zero_variance_note <- function(diagnostic, parameters,
                                           variance = "LOO posterior") {

  parameters <- paste(parameters, collapse = ", ")
  return(paste0(
    diagnostic,
    " could not be computed for parameter(s) ",
    parameters,
    " because the ",
    variance,
    " variance is zero; values are reported as NaN."
  ))
}


.diagnostic_with_note <- function(x, class, note) {

  attr(x, "note") <- note
  class(x)        <- c(class, class(x))

  return(x)
}


.print_diagnostic_note <- function(note) {

  if (!is.null(note) && nzchar(note)) {
    cat("\nNote: ", note, "\n", sep = "")
  }

  return(invisible(NULL))
}


.diagnostic_psis_context <- function(model, context = NULL) {

  if (!is.null(context)) {
    return(context)
  }

  loo_result   <- loo.brma(model, unit = "estimate")
  psis_object  <- loo_result[["psis_object"]]
  psis_weights <- loo::weights.importance_sampling(
    psis_object,
    log       = FALSE,
    normalize = TRUE
  )

  return(list(
    loo_result   = loo_result,
    psis_object  = psis_object,
    psis_weights = psis_weights
  ))
}


.diagnostic_psis_weights <- function(model, weights = NULL) {

  if (!is.null(weights)) {
    return(weights)
  }

  return(loo_weights(model, unit = "estimate"))
}


.diagnostic_check_loo <- function(model, context = NULL, unit = "estimate") {

  if (is.null(context)) {
    return(check_loo(model, unit = unit))
  }

  loo_result <- context[["loo_result"]]
  pareto_k   <- loo_result[["diagnostics"]][["pareto_k"]]
  bad_k      <- which(pareto_k > 0.7)

  if (length(bad_k) > 0L) {
    warning(
      "Some Pareto k values are high (> 0.7), indicating potentially unreliable ",
      "LOO diagnostics for ", unit, "s: ", paste(bad_k, collapse = ", "), ". ",
      "Inspect the loo fit by using loo(object).",
      call. = FALSE
    )
  }

  return(invisible(NULL))
}


.diagnostic_study_labels <- function(model) {

  labels <- as.character(.get_estimate_labels(model))
  empty  <- is.na(labels) | !nzchar(labels)

  if (any(empty)) {
    labels[empty] <- as.character(which(empty))
  }
  if (anyDuplicated(labels)) {
    labels <- make.unique(labels)
  }

  return(labels)
}


.diagnostic_set_rownames <- function(x, model) {

  labels <- .diagnostic_study_labels(model)
  if (NROW(x) == length(labels)) {
    rownames(x) <- labels
  }

  return(x)
}


.diagnostic_set_names <- function(x, model) {

  labels <- .diagnostic_study_labels(model)
  if (length(x) == length(labels)) {
    names(x) <- labels
  }

  return(x)
}


.diagnostic_collect_notes <- function(...) {

  notes <- unlist(list(...), use.names = FALSE)
  notes <- notes[!is.na(notes) & nzchar(notes)]
  notes <- unique(notes)

  if (length(notes) == 0L) {
    return(NULL)
  }

  return(paste(notes, collapse = " "))
}
