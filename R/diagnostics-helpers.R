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


.diagnostic_excluded_zero_variance_note <- function(diagnostic, parameters,
                                                    variance = "posterior") {

  parameters <- paste(parameters, collapse = ", ")
  return(paste0(
    diagnostic,
    " excluded parameter(s) ",
    parameters,
    " because the ",
    variance,
    " variance is zero."
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


.check_fixed_location_influence_available <- function(model, caller) {

  if (.outcome_type(model) != "norm") {
    stop(caller, " is only available for normal outcome models.",
         call. = FALSE)
  }
  if (.is_weightfunction(model)) {
    stop(caller, " is not available for selection models (weightfunction).",
         call. = FALSE)
  }

  return(invisible(TRUE))
}


.diagnostic_location_parameter_samples <- function(model,
                                                   standardized_coefficients = FALSE,
                                                   transform_factors = TRUE) {

  if (.is_mods(model) || .is_random(model)) {
    # a fixed-zero intercept is not a parameter of the posterior; the table
    # would otherwise append it as a constant (or dependent) column (as summary)
    remove_parameters <- if (.location_omit_fixed_zero_intercept(
      model, transform_scaled = !standardized_coefficients
    )) {
      "mu_intercept"
    } else {
      NULL
    }
    samples <- BayesTools::JAGS_estimates_table(
      fit                    = model[["fit"]],
      keep_formulas          = "mu",
      remove_inclusion       = TRUE,
      remove_parameters      = remove_parameters,
      random_effects_summary = "none",
      remove_diagnostics     = TRUE,
      transform_factors      = transform_factors,
      transform_scaled       = !standardized_coefficients,
      return_samples         = TRUE
    )

    return(.diagnostic_fixed_location_coefficient_samples(as.matrix(samples)))
  }

  samples <- BayesTools::JAGS_estimates_table(
    fit                = model[["fit"]],
    keep_parameters    = "mu",
    remove_inclusion   = TRUE,
    remove_diagnostics = TRUE,
    transform_factors  = transform_factors,
    transform_scaled   = !standardized_coefficients,
    return_samples     = TRUE
  )

  return(as.matrix(samples))
}

.diagnostic_parameter_component <- function(type, component,
                                            type_supplied = FALSE,
                                            allow_bias = FALSE) {

  type_component <- switch(
    type,
    mods  = "mods",
    scale = "scale",
    bias  = "bias"
  )
  if (is.null(component)) {
    return(type_component)
  }

  component <- .parameter_component_normalize(component)
  if (type_supplied && !identical(component, type_component)) {
    stop(
      "Explicit 'type' and 'component' select different parameter namespaces.",
      call. = FALSE
    )
  }
  if (identical(component, "auto")) {
    stop("'component' must select a parameter namespace.", call. = FALSE)
  }
  if (identical(component, "bias") && !allow_bias) {
    stop("component = 'bias' is not available for this diagnostic.",
         call. = FALSE)
  }

  component
}

.diagnostic_random_parameter_samples <- function(
    model, parameter = NULL, standardized_coefficients = FALSE) {

  if (!is.null(parameter)) {
    selected <- .brma_random_parameter_select(
      object                    = model,
      parameter                 = parameter,
      standardized_coefficients = standardized_coefficients
    )
    return(selected[["samples"]])
  }

  bundle <- .brma_random_parameter_bundle(
    object                    = model,
    standardized_coefficients = standardized_coefficients
  )
  if (ncol(bundle[["samples"]]) == 0L) {
    stop("No random-effect quantities are available for this model.",
         call. = FALSE)
  }

  return(bundle[["samples"]])
}


.diagnostic_fixed_location_coefficient_samples <- function(samples_mat) {

  column_names <- colnames(samples_mat)
  keep <- if (is.null(column_names)) logical() else startsWith(column_names, "(mu) ")
  if (!any(keep)) {
    stop("Fixed location coefficient columns could not be identified.",
         call. = FALSE)
  }

  return(samples_mat[, keep, drop = FALSE])
}


# Influence retains the complete posterior request. A semantic correlation
# declared undefined at zero SD is unavailable, rather than a malformed draw.
.diagnostic_check_random_defined <- function(samples, diagnostic) {

  declared <- BayesTools::posterior_metadata(samples, "undefined_draws")
  if (is.null(declared)) return(invisible(NULL))
  samples <- as.matrix(samples)
  correlation <- colnames(samples) %in% names(declared)[declared == "correlation"]
  missing <- is.na(samples)
  allowed <- matrix(correlation, nrow(samples), ncol(samples), byrow = TRUE)
  # Leave undeclared missing or infinite values to the ordinary validator.
  if (any(!is.finite(samples) & !(missing & allowed))) return(invisible(NULL))
  affected <- which(correlation & colSums(missing) > 0L)
  if (length(affected) == 0L) return(invisible(NULL))
  quantities <- vapply(affected, function(j) {
    paste0("'", colnames(samples)[j], "' (", sum(!missing[, j]), " of ",
      nrow(samples), " posterior draws defined)")
  }, character(1))
  stop(
    diagnostic, " is unavailable for random-effect quantity ",
    paste(quantities, collapse = ", "),
    " because the correlation is undefined when either SD is zero. ",
    "Use 'component = \"random\"' and 'parameter' to explicitly select an ",
    "everywhere-defined quantity.",
    call. = FALSE
  )
}


.diagnostic_check_loo <- function(model, context = NULL, unit = "estimate") {

  if (is.null(context)) {
    return(check_loo(model, unit = unit))
  }

  .loo_warn_pareto_k(context[["loo_result"]], unit)

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
