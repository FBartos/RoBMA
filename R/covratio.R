# ============================================================================ #
# brma.covratio.R
# ============================================================================ #
#
# This file contains the covratio method for brma model fits. COVRATIO
# measures the impact of an observation on the precision of the estimates.
#
# ============================================================================ #


#' @export
covratio <- function(model, ...) UseMethod("covratio")

#' @title COVRATIO for brma Objects
#'
#' @description Computes COVRATIO for a fitted brma object. COVRATIO measures
#' the change in the determinant of the covariance matrix of the estimates
#' when observation \eqn{i} is removed.
#'
#' @param model a fitted brma object.
#' @param type type of parameters to be summarized. Defaults to \code{"mods"}
#' (for the effect size and meta-regression coefficients). Use \code{"scale"}
#' for heterogeneity and scale-regression coefficients.
#' @param component optional parameter namespace. Use \code{"random"} with one
#'   explicitly selected semantic random-effect quantity.
#' @param parameter semantic random-effect quantity used when
#'   \code{component = "random"}.
#' @param ... additional arguments. The internal \code{.weights} argument can
#' supply precomputed PSIS weights for callers that already extracted them.
#'
#' @details
#' COVRATIO is computed using importance sampling weights to approximate the
#' leave-one-out covariance matrices without refitting the model.
#' Estimate-unit LOO must first be computed with
#' \code{model <- add_loo(model, unit = "estimate")}, unless internal weights
#' are supplied.
#'
#' \deqn{COVRATIO_i = \frac{\det(Cov(\beta)_{-i})}{\det(Cov(\beta))}}
#'
#' Values > 1 indicate that the observation improves precision (decreases
#' variance), while values < 1 indicate that the observation decreases precision
#' (increases variance).
#' For \code{brma.mv()} objects, COVRATIO uses estimate-unit PSIS weights. With
#' correlated known-\code{V} data, this is influence under conditional estimate
#' deletion, not independent-study deletion.
#' Parameters with zero posterior variance are excluded from the covariance
#' determinant and reported in a printed note when available. If no estimable
#' parameters remain, or if a full or LOO covariance determinant is zero or
#' non-finite after exclusion, COVRATIO is undefined and values are reported as
#' \code{NaN}.
#'
#' @return A named numeric vector of COVRATIO values, one for each observation.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit <- add_loo(fit)
#'
#'   covratio(fit)
#' }
#' }
#'
#' @seealso \code{\link{influence.brma}}, \code{\link{dffits.brma}}, \code{\link{cooks.distance.brma}}
#' @aliases covratio
#' @exportS3Method
covratio.brma <- function(model, type = "mods", component = NULL,
                          parameter = NULL, ...) {

  dots <- list(...)
  .weights <- dots[[".weights"]]
  BayesTools::check_char(type, "type", allow_values = c("mods", "scale"))
  BayesTools::check_char(component, "component", check_length = 1,
                         allow_NULL = TRUE)
  BayesTools::check_char(parameter, "parameter", check_length = 1,
                         allow_NULL = TRUE)
  component <- .diagnostic_parameter_component(
    type          = type,
    component     = component,
    type_supplied = !missing(type),
    allow_bias    = FALSE
  )
  if (identical(component, "random") && is.null(parameter)) {
    stop(
      "COVRATIO for random-effect quantities requires one explicit 'parameter'.",
      call. = FALSE
    )
  }
  if (!identical(component, "random") && !is.null(parameter)) {
    stop("'parameter' is currently available only with component = 'random'.",
         call. = FALSE)
  }

  if (is.null(.weights)) {
    psis_context <- .diagnostic_psis_context(model)
    .diagnostic_check_loo(model, context = psis_context, unit = "estimate")
    weights <- psis_context[["psis_weights"]]
  } else {
    weights <- .diagnostic_psis_weights(model, .weights)
  }

  # determine whether to extract formula (for meta-regression) or parameter (for intercept-only)
  is_scale <- .is_scale(model)

  # We follow dfbetas logic here:
  if (component == "mods") {
    samples_table <- .diagnostic_location_parameter_samples(model)
  } else if (component == "scale") {
    if (is_scale) {
      # scale-regression: tau is modeled via log_tau formula
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                    = model[["fit"]],
        keep_formulas          = "log_tau",
        random_effects_summary = "none",
        remove_diagnostics     = TRUE,
        return_samples         = TRUE
      )
    } else {
      # random/fixed effects: tau is a parameter (intercept)
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                    = model[["fit"]],
        keep_parameters        = "tau",
        random_effects_summary = "none",
        remove_diagnostics     = TRUE,
        return_samples         = TRUE
      )
    }
  } else if (component == "random") {
    samples_table <- .diagnostic_random_parameter_samples(
      model     = model,
      parameter = parameter
    )
  }

  # Ensure matrix (S x P)
  beta_samples <- as.matrix(samples_table)

  result <- .covratio_internal(beta_samples, weights)
  K      <- length(result[["values"]])
  note   <- NULL

  if (any(result[["excluded"]])) {
    zero_note <- .diagnostic_excluded_zero_variance_note(
      diagnostic = "COVRATIO",
      parameters = colnames(beta_samples)[result[["excluded"]]],
      variance   = "posterior"
    )
    if (all(result[["excluded"]])) {
      note <- .diagnostic_collect_notes(
        zero_note,
        "COVRATIO could not be computed because no parameters with non-zero posterior variance remain; values are reported as NaN."
      )
      return(.diagnostic_with_note(
        .diagnostic_set_names(rep(NaN, K), model),
        class = "covratio.brma",
        note  = note
      ))
    }
    note <- .diagnostic_collect_notes(note, zero_note)
  }

  if (!result[["full_defined"]]) {
    note <- .diagnostic_collect_notes(
      note,
      "COVRATIO could not be computed because the full posterior covariance determinant is zero or non-finite; values are reported as NaN."
    )
    return(.diagnostic_with_note(
      .diagnostic_set_names(rep(NaN, K), model),
      class = "covratio.brma",
      note  = note
    ))
  }

  out <- result[["values"]]
  out <- .diagnostic_set_names(out, model)
  if (any(!result[["loo_defined"]])) {
    note <- .diagnostic_collect_notes(
      note,
      "COVRATIO could not be computed for one or more observations because the LOO covariance determinant is zero or non-finite; affected values are reported as NaN."
    )
  }
  if (!is.null(note)) {
    out <- .diagnostic_with_note(
      out,
      class = "covratio.brma",
      note  = note
    )
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .covratio_internal
# ---------------------------------------------------------------------------- #
#
# Compute covariance determinant ratios in affine-standardized coordinates.
#
# ---------------------------------------------------------------------------- #
.covratio_internal <- function(samples, weights) {

  coordinates <- .influence_sample_coordinates(samples)
  weights     <- .influence_normalize_weights(weights, nrow(samples))
  excluded    <- !coordinates[["variable"]]
  K           <- ncol(weights)
  values      <- rep(NaN, K)
  loo_defined <- rep(FALSE, K)

  if (all(excluded)) {
    return(list(
      values       = values,
      excluded     = excluded,
      full_defined = FALSE,
      loo_defined  = loo_defined
    ))
  }

  samples  <- coordinates[["samples"]][, !excluded, drop = FALSE]
  full_cov <- stats::cov.wt(
    samples,
    wt     = rep(1 / nrow(samples), nrow(samples)),
    method = "ML"
  )[["cov"]]
  full_log_det <- .positive_log_determinant(full_cov)
  if (!is.finite(full_log_det)) {
    return(list(
      values       = values,
      excluded     = excluded,
      full_defined = FALSE,
      loo_defined  = loo_defined
    ))
  }

  for (i in seq_len(K)) {
    loo_cov     <- stats::cov.wt(
      samples,
      wt     = weights[, i],
      method = "ML"
    )[["cov"]]
    loo_log_det <- .positive_log_determinant(loo_cov)
    if (is.finite(loo_log_det)) {
      values[[i]]      <- exp(loo_log_det - full_log_det)
      loo_defined[[i]] <- TRUE
    }
  }

  return(list(
    values       = values,
    excluded     = excluded,
    full_defined = TRUE,
    loo_defined  = loo_defined
  ))
}


# ---------------------------------------------------------------------------- #
# .positive_log_determinant
# ---------------------------------------------------------------------------- #
#
# Return a finite log determinant only for positive-definite matrices.
#
# ---------------------------------------------------------------------------- #
.positive_log_determinant <- function(x) {

  value <- determinant(x, logarithm = TRUE)
  if (value[["sign"]] != 1 || !is.finite(value[["modulus"]])) {
    return(NA_real_)
  }

  return(as.numeric(value[["modulus"]]))
}


#' @exportS3Method
print.covratio.brma <- function(x, ...) {

  note <- attr(x, "note")
  attr(x, "note") <- NULL
  class(x) <- NULL
  print(x, ...)
  .print_diagnostic_note(note)

  return(invisible(x))
}


#' @title Convert COVRATIO Results to a Data Frame
#'
#' @description Converts a \code{covratio.brma} vector to a component-aware
#' long data frame.
#'
#' @param x a \code{covratio.brma} object.
#' @param row.names \code{NULL} or a character vector giving the row names.
#' @param optional logical; passed to the final data-frame coercion.
#' @param stringsAsFactors accepted for compatibility with \code{data.frame()}.
#' @param ... unused additional arguments.
#'
#' @return A plain \code{data.frame} with leading \code{component} and
#' \code{parameter} columns and a numeric \code{value} column.
#'
#' @export
as.data.frame.covratio.brma <- function(
    x, row.names = NULL, optional = FALSE, stringsAsFactors = FALSE, ...) {

  table <- data.frame(
    value       = as.numeric(x),
    row.names   = names(x),
    check.names = FALSE
  )
  output <- .output_table_as_long_data_frame(
    table            = table,
    component        = "covratio",
    row.names        = row.names,
    optional         = optional,
    stringsAsFactors = stringsAsFactors
  )

  return(output)
}
