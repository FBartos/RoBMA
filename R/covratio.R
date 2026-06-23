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
covratio.brma <- function(model, type = "mods", ...) {

  dots <- list(...)
  .weights <- dots[[".weights"]]
  BayesTools::check_char(type, "type", allow_values = c("mods", "scale"))

  # Get PSIS weights (S x K matrix)
  weights <- .diagnostic_psis_weights(model, .weights)

  # determine whether to extract formula (for meta-regression) or parameter (for intercept-only)
  is_scale <- .is_scale(model)

  # We follow dfbetas logic here:
  if (type == "mods") {
    samples_table <- .diagnostic_location_parameter_samples(model)
  } else if (type == "scale") {
    if (is_scale) {
      # scale-regression: tau is modeled via log_tau formula
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                = model[["fit"]],
        keep_formulas      = "log_tau",
        remove_diagnostics = TRUE,
        return_samples     = TRUE
      )
    } else {
      # random/fixed effects: tau is a parameter (intercept)
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                = model[["fit"]],
        keep_parameters    = "tau",
        remove_diagnostics = TRUE,
        return_samples     = TRUE
      )
    }
  }

  # Ensure matrix (S x P)
  beta_samples <- as.matrix(samples_table)

  S <- nrow(beta_samples)
  K <- ncol(weights)

  note <- NULL
  undefined <- apply(
    beta_samples,
    2,
    function(x) diff(range(x, na.rm = TRUE)) <= sqrt(.Machine$double.eps)
  )
  if (any(undefined)) {
    zero_note <- .diagnostic_excluded_zero_variance_note(
      diagnostic = "COVRATIO",
      parameters = colnames(beta_samples)[undefined],
      variance   = "posterior"
    )
    if (all(undefined)) {
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
    beta_samples <- beta_samples[, !undefined, drop = FALSE]
  }

  # 1. Full Covariance (Method = "ML" for consistency)
  # We construct uniform weights 1/S to use cov.wt with "ML"
  # This avoids the (S-1)/S mismatch with the LOO calculation
  w_full <- rep(1/S, S)
  cov_full_res <- stats::cov.wt(beta_samples, wt = w_full, method = "ML")

  # Log-Determinant of Full Covariance
  # determinant() returns $modulus which is the log-abs-determinant
  val_full <- as.numeric(determinant(cov_full_res$cov, logarithm = TRUE)$modulus)
  if (!is.finite(val_full)) {
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

  out <- numeric(K)

  # 2. Loop over Studies for LOO Covariance
  for (i in seq_len(K)) {

    # Extract weights for observation i
    # (These should be pre-normalized, but cov.wt handles normalization too)
    w_i <- weights[, i]

    # STABILITY FIX: Use method = "ML"
    # "unbiased" uses factor 1/(1 - sum(w^2)), which explodes if ESS is low.
    # "ML" uses factor 1/sum(w) = 1, which is stable for distribution variance.
    cov_loo_res <- stats::cov.wt(beta_samples, wt = w_i, method = "ML")
    val_loo     <- as.numeric(determinant(cov_loo_res$cov, logarithm = TRUE)$modulus)

    # Compute Ratio
    out[i] <- if (is.finite(val_loo)) exp(val_loo - val_full) else NaN
  }

  out <- .diagnostic_set_names(out, model)
  if (any(is.nan(out))) {
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


#' @exportS3Method
print.covratio.brma <- function(x, ...) {

  note <- attr(x, "note")
  attr(x, "note") <- NULL
  class(x) <- NULL
  print(x, ...)
  .print_diagnostic_note(note)

  return(invisible(x))
}
