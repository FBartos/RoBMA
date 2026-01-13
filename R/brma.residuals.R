# ============================================================================ #
# brma.residuals.R
# ============================================================================ #
#
# This file contains the residuals method for brma model fits. Residuals are
# computed as the difference between observed effect sizes and fitted values.
#
# Design principles:
# - Consistent with predict.brma: use same underlying machinery
# - Simple interface: minimal required arguments
# - Supports both summary statistics and posterior samples
#
# ============================================================================ #


#' @title Residuals for brma Objects
#'
#' @description Computes residuals (observed minus fitted values) from a
#' fitted brma object.
#'
#' @param object a fitted brma object
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param as_samples whether posterior samples should be returned instead of
#' a summary table. Defaults to \code{FALSE}.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' Residuals are computed as:
#' \deqn{r_i = y_i - \hat{\mu}_i}
#' where \eqn{y_i} is the observed effect size and \eqn{\hat{\mu}_i} is the
#' fitted value (prediction from the fixed effects).
#'
#' For meta-regression models, fitted values incorporate moderator effects.
#' For models without moderators, all fitted values equal the pooled effect.
#'
#' For GLMM models (binomial or Poisson), observed effect sizes are computed
#' from the raw frequency data using the same formulas as \code{metafor::escalc}
#' with the default zero-cell adjustment (adding 0.5 to all cells when any cell
#' is zero). Specifically:
#' \itemize{
#'   \item For binomial models (\code{measure = "OR"}): log odds ratio is computed as
#'     \code{log((ai + 0.5) * (di + 0.5) / ((bi + 0.5) * (ci + 0.5)))} when any cell is zero,
#'     otherwise \code{log(ai * di / (bi * ci))}.
#'   \item For Poisson models (\code{measure = "IRR"}): log incidence rate ratio is computed as
#'     \code{log((x1i / t1i) / (x2i / t2i))} with adjustment when any count is zero.
#' }
#'
#' @return If \code{as_samples = FALSE}, returns a \code{brma.residuals} object
#' containing a summary table with one row per observation showing the residual
#' estimate and credible interval. If \code{as_samples = TRUE}, returns a matrix
#' of posterior samples with one column per observation.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get residuals
#' residuals(fit)
#'
#' # get posterior samples
#' samples <- residuals(fit, as_samples = TRUE)
#' }
#'
#' @seealso [predict.brma()], [blup.brma()], [pooled_effect()]
#' @exportS3Method
residuals.brma <- function(object, probs = c(.025, .975), as_samples = FALSE, ...) {

  # input validation
  BayesTools::check_bool(as_samples, "as_samples")

  # get model type

  outcome_type <- .outcome_type(object)

  # get observed effect sizes based on outcome type
  yi <- .get_observed_yi(object, outcome_type)
  K  <- length(yi)

  # get fitted values (mu predictions for each observation)
  # type = "terms" gives fixed effects predictions without random effects
  fitted_samples <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = "terms",
    bias_adjusted = FALSE,
    as_samples    = TRUE,
    quiet         = TRUE
  )

  # compute residuals: yi - fitted
  # fitted_samples is S x K matrix, yi is vector of length K
  # replicate yi across samples for vectorized subtraction
  S      <- nrow(fitted_samples)
  yi_mat <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)
  resid_samples <- yi_mat - fitted_samples

  # name columns
  colnames(resid_samples) <- paste0("resid[", seq_len(K), "]")

  # return samples if requested
  if (as_samples) {
    return(resid_samples)
  }

  # create summary table
  resid_table <- BayesTools::ensemble_estimates_table(
    samples    = asplit(resid_samples, 2),
    parameters = colnames(resid_samples),
    probs      = probs,
    title      = "Residuals:"
  )

  out <- list(
    summary = resid_table,
    data    = object[["data"]]
  )
  class(out) <- "brma.residuals"

  return(out)
}


# ---------------------------------------------------------------------------- #
# .get_observed_yi
# ---------------------------------------------------------------------------- #
#
# Extract or compute observed effect sizes from brma object.
#
# For normal models, yi is directly available in the data.
# For GLMM models, yi is computed from raw frequency data using formulas
# equivalent to metafor::escalc with default zero-cell adjustment.
#
# @param object       brma object
# @param outcome_type character; "norm", "bin", or "pois"
#
# @return numeric vector of observed effect sizes
#
# ---------------------------------------------------------------------------- #
.get_observed_yi <- function(object, outcome_type) {

  outcome_data <- object[["data"]][["outcome"]]

  if (outcome_type == "norm") {

    # normal models: yi is directly available
    return(outcome_data[["yi"]])

  } else if (outcome_type == "bin") {

    # binomial models: compute log odds ratio from 2x2 table
    # using metafor::escalc formula with default zero-cell adjustment (add = 0.5)
    return(.compute_logOR(
      ai  = outcome_data[["ai"]],
      bi  = outcome_data[["n1i"]] - outcome_data[["ai"]],  # bi = n1i - ai
      ci  = outcome_data[["ci"]],
      di  = outcome_data[["n2i"]] - outcome_data[["ci"]],  # di = n2i - ci
      add = 0.5
    ))

  } else if (outcome_type == "pois") {

    # Poisson models: compute log incidence rate ratio
    return(.compute_logIRR(
      x1i = outcome_data[["x1i"]],
      t1i = outcome_data[["t1i"]],
      x2i = outcome_data[["x2i"]],
      t2i = outcome_data[["t2i"]],
      add = 0.5
    ))

  } else {
    stop("Unknown outcome type: ", outcome_type, call. = FALSE)
  }
}


# ---------------------------------------------------------------------------- #
# .compute_logOR
# ---------------------------------------------------------------------------- #
#
# Compute log odds ratio from 2x2 table with zero-cell adjustment.
#
# Uses the same formula as metafor::escalc(measure = "OR") with default
# add = 0.5 and to = "only0" (add only when any cell is zero).
#
# Formula: logOR = log((ai * di) / (bi * ci))
# With adjustment: logOR = log(((ai + add) * (di + add)) / ((bi + add) * (ci + add)))
#
# @param ai  integer vector; events in treatment group
# @param bi  integer vector; non-events in treatment group
# @param ci  integer vector; events in control group
# @param di  integer vector; non-events in control group
# @param add numeric; constant to add to cells when any cell is zero (default 0.5)
#
# @return numeric vector of log odds ratios
#
# ---------------------------------------------------------------------------- #
.compute_logOR <- function(ai, bi, ci, di, add = 0.5) {

  K <- length(ai)
  logOR <- numeric(K)

  for (k in seq_len(K)) {
    # check if any cell is zero
    if (ai[k] == 0 || bi[k] == 0 || ci[k] == 0 || di[k] == 0) {
      # apply zero-cell adjustment
      logOR[k] <- log(((ai[k] + add) * (di[k] + add)) / ((bi[k] + add) * (ci[k] + add)))
    } else {
      # no adjustment needed
      logOR[k] <- log((ai[k] * di[k]) / (bi[k] * ci[k]))
    }
  }

  return(logOR)
}


# ---------------------------------------------------------------------------- #
# .compute_logIRR
# ---------------------------------------------------------------------------- #
#
# Compute log incidence rate ratio from count and person-time data.
#
# Uses the same formula as metafor::escalc(measure = "IRR") with adjustment
# when any count is zero.
#
# Formula: logIRR = log((x1i / t1i) / (x2i / t2i)) = log(x1i * t2i) - log(x2i * t1i)
# With adjustment when x1i = 0 or x2i = 0: add constant to counts
#
# @param x1i integer vector; events in treatment group
# @param t1i numeric vector; person-time in treatment group
# @param x2i integer vector; events in control group
# @param t2i numeric vector; person-time in control group
# @param add numeric; constant to add to counts when any count is zero (default 0.5)
#
# @return numeric vector of log incidence rate ratios
#
# ---------------------------------------------------------------------------- #
.compute_logIRR <- function(x1i, t1i, x2i, t2i, add = 0.5) {

  K <- length(x1i)
  logIRR <- numeric(K)

  for (k in seq_len(K)) {
    # check if any count is zero
    if (x1i[k] == 0 || x2i[k] == 0) {
      # apply adjustment to counts
      logIRR[k] <- log(((x1i[k] + add) / t1i[k]) / ((x2i[k] + add) / t2i[k]))
    } else {
      # no adjustment needed
      logIRR[k] <- log((x1i[k] / t1i[k]) / (x2i[k] / t2i[k]))
    }
  }

  return(logIRR)
}


#' @exportS3Method
summary.brma.residuals <- function(object, ...) {
  print(object, ...)
}


#' @exportS3Method
print.brma.residuals <- function(x, ...) {

  cat("\n")
  print(x[["summary"]])
  cat("\n")

  return(invisible(x))
}
