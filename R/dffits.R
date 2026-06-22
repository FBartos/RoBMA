# ============================================================================ #
# brma.dffits.R
# ============================================================================ #
#
# This file contains the dffits method for brma model fits. DFFITS measures
# the number of standard deviations that the fitted value changes if observation
# i is removed.
#
# ============================================================================ #


#' @export
dffits <- function(model, ...) UseMethod("dffits")

#' @title DFFITS for brma Objects
#'
#' @description Computes DFFITS (Difference in FITS, standardized) for a
#' fitted brma object. DFFITS measures how much the fitted value for observation
#' \eqn{i} changes if observation \eqn{i} is removed, standardized by the
#' estimated standard error of the fit.
#'
#' @param model a fitted normal-outcome \code{brma} object without a
#' weightfunction component.
#' @param ... additional arguments (currently ignored).
#'
#' @details
#' DFFITS values are computed as a PSIS leave-one-out deletion diagnostic. For
#' each observation \eqn{i}, the leave-one-out posterior mean fitted value at
#' that observation is estimated with normalized PSIS weights and compared to
#' the full-posterior fitted value:
#' \deqn{DFFITS_i =
#'   \frac{\hat{\mu}_i - \hat{\mu}_{i(-i)}}{SD_{(-i)}(\mu_i)}}
#'
#' This targets deletion influence on fitted values directly. It does not use
#' LOO-PIT residuals, which are predictive outlier diagnostics rather than
#' fitted-value deletion diagnostics.
#' For \code{brma.mv()} known-\code{V} models, DFFITS uses estimate-unit PSIS
#' weights. With correlated known-\code{V}, deletion is conditional estimate
#' deletion and the reported fitted-value target is the fixed-location mean
#' \eqn{\mu = X\beta}; sampled or marginalized random effects are not included
#' in the reported fitted value.
#'
#' Estimate-unit LOO must first be computed with
#' \code{model <- add_loo(model, unit = "estimate")}. If the leave-one-out
#' posterior SD of a fitted value is near zero, the corresponding DFFITS value
#' is returned as \code{NA}.
#'
#' @return A named numeric vector of DFFITS values, one for each observation.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit <- add_loo(fit)
#'
#'   dffits(fit)
#' }
#' }
#'
#' @seealso \code{\link{influence.brma}}, \code{\link{cooks.distance.brma}}, \code{\link{hatvalues.brma}}
#' @aliases dffits
#' @exportS3Method
dffits.brma <- function(model, ...) {

    # the function relies on the normal-normal hat matrix
    outcome_type      <- .outcome_type(model)
    is_weightfunction <- .is_weightfunction(model)

    if (outcome_type != "norm") {
      stop("dffits is only available for normal outcome models.", call. = FALSE)
    }
    if (is_weightfunction) {
      stop("dffits is not available for selection models (weightfunction).", call. = FALSE)
    }

    fit_samples <- .influence_fit_samples(model)
    weights     <- .diagnostic_psis_weights(model)
    dffits_vec  <- .dffits_internal(fit_samples, weights)
    dffits_vec  <- .diagnostic_set_names(dffits_vec, model)
    if (inherits(model, "brma.mv")) {
      dffits_vec <- .brma_mv_attach_target_metadata(dffits_vec, "dffits()")
    }

    return(dffits_vec)
}

.dffits_internal <- function(fit_samples, weights) {

  summary <- .psis_fit_influence_summary(fit_samples, weights)

  delta <- summary[["full_fit"]] - diag(summary[["loo_fit"]])
  se    <- sqrt(diag(summary[["loo_var"]]))
  out   <- delta / se
  out[se <= sqrt(.Machine$double.eps)] <- NA_real_

  names(out) <- colnames(fit_samples)

  return(out)
}


# ---------------------------------------------------------------------------- #
# .influence_fit_samples
# ---------------------------------------------------------------------------- #
#
# Posterior fitted-value samples for fixed location terms used by PSIS deletion
# diagnostics.
#
# ---------------------------------------------------------------------------- #
.influence_fit_samples <- function(model) {

  data              <- model[["data"]]
  priors            <- model[["priors"]]
  posterior_samples <- .get_posterior_samples(model[["fit"]])
  K                 <- nrow(data[["outcome"]])

  fit_samples <- .evaluate.brma.mu(
    fit               = model[["fit"]],
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (.is_mods(model)) .create_fit_formula_list(data = data, "mods") else NULL,
    mods_priors       = if (.is_random(model)) priors[["location"]] else priors[["mods"]],
    is_mods           = .is_mods(model),
    is_PET            = .is_PET(model),
    is_PEESE          = .is_PEESE(model),
    effect_direction  = .effect_direction(model),
    bias_adjusted     = FALSE,
    K                 = K,
    posterior_samples = posterior_samples
  )

  colnames(fit_samples) <- .get_estimate_labels(model)

  return(fit_samples)
}


# ---------------------------------------------------------------------------- #
# .psis_fit_influence_summary
# ---------------------------------------------------------------------------- #
#
# Full and PSIS leave-one-out fitted-value moments.
#
# ---------------------------------------------------------------------------- #
.psis_fit_influence_summary <- function(fit_samples, weights) {

  if (!is.matrix(fit_samples)) {
    fit_samples <- as.matrix(fit_samples)
  }
  if (!is.matrix(weights)) {
    weights <- as.matrix(weights)
  }
  if (nrow(fit_samples) != nrow(weights) ||
      ncol(fit_samples) != ncol(weights)) {
    stop("'fit_samples' and 'weights' must have the same dimensions.",
         call. = FALSE)
  }

  full_fit <- colMeans(fit_samples)
  loo_fit  <- crossprod(weights, fit_samples)
  loo_m2   <- crossprod(weights, fit_samples^2)
  loo_var  <- pmax(loo_m2 - loo_fit^2, 0)

  return(list(
    full_fit = full_fit,
    loo_fit  = loo_fit,
    loo_var  = loo_var
  ))
}
