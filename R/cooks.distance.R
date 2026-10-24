# ============================================================================ #
# brma.cooks.distance.R
# ============================================================================ #
#
# This file contains the cooks.distance method for brma model fits. Cook's
# distance measures the aggregate impact of an observation on all model
# coefficients.
#
# ============================================================================ #


#' @title Cook's Distance for brma Objects
#'
#' @description Computes Cook's distance for a fitted brma object. Cook's
#' distance measures the aggregate influence of each observation on the
#' model coefficients.
#'
#' @param model a fitted normal-outcome \code{brma} object without a
#' weightfunction component.
#' @param ... additional arguments (currently ignored).
#'
#' @details
#' Cook's distance is computed as a PSIS leave-one-out deletion diagnostic. For
#' each observation \eqn{i}, normalized PSIS weights estimate the fitted values
#' under the leave-one-out posterior. The distance is the posterior Mahalanobis
#' distance between the full-data and leave-one-out fitted-value vectors,
#' divided by the rank of the fixed-effect model matrix:
#' \deqn{D_i = \frac{\Delta_i' V_\mu^+ \Delta_i}{P}}
#'
#' where \eqn{\Delta_i = \hat{\mu} - \hat{\mu}_{(-i)}}, \eqn{V_\mu^+} is the
#' generalized inverse of the full-posterior fitted-value covariance, and
#' \eqn{P} is the rank of the fixed-effect model matrix (the convention of
#' Cook's distance for linear models; the unscaled squared Mahalanobis distance
#' used for meta-analysis by \insertCite{viechtbauer2010outlier;textual}{RoBMA}
#' is \eqn{P} times this value).
#' For \code{brma.mv()} known-\code{V} models, Cook's distance uses
#' estimate-unit PSIS weights. With correlated known-\code{V}, deletion is
#' conditional estimate deletion and the reported fitted-value target is the
#' fixed-location mean \eqn{\mu = X\beta}; sampled or marginalized random
#' effects are not included in the reported fitted value.
#' Selection \code{brma.mv()} influence is currently unavailable. The PSIS
#' fixed-location moments do not require the inverse marginal covariance used
#' by \code{hatvalues()}.
#'
#' @return A numeric vector of Cook's distance values, one for each observation.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit <- add_loo(fit)
#'
#'   cooks.distance(fit)
#' }
#' }
#'
#' @references
#' \insertCite{viechtbauer2010outlier}{RoBMA}
#'
#' @seealso \code{\link{influence.brma}}, \code{\link{dffits.brma}}, \code{\link{hatvalues.brma}}
#' @importFrom stats cooks.distance
#' @exportS3Method
cooks.distance.brma <- function(model, ...) {

  .check_fixed_location_influence_available(model, "cooks.distance")

  psis_context <- .diagnostic_psis_context(model)
  .diagnostic_check_loo(model, context = psis_context, unit = "estimate")

  fit_samples <- .influence_fit_samples(model)
  weights     <- psis_context[["psis_weights"]]
  P           <- .cooks_distance_rank(model)
  d_vec       <- .cooks.distance_internal(fit_samples, weights, P = P)
  d_vec       <- .diagnostic_set_names(d_vec, model)
  if (inherits(model, "brma.mv")) {
    d_vec <- .brma_mv_attach_target_metadata(d_vec, "cooks.distance()")
  }

  return(d_vec)
}

# The rank of the fixed-effect model matrix that scales Cook's distance.
.cooks_distance_rank <- function(model) {

  return(qr(.get_model_matrix(model))[["rank"]])
}

.cooks.distance_internal <- function(fit_samples, weights, summary = NULL,
                                     P = 1L) {

  if (is.null(summary)) {
    summary <- .psis_influence_summary(
      samples     = fit_samples,
      weights     = weights,
      fit_moments = "all",
      variance    = "none"
    )
  }
  variable <- summary[["variable"]]
  if (!any(variable)) {
    out <- rep(0, ncol(weights))
    names(out) <- colnames(fit_samples)
    return(out)
  }

  delta   <- sweep(summary[["loo_fit"]], 2, summary[["full_fit"]], "-")
  delta   <- -delta[, variable, drop = FALSE]

  vcov_fit <- stats::cov(summary[["samples"]][, variable, drop = FALSE])
  vcov_inv <- .symmetric_ginv(vcov_fit)

  d_vec <- rowSums((delta %*% vcov_inv) * delta) / max(P, 1L)
  names(d_vec) <- colnames(fit_samples)

  return(d_vec)
}


# ---------------------------------------------------------------------------- #
# .symmetric_ginv
# ---------------------------------------------------------------------------- #
#
# Generalized inverse for symmetric positive semi-definite covariance matrices.
#
# ---------------------------------------------------------------------------- #
.symmetric_ginv <- function(x) {

  x          <- (x + t(x)) / 2
  components <- svd(x)

  scale <- max(components[["d"]])
  if (scale == 0) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  }

  tolerance <- max(dim(x)) * .Machine$double.eps * scale
  keep      <- components[["d"]] > tolerance
  if (!any(keep)) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  }

  left   <- components[["u"]][, keep, drop = FALSE]
  right  <- components[["v"]][, keep, drop = FALSE]
  values <- components[["d"]][keep]

  return(right %*% diag(1 / values, nrow = length(values)) %*% t(left))
}
