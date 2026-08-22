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
#' posterior SD of a fitted value is zero, the corresponding DFFITS value
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

  .check_fixed_location_influence_available(model, "dffits")

  psis_context <- .diagnostic_psis_context(model)
  .diagnostic_check_loo(model, context = psis_context, unit = "estimate")

  fit_samples <- .influence_fit_samples(model)
  weights     <- psis_context[["psis_weights"]]
  dffits_vec  <- .dffits_internal(fit_samples, weights)
  dffits_vec  <- .diagnostic_set_names(dffits_vec, model)
  if (inherits(model, "brma.mv")) {
    dffits_vec <- .brma_mv_attach_target_metadata(dffits_vec, "dffits()")
  }

  return(dffits_vec)
}

.dffits_internal <- function(fit_samples, weights, summary = NULL) {

  if (ncol(as.matrix(fit_samples)) != ncol(as.matrix(weights))) {
    stop("'fit_samples' and 'weights' must have the same number of columns.",
         call. = FALSE)
  }

  if (is.null(summary)) {
    summary <- .psis_influence_summary(
      samples     = fit_samples,
      weights     = weights,
      fit_moments = "matching",
      variance    = "matching"
    )
  }

  loo_fit <- summary[["loo_fit"]]
  if (is.matrix(loo_fit)) {
    loo_fit <- diag(loo_fit)
  }
  loo_var <- summary[["loo_var"]]
  if (is.matrix(loo_var)) {
    loo_var <- diag(loo_var)
  }
  delta <- summary[["full_fit"]] - loo_fit
  se    <- sqrt(loo_var)
  out   <- rep(NA_real_, length(se))
  valid <- se > 0
  out[valid] <- delta[valid] / se[valid]

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
    posterior_samples = posterior_samples,
    priors            = priors
  )

  colnames(fit_samples) <- .get_estimate_labels(model)

  return(fit_samples)
}


# ---------------------------------------------------------------------------- #
# .influence_sample_coordinates
# ---------------------------------------------------------------------------- #
#
# Express posterior samples in dimensionless, affine-standardized coordinates.
# This keeps standardized influence diagnostics invariant to parameter units
# and avoids underflow for small but non-zero posterior variation.
#
# ---------------------------------------------------------------------------- #
.influence_sample_coordinates <- function(samples) {

  if (!is.matrix(samples)) {
    samples <- as.matrix(samples)
  }
  if (nrow(samples) == 0L || ncol(samples) == 0L || any(!is.finite(samples))) {
    stop("'samples' must be a non-empty matrix of finite values.",
         call. = FALSE)
  }

  origin  <- samples[1L, ]
  shifted <- sweep(samples, 2, origin, "-")
  if (any(!is.finite(shifted))) {
    stop("The posterior sample range exceeds finite arithmetic.",
         call. = FALSE)
  }

  scale    <- apply(abs(shifted), 2, max)
  variable <- scale > 0
  shifted[, variable] <- sweep(
    shifted[, variable, drop = FALSE],
    2,
    scale[variable],
    "/"
  )
  shifted[, !variable] <- 0

  return(list(
    samples  = shifted,
    origin   = origin,
    scale    = scale,
    variable = variable
  ))
}


# ---------------------------------------------------------------------------- #
# .influence_normalize_weights
# ---------------------------------------------------------------------------- #
#
# Validate and normalize observation-specific importance weights.
#
# ---------------------------------------------------------------------------- #
.influence_normalize_weights <- function(weights, n_samples) {

  if (!is.matrix(weights)) {
    weights <- as.matrix(weights)
  }
  if (nrow(weights) != n_samples || ncol(weights) == 0L ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must be a non-negative finite matrix with one row per sample.",
         call. = FALSE)
  }

  weight_sums <- colSums(weights)
  if (any(!is.finite(weight_sums)) || any(weight_sums <= 0)) {
    stop("Each column of 'weights' must have a positive finite sum.",
         call. = FALSE)
  }

  return(sweep(weights, 2, weight_sums, "/"))
}


# ---------------------------------------------------------------------------- #
# .psis_influence_summary
# ---------------------------------------------------------------------------- #
#
# Full and PSIS leave-one-out moments in affine-standardized coordinates.
#
# ---------------------------------------------------------------------------- #
.psis_influence_summary <- function(
    samples, weights, fit_moments = c("all", "matching"),
    variance = c("all", "matching", "none")) {

  fit_moments <- match.arg(fit_moments)
  variance    <- match.arg(variance)
  coordinates <- .influence_sample_coordinates(samples)
  samples     <- coordinates[["samples"]]
  weights     <- .influence_normalize_weights(weights, nrow(samples))
  matching <- ncol(samples) == ncol(weights)
  if ((fit_moments == "matching" || variance == "matching") && !matching) {
    stop(
      "Matching PSIS influence moments require the same number of sample and weight columns.",
      call. = FALSE
    )
  }
  if (fit_moments == "matching" && variance == "all") {
    stop(
      "All PSIS influence variances require all leave-one-out fit moments.",
      call. = FALSE
    )
  }

  full_fit <- colMeans(samples)
  loo_fit <- if (fit_moments == "all") {
    crossprod(weights, samples)
  } else {
    colSums(weights * samples)
  }
  loo_var <- NULL
  if (variance == "matching") {
    matching_fit <- if (is.matrix(loo_fit)) diag(loo_fit) else loo_fit
    centered     <- sweep(samples, 2L, matching_fit, "-")
    loo_var      <- colSums(weights * centered^2)
  } else if (variance == "all") {
    loo_var <- matrix(
      0,
      nrow     = ncol(weights),
      ncol     = ncol(samples),
      dimnames = list(colnames(weights), colnames(samples))
    )
    for (j in seq_len(ncol(samples))) {
      centered     <- outer(samples[, j], loo_fit[, j], "-")
      loo_var[, j] <- colSums(weights * centered^2)
    }
  }

  return(list(
    full_fit    = full_fit,
    loo_fit     = loo_fit,
    loo_var     = loo_var,
    samples     = samples,
    origin      = coordinates[["origin"]],
    scale       = coordinates[["scale"]],
    variable    = coordinates[["variable"]],
    fit_moments = fit_moments,
    variance    = variance
  ))
}
