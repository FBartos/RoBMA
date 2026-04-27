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
#' fitted brma object, with options for standardization.
#'
#' @param object a fitted brma object
#' @param type the type of residuals to compute. Options are:
#' \itemize{
#'   \item \code{"outcome"} (default): Raw residuals (observed - fitted).
#'     Available for all model types.
#'   \item \code{"pearson"}: Pearson (semi-standardized) residuals, dividing
#'     raw residuals by the marginal standard error \eqn{\sqrt{v_i + \tau^2}}.
#'     Only available for normal outcome models without selection (weightfunction)
#'     bias adjustment.
#'   \item \code{"rstandard"}: Internally standardized residuals, dividing
#'     raw residuals by their standard errors computed using the hat matrix.
#'     Only available for normal outcome models without selection (weightfunction)
#'     bias adjustment.
#'   \item \code{"LOO-PIT" (alias: \code{"rstudent"}}: Leave-one-out probability
#'     integral transform (PIT) residuals computed via Pareto smoothed importance
#'     sampling. The LOO-CDF value for each observation is computed and transformed
#'     to standard normal quantiles via \eqn{\Phi^{-1}(u_i)}. Under a correctly
#'     specified model, these residuals should follow a standard normal distribution.
#'     This is the recommended standardized residual for Bayesian models as it properly
#'     accounts for estimation uncertainty and leverage. Available for all model
#'     types. Note: This requires that the loo has been computed previously (see
#'     [add_loo()] function).
#' }
#' @param unit output unit. Only \code{"estimate"} is implemented currently..
#' @param conditioning_depth conditioning depth for non-LOO residuals. \code{"marginal"}
#' uses fixed effects only, \code{"cluster"} conditions on cluster-level random
#' effects, and \code{"estimate"} conditions on the full estimate-level fitted
#' value. LOO-PIT residuals always use the estimate-unit LOO target.
#' @param bias_adjusted whether residuals should be computed from bias-adjusted
#' fitted values. Defaults to \code{FALSE}, which means residuals are computed
#' as the difference between observed values and raw (biased) predictions
#' including PET/PEESE terms. Set to \code{TRUE} to compute residuals from
#' bias-corrected fitted values. Only applies to \code{type = "outcome"}. Note
#' that the bias adjustment residuals are not residuals in the traditional sense.
#' @param ... additional arguments.
#'
#' @details
#' Raw residuals (\code{type = "outcome"}) are computed as:
#' \deqn{e_i = y_i - \hat{\mu}_i}
#' where \eqn{y_i} is the observed effect size and \eqn{\hat{\mu}_i} is the
#' fitted value (prediction from the fixed effects).
#'
#' Pearson residuals (\code{type = "pearson"}) divide raw residuals by the
#' marginal standard error:
#' \deqn{r_i^{Pearson} = \frac{e_i}{\sqrt{v_i + \tau^2}}}
#' where \eqn{v_i} is the sampling variance and \eqn{\tau^2} is the
#' relevant heterogeneity component. Only available for normal outcome models
#' without selection (weightfunction) bias adjustment.
#'
#' Standardized residuals (\code{type = "rstandard"}) use the hat matrix to
#' compute residual standard errors that account for the uncertainty in
#' estimated coefficients:
#' \deqn{z_i = \frac{e_i}{\sqrt{[(I-H)M(I-H)']_{ii}}}}
#' where \eqn{H} is the hat matrix and \eqn{M} is the marginal variance-covariance
#' matrix. For models without moderators, this simplifies to the Pearson formula.
#' Only available for normal outcome models without selection (weightfunction)
#' bias adjustment.
#'
#' LOO-PIT residuals (\code{type = "LOO-PIT"}) are the Bayesian equivalent of
#' studentized deleted residuals \insertCite{vehtari2017practical}{RoBMA}. They
#' are computed via leave-one-out probability integral transformation:
#' \deqn{r_i = \Phi^{-1}(u_i)}
#' where \eqn{u_i = \sum_s w_{is} F(y_i | \theta^{(s)})} is the LOO-weighted CDF
#' value, \eqn{w_{is}} are the normalized PSIS weights, and \eqn{F} is the
#' cumulative distribution function of the estimate-unit predictive
#' distribution used by LOO. Under a correctly specified model, LOO-PIT
#' residuals should follow a standard normal distribution. Unlike traditional
#' standardized residuals, LOO-PIT residuals properly account for estimation
#' uncertainty and leverage without requiring a hat matrix. This is the
#' recommended method for standardized residuals in Bayesian meta-analysis.
#'
#' For meta-regression models, fitted values incorporate moderator effects.
#' For models without moderators, all fitted values equal the pooled effect.
#'
#' For GLMM models (binomial or Poisson), observed effect sizes and their
#' sampling variances are computed from the raw frequency data using the
#' same formulas as \code{metafor::escalc} with the default zero-cell
#' adjustment (adding 0.5 to all cells when any cell is zero).
#'
#' The residuals are computed separately for each posterior sample,
#' naturally propagating uncertainty in model parameters to the residuals.
#'
#' @return A numeric vector of residual means, one per estimate.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get raw residuals (default)
#' residuals(fit)
#'
#' # get Pearson (semi-standardized) residuals (normal outcome only)
#' residuals(fit, type = "pearson")
#'
#' # get internally standardized residuals (normal outcome only)
#' residuals(fit, type = "rstandard")
#'
#' # get LOO-PIT residuals (recommended for standardized residuals)
#' # should be ~N(0,1) if model is correctly specified
#' residuals(fit, type = "LOO-PIT")
#'
#' # get residuals from bias-adjusted predictions
#' residuals(fit, bias_adjusted = TRUE)
#'
#' # check LOO diagnostics before using LOO-PIT residuals
#' loo_result <- loo(fit)
#' plot(loo_result) # check Pareto k values
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [predict.brma()], [blup.brma()], [pooled_effect()], [rstandard.brma()], [rstudent.brma()]
#' @exportS3Method
residuals.brma <- function(object, type = "outcome", unit = "estimate",
                           conditioning_depth = "marginal",
                           bias_adjusted = FALSE, ...) {

  # input validation
  dots                         <- list(...)
  conditioning_depth_specified <- !missing(conditioning_depth)
  .check_legacy_level_arg(dots, "residuals()")

  BayesTools::check_char(type, "type", allow_values = c("outcome", "pearson", "rstandard", "rstudent", "LOO-PIT"))
  unit                         <- .normalize_unit(unit)
  conditioning_depth           <- .normalize_conditioning_depth(conditioning_depth)
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")

  if (is.element(type, c("LOO-PIT", "rstudent")) && conditioning_depth_specified) {
    stop(
      "LOO-PIT residuals use the estimate-unit LOO target; ",
      "do not set 'conditioning_depth'.",
      call. = FALSE
    )
  }

  .check_unit_conditioning_depth(
    object             = object,
    unit               = unit,
    conditioning_depth = conditioning_depth,
    caller             = "residuals()"
  )

  if (unit == "cluster") {
    .check_cluster_unit_deferred("residuals()")
  }

  # extract model characteristics for error checking
  outcome_type      <- .outcome_type(object)
  is_weightfunction <- .is_weightfunction(object)

  # check for unsupported type/model combinations
  .check_residual_type_availability(type, outcome_type, is_weightfunction)

  out <- .residuals_estimate.brma(
    object             = object,
    type               = type,
    conditioning_depth = conditioning_depth,
    bias_adjusted      = bias_adjusted
  )

  # clean names
  if (!is.data.frame(out)) {
    names(out) <- NULL
  }
  return(out)
}


# ---------------------------------------------------------------------------- #
# .residuals_estimate.brma
# ---------------------------------------------------------------------------- #
#
# Estimate-unit residuals, one value per effect-size estimate.
#
# @param object        brma object.
# @param type          character; residual type.
# @param conditioning_depth character; conditioning depth.
# @param bias_adjusted      logical; whether PET/PEESE bias adjustments are removed.
#
# @return numeric vector of residuals or standardized residuals.
#
# ---------------------------------------------------------------------------- #
.residuals_estimate.brma <- function(object, type, conditioning_depth,
                                     bias_adjusted) {

  if (is.element(type, c("LOO-PIT", "rstudent"))) {
    return(.standardized_residuals_loopit(object))
  }

  if (type == "rstandard") {
    return(rstandard.brma(
      model              = object,
      unit               = "estimate",
      conditioning_depth = conditioning_depth
    )[["z"]])
  }

  pred_type <- switch(conditioning_depth,
    "marginal" = "terms",
    "cluster"  = "cluster",
    "estimate" = "estimate"
  )

  fitted_samples <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = pred_type,
    bias_adjusted = bias_adjusted,
    quiet         = TRUE
  )

  yi <- .outcome_data_yi(object)
  K  <- length(yi)
  S  <- nrow(fitted_samples)

  yi_mat        <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)
  resid_samples <- yi_mat - fitted_samples

  if (type == "pearson") {
    se_samples    <- .pearson_residual_se_samples(
      object             = object,
      conditioning_depth = conditioning_depth
    )
    resid_samples <- resid_samples / se_samples
  }

  return(colMeans(resid_samples))
}


# ---------------------------------------------------------------------------- #
# .pearson_residual_se_samples
# ---------------------------------------------------------------------------- #
#
# Compute posterior samples of Pearson residual standard errors.
#
# @param object brma object.
# @param conditioning_depth character; conditioning depth.
#
# @return S x K matrix of standard errors.
#
# ---------------------------------------------------------------------------- #
.pearson_residual_se_samples <- function(object, conditioning_depth) {

  is_multilevel <- .is_multilevel(object)
  is_scale      <- .is_scale(object)
  priors        <- object[["priors"]]
  vi            <- .outcome_data_vi(object)
  K             <- length(vi)

  tau_result <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = object[["data"]][["scale"]],
    scale_formula = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors  = priors[["scale"]],
    is_scale      = is_scale,
    is_multilevel = is_multilevel,
    K             = K
  )

  tau_samples <- switch(conditioning_depth,
    "marginal" = tau_result[["tau_total"]],
    "cluster"  = tau_result[["tau_within"]],
    "estimate" = matrix(0, nrow = nrow(tau_result[["tau_total"]]), ncol = K)
  )

  vi_mat <- matrix(vi, nrow = nrow(tau_samples), ncol = K, byrow = TRUE)

  return(sqrt(vi_mat + tau_samples^2))
}


#' @title Internally Standardized Residuals for brma Objects
#'
#' @description Computes internally standardized residuals from a fitted brma
#' object using the hat matrix. Returns a data frame with raw residuals,
#' standard errors, and standardized residuals (z-values). Available for normal
#' outcome models only.
#'
#' @param model a fitted brma object.
#' @param unit output unit. Only \code{"estimate"} is implemented currently.
#' @param conditioning_depth conditioning depth. Options are:
#' \itemize{
#'   \item \code{"marginal"} (default): Residuals from fixed effects predictions
#'     (observed - Xβ).
#'   \item \code{"cluster"}: Residuals from cluster-level predictions (observed -
#'     (Xβ + gamma)). Only available for multilevel (3-level) models.
#'   \item \code{"estimate"}: Residuals from BLUPs,
#'     i.e., deviations of the observed effect sizes from the best linear
#'     unbiased predictions of the estimate-specific true effects (observed - theta).
#' }
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' This function returns a data frame with three columns matching the output
#' of \code{metafor::rstandard}:
#' \itemize{
#'   \item \code{resid}: Raw residuals (observed - fitted values)
#'   \item \code{se}: Standard errors of the residuals
#'   \item \code{z}: Standardized residuals (resid / se)
#' }
#'
#' Internally standardized residuals divide the observed residuals by their
#' corresponding standard errors computed using the hat matrix. For a correctly
#' specified model, these residuals should approximately follow a standard
#' normal distribution.
#'
#' This function is only available for normal outcome models without selection
#' (weightfunction) bias adjustment. For other model types, use
#' \code{\link{rstudent.brma}} which uses LOO-PIT.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{resid}: Raw residuals
#'   \item \code{se}: Standard errors of the residuals
#'   \item \code{z}: Standardized residuals
#' }
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get marginal internally standardized residuals (default)
#' rstandard(fit)
#'
#' # get estimate-level (BLUP-based) residuals
#' rstandard(fit, conditioning_depth = "estimate")
#' }
#'
#' @seealso [rstudent.brma()], [residuals.brma()], [blup.brma()], [predict.brma()]
#' @exportS3Method
rstandard.brma <- function(model, unit = "estimate",
                           conditioning_depth = "marginal", ...) {

  dots <- list(...)
  .check_legacy_level_arg(dots, "rstandard()")

  unit               <- .normalize_unit(unit)
  conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)
  .check_unit_conditioning_depth(
    object             = model,
    unit               = unit,
    conditioning_depth = conditioning_depth,
    caller             = "rstandard()"
  )

  if (unit == "cluster") {
    .check_cluster_unit_deferred("rstandard()")
  }

  # check model type availability
  outcome_type      <- .outcome_type(model)
  is_weightfunction <- .is_weightfunction(model)

  if (outcome_type != "norm") {
    stop(
      "rstandard is only available for normal outcome models. ",
      "Use rstudent() for GLMM models.",
      call. = FALSE
    )
  }
  if (is_weightfunction) {
    stop(
      "rstandard is not available for selection models (weightfunction). ",
      "Use rstudent() for selection models.",
      call. = FALSE
    )
  }
  # compute residuals and standard errors using the same GLS projection
  hat_res <- .compute_hat_matrix_samples(
    object             = model,
    conditioning_depth = conditioning_depth,
    return_se          = TRUE,
    return_resid       = TRUE,
    summarize          = TRUE
  )

  # construct output data frame matching metafor::rstandard format
  out <- data.frame(
    resid = hat_res[["resid"]],
    se    = hat_res[["se"]],
    z     = hat_res[["z"]]
  )
  rownames(out) <- NULL

  return(out)
}


#' @title Externally Standardized (Studentized) Residuals for brma Objects
#'
#' @description Computes externally standardized residuals (also called
#' studentized residuals or standardized deleted residuals) from a fitted brma
#' object using LOO-PIT (Leave-One-Out Probability Integral Transform).
#' Returns a data frame with raw residuals, standard errors, and standardized
#' residuals (z-values).
#'
#' @param model a fitted brma object.
#' @param unit output unit. Only \code{"estimate"} is available for LOO-PIT
#' residuals.
#' @param conditioning_depth unused for LOO-PIT residuals. LOO-PIT residuals
#' always use the estimate-unit LOO target.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' This function returns a data frame with three columns matching the output
#' of \code{metafor::rstudent}:
#' \itemize{
#'   \item \code{resid}: LOO predictive residuals (observed - fitted values)
#'   \item \code{se}: LOO predictive standard errors when available
#'   \item \code{z}: Externally standardized residuals (LOO-PIT transformed)
#' }
#'
#' LOO-PIT residuals are the Bayesian equivalent of studentized deleted
#' residuals. They are computed via leave-one-out probability integral
#' transformation using Pareto smoothed importance sampling. For each
#' observation, the LOO-weighted CDF is computed and transformed to a
#' standard normal quantile.
#'
#' Under a correctly specified model, LOO-PIT residuals should follow a
#' standard normal distribution. Large absolute values may indicate outliers
#' or model misspecification.
#'
#' The \code{z} column is the primary standardized diagnostic. The \code{resid}
#' and \code{se} columns are raw-scale companions computed from LOO predictive
#' moments using the normalized PSIS weights. For selection models, these moments
#' are computed from the fitted weighted-normal predictive distribution. For
#' GLMMs, they are computed on the approximate effect-size scale used by the
#' LOO-PIT diagnostic.
#'
#' Unlike \code{\link{rstandard.brma}} (which uses the hat matrix), LOO-PIT
#' residuals properly account for estimation uncertainty and leverage without
#' requiring explicit hat matrix computation. This makes \code{rstudent.brma}
#' suitable for all model types including selection models and GLMMs.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{resid}: Raw residuals
#'   \item \code{se}: Standard errors of the residuals
#'   \item \code{z}: Externally standardized residuals (LOO-PIT)
#' }
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#' fit <- add_loo(fit)
#'
#' # get externally standardized residuals
#' rstudent(fit)
#'
#' # check LOO diagnostics
#' loo_obj <- loo(fit)
#' plot(loo_obj) # check Pareto k values
#' }
#'
#' @seealso [rstandard.brma()], [residuals.brma()], [loo.brma()], [blup.brma()]
#' @exportS3Method
rstudent.brma <- function(model, unit = "estimate",
                          conditioning_depth = "marginal", ...) {

  dots                         <- list(...)
  conditioning_depth_specified <- !missing(conditioning_depth)
  .check_legacy_level_arg(dots, "rstudent()")

  unit <- .normalize_unit(unit)

  if (unit == "cluster") {
    stop(
      "Cluster-unit rstudent residuals are not available because multivariate ",
      "LOO-PIT residuals are not unique. Use unit = 'estimate'.",
      call. = FALSE
    )
  }

  if (conditioning_depth_specified) {
    stop(
      "LOO-PIT residuals use the estimate-unit LOO target; ",
      "do not set 'conditioning_depth'.",
      call. = FALSE
    )
  }

  .check_unit_conditioning_depth(
    object             = model,
    unit               = unit,
    conditioning_depth = "estimate",
    caller             = "rstudent()"
  )

  setup <- .estimate_likelihood_setup.brma(model)

  # extract PSIS object once and reuse it for PIT and LOO expectations
  loo_result   <- loo.brma(model, unit = "estimate")
  psis_object  <- loo_result[["psis_object"]]
  psis_weights <- loo::weights.importance_sampling(
    psis_object,
    log       = FALSE,
    normalize = TRUE
  )

  # get LOO-PIT z values using the estimate-unit LOO target
  check_loo(model, unit = "estimate")
  z <- .standardized_residuals_loopit(
    object       = model,
    psis_weights = psis_weights,
    check        = FALSE
  )

  moments <- .loo_predictive_moments_estimate(
    object       = model,
    setup        = setup,
    psis_object  = psis_object,
    psis_weights = psis_weights
  )

  resid <- moments[["resid"]]
  se    <- moments[["se"]]

  # construct output data frame matching metafor::rstudent format
  out <- data.frame(
    resid = resid,
    se    = se,
    z     = z
  )
  rownames(out) <- NULL

  return(out)
}


# ---------------------------------------------------------------------------- #
# .loo_predictive_moments_estimate
# ---------------------------------------------------------------------------- #
#
# Estimate-unit LOO predictive residual mean and standard deviation.
#
# The LOO-PIT z-value is the primary standardized residual. These moments provide
# raw-scale companions for the output table and residual funnel. The aggregation
# uses normalized PSIS weights on analytic predictive first and second moments,
# matching loo::E_loo(..., type = "mean") values without recomputing Pareto k.
#
# @param object       brma object.
# @param setup        output from .estimate_likelihood_setup.brma().
# @param psis_object  PSIS object from loo.
# @param psis_weights optional normalized PSIS weights.
#
# @return list with numeric vectors resid and se.
#
# ---------------------------------------------------------------------------- #
.loo_predictive_moments_estimate <- function(object, setup, psis_object,
                                             psis_weights = NULL) {

  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  K                 <- setup[["K"]]
  mu_samples        <- setup[["mu"]]
  tau_within        <- setup[["tau_within"]]
  outcome_type      <- setup[["outcome_type"]]
  is_weightfunction <- setup[["is_weightfunction"]]
  effect_direction  <- setup[["effect_direction"]]

  sei_mat <- matrix(sei, nrow = nrow(mu_samples), ncol = K, byrow = TRUE)

  if (outcome_type == "norm" && is_weightfunction) {
    posterior_samples <- setup[["posterior_samples"]]
    omega_samples     <- .extract_omega_samples(posterior_samples)

    if (effect_direction == "negative") {
      mu_moments <- -mu_samples
    } else {
      mu_moments <- mu_samples
    }

    moments <- .outcome_moments.wnorm(
      mu_samples          = mu_moments,
      tau_within          = tau_within,
      sei                 = sei,
      omega               = omega_samples,
      crit_yi             = setup[["fit_data"]][["crit_yi"]],
      use_normal          = .extract_use_normal(object, posterior_samples = posterior_samples),
      bias_indicator      = .extract_bias_indicator(object, posterior_samples = posterior_samples),
      crit_yi_mapping     = setup[["fit_data"]][["crit_yi_mapping"]],
      crit_yi_mapping_max = setup[["fit_data"]][["crit_yi_mapping_max"]]
    )

    mean_samples   <- moments[["mean"]]
    second_samples <- moments[["second"]]

    if (effect_direction == "negative") {
      mean_samples <- -mean_samples
    }
  } else {
    pred_var       <- tau_within^2 + sei_mat^2
    mean_samples   <- mu_samples
    second_samples <- pred_var + mu_samples^2
  }

  mean_samples   <- as.matrix(mean_samples)
  second_samples <- as.matrix(second_samples)

  if (is.null(psis_weights)) {
    psis_weights <- loo::weights.importance_sampling(
      psis_object,
      log       = FALSE,
      normalize = TRUE
    )
  }

  pred_mean   <- colSums(psis_weights * mean_samples)
  pred_second <- colSums(psis_weights * second_samples)

  resid <- yi - pred_mean
  se    <- sqrt(pmax(pred_second - pred_mean^2, 0))

  return(list(
    resid = resid,
    se    = se
  ))
}


# ---------------------------------------------------------------------------- #
# .outcome_moments.wnorm
# ---------------------------------------------------------------------------- #
#
# Per-posterior-sample first and second predictive moments for selection models.
#
# @param mu_samples       S x K matrix of location samples in selection space.
# @param tau_within       S x K matrix of estimate-level heterogeneity samples.
# @param sei              numeric vector of length K.
# @param omega            S x W matrix of omega samples.
# @param crit_yi          W - 1 x K matrix of selection cutoffs.
# @param use_normal       optional logical vector for normal fast path.
# @param bias_indicator   optional posterior branch indicator.
# @param crit_yi_mapping  optional branch-specific cutoff mapping.
# @param crit_yi_mapping_max optional active cutoff count per branch.
#
# @return list with S x K matrices mean and second.
#
# ---------------------------------------------------------------------------- #
.outcome_moments.wnorm <- function(mu_samples, tau_within, sei, omega, crit_yi,
                                   use_normal = NULL, bias_indicator = NULL,
                                   crit_yi_mapping = NULL,
                                   crit_yi_mapping_max = NULL) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  mean_samples   <- matrix(NA_real_, nrow = S, ncol = K)
  second_samples <- matrix(NA_real_, nrow = S, ncol = K)

  has_mapping <- .selection_has_mapping(
    bias_indicator      = bias_indicator,
    crit_yi_mapping     = crit_yi_mapping,
    crit_yi_mapping_max = crit_yi_mapping_max
  )

  if (has_mapping) {
    return(.wnorm_mix_moments_matrix(
      mean                = mu_samples,
      sd                  = total_sd,
      omega               = omega,
      crit_yi             = crit_yi,
      bias_indicator      = bias_indicator,
      crit_yi_mapping     = crit_yi_mapping,
      crit_yi_mapping_max = crit_yi_mapping_max
    ))
  }

  for (k in seq_len(K)) {
    moments <- .wnorm_moments_fast.ss.matrix(
      mean       = mu_samples[, k],
      sd         = total_sd[, k],
      omega      = omega,
      crit_x     = crit_yi[, k],
      use_normal = use_normal
    )
    mean_samples[, k]   <- moments[["mean"]]
    second_samples[, k] <- moments[["second"]]
  }

  return(list(
    mean   = mean_samples,
    second = second_samples
  ))
}


# ---------------------------------------------------------------------------- #
# .residuals_cluster.brma
# ---------------------------------------------------------------------------- #
#
# Deferred cluster-unit residual diagnostics.
#
# @param object brma object.
# @param type   character; residual type requested by the public caller.
# @param conditioning_depth character; conditioning depth.
#
# @return stops with a deferred-design error.
#
# ---------------------------------------------------------------------------- #
.residuals_cluster.brma <- function(object, type, conditioning_depth) {

  .check_cluster_unit_deferred("residuals()")
}


# ---------------------------------------------------------------------------- #
# .outcome_data_yi
# ---------------------------------------------------------------------------- #
#
# Extract or compute observed effect sizes from brma object.
#
# For normal models, yi is directly available in the data.
# For GLMM models, yi is computed from raw frequency data using formulas
# equivalent to metafor::escalc with default zero-cell adjustment (add = 0.5).
#
# @param object brma object
#
# @return numeric vector of observed effect sizes
#
# ---------------------------------------------------------------------------- #
.outcome_data_yi <- function(object) {

  outcome_type <- .outcome_type(object)
  outcome_data <- object[["data"]][["outcome"]]

  if (outcome_type == "norm") {
    # normal models: yi is directly available
    return(outcome_data[["yi"]])
  } else if (outcome_type == "bin") {
    # binomial models: compute log odds ratio from 2x2 table
    # using metafor::escalc formula with default zero-cell adjustment (add = 0.5)
    ai <- outcome_data[["ai"]]
    bi <- outcome_data[["n1i"]] - outcome_data[["ai"]] # bi = n1i - ai
    ci <- outcome_data[["ci"]]
    di <- outcome_data[["n2i"]] - outcome_data[["ci"]] # di = n2i - ci
    add <- 0.5

    # vectorized computation with conditional zero-cell adjustment
    needs_adj <- (ai == 0 | bi == 0 | ci == 0 | di == 0)
    ai_adj <- ai + add * needs_adj
    bi_adj <- bi + add * needs_adj
    ci_adj <- ci + add * needs_adj
    di_adj <- di + add * needs_adj

    return(log((ai_adj * di_adj) / (bi_adj * ci_adj)))
  } else if (outcome_type == "pois") {
    # Poisson models: compute log incidence rate ratio
    x1i <- outcome_data[["x1i"]]
    t1i <- outcome_data[["t1i"]]
    x2i <- outcome_data[["x2i"]]
    t2i <- outcome_data[["t2i"]]
    add <- 0.5

    # vectorized computation with conditional zero-cell adjustment
    needs_adj <- (x1i == 0 | x2i == 0)
    x1i_adj <- x1i + add * needs_adj
    x2i_adj <- x2i + add * needs_adj

    return(log((x1i_adj / t1i) / (x2i_adj / t2i)))
  }
}


# ---------------------------------------------------------------------------- #
# .outcome_data_sei
# ---------------------------------------------------------------------------- #
#
# Extract or compute sampling standard errors from brma object.
#
# For normal models, sei is directly available in the data.
# For GLMM models, sei is computed from raw frequency data using formulas
# equivalent to metafor::escalc with default zero-cell adjustment (add = 0.5).
#
# @param object brma object
#
# @return numeric vector of sampling standard errors
#
# ---------------------------------------------------------------------------- #
.outcome_data_sei <- function(object) {
  outcome_type <- .outcome_type(object)
  outcome_data <- object[["data"]][["outcome"]]

  if (outcome_type == "norm") {
    # normal models: sei is directly available
    return(outcome_data[["sei"]])
  } else if (outcome_type == "bin") {
    # binomial models: compute approximate SE for log odds ratio
    # SE(logOR) = sqrt(1/ai + 1/bi + 1/ci + 1/di) with zero-cell adjustment
    ai <- outcome_data[["ai"]]
    bi <- outcome_data[["n1i"]] - outcome_data[["ai"]]
    ci <- outcome_data[["ci"]]
    di <- outcome_data[["n2i"]] - outcome_data[["ci"]]
    add <- 0.5

    # vectorized computation with conditional zero-cell adjustment
    needs_adj <- (ai == 0 | bi == 0 | ci == 0 | di == 0)
    ai_adj <- ai + add * needs_adj
    bi_adj <- bi + add * needs_adj
    ci_adj <- ci + add * needs_adj
    di_adj <- di + add * needs_adj

    return(sqrt(1 / ai_adj + 1 / bi_adj + 1 / ci_adj + 1 / di_adj))
  } else if (outcome_type == "pois") {
    # Poisson models: compute approximate SE for log IRR
    # SE(logIRR) = sqrt(1/x1i + 1/x2i)
    x1i <- outcome_data[["x1i"]]
    x2i <- outcome_data[["x2i"]]
    add <- 0.5

    # vectorized computation with conditional zero-cell adjustment
    needs_adj <- (x1i == 0 | x2i == 0)
    x1i_adj <- x1i + add * needs_adj
    x2i_adj <- x2i + add * needs_adj

    return(sqrt(1 / x1i_adj + 1 / x2i_adj))
  }
}


# ---------------------------------------------------------------------------- #
# .outcome_data_vi
# ---------------------------------------------------------------------------- #
#
# Extract or compute sampling variances from brma object.
#
# Convenience wrapper that returns sei^2 from .outcome_data_sei().
#
# @param object brma object
#
# @return numeric vector of sampling variances
#
# ---------------------------------------------------------------------------- #
.outcome_data_vi <- function(object) {
  return(.outcome_data_sei(object)^2)
}


# ---------------------------------------------------------------------------- #
# .standardized_residuals_loopit
# ---------------------------------------------------------------------------- #
#
# Compute LOO-PIT residuals.
#
# This is the Bayesian equivalent of studentized deleted residuals. It computes
# leave-one-out probability integral transform (PIT) residuals using PSIS weights
# to reweight the posterior samples.
#
# The algorithm:
# 1. Compute LOO via PSIS (reusing loo.brma)
# 2. Extract normalized PSIS weights (S x K matrix)
# 3. Compute CDF matrix F(yi | parameters^(s)) (S x K matrix) for the LOO target
# 4. For each observation i, compute LOO-weighted CDF:
#    u_i = sum_s w_{is} * F(yi | theta^(s))
# 5. Transform to standard normal: r_i = qnorm(u_i)
#
# Unlike traditional standardized residuals, LOO-PIT residuals:
# - Account for estimation uncertainty (integrate over posterior)
# - Account for leverage (PSIS effectively removes yi from posterior)
# - Work for all model types (selection models, GLMMs)
#
# @param object brma object.
#
# @return If loo_only = TRUE, returns the psis_loo object.
#         Otherwise, returns numeric vector of LOO-PIT residuals.
#
# ---------------------------------------------------------------------------- #
.standardized_residuals_loopit <- function(object, psis_weights = NULL, check = TRUE) {

  # extract PSIS object and get normalized weights
  if (is.null(psis_weights)) {
    psis_weights <- loo_weights(object, unit = "estimate")
  }

  # check Pareto k diagnostics and warn if unreliable
  if (check) {
    check_loo(object, unit = "estimate")
  }

  # compute CDF matrix (S x K) for the estimate-unit LOO target
  cdf_matrix <- .cdf_lik_estimate.brma(object)

  # compute LOO-weighted CDF for each observation
# u_i = sum_s w_{is} * F(yi | parameters^(s))
  # this is a weighted average across posterior samples
  u_values <- colSums(psis_weights * cdf_matrix)

  # clamp to avoid infinite quantiles
  u_values <- pmax(pmin(u_values, 1 - 1e-10), 1e-10)

  # transform to standard normal quantiles: r_i = qnorm(u_i)
  resid_point <- stats::qnorm(u_values)

  return(resid_point)
}


# ---------------------------------------------------------------------------- #
# .check_residual_type_availability
# ---------------------------------------------------------------------------- #
#
# Helper function to validate that the requested residual type is compatible
# with the model's outcome type and presence of selection (weightfunction).
#
# @param type Residual type requested
# @param outcome_type Model outcome type
# @param is_weightfunction Whether model includes weightfunction
#
# @return NULL (throws error if invalid combination)
# @keywords internal
#
# ---------------------------------------------------------------------------- #
.check_residual_type_availability <- function(type, outcome_type, is_weightfunction) {
  if (type == "pearson") {
    if (outcome_type != "norm") {
      stop(
        "Pearson residuals are only available for normal outcome models. ",
        "Use type = 'LOO-PIT' for GLMM models.",
        call. = FALSE
      )
    }
    if (is_weightfunction) {
      stop(
        "Pearson residuals are not available for selection models (weightfunction). ",
        "Use type = 'LOO-PIT' for selection models.",
        call. = FALSE
      )
    }
  }

  if (type == "rstandard") {
    if (outcome_type != "norm") {
      stop(
        "Standardized residuals (rstandard) are only available for normal outcome models. ",
        "Use type = 'LOO-PIT' for GLMM models.",
        call. = FALSE
      )
    }
    if (is_weightfunction) {
      stop(
        "Standardized residuals (rstandard) are not available for selection models (weightfunction). ",
        "Use type = 'LOO-PIT' for selection models.",
        call. = FALSE
      )
    }
  }

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .get_model_matrix
# ---------------------------------------------------------------------------- #
#
# Extract the model matrix X from a brma object.
#
# @param object brma object
#
# @return numeric matrix with K rows (one per observation)
#
# ---------------------------------------------------------------------------- #
.get_model_matrix <- function(object) {
  # check if there are moderators
  is_mods  <- .is_mods(object)
  is_PET   <- .is_PET(object)
  is_PEESE <- .is_PEESE(object)

  if (!is_mods) {
    # intercept-only model: X is a column of 1s
    K <- nrow(object[["data"]][["outcome"]])
    X <- matrix(1, nrow = K, ncol = 1)
    colnames(X) <- "(Intercept)"
  } else {
    # for models with moderators, reconstruct the model matrix
    # from the stored data and formula

    # TODO: use BayesTools to obtain the actual adjusted model matrix
    # (i.e., for handling different factor contrasts, etc.)
    mods_data    <- object[["data"]][["mods"]]
    mods_formula <- attr(mods_data, "formula")

    # create the model matrix
    X <- model.matrix(mods_formula, data = mods_data)
  }

  if (is_PET || is_PEESE) {
    outcome_data <- object[["data"]][["outcome"]]
    direction    <- if (.effect_direction(object) == "negative") -1 else 1

    if (is_PET) {
      X <- cbind(X, PET = direction * outcome_data[["sei"]])
    }
    if (is_PEESE) {
      X <- cbind(X, PEESE = direction * outcome_data[["sei"]]^2)
    }
  }

  return(X)
}
