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
#'   \item \code{"LOO-PIT"} (alias: \code{"rstudent"}): Leave-one-out probability
#'     integral transform (PIT) residuals computed via Pareto smoothed importance
#'     sampling. The LOO-CDF value for each observation is computed and transformed
#'     to standard normal quantiles via \eqn{\Phi^{-1}(u_i)}. Under a correctly
#'     specified model, these residuals should follow a standard normal distribution.
#'     This is the recommended standardized residual for Bayesian models as it properly
#'     accounts for estimation uncertainty and leverage. Available for normal
#'     outcome models. Binomial and Poisson models require an explicit discrete
#'     PIT convention, which is not implemented. Note: This requires that the loo
#'     has been computed previously (see \code{add_loo()} function).
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
#' bias-corrected fitted values. Applies to outcome and Pearson residuals. Note
#' that bias-adjusted residuals are not residuals in the traditional sense.
#' @param max_samples maximum posterior draws used for known-\code{V}
#' \code{brma.mv()} marginal covariance computations in Pearson and internally
#' standardized residuals. Defaults to \code{Inf}. Finite values
#' deterministically thin draws across posterior row order.
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
#' distribution used by LOO. Lower and upper probabilities are accumulated
#' separately on the log scale, and the smaller tail is passed directly to the
#' normal quantile function. Therefore extreme residuals are not truncated by
#' probability clipping. Under a correctly specified model, LOO-PIT
#' residuals should follow a standard normal distribution. Unlike traditional
#' standardized residuals, LOO-PIT residuals properly account for estimation
#' uncertainty and leverage without requiring a hat matrix. This is the
#' recommended method for standardized residuals in Bayesian meta-analysis.
#'
#' For meta-regression models, fitted values incorporate moderator effects.
#' For models without moderators, all fitted values equal the pooled effect.
#'
#' LOO-PIT is not available for binomial or Poisson GLMMs. Their fitted
#' predictive distributions are discrete, so a PIT diagnostic requires an
#' explicit randomized, mid-P, or other discrete convention. The former normal
#' approximation on a derived effect-size scale is not a PIT for the fitted
#' count likelihood.
#'
#' For correlated known-\code{V} \code{brma.mv()} models, LOO-PIT residuals use
#' the same conditional estimate-unit target as LOO,
#' \eqn{p(y_i \mid y_{-i}, \theta)}. Random-formula \code{brma.mv()} models
#' support marginal Pearson and internally standardized residuals through the
#' marginal covariance \eqn{V + ZGZ'} and estimate-depth residuals through the
#' fitted existing-level random-effect target.
#'
#' The residuals are computed separately for each posterior sample,
#' naturally propagating uncertainty in model parameters to the residuals.
#'
#' @return A numeric vector of residual means, one per estimate.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit <- bPET(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   # raw residuals (default)
#'   residuals(fit)
#'
#'   # Pearson and internally standardized residuals
#'   residuals(fit, type = "pearson")
#'   residuals(fit, type = "rstandard")
#'
#'   # LOO-PIT residuals require stored LOO
#'   fit <- add_loo(fit)
#'   residuals(fit, type = "LOO-PIT")
#'
#'   # residuals from bias-adjusted predictions
#'   residuals(fit, bias_adjusted = TRUE)
#'
#'   # check Pareto k diagnostics
#'   plot(loo(fit))
#' }
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [predict.brma()], [blup.brma()], [pooled_effect()], [rstandard.brma()], [rstudent.brma()]
#' @exportS3Method
residuals.brma <- function(object, type = "outcome", unit = "estimate",
                           conditioning_depth = "marginal",
                           bias_adjusted = FALSE, max_samples = Inf, ...) {

  # input validation
  dots                         <- list(...)
  conditioning_depth_specified <- !missing(conditioning_depth)
  .check_legacy_level_arg(dots, "residuals()")

  BayesTools::check_char(type, "type", allow_values = c("outcome", "pearson", "rstandard", "rstudent", "LOO-PIT"))
  unit                         <- .normalize_unit(unit)
  conditioning_depth           <- .normalize_conditioning_depth(conditioning_depth)
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  max_samples                  <- .normalize_max_samples(max_samples, "max_samples")

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
  .check_residual_random_formula_availability(
    object             = object,
    type               = type,
    conditioning_depth = conditioning_depth,
    caller             = "residuals()"
  )

  # extract model characteristics for error checking
  outcome_type      <- .outcome_type(object)
  is_weightfunction <- .is_weightfunction(object)

  # check for unsupported type/model combinations
  .check_residual_type_availability(type, outcome_type, is_weightfunction)

  out <- .residuals_estimate.brma(
    object             = object,
    type               = type,
    conditioning_depth = conditioning_depth,
    bias_adjusted      = bias_adjusted,
    max_samples        = max_samples
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
                                     bias_adjusted, max_samples = Inf) {

  if (is.element(type, c("LOO-PIT", "rstudent"))) {
    return(.standardized_residuals_loopit(object))
  }

  if (type == "rstandard") {
    standard <- rstandard.brma(
      model              = object,
      unit               = "estimate",
      conditioning_depth = conditioning_depth,
      max_samples        = max_samples
    )
    out <- standard[["z"]]
    if (!is.null(attr(standard, "known_v_diagnostic"))) {
      out <- .known_v_attach_diagnostic_metadata(
        out,
        attr(standard, "known_v_diagnostic")
      )
    }
    return(out)
  }

  known_v_sample_info <- NULL
  known_v_samples     <- NULL
  if (type == "pearson" && .is_data_known_v(object[["data"]])) {
    known_v_sample_info <- .known_v_diagnostic_posterior_samples(
      object      = object,
      max_samples = max_samples,
      caller      = "known-V Pearson residual diagnostics"
    )
    known_v_samples <- known_v_sample_info[["posterior_samples"]]
  }

  if (.is_data_known_v(object[["data"]]) && conditioning_depth == "estimate") {
    setup          <- .estimate_likelihood_setup.brma(
      object            = object,
      bias_adjusted     = bias_adjusted,
      posterior_samples = known_v_samples
    )
    fitted_samples <- .known_v_estimate_blup_from_setup(setup)
  } else {
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
      quiet         = TRUE,
      .posterior_samples = known_v_samples
    )
  }

  yi <- .outcome_data_yi(object)
  K  <- length(yi)
  S  <- nrow(fitted_samples)

  yi_mat        <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)
  resid_samples <- yi_mat - fitted_samples

  if (type == "pearson") {
    se_samples    <- .pearson_residual_se_samples(
      object             = object,
      conditioning_depth = conditioning_depth,
      max_samples        = max_samples,
      posterior_samples  = known_v_samples,
      sample_info        = known_v_sample_info
    )
    resid_samples <- resid_samples / se_samples
    known_v_metadata <- attr(se_samples, "known_v_diagnostic")
  }

  out <- colMeans(resid_samples)
  if (exists("known_v_metadata", inherits = FALSE) &&
      !is.null(known_v_metadata)) {
    out <- .known_v_attach_diagnostic_metadata(out, known_v_metadata)
  }

  return(out)
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
.pearson_residual_se_samples <- function(object, conditioning_depth,
                                         max_samples = Inf,
                                         posterior_samples = NULL,
                                         sample_info = NULL) {

  if (.is_data_known_v(object[["data"]])) {
    return(.known_v_pearson_residual_se_samples(
      object             = object,
      conditioning_depth = conditioning_depth,
      max_samples        = max_samples,
      posterior_samples  = posterior_samples,
      sample_info        = sample_info
    ))
  }

  is_multilevel <- .is_multilevel(object)
  is_scale      <- .is_scale(object)
  priors        <- object[["priors"]]
  vi            <- .outcome_data_vi(object)
  weights       <- .outcome_data_weights(object)
  K             <- length(vi)

  tau_result <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = object[["data"]][["scale"]],
    scale_formula = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors  = priors[["scale"]],
    is_scale      = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    fixed_tau         = .fixed_tau_prior_value(priors),
    fixed_rho         = .fixed_rho_prior_value(priors)
  )

  tau_within  <- tau_result[["tau_within"]]
  tau_between <- tau_result[["tau_between"]]
  if (is.null(tau_between)) {
    tau_between <- matrix(0, nrow = nrow(tau_within), ncol = K)
  }

  vi_mat      <- matrix(vi, nrow = nrow(tau_within), ncol = K, byrow = TRUE)
  weights_mat <- matrix(weights, nrow = nrow(tau_within), ncol = K, byrow = TRUE)

  se2 <- switch(conditioning_depth,
    "marginal" = (vi_mat + tau_within^2) / weights_mat + tau_between^2,
    "cluster"  = (vi_mat + tau_within^2) / weights_mat,
    "estimate" = vi_mat / weights_mat
  )

  return(sqrt(se2))
}


# ---------------------------------------------------------------------------- #
# .known_v_pearson_residual_se_samples
# ---------------------------------------------------------------------------- #
#
# Posterior standard errors for known-V Pearson residuals.
#
# ---------------------------------------------------------------------------- #
.known_v_pearson_residual_se_samples <- function(object, conditioning_depth,
                                                 max_samples = Inf,
                                                 posterior_samples = NULL,
                                                 sample_info = NULL) {

  if (is.null(sample_info)) {
    sample_info <- .known_v_diagnostic_posterior_samples(
      object            = object,
      posterior_samples = posterior_samples,
      max_samples       = max_samples,
      caller            = "known-V Pearson residual diagnostics"
    )
  }
  setup           <- .estimate_likelihood_setup.brma(
    object            = object,
    posterior_samples = sample_info[["posterior_samples"]]
  )
  known_V         <- .data_known_v_data(setup[["data"]])
  sampling_diag   <- .known_v_diagonal(known_V)
  sampling_matrix <- matrix(
    sampling_diag,
    nrow  = setup[["S"]],
    ncol  = setup[["K"]],
    byrow = TRUE
  )

  se2 <- if (conditioning_depth == "estimate") {
    sampling_matrix
  } else if (conditioning_depth == "marginal") {
    out        <- matrix(NA_real_, nrow = setup[["S"]], ncol = setup[["K"]])
    chunk_info <- .known_v_apply_marginal_covariance_chunks(
      object            = object,
      posterior_samples = setup[["posterior_samples"]],
      FUN               = function(covariance_samples, rows) {
        for (j in seq_along(rows)) {
          out[rows[j], ] <<- diag(covariance_samples[j, , ])
        }
      }
    )
    out
  } else {
    .check_cluster_unit_deferred("residuals()", argument = "conditioning_depth")
  }

  out <- sqrt(se2)
  out <- .known_v_attach_diagnostic_metadata(
    out,
    .known_v_diagnostic_metadata(
      sample_info = sample_info,
      chunk_info  = if (exists("chunk_info", inherits = FALSE)) chunk_info else NULL
    )
  )

  return(out)
}


#' @title Internally Standardized Residuals for brma Objects
#'
#' @description Computes internally standardized residuals from a fitted brma
#' object using the hat matrix or the corresponding known-\code{V} GLS/BLUP
#' projection. Returns a data frame with raw residuals, standard errors, and
#' standardized residuals (z-values). Available for normal outcome models only.
#'
#' @param model a fitted brma object.
#' @param unit output unit. Only \code{"estimate"} is implemented currently.
#' @param conditioning_depth conditioning depth. Options are:
#' \itemize{
#'   \item \code{"marginal"} (default): Residuals from fixed effects predictions
#'     (\eqn{observed - X\beta}).
#'   \item \code{"cluster"}: Residuals from cluster-level predictions (observed -
#'     (\eqn{X\beta + gamma})). Only available for multilevel (3-level) models.
#'   \item \code{"estimate"}: Residuals from BLUPs,
#'     i.e., deviations of the observed effect sizes from the best linear
#'     unbiased predictions of the estimate-specific true effects (observed - theta).
#' }
#' @param ... additional arguments (currently ignored)
#' @param max_samples maximum posterior draws used for known-\code{V}
#' \code{brma.mv()} marginal covariance computations. Defaults to \code{Inf}.
#' Finite values deterministically thin draws across posterior row order.
#'
#' @details
#' This function returns a data frame with three columns matching the output
#' of \code{metafor::rstandard}:
#' \itemize{
#'   \item \code{resid}: Raw residuals (observed - fitted values)
#'   \item \code{se}: Standard errors of the residuals
#'   \item \code{z}: Posterior mean of draw-wise standardized residuals
#' }
#'
#' Internally standardized residuals are computed draw by draw by dividing the
#' observed residual by its corresponding residual standard error from the hat
#' matrix. The returned \code{z} column is the posterior mean of those draw-wise
#' ratios, so it is not generally equal to the displayed \code{resid / se}
#' summaries. For correlated known-\code{V} models, marginal residuals use the
#' GLS residual covariance and estimate-depth residuals use the BLUP-style
#' sampling residual projection. For a correctly specified model, these
#' residuals should approximately follow a standard normal distribution.
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
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   # marginal internally standardized residuals (default)
#'   rstandard(fit)
#'
#'   # estimate-level (BLUP-based) residuals
#'   rstandard(fit, conditioning_depth = "estimate")
#' }
#' }
#'
#' @seealso [rstudent.brma()], [residuals.brma()], [blup.brma()], [predict.brma()]
#' @exportS3Method
rstandard.brma <- function(model, unit = "estimate",
                           conditioning_depth = "marginal",
                           max_samples = Inf, ...) {

  dots <- list(...)
  .check_legacy_level_arg(dots, "rstandard()")

  unit               <- .normalize_unit(unit)
  conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)
  max_samples        <- .normalize_max_samples(max_samples, "max_samples")
  .check_unit_conditioning_depth(
    object             = model,
    unit               = unit,
    conditioning_depth = conditioning_depth,
    caller             = "rstandard()"
  )

  if (unit == "cluster") {
    .check_cluster_unit_deferred("rstandard()")
  }
  .check_residual_random_formula_availability(
    object             = model,
    type               = "rstandard",
    conditioning_depth = conditioning_depth,
    caller             = "rstandard()"
  )

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
    summarize          = TRUE,
    max_samples        = max_samples
  )

  # construct output data frame matching metafor::rstandard format
  out <- data.frame(
    resid = hat_res[["resid"]],
    se    = hat_res[["se"]],
    z     = hat_res[["z"]]
  )
  out <- .diagnostic_set_rownames(out, model)
  if (!is.null(hat_res[["known_v_diagnostic"]])) {
    out <- .known_v_attach_diagnostic_metadata(
      out,
      hat_res[["known_v_diagnostic"]]
    )
  }

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
#' are computed from the fitted selected-normal predictive distribution.
#'
#' Unlike \code{\link{rstandard.brma}} (which uses the hat matrix), LOO-PIT
#' residuals properly account for estimation uncertainty and leverage without
#' requiring explicit hat matrix computation. Binomial and Poisson GLMMs are
#' unavailable until a discrete PIT convention is defined.
#' For correlated known-\code{V} \code{brma.mv()} models, the PIT and companion
#' residual moments are computed from the Schur-complement conditional
#' predictive distribution used by estimate-unit LOO.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{resid}: Raw residuals
#'   \item \code{se}: Standard errors of the residuals
#'   \item \code{z}: Externally standardized residuals (LOO-PIT)
#' }
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit <- add_loo(fit)
#'
#'   # externally standardized residuals
#'   rstudent(fit)
#'
#'   # check Pareto k values
#'   plot(loo(fit))
#' }
#' }
#'
#' @seealso [rstandard.brma()], [residuals.brma()], [loo.brma()], [blup.brma()]
#' @exportS3Method
rstudent.brma <- function(model, unit = "estimate",
                          conditioning_depth = "marginal", ...) {

  dots                         <- list(...)
  .psis_context                <- dots[[".psis_context"]]
  dots[[".psis_context"]]      <- NULL
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
  .check_residual_random_formula_availability(
    object             = model,
    type               = "rstudent",
    conditioning_depth = "estimate",
    caller             = "rstudent()"
  )

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
  .check_residual_type_availability(
    type              = "rstudent",
    outcome_type      = setup[["outcome_type"]],
    is_weightfunction = setup[["is_weightfunction"]]
  )

  # extract PSIS object once and reuse it for PIT and LOO expectations
  psis_context <- .diagnostic_psis_context(model, .psis_context)
  psis_object  <- psis_context[["psis_object"]]
  psis_weights <- psis_context[["psis_weights"]]

  # get LOO-PIT z values using the estimate-unit LOO target
  .diagnostic_check_loo(model, psis_context, unit = "estimate")

  if (setup[["outcome_type"]] == "norm" && setup[["is_weightfunction"]]) {
    summary <- .loo_predictive_selnorm_summary_estimate(
      object       = model,
      setup        = setup,
      psis_weights = psis_weights
    )

    out <- data.frame(
      resid = summary[["resid"]],
      se    = summary[["se"]],
      z     = summary[["z"]]
    )
    out <- .diagnostic_set_rownames(out, model)

    return(out)
  }

  z <- .standardized_residuals_loopit(
    object       = model,
    psis_weights = psis_weights,
    check        = FALSE,
    setup        = setup
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
  out <- .diagnostic_set_rownames(out, model)

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
.weighted_predictive_moments <- function(mean, variance, weights) {

  mean     <- as.matrix(mean)
  variance <- as.matrix(variance)
  weights  <- as.matrix(weights)
  if (!identical(dim(mean), dim(variance)) ||
      !identical(dim(mean), dim(weights)) ||
      any(!is.finite(mean)) || any(!is.finite(variance)) ||
      any(variance < 0) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("Predictive moments and weights are invalid.", call. = FALSE)
  }

  weight_sum <- colSums(weights)
  if (any(!is.finite(weight_sum)) || any(weight_sum <= 0)) {
    stop("Predictive weights must have positive finite column sums.",
         call. = FALSE)
  }
  weights  <- sweep(weights, 2L, weight_sum, "/")
  location <- colSums(weights * mean)
  centered <- sweep(mean, 2L, location, "-")
  variance <- colSums(weights * (variance + centered^2))
  if (any(!is.finite(variance)) || any(variance < 0)) {
    stop("Predictive variance is not finite and non-negative.", call. = FALSE)
  }

  return(list(
    mean     = location,
    variance = variance,
    se       = sqrt(variance)
  ))
}


.loo_predictive_moments_estimate <- function(object, setup, psis_object,
                                             psis_weights = NULL) {

  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  K                 <- setup[["K"]]
  mu_samples        <- setup[["mu"]]
  tau_within        <- setup[["tau_within"]]
  outcome_type      <- setup[["outcome_type"]]
  is_weightfunction <- setup[["is_weightfunction"]]

  if (is.null(psis_weights)) {
    psis_weights <- loo::weights.importance_sampling(
      psis_object,
      log       = FALSE,
      normalize = TRUE
    )
  }

  if (.known_v_estimate_target_uses_backend(setup[["data"]])) {
    summary <- .known_v_estimate_target_summary_from_setup(setup)
    moments <- .weighted_predictive_moments(
      mean     = summary[["mean"]],
      variance = summary[["variance"]],
      weights  = psis_weights
    )

    return(list(
      resid = yi - moments[["mean"]],
      se    = moments[["se"]]
    ))
  }

  sei_mat <- matrix(sei, nrow = nrow(mu_samples), ncol = K, byrow = TRUE)

  if (outcome_type == "norm" && is_weightfunction) {
    summary <- .loo_predictive_selnorm_summary_estimate(
      object       = object,
      setup        = setup,
      psis_weights = psis_weights
    )

    return(list(
      resid = summary[["resid"]],
      se    = summary[["se"]]
    ))
  } else {
    variance_samples <- tau_within^2 + sei_mat^2
    mean_samples     <- mu_samples
  }

  moments <- .weighted_predictive_moments(
    mean     = mean_samples,
    variance = variance_samples,
    weights  = psis_weights
  )

  return(list(
    resid = yi - moments[["mean"]],
    se    = moments[["se"]]
  ))
}


# ---------------------------------------------------------------------------- #
# .loo_predictive_selnorm_summary_estimate
# ---------------------------------------------------------------------------- #
#
# Fused selected-normal LOO-PIT and predictive moment reducer.
#
# ---------------------------------------------------------------------------- #
.loo_predictive_selnorm_summary_estimate <- function(object, setup,
                                                     psis_weights) {

  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  K                 <- setup[["K"]]
  mu_samples        <- setup[["mu"]]
  tau_within        <- setup[["tau_within"]]
  posterior_samples <- setup[["posterior_samples"]]
  S                 <- nrow(mu_samples)
  sei_mat           <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd          <- sqrt(tau_within^2 + sei_mat^2)

  selection_context <- .selection_context(
    object            = object,
    posterior_samples = posterior_samples
  )
  summary <- .selection_step_weighted_summary(
    yi                = yi,
    mean              = mu_samples,
    sd                = total_sd,
    sei               = sei,
    psis_weights      = psis_weights,
    selection_context = selection_context
  )

  resid <- yi - summary[["mean"]]

  return(list(
    resid = resid,
    se    = sqrt(summary[["variance"]]),
    z     = .loo_pit_z_from_log_probabilities(
      summary[["log_lower"]],
      summary[["log_upper"]]
    )
  ))
}


# ---------------------------------------------------------------------------- #
# .outcome_moments.selnorm
# ---------------------------------------------------------------------------- #
#
# Per-posterior-sample first and second moments for selected-normal kernels.
#
# ---------------------------------------------------------------------------- #
.outcome_moments.selnorm <- function(mu_samples, tau_within, sei,
                                     selection_context) {

  S        <- nrow(mu_samples)
  K        <- ncol(mu_samples)
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  return(.selection_step_moments_matrix(
    mean              = mu_samples,
    sd                = total_sd,
    sei               = sei,
    selection_context = selection_context
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
.standardized_residuals_loopit <- function(object, psis_weights = NULL,
                                           check = TRUE, setup = NULL) {

  # extract PSIS object and get normalized weights
  if (is.null(psis_weights)) {
    psis_weights <- loo_weights(object, unit = "estimate")
  }

  # check Pareto k diagnostics and warn if unreliable
  if (check) {
    check_loo(object, unit = "estimate")
  }

  if (is.null(setup)) {
    setup <- .estimate_likelihood_setup.brma(object)
  }
  if (setup[["outcome_type"]] == "norm" && setup[["is_weightfunction"]]) {
    summary <- .loo_predictive_selnorm_summary_estimate(
      object       = object,
      setup        = setup,
      psis_weights = psis_weights
    )
    return(summary[["z"]])
  }

  log_tails <- .loo_predictive_log_tails_estimate(setup)
  return(.loo_pit_z_from_log_tails(
    log_lower    = log_tails[["log_lower"]],
    log_upper    = log_tails[["log_upper"]],
    psis_weights = psis_weights
  ))
}


# Compute paired log tails for the normal estimate-unit LOO target.
.loo_predictive_log_tails_estimate <- function(setup) {

  if (!identical(setup[["outcome_type"]], "norm")) {
    stop(
      "LOO-PIT residuals are not available for binomial or Poisson GLMMs ",
      "because a discrete PIT convention has not been defined.",
      call. = FALSE
    )
  }

  if (.known_v_estimate_target_uses_backend(setup[["data"]])) {
    summary <- .known_v_estimate_target_summary_from_setup(setup)
    return(summary[c("log_lower", "log_upper")])
  }

  yi_mat <- matrix(
    setup[["yi"]],
    nrow  = setup[["S"]],
    ncol  = setup[["K"]],
    byrow = TRUE
  )
  sei_mat <- matrix(
    setup[["sei"]],
    nrow  = setup[["S"]],
    ncol  = setup[["K"]],
    byrow = TRUE
  )
  total_sd <- sqrt(setup[["tau_within"]]^2 + sei_mat^2)

  return(list(
    log_lower = stats::pnorm(
      yi_mat,
      mean  = setup[["mu"]],
      sd    = total_sd,
      log.p = TRUE
    ),
    log_upper = stats::pnorm(
      yi_mat,
      mean       = setup[["mu"]],
      sd         = total_sd,
      lower.tail = FALSE,
      log.p      = TRUE
    )
  ))
}


# Transform paired log probabilities without materializing the larger tail.
.loo_pit_z_from_log_probabilities <- function(log_lower, log_upper) {

  if (!is.numeric(log_lower) || !is.numeric(log_upper) ||
      length(log_lower) != length(log_upper) ||
      anyNA(log_lower) || anyNA(log_upper) ||
      any(log_lower > 0) || any(log_upper > 0)) {
    stop("Internal error: invalid paired log PIT probabilities.", call. = FALSE)
  }

  use_lower <- log_lower <= log_upper
  z         <- numeric(length(log_lower))
  z[use_lower] <- stats::qnorm(
    log_lower[use_lower],
    log.p = TRUE
  )
  z[!use_lower] <- stats::qnorm(
    log_upper[!use_lower],
    lower.tail = FALSE,
    log.p      = TRUE
  )
  return(z)
}


# PSIS-average paired per-draw log tails and transform them to z scores.
.loo_pit_z_from_log_tails <- function(log_lower, log_upper, psis_weights) {

  log_lower    <- as.matrix(log_lower)
  log_upper    <- as.matrix(log_upper)
  psis_weights <- as.matrix(psis_weights)
  if (!identical(dim(log_lower), dim(log_upper)) ||
      !identical(dim(log_lower), dim(psis_weights)) ||
      anyNA(log_lower) || anyNA(log_upper) ||
      any(log_lower > 0) || any(log_upper > 0) ||
      anyNA(psis_weights) || any(!is.finite(psis_weights)) ||
      any(psis_weights < 0) || any(colSums(psis_weights) <= 0)) {
    stop("Internal error: invalid LOO-PIT tails or PSIS weights.", call. = FALSE)
  }

  log_weights <- log(psis_weights)
  log_average <- function(log_values) {

    terms     <- log_weights + log_values
    maxima    <- apply(terms, 2L, max)
    out       <- rep(-Inf, ncol(terms))
    is_finite <- is.finite(maxima)
    out[is_finite] <- maxima[is_finite] + log(colSums(exp(sweep(
      terms[, is_finite, drop = FALSE],
      2L,
      maxima[is_finite],
      "-"
    ))))
    return(out)
  }

  return(.loo_pit_z_from_log_probabilities(
    log_average(log_lower),
    log_average(log_upper)
  ))
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
  if (is.element(type, c("LOO-PIT", "rstudent")) &&
      outcome_type != "norm") {
    stop(
      "LOO-PIT residuals are not available for binomial or Poisson GLMMs ",
      "because a discrete PIT convention has not been defined.",
      call. = FALSE
    )
  }

  if (type == "pearson") {
    if (outcome_type != "norm") {
      stop(
        "Pearson residuals are only available for normal outcome models.",
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
        "Standardized residuals (rstandard) are only available for normal outcome models.",
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
# .check_residual_random_formula_availability
# ---------------------------------------------------------------------------- #
#
# Random-formula post-fit residual diagnostics are implemented only for the
# known-V brma.mv() estimate-unit targets where the likelihood target is defined.
#
# ---------------------------------------------------------------------------- #
.check_residual_random_formula_availability <- function(object, type,
                                                        conditioning_depth,
                                                        caller) {

  if (!.is_random(object)) {
    return(invisible(TRUE))
  }

  if (!(inherits(object, "brma.mv") && .is_data_known_v(object[["data"]]))) {
    .check_random_formula_postfit_deferred(object, caller)
  }

  if (type == "outcome") {
    return(invisible(TRUE))
  }

  if (is.element(type, c("LOO-PIT", "rstudent"))) {
    return(invisible(TRUE))
  }

  if (conditioning_depth == "cluster") {
    stop(
      caller,
      " for brma.mv() random-formula known-V models is available only with ",
      "conditioning_depth = 'marginal' or 'estimate'. Cluster-depth ",
      "standardization is not implemented yet.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}
