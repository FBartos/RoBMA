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
#'   \item \code{"outcome"} (default): Raw residuals (observed - fitted)
#'   \item \code{"pearson"}: Pearson (semi-standardized) residuals, dividing
#'     raw residuals by the marginal standard error \eqn{\sqrt{v_i + \tau^2}}
#'   \item \code{"rstandard"}: Internally standardized residuals, dividing
#'     raw residuals by their standard errors computed using the hat matrix
#'   \item \code{"quantile"}: Quantile residuals computed via inverse probability
#'     transformation. The CDF value \eqn{F(y_i | \mu_i, \sigma_i)} is computed
#'     for each observation and transformed to standard normal quantiles via
#'     \eqn{\Phi^{-1}(F(y_i))}. For selection models, the weighted normal CDF
#'     is used. Under a correctly specified model, quantile residuals should
#'     follow a standard normal distribution.
#' }
#' @param bias_adjusted whether residuals should be computed from bias-adjusted
#' fitted values. Defaults to \code{FALSE}, which means residuals are computed
#' as the difference between observed values and raw (biased) predictions
#' including PET/PEESE terms. Set to \code{TRUE} to compute residuals from
#' bias-corrected fitted values.
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param as_samples whether posterior samples should be returned instead of
#' a summary table. Defaults to \code{FALSE}.
#' @param ... additional arguments (currently ignored)
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
#' between-study heterogeneity.
#'
#' Standardized residuals (\code{type = "rstandard"}) use the hat matrix to
#' compute residual standard errors that account for the uncertainty in
#' estimated coefficients:
#' \deqn{z_i = \frac{e_i}{\sqrt{[(I-H)M(I-H)']_{ii}}}}
#' where \eqn{H} is the hat matrix and \eqn{M} is the marginal variance-covariance
#' matrix. For models without moderators, this simplifies to the Pearson formula.
#'
#' Quantile residuals (\code{type = "quantile"}) are computed via inverse
#' probability transformation \insertCite{cox1968general}{RoBMA}:
#' \deqn{q_i = \Phi^{-1}(F(y_i | \mu_i, \sigma_i))}
#' where \eqn{F} is the cumulative distribution function of the marginal
#' distribution, \eqn{\sigma_i = \sqrt{v_i + \tau^2}}, and \eqn{\Phi^{-1}}
#' is the standard normal quantile function. For selection models
#' (weightfunction), the weighted normal CDF is used. Under a correctly
#' specified model, quantile residuals should follow a standard normal
#' distribution, making them useful for model diagnostics.
#'
#' For meta-regression models, fitted values incorporate moderator effects.
#' For models without moderators, all fitted values equal the pooled effect.
#'
#' When \code{bias_adjusted = FALSE} (default), the fitted values include
#' PET/PEESE bias terms, so residuals represent deviation from what we
#' expect to observe given publication bias. When \code{bias_adjusted = TRUE},
#' fitted values are bias-corrected, so residuals represent deviation from
#' the estimated true effect.
#'
#' For GLMM models (binomial or Poisson), observed effect sizes and their
#' sampling variances are computed from the raw frequency data using the
#' same formulas as \code{metafor::escalc} with the default zero-cell
#' adjustment (adding 0.5 to all cells when any cell is zero).
#'
#' The standardized residuals are computed separately for each posterior
#' sample of \eqn{\tau}, naturally propagating uncertainty in heterogeneity
#' to the standardized residuals.
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
#' # get raw residuals (default)
#' residuals(fit)
#'
#' # get Pearson (semi-standardized) residuals
#' residuals(fit, type = "pearson")
#'
#' # get internally standardized residuals
#' residuals(fit, type = "rstandard")
#'
#' # get quantile residuals (should be ~N(0,1) if model is correct)
#' residuals(fit, type = "quantile")
#'
#' # get residuals from bias-adjusted predictions
#' residuals(fit, bias_adjusted = TRUE)
#'
#' # get posterior samples
#' samples <- residuals(fit, as_samples = TRUE)
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [predict.brma()], [blup.brma()], [pooled_effect()], [rstandard.brma()]
#' @exportS3Method
residuals.brma <- function(object, type = "outcome", bias_adjusted = FALSE,
                           probs = c(.025, .975), as_samples = FALSE, ...) {

  # input validation
  BayesTools::check_char(type, "type", allow_values = c("outcome", "pearson", "rstandard", "quantile"))
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  BayesTools::check_bool(as_samples, "as_samples")

  # get observed effect sizes
  yi <- .outcome_data_yi(object)
  K  <- length(yi)

  # get fitted values (mu predictions for each observation)
  # type = "terms" gives fixed effects predictions without random effects
  fitted_samples <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = "terms",
    bias_adjusted = bias_adjusted,
    as_samples    = TRUE,
    quiet         = TRUE
  )

  # compute raw residuals: yi - fitted
  # fitted_samples is S x K matrix, yi is vector of length K
  # replicate yi across samples for vectorized subtraction
  S      <- nrow(fitted_samples)
  yi_mat <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)
  resid_samples <- yi_mat - fitted_samples

  # standardize residuals if requested
 if (type == "pearson") {

    resid_samples <- .standardize_residuals_pearson(
      resid_samples = resid_samples,
      object        = object
    )

  } else if (type == "rstandard") {

    resid_samples <- .standardize_residuals_rstandard(
      resid_samples = resid_samples,
      object        = object
    )

  } else if (type == "quantile") {

    resid_samples <- .standardize_residuals_quantile(
      object        = object,
      bias_adjusted = bias_adjusted
    )

  }

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
    title      = switch(
      type,
      "outcome"   = "Residuals:",
      "pearson"   = "Pearson Residuals:",
      "rstandard" = "Standardized Residuals:",
      "quantile"  = "Quantile Residuals:"
    )
  )

  out <- list(
    summary = resid_table,
    data    = object[["data"]]
  )
  class(out) <- "brma.residuals"

  return(out)
}


#' @title Standardized Residuals for brma Objects
#'
#' @description Computes internally standardized residuals from a fitted brma
#' object. This is a convenience wrapper for \code{residuals(object, type = "rstandard", ...)}.
#'
#' @param model a fitted brma object
#' @param ... additional arguments passed to \code{\link{residuals.brma}}
#'
#' @details
#' Standardized residuals use the hat matrix to compute residual standard errors
#' that account for the uncertainty in estimated coefficients:
#' \deqn{z_i = \frac{e_i}{\sqrt{[(I-H)M(I-H)']_{ii}}}}
#' where \eqn{H} is the hat matrix and \eqn{M} is the marginal variance-covariance
#' matrix. For models without moderators, this simplifies to Pearson residuals.
#'
#' @return A \code{brma.residuals} object containing standardized residuals,
#' or a matrix of posterior samples if \code{as_samples = TRUE}.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get standardized residuals
#' rstandard(fit)
#'
#' # equivalent to:
#' residuals(fit, type = "rstandard")
#' }
#'
#' @seealso [residuals.brma()], [predict.brma()]
#' @exportS3Method
rstandard.brma <- function(model, ...) {
  residuals.brma(object = model, type = "rstandard", ...)
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
    ai  <- outcome_data[["ai"]]
    bi  <- outcome_data[["n1i"]] - outcome_data[["ai"]]  # bi = n1i - ai
    ci  <- outcome_data[["ci"]]
    di  <- outcome_data[["n2i"]] - outcome_data[["ci"]]  # di = n2i - ci
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
    ai  <- outcome_data[["ai"]]
    bi  <- outcome_data[["n1i"]] - outcome_data[["ai"]]
    ci  <- outcome_data[["ci"]]
    di  <- outcome_data[["n2i"]] - outcome_data[["ci"]]
    add <- 0.5

    # vectorized computation with conditional zero-cell adjustment
    needs_adj <- (ai == 0 | bi == 0 | ci == 0 | di == 0)
    ai_adj <- ai + add * needs_adj
    bi_adj <- bi + add * needs_adj
    ci_adj <- ci + add * needs_adj
    di_adj <- di + add * needs_adj

    return(sqrt(1/ai_adj + 1/bi_adj + 1/ci_adj + 1/di_adj))

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

    return(sqrt(1/x1i_adj + 1/x2i_adj))

  }
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
# .standardize_residuals_quantile
# ---------------------------------------------------------------------------- #
#
# Compute quantile residuals via inverse probability transformation.
#
# This is the main dispatcher for quantile residuals. It:
# 1. Extracts observed effect sizes and standard errors
# 2. Gets fitted values (mu) and heterogeneity samples (tau)
# 3. Adds study-level random effects for multilevel models
# 4. Dispatches to normal or weighted normal CDF based on model type
#
# Quantile residuals are computed as:
#   q_i = Phi^{-1}(F(y_i | mu_i + gamma_i*tau_between, sigma_i))
# where sigma_i = sqrt(sei^2 + tau_within^2) is the marginal SD at the
# observation level (conditional on study effects).
#
# For selection models with bias_adjusted = FALSE, uses weighted normal CDF.
# Otherwise uses standard normal CDF.
#
# @param object        brma object
# @param bias_adjusted whether to use bias-adjusted fitted values
#
# @return S x K matrix of quantile residual samples
#
# ---------------------------------------------------------------------------- #
.standardize_residuals_quantile <- function(object, bias_adjusted) {

  # get observed effect sizes and standard errors
  yi  <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)
  K   <- length(yi)

  # extract model characteristics
  priors            <- object[["priors"]]
  is_multilevel     <- .is_multilevel(object)
  is_scale          <- .is_scale(object)
  is_weightfunction <- .is_weightfunction(object)
  effect_direction  <- .effect_direction(object)

  # get fitted values (mu predictions for each observation)
  # type = "terms" gives fixed effects predictions without random effects
  fitted_samples <- predict.brma(
    object        = object,
    newdata       = NULL,
    type          = "terms",
    bias_adjusted = bias_adjusted,
    as_samples    = TRUE,
    quiet         = TRUE
  )
  S <- nrow(fitted_samples)

  # get tau samples (heterogeneity) using helper function
  # returns list(tau_within, tau_between) - both S x K matrices
  tau_result <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = object[["data"]][["scale"]],
    scale_formula = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors  = priors[["scale"]],
    is_scale      = is_scale,
    is_multilevel = is_multilevel,
    K             = K
  )
  tau_within_samples  <- tau_result[["tau_within"]]
  tau_between_samples <- tau_result[["tau_between"]]

  # add study-level random effects for multilevel models
  # this makes residuals conditional on study effects (gamma * tau_between)
  if (is_multilevel) {
    fit_data <- .create_fit_data(data = object[["data"]], priors = priors)
    study_contribution <- .evaluate.brma.study_effects(
      fit              = object[["fit"]],
      tau_between      = tau_between_samples,
      study_ids        = fit_data[["study_ids"]],
      same_data        = TRUE,
      effect_direction = effect_direction
    )
    fitted_samples <- fitted_samples + study_contribution
  }

  # compute marginal SD at observation level: sqrt(sei^2 + tau_within^2)
  # note: tau_within is used here because study effects are already in mu
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(sei_mat^2 + tau_within_samples^2)

  # dispatch to appropriate CDF function
  if (!bias_adjusted && is_weightfunction) {
    # use weighted normal CDF for selection models
    quantile_resid <- .standardize_residuals_quantile.wnorm(
      yi               = yi,
      mu_samples       = fitted_samples,
      total_sd         = total_sd,
      object           = object,
      effect_direction = effect_direction
    )
  } else {
    # use standard normal CDF
    quantile_resid <- .standardize_residuals_quantile.norm(
      yi         = yi,
      mu_samples = fitted_samples,
      total_sd   = total_sd
    )
  }

  return(quantile_resid)
}


# ---------------------------------------------------------------------------- #
# .standardize_residuals_quantile.norm
# ---------------------------------------------------------------------------- #
#
# Compute quantile residuals using standard normal CDF.
#
# For each observation, computes:
#   p_i = pnorm(y_i, mu_i, sigma_i)
#   q_i = qnorm(p_i)
#
# @param yi         numeric vector of observed effect sizes
# @param mu_samples S x K matrix of location samples (includes study effects)
# @param total_sd   S x K matrix of marginal SD samples
#
# @return S x K matrix of quantile residual samples
#
# ---------------------------------------------------------------------------- #
.standardize_residuals_quantile.norm <- function(yi, mu_samples, total_sd) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate yi across samples
  yi_mat <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)

  # compute CDF values: F(yi | mu, total_sd)
  p_values <- stats::pnorm(yi_mat, mean = mu_samples, sd = total_sd)

  # transform to standard normal quantiles
  quantile_resid <- stats::qnorm(p_values)

  return(quantile_resid)
}


# ---------------------------------------------------------------------------- #
# .standardize_residuals_quantile.wnorm
# ---------------------------------------------------------------------------- #
#
# Compute quantile residuals using weighted normal CDF for selection models.
#
# For selection models, the observed effect sizes follow a weighted normal
# distribution that accounts for publication bias. The CDF is computed using
# the fast spike-and-slab implementation.
#
# Handles effect direction flipping to match the internal representation where
# selection models operate in "positive" space.
#
# @param yi               numeric vector of observed effect sizes
# @param mu_samples       S x K matrix of location samples (includes study effects)
# @param total_sd         S x K matrix of marginal SD samples
# @param object           brma object
# @param effect_direction "positive" or "negative"
#
# @return S x K matrix of quantile residual samples
#
# ---------------------------------------------------------------------------- #
.standardize_residuals_quantile.wnorm <- function(yi, mu_samples, total_sd, object, effect_direction) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # extract omega samples for weight function
  posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
  omega_samples     <- posterior_samples[, grep("omega", colnames(posterior_samples)), drop = FALSE]

  # get fit_data which contains crit_yi (critical values in "positive" space)
  fit_data <- .create_fit_data(data = object[["data"]], priors = object[["priors"]])

  # flip to positive space if negative effect direction
  if (effect_direction == "negative") {
    yi_for_cdf <- -yi
    mu_for_cdf <- -mu_samples
  } else {
    yi_for_cdf <- yi
    mu_for_cdf <- mu_samples
  }

  # compute CDF values using weighted normal
  # need to loop over observations due to .pwnorm_fast.ss interface
  p_values <- matrix(NA_real_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    p_values[, k] <- .pwnorm_fast.ss(
      q      = rep(yi_for_cdf[k], S),
      mean   = mu_for_cdf[, k],
      sd     = total_sd[, k],
      omega  = omega_samples,
      crit_x = fit_data$crit_yi[, k]
    )
  }

  # transform to standard normal quantiles
  quantile_resid <- stats::qnorm(p_values)

  return(quantile_resid)
}


# ---------------------------------------------------------------------------- #
# .standardize_residuals_pearson
# ---------------------------------------------------------------------------- #
#
# Compute Pearson (semi-standardized) residuals.
#
# Divides raw residuals by the marginal standard error sqrt(vi + tau^2),
# where vi is the sampling variance and tau^2 is the total heterogeneity.
#
# This is computed per posterior sample, naturally propagating uncertainty
# in tau to the standardized residuals.
#
# @param resid_samples S x K matrix of raw residual samples
# @param object        brma object
# @param outcome_type  character; "norm", "bin", or "pois"
#
# @return S x K matrix of Pearson residual samples
#
# ---------------------------------------------------------------------------- #
.standardize_residuals_pearson <- function(resid_samples, object) {

  S  <- nrow(resid_samples)
  K  <- ncol(resid_samples)

  # get sampling variances
  vi <- .outcome_data_vi(object)

  # get tau samples (heterogeneity)
  # use predict.brma with type = "terms.scale" to get tau samples
  tau_samples <- predict.brma(
    object     = object,
    newdata    = NULL,
    type       = "terms.scale",
    as_samples = TRUE,
    quiet      = TRUE
  )

  # tau_samples is S x K matrix (one column per observation for location-scale models,
  # or all columns identical for non-scale models)

  # compute marginal variance: vi + tau^2
  # vi is vector of length K, tau_samples is S x K matrix
  vi_mat <- matrix(vi, nrow = S, ncol = K, byrow = TRUE)
  marginal_var <- vi_mat + tau_samples^2

  # compute Pearson residuals: ei / sqrt(vi + tau^2)
  marginal_se <- sqrt(marginal_var)
  pearson_resid <- resid_samples / marginal_se

  return(pearson_resid)
}


# ---------------------------------------------------------------------------- #
# .standardize_residuals_rstandard
# ---------------------------------------------------------------------------- #
#
# Compute internally standardized residuals using the hat matrix.
#
# Following metafor::rstandard.rma.uni:
#   H = X (X'WX)^{-1} X' W     (hat matrix)
#   ImH = I - H
#   ve = ImH M ImH'            (variance-covariance of residuals)
#   SE_i = sqrt(ve_ii)         (standard error of each residual)
#   z_i = e_i / SE_i           (standardized residual)
#
# where M is the marginal variance-covariance matrix (diagonal with vi + tau^2).
#
# For models without moderators, this simplifies to Pearson residuals since
# the hat matrix correction has minimal effect on the diagonal.
#
# This is computed per posterior sample, naturally propagating uncertainty
# in tau to the standardized residuals.
#
# @param resid_samples S x K matrix of raw residual samples
# @param object        brma object
# @param outcome_type  character; "norm", "bin", or "pois"
#
# @return S x K matrix of standardized residual samples
#
# ---------------------------------------------------------------------------- #
.standardize_residuals_rstandard <- function(resid_samples, object) {

  S  <- nrow(resid_samples)
  K  <- ncol(resid_samples)

  # get sampling variances
  vi <- .outcome_data_vi(object)

  # get tau samples (heterogeneity)
  tau_samples <- predict.brma(
    object     = object,
    newdata    = NULL,
    type       = "terms.scale",
    as_samples = TRUE,
    quiet      = TRUE
  )

  # get model matrix X
  X <- .get_model_matrix(object)

  # for intercept-only models, use Pearson formula (hat matrix correction is minimal)
  if (ncol(X) == 1 && all(X == 1)) {
    return(.standardize_residuals_pearson(resid_samples, object))
  }

  # for models with moderators, compute full hat matrix correction per sample
  # this is more expensive but accounts for the projection
  rstandard_samples <- matrix(0, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    # get tau^2 for this sample (use first column, or column-specific for scale models)
    tau2_s <- tau_samples[s, ]^2

    # construct marginal variance M = diag(vi + tau^2)
    M <- diag(vi + tau2_s)

    # weight matrix W = M^{-1}
    W <- diag(1 / (vi + tau2_s))

    # hat matrix H = X (X'WX)^{-1} X' W
    XtWX <- crossprod(X, W %*% X)
    XtWX_inv <- tryCatch(
      solve(XtWX),
      error = function(e) MASS::ginv(XtWX)
    )
    H <- X %*% XtWX_inv %*% crossprod(X, W)

    # I - H
    ImH <- diag(K) - H

    # variance-covariance of residuals: (I-H) M (I-H)'
    ve <- ImH %*% tcrossprod(M, ImH)

    # standard errors of residuals
    se_resid <- sqrt(diag(ve))

    # standardized residuals
    rstandard_samples[s, ] <- resid_samples[s, ] / se_resid
  }

  return(rstandard_samples)
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
  is_mods <- .is_mods(object)

  if (!is_mods) {
    # intercept-only model: X is a column of 1s
    K <- nrow(object[["data"]][["outcome"]])
    return(matrix(1, nrow = K, ncol = 1))
  }

  # for models with moderators, reconstruct the model matrix
  # from the stored data and formula

  # TODO: use BayesTools to obtain the actual adjusted model matrix
  # (i.e., for handling different factor contrasts, etc.)
  mods_data    <- object[["data"]][["mods"]]
  mods_formula <- attr(mods_data, "formula")

  # create the model matrix
  X <- model.matrix(mods_formula, data = mods_data)

  return(X)
}
