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
#' between-study heterogeneity. Only available for normal outcome models
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
#' cumulative distribution function of the marginal distribution. Under a
#' correctly specified model, LOO-PIT residuals should follow a standard
#' normal distribution. Unlike traditional standardized residuals, LOO-PIT
#' residuals properly account for estimation uncertainty and leverage without
#' requiring a hat matrix. This is the recommended method for standardized
#' residuals in Bayesian meta-analysis.
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
#' @return A numeric vector of residual means, one per observation.
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
residuals.brma <- function(object, type = "outcome", bias_adjusted = FALSE, ...) {

  # input validation
  BayesTools::check_char(type, "type", allow_values = c("outcome", "pearson", "rstandard", "rstudent", "LOO-PIT"))
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")

  # extract model characteristics for error checking
  outcome_type      <- .outcome_type(object)
  is_weightfunction <- .is_weightfunction(object)

  # check for unsupported type/model combinations
  .check_residual_type_availability(type, outcome_type, is_weightfunction)

  # dispatch to LOO-PIT method (separate path due to different computation)
  if (is.element(type, c("LOO-PIT", "rstudent"))) {
    out <- .standardized_residuals_loopit(object)
  } else {
    # get fitted values (mu predictions for each observation)
    # type = "terms" gives fixed effects predictions estimate-level random effects
    fitted_samples <- predict.brma(
      object        = object,
      newdata       = NULL,
      type          = "terms",
      bias_adjusted = bias_adjusted,
      quiet         = TRUE
    )

    # get observed effect sizes
    yi <- .outcome_data_yi(object)
    K  <- length(yi)

    # compute raw residuals: yi - fitted
    # fitted_samples is S x K matrix, yi is vector of length K
    # replicate yi across samples for vectorized subtraction
    S <- nrow(fitted_samples)
    yi_mat        <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)
    resid_samples <- yi_mat - fitted_samples

    # standardize residuals if requested
    if (type == "pearson" || type == "rstandard") {

      if (type == "rstandard") {
        # use hat matrix SE
        hat_res       <- .compute_hat_matrix_samples(object, type = "marginal", return_se = TRUE)
        se_samples    <- hat_res[["se"]]
        resid_samples <- resid_samples / se_samples
      } else if (type == "pearson") {
        # use marginal SE without hat matrix
        # SE = sqrt(vi + tau_total^2)
        vi <- .outcome_data_vi(object)
        tau_samples <- predict.brma(
          object     = object,
          newdata    = NULL,
          type       = "terms.scale",
          quiet      = TRUE
        )
        vi_mat        <- matrix(vi, nrow = S, ncol = K, byrow = TRUE)
        se_samples    <- sqrt(vi_mat + tau_samples^2)
        resid_samples <- resid_samples / se_samples
      }
    }

    # compute posterior mean residuals
    out <- colMeans(resid_samples)
  }

  # clean names
  names(out) <- NULL
  return(out)
}


#' @title Internally Standardized Residuals for brma Objects
#'
#' @description Computes internally standardized residuals from a fitted brma
#' object using the hat matrix. Returns a data frame with raw residuals,
#' standard errors, and standardized residuals (z-values). Available for normal
#' outcome models only.
#'
#' @param model a fitted brma object
#' @param type the type of residuals to compute. Options are:
#' \itemize{
#'   \item \code{"marginal"} (default): Residuals from fixed effects predictions
#'     (observed - Xβ).
#'   \item \code{"study"}: Residuals from study-level predictions (observed -
#'     (Xβ + gamma)). Only available for multilevel (3-level) models.
#'   \item \code{"estimate"} (alias: \code{"conditional"}): Residuals from BLUPs,
#'     i.e., deviations of the observed effect sizes from the best linear
#'     unbiased predictions of the study-specific true effects (observed - theta).
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
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get marginal internally standardized residuals (default)
#' rstandard(fit)
#'
#' # get estimate-level (BLUP-based) residuals
#' rstandard(fit, type = "estimate")
#' }
#'
#' @seealso [rstudent.brma()], [residuals.brma()], [blup.brma()], [predict.brma()]
#' @exportS3Method
rstandard.brma <- function(model, type = "marginal", ...) {

  # normalize type alias
  type <- match.arg(type, c("marginal", "study", "estimate", "conditional"))
  if (type == "conditional") {
    type <- "estimate"
  }

  # check model type availability
  outcome_type <- .outcome_type(model)
  is_weightfunction <- .is_weightfunction(model)
  is_multilevel <- .is_multilevel(model)

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
  if (type == "study" && !is_multilevel) {
    stop(
      "type = 'study' is only available for multilevel (3-level) models.",
      call. = FALSE
    )
  }

  # get observed effect sizes
  yi <- .outcome_data_yi(model)
  K  <- length(yi)

  # map type to predict.brma type
  pred_type <- switch(type,
    "marginal" = "terms",
    "study"    = "study",
    "estimate" = "estimate"
  )

  # get fitted values via predict.brma
  fitted_samples <- predict.brma(
    object  = model,
    newdata = NULL,
    type    = pred_type,
    quiet   = TRUE
  )

  # compute raw residuals: yi - fitted
  S <- nrow(fitted_samples)
  yi_mat <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)
  resid_samples <- yi_mat - fitted_samples

  # compute standard errors using hat matrix helper
  hat_res <- .compute_hat_matrix_samples(
    object    = model,
    type      = type,
    return_se = TRUE
  )
  se_samples <- hat_res[["se"]]

  # compute standardized residuals
  z_samples <- resid_samples / se_samples

  resid <- colMeans(resid_samples)
  se    <- colMeans(se_samples)

  # construct output data frame matching metafor::rstandard format
  out <- data.frame(
    resid = resid,
    se    = se,
    z     = colMeans(z_samples)
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
#' @param model a fitted brma object
#' @param type the type of residuals to compute. Options are:
#' \itemize{
#'   \item \code{"marginal"} (default): Residuals from marginal (fixed effects)
#'     predictions (observed - Xβ).
#'   \item \code{"study"}: Residuals from study-level predictions (observed -
#'     (Xβ + gamma)). Only available for multilevel (3-level) models.
#'   \item \code{"estimate"} (alias: \code{"conditional"}): Residuals from BLUPs
#'     (observed - theta). Note that these residuals are conceptually
#'     questionable since it would be impossible to retain the true study effect
#'     had the study been ommited.
#' }
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' This function returns a data frame with three columns matching the output
#' of \code{metafor::rstudent}:
#' \itemize{
#'   \item \code{resid}: Raw residuals (observed - fitted values)
#'   \item \code{se}: Standard errors of the residuals (back-computed from z)
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
#' # get marginal externally standardized residuals (default)
#' rstudent(fit)
#'
#' # get estimate-level (BLUP-based) residuals
#' rstudent(fit, type = "estimate")
#'
#' # check LOO diagnostics
#' loo_obj <- loo(fit)
#' plot(loo_obj) # check Pareto k values
#' }
#'
#' @seealso [rstandard.brma()], [residuals.brma()], [loo.brma()], [blup.brma()]
#' @exportS3Method
rstudent.brma <- function(model, type = "marginal", ...) {
  # normalize type alias
  type <- match.arg(type, c("marginal", "study", "estimate", "conditional"))
  if (type == "conditional") type <- "estimate"

  # check model type availability
  is_multilevel <- .is_multilevel(model)

  # study type requires multilevel model
  if (type == "study" && !is_multilevel) {
    stop("type = 'study' is only available for multilevel (3-level) models.", call. = FALSE)
  }

  # get observed effect sizes
  yi <- .outcome_data_yi(model)
  K  <- length(yi)

  # map type to predict.brma type
  pred_type <- switch(type,
    "marginal" = "terms",
    "study"    = "study",
    "estimate" = "estimate"
  )

  # get fitted values via predict.brma
  fitted_samples <- predict.brma(
    object  = model,
    newdata = NULL,
    type    = pred_type,
    quiet   = TRUE
  )

  # compute raw residuals: yi - fitted
  S <- nrow(fitted_samples)
  yi_mat <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)
  resid_samples <- yi_mat - fitted_samples

  # extract PSIS object and get normalized weights
  # (PSIS warning is done within .standardized_residuals_loopit)
  psis_weights <- loo_weights(model)

  # compute LOO-weighted residuals for each observation
  resid <- colSums(psis_weights * resid_samples)

  # get LOO-PIT z values using the specified type
  z <- .standardized_residuals_loopit(model, type = type)

  # back-compute SE from z and resid: se = resid / z
  # handle zero z values to avoid division by zero
  se <- abs(resid / z)
  se[!is.finite(se)] <- NA
  # TODO: implement a workaround by using posterior sd?

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
# 3. Compute CDF matrix F(yi | theta^(s)) (S x K matrix) for the specified type
# 4. For each observation i, compute LOO-weighted CDF:
#    u_i = sum_s w_{is} * F(yi | theta^(s))
# 5. Transform to standard normal: r_i = qnorm(u_i)
#
# Unlike traditional standardized residuals, LOO-PIT residuals:
# - Account for estimation uncertainty (integrate over posterior)
# - Account for leverage (PSIS effectively removes yi from posterior)
# - Work for all model types (selection models, GLMMs)
#
# @param object   brma object
# @param type     character; type of CDF to use:
#                 - "marginal" (default): marginal CDF yi ~ N(mu_i, tau^2 + se_i^2)
#                 - "study": conditional CDF yi ~ N(theta_i + gamma_j, se_i^2)
#                 - "estimate": conditional CDF yi ~ N(theta_i, se_i^2)
#
# @return If loo_only = TRUE, returns the psis_loo object.
#         Otherwise, returns numeric vector of LOO-PIT residuals.
#
# ---------------------------------------------------------------------------- #
.standardized_residuals_loopit <- function(object, type = "marginal") {
  # extract PSIS object and get normalized weights
  psis_weights <- loo_weights(object)

  # check Pareto k diagnostics and warn if unreliable
  check_loo(object)

  # compute CDF matrix (S x K) for the specified type
  cdf_matrix <- .cdf.brma(object, type = type)

  # compute LOO-weighted CDF for each observation
  # u_i = sum_s w_{is} * F(yi | theta^(s))
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
  mods_data <- object[["data"]][["mods"]]
  mods_formula <- attr(mods_data, "formula")

  # create the model matrix
  X <- model.matrix(mods_formula, data = mods_data)

  return(X)
}
