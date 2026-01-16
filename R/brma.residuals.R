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
#'   \item \code{"LOO-PIT"}: Leave-one-out probability integral transform (PIT)
#'     residuals computed via Pareto smoothed importance sampling. The LOO-CDF
#'     value for each observation is computed and transformed to standard normal
#'     quantiles via \eqn{\Phi^{-1}(u_i)}. Under a correctly specified model,
#'     these residuals should follow a standard normal distribution. This is the
#'     recommended standardized residual for Bayesian models as it properly
#'     accounts for estimation uncertainty and leverage. Available for all model
#'     types. Note: This is computationally more expensive than other methods.
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
#' plot(loo_result)  # check Pareto k values
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [predict.brma()], [blup.brma()], [pooled_effect()], [rstandard.brma()]
#' @exportS3Method
residuals.brma <- function(object, type = "outcome", bias_adjusted = FALSE, ...) {

  # input validation
  BayesTools::check_char(type, "type", allow_values = c("outcome", "pearson", "rstandard", "LOO-PIT"))
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")

  # extract model characteristics for error checking
  outcome_type      <- .outcome_type(object)
  is_weightfunction <- .is_weightfunction(object)

  # check for unsupported type/model combinations
  .check_residual_type_availability(type, outcome_type, is_weightfunction)

  # dispatch to LOO-PIT method (separate path due to different computation)
  if (type == "LOO-PIT") {

    out <- .standardize_residuals_loopit(object)

  } else {

    # get fitted values (mu predictions for each observation)
    # type = "terms" gives fixed effects predictions estimate-level random effects
    fitted_samples <- predict.brma(
      object        = object,
      newdata       = NULL,
      type          = "terms",
      bias_adjusted = bias_adjusted,
      as_samples    = TRUE,
      quiet         = TRUE
    )

    # get observed effect sizes
    yi <- .outcome_data_yi(object)
    K  <- length(yi)

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

    }

    # compute posterior mean residuals
    out <- colMeans(resid_samples)

  }

  # clean names
  names(out) <- NULL
  return(out)
}


#' @title Standardized Residuals for brma Objects
#'
#' @description Computes standardized residuals from a fitted brma object using
#' LOO-PIT (Leave-One-Out Probability Integral Transform). This is a
#' convenience wrapper for \code{residuals(object, type = "LOO-PIT", ...)}.
#'
#' @param model a fitted brma object
#' @param ... additional arguments passed to \code{\link{residuals.brma}}
#'
#' @details
#' LOO-PIT residuals are the Bayesian equivalent of studentized deleted
#' residuals. They are computed via leave-one-out probability integral
#' transformation using Pareto smoothed importance sampling. Under a correctly
#' specified model, LOO-PIT residuals should follow a standard normal
#' distribution.
#'
#' Unlike traditional standardized residuals (which use the hat matrix),
#' LOO-PIT residuals properly account for estimation uncertainty and leverage
#' without requiring explicit hat matrix computation. This makes them suitable
#' for all model types including selection models and GLMMs.
#'
#' Note: LOO-PIT residuals are computationally more expensive than traditional
#' standardized residuals as they require importance sampling.
#'
#' @return A numeric vector of LOO-PIT standardized residual means.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get standardized residuals (LOO-PIT)
#' rstandard(fit)
#'
#' # equivalent to:
#' residuals(fit, type = "LOO-PIT")
#'
#' # check LOO diagnostics
#' loo_obj <- residuals(fit, type = "LOO-PIT", loo_only = TRUE)
#' plot(loo_obj)  # check Pareto k values
#' }
#'
#' @seealso [residuals.brma()], [loo.brma()], [predict.brma()]
#' @exportS3Method
rstandard.brma <- function(model, type = "LOO-PIT", ...) {
  residuals.brma(object = model, type = type, ...)
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
# .standardize_residuals_loopit
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
# 3. Compute full CDF matrix F(yi | theta^(s)) (S x K matrix)
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
# @param loo_only if TRUE, return only the LOO object (for diagnostics)
# @param probs    quantiles for summary (not used if loo_only = TRUE)
#
# @return If loo_only = TRUE, returns the psis_loo object.
#         Otherwise, returns list with:
#         - samples: S x K matrix of LOO-PIT residual samples
#         - loo: the psis_loo object (for diagnostics)
#
# ---------------------------------------------------------------------------- #
.standardize_residuals_loopit <- function(object, loo_only = FALSE, probs = c(.025, .975)) {

  # check that loo package is available
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Package 'loo' is required for LOO-PIT residuals. ",
         "Please install it with: install.packages('loo')", call. = FALSE)
  }

  # compute LOO via PSIS
  loo_result <- loo.brma(object, save_psis = TRUE)

  # early exit if only LOO object requested
  if (loo_only) {
    return(loo_result)
  }

  # extract PSIS object and get normalized weights
  psis_object <- loo_result[["psis_object"]]
  # weights() returns S x K matrix of normalized (not log) weights
  psis_weights <- loo::weights.importance_sampling(psis_object, log = FALSE, normalize = TRUE)

  # check Pareto k diagnostics and warn if unreliable
  pareto_k <- loo_result[["diagnostics"]][["pareto_k"]]
  bad_k <- which(pareto_k > 0.7)
  if (length(bad_k) > 0) {
    warning(
      "Some Pareto k values are high (> 0.7), indicating potentially unreliable ",
      "LOO-PIT residuals for observations: ", paste(bad_k, collapse = ", "), ". ",
      "Consider using loo_only = TRUE to inspect diagnostics.",
      call. = FALSE
    )
  }

  # compute full CDF matrix (S x K)
  cdf_matrix <- .cdf.brma(object)

  # get dimensions
  S <- nrow(cdf_matrix)
  K <- ncol(cdf_matrix)

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

    # hat matrix H = X (X' W X)^{-1} X' W
    H <- .hat_matrix(X, W)

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




#' Check if a residual type is supported for the given model
#'
#' Helper function to validate that the requested residual type is compatible
#' with the model's outcome type and presence of selection (weightfunction).
#'
#' @param type Residual type requested
#' @param outcome_type Model outcome type
#' @param is_weightfunction Whether model includes weightfunction
#'
#' @return NULL (throws error if invalid combination)
#' @keywords internal
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
  mods_data    <- object[["data"]][["mods"]]
  mods_formula <- attr(mods_data, "formula")

  # create the model matrix
  X <- model.matrix(mods_formula, data = mods_data)

  return(X)
}
