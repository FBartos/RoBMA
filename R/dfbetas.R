# ============================================================================ #
# brma.dfbetas.R
# ============================================================================ #
#
# This file contains the dfbetas method for brma model fits. DFBETAS values
# measure the influence of each observation on the estimated model coefficients.
#
# Instead of computationally expensive refitting (leave-one-out), we use
# Pareto Smoothed Importance Sampling (PSIS) weights from the LOO-CV
# approximation to estimate the leave-one-out coefficients.
#
# ============================================================================ #


#' @title DFBETAS for brma Objects
#'
#' @description Computes DFBETAS (Difference in FITS, standardized) for a
#' fitted brma object. DFBETAS measures the influence of each observation on
#' the estimated model coefficients. Positive values indicate that deleting the
#' observation yields a smaller estimate, negative values indicate that deleting
#' the observation yields a larger estimate.
#'
#' @param model a fitted brma object.
#' @param type type of parameters to be summarized. Defaults to \code{"mods"}
#' (for the effect size and meta-regression coefficients). The other option is
#' \code{"scale"} (for the heterogeneity and scale-regression coefficients).
#' @param standardized_coefficients whether to show standardized meta-regression coefficients.
#' Defaults to \code{FALSE}. When set to \code{TRUE}, standardized meta-regression
#' coefficients are returned for the intercept and continuous predictors. These coefficients
#' correspond to the standardized scale on which prior distributions are specified by default
#' (i.e., `standardize_continuous_predictors = TRUE`).
#' @param transform_factors whether to transform factors to their original names.
#' Defaults to \code{TRUE}.
#' @param ... additional arguments (currently ignored).
#'
#' @details
#' This function computes DFBETAS values using the Leave-One-Out (LOO)
#' approximation based on Pareto Smoothed Importance Sampling (PSIS) weights.
#' Ideally, DFBETAS is defined as:
#' \deqn{DFBETAS_{ij} = \frac{\hat{\beta}_j - \hat{\beta}_{j(-i)}}{SE(\hat{\beta}_{j(-i)})}}
#' where \eqn{\hat{\beta}_j} is the estimate of the \eqn{j}-th coefficient using
#' the full data, \eqn{\hat{\beta}_{j(-i)}} is the estimate when observation \eqn{i}
#' is omitted, and \eqn{SE(\hat{\beta}_{j(-i)})} is the standard error of the
#' coefficient when observation \eqn{i} is omitted.
#'
#' In the Bayesian context using LOO approximation:
#' \itemize{
#'   \item \eqn{\hat{\beta}_{j(-i)}} is estimated as the importance sampling
#'     weighted mean of the posterior samples, using PSIS weights \eqn{w_{is}}.
#'   \item \eqn{SE(\hat{\beta}_{j(-i)})} is estimated as the importance sampling
#'     weighted standard deviation of the posterior samples.
#' }
#'
#' This approximation allows computing influence statistics without refitting
#' the model \eqn{K} times, making it computationally efficient.
#'
#' Note: This function requires that LOO-CV has been computed for the model
#' using \code{\link{add_loo}}.
#'
#' @return A data frame with \eqn{K} rows (observations) and \eqn{P} columns
#' (coefficients), containing the DFBETAS values. Row names correspond to
#' study labels (if available) or indices.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # compute LOO (required)
#' fit <- add_loo(fit)
#'
#' # compute DFBETAS
#' inf <- dfbetas(fit)
#' plot(inf[, 1], type = "h")
#' }
#'
#' @seealso \code{\link{add_loo}}, \code{\link{loo_weights.brma}}
#' @exportS3Method
dfbetas.brma <- function(model, type = "mods", standardized_coefficients = FALSE, transform_factors = TRUE, return_loo_estimates = FALSE, ...) {

  BayesTools::check_char(type, "type", allow_values = c("mods", "scale"))

  # get PSIS weights (S x K matrix)
  weights <- loo_weights(model)

  # determine whether to extract formula (for meta-regression) or parameter (for intercept-only)
  is_mods  <- .is_mods(model)
  is_scale <- .is_scale(model)

  if (type == "mods") {
    if (is_mods) {
      # meta-regression: means mu is a formula
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                = model[["fit"]],
        keep_formulas      = "mu",
        remove_diagnostics = TRUE,
        transform_factors  = transform_factors,
        transform_scaled   = !standardized_coefficients,
        return_samples     = TRUE
      )
    } else {
      # random/fixed effects: mu is a parameter (intercept)
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                = model[["fit"]],
        keep_parameters    = "mu",
        remove_diagnostics = TRUE,
        transform_factors  = transform_factors,
        transform_scaled   = !standardized_coefficients,
        return_samples     = TRUE
      )
    }
  } else if (type == "scale") {
    if (is_scale) {
      # scale-regression: tau is modeled via log_tau formula
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                = model[["fit"]],
        keep_formulas      = "log_tau",
        remove_diagnostics = TRUE,
        transform_factors  = transform_factors,
        return_samples     = TRUE
      )
    } else {
      # random/fixed effects: tau is a parameter (intercept)
      samples_table <- BayesTools::JAGS_estimates_table(
        fit                = model[["fit"]],
        keep_parameters    = "tau",
        remove_diagnostics = TRUE,
        transform_factors  = transform_factors,
        return_samples     = TRUE
      )
    }
  }


  # samples_table is a flattened data frame (S rows x P columns)
  # ensure it's a matrix for computation
  samples_mat <- as.matrix(samples_table)

  # dimensions
  S <- nrow(samples_mat) # number of samples
  P <- ncol(samples_mat) # number of coefficients
  K <- ncol(weights)     # number of observations

  # 1. Compute LOO-weighted means for each observation i and parameter j
  # beta_loo[i, j] = sum_s (w_{is} * beta_{js})
  # weights is S x K -> t(weights) is K x S
  # samples_mat is S x P
  # K x S %*% S x P -> K x P
  beta_loo <- crossprod(weights, samples_mat)

  # 2. Compute Robust Weighted SD (SE LOO)
  # Uses centered moment calculation: sum(w * (x - mu)^2)
  # This avoids catastrophic cancellation issues with E[x^2] - E[x]^2
  se_loo <- matrix(NA, nrow = K, ncol = P)

  for (j in seq_len(P)) {
    # Get samples for parameter j (vector length S)
    beta_s <- samples_mat[, j]

    # Get LOO means for parameter j (vector length K)
    mu_k <- beta_loo[, j]

    # Efficient Centering:
    # Use outer() to create an (S x K) matrix of differences
    # element [s, k] = beta_s[s] - mu_k[k]
    # This vectorizes the subtraction of every sample from every LOO mean
    diff_mat <- outer(beta_s, mu_k, "-")

    # Square differences
    diff_sq <- diff_mat^2

    # Weighted Sum of Squares (Variance)
    # colSums of (S x K) * (S x K) -> Vector length K
    # This is the "Population Variance" of the LOO posterior distribution
    var_k <- colSums(weights * diff_sq)

    # Store SE
    se_loo[, j] <- sqrt(var_k)
  }

  # 3. Compute full posterior means (unweighted / equal weights)
  # beta_full[j] = mean(beta_{js})
  beta_full <- colMeans(samples_mat)

  # expand beta_full to K x P matrix for vectorized subtraction
  beta_full_mat <- matrix(beta_full, nrow = K, ncol = P, byrow = TRUE)

  # 4. Return LOO estimates if requested
  if (return_loo_estimates) {
    beta_loo_df <- as.data.frame(beta_loo)
    colnames(beta_loo_df) <- colnames(samples_mat)
    return(beta_loo_df)
  }

  # 4. Compute DFBETAS
  # (beta_full - beta_loo) / se_loo
  dfbetas_val <- (beta_full_mat - beta_loo) / se_loo

  # convert to data frame
  dfbetas_df <- as.data.frame(dfbetas_val)
  colnames(dfbetas_df) <- colnames(samples_mat)

  return(dfbetas_df)
}
