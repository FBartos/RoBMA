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
#' @description Computes DFBETAS (Difference in BETAS, standardized) for a
#' fitted brma object. DFBETAS measures the influence of each observation on
#' the estimated model coefficients. Positive values indicate that deleting the
#' observation yields a smaller estimate, negative values indicate that deleting
#' the observation yields a larger estimate.
#'
#' @param model a fitted brma object.
#' @param type type of parameters to be summarized. Defaults to \code{"mods"}
#' (for the effect size and meta-regression coefficients). The other options are
#' \code{"scale"} (for the heterogeneity and scale-regression coefficients) and
#' \code{"bias"} (for omega, PET, and PEESE publication-bias parameters).
#' @param standardized_coefficients whether to show standardized meta-regression coefficients.
#' Defaults to \code{FALSE}. When set to \code{TRUE}, standardized meta-regression
#' coefficients are returned for the intercept and continuous predictors. These coefficients
#' correspond to the standardized scale on which prior distributions are specified by default
#' (i.e., `standardize_continuous_predictors = TRUE`).
#' @param transform_factors whether to transform factors to their original names.
#' Defaults to \code{TRUE}.
#' @param return_loo_estimates whether to return the leave-one-out coefficient
#' estimates used to compute DFBETAS instead of standardized DFBETAS values.
#' Defaults to \code{FALSE}.
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
#' For \code{brma.mv()} objects, DFBETAS use estimate-unit PSIS weights. With
#' correlated known-\code{V} data, this is influence under conditional estimate
#' deletion, not independent-study deletion.
#' For \code{type = "bias"}, fixed identification parameters (e.g., the reference
#' \eqn{\omega = 1} interval) can have zero LOO posterior standard deviation.
#' These parameters are retained in the output, but their DFBETAS values are
#' reported as \code{NaN} because the standardized diagnostic is undefined.
#'
#' Note: This function requires that LOO-CV has been computed for the model
#' using \code{\link{add_loo}}.
#'
#' @return If `return_loo_estimates = FALSE`, a data frame with \eqn{K} rows
#' (observations) and \eqn{P} columns (parameters), containing DFBETAS values.
#' If `return_loo_estimates = TRUE`, returns the corresponding leave-one-out
#' coefficient estimates. Row names correspond to study labels (if available)
#' or indices.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit <- add_loo(fit)
#'
#'   inf <- dfbetas(fit)
#'   plot(inf[, 1], type = "h")
#' }
#' }
#'
#' @seealso \code{\link{add_loo}}, \code{\link{loo_weights.brma}}
#' @importFrom stats dfbetas
#' @exportS3Method
dfbetas.brma <- function(model, type = "mods", standardized_coefficients = FALSE,
                         transform_factors = TRUE,
                         return_loo_estimates = FALSE, ...) {

  dots <- list(...)
  .weights <- dots[[".weights"]]
  BayesTools::check_char(type, "type", allow_values = c("mods", "scale", "bias"))

  # get PSIS weights (S x K matrix)
  weights <- .diagnostic_psis_weights(model, .weights)

  # determine whether to extract formula (for meta-regression) or parameter (for intercept-only)
  is_scale <- .is_scale(model)
  is_bias  <- .is_bias(model)

  if (type == "mods") {
    samples_table <- .diagnostic_location_parameter_samples(
      model                     = model,
      standardized_coefficients = standardized_coefficients,
      transform_factors         = transform_factors
    )
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
  } else if (type == "bias") {
    if (!is_bias) {
      stop("type = 'bias' is only available for models with publication bias adjustment.", call. = FALSE)
    }

    samples_table <- .dfbetas_bias_samples(model)
  }


  # samples_table is a flattened data frame (S rows x P columns)
  # ensure it's a matrix for computation
  samples_mat <- as.matrix(samples_table)

  if (ncol(samples_mat) == 0) {
    stop("No parameters available for DFBETAS with the requested type.", call. = FALSE)
  }

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
    beta_loo_df <- .diagnostic_set_rownames(beta_loo_df, model)
    return(beta_loo_df)
  }

  # 4. Compute DFBETAS
  # (beta_full - beta_loo) / se_loo
  dfbetas_val <- (beta_full_mat - beta_loo) / se_loo
  undefined   <- apply(se_loo <= sqrt(.Machine$double.eps), 2, any)
  note        <- NULL
  if (any(undefined)) {
    dfbetas_val[, undefined] <- NaN
    note <- .diagnostic_zero_variance_note(
      diagnostic = "DFBETAS",
      parameters = colnames(samples_mat)[undefined]
    )
  }
  # convert to data frame
  dfbetas_df <- as.data.frame(dfbetas_val)
  colnames(dfbetas_df) <- colnames(samples_mat)
  dfbetas_df <- .diagnostic_set_rownames(dfbetas_df, model)

  if (!is.null(note)) {
    attr(dfbetas_df, "note") <- note
  }
  class(dfbetas_df) <- c("dfbetas.brma", class(dfbetas_df))

  return(dfbetas_df)
}


#' @exportS3Method
print.dfbetas.brma <- function(x, ...) {

  note <- attr(x, "note")
  class(x) <- "data.frame"
  print(x, ...)
  .print_diagnostic_note(note)

  return(invisible(x))
}


# ---------------------------------------------------------------------------- #
# .dfbetas_bias_samples
# ---------------------------------------------------------------------------- #
#
# Extract publication-bias posterior samples without using JAGS_estimates_table.
# BayesTools summary tables always run formula/factor post-processing, which can
# fail for RoBMA mixture-prior objects when we only need raw omega/PET/PEESE draws.
#
# ---------------------------------------------------------------------------- #
.dfbetas_bias_samples <- function(model) {

  posterior_samples <- .get_posterior_samples(model[["fit"]])
  bias_samples      <- list()
  omega_cols        <- grep("^omega(\\[|$)", colnames(posterior_samples), value = TRUE)

  if (length(omega_cols) > 0L) {
    bias_samples[["omega"]] <- .extract_omega_samples(posterior_samples)
  }

  for (par in c("PET", "PEESE")) {
    if (par %in% colnames(posterior_samples)) {
      bias_samples[[par]] <- as.matrix(posterior_samples[, par, drop = FALSE])
    }
  }

  if (length(bias_samples) == 0L) {
    return(matrix(nrow = nrow(posterior_samples), ncol = 0))
  }

  return(do.call(cbind, bias_samples))
}
