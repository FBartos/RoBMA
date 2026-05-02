# ============================================================================ #
# brma.vif.R
# ============================================================================ #
#
# Variance Inflation Factor (VIF) and Generalized VIF (GVIF) for brma objects.
#
# Two complementary diagnostics for multicollinearity:
# 1. Classical VIF/GVIF from cov2cor(vcov) (matches metafor)
# 2. Posterior correlation of regression coefficients (Bayesian diagnostic)
#
# VIF is computed from the correlation matrix of the coefficient
# variance-covariance matrix: vcov = (X'WX)^{-1}. For simple models, W uses
# diag(weight_i/(vi + tau^2)); scale and multilevel models average this
# coefficient covariance over posterior heterogeneity draws.
#
# GVIF formula (Fox & Monette, 1992, for multi-df terms):
#   GVIF = det(R_11) * det(R_22) / det(R)
#
# GSIF (comparable across terms with different df):
#   GSIF = GVIF^(1/(2*m))
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# vif generic
# ---------------------------------------------------------------------------- #

#' @title Variance Inflation Factors
#'
#' @description Computes variance inflation factors (VIF) for a fitted model
#' with moderators, to assess multicollinearity among predictors.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value containing VIF diagnostics.
#'
#' @references
#' \insertCite{fox1992generalized}{RoBMA}
#'
#' @seealso [regplot()], [summary.brma()]
#' @export
vif <- function(object, ...) {

  UseMethod("vif")
}


# ---------------------------------------------------------------------------- #
# vif.brma
# ---------------------------------------------------------------------------- #

#' @title Variance Inflation Factors for brma Objects
#'
#' @description Computes variance inflation factors (VIF) and generalized
#' VIF (GVIF) for a fitted brma meta-regression model. Also optionally
#' returns the posterior correlation matrix of regression coefficients.
#'
#' @param object a fitted brma object with moderators
#' @param posterior_correlation logical; whether to also compute and return
#' the posterior correlation matrix of regression coefficients. Defaults
#' to \code{TRUE}.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' VIF is computed from the correlation matrix derived from the
#' coefficient variance-covariance matrix \eqn{(X'WX)^{-1}}. For standard
#' meta-regression models, \eqn{W = \mathrm{diag}(w_i/(v_i + \hat\tau^2))}
#' with \eqn{\hat\tau} equal to the posterior mean heterogeneity. For scale
#' and multilevel models, the coefficient covariance is averaged across
#' posterior heterogeneity draws, using observation-specific \eqn{\tau_i}
#' and block-structured multilevel covariance where applicable.
#'
#' A VIF of 1 indicates no collinearity; values above 5 or 10 are
#' commonly considered problematic.
#'
#' For categorical moderators (factors) that produce multiple dummy variables,
#' the Generalized VIF (GVIF) of \insertCite{fox1992generalized;textual}{RoBMA}
#' is reported. GVIF captures the joint inflation for all coefficients
#' belonging to the same term. To enable comparison across terms with different
#' degrees of freedom, \eqn{GVIF^{1/(2 \cdot df)}} is also reported; this
#' value can be compared against the usual VIF thresholds (after squaring).
#'
#' When \code{posterior_correlation = TRUE}, the function also returns the
#' posterior correlation matrix of the regression coefficients. This
#' Bayesian diagnostic complements VIF: while VIF diagnoses the
#' \emph{potential} for collinearity problems (a data property), the
#' posterior correlation shows the \emph{realized} identification
#' given the data and priors. Informative priors can mitigate
#' collinearity, reducing posterior correlations even when VIF is high.
#'
#' @return An object of class \code{vif.brma} containing:
#' \item{vif}{A data frame with columns \code{term}, \code{df}, \code{GVIF},
#' and \code{GVIF_df} (= \eqn{GVIF^{1/(2 \cdot df)}}). For single-df terms,
#' GVIF equals the standard VIF.}
#' \item{posterior_correlation}{(If requested) A correlation matrix of
#' posterior regression coefficient samples.}
#'
#' @examples \dontrun{
#' # fit a meta-regression
#' fit <- brma(yi ~ ablat + year, sei = sei, data = dat)
#'
#' # compute VIF
#' vif(fit)
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [regplot()], [summary.brma()]
#' @exportS3Method
vif.brma <- function(object, posterior_correlation = TRUE, ...) {

  BayesTools::check_bool(posterior_correlation, "posterior_correlation")

  # require moderators
  if (!.is_mods(object)) {
    stop("VIF is only meaningful for models with moderators (meta-regression).", call. = FALSE)
  }

  # compute VIF from vcov correlation matrix
  vif_df <- .compute_vif(object)

  # posterior correlation of regression coefficients
  post_cor <- NULL
  if (posterior_correlation) {
    post_cor <- .compute_posterior_correlation(object)
  }

  output <- list(
    vif                  = vif_df,
    posterior_correlation = post_cor
  )

  class(output) <- "vif.brma"

  return(output)
}


# ---------------------------------------------------------------------------- #
# .compute_vif
# ---------------------------------------------------------------------------- #
#
# Compute VIF/GVIF from cov2cor(solve(X'WX)). Simple models use the posterior
# mean tau plug-in. Scale and multilevel models average solve(X'WX) over
# posterior heterogeneity draws so observation-specific and block covariances
# are preserved.
#
# @param object brma object with moderators
#
# @return data frame with columns: term, df, GVIF, GVIF^(1/(2*df))
#
# ---------------------------------------------------------------------------- #
.compute_vif <- function(object) {

  # extract model matrix
  X      <- .get_model_matrix(object)
  assign <- attr(X, "assign")

  # extract sampling variances
  vi      <- .outcome_data_vi(object)
  weights <- .outcome_data_weights(object)
  K       <- length(vi)

  # extract posterior heterogeneity
  # for simple models: use summary["tau", "Mean"] to preserve metafor-style VIF
  # for scale/multilevel models: use full posterior tau matrices
  is_scale      <- .is_scale(object)
  is_multilevel <- .is_multilevel(object)

  if (!is_scale && !is_multilevel) {
    tau_within_samples  <- matrix(object[["summary"]]["tau", "Mean"], nrow = 1, ncol = K)
    tau_between_samples <- matrix(0, nrow = 1, ncol = K)
  } else {
    tau_result <- .evaluate.brma.tau(
      fit           = object[["fit"]],
      scale_data    = object[["data"]][["scale"]],
      scale_formula = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
      scale_priors  = object[["priors"]][["scale"]],
      is_scale      = is_scale,
      is_multilevel = is_multilevel,
      K             = K
    )
    tau_within_samples  <- tau_result[["tau_within"]]
    tau_between_samples <- tau_result[["tau_between"]]
  }

  # compute posterior-averaged vcov = (X'WX)^{-1} using meta-analytic weights
  vcov <- .vif_vcov_from_tau_samples(
    X                   = X,
    vi                  = vi,
    weights             = weights,
    tau_within_samples  = tau_within_samples,
    tau_between_samples = tau_between_samples,
    cluster             = if (is_multilevel) object[["data"]][["outcome"]][["cluster"]] else NULL
  )

  # identify and remove intercept (assign == 0)
  has_intercept <- 0 %in% assign
  if (has_intercept) {
    keep   <- assign != 0
    vcov   <- vcov[keep, keep, drop = FALSE]
    assign <- assign[keep]
  }

  p <- ncol(vcov)

  # term labels from formula
  mods_data    <- object[["data"]][["mods"]]
  mods_formula <- attr(mods_data, "formula")
  term_labels  <- attr(terms(mods_formula), "term.labels")

  # compute GVIF per term from vcov correlation matrix
  terms_unique <- sort(unique(assign))
  n_terms      <- length(terms_unique)

  gvif_values    <- numeric(n_terms)
  gvif_df_values <- numeric(n_terms)
  df_values      <- integer(n_terms)
  term_names     <- character(n_terms)

  if (p > 1) {
    R    <- stats::cov2cor(vcov)
    detR <- det(R)
  }

  for (k in seq_len(n_terms)) {

    term_index <- terms_unique[k]
    cols       <- which(assign == term_index)
    m          <- length(cols)

    df_values[k]  <- m
    term_names[k] <- if (term_index <= length(term_labels)) term_labels[term_index] else colnames(vcov)[cols[1]]

    if (p <= 1 || m == p) {
      # single predictor or all columns belong to one term
      gvif_values[k]    <- 1
      gvif_df_values[k] <- 1
    } else {
      # Fox & Monette (1992) formula
      gvif <- det(R[cols, cols, drop = FALSE]) *
              det(R[-cols, -cols, drop = FALSE]) / detR

      gvif_values[k]    <- gvif
      gvif_df_values[k] <- gvif^(1 / (2 * m))
    }
  }

  vif_df <- data.frame(
    term               = term_names,
    df                 = df_values,
    GVIF               = gvif_values,
    "GVIF^(1/(2*df))"  = gvif_df_values,
    stringsAsFactors   = FALSE,
    check.names        = FALSE
  )

  return(vif_df)
}


# ---------------------------------------------------------------------------- #
# .vif_vcov_from_tau_samples
# ---------------------------------------------------------------------------- #
#
# Compute the coefficient covariance used by VIF from one or more posterior
# heterogeneity draws.
#
# @param X design matrix.
# @param vi sampling variances.
# @param weights likelihood weights.
# @param tau_within_samples S x K matrix of estimate-level heterogeneity SDs.
# @param tau_between_samples optional S x K matrix of cluster-level SDs.
# @param cluster optional cluster identifiers for multilevel block covariance.
#
# @return posterior mean of solve(X'WX).
#
# ---------------------------------------------------------------------------- #
.vif_vcov_from_tau_samples <- function(X, vi, weights, tau_within_samples,
                                       tau_between_samples = NULL,
                                       cluster = NULL) {

  K <- length(vi)
  P <- ncol(X)

  tau_within_samples <- .vif_tau_matrix(tau_within_samples, K, "tau_within_samples")

  if (is.null(tau_between_samples)) {
    tau_between_samples <- matrix(0, nrow = nrow(tau_within_samples), ncol = K)
  } else {
    tau_between_samples <- .vif_tau_matrix(tau_between_samples, K, "tau_between_samples")
  }

  if (nrow(tau_between_samples) != nrow(tau_within_samples)) {
    stop("'tau_within_samples' and 'tau_between_samples' must have the same number of rows.",
         call. = FALSE)
  }

  block_indices <- list()
  if (!is.null(cluster)) {
    block_indices <- .get_multilevel_block_indices(cluster)
  }

  S        <- nrow(tau_within_samples)
  vcov_sum <- matrix(0, nrow = P, ncol = P)

  for (s in seq_len(S)) {
    diagonal_s <- (vi + tau_within_samples[s, ]^2) / weights
    WX <- .hat_apply_precision(
      x             = X,
      diagonal      = diagonal_s,
      rank_one      = if (!is.null(cluster)) tau_between_samples[s, ] else NULL,
      block_indices = block_indices
    )

    vcov_sum <- vcov_sum + .hat_solve_crossprod(crossprod(X, WX))
  }

  vcov <- vcov_sum / S
  if (!is.null(colnames(X))) {
    dimnames(vcov) <- list(colnames(X), colnames(X))
  }

  return(vcov)
}


# ---------------------------------------------------------------------------- #
# .vif_tau_matrix
# ---------------------------------------------------------------------------- #
#
# Normalize heterogeneity input to an S x K matrix.
#
# ---------------------------------------------------------------------------- #
.vif_tau_matrix <- function(x, K, name) {

  if (is.null(dim(x))) {
    if (length(x) == 1L) {
      return(matrix(x, nrow = 1, ncol = K))
    }
    if (length(x) == K) {
      return(matrix(x, nrow = 1, ncol = K))
    }
    stop("'", name, "' must have length 1, length K, or be an S x K matrix.",
         call. = FALSE)
  }

  x <- as.matrix(x)

  if (ncol(x) == 1L && K > 1L) {
    x <- matrix(x[, 1], nrow = nrow(x), ncol = K)
  }

  if (ncol(x) != K) {
    stop("'", name, "' must have K columns.", call. = FALSE)
  }

  return(x)
}


# ---------------------------------------------------------------------------- #
# .compute_posterior_correlation
# ---------------------------------------------------------------------------- #
#
# Extract posterior samples of regression coefficients and compute their
# correlation matrix. Caller must ensure the model has moderators.
#
# @param object brma object with moderators
#
# @return correlation matrix, or NULL if fewer than 2 coefficients
#
# ---------------------------------------------------------------------------- #
.compute_posterior_correlation <- function(object) {

  # extract coefficient samples via BayesTools
  # keep_formulas = "mu" extracts all mu-formula coefficients
  # return_samples = TRUE returns S x P matrix instead of summary
  samples_mat <- as.matrix(BayesTools::JAGS_estimates_table(
    fit                = object[["fit"]],
    keep_formulas      = "mu",
    remove_diagnostics = TRUE,
    transform_factors  = TRUE,
    transform_scaled   = TRUE,
    return_samples     = TRUE
  ))

  # need at least 2 coefficients for correlation
  if (ncol(samples_mat) < 2) {
    return(NULL)
  }

  return(stats::cor(samples_mat))
}


# ---------------------------------------------------------------------------- #
# print.vif.brma
# ---------------------------------------------------------------------------- #

#' @title Print VIF Results
#'
#' @description Prints variance inflation factors and optional posterior
#' correlation matrix.
#'
#' @param x a \code{vif.brma} object
#' @param digits number of decimal places. Defaults to 3.
#' @param ... additional arguments (currently ignored)
#'
#' @return Returns \code{x} invisibly.
#'
#' @exportS3Method
print.vif.brma <- function(x, digits = 3, ...) {

  vif_df <- x$vif

  cat("\nVariance Inflation Factors:\n")

  if (all(vif_df$df == 1)) {
    # simple VIF display (no GVIF columns needed)
    out        <- vif_df$GVIF
    names(out) <- vif_df$term
    print(round(out, digits))
  } else {
    # full GVIF display with term as rownames
    print_df           <- vif_df[, -1, drop = FALSE]
    rownames(print_df) <- vif_df$term
    print(round(print_df, digits))
  }

  if (!is.null(x$posterior_correlation)) {
    cat("\nPosterior Correlation of Coefficients:\n")
    print(round(x$posterior_correlation, digits))
  }

  cat("\n")

  return(invisible(x))
}
