# ============================================================================ #
# brma.vif.R
# ============================================================================ #
#
# Variance Inflation Factor (VIF) and Generalized VIF (GVIF) for brma objects.
#
# Two complementary diagnostics for multicollinearity:
# 1. Posterior-averaged working VIF/GVIF from cov2cor(vcov)
# 2. Posterior correlation of regression coefficients (Bayesian diagnostic)
#
# VIF is computed from the correlation matrix of the posterior mean coefficient
# variance-covariance matrix: vcov = mean_s((X'W_sX)^{-1}). brma.mv known-V
# models use the full supplied sampling covariance plus any marginal
# heterogeneity or formula random-effect covariance in each posterior draw.
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
#' @param max_samples maximum posterior draws used for known-\code{V}
#' \code{brma.mv()} marginal GLS covariance computations. Defaults to
#' \code{Inf}. Finite values deterministically thin draws across posterior row
#' order.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' VIF is computed from the correlation matrix derived from the
#' coefficient variance-covariance matrix. The package computes
#' \eqn{\mathrm{E}[(X'W_sX)^{-1} \mid y]} by averaging coefficient covariance
#' across posterior heterogeneity draws. Standard models use
#' \eqn{W_s = \mathrm{diag}(w_i/(v_i + \tau_s^2))}; scale and multilevel models
#' use observation-specific \eqn{\tau_{is}} and block-structured multilevel
#' covariance where applicable. Averaging occurs before conversion to a
#' correlation matrix and computation of VIF/GVIF.
#'
#' A VIF of 1 indicates no collinearity; values above 5 or 10 are
#' commonly considered problematic.
#'
#' For multi-column terms, such as factor contrasts,
#' the Generalized VIF (GVIF) of \insertCite{fox1992generalized;textual}{RoBMA}
#' is reported. GVIF captures the joint inflation for all coefficients
#' belonging to the same term. To enable comparison across terms with different
#' degrees of freedom, \eqn{GVIF^{1/(2 \cdot df)}} is also reported; this
#' value can be compared against the usual VIF thresholds (after squaring).
#'
#' When \code{posterior_correlation = TRUE}, the function also returns the
#' posterior correlation matrix of the regression coefficients. This
#' Bayesian diagnostic complements VIF: VIF summarizes collinearity under the
#' posterior distribution of the marginal working covariance, whereas the
#' posterior correlation shows the realized joint identification of the
#' coefficients given the data and priors. Unlike a fixed-covariance design
#' VIF, posterior-averaged VIF can therefore depend on the outcome and the
#' heterogeneity prior. Informative coefficient priors can also reduce posterior
#' correlations even when VIF is high.
#'
#' Data-only objects created with \code{only_data = TRUE} are rejected because
#' VIF depends on fitted heterogeneity information, not just the parsed data.
#'
#' For \code{brma.mv()} known-\code{V} models, VIF is computed from the full
#' marginal GLS covariance. Formula random effects are marginalized through
#' BayesTools' random-effect covariance metadata rather than treated as
#' conditioned fitted effects.
#'
#' For generalized linear mixed models, this is a working effect-size-scale
#' diagnostic based on the fitted model's approximate sampling variances. It is
#' not a VIF derived directly from the binomial or Poisson information matrix.
#' When the marginal covariance is fixed, the calculation reduces to the
#' conventional fixed-covariance VIF used by \pkg{metafor}.
#'
#' @return An object of class \code{vif.brma} containing:
#' \item{vif}{A data frame with columns \code{term}, \code{df}, \code{GVIF},
#' and \code{GVIF^(1/(2*df))} (= \eqn{GVIF^{1/(2 \cdot df)}}). For single-df
#' terms, GVIF equals the standard VIF.}
#' \item{posterior_correlation}{A correlation matrix of posterior regression
#' coefficient samples when requested and available; otherwise \code{NULL}.}
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE) &&
#'     requireNamespace("metafor", quietly = TRUE)) {
#'   data(dat.bcg, package = "metadat")
#'   dat <- metafor::escalc(
#'     measure = "RR",
#'     ai      = tpos,
#'     bi      = tneg,
#'     ci      = cpos,
#'     di      = cneg,
#'     data    = dat.bcg
#'   )
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     mods    = ~ ablat + year,
#'     data    = dat,
#'     measure = "RR",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   vif(fit)
#' }
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [regplot()], [summary.brma()]
#' @exportS3Method
vif.brma <- function(object, posterior_correlation = TRUE,
                     max_samples = Inf, ...) {

  BayesTools::check_bool(posterior_correlation, "posterior_correlation")
  max_samples <- .normalize_max_samples(max_samples, "max_samples")

  # require moderators
  if (is.null(object[["priors"]]) && is.null(object[["fit"]])) {
    stop("VIF is not available for data-only objects.", call. = FALSE)
  }
  if (is.null(object[["fit"]])) {
    stop("VIF is available only for fitted objects.", call. = FALSE)
  }
  if (!.is_mods(object)) {
    stop("VIF is only meaningful for models with moderators (meta-regression).", call. = FALSE)
  }

  # compute VIF from vcov correlation matrix
  vif_df <- .compute_vif(object, max_samples = max_samples)

  # posterior correlation of regression coefficients
  post_cor <- NULL
  if (posterior_correlation) {
    post_cor <- .compute_posterior_correlation(object)
  }

  output <- list(
    vif                  = vif_df,
    posterior_correlation = post_cor
  )
  if (!is.null(attr(vif_df, "known_v_diagnostic"))) {
    output <- .known_v_attach_diagnostic_metadata(
      output,
      attr(vif_df, "known_v_diagnostic")
    )
  }

  class(output) <- "vif.brma"

  return(output)
}


# ---------------------------------------------------------------------------- #
# .compute_vif
# ---------------------------------------------------------------------------- #
#
# Compute VIF/GVIF after averaging solve(X'WX) over posterior heterogeneity
# draws so scalar, observation-specific, and block covariances share one
# package-wide estimand.
#
# @param object brma object with moderators
#
# @return data frame with columns: term, df, GVIF, GVIF^(1/(2*df))
#
# ---------------------------------------------------------------------------- #
.compute_vif <- function(object, max_samples = Inf) {

  # extract model matrix
  X      <- .get_model_matrix(object)
  assign <- attr(X, "assign")

  if (inherits(object, "brma.mv")) {
    vcov <- .vif_vcov_brma_mv(
      object      = object,
      X           = X,
      max_samples = max_samples
    )
  } else {
    vcov <- .vif_vcov_brma(object, X)
  }
  known_v_metadata <- attr(vcov, "known_v_diagnostic")

  vif_df <- .vif_table_from_vcov(
    vcov   = vcov,
    assign = assign,
    object = object,
    X      = X
  )
  if (!is.null(known_v_metadata)) {
    vif_df <- .known_v_attach_diagnostic_metadata(vif_df, known_v_metadata)
  }

  return(vif_df)
}


# ---------------------------------------------------------------------------- #
# .vif_vcov_brma
# ---------------------------------------------------------------------------- #
#
# Coefficient covariance for standard brma objects.
#
# ---------------------------------------------------------------------------- #
.vif_vcov_brma <- function(object, X) {

  # extract sampling variances
  vi      <- .outcome_data_vi(object)
  weights <- .outcome_data_weights(object)
  K       <- length(vi)

  # extract posterior heterogeneity
  is_scale      <- .is_scale(object)
  is_multilevel <- .is_multilevel(object)
  tau_K         <- if (!is_scale && !is_multilevel) 1L else K

  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(
      data = object[["data"]],
      "scale"
    ) else NULL,
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = tau_K,
    allow_missing_tau = .fixed_tau_prior_value(object[["priors"]])
  )
  tau_within_samples  <- tau_result[["tau_within"]]
  tau_between_samples <- tau_result[["tau_between"]]

  # compute posterior-averaged vcov = (X'WX)^{-1} using meta-analytic weights
  vcov <- .vif_vcov_from_tau_samples(
    X                   = X,
    vi                  = vi,
    weights             = weights,
    tau_within_samples  = tau_within_samples,
    tau_between_samples = tau_between_samples,
    cluster             = if (is_multilevel) object[["data"]][["outcome"]][["cluster"]] else NULL
  )

  return(vcov)
}


# ---------------------------------------------------------------------------- #
# .vif_vcov_brma_mv
# ---------------------------------------------------------------------------- #
#
# Coefficient covariance for brma.mv() known-V models. This is a marginal GLS
# design diagnostic: sampled and marginalized formula random effects contribute
# through ZGZ', not through fitted random-effect means.
#
# ---------------------------------------------------------------------------- #
.vif_vcov_brma_mv <- function(object, X, max_samples = Inf) {

  if (!.is_data_known_v(object[["data"]])) {
    stop("VIF for brma.mv() requires known-V covariance metadata.",
         call. = FALSE)
  }
  if (!identical(.outcome_type(object), "norm")) {
    stop("VIF for brma.mv() is available only for normal outcome models.",
         call. = FALSE)
  }

  sample_info <- .known_v_diagnostic_posterior_samples(
    object      = object,
    max_samples = max_samples,
    caller      = "known-V VIF diagnostics"
  )
  vcov <- .vif_vcov_from_covariance_chunks(
    object            = object,
    X                 = X,
    posterior_samples = sample_info[["posterior_samples"]]
  )
  vcov <- .known_v_attach_diagnostic_metadata(
    vcov,
    .known_v_diagnostic_metadata(
      sample_info = sample_info,
      chunk_info  = attr(vcov, "known_v_covariance_chunks")
    )
  )
  attr(vcov, "known_v_covariance_chunks") <- NULL

  return(vcov)
}


# ---------------------------------------------------------------------------- #
# .vif_table_from_vcov
# ---------------------------------------------------------------------------- #
#
# Compute GVIF/GVIF-adjusted table from the fixed-effect covariance.
#
# ---------------------------------------------------------------------------- #
.vif_table_from_vcov <- function(vcov, assign, object, X) {

  # identify and remove intercept (assign == 0)
  has_intercept <- 0 %in% assign
  if (has_intercept) {
    keep   <- assign != 0
    vcov   <- vcov[keep, keep, drop = FALSE]
    assign <- assign[keep]
  }

  p <- ncol(vcov)

  term_labels <- attr(X, "term_labels")
  if (is.null(term_labels)) {
    term_labels <- .fitted_formula_terms(
      object            = object,
      parameter         = "mu",
      include_intercept = FALSE,
      display           = TRUE
    )
  }

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
# .vif_vcov_from_covariance_samples
# ---------------------------------------------------------------------------- #
#
# Compute mean_s solve(X' M_s^{-1} X) from full marginal covariance samples.
#
# ---------------------------------------------------------------------------- #
.vif_vcov_from_covariance_samples <- function(X, covariance_samples,
                                              average = TRUE) {

  K <- nrow(X)
  P <- ncol(X)

  if (is.matrix(covariance_samples)) {
    covariance_samples <- array(covariance_samples, dim = c(1L, K, K))
  }
  if (length(dim(covariance_samples)) != 3L ||
      dim(covariance_samples)[2L] != K ||
      dim(covariance_samples)[3L] != K) {
    stop("VIF covariance samples must have dimensions draw x row x row.",
         call. = FALSE)
  }

  S        <- dim(covariance_samples)[1L]
  vcov_sum <- matrix(0, nrow = P, ncol = P)

  for (s in seq_len(S)) {
    covariance <- covariance_samples[s, , , drop = FALSE]
    covariance <- matrix(covariance, nrow = K, ncol = K)
    covariance <- (covariance + t(covariance)) / 2

    chol_covariance <- tryCatch(
      chol(covariance),
      error = function(e) NULL
    )
    if (is.null(chol_covariance)) {
      stop("VIF marginal covariance is not positive definite.",
           call. = FALSE)
    }

    W <- chol2inv(chol_covariance)
    vcov_sum <- vcov_sum + .hat_solve_crossprod(crossprod(X, W %*% X))
  }

  vcov <- if (average) vcov_sum / S else vcov_sum
  if (!is.null(colnames(X))) {
    dimnames(vcov) <- list(colnames(X), colnames(X))
  }

  return(vcov)
}


.vif_vcov_from_covariance_chunks <- function(object, X, posterior_samples) {

  P        <- ncol(X)
  vcov_sum <- matrix(0, nrow = P, ncol = P)

  chunk_info <- .known_v_apply_marginal_covariance_chunks(
    object            = object,
    posterior_samples = posterior_samples,
    FUN               = function(covariance_samples, rows) {
      vcov_sum <<- vcov_sum + .vif_vcov_from_covariance_samples(
        X                  = X,
        covariance_samples = covariance_samples,
        average            = FALSE
      )
    }
  )

  vcov <- vcov_sum / nrow(posterior_samples)
  if (!is.null(colnames(X))) {
    dimnames(vcov) <- list(colnames(X), colnames(X))
  }
  attr(vcov, "known_v_covariance_chunks") <- chunk_info

  return(vcov)
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
# @param tau_within_samples S x 1 or S x K matrix of estimate-level SDs.
# @param tau_between_samples optional S x 1 or S x K matrix of cluster-level SDs.
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
    tau_between_samples <- matrix(
      0,
      nrow = nrow(tau_within_samples),
      ncol = ncol(tau_within_samples)
    )
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
    tau_within_s  <- tau_within_samples[s, ]
    tau_between_s <- tau_between_samples[s, ]
    if (!is.null(cluster) && length(tau_between_s) == 1L) {
      tau_between_s <- rep(tau_between_s, K)
    }
    diagonal_s    <- (vi + tau_within_s^2) / weights
    WX <- .hat_apply_precision(
      x             = X,
      diagonal      = diagonal_s,
      rank_one      = if (!is.null(cluster)) tau_between_s else NULL,
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
# Normalize heterogeneity input to a compact S x 1 or row-specific S x K matrix.
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

  if (!ncol(x) %in% c(1L, K)) {
    stop("'", name, "' must have one or K columns.", call. = FALSE)
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
  samples_mat <- .diagnostic_fixed_location_coefficient_samples(
    samples_mat,
    require_mu_columns = TRUE
  )

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
