# ============================================================================ #
# brma.vif.R
# ============================================================================ #
#
# Variance Inflation Factor (VIF) and Generalized VIF (GVIF) for brma objects.
#
# Two complementary diagnostics for multicollinearity:
# 1. Posterior summaries of the design-based VIF/GVIF from cov2cor(vcov)
# 2. Posterior correlation of regression coefficients (Bayesian diagnostic)
#
# VIF is computed from the correlation matrix of the coefficient
# variance-covariance matrix: vcov = (X'WX)^{-1}, calculated separately in
# each posterior draw before computing and summarizing VIF/GVIF. brma.mv
# known-V models use the full supplied sampling covariance plus
# any marginal heterogeneity or formula random-effect covariance in each
# posterior draw.
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
#' @param max_samples maximum posterior draws used for the VIF/GVIF
#' computation in every model class. Defaults to
#' \code{Inf}. Finite values deterministically thin draws across posterior row
#' order. The optional posterior coefficient correlation still uses all draws.
#' @param probs optional credible-interval probabilities for GVIF (ordinary
#' VIF for one-df terms). Defaults to \code{NULL}, which retains the point-summary
#' table. Requested empirical quantiles appear as \code{CI_<probability>} columns.
#' These intervals describe GVIF, not its adjusted root.
#' @param include_samples whether to include a \code{samples} data frame with
#' posterior \code{draw} row indices before thinning, \code{term_index} (the row
#' of the \code{vif} table, also distinguishing repeated display labels),
#' \code{term}, \code{df}, \code{GVIF},
#' and \code{GVIF^(1/(2*df))}. Defaults to \code{FALSE}.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' VIF is computed from the correlation matrix derived from the
#' coefficient variance-covariance matrix \eqn{(X'WX)^{-1}}. For standard
#' meta-regression models, \eqn{W = \mathrm{diag}(w_i/(v_i + \tau^2))}
#' at each posterior heterogeneity draw. Scale and multilevel models use
#' observation-specific \eqn{\tau_i} and block-structured covariance where
#' applicable. VIF/GVIF and its adjusted root are computed in each draw and
#' averaged separately. These posterior means are generally different from
#' applying VIF to an averaged covariance or to a posterior-mean heterogeneity.
#' The mean adjusted GVIF is likewise generally different from transforming
#' the mean GVIF. The diagnostic retains its GLS design-based meaning; it is
#' not calculated from the empirical posterior covariance of the coefficients.
#'
#' A VIF of 1 indicates no collinearity; values above 5 or 10 are
#' commonly considered problematic.
#'
#' The members of a \code{RoBMA()} ensemble can differ in their bias terms (PET
#' adds a \code{PET} column, PEESE a \code{PEESE} column, and the null and
#' weight-function members add none), so each posterior draw uses the design of
#' its member. A location term summarizes all selected posterior draws;
#' \code{PET} and \code{PEESE} summarize only draws in their own class
#' (\code{NA} when no selected draw belongs to it). Their samples are \code{NA}
#' outside that class, not zero, and their credible intervals are conditional
#' on class membership. Location terms keep the specified location design,
#' including terms with spike-and-slab priors.
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
#' Bayesian diagnostic complements VIF: while VIF diagnoses the
#' \emph{potential} for collinearity problems (a data property), the
#' posterior correlation shows the \emph{realized} identification
#' given the data and priors. Informative priors can mitigate
#' collinearity, reducing posterior correlations even when VIF is high.
#'
#' Data-only objects created with \code{only_data = TRUE} are rejected because
#' VIF depends on fitted heterogeneity information, not just the parsed data.
#'
#' For \code{brma.mv()} known-\code{V} models, VIF is computed from the full
#' marginal GLS covariance. Formula random effects are marginalized through
#' BayesTools' random-effect covariance metadata rather than treated as
#' conditioned fitted effects.
#'
#' For known-\code{V} models, the VIF requires applying the inverse of the
#' marginal covariance. A covariance that is too ill-conditioned for that (a
#' reciprocal condition number of about \code{2e-10} or less in any posterior
#' draw used) is unavailable: the call stops with an error of class
#' \code{RoBMA_known_v_precision_unavailable}, whose parent class is
#' \code{RoBMA_numerical_unavailable}, instead of returning a rounded number.
#' The condition carries the draw positions (\code{rows}) and their
#' \code{reciprocal_condition}.
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
#' \item{samples}{The optional long data frame of draw-wise diagnostics when
#' \code{include_samples = TRUE}. Printing and data-frame coercion of the result
#' retain the summary-table view; extract this element explicitly for draws.}
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
                     max_samples = Inf, ..., probs = NULL, include_samples = FALSE) {

  BayesTools::check_bool(posterior_correlation, "posterior_correlation")
  max_samples <- .normalize_max_samples(max_samples, "max_samples")
  BayesTools::check_bool(include_samples, "include_samples")
  BayesTools::check_real(probs, "probs", lower = 0, upper = 1,
                         allow_NULL = TRUE, allow_NA = FALSE, check_length = 0)
  if (anyDuplicated(probs)) {
    stop("'probs' must not contain duplicate probabilities.", call. = FALSE)
  }

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
  vif_df <- .compute_vif(object, max_samples = max_samples, probs = probs,
                         include_samples = include_samples)

  # posterior correlation of regression coefficients
  post_cor <- NULL
  if (posterior_correlation) {
    post_cor <- .compute_posterior_correlation(object)
  }

  output <- list(
    vif                  = vif_df,
    posterior_correlation = post_cor
  )
  if (include_samples) {
    output[["samples"]] <- attr(vif_df, "posterior_samples")
    attr(output[["vif"]], "posterior_samples") <- NULL
  }
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
# Compute VIF/GVIF from cov2cor(solve(X'WX)) in each posterior draw, preserving
# its design class and covariance, then summarize the diagnostic values.
#
# @param object brma object with moderators
#
# @return data frame with columns: term, df, GVIF, GVIF^(1/(2*df))
#
# ---------------------------------------------------------------------------- #
.compute_vif <- function(object, max_samples = Inf, probs = NULL,
                          include_samples = FALSE) {

  X <- .get_model_matrix(object)
  is_known_v <- inherits(object, "brma.mv")
  sample_info <- if (is_known_v) {
    .vif_known_v_sample_info(object, max_samples)
  } else {
    .known_v_diagnostic_posterior_samples(object, max_samples = max_samples,
                                          caller = "VIF")
  }
  draw_rows <- .thin_sample_rows(sample_info[["n_total"]], max_samples)
  if (is.null(draw_rows)) draw_rows <- seq_len(sample_info[["n_total"]])
  bias_classes <- .bias_design_classes(object, X = X,
    posterior_samples = sample_info[["posterior_samples"]])

  values_of <- function(design, rows) {
    assign <- attr(design, "assign")
    keep <- assign != 0L
    groups <- split(seq_len(sum(keep)), assign[keep])
    transform <- function(covariance) {
      .vif_gvif(covariance[keep, keep, drop = FALSE], groups)
    }
    if (is_known_v) {
      return(.vif_vcov_brma_mv(object, design, sample_info, rows, transform))
    }
    .vif_vcov_brma(object, design, rows = rows, transform = transform,
                   posterior_samples = sample_info[["posterior_samples"]])
  }
  template_of <- function(design) {
    if (is.null(design)) return(NULL)
    .vif_table_from_vcov(diag(ncol(design)), attr(design, "assign"), object, design)
  }

  if (is.null(bias_classes)) {
    template <- template_of(X)
    values <- values_of(X, seq_along(draw_rows))
    included <- matrix(TRUE, nrow(values), ncol(values))
  } else {
    templates <- lapply(bias_classes[["designs"]], template_of)
    first_class <- which(!vapply(templates, is.null, logical(1L)))[1L]
    n_location <- nrow(templates[[first_class]]) - (first_class != 1L)
    bias_terms <- which(!vapply(templates, is.null, logical(1L)) & seq_len(3L) != 1L)
    template <- templates[[first_class]][seq_len(n_location), , drop = FALSE]
    for (design_class in bias_terms) {
      template <- rbind(template, templates[[design_class]][n_location + 1L, , drop = FALSE])
    }
    values <- matrix(NA_real_, nrow = length(draw_rows), ncol = nrow(template))
    included <- matrix(FALSE, nrow = length(draw_rows), ncol = nrow(template))
    for (design_class in which(bias_classes[["n_draws"]] > 0L)) {
      design <- bias_classes[["designs"]][[design_class]]
      rows <- which(bias_classes[["class"]] == design_class)
      columns <- c(seq_len(n_location),
        if (design_class != 1L) n_location + match(design_class, bias_terms))
      values[rows, columns] <- values_of(design, rows)
      included[rows, columns] <- TRUE
    }
  }
  vif_df <- .vif_summarize_draws(values, template, draw_rows, probs, include_samples, included)
  if (is_known_v) {
    vif_df <- .known_v_attach_diagnostic_metadata(vif_df,
      .known_v_diagnostic_metadata(sample_info = sample_info))
  }
  vif_df
}


.vif_summarize_draws <- function(values, template, draw_rows, probs = NULL,
                                 include_samples = FALSE,
                                 included = matrix(TRUE, nrow(values), ncol(values))) {

  adjusted <- values^rep(1 / (2 * template[["df"]]), each = nrow(values))
  summarize <- function(x, fun) vapply(seq_len(ncol(x)), function(j) {
    selected <- x[included[, j], j]
    if (!length(selected)) NA_real_ else fun(selected)
  }, numeric(1L))
  template[["GVIF"]] <- summarize(values, mean)
  template[["GVIF^(1/(2*df))"]] <- summarize(adjusted, mean)
  for (prob in probs) {
    template[[paste0("CI_", prob)]] <- summarize(values, function(x) {
      if (anyNA(x)) return(NA_real_)
      stats::quantile(x, probs = prob, names = FALSE)
    })
  }
  rownames(template) <- NULL
  if (include_samples) {
    attr(template, "posterior_samples") <- data.frame(
      draw = rep(draw_rows, times = ncol(values)),
      term_index = rep(seq_len(ncol(values)), each = nrow(values)),
      term = rep(template[["term"]], each = nrow(values)),
      df = rep(template[["df"]], each = nrow(values)),
      GVIF = as.vector(values),
      "GVIF^(1/(2*df))" = as.vector(adjusted), check.names = FALSE
    )
  }
  template
}


# ---------------------------------------------------------------------------- #
# .vif_vcov_brma
# ---------------------------------------------------------------------------- #
#
# Coefficient covariance for standard brma objects.
#
# @param rows optional posterior rows (default: all draws).
# @param transform optional scalar/vector functional evaluated before averaging;
#   returns one row of functional values per selected posterior draw.
#
# ---------------------------------------------------------------------------- #
.vif_vcov_brma <- function(object, X, rows = NULL, transform = NULL,
                            posterior_samples = NULL) {

  # extract sampling variances
  vi      <- .outcome_data_vi(object)
  weights <- .outcome_data_weights(object)
  K       <- length(vi)

  # extract posterior heterogeneity
  is_scale      <- .is_scale(object)
  is_multilevel <- .is_multilevel(object)
  tau_K         <- if (!is_scale && !is_multilevel) 1L else K
  posterior_samples <- .get_posterior_samples(object[["fit"]], posterior_samples)
  if (!is.null(rows)) posterior_samples <- posterior_samples[rows, , drop = FALSE]

  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    posterior_samples = posterior_samples,
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(
      data = object[["data"]],
      "scale"
    ) else NULL,
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    is_random         = .is_random(object),
    K                 = tau_K,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
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
    cluster             = if (is_multilevel) object[["data"]][["outcome"]][["cluster"]] else NULL,
    transform           = transform
  )

  return(vcov)
}


# ---------------------------------------------------------------------------- #
# .vif_known_v_sample_info
# ---------------------------------------------------------------------------- #
#
# Validate a brma.mv() VIF request and select (thin) the posterior draws.
#
# ---------------------------------------------------------------------------- #
.vif_known_v_sample_info <- function(object, max_samples = Inf) {

  if (!.is_data_known_v(object[["data"]])) {
    stop("VIF for brma.mv() requires known-V covariance metadata.",
         call. = FALSE)
  }
  if (!identical(.outcome_type(object), "norm")) {
    stop("VIF for brma.mv() is available only for normal outcome models.",
         call. = FALSE)
  }

  return(.known_v_diagnostic_posterior_samples(
    object      = object,
    max_samples = max_samples,
    caller      = "known-V VIF diagnostics"
  ))
}


# ---------------------------------------------------------------------------- #
# .vif_vcov_brma_mv
# ---------------------------------------------------------------------------- #
#
# Coefficient covariance for brma.mv() known-V models. This is a marginal GLS
# design diagnostic: sampled and marginalized formula random effects contribute
# through ZGZ', not through fitted random-effect means.
#
# @param sample_info selected posterior draws (.vif_known_v_sample_info()).
# @param rows optional rows of the selected draws to average over.
#
# ---------------------------------------------------------------------------- #
.vif_vcov_brma_mv <- function(object, X, sample_info, rows = NULL, transform = NULL) {

  posterior_samples <- sample_info[["posterior_samples"]]
  if (!is.null(rows)) {
    posterior_samples <- posterior_samples[rows, , drop = FALSE]
  }
  known_V <- .data_known_v_data(object[["data"]])
  vcov <- .known_v_with_request_rows(
    .vif_vcov_from_factor_plan(
      object            = object,
      X                 = X,
      posterior_samples = posterior_samples,
      known_V           = known_V,
      transform         = transform
    ),
    rows = rows
  )
  vcov <- .known_v_attach_diagnostic_metadata(
    vcov,
    .known_v_diagnostic_metadata(sample_info = sample_info)
  )

  return(vcov)
}


# ---------------------------------------------------------------------------- #
# .vif_table_from_vcov
# ---------------------------------------------------------------------------- #
#
# Compute GVIF/GVIF-adjusted table from the fixed-effect covariance.
#
# ---------------------------------------------------------------------------- #
.vif_gvif <- function(vcov, groups) {

  p <- ncol(vcov)
  if (p <= 1L) return(rep(1, length(groups)))
  correlation <- stats::cov2cor(vcov)
  determinant <- det(correlation)
  vapply(groups, function(columns) {
    if (length(columns) == p) return(1)
    det(correlation[columns, columns, drop = FALSE]) *
      det(correlation[-columns, -columns, drop = FALSE]) / determinant
  }, numeric(1L), USE.NAMES = FALSE)
}


.vif_table_from_vcov <- function(vcov, assign, object, X) {

  # identify and remove intercept (assign == 0)
  has_intercept <- 0 %in% assign
  if (has_intercept) {
    keep   <- assign != 0
    vcov   <- vcov[keep, keep, drop = FALSE]
    assign <- assign[keep]
  }

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

  groups <- lapply(terms_unique, function(term) which(assign == term))
  gvif_values    <- .vif_gvif(vcov, groups)
  gvif_df_values <- numeric(n_terms)
  df_values      <- integer(n_terms)
  term_names     <- character(n_terms)

  for (k in seq_len(n_terms)) {

    term_index <- terms_unique[k]
    cols       <- which(assign == term_index)
    m          <- length(cols)

    df_values[k]  <- m
    term_names[k] <- if (term_index <= length(term_labels)) term_labels[term_index] else colnames(vcov)[cols[1]]

    gvif_df_values[k] <- gvif_values[k]^(1 / (2 * m))
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


.vif_vcov_from_factor_plan <- function(object, X, posterior_samples, known_V,
                                       transform = NULL) {

  S <- nrow(posterior_samples)
  K <- nrow(X)
  P <- ncol(X)
  plan_data <- .known_v_marginal_factor_plan(
    object            = object,
    posterior_samples = posterior_samples,
    known_V           = known_V
  )
  precision_X <- .known_v_covariance_plan_precision_rhs_batch(
    plan_data = plan_data,
    rhs       = X
  )

  vcov_sum <- matrix(0, nrow = P, ncol = P)
  values <- NULL
  for (s in seq_len(S)) {
    W_X <- matrix(precision_X[s, , ], nrow = K, ncol = P)
    covariance <- .hat_solve_crossprod(
      crossprod(X, W_X),
      rank = attr(X, "rank")
    )
    if (is.null(transform)) {
      vcov_sum <- vcov_sum + covariance
    } else {
      value <- transform(covariance)
      if (is.null(values)) values <- matrix(NA_real_, nrow = S, ncol = length(value))
      values[s, ] <- value
    }
  }

  if (!is.null(transform)) return(values)

  vcov <- vcov_sum / S
  if (!is.null(colnames(X))) {
    dimnames(vcov) <- list(colnames(X), colnames(X))
  }

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
                                       cluster = NULL, transform = NULL) {

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
  values <- NULL

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

    covariance <- .hat_solve_crossprod(
      crossprod(X, WX),
      rank = attr(X, "rank")
    )
    if (is.null(transform)) {
      vcov_sum <- vcov_sum + covariance
    } else {
      value <- transform(covariance)
      if (is.null(values)) values <- matrix(NA_real_, nrow = S, ncol = length(value))
      values[s, ] <- value
    }
  }

  if (!is.null(transform)) return(values)

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

  # the S x P matrix of the mu-formula coefficients on the original scale,
  # without a fixed-zero intercept (as summary())
  samples_mat <- .diagnostic_location_parameter_samples(object)

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

  intervals <- grep("^CI_", names(vif_df), value = TRUE)
  if (all(vif_df$df == 1) && !length(intervals)) {
    # simple VIF display (no GVIF columns needed)
    out        <- vif_df$GVIF
    names(out) <- vif_df$term
    print(round(out, digits))
  } else {
    # full GVIF display with term as rownames
    print_df <- if (all(vif_df$df == 1)) {
      out <- vif_df[, c("GVIF", intervals), drop = FALSE]
      names(out)[1L] <- "VIF"
      out
    } else {
      vif_df[, -1, drop = FALSE]
    }
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


#' @title Convert VIF Results to a Data Frame
#'
#' @description Converts every table displayed by a \code{vif.brma} object to
#' one component-aware long data frame.
#'
#' @param x a \code{vif.brma} object.
#' @param row.names \code{NULL} or a character vector giving the row names.
#' @param optional logical; passed to the final data-frame coercion.
#' @param stringsAsFactors accepted for compatibility with \code{data.frame()}.
#' @param ... unused additional arguments.
#'
#' @return A plain \code{data.frame} with leading \code{component} and
#' \code{parameter} columns.
#'
#' @export
as.data.frame.vif.brma <- function(
    x, row.names = NULL, optional = FALSE, stringsAsFactors = FALSE, ...) {

  vif_table <- x[["vif"]]
  vif_parameter <- vif_table[["term"]]
  if (all(vif_table[["df"]] == 1)) {
    vif_table <- vif_table[, c("GVIF", grep("^CI_", names(vif_table), value = TRUE)), drop = FALSE]
    names(vif_table)[1L] <- "VIF"
  } else {
    vif_table[["term"]] <- NULL
  }
  tables <- list(.output_table_as_long_data_frame(
    table            = vif_table,
    component        = "vif",
    parameter        = vif_parameter,
    stringsAsFactors = stringsAsFactors
  ))

  if (!is.null(x[["posterior_correlation"]])) {
    posterior_correlation <- x[["posterior_correlation"]]
    tables[[length(tables) + 1L]] <- .output_table_as_long_data_frame(
      table            = posterior_correlation,
      component        = "posterior correlation",
      parameter        = rownames(posterior_correlation),
      stringsAsFactors = stringsAsFactors
    )
  }

  output <- .output_bind_long_data_frames(
    tables    = tables,
    row.names = row.names,
    optional  = optional
  )

  return(output)
}
