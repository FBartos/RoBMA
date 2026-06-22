# ============================================================================ #
# brma.summary_heterogeneity.R
# ============================================================================ #
#
# This file contains the summary_heterogeneity method for brma model fits.
# It computes absolute and relative heterogeneity measures:
# - tau, tau^2: Absolute heterogeneity (SD and variance)
# - I^2: Percentage of total variance due to heterogeneity
# - H^2: Ratio of total variance to sampling variance
#
# For multilevel (3-level) models, heterogeneity is partitioned into:
# - tau [within]: Estimate-level heterogeneity
# - tau [between]: Cluster-level heterogeneity
# - I^2 [within], I^2 [between]: Partitioned relative heterogeneity
#
# Formulas follow metafor documentation and Higgins & Thompson (2002).
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# summary_heterogeneity generic
# ---------------------------------------------------------------------------- #

#' @title Summary of Heterogeneity
#'
#' @description Computes the absolute heterogeneity (tau, tau^2) and
#' relative measures of heterogeneity (I^2, H^2) for a fitted model.
#'
#' @param object a fitted model object
#' @param ... additional arguments passed to methods
#'
#' @return Method-specific return value containing heterogeneity estimates.
#'
#' @seealso [pooled_heterogeneity()], [summary.brma()]
#' @export
summary_heterogeneity <- function(object, ...) {

  UseMethod("summary_heterogeneity")
}


# ---------------------------------------------------------------------------- #
# summary_heterogeneity.brma
# ---------------------------------------------------------------------------- #

#' @title Summary of Heterogeneity for brma Objects
#'
#' @description Computes the absolute heterogeneity (tau, tau^2) and
#' relative measures of heterogeneity (I^2, H^2) for a fitted brma object.
#'
#' @param object a fitted brma object
#' @param probs quantiles of the posterior distribution to be displayed.
#' Defaults to \code{c(.025, .975)} for 95% credible intervals.
#' @param component heterogeneity component to return for \code{brma.mv()}
#' models. Defaults to \code{"all"}. Use \code{"total"} for the
#' variance-additive total heterogeneity. Random-formula allocation nodes can
#' be selected by their displayed name, such as \code{"study/esid"}.
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' For standard (2-level) random-effects models, the function reports:
#' \itemize{
#'   \item \code{tau}: Between-study standard deviation
#'   \item \code{tau2}: Between-study variance
#'   \item \code{I2}: Percentage of total variance due to heterogeneity
#'   \item \code{H2}: Ratio of total to sampling variance
#' }
#'
#' For multilevel (3-level) models with nested effects, the function additionally
#' partitions heterogeneity into estimate-level and cluster-level components:
#' \itemize{
#'   \item \code{rho}: Proportion of heterogeneity variance allocated to clusters
#'   \item \code{tau [within]}: Estimate-level standard deviation
#'   \item \code{tau [between]}: Cluster-level standard deviation
#'   \item \code{tau2 [within]}: Estimate-level variance
#'   \item \code{tau2 [between]}: Cluster-level variance
#'   \item \code{I2 [within]}: Percentage of variance due to estimate-level heterogeneity
#'   \item \code{I2 [between]}: Percentage of variance due to cluster-level heterogeneity
#' }
#'
#' For location-scale models, \code{tau2} aggregates the observation-specific
#' heterogeneity variances \eqn{\tau_i^2}; the corresponding \code{tau} summary
#' is the square root of this aggregate variance. The relative \eqn{I^2} and
#' \eqn{H^2} measures average the observation-specific indices.
#'
#' The I^2 and H^2 statistics are computed following the metafor package
#' implementation, using the "typical" sampling variance formula from
#' \insertCite{higgins2002quantifying;textual}{RoBMA}. For multilevel models,
#' the partitioned I^2 follows the approach described in the metafor documentation.
#'
#' For \code{brma.mv()} models, the method reports absolute heterogeneity
#' summaries (\code{tau}, \code{tau2}) for selected random-effect components.
#' Relative \eqn{I^2} and \eqn{H^2} summaries are not reported for general
#' known-V covariance structures. For \code{component = "total"}, independent
#' component variances are summed before summarizing \code{tau2}; \code{tau}
#' is the square root of this variance draw.
#' When a random-formula model uses a shared total-SD plus variance-allocation
#' node, \code{component = "all"} also includes an allocation-node table with
#' \code{tau}, \code{tau2}, and \code{var_frac(...)} rows. For nested formulas
#' such as \code{random = ~ 1 | study / esid}, this table is displayed under
#' the user-facing path \code{"study/esid"}.
#' For heterogeneous structured random effects, such as
#' \code{random = ~ har(time | study)}, the allocation table reports
#' \code{sd_total(...)}, level-specific \code{sd(...)} rows, and
#' \code{var_ratio(...)} rows for the heterogeneous standard-deviation
#' components.
#'
#' @return A list of class \code{summary_heterogeneity.brma} containing:
#' \itemize{
#'   \item \code{estimates}: A \code{BayesTools_table} with heterogeneity statistics
#' }
#' For decomposed \code{brma.mv()} models with multiple selected components,
#' a named list of such summary objects is returned.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     data    = dat.lehmann2018,
#'     measure = "SMD",
#'     seed    = 1,
#'     silent  = TRUE
#'   )
#'
#'   summary_heterogeneity(fit)
#' }
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [pooled_heterogeneity()], [summary.brma()]
#' @export
summary_heterogeneity.brma <- function(object, probs = c(.025, .975),
                                       component = "all", ...) {

  # input validation
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)
  if (inherits(object, "brma.mv")) {
    return(
      .summary_heterogeneity_brma_mv(
        object    = object,
        probs     = probs,
        component = component,
        ...
      )
    )
  }

  .check_univariate_heterogeneity_component(component)

  # extract model characteristics
  is_multilevel <- .is_multilevel(object)
  is_scale      <- .is_scale(object)

  # extract sampling variances using helper
  vi <- .outcome_data_vi(object)
  K  <- length(vi)

  # compute typical sampling variance (v_tilde)
  # following metafor formula (Higgins & Thompson, 2002, eq. 9)
  # using the generalized formula with projection matrix for meta-regression
  X <- .get_model_matrix(object)
  p <- ncol(X)
  W <- diag(1/vi, nrow = K)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  v_tilde <- (K - p) / sum(diag(P))

  posterior_samples <- .get_posterior_samples(object[["fit"]])

  # extract tau samples using the evaluate helper
  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = posterior_samples
  )

  samples_list <- .summary_heterogeneity_samples(
    tau_within_samples  = tau_result[["tau_within"]],
    tau_between_samples = tau_result[["tau_between"]],
    rho_samples         = if (is_multilevel) posterior_samples[, "rho"] else NULL,
    v_tilde             = v_tilde,
    is_multilevel       = is_multilevel
  )

  # generate summary table
  estimates <- BayesTools::ensemble_estimates_table(
    samples    = samples_list,
    parameters = names(samples_list),
    probs      = probs,
    title      = "Heterogeneity Estimates:"
  )

  # create output object
  output <- list(
    estimates = estimates
  )

  class(output) <- "summary_heterogeneity.brma"

  return(output)
}


# ---------------------------------------------------------------------------- #
# .summary_heterogeneity_samples
# ---------------------------------------------------------------------------- #
#
# Compute heterogeneity summaries from observation-level tau samples.
#
.summary_heterogeneity_samples <- function(tau_within_samples, tau_between_samples,
                                           v_tilde, is_multilevel,
                                           rho_samples = NULL) {

  if (is.null(tau_between_samples)) {
    tau_between_samples <- matrix(0, nrow = nrow(tau_within_samples),
                                  ncol = ncol(tau_within_samples))
  }

  if (!is.matrix(tau_within_samples) || !is.matrix(tau_between_samples)) {
    stop("'tau_within_samples' and 'tau_between_samples' must be matrices.",
         call. = FALSE)
  }
  if (!identical(dim(tau_within_samples), dim(tau_between_samples))) {
    stop("'tau_within_samples' and 'tau_between_samples' must have matching dimensions.",
         call. = FALSE)
  }
  if (!is.numeric(v_tilde) || length(v_tilde) != 1 || !is.finite(v_tilde) || v_tilde <= 0) {
    stop("'v_tilde' must be a positive finite number.", call. = FALSE)
  }
  if (is_multilevel && !is.null(rho_samples) &&
      (!is.numeric(rho_samples) || length(rho_samples) != nrow(tau_within_samples))) {
    stop("'rho_samples' must be a numeric vector matching posterior sample rows.",
         call. = FALSE)
  }

  sigma2_within_matrix  <- tau_within_samples^2
  sigma2_between_matrix <- tau_between_samples^2
  sigma2_total_matrix   <- sigma2_within_matrix + sigma2_between_matrix

  sigma2_within  <- rowMeans(sigma2_within_matrix)
  sigma2_between <- rowMeans(sigma2_between_matrix)
  sigma2_total   <- sigma2_within + sigma2_between

  denominator_matrix <- sigma2_total_matrix + v_tilde

  I2_total <- rowMeans(100 * sigma2_total_matrix / denominator_matrix)
  H2       <- rowMeans(denominator_matrix / v_tilde)

  if (is_multilevel) {

    I2_within  <- rowMeans(100 * sigma2_within_matrix / denominator_matrix)
    I2_between <- rowMeans(100 * sigma2_between_matrix / denominator_matrix)
    if (is.null(rho_samples)) {
      rho_samples <- ifelse(sigma2_total > 0, sigma2_between / sigma2_total, 0)
    } else {
      rho_samples <- pmin(pmax(rho_samples, 0), 1)
    }

    return(list(
      "tau"             = sqrt(sigma2_total),
      "rho"             = rho_samples,
      "tau [within]"    = sqrt(sigma2_within),
      "tau [between]"   = sqrt(sigma2_between),
      "tau2"            = sigma2_total,
      "tau2 [within]"   = sigma2_within,
      "tau2 [between]"  = sigma2_between,
      "I2"              = I2_total,
      "I2 [within]"     = I2_within,
      "I2 [between]"    = I2_between,
      "H2"              = H2
    ))
  }

  return(list(
    "tau"  = sqrt(sigma2_total),
    "tau2" = sigma2_total,
    "I2"   = I2_total,
    "H2"   = H2
  ))
}


# ---------------------------------------------------------------------------- #
# print method for summary.brma_heterogeneity
# ---------------------------------------------------------------------------- #

#' @title Print Summary of Heterogeneity
#'
#' @description Prints the heterogeneity summary table.
#'
#' @param x a \code{summary_heterogeneity.brma} object
#' @param ... additional arguments (currently ignored)
#'
#' @return Returns \code{x} invisibly.
#'
#' @export
print.summary_heterogeneity.brma <- function(x, ...) {

  cat("\n")
  print(x$estimates)
  cat("\n")

  return(invisible(x))
}


#' @rdname summary_heterogeneity.brma
#' @export
print.summary_heterogeneity.brma_list <- function(x, ...) {

  cat("\n")
  for (i in seq_along(x)) {
    print(x[[i]][["estimates"]])
    cat("\n")
  }

  return(invisible(x))
}
