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
# - tau [within]: Estimate-level (within-study) heterogeneity
# - tau [between]: Study-level (between-study) heterogeneity
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
#' partitions heterogeneity into within-study and between-study components:
#' \itemize{
#'   \item \code{tau [within]}: Estimate-level (within-study) standard deviation
#'   \item \code{tau [between]}: Study-level (between-study) standard deviation
#'   \item \code{I2 [within]}: Percentage of variance due to within-study heterogeneity
#'   \item \code{I2 [between]}: Percentage of variance due to between-study heterogeneity
#' }
#'
#' The I^2 and H^2 statistics are computed following the metafor package
#' implementation, using the "typical" sampling variance formula from
#' \insertCite{higgins2002quantifying;textual}{RoBMA} For multilevel models,
#' the partitioned I^2 follows the approach described in the metafor documentation.
#'
#' @return A list of class \code{summary_heterogeneity.brma} containing:
#' \itemize{
#'   \item \code{estimates}: A \code{BayesTools_table} with heterogeneity statistics
#' }
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get heterogeneity summary
#' summary_heterogeneity(fit)
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [pooled_heterogeneity()], [summary.brma()]
#' @export
summary_heterogeneity.brma <- function(object, probs = c(.025, .975), ...) {

  # input validation
  BayesTools::check_real(probs, "probs", allow_NULL = TRUE, check_length = 0)

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

  # extract tau samples using the evaluate helper
  tau_result <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = object[["data"]][["scale"]],
    scale_formula = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors  = object[["priors"]][["scale"]],
    is_scale      = is_scale,
    is_multilevel = is_multilevel,
    K             = K
  )

  # aggregate tau samples across observations (for scale models)
  tau_within_samples  <- rowMeans(tau_result[["tau_within"]])
  tau_between_samples <- rowMeans(tau_result[["tau_between"]])

  # compute heterogeneity indices
  if (is_multilevel) {

    # multilevel model: partition heterogeneity
    sigma2_within  <- tau_within_samples^2
    sigma2_between <- tau_between_samples^2
    sigma2_total   <- sigma2_within + sigma2_between

    tau_total <- sqrt(sigma2_total)
    tau2      <- sigma2_total

    # partitioned I^2
    I2_total   <- 100 * sigma2_total   / (sigma2_total + v_tilde)
    I2_within  <- 100 * sigma2_within  / (sigma2_total + v_tilde)
    I2_between <- 100 * sigma2_between / (sigma2_total + v_tilde)

    # H^2
    H2 <- (sigma2_total + v_tilde) / v_tilde

    # build samples list for multilevel output
    samples_list <- list(
      "tau"             = tau_total,
      "tau [within]"    = sqrt(sigma2_within),
      "tau [between]"   = sqrt(sigma2_between),
      "tau2"            = tau2,
      "tau2 [within]"   = sigma2_within,
      "tau2 [between]"  = sigma2_between,
      "I2"              = I2_total,
      "I2 [within]"     = I2_within,
      "I2 [between]"    = I2_between,
      "H2"              = H2
    )

  } else {

    # standard (non-multilevel) model
    tau_samples  <- tau_within_samples
    tau2_samples <- tau_samples^2

    # I^2 and H^2
    I2 <- 100 * tau2_samples / (tau2_samples + v_tilde)
    H2 <- (tau2_samples + v_tilde) / v_tilde

    # build samples list
    samples_list <- list(
      "tau"  = tau_samples,
      "tau2" = tau2_samples,
      "I2"   = I2,
      "H2"   = H2
    )
  }

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
