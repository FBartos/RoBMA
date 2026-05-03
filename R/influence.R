# ============================================================================ #
# brma.influence.R
# ============================================================================ #
#
# This file contains the influence method for brma model fits. It computes
# various diagnostics (DFFITS, Cook's distance, COVRATIO, etc.) to identify
# influential observations.
#
# ============================================================================ #


#' @title Measure Influence for brma Objects
#'
#' @description Computes DFFITS, Cook's distance, COVRATIO, and other
#' influence diagnostics for a fitted brma object.
#'
#' @param model a fitted brma object.
#' @param ... additional arguments (currently ignored).
#'
#' @return An object of class \code{"infl.brma"}, which corresponds to the
#' structure of \code{metafor::influence} objects. It is a list containing:
#' \item{inf}{A data frame with columns: \code{rstudent}, \code{dffits},
#' \code{cook.d}, \code{cov.r}, \code{tau.del}, and \code{hat}, where available.
#' The \code{tau.del} column is the PSIS leave-one-out posterior mean of the
#' aggregate heterogeneity. For scale models, aggregation is over the remaining
#' scale-model rows after each deletion.}
#' \item{dfbs}{A data frame with DFBETAS values for the model coefficients.}
#' Undefined determinant- or variance-standardized diagnostics are reported as
#' \code{NaN} and printed with an explanatory note.
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
#'   fit <- add_loo(fit)
#'
#'   inf <- influence(fit)
#'   print(inf)
#' }
#' }
#'
#' @seealso \code{\link{dffits.brma}}, \code{\link{cooks.distance.brma}},
#' \code{\link{covratio.brma}}, \code{\link{dfbetas.brma}}, \code{\link{hatvalues.brma}}
#' @exportS3Method
influence.brma <- function(model, ...) {

  # object information
  # hatvalues, dffit, and cooks distance are possible only for normal-normal models
  outcome_type       <- .outcome_type(model)
  is_weightfunction  <- .is_weightfunction(model)


  ### Precomputed Shared Components
  # rstudent (LOO-PIT residuals)
  rstud_df     <- rstudent(model)
  rstudent_val <- rstud_df[["z"]]

  # Hat Matrix Samples (S x K)
  if (outcome_type == "norm" && !is_weightfunction) {
    hat_res <- .compute_hat_matrix_samples(
      object             = model,
      conditioning_depth = "marginal",
      return_full_H      = FALSE,
      return_se          = FALSE
    )
    hat_samples <- hat_res[["H_diag"]]

    # Hat Values (Posterior Means)
    hat_val <- colMeans(hat_samples)

    fit_samples <- .influence_fit_samples(model)
    loo_wts     <- loo_weights(model)
  }

  ### Compute Individual Diagnostics using Shared Components
  # DFFITS
  if (outcome_type == "norm" && !is_weightfunction) {
    dffits_val <- .dffits_internal(fit_samples, loo_wts)
  }

  # Cook's Distance
  # Need P (number of coefficients)
  if (outcome_type == "norm" && !is_weightfunction) {
    X        <- .get_model_matrix(model)
    P        <- qr(X)[["rank"]]
    cook_val <- .cooks.distance_internal(fit_samples, loo_wts, P)
  }

  # COVRATIO
  # covratio uses weights and beta samples, not hat_samples directly.
  cov_val <- covratio(model)
  cov_note <- attr(cov_val, "note")
  cov_val <- as.numeric(cov_val)

  # DFBETAS
  # dfbetas uses weights and beta samples.
  dfb_val <- dfbetas(model)

  # Compute tau.del (LOO aggregate heterogeneity estimate)
  tau_del_val <- .influence_tau_del(model)

  # Construct 'inf' data frame
  if (outcome_type == "norm" && !is_weightfunction) {
    inf_df <- data.frame(
      rstudent = rstudent_val,
      dffits   = dffits_val,
      cook.d   = cook_val,
      cov.r    = cov_val,
      tau.del  = tau_del_val,
      hat      = hat_val
    )
  } else {
    inf_df <- data.frame(
      rstudent = rstudent_val,
      cov.r    = cov_val,
      tau.del  = tau_del_val
    )
  }


  ### Return list structure matching metafor
  res <- list(
    inf  = inf_df,
    dfbs = dfb_val
  )
  if (!is.null(cov_note)) {
    attr(res, "note") <- cov_note
  }

  class(res) <- "infl.brma"
  return(res)
}


# ---------------------------------------------------------------------------- #
# .influence_tau_del
# ---------------------------------------------------------------------------- #
#
# PSIS leave-one-out aggregate tau estimates for influence diagnostics.
#
.influence_tau_del <- function(model) {

  weights           <- loo_weights(model)
  is_scale          <- .is_scale(model)
  is_multilevel     <- .is_multilevel(model)
  K                 <- ncol(weights)
  posterior_samples <- .get_posterior_samples(model[["fit"]])

  tau_result <- .evaluate.brma.tau(
    fit               = model[["fit"]],
    scale_data        = model[["data"]][["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = model[["data"]], "scale") else NULL,
    scale_priors      = model[["priors"]][["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = posterior_samples
  )
  tau_total <- sqrt(tau_result[["tau_within"]]^2 + tau_result[["tau_between"]]^2)

  return(.influence_tau_del_from_samples(tau_total, weights))
}


# ---------------------------------------------------------------------------- #
# .influence_tau_del_from_samples
# ---------------------------------------------------------------------------- #
#
# @param tau_samples S x K matrix of observation-level total tau samples.
# @param weights S x K matrix of normalized PSIS weights.
#
.influence_tau_del_from_samples <- function(tau_samples, weights) {

  if (!is.matrix(tau_samples)) {
    tau_samples <- as.matrix(tau_samples)
  }
  if (!is.matrix(weights)) {
    weights <- as.matrix(weights)
  }
  if (!identical(dim(tau_samples), dim(weights))) {
    stop("'tau_samples' and 'weights' must have the same dimensions.",
         call. = FALSE)
  }

  K <- ncol(tau_samples)
  if (K > 1L) {
    tau_deleted <- (rowSums(tau_samples) - tau_samples) / (K - 1L)
  } else {
    tau_deleted <- tau_samples
  }

  return(unname(colSums(weights * tau_deleted)))
}

#' @exportS3Method
print.infl.brma <- function(x, digits = 3, ...) {
  cat("Influence diagnostics:\n")
  print(x$inf, digits = digits, ...)
  cat("\nDFBETAS:\n")
  print(x$dfbs, digits = digits, ...)
  .print_diagnostic_note(attr(x, "note"))
  invisible(x)
}
