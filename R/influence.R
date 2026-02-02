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
#' \code{cook.d}, \code{cov.r}, \code{hat}, and weight (normalized PSIS weight mean).}
#' \item{dfbs}{A data frame with DFBETAS values for the model coefficients.}
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#' fit <- add_loo(fit)
#'
#' # compute influence diagnostics
#' inf <- influence(fit)
#' print(inf)
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
      object        = model,
      type          = "marginal",
      return_full_H = FALSE,
      return_se     = FALSE
    )
    hat_samples <- hat_res[["H_diag"]]

    # Hat Values (Posterior Means)
    hat_val <- colMeans(hat_samples)
  }

  ### Compute Individual Diagnostics using Shared Components
  # DFFITS
  if (outcome_type == "norm" && !is_weightfunction) {
    dffits_val <- .dffits_internal(rstudent_val, hat_samples)
  }

  # Cook's Distance
  # Need P (number of coefficients)
  if (outcome_type == "norm" && !is_weightfunction) {
    X <- .get_model_matrix(model)
    P <- ncol(X)
    cook_val <- .cooks.distance_internal(rstudent_val, hat_samples, P)
  }

  # COVRATIO
  # covratio uses weights and beta samples, not hat_samples directly.
  cov_val <- covratio(model)

  # DFBETAS
  # dfbetas uses weights and beta samples.
  dfb_val <- dfbetas(model)

  # Compute tau.del (LOO tau estimate)
  tau_del_val <- dfbetas(model, type = "scale", return_loo_estimates = TRUE)[[1]]

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

  class(res) <- "infl.brma"
  return(res)
}

#' @exportS3Method
print.infl.brma <- function(x, digits = 3, ...) {
  cat("Influence diagnostics:\n")
  print(x$inf, digits = digits, ...)
  cat("\nDFBETAS:\n")
  print(x$dfbs, digits = digits, ...)
  invisible(x)
}
