#' @title Hat Values for brma Objects
#'
#' @description Computes hat values (leverages) from a fitted brma object.
#' Returns a matrix of hat values, one for each observation and each posterior sample.
#'
#' @param model a fitted brma object
#' @param ... additional arguments (currently ignored)
#'
#' @details
#' Hat values (leverages) measure the influence of each observation on the
#' fitted values. In a Bayesian meta-analysis, the random effects variance \eqn{\tau^2}
#' is uncertain, so the hat matrix depends on the posterior samples of \eqn{\tau^2}.
#'
#' This function computes the diagonal elements of the hat matrix:
#' \deqn{h_i = (X (X^T W X)^{-1} X^T W)_{ii}}
#' where \eqn{W} is the weight matrix inverse to the marginal variance matrix.
#'
#' The result is a matrix of coordinates, where rows correspond to posterior samples
#' and columns correspond to observations.
#'
#' @return A matrix of hat values with dimensions \code{S x K}, where \code{S} is the number
#' of posterior samples and \code{K} is the number of observations.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # get hat values
#' hv <- hatvalues(fit)
#'
#' # compute posterior mean leverage
#' colMeans(hv)
#' }
#'
#' @exportS3Method
hatvalues.brma <- function(model, ...) {

  # the function relies on normal-normal approximations
  # as such it is sensible only sensible for normal models
  outcome_type       <- .outcome_type(model)
  is_weightfunction  <- .is_weightfunction(model)

  if (outcome_type != "norm") {
    stop("hatvalues is only available for normal outcome models.", call. = FALSE)
  }
  if (is_weightfunction) {
    stop("hatvalues is not available for selection models (weightfunction).", call. = FALSE)
  }

  # compute hat matrix samples
  # returns list(H_diag, M_diag, ...)
  # we only need H_diag
  res <- .compute_hat_matrix_samples(
    object        = model,
    type          = "marginal",
    return_full_H = FALSE,
    return_se     = FALSE
  )

  # extract the diagonal
  res <- res[["H_diag"]]

  # return only the posterior mean for consistency with metafor
  res <- colMeans(res)

  return(res)
}
