#' @title Hat Values for brma Objects
#'
#' @description Computes hat values (leverages) from a fitted brma object.
#' Returns posterior mean leverages, one for each observation.
#'
#' @param model a fitted brma object
#' @param max_samples maximum posterior draws used for known-\code{V}
#' \code{brma.mv()} marginal covariance computations. Defaults to \code{Inf}.
#' Finite values deterministically thin draws across posterior row order.
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
#' The hat matrix is computed for each posterior draw and then averaged over
#' draws, matching the vector output shape used by `metafor`.
#' For known-\code{V} \code{brma.mv()} random-formula models, leverages use the
#' marginal GLS covariance \eqn{V + ZGZ'}.
#'
#' This method is available only for normal outcome models without
#' weightfunction selection.
#'
#' @return A numeric vector of length `K`, one posterior mean leverage per
#' observation.
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
#'
#'   fit <- brma(yi = yi, vi = vi, data = dat, measure = "RR")
#'   hatvalues(fit)
#' }
#' }
#'
#' @exportS3Method
hatvalues.brma <- function(model, max_samples = Inf, ...) {

  max_samples <- .normalize_max_samples(max_samples, "max_samples")

  if (inherits(model, "only_priors.brma") || is.null(model[["fit"]])) {
    stop("hatvalues require a fitted brma object with posterior samples.",
         call. = FALSE)
  }

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
  if (inherits(model, "brma.mv") && !.is_data_known_v(model[["data"]])) {
    stop("hatvalues for brma.mv() require known-V covariance metadata.",
         call. = FALSE)
  }
  if (!inherits(model, "brma.mv")) {
    .check_random_formula_postfit_deferred(model, "hatvalues()")
  }

  # compute hat matrix samples
  # returns list(H_diag, M_diag, ...)
  # we only need H_diag
  res <- .compute_hat_matrix_samples(
    object             = model,
    conditioning_depth = "marginal",
    return_full_H      = FALSE,
    return_se          = FALSE,
    max_samples        = max_samples
  )

  # extract the diagonal
  known_v_metadata <- res[["known_v_diagnostic"]]
  res              <- res[["H_diag"]]

  # return only the posterior mean for consistency with metafor
  res <- colMeans(res)
  res <- .diagnostic_set_names(res, model)
  if (!is.null(known_v_metadata)) {
    res <- .known_v_attach_diagnostic_metadata(res, known_v_metadata)
  }

  return(res)
}
