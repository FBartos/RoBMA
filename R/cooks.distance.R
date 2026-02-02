# ============================================================================ #
# brma.cooks.distance.R
# ============================================================================ #
#
# This file contains the cooks.distance method for brma model fits. Cook's
# distance measures the aggregate impact of an observation on all model
# coefficients.
#
# ============================================================================ #


#' @title Cook's Distance for brma Objects
#'
#' @description Computes Cook's distance for a fitted brma object. Cook's
#' distance measures the aggregate influence of each observation on the
#' model coefficients.
#'
#' @param model a fitted brma object.
#' @param ... additional arguments (currently ignored).
#'
#' @details
#' Cook's distance is computed using a hybrid Bayesian approach. Like DFFITS,
#' we calculate the distribution sample-wise to capture the non-linear relationship
#' between leverage and influence.
#'
#' \deqn{D_i = \frac{r_i^2}{P} \times \frac{h_i}{1 - h_i}}
#'
#' where \eqn{P} is the number of regression coefficients (beta parameters),
#' \eqn{r_i} is the LOO-PIT residual, and \eqn{h_i} is the hat value (leverage).
#'
#' @return A numeric vector of Cook's distance values, one for each observation.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#' fit <- add_loo(fit)
#'
#' # compute Cook's distance
#' cooks.distance(fit)
#' }
#'
#' @seealso \code{\link{influence.brma}}, \code{\link{dffits.brma}}, \code{\link{hatvalues.brma}}
#' @exportS3Method
cooks.distance.brma <- function(model, ...) {

  # the function relies on hatvalues
  # as such it is sensible only sensible for normal models
  outcome_type       <- .outcome_type(model)
  is_weightfunction  <- .is_weightfunction(model)

  if (outcome_type != "norm") {
    stop("cooks.distance is only available for normal outcome models.", call. = FALSE)
  }
  if (is_weightfunction) {
    stop("cooks.distance is not available for selection models (weightfunction).", call. = FALSE)
  }

  # 1. Get rstudent (LOO-PIT residuals) - Vector of length K
  r_res        <- rstudent(model)
  rstudent_vec <- r_res[["z"]]

  # 2. Identify P (number of coefficients)
  X <- .get_model_matrix(model)
  P <- ncol(X)

  # 3. Get hat matrix samples (S x K)
  hat_res <- .compute_hat_matrix_samples(
    object    = model,
    type      = "marginal",
    return_full_H = FALSE,
    return_se = FALSE
  )
  hat_samples <- hat_res[["H_diag"]]

  # 4. Call internal function
  d_vec <- .cooks.distance_internal(rstudent_vec, hat_samples, P)

  return(d_vec)
}

.cooks.distance_internal <- function(rstudent_vec, hat_samples, P) {

  K <- length(rstudent_vec)
  S <- nrow(hat_samples)

  # 4. Calculate Cook's Distance Samples
  # Broadcast rstudent squared
  rstudent_sq_mat <- matrix(rstudent_vec^2, nrow = S, ncol = K, byrow = TRUE)

  # Safety Clip Hat Values
  hat_samples_safe <- pmin(hat_samples, 0.9999)

  # Calculate D samples:
  # D_{s,i} = (rstudent_i^2 / P) * (h_{s,i} / (1 - h_{s,i}))
  factor_mat <- hat_samples_safe / (1 - hat_samples_safe)

  d_samples <- (rstudent_sq_mat / P) * factor_mat

  # 5. Aggregate: Compute column means
  d_vec <- colMeans(d_samples)
  names(d_vec) <- names(rstudent_vec)

  return(d_vec)
}
