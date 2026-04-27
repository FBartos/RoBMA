# ============================================================================ #
# brma.dffits.R
# ============================================================================ #
#
# This file contains the dffits method for brma model fits. DFFITS measures
# the number of standard deviations that the fitted value changes if observation
# i is removed.
#
# ============================================================================ #


#' @export
dffits <- function(model, ...) UseMethod("dffits")

#' @title DFFITS for brma Objects
#'
#' @description Computes DFFITS (Difference in FITS, standardized) for a
#' fitted brma object. DFFITS measures how much the fitted value for observation
#' \eqn{i} changes if observation \eqn{i} is removed, standardized by the
#' estimated standard error of the fit.
#'
#' @param model a fitted brma object.
#' @param ... additional arguments (currently ignored).
#'
#' @details
#' DFFITS values are computed using a hybrid Bayesian approach. Since the leverage
#' (hat values) is uncertain in Bayesian models (depending on \eqn{\tau^2}),
#' we calculate the distribution of DFFITS values across posterior samples and
#' report the posterior mean.
#'
#' The computation uses the LOO-PIT residuals (which account for estimation
#' uncertainty and leverage) combined with the full posterior distribution of
#' the hat matrix diagonals.
#'
#' \deqn{DFFITS_i = r_i \times \sqrt{\frac{h_i}{1 - h_i}}}
#'
#' where \eqn{r_i} is the LOO-PIT residual (Studentized residual) and \eqn{h_i}
#' is the hat value (leverage). In the Bayesian implementation, this is averaged
#' over the posterior samples of \eqn{h_i}.
#'
#' @return A numeric vector of DFFITS values, one for each observation.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#' fit <- add_loo(fit)
#'
#' # compute DFFITS
#' dffits(fit)
#' }
#'
#' @seealso \code{\link{influence.brma}}, \code{\link{cooks.distance.brma}}, \code{\link{hatvalues.brma}}
#' @exportS3Method
dffits.brma <- function(model, ...) {

    # the function relies on the normal-normal hat matrix
    outcome_type      <- .outcome_type(model)
    is_weightfunction <- .is_weightfunction(model)

    if (outcome_type != "norm") {
      stop("dffits is only available for normal outcome models.", call. = FALSE)
    }
    if (is_weightfunction) {
      stop("dffits is not available for selection models (weightfunction).", call. = FALSE)
    }

    # 1. Get rstudent (LOO-PIT residuals) - Vector of length K
    # We use the "z" column (standardized residuals)
    r_res <- rstudent(model)
    rstudent_vec <- r_res[["z"]]

    # 2. Get hat matrix samples (S x K)
    # returns list(H_diag, ...)
    hat_res <- .compute_hat_matrix_samples(
        object             = model,
        conditioning_depth = "marginal",
        return_full_H      = FALSE,
        return_se          = FALSE
    )
    hat_samples <- hat_res[["H_diag"]]
    
    # 3. Call internal function
    dffits_vec <- .dffits_internal(rstudent_vec, hat_samples)

    return(dffits_vec)
}

.dffits_internal <- function(rstudent_vec, hat_samples) {
    
    K <- length(rstudent_vec)
    S <- nrow(hat_samples)

    # Broadcast: Expand the rstudent vector (Length K) into a matrix (Dimensions S x K)
    # to match hat_samples.
    rstudent_mat <- matrix(rstudent_vec, nrow = S, ncol = K, byrow = TRUE)

    # Leverage-one cases do not have a defined linear-model deletion diagnostic:
    # deleting the point removes the information that identifies its fitted cell.
    hat_samples_safe <- pmax(hat_samples, 0)
    hat_samples_safe[hat_samples_safe >= 1 - sqrt(.Machine$double.eps)] <- NA_real_

    # Calculate Leverage Factor
    # factor_{s,i} = sqrt(h_{s,i} / (1 - h_{s,i}))
    factor_mat <- sqrt(hat_samples_safe / (1 - hat_samples_safe))

    # Compute Samples
    # DFFITS_{s,i} = rstudent_i * factor_{s,i}
    dffits_samples <- rstudent_mat * factor_mat

    # Aggregate: Compute column means to return a vector of length K
    dffits_vec <- colMeans(dffits_samples, na.rm = FALSE)
    names(dffits_vec) <- names(rstudent_vec)
    
    return(dffits_vec)
}
