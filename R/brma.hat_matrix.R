# TODO:
# Assess whether this separate file is actually needed once all function is
# implemented (the functionality might be merged into residual directly?)
# (unless it is reused in other influence diagnostics)


# ---------------------------------------------------------------------------- #
# brma.hat_matrix.R
# ---------------------------------------------------------------------------- #
#
# Utilities to compute the hat matrix for brma models. The core function
# `.hat_matrix()` computes the hat matrix H = X (X' W X)^{-1} X' W given a
# model matrix X and a weight matrix W. A convenience wrapper
# `.compute_hat_matrix()` builds W from sampling variances and tau^2 and uses
# `.get_model_matrix()` to extract X from a brma object.
#
# This file intentionally contains only hat-matrix related helpers so we can
# later swap or extend the computation depending on model / likelihood type.
#
# ---------------------------------------------------------------------------- #

#' Compute hat matrix from model matrix and weight matrix
#'
#' @param X numeric matrix; model matrix with K rows (observations) and p cols
#' @param W numeric matrix; weight matrix (usually diagonal) of dimension K x K
#'
#' @return numeric matrix H (K x K) hat matrix
#' @keywords internal
.hat_matrix <- function(X, W) {

  # Xt W X
  XtWX <- crossprod(X, W %*% X)

  # invert robustly (use generalized inverse on failure)
  XtWX_inv <- tryCatch(
    solve(XtWX),
    error = function(e) MASS::ginv(XtWX)
  )

  # H = X (X' W X)^{-1} X' W
  H <- X %*% XtWX_inv %*% crossprod(X, W)

  return(H)
}


#' Convenience wrapper to compute hat matrix for a brma object and a sample of tau^2
#'
#' Builds the weight matrix W = diag(1 / (vi + tau2_s)) where `vi` are the
#' sampling variances returned by `.outcome_data_vi()` and `tau2_s` is a length-
#' K numeric vector of tau^2 values for one posterior sample. The function
#' extracts the model matrix using `.get_model_matrix()` and returns H.
#'
#' @param object brma object
#' @param tau2_s numeric vector of length K with tau^2 values for a single sample
#'
#' @return numeric matrix H (K x K)
#' @keywords internal
.compute_hat_matrix <- function(object, tau2_s) {

  # model matrix
  X <- .get_model_matrix(object)

  # sampling variances
  vi <- .outcome_data_vi(object)

  # construct weight matrix W = diag(1 / (vi + tau2_s))
  W <- diag(1 / (vi + tau2_s))

  H <- .hat_matrix(X, W)

  return(H)
}
