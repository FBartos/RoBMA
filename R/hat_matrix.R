# ============================================================================ #
# brma.hat_matrix.R
# ============================================================================ #
#
# This file contains internal functions for computing the hat matrix and related
# quantities for brma objects. These functions are shared between residuals()
# (for rstandard) and hatvalues().
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .compute_hat_matrix_samples
# ---------------------------------------------------------------------------- #
#
# Computes hat matrix components for each posterior sample.
#
# Returns a list containing:
# - H_diag: S x K matrix of hat matrix diagonals (leverages)
# - H: S x K x K array of full hat matrices (only if return_full_H = TRUE)
# - se_components: S x K matrix of standard errors components for residuals (if return_se = TRUE)
# - resid_components: S x K matrix of GLS residuals (if return_resid = TRUE)
# - M_diag: S x K matrix of marginal variance diagonals
#
# @param object brma object
# @param conditioning_depth character; conditioning depth for residual SEs:
#                           - "marginal" (default): Uses marginal variance M
#                           - "cluster": Uses within-cluster variance M_within
#                           - "estimate": Uses sampling variance V
# @param return_full_H logical; whether to return the full hat matrix for each sample
# @param return_se logical; whether to compute and return standard error components
# @param return_resid logical; whether to compute and return GLS residual components
#
# ---------------------------------------------------------------------------- #
.compute_hat_matrix_samples <- function(object, conditioning_depth = "marginal",
                                        return_full_H = FALSE,
                                        return_se = FALSE, return_resid = FALSE,
                                        summarize = FALSE) {
  # check inputs
  conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)

  # get model characteristics
  is_multilevel <- .is_multilevel(object)
  is_scale      <- .is_scale(object)
  priors        <- object[["priors"]]

  # get observed data
  outcome_data <- object[["data"]][["outcome"]]
  yi           <- .outcome_data_yi(object)
  vi           <- .outcome_data_vi(object)
  K            <- length(yi)

  # get design matrix X
  X <- .get_model_matrix(object)

  # get tau samples (heterogeneity)
  # tau_within: estimate-level heterogeneity (used for cluster and estimate residuals)
  # tau_between: cluster-level heterogeneity (only for multilevel models)
  tau_result <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = object[["data"]][["scale"]],
    scale_formula = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors  = priors[["scale"]],
    is_scale      = is_scale,
    is_multilevel = is_multilevel,
    K             = K
  )
  tau_within_samples  <- tau_result[["tau_within"]]
  tau_between_samples <- tau_result[["tau_between"]]

  S <- nrow(tau_within_samples)

  # get cluster for multilevel structure
  ids <- NULL
  if (is_multilevel) {
    ids <- outcome_data[["cluster"]]
  }

  # initialize outputs
  H_diag_samples <- matrix(0, nrow = S, ncol = K)
  H_samples      <- if (return_full_H) array(0, dim = c(S, K, K)) else NULL
  se_samples     <- if (return_se && !summarize) matrix(0, nrow = S, ncol = K) else NULL
  resid_samples  <- if (return_resid && !summarize) matrix(0, nrow = S, ncol = K) else NULL
  se_sum         <- if (return_se && summarize) rep(0, K) else NULL
  resid_sum      <- if (return_resid && summarize) rep(0, K) else NULL
  z_sum          <- if (return_se && return_resid && summarize) rep(0, K) else NULL
  M_diag_samples <- matrix(0, nrow = S, ncol = K) # useful debug/checking
  V              <- diag(vi, nrow = K, ncol = K)
  I_K            <- diag(K)
  residual_tol   <- 100 * .Machine$double.eps * max(1, max(abs(yi)))

  # Pre-calculate indices for multilevel blocks to avoid repeating inside loop
  block_indices <- list()
  if (is_multilevel) {
    block_indices <- .get_multilevel_block_indices(ids)
  }

  for (s in seq_len(S)) {
    tau_w_s    <- tau_within_samples[s, ]
    tau_b_s    <- tau_between_samples[s, ]
    diagonal_s <- vi + tau_w_s^2
    M_diag_s   <- diagonal_s

    if (is_multilevel) {
      M_diag_s <- M_diag_s + tau_b_s^2
    }

    M_diag_samples[s, ] <- M_diag_s

    WX <- .hat_apply_precision(
      x             = X,
      diagonal      = diagonal_s,
      rank_one      = if (is_multilevel) tau_b_s else NULL,
      block_indices = block_indices
    )
    Wy <- .hat_apply_precision(
      x             = yi,
      diagonal      = diagonal_s,
      rank_one      = if (is_multilevel) tau_b_s else NULL,
      block_indices = block_indices
    )

    XtWX     <- crossprod(X, WX)
    XtWX_inv <- .hat_solve_crossprod(XtWX)
    XB       <- X %*% XtWX_inv
    Q_diag   <- rowSums(XB * X)

    H_diag_samples[s, ] <- rowSums(XB * WX)

    if (return_full_H) {
      H_samples[s, , ] <- X %*% XtWX_inv %*% t(WX)
    }

    if (return_se || return_resid) {
      residual_s <- NULL
      se_s       <- NULL

      beta_hat <- as.vector(XtWX_inv %*% crossprod(X, Wy))
      residual <- yi - as.vector(X %*% beta_hat)

      if (conditioning_depth == "estimate") {
        weighted_residual <- .hat_apply_precision(
          x             = residual,
          diagonal      = diagonal_s,
          rank_one      = if (is_multilevel) tau_b_s else NULL,
          block_indices = block_indices
        )
        residual <- vi * weighted_residual

      } else if (conditioning_depth == "cluster" && is_multilevel) {
        weighted_residual <- .hat_apply_precision(
          x             = residual,
          diagonal      = diagonal_s,
          rank_one      = tau_b_s,
          block_indices = block_indices
        )
        cluster_adjust <- rep(0, K)
        for (idx in block_indices) {
          cluster_adjust[idx] <- tau_b_s[idx] * sum(tau_b_s[idx] * weighted_residual[idx])
        }
        residual <- residual - cluster_adjust
      }

      if (return_resid) {
        residual[abs(residual) < residual_tol] <- 0
        residual_s <- residual
        if (summarize) {
          resid_sum <- resid_sum + residual
        } else {
          resid_samples[s, ] <- residual
        }
      }

      if (return_se) {
        if (conditioning_depth == "marginal" ||
            (conditioning_depth == "cluster" && !is_multilevel)) {
          se2 <- M_diag_s - Q_diag

        } else if (conditioning_depth == "estimate") {
          W_diag  <- .hat_precision_diag(
            diagonal      = diagonal_s,
            rank_one      = if (is_multilevel) tau_b_s else NULL,
            block_indices = block_indices
          )
          QW_diag <- rowSums((WX %*% XtWX_inv) * WX)
          se2     <- vi^2 * (W_diag - QW_diag)

        } else {
          se2 <- .hat_cluster_se2(
            X             = X,
            XtWX_inv      = XtWX_inv,
            diagonal      = diagonal_s,
            tau_between   = tau_b_s,
            vi            = vi,
            block_indices = block_indices,
            I_K           = I_K
          )
        }

        se_s <- sqrt(pmax(se2, 0))
        if (summarize) {
          se_sum <- se_sum + se_s
        } else {
          se_samples[s, ] <- se_s
        }
      }

      if (summarize && return_se && return_resid) {
        z_s <- residual_s / se_s
        z_s[residual_s == 0 & se_s == 0] <- 0
        z_sum <- z_sum + z_s
      }
    }
  }

  result <- list(
    H_diag = H_diag_samples,
    M_diag = M_diag_samples
  )

  if (return_full_H) {
    result[["H"]] <- H_samples
  }
  if (return_se) {
    result[["se"]] <- if (summarize) se_sum / S else se_samples
  }
  if (return_resid) {
    result[["resid"]] <- if (summarize) resid_sum / S else resid_samples
  }
  if (summarize && return_se && return_resid) {
    result[["z"]] <- z_sum / S
  }

  return(result)
}


# ---------------------------------------------------------------------------- #
# .hat_apply_precision
# ---------------------------------------------------------------------------- #
#
# Apply the inverse marginal covariance to a vector or matrix.
#
# ---------------------------------------------------------------------------- #
.hat_apply_precision <- function(x, diagonal, rank_one = NULL, block_indices = list()) {

  was_vector <- is.null(dim(x))
  if (was_vector) {
    x <- matrix(x, ncol = 1)
  }

  out <- matrix(0, nrow = nrow(x), ncol = ncol(x))

  if (is.null(rank_one)) {
    out <- x / diagonal
  } else {
    for (idx in block_indices) {
      inv_diag <- 1 / diagonal[idx]
      inv_x    <- x[idx, , drop = FALSE] * inv_diag
      inv_u    <- rank_one[idx] * inv_diag
      denom    <- 1 + sum(rank_one[idx] * inv_u)

      out[idx, ] <- inv_x -
        tcrossprod(inv_u, colSums(rank_one[idx] * inv_x) / denom)
    }
  }

  if (was_vector) {
    return(as.vector(out))
  }
  return(out)
}


# ---------------------------------------------------------------------------- #
# .hat_precision_diag
# ---------------------------------------------------------------------------- #
#
# Diagonal of the inverse marginal covariance.
#
# ---------------------------------------------------------------------------- #
.hat_precision_diag <- function(diagonal, rank_one = NULL, block_indices = list()) {

  if (is.null(rank_one)) {
    return(1 / diagonal)
  }

  out <- rep(NA_real_, length(diagonal))
  for (idx in block_indices) {
    inv_diag <- 1 / diagonal[idx]
    inv_u    <- rank_one[idx] * inv_diag
    denom    <- 1 + sum(rank_one[idx] * inv_u)

    out[idx] <- inv_diag - inv_u^2 / denom
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .hat_precision_matrix
# ---------------------------------------------------------------------------- #
#
# Build the full inverse marginal covariance only for rare full-matrix paths.
#
# ---------------------------------------------------------------------------- #
.hat_precision_matrix <- function(diagonal, rank_one = NULL, block_indices = list()) {

  K <- length(diagonal)
  W <- matrix(0, nrow = K, ncol = K)

  if (is.null(rank_one)) {
    diag(W) <- 1 / diagonal
    return(W)
  }

  for (idx in block_indices) {
    inv_diag <- 1 / diagonal[idx]
    inv_u    <- rank_one[idx] * inv_diag
    denom    <- 1 + sum(rank_one[idx] * inv_u)

    W[idx, idx] <- diag(inv_diag, nrow = length(idx), ncol = length(idx)) -
      tcrossprod(inv_u) / denom
  }

  return(W)
}


# ---------------------------------------------------------------------------- #
# .hat_solve_crossprod
# ---------------------------------------------------------------------------- #
#
# Stable inverse for the small fixed-effect crossproduct.
#
# ---------------------------------------------------------------------------- #
.hat_solve_crossprod <- function(x) {

  chk <- try(chol(x), silent = TRUE)
  if (!inherits(chk, "try-error")) {
    return(chol2inv(chk))
  }

  return(tryCatch(solve(x), error = function(e) MASS::ginv(x)))
}


# ---------------------------------------------------------------------------- #
# .hat_cluster_se2
# ---------------------------------------------------------------------------- #
#
# Fallback cluster-level residual variance using block-structured W and M.
#
# ---------------------------------------------------------------------------- #
.hat_cluster_se2 <- function(X, XtWX_inv, diagonal, tau_between, vi,
                             block_indices, I_K) {

  K <- length(vi)
  W <- .hat_precision_matrix(
    diagonal      = diagonal,
    rank_one      = tau_between,
    block_indices = block_indices
  )
  M <- .build_multilevel_marginal_covariance(
    tau_within    = sqrt(pmax(diagonal - vi, 0)),
    tau_between   = tau_between,
    vi            = vi,
    block_indices = block_indices
  )

  between_cov <- matrix(0, nrow = K, ncol = K)
  for (idx in block_indices) {
    between_cov[idx, idx] <- tcrossprod(tau_between[idx])
  }

  Q <- X %*% XtWX_inv %*% t(X)
  A <- (I_K - between_cov %*% W)

  return(diag(A %*% (M - Q) %*% t(A)))
}
