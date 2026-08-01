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
                                        summarize = FALSE,
                                        max_samples = Inf) {
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
  weights      <- .outcome_data_weights(object)
  K            <- length(yi)

  # get design matrix X
  X <- .get_model_matrix(object)

  if (.is_data_known_v(object[["data"]])) {
    return(.compute_known_v_hat_matrix_samples(
      object             = object,
      conditioning_depth = conditioning_depth,
      return_full_H      = return_full_H,
      return_se          = return_se,
      return_resid       = return_resid,
      summarize          = summarize,
      X                  = X,
      max_samples        = max_samples
    ))
  }

  # get tau samples (heterogeneity)
  # tau_within: estimate-level heterogeneity (used for cluster and estimate residuals)
  # tau_between: cluster-level heterogeneity (only for multilevel models)
  tau_result <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = object[["data"]][["scale"]],
    scale_formula = if (is_scale) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors  = priors[["scale"]],
    is_scale      = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    fixed_tau         = .fixed_tau_prior_value(priors),
    fixed_rho         = .fixed_rho_prior_value(priors)
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
  H_diag_samples    <- matrix(0, nrow = S, ncol = K)
  H_samples         <- if (return_full_H) array(0, dim = c(S, K, K)) else NULL
  se_samples        <- if (return_se && !summarize) matrix(0, nrow = S, ncol = K) else NULL
  resid_samples     <- if (return_resid && !summarize) matrix(0, nrow = S, ncol = K) else NULL
  se_sum            <- if (return_se && summarize) rep(0, K) else NULL
  resid_sum         <- if (return_resid && summarize) rep(0, K) else NULL
  z_sum             <- if (return_se && return_resid && summarize) rep(0, K) else NULL
  M_diag_samples    <- matrix(0, nrow = S, ncol = K) # useful debug/checking
  I_K               <- diag(K)
  zero_residual_rows <- if (return_se || return_resid) {
    .hat_zero_residual_rows(X)
  } else {
    integer()
  }

  # Pre-calculate indices for multilevel blocks to avoid repeating inside loop
  block_indices <- list()
  if (is_multilevel) {
    block_indices <- .get_multilevel_block_indices(ids)
  }

  for (s in seq_len(S)) {
    tau_w_s             <- tau_within_samples[s, ]
    tau_b_s             <- tau_between_samples[s, ]
    sampling_diagonal_s <- vi / weights
    diagonal_s          <- (vi + tau_w_s^2) / weights
    M_diag_s            <- diagonal_s

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
    variance_transform <- NULL
    if (return_se) {
      residual_projection <- I_K - XB %*% t(WX)
      variance_transform   <- residual_projection

      if (conditioning_depth == "estimate") {
        W <- .hat_precision_matrix(
          diagonal      = diagonal_s,
          rank_one      = if (is_multilevel) tau_b_s else NULL,
          block_indices = block_indices
        )
        variance_transform <- sweep(
          W - WX %*% XtWX_inv %*% t(WX),
          1L,
          sampling_diagonal_s,
          "*"
        )

      } else if (conditioning_depth == "cluster" && is_multilevel) {
        W <- .hat_precision_matrix(
          diagonal      = diagonal_s,
          rank_one      = tau_b_s,
          block_indices = block_indices
        )
        between_covariance <- matrix(0, nrow = K, ncol = K)
        for (idx in block_indices) {
          between_covariance[idx, idx] <- tcrossprod(tau_b_s[idx])
        }
        variance_transform <-
          (I_K - between_covariance %*% W) %*% residual_projection
      }
    }

    H_diag_samples[s, ] <- rowSums(XB * WX)

    if (return_full_H) {
      H_samples[s, , ] <- X %*% XtWX_inv %*% t(WX)
    }

    if (return_se || return_resid) {
      residual_s <- NULL
      se_s       <- NULL

      beta_hat <- as.vector(XtWX_inv %*% crossprod(X, Wy))
      residual <- yi - as.vector(X %*% beta_hat)

      if (!is.null(variance_transform)) {
        residual <- as.vector(variance_transform %*% yi)

      } else if (conditioning_depth == "estimate") {
        weighted_residual <- .hat_apply_precision(
          x             = residual,
          diagonal      = diagonal_s,
          rank_one      = if (is_multilevel) tau_b_s else NULL,
          block_indices = block_indices
        )
        residual <- sampling_diagonal_s * weighted_residual

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
      residual[zero_residual_rows] <- 0

      if (return_resid) {
        residual_s <- residual
        if (summarize) {
          resid_sum <- resid_sum + residual
        } else {
          resid_samples[s, ] <- residual
        }
      }

      if (return_se) {
        se2 <- .hat_transformed_covariance_diag(
          transform     = variance_transform,
          diagonal      = diagonal_s,
          rank_one      = if (is_multilevel) tau_b_s else NULL,
          block_indices = block_indices
        )
        se2[zero_residual_rows] <- 0

        se_s <- .hat_variance_sd(se2, "Residual variance")
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
# .compute_known_v_hat_matrix_samples
# ---------------------------------------------------------------------------- #
#
# Known-V residual projection for correlated sampling covariance targets.
#
# ---------------------------------------------------------------------------- #
.compute_known_v_hat_matrix_samples <- function(object, conditioning_depth,
                                                return_full_H, return_se,
                                                return_resid, summarize, X,
                                                max_samples = Inf) {

  sample_info <- .known_v_diagnostic_posterior_samples(
    object      = object,
    max_samples = max_samples,
    caller      = "known-V hat/residual diagnostics"
  )
  setup <- .estimate_likelihood_setup.brma(
    object            = object,
    posterior_samples = sample_info[["posterior_samples"]]
  )
  known_V            <- .data_known_v_data(setup[["data"]])
  yi                 <- setup[["yi"]]
  S                  <- setup[["S"]]
  K                  <- setup[["K"]]
  zero_residual_rows <- if (return_se || return_resid) {
    .hat_zero_residual_rows(X)
  } else {
    integer()
  }

  if (conditioning_depth == "estimate") {
    extra_variance <- .known_v_extra_variance_from_setup(setup)
  }

  H_diag_samples <- matrix(0, nrow = S, ncol = K)
  H_samples      <- if (return_full_H) array(0, dim = c(S, K, K)) else NULL
  se_samples     <- if (return_se && !summarize) matrix(0, nrow = S, ncol = K) else NULL
  resid_samples  <- if (return_resid && !summarize) matrix(0, nrow = S, ncol = K) else NULL
  se_sum         <- if (return_se && summarize) rep(0, K) else NULL
  resid_sum      <- if (return_resid && summarize) rep(0, K) else NULL
  z_sum          <- if (return_se && return_resid && summarize) rep(0, K) else NULL
  M_diag_samples <- matrix(0, nrow = S, ncol = K)
  chunk_info     <- NULL

  if (conditioning_depth == "estimate") {
    offset_samples <- setup[["mu_random"]]
    if (is.null(offset_samples)) {
      offset_samples <- matrix(0, nrow = S, ncol = K)
    }

    for (s in seq_len(S)) {
      y_offset  <- yi - offset_samples[s, ]
      projection <- .known_v_gls_projection_blocks(
        X              = X,
        y              = y_offset,
        known_V        = known_V,
        extra_variance = extra_variance[s, ],
        return_full_H  = return_full_H
      )
      residual   <- projection[["sampling_residual"]]
      residual[zero_residual_rows] <- 0
      projection[["sampling_residual_variance"]][zero_residual_rows] <- 0
      se         <- .hat_variance_sd(
        projection[["sampling_residual_variance"]],
        "Known-V sampling residual variance"
      )

      H_diag_samples[s, ] <- projection[["H_diag"]]
      M_diag_samples[s, ] <- .known_v_diagonal(known_V) + extra_variance[s, ]

      if (return_full_H) {
        H_samples[s, , ] <- projection[["H"]]
      }
      if (return_resid) {
        if (summarize) {
          resid_sum <- resid_sum + residual
        } else {
          resid_samples[s, ] <- residual
        }
      }
      if (return_se) {
        if (summarize) {
          se_sum <- se_sum + se
        } else {
          se_samples[s, ] <- se
        }
      }
      if (summarize && return_se && return_resid) {
        z_s <- residual / se
        z_s[residual == 0 & se == 0] <- 0
        z_sum <- z_sum + z_s
      }
    }

  } else if (conditioning_depth == "marginal") {
    I_K <- diag(K)
    chunk_info <- .known_v_apply_marginal_covariance_chunks(
      object            = object,
      posterior_samples = setup[["posterior_samples"]],
      FUN               = function(covariance_samples, rows) {
        for (j in seq_along(rows)) {
          s          <- rows[j]
          covariance <- covariance_samples[j, , ]
          projection <- .known_v_gls_projection(
            X          = X,
            y          = yi,
            covariance = covariance
          )
          residual <- projection[["residual"]]
          A        <- I_K - projection[["H"]]
          residual <- as.vector(A %*% yi)
          residual[zero_residual_rows] <- 0
          residual_factor <- A %*% projection[["covariance_factor"]]
          se2      <- rowSums(residual_factor^2)
          se2[zero_residual_rows] <- 0
          se       <- .hat_variance_sd(
            se2,
            "Known-V marginal residual variance"
          )

          H_diag_samples[s, ] <<- diag(projection[["H"]])
          M_diag_samples[s, ] <<- diag(projection[["covariance"]])

          if (return_full_H) {
            H_samples[s, , ] <<- projection[["H"]]
          }
          if (return_resid) {
            if (summarize) {
              resid_sum <<- resid_sum + residual
            } else {
              resid_samples[s, ] <<- residual
            }
          }
          if (return_se) {
            if (summarize) {
              se_sum <<- se_sum + se
            } else {
              se_samples[s, ] <<- se
            }
          }
          if (summarize && return_se && return_resid) {
            z_s <- residual / se
            z_s[residual == 0 & se == 0] <- 0
            z_sum <<- z_sum + z_s
          }
        }
      }
    )

  } else {
    .check_cluster_unit_deferred("rstandard()", argument = "conditioning_depth")
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
  result[["known_v_diagnostic"]] <- .known_v_diagnostic_metadata(
    sample_info = sample_info,
    chunk_info  = chunk_info
  )

  return(result)
}


.known_v_gls_projection_blocks <- function(X, y, known_V, extra_variance,
                                           return_full_H = FALSE) {

  K <- .known_v_nrow(known_V)
  p <- ncol(X)
  if (nrow(X) != K || length(y) != K || length(extra_variance) != K) {
    stop("Known-V block GLS inputs have inconsistent dimensions.", call. = FALSE)
  }

  if (any(!is.finite(extra_variance)) || any(extra_variance < 0)) {
    stop("Known-V residual covariance is not positive definite.", call. = FALSE)
  }

  W_X               <- matrix(0, nrow = K, ncol = p)
  W_y               <- numeric(K)
  covariance_blocks <- .known_v_blocks(known_V)
  projection_blocks <- vector("list", length(covariance_blocks))

  for (b in seq_along(covariance_blocks)) {
    block        <- covariance_blocks[[b]]
    index        <- block[["index"]]
    V_block      <- block[["covariance"]]
    block_extra  <- extra_variance[index]
    rank_one     <- .covariance_exact_rank_one_factor(V_block)

    if (!is.null(rank_one) && all(block_extra > 0)) {
      local_index <- seq_along(index)
      W_block <- .hat_precision_matrix(
        diagonal      = block_extra,
        rank_one      = rank_one,
        block_indices = list(local_index)
      )
      covariance_factor <- cbind(
        diag(sqrt(block_extra), nrow = length(index)),
        rank_one
      )
    } else {
      covariance <- V_block
      diag(covariance) <- diag(covariance) + block_extra
      chol_covariance <- .covariance_cholesky(
        .covariance_factorization(covariance)
      )
      if (is.null(chol_covariance)) {
        stop("Known-V residual covariance is not positive definite.",
             call. = FALSE)
      }
      W_block          <- chol2inv(chol_covariance)
      covariance_factor <- t(chol_covariance)
    }

    W_X[index, ] <- W_block %*% X[index, , drop = FALSE]
    W_y[index]   <- as.vector(W_block %*% y[index])
    projection_blocks[[b]] <- list(
      index             = index,
      V                 = V_block,
      W                 = W_block,
      covariance_factor = covariance_factor
    )
  }

  XtWX_inv <- .hat_solve_crossprod(crossprod(X, W_X))
  beta_hat <- as.vector(XtWX_inv %*% crossprod(X, W_y))
  residual <- y - as.vector(X %*% beta_hat)

  sampling_residual <- numeric(K)
  V_W_X             <- matrix(0, nrow = K, ncol = p)
  for (block in projection_blocks) {
    index   <- block[["index"]]
    W_r     <- as.vector(block[["W"]] %*% residual[index])
    V_W_X[index, ] <- block[["V"]] %*% W_X[index, , drop = FALSE]
    sampling_residual[index] <- as.vector(block[["V"]] %*% W_r)
  }

  V_W_X_B                  <- V_W_X %*% XtWX_inv
  sampling_residual_variance <- numeric(K)
  for (block in projection_blocks) {
    index <- block[["index"]]
    W_L   <- block[["W"]] %*% block[["covariance_factor"]]
    residual_factor <- -V_W_X_B %*%
      crossprod(X[index, , drop = FALSE], W_L)
    residual_factor[index, ] <- residual_factor[index, , drop = FALSE] +
      block[["V"]] %*% W_L
    sampling_residual_variance <- sampling_residual_variance +
      rowSums(residual_factor^2)
  }

  H_diag <- rowSums((X %*% XtWX_inv) * W_X)

  out <- list(
    H_diag                       = H_diag,
    sampling_residual            = sampling_residual,
    sampling_residual_variance   = sampling_residual_variance
  )
  if (isTRUE(return_full_H)) {
    out[["H"]] <- X %*% XtWX_inv %*% t(W_X)
  }
  out
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
# .hat_transformed_covariance_diag
# ---------------------------------------------------------------------------- #
#
# Compute diag(T M T') from a diagonal-plus-rank-one covariance factor. This
# avoids subtracting algebraically equivalent positive-semidefinite matrices.
#
# ---------------------------------------------------------------------------- #
.hat_transformed_covariance_diag <- function(transform, diagonal, rank_one,
                                             block_indices) {

  diagonal_factor <- sweep(transform, 2L, sqrt(diagonal), "*")
  variance        <- rowSums(diagonal_factor^2)

  if (!is.null(rank_one)) {
    for (idx in block_indices) {
      block_factor <- as.vector(
        transform[, idx, drop = FALSE] %*% rank_one[idx]
      )
      variance <- variance + block_factor^2
    }
  }

  return(variance)
}


# ---------------------------------------------------------------------------- #
# .hat_zero_residual_rows
# ---------------------------------------------------------------------------- #
#
# Rows whose removal lowers the design rank are fitted exactly by every
# weighted projection: e_i belongs to the column space of X, so both their
# residual and residual variance are identically zero.
#
# ---------------------------------------------------------------------------- #
.hat_zero_residual_rows <- function(X) {

  zero_tolerance <- max(dim(X)) * .Machine$double.eps
  full_rank      <- qr(X, tol = zero_tolerance)[["rank"]]
  return(which(vapply(
    seq_len(nrow(X)),
    function(i) {
      qr(X[-i, , drop = FALSE], tol = zero_tolerance)[["rank"]] < full_rank
    },
    logical(1L)
  )))
}


.hat_variance_sd <- function(variance, context) {

  if (!is.numeric(variance) || any(!is.finite(variance)) ||
      any(variance < 0)) {
    stop(context, " must be finite and non-negative.", call. = FALSE)
  }

  return(sqrt(variance))
}
