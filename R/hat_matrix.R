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
# - M_diag: S x K matrix of marginal variance diagonals
#
# @param object brma object
# @param type character; type of residuals/hatvalues to compute:
#               - "marginal" (default): Uses marginal variance M
#               - "study": Uses within-study variance M_within
#               - "estimate": Uses sampling variance V
# @param return_full_H logical; whether to return the full hat matrix for each sample
# @param return_se logical; whether to compute and return standard error components
#
# ---------------------------------------------------------------------------- #
.compute_hat_matrix_samples <- function(object, type = "marginal", return_full_H = FALSE, return_se = FALSE) {
  # check inputs
  type <- match.arg(type, c("marginal", "study", "estimate"))

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
  # tau_within: estimate-level heterogeneity (used for study and estimate residuals)
  # tau_between: study-level heterogeneity (only for multilevel models)
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

  # get study_ids for multilevel structure
  ids <- NULL
  if (is_multilevel) {
    ids <- outcome_data[["study_ids"]]
  }

  # initialize outputs
  H_diag_samples <- matrix(0, nrow = S, ncol = K)
  H_samples      <- if (return_full_H) array(0, dim = c(S, K, K)) else NULL
  se_samples     <- if (return_se) matrix(0, nrow = S, ncol = K) else NULL
  M_diag_samples <- matrix(0, nrow = S, ncol = K) # useful debug/checking
  V              <- diag(vi, nrow = K, ncol = K)
  I_K            <- diag(K)

  # Pre-calculate indices for multilevel blocks to avoid repeating inside loop
  block_indices <- list()
  if (is_multilevel) {
    block_indices <- .get_multilevel_block_indices(ids)
  }

  for (s in seq_len(S)) {
    # 1. Construct full Marginal Variance Matrix M
    tau2_w_s <- tau_within_samples[s, ]^2
    tau_b_s  <- tau_between_samples[s, ]

    if (is_multilevel) {
      M <- .build_multilevel_marginal_covariance(
        tau_within    = tau_within_samples[s, ],
        tau_between   = tau_b_s,
        vi            = vi,
        block_indices = block_indices
      )
    } else {
      M <- diag(vi + tau2_w_s, nrow = K, ncol = K)
    }

    M_diag_samples[s, ] <- diag(M)

    # 2. Compute Weight matrix W and Hat matrix H
    # W = M^{-1}
    # Use Cholesky for stability and validity
    chk <- try(chol(M), silent = TRUE)
    if (inherits(chk, "try-error")) {
      # If Cholesky fails, M is likely not positive definite -> invalid sample
      # For now, we return NA which will likely propagate to NA in outputs
      W <- matrix(NA, nrow = K, ncol = K)
    } else {
      W <- chol2inv(chk)
    }

    # H = X (X' W X)^{-1} X' W
    # Note: X' W is K x K if computed directly, but X is K x p.
    # XtW should be p x K.
    XtW <- crossprod(X, W) # X' W
    XtWX <- XtW %*% X # X' W X

    # Invert XtWX (p x p)
    XtWX_inv <- tryCatch(solve(XtWX), error = function(e) MASS::ginv(XtWX))

    # H = X * (XtWX_inv) * XtW
    H <- X %*% XtWX_inv %*% XtW

    # store leverage (diagonal)
    H_diag_samples[s, ] <- diag(H)

    # store full H if requested
    if (return_full_H) {
      H_samples[s, , ] <- H
    }

    # 3. Compute Standard Errors if requested
    if (return_se) {
      ImH <- I_K - H

      if (type == "marginal") {
        # Marginal: Var(e) = (I-H) M (I-H)'
        ve <- ImH %*% M %*% t(ImH)
        se_samples[s, ] <- sqrt(pmax(diag(ve), 0))
      } else if (type == "study") {
        # Study: Var(e) = (I-H) M_within (I-H)'
        # M_within = diag(vi + tau_within^2)
        M_within <- diag(vi + tau2_w_s, nrow = K, ncol = K)

        ve <- ImH %*% M_within %*% t(ImH)
        se_samples[s, ] <- sqrt(pmax(diag(ve), 0))
      } else if (type == "estimate") {
        # Estimate/Conditional:
        # Metafor's multilevel conditional residual variance is:
        # V W (I-H) M (I-H)' W V, where V is the sampling covariance.
        residual_v <- V %*% W %*% ImH %*% M %*% t(ImH) %*% W %*% V
        se_samples[s, ] <- sqrt(pmax(diag(residual_v), 0))
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
    result[["se"]] <- se_samples
  }

  return(result)
}
