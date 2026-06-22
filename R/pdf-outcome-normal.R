# ============================================================================ #
# Normal Outcome Log-Likelihoods
# ============================================================================ #



# ---------------------------------------------------------------------------- #
# .outcome_pdf.norm
# ---------------------------------------------------------------------------- #
#
# Compute pointwise log-likelihoods for normal outcome models.
#
# For normal outcome models, the observed effect y_i has likelihood:
#   y_i ~ N(mu_i, tau_within^2 + se_i^2)
#
# This function computes the log-density for each observation at each posterior
# sample.
#
# @param yi               numeric vector of length K; observed effect sizes
# @param mu_samples       S x K matrix of location samples (with cluster effects if multilevel)
# @param tau_within       S x K matrix of estimate-level heterogeneity samples
# @param sei              numeric vector of length K; standard errors
#
# @return S x K matrix of log-likelihood values
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.norm <- function(yi, mu_samples, tau_within, sei) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate yi and sei across samples: matrix(vec, S, K, byrow = TRUE)
  yi_mat  <- matrix(yi,  nrow = S, ncol = K, byrow = TRUE)
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)

  # compute total SD: sqrt(tau^2 + se^2)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  # compute log-likelihood for each cell
  log_lik <- stats::dnorm(yi_mat, mean = mu_samples, sd = total_sd, log = TRUE)

  return(log_lik)
}



# ---------------------------------------------------------------------------- #
# .outcome_pdf.selnorm
# ---------------------------------------------------------------------------- #
#
# Selected-normal likelihood using the compiled p-bin kernel.
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.selnorm <- function(yi, mu_samples, tau_within, sei,
                                 selection_context, weights = NULL,
                                 selection_sei = sei) {

  S        <- nrow(mu_samples)
  K        <- ncol(mu_samples)
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  return(.selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu_samples,
    sigma_num      = total_sd,
    mu_norm        = mu_samples,
    sigma_norm     = total_sd,
    sei            = selection_sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]],
    weights        = weights
  ))
}


.outcome_pdf_sum.norm <- function(yi, mu_samples, tau_within, sei,
                                  weights = NULL) {

  if (.has_native_norm_loglik_row_sum(selection = FALSE, cluster = FALSE)) {
    weights_arg <- if (is.null(weights)) NULL else .native_numeric_vector(weights)
    return(.Call(
      "RoBMA_norm_loglik_row_sum",
      .native_numeric_vector(yi),
      .native_numeric_matrix(mu_samples),
      .native_numeric_matrix(tau_within),
      .native_numeric_vector(sei),
      weights_arg,
      PACKAGE = "RoBMA"
    ))
  }

  log_lik <- .outcome_pdf.norm(
    yi         = yi,
    mu_samples = mu_samples,
    tau_within = tau_within,
    sei        = sei
  )

  return(rowSums(.apply_log_lik_weights(log_lik, weights)))
}


.outcome_pdf_sum.selnorm <- function(yi, mu_samples, tau_within, sei,
                                     selection_context, weights = NULL,
                                     selection_sei = sei) {

  S        <- nrow(mu_samples)
  K        <- ncol(mu_samples)
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  return(.selnorm_kernel_loglik_row_sum(
    yi             = yi,
    mu_num         = mu_samples,
    sigma_num      = total_sd,
    mu_norm        = mu_samples,
    sigma_norm     = total_sd,
    sei            = selection_sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]],
    weights        = weights
  ))
}
