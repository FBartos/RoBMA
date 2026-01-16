# ============================================================================ #
# Outcome CDF Functions (Cumulative Distribution)
# ============================================================================ #
#
# These functions compute pointwise CDF values F(yi | theta) for each
# observation and posterior sample. They are used for LOO-PIT residuals
# via probability integral transformation.
#
# Parallels the structure of brma.pdf.R but returns CDF values instead of
# density values.
#
# Note: For binomial and Poisson models, each "observation" consists of a pair
# of data points (ai+ci or x1i+x2i) that together define a single effect size
# estimate. The CDF is computed using the implied normal approximation for the
# effect size (log-OR or log-IRR).
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .outcome_cdf.norm
# ---------------------------------------------------------------------------- #
#
# Compute pointwise CDF values for normal outcome models.
#
# For normal outcome models, the observed effect y_i has distribution:
#   y_i ~ N(mu_i, tau_within^2 + se_i^2)
#
# This function computes F(y_i | mu_i, sigma_i) for each observation at each
# posterior sample.
#
# @param yi               numeric vector of length K; observed effect sizes
# @param mu_samples       S x K matrix of location samples (with study effects if multilevel)
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param sei              numeric vector of length K; standard errors
#
# @return S x K matrix of CDF values in (0, 1)
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.norm <- function(yi, mu_samples, tau_within, sei) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate yi and sei across samples: matrix(vec, S, K, byrow = TRUE)
  yi_mat  <- matrix(yi,  nrow = S, ncol = K, byrow = TRUE)
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)

  # compute total SD: sqrt(tau^2 + se^2)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  # compute CDF value for each cell
  cdf_vals <- stats::pnorm(yi_mat, mean = mu_samples, sd = total_sd)

  return(cdf_vals)
}


# ---------------------------------------------------------------------------- #
# .outcome_cdf.wnorm
# ---------------------------------------------------------------------------- #
#
# Compute pointwise CDF values for weighted normal distribution (selection models).
#
# For selection models, the observed effect y_i follows a weighted normal:
#   y_i ~ f(y) * omega(y) / Z
# where f(y) is normal density and omega(y) is the selection weight function.
#
# This uses the .pwnorm_fast.ss function that handles the weighted CDF.
#
# @param yi               numeric vector of length K; observed effect sizes
# @param mu_samples       S x K matrix of location samples
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param sei              numeric vector of length K; standard errors
# @param omega            S x W matrix of omega (weight) samples
# @param crit_yi          W x K matrix of critical values for each observation
#
# @return S x K matrix of CDF values from weighted distribution
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.wnorm <- function(yi, mu_samples, tau_within, sei, omega, crit_yi) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # pre-compute total SD for each observation
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  # compute CDF values for each observation using weighted normal
  # (loop is necessary due to observation-specific crit_yi values)
  cdf_vals <- matrix(NA_real_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    cdf_vals[, k] <- .pwnorm_fast.ss(
      q      = rep(yi[k], S),
      mean   = mu_samples[, k],
      sd     = total_sd[, k],
      omega  = omega,
      crit_x = crit_yi[, k]
    )
  }

  return(cdf_vals)
}


# ---------------------------------------------------------------------------- #
# .outcome_cdf.binom
# ---------------------------------------------------------------------------- #
#
# Compute pointwise CDF values for binomial outcome models.
#
# For binomial outcome models, we use a normal approximation based on the
# implied log-odds ratio effect size and its approximate sampling variance.
# This approach is consistent with how metafor computes residuals for GLMM.
#
# The marginal distribution of the effect size approximation is:
#   y_i ~ N(mu_i + theta_i * tau, sigma_i)
# where sigma_i is the approximate sampling SE from the cell counts.
#
# @param yi               numeric vector of length K; approximate log-OR effect sizes
# @param sei              numeric vector of length K; approximate sampling SEs
# @param mu_samples       S x K matrix of log-odds ratio samples (includes theta contribution)
#
# @return S x K matrix of CDF values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.binom <- function(yi, sei, mu_samples) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate yi and sei across samples
  yi_mat  <- matrix(yi,  nrow = S, ncol = K, byrow = TRUE)
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)

  # compute CDF using normal approximation
  # note: tau contribution is already in mu_samples via theta
  cdf_vals <- stats::pnorm(yi_mat, mean = mu_samples, sd = sei_mat)

  return(cdf_vals)
}


# ---------------------------------------------------------------------------- #
# .outcome_cdf.pois
# ---------------------------------------------------------------------------- #
#
# Compute pointwise CDF values for Poisson outcome models.
#
# For Poisson outcome models, we use a normal approximation based on the
# implied log incidence rate ratio effect size and its approximate sampling
# variance. This is consistent with how metafor computes residuals for GLMM.
#
# The marginal distribution of the effect size approximation is:
#   y_i ~ N(mu_i + theta_i * tau, sigma_i)
# where sigma_i is the approximate sampling SE from the counts.
#
# @param yi               numeric vector of length K; approximate log-IRR effect sizes
# @param sei              numeric vector of length K; approximate sampling SEs
# @param mu_samples       S x K matrix of log-IRR samples (includes theta contribution)
#
# @return S x K matrix of CDF values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.pois <- function(yi, sei, mu_samples) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate yi and sei across samples
  yi_mat  <- matrix(yi,  nrow = S, ncol = K, byrow = TRUE)
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)

  # compute CDF using normal approximation
  # note: tau contribution is already in mu_samples via theta
  cdf_vals <- stats::pnorm(yi_mat, mean = mu_samples, sd = sei_mat)

  return(cdf_vals)
}


# ---------------------------------------------------------------------------- #
# .cdf.brma
# ---------------------------------------------------------------------------- #
#
# Compute the full CDF matrix for a brma object.
#
# This is the main dispatcher function that coordinates extraction of posterior
# samples and computation of pointwise CDF values using the appropriate
# CDF function for the outcome type.
#
# The CDF is computed at the estimate level: for multilevel models, we condition
# on the fitted study-level random effects (gamma); for selection models, we
# condition on the omega samples.
#
# @param object           brma object
#
# @return S x K matrix of CDF values in (0, 1)
#
# ---------------------------------------------------------------------------- #
.cdf.brma <- function(object) {

  ### extract priors and structural information about the model
  priors            <- object[["priors"]]
  data              <- object[["data"]]
  is_mods           <- .is_mods(object)
  is_multilevel     <- .is_multilevel(object)
  is_scale          <- .is_scale(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)

  ### extract outcome data and fit data
  outcome_data <- data[["outcome"]]
  fit_data     <- .create_fit_data(data = data, priors = priors)
  K            <- nrow(outcome_data)

  ### obtain tau samples using helper function
  tau_result          <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = data[["scale"]],
    scale_formula = if (is_scale) .create_fit_formula_list(data = data, "scale") else NULL,
    scale_priors  = priors[["scale"]],
    is_scale      = is_scale,
    is_multilevel = is_multilevel,
    K             = K
  )
  tau_within_samples  <- tau_result[["tau_within"]]
  tau_between_samples <- tau_result[["tau_between"]]

  ### get the base mu samples using helper function
  # for CDF evaluation, we do NOT adjust for publication bias
  # (we want to evaluate CDF under the assumed model, which includes bias)
  mu_samples <- .evaluate.brma.mu(
    fit               = object[["fit"]],
    outcome_data      = outcome_data,
    mods_data         = data[["mods"]],
    mods_formula      = if (is_mods) .create_fit_formula_list(data = data, "mods") else NULL,
    mods_priors       = priors[["mods"]],
    is_mods           = is_mods,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction,
    bias_adjusted     = FALSE,  # include PET/PEESE terms in CDF
    K                 = K
  )

  ### include 3-level (study-level) random-effects for multilevel models
  if (is_multilevel) {
    study_contribution <- .evaluate.brma.study_effects(
      fit              = object[["fit"]],
      tau_between      = tau_between_samples,
      study_ids        = fit_data[["study_ids"]],
      same_data        = TRUE,  # use fitted gamma values
      effect_direction = effect_direction
    )
    mu_samples <- mu_samples + study_contribution
  }

  ### obtain outcome data: yi and sei
  yi  <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)

  ### dispatch to appropriate CDF function based on outcome type
  if (outcome_type == "norm") {

    # for PET, PEESE, and selection models, the outcome and crit_yi is computed in "positive" space
    # (yi flipped for negative effect direction in .create_fit_data)
    # so we need to flip mu_samples to match
    if (effect_direction == "negative") {
      mu_samples <- -mu_samples
      yi         <- -yi
    }

    if (is_weightfunction) {

      # extract omega samples for weight function
      posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
      omega_samples     <- posterior_samples[, grep("omega", colnames(posterior_samples)), drop = FALSE]

      cdf_vals <- .outcome_cdf.wnorm(
        yi         = yi,
        mu_samples = mu_samples,
        tau_within = tau_within_samples,
        sei        = sei,
        omega      = omega_samples,
        crit_yi    = fit_data$crit_yi
      )

    } else {

      # standard normal CDF
      cdf_vals <- .outcome_cdf.norm(
        yi         = yi,
        mu_samples = mu_samples,
        tau_within = tau_within_samples,
        sei        = sei
      )

    }

    # flip CDF for negative effect direction
    # since F(-y) under flipped model = 1 - F(y) under original
    if (effect_direction == "negative") {
      cdf_vals <- 1 - cdf_vals
    }

  } else if (outcome_type == "bin") {

    # add GLMM random effects (theta * tau_within)
    theta_contribution <- .evaluate.brma.theta.glmm(
      fit        = object[["fit"]],
      tau_within = tau_within_samples,
      same_data  = TRUE,  # use fitted theta values
      K          = K
    )
    mu_samples <- mu_samples + theta_contribution

    cdf_vals <- .outcome_cdf.binom(
      yi         = yi,
      sei        = sei,
      mu_samples = mu_samples
    )

  } else if (outcome_type == "pois") {

    # add GLMM random effects (theta * tau_within)
    theta_contribution <- .evaluate.brma.theta.glmm(
      fit        = object[["fit"]],
      tau_within = tau_within_samples,
      same_data  = TRUE,  # use fitted theta values
      K          = K
    )
    mu_samples <- mu_samples + theta_contribution

    cdf_vals <- .outcome_cdf.pois(
      yi         = yi,
      sei        = sei,
      mu_samples = mu_samples
    )

  }

  # add column names for observations
  colnames(cdf_vals) <- paste0("cdf[", seq_len(K), "]")

  return(cdf_vals)
}
