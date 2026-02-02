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
# Performance optimization: When use_normal is provided, samples with
# use_normal[i] = TRUE use the fast normal path (pnorm), skipping the
# expensive weighted CDF computation. This is correct because these samples
# have omega = all 1s (set by JAGS for non-weightfunction models).
#
# @param yi               numeric vector of length K; observed effect sizes
# @param mu_samples       S x K matrix of location samples
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param sei              numeric vector of length K; standard errors
# @param omega            S x W matrix of omega (weight) samples
# @param crit_yi          W x K matrix of critical values for each observation
# @param use_normal       optional logical vector of length S; TRUE if sample should
#                         use fast normal path (for samples with omega = all 1s)
#
# @return S x K matrix of CDF values from weighted distribution
#
# ---------------------------------------------------------------------------- #
.outcome_cdf.wnorm <- function(yi, mu_samples, tau_within, sei, omega, crit_yi,
                               use_normal = NULL) {

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
      q          = rep(yi[k], S),
      mean       = mu_samples[, k],
      sd         = total_sd[, k],
      omega      = omega,
      crit_x     = crit_yi[, k],
      use_normal = use_normal  # <-- Pass through for internal subdispatch
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
# This function computes pointwise CDF values F(yi | mu, tau) for each
# observation and posterior sample. It uses predict.brma to obtain the
# appropriate mu samples at the requested level.
#
# @param object brma object
# @param type   character; level at which to compute CDF. Options are:
#               - "marginal" (default): Fixed effects only (mu).
#                 CDF: yi ~ N(mu_i, tau^2 + sei^2)
#               - "study": Fixed effects + study-level random effects (mu + gamma).
#                 Only available for multilevel models.
#                 CDF: yi ~ N(mu_i + gamma_j, tau_within^2 + sei^2)
#               - "estimate": True study effects (mu + gamma + theta).
#                 CDF: yi ~ N(theta_i, sei^2)
#
# @return S x K matrix of CDF values in (0, 1)
#
# ---------------------------------------------------------------------------- #
.cdf.brma <- function(object, type = "marginal") {

  ### input validation
  type <- match.arg(type, c("marginal", "study", "estimate"))

  ### extract structural information about the model
  priors            <- object[["priors"]]
  data              <- object[["data"]]
  is_multilevel     <- .is_multilevel(object)
  is_scale          <- .is_scale(object)
  is_weightfunction <- .is_weightfunction(object)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)

  # check: study type requires multilevel model
  if (type == "study" && !is_multilevel) {
    stop("type = 'study' is only available for multilevel (3-level) models. ",
         "Use type = 'marginal' for non-multilevel models.", call. = FALSE)
  }

  ### obtain observed effect sizes and sampling SEs
  yi  <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)
  K   <- length(yi)

  ### get mu samples at the appropriate level using predict.brma
  # map CDF types to predict.brma types
  predict_type <- switch(type,
    "marginal"  = "terms",
    "study"     = "study",
    "estimate"  = "estimate"
  )

  mu_samples <- predict.brma(
    object  = object,
    newdata = NULL,
    type    = predict_type,
    quiet   = TRUE
  )
  S <- nrow(mu_samples)

  ### determine tau for CDF computation
  # for "estimate" type, tau = 0 (conditional on true effect, only sampling variance)
  # for "marginal" and "study" types, tau = tau_within (from scale model)
  if (type == "estimate") {
    tau_within_samples <- matrix(0, nrow = S, ncol = K)
  } else {
    # get tau samples from predict.brma
    tau_within_samples <- predict.brma(
      object  = object,
      newdata = NULL,
      type    = "terms.scale",
      quiet   = TRUE
    )
  }

  ### compute CDF based on outcome type
  if (outcome_type == "norm") {

    # flip for negative effect direction (applies to normal and weighted normal)
    if (effect_direction == "negative") {
      mu_samples_cdf <- -mu_samples
      yi_cdf         <- -yi
    } else {
      mu_samples_cdf <- mu_samples
      yi_cdf         <- yi
    }

    # dispatch between weighted and standard normal
    if (is_weightfunction) {

      # extract omega samples for weight function
      posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
      omega_samples     <- posterior_samples[, grep("omega", colnames(posterior_samples)), drop = FALSE]

      # get fit_data for crit_yi
      fit_data <- .create_fit_data(data = data, priors = priors)

      # compute use_normal indicator for performance optimization
      # this identifies which samples come from non-weightfunction bias models
      use_normal <- .extract_use_normal(object)

      cdf_vals <- .outcome_cdf.wnorm(
        yi         = yi_cdf,
        mu_samples = mu_samples_cdf,
        tau_within = tau_within_samples,
        sei        = sei,
        omega      = omega_samples,
        crit_yi    = fit_data$crit_yi,
        use_normal = use_normal
      )

    } else {

      # standard normal CDF
      cdf_vals <- .outcome_cdf.norm(
        yi         = yi_cdf,
        mu_samples = mu_samples_cdf,
        tau_within = tau_within_samples,
        sei        = sei
      )

    }

    # flip CDF for negative effect direction
    if (effect_direction == "negative") {
      cdf_vals <- 1 - cdf_vals
    }

  } else if (outcome_type == "bin") {

    # binomial CDF using normal approximation
    cdf_vals <- .outcome_cdf.binom(
      yi         = yi,
      sei        = sei,
      mu_samples = mu_samples
    )

  } else if (outcome_type == "pois") {

    # Poisson CDF using normal approximation
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


