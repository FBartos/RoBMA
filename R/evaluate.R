# ============================================================================ #
# brma.evaluate.R
# ============================================================================ #
#
# This file contains modular helper functions for evaluating posterior samples
# from brma model fits. These functions extract and transform MCMC samples for:
# - heterogeneity (tau) parameters
# - location (mu) parameters
# - true study effects (theta)
# - observed responses
# - GLMM-specific parameters (baserate, lograte)
#
# The functions are designed for:
# - Prediction (predict.brma)
# - Simulation from the posterior predictive distribution
# - Likelihood evaluation
# - Other downstream tasks requiring posterior sample manipulation
#
# Design principles:
# - Decomposed parameters: minimize memory by passing only required components
# - Vectorized operations: use outer(), sweep(), matrix ops instead of loops
# - Consistent return structures: always return lists with predictable components
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .evaluate.brma.tau
# ---------------------------------------------------------------------------- #
#
# Extract and compute heterogeneity (tau) posterior samples from a brma fit.
#
# This function handles:
# - Scale regression models: evaluates log_tau formula, then exponentiates
# - Simple models: extracts tau column and replicates to K columns
# - Multilevel models: splits tau into within/between components via rho
#
# @param fit         runjags fit object containing posterior samples
# @param scale_data  data.frame of scale predictors for formula evaluation
#                    (NULL if not a scale regression model)
# @param scale_formula formula object for scale regression (from attr(data$scale, "formula"))
# @param scale_priors list of priors for scale parameters
# @param is_scale    logical; whether model uses scale regression
# @param is_multilevel logical; whether model is 3-level (study-level clustering)
# @param K           integer; number of observations (determines output columns)
#
# @return A list with two components (all S x K matrices):
#   - tau_within: estimate-level (within-study) heterogeneity component
#   - tau_between: study-level (between-study) heterogeneity component
#   For non-multilevel models: tau_within = tau (total), tau_between = 0 matrix
#   Total tau can be reconstructed as: sqrt(tau_within^2 + tau_between^2)
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.tau <- function(fit, scale_data, scale_formula, scale_priors,
                               is_scale, is_multilevel, K) {

  # extract posterior samples from JAGS fit
  # suppressWarnings: coda may warn about thinning or chain length

  posterior_samples <- suppressWarnings(coda::as.mcmc(fit))
  S <- nrow(posterior_samples)  # number of posterior samples

  ### compute tau samples based on model type
  if (is_scale) {

    # scale regression: evaluate log_tau formula then exponentiate
    # BayesTools::JAGS_evaluate_formula returns K x S matrix, we need S x K
    log_tau_samples <- t(BayesTools::JAGS_evaluate_formula(
      fit        = fit,
      formula    = scale_formula,
      parameter  = "log_tau",
      data       = scale_data,
      prior_list = scale_priors
    ))
    tau_samples <- exp(log_tau_samples)

  } else {

    # simple model: extract tau column and replicate to K columns
    # matrix(vec, nrow = S, ncol = K) replicates vec across columns
    # equivalent to: for (k in 1:K) result[, k] <- posterior_samples[, "tau"]
    tau_samples <- matrix(posterior_samples[, "tau"], nrow = S, ncol = K)

  }

  ### split tau into within/between components for multilevel models
  if (is_multilevel) {

    # extract rho (proportion of variance at estimate-level)
    rho_samples <- posterior_samples[, "rho"]

    # clamp rho to [0, 1] to handle JAGS numerical precision issues
    # pmin/pmax are vectorized min/max: faster than rho[rho > 1] <- 1
    rho_samples <- pmin(pmax(rho_samples, 0), 1)

    # tau_within = tau * sqrt(rho)       (estimate-level heterogeneity)
    # tau_between = tau * sqrt(1 - rho)  (study-level heterogeneity)
    # multiplication by vector rho_samples broadcasts across columns
    tau_within_samples  <- tau_samples * sqrt(rho_samples)
    tau_between_samples <- tau_samples * sqrt(1 - rho_samples)

  } else {

    # non-multilevel: all heterogeneity is at estimate-level
    # tau_between is zero matrix for consistent interface
    tau_within_samples  <- tau_samples
    tau_between_samples <- matrix(0, nrow = S, ncol = K)

  }

  return(list(
    tau_within  = tau_within_samples,
    tau_between = tau_between_samples
  ))
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.mu
# ---------------------------------------------------------------------------- #
#
# Extract and compute location (mu) posterior samples from a brma fit.
#
# This function handles:
# - Meta-regression models: evaluates mu formula with moderators
# - Simple models: extracts mu column and replicates to K columns
# - Effect direction flipping (for models fit with negative direction)
# - PET/PEESE bias adjustments (added when bias_adjusted = FALSE)
#
# @param fit              runjags fit object containing posterior samples
# @param outcome_data     data.frame with outcome info (must contain 'sei')
# @param mods_data        data.frame of moderators for formula evaluation
#                         (NULL if not a meta-regression model)
# @param mods_formula     formula object for meta-regression
# @param mods_priors      list of priors for moderator parameters
# @param is_mods          logical; whether model is meta-regression
# @param is_PET           logical; whether model includes PET adjustment
# @param is_PEESE         logical; whether model includes PEESE adjustment
# @param effect_direction character; "positive" or "negative" - direction of true effect
# @param bias_adjusted    logical; if TRUE, PET/PEESE adjustments are skipped
#                         (returns bias-adjusted mu); if FALSE, adds bias terms
# @param K                integer; number of observations
#
# @return S x K matrix of mu (location) posterior samples
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.mu <- function(fit, outcome_data, mods_data, mods_formula, mods_priors,
                              is_mods, is_PET, is_PEESE, effect_direction,
                              bias_adjusted, K) {

  # extract posterior samples from JAGS fit
  posterior_samples <- suppressWarnings(coda::as.mcmc(fit))
  S <- nrow(posterior_samples)

  ### compute base mu samples
  if (is_mods) {

    # meta-regression: evaluate mu formula with moderators
    # returns K x S, transpose to S x K
    mu_samples <- t(BayesTools::JAGS_evaluate_formula(
      fit        = fit,
      formula    = mods_formula,
      parameter  = "mu",
      data       = mods_data,
      prior_list = mods_priors
    ))

  } else {

    # simple model: replicate mu column to K columns
    mu_samples <- matrix(posterior_samples[, "mu"], nrow = S, ncol = K)

  }

  # NOTE: No effect direction flipping needed here!
  # The JAGS model uses: ifelse(effect_direction == "negative", "-mu", "mu") in the likelihood
  # This means mu in the posterior already represents the true effect in its original scale
  # (e.g., if true effect is -0.13, mu ≈ -0.13 in the posterior)
  # The data flip in .create_fit_data() is matched by -mu in the likelihood

  ### add PET adjustment when NOT incorporating publication bias adjustment
  # (i.e., when we want to show the biased predictions)
  # PET model in JAGS: yi_flipped ~ N(-mu + PET * sei, tau) for negative effects
  #                    yi ~ N(mu + PET * sei, tau) for positive effects
  # To get biased effect in original scale:
  #   - positive: E[yi] = mu + PET * sei
  #   - negative: E[yi_original] = -E[yi_flipped] = -(-mu + PET * sei) = mu - PET * sei
  # So the sign of PET adjustment depends on effect_direction
  if (is_PET && !bias_adjusted) {

    PET_samples <- posterior_samples[, "PET"]
    sei_vec     <- outcome_data[["sei"]]

    # vectorized: outer(PET_samples, sei_vec) creates S x K matrix
    # outer(a, b) computes a[i] * b[j] for all i,j pairs
    # direction multiplier: +1 for positive, -1 for negative
    direction <- ifelse(effect_direction == "negative", -1, 1)
    mu_samples <- mu_samples + direction * outer(PET_samples, sei_vec)

  }

  ### add PEESE adjustment when NOT incorporating publication bias adjustment
  # PEESE model: Same logic as PET but with sei^2
  if (is_PEESE && !bias_adjusted) {

    PEESE_samples <- posterior_samples[, "PEESE"]
    sei_sq_vec    <- outcome_data[["sei"]]^2

    # direction multiplier: +1 for positive, -1 for negative
    direction <- ifelse(effect_direction == "negative", -1, 1)
    mu_samples <- mu_samples + direction * outer(PEESE_samples, sei_sq_vec)

  }

  return(mu_samples)
}


# ---------------------------------------------------------------------------- #
# .get_multilevel_block_indices
# ---------------------------------------------------------------------------- #
#
# Create a reusable block index list for multilevel covariance calculations.
#
# @param study_ids integer vector of study identifiers
#
# @return A named list of integer index vectors, one per study.
#
# ---------------------------------------------------------------------------- #
.get_multilevel_block_indices <- function(study_ids) {

  split(seq_along(study_ids), study_ids)
}


# ---------------------------------------------------------------------------- #
# .build_multilevel_marginal_covariance
# ---------------------------------------------------------------------------- #
#
# Construct the marginal covariance matrix for a 3-level normal model.
#
# The covariance decomposes into:
# - sampling variance: diag(vi)
# - estimate-level heterogeneity: diag(tau_within^2)
# - study-level heterogeneity: block-wise tcrossprod(tau_between)
#
# @param tau_within    numeric vector of length K with within-study SDs
# @param tau_between   numeric vector of length K with between-study SDs
# @param vi            numeric vector of length K with sampling variances
# @param block_indices list of observation indices for each study
#
# @return A K x K marginal covariance matrix.
#
# ---------------------------------------------------------------------------- #
.build_multilevel_marginal_covariance <- function(tau_within, tau_between, vi,
                                                  block_indices) {

  K <- length(vi)

  marginal_covariance <- diag(vi + tau_within^2, nrow = K, ncol = K)

  for (idx in block_indices) {
    marginal_covariance[idx, idx] <- marginal_covariance[idx, idx] +
      tcrossprod(tau_between[idx])
  }

  return(marginal_covariance)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.multilevel_blup.norm
# ---------------------------------------------------------------------------- #
#
# Compute exact same-data BLUP contributions for 3-level normal models.
#
# The conditional predictor is obtained from the marginal mixed-model identity:
#   b_hat = G M^{-1} (y - X beta)
# where G is decomposed into study-level and estimate-level components.
#
# @param mu_samples   S x K matrix of fixed-effect predictions
# @param tau_within   S x K matrix of within-study SDs
# @param tau_between  S x K matrix of between-study SDs
# @param yi           numeric vector of observed outcomes
# @param vi           numeric vector of sampling variances
# @param study_ids    integer vector of study identifiers
#
# @return A list with S x K matrices `study` and `estimate`.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.multilevel_blup.norm <- function(mu_samples, tau_within, tau_between,
                                                yi, vi, study_ids) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  block_indices <- .get_multilevel_block_indices(study_ids)

  study_contribution    <- matrix(0, nrow = S, ncol = K)
  estimate_contribution <- matrix(0, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    marginal_covariance <- .build_multilevel_marginal_covariance(
      tau_within    = tau_within[s, ],
      tau_between   = tau_between[s, ],
      vi            = vi,
      block_indices = block_indices
    )

    chol_m <- try(chol(marginal_covariance), silent = TRUE)

    if (inherits(chol_m, "try-error")) {
      weights <- tryCatch(solve(marginal_covariance), error = function(e) MASS::ginv(marginal_covariance))
    } else {
      weights <- chol2inv(chol_m)
    }

    weighted_residual <- as.vector(weights %*% (yi - mu_samples[s, ]))

    estimate_contribution[s, ] <- tau_within[s, ]^2 * weighted_residual

    for (idx in block_indices) {
      study_scale <- sum(tau_between[s, idx] * weighted_residual[idx])
      study_contribution[s, idx] <- tau_between[s, idx] * study_scale
    }
  }

  return(list(
    study    = study_contribution,
    estimate = estimate_contribution
  ))
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.study_effects
# ---------------------------------------------------------------------------- #
#
# Extract or sample study-level (gamma) random effects for multilevel models.
#
# For multilevel models, gamma[study_id] represents the standardized study-level
# random effect (i.e., gamma ~ N(0, 1)). The actual contribution to mu is
# gamma * tau_between.
#
# This function handles:
# - Same data: extracts fitted gamma samples from posterior
# - New data: samples new gamma from N(0, 1) (marginalizes over study effects)
#
# @param fit              runjags fit object (needed to extract gamma if same_data)
# @param tau_between      S x K matrix of between-study heterogeneity samples
# @param study_ids        integer vector of length K; study ID for each observation
# @param same_data        logical; TRUE if predicting on original data (use fitted gamma)
# @param effect_direction character; "positive" or "negative" (currently unused but kept for interface)
#
# @return S x K matrix: contribution from study-level random effects
#         (gamma[study_ids[k]] * tau_between[,k])
#         Can be added directly to mu samples.
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.study_effects <- function(fit, tau_between, study_ids,
                                         same_data, effect_direction) {

  S <- nrow(tau_between)
  K <- ncol(tau_between)

  # NOTE: No direction flipping needed for study effects!
  # The JAGS model uses: "-gamma*tau_between" for negative effects, "+gamma*tau_between" for positive
  # But when converting to original scale:
  # E[yi_original] = -E[yi_flipped] = -(-mu - gamma*tau_between) = mu + gamma*tau_between
  # So the contribution to the original-scale effect is always +gamma*tau_between

  if (same_data) {

    # extract fitted gamma samples from posterior
    posterior_samples <- suppressWarnings(coda::as.mcmc(fit))
    n_studies         <- max(study_ids)

    # extract all gamma columns at once: S x n_studies matrix
    gamma_col_names <- paste0("gamma[", seq_len(n_studies), "]")
    gamma_samples   <- posterior_samples[, gamma_col_names, drop = FALSE]

    # gamma_samples[, study_ids] reorders columns to match observations
    # this is S x K after reordering, element-wise multiply with tau_between
    study_contribution <- gamma_samples[, study_ids] * tau_between

  } else {

    # new data: sample fresh gamma ~ N(0, 1) for each observation
    # each observation is treated as from a new study (marginalize over gamma)
    # matrix(rnorm(S*K), S, K) generates S x K matrix of standard normals
    gamma_new          <- matrix(stats::rnorm(S * K), nrow = S, ncol = K)
    study_contribution <- gamma_new * tau_between

  }

  return(study_contribution)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.true_effects.norm
# ---------------------------------------------------------------------------- #
#
# Compute posterior samples of true study effects (theta) for normal models.
# (applies to selection models too)
#
# For same_data = TRUE: Uses empirical Bayes shrinkage (BLUP) to estimate true effects:
#   theta_i = lambda_i * y_i + (1 - lambda_i) * mu_i
# where:
#   lambda_i = tau_within^2 / (tau_within^2 + se_i^2)
#
# For same_data = FALSE: Samples from the marginal distribution of true effects:
#   theta_i = mu_i + epsilon_i * tau_within_i, where epsilon_i ~ N(0, 1)
#
# IMPORTANT: For multilevel models, mu_samples must already include the
# study-level contribution (gamma * tau_between) before calling this function.
#
# @param mu_samples       S x K matrix of location samples (must include
#                         gamma * tau_between contribution for multilevel models)
# @param tau_within       S x K matrix of within-study (estimate-level) heterogeneity
# @param yi               numeric vector of length K; observed effect sizes (used only if same_data = TRUE)
# @param sei              numeric vector of length K; standard errors (used only if same_data = TRUE)
# @param same_data        logical; TRUE for BLUP estimates, FALSE for marginal sampling
#
# @return S x K matrix of true effect (theta) posterior samples
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.true_effects.norm <- function(mu_samples, tau_within, yi, sei, same_data) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  if (same_data) {

    # BLUP: empirical Bayes shrinkage estimates for existing observations
    # replicate yi and sei across samples for vectorized computation
    yi_mat  <- matrix(yi,  nrow = S, ncol = K, byrow = TRUE)
    sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)

    # compute shrinkage factor lambda (S x K matrix)
    # lambda = tau^2 / (tau^2 + se^2) ranges from 0 (strong shrinkage) to 1 (weak)
    tau_sq <- tau_within^2
    se_sq  <- sei_mat^2
    lambda <- tau_sq / (tau_sq + se_sq)

    # BLUP: weighted average of observed effect and model prediction
    # high tau -> lambda -> 1 -> trust data more
    # low tau -> lambda -> 0 -> trust model more
    true_effects_samples <- lambda * yi_mat + (1 - lambda) * mu_samples

  } else {

    # Marginal sampling: sample new theta from N(mu, tau_within)
    # epsilon ~ N(0, 1), then theta = mu + epsilon * tau_within
    epsilon <- matrix(stats::rnorm(S * K), nrow = S, ncol = K)
    true_effects_samples <- mu_samples + epsilon * tau_within

  }

  return(true_effects_samples)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.true_effects.glmm
# ---------------------------------------------------------------------------- #
#
# Compute posterior samples of true study effects (theta) for GLMM models.
#
# For GLMM models (binomial or Poisson), the estimate-level random effects
# (theta) are directly sampled in JAGS (not marginalized as in normal models).
# The true effect is:
#   true_effect_i = mu_i + theta_i * tau_within_i
#
# For same_data: extracts theta[k] from posterior samples
# For new_data: samples new theta ~ N(0, 1) to marginalize over random effects
#
# IMPORTANT: For multilevel models, mu_samples must already include the
# study-level contribution (gamma * tau_between) before calling this function.
# This is handled by .evaluate.brma.study_effects() in brma.predict.R.
#
# @param fit              runjags fit object (needed to extract theta if same_data = TRUE)
# @param mu_samples       S x K matrix of location samples (must include
#                         gamma * tau_between contribution for multilevel models)
# @param tau_within       S x K matrix of within-study (estimate-level) heterogeneity
# @param same_data        logical; TRUE if predicting on original data
# @param K                integer; number of observations
#
# @return S x K matrix of true effect (theta) posterior samples
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.true_effects.glmm <- function(fit, mu_samples, tau_within, same_data, K) {

  # add the estimate-level random effects (theta * tau_within) to mu
  theta_contribution <- .evaluate.brma.theta.glmm(
    fit        = fit,
    tau_within = tau_within,
    same_data  = same_data,
    K          = K
  )
  true_effects_samples <- mu_samples + theta_contribution

  return(true_effects_samples)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.baserate
# ---------------------------------------------------------------------------- #
#
# Extract base-rate (pi) samples for binomial GLMM models.
#
# For binomial outcomes, the base-rate pi[i] represents the probability of
# success in the control/reference condition. The log-odds ratio effect size
# is then applied relative to this baseline.
#
# @param fit              runjags fit object containing posterior samples
# @param K                integer; number of observations
#
# @return A matrix of logit(pi) samples for each observation
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.baserate <- function(fit, K) {

  posterior_samples <- suppressWarnings(coda::as.mcmc(fit))
  S <- nrow(posterior_samples)

  # extract pi[i] for each observation and compute logit transform
  logit_baserate <- matrix(NA_real_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    pi_k             <- posterior_samples[, paste0("pi[", k, "]")]
    logit_baserate[, k] <- .logit(pi_k)
  }

  return(logit_baserate)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.lograte
# ---------------------------------------------------------------------------- #
#
# Extract log-rate (phi) samples for Poisson GLMM models.
#
# For Poisson outcomes, phi[i] represents the log of the baseline event rate.
# The incidence rate ratio effect size is then applied relative to this baseline.
#
# @param fit              runjags fit object containing posterior samples
# @param K                integer; number of observations
#
# @return A matrix of log-rate (phi) samples for each observation
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.lograte <- function(fit, K) {

  posterior_samples <- suppressWarnings(coda::as.mcmc(fit))
  S <- nrow(posterior_samples)

  # extract phi[i] for each observation
  log_rate <- matrix(NA_real_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    log_rate[, k]  <- posterior_samples[, paste0("phi[", k, "]")]
  }

  return(log_rate)
}


# ---------------------------------------------------------------------------- #
# .evaluate.brma.theta.glmm
# ---------------------------------------------------------------------------- #
#
# Extract or sample estimate-level random effects (theta) for GLMM models.
#
# For GLMMs, theta[i] represents the standardized estimate-level random effect
# (i.e., theta ~ N(0, 1)). The actual random effect is theta * tau_within.
#
# @param fit              runjags fit object (needed to extract theta if same_data)
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param same_data        logical; TRUE if predicting on original data
# @param K                integer; number of observations
#
# @return S x K matrix: mu contribution from estimate-level random effects
#         (theta[k] * tau_within[,k])
#
# ---------------------------------------------------------------------------- #
.evaluate.brma.theta.glmm <- function(fit, tau_within, same_data, K) {

  S <- nrow(tau_within)

  if (same_data) {

    # extract fitted theta samples from posterior
    posterior_samples <- suppressWarnings(coda::as.mcmc(fit))

    theta_contribution <- matrix(NA_real_, nrow = S, ncol = K)
    for (k in seq_len(K)) {
      theta_k                 <- posterior_samples[, paste0("theta[", k, "]")]
      theta_contribution[, k] <- theta_k * tau_within[, k]
    }

  } else {

    # new data: sample fresh theta ~ N(0, 1) for each observation
    theta_new          <- matrix(stats::rnorm(S * K), nrow = S, ncol = K)
    theta_contribution <- theta_new * tau_within

  }

  return(theta_contribution)
}


.logit <- function(p) {
  return(log(p / (1 - p)))
}
.inv_logit <- function(x) {
  return(1 / (1 + exp(-x)))
}


# ---------------------------------------------------------------------------- #
# .extract_use_normal
# ---------------------------------------------------------------------------- #
#
# Extract bias indicator and compute use_normal logical vector.
#
# For RoBMA objects with mixture bias priors, `bias_indicator` tracks which
# bias model generated each posterior sample. For brma objects with a single
# bias prior, no indicator column exists (all samples use the same model).
#
# This function is used for performance optimization in PDF/CDF/RNG functions.
# When `use_normal[i] = TRUE`, the sample uses a non-weightfunction bias model
# (PET, PEESE, or no bias), so the fast normal path can be used instead of
# the expensive weighted normal computation.
#
# @param object    brma or RoBMA object
#
# @return A logical vector of length S (number of posterior samples):
#   - TRUE if sample uses non-weightfunction bias model (PET, PEESE, or no bias)
#   - FALSE if sample uses weightfunction bias model
#
# ---------------------------------------------------------------------------- #
.extract_use_normal <- function(object) {

  priors <- object[["priors"]]
  fit    <- object[["fit"]]

  # extract posterior samples
  posterior_samples <- suppressWarnings(coda::as.mcmc(fit))
  S <- nrow(posterior_samples)

  # check if bias_indicator column exists (RoBMA with mixture priors)
  has_bias_indicator <- "bias_indicator" %in% colnames(posterior_samples)

  if (has_bias_indicator) {

    # RoBMA case: multiple bias priors in mixture
    bias_indicator <- as.integer(posterior_samples[, "bias_indicator"])

    # extract bias priors and ensure list format
    priors_bias <- priors[["outcome"]][["bias"]]
    if (!BayesTools::is.prior.mixture(priors_bias)) {
      priors_bias <- list(priors_bias)
    }

    # identify which bias priors are weightfunctions
    weightfunction_indices <- which(sapply(priors_bias, BayesTools::is.prior.weightfunction))

    # use_normal = TRUE for samples NOT from weightfunction priors
    use_normal <- !(bias_indicator %in% weightfunction_indices)

  } else {

    # brma case: single bias prior (or no bias)
    # check if the single prior is a weightfunction
    is_weightfunction <- .is_priors_weightfunction(priors)

    # if weightfunction, all samples use weighted path; otherwise all use normal
    use_normal <- rep(!is_weightfunction, S)

  }

  return(use_normal)
}
