# ============================================================================ #
# Outcome PDF Functions (Log-Likelihood)
# ============================================================================ #
#
# These functions compute pointwise log-likelihoods for each observation and
# posterior sample. They are used for LOO-PSIS diagnostics and model comparison.
#
# Parallels the structure of brma.rng.R but returns density values instead of
# random samples.
#
# Note: For binomial and Poisson models, each "observation" consists of a pair
# of data points (ai+ci or x1i+x2i) that together define a single effect size
# estimate. The log-likelihood for the observation is the sum of log-likelihoods
# for both components.
#
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
# @param mu_samples       S x K matrix of location samples (with study effects if multilevel)
# @param tau_within       S x K matrix of within-study heterogeneity samples
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
# .outcome_pdf.wnorm
# ---------------------------------------------------------------------------- #
#
# Compute pointwise log-likelihoods for weighted normal distribution (selection models).
#
# For selection models, the observed effect y_i follows a weighted normal:
#   y_i ~ f(y) * omega(y) / Z
# where f(y) is normal density and omega(y) is the selection weight function.
#
# This uses the extended weighted normal density function that handles
# matrix omega inputs (one row per posterior sample).
#
# @param yi               numeric vector of length K; observed effect sizes
# @param mu_samples       S x K matrix of location samples
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param sei              numeric vector of length K; standard errors
# @param omega            S x W matrix of omega (weight) samples
# @param crit_yi          W x K matrix of critical values for each observation
#
# @return S x K matrix of log-likelihood values from weighted distribution
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.wnorm <- function(yi, mu_samples, tau_within, sei, omega, crit_yi) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # pre-compute total SD for each observation
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  # compute log-likelihood for each observation using extended weighted normal
  # (loop is necessary due to observation-specific crit_yi values)
  log_lik <- matrix(NA_real_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    log_lik[, k] <- .dwnorm_fast.ss.matrix(
      x      = yi[k],
      mean   = mu_samples[, k],
      sd     = total_sd[, k],
      omega  = omega,
      crit_x = crit_yi[, k],
      log    = TRUE
    )
  }

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .outcome_pdf.binom
# ---------------------------------------------------------------------------- #
#
# Compute pointwise log-likelihoods for binomial outcome models.
#
# For binomial outcome models (log-odds ratio), the observed counts have likelihood:
#   ai ~ Binom(n1i, p1) where logit(p1) = logit(pi) + 0.5 * mu (treatment)
#   ci ~ Binom(n2i, p2) where logit(p2) = logit(pi) - 0.5 * mu (control)
#
# The log-likelihood for the observation is the sum of log-likelihoods for
# both groups since they jointly define the effect size estimate.
#
# @param ai               integer vector of length K; events in treatment group
# @param ci               integer vector of length K; events in control group
# @param n1i              integer vector of length K; treatment group sizes
# @param n2i              integer vector of length K; control group sizes
# @param mu_samples       S x K matrix of log-odds ratio samples
# @param logit_baserate   S x K matrix of logit(pi) base-rate samples
#
# @return S x K matrix of log-likelihood values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.binom <- function(ai, ci, n1i, n2i, mu_samples, logit_baserate) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # compute logit probabilities for each group
  # group 1: logit(p1) = logit(pi) + 0.5 * mu (treatment/exposed)
  # group 2: logit(p2) = logit(pi) - 0.5 * mu (control/unexposed)
  logit_p1 <- logit_baserate + 0.5 * mu_samples
  logit_p2 <- logit_baserate - 0.5 * mu_samples

  # convert to probabilities
  p1 <- .inv_logit(logit_p1)
  p2 <- .inv_logit(logit_p2)

  # replicate counts across samples
  ai_mat  <- matrix(ai,  nrow = S, ncol = K, byrow = TRUE)
  ci_mat  <- matrix(ci,  nrow = S, ncol = K, byrow = TRUE)
  n1i_mat <- matrix(n1i, nrow = S, ncol = K, byrow = TRUE)
  n2i_mat <- matrix(n2i, nrow = S, ncol = K, byrow = TRUE)

  # compute log-likelihood as sum of both binomial components
  log_lik_ai <- stats::dbinom(ai_mat, size = n1i_mat, prob = p1, log = TRUE)
  log_lik_ci <- stats::dbinom(ci_mat, size = n2i_mat, prob = p2, log = TRUE)

  log_lik <- log_lik_ai + log_lik_ci

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .outcome_pdf.pois
# ---------------------------------------------------------------------------- #
#
# Compute pointwise log-likelihoods for Poisson outcome models.
#
# For Poisson outcome models (log incidence rate ratio), the observed counts are:
#   x1i ~ Pois(t1i * exp(phi + 0.5 * mu))  (treatment/exposed)
#   x2i ~ Pois(t2i * exp(phi - 0.5 * mu))  (control/unexposed)
#
# The log-likelihood for the observation is the sum of log-likelihoods for
# both groups since they jointly define the effect size estimate.
#
# @param x1i              integer vector of length K; events in treatment group
# @param x2i              integer vector of length K; events in control group
# @param t1i              numeric vector of length K; treatment exposure times
# @param t2i              numeric vector of length K; control exposure times
# @param mu_samples       S x K matrix of log incidence rate ratio samples
# @param log_phi          S x K matrix of log baseline rate samples
#
# @return S x K matrix of log-likelihood values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.pois <- function(x1i, x2i, t1i, t2i, mu_samples, log_phi) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate exposure times to S x K matrices for vectorized ops
  log_t1i <- matrix(log(t1i), nrow = S, ncol = K, byrow = TRUE)
  log_t2i <- matrix(log(t2i), nrow = S, ncol = K, byrow = TRUE)

  # compute log rates for each group
  # group 1: log(lambda1) = phi + 0.5 * mu + log(t1i)
  # group 2: log(lambda2) = phi - 0.5 * mu + log(t2i)
  log_lambda1 <- log_phi + 0.5 * mu_samples + log_t1i
  log_lambda2 <- log_phi - 0.5 * mu_samples + log_t2i

  # replicate counts across samples
  x1i_mat <- matrix(x1i, nrow = S, ncol = K, byrow = TRUE)
  x2i_mat <- matrix(x2i, nrow = S, ncol = K, byrow = TRUE)

  # compute log-likelihood as sum of both Poisson components
  log_lik_x1i <- stats::dpois(x1i_mat, lambda = exp(log_lambda1), log = TRUE)
  log_lik_x2i <- stats::dpois(x2i_mat, lambda = exp(log_lambda2), log = TRUE)

  log_lik <- log_lik_x1i + log_lik_x2i

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .pdf.brma
# ---------------------------------------------------------------------------- #
#
# Compute the full log-likelihood matrix for a brma object.
#
# This is the main dispatcher function that coordinates extraction of posterior
# samples and computation of pointwise log-likelihoods using the appropriate
# PDF function for the outcome type.
#
# The LOO is computed at the estimate level: for multilevel models, we condition
# on the fitted study-level random effects (gamma); for selection models, we
# condition on the omega samples.
#
# @param object           brma object
#
# @return S x K matrix of log-likelihood values
#
# ---------------------------------------------------------------------------- #
.pdf.brma <- function(object) {

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
  # for likelihood evaluation, we do NOT adjust for publication bias
  # (we want to evaluate likelihood under the assumed model, which includes bias)
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
    bias_adjusted     = FALSE,  # include PET/PEESE terms in likelihood
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

  ### dispatch to appropriate PDF function based on outcome type
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

      log_lik <- .outcome_pdf.wnorm(
        yi         = yi,
        mu_samples = mu_samples,
        tau_within = tau_within_samples,
        sei        = outcome_data[["sei"]],
        omega      = omega_samples,
        crit_yi    = fit_data$crit_yi
      )

    } else {

      # standard normal likelihood
      log_lik <- .outcome_pdf.norm(
        yi         = yi,
        mu_samples = mu_samples,
        tau_within = tau_within_samples,
        sei        = sei
      )

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

    # get baserate samples
    logit_baserate <- .evaluate.brma.baserate(fit = object[["fit"]], K = K)

    log_lik <- .outcome_pdf.binom(
      ai             = outcome_data[["ai"]],
      ci             = outcome_data[["ci"]],
      n1i            = outcome_data[["n1i"]],
      n2i            = outcome_data[["n2i"]],
      mu_samples     = mu_samples,
      logit_baserate = logit_baserate
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

    # get log-rate samples
    log_phi <- .evaluate.brma.lograte(fit = object[["fit"]], K = K)

    log_lik <- .outcome_pdf.pois(
      x1i        = outcome_data[["x1i"]],
      x2i        = outcome_data[["x2i"]],
      t1i        = outcome_data[["t1i"]],
      t2i        = outcome_data[["t2i"]],
      mu_samples = mu_samples,
      log_phi    = log_phi
    )

  }

  # add column names for observations
  colnames(log_lik) <- paste0("log_lik[", seq_len(K), "]")

  return(log_lik)
}
