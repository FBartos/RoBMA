# ============================================================================ #
# Outcome RNG Functions
# ============================================================================ #
#
# These are prediction-specific RNG functions for sampling from the posterior
# predictive distribution of observed effect sizes. They are kept separate from
# the .evaluate.brma.* functions because they involve random number generation
# and are specific to the prediction task.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .outcome_rng.norm
# ---------------------------------------------------------------------------- #
#
# Sample observed effect sizes from normal posterior predictive distribution.
#
# For normal outcome models, the observed effect y_i is:
#   y_i ~ N(mu_i, tau_within^2 + se_i^2)
#
# This function samples from this distribution for each posterior sample.
#
# @param mu_samples       S x K matrix of location samples (with study effects if multilevel)
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param sei              numeric vector of length K; standard errors
#
# @return S x K matrix of sampled observed effect sizes
#
# ---------------------------------------------------------------------------- #
.outcome_rng.norm <- function(mu_samples, tau_within, sei) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # replicate sei across samples: matrix(sei, S, K, byrow = TRUE)
  # each row contains the same se values
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)

  # compute total SD: sqrt(tau^2 + se^2)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  # sample from N(mu, total_sd) for each cell
  # matrix(rnorm(S*K), S, K) * total_sd + mu is vectorized sampling
  # equivalent: for each s,k: y[s,k] ~ N(mu[s,k], total_sd[s,k])
  response_samples <- mu_samples + matrix(stats::rnorm(S * K), nrow = S, ncol = K) * total_sd

  return(response_samples)
}


# ---------------------------------------------------------------------------- #
# .outcome_rng.wnorm
# ---------------------------------------------------------------------------- #
#
# Sample observed effect sizes from weighted normal distribution (selection models).
#
# For selection models, the observed effect y_i follows a weighted normal:
#   y_i ~ f(y) * omega(y) / Z
# where f(y) is normal density and omega(y) is the selection weight function.
#
# This uses the fast spike-and-slab weighted normal sampler.
#
# @param mu_samples       S x K matrix of location samples
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param sei              numeric vector of length K; standard errors
# @param omega            S x W matrix of omega (weight) samples
# @param crit_yi          W x K matrix of critical values for each observation
#
# @return S x K matrix of sampled observed effect sizes from weighted distribution
#
# ---------------------------------------------------------------------------- #
.outcome_rng.wnorm <- function(mu_samples, tau_within, sei, omega, crit_yi) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # pre-compute total SD for each observation
  sei_mat  <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  # sample from weighted normal for each observation
  # (loop is necessary due to .rwnorm_fast.ss interface)
  response_samples <- matrix(NA_real_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    response_samples[, k] <- .rwnorm_fast.ss(
      mean   = mu_samples[, k],
      sd     = total_sd[, k],
      omega  = omega,
      crit_x = crit_yi[, k]
    )
  }

  return(response_samples)
}


# ---------------------------------------------------------------------------- #
# .outcome_rng.binom
# ---------------------------------------------------------------------------- #
#
# Sample observed counts from binomial posterior predictive distribution.
#
# For binomial outcome models (log-odds ratio), the observed counts are:
#   ai ~ Binom(n1i, p1) where logit(p1) = logit(pi) + 0.5 * mu (treatment)
#   ci ~ Binom(n2i, p2) where logit(p2) = logit(pi) - 0.5 * mu (control)
#
# The 0.5 factor arises from the log-OR parameterization where mu represents
# the log-odds ratio between treatment and control groups.
#
# @param mu_samples       S x K matrix of log-odds ratio samples
# @param logit_baserate   S x K matrix of logit(pi) base-rate samples
# @param n1i              numeric vector of length K; treatment group sizes
# @param n2i              numeric vector of length K; control group sizes
#
# @return S x (2*K) matrix with columns interleaved: ai[1], ci[1], ai[2], ci[2], ...
#
# ---------------------------------------------------------------------------- #
.outcome_rng.binom <- function(mu_samples, logit_baserate, n1i, n2i) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  # compute logit probabilities for each group
  # group 1: logit(p1) = logit(pi) + 0.5 * mu (treatment/exposed)
  # group 2: logit(p2) = logit(pi) - 0.5 * mu (control/unexposed)
  logit_p1 <- logit_baserate + 0.5 * mu_samples
  logit_p2 <- logit_baserate - 0.5 * mu_samples

  # sample from binomial for each group
  outcome_samples_ai <- matrix(NA_integer_, nrow = S, ncol = K)
  outcome_samples_ci <- matrix(NA_integer_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    outcome_samples_ai[, k] <- stats::rbinom(n = S, size = n1i[k], prob = .inv_logit(logit_p1[, k]))
    outcome_samples_ci[, k] <- stats::rbinom(n = S, size = n2i[k], prob = .inv_logit(logit_p2[, k]))
  }

  # name columns
  colnames(outcome_samples_ai) <- paste0("ai[", seq_len(K), "]")
  colnames(outcome_samples_ci) <- paste0("ci[", seq_len(K), "]")

  # merge samples with interleaved columns: ai[1], ci[1], ai[2], ci[2], ...
  outcome_samples <- matrix(NA_integer_, nrow = S, ncol = 2 * K)
  outcome_samples[, seq(1, 2 * K, by = 2)] <- outcome_samples_ai
  outcome_samples[, seq(2, 2 * K, by = 2)] <- outcome_samples_ci
  colnames(outcome_samples) <- as.vector(rbind(
    colnames(outcome_samples_ai),
    colnames(outcome_samples_ci)
  ))

  return(outcome_samples)
}


# ---------------------------------------------------------------------------- #
# .outcome_rng.pois
# ---------------------------------------------------------------------------- #
#
# Sample observed counts from Poisson posterior predictive distribution.
#
# For Poisson outcome models (log incidence rate ratio), the observed counts are:
#   x1i ~ Pois(t1i * exp(phi + 0.5 * mu))  (treatment/exposed)
#   x2i ~ Pois(t2i * exp(phi - 0.5 * mu))  (control/unexposed)
#
# where phi is the log baseline rate and t1i, t2i are exposure times.
#
# @param mu_samples       S x K matrix of log incidence rate ratio samples
# @param log_phi          S x K matrix of log baseline rate samples
# @param t1i              numeric vector of length K; treatment exposure times
# @param t2i              numeric vector of length K; control exposure times
#
# @return S x (2*K) matrix with columns interleaved: x1i[1], x2i[1], x1i[2], x2i[2], ...
#
# ---------------------------------------------------------------------------- #
.outcome_rng.pois <- function(mu_samples, log_phi, t1i, t2i) {

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

  # sample from Poisson for each group
  outcome_samples_x1i <- matrix(NA_integer_, nrow = S, ncol = K)
  outcome_samples_x2i <- matrix(NA_integer_, nrow = S, ncol = K)

  for (k in seq_len(K)) {
    outcome_samples_x1i[, k] <- stats::rpois(n = S, lambda = exp(log_lambda1[, k]))
    outcome_samples_x2i[, k] <- stats::rpois(n = S, lambda = exp(log_lambda2[, k]))
  }

  # name columns
  colnames(outcome_samples_x1i) <- paste0("x1i[", seq_len(K), "]")
  colnames(outcome_samples_x2i) <- paste0("x2i[", seq_len(K), "]")

  # merge samples with interleaved columns: x1i[1], x2i[1], x1i[2], x2i[2], ...
  outcome_samples <- matrix(NA_integer_, nrow = S, ncol = 2 * K)
  outcome_samples[, seq(1, 2 * K, by = 2)] <- outcome_samples_x1i
  outcome_samples[, seq(2, 2 * K, by = 2)] <- outcome_samples_x2i
  colnames(outcome_samples) <- as.vector(rbind(
    colnames(outcome_samples_x1i),
    colnames(outcome_samples_x2i)
  ))

  return(outcome_samples)
}