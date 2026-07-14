# ============================================================================ #
# GLMM Outcome Log-Likelihoods
# ============================================================================ #



# ---------------------------------------------------------------------------- #
# .outcome_pdf.binom
# ---------------------------------------------------------------------------- #
#
# Native wrapper for binomial GLMM marginal log-likelihoods.
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.binom <- function(ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi,
                               weights = NULL, n_theta = 15, n_pi = 30) {

  use_aghq   <- missing(n_theta) && missing(n_pi)
  prior_spec <- .glmm_aghq_prior_spec(prior_pi, "bin")
  if (use_aghq) {
    .glmm_aghq_require_default(prior_spec, "bin")
    out <- .glmm_binom_aghq(
      ai          = ai,
      ci          = ci,
      n1i         = n1i,
      n2i         = n2i,
      mu_samples  = mu_samples,
      tau_within  = tau_within,
      weights     = weights,
      prior_spec  = prior_spec
    )
    return(.glmm_aghq_value(out))
  }

  if (!.has_native_glmm("bin")) {
    return(.outcome_pdf.binom_r(
      ai         = ai,
      ci         = ci,
      n1i        = n1i,
      n2i        = n2i,
      mu_samples = mu_samples,
      tau_within = tau_within,
      prior_pi   = prior_pi,
      weights    = weights,
      n_theta    = n_theta,
      n_pi       = n_pi
    ))
  }

  gh      <- .gauss_hermite_nodes(n_theta)
  pi_grid <- .glmm_binom_logit_pi_grid(
    ai       = ai,
    ci       = ci,
    n1i      = n1i,
    n2i      = n2i,
    prior_pi = prior_pi,
    n_pi     = n_pi
  )
  weights_arg <- if (is.null(weights)) NULL else .native_numeric_vector(
    .glmm_likelihood_weights(weights, ncol(mu_samples))
  )

  return(.Call(
    "RoBMA_glmm_binom_marginal_loglik",
    .native_integer_vector(ai),
    .native_integer_vector(ci),
    .native_integer_vector(n1i),
    .native_integer_vector(n2i),
    .native_numeric_matrix(mu_samples),
    .native_numeric_matrix(tau_within),
    weights_arg,
    .native_numeric_vector(gh[["nodes"]]),
    .native_numeric_vector(gh[["log_weights"]]),
    .native_numeric_matrix(pi_grid[["grid"]]),
    .native_numeric_matrix(pi_grid[["log_weights"]]),
    PACKAGE = "RoBMA"
  ))
}



# ---------------------------------------------------------------------------- #
# .outcome_pdf.binom
# ---------------------------------------------------------------------------- #
#
# Compute pointwise marginal log-likelihoods for binomial outcome models.
#
# For binomial GLMM models, we need to compute the marginal likelihood by
# integrating over the latent parameters (theta, pi):
#   log p(data) = log ∫∫ p(ai, ci | theta, pi, mu) p(theta) p(pi) d(theta) d(pi)
#
# This function uses Gauss-Hermite quadrature for the theta dimension (N(0,1)
# prior) which requires far fewer points than quantile spacing for similar
# accuracy. The pi dimension uses Gauss-Legendre quadrature on the prior CDF
# scale, so the integration covers the full baserate prior support.
#
# The integration over samples is vectorized for efficiency: instead of looping
# over S samples, we compute all samples simultaneously using matrix operations.
#
# Performance optimizations:
# - Gauss-Hermite quadrature for theta (15 points vs 30 for same accuracy)
# - Pre-computed log-binomial coefficients (avoids repeated lgamma calls)
# - Direct log-likelihood calculation without dbinom overhead
# - Log-sum-exp trick for numerical stability
# - Minimized memory allocations in inner loop
#
# @param ai               integer vector of length K; events in treatment group
# @param ci               integer vector of length K; events in control group
# @param n1i              integer vector of length K; treatment group sizes
# @param n2i              integer vector of length K; control group sizes
# @param mu_samples       S x K matrix of log-odds ratio samples (without theta contribution)
# @param tau_within       S x K matrix of estimate-level heterogeneity samples
# @param prior_pi         BayesTools prior object for pi
# @param n_theta          integer; number of Gauss-Hermite points for theta (default: 15)
# @param n_pi             integer; number of grid points for pi dimension (default: 30)
#
# @return S x K matrix of marginal log-likelihood values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.binom_r <- function(ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi,
                                 weights = NULL, n_theta = 15, n_pi = 30) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)
  weights <- .glmm_likelihood_weights(weights, K)

  # initialize output matrix
  log_lik <- matrix(NA_real_, nrow = S, ncol = K)

  # --- 1. SETUP THETA GRID (Gauss-Hermite Quadrature) ---
  gh <- .gauss_hermite_nodes(n_theta)
  theta_grid        <- gh$nodes
  log_theta_weights <- gh$log_weights

  pi_quadrature <- .glmm_binom_logit_pi_grid(
    ai       = ai,
    ci       = ci,
    n1i      = n1i,
    n2i      = n2i,
    prior_pi = prior_pi,
    n_pi     = n_pi
  )

  # Expand theta grid for vectorized computation
  theta_expanded <- rep(theta_grid, times = n_pi)  # length G

  # --- 2. PRE-COMPUTE BINOMIAL COEFFICIENTS (Outside Loop) ---
  log_binom_coef_a <- lgamma(n1i + 1) - lgamma(ai + 1) - lgamma(n1i - ai + 1)
  log_binom_coef_c <- lgamma(n2i + 1) - lgamma(ci + 1) - lgamma(n2i - ci + 1)

  # --- 3. LOOP OVER OBSERVATIONS ---
  for (k in seq_len(K)) {

    lpi_grid       <- pi_quadrature[["grid"]][, k]
    log_pi_weights <- pi_quadrature[["log_weights"]][, k]

    # Create log-weight vector for all grid points
    log_weights <- as.vector(outer(log_theta_weights, log_pi_weights, "+"))

    # Expand pi grid
    lpi_expanded <- rep(lpi_grid, each = n_theta)

    # --- 4. VECTORIZED INTEGRATION ---
    effect_mat <- mu_samples[, k] + outer(tau_within[, k], theta_expanded)

    half_effect <- 0.5 * effect_mat
    logit_p1    <- sweep(half_effect, 2, lpi_expanded, "+")
    logit_p2    <- sweep(-half_effect, 2, lpi_expanded, "+")

    log_p1   <- plogis(logit_p1, log.p = TRUE)
    log_1mp1 <- plogis(logit_p1, lower.tail = FALSE, log.p = TRUE)
    log_p2   <- plogis(logit_p2, log.p = TRUE)
    log_1mp2 <- plogis(logit_p2, lower.tail = FALSE, log.p = TRUE)

    log_lik_grid <- log_binom_coef_a[k] + ai[k] * log_p1 + (n1i[k] - ai[k]) * log_1mp1 +
                    log_binom_coef_c[k] + ci[k] * log_p2 + (n2i[k] - ci[k]) * log_1mp2
    log_lik_grid <- weights[k] * log_lik_grid

    log_terms <- sweep(log_lik_grid, 2, log_weights, "+")
    log_lik[, k] <- .rowLogSumExps(log_terms)
  }

  return(log_lik)
}



# ---------------------------------------------------------------------------- #
# .outcome_pdf.binom_conditional
# ---------------------------------------------------------------------------- #
#
# Compute conditional log-likelihoods for binomial outcome models.
#
# This path is used by bridge sampling where pi and theta are sampled model
# parameters and should not be integrated out.
#
# @param ai             integer vector of length K; events in treatment group
# @param ci             integer vector of length K; events in control group
# @param n1i            integer vector of length K; treatment group sizes
# @param n2i            integer vector of length K; control group sizes
# @param mu_samples     S x K matrix of log-odds ratio samples
# @param logit_baserate S x K matrix of logit base-rate samples
#
# @return S x K matrix of conditional log-likelihood values
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.binom_conditional <- function(ai, ci, n1i, n2i, mu_samples,
                                           logit_baserate) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  logit_p1 <- logit_baserate + 0.5 * mu_samples
  logit_p2 <- logit_baserate - 0.5 * mu_samples

  log_lik <- matrix(
    stats::dbinom(
      x    = rep(ai, each = S),
      size = rep(n1i, each = S),
      prob = as.vector(.inv_logit(logit_p1)),
      log  = TRUE
    ) + stats::dbinom(
      x    = rep(ci, each = S),
      size = rep(n2i, each = S),
      prob = as.vector(.inv_logit(logit_p2)),
      log  = TRUE
    ),
    nrow = S,
    ncol = K
  )

  return(log_lik)
}


.outcome_pdf_sum.binom_conditional <- function(ai, ci, n1i, n2i,
                                               mu_samples, logit_baserate,
                                               weights = NULL) {

  if (.has_native_glmm_row_sum("bin", conditional = TRUE)) {
    weights_arg <- if (is.null(weights)) NULL else .native_numeric_vector(
      .glmm_likelihood_weights(weights, ncol(mu_samples))
    )
    return(.Call(
      "RoBMA_glmm_binom_conditional_loglik_sum",
      .native_integer_vector(ai),
      .native_integer_vector(ci),
      .native_integer_vector(n1i),
      .native_integer_vector(n2i),
      .native_numeric_matrix(mu_samples),
      .native_numeric_matrix(logit_baserate),
      weights_arg,
      PACKAGE = "RoBMA"
    ))
  }

  log_lik <- .outcome_pdf.binom_conditional(
    ai             = ai,
    ci             = ci,
    n1i            = n1i,
    n2i            = n2i,
    mu_samples     = mu_samples,
    logit_baserate = logit_baserate
  )

  return(rowSums(.apply_log_lik_weights(log_lik, weights)))
}


.outcome_pdf_sum.binom <- function(ai, ci, n1i, n2i, mu_samples, tau_within,
                                   prior_pi, weights = NULL,
                                   n_theta = 15, n_pi = 30) {

  use_aghq   <- missing(n_theta) && missing(n_pi)
  prior_spec <- .glmm_aghq_prior_spec(prior_pi, "bin")
  if (use_aghq) {
    .glmm_aghq_require_default(prior_spec, "bin", row_sum = TRUE)
    out <- .glmm_binom_aghq(
      ai          = ai,
      ci          = ci,
      n1i         = n1i,
      n2i         = n2i,
      mu_samples  = mu_samples,
      tau_within  = tau_within,
      weights     = weights,
      prior_spec  = prior_spec,
      row_sum     = TRUE
    )
    return(.glmm_aghq_value(out))
  }

  if (!.has_native_glmm_row_sum("bin")) {
    return(rowSums(.outcome_pdf.binom(
      ai         = ai,
      ci         = ci,
      n1i        = n1i,
      n2i        = n2i,
      mu_samples = mu_samples,
      tau_within = tau_within,
      prior_pi   = prior_pi,
      weights    = weights,
      n_theta    = n_theta,
      n_pi       = n_pi
    )))
  }

  gh      <- .gauss_hermite_nodes(n_theta)
  pi_grid <- .glmm_binom_logit_pi_grid(
    ai       = ai,
    ci       = ci,
    n1i      = n1i,
    n2i      = n2i,
    prior_pi = prior_pi,
    n_pi     = n_pi
  )
  weights_arg <- if (is.null(weights)) NULL else .native_numeric_vector(
    .glmm_likelihood_weights(weights, ncol(mu_samples))
  )

  return(.Call(
    "RoBMA_glmm_binom_marginal_loglik_row_sum",
    .native_integer_vector(ai),
    .native_integer_vector(ci),
    .native_integer_vector(n1i),
    .native_integer_vector(n2i),
    .native_numeric_matrix(mu_samples),
    .native_numeric_matrix(tau_within),
    weights_arg,
    .native_numeric_vector(gh[["nodes"]]),
    .native_numeric_vector(gh[["log_weights"]]),
    .native_numeric_matrix(pi_grid[["grid"]]),
    .native_numeric_matrix(pi_grid[["log_weights"]]),
    PACKAGE = "RoBMA"
  ))
}



# ---------------------------------------------------------------------------- #
# .outcome_pdf.pois
# ---------------------------------------------------------------------------- #
#
# Native wrapper for Poisson GLMM marginal log-likelihoods.
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.pois <- function(x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi,
                              weights = NULL, n_theta = 15, n_phi = 30) {

  use_aghq   <- missing(n_theta) && missing(n_phi)
  prior_spec <- .glmm_aghq_prior_spec(prior_phi, "pois")
  if (use_aghq) {
    .glmm_aghq_require_default(prior_spec, "pois")
    out <- .glmm_pois_aghq(
      x1i        = x1i,
      x2i        = x2i,
      t1i        = t1i,
      t2i        = t2i,
      mu_samples = mu_samples,
      tau_within = tau_within,
      weights    = weights,
      prior_spec = prior_spec
    )
    return(.glmm_aghq_value(out))
  }

  if (!.has_native_glmm("pois")) {
    return(.outcome_pdf.pois_r(
      x1i        = x1i,
      x2i        = x2i,
      t1i        = t1i,
      t2i        = t2i,
      mu_samples = mu_samples,
      tau_within = tau_within,
      prior_phi  = prior_phi,
      weights    = weights,
      n_theta    = n_theta,
      n_phi      = n_phi
    ))
  }

  gh       <- .gauss_hermite_nodes(n_theta)
  phi_grid <- .glmm_pois_log_phi_grid(
    x1i      = x1i,
    x2i      = x2i,
    t1i      = t1i,
    t2i      = t2i,
    prior_phi = prior_phi,
    n_phi    = n_phi
  )
  weights_arg <- if (is.null(weights)) NULL else .native_numeric_vector(
    .glmm_likelihood_weights(weights, ncol(mu_samples))
  )

  return(.Call(
    "RoBMA_glmm_pois_marginal_loglik",
    .native_integer_vector(x1i),
    .native_integer_vector(x2i),
    .native_numeric_vector(t1i),
    .native_numeric_vector(t2i),
    .native_numeric_matrix(mu_samples),
    .native_numeric_matrix(tau_within),
    weights_arg,
    .native_numeric_vector(gh[["nodes"]]),
    .native_numeric_vector(gh[["log_weights"]]),
    .native_numeric_matrix(phi_grid[["grid"]]),
    .native_numeric_matrix(phi_grid[["log_weights"]]),
    PACKAGE = "RoBMA"
  ))
}


.outcome_pdf_sum.pois <- function(x1i, x2i, t1i, t2i, mu_samples,
                                  tau_within, prior_phi, weights = NULL,
                                  n_theta = 15, n_phi = 30) {

  use_aghq   <- missing(n_theta) && missing(n_phi)
  prior_spec <- .glmm_aghq_prior_spec(prior_phi, "pois")
  if (use_aghq) {
    .glmm_aghq_require_default(prior_spec, "pois", row_sum = TRUE)
    out <- .glmm_pois_aghq(
      x1i        = x1i,
      x2i        = x2i,
      t1i        = t1i,
      t2i        = t2i,
      mu_samples = mu_samples,
      tau_within = tau_within,
      weights    = weights,
      prior_spec = prior_spec,
      row_sum    = TRUE
    )
    return(.glmm_aghq_value(out))
  }

  if (!.has_native_glmm_row_sum("pois")) {
    return(rowSums(.outcome_pdf.pois(
      x1i        = x1i,
      x2i        = x2i,
      t1i        = t1i,
      t2i        = t2i,
      mu_samples = mu_samples,
      tau_within = tau_within,
      prior_phi  = prior_phi,
      weights    = weights,
      n_theta    = n_theta,
      n_phi      = n_phi
    )))
  }

  gh       <- .gauss_hermite_nodes(n_theta)
  phi_grid <- .glmm_pois_log_phi_grid(
    x1i      = x1i,
    x2i      = x2i,
    t1i      = t1i,
    t2i      = t2i,
    prior_phi = prior_phi,
    n_phi    = n_phi
  )
  weights_arg <- if (is.null(weights)) NULL else .native_numeric_vector(
    .glmm_likelihood_weights(weights, ncol(mu_samples))
  )

  return(.Call(
    "RoBMA_glmm_pois_marginal_loglik_row_sum",
    .native_integer_vector(x1i),
    .native_integer_vector(x2i),
    .native_numeric_vector(t1i),
    .native_numeric_vector(t2i),
    .native_numeric_matrix(mu_samples),
    .native_numeric_matrix(tau_within),
    weights_arg,
    .native_numeric_vector(gh[["nodes"]]),
    .native_numeric_vector(gh[["log_weights"]]),
    .native_numeric_matrix(phi_grid[["grid"]]),
    .native_numeric_matrix(phi_grid[["log_weights"]]),
    PACKAGE = "RoBMA"
  ))
}



# ---------------------------------------------------------------------------- #
# .outcome_pdf.pois
# ---------------------------------------------------------------------------- #
#
# Compute pointwise marginal log-likelihoods for Poisson outcome models.
#
# For Poisson GLMM models, we need to compute the marginal likelihood by
# integrating over the latent parameters (theta, phi):
#   log p(data) = log ∫∫ p(x1i, x2i | theta, phi, mu) p(theta) p(phi) d(theta) d(phi)
#
# This function uses Gauss-Hermite quadrature for the theta dimension (N(0,1)
# prior) which requires far fewer points than quantile spacing for similar
# accuracy. The phi dimension uses Gauss-Legendre quadrature on the prior CDF
# scale, so the integration covers the full log-rate prior support.
#
# The integration over samples is vectorized for efficiency: instead of looping
# over S samples, we compute all samples simultaneously using matrix operations.
#
# Performance optimizations:
# - Gauss-Hermite quadrature for theta (15 points vs 30 for same accuracy)
# - Pre-computed log-factorial terms (avoids repeated lgamma calls)
# - Direct log-likelihood calculation without dpois overhead
# - Log-sum-exp trick for numerical stability
# - Minimized memory allocations in inner loop
#
# @param x1i              integer vector of length K; events in treatment group
# @param x2i              integer vector of length K; events in control group
# @param t1i              numeric vector of length K; treatment exposure times
# @param t2i              numeric vector of length K; control exposure times
# @param mu_samples       S x K matrix of log incidence rate ratio samples (without theta contribution)
# @param tau_within       S x K matrix of estimate-level heterogeneity samples
# @param prior_phi        BayesTools prior object for phi
# @param n_theta          integer; number of Gauss-Hermite points for theta (default: 15)
# @param n_phi            integer; number of grid points for phi dimension (default: 30)
#
# @return S x K matrix of marginal log-likelihood values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.pois_r <- function(x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi,
                                weights = NULL, n_theta = 15, n_phi = 30) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)
  weights <- .glmm_likelihood_weights(weights, K)

  # initialize output matrix
  log_lik <- matrix(NA_real_, nrow = S, ncol = K)

  # --- 1. SETUP THETA GRID (Gauss-Hermite Quadrature) ---
  gh <- .gauss_hermite_nodes(n_theta)
  theta_grid        <- gh$nodes
  log_theta_weights <- gh$log_weights

  phi_quadrature <- .glmm_pois_log_phi_grid(
    x1i       = x1i,
    x2i       = x2i,
    t1i       = t1i,
    t2i       = t2i,
    prior_phi = prior_phi,
    n_phi     = n_phi
  )

  # Expand theta grid for vectorized computation
  theta_expanded <- rep(theta_grid, times = n_phi)  # length G

  # --- 2. PRE-COMPUTE LOG-FACTORIAL TERMS (Outside Loop) ---
  log_factorial_x1 <- lgamma(x1i + 1)
  log_factorial_x2 <- lgamma(x2i + 1)
  log_t1i <- log(t1i)
  log_t2i <- log(t2i)

  # --- 3. LOOP OVER OBSERVATIONS ---
  for (k in seq_len(K)) {

    phi_grid        <- phi_quadrature[["grid"]][, k]
    log_phi_weights <- phi_quadrature[["log_weights"]][, k]

    # Create log-weight vector for all grid points
    log_weights <- as.vector(outer(log_theta_weights, log_phi_weights, "+"))

    # Expand phi grid
    phi_expanded <- rep(phi_grid, each = n_theta)

    # --- 4. VECTORIZED INTEGRATION ---
    effect_mat <- mu_samples[, k] + outer(tau_within[, k], theta_expanded)

    half_effect <- 0.5 * effect_mat
    log_lambda1 <- sweep(half_effect, 2, phi_expanded + log_t1i[k], "+")
    log_lambda2 <- sweep(-half_effect, 2, phi_expanded + log_t2i[k], "+")

    log_lik_grid <- x1i[k] * log_lambda1 - exp(log_lambda1) - log_factorial_x1[k] +
                    x2i[k] * log_lambda2 - exp(log_lambda2) - log_factorial_x2[k]
    log_lik_grid <- weights[k] * log_lik_grid

    log_terms <- sweep(log_lik_grid, 2, log_weights, "+")
    log_lik[, k] <- .rowLogSumExps(log_terms)
  }

  return(log_lik)
}



# ---------------------------------------------------------------------------- #
# .outcome_pdf.pois_conditional
# ---------------------------------------------------------------------------- #
#
# Compute conditional log-likelihoods for Poisson outcome models.
#
# This path is used by bridge sampling where phi and theta are sampled model
# parameters and should not be integrated out.
#
# @param x1i        integer vector of length K; treatment events
# @param x2i        integer vector of length K; control events
# @param t1i        numeric vector of length K; treatment exposure times
# @param t2i        numeric vector of length K; control exposure times
# @param mu_samples S x K matrix of log incidence rate ratio samples
# @param log_phi    S x K matrix of midpoint log-rate samples
#
# @return S x K matrix of conditional log-likelihood values
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.pois_conditional <- function(x1i, x2i, t1i, t2i, mu_samples,
                                          log_phi) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  log_t1i <- matrix(log(t1i), nrow = S, ncol = K, byrow = TRUE)
  log_t2i <- matrix(log(t2i), nrow = S, ncol = K, byrow = TRUE)

  log_lambda1 <- log_phi + 0.5 * mu_samples + log_t1i
  log_lambda2 <- log_phi - 0.5 * mu_samples + log_t2i

  log_lik <- matrix(
    stats::dpois(
      x      = rep(x1i, each = S),
      lambda = as.vector(exp(log_lambda1)),
      log    = TRUE
    ) + stats::dpois(
      x      = rep(x2i, each = S),
      lambda = as.vector(exp(log_lambda2)),
      log    = TRUE
    ),
    nrow = S,
    ncol = K
  )

  return(log_lik)
}


.outcome_pdf_sum.pois_conditional <- function(x1i, x2i, t1i, t2i,
                                              mu_samples, log_phi,
                                              weights = NULL) {

  if (.has_native_glmm_row_sum("pois", conditional = TRUE)) {
    weights_arg <- if (is.null(weights)) NULL else .native_numeric_vector(
      .glmm_likelihood_weights(weights, ncol(mu_samples))
    )
    return(.Call(
      "RoBMA_glmm_pois_conditional_loglik_sum",
      .native_integer_vector(x1i),
      .native_integer_vector(x2i),
      .native_numeric_vector(t1i),
      .native_numeric_vector(t2i),
      .native_numeric_matrix(mu_samples),
      .native_numeric_matrix(log_phi),
      weights_arg,
      PACKAGE = "RoBMA"
    ))
  }

  log_lik <- .outcome_pdf.pois_conditional(
    x1i        = x1i,
    x2i        = x2i,
    t1i        = t1i,
    t2i        = t2i,
    mu_samples = mu_samples,
    log_phi    = log_phi
  )

  return(rowSums(.apply_log_lik_weights(log_lik, weights)))
}
