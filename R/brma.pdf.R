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
# .rowLogSumExps
# ---------------------------------------------------------------------------- #
#
# Compute log(sum(exp(x))) for each row of a matrix, using the log-sum-exp
# trick for numerical stability.
#
# @param lx  numeric matrix of log-values
#
# @return numeric vector of length nrow(lx)
#
# ---------------------------------------------------------------------------- #
.rowLogSumExps <- function(lx) {
  row_max <- do.call(pmax, c(as.data.frame(lx), na.rm = TRUE))
  row_max + log(rowSums(exp(lx - row_max)))
}


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
# Compute pointwise marginal log-likelihoods for binomial outcome models.
#
# For binomial GLMM models, we need to compute the marginal likelihood by
# integrating over the latent parameters (theta, pi):
#   log p(data) = log ∫∫ p(ai, ci | theta, pi, mu) p(theta) p(pi) d(theta) d(pi)
#
# This function uses Gauss-Hermite quadrature for the theta dimension (N(0,1)
# prior) which requires far fewer points than quantile spacing for similar
# accuracy. The pi dimension uses adaptive data-guided linear spacing.
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
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param prior_pi         BayesTools prior object for pi
# @param n_theta          integer; number of Gauss-Hermite points for theta (default: 15)
# @param n_pi             integer; number of grid points for pi dimension (default: 30)
#
# @return S x K matrix of marginal log-likelihood values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.binom <- function(ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi,
                                n_theta = 15, n_pi = 30) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)
  G <- n_theta * n_pi  # total number of grid points

  # initialize output matrix
  log_lik <- matrix(NA_real_, nrow = S, ncol = K)

  # --- 1. SETUP THETA GRID (Gauss-Hermite Quadrature) ---
  gh <- .gauss_hermite_nodes(n_theta)
  theta_grid        <- gh$nodes
  log_theta_weights <- gh$log_weights

  # Expand theta grid for vectorized computation
  theta_expanded <- rep(theta_grid, times = n_pi)  # length G

  # --- 2. PRE-COMPUTE BINOMIAL COEFFICIENTS (Outside Loop) ---
  log_binom_coef_a <- lgamma(n1i + 1) - lgamma(ai + 1) - lgamma(n1i - ai + 1)
  log_binom_coef_c <- lgamma(n2i + 1) - lgamma(ci + 1) - lgamma(n2i - ci + 1)

  # --- 3. LOOP OVER OBSERVATIONS (Data-Guided Centering for Pi) ---
  for (k in seq_len(K)) {

    # --- ADAPTIVE STEP: Center Pi Grid on Empirical Data ---
    p_hat <- if (ai[k] == 0 || ci[k] == 0) {
      (ai[k] + ci[k] + 0.5) / (n1i[k] + n2i[k] + 1)
    } else {
      (ai[k] + ci[k]) / (n1i[k] + n2i[k])
    }
    logit_center <- qlogis(p_hat)

    # Define range around center (+/- 5 on logit scale)
    lpi_grid <- seq(logit_center - 5, logit_center + 5, length.out = n_pi)
    pi_grid  <- plogis(lpi_grid)

    # Calculate log-weights with Jacobian correction
    lpi_step       <- lpi_grid[2] - lpi_grid[1]
    log_jacobian   <- stats::dlogis(lpi_grid, log = TRUE)
    log_prior      <- BayesTools::lpdf(prior_pi, pi_grid)
    log_pi_weights <- log_prior + log_jacobian + log(lpi_step)

    # Create log-weight vector for all G grid points
    log_weights <- as.vector(outer(log_theta_weights, log_pi_weights, "+"))  # length G

    # Expand pi grid
    lpi_expanded <- rep(lpi_grid, each = n_theta)  # length G

    # --- 4. VECTORIZED INTEGRATION ---
    effect_mat <- mu_samples[, k] + outer(tau_within[, k], theta_expanded)  # S x G

    half_effect <- 0.5 * effect_mat
    logit_p1    <- sweep(half_effect, 2, lpi_expanded, "+")
    logit_p2    <- sweep(-half_effect, 2, lpi_expanded, "+")

    log_p1   <- plogis(logit_p1, log.p = TRUE)
    log_1mp1 <- plogis(logit_p1, lower.tail = FALSE, log.p = TRUE)
    log_p2   <- plogis(logit_p2, log.p = TRUE)
    log_1mp2 <- plogis(logit_p2, lower.tail = FALSE, log.p = TRUE)

    log_lik_grid <- log_binom_coef_a[k] + ai[k] * log_p1 + (n1i[k] - ai[k]) * log_1mp1 +
                    log_binom_coef_c[k] + ci[k] * log_p2 + (n2i[k] - ci[k]) * log_1mp2

    log_terms <- sweep(log_lik_grid, 2, log_weights, "+")
    log_lik[, k] <- .rowLogSumExps(log_terms)
  }

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .gauss_hermite_nodes
# ---------------------------------------------------------------------------- #
#
# Compute Gauss-Hermite quadrature nodes and weights for standard normal.
#
# Gauss-Hermite quadrature integrates ∫ f(x) exp(-x²) dx exactly for polynomials
# up to degree 2n-1. For the standard normal N(0,1), we transform:
#   ∫ f(z) * (1/√(2π)) * exp(-z²/2) dz
#
# Using substitution u = z/√2, this becomes:
#   (1/√π) ∫ f(√2 * u) * exp(-u²) du
#
# So for N(0,1): nodes = √2 * GH_nodes, and we need to normalize weights
# so they sum to 1 (representing probability mass).
#
# This function uses the Golub-Welsch algorithm (eigenvalue method).
#
# @param n integer; number of quadrature points
#
# @return list with components:
#   - nodes: n quadrature points (transformed for N(0,1))
#   - weights: n quadrature weights (sum to 1 for N(0,1))
#   - log_weights: log of weights for numerical stability
#
# ---------------------------------------------------------------------------- #
.gauss_hermite_nodes <- function(n) {

  # Golub-Welsch algorithm: eigenvalues of symmetric tridiagonal matrix

  # For Gauss-Hermite (physicist's), the recurrence coefficients are:
  #   α_k = 0, β_k = k/2 for k = 1, ..., n-1
  # The Jacobi matrix has off-diagonal elements √β_k = √(k/2)

  i <- seq_len(n - 1)
  b <- sqrt(i / 2)

  # Build symmetric tridiagonal matrix
  H <- diag(0, n)
  H[cbind(i, i + 1)] <- b
  H[cbind(i + 1, i)] <- b

  # Eigendecomposition
  eig <- eigen(H, symmetric = TRUE)

  # Nodes are eigenvalues
  # Weights: w_i = μ₀ * (first component of eigenvector i)²
  # where μ₀ = ∫ exp(-x²) dx = √π
  nodes_raw   <- eig$values
  weights_raw <- sqrt(pi) * eig$vectors[1, ]^2

  # Transform for standard normal N(0,1):
  # ∫ f(z) φ(z) dz = (1/√π) ∫ f(√2 u) exp(-u²) du ≈ (1/√π) Σ w_i f(√2 u_i)
  # So: nodes = √2 * raw_nodes, weights = raw_weights / √π
  nodes   <- sqrt(2) * nodes_raw
  weights <- weights_raw / sqrt(pi)

  # Verify: weights should sum to 1
  # (The raw weights sum to √π because μ₀ = √π and Σ v_{1i}² = 1 for orthonormal eigenvectors)

  # Sort by nodes (eigen returns in decreasing order of eigenvalue)
  ord <- order(nodes)

  list(
    nodes       = nodes[ord],
    weights     = weights[ord],
    log_weights = log(weights[ord])
  )
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
# accuracy. The phi dimension uses adaptive data-guided linear spacing.
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
# @param tau_within       S x K matrix of within-study heterogeneity samples
# @param prior_phi        BayesTools prior object for phi
# @param n_theta          integer; number of Gauss-Hermite points for theta (default: 15)
# @param n_phi            integer; number of grid points for phi dimension (default: 30)
#
# @return S x K matrix of marginal log-likelihood values (one per estimate)
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.pois <- function(x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi,
                               n_theta = 15, n_phi = 30) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)
  G <- n_theta * n_phi  # total number of grid points

  # initialize output matrix
  log_lik <- matrix(NA_real_, nrow = S, ncol = K)

  # --- 1. SETUP THETA GRID (Gauss-Hermite Quadrature) ---
  gh <- .gauss_hermite_nodes(n_theta)
  theta_grid        <- gh$nodes
  log_theta_weights <- gh$log_weights

  # Expand theta grid for vectorized computation
  theta_expanded <- rep(theta_grid, times = n_phi)  # length G

  # --- 2. PRE-COMPUTE LOG-FACTORIAL TERMS (Outside Loop) ---
  log_factorial_x1 <- lgamma(x1i + 1)
  log_factorial_x2 <- lgamma(x2i + 1)
  log_t1i <- log(t1i)
  log_t2i <- log(t2i)

  # --- 3. LOOP OVER OBSERVATIONS (Data-Guided Centering for Phi) ---
  for (k in seq_len(K)) {

    # --- ADAPTIVE STEP: Center Phi Grid on Empirical Data ---
    if (x1i[k] == 0 || x2i[k] == 0) {
      empirical_rate <- (x1i[k] + x2i[k] + 0.5) / (t1i[k] + t2i[k])
    } else {
      empirical_rate <- (x1i[k] + x2i[k]) / (t1i[k] + t2i[k])
    }
    log_rate_center <- log(empirical_rate)

    # Define range around center (+/- 5 on log scale)
    phi_grid <- seq(log_rate_center - 5, log_rate_center + 5, length.out = n_phi)

    # Calculate log-weights (no Jacobian needed - phi is on log scale)
    phi_step        <- phi_grid[2] - phi_grid[1]
    log_phi_weights <- BayesTools::lpdf(prior_phi, phi_grid) + log(phi_step)

    # Create log-weight vector for all G grid points
    log_weights <- as.vector(outer(log_theta_weights, log_phi_weights, "+"))  # length G

    # Expand phi grid
    phi_expanded <- rep(phi_grid, each = n_theta)  # length G

    # --- 4. VECTORIZED INTEGRATION ---
    effect_mat <- mu_samples[, k] + outer(tau_within[, k], theta_expanded)  # S x G

    half_effect <- 0.5 * effect_mat
    log_lambda1 <- sweep(half_effect, 2, phi_expanded + log_t1i[k], "+")
    log_lambda2 <- sweep(-half_effect, 2, phi_expanded + log_t2i[k], "+")

    log_lik_grid <- x1i[k] * log_lambda1 - exp(log_lambda1) - log_factorial_x1[k] +
                    x2i[k] * log_lambda2 - exp(log_lambda2) - log_factorial_x2[k]

    log_terms <- sweep(log_lik_grid, 2, log_weights, "+")
    log_lik[, k] <- .rowLogSumExps(log_terms)
  }

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .pdf.brma
# ---------------------------------------------------------------------------- #
#
# Compute the full log-likelihood matrix for a brma object.
#
# This function computes pointwise log-likelihoods f(yi | mu, tau) for each
# observation and posterior sample. It uses predict.brma to obtain the
# appropriate mu samples.
#
# For GLMMs, the log-likelihood includes p(theta_z) where theta_z ~ N(0,1) is
# the standardized estimate-level random effect. This is needed for proper LOO
# computation since theta_z is a latent variable in the model.
#
# @param object brma object
#
# @return S x K matrix of log-likelihood values
#
# ---------------------------------------------------------------------------- #
.pdf.brma <- function(object) {

  ### extract structural information about the model
  priors            <- object[["priors"]]
  data              <- object[["data"]]
  is_multilevel     <- .is_multilevel(object)
  is_weightfunction <- .is_weightfunction(object)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)

  ### obtain observed effect sizes and sampling SEs
  yi  <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)
  K   <- length(yi)

  ### get mu samples at appropriate level using predict.brma
  # for likelihood, we use mu + gamma (study level for multilevel), terms for non-multilevel
  # the theta contribution is handled separately for GLMMs
  mu_type <- if (is_multilevel) "study" else "terms"
  mu_samples <- predict.brma(
    object  = object,
    newdata = NULL,
    type    = mu_type,
    quiet   = TRUE
  )
  S <- nrow(mu_samples)

  ### get tau samples from predict.brma
  tau_within_samples <- predict.brma(
    object  = object,
    newdata = NULL,
    type    = "terms.scale",
    quiet   = TRUE
  )

  ### extract posterior samples once if needed
  # (only for weightfunction or GLMMs)
  if (is_weightfunction || outcome_type %in% c("bin", "pois")) {
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
  }

  ### compute log-likelihood based on outcome type
  if (outcome_type == "norm") {

    # flip for negative effect direction (applies to normal and weighted normal)
    if (effect_direction == "negative") {
      mu_samples_pdf <- -mu_samples
      yi_pdf         <- -yi
    } else {
      mu_samples_pdf <- mu_samples
      yi_pdf         <- yi
    }

    # dispatch between weighted and standard normal
    if (is_weightfunction) {

      # extract omega samples for weight function
      omega_samples <- posterior_samples[, grep("omega", colnames(posterior_samples)), drop = FALSE]

      # get fit_data for crit_yi
      fit_data <- .create_fit_data(data = data, priors = priors)

      log_lik <- .outcome_pdf.wnorm(
        yi         = yi_pdf,
        mu_samples = mu_samples_pdf,
        tau_within = tau_within_samples,
        sei        = sei,
        omega      = omega_samples,
        crit_yi    = fit_data$crit_yi
      )

    } else {

      # standard normal likelihood
      log_lik <- .outcome_pdf.norm(
        yi         = yi_pdf,
        mu_samples = mu_samples_pdf,
        tau_within = tau_within_samples,
        sei        = sei
      )

    }

  } else if (outcome_type == "bin") {

    # compute marginal log-likelihood via integration over theta and pi
    # the new .outcome_pdf.binom handles integration internally
    log_lik <- .outcome_pdf.binom(
      ai         = data[["outcome"]][["ai"]],
      ci         = data[["outcome"]][["ci"]],
      n1i        = data[["outcome"]][["n1i"]],
      n2i        = data[["outcome"]][["n2i"]],
      mu_samples = mu_samples,
      tau_within = tau_within_samples,
      prior_pi   = priors[["outcome"]][["pi"]]
    )


  } else if (outcome_type == "pois") {

    # compute marginal log-likelihood via integration over theta and phi
    # the new .outcome_pdf.pois handles integration internally
    log_lik <- .outcome_pdf.pois(
      x1i        = data[["outcome"]][["x1i"]],
      x2i        = data[["outcome"]][["x2i"]],
      t1i        = data[["outcome"]][["t1i"]],
      t2i        = data[["outcome"]][["t2i"]],
      mu_samples = mu_samples,
      tau_within = tau_within_samples,
      prior_phi  = priors[["outcome"]][["phi"]]
    )


  }

  # add column names for observations
  colnames(log_lik) <- paste0("log_lik[", seq_len(K), "]")

  return(log_lik)
}
