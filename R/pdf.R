# ============================================================================ #
# Outcome PDF Functions (Log-Likelihood)
# ============================================================================ #
#
# These functions compute pointwise log-likelihoods for each observation and
# posterior sample. They are used for LOO-PSIS diagnostics and model comparison.
#
# Parallels the structure of rng.R but returns density values instead of
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
  if (!is.matrix(lx)) {
    lx <- as.matrix(lx)
  }

  if (!anyNA(lx)) {
    row_max <- lx[cbind(seq_len(nrow(lx)), max.col(lx, ties.method = "first"))]
    out     <- row_max

    finite_rows <- is.finite(row_max)
    if (any(finite_rows)) {
      out[finite_rows] <- row_max[finite_rows] +
        log(rowSums(exp(lx[finite_rows, , drop = FALSE] - row_max[finite_rows])))
    }

    return(out)
  }

  has_value <- rowSums(!is.na(lx)) > 0
  lx[is.na(lx)] <- -Inf

  row_max <- lx[cbind(seq_len(nrow(lx)), max.col(lx, ties.method = "first"))]
  out     <- rep(NA_real_, nrow(lx))

  infinite_rows <- has_value & is.infinite(row_max)
  out[infinite_rows] <- row_max[infinite_rows]

  finite_rows <- has_value & is.finite(row_max)
  if (any(finite_rows)) {
    out[finite_rows] <- row_max[finite_rows] +
      log(rowSums(exp(lx[finite_rows, , drop = FALSE] - row_max[finite_rows])))
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .apply_log_lik_weights
# ---------------------------------------------------------------------------- #
#
# Match JAGS weighted-density semantics by multiplying log-likelihood
# contributions by the user-supplied data weights.
#
# @param log_lik S x K matrix of log-likelihood values.
# @param weights numeric vector of length K, or NULL.
#
# @return weighted log-likelihood matrix.
#
# ---------------------------------------------------------------------------- #
.apply_log_lik_weights <- function(log_lik, weights) {

  if (is.null(weights)) {
    return(log_lik)
  }

  if (length(weights) != ncol(log_lik)) {
    stop("Data weights length does not match the log-likelihood matrix.",
         call. = FALSE)
  }

  return(sweep(log_lik, 2, weights, "*"))
}


# ---------------------------------------------------------------------------- #
# .get_log_lik_data_weights
# ---------------------------------------------------------------------------- #
#
# @param object brma object.
#
# @return numeric vector of data weights, or NULL.
#
# ---------------------------------------------------------------------------- #
.get_log_lik_data_weights <- function(object) {

  if (!.is_weights(object)) {
    return(NULL)
  }

  return(object[["data"]][["outcome"]][["weights"]])
}


# ---------------------------------------------------------------------------- #
# .glmm_likelihood_weights
# ---------------------------------------------------------------------------- #
#
# Validate and default GLMM pair likelihood weights.
#
# ---------------------------------------------------------------------------- #
.glmm_likelihood_weights <- function(weights, K) {

  if (is.null(weights)) {
    return(rep(1, K))
  }

  BayesTools::check_real(
    weights, "weights",
    check_length = K, allow_NULL = FALSE, allow_NA = FALSE,
    lower = 0, allow_bound = FALSE
  )

  return(as.numeric(weights))
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
                                 selection_context, weights = NULL) {

  S       <- nrow(mu_samples)
  K       <- ncol(mu_samples)
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd <- sqrt(tau_within^2 + sei_mat^2)

  return(.selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = mu_samples,
    sigma_num      = total_sd,
    mu_norm        = mu_samples,
    sigma_norm     = total_sd,
    sei            = sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]],
    weights        = weights
  ))
}


# ---------------------------------------------------------------------------- #
# .has_native_glmm
# ---------------------------------------------------------------------------- #
#
# @return TRUE when native GLMM likelihood kernels are available.
#
# ---------------------------------------------------------------------------- #
.has_native_glmm <- function() {

  return(is.loaded("RoBMA_glmm_binom_marginal_loglik", PACKAGE = "RoBMA"))
}

.has_native_glmm_cluster <- function() {

  return(
    is.loaded("RoBMA_glmm_binom_cluster_loglik", PACKAGE = "RoBMA") &&
      is.loaded("RoBMA_glmm_pois_cluster_loglik",  PACKAGE = "RoBMA")
  )
}


# ---------------------------------------------------------------------------- #
# .cluster_indices_flatten
# ---------------------------------------------------------------------------- #
#
# Flatten a cluster index list while preserving list order.
#
# ---------------------------------------------------------------------------- #
.cluster_indices_flatten <- function(cluster_indices) {

  return(list(
    index = as.integer(unlist(cluster_indices, use.names = FALSE)),
    size  = as.integer(lengths(cluster_indices))
  ))
}


# ---------------------------------------------------------------------------- #
# .gauss_legendre_nodes
# ---------------------------------------------------------------------------- #
#
# Compute Gauss-Legendre quadrature nodes and weights on (0, 1).
#
# ---------------------------------------------------------------------------- #
.gauss_legendre_nodes_cache <- new.env(parent = emptyenv())

.gauss_legendre_nodes <- function(n) {

  key <- as.character(n)
  if (exists(key, envir = .gauss_legendre_nodes_cache, inherits = FALSE)) {
    return(get(key, envir = .gauss_legendre_nodes_cache, inherits = FALSE))
  }

  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)

  J <- diag(0, n)
  J[cbind(i, i + 1)] <- b
  J[cbind(i + 1, i)] <- b

  eig <- eigen(J, symmetric = TRUE)

  nodes_raw   <- eig$values
  weights_raw <- 2 * eig$vectors[1, ]^2

  nodes   <- (nodes_raw + 1) / 2
  weights <- weights_raw / 2

  ord <- order(nodes)

  out <- list(
    nodes       = nodes[ord],
    weights     = weights[ord],
    log_weights = log(weights[ord])
  )
  assign(key, out, envir = .gauss_legendre_nodes_cache)

  return(out)
}


# ---------------------------------------------------------------------------- #
# .glmm_prior_quantile_grid
# ---------------------------------------------------------------------------- #
#
# Construct quadrature nodes under a prior using the probability integral
# transform. The returned weights integrate with respect to the prior measure.
#
# ---------------------------------------------------------------------------- #
.glmm_prior_quantile_grid <- function(prior, n) {

  gl      <- .gauss_legendre_nodes(n)
  q_grid  <- as.numeric(BayesTools::quant(prior, gl[["nodes"]]))

  if (anyNA(q_grid) || any(!is.finite(q_grid))) {
    stop("The GLMM nuisance prior produced non-finite quadrature nodes.",
         call. = FALSE)
  }

  return(list(
    grid        = q_grid,
    log_weights = gl[["log_weights"]]
  ))
}


# ---------------------------------------------------------------------------- #
# .glmm_binom_logit_pi_grid
# ---------------------------------------------------------------------------- #
#
# Construct prior-scale nuisance grids for binomial GLMM likelihoods.
#
# ---------------------------------------------------------------------------- #
.glmm_binom_logit_pi_grid <- function(ai, ci, n1i, n2i, prior_pi, n_pi) {

  K       <- length(ai)
  pi_grid <- .glmm_prior_quantile_grid(prior = prior_pi, n = n_pi)

  if (any(pi_grid[["grid"]] < 0 | pi_grid[["grid"]] > 1)) {
    stop("The binomial GLMM baserate prior produced nodes outside [0, 1].",
         call. = FALSE)
  }

  pi_nodes <- pmin(pmax(pi_grid[["grid"]], .Machine$double.eps),
                   1 - .Machine$double.eps)

  logit_pi_grid  <- matrix(qlogis(pi_nodes), nrow = n_pi, ncol = K)
  log_pi_weights <- matrix(pi_grid[["log_weights"]], nrow = n_pi, ncol = K)

  return(list(
    grid        = logit_pi_grid,
    log_weights = log_pi_weights
  ))
}


# ---------------------------------------------------------------------------- #
# .glmm_pois_log_phi_grid
# ---------------------------------------------------------------------------- #
#
# Construct prior-scale nuisance grids for Poisson GLMM likelihoods.
#
# ---------------------------------------------------------------------------- #
.glmm_pois_log_phi_grid <- function(x1i, x2i, t1i, t2i, prior_phi, n_phi) {

  K        <- length(x1i)
  phi_grid <- .glmm_prior_quantile_grid(prior = prior_phi, n = n_phi)

  log_phi_grid    <- matrix(phi_grid[["grid"]], nrow = n_phi, ncol = K)
  log_phi_weights <- matrix(phi_grid[["log_weights"]], nrow = n_phi, ncol = K)

  return(list(
    grid        = log_phi_grid,
    log_weights = log_phi_weights
  ))
}


# ---------------------------------------------------------------------------- #
# .outcome_pdf.binom
# ---------------------------------------------------------------------------- #
#
# Native wrapper for binomial GLMM marginal log-likelihoods.
#
# ---------------------------------------------------------------------------- #
.outcome_pdf.binom <- function(ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi,
                               weights = NULL, n_theta = 15, n_pi = 30) {

  if (!.has_native_glmm()) {
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
.gauss_hermite_nodes_cache <- new.env(parent = emptyenv())

.gauss_hermite_nodes <- function(n) {

  key <- as.character(n)
  if (exists(key, envir = .gauss_hermite_nodes_cache, inherits = FALSE)) {
    return(get(key, envir = .gauss_hermite_nodes_cache, inherits = FALSE))
  }

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

  out <- list(
    nodes       = nodes[ord],
    weights     = weights[ord],
    log_weights = log(weights[ord])
  )
  assign(key, out, envir = .gauss_hermite_nodes_cache)

  return(out)
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

  if (!.has_native_glmm()) {
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
# @param log_phi    S x K matrix of log baseline-rate samples
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


# ---------------------------------------------------------------------------- #
# .get_cluster_likelihood_n_gamma
# ---------------------------------------------------------------------------- #
#
# @return integer; number of Gauss-Hermite nodes for cluster likelihoods.
#
# ---------------------------------------------------------------------------- #
.get_cluster_likelihood_n_gamma <- function() {

  n_gamma <- RoBMA.get_option("cluster_likelihood.n_gamma")
  BayesTools::check_int(n_gamma, "cluster_likelihood.n_gamma", lower = 3)

  return(as.integer(n_gamma))
}


# ---------------------------------------------------------------------------- #
# .log_lik_cluster_gamma_quadrature
# ---------------------------------------------------------------------------- #
#
# Integrate the held-out cluster effect gamma_g with Gauss-Hermite quadrature.
#
# @param cluster_indices     list of cluster index vectors.
# @param mu_samples          S x K matrix of fixed-effect means.
# @param tau_between_samples S x K matrix of cluster-level SDs.
# @param log_lik_fun         function(idx, mu_node) returning S x length(idx)
#                            conditional log-likelihood matrix.
# @param weights             optional data weights.
# @param n_gamma             number of Gauss-Hermite nodes.
#
# @return S x G cluster-unit log-likelihood matrix.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_gamma_quadrature <- function(cluster_indices, mu_samples,
                                              tau_between_samples, log_lik_fun,
                                              weights = NULL,
                                              n_gamma = .get_cluster_likelihood_n_gamma()) {

  gh      <- .gauss_hermite_nodes(n_gamma)
  S       <- nrow(mu_samples)
  G       <- length(cluster_indices)
  log_lik <- matrix(NA_real_, nrow = S, ncol = G)

  for (g in seq_along(cluster_indices)) {
    idx       <- cluster_indices[[g]]
    log_terms <- matrix(NA_real_, nrow = S, ncol = n_gamma)

    for (j in seq_len(n_gamma)) {
      mu_node <- mu_samples[, idx, drop = FALSE] +
        gh$nodes[j] * tau_between_samples[, idx, drop = FALSE]

      point_log_lik <- log_lik_fun(idx, mu_node)
      point_log_lik <- .apply_log_lik_weights(point_log_lik, weights[idx])

      log_terms[, j] <- rowSums(point_log_lik) + gh$log_weights[j]
    }

    log_lik[, g] <- .rowLogSumExps(log_terms)
  }

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .log_lik.brma
# ---------------------------------------------------------------------------- #
#
# Compute the log-likelihood matrix for a brma object.
#
# @param object brma object.
# @param unit   character; output/deletion unit.
#
# @return S x K or S x G matrix of log-likelihood values.
#
# ---------------------------------------------------------------------------- #
.log_lik.brma <- function(object, unit = "estimate") {

  unit <- .normalize_unit(unit)

  if (unit == "estimate") {
    return(.log_lik_estimate.brma(object))
  } else {
    return(.log_lik_cluster.brma(object))
  }
}


# ---------------------------------------------------------------------------- #
# .pdf.brma
# ---------------------------------------------------------------------------- #
#
# Back-end wrapper kept for internal callers that use the old helper name.
#
# @param object brma object.
# @param unit   character; output/deletion unit.
#
# @return S x K or S x G matrix of log-likelihood values.
#
# ---------------------------------------------------------------------------- #
.pdf.brma <- function(object, unit = "estimate") {

  return(.log_lik.brma(object = object, unit = unit))
}


# ---------------------------------------------------------------------------- #
# .estimate_likelihood_setup.brma
# ---------------------------------------------------------------------------- #
#
# Extract common posterior quantities for estimate-unit predictive likelihoods.
#
# The estimate-unit target conditions on fitted cluster effects for multilevel
# models and integrates over estimate-level heterogeneity.
#
# @param object brma object.
#
# @return list with model metadata and posterior matrices.
#
# ---------------------------------------------------------------------------- #
.estimate_likelihood_setup.brma <- function(object) {

  priors            <- object[["priors"]]
  data              <- object[["data"]]
  is_multilevel     <- .is_multilevel(object)
  is_mods           <- .is_mods(object)
  is_scale          <- .is_scale(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)
  posterior_samples <- .get_posterior_samples(object[["fit"]])

  yi  <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)
  K   <- length(yi)

  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = data[["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = data, "scale") else NULL,
    scale_priors      = priors[["scale"]],
    is_scale          = is_scale,
    is_multilevel     = is_multilevel,
    K                 = K,
    posterior_samples = posterior_samples
  )

  mu_samples <- .evaluate.brma.mu(
    fit               = object[["fit"]],
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (is_mods) .create_fit_formula_list(data = data, "mods") else NULL,
    mods_priors       = priors[["mods"]],
    is_mods           = is_mods,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction,
    bias_adjusted     = FALSE,
    K                 = K,
    posterior_samples = posterior_samples
  )

  fit_data <- NULL
  if (is_multilevel || is_weightfunction) {
    fit_data <- .create_fit_data(data = data, priors = priors)
  }

  if (is_multilevel) {
    mu_samples <- mu_samples + .evaluate.brma.cluster_effects(
      fit               = object[["fit"]],
      tau_between       = tau_result[["tau_between"]],
      cluster           = fit_data[["cluster"]],
      same_data         = TRUE,
      effect_direction  = effect_direction,
      posterior_samples = posterior_samples
    )
  }

  return(list(
    priors            = priors,
    data              = data,
    fit_data          = fit_data,
    yi                = yi,
    sei               = sei,
    K                 = K,
    S                 = nrow(mu_samples),
    mu                = mu_samples,
    tau_within        = tau_result[["tau_within"]],
    tau_between       = tau_result[["tau_between"]],
    is_weightfunction = is_weightfunction,
    outcome_type      = outcome_type,
    effect_direction  = effect_direction,
    posterior_samples = posterior_samples
  ))
}


# ---------------------------------------------------------------------------- #
# .log_lik_estimate.brma
# ---------------------------------------------------------------------------- #
#
# Estimate-unit likelihood, one contribution per effect-size estimate.
#
# For normal multilevel models this is the model likelihood factor conditional
# on the cluster effect and marginal over estimate-level heterogeneity.
#
# @param object brma object.
#
# @return S x K matrix of log-likelihood values.
#
# ---------------------------------------------------------------------------- #
.log_lik_estimate.brma <- function(object) {

  setup <- .estimate_likelihood_setup.brma(object)

  priors            <- setup[["priors"]]
  data              <- setup[["data"]]
  yi                <- setup[["yi"]]
  sei               <- setup[["sei"]]
  K                 <- setup[["K"]]
  mu_samples         <- setup[["mu"]]
  tau_within_samples <- setup[["tau_within"]]
  is_weightfunction <- setup[["is_weightfunction"]]
  outcome_type      <- setup[["outcome_type"]]
  effect_direction  <- setup[["effect_direction"]]
  data_weights      <- .get_log_lik_data_weights(object)

  ### extract posterior samples once if needed
  if (is_weightfunction) {
    posterior_samples <- setup[["posterior_samples"]]
  }

  ### compute log-likelihood based on outcome type
  if (outcome_type == "norm") {

    # flip for negative effect direction
    if (effect_direction == "negative") {
      mu_samples_pdf <- -mu_samples
      yi_pdf         <- -yi
    } else {
      mu_samples_pdf <- mu_samples
      yi_pdf         <- yi
    }

    # dispatch between weighted and standard normal
    if (is_weightfunction) {
      selection_context <- .selection_context(
        object            = object,
        posterior_samples = posterior_samples
      )

      log_lik <- .outcome_pdf.selnorm(
        yi                = yi,
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        sei               = sei,
        selection_context = selection_context
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
      prior_pi   = priors[["outcome"]][["pi"]],
      weights    = data_weights
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
      prior_phi  = priors[["outcome"]][["phi"]],
      weights    = data_weights
    )

  } else {

    stop("Unsupported outcome type for estimate-unit log-likelihood.",
         call. = FALSE)

  }

  # add column names and target metadata
  if (outcome_type == "norm") {
    log_lik <- .apply_log_lik_weights(log_lik, data_weights)
  }

  colnames(log_lik) <- paste0("log_lik[", seq_len(K), "]")
  attr(log_lik, "RoBMA_target") <- list(
    unit               = "estimate",
    conditioning_depth = "estimate",
    n                  = K,
    targets            = seq_len(K),
    data_hash          = .get_outcome_hash(object)
  )

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .log_lik_cluster.brma
# ---------------------------------------------------------------------------- #
#
# Cluster-unit likelihood for multilevel models.
#
# Each column is the joint held-out-cluster log-likelihood contribution.
#
# @param object brma object.
#
# @return S x G matrix of log-likelihood values.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster.brma <- function(object) {

  if (!.is_multilevel(object)) {
    stop("Cluster-unit log-likelihood is only available for multilevel models.",
         call. = FALSE)
  }

  outcome_type      <- .outcome_type(object)
  is_weightfunction <- .is_weightfunction(object)
  is_weights        <- .is_weights(object)

  if (outcome_type == "norm" && !is_weightfunction && !is_weights) {
    return(.log_lik_cluster_norm.brma(object))
  } else if (outcome_type == "norm") {
    return(.log_lik_cluster_norm_quadrature.brma(object))
  } else if (outcome_type %in% c("bin", "pois")) {
    return(.log_lik_cluster_glmm.brma(object))
  } else {
    stop("Unsupported outcome type for cluster-unit log-likelihood.",
         call. = FALSE)
  }
}


# ---------------------------------------------------------------------------- #
# .log_lik_cluster_setup.brma
# ---------------------------------------------------------------------------- #
#
# Extract common posterior matrices for cluster-unit likelihoods.
#
# @param object brma object.
#
# @return list with model metadata and posterior matrices.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_setup.brma <- function(object) {

  priors           <- object[["priors"]]
  data             <- object[["data"]]
  is_mods          <- .is_mods(object)
  is_scale          <- .is_scale(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  effect_direction  <- .effect_direction(object)
  K                 <- nrow(data[["outcome"]])
  posterior_samples <- .get_posterior_samples(object[["fit"]])

  mu_samples <- .evaluate.brma.mu(
    fit               = object[["fit"]],
    outcome_data      = data[["outcome"]],
    mods_data         = data[["mods"]],
    mods_formula      = if (is_mods) .create_fit_formula_list(data = data, "mods") else NULL,
    mods_priors       = priors[["mods"]],
    is_mods           = is_mods,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction,
    bias_adjusted     = FALSE,
    K                 = K,
    posterior_samples = posterior_samples
  )

  tau_result <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = data[["scale"]],
    scale_formula     = if (is_scale) .create_fit_formula_list(data = data, "scale") else NULL,
    scale_priors      = priors[["scale"]],
    is_scale          = is_scale,
    is_multilevel     = TRUE,
    K                 = K,
    posterior_samples = posterior_samples
  )

  return(list(
    priors            = priors,
    data              = data,
    K                 = K,
    S                 = nrow(mu_samples),
    mu                = mu_samples,
    tau_within        = tau_result[["tau_within"]],
    tau_between       = tau_result[["tau_between"]],
    cluster           = .get_cluster_indices(object),
    weights           = .get_log_lik_data_weights(object),
    data_hash         = .get_outcome_hash(object),
    effect_direction  = effect_direction,
    posterior_samples = posterior_samples
  ))
}


# ---------------------------------------------------------------------------- #
# .add_cluster_log_lik_metadata
# ---------------------------------------------------------------------------- #
#
# @param log_lik         S x G log-likelihood matrix.
# @param cluster_indices named list of cluster index vectors.
# @param data_hash       character; hash of the yi/sei target.
#
# @return log-likelihood matrix with names and metadata.
#
# ---------------------------------------------------------------------------- #
.add_cluster_log_lik_metadata <- function(log_lik, cluster_indices, data_hash) {

  cluster_labels    <- names(cluster_indices)
  colnames(log_lik) <- paste0("log_lik_cluster[", cluster_labels, "]")
  attr(log_lik, "RoBMA_target") <- list(
    unit               = "cluster",
    conditioning_depth = "cluster",
    n                  = length(cluster_indices),
    targets            = cluster_labels,
    data_hash          = data_hash
  )

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm.brma
# ---------------------------------------------------------------------------- #
#
# Analytic cluster-unit likelihood for unweighted normal multilevel models
# without selection. The normal cluster effect is integrated by the block
# covariance.
#
# @param object brma object.
#
# @return S x G cluster-unit log-likelihood matrix.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm.brma <- function(object) {

  setup <- .log_lik_cluster_setup.brma(object)
  yi    <- .outcome_data_yi(object)
  vi    <- .outcome_data_vi(object)

  if (setup[["effect_direction"]] == "negative") {
    mu_samples <- -setup[["mu"]]
    yi         <- -yi
  } else {
    mu_samples <- setup[["mu"]]
  }

  cluster_indices <- setup[["cluster"]]
  log_lik         <- matrix(NA_real_, nrow = setup[["S"]], ncol = length(cluster_indices))

  for (g in seq_along(cluster_indices)) {
    idx <- cluster_indices[[g]]

    for (s in seq_len(setup[["S"]])) {
      log_lik[s, g] <- .log_dmvnorm_diag_rank_one(
        x        = yi[idx],
        mean     = mu_samples[s, idx],
        diagonal = vi[idx] + setup[["tau_within"]][s, idx]^2,
        rank_one = setup[["tau_between"]][s, idx]
      )
    }
  }

  return(.add_cluster_log_lik_metadata(log_lik, cluster_indices, setup[["data_hash"]]))
}


# ---------------------------------------------------------------------------- #
# .log_dmvnorm_diag_rank_one
# ---------------------------------------------------------------------------- #
#
# Log-density for N(mean, diag(diagonal) + rank_one %*% t(rank_one)).
#
# ---------------------------------------------------------------------------- #
.log_dmvnorm_diag_rank_one <- function(x, mean, diagonal, rank_one) {

  residual <- x - mean

  if (!all(is.finite(diagonal)) || any(diagonal <= 0)) {
    covariance <- diag(diagonal, nrow = length(diagonal), ncol = length(diagonal)) +
      tcrossprod(rank_one)
    return(mvtnorm::dmvnorm(
      x     = x,
      mean  = mean,
      sigma = covariance,
      log   = TRUE
    ))
  }

  inv_diag <- 1 / diagonal
  inv_rank <- rank_one * inv_diag
  denom    <- 1 + sum(rank_one * inv_rank)

  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    covariance <- diag(diagonal, nrow = length(diagonal), ncol = length(diagonal)) +
      tcrossprod(rank_one)
    return(mvtnorm::dmvnorm(
      x     = x,
      mean  = mean,
      sigma = covariance,
      log   = TRUE
    ))
  }

  log_det <- sum(log(diagonal)) + log(denom)
  quad    <- sum(residual^2 * inv_diag) -
    sum(rank_one * residual * inv_diag)^2 / denom

  return(-0.5 * (length(x) * log(2 * pi) + log_det + quad))
}


# ---------------------------------------------------------------------------- #
# .log_lik_cluster_norm_quadrature.brma
# ---------------------------------------------------------------------------- #
#
# Gamma-quadrature cluster-unit likelihood for selected-normal models and
# data-weighted normal models. Conditional on gamma, per-estimate likelihood
# contributions factorize.
#
# @param object brma object.
#
# @return S x G cluster-unit log-likelihood matrix.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_norm_quadrature.brma <- function(object) {

  setup             <- .log_lik_cluster_setup.brma(object)
  is_weightfunction <- .is_weightfunction(object)
  yi                <- .outcome_data_yi(object)
  sei               <- .outcome_data_sei(object)

  if (setup[["effect_direction"]] == "negative" && !is_weightfunction) {
    mu_samples <- -setup[["mu"]]
    yi         <- -yi
  } else {
    mu_samples <- setup[["mu"]]
  }

  if (is_weightfunction) {
    posterior_samples <- setup[["posterior_samples"]]
    selection_context <- .selection_context(
      object            = object,
      posterior_samples = posterior_samples
    )
  }

  log_lik_fun <- function(idx, mu_node) {
    if (is_weightfunction) {
      local_context <- selection_context
      local_context[["obs_bin"]] <- selection_context[["obs_bin"]][idx]

      return(.outcome_pdf.selnorm(
        yi                = yi[idx],
        mu_samples        = mu_node,
        tau_within        = setup[["tau_within"]][, idx, drop = FALSE],
        sei               = sei[idx],
        selection_context = local_context
      ))
    } else {
      return(.outcome_pdf.norm(
        yi         = yi[idx],
        mu_samples = mu_node,
        tau_within = setup[["tau_within"]][, idx, drop = FALSE],
        sei        = sei[idx]
      ))
    }
  }

  log_lik <- .log_lik_cluster_gamma_quadrature(
    cluster_indices     = setup[["cluster"]],
    mu_samples          = mu_samples,
    tau_between_samples = setup[["tau_between"]],
    log_lik_fun         = log_lik_fun,
    weights             = setup[["weights"]]
  )

  return(.add_cluster_log_lik_metadata(log_lik, setup[["cluster"]], setup[["data_hash"]]))
}


# ---------------------------------------------------------------------------- #
# .log_lik_cluster_glmm.brma
# ---------------------------------------------------------------------------- #
#
# Gamma-quadrature cluster-unit likelihood for binomial and Poisson multilevel
# models. Conditional on gamma, existing estimate-level GLMM quadrature helpers
# are reused.
#
# @param object brma object.
#
# @return S x G cluster-unit log-likelihood matrix.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_glmm.brma <- function(object) {

  setup       <- .log_lik_cluster_setup.brma(object)
  data        <- setup[["data"]]
  priors      <- setup[["priors"]]
  outcome_type <- .outcome_type(object)

  if (.has_native_glmm_cluster()) {
    log_lik <- .log_lik_cluster_glmm_native(
      setup        = setup,
      data         = data,
      priors       = priors,
      outcome_type = outcome_type
    )
  } else {
    log_lik <- .log_lik_cluster_glmm_r(
      setup        = setup,
      data         = data,
      priors       = priors,
      outcome_type = outcome_type
    )
  }

  return(.add_cluster_log_lik_metadata(log_lik, setup[["cluster"]], setup[["data_hash"]]))
}


# ---------------------------------------------------------------------------- #
# .log_lik_cluster_glmm_r
# ---------------------------------------------------------------------------- #
#
# R-composed reference implementation for GLMM cluster likelihoods.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_glmm_r <- function(setup, data, priors, outcome_type,
                                    n_theta = 15,
                                    n_gamma = .get_cluster_likelihood_n_gamma(),
                                    n_pi = 30, n_phi = 30) {

  log_lik_fun <- switch(
    outcome_type,
    "bin" = function(idx, mu_node) {
      .outcome_pdf.binom(
        ai         = data[["outcome"]][["ai"]][idx],
        ci         = data[["outcome"]][["ci"]][idx],
        n1i        = data[["outcome"]][["n1i"]][idx],
        n2i        = data[["outcome"]][["n2i"]][idx],
        mu_samples = mu_node,
        tau_within = setup[["tau_within"]][, idx, drop = FALSE],
        prior_pi   = priors[["outcome"]][["pi"]],
        weights    = setup[["weights"]][idx],
        n_theta    = n_theta,
        n_pi       = n_pi
      )
    },
    "pois" = function(idx, mu_node) {
      .outcome_pdf.pois(
        x1i        = data[["outcome"]][["x1i"]][idx],
        x2i        = data[["outcome"]][["x2i"]][idx],
        t1i        = data[["outcome"]][["t1i"]][idx],
        t2i        = data[["outcome"]][["t2i"]][idx],
        mu_samples = mu_node,
        tau_within = setup[["tau_within"]][, idx, drop = FALSE],
        prior_phi  = priors[["outcome"]][["phi"]],
        weights    = setup[["weights"]][idx],
        n_theta    = n_theta,
        n_phi      = n_phi
      )
    }
  )

  log_lik <- .log_lik_cluster_gamma_quadrature(
    cluster_indices     = setup[["cluster"]],
    mu_samples          = setup[["mu"]],
    tau_between_samples = setup[["tau_between"]],
    log_lik_fun         = log_lik_fun,
    n_gamma             = n_gamma
  )

  return(log_lik)
}


# ---------------------------------------------------------------------------- #
# .log_lik_cluster_glmm_native
# ---------------------------------------------------------------------------- #
#
# Native batched implementation for GLMM cluster likelihoods.
#
# ---------------------------------------------------------------------------- #
.log_lik_cluster_glmm_native <- function(setup, data, priors, outcome_type,
                                         n_theta = 15,
                                         n_gamma = .get_cluster_likelihood_n_gamma(),
                                         n_pi = 30, n_phi = 30) {

  gh_theta <- .gauss_hermite_nodes(n_theta)
  gh_gamma <- .gauss_hermite_nodes(n_gamma)
  cluster  <- .cluster_indices_flatten(setup[["cluster"]])
  weights  <- if (is.null(setup[["weights"]])) NULL else .native_numeric_vector(setup[["weights"]])

  if (outcome_type == "bin") {
    pi_grid <- .glmm_binom_logit_pi_grid(
      ai       = data[["outcome"]][["ai"]],
      ci       = data[["outcome"]][["ci"]],
      n1i      = data[["outcome"]][["n1i"]],
      n2i      = data[["outcome"]][["n2i"]],
      prior_pi = priors[["outcome"]][["pi"]],
      n_pi     = n_pi
    )

    return(.Call(
      "RoBMA_glmm_binom_cluster_loglik",
      .native_integer_vector(data[["outcome"]][["ai"]]),
      .native_integer_vector(data[["outcome"]][["ci"]]),
      .native_integer_vector(data[["outcome"]][["n1i"]]),
      .native_integer_vector(data[["outcome"]][["n2i"]]),
      .native_numeric_matrix(setup[["mu"]]),
      .native_numeric_matrix(setup[["tau_within"]]),
      .native_numeric_matrix(setup[["tau_between"]]),
      cluster[["index"]],
      cluster[["size"]],
      weights,
      .native_numeric_vector(gh_theta[["nodes"]]),
      .native_numeric_vector(gh_theta[["log_weights"]]),
      .native_numeric_matrix(pi_grid[["grid"]]),
      .native_numeric_matrix(pi_grid[["log_weights"]]),
      .native_numeric_vector(gh_gamma[["nodes"]]),
      .native_numeric_vector(gh_gamma[["log_weights"]]),
      PACKAGE = "RoBMA"
    ))
  }

  if (outcome_type == "pois") {
    phi_grid <- .glmm_pois_log_phi_grid(
      x1i       = data[["outcome"]][["x1i"]],
      x2i       = data[["outcome"]][["x2i"]],
      t1i       = data[["outcome"]][["t1i"]],
      t2i       = data[["outcome"]][["t2i"]],
      prior_phi = priors[["outcome"]][["phi"]],
      n_phi     = n_phi
    )

    return(.Call(
      "RoBMA_glmm_pois_cluster_loglik",
      .native_integer_vector(data[["outcome"]][["x1i"]]),
      .native_integer_vector(data[["outcome"]][["x2i"]]),
      .native_numeric_vector(data[["outcome"]][["t1i"]]),
      .native_numeric_vector(data[["outcome"]][["t2i"]]),
      .native_numeric_matrix(setup[["mu"]]),
      .native_numeric_matrix(setup[["tau_within"]]),
      .native_numeric_matrix(setup[["tau_between"]]),
      cluster[["index"]],
      cluster[["size"]],
      weights,
      .native_numeric_vector(gh_theta[["nodes"]]),
      .native_numeric_vector(gh_theta[["log_weights"]]),
      .native_numeric_matrix(phi_grid[["grid"]]),
      .native_numeric_matrix(phi_grid[["log_weights"]]),
      .native_numeric_vector(gh_gamma[["nodes"]]),
      .native_numeric_vector(gh_gamma[["log_weights"]]),
      PACKAGE = "RoBMA"
    ))
  }

  stop("Unsupported outcome type for GLMM cluster-unit log-likelihood.",
       call. = FALSE)
}
