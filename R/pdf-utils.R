# ============================================================================ #
# Log-Likelihood Utilities
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

  valid_rows <- rowSums(is.na(lx)) == 0
  if (!any(valid_rows)) {
    return(rep(NA_real_, nrow(lx)))
  }

  valid_index <- which(valid_rows)
  lx_valid    <- lx[valid_rows, , drop = FALSE]
  row_max     <- lx_valid[cbind(seq_len(nrow(lx_valid)), max.col(lx_valid, ties.method = "first"))]
  out         <- rep(NA_real_, nrow(lx))

  infinite_rows <- is.infinite(row_max)
  out[valid_index[infinite_rows]] <- row_max[infinite_rows]

  finite_rows <- is.finite(row_max)
  if (any(finite_rows)) {
    out[valid_index[finite_rows]] <- row_max[finite_rows] +
      log(rowSums(exp(lx_valid[finite_rows, , drop = FALSE] - row_max[finite_rows])))
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
# .has_native_glmm
# ---------------------------------------------------------------------------- #
#
# @param outcome_type character; GLMM outcome type, either "bin" or "pois".
#
# @return TRUE when native GLMM likelihood kernels are available for the
# requested outcome type.
#
# ---------------------------------------------------------------------------- #
.has_native_glmm <- function(outcome_type) {

  if (outcome_type == "bin") {
    return(is.loaded("RoBMA_glmm_binom_marginal_loglik", PACKAGE = "RoBMA"))
  }
  if (outcome_type == "pois") {
    return(is.loaded("RoBMA_glmm_pois_marginal_loglik", PACKAGE = "RoBMA"))
  }

  stop("Unknown GLMM outcome type.", call. = FALSE)
}


.has_native_glmm_cluster <- function() {

  return(
    is.loaded("RoBMA_glmm_binom_cluster_loglik", PACKAGE = "RoBMA") &&
      is.loaded("RoBMA_glmm_pois_cluster_loglik",  PACKAGE = "RoBMA")
  )
}


.has_native_glmm_row_sum <- function(outcome_type, cluster = FALSE,
                                     conditional = FALSE) {

  suffix <- if (conditional) {
    "conditional_loglik_sum"
  } else if (cluster) {
    "cluster_loglik_row_sum"
  } else {
    "marginal_loglik_row_sum"
  }

  if (outcome_type == "bin") {
    return(is.loaded(
      paste("RoBMA_glmm_binom", suffix, sep = "_"),
      PACKAGE = "RoBMA"
    ))
  }
  if (outcome_type == "pois") {
    return(is.loaded(
      paste("RoBMA_glmm_pois", suffix, sep = "_"),
      PACKAGE = "RoBMA"
    ))
  }

  stop("Unknown GLMM outcome type.", call. = FALSE)
}



.has_native_norm_cluster_quadrature <- function(selection = FALSE) {

  if (isTRUE(selection)) {
    return(is.loaded("RoBMA_selnorm_cluster_loglik", PACKAGE = "RoBMA"))
  }

  return(is.loaded("RoBMA_norm_cluster_loglik", PACKAGE = "RoBMA"))
}


.has_native_norm_loglik_row_sum <- function(selection = FALSE,
                                            cluster = FALSE) {

  if (isTRUE(cluster)) {
    if (isTRUE(selection)) {
      return(is.loaded("RoBMA_selnorm_cluster_loglik_row_sum", PACKAGE = "RoBMA"))
    }
    return(is.loaded("RoBMA_norm_cluster_loglik_row_sum", PACKAGE = "RoBMA"))
  }

  if (isTRUE(selection)) {
    return(is.loaded("RoBMA_selnorm_kernel_loglik_row_sum", PACKAGE = "RoBMA"))
  }

  return(is.loaded("RoBMA_norm_loglik_row_sum", PACKAGE = "RoBMA"))
}


.has_native_selnorm_cluster_location_grid <- function() {

  return(is.loaded("RoBMA_selnorm_cluster_location_grid", PACKAGE = "RoBMA"))
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
.glmm_prior_quantile_grid_cache <- new.env(parent = emptyenv())

.glmm_binom_logit_pi_grid_cache <- new.env(parent = emptyenv())

.glmm_pois_log_phi_grid_cache   <- new.env(parent = emptyenv())

.glmm_prior_grid_hash_cache     <- new.env(parent = emptyenv())


.glmm_prior_grid_hash_compute <- function(prior) {

  bytes <- as.integer(serialize(prior, NULL, version = 3))
  hash1 <- 5381
  hash2 <- 0

  for (byte in bytes) {
    hash1 <- (hash1 * 33 + byte) %% 2147483647
    hash2 <- (hash2 * 65599 + byte) %% 2147483629
  }

  return(c(
    sprintf("%08x", as.integer(hash1)),
    sprintf("%08x", as.integer(hash2))
  ))
}


.glmm_prior_grid_hash <- function(prior) {

  cached_prior <- get0(
    "prior",
    envir      = .glmm_prior_grid_hash_cache,
    inherits   = FALSE,
    ifnotfound = NULL
  )
  if (!is.null(cached_prior) && identical(prior, cached_prior)) {
    return(get("hash", envir = .glmm_prior_grid_hash_cache, inherits = FALSE))
  }

  hash <- .glmm_prior_grid_hash_compute(prior)
  assign("prior", prior, envir = .glmm_prior_grid_hash_cache)
  assign("hash", hash, envir = .glmm_prior_grid_hash_cache)

  return(hash)
}


.glmm_prior_grid_key <- function(prior, n, suffix = NULL) {

  hash <- .glmm_prior_grid_hash(prior)

  return(paste(
    n,
    hash[[1L]],
    hash[[2L]],
    suffix,
    sep = "|"
  ))
}


.glmm_prior_quantile_grid <- function(prior, n) {

  key <- .glmm_prior_grid_key(prior, n)
  if (exists(key, envir = .glmm_prior_quantile_grid_cache, inherits = FALSE)) {
    return(get(key, envir = .glmm_prior_quantile_grid_cache, inherits = FALSE))
  }

  gl      <- .gauss_legendre_nodes(n)
  q_grid  <- as.numeric(BayesTools::quant(prior, gl[["nodes"]]))

  if (anyNA(q_grid) || any(!is.finite(q_grid))) {
    stop("The GLMM nuisance prior produced non-finite quadrature nodes.",
         call. = FALSE)
  }

  out <- list(
    grid        = q_grid,
    log_weights = gl[["log_weights"]]
  )

  assign(key, out, envir = .glmm_prior_quantile_grid_cache)
  return(out)
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
  key     <- .glmm_prior_grid_key(prior_pi, n_pi, paste("bin", K, sep = "|"))
  if (exists(key, envir = .glmm_binom_logit_pi_grid_cache, inherits = FALSE)) {
    return(get(key, envir = .glmm_binom_logit_pi_grid_cache, inherits = FALSE))
  }

  pi_grid <- .glmm_prior_quantile_grid(prior = prior_pi, n = n_pi)

  pi_nodes <- pi_grid[["grid"]]
  if (any(pi_nodes <= 0 | pi_nodes >= 1)) {
    stop(
      "The binomial GLMM baserate prior must produce quadrature nodes ",
      "strictly inside (0, 1).",
      call. = FALSE
    )
  }

  logit_pi_grid  <- matrix(qlogis(pi_nodes), nrow = n_pi, ncol = K)
  log_pi_weights <- matrix(pi_grid[["log_weights"]], nrow = n_pi, ncol = K)

  out <- list(
    grid        = logit_pi_grid,
    log_weights = log_pi_weights
  )

  assign(key, out, envir = .glmm_binom_logit_pi_grid_cache)
  return(out)
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
  key      <- .glmm_prior_grid_key(prior_phi, n_phi, paste("pois", K, sep = "|"))
  if (exists(key, envir = .glmm_pois_log_phi_grid_cache, inherits = FALSE)) {
    return(get(key, envir = .glmm_pois_log_phi_grid_cache, inherits = FALSE))
  }

  phi_grid <- .glmm_prior_quantile_grid(prior = prior_phi, n = n_phi)

  log_phi_grid    <- matrix(phi_grid[["grid"]], nrow = n_phi, ncol = K)
  log_phi_weights <- matrix(phi_grid[["log_weights"]], nrow = n_phi, ncol = K)

  out <- list(
    grid        = log_phi_grid,
    log_weights = log_phi_weights
  )

  assign(key, out, envir = .glmm_pois_log_phi_grid_cache)
  return(out)
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
