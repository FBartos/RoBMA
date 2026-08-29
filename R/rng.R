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
# @param mu_samples       S x K matrix of location samples (with cluster effects if multilevel)
# @param tau_within       S x K matrix of within-cluster heterogeneity samples
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
  total_sd <- .root_sum_squares(tau_within, sei_mat)

  # sample from N(mu, total_sd) for each cell
  # matrix(rnorm(S*K), S, K) * total_sd + mu is vectorized sampling
  # equivalent: for each s,k: y[s,k] ~ N(mu[s,k], total_sd[s,k])
  response_samples <- mu_samples + matrix(stats::rnorm(S * K), nrow = S, ncol = K) * total_sd

  return(response_samples)
}


.outcome_rng.norm_known_v <- function(mu_samples, tau_within, known_V) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  sampling_noise      <- .known_v_sampling_noise(known_V, S = S, K = K)
  heterogeneity_noise <- matrix(stats::rnorm(S * K), nrow = S, ncol = K) *
    tau_within
  response_samples    <- mu_samples + sampling_noise + heterogeneity_noise

  return(response_samples)
}


# Draw known-V sampling noise without materializing a global covariance factor.
.known_v_sampling_noise <- function(known_V, S, K) {

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V covariance dimensions do not match prediction rows.",
         call. = FALSE)
  }

  sampling_noise <- matrix(0, nrow = S, ncol = K)
  independent    <- .known_v_independent_indices(known_V)
  if (length(independent) > 0L) {
    independent_noise <- matrix(
      stats::rnorm(S * length(independent)),
      nrow = S,
      ncol = length(independent)
    )
    sampling_noise[, independent] <- sweep(
      independent_noise,
      MARGIN = 2L,
      STATS  = sqrt(.known_v_diagonal(known_V)[independent]),
      FUN    = "*"
    )
  }

  for (block in .known_v_correlated_blocks(known_V)) {
    index  <- block[["index"]]
    factor <- .known_v_sampling_factor(block[["covariance"]])
    sampling_noise[, index] <- matrix(
      stats::rnorm(S * length(index)),
      nrow = S,
      ncol = length(index)
    ) %*% factor
  }

  return(sampling_noise)
}


.outcome_rng.norm_known_v_covariance <- function(mu_samples,
                                                 covariance_samples) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)

  if (length(dim(covariance_samples)) != 3L ||
      dim(covariance_samples)[1L] != S ||
      dim(covariance_samples)[2L] != K ||
      dim(covariance_samples)[3L] != K) {
    stop("Known-V response covariance samples have inconsistent dimensions.",
         call. = FALSE)
  }

  response_samples <- mu_samples
  for (s in seq_len(S)) {
    covariance <- matrix(
      covariance_samples[s, , ],
      nrow = K,
      ncol = K
    )
    factorization <- .covariance_factorization(covariance)
    factor        <- .covariance_sampling_factor(factorization)
    if (is.null(factor)) {
      stop("Known-V response covariance is not positive semidefinite.",
           call. = FALSE)
    }

    response_samples[s, ] <- mu_samples[s, ] +
      as.vector(stats::rnorm(K) %*% factor)
  }

  return(response_samples)
}


# ---------------------------------------------------------------------------- #
# .outcome_rng.selnorm
# ---------------------------------------------------------------------------- #
#
# Sample observed effect sizes from the selected-normal kernel.
#
# ---------------------------------------------------------------------------- #
.outcome_rng.selnorm <- function(mu_samples, tau_within, sei,
                                 selection_context) {

  S          <- nrow(mu_samples)
  K          <- ncol(mu_samples)
  sei_mat    <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd   <- .root_sum_squares(tau_within, sei_mat)
  selection_context <- BayesTools::selection_context_validate(
    context   = selection_context,
    n_samples = S,
    required  = "use_normal"
  )
  use_normal <- selection_context[["use_normal"]]

  if (all(use_normal)) {
    return(.outcome_rng.norm(
      mu_samples = mu_samples,
      tau_within = tau_within,
      sei        = sei
    ))
  }

  if (any(use_normal)) {
    out <- matrix(NA_real_, nrow = S, ncol = K)

    normal_rows <- which(use_normal)
    step_rows   <- which(!use_normal)

    out[normal_rows, ] <- .outcome_rng.norm(
      mu_samples = mu_samples[normal_rows, , drop = FALSE],
      tau_within = tau_within[normal_rows, , drop = FALSE],
      sei        = sei
    )

    out[step_rows, ] <- .selection_step_rng_matrix(
      mean              = mu_samples[step_rows, , drop = FALSE],
      sd                = total_sd[step_rows, , drop = FALSE],
      sei               = sei,
      selection_context = BayesTools::selection_context_subset_rows(
        context = selection_context,
        rows    = step_rows
      )
    )

    return(out)
  }

  return(.selection_step_rng_matrix(
    mean              = mu_samples,
    sd                = total_sd,
    sei               = sei,
    selection_context = selection_context
  ))
}


# Draw one finite-vector product-selected Gaussian response per posterior row.
# Rejection is performed independently by structural covariance block. The
# native sampler normalizes each row's nonnegative relative weights by their
# maximum, so accepting with their product gives the exact selected density.
.outcome_rng.selnorm_mvn <- function(
    mu_samples, covariance_samples, sei, selection_context,
    dependency_blocks, max_attempts = 100000L) {

  S <- nrow(mu_samples)
  K <- ncol(mu_samples)
  if (!identical(dim(covariance_samples), c(S, K, K)) ||
      any(!is.finite(covariance_samples))) {
    stop("Selected response covariance samples have invalid dimensions.",
         call. = FALSE)
  }
  if (!is.numeric(sei) || length(sei) != K || anyNA(sei) ||
      any(!is.finite(sei)) || any(sei <= 0)) {
    stop("Selected response standard errors must be finite and positive.",
         call. = FALSE)
  }
  BayesTools::check_int(
    max_attempts,
    "max_attempts",
    lower = 1L,
    check_length = 1L,
    allow_NA = FALSE
  )
  selection_context <- BayesTools::selection_context_validate(
    context   = selection_context,
    n_samples = S,
    required  = c("omega", "kernel_mode", "use_normal")
  )
  if (any(!selection_context[["kernel_mode"]] %in%
          c(SELKERNEL_NORMAL, SELKERNEL_STEP))) {
    stop("Exact selected response simulation requires a step selection kernel.",
         call. = FALSE)
  }
  omega <- selection_context[["omega"]]
  if (any(!is.finite(omega)) || any(omega < 0)) {
    stop("Exact selected response simulation requires nonnegative weights.",
         call. = FALSE)
  }
  .known_v_validate_dependency_blocks(dependency_blocks, K)

  p_cuts  <- .selection_assert_p_cuts(selection_context[["p_cuts"]])
  z_lower <- stats::qnorm(p_cuts[-1L], lower.tail = FALSE)
  z_upper <- c(Inf, z_lower[-length(z_lower)])
  result <- .Call(
    "RoBMA_selnorm_mnorm_step_rng_batch",
    .native_numeric_matrix(mu_samples),
    covariance_samples,
    .native_numeric_vector(sei),
    .native_numeric_matrix(omega),
    .native_numeric_vector(z_lower),
    .native_numeric_vector(z_upper),
    .native_integer_vector(selection_context[["sign"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    dependency_blocks,
    .native_integer_vector(max_attempts),
    PACKAGE = "RoBMA"
  )
  if (!is.list(result) ||
      !identical(names(result), c(
        "draws", "failure_code", "failure_size"
      ))) {
    stop("The exact selected response native kernel returned invalid output.",
         call. = FALSE)
  }
  if (result[["failure_code"]] == 1L) {
    stop("Selected response covariance is not positive semidefinite.",
         call. = FALSE)
  }
  if (result[["failure_code"]] == 2L) {
    stop(
      "Exact selected response RNG was rejected by diagnostics: no ",
      "proposal was accepted in ", max_attempts, " attempts for a ",
      "dependency block of size ", result[["failure_size"]], ". Use ",
      "'bias_adjusted = TRUE' to draw responses before selection.",
      call. = FALSE
    )
  }
  if (result[["failure_code"]] == 3L) {
    stop("Exact selected response simulation requires a positive weight.",
         call. = FALSE)
  }
  if (result[["failure_code"]] == 4L) {
    stop("Selected response covariance must be symmetric.", call. = FALSE)
  }
  if (result[["failure_code"]] != 0L ||
      !identical(dim(result[["draws"]]), c(S, K)) ||
      any(!is.finite(result[["draws"]]))) {
    stop("The exact selected response native kernel returned invalid output.",
         call. = FALSE)
  }

  return(result[["draws"]])
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
    outcome_samples_ai[, k] <- stats::rbinom(
      n    = S,
      size = n1i[k],
      prob = .inv_logit(logit_p1[, k])
    )
    outcome_samples_ci[, k] <- stats::rbinom(
      n    = S,
      size = n2i[k],
      prob = .inv_logit(logit_p2[, k])
    )
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
# where phi is the midpoint log rate and t1i, t2i are exposure times.
#
# @param mu_samples       S x K matrix of log incidence rate ratio samples
# @param log_phi          S x K matrix of midpoint log-rate samples
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
    outcome_samples_x1i[, k] <- stats::rpois(
      n      = S,
      lambda = exp(log_lambda1[, k])
    )
    outcome_samples_x2i[, k] <- stats::rpois(
      n      = S,
      lambda = exp(log_lambda2[, k])
    )
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
