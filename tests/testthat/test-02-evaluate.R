context("Evaluate posterior samples")
# ============================================================================ #
# test-02-evaluate.R
# ============================================================================ #
#
# Unit and integration tests for evaluate.R helper functions
# defined in R/evaluate.R
#
# Structure:
# 1. Unit tests: Mock posterior matrices (no JAGS fit needed)
# 2. Integration tests: Pre-fitted models from list_fits()
#
# ============================================================================ #

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "evaluate")

# ============================================================================ #
# SECTION 1: Unit Tests with Mock Data
# ============================================================================ #
# These tests verify the computational logic of helper functions
# using mock posterior sample matrices (no JAGS fit required)
# ============================================================================ #

test_that(".evaluate.brma.true_effects.norm computes correct BLUPs", {

  # create mock data
  S <- 100  # samples
  K <- 4    # observations

  set.seed(1234)
  mu_samples <- matrix(rnorm(S * K, mean = 0.5, sd = 0.1), nrow = S, ncol = K)
  tau_within <- matrix(abs(rnorm(S * K, mean = 0.2, sd = 0.01)), nrow = S, ncol = K)
  yi  <- c(0.3, 0.5, 0.4, 0.6)
  sei <- c(0.1, 0.15, 0.25, 0.30)

  # compute true effects with same_data = TRUE (BLUPs)
  theta <- .evaluate.brma.true_effects.norm(
    mu_samples = mu_samples,
    tau_within = tau_within,
    yi         = yi,
    sei        = sei,
    same_data  = TRUE
  )

  # verify dimensions
  expect_equal(dim(theta), c(S, K))

  # verify shrinkage behavior: theta should be between yi and mu
  # (on average across samples)
  mean_theta <- colMeans(theta)
  mean_mu    <- colMeans(mu_samples)

  for (k in seq_len(K)) {
    # theta should be closer to data when tau is large relative to se
    # and closer to model when tau is small relative to se
    if (mean(tau_within[, k]) > sei[k]) {
      # high tau: theta closer to yi
      expect_true(abs(mean_theta[k] - yi[k]) < abs(mean_theta[k] - mean_mu[k]),
                  info = paste("Observation", k, "shrinks toward data with high tau"))
    } else {
      expect_true(abs(mean_theta[k] - yi[k]) > abs(mean_theta[k] - mean_mu[k]),
                  info = paste("Observation", k, "shrinks toward mean with high sei"))
    }
  }

  # verify lambda formula against a closed-form calculation
  # lambda = tau^2 / (tau^2 + se^2)
  # theta = lambda * yi + (1 - lambda) * mu
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  lambda  <- tau_within^2 / (tau_within^2 + sei_mat^2)
  yi_mat  <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)
  expected_theta <- lambda * yi_mat + (1 - lambda) * mu_samples

  expect_equal(theta, expected_theta, tolerance = 1e-10)
})

test_that(".evaluate.brma.true_effects.norm returns BLUP means for same data", {

  S <- 25
  K <- 2

  mu_samples <- matrix(c(0.2, -0.1), nrow = S, ncol = K, byrow = TRUE)
  tau_within <- matrix(c(0.5, 0.25), nrow = S, ncol = K, byrow = TRUE)
  yi         <- c(1.0, -0.7)
  sei        <- c(0.2, 0.4)

  theta <- .evaluate.brma.true_effects.norm(
    mu_samples = mu_samples,
    tau_within = tau_within,
    yi         = yi,
    sei        = sei,
    same_data  = TRUE
  )

  lambda   <- tau_within[1, ]^2 / (tau_within[1, ]^2 + sei^2)
  expected <- mu_samples[1, ] + lambda * (yi - mu_samples[1, ])

  expect_equal(theta[1, ], expected, tolerance = 1e-12)
  expect_equal(apply(theta, 2, stats::var), c(0, 0), tolerance = 1e-14)
})

test_that(".evaluate.brma.true_effects.norm subtracts posterior-row bias offsets", {

  mu_samples <- matrix(
    c(
      0.0,  0.0,
      0.1, -0.1
    ),
    nrow = 2,
    byrow = TRUE
  )
  tau_within <- matrix(1, nrow = 2, ncol = 2)
  yi          <- c(0.4, -0.2)
  sei         <- c(1, 2)
  bias_offset <- matrix(
    c(
       0.2, -0.1,
      -0.3,  0.4
    ),
    nrow = 2,
    byrow = TRUE
  )

  theta <- .evaluate.brma.true_effects.norm(
    mu_samples  = mu_samples,
    tau_within  = tau_within,
    yi          = yi,
    sei         = sei,
    same_data   = TRUE,
    bias_offset = bias_offset
  )

  yi_samples     <- matrix(yi, nrow = 2, ncol = 2, byrow = TRUE)
  lambda         <- tau_within^2 / sweep(tau_within^2, 2, sei^2, "+")
  expected_theta <- mu_samples + lambda * (yi_samples - bias_offset - mu_samples)

  expect_equal(theta, expected_theta, tolerance = 1e-12)
})

test_that(".evaluate.brma.true_effects.norm samples from marginal with same_data = FALSE", {

  # test that same_data = FALSE samples from N(mu, tau) marginal distribution

  S <- 10000  # many samples for stable variance estimation
  K <- 3

  set.seed(9999)
  mu_samples <- matrix(0.4, nrow = S, ncol = K)  # constant mu
  tau_within <- matrix(0.25, nrow = S, ncol = K)  # constant tau

  # compute true effects with same_data = FALSE (marginal sampling)
  theta <- .evaluate.brma.true_effects.norm(
    mu_samples = mu_samples,
    tau_within = tau_within,
    yi         = NULL,  # not used when same_data = FALSE
    sei        = NULL,  # not used when same_data = FALSE
    same_data  = FALSE
  )

  # verify dimensions
  expect_equal(dim(theta), c(S, K))

  # verify mean is approximately mu
  for (k in seq_len(K)) {
    expect_equal(mean(theta[, k]), 0.4, tolerance = 0.02,
                 info = paste("Mean for observation", k))
  }

  # verify variance is approximately tau^2
  for (k in seq_len(K)) {
    expect_equal(var(theta[, k]), 0.25^2, tolerance = 0.005,
                 info = paste("Variance for observation", k))
  }
})

test_that(".outcome_rng.norm has correct sampling variance", {

  # create mock data
  S <- 10000  # many samples for stable variance estimation
  K <- 3

  set.seed(5678)
  mu_samples <- matrix(0.5, nrow = S, ncol = K)  # constant mu for variance test
  tau_within <- matrix(0.2, nrow = S, ncol = K)  # constant tau
  sei <- c(0.1, 0.2, 0.3)

  # compute response samples
  response <- .outcome_rng.norm(mu_samples, tau_within, sei)

  # verify dimensions
  expect_equal(dim(response), c(S, K))

  # verify sampling variance matches theoretical: Var = tau^2 + se^2
  for (k in seq_len(K)) {
    expected_var <- 0.2^2 + sei[k]^2
    observed_var <- var(response[, k])
    # use tolerance appropriate for 10k samples
    expect_equal(observed_var, expected_var, tolerance = 0.01,
                 info = paste("Variance for observation", k))
  }

  # verify means are close to mu
  for (k in seq_len(K)) {
    expect_equal(mean(response[, k]), 0.5, tolerance = 0.05,
                 info = paste("Mean for observation", k))
  }
})

test_that(".evaluate.brma.cluster_effects returns contribution matrix for new data", {

  S <- 10000  # many samples
  K <- 3

  set.seed(222)
  tau_between <- matrix(0.3, nrow = S, ncol = K)  # constant tau_between

  # for new data, gamma ~ N(0,1) is sampled per cluster, so contribution = gamma * tau_between
  # expected variance = tau_between^2 = 0.09
  contribution <- .evaluate.brma.cluster_effects(
    fit              = NULL,
    tau_between      = tau_between,
    cluster          = c(1, 1, 2),
    same_data        = FALSE,
    effect_direction = "positive"
  )

  expect_equal(dim(contribution), c(S, K))

  # variance should be approximately tau_between^2
  for (k in seq_len(K)) {
    expect_equal(var(contribution[, k]), 0.3^2, tolerance = 0.01,
                 info = paste("Variance for observation", k))
  }

  expect_equal(contribution[, 1], contribution[, 2])
})

test_that(".evaluate.brma.multilevel_blup.norm matches full covariance solve", {

  S       <- 8
  K       <- 5
  cluster <- c(1, 1, 2, 2, 2)

  set.seed(333)
  mu_samples  <- matrix(rnorm(S * K, mean = 0.2, sd = 0.1), nrow = S, ncol = K)
  tau_within  <- matrix(runif(S * K, min = 0.05, max = 0.25), nrow = S, ncol = K)
  tau_between <- matrix(runif(S * K, min = 0.05, max = 0.20), nrow = S, ncol = K)
  yi          <- c(0.10, 0.25, -0.10, 0.05, 0.30)
  vi          <- c(0.02, 0.03, 0.01, 0.04, 0.02)

  result <- .evaluate.brma.multilevel_blup.norm(
    mu_samples  = mu_samples,
    tau_within  = tau_within,
    tau_between = tau_between,
    yi          = yi,
    vi          = vi,
    cluster     = cluster
  )

  block_indices     <- .get_multilevel_block_indices(cluster)
  expected_cluster  <- matrix(0, nrow = S, ncol = K)
  expected_estimate <- matrix(0, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    covariance <- .build_multilevel_marginal_covariance(
      tau_within    = tau_within[s, ],
      tau_between   = tau_between[s, ],
      vi            = vi,
      block_indices = block_indices
    )
    weighted_residual <- as.vector(chol2inv(chol(covariance)) %*%
      (yi - mu_samples[s, ]))

    expected_estimate[s, ] <- tau_within[s, ]^2 * weighted_residual

    for (idx in block_indices) {
      cluster_scale <- sum(tau_between[s, idx] * weighted_residual[idx])
      expected_cluster[s, idx] <- tau_between[s, idx] * cluster_scale
    }
  }

  expect_equal(result[["cluster"]], expected_cluster, tolerance = 1e-12)
  expect_equal(result[["estimate"]], expected_estimate, tolerance = 1e-12)
})

test_that(".evaluate.brma.multilevel_blup.norm subtracts posterior-row bias offsets", {

  S       <- 3
  K       <- 4
  cluster <- c(1, 1, 2, 2)

  mu_samples <- matrix(
    c(
      0.00, 0.10, 0.20, 0.30,
      0.05, 0.15, 0.25, 0.35,
      0.10, 0.20, 0.30, 0.40
    ),
    nrow = S,
    byrow = TRUE
  )
  tau_within  <- matrix(0.20, nrow = S, ncol = K)
  tau_between <- matrix(0.15, nrow = S, ncol = K)
  yi          <- c(0.30, 0.10, 0.50, 0.20)
  vi          <- c(0.02, 0.03, 0.02, 0.04)
  bias_offset <- matrix(
    c(
       0.05,  0.00, 0.10, 0.00,
      -0.05,  0.03, 0.00, 0.08,
       0.02, -0.04, 0.05, 0.01
    ),
    nrow = S,
    byrow = TRUE
  )

  result <- .evaluate.brma.multilevel_blup.norm(
    mu_samples  = mu_samples,
    tau_within  = tau_within,
    tau_between = tau_between,
    yi          = yi,
    vi          = vi,
    cluster     = cluster,
    bias_offset = bias_offset
  )

  block_indices     <- .get_multilevel_block_indices(cluster)
  expected_cluster  <- matrix(0, nrow = S, ncol = K)
  expected_estimate <- matrix(0, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    covariance <- .build_multilevel_marginal_covariance(
      tau_within    = tau_within[s, ],
      tau_between   = tau_between[s, ],
      vi            = vi,
      block_indices = block_indices
    )
    residual <- yi - bias_offset[s, ] - mu_samples[s, ]
    weighted_residual <- as.vector(chol2inv(chol(covariance)) %*% residual)

    expected_estimate[s, ] <- tau_within[s, ]^2 * weighted_residual

    for (idx in block_indices) {
      cluster_scale <- sum(tau_between[s, idx] * weighted_residual[idx])
      expected_cluster[s, idx] <- tau_between[s, idx] * cluster_scale
    }
  }

  expect_equal(result[["cluster"]], expected_cluster, tolerance = 1e-12)
  expect_equal(result[["estimate"]], expected_estimate, tolerance = 1e-12)
})

test_that("variance ordering: mu < theta < response", {

  # this tests the theoretical relationship:
  # mu (fixed only) has least variance
  # theta (includes tau) has more variance
  # response (includes tau + se) has most variance

  S <- 5000
  K <- 4

  set.seed(444)
  # create realistic mock data
  mu_samples <- matrix(rnorm(S * K, mean = 0.33, sd = 0.05), nrow = S, ncol = K)
  tau_within <- matrix(abs(rnorm(S * K, mean = 0.25, sd = 0.05)), nrow = S, ncol = K)
  yi  <- c(0.30, 0.2, 0.5, 0.45)
  sei <- c(0.15, 0.1, 0.2, 0.12)

  # compute true effects and response
  theta <- .evaluate.brma.true_effects.norm(
    mu_samples = mu_samples,
    tau_within = tau_within,
    yi         = yi,
    sei        = sei,
    same_data  = TRUE
  )

  response <- .outcome_rng.norm(mu_samples, tau_within, sei)

  # compute variances across the distributions
  sd_mu       <- median(apply(mu_samples, 1, sd))
  sd_theta    <- median(apply(theta, 1, sd))
  sd_response <- median(apply(response, 1, sd))

  # theta variance should be less than or equal to mu variance
  # (shrinkage reduces variance toward the mean)
  # response variance should be larger than theta (adds sampling error)
  expect_true(sd_response > sd_theta,
              info = "response has more variance than theta")
  expect_true(sd_theta > sd_mu,
              info = "theta has more variance than terms")
})


# ============================================================================ #
# SECTION 2: Integration Tests with Pre-fitted Models
# ============================================================================ #
# These tests verify helper functions work correctly with actual brma fits
# ============================================================================ #

skip_if_no_fits()
fit_names     <- list_fits(tier = test_tier())
all_fit_names <- list_fits()
fits          <- lazy_fits(fit_names, validate = FALSE)
all_fits      <- lazy_fits(all_fit_names, validate = FALSE)

test_that(".evaluate.brma.tau returns correct structure", {

  for (name in names(fits)) {

    object <- fits[[name]]

    if (!inherits(object, "brma")) next

    priors       <- object[["priors"]]
    is_scale     <- .is_scale(object)
    is_multilevel <- .is_multilevel(object)
    K            <- nrow(object[["data"]][["outcome"]])

    # Direct tau extraction for scale formulas is exercised through predict().
    if (is_scale) next

    result <- .evaluate.brma.tau(
      fit           = object[["fit"]],
      scale_data    = object[["data"]][["scale"]],
      scale_formula = if (is_scale) attr(object[["data"]][["scale"]], "formula") else NULL,
      scale_priors  = priors[["scale"]],
      is_scale      = is_scale,
      is_multilevel = is_multilevel,
      K             = K
    )

    # verify structure
    expect_true(is.list(result), info = paste(name, ": returns list"))
    expect_true(all(c("tau_within", "tau_between") %in% names(result)),
                info = paste(name, ": has tau_within and tau_between components"))

    # verify dimensions
    S <- nrow(result[["tau_within"]])
    expect_equal(dim(result[["tau_within"]]), c(S, K), info = paste(name, ": tau_within dimensions"))
    expect_equal(dim(result[["tau_between"]]), c(S, K), info = paste(name, ": tau_between dimensions"))

    # verify positivity
    expect_true(all(result[["tau_within"]] >= 0), info = paste(name, ": tau_within is non-negative"))
    expect_true(all(result[["tau_between"]] >= 0), info = paste(name, ": tau_between is non-negative"))

    # verify relationship for non-multilevel: tau_between = 0
    if (!is_multilevel) {
      expect_true(all(result[["tau_between"]] == 0),
                  info = paste(name, ": non-multilevel has tau_between = 0"))
    }

    # verify total tau can be reconstructed: tau = sqrt(tau_within^2 + tau_between^2)
    tau_reconstructed <- sqrt(result[["tau_within"]]^2 + result[["tau_between"]]^2)
    expect_true(all(tau_reconstructed >= 0), info = paste(name, ": reconstructed tau is non-negative"))
  }
})

test_that(".evaluate.brma.mu returns correct dimensions", {

  for (name in names(fits)) {

    object <- fits[[name]]

    if (!inherits(object, "brma")) next

    priors           <- object[["priors"]]
    is_mods          <- .is_mods(object)
    is_PET           <- .is_PET(object)
    is_PEESE         <- .is_PEESE(object)
    effect_direction <- .effect_direction(object)
    outcome_data     <- object[["data"]][["outcome"]]
    K                <- nrow(outcome_data)

    # Direct moderator design handling is exercised through predict().
    if (is_mods) next

    mu_samples <- .evaluate.brma.mu(
      fit               = object[["fit"]],
      outcome_data      = outcome_data,
      mods_data         = object[["data"]][["mods"]],
      mods_formula      = if (is_mods) attr(object[["data"]][["mods"]], "formula") else NULL,
      mods_priors       = priors[["mods"]],
      is_mods           = is_mods,
      is_PET            = is_PET,
      is_PEESE          = is_PEESE,
      effect_direction  = effect_direction,
      bias_adjusted     = TRUE,
      K                 = K
    )

    # verify dimensions
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    S <- nrow(posterior_samples)
    expect_equal(dim(mu_samples), c(S, K),
                 info = paste(name, ": mu dimensions are S x K"))
  }
})

test_that("PET/PEESE bias offsets match column-wise algebra", {

  S   <- 100
  K   <- 5
  sei <- c(0.1, 0.2, 0.15, 0.25, 0.12)

  set.seed(666)
  mu_samples  <- matrix(rnorm(S * K), nrow = S, ncol = K)
  PET_samples <- rnorm(S, mean = 0.5, sd = 0.1)

  mu_loop <- mu_samples
  for (i in seq_len(K)) {
    mu_loop[, i] <- mu_loop[, i] + PET_samples * sei[i]
  }

  mu_vec <- mu_samples + outer(PET_samples, sei)

  expect_equal(mu_loop, mu_vec, tolerance = 1e-14,
               info = "Vectorized outer() matches loop for PET")

  # same test for PEESE (se^2)
  set.seed(777)
  mu_samples    <- matrix(rnorm(S * K), nrow = S, ncol = K)
  PEESE_samples <- rnorm(S, mean = 0.3, sd = 0.05)
  sei_sq        <- sei^2

  mu_loop <- mu_samples
  for (i in seq_len(K)) {
    mu_loop[, i] <- mu_loop[, i] + PEESE_samples * sei_sq[i]
  }

  mu_vec <- mu_samples + outer(PEESE_samples, sei_sq)

  expect_equal(mu_loop, mu_vec, tolerance = 1e-14,
               info = "Vectorized outer() matches loop for PEESE")
})

test_that(".evaluate.brma.bias_offset handles PET/PEESE and effect direction", {

  posterior_samples <- matrix(
    c(
      0.50, 1.00,
      0.25, 2.00
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(posterior_samples) <- c("PET", "PEESE")

  outcome_data <- data.frame(sei = c(0.10, 0.20, 0.30))
  expected_positive <- outer(posterior_samples[, "PET"], outcome_data[["sei"]]) +
    outer(posterior_samples[, "PEESE"], outcome_data[["sei"]]^2)

  offset_positive <- .evaluate.brma.bias_offset(
    fit               = NULL,
    outcome_data      = outcome_data,
    is_PET            = TRUE,
    is_PEESE          = TRUE,
    effect_direction  = "positive",
    K                 = 3,
    posterior_samples = posterior_samples
  )
  offset_negative <- .evaluate.brma.bias_offset(
    fit               = NULL,
    outcome_data      = outcome_data,
    is_PET            = TRUE,
    is_PEESE          = TRUE,
    effect_direction  = "negative",
    K                 = 3,
    posterior_samples = posterior_samples
  )

  expect_equal(offset_positive, expected_positive, tolerance = 1e-12)
  expect_equal(offset_negative, -expected_positive, tolerance = 1e-12)
})

test_that("rho clamping handles boundary values", {

  # test that rho values outside [0, 1] are properly clamped

  # create mock rho with edge cases
  rho <- c(-0.01, 0, 0.5, 1, 1.001)
  rho_clamped <- pmin(pmax(rho, 0), 1)

  expect_equal(rho_clamped, c(0, 0, 0.5, 1, 1))

  # verify tau decomposition is valid after clamping
  tau <- 0.3
  tau_within  <- tau * sqrt(1 - rho_clamped)
  tau_between <- tau * sqrt(rho_clamped)

  # all values should be non-negative
  expect_true(all(tau_within >= 0))
  expect_true(all(tau_between >= 0))

  # verify Pythagorean relationship: tau^2 = tau_within^2 + tau_between^2
  tau_reconstructed <- sqrt(tau_within^2 + tau_between^2)
  expect_equal(tau_reconstructed, rep(tau, length(rho)), tolerance = 1e-10)
})

test_that("GLMM posterior extraction helpers are vectorized", {

  posterior_samples <- matrix(
    c(
      0.25, 0.75, -1.0, 1.0, 0.10, -0.20,
      0.40, 0.60, -0.5, 0.5, 0.30,  0.40
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(posterior_samples) <- c(
    "pi[1]", "pi[2]", "phi[1]", "phi[2]", "theta[1]", "theta[2]"
  )

  tau_within <- matrix(c(0.2, 0.3, 0.4, 0.5), nrow = 2, byrow = TRUE)

  expect_equal(
    .evaluate.brma.baserate(fit = NULL, K = 2, posterior_samples = posterior_samples),
    stats::qlogis(posterior_samples[, c("pi[1]", "pi[2]")])
  )
  expect_equal(
    .evaluate.brma.lograte(fit = NULL, K = 2, posterior_samples = posterior_samples),
    posterior_samples[, c("phi[1]", "phi[2]")]
  )
  expect_equal(
    .evaluate.brma.theta.glmm(
      fit               = NULL,
      tau_within        = tau_within,
      same_data         = TRUE,
      K                 = 2,
      posterior_samples = posterior_samples
    ),
    posterior_samples[, c("theta[1]", "theta[2]")] * tau_within
  )
})

test_that("matrix replication patterns preserve dimensions", {

  # verify that matrix(vec, S, K, byrow = TRUE) works as expected
  # this pattern is used throughout the helpers

  S   <- 3
  K   <- 4
  vec <- c(0.1, 0.2, 0.3, 0.4)

  # byrow = TRUE: each row is the vector
  mat <- matrix(vec, nrow = S, ncol = K, byrow = TRUE)

  expect_equal(nrow(mat), S)
  expect_equal(ncol(mat), K)

  # each row should equal vec
  for (s in seq_len(S)) {
    expect_equal(mat[s, ], vec)
  }

  # byrow = FALSE: each column is the vector (recycled)
  vec2 <- c(1, 2, 3)
  mat2 <- matrix(vec2, nrow = S, ncol = K, byrow = FALSE)

  # each column should equal vec2
  for (k in seq_len(K)) {
    expect_equal(mat2[, k], vec2)
  }
})

# ============================================================================ #
# SECTION 3: Unit Tests for Aggregated Predictions (newdata = TRUE)
# ============================================================================ #
# These tests verify the aggregation logic for predict.brma with newdata = TRUE
# ============================================================================ #

test_that("aggregation matches rowMeans", {

  # test that rowMeans aggregation works as expected
  S <- 100
  K <- 5

  set.seed(888)
  mu_samples <- matrix(rnorm(S * K, mean = 0.3, sd = 0.1), nrow = S, ncol = K)

  # aggregation: rowMeans across K observations

  mu_aggregated <- matrix(rowMeans(mu_samples), ncol = 1)

  # verify dimensions
  expect_equal(dim(mu_aggregated), c(S, 1))

  # verify mean is preserved (on average)
  expect_equal(mean(mu_aggregated), mean(mu_samples), tolerance = 0.01)

  # verify aggregation against a direct matrix calculation
  expected_aggregated <- matrix(apply(mu_samples, 1, mean), ncol = 1)
  expect_equal(mu_aggregated, expected_aggregated, tolerance = 1e-14)
})

test_that("aggregation is no-op for identical columns (non-mods/scale models)", {

  # for models without moderators/scale, all columns are identical
  # aggregation should return the same value
  S <- 100
  K <- 5

  # create mock data where all columns are identical
  mu_values <- rnorm(S, mean = 0.5, sd = 0.1)
  mu_samples <- matrix(mu_values, nrow = S, ncol = K)

  # aggregation
  mu_aggregated <- matrix(rowMeans(mu_samples), ncol = 1)

  # should be identical to original column
  expect_equal(mu_aggregated[, 1], mu_values, tolerance = 1e-14)
})

test_that("aggregated true effects have expected variance", {

  # for aggregated effect predictions:
  # theta ~ N(mu, tau) where mu and tau are aggregated
  # variance of theta should be approximately var(mu) + mean(tau)^2

  S <- 10000  # many samples for stable variance estimation

  set.seed(999)
  mu_aggregated  <- matrix(rnorm(S, mean = 0.4, sd = 0.05), ncol = 1)
  tau_aggregated <- matrix(abs(rnorm(S, mean = 0.25, sd = 0.02)), ncol = 1)

  # sample true effects: theta = mu + rnorm(S) * tau
  theta_samples <- mu_aggregated + rnorm(S) * tau_aggregated

  # expected variance: Var(mu) + E[tau^2] = Var(mu) + Var(tau) + E[tau]^2
  expected_var <- var(mu_aggregated) + var(tau_aggregated) + mean(tau_aggregated)^2

  # observed variance should be close
  observed_var <- var(theta_samples)
  expect_equal(observed_var, expected_var, tolerance = 0.02,
               info = "Aggregated theta variance matches theory")
})


# ============================================================================ #
# SECTION 4: Integration Tests for predict.brma with newdata = TRUE
# ============================================================================ #
# These tests verify predict.brma with aggregated predictions using cached fits
# ============================================================================ #

skip_if_no_fits()

test_that("predict.brma with newdata = TRUE returns single aggregated prediction", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects
    if (!inherits(object, "brma")) next

    # test type = "terms"
    result_terms <- predict(object, newdata = TRUE, type = "terms")

    expect_s3_class(result_terms, "brma_samples")
    expect_null(attr(result_terms, "data"), info = paste(name, ": data is NULL for aggregate"))
    expect_equal(nrow(summary(result_terms)), 1,
                 info = paste(name, ": has single row for aggregate terms"))

    # test type = "terms.scale"
    result_scale <- predict(object, newdata = TRUE, type = "terms.scale")

    expect_s3_class(result_scale, "brma_samples")
    expect_null(attr(result_scale, "data"), info = paste(name, ": data is NULL for aggregate scale"))
    expect_equal(nrow(summary(result_scale)), 1,
                 info = paste(name, ": has single row for aggregate scale"))

    # test type = "effect"
    result_effect <- predict(object, newdata = TRUE, type = "effect")

    expect_s3_class(result_effect, "brma_samples")
    expect_null(attr(result_effect, "data"), info = paste(name, ": data is NULL for aggregate effect"))
    expect_equal(nrow(summary(result_effect)), 1,
                 info = paste(name, ": has single row for aggregate effect"))
  }
})

test_that("predict.brma with newdata = TRUE returns S x 1 matrix", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects
    if (!inherits(object, "brma")) next

    # get expected sample count
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    S <- nrow(posterior_samples)

    # test type = "terms"
    samples_terms <- predict(object, newdata = TRUE, type = "terms")
    expect_equal(dim(samples_terms), c(S, 1),
                 info = paste(name, ": terms samples are S x 1"))
    expect_equal(colnames(samples_terms), "mu",
                 info = paste(name, ": terms samples have 'mu' column name"))

    # test type = "terms.scale"
    samples_scale <- predict(object, newdata = TRUE, type = "terms.scale")
    expect_equal(dim(samples_scale), c(S, 1),
                 info = paste(name, ": scale samples are S x 1"))
    expect_equal(colnames(samples_scale), "tau",
                 info = paste(name, ": scale samples have 'tau' column name"))

    # test type = "effect"
    samples_effect <- predict(object, newdata = TRUE, type = "effect")
    expect_equal(dim(samples_effect), c(S, 1),
                 info = paste(name, ": effect samples are S x 1"))
    expect_equal(colnames(samples_effect), "theta",
                 info = paste(name, ": effect samples have 'theta' column name"))
  }
})

test_that("predict.brma rejects response-scale newdata predictions", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects
    if (!inherits(object, "brma")) next

    # should throw error for type = "response"
    expect_error(
      predict(object, newdata = TRUE, type = "response"),
      "Aggregated predictions.*not available for type = 'response'",
      info = paste(name, ": aggregate + response is rejected")
    )
  }
})

test_that("aggregated mu equals rowMeans of non-aggregated mu", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects
    if (!inherits(object, "brma")) next

    # get non-aggregated samples (as plain matrix)
    samples_full <- as.matrix(predict(object, newdata = NULL, type = "terms"))

    # get aggregated samples (as plain matrix)
    samples_agg <- as.matrix(predict(object, newdata = TRUE, type = "terms"))

    # aggregated should equal rowMeans of full
    expected_agg <- matrix(rowMeans(samples_full), ncol = 1)
    colnames(expected_agg) <- "mu"

    expect_equal(samples_agg, expected_agg, tolerance = 1e-10,
                 info = paste(name, ": aggregated mu equals rowMeans of full mu"))
  }
})

test_that("aggregated tau equals rowMeans of non-aggregated tau", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects
    if (!inherits(object, "brma")) next

    # get non-aggregated samples (as plain matrix)
    samples_full <- as.matrix(predict(object, newdata = NULL, type = "terms.scale"))

    # get aggregated samples (as plain matrix)
    samples_agg <- as.matrix(predict(object, newdata = TRUE, type = "terms.scale"))

    # aggregated should equal rowMeans of full
    expected_agg <- matrix(rowMeans(samples_full), ncol = 1)
    colnames(expected_agg) <- "tau"

    expect_equal(samples_agg, expected_agg, tolerance = 1e-10,
                 info = paste(name, ": aggregated tau equals rowMeans of full tau"))
  }
})


# ============================================================================ #
# SECTION 3: Tests for .extract_use_normal()
# ============================================================================ #
# These tests verify the bias indicator extraction helper function
# that identifies which posterior samples use weightfunction vs normal path
# ============================================================================ #

test_that(".extract_omega_samples orders indexed omega columns", {

  posterior_samples <- matrix(
    seq_len(12),
    nrow = 3,
    ncol = 4
  )
  colnames(posterior_samples) <- c("omega[10]", "mu", "omega[2]", "omega[1]")

  omega <- .extract_omega_samples(posterior_samples)

  expect_equal(colnames(omega), c("omega[1]", "omega[2]", "omega[10]"))
  expect_equal(omega[, 1], posterior_samples[, "omega[1]"])
  expect_equal(omega[, 2], posterior_samples[, "omega[2]"])
  expect_equal(omega[, 3], posterior_samples[, "omega[10]"])
})

test_that(".extract_use_normal returns correct structure for brma without bias", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects or objects with bias priors
    if (!inherits(object, "brma")) next
    if (.is_bias(object)) next

    # test the function
    use_normal <- .extract_use_normal(object)

    # get expected sample count
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    S <- nrow(posterior_samples)

    # verify structure
    expect_type(use_normal, "logical")
    expect_length(use_normal, S)

    # all samples should use normal path (no weightfunction)
    expect_true(all(use_normal),
                info = paste(name, ": brma without bias has all use_normal = TRUE"))
  }
})

test_that(".extract_use_normal returns correct structure for brma with weightfunction", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects or objects without weightfunction priors
    if (!inherits(object, "brma")) next
    if (!.is_weightfunction(object)) next
    if (inherits(object, "RoBMA")) next  # RoBMA has mixture priors, test separately

    # test the function
    use_normal <- .extract_use_normal(object)

    # get expected sample count
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    S <- nrow(posterior_samples)

    # verify structure
    expect_type(use_normal, "logical")
    expect_length(use_normal, S)

    # all samples should use weighted path (all from weightfunction)
    expect_true(all(!use_normal),
                info = paste(name, ": brma with single weightfunction has all use_normal = FALSE"))
  }
})

test_that(".extract_use_normal returns correct structure for brma with PET", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects or objects without PET priors
    if (!inherits(object, "brma")) next
    if (!.is_PET(object)) next
    if (.is_weightfunction(object)) next  # skip if also has weightfunction
    if (inherits(object, "RoBMA")) next

    # test the function
    use_normal <- .extract_use_normal(object)

    # get expected sample count
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    S <- nrow(posterior_samples)

    # verify structure
    expect_type(use_normal, "logical")
    expect_length(use_normal, S)

    # all samples should use normal path (PET is not a weightfunction)
    expect_true(all(use_normal),
                info = paste(name, ": brma with PET has all use_normal = TRUE"))
  }
})

test_that(".extract_use_normal returns correct structure for RoBMA", {

  robma_names <- intersect(names(fits), catalog_fits(class = "RoBMA", tier = test_tier()))
  if (length(robma_names) == 0L) {
    skip("No cached RoBMA fits available for use_normal branch checks.")
  }

  for (name in robma_names) {

    object <- fits[[name]]

    # test the function
    use_normal <- .extract_use_normal(object)

    # get expected sample count
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    S <- nrow(posterior_samples)

    # verify structure
    expect_type(use_normal, "logical")
    expect_length(use_normal, S)

    # RoBMA should have a mix of TRUE and FALSE (different bias models)
    # get bias_indicator to verify
    if ("bias_indicator" %in% colnames(posterior_samples)) {
      bias_indicator <- as.integer(posterior_samples[, "bias_indicator"])

      # extract bias priors
      priors_bias <- object[["priors"]][["outcome"]][["bias"]]
      if (!BayesTools::is.prior.mixture(priors_bias)) {
        priors_bias <- list(priors_bias)
      }

      # identify selected-normal kernel indices
      wf_indices <- which(sapply(priors_bias, .prior_is_selection_kernel))

      # verify use_normal is correctly computed
      expected_use_normal <- !(bias_indicator %in% wf_indices)
      expect_equal(use_normal, expected_use_normal,
                   info = paste(name, ": use_normal matches bias_indicator logic"))

      # verify we have both types of samples (unless all priors are one type)
      if (length(wf_indices) > 0 && length(wf_indices) < length(priors_bias)) {
        expect_true(any(use_normal) && any(!use_normal),
                    info = paste(name, ": RoBMA has mix of normal and weighted samples"))
      }
    }
  }
})

test_that("RoBMA mixed-bias branch PDF, CDF, and RNG use selected-normal branches", {

  name <- "dat.lehmann2018_RoBMA"
  skip_if_missing_fits(name)

  object            <- fits[[name]]
  setup             <- .estimate_likelihood_setup.brma(object)
  posterior_samples <- setup[["posterior_samples"]]
  bias_indicator    <- .extract_bias_indicator(object, posterior_samples = posterior_samples)
  use_normal        <- .extract_use_normal(object, posterior_samples = posterior_samples)
  fit_data          <- setup[["fit_data"]]

  normal_rows   <- head(which(use_normal), 100)
  weighted_rows <- head(which(!use_normal), 100)
  if (length(normal_rows) == 0L || length(weighted_rows) == 0L) {
    skip("Cached RoBMA fit lacks both normal and weighted publication-bias branches.")
  }

  selected_rows <- sort(c(normal_rows, weighted_rows))
  selection_context <- .selection_context(
    object            = object,
    posterior_samples = posterior_samples
  )
  selected_context <- .selection_context_subset_rows(selection_context, selected_rows)

  log_lik <- .outcome_pdf.selnorm(
    yi                = setup[["yi"]],
    mu_samples        = setup[["mu"]][selected_rows, , drop = FALSE],
    tau_within        = setup[["tau_within"]][selected_rows, , drop = FALSE],
    sei               = setup[["sei"]],
    selection_context = selected_context
  )
  cdf_vals <- .outcome_cdf.selnorm(
    yi                = setup[["yi"]],
    mu_samples        = setup[["mu"]][selected_rows, , drop = FALSE],
    tau_within        = setup[["tau_within"]][selected_rows, , drop = FALSE],
    sei               = setup[["sei"]],
    selection_context = selected_context
  )

  expect_equal(dim(log_lik), c(length(selected_rows), setup[["K"]]))
  expect_equal(dim(cdf_vals), c(length(selected_rows), setup[["K"]]))
  expect_true(all(is.finite(log_lik)))
  expect_true(all(is.finite(cdf_vals)))
  expect_true(all(cdf_vals > 0 & cdf_vals < 1))

  normal_branch   <- sort(unique(bias_indicator[normal_rows]))[1]
  weighted_branch <- sort(unique(bias_indicator[weighted_rows]))[1]
  normal_branch_rows   <- normal_rows[bias_indicator[normal_rows] == normal_branch]
  weighted_branch_rows <- weighted_rows[bias_indicator[weighted_rows] == weighted_branch]

  normal_pdf <- .outcome_pdf.norm(
    yi         = if (setup[["effect_direction"]] == "negative") -setup[["yi"]] else setup[["yi"]],
    mu_samples = if (setup[["effect_direction"]] == "negative") {
      -setup[["mu"]][normal_branch_rows, , drop = FALSE]
    } else {
      setup[["mu"]][normal_branch_rows, , drop = FALSE]
    },
    tau_within = setup[["tau_within"]][normal_branch_rows, , drop = FALSE],
    sei        = setup[["sei"]]
  )
  normal_cdf <- .outcome_cdf.norm(
    yi         = if (setup[["effect_direction"]] == "negative") -setup[["yi"]] else setup[["yi"]],
    mu_samples = if (setup[["effect_direction"]] == "negative") {
      -setup[["mu"]][normal_branch_rows, , drop = FALSE]
    } else {
      setup[["mu"]][normal_branch_rows, , drop = FALSE]
    },
    tau_within = setup[["tau_within"]][normal_branch_rows, , drop = FALSE],
    sei        = setup[["sei"]],
    lower.tail = if (setup[["effect_direction"]] == "negative") FALSE else TRUE
  )

  expect_equal(
    log_lik[match(normal_branch_rows, selected_rows), , drop = FALSE],
    normal_pdf,
    tolerance = 1e-10
  )
  expect_equal(
    cdf_vals[match(normal_branch_rows, selected_rows), , drop = FALSE],
    normal_cdf,
    tolerance = 1e-10
  )

  weighted_pdf_norm <- .outcome_pdf.norm(
    yi         = if (setup[["effect_direction"]] == "negative") -setup[["yi"]] else setup[["yi"]],
    mu_samples = if (setup[["effect_direction"]] == "negative") {
      -setup[["mu"]][weighted_branch_rows, , drop = FALSE]
    } else {
      setup[["mu"]][weighted_branch_rows, , drop = FALSE]
    },
    tau_within = setup[["tau_within"]][weighted_branch_rows, , drop = FALSE],
    sei        = setup[["sei"]]
  )
  expect_false(isTRUE(all.equal(
    log_lik[match(weighted_branch_rows, selected_rows), , drop = FALSE],
    weighted_pdf_norm,
    tolerance = 1e-8
  )))

  set.seed(11)
  rng_normal <- .outcome_rng.selnorm(
    mu_samples        = setup[["mu"]][normal_branch_rows, , drop = FALSE],
    tau_within        = setup[["tau_within"]][normal_branch_rows, , drop = FALSE],
    sei               = setup[["sei"]],
    selection_context = .selection_context_subset_rows(selection_context, normal_branch_rows)
  )
  set.seed(11)
  rng_normal_expected <- .outcome_rng.norm(
    mu_samples = setup[["mu"]][normal_branch_rows, , drop = FALSE],
    tau_within = setup[["tau_within"]][normal_branch_rows, , drop = FALSE],
    sei        = setup[["sei"]]
  )
  expect_equal(dim(rng_normal), dim(rng_normal_expected))
  expect_equal(as.vector(rng_normal), as.vector(rng_normal_expected),
               tolerance = 1e-12)

  set.seed(12)
  rng_mixed <- .outcome_rng.selnorm(
    mu_samples        = setup[["mu"]][selected_rows, , drop = FALSE],
    tau_within        = setup[["tau_within"]][selected_rows, , drop = FALSE],
    sei               = setup[["sei"]],
    selection_context = selected_context
  )

  expect_equal(dim(rng_mixed), c(length(selected_rows), setup[["K"]]))
  expect_true(all(is.finite(rng_mixed)))
})

.expect_indicator_for_prior <- function(posterior_samples, prior, column, info) {

  if (is.null(prior) || !BayesTools::is.prior.mixture(prior)) {
    return(invisible(TRUE))
  }

  expect_true(column %in% colnames(posterior_samples),
              info = paste(info, "missing posterior indicator", column))
  if (!column %in% colnames(posterior_samples)) {
    return(invisible(FALSE))
  }

  values <- posterior_samples[, column]

  expect_true(all(is.finite(values)),
              info = paste(info, column, "contains non-finite values"))
  expect_true(all(values == as.integer(values)),
              info = paste(info, column, "contains non-integer values"))
  expect_true(all(values >= 1 & values <= length(prior)),
              info = paste(info, column, "outside prior component range"))
}

test_that("product-space posterior indicators have expected columns and valid ranges", {

  product_names <- c(
    "dat.lehmann2018_BMA.norm",
    "dat.lehmann2018_BMA.norm_mods",
    "dat.lehmann2018_BMA.norm_scale",
    "bcg_BMA.glmm",
    "bcg_BMA.glmm_3lvl_location_scale",
    "dat.lehmann2018_RoBMA",
    "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  )
  skip_if_missing_fits(product_names)

  for (name in product_names) {

    object            <- all_fits[[name]]
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    priors            <- object[["priors"]]

    .expect_indicator_for_prior(
      posterior_samples = posterior_samples,
      prior             = priors[["outcome"]][["mu"]],
      column            = "mu_indicator",
      info              = name
    )
    .expect_indicator_for_prior(
      posterior_samples = posterior_samples,
      prior             = priors[["outcome"]][["tau"]],
      column            = "tau_indicator",
      info              = name
    )
    .expect_indicator_for_prior(
      posterior_samples = posterior_samples,
      prior             = priors[["outcome"]][["rho"]],
      column            = "rho_indicator",
      info              = name
    )
    .expect_indicator_for_prior(
      posterior_samples = posterior_samples,
      prior             = priors[["outcome"]][["bias"]],
      column            = "bias_indicator",
      info              = name
    )

    for (term in names(priors[["mods"]])) {
      column <- paste0(
        BayesTools::JAGS_parameter_names(term, formula_parameter = "mu"),
        "_indicator"
      )
      .expect_indicator_for_prior(
        posterior_samples = posterior_samples,
        prior             = priors[["mods"]][[term]],
        column            = column,
        info              = name
      )
    }

    for (term in names(priors[["scale"]])) {
      column <- paste0(
        BayesTools::JAGS_parameter_names(term, formula_parameter = "log_tau"),
        "_indicator"
      )
      .expect_indicator_for_prior(
        posterior_samples = posterior_samples,
        prior             = priors[["scale"]][[term]],
        column            = column,
        info              = name
      )
    }
  }
})

test_that("RoBMA inactive bias branches expose zeros or neutral weights", {

  product_names <- c("dat.lehmann2018_RoBMA", "dat.lehmann2018_RoBMA_3lvl_mods_scale")
  skip_if_missing_fits(product_names)

  for (name in product_names) {

    object            <- all_fits[[name]]
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    priors_bias       <- object[["priors"]][["outcome"]][["bias"]]
    if (is.null(priors_bias) || !"bias_indicator" %in% colnames(posterior_samples)) {
      next
    }
    if (!BayesTools::is.prior.mixture(priors_bias)) {
      priors_bias <- list(priors_bias)
    }

    bias_indicator      <- as.integer(posterior_samples[, "bias_indicator"])
    branch_is_PET       <- vapply(priors_bias, BayesTools::is.prior.PET, logical(1))
    branch_is_PEESE     <- vapply(priors_bias, BayesTools::is.prior.PEESE, logical(1))
    branch_is_weightfun <- vapply(priors_bias, BayesTools::is.prior.weightfunction, logical(1))

    if ("PET" %in% colnames(posterior_samples)) {
      inactive <- which(!branch_is_PET[bias_indicator])
      expect_equal(as.numeric(posterior_samples[inactive, "PET"]), rep(0, length(inactive)),
                   tolerance = 1e-12, info = paste(name, "inactive PET branch"))
    }
    if ("PEESE" %in% colnames(posterior_samples)) {
      inactive <- which(!branch_is_PEESE[bias_indicator])
      expect_equal(as.numeric(posterior_samples[inactive, "PEESE"]), rep(0, length(inactive)),
                   tolerance = 1e-12, info = paste(name, "inactive PEESE branch"))
    }

    omega_cols <- grep("^omega(\\[|$)", colnames(posterior_samples), value = TRUE)
    if (length(omega_cols) > 0L) {
      inactive <- which(!branch_is_weightfun[bias_indicator])
      omega_inactive <- as.matrix(posterior_samples[inactive, omega_cols, drop = FALSE])
      expect_equal(dim(omega_inactive), c(length(inactive), length(omega_cols)),
                   info = paste(name, "inactive weightfunction branch"))
      expect_equal(as.vector(omega_inactive), rep(1, length(omega_inactive)),
                   tolerance = 1e-12,
                   info = paste(name, "inactive weightfunction branch"))
    }
  }
})


# ============================================================================ #
# SECTION 4: Integration Tests for .pdf.brma() and .cdf.brma() with use_normal
# ============================================================================ #
# These tests verify that PDF and CDF functions work correctly with the
# use_normal fast-path optimization
# ============================================================================ #

test_that(".pdf.brma returns finite log-likelihoods for weightfunction models", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects or objects without weightfunction priors
    if (!inherits(object, "brma")) next
    if (!.is_weightfunction(object)) next

    # compute PDF (this internally uses use_normal optimization)
    log_lik <- .pdf.brma(object)

    # verify structure
    expect_true(is.matrix(log_lik),
                info = paste(name, ": log_lik is a matrix"))
    expect_true(all(is.finite(log_lik)),
                info = paste(name, ": all log-likelihoods are finite"))

    # verify dimensions match
    K <- length(.outcome_data_yi(object))
    expect_equal(ncol(log_lik), K,
                 info = paste(name, ": ncol matches number of observations"))
  }
})

test_that(".cdf.brma returns valid CDF values for weightfunction models", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects or objects without weightfunction priors
    if (!inherits(object, "brma")) next
    if (!.is_weightfunction(object)) next

    # compute CDF (this internally uses use_normal optimization)
    cdf_vals <- .cdf.brma(object)

    # verify structure
    expect_true(is.matrix(cdf_vals),
                info = paste(name, ": cdf_vals is a matrix"))
    expect_true(all(is.finite(cdf_vals)),
                info = paste(name, ": all CDF values are finite"))
    expect_true(all(cdf_vals > 0 & cdf_vals < 1),
                info = paste(name, ": all CDF values are in (0, 1)"))

    # verify dimensions match
    K <- length(.outcome_data_yi(object))
    expect_equal(ncol(cdf_vals), K,
                 info = paste(name, ": ncol matches number of observations"))
  }
})
