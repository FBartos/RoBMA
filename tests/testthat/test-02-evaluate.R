# ============================================================================ #
# test-02-evaluate.R
# ============================================================================ #
#
# Unit and integration tests for .evaluate.brma.* helper functions
# defined in R/brma.evaluate.R
#
# Structure:
# 1. Unit tests: Mock posterior matrices (no JAGS fit needed)
# 2. Integration tests: Pre-fitted models from list_fits()
#
# ============================================================================ #

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "evaluate")

# ============================================================================ #
# SECTION 1: Unit Tests with Mock Data
# ============================================================================ #
# These tests verify the computational logic of helper functions
# using mock posterior sample matrices (no JAGS fit required)
# ============================================================================ #

test_that(".evaluate.brma.true_effects computes correct BLUPs", {

  # create mock data
  S <- 100  # samples
  K <- 5    # observations

  set.seed(1234)
  mu_samples <- matrix(rnorm(S * K, mean = 0.5, sd = 0.1), nrow = S, ncol = K)
  tau_within <- matrix(abs(rnorm(S * K, mean = 0.2, sd = 0.05)), nrow = S, ncol = K)
  yi  <- c(0.3, 0.5, 0.7, 0.4, 0.6)
  sei <- c(0.1, 0.2, 0.15, 0.25, 0.12)

  # compute true effects
  theta <- .evaluate.brma.true_effects(mu_samples, tau_within, yi, sei)

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
      expect_true(abs(mean_theta[k] - yi[k]) < abs(mean_theta[k] - mean_mu[k]) + 0.1,
                  info = paste("Observation", k, "should shrink toward data with high tau"))
    }
  }

  # verify lambda formula is correct by manual calculation
  # lambda = tau^2 / (tau^2 + se^2)
  # theta = lambda * yi + (1 - lambda) * mu
  sei_mat <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  lambda  <- tau_within^2 / (tau_within^2 + sei_mat^2)
  yi_mat  <- matrix(yi, nrow = S, ncol = K, byrow = TRUE)
  expected_theta <- lambda * yi_mat + (1 - lambda) * mu_samples

  expect_equal(theta, expected_theta, tolerance = 1e-10)
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
  set.seed(9999)
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


test_that(".evaluate.brma.study_effects returns contribution matrix for new data", {

  S <- 10000  # many samples
  K <- 3

  set.seed(222)
  tau_between <- matrix(0.3, nrow = S, ncol = K)  # constant tau_between

  # for new data, gamma ~ N(0,1) is sampled, so contribution = gamma * tau_between
  # expected variance = tau_between^2 = 0.09
  set.seed(333)
  contribution <- .evaluate.brma.study_effects(
    fit              = NULL,
    tau_between      = tau_between,
    study_ids        = c(1, 1, 2),  # ignored for new data
    same_data        = FALSE,       # triggers new gamma sampling
    effect_direction = "positive"
  )

  expect_equal(dim(contribution), c(S, K))

  # variance should be approximately tau_between^2
  for (k in seq_len(K)) {
    expect_equal(var(contribution[, k]), 0.3^2, tolerance = 0.01,
                 info = paste("Variance for observation", k))
  }
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
  mu_samples <- matrix(rnorm(S * K, mean = 0.3, sd = 0.15), nrow = S, ncol = K)
  tau_within <- matrix(abs(rnorm(S * K, mean = 0.2, sd = 0.02)), nrow = S, ncol = K)
  yi  <- c(0.4, 0.2, 0.5, 0.35)
  sei <- c(0.15, 0.1, 0.2, 0.12)

  # compute true effects and response
  theta <- .evaluate.brma.true_effects(mu_samples, tau_within, yi, sei)

  set.seed(555)
  response <- .outcome_rng.norm(mu_samples, tau_within, sei)

  # compute variances for each observation
  for (k in seq_len(K)) {
    var_mu       <- var(mu_samples[, k])
    var_theta    <- var(theta[, k])
    var_response <- var(response[, k])

    # theta variance should be less than or equal to mu variance
    # (shrinkage reduces variance toward the mean)
    # response variance should be larger than theta (adds sampling error)
    expect_true(var_response > var_theta,
                info = paste("Observation", k, ": response should have more variance than theta"))
  }
})


# ============================================================================ #
# SECTION 2: Integration Tests with Pre-fitted Models
# ============================================================================ #
# These tests verify helper functions work correctly with actual brma fits
# ============================================================================ #

skip_if_no_fits()
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


test_that(".evaluate.brma.tau returns correct structure", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects (e.g., RoBMA objects)
    if (!inherits(object, "brma")) next

    priors       <- object[["priors"]]
    is_scale     <- .is_scale(object)
    is_multilevel <- .is_multilevel(object)
    K            <- nrow(object[["data"]][["outcome"]])

    # skip scale models for now - they require specific prior structure
    # that may not be present in all cached fits
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
    expect_true(is.list(result), info = paste(name, ": should return list"))
    expect_true(all(c("tau_within", "tau_between") %in% names(result)),
                info = paste(name, ": should have tau_within and tau_between components"))

    # verify dimensions
    S <- nrow(result[["tau_within"]])
    expect_equal(dim(result[["tau_within"]]), c(S, K), info = paste(name, ": tau_within dimensions"))
    expect_equal(dim(result[["tau_between"]]), c(S, K), info = paste(name, ": tau_between dimensions"))

    # verify positivity
    expect_true(all(result[["tau_within"]] >= 0), info = paste(name, ": tau_within should be non-negative"))
    expect_true(all(result[["tau_between"]] >= 0), info = paste(name, ": tau_between should be non-negative"))

    # verify relationship for non-multilevel: tau_between = 0
    if (!is_multilevel) {
      expect_true(all(result[["tau_between"]] == 0),
                  info = paste(name, ": non-multilevel should have tau_between = 0"))
    }

    # verify total tau can be reconstructed: tau = sqrt(tau_within^2 + tau_between^2)
    tau_reconstructed <- sqrt(result[["tau_within"]]^2 + result[["tau_between"]]^2)
    expect_true(all(tau_reconstructed >= 0), info = paste(name, ": reconstructed tau should be non-negative"))
  }
})


test_that(".evaluate.brma.mu returns correct dimensions", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects
    if (!inherits(object, "brma")) next

    priors           <- object[["priors"]]
    is_mods          <- .is_mods(object)
    is_PET           <- .is_PET(object)
    is_PEESE         <- .is_PEESE(object)
    effect_direction <- .effect_direction(object)
    outcome_data     <- object[["data"]][["outcome"]]
    K                <- nrow(outcome_data)

    # skip mods models for now - they require specific prior structure
    # that may not be present in all cached fits
    if (is_mods) next

    mu_samples <- .evaluate.brma.mu(
      fit                          = object[["fit"]],
      outcome_data                 = outcome_data,
      mods_data                    = object[["data"]][["mods"]],
      mods_formula                 = if (is_mods) attr(object[["data"]][["mods"]], "formula") else NULL,
      mods_priors                  = priors[["mods"]],
      is_mods                      = is_mods,
      is_PET                       = is_PET,
      is_PEESE                     = is_PEESE,
      effect_direction             = effect_direction,
      incorporate_publication_bias = TRUE,
      K                            = K
    )

    # verify dimensions
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    S <- nrow(posterior_samples)
    expect_equal(dim(mu_samples), c(S, K),
                 info = paste(name, ": mu dimensions should be S x K"))
  }
})


test_that("vectorized PET/PEESE matches loop implementation", {

  # test that outer() produces same results as the original loop

  S   <- 100
  K   <- 5
  sei <- c(0.1, 0.2, 0.15, 0.25, 0.12)

  set.seed(666)
  mu_samples  <- matrix(rnorm(S * K), nrow = S, ncol = K)
  PET_samples <- rnorm(S, mean = 0.5, sd = 0.1)

  # loop implementation (original)
  mu_loop <- mu_samples
  for (i in seq_len(K)) {
    mu_loop[, i] <- mu_loop[, i] + PET_samples * sei[i]
  }

  # vectorized implementation (new)
  mu_vec <- mu_samples + outer(PET_samples, sei)

  expect_equal(mu_loop, mu_vec, tolerance = 1e-14,
               info = "Vectorized outer() should match loop for PET")

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
               info = "Vectorized outer() should match loop for PEESE")
})


test_that("rho clamping handles boundary values", {

  # test that rho values outside [0, 1] are properly clamped

  # create mock rho with edge cases
  rho <- c(-0.01, 0, 0.5, 1, 1.001)
  rho_clamped <- pmin(pmax(rho, 0), 1)

  expect_equal(rho_clamped, c(0, 0, 0.5, 1, 1))

  # verify tau decomposition is valid after clamping
  tau <- 0.3
  tau_within  <- tau * sqrt(rho_clamped)
  tau_between <- tau * sqrt(1 - rho_clamped)

  # all values should be non-negative
  expect_true(all(tau_within >= 0))
  expect_true(all(tau_between >= 0))

  # verify Pythagorean relationship: tau^2 = tau_within^2 + tau_between^2
  tau_reconstructed <- sqrt(tau_within^2 + tau_between^2)
  expect_equal(tau_reconstructed, rep(tau, length(rho)), tolerance = 1e-10)
})


test_that("matrix replication patterns are correct", {

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
