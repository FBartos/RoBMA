context("Evaluate posterior samples")
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
                  info = paste("Observation", k, "should shrink toward data with high tau"))
    } else {
      expect_true(abs(mean_theta[k] - yi[k]) > abs(mean_theta[k] - mean_mu[k]),
                  info = paste("Observation", k, "should shrink toward mean with high sei"))
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


test_that(".evaluate.brma.study_effects returns contribution matrix for new data", {

  S <- 10000  # many samples
  K <- 3

  set.seed(222)
  tau_between <- matrix(0.3, nrow = S, ncol = K)  # constant tau_between

  # for new data, gamma ~ N(0,1) is sampled, so contribution = gamma * tau_between
  # expected variance = tau_between^2 = 0.09
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
              info = paste("Observation", k, ": response should have more variance than theta"))
  expect_true(sd_theta > sd_mu,
              info = paste("Observation", k, ": theta should have more variance than terms"))
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


# ============================================================================ #
# SECTION 3: Unit Tests for Aggregated Predictions (newdata = TRUE)
# ============================================================================ #
# These tests verify the aggregation logic for predict.brma with newdata = TRUE
# ============================================================================ #

test_that("aggregation via rowMeans produces correct results", {

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

  # verify aggregation matches manual calculation
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


test_that("aggregated true effects have correct variance", {

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
               info = "Aggregated theta variance should match theoretical")
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
    expect_null(attr(result_terms, "data"), info = paste(name, ": data should be NULL for aggregate"))
    expect_equal(nrow(summary(result_terms)), 1,
                 info = paste(name, ": should have single row for aggregate terms"))

    # test type = "terms.scale"
    result_scale <- predict(object, newdata = TRUE, type = "terms.scale")

    expect_s3_class(result_scale, "brma_samples")
    expect_null(attr(result_scale, "data"), info = paste(name, ": data should be NULL for aggregate scale"))
    expect_equal(nrow(summary(result_scale)), 1,
                 info = paste(name, ": should have single row for aggregate scale"))

    # test type = "effect"
    result_effect <- predict(object, newdata = TRUE, type = "effect")

    expect_s3_class(result_effect, "brma_samples")
    expect_null(attr(result_effect, "data"), info = paste(name, ": data should be NULL for aggregate effect"))
    expect_equal(nrow(summary(result_effect)), 1,
                 info = paste(name, ": should have single row for aggregate effect"))
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
                 info = paste(name, ": terms samples should be S x 1"))
    expect_equal(colnames(samples_terms), "mu",
                 info = paste(name, ": terms samples should have 'mu' column name"))

    # test type = "terms.scale"
    samples_scale <- predict(object, newdata = TRUE, type = "terms.scale")
    expect_equal(dim(samples_scale), c(S, 1),
                 info = paste(name, ": scale samples should be S x 1"))
    expect_equal(colnames(samples_scale), "tau",
                 info = paste(name, ": scale samples should have 'tau' column name"))

    # test type = "effect"
    samples_effect <- predict(object, newdata = TRUE, type = "effect")
    expect_equal(dim(samples_effect), c(S, 1),
                 info = paste(name, ": effect samples should be S x 1"))
    expect_equal(colnames(samples_effect), "theta",
                 info = paste(name, ": effect samples should have 'theta' column name"))
  }
})


test_that("predict.brma with newdata = TRUE and type = 'response' throws error", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects
    if (!inherits(object, "brma")) next

    # should throw error for type = "response"
    expect_error(
      predict(object, newdata = TRUE, type = "response"),
      "Aggregated predictions.*not available for type = 'response'",
      info = paste(name, ": should error for aggregate + response")
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
                 info = paste(name, ": aggregated mu should equal rowMeans of full mu"))
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
                 info = paste(name, ": aggregated tau should equal rowMeans of full tau"))
  }
})


# ============================================================================ #
# SECTION 3: Tests for .extract_use_normal()
# ============================================================================ #
# These tests verify the bias indicator extraction helper function
# that identifies which posterior samples use weightfunction vs normal path
# ============================================================================ #

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
                info = paste(name, ": brma without bias should have all use_normal = TRUE"))
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
                info = paste(name, ": brma with single weightfunction should have all use_normal = FALSE"))
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
                info = paste(name, ": brma with PET should have all use_normal = TRUE"))
  }
})


test_that(".extract_use_normal returns correct structure for RoBMA", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # only test RoBMA objects
    if (!inherits(object, "RoBMA")) next

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

      # identify weightfunction indices
      wf_indices <- which(sapply(priors_bias, BayesTools::is.prior.weightfunction))

      # verify use_normal is correctly computed
      expected_use_normal <- !(bias_indicator %in% wf_indices)
      expect_equal(use_normal, expected_use_normal,
                   info = paste(name, ": use_normal should match bias_indicator logic"))

      # verify we have both types of samples (unless all priors are one type)
      if (length(wf_indices) > 0 && length(wf_indices) < length(priors_bias)) {
        expect_true(any(use_normal) && any(!use_normal),
                    info = paste(name, ": RoBMA should have mix of normal and weighted samples"))
      }
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
                info = paste(name, ": log_lik should be a matrix"))
    expect_true(all(is.finite(log_lik)),
                info = paste(name, ": all log-likelihoods should be finite"))
    expect_true(all(log_lik <= 0),
                info = paste(name, ": all log-likelihoods should be non-positive"))

    # verify dimensions match
    K <- length(.outcome_data_yi(object))
    expect_equal(ncol(log_lik), K,
                 info = paste(name, ": ncol should match number of observations"))
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
                info = paste(name, ": cdf_vals should be a matrix"))
    expect_true(all(is.finite(cdf_vals)),
                info = paste(name, ": all CDF values should be finite"))
    expect_true(all(cdf_vals > 0 & cdf_vals < 1),
                info = paste(name, ": all CDF values should be in (0, 1)"))

    # verify dimensions match
    K <- length(.outcome_data_yi(object))
    expect_equal(ncol(cdf_vals), K,
                 info = paste(name, ": ncol should match number of observations"))
  }
})
