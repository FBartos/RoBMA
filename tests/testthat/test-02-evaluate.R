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

test_that("known-V BLUP block solver fails on invalid covariance blocks", {

  expect_error(
    .solve_diagonal_rank_one_block(
      diagonal = c(1, -2),
      rank_one = c(0, 0),
      residual = c(1, 1)
    ),
    "not positive definite"
  )
  expect_error(
    .solve_diagonal_rank_one_block(
      diagonal = c(1, NA_real_),
      rank_one = c(0, 0),
      residual = c(1, 1)
    ),
    "non-finite"
  )
})

test_that("known-V BLUP uses full covariance blocks", {

  mu_samples <- matrix(c(0.2, -0.1), nrow = 1)
  tau_within <- matrix(c(0.5, 0.25), nrow = 1)
  yi         <- c(1.0, -0.7)
  V          <- matrix(c(0.04, 0.02, 0.02, 0.09), nrow = 2)
  known_V    <- .known_v_newdata_prepare(V, k = length(yi))

  theta <- .evaluate.brma.known_v_blup.norm(
    mu_samples = mu_samples,
    tau_within = tau_within,
    yi         = yi,
    known_V    = known_V
  )

  T_block  <- diag(tau_within[1, ]^2)
  expected <- mu_samples[1, ] +
    as.vector(T_block %*% solve(T_block + V) %*% (yi - mu_samples[1, ]))
  diagonal <- mu_samples[1, ] +
    tau_within[1, ]^2 / (tau_within[1, ]^2 + diag(V)) * (yi - mu_samples[1, ])

  expect_equal(theta[1, ], expected, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(expected, diagonal, tolerance = 1e-8)))
})

test_that("known-V BLUP scalar-tau blocks match inverse oracle", {

  mu_samples <- matrix(
    c(.10, .20, .30,
      .25, .35, .45,
      .40, .50, .60),
    nrow = 3,
    byrow = TRUE
  )
  tau_within <- matrix(
    c(.20, .20, .20,
      .35, .35, .35,
      .50, .50, .50),
    nrow = 3,
    byrow = TRUE
  )
  yi <- c(.50, -.10, .80)
  V <- matrix(
    c(.09, .03, .01,
      .03, .16, .02,
      .01, .02, .25),
    nrow = 3
  )
  bias_offset <- matrix(
    c(.02, .01, -.03,
      .00, .04,  .01,
      .03, .02,  .00),
    nrow = 3,
    byrow = TRUE
  )

  theta <- .evaluate.brma.known_v_blup.norm(
    mu_samples  = mu_samples,
    tau_within  = tau_within,
    yi          = yi,
    known_V     = .known_v_newdata_prepare(V, k = length(yi)),
    bias_offset = bias_offset
  )
  expected <- mu_samples
  for (s in seq_len(nrow(mu_samples))) {
    T_block <- diag(tau_within[s, ]^2)
    expected[s, ] <- mu_samples[s, ] +
      as.vector(T_block %*% solve(T_block + V) %*%
                  (yi - bias_offset[s, ] - mu_samples[s, ]))
  }

  expect_equal(theta, expected, tolerance = 1e-12)
})


test_that("known-V BLUP row-varying tau blocks match inverse oracle", {

  mu_samples <- matrix(
    c(.10, .20, .30,
      .25, .35, .45),
    nrow = 2,
    byrow = TRUE
  )
  tau_within <- matrix(
    c(.20, .35, .50,
      .45, .25, .30),
    nrow = 2,
    byrow = TRUE
  )
  yi <- c(.50, -.10, .80)
  V <- matrix(
    c(.09, .03, .01,
      .03, .16, .02,
      .01, .02, .25),
    nrow = 3
  )

  theta <- .evaluate.brma.known_v_blup.norm(
    mu_samples = mu_samples,
    tau_within = tau_within,
    yi         = yi,
    known_V    = .known_v_newdata_prepare(V, k = length(yi))
  )
  expected <- mu_samples
  for (s in seq_len(nrow(mu_samples))) {
    T_block <- diag(tau_within[s, ]^2)
    expected[s, ] <- mu_samples[s, ] +
      as.vector(T_block %*% solve(T_block + V) %*%
                  (yi - mu_samples[s, ]))
  }

  expect_equal(theta, expected, tolerance = 1e-12)
})


test_that("known-V BLUP singleton blocks use scalar shrinkage", {

  mu_samples <- matrix(c(.10, .30), ncol = 1)
  tau_within <- matrix(c(.20, .50), ncol = 1)
  yi         <- .70
  V          <- matrix(.09, nrow = 1)
  bias_offset <- matrix(c(.05, -.10), ncol = 1)

  theta <- .evaluate.brma.known_v_blup.norm(
    mu_samples  = mu_samples,
    tau_within  = tau_within,
    yi          = yi,
    known_V     = .known_v_newdata_prepare(V, k = length(yi)),
    bias_offset = bias_offset
  )
  expected <- mu_samples[, 1L] +
    tau_within[, 1L]^2 / (tau_within[, 1L]^2 + V[1L, 1L]) *
      (yi - bias_offset[, 1L] - mu_samples[, 1L])

  expect_equal(theta[, 1L], expected, tolerance = 1e-12)
})


test_that("known-V BLUP helper rejects unsolvable covariance blocks", {

  invalid_known_V <- .new_known_v(list(
    version       = 2L,
    storage       = "dense",
    K             = 1L,
    diagonal      = 0,
    V             = matrix(-.5, nrow = 1),
    blocks        = NULL,
    block_indices = list(1L)
  ))

  expect_error(
    .evaluate.brma.known_v_blup.norm(
      mu_samples = matrix(0, nrow = 1, ncol = 1),
      tau_within = matrix(.1, nrow = 1, ncol = 1),
      yi         = 0,
      known_V    = invalid_known_V
    ),
    "positive definite"
  )
})

test_that("LOO chain IDs require retained MCMC draws", {

  expect_error(
    .loo_chain_id(list(sample = 10), n_samples = 10),
    "fitted MCMC draws are missing"
  )
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

test_that(".outcome_rng.norm_known_v uses Cholesky orientation for full covariance", {

  S <- 60000L
  K <- 3L
  V <- matrix(
    c(
      .040, .018, .012,
      .018, .090, .026,
      .012, .026, .160
    ),
    nrow  = K,
    byrow = TRUE
  )
  tau <- .15
  known_V         <- .known_v_canonicalize(V)
  sampling_factor <- .known_v_sampling_factor(V)

  expect_equal(
    t(sampling_factor) %*% sampling_factor,
    V,
    tolerance = 1e-14
  )

  set.seed(20260615)
  response <- .outcome_rng.norm_known_v(
    mu_samples = matrix(0, nrow = S, ncol = K),
    tau_within = matrix(tau, nrow = S, ncol = K),
    known_V    = known_V
  )
  expected_cov <- V + diag(tau^2, nrow = K)

  expect_equal(dim(response), c(S, K))
  expect_equal(colMeans(response), rep(0, K), tolerance = .01)
  expect_equal(stats::cov(response), expected_cov, tolerance = .01)
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

.multilevel_reference_covariance <- function(tau_within, tau_between, vi,
                                             block_indices) {

  covariance <- diag(
    vi + tau_within^2,
    nrow = length(vi),
    ncol = length(vi)
  )
  for (indices in block_indices) {
    covariance[indices, indices] <- covariance[indices, indices] +
      tcrossprod(tau_between[indices])
  }

  return(covariance)
}


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
    covariance <- .multilevel_reference_covariance(
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
    covariance <- .multilevel_reference_covariance(
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

test_that("selection row routing validates posterior bias indicators", {

  bias <- BayesTools::prior_mixture(list(
    BayesTools::prior_none(),
    BayesTools::prior_weightfunction(
      side    = "one-sided",
      steps   = c(.025),
      weights = BayesTools::wf_fixed(c(1, .5))
    )
  ))
  object <- list(
    fit    = NULL,
    priors = list(outcome = list(bias = bias))
  )

  posterior_samples <- matrix(c(1, 2, 1, 2), ncol = 1)
  colnames(posterior_samples) <- "bias_indicator"

  expect_equal(
    .extract_bias_indicator(object, posterior_samples = posterior_samples),
    c(1L, 2L, 1L, 2L)
  )
  expect_equal(
    .extract_use_normal(object, posterior_samples = posterior_samples),
    c(TRUE, FALSE, TRUE, FALSE)
  )

  invalid <- posterior_samples
  invalid[, "bias_indicator"] <- c(1, 0, 1, 2)
  expect_error(
    .extract_bias_indicator(object, posterior_samples = invalid),
    "Invalid posterior model indicator range"
  )

  invalid[, "bias_indicator"] <- c(1, NA, 1, 2)
  expect_error(
    .extract_use_normal(object, posterior_samples = invalid),
    "Invalid posterior model indicator"
  )

  invalid[, "bias_indicator"] <- c(1, 3, 1, 2)
  expect_error(
    .selection_row_routing(
      priors               = object[["priors"]],
      posterior_samples    = invalid
    ),
    "Invalid posterior model indicator range"
  )
})

test_that("selected-normal RNG requires explicit row routing", {

  expect_error(
    .outcome_rng.selnorm(
      mu_samples        = matrix(0, nrow = 2, ncol = 1),
      tau_within        = matrix(0, nrow = 2, ncol = 1),
      sei               = 1,
      selection_context = list(kernel_mode = c(0L, 0L))
    ),
    "use_normal"
  )
})

test_that("outcome CDF values retain exact probability endpoints", {

  cdf_vals <- .outcome_cdf.norm(
    yi         = c(-100, 100),
    mu_samples = matrix(0, nrow = 1L, ncol = 2L),
    tau_within = matrix(0, nrow = 1L, ncol = 2L),
    sei        = c(1, 1)
  )

  expect_identical(cdf_vals, matrix(c(0, 1), nrow = 1L))

  interior_tail <- .outcome_cdf.norm(
    yi         = c(-38, 8.2),
    mu_samples = matrix(0, nrow = 1L, ncol = 2L),
    tau_within = matrix(0, nrow = 1L, ncol = 2L),
    sei        = c(1, 1)
  )
  expect_equal(as.numeric(interior_tail), stats::pnorm(c(-38, 8.2)))
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

    posterior_samples <- as.matrix(object[["fit"]][["mcmc"]])
    has_tau_samples <- !is.null(.extract_indexed_parameter_samples(
      posterior_samples = posterior_samples,
      parameter         = "tau",
      required          = FALSE
    ))

    evaluate_tau <- function(fixed_tau = NULL) {
      .evaluate.brma.tau(
        fit               = object[["fit"]],
        scale_data        = object[["data"]][["scale"]],
        scale_formula     = if (is_scale) attr(object[["data"]][["scale"]], "formula") else NULL,
        scale_priors      = priors[["scale"]],
        is_scale          = is_scale,
        is_multilevel     = is_multilevel,
        K                 = K,
        posterior_samples = posterior_samples,
        fixed_tau         = fixed_tau,
        fixed_rho         = .fixed_rho_prior_value(priors)
      )
    }

    fixed_tau <- .fixed_tau_prior_value(priors)
    if (!has_tau_samples && !is.null(fixed_tau)) {
      expect_error(
        evaluate_tau(),
        "Missing posterior tau columns",
        info = paste(name, ": missing tau requires explicit allowance")
      )
      result <- evaluate_tau(fixed_tau = fixed_tau)
      expect_equal(
        result[["tau_total"]],
        matrix(fixed_tau, nrow = nrow(posterior_samples), ncol = K),
        info = paste(name, ": fixed tau is reconstructed from the prior")
      )
    } else {
      result <- evaluate_tau()
    }

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

    # Direct moderator/random-design handling is exercised through predict().
    if (is_mods || .is_random(object)) next

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
      K                 = K,
      priors            = priors
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


test_that("Binomial baserate evaluation preserves exact endpoints", {

  posterior_samples <- cbind(mu = c(0.1, 0.2), "pi[1]" = c(0, 1))

  expect_identical(
    as.vector(.evaluate.brma.baserate(
      fit               = NULL,
      K                 = 1L,
      posterior_samples = posterior_samples
    )),
    c(-Inf, Inf)
  )
  expect_identical(
    as.vector(.evaluate.brma.baserate_newdata(
      prior_pi = BayesTools::prior("spike", parameters = list(location = 0)),
      S        = 2L,
      K        = 1L
    )),
    rep(-Inf, 2L)
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

# ============================================================================ #
# SECTION 4: Integration Tests for .log_lik.brma() and .cdf.brma() with use_normal
# ============================================================================ #
# These tests verify that PDF and CDF functions work correctly with the
# use_normal fast-path optimization
# ============================================================================ #

test_that(".log_lik.brma returns finite log-likelihoods for weightfunction models", {

  for (name in names(fits)) {

    object <- fits[[name]]

    # skip non-brma objects or objects without weightfunction priors
    if (!inherits(object, "brma")) next
    if (!.is_weightfunction(object)) next

    # compute PDF (this internally uses use_normal optimization)
    log_lik <- .log_lik.brma(object)

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
    expect_true(all(cdf_vals >= 0 & cdf_vals <= 1),
                info = paste(name, ": all CDF values are in [0, 1]"))

    # verify dimensions match
    K <- length(.outcome_data_yi(object))
    expect_equal(ncol(cdf_vals), K,
                 info = paste(name, ": ncol matches number of observations"))
  }
})
