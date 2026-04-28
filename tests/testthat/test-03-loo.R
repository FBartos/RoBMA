# ============================================================================ #
#
# Test LOO-PSIS functionality for brma objects
#
# ============================================================================
context("loo methods for brma objects")

# The LOO-PSIS functionality is necessary for the residuals and funnel plot
# functionality. Due to the computational costs (and the possibility to test)
# against other metafor output) it is primary tested therein.

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_no_fits()
skip_if_not_installed("loo")

# list cached fits lazily
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)


# ---------------------------------------------------------------------------- #
# logLik simple function test
# ---------------------------------------------------------------------------- #

test_that("logLik computes log-likelihood matrix with correct dimensions", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  name2        <- "bcg_meta-regression"
  fit_metafor2 <- info[[name2]][["metafor"]]
  fit_brma2    <- fits[[name2]]

  # the loglik comparison is meaningless, however, the average logLik should be
  # at least in the same ballpark as frequentist logLik + test for consistency

  ll_metafor  <- logLik(fit_metafor)
  ll_metafor2 <- logLik(fit_metafor2)

  ll_brma  <- logLik(fit_brma)
  ll_brma2 <- logLik(fit_brma2)

  expect_equal(mean(apply(ll_brma,  1, sum)), -13.60, tolerance = 0.01) # metafor: 'log Lik.' -12.20237 (df=2)
  expect_equal(mean(apply(ll_brma2, 1, sum)), -10.58, tolerance = 0.01) # metafor: 'log Lik.' -8.106874 (df=4)

  expect_equal(ncol(ll_brma),  nrow(fit_brma$data$outcome))
  expect_equal(ncol(ll_brma2), nrow(fit_brma2$data$outcome))

})


# ---------------------------------------------------------------------------- #
# loo/WAIC simple function test
# ---------------------------------------------------------------------------- #

test_that("loo/WAIC computes roughly matches AIC", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  name2        <- "bcg_meta-regression"
  fit_metafor2 <- info[[name2]][["metafor"]]
  fit_brma2    <- fits[[name2]]

  ### the AIC - LOO/WAIC comparison is not exact but there should be some mapping

  AIC_metafor  <- AIC(fit_metafor)
  AIC_metafor2 <- AIC(fit_metafor2)

  ## loo
  loo_brma  <- loo(fit_brma)
  loo_brma2 <- loo(fit_brma2)

  # simulate loo not being computed to test the error
  fit_brma[["loo"]] <- NULL
  expect_error(loo(fit_brma), "LOO has not been computed")

  ## waic (not precompted for test objects)
  expect_error(waic(fit_brma), "WAIC has not been computed")
  fit_brma   <- add_waic(fit_brma)
  fit_brma2  <- suppressWarnings(add_waic(fit_brma2))
  waic_brma  <- waic(fit_brma)
  waic_brma2 <- waic(fit_brma2)

  expect_equal(loo_brma$estimates["looic", "Estimate"],  28.71, tolerance = 0.01) # metafor: 28.40474
  expect_equal(loo_brma2$estimates["looic", "Estimate"], 26.00, tolerance = 0.01) # metafor: 24.21375

  expect_equal(waic_brma$estimates["waic", "Estimate"],  28.62, tolerance = 0.01) # metafor: 28.40474
  expect_equal(waic_brma2$estimates["waic", "Estimate"], 25.12, tolerance = 0.01) # metafor: 24.21375
})

# ---------------------------------------------------------------------------- #
# loo_compare simple function test
# ---------------------------------------------------------------------------- #

test_that("loo_compare compares two brma models", {

  # get two brma fits
  fit_brma  <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  # compare
  out <- suppressWarnings(loo_compare(fit_brma, fit_brma2))

  # check structure
  expect_true(is.matrix(out))
  expect_equal(nrow(out), 2)
  expect_true("elpd_diff" %in% colnames(out))
  expect_true("se_diff" %in% colnames(out))
})

test_that("loo_compare works with loo objects", {

  # get two brma fits
  fit_brma  <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  # get two loo brma fits
  loo_brma  <- loo(fits[["bcg_meta-analysis"]])
  loo_brma2 <- suppressWarnings(loo(fits[["bcg_meta-regression"]]))

  # compare
  out <- loo_compare(loo_brma, loo_brma2)

  # check structure
  expect_true(is.matrix(out))
  expect_equal(nrow(out), 2)
  expect_true("elpd_diff" %in% colnames(out))
  expect_true("se_diff" %in% colnames(out))
})

test_that("loo_compare errors with < 2 models", {

  # get one brma fit
  fit_brma <- fits[["bcg_meta-analysis"]]

  expect_error(loo_compare(fit_brma, "At least two models"))
})

test_that("logLik, LOO, weights, diagnostics, and WAIC work for product-space fits", {

  product_names <- c(
    "dat.lehmann2018_BMA.norm",
    "bcg_BMA.glmm",
    "dat.lehmann2018_RoBMA"
  )

  for (name in product_names) {

    fit_brma <- fits[[name]]

    log_lik <- logLik(fit_brma)
    expect_s3_class(log_lik, "logLik.brma")
    expect_true(is.matrix(log_lik), info = name)
    expect_equal(ncol(log_lik), nobs(fit_brma), info = name)
    expect_true(all(is.finite(log_lik)), info = name)

    loo_result <- loo(fit_brma)
    expect_s3_class(loo_result, "loo")

    weights <- loo_weights(fit_brma)
    expect_true(is.matrix(weights), info = name)
    expect_equal(dim(weights), dim(log_lik), info = name)
    expect_equal(colSums(weights), rep(1, ncol(weights)), tolerance = 1e-10,
                 info = name)

    expect_no_error(suppressWarnings(check_loo(fit_brma)))

    fit_waic <- fit_brma
    fit_waic[["waic"]] <- NULL
    fit_waic <- suppressWarnings(add_waic(fit_waic))
    waic_result <- waic(fit_waic)
    expect_s3_class(waic_result, "waic")
  }
})

test_that("loo_compare compares BMA and RoBMA product-space fits on the same data", {

  product_names <- c("dat.lehmann2018_BMA.norm", "dat.lehmann2018_RoBMA")

  out <- suppressWarnings(loo_compare(
    fits[["dat.lehmann2018_BMA.norm"]],
    fits[["dat.lehmann2018_RoBMA"]]
  ))

  expect_true(is.matrix(out))
  expect_equal(nrow(out), 2)
  expect_true("elpd_diff" %in% colnames(out))
  expect_true("se_diff" %in% colnames(out))
})


# ---------------------------------------------------------------------------- #
# Outcome type specific tests
# ---------------------------------------------------------------------------- #

test_that(".outcome_pdf.norm computes correct log-likelihood", {

  # manual test with known values
  set.seed(123)

  yi  <- c(0.1, 0.2, 0.3)
  sei <- c(0.1, 0.1, 0.1)
  K   <- length(yi)
  S   <- 10

  mu_samples <- matrix(0.15, nrow = S, ncol = K)
  tau_within <- matrix(0.05, nrow = S, ncol = K)

  log_lik <- .outcome_pdf.norm(yi, mu_samples, tau_within, sei)

  expect_equal(dim(log_lik), c(S, K))

  # compute expected log-likelihood manually for first observation
  total_sd <- sqrt(0.05^2 + 0.1^2)
  expected_ll <- dnorm(0.1, mean = 0.15, sd = total_sd, log = TRUE)
  expect_equal(log_lik[1, 1], expected_ll, tolerance = 1e-10)
})

test_that(".outcome_cdf.norm keeps negative-direction tail precision", {

  yi         <- 10
  sei        <- 1
  mu_samples <- matrix(0, nrow = 1, ncol = 1)
  tau_within <- matrix(0, nrow = 1, ncol = 1)

  cdf_vals <- .outcome_cdf.norm(
    yi         = yi,
    mu_samples = mu_samples,
    tau_within = tau_within,
    sei        = sei,
    lower.tail = FALSE
  )

  expect_equal(cdf_vals[1, 1], stats::pnorm(yi, lower.tail = FALSE))
  expect_gt(cdf_vals[1, 1], 0)
})

test_that(".outcome_cdf.wnorm forwards lower.tail to normal fast path", {

  yi         <- 10
  sei        <- 1
  mu_samples <- matrix(0, nrow = 1, ncol = 1)
  tau_within <- matrix(0, nrow = 1, ncol = 1)
  omega      <- matrix(c(1, 1), nrow = 1)
  crit_yi    <- matrix(0, nrow = 1, ncol = 1)

  cdf_vals <- .outcome_cdf.wnorm(
    yi         = yi,
    mu_samples = mu_samples,
    tau_within = tau_within,
    sei        = sei,
    omega      = omega,
    crit_yi    = crit_yi,
    lower.tail = FALSE,
    use_normal = TRUE
  )

  expect_equal(cdf_vals[1, 1], stats::pnorm(yi, lower.tail = FALSE))
  expect_gt(cdf_vals[1, 1], 0)
})

test_that(".outcome_pdf.binom computes correct log-likelihood", {

  # manual test with known values
  ai  <- c(10, 20)
  ci  <- c(15, 25)
  n1i <- c(100, 100)
  n2i <- c(100, 100)
  K   <- length(ai)
  S   <- 5

  mu_samples <- matrix(0.0, nrow = S, ncol = K) # log-OR = 0
  tau_within <- matrix(0.0, nrow = S, ncol = K)
  prior_pi   <- BayesTools::prior("beta", list(1, 1))

  log_lik <- .outcome_pdf.binom(ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi)

  expect_equal(dim(log_lik), c(S, K))

  # with log-OR = 0, both groups should have same probability
  # value updated to match the marginal likelihood calculation
  expect_equal(log_lik[1, 1], -7.64, tolerance = 0.01)
})

test_that(".outcome_pdf.binom handles boundary cell studies", {

  prior_pi   <- BayesTools::prior("beta", list(1, 1))
  mu_samples <- matrix(0, nrow = 2, ncol = 4)
  tau_within <- matrix(0, nrow = 2, ncol = 4)

  log_lik <- .outcome_pdf.binom(
    ai         = c(10, 0, 10, 0),
    ci         = c(10, 0, 0, 10),
    n1i        = c(10, 10, 10, 10),
    n2i        = c(10, 10, 10, 10),
    mu_samples = mu_samples,
    tau_within = tau_within,
    prior_pi   = prior_pi
  )

  expect_equal(dim(log_lik), c(2, 4))
  expect_true(all(is.finite(log_lik)))
})

test_that(".outcome_pdf.binom matches R reference", {

  prior_pi <- BayesTools::prior("beta", list(1.5, 2.5))

  ai  <- c(3L, 0L, 8L, 12L)
  ci  <- c(2L, 4L, 0L, 10L)
  n1i <- c(20L, 18L, 15L, 30L)
  n2i <- c(22L, 17L, 16L, 31L)

  mu_samples <- matrix(
    c(-0.4, 0.1, 0.7, 0.2,
      -0.2, 0.3, 0.5, 0.0,
       0.1, 0.5, 0.3, -0.2),
    nrow = 3,
    ncol = 4
  )
  tau_within <- matrix(
    c(0.0, 0.1, 0.4, 0.2,
      0.2, 0.3, 0.1, 0.5,
      0.4, 0.2, 0.0, 0.3),
    nrow = 3,
    ncol = 4
  )

  expect_equal(
    .outcome_pdf.binom(ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi, n_theta = 7, n_pi = 9),
    .outcome_pdf.binom_r(ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi, n_theta = 7, n_pi = 9),
    tolerance = 1e-10
  )
})

test_that(".outcome_pdf.pois computes correct log-likelihood", {

  # manual test with known values
  x1i <- c(10, 20)
  x2i <- c(15, 25)
  t1i <- c(100, 100)
  t2i <- c(100, 100)
  K   <- length(x1i)
  S   <- 5

  mu_samples <- matrix(0.0, nrow = S, ncol = K) # log-IRR = 0
  tau_within <- matrix(0.0, nrow = S, ncol = K)
  prior_phi  <- BayesTools::prior("normal", list(0, 1))

  log_lik <- .outcome_pdf.pois(x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi)

  expect_equal(dim(log_lik), c(S, K))

  # with log-IRR = 0, both groups should have same rate
  # value updated to match the marginal likelihood calculation
  expect_equal(log_lik[1, 1], -8.60, tolerance = 0.01)
})

test_that(".outcome_pdf.pois matches R reference", {

  prior_phi <- BayesTools::prior("normal", list(-1, 1.5))

  x1i <- c(0L, 3L, 10L, 12L)
  x2i <- c(1L, 0L, 8L, 9L)
  t1i <- c(12, 30, 45, 50)
  t2i <- c(10, 28, 43, 48)

  mu_samples <- matrix(
    c(-0.4, 0.1, 0.7, 0.2,
      -0.2, 0.3, 0.5, 0.0,
       0.1, 0.5, 0.3, -0.2),
    nrow = 3,
    ncol = 4
  )
  tau_within <- matrix(
    c(0.0, 0.1, 0.4, 0.2,
      0.2, 0.3, 0.1, 0.5,
      0.4, 0.2, 0.0, 0.3),
    nrow = 3,
    ncol = 4
  )

  expect_equal(
    .outcome_pdf.pois(x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi, n_theta = 7, n_phi = 9),
    .outcome_pdf.pois_r(x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi, n_theta = 7, n_phi = 9),
    tolerance = 1e-10
  )
})

test_that("native GLMM cluster likelihood matches R composition", {

  skip_if_not(.has_native_glmm_cluster(), "Native GLMM cluster kernels unavailable.")

  set.seed(2024)
  S <- 5
  K <- 5

  setup <- list(
    mu          = matrix(rnorm(S * K, 0, 0.25), nrow = S, ncol = K),
    tau_within  = matrix(runif(S * K, 0.05, 0.25), nrow = S, ncol = K),
    tau_between = matrix(runif(S * K, 0.02, 0.18), nrow = S, ncol = K),
    cluster     = list(a = c(1L, 3L), b = c(2L, 4L, 5L)),
    weights     = c(1, 0.5, 1.25, 2, 0.75)
  )

  bin_data <- list(outcome = data.frame(
    ai  = c(3L, 0L, 8L, 12L, 2L),
    ci  = c(2L, 4L, 0L, 10L, 1L),
    n1i = c(20L, 18L, 15L, 30L, 22L),
    n2i = c(22L, 17L, 16L, 31L, 21L)
  ))
  bin_priors <- list(outcome = list(
    pi = BayesTools::prior("beta", list(1.5, 2.5))
  ))

  expect_equal(
    .log_lik_cluster_glmm_native(
      setup        = setup,
      data         = bin_data,
      priors       = bin_priors,
      outcome_type = "bin",
      n_theta      = 5,
      n_gamma      = 5,
      n_pi         = 7
    ),
    .log_lik_cluster_glmm_r(
      setup        = setup,
      data         = bin_data,
      priors       = bin_priors,
      outcome_type = "bin",
      n_theta      = 5,
      n_gamma      = 5,
      n_pi         = 7
    ),
    tolerance = 1e-10
  )

  pois_data <- list(outcome = data.frame(
    x1i = c(0L, 3L, 10L, 12L, 1L),
    x2i = c(1L, 0L, 8L, 9L, 2L),
    t1i = c(12, 30, 45, 50, 25),
    t2i = c(10, 28, 43, 48, 24)
  ))
  pois_priors <- list(outcome = list(
    phi = BayesTools::prior("normal", list(-1, 1.5))
  ))

  expect_equal(
    .log_lik_cluster_glmm_native(
      setup        = setup,
      data         = pois_data,
      priors       = pois_priors,
      outcome_type = "pois",
      n_theta      = 5,
      n_gamma      = 5,
      n_phi        = 7
    ),
    .log_lik_cluster_glmm_r(
      setup        = setup,
      data         = pois_data,
      priors       = pois_priors,
      outcome_type = "pois",
      n_theta      = 5,
      n_gamma      = 5,
      n_phi        = 7
    ),
    tolerance = 1e-10
  )
})


# ---------------------------------------------------------------------------- #
# loo_weights and check_loo S3 tests
# ---------------------------------------------------------------------------- #

test_that("loo_weights and check_loo work correctly", {

  fit_brma <- fits[["bcg_meta-analysis"]]

  # check loo_weights
  weights <- loo_weights(fit_brma)
  expect_true(is.matrix(weights))
  expect_equal(dim(weights), dim(logLik(fit_brma)))
  expect_equal(colSums(weights), rep(1, ncol(weights)), tolerance = 1e-10)

  # check check_loo (should not throw anything for this clean fit)
  expect_silent(check_loo(fit_brma))

  # simulate bad k
  fit_bad <- fit_brma
  fit_bad[["loo"]][["estimate"]][["diagnostics"]][["pareto_k"]][1] <- 0.8
  expect_warning(check_loo(fit_bad), "Some Pareto k values are high")
})

