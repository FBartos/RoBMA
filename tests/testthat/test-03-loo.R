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

# load cached fits
fits <- list()
info <- list()
for (name in list_fits()) {
  fits[[name]] <- load_fit(name)
  info[[name]] <- load_info(name)
}


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
  expect_equal(mean(apply(ll_brma2, 1, sum)), -10.05, tolerance = 0.01) # metafor: 'log Lik.' -8.106874 (df=4)

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

  # the AIC - LOO comparison is not exact but there should be some mapping

  AIC_metafor  <- AIC(fit_metafor)
  AIC_metafor2 <- AIC(fit_metafor2)

  loo_brma  <- loo(fit_brma)
  loo_brma2 <- suppressWarnings(loo(fit_brma2))

  waic_brma  <- waic(fit_brma)
  waic_brma2 <- suppressWarnings(waic(fit_brma2))

  expect_equal(loo_brma$estimates["looic", "Estimate"],  28.71, tolerance = 0.01) # metafor: 28.40474
  expect_equal(loo_brma2$estimates["looic", "Estimate"], 24.93, tolerance = 0.01) # metafor: 24.21375

  expect_equal(waic_brma$estimates["waic", "Estimate"],  28.62, tolerance = 0.01) # metafor: 28.40474
  expect_equal(waic_brma2$estimates["waic", "Estimate"], 23.93, tolerance = 0.01) # metafor: 24.21375
})

# ---------------------------------------------------------------------------- #
# loo_compare simple function test
# ---------------------------------------------------------------------------- #

test_that("loo_compare compares two brma models", {

  # get two brma fits
  fit_brma    <- fits[["bcg_meta-analysis"]]
  fit_brma2   <- fits[["bcg_meta-regression"]]

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
  fit_brma    <- fits[["bcg_meta-analysis"]]
  fit_brma2   <- fits[["bcg_meta-regression"]]

  # get two loo brma fits
  loo_brma    <- loo(fits[["bcg_meta-analysis"]])
  loo_brma2   <- suppressWarnings(loo(fits[["bcg_meta-regression"]]))

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

  mu_samples   <- matrix(0.15, nrow = S, ncol = K)
  tau_within   <- matrix(0.05, nrow = S, ncol = K)

  log_lik <- .outcome_pdf.norm(yi, mu_samples, tau_within, sei)

  expect_equal(dim(log_lik), c(S, K))

  # compute expected log-likelihood manually for first observation
  total_sd <- sqrt(0.05^2 + 0.1^2)
  expected_ll <- dnorm(0.1, mean = 0.15, sd = total_sd, log = TRUE)
  expect_equal(log_lik[1, 1], expected_ll, tolerance = 1e-10)
})

test_that(".outcome_pdf.binom computes correct log-likelihood", {

  # manual test with known values
  ai  <- c(10, 20)
  ci  <- c(15, 25)
  n1i <- c(100, 100)
  n2i <- c(100, 100)
  K   <- length(ai)
  S   <- 5

  mu_samples     <- matrix(0.0, nrow = S, ncol = K)  # log-OR = 0
  logit_baserate <- matrix(qlogis(0.15), nrow = S, ncol = K)  # baseline ~15%

  log_lik <- .outcome_pdf.binom(ai, ci, n1i, n2i, mu_samples, logit_baserate)

  expect_equal(dim(log_lik), c(S, K))

  # with log-OR = 0, both groups should have same probability
  p <- plogis(qlogis(0.15))
  expected_ll <- dbinom(10, 100, p, log = TRUE) + dbinom(15, 100, p, log = TRUE)
  expect_equal(log_lik[1, 1], expected_ll, tolerance = 1e-10)
})

test_that(".outcome_pdf.pois computes correct log-likelihood", {

  # manual test with known values
  x1i <- c(10, 20)
  x2i <- c(15, 25)
  t1i <- c(100, 100)
  t2i <- c(100, 100)
  K   <- length(x1i)
  S   <- 5

  mu_samples <- matrix(0.0, nrow = S, ncol = K)   # log-IRR = 0
  log_phi    <- matrix(log(0.1), nrow = S, ncol = K)  # baseline rate = 0.1

  log_lik <- .outcome_pdf.pois(x1i, x2i, t1i, t2i, mu_samples, log_phi)

  expect_equal(dim(log_lik), c(S, K))

  # with log-IRR = 0, both groups should have same rate
  lambda <- 0.1 * 100  # rate * time
  expected_ll <- dpois(10, lambda, log = TRUE) + dpois(15, lambda, log = TRUE)
  expect_equal(log_lik[1, 1], expected_ll, tolerance = 1e-10)
})
