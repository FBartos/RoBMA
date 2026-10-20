context("Direct arithmetic of the multilevel and true-effect evaluators")
skip_on_cran()

.direct_arithmetic_fixture <- function(scale = 1, S = 40L, K = 12L) {

  set.seed(20261020)
  list(
    mu = matrix(rnorm(S * K), S, K) * scale,
    tau_within = matrix(abs(rnorm(S * K, .3, .1)), S, K) * scale,
    tau_between = matrix(abs(rnorm(S * K, .2, .1)), S, K) * scale,
    bias = matrix(rnorm(S * K, 0, .05), S, K) * scale,
    yi = rnorm(K) * scale,
    vi = runif(K, .01, .1) * scale^2,
    cluster = rep(seq_len(K / 3L), each = 3L)
  )
}

.direct_arithmetic_runs <- function(x, bias = TRUE) {

  offset <- if (bias) x$bias else NULL
  list(
    blup = function() .evaluate.brma.multilevel_blup.norm(
      x$mu, x$tau_within, x$tau_between, x$yi, x$vi, x$cluster, offset),
    total = function() {
      set.seed(2)
      .evaluate.brma.multilevel_posterior.norm(
        x$mu, x$tau_within, x$tau_between, x$yi, x$vi, x$cluster, offset, "total")
    },
    cluster = function() {
      set.seed(2)
      .evaluate.brma.multilevel_posterior.norm(
        x$mu, x$tau_within, x$tau_between, x$yi, x$vi, x$cluster, offset, "cluster")
    },
    true_effects = function() .evaluate.brma.true_effects.norm(
      x$mu, x$tau_within, x$yi, sqrt(x$vi), TRUE, offset)
  )
}

# The direct helpers are spied so that the route taken is observable: a call
# counts when a helper returns a result (the true-effect helper returns NULL for
# inputs outside its range).
.direct_arithmetic_spy <- function(env = parent.frame()) {

  force(env)
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  for (name in c(".multilevel_blup_direct", ".multilevel_posterior_direct", ".true_effects_direct")) {
    local({
      original <- get(name, envir = asNamespace("RoBMA"))
      spy <- function(...) {
        result <- original(...)
        if (!is.null(result)) calls$n <- calls$n + 1L
        result
      }
      testthat::local_mocked_bindings(.package = "RoBMA", .env = env, !!name := spy)
    })
  }
  calls
}


test_that("in-range inputs take the direct route and agree with the checked route", {

  x <- .direct_arithmetic_fixture()
  calls <- .direct_arithmetic_spy()
  for (bias in c(FALSE, TRUE)) {
    runs <- .direct_arithmetic_runs(x, bias)
    for (name in names(runs)) {
      calls$n <- 0L
      direct <- runs[[name]]()
      expect_gt(calls$n, 0L, label = paste("direct calls of", name, "with bias", bias))
      testthat::with_mocked_bindings(
        checked <- runs[[name]](),
        .evaluate_direct_arithmetic = function(...) FALSE,
        .package = "RoBMA"
      )
      direct <- unlist(direct, use.names = FALSE)
      checked <- unlist(checked, use.names = FALSE)
      # Both routes round at the scale of their inputs: compare against the
      # largest magnitude, not elementwise near cancellations.
      expect_lt(max(abs(direct - checked)) / max(abs(checked)), 1e-13,
        label = paste(name, "with bias", bias))
    }
  }
})


test_that("inputs outside the direct range take the checked route", {

  x <- .direct_arithmetic_fixture()
  # Scaling all SDs, variances and residuals by c scales every contribution by c.
  scale <- 2^150
  big <- .direct_arithmetic_fixture(scale)
  calls <- .direct_arithmetic_spy()
  reference <- .direct_arithmetic_runs(x, TRUE)
  scaled <- .direct_arithmetic_runs(big, TRUE)
  for (name in names(scaled)) {
    calls$n <- 0L
    actual <- unlist(scaled[[name]](), use.names = FALSE)
    expect_identical(calls$n, 0L, label = paste(name, "outside the range"))
    calls$n <- 0L
    expected <- unlist(reference[[name]](), use.names = FALSE) * scale
    expect_gt(calls$n, 0L)
    expect_lt(max(abs(actual - expected)) / max(abs(expected)), 1e-12, label = name)
  }

  # One tiny SD, one huge SD, one huge value, one tiny residual and one zero
  # sampling SD each leave the range.
  cases <- list(
    tau_within = function(y) { y$tau_within[3L, 2L] <- 2^-120; y },
    tau_between = function(y) { y$tau_between[5L, 4L] <- 2^101; y },
    mu = function(y) { y$mu[7L, 1L] <- 2^105; y },
    residual = function(y) { y$yi[2L] <- 0; y$mu[, 2L] <- 2^-110; y },
    sampling = function(y) { y$vi[6L] <- 0; y$tau_within[, 6L] <- abs(y$tau_within[, 6L]) + 1; y }
  )
  for (name in names(cases)) {
    y <- cases[[name]](.direct_arithmetic_fixture())
    runs <- .direct_arithmetic_runs(y, FALSE)
    for (run in names(runs)) {
      # The true-effect evaluator reads neither the cluster SDs nor the clusters.
      if (identical(name, "tau_between") && identical(run, "true_effects")) next
      calls$n <- 0L
      suppressWarnings(try(runs[[run]](), silent = TRUE))
      expect_identical(calls$n, 0L, label = paste(name, run))
    }
  }
})


test_that("the direct range is zero or an absolute value within the bound", {

  bound <- EVALUATE_DIRECT_BOUND
  expect_true(.within_direct_range(c(0, 1, bound, 1 / bound), nonnegative = TRUE))
  expect_true(.within_direct_range(c(-bound, 0, -1 / bound, 3)))
  expect_false(.within_direct_range(c(0, bound * 2), nonnegative = TRUE))
  expect_false(.within_direct_range(c(0, 1 / bound / 2), nonnegative = TRUE))
  expect_false(.within_direct_range(c(-1, 1), nonnegative = TRUE))
  expect_false(.within_direct_range(c(1, -bound * 2)))
  expect_false(.within_direct_range(c(1, NA)))
  expect_false(.within_direct_range(c(1, Inf)))
  expect_false(.within_direct_range(numeric()))
  expect_false(.evaluate_direct_arithmetic(list(matrix(1, 2L, 3L)), c(1, 0, 1), list(matrix(1, 2L, 3L))))
  expect_false(.evaluate_direct_arithmetic(list(matrix(1, 2L, 3L)), c(1, 1), list(matrix(1, 2L, 3L))))
  expect_false(.evaluate_direct_arithmetic(list(matrix(1, 2L, 2L)), c(1, 1, 1), list(matrix(1, 2L, 3L))))
  expect_true(.evaluate_direct_arithmetic(list(matrix(1, 2L, 3L)), c(1, 1, 1), list(matrix(1, 2L, 3L))))
})


test_that("a direct product below the double range rounds to zero where the checked route refuses", {

  # One cluster of three estimates with every input inside the direct range.
  # tau_within = 2^100 with vi = 2^-200 gives a total contribution of about
  # 2^-1300 in the second estimate, which is not representable.
  mu_samples  <- matrix(0, 1L, 3L)
  tau_within  <- matrix(c(0, 2^100, 2^-100), 1L, 3L)
  tau_between <- matrix(c(2^100, 2^-100, 2^-100), 1L, 3L)
  vi          <- c(2^-200, 2^-200, 2^200)
  yi          <- c(0, 0, 2^-100)
  cluster     <- c(1L, 1L, 1L)
  evaluate <- function() .evaluate.brma.multilevel_blup.norm(
    mu_samples, tau_within, tau_between, yi, vi, cluster)

  direct <- evaluate()
  expect_true(all(is.finite(unlist(direct))))
  expect_identical(direct[["total"]][[2L]], 0)
  testthat::with_mocked_bindings(
    expect_error(evaluate(), "numerically unavailable"),
    .evaluate_direct_arithmetic = function(...) FALSE,
    .package = "RoBMA"
  )
})
