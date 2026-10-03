.known_v_score_range_values <- function(y, covariance, blocks) {

  for (symbol in c("RoBMA_known_v_covariance_plan_create",
    "RoBMA_known_v_covariance_plan_conditional_loglik",
    "RoBMA_known_v_covariance_plan_conditional_loglik_batch",
    "RoBMA_known_v_covariance_plan_conditional_summary_batch",
    "RoBMA_known_v_covariance_plan_precision_residual_batch")) {
    expect_s3_class(getNativeSymbolInfo(symbol, PACKAGE = "RoBMA"),
      "NativeSymbolInfo")
  }
  K <- length(y)
  plan <- .Call("RoBMA_known_v_covariance_plan_create", as.numeric(y),
    covariance, list(), blocks, PACKAGE = "RoBMA")
  mean <- numeric(K)
  extra <- numeric(K)
  list(
    scalar = .Call("RoBMA_known_v_covariance_plan_conditional_loglik",
      plan, mean, list(), extra, PACKAGE = "RoBMA"),
    batch = .Call("RoBMA_known_v_covariance_plan_conditional_loglik_batch",
      plan, matrix(mean, K), list(list()), matrix(extra, K), PACKAGE = "RoBMA"),
    summary = .Call("RoBMA_known_v_covariance_plan_conditional_summary_batch",
      plan, matrix(mean, K), list(list()), matrix(extra, K), PACKAGE = "RoBMA"),
    precision_residual = .Call("RoBMA_known_v_covariance_plan_precision_residual_batch",
      plan, matrix(mean, K), list(list()), matrix(extra, K), PACKAGE = "RoBMA")
  )
}

test_that("conditional scores standardize finite precision residuals before squaring", {

  covariance <- 1e-200 * matrix(c(1, .5, .5, 1), 2L)
  actual <- .known_v_score_range_values(c(1, 0), covariance, list(1:2))
  residual <- c(1, -.5)
  variance <- rep(.75e-200, 2L)
  reference <- stats::dnorm(residual, sd = sqrt(variance), log = TRUE)
  expect_true(all(is.finite(actual$scalar)))
  expect_lte(max(abs(actual$scalar / reference - 1)), 1e-12)
  expect_equal(as.numeric(actual$batch), actual$scalar, tolerance = 0)
  expect_equal(as.numeric(actual$summary$residual), residual, tolerance = 1e-14)
  expect_lte(max(abs(as.numeric(actual$summary$variance) / variance - 1)), 1e-14)
  expected_precision_residual <- residual / variance
  expect_lte(max(abs(as.numeric(actual$precision_residual) /
    expected_precision_residual - 1)), 1e-14)
})

test_that("singleton scores separate Gaussian constants at huge finite variance", {

  actual <- .known_v_score_range_values(0, matrix(1e308, 1L), list(1L))
  reference <- stats::dnorm(0, sd = 1e154, log = TRUE)
  expect_true(is.finite(actual$scalar))
  expect_equal(actual$scalar, reference, tolerance = 1e-14)
  expect_equal(as.numeric(actual$batch), reference, tolerance = 1e-14)
  expect_identical(as.numeric(actual$summary$residual), 0)
  expect_identical(as.numeric(actual$summary$variance), 1e308)
  expect_identical(as.numeric(actual$precision_residual), 0)
})

test_that("finite Gaussian half squares survive an overflowing unhalved square", {

  z <- 1.5e154
  # Explicit left association keeps the representable half square finite.
  reference <- -.5 * log(2 * pi) - (.5 * z) * z
  expect_true(is.finite(reference))
  expect_equal(reference, stats::dnorm(z, log = TRUE), tolerance = 1e-14)
  actual <- .known_v_score_range_values(z, matrix(1, 1L), list(1L))
  expect_true(is.finite(actual$scalar))
  expect_lte(abs(actual$scalar / reference - 1), 1e-14)
  expect_equal(as.numeric(actual$batch), actual$scalar, tolerance = 0)
  expect_identical(as.numeric(actual$summary$residual), z)
  expect_identical(as.numeric(actual$summary$variance), 1)
  expect_identical(as.numeric(actual$precision_residual), z)

  overflow <- .known_v_score_range_values(2e154, matrix(1, 1L), list(1L))
  expect_identical(overflow$scalar, -Inf)
})

test_that("rank-zero dense covariance scores and ordinary controls preserve their laws", {

  # Factory rank zero cannot enter diagonal_low_rank_conditional_values;
  # this admitted multirow plan exercises its dense covariance fallback.
  actual <- .known_v_score_range_values(c(0, 0), diag(1e308, 2L), list(1:2))
  reference <- rep(stats::dnorm(0, sd = 1e154, log = TRUE), 2L)
  expect_equal(actual$scalar, reference, tolerance = 1e-14)
  expect_equal(as.numeric(actual$summary$variance), rep(1e308, 2L), tolerance = 1e-14)

  covariance <- matrix(c(1, .5, .5, 1), 2L)
  actual <- .known_v_score_range_values(c(1, 0), covariance, list(1:2))
  expect_equal(actual$scalar, stats::dnorm(c(1, -.5), sd = sqrt(.75), log = TRUE),
    tolerance = 1e-14)
  expect_equal(as.numeric(actual$summary$residual), c(1, -.5), tolerance = 1e-14)
  expect_equal(as.numeric(actual$summary$variance), rep(.75, 2L), tolerance = 1e-14)
})
