# ============================================================================ #
# Current loo comparison contract
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

test_that("loo_compare returns the current mixed-type comparison table", {

  set.seed(11)
  log_lik_a <- matrix(stats::rnorm(2000, mean = -2, sd = 0.3), ncol = 10)
  log_lik_b <- log_lik_a + rep(seq(-0.02, 0.02, length.out = 10), each = 200)

  loo_a <- suppressWarnings(loo::loo(log_lik_a))
  loo_b <- suppressWarnings(loo::loo(log_lik_b))
  target <- list(
    unit               = "estimate",
    conditioning_depth = "estimate",
    data_hash          = "loo-contract",
    target             = "factorized_estimate"
  )
  attr(loo_a, "RoBMA_target") <- target
  attr(loo_b, "RoBMA_target") <- target

  out <- loo_compare(loo_a, loo_b)

  expect_s3_class(out, "compare.loo")
  expect_s3_class(out, "data.frame")
  expect_true(all(c(
    "model", "elpd_diff", "se_diff", "p_worse", "diag_diff", "diag_elpd"
  ) %in% colnames(out)))
  expect_true(all(is.finite(out[["elpd_diff"]])))
  expect_true(all(is.finite(out[["se_diff"]])))
  expect_true(is.na(out[["p_worse"]][[1L]]))
  expect_type(out[["diag_diff"]], "character")
  expect_type(out[["diag_elpd"]], "character")
})
