# ============================================================================ #
# Current loo comparison contract
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

test_that("loo_compare preserves the released numeric comparison table", {

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

  out      <- loo_compare(loo_a, loo_b)
  upstream <- do.call(
    get("loo_compare.default", envir = asNamespace("loo"), inherits = FALSE),
    list(loo_a, loo_b)
  )
  legacy_columns <- c(
    "elpd_diff", "se_diff", "elpd_loo", "se_elpd_loo",
    "p_loo", "se_p_loo", "looic", "se_looic"
  )

  expect_identical(class(out), c("compare.loo", "matrix", "array"))
  expect_true(is.numeric(out))
  expect_identical(colnames(out), legacy_columns)
  expect_identical(rownames(out), as.character(upstream[["model"]]))
  expect_equal(
    as.vector(out),
    as.vector(as.matrix(upstream[, legacy_columns, drop = FALSE]))
  )
  expect_true(all(is.finite(out)))
  expect_true("candidate" %in% rownames(loo_compare(loo_a, candidate = loo_b)))

  simplified     <- capture.output(print(out))
  full           <- capture.output(print(out, simplify = FALSE))
  upstream_print <- capture.output(print(upstream))
  expect_match(paste(simplified, collapse = "\n"), "elpd_diff")
  expect_false(any(grepl("elpd_loo", simplified, fixed = TRUE)))
  expect_true(any(grepl("elpd_loo", full, fixed = TRUE)))
  expect_true(any(grepl("model", upstream_print, fixed = TRUE)))
})


test_that("loo_compare preserves the released numeric WAIC table", {

  set.seed(12)
  log_lik_a <- matrix(stats::rnorm(2000, mean = -2, sd = 0.3), ncol = 10)
  log_lik_b <- log_lik_a + rep(seq(-0.02, 0.02, length.out = 10), each = 200)

  waic_a <- suppressWarnings(loo::waic(log_lik_a))
  waic_b <- suppressWarnings(loo::waic(log_lik_b))
  target <- list(
    unit               = "estimate",
    conditioning_depth = "estimate",
    data_hash          = "waic-contract",
    target             = "factorized_estimate"
  )
  attr(waic_a, "RoBMA_target") <- target
  attr(waic_b, "RoBMA_target") <- target

  out      <- loo_compare(waic_a, waic_b)
  upstream <- do.call(
    get("loo_compare.default", envir = asNamespace("loo"), inherits = FALSE),
    list(waic_a, waic_b)
  )
  legacy_columns <- c(
    "elpd_diff", "se_diff", "elpd_waic", "se_elpd_waic",
    "p_waic", "se_p_waic", "waic", "se_waic"
  )

  expect_identical(class(out), c("compare.loo", "matrix", "array"))
  expect_true(is.numeric(out))
  expect_identical(colnames(out), legacy_columns)
  expect_identical(rownames(out), as.character(upstream[["model"]]))
  expect_equal(
    as.vector(out),
    as.vector(as.matrix(upstream[, legacy_columns, drop = FALSE]))
  )
  expect_true(all(is.finite(out)))

  simplified <- capture.output(print(out))
  full       <- capture.output(print(out, simplify = FALSE))
  expect_false(any(grepl("elpd_waic", simplified, fixed = TRUE)))
  expect_true(any(grepl("elpd_waic", full, fixed = TRUE)))
})
