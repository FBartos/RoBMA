test_that("LOO comparison retains model identities after ranking", {

  base_scores <- matrix(seq(-2, -1, length.out = 200L), 40L, 5L)
  fit_low <- loo::waic(base_scores)
  fit_high <- loo::waic(base_scores + 1)
  metadata <- list(unit = "estimate", retained_context = "remaining_data",
                   target = "estimate_log_score", data_hash = "same-data")
  attr(fit_low, "RoBMA_target") <- metadata
  attr(fit_high, "RoBMA_target") <- metadata
  result <- loo_compare.loo(fit_low, fit_high)
  expect_identical(rownames(result), c("fit_high", "fit_low"))
  expect_equal(unname(result["fit_low", "elpd_diff"]), -5)
  expect_identical(rownames(loo_compare.loo(fit_low, better = fit_high)),
                   c("better", "fit_low"))
})

test_that("deterministic LOO merge matches pointwise columns by name", {

  scores <- matrix(seq(-2, -1, length.out = 400L), 100L, 4L)
  scores[, 1L] <- -2
  variable <- loo::loo(scores[, -1L], r_eff = 1, save_psis = TRUE)
  baseline <- .loo_combine_deterministic_columns(variable, scores, rep(1, 4L),
                                                 c(TRUE, FALSE, FALSE, FALSE))
  variable$pointwise <- variable$pointwise[, rev(seq_len(ncol(variable$pointwise)))]
  actual <- .loo_combine_deterministic_columns(variable, scores, rep(1, 4L),
                                               c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(actual$pointwise, baseline$pointwise)
})

test_that("LOO model labels remain bounded for expressions and do.call values", {

  scores <- matrix(seq(-2, -1, length.out = 200L), 40L, 5L)
  low <- loo::waic(scores)
  high <- loo::waic(scores + 1)
  metadata <- list(unit = "estimate", retained_context = "remaining_data",
    target = "estimate_log_score", data_hash = "same-data")
  attr(low, "RoBMA_target") <- attr(high, "RoBMA_target") <- metadata
  expect_identical(rownames(do.call(loo_compare.loo, list(low, high))), c("model2", "model1"))
  expect_identical(rownames(do.call(loo_compare.loo, list(x = low, second = high))), c("second", "x"))
  expect_identical(rownames(loo_compare.loo(identity(low), identity(high))), c("model2", "model1"))
  expect_identical(rownames(loo_compare.loo(low, low = high)), c("low.1", "low"))
})
