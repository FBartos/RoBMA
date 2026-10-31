test_that("LOO comparison retains model identities after ranking", {

  base_scores <- matrix(seq(-2, -1, length.out = 200L), 40L, 5L)
  fit_low <- loo::waic(base_scores)
  fit_high <- loo::waic(base_scores + 1)
  metadata <- list(unit = "estimate", retained_context = "remaining_data",
                   target = "estimate_log_score", data_hash = "same-data")
  attr(fit_low, "RoBMA_target") <- metadata
  attr(fit_high, "RoBMA_target") <- metadata
  result <- loo_compare.loo(fit_low, fit_high)
  expect_identical(as.character(result[["model"]]), c("model2", "model1"))
  expect_equal(result[["elpd_diff"]][result[["model"]] == "model1"], -5)
  expect_identical(as.character(loo_compare.loo(fit_low, better = fit_high)[["model"]]),
                   c("better", "model1"))
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

test_that("LOO comparison labels follow loo for unnamed and named models", {

  scores <- matrix(seq(-2, -1, length.out = 200L), 40L, 5L)
  low <- loo::waic(scores)
  high <- loo::waic(scores + 1)
  metadata <- list(unit = "estimate", retained_context = "remaining_data",
    target = "estimate_log_score", data_hash = "same-data")
  attr(low, "RoBMA_target") <- attr(high, "RoBMA_target") <- metadata
  labels <- function(x) as.character(x[["model"]])
  expect_identical(labels(do.call(loo_compare.loo, list(low, high))), c("model2", "model1"))
  expect_identical(labels(do.call(loo_compare.loo, list(x = low, second = high))), c("second", "model1"))
  expect_identical(labels(loo_compare.loo(identity(low), identity(high))), c("model2", "model1"))
  expect_identical(labels(loo_compare.loo(low, low = high)), c("low", "model1"))
})
