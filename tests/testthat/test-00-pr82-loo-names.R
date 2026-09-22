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
