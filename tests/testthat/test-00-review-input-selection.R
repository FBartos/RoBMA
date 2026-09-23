test_that("sensitivity rejects partitions that split normalization events", {

  expect_error(.selection_sensitivity_run(
    means = matrix(0, 1L, 2L), covariance = .block_covariance_zero(1L, 2L),
    context_covariance = .block_covariance_zero(1L, 2L),
    selection_sei = c(1, 1), selection_context = list(use_normal = TRUE),
    normalization_units = list(1:2), row_blocks = list(1L, 2L),
    execution_plan = list(), latent_samples = 5L
  ), "^Selection sensitivity blocks must contain complete normalization events\\.$")
})

