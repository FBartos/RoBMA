test_that("generic brma_samples draws preserve chains and iteration identity", {

  samples <- .new_brma_samples(matrix(c(1, 2, 10, 20, 3, 4, 30, 40), 4L,
    dimnames = list(NULL, c("a", "b"))), n_chains = 2L, n_iter = 2L,
    title = "Draws", component = "location", probs = c(0.025, 0.975), data = list())
  for (converter in list(RoBMA::as_draws, posterior::as_draws)) {
    draws <- converter(samples)
    expect_s3_class(draws, "draws_array")
    expect_identical(posterior::nchains(draws), 2L)
    expect_identical(posterior::niterations(draws), 2L)
    expect_identical(as.numeric(draws[, 1L, "a"]), c(1, 2))
    expect_identical(as.numeric(draws[, 2L, "a"]), c(10, 20))
    expect_identical(as.numeric(draws[, 2L, "b"]), c(30, 40))
  }
  flat <- as_draws_matrix(samples)
  expect_s3_class(flat, "draws_matrix")
  expect_identical(posterior::nchains(flat), 2L)
  expect_identical(as.numeric(flat[, "a"]), c(1, 2, 10, 20))
})
