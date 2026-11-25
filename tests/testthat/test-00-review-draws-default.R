test_that("generic brma_samples draws preserve chains and iteration identity", {

  samples <- .new_brma_samples(matrix(c(1, 2, 10, 20, 3, 4, 30, 40), 4L,
    dimnames = list(NULL, c("a", "b"))), n_chains = 2L, n_iter = 2L,
    title = "Draws", component = "location", probs = c(0.025, 0.975), data = list())
  # the generic returns the draws_matrix of 4.0.0, with the chain counts
  for (converter in list(RoBMA::as_draws, posterior::as_draws, as_draws_matrix)) {
    flat <- converter(samples)
    expect_s3_class(flat, "draws_matrix")
    expect_identical(posterior::nchains(flat), 2L)
    expect_equal(posterior::niterations(flat), 2L)
    expect_identical(as.numeric(flat[, "a"]), c(1, 2, 10, 20))
  }
  draws <- as_draws_array(samples)
  expect_s3_class(draws, "draws_array")
  expect_identical(posterior::nchains(draws), 2L)
  expect_identical(posterior::niterations(draws), 2L)
  expect_identical(as.numeric(draws[, 1L, "a"]), c(1, 2))
  expect_identical(as.numeric(draws[, 2L, "a"]), c(10, 20))
  expect_identical(as.numeric(draws[, 2L, "b"]), c(30, 40))
})


test_that("non-brma draw conversion retains posterior generic dispatch", {

  values <- matrix(seq_len(24), 8, dimnames = list(NULL,
    c("scalar", "ordinary[1]", "ordinary[2]")))
  input <- coda::mcmc.list(coda::mcmc(values, start = 5, thin = 3),
    coda::mcmc(values + 100, start = 5, thin = 3))
  for (generic in c("as_draws", "as_draws_array", "as_draws_df",
      "as_draws_list", "as_draws_matrix", "as_draws_rvars")) {
    actual <- getExportedValue("RoBMA", generic)(input)
    expected <- getExportedValue("posterior", generic)(input)
    expect_identical(class(actual), class(expected))
    expect_identical(posterior::as_draws_array(actual), posterior::as_draws_array(expected))
  }
})
