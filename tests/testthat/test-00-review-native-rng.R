context("Native selection random-number ordering")

test_that("kernel RNG consumes the bin uniform before the interval uniform", {

  spec <- .test_step_spec(0, 1)
  withr::local_seed(409L)
  uniforms <- runif(12L)
  expected <- matrix(qnorm(uniforms[seq.int(2L, 12L, by = 2L)]), 6L)
  set.seed(409L)
  actual <- .selnorm_kernel_rng_matrix(matrix(0, 6L), matrix(1, 6L), 1,
    matrix(1, 6L, spec$n_bins), spec, kernel_mode = rep(0L, 6L))
  expect_equal(actual, expected, tolerance = 1e-14)
})
