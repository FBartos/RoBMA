test_that("conditional Gaussian variances preserve the scalar harmonic law", {

  for (variances in list(c(1, 1e-20), c(1e-20, 1), c(1, 1), c(1e-20, 1e-20))) {
    latent <- variances[[1L]]
    sampling <- variances[[2L]]
    expected <- latent * sampling / (latent + sampling)
    observed <- .zplot_gaussian_conditional_variance(matrix(latent), matrix(sampling))
    expect_gt(observed, 0)
    expect_equal(observed / expected, 1, tolerance = 1e-12)
  }
})

test_that("conditional Gaussian variance agrees with independent precision addition", {

  latent <- matrix(c(1, .5, .5, 2), 2L)
  for (sampling in list(matrix(c(.4, .1, .1, .3), 2L), diag(c(1e-20, 3e-20)))) {
    expected <- diag(solve(solve(latent) + solve(sampling)))
    observed <- .zplot_gaussian_conditional_variance(latent, sampling)
    expect_equal(observed / expected, c(1, 1), tolerance = 1e-12)
  }

  # The first latent coordinate is observed without noise. The remaining
  # latent conditional variance is 2 - .5^2; combine it with noise variance .3.
  sampling <- diag(c(0, .3))
  observed <- .zplot_gaussian_conditional_variance(latent, sampling)
  expect_identical(observed[[1L]], 0)
  expect_equal(observed[[2L]], 1.75 * .3 / (1.75 + .3), tolerance = 1e-12)
  expect_identical(.zplot_gaussian_conditional_variance(diag(2L), matrix(0, 2L, 2L)), c(0, 0))
  expect_identical(.zplot_gaussian_conditional_variance(matrix(0, 2L, 2L), diag(2L)), c(0, 0))
})

test_that("conditional Gaussian variance keeps invalid covariance failures explicit", {

  expect_error(.zplot_gaussian_conditional_variance(matrix(0), matrix(0)),
    "The estimate-depth Gaussian conditional variance is unavailable because positive definiteness of the marginal covariance cannot be resolved at working precision.", fixed = TRUE)
  expect_error(.zplot_gaussian_conditional_variance(matrix(-1), matrix(2)),
    "Latent and sampling covariance matrices must be positive semidefinite.", fixed = TRUE)
})
