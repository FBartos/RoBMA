test_that("known-V collinearity diagnostics are separate from physical variance scale", {

  for (covariance in list(diag(c(1, 1e-20)),
      matrix(c(1, .5e-10, .5e-10, 1e-20), 2L))) {
    expect_false(.known_v_covariance_classification(covariance)$singular)
    expect_silent(value <- .known_v_as_matrix(covariance))
    expect_identical(value, covariance)
  }

  # Exact eigenvalues are 1+rho and 1-rho > 0. Under the existing
  # numerical PSD policy the second direction is unresolved at this scale;
  # it cannot authorize a positive-definite conditional solve.
  rho <- 1 - .Machine$double.eps
  covariance <- matrix(c(1, rho, rho, 1), 2L)
  factorization <- .covariance_factorization(covariance)
  expect_identical(factorization$status, "positive_semidefinite")
  expect_null(.covariance_cholesky(factorization))
  expect_true(.known_v_covariance_classification(covariance)$singular)
  expect_warning(value <- .known_v_as_matrix(covariance), "positive semidefinite")
  expect_identical(value, covariance)
})
