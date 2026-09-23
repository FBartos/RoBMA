test_that("small positive covariance directions survive structural singularity", {

  small <- 1e-20
  cases <- list(diag(c(1, small, 0)),
    matrix(c(1, 1, 0, 1, 1, 0, 0, 0, small), 3L))
  coordinates <- c(2L, 3L)
  for (index in seq_along(cases)) {
    covariance <- cases[[index]]
    value <- .covariance_factorization(covariance)
    factor <- .covariance_sampling_factor(value)
    expect_false(is.null(factor))
    expect_identical(dim(factor), dim(covariance))
    observed <- crossprod(factor)[coordinates[index], coordinates[index]]
    # Relative to this coordinate's own variance: a global absolute matrix
    # tolerance would incorrectly accept erasing the whole small direction.
    expect_gt(observed, 0)
    expect_equal(observed / small, 1, tolerance = 1e-12)
  }
  expect_true(.known_v_nullspace_is_regularized(diag(c(1, small, 0)), c(FALSE, FALSE, TRUE)))
  expect_false(.known_v_nullspace_is_regularized(diag(c(1, small, 0)), c(FALSE, TRUE, FALSE)))
})

test_that("accepted PSD support preserves exact rank and a separate small source", {

  # This dyadic Gram matrix is exactly rank two: its two factor columns are
  # orthogonal, with squared lengths 4 and 1. Its remaining directions are null.
  source <- cbind(rep(1, 4L), c(.5, -.5, .5, -.5))
  covariance <- tcrossprod(source)
  value <- .covariance_factorization(covariance)
  factor <- .covariance_sampling_factor(value)
  expect_identical(sum(value$spectral_values == 0), 2L)
  expect_equal(value$spectral_values[1:2], c(4, 1), tolerance = 1e-12)
  null <- cbind(c(1, 0, -1, 0), c(0, 1, 0, -1))
  expect_true(all(colSums((factor %*% null)^2) < 1e-26))
  expect_equal(crossprod(factor), covariance, tolerance = 1e-12)

  # Independent coordinate three is unaffected by observing the shared pair.
  covariance <- matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1e-20), 3L)
  conditional <- .selection_deleted_gaussian_block(covariance, c(.3, .3, 0), 3L)
  expect_equal(conditional$mean, 0)
  expect_equal(conditional$covariance[1L, 1L] / 1e-20, 1, tolerance = 1e-12)
})

test_that("successful Cholesky defines positive definiteness even for small variances", {

  covariance <- diag(c(1, 1e-20))
  value <- .covariance_factorization(covariance)
  expect_identical(value$status, "positive_definite")
  expect_false(value$singular)
  expect_identical(.covariance_cholesky(value), chol(covariance))
  expect_identical(.covariance_sampling_factor(value), chol(covariance))
  expect_true(.known_v_nullspace_is_regularized(covariance, c(FALSE, FALSE)))
})

test_that("covariance factorization does not repair asymmetric inputs", {

  covariance <- matrix(c(1, .4, .4, 1), 2L)
  covariance[1L, 2L] <- covariance[1L, 2L] + .Machine$double.eps
  expect_error(.covariance_factorization(covariance),
    "Covariance must be symmetric.", fixed = TRUE)
})

test_that("exact PSD variance necessities precede the eigenvalue tolerance", {

  withr::local_seed(123)
  invalid <- list(diag(c(1, -1e-20, 0)), matrix(c(0, 1e-20, 1e-20, 1), 2L))
  for (covariance in invalid) {
    value <- .covariance_factorization(covariance)
    expect_identical(value$status, "indefinite")
    expect_identical(value$covariance, covariance)
    expect_false(.covariance_is_positive_semidefinite(value))
    expect_null(.covariance_sampling_factor(value))
  }
})
