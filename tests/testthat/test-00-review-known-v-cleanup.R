test_that("null-space coverage refuses an indefinite covariance factorization", {

  indefinite <- .covariance_factorization(matrix(c(1, 2, 2, 1), 2L))
  expect_false(.covariance_is_positive_semidefinite(indefinite))
  expect_false(.known_v_nullspace_is_covered(indefinite, diag(2L)))
  positive <- .covariance_factorization(diag(2L))
  expect_true(.known_v_nullspace_is_covered(positive, diag(2L)))
  singular <- .covariance_factorization(matrix(1, 2L, 2L))
  expect_true(.known_v_nullspace_is_covered(singular, diag(2L)))
})
