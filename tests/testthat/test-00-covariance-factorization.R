test_that("covariance policy classifies positive semidefinite matrices", {

  positive_definite <- .covariance_factorization(
    matrix(c(2, 0.25, 0.25, 1), nrow = 2L)
  )
  singular <- .covariance_factorization(
    matrix(c(1, 1, 1, 1), nrow = 2L)
  )
  indefinite <- .covariance_factorization(
    matrix(c(1, 2, 2, 1), nrow = 2L)
  )

  expect_identical(positive_definite[["status"]], "positive_definite")
  expect_true(.covariance_is_numerically_positive_definite(positive_definite))
  expect_identical(singular[["status"]], "positive_semidefinite")
  expect_true(.covariance_is_positive_semidefinite(singular))
  expect_identical(indefinite[["status"]], "indefinite")
  expect_false(.covariance_is_positive_semidefinite(indefinite))
})


test_that("covariance policy rejects materially invalid correlations", {

  correlation <- matrix(
    c(1, 1 + 1e-9, 1 + 1e-9, 1),
    nrow = 2L
  )
  factorization <- .covariance_factorization(correlation)

  expect_identical(factorization[["status"]], "indefinite")
  expect_false(.covariance_is_positive_semidefinite(factorization))
})


test_that("strict covariance policy rejects adjacent-above-one correlations", {

  correlation <- matrix(
    c(1, 1 + .Machine$double.eps, 1 + .Machine$double.eps, 1),
    nrow = 2L
  )
  factorization <- .covariance_factorization(correlation, strict = TRUE)

  expect_identical(factorization[["status"]], "indefinite")
  expect_false(.covariance_is_positive_semidefinite(factorization))
  expect_null(.covariance_sampling_factor(factorization))

  rho <- -0.5 - .Machine$double.eps
  near_indefinite <- matrix(rho, nrow = 3L, ncol = 3L)
  diag(near_indefinite) <- 1
  expect_identical(
    .covariance_factorization(near_indefinite, strict = TRUE)[["status"]],
    "indefinite"
  )
})


test_that("covariance sampling factor preserves singular covariance", {

  covariance    <- matrix(c(1, 1, 1, 1), nrow = 2L)
  factorization <- .covariance_factorization(covariance)
  factor        <- .covariance_sampling_factor(factorization)

  expect_equal(crossprod(factor), covariance, tolerance = 1e-12)
  expect_null(.covariance_cholesky(factorization))
})


test_that("covariance sampling factor preserves exact rank-one covariance", {

  factor          <- c(2^-40, 1, 2^40)
  covariance      <- tcrossprod(factor)
  factorization   <- .covariance_factorization(covariance, strict = TRUE)
  sampling_factor <- .covariance_sampling_factor(factorization)

  expect_identical(factorization[["status"]], "positive_semidefinite")
  expect_identical(crossprod(sampling_factor), covariance)
})


test_that("draw-specific covariance sampling preserves singleton dimensions", {

  mu_samples         <- matrix(c(1, 2), ncol = 1L)
  covariance_samples <- array(c(0, 4), dim = c(2L, 1L, 1L))

  set.seed(127)
  draws <- .outcome_rng.norm_known_v_covariance(
    mu_samples         = mu_samples,
    covariance_samples = covariance_samples
  )

  expect_identical(dim(draws), c(2L, 1L))
  expect_identical(draws[1L, 1L], 1)
  expect_true(is.finite(draws[2L, 1L]))
})
