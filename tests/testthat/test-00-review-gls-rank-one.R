.review_rank_one_gls <- function(X, y, u, D, full = TRUE) {

  known <- .known_v_prepare(known_v_factor(rep(0, length(u)), matrix(u, ncol = 1)),
    rep(TRUE, length(u)), "whitened", warn_singular = FALSE)
  .known_v_gls_projection_blocks(X, y, known, D, return_full_H = full)
}

test_that("rank-one sampling residuals retain tiny independent variance", {

  for (K in c(2L, 4L)) {
    X <- if (K == 2L) matrix(c(1, 2), 2L) else cbind(1, 1:4)
    u <- if (K == 2L) c(1, 1) else c(1, 1, -1, -1)
    for (d in c(2^-60, .1)) {
      actual <- .review_rank_one_gls(X, u, u, rep(d, K))
      multiplier <- if (K == 2L) 1 / (1 + 5 * d) else 4 / (4 + 5 * d)
      expect_equal(actual$sampling_residual, u * multiplier, tolerance = 1e-12)
      expect_equal(actual$sampling_residual_variance, rep(multiplier, K), tolerance = 1e-12)
      expect_equal(actual$H_diag, diag(actual$H), tolerance = 1e-12)
      expect_equal(.review_rank_one_gls(X, u, u, rep(d, K), full = FALSE)$H_diag,
        actual$H_diag, tolerance = 1e-12)
    }
  }
})

test_that("sampling residual factor products retain a global mixed-block projection", {

  V <- matrix(0, 6L, 6L)
  V[1:4, 1:4] <- tcrossprod(c(1, 1, -1, -1))
  V[5:6, 5:6] <- matrix(c(1, .2, .2, 2), 2L)
  known <- .known_v_prepare(V, rep(TRUE, 6L), "whitened", warn_singular = FALSE)
  X <- cbind(1, seq_len(6L))
  y <- c(1, 1, -1, -1, .3, -.2)
  D <- c(.1, .2, .3, .4, .15, .25)
  actual <- .known_v_gls_projection_blocks(X, y, known, D, TRUE)
  M <- V + diag(D)
  W <- chol2inv(chol(M))
  B <- chol2inv(chol(crossprod(X, W %*% X)))
  H <- X %*% B %*% crossprod(X, W)
  A <- V %*% W %*% (diag(6L) - H)
  expect_equal(actual$sampling_residual, as.vector(A %*% y), tolerance = 1e-12)
  expect_equal(actual$sampling_residual_variance, diag(A %*% M %*% t(A)), tolerance = 1e-12)
  expect_equal(actual$H, H, tolerance = 1e-12)
  expect_equal(actual$H_diag, diag(H), tolerance = 1e-12)
})

test_that("rank-one products avoid a subtractive precision matrix and retain D0 controls", {

  testthat::local_mocked_bindings(.hat_precision_matrix = function(...) stop("subtractive precision"),
    .package = "RoBMA")
  expect_true(all(is.finite(.review_rank_one_gls(matrix(c(1, 2), 2L), c(1, 1),
    c(1, 1), c(.1, .2))$sampling_residual)))
  known <- .known_v_prepare(known_v_factor(c(0, 0), matrix(1, 2L, 1L)),
    c(TRUE, TRUE), "whitened", warn_singular = FALSE)
  actual <- .known_v_gls_projection_blocks(matrix(1, 2L), c(1, -1), known, c(0, .1))
  expect_equal(actual$sampling_residual, c(0, 0), tolerance = 1e-12)
  expect_equal(actual$sampling_residual_variance, c(0, 0), tolerance = 1e-12)
  zero <- .known_v_prepare(known_v_factor(c(1, 1), matrix(0, 2L, 1L)),
    c(TRUE, TRUE), "whitened")
  actual <- .known_v_gls_projection_blocks(matrix(1, 2L), c(1, -1), zero, c(.1, .2))
  expect_true(all(is.finite(actual$sampling_residual_variance)))
})
