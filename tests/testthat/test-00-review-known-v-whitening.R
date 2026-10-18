.review_whitened_known_v <- function(mixed = FALSE) {

  L <- rbind(c(1, 0), c(0, 1), c(1, 1))
  V <- tcrossprod(L)
  if (mixed) {
    full <- diag(c(1, 1, 2, 2))
    full[1:3, 1:3] <- V
    V <- full
  }
  .known_v_prepare(V, rep(TRUE, nrow(V)), "whitened", warn_singular = FALSE)
}

test_that("accepted whitening zeros reach the sampling factor adapter", {

  for (mixed in c(FALSE, TRUE)) {
    known <- .review_whitened_known_v(mixed)
    bytes <- serialize(known, NULL)
    factor <- .known_v_sampling_factor_plan(known)
    expect_false(is.null(factor))
    if (is.null(factor)) next
    expect_identical(serialize(known, NULL), bytes)
    expect_equal(ncol(factor$factor_plan$model_matrix), 2L)
    reconstruction <- factor$sampling_covariance + tcrossprod(factor$factor_plan$model_matrix)
    expect_lt(max(abs(reconstruction - .known_v_covariance_matrix(known))), 5e-15)
    if (mixed) expect_identical(factor$sampling_covariance[4, 4], 2)
    K <- .known_v_nrow(known)
    for (y in list(rep(0, K), c(1e-8 * c(1, 1, -1), if (mixed) 0))) {
      actual <- .marglik_covariance_plan_loglik(NULL, y, rep(0, K),
        factor$sampling_covariance, list(factor$factor_plan), list(factor$factor_state),
        .known_v_block_indices(.known_v_covariance_matrix(known)), rep(1e-18, K))
      expected <- if (all(y == 0)) 17.417144092998338 else -132.58285590700166
      if (mixed) expected <- expected + dnorm(0, 0, sqrt(2), log = TRUE)
      expect_lt(abs(actual - expected), 1e-11)
    }
  }
  pd <- .known_v_prepare(matrix(c(1, .25, .25, 1), 2), rep(TRUE, 2), "whitened")
  expect_null(.known_v_sampling_factor_plan(pd))
  damaged <- .review_whitened_known_v()
  damaged$whitening_blocks[[1]]$index <- c(1L, 1L, 3L)
  expect_error(.known_v_sampling_factor_plan(damaged), "whitening")
})
