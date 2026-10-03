test_that("integer Gram support does not depend on rounded Cholesky pivots", {

  # C = U U' has exact integer entries and rank two: U's two columns are
  # independent and row1=row2+row3. Its nonzero eigenvalues are exactly 3,1.
  loading <- rbind(c(1, 1), c(1, 0), c(0, 1))
  null <- c(1, -1, -1)
  permutations <- list(c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
    c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L))
  for (index in permutations) {
    covariance <- tcrossprod(loading[index, , drop = FALSE])
    value <- .covariance_factorization(covariance)
    root <- .covariance_sampling_factor(value)
    expect_identical(value$status, "positive_semidefinite")
    expect_identical(sum(value$spectral_values > 0), 2L)
    expect_equal(value$spectral_values[1:2], c(3, 1), tolerance = 1e-12)
    expect_null(.covariance_cholesky(value))
    expect_true(sum((root %*% null[index])^2) < 1e-26)
    expect_equal(crossprod(root), covariance, tolerance = 1e-12)
  }

  # The same exact singular total must not enter a PD conditional solve.
  latent <- tcrossprod(loading[, 1L])
  sampling <- tcrossprod(loading[, 2L])
  expect_error(.zplot_gaussian_conditional_variance(latent, sampling),
    "The estimate-depth Gaussian conditional variance is unavailable because positive definiteness of the marginal covariance cannot be resolved at working precision.", fixed = TRUE)
})

test_that("retained integer Gram sampling uses the accepted compact support", {

  loading <- rbind(c(1, 1), c(1, 0), c(0, 1))
  null <- c(1, -1, -1)
  permutations <- list(c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
    c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L))
  for (index in permutations) {
    covariance <- tcrossprod(loading[index, , drop = FALSE])
    accepted <- .covariance_factorization(covariance)$sampling_factor
    known_V <- .known_v_prepare(covariance, rep(TRUE, 3L), "block_mvn",
      warn_singular = FALSE)
    sampling <- .known_v_check_selection_structure(
      .known_v_resolve_selection_structure(known_V))
    block <- sampling$latent_blocks[[1L]]
    expect_identical(sampling$rank, 2L)
    expect_identical(block$rank, 2L)
    expect_identical(dim(block$B), c(3L, 2L))
    expect_identical(block$z_start, 1L)
    expect_identical(block$z_end, 2L)
    expect_equal(block$B, unname(t(accepted)), tolerance = 0)
    expect_true(sum((crossprod(block$B, null[index]))^2) < 1e-26)
    expect_equal(tcrossprod(block$B), covariance, tolerance = 1e-12)
    expect_identical(sampling$residual_variance, numeric(3L))
    expect_identical(sampling$provenance$policy, "whole_sampling_error")
  }
})
