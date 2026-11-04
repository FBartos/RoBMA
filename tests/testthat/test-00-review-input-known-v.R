test_that("zero-rank latent blocks preserve later sampling-factor columns", {

  covariance <- matrix(0, 3L, 3L)
  covariance[1L, 1L] <- 0.1
  covariance[2:3, 2:3] <- tcrossprod(c(0.2, 0.3))
  known_V <- .known_v_prepare(
    covariance, keep_rows = rep(TRUE, 3L),
    known_v_parameterization = "latent", warn_singular = FALSE
  )
  factor <- .known_v_sampling_factor_plan(known_V)
  reconstructed <- factor[["sampling_covariance"]] +
    tcrossprod(factor[["factor_plan"]][["model_matrix"]])
  expect_equal(reconstructed, covariance, tolerance = 1e-15)
})

test_that("precision right-hand sides match the covariance-plan dimension", {

  plan <- list(extra_variances = matrix(1, nrow = 2L, ncol = 3L))
  for (rhs in list(matrix(1, 2L, 1L), matrix(1, 4L, 1L), 1:3,
                   matrix(1, 3L, 0L), matrix(NA_real_, 3L, 1L))) {
    expect_error(.known_v_covariance_plan_precision_rhs_batch(plan, rhs),
                 "^Known-V precision right-hand sides are invalid\\.$")
  }
})
