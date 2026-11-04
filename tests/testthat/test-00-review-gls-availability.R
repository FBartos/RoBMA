test_that("estimate projections use the existing precision availability owner", {

  withr::local_options(scipen = 0)
  rho <- 1 - 1e-11
  V <- matrix(c(1, rho, rho, 1), 2L)
  known <- .known_v_prepare(V, rep(TRUE, 2L), "whitened", warn_singular = FALSE)
  for (extra in list(c(0, 0), c(0, 1e-30))) {
    outcome <- tryCatch(.known_v_gls_projection_blocks(matrix(1, 2L), c(1, -1), known, extra),
      RoBMA_known_v_precision_unavailable = identity)
    expect_s3_class(outcome, "RoBMA_known_v_precision_unavailable")
    if (inherits(outcome, "condition")) {
      expect_identical(outcome$rows, 1L)
      expect_true(outcome$reciprocal_condition < 1e-10)
      expect_identical(conditionMessage(outcome), paste0(
        "The precision application of the known covariance is unavailable: ",
        "the known sampling covariance with its random effects is too ",
        "ill-conditioned at draw 1 among the posterior draws it was applied to ",
        "(reciprocal condition number 5e-12). ",
        "Inspect the known covariance for near-singular rows."))
    }
  }
  expect_error(.known_v_gls_projection_blocks(matrix(1, 2L), c(1, -1), known, c(-1, 0)),
    "'extra_variance' must contain finite, non-negative variances.", fixed = TRUE)
  expect_error(.known_v_gls_projection_blocks(matrix(1, 2L), c(1, -1), known, c(Inf, 0)),
    "'extra_variance' must contain finite, non-negative variances.", fixed = TRUE)
  ordinary <- .known_v_gls_projection_blocks(matrix(1, 2L), c(1, -1), known, c(.1, .2))
  expect_true(all(is.finite(ordinary$sampling_residual)))
})
