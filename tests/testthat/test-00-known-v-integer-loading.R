test_that("integer sampling factors share double identity and native precision solves", {

  integer_loading <- matrix(1L, 2L, 1L, dimnames = list(c("a", "b"), "factor"))
  double_loading <- integer_loading
  storage.mode(double_loading) <- "double"
  integer_factor <- known_v_factor(c(.5, .5), integer_loading)
  double_factor <- known_v_factor(c(.5, .5), double_loading)
  expect_identical(integer_factor, double_factor)
  expect_type(integer_factor$loading, "double")
  expect_identical(dim(integer_factor$loading), c(2L, 1L))
  expect_null(dimnames(integer_factor$loading))
  integer_known <- .known_v_prepare(integer_factor, keep_rows = rep(TRUE, 2L),
    known_v_parameterization = "auto", warn_singular = FALSE)
  double_known <- .known_v_prepare(double_factor, keep_rows = rep(TRUE, 2L),
    known_v_parameterization = "auto", warn_singular = FALSE)
  expect_identical(integer_known, double_known)
  rhs <- matrix(c(.3, -.5), 2L, 1L)
  for (known in list(double_known, integer_known)) {
    factor <- .known_v_sampling_factor_plan(known)
    plan <- list(sampling_covariance = factor$sampling_covariance,
      random_covariance_plans = list(factor$factor_plan),
      random_covariance_states = list(list(factor$factor_state)),
      block_indices = list(1:2), extra_variances = matrix(0, 1L, 2L))
    actual <- .known_v_covariance_plan_precision_rhs_batch(plan, rhs)
    # Independent dense solve: V = .5 I + 11'; expected (.76, -.84).
    reference <- solve(diag(.5, 2L) + matrix(1, 2L, 2L), rhs)
    expect_equal(as.numeric(actual), as.numeric(reference), tolerance = 1e-14)
    expect_equal(as.numeric(actual), c(.76, -.84), tolerance = 1e-14)
  }
  expect_identical(known_v_factor(c(0, .5), matrix(integer(), 2L, 0L)),
    known_v_factor(c(0, .5), matrix(numeric(), 2L, 0L)))
  expect_identical(known_v_factor(c(0, .5), matrix(0L, 2L, 1L)),
    known_v_factor(c(0, .5), matrix(0, 2L, 1L)))
})

test_that("sampling loading validation still rejects unsupported matrices", {

  for (loading in list(c(1, 1), matrix("1", 2L, 1L), matrix(NA_real_, 2L, 1L),
    matrix(Inf, 2L, 1L), matrix(1, 3L, 1L))) {
    error <- tryCatch(known_v_factor(c(.5, .5), loading), error = identity)
    expect_s3_class(error, "error")
    expect_identical(conditionMessage(error),
      "'loading' must be a finite numeric matrix with one row per element of 'diagonal'.")
  }
})
