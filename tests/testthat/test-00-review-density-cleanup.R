test_that("density aggregation refuses an undersized denominator", {

  expect_error(.iwmde_density_aggregate(matrix(0, 1L, 2L), 1, 1L),
    "IWMDE density aggregation requires a finite denominator at least as large as the number of terms.",
    fixed = TRUE)
})


test_that("construction and aggregation failures describe their actual evidence", {

  error <- tryCatch(.iwmde_stop_construction_failure(NULL, "mu", 2L,
    "baseline evaluation", "the baseline was unavailable"), error = identity)
  expect_s3_class(error, "iwmde_construction_error")
  expect_identical(conditionMessage(error), paste0(
    "Posterior density construction failed for target 'mu' at posterior row 2 ",
    "during baseline evaluation: the baseline was unavailable."))
  error <- tryCatch(.iwmde_stop_ordinate_numerical_failure(NA_real_, NA_real_),
    error = identity)
  expect_s3_class(error, "iwmde_ordinate_numerical_error")
  expect_identical(conditionMessage(error),
    "Posterior density aggregation failed on the ordinary scale.")
  expect_true(is.na(error[["log_ordinate"]]))
})


test_that("invalid conditional supports are distinguished from zero proposal density", {

  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows, ...) {
      matrix(c(1, 0), nrow = length(rows), ncol = 2L, byrow = TRUE)
    },
    .package = "RoBMA"
  )
  error <- tryCatch(.iwmde_chen_log_weight(list(), "mu", list(type = "primitive"),
    1L, .5, 1L, .5), error = identity)
  expect_s3_class(error, "iwmde_construction_error")
  if (inherits(error, "error")) {
    expect_identical(error[["stage"]], "proposal-support construction")
    expect_identical(error[["posterior_rows"]], 1L)
    expect_identical(conditionMessage(error), paste0(
      "IWMDE construction failed for target 'mu' at posterior row 1 during ",
      "proposal-support construction: the declared conditional target support ",
      "has invalid bounds."))
  }
})


test_that("grouped baseline failures retain the construction family and stage", {

  testthat::local_mocked_bindings(
    .iwmde_retained_location_row_states = function(...) NULL,
    .iwmde_row_states_grouped_marginal = function(...) stop("group evaluation failed"),
    .package = "RoBMA"
  )
  plan <- list(method = "iwmde", execution_spec = list(type = "primitive"),
    target = list(parameter = "mu"))
  error <- tryCatch(.iwmde_plan_baseline_contract(list(), plan, 1:2, c(0, 1)),
    error = identity)
  expect_s3_class(error, "iwmde_construction_error")
  expect_identical(conditionMessage(error), paste0(
    "IWMDE construction failed for target 'mu' at posterior rows 1, 2 during ",
    "grouped baseline joint-density evaluation: group evaluation failed."))
  expect_identical(error[["posterior_rows"]], 1:2)
})
