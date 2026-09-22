context("Native covariance plan ownership")

test_that("partially initialized covariance plans can be finalized after rejection", {

  for (iteration in 1:3) {
    expect_error(.Call("RoBMA_known_v_covariance_plan_create", c(0, 0),
      diag(2L), list(list(type = "group")), list(1:2), PACKAGE = "RoBMA"),
      "missing 'model_matrix'")
    expect_error(.Call("RoBMA_known_v_covariance_plan_create", c(0, 0),
      diag(2L), list(), list(1L, 1:2), PACKAGE = "RoBMA"),
      "must not overlap")
    # The external-pointer finalizer owns rejected partial plans too. GC
    # exercises this path, while the successful plan checks later calls work.
    gc()
  }
  plan <- .Call("RoBMA_known_v_covariance_plan_create", c(0, 0), diag(2L),
    list(), list(1:2), PACKAGE = "RoBMA")
  actual <- .Call("RoBMA_known_v_covariance_plan_loglik", plan, c(0, 0),
    list(), c(0, 0), PACKAGE = "RoBMA")
  expect_equal(actual, 2 * dnorm(0, log = TRUE), tolerance = 1e-14)
})

test_that("factor sweeps preserve declared diagonal coefficient structure", {

  evaluate <- function(root) {

    .marglik_covariance_plan_factor_grid_loglik(
      cache = new.env(parent = emptyenv()), y = c(0, .1),
      means = matrix(c(0, 0), 1L), sampling_covariance = diag(2L),
      random_covariance_plans = list(group = list(
        type = "group", model_matrix = diag(2L), group_map = c(1L, 1L),
        coefficient_structure = "diagonal")),
      random_covariance_states = list(list(list(coefficient_factor = diag(2L)))),
      block_indices = list(1:2), extra_variances = matrix(0, 1L, 2L),
      update_grid = list(family = "factor", factor_index = 1L,
        component_index = 1L, coefficient_scale = matrix(c(1, 1), 1L),
        coefficient_cholesky = array(root, c(1L, 2L, 2L)),
        candidate_scale = matrix(.5)))
  }
  expect_true(all(is.finite(evaluate(diag(2L)))))
  expect_error(evaluate(matrix(c(1, .4, 0, 1), 2L)),
               "must remain diagonal")
})
