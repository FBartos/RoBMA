test_that("fixed-group precision is refined while location shortcuts decline", {

  groups <- c(1L, 1L, 2L, 2L)
  G <- matrix(c(1, .5, .5, 1), 2)
  factor_plan <- list(type = "known_group", model_matrix = matrix(1, 4, 1),
    group_map = groups, group_covariance = G, coefficient_structure = "diagonal")
  args <- list(cache = new.env(parent = emptyenv()), y = c(1, 1, -.5, -.5),
    means = matrix(0, 1, 4), sampling_covariance = diag(1e-14, 4),
    random_covariance_plans = list(factor_plan),
    random_covariance_states = list(list(list(coefficient_factor = matrix(1, 1, 1)))),
    block_indices = list(1:4), extra_variances = matrix(0, 1, 4))
  reference <- rep(as.numeric(solve(G + diag(1e-14 / 2, 2), c(1, -.5))) / 2, each = 2)
  actual <- do.call(.marglik_covariance_plan_precision_residual_batch, args)
  expect_lt(max(abs(actual - reference)), 1e-12)
  location <- do.call(.marglik_covariance_plan_location_quadratic_batch,
    c(args, list(bases = matrix(c(1, -1, 0, 0), 1))))
  expect_null(location)
  args$sampling_covariance <- diag(.1, 4)
  args$cache <- new.env(parent = emptyenv())
  ordinary <- do.call(.marglik_covariance_plan_location_quadratic_batch,
    c(args, list(bases = matrix(1, 1, 4))))
  covariance <- G[groups, groups] + diag(.1, 4)
  expect_equal(ordinary$linear, sum(solve(covariance, args$y)), tolerance = 1e-12)
  expect_equal(ordinary$quadratic, sum(solve(covariance, rep(1, 4))), tolerance = 1e-12)
})

test_that("spectral precision requests use the existing dense availability gate", {

  V <- matrix(c(1, 1 - 1e-12, 1 - 1e-12, 1), 2)
  args <- list(cache = NULL, y = c(0, 1), means = matrix(0, 1, 2),
    sampling_covariance = V, random_covariance_plans = list(),
    random_covariance_states = list(list()), block_indices = list(1:2),
    extra_variances = matrix(0, 1, 2))
  expect_error(do.call(.marglik_covariance_plan_precision_residual_batch, args),
    class = "RoBMA_known_v_precision_unavailable")
  score <- do.call(.marglik_covariance_plan_conditional_loglik_batch, args)
  expect_true(all(is.finite(score)))
})
