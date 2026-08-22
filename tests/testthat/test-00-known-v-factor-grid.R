.factor_grid_dense_loglik <- function(outcome, mean, covariance) {

  root <- chol(covariance)
  residual <- outcome - mean
  solved <- backsolve(root, forwardsolve(t(root), residual))
  -0.5 * (
    length(outcome) * log(2 * pi) +
      2 * sum(log(diag(root))) +
      sum(residual * solved)
  )
}


test_that("compact factor and Markov grids match dense covariance likelihoods", {

  outcome  <- c(0.2, -0.4, 0.5, 0.1)
  means    <- rbind(c(0.1, -0.2, 0.3, 0), c(-0.1, 0.2, 0.1, 0.15))
  sampling <- matrix(c(
    0.50, 0.08, 0.00, 0.00,
    0.08, 0.60, 0.04, 0.00,
    0.00, 0.04, 0.55, 0.06,
    0.00, 0.00, 0.06, 0.65
  ), nrow = 4L, byrow = TRUE)
  plan_cache <- new.env(parent = emptyenv())
  scales <- rbind(
    c(0.35, 0.45, 0.55, 0.65),
    c(0.50, 0.40, 0.60, 0.30)
  )
  factor_plan <- list(
    type = "group",
    model_matrix = diag(4L),
    group_map = rep(1L, 4L),
    coefficient_structure = "markov"
  )
  baseline_states <- lapply(seq_len(nrow(scales)), function(draw) {
    list(list(
      coefficient_factor = diag(scales[draw, ]),
      coefficient_scale = scales[draw, ],
      markov_transition = rep(0, 3L),
      markov_innovation_variance = rep(1, 3L)
    ))
  })
  rho <- c(-0.25, 0.2, 0.65)
  correlation <- lapply(rho, function(value) {
    value^abs(outer(seq_len(4L), seq_len(4L), "-"))
  })
  cholesky <- array(NA_real_, dim = c(length(rho), 4L, 4L))
  transition <- innovation <- matrix(NA_real_, length(rho), 3L)
  for (grid in seq_along(rho)) {
    current <- t(chol(correlation[[grid]]))
    cholesky[grid, , ] <- current
    transition[grid, ] <-
      current[cbind(2:4, 1:3)] / diag(current)[1:3]
    innovation[grid, ] <- diag(current)[2:4]^2
  }
  update <- list(
    family = "markov",
    factor_index = 1L,
    coefficient_scale = scales,
    candidate_cholesky = cholesky,
    candidate_transition = transition,
    candidate_innovation_variance = innovation
  )

  observed <- .marglik_covariance_plan_factor_grid_loglik(
    cache = plan_cache,
    y = outcome,
    means = means,
    sampling_covariance = sampling,
    random_covariance_plans = list(time = factor_plan),
    random_covariance_states = baseline_states,
    block_indices = list(1:4),
    extra_variances = matrix(0, nrow(means), ncol(means)),
    update_grid = update
  )
  expected <- matrix(NA_real_, length(rho), nrow(means))
  for (grid in seq_along(rho)) {
    for (draw in seq_len(nrow(means))) {
      covariance <- sampling +
        correlation[[grid]] * tcrossprod(scales[draw, ])
      expected[grid, draw] <- .factor_grid_dense_loglik(
        outcome,
        means[draw, ],
        covariance
      )
    }
  }
  expect_equal(observed, expected, tolerance = 1e-12)

  correlation_draws <- array(NA_real_, dim = c(2L, 4L, 4L))
  for (draw in 1:2) {
    correlation_draws[draw, , ] <- t(chol(
      matrix(0.15 * draw, 4L, 4L) + diag(1 - 0.15 * draw, 4L)
    ))
  }
  candidate_scale <- rbind(c(0.2, 0.3), c(0.7, 0.8), c(1.1, 0.9))
  factor_update <- list(
    family = "factor",
    factor_index = 1L,
    component_index = 2L,
    coefficient_scale = scales,
    coefficient_cholesky = correlation_draws,
    candidate_scale = candidate_scale
  )
  factor_plan$coefficient_structure <- "dense"
  baseline_states <- lapply(seq_len(nrow(scales)), function(draw) {
    list(list(coefficient_factor = diag(scales[draw, ])))
  })
  observed <- .marglik_covariance_plan_factor_grid_loglik(
    cache = new.env(parent = emptyenv()),
    y = outcome,
    means = means,
    sampling_covariance = sampling,
    random_covariance_plans = list(time = factor_plan),
    random_covariance_states = baseline_states,
    block_indices = list(1:4),
    extra_variances = matrix(0, nrow(means), ncol(means)),
    update_grid = factor_update
  )
  for (grid in seq_len(nrow(candidate_scale))) {
    for (draw in seq_len(nrow(means))) {
      scale <- scales[draw, ]
      scale[[2L]] <- candidate_scale[grid, draw]
      root <- matrix(correlation_draws[draw, , ], 4L, 4L)
      covariance <- sampling + tcrossprod(root * scale)
      expected[grid, draw] <- .factor_grid_dense_loglik(
        outcome,
        means[draw, ],
        covariance
      )
    }
  }
  expect_equal(observed, expected, tolerance = 1e-12)
})
