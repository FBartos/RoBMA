test_that("affine spectral likelihood matches direct dense evaluation", {

  base <- matrix(c(
    2.0, 0.4, 0.3,
    0.4, 1.5, 0.2,
    0.3, 0.2, 1.2
  ), nrow = 3L, byrow = TRUE)
  update <- matrix(c(
    0.7,  0.2, -0.1,
    0.2,  0.4,  0.1,
   -0.1,  0.1,  0.6
  ), nrow = 3L, byrow = TRUE)
  reference <- 0.35
  covariance_reference <- base + reference * update
  residual <- c(-0.4, 0.8, 0.2)
  coefficients <- c(0, 0.2, 0.9, 1.4)

  plan <- .known_v_affine_spectral_plan(
    base_covariance      = covariance_reference,
    update_covariance    = update,
    reference_coefficient = reference
  )
  observed <- .known_v_affine_log_likelihood(
    plan         = plan,
    residual     = residual,
    coefficients = coefficients
  )
  expected <- vapply(coefficients, function(coefficient) {
    covariance <- base + coefficient * update
    root <- chol(covariance)
    solved <- backsolve(root, forwardsolve(t(root), residual))
    -0.5 * (
      length(residual) * log(2 * pi) +
        2 * sum(log(diag(root))) +
        sum(residual * solved)
    )
  }, numeric(1))

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("affine spectral GLS matches direct dense projection", {

  base <- matrix(c(
    1.8, 0.5, 0.2, 0.1,
    0.5, 1.7, 0.4, 0.2,
    0.2, 0.4, 1.5, 0.3,
    0.1, 0.2, 0.3, 1.3
  ), nrow = 4L, byrow = TRUE)
  update <- tcrossprod(c(0.6, -0.3, 0.4, 0.8)) + diag(c(0.1, 0.2, 0.3, 0.4))
  coefficient <- 0.65
  X <- cbind(1, c(-1, -0.2, 0.5, 1.4))
  y <- c(0.2, -0.4, 0.9, 1.1)

  plan <- .known_v_affine_spectral_plan(
    base_covariance   = base,
    update_covariance = update
  )
  observed <- .known_v_affine_gls_projection(
    plan          = plan,
    X             = X,
    y             = y,
    coefficient   = coefficient,
    return_full_H = TRUE
  )

  covariance <- base + coefficient * update
  W <- solve(covariance)
  XtWX_inv <- solve(crossprod(X, W %*% X))
  beta_hat <- as.vector(XtWX_inv %*% crossprod(X, W %*% y))
  H <- X %*% XtWX_inv %*% crossprod(X, W)
  residual_variance <- diag(covariance - X %*% XtWX_inv %*% t(X))

  expect_equal(observed$WX, W %*% X, tolerance = 1e-12)
  expect_equal(observed$XtWX_inv, XtWX_inv, tolerance = 1e-12)
  expect_equal(observed$beta_hat, beta_hat, tolerance = 1e-12)
  expect_equal(
    observed$residual,
    as.vector(y - X %*% beta_hat),
    tolerance = 1e-12
  )
  expect_equal(observed$H, H, tolerance = 1e-12)
  expect_equal(observed$H_diag, diag(H), tolerance = 1e-12)
  expect_equal(
    observed$residual_variance,
    residual_variance,
    tolerance = 1e-12
  )
})


test_that("affine spectral engine reports invalid covariance candidates", {

  plan <- .known_v_affine_spectral_plan(
    base_covariance   = diag(2),
    update_covariance = -diag(2)
  )

  expect_null(.known_v_affine_log_likelihood(
    plan         = plan,
    residual     = c(0, 0),
    coefficients = c(0, 1)
  ))
  expect_null(.known_v_affine_gls_projection(
    plan        = plan,
    X           = matrix(1, nrow = 2L),
    y           = c(0, 0),
    coefficient = 1
  ))
})


test_that("draw-specific affine covariance chunks match dense block likelihoods", {

  reference <- 0.3
  coefficients <- c(0, 0.2, 0.8)
  base <- array(0, dim = c(2L, 4L, 4L))
  update <- array(0, dim = c(2L, 4L, 4L))
  base[1L, , ] <- diag(c(1.2, 1.4, 1.1, 1.3)) +
    reference * diag(c(0.3, 0.2, 0.4, 0.1))
  base[2L, , ] <- diag(c(1.5, 1.1, 1.6, 1.2)) +
    reference * diag(c(0.2, 0.5, 0.1, 0.3))
  update[1L, , ] <- diag(c(0.3, 0.2, 0.4, 0.1))
  update[2L, , ] <- diag(c(0.2, 0.5, 0.1, 0.3))
  update[, 1L, 2L] <- update[, 2L, 1L] <- c(0.05, -0.04)
  update[, 3L, 4L] <- update[, 4L, 3L] <- c(-0.03, 0.02)
  means <- rbind(c(0.1, -0.2, 0.4, 0.3), c(-0.1, 0.3, 0.2, -0.2))
  outcome <- c(0.2, -0.5, 0.7, 0.1)
  blocks <- list(1:2, 3:4)

  observed <- .iwmde_random_affine_log_likelihood_chunk(
    base_covariance       = base,
    update_covariance     = update,
    reference_coefficient = reference,
    coefficients          = coefficients,
    means                 = means,
    outcome               = outcome,
    blocks                = blocks
  )
  expected <- matrix(0, nrow = length(coefficients), ncol = 2L)
  for (s in seq_len(2L)) {
    for (g in seq_along(coefficients)) {
      covariance <- base[s, , ] +
        (coefficients[[g]] - reference) * update[s, , ]
      for (block in blocks) {
        block_covariance <- covariance[block, block, drop = FALSE]
        residual <- outcome[block] - means[s, block]
        root <- chol(block_covariance)
        solved <- backsolve(root, forwardsolve(t(root), residual))
        expected[g, s] <- expected[g, s] - 0.5 * (
          length(block) * log(2 * pi) +
            2 * sum(log(diag(root))) +
            sum(residual * solved)
        )
      }
    }
  }

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("single unit random component receives the exact IID grid plan", {

  update <- structure(
    list(
      family = "affine",
      update = "scale",
      source_parameter = "tau",
      coefficient_transform = list(type = "square"),
      coefficient_input = "source",
      invariant_covariance = list(
        reference_coefficient = 0,
        base_covariance = matrix(0, 3L, 3L),
        update_covariance = diag(3L)
      )
    ),
    class = c(
      "BayesTools_random_effects_marginal_update_plan",
      "list"
    )
  )
  term <- list(
    n_columns = 1L,
    model_matrix = matrix(1, nrow = 3L),
    group_map = 1:3,
    sd_parameter_names = "tau",
    group_covariance = NULL
  )
  replacement <- list(
    type = "scalar",
    covariance_update = update
  )

  plan <- .iwmde_known_v_random_single_iid_plan(
    term = term,
    parameter = "tau",
    replacement = replacement,
    K = 3L
  )
  expect_identical(plan$target_mode, "direct_sd")
  expect_true(plan$unique_groups)
  expect_identical(plan$coefficient_transform, list(type = "square"))

  term$model_matrix[2L, 1L] <- 2
  expect_null(.iwmde_known_v_random_single_iid_plan(
    term = term,
    parameter = "tau",
    replacement = replacement,
    K = 3L
  ))
})


test_that("direct grouped-scale grid matches dense Gaussian likelihoods", {

  yi <- c(0.1, -0.2, 0.3, 0.05)
  vi <- c(0.04, 0.05, 0.06, 0.07)
  group_map <- c(1L, 1L, 2L, 2L)
  posterior_samples <- matrix(
    c(0.2, 0.45),
    ncol = 1L,
    dimnames = list(NULL, "tau")
  )
  mu_samples <- rbind(
    c(0.02, 0.02, -0.01, -0.01),
    c(-0.03, 0.01, 0.04, 0.02)
  )
  values <- c(0, 0.4, 0.9)
  row_states <- list(list(row_index = 1L), list(row_index = 2L))
  context <- list(
    posterior_samples = posterior_samples,
    data = list(outcome = data.frame(yi = yi))
  )
  testthat::local_mocked_bindings(
    .iwmde_known_v_random_group_iid_plan = function(...) {
      list(
        source_parameter      = "tau",
        coefficient_transform = list(type = "square"),
        target_mode           = "direct_sd",
        group_maps            = list(group_map),
        unique_groups         = FALSE
      )
    },
    .iwmde_known_v_random_marginal_setup = function(...) {
      list(
        sampling_covariance = diag(vi),
        dependency_blocks   = unname(split(seq_along(group_map), group_map)),
        group_iid_plan_cache = new.env(parent = emptyenv())
      )
    },
    .iwmde_predictor_evaluate_fixed_mu = function(...) mu_samples,
    .iwmde_predictor_log_prior = function(...) {
      numeric(length(values) * length(row_states))
    },
    .data_effect_direction = function(...) "positive",
    .package = "RoBMA"
  )

  observed <- .iwmde_log_q_grid_known_v_random_group_iid(
    context      = context,
    parameter    = "tau",
    values       = values,
    row_states   = row_states,
    replacement  = list(),
    active_setup = list()
  )
  expected <- matrix(NA_real_, nrow = length(values), ncol = nrow(mu_samples))
  same_group <- outer(group_map, group_map, "==") * 1
  for (value_i in seq_along(values)) {
    covariance <- diag(vi) + values[[value_i]]^2 * same_group
    for (draw in seq_len(nrow(mu_samples))) {
      expected[value_i, draw] <- .marglik_mvn_log_density(
        y          = yi,
        mean       = mu_samples[draw, ],
        covariance = covariance
      )
    }
  }

  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("invariant affine covariance uses the declared basis directly", {

  sampling <- matrix(c(
    0.5, 0.1, 0.0,
    0.1, 0.6, 0.2,
    0.0, 0.2, 0.7
  ), nrow = 3L, byrow = TRUE)
  basis <- tcrossprod(c(1, 0.5, -0.25))
  update <- list(invariant_covariance = list(
    reference_coefficient = 0,
    base_covariance       = matrix(0, 3L, 3L),
    update_covariance     = basis
  ))
  invariant <- .iwmde_random_affine_invariant_covariance(
    update = update,
    setup  = list(sampling_covariance = sampling),
    K      = 3L
  )
  coefficients <- c(0, 0.2, 0.8)
  means <- rbind(c(0.1, -0.1, 0.2), c(-0.2, 0.3, 0.1))
  outcome <- c(0.3, -0.4, 0.6)
  observed <- .iwmde_random_affine_log_likelihood_invariant(
    base_covariance       = invariant$base_covariance,
    update_covariance     = invariant$update_covariance,
    reference_coefficient = invariant$reference_coefficient,
    coefficients          = coefficients,
    means                 = means,
    outcome               = outcome,
    blocks                = list(1:3)
  )
  expected <- matrix(NA_real_, nrow = length(coefficients), ncol = nrow(means))
  for (g in seq_along(coefficients)) {
    covariance <- sampling + coefficients[[g]] * basis
    for (s in seq_len(nrow(means))) {
      expected[g, s] <- .marglik_mvn_log_density(
        y          = outcome,
        mean       = means[s, ],
        covariance = covariance
      )
    }
  }

  expect_equal(invariant$base_covariance, sampling)
  expect_identical(invariant$update_covariance, basis)
  expect_equal(observed, expected, tolerance = 1e-12)
})


test_that("generalized affine routing requires sampling dependence", {

  expect_false(.iwmde_random_affine_has_sampling_dependence(list(
    sampling_covariance = diag(c(0.1, 0.2, 0.3))
  )))

  covariance <- diag(c(0.1, 0.2, 0.3))
  covariance[1L, 2L] <- covariance[2L, 1L] <- 0.04
  expect_true(.iwmde_random_affine_has_sampling_dependence(list(
    sampling_covariance = covariance
  )))
})
