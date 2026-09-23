test_that("normal cluster quadrature fallback resolves its omitted order", {

  setup <- list(S = 2L, cluster = list(1:2), mu = matrix(0, 2L, 2L),
    tau_within = matrix(0.2, 2L, 2L), tau_between = matrix(0, 2L, 2L),
    weights = c(2, 0.5))
  local_mocked_bindings(.has_native_norm_cluster_quadrature = function(...) FALSE,
                        .package = "RoBMA")
  actual <- .log_lik_cluster_norm_quadrature(setup, c(0.1, -0.3), c(0.3, 0.4), FALSE)
  expected <- sum(dnorm(c(0.1, -0.3), 0, sqrt(c(0.3, 0.4)^2 + 0.2^2),
                        log = TRUE) * c(2, 0.5))
  expect_equal(actual, matrix(expected, 2L, 1L), tolerance = 1e-12)
})

test_that("selection deletion rejects a nonfinite full event density", {

  local_mocked_bindings(
    .selection_retains_sampling = function(...) FALSE,
    .data_selection_execution_plan = function(...) list(row_blocks = list(1:2)),
    .estimate_normal_covariance_target_location_from_setup = function(...) {
      list(means = matrix(0, 1L, 2L), y = c(0, 0))
    },
    .selection_joint_signed_context = function(...) list(),
    .selection_joint_random_factor_samples = function(...) list(),
    .selection_joint_covariance_lower = function(...) matrix(c(1, 0, 1), 1L),
    .selection_joint_event_numerator = function(...) list(log_numerator = -Inf),
    .package = "RoBMA"
  )
  local_mocked_bindings(selection_context_subset_observations = function(...) list(),
                        .package = "BayesTools")
  expect_error(.selection_joint_deletion_loglik_from_setup(
    list(S = 1L, data = list(), selection_sei = c(1, 1)), list(1L)),
    "observed outcomes have zero or invalid event density.", fixed = TRUE)
})

test_that("deletion fallback preserves coordinates resolved by the fast path", {

  local_mocked_bindings(
    .data_selection_execution_plan = function(...) list(row_blocks = list(1L, 2L),
      relative_tolerance = 1e-6, max_points_per_scramble = 128L),
    .data_selection_model = function(...) list(groups = list(row_blocks = list(1:2))),
    .selection_joint_signed_context = function(...) list(),
    .selection_conditioned_sampling_independent_targets = function(...) {
      list(mean = matrix(c(5, NA_real_), 1L))
    },
    .selection_conditioned_sampling_state = function(...) list(
      baseline_mu = matrix(0, 1L, 2L), e = matrix(0, 1L, 2L),
      sampling_covariance = diag(2), integrated_covariance = array(0, c(1L, 2L, 2L)),
      total_covariance = array(diag(c(0, 1)), c(1L, 2L, 2L))),
    .selection_deleted_gaussian_coordinates = function(...) {
      list(mean = c(0, 0), variance = c(0, 1))
    }, .package = "RoBMA"
  )
  out <- .selection_conditioned_sampling_estimate_targets(
    list(S = 1L, K = 2L, yi = c(0, 0), mu = matrix(c(2, 3), 1L)), "mean")
  expect_identical(out$mean, matrix(c(5, 3), 1L))
})
