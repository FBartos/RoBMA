test_that("mixed-state uncertainty accepts the plan's logical active mass", {

  rows <- seq_len(3000L)
  for (count in c(5L, 7L, 10L, 23L)) {
    active <- seq_len(count)
    mass <- mean(rows %in% active)
    contributions <- matrix(mass, nrow = 1L, ncol = count)
    attr(contributions, "contribution_rows") <- active
    attr(contributions, "expected_chain_ids") <- 1L
    attr(contributions, "target") <- mass
    result <- .iwmde_mcmc_contributions(
      contributions, active, rows, rep(1L, length(rows)), mass
    )
    expect_equal(as.numeric(result), as.numeric(rows %in% active))
    error <- .iwmde_active_mass_mcse(
      active, rows, rep(1L, length(rows)), 1L, mass
    )
    expect_true(is.finite(error[["mcse"]]))
  }
  expect_error(.iwmde_active_mass_mcse(
    1:5, rows, rep(1L, length(rows)), 1L, .5
  ), "Inconsistent active mass for IWMDE contributions.", fixed = TRUE)
})

test_that("two qCMDE grids retain their nonzero ordinate disagreement", {

  log_q <- matrix(0, nrow = 1L, ncol = 4L)
  normalizers <- list(rep(0, 4L), rep(log(2), 4L))
  selected <- .iwmde_qcmde_select_refinement(log_q, normalizers, 1, 4L)
  final <- .iwmde_qcmde_density_from_normalizer(
    log_q, normalizers[[selected[["final_index"]]]], 1, 4L
  )
  validation <- .iwmde_qcmde_density_from_normalizer(
    log_q, normalizers[[selected[["validation_index"]]]], 1, 4L
  )
  expect_equal(final, .5)
  expect_equal(validation, 1)
  change <- .iwmde_qcmde_ordinate_change(final, validation)
  expect_gt(change[["relative"]], .1)
})

test_that("ordinate-only qCMDE does not infer bulk mass from requested points", {

  testthat::local_mocked_bindings(
    .selection_normalizer_grid = function(...) NULL,
    .selection_covariance_grid = function(...) NULL,
    .iwmde_log_q_grid = function(context, parameter, values, ...) {
      matrix(stats::dnorm(values, log = TRUE), ncol = 1L)
    },
    .iwmde_qcmde_refinement_pair_converged = function(...) FALSE,
    .iwmde_qcmde_pilot_bulk_ess = function(...) stop("not a density grid"),
    .package = "RoBMA"
  )
  grid <- list(x = c(-1, 0, 1), z = c(-1, 0, 1),
               log_jacobian = rep(0, 3L), all_index = 1:3)
  result <- .iwmde_qcmde_evaluate_grid_sequence(
    context = list(), parameter = "mu", display_grid = c(0, 2, 50),
    normalizer_plan = list(grid_sequence = rep(list(grid), 3L), all_grid = grid),
    row_states = list(list(row_index = 1L)), replacement = list(),
    estimator_rows = 1L, active_mass = 1, denominator = 1,
    density_output = FALSE
  )
  expect_false(result[["pilot_gate_stopped"]])
  expect_length(result[["log_normalizer_sequence"]], 3L)
})

test_that("density mass remedies use reported row sampling uncertainty", {

  diagnostic <- list(max_sampling_relative_mcse = .2,
                     target_relative_mcse = .05, all_rows_used = FALSE)
  expect_identical(.iwmde_diagnostics_mass_failure_action(diagnostic, "iwmde"),
    paste0("Try increasing 'samples' in the 'density_control' argument or using ",
           "'samples = Inf' for the eligible-row census"))
  diagnostic[["all_rows_used"]] <- TRUE
  expect_match(.iwmde_diagnostics_mass_failure_action(diagnostic, "iwmde"),
               "normalization_points", fixed = TRUE)
})

test_that("scalar, indexed, and fallback replacements retain covariance plans", {

  update <- structure(list(family = "affine", coefficient_input = "source"),
    class = "BayesTools_random_effects_marginal_update_plan")
  testthat::local_mocked_bindings(
    .iwmde_predictor_formula_parameter = function(...) NULL,
    .package = "RoBMA"
  )
  expected <- c(mu = "scalar", `gamma[1]` = "indexed", custom = "fallback")
  for (parameter in names(expected)) {
    replacement <- .iwmde_replacement_spec(
      list(), parameter, list(type = "primitive", covariance_update = update)
    )
    expect_identical(replacement[["type"]], expected[[parameter]])
    expect_identical(replacement[["covariance_update"]], update)
  }
})

test_that("equal selection weights never reuse a varying Gaussian likelihood", {

  base <- list(covariance = diag(2), coefficient = 0, direction = c(1, 1))
  for (family in c("mean", "variance")) {
    shared <- list(family = family, loading = matrix(1, 2, 1))
    expect_null(.selection_covariance_grid_binding(shared, 1:2, base, c(1, 1)))
  }
  # With the selection weight cancelled, the remaining exact likelihood varies.
  expect_false(isTRUE(all.equal(sum(dnorm(c(0, 0), 0, 1, log = TRUE)),
                               sum(dnorm(c(0, 0), 1, 1, log = TRUE)))))
})

test_that("factor covariance grids require the replacement's coefficient scale", {

  testthat::local_mocked_bindings(
    .iwmde_uses_known_v_random_marginal_likelihood = function(...) TRUE,
    .is_data_weights = function(...) FALSE,
    .iwmde_known_v_random_marginal_setup = function(...) stop("wrong scale routed"),
    .package = "RoBMA"
  )
  for (type in c("scalar", "random_component_sd")) {
    update <- structure(list(family = "factor", coefficient_input =
      if (type == "scalar") "quantity" else "source"),
      class = "BayesTools_random_effects_marginal_update_plan")
    expect_null(.iwmde_log_q_grid_known_v_random_factor(
      list(), "tau", c(.1, .2), list(),
      list(type = type, covariance_update = update), list()
    ))
  }
})

test_that("unsupported factor update grids decline to the generic evaluator", {

  testthat::local_mocked_bindings(
    .iwmde_uses_known_v_random_marginal_likelihood = function(...) TRUE,
    .is_data_weights = function(...) FALSE,
    .iwmde_known_v_random_marginal_setup = function(...) list(),
    .iwmde_predictor_evaluate_fixed_mu = function(...) matrix(0, 1L, 2L),
    .data_effect_direction = function(...) "positive",
    .brma_mv_random_effects_marginal_inputs = function(...) list(),
    .brma_mv_random_effects_marginal_factor_states = function(...) list(),
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    random_effects_marginal_update_grid = function(...) stop("source unavailable"),
    .package = "BayesTools"
  )
  update <- structure(list(family = "factor", coefficient_input = "source"),
    class = "BayesTools_random_effects_marginal_update_plan")
  expect_null(.iwmde_log_q_grid_known_v_random_factor(
    list(posterior_samples = matrix(1, 1L, 1L), data = list(outcome = list(yi = c(0, 1)))),
    "tau", c(.1, .2), list(list(row_index = 1L)),
    list(type = "scalar", covariance_update = update), list()
  ))
})

test_that("asymmetric affine anchors decline before the native symmetric kernel", {

  testthat::local_mocked_bindings(
    .iwmde_uses_known_v_random_marginal_likelihood = function(...) TRUE,
    .is_data_weights = function(...) FALSE,
    .iwmde_known_v_random_marginal_setup = function(...) list(sampling_covariance = diag(2)),
    .iwmde_predictor_evaluate_fixed_mu = function(...) matrix(0, 1L, 2L),
    .iwmde_random_affine_coefficients = function(values, ...) values,
    .data_effect_direction = function(...) "positive",
    .iwmde_build_replacement_samples = function(...) list(valid = TRUE, samples = matrix(1)),
    .brma_mv_random_effects_marginal_vcov = function(...) {
      list(samples = array(c(1, .1, .2, 1), c(1L, 2L, 2L)))
    },
    .iwmde_random_affine_log_likelihood_chunk = function(...) stop("asymmetric kernel reached"),
    .package = "RoBMA"
  )
  update <- structure(list(family = "affine", coefficient_input = "source"),
    class = "BayesTools_random_effects_marginal_update_plan")
  expect_null(.iwmde_log_q_grid_known_v_random_affine(
    list(posterior_samples = matrix(1, 1L, 1L), data = list(outcome = list(yi = c(0, 1)))),
    "tau", c(.1, .2), list(list(row_index = 1L)),
    list(type = "scalar", covariance_update = update), list()
  ))
})

test_that("density windows use the continuous branch of spike-heavy posteriors", {

  values <- c(rep(0, 6000L), seq(10, 20, length.out = 30L))
  active <- values > 0
  testthat::local_mocked_bindings(
    .iwmde_context_unavailable_reason = function(...) NULL,
    .iwmde_known_v_tau_zero_boundary_reason = function(...) NULL,
    .iwmde_parameter_values = function(...) values,
    .iwmde_parameter_condition_rows = function(...) rep(TRUE, length(values)),
    .iwmde_parameter_components = function(...) list(active = active,
      point_masses = data.frame(x = 0, mass = mean(!active))),
    .iwmde_parameter_support = function(...) c(-Inf, Inf),
    .iwmde_plan_baseline_contract = function(context, plan, candidate_rows, candidate_values) {
      list(estimator_rows = candidate_rows, estimator_values = candidate_values)
    },
    .iwmde_predictor_formula_parameter = function(...) NULL,
    .package = "RoBMA"
  )
  plan <- .iwmde_plan_prepare_contract(
    list(posterior_samples = matrix(values, ncol = 1L), chain_id = rep(1L, length(values))),
    list(target = list(parameter = "mu"), parameter_spec = list(type = "primitive"),
      row_budget = Inf, method = "iwmde",
      outputs = list(need_density = TRUE, requested_values = numeric()),
      control = list(n_points = 21L, display_grid = "uniform",
                     normalization_points = 21L, normalization_prob = .99))
  )
  expect_identical(plan[["status"]], "ok")
  expect_equal(plan[["support"]][["xlim"]], .iwmde_plot_range(values[active], c(-Inf, Inf)))
  expect_gt(plan[["support"]][["xlim"]][1L], 0)
  expect_gt(plan[["support"]][["xlim"]][2L], max(values))
})

test_that("histograms reject degenerate windows and retain narrow finite windows", {

  expect_error(.iwmde_histogram(1, c(1, 1)),
    "IWMDE histogram range must contain two increasing finite values.", fixed = TRUE)
  values <- c(1, 1 + .Machine$double.eps)
  histogram <- .iwmde_histogram(values, range(values))
  expect_true(all(is.finite(histogram[["density"]])))
  expect_equal(sum(histogram[["density"]] * diff(histogram[["breaks"]])), 1)
})

test_that("computed ordinates with coincident displayed diagnostics retain rows", {

  entry <- function(height) BayesTools::posterior_ordinate_attribute(
    value = 0, ordinate = height, method = "iwmde", density_method = "IWMDE",
    diagnostics = list(estimator = "iwmde", normalization_relative_error = 0),
    parameter = "mu"
  )
  posterior <- list(structure(1:3, posterior_ordinate = entry(1)),
                    structure(2:4, posterior_ordinate = entry(2)))
  diagnostics <- .iwmde_collect_public_density_diagnostics(posterior)
  expect_equal(nrow(diagnostics), 2L)
  expect_identical(diagnostics[["parameter"]], c("mu", "mu"))
})
