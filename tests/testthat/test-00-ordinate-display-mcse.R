test_that("display Jacobians scale every absolute mixture density MCSE component", {

  # The partial mixed-chain aggregation fixture from the estimator tests:
  # selected active rows have a conservative active-branch plus mass bound.
  withr::local_seed(20260806)
  conditioned_rows <- seq_len(4000L)
  chain_id <- rep(1:4, each = 1000L)
  active_rows <- unlist(lapply(0:3, function(chain) {
    chain * 1000L + sort(sample.int(1000L, 50L))
  }))
  selected <- active_rows[c(seq_len(25L), 50L + seq_len(25L),
    100L + seq_len(25L), 150L + seq_len(25L))]
  density <- .iwmde_density_aggregate(
    log_terms = matrix(log(rep(c(.5, 1, 1.5, 2), 25L)), nrow = 1L),
    active_mass = .05, denominator = 100L,
    contribution_rows = selected, sampling_population_rows = active_rows,
    chain_id = chain_id[selected], expected_chain_ids = 1:4,
    conditioned_rows = conditioned_rows, conditioned_chain_id = chain_id)
  mcse <- .iwmde_mixture_mcse(density$contributions, density$mcmc_contributions,
    density$active_mass_error, active_mass = .05)
  density[names(mcse)] <- mcse
  density$x <- .5
  density$evaluation_x <- .5
  expect_gt(density$active_branch_mcse, 0)
  expect_gt(density$active_mass_component_mcse, 0)
  expect_equal(density$active_mass_component_mcse,
    abs(density$y / .05) * density$active_mass_mcse, tolerance = 1e-14)
  expect_equal(density$mcse,
    density$active_branch_mcse + density$active_mass_component_mcse, tolerance = 1e-14)
  diagnostics <- .iwmde_density_bf_diagnostics(density, .5)
  diagnostics$estimator <- "q_grid_cmde"
  diagnostics$active_mass <- .05
  for (field in c("active_mass_mcse", "active_mass_relative_mcse", "mixture_mcse_type",
    "sampling_fraction", "sampling_uncertainty_type")) {
    diagnostics[[field]] <- density[[field]]
  }
  diagnostics$mcmc_uncertainty_scope <- mcse$uncertainty_scope
  diagnostics$mcmc_uncertainty_status <- mcse$uncertainty_status
  source <- .iwmde_posterior_ordinate_attribute(
    list(status = "ok", diagnostics = diagnostics), "qCMDE",
    metadata = list(parameter = "mu", conditional = FALSE, point_mass_total = .95),
    allow_rejected = TRUE)
  expect_s3_class(source, "BayesTools_posterior_ordinate")
  transformed <- .hypothesis_brma_transform_iwmde_ordinate(source, list(type = "tanh"))
  jacobian <- 1 - tanh(.5)^2
  expect_equal(transformed$ordinate, density$y / jacobian, tolerance = 1e-14)
  for (field in c("mcse", "active_branch_mcse", "active_mass_component_mcse", "sampling_mcse")) {
    expect_equal(transformed$diagnostics[[field]], source$diagnostics[[field]] / jacobian,
      tolerance = 1e-14, info = field)
  }
  expect_equal(transformed$diagnostics$mcse,
    transformed$diagnostics$active_branch_mcse +
      transformed$diagnostics$active_mass_component_mcse, tolerance = 1e-14)
  for (field in c("active_mass", "active_mass_mcse", "active_mass_relative_mcse",
    "relative_mcse", "active_branch_relative_mcse", "sampling_relative_mcse", "sampling_fraction")) {
    expect_identical(transformed$diagnostics[[field]], source$diagnostics[[field]], info = field)
  }
  expect_identical(transformed[c("parameter", "conditional", "point_mass_total")],
    source[c("parameter", "conditional", "point_mass_total")])
  expect_equal(transformed$value, tanh(source$value))
  expect_equal((stats::dnorm(.5) / jacobian) / transformed$ordinate,
    stats::dnorm(.5) / source$ordinate, tolerance = 1e-14)
})
