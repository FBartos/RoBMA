test_that("gate-only allocations remain independently addressable", {

  first <- list(label = "first", weight_name = NULL, n_targets = 2L)
  second <- list(label = "second", weight_name = NULL, n_targets = 3L)
  design <- list(mu = list(parameter = "mu", random_allocations = list(first, second)))
  expect_identical(.brma_random_parameter_design_allocation(
    design, list(allocation = "second", formula_parameter = "mu")
  ), second)
})

test_that("allocation prior gates reject missing and fractional component indices", {

  for (index in list(NA_integer_, NULL, 0L, 3L, 1.5)) {
    metadata <- list(quantity = "var_prop", index = index,
                     component_indicators = c(NA_character_, NA_character_),
                     parent_indicators = character())
    expect_error(.brma_random_parameter_allocation_gate_prior(list(), metadata),
      "Variance-proportion gate metadata have no valid component index.", fixed = TRUE)
  }
})

test_that("incomplete allocation scale metadata has no boundary alternative", {

  term <- list(sd_binding = list(true_allocation = TRUE,
    allocations = list(list(scale = NULL, target = "block"))))
  testthat::local_mocked_bindings(
    .brma_random_parameter_design_term = function(...) term,
    .package = "RoBMA"
  )
  selected <- list(spec = list(quantity = "sd", allocation_derived = TRUE))
  expect_null(.brma_random_parameter_zero_boundary_alternative(list(), selected))
})

test_that("gate-only proportions decline unsupported density coordinates", {

  selected <- list(spec = list(quantity = "var_prop", source_type = "composite",
    source_parameter = "", source_transform = "var_prop", evaluator = "allocation",
    allocation_index = 1L, label = "tau2_prop(study)"),
    allocation_definition = list(scale = "total_variance", n_targets = 1L,
      parent_factors = list(list(inclusion_name = "gate"))))
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .get_posterior_samples = function(...) cbind(gate = c(0, 1)),
    .brma_random_parameter_simplex_exclusions = function(...) character(),
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    random_effects_marginal_update_plan = function(...) NULL,
    .package = "BayesTools"
  )
  target <- .brma_random_parameter_density_target(list(fit = list()), "tau2_prop(study)")
  expect_match(target[["reason"]], "no supported scalar random-component coordinate",
               fixed = TRUE)
})
