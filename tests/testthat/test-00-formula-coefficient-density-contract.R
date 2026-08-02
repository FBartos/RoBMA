context("BayesTools formula coefficient density contract")


test_that("transformed coefficient hypotheses use exact structural weights", {

  transform <- list(
    schema_version             = 1L,
    formula_design_version     = 3L,
    parameter_registry_version = 3L,
    parameter                  = "mu",
    target_scale               = "original",
    source_names               = c("mu_intercept", "mu_x"),
    target_names               = c("mu_intercept", "mu_x"),
    matrix = rbind(
      mu_intercept = c(mu_intercept = 1, mu_x = -2),
      mu_x         = c(mu_intercept = 0, mu_x = 0.5)
    ),
    source_transforms = c(mu_intercept = "identity", mu_x = "identity"),
    output_transforms = c(mu_intercept = "identity", mu_x = "identity"),
    dependencies = data.frame(
      target      = c("mu_intercept", "mu_intercept", "mu_x"),
      source      = c("mu_intercept", "mu_x", "mu_x"),
      coefficient = c(1, -2, 0.5),
      stringsAsFactors = FALSE
    ),
    sources = data.frame(),
    targets = data.frame(
      target            = c("mu_intercept", "mu_x"),
      structural_status = "dependent",
      fixed_value       = NA_real_,
      reason            = "",
      stringsAsFactors  = FALSE
    )
  )
  class(transform) <- c("BayesTools_formula_coefficient_transform", "list")
  density <- structure(list(sentinel = TRUE), class = "prior_linear_density")
  fit <- structure(list(sentinel = TRUE), class = "BayesTools_fit")
  object <- list(fit = fit)
  selected <- list(
    parameter = "mu_intercept",
    component = "mods",
    entry     = list(formula_parameter = "mu")
  )
  context <- structure(list(sentinel = TRUE), class = "prior_density_context")
  samples <- structure(list(), prior_density_context = context)
  observed_context <- NULL
  observed_values  <- numeric()
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    JAGS_formula_prior_density = function(..., context) {
      observed_context <<- context
      density
    },
    prior_density_ordinate = function(density, value) {
      observed_values <<- c(observed_values, value)
      list(behavior = "regular", exact = TRUE)
    },
    .package = "BayesTools"
  )

  target <- .hypothesis_brma_formula_coefficient_target(
    object                    = object,
    selected                  = selected,
    standardized_coefficients = FALSE
  )
  target <- .hypothesis_brma_formula_prior_target(
    object      = object,
    samples     = samples,
    hypothesis  = BayesTools::hypothesis_parse("mu_intercept = 3"),
    target_info = target
  )

  expect_identical(observed_context, context)
  expect_identical(observed_values, 3)
  expect_identical(target[["prior_density"]], density)
  expect_identical(target[["parameter_spec"]][["type"]], "linear")
  expect_identical(
    target[["parameter_spec"]][["weights"]],
    c(mu_intercept = 1, mu_x = -2)
  )
})


test_that("nonlinear joint coefficient transforms fail qCMDE/IWMDE closed", {

  target <- list(
    formula_parameter = "log_tau",
    target            = "log_tau_intercept",
    target_i          = 1L,
    transform = list(
      matrix = matrix(
        c(1, -2),
        nrow = 1L,
        dimnames = list(
          "log_tau_intercept",
          c("log_tau_intercept", "log_tau_x")
        )
      ),
      source_transforms = c(
        log_tau_intercept = "log",
        log_tau_x         = "identity"
      ),
      output_transforms = c(log_tau_intercept = "exp")
    )
  )
  density <- structure(list(sentinel = TRUE), class = "prior_linear_density")
  testthat::local_mocked_bindings(
    JAGS_formula_prior_density = function(...) density,
    prior_density_ordinate = function(...) list(exact = TRUE),
    .package = "BayesTools"
  )

  observed <- .hypothesis_brma_formula_prior_target(
    object      = list(fit = structure(list(), class = "BayesTools_fit")),
    samples     = list(),
    hypothesis  = BayesTools::hypothesis_parse("log_tau_intercept = 1"),
    target_info = target
  )

  expect_identical(
    observed[["parameter_spec"]][["type"]],
    "unsupported_formula_transform"
  )
  expect_match(observed[["parameter_spec"]][["reason"]], "nonlinear joint")
})
