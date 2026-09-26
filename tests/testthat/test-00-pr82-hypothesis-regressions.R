test_that("coefficient transform routes follow the declared BayesTools map type", {

  transform <- function(targets) {
    structure(list(target_scale = "original", target_names = "mu",
      matrix = matrix(c(1, 1), 1L, dimnames = list("mu", c("a", "b"))),
      source_transforms = c(a = "identity", b = "log"),
      output_transforms = c(mu = "identity"), targets = targets),
      class = "BayesTools_formula_coefficient_transform")
  }
  route <- function(targets) {
    testthat::local_mocked_bindings(
      JAGS_formula_coefficient_transform = function(...) {
        out <- transform(targets)
        out[["schema_version"]] <- 2L
        out
      },
      .package = "BayesTools"
    )
    .brma_formula_coefficient_route(
      object   = list(fit = structure(list(), class = "BayesTools_fit")),
      selected = list(parameter = "mu", component = "mods",
                      entry = list(formula_parameter = "mu"))
    )
  }

  unsupported <- route(formula_transform_targets(
    c(mu = "unsupported"), c(mu = "identity")
  ))
  expect_identical(unsupported[["type"]], "unsupported")
  expect_match(unsupported[["reason"]], "not supported by hypothesis()", fixed = TRUE)

  affine <- route(formula_transform_targets(c(mu = "affine"), c(mu = "identity")))
  expect_identical(affine[["type"]], "affine")
  expect_identical(affine[["weights"]], c(a = 1, b = 1))
  expect_identical(affine[["support"]], c(-Inf, Inf))

  # Transforms without the declared map types are not certified.
  for (targets in list(NULL, data.frame(target = "mu"))) {
    uncertified <- route(targets)
    expect_identical(uncertified[["type"]], "unsupported")
    expect_match(uncertified[["reason"]], "lacks the certified structural metadata",
                 fixed = TRUE)
  }
})

test_that("cross-level hypotheses keep all rows and build one shared context", {

  calls    <- 0L
  attached <- NULL
  testthat::local_mocked_bindings(
    hypothesis_BF = function(...) data.frame(BF = c(2, 3)),
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .iwmde_check_point_ordinate_supported = function(...) NULL,
    .iwmde_context = function(...) { calls <<- calls + 1L; list() },
    .hypothesis_brma_attach_iwmde_scalar = function(posterior, value, parameter_spec, ...) {
      attached <<- list(value = value, spec = parameter_spec)
      posterior
    },
    .hypothesis_brma_append_iwmde_warnings = function(table, ...) table,
    .package = "RoBMA"
  )
  spec <- list(type = "linear", weights = c(a = 1, b = -1), prior_density = NULL)
  plan <- list(
    parameter     = "mu_group",
    label         = "group",
    point         = TRUE,
    conditional   = FALSE,
    linear_target = list(
      posterior = 1:4, hypothesis = c("contrast = 0", "contrast = 1"),
      parameter = "contrast", weights = c(a = 1, b = -1)
    ),
    targets = list(
      list(value = 0, spec = spec),
      list(value = 1, spec = spec)
    )
  )
  result <- .hypothesis_plan_execute_combination(
    plans = list(plan), hypothesis = NULL, object = structure(list(fit = list(TRUE)), class = "brma"),
    logBF = FALSE, BF01 = FALSE, seed = NULL, density_method = "qCMDE",
    density_control = list(n_points = 50L), columns = NULL
  )
  expect_identical(result[["BF"]], c(2, 3))
  expect_identical(rownames(result), c("mu_group (1)", "mu_group (2)"))
  expect_identical(calls, 1L)
  # Both values share one estimate of the linear target with the plan's
  # fitted-coordinate weights.
  expect_identical(attached[["value"]], c(0, 1))
  expect_identical(attached[["spec"]], spec)
})

test_that("empty resolved hypothesis quantities fail with a metadata message", {

  testthat::local_mocked_bindings(
    hypothesis_resolve = function(...) list(occurrences = data.frame(quantity_id = character())),
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .brma_parameter_catalog_entries_for_quantities = function(...) data.frame(),
    .package = "RoBMA"
  )
  expect_error(.hypothesis_brma_select_parameter(
    list(), "mu = 0", "auto", list(catalog = NULL, entries = data.frame())
  ), paste0("Resolved hypothesis metadata are unavailable. Refit the model with ",
            "the current RoBMA/BayesTools build."), fixed = TRUE,
  class = "RoBMA_refit_required")
})

test_that("BayesTools resolution refusals are statement errors by their class", {

  metadata <- list(catalog = NULL, entries = data.frame())
  select   <- function(hypothesis) {
    tryCatch(
      .hypothesis_brma_select_parameter(list(), hypothesis, "auto", metadata),
      error = identity
    )
  }
  resolution <- function(class, ...) {
    structure(
      class = c(class, "BayesTools_parameter_resolution_error", "error",
                "condition"),
      list(message = "The statement was refused.", call = NULL, ...)
    )
  }
  # A statement without parameter symbols and a level of another component
  # are statement errors by their BayesTools class, with BayesTools' message
  # and fields, whatever the statement's form.
  refusals <- list(
    no_parameters      = resolution("BayesTools_hypothesis_no_parameters"),
    component_mismatch = resolution(
      "BayesTools_hypothesis_component_mismatch",
      symbol = "mu_g[a]", component = "mods"
    )
  )
  for (name in names(refusals)) {
    testthat::local_mocked_bindings(
      hypothesis_resolve = function(...) stop(refusals[[name]]),
      .package = "BayesTools"
    )
    expected <- refusals[[name]]
    class(expected) <- c("RoBMA_hypothesis_statement", class(expected))
    for (statement in c("1 > 0", "mu = 0")) {
      expect_identical(select(statement), expected,
                       info = paste(name, statement))
    }
  }
  # Other errors of the resolution are no statement refusals.
  unclassed <- simpleError("The resolution failed.")
  testthat::local_mocked_bindings(
    hypothesis_resolve = function(...) stop(unclassed),
    .package = "BayesTools"
  )
  expect_identical(select("1 > 0"), unclassed)
})

test_that("stale fitted metadata of hypothesis targets require a refit", {

  # Every hypothesis() stop on missing or unsupported fitted metadata has
  # the class "RoBMA_refit_required" with its parent
  # "BayesTools_refit_required", and no other class.
  refit <- c("RoBMA_refit_required", "BayesTools_refit_required",
             "error", "condition")
  selected <- list(
    parameter = "mu_x",
    component = "mods",
    entry     = list(formula_parameter = "mu", role = "fixed_coefficient")
  )
  transform <- structure(
    list(schema_version = 1L, target_names = "mu_x"),
    class = "BayesTools_formula_coefficient_transform"
  )
  testthat::local_mocked_bindings(
    JAGS_formula_coefficient_transform = function(...) transform,
    parameter_coordinates = function(...) data.frame(
      coordinate_name = "mu_z", formula_parameter = "mu", term = "z"
    ),
    .package = "BayesTools"
  )
  # A coefficient transform of an unsupported schema.
  error <- tryCatch(
    .hypothesis_brma_formula_coefficient_target(list(fit = list()), selected),
    error = identity
  )
  expect_identical(class(error), refit)
  expect_match(conditionMessage(error), "Refit the model", fixed = TRUE)
  # A resolved coefficient absent from the fitted transform.
  transform[["schema_version"]] <- 2L
  transform[["target_names"]]   <- "mu_z"
  error <- tryCatch(
    .hypothesis_brma_formula_coefficient_target(list(fit = list()), selected),
    error = identity
  )
  expect_identical(class(error), refit)
  # Weights on coordinates absent from the fitted coordinate table.
  error <- tryCatch(
    .hypothesis_brma_target_prior_parameters(list(fit = list()), c(mu_x = 1)),
    error = identity
  )
  expect_identical(class(error), refit)
  expect_match(conditionMessage(error), "Refit the model", fixed = TRUE)
  # A factor level without fitted linear weights is stale fitted metadata,
  # not a target refusal of a current fit.
  testthat::local_mocked_bindings(
    .iwmde_linear_weights = function(...) numeric(),
    .package = "RoBMA"
  )
  error <- tryCatch(
    .hypothesis_plan_level_target(
      plan   = list(draws = list(posterior = list(a = NULL)), label = "g"),
      object = NULL,
      level  = "a",
      value  = 0,
      label  = "g[a]"
    ),
    error = identity
  )
  expect_identical(class(error), refit)
  expect_identical(
    conditionMessage(error),
    paste0(
      "The linear weights of factor level 'g[a]' on the fitted coefficients ",
      "are unavailable. Refit the model with the current RoBMA/BayesTools ",
      "build."
    )
  )
})

test_that("BayesTools refusals are matched by their condition class", {

  status <- function(condition, reason) {
    data.frame(
      value = 0, eligible = FALSE, condition = condition, reason = reason,
      continuous_behavior = "regular", stringsAsFactors = FALSE
    )
  }
  robma <- structure(list(), class = c("RoBMA", "brma"))

  # A prior point mass at the value names the inclusion Bayes factor for
  # model-averaged objects whatever the wording of the BayesTools message;
  # the refusal keeps the BayesTools class.
  density <- BayesTools::prior("normal", list(0, 1))
  target <- list(
    prior_density     = density,
    status            = status("BayesTools_point_mass_at_null", "Reworded refusal."),
    point_mass_reason = .hypothesis_plan_point_mass_reason(robma)
  )
  refusal <- .hypothesis_plan_target_refusal(target)
  expect_identical(
    refusal[["class"]],
    c("BayesTools_point_mass_at_null", "BayesTools_hypothesis_ordinate")
  )
  expect_match(
    refusal[["reason"]],
    paste0(
      "This parameter has a null component, so its evidence against the ",
      "null is the inclusion Bayes factor reported by 'summary()' and ",
      "'summary_models()'."
    ),
    fixed = TRUE
  )
  expect_true(refusal[["value_specific"]])
  expect_error(.hypothesis_stop(refusal), class = "BayesTools_point_mass_at_null")

  # Other classes, and point masses of single models, keep the BayesTools
  # message.
  infinite <- .hypothesis_plan_target_refusal(list(
    prior_density     = density,
    status            = status("BayesTools_infinite_ordinate", "Infinite ordinate."),
    point_mass_reason = .hypothesis_plan_point_mass_reason(robma)
  ))
  expect_identical(infinite[["reason"]], "Infinite ordinate.")
  expect_identical(infinite[["class"]][[1L]], "BayesTools_infinite_ordinate")
  single <- .hypothesis_plan_target_refusal(list(
    prior_density     = density,
    status            = status("BayesTools_point_mass_at_null", "Point mass."),
    point_mass_reason = .hypothesis_plan_point_mass_reason(structure(list(), class = "brma"))
  ))
  expect_identical(single[["reason"]], "Point mass.")

  # Missing monitored columns of a prior are recognized by class only.
  missing_columns <- structure(
    class = c(
      "BayesTools_missing_monitored_columns", "BayesTools_marglik_input",
      "error", "condition"
    ),
    list(message = "Reworded missing columns.", call = NULL)
  )
  expect_true(.iwmde_marglik_parameters_missing(missing_columns))
  expect_false(.iwmde_marglik_parameters_missing(simpleError(
    "'samples' does not contain all monitored parameters."
  )))
})

test_that("qCMDE and IWMDE are RoBMA's precomputed density methods", {

  expect_true(.density_method_uses_precomputed("qCMDE"))
  expect_true(.density_method_uses_precomputed("iwmde"))
  expect_false(.density_method_uses_precomputed("KDE"))
  expect_false(.density_method_uses_precomputed("normal", allow_normal = TRUE))
})
