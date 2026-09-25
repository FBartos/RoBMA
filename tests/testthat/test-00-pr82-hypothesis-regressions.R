test_that("coefficient transform routes follow the declared BayesTools map type", {

  transform <- function(targets) {
    structure(list(target_scale = "original",
      matrix = matrix(c(1, 1), 1L, dimnames = list("mu", c("a", "b"))),
      source_transforms = c(a = "identity", b = "log"),
      output_transforms = c(mu = "identity"), targets = targets),
      class = "BayesTools_formula_coefficient_transform")
  }
  route <- function(targets) .hypothesis_brma_formula_transform_route(
    list(transform = transform(targets), target = "mu", target_i = 1L)
  )

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

  calls <- 0L
  testthat::local_mocked_bindings(
    hypothesis_linear_target = function(...) list(
      posterior = 1:4, hypothesis = c("contrast = 0", "contrast = 1"),
      parameter = "contrast", weights = c(a = 1, b = -1)),
    hypothesis_BF = function(...) data.frame(BF = c(2, 3)),
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .check_iwmde_available = function(...) NULL,
    .iwmde_check_point_ordinate_supported = function(...) NULL,
    .iwmde_context = function(...) { calls <<- calls + 1L; list() },
    .hypothesis_brma_attach_iwmde_scalar = function(posterior, ...) posterior,
    .hypothesis_brma_append_iwmde_warnings = function(table, ...) table,
    .package = "RoBMA"
  )
  result <- .hypothesis_brma_level_contrast_BF(
    object = structure(list(fit = list(TRUE)), class = "brma"),
    posterior = NULL, hypothesis = "a - b = 0", parameter = "mu_group",
    density_method = "qCMDE", density_control = list(n_points = 50L),
    logBF = FALSE, BF01 = FALSE, seed = NULL, columns = NULL
  )
  expect_identical(result[["BF"]], c(2, 3))
  expect_identical(rownames(result), c("mu_group (1)", "mu_group (2)"))
  expect_identical(calls, 1L)
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
    list(), "mu = 0", "auto", list(catalog = list(), entries = data.frame())
  ), paste0("Resolved hypothesis metadata are unavailable. Refit the model with ",
            "the current RoBMA/BayesTools build."), fixed = TRUE)
})

test_that("BayesTools refusals are matched by their condition class", {

  classed <- function(class, message) {
    structure(
      class = c(class, "BayesTools_hypothesis_ordinate", "error", "condition"),
      list(message = message, call = NULL)
    )
  }
  robma <- structure(list(), class = c("RoBMA", "brma"))

  # Prior and posterior point masses at the null name the inclusion Bayes
  # factor whatever the wording of the BayesTools message.
  for (class in c("BayesTools_point_mass_at_null",
                  "BayesTools_posterior_point_mass_at_null")) {
    expect_error(
      .hypothesis_brma_stop_point_mass(robma, classed(class, "Reworded refusal.")),
      paste0(
        "Reworded refusal. This parameter has a null component, so its ",
        "evidence against the null is the inclusion Bayes factor reported by ",
        "'summary()' and 'summary_models()'."
      ),
      fixed = TRUE,
      info = class
    )
  }

  # Other conditions, and point-mass messages without the class, are
  # rethrown unchanged.
  unclassed <- simpleError(paste0(
    "There is a point mass in the prior at the exact null hypothesis value. ",
    "The Savage-Dickey density ratio is invalid."
  ))
  expect_identical(
    tryCatch(.hypothesis_brma_stop_point_mass(robma, unclassed), error = identity),
    unclassed
  )
  infinite <- classed("BayesTools_infinite_ordinate", "Infinite ordinate.")
  expect_identical(
    tryCatch(.hypothesis_brma_stop_point_mass(robma, infinite), error = identity),
    infinite
  )
  point_mass <- classed("BayesTools_point_mass_at_null", "Point mass.")
  expect_identical(
    tryCatch(
      .hypothesis_brma_stop_point_mass(structure(list(), class = "brma"), point_mass),
      error = identity
    ),
    point_mass
  )

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
