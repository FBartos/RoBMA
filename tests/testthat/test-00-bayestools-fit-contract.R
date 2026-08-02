context("BayesTools fitted metadata contract")


test_that("fitted formula identities come from the versioned name map", {

  fit <- structure(list(), class = "BayesTools_fit")
  object <- list(fit = fit)
  checked <- NULL
  name_map <- data.frame(
    encoded_name      = c("BT1_intercept", "BT1_first", "BT1_second"),
    jags_name         = c("mu_intercept", "opaque_backend_1", "opaque_backend_2"),
    kind              = "fixed",
    formula_parameter = "mu",
    term              = c("intercept", "x:a_b", "x__xXx__a_b"),
    role              = "coefficient",
    stringsAsFactors  = FALSE
  )
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(fit, requires) {
      checked <<- requires
      invisible(TRUE)
    },
    JAGS_formula_name_map = function(fit, parameter) name_map,
    JAGS_parameter_names = function(...) {
      stop("backend names must not be reconstructed")
    },
    .package = "BayesTools"
  )

  observed <- list()
  add <- function(parameter, component, term, aliases,
                  source, formula_parameter) {
    observed[[length(observed) + 1L]] <<- list(
      parameter         = parameter,
      component         = component,
      term              = term,
      aliases           = aliases,
      source            = source,
      formula_parameter = formula_parameter
    )
  }
  expect_silent(.brma_parameter_catalog_terms(
    object            = object,
    model_parameter   = "mu",
    component         = "mods",
    formula_parameter = "mu",
    source            = "mods",
    add               = add
  ))

  expect_setequal(
    checked,
    c(
      "name_encoding",
      "formula_name_map",
      "formula_design",
      "parameter_registry"
    )
  )
  expect_identical(
    vapply(observed, `[[`, character(1), "parameter"),
    c("opaque_backend_1", "opaque_backend_2")
  )
  expect_identical(
    vapply(observed, `[[`, character(1), "term"),
    c("x:a_b", "x__xXx__a_b")
  )
})


test_that("fitted metadata consumers fail closed on stale contracts", {

  object <- list(fit = structure(list(), class = "BayesTools_fit"))
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(...) {
      stop("Refit the model with this version of BayesTools.", call. = FALSE)
    },
    .package = "BayesTools"
  )

  expect_error(
    .fitted_formula_name_map(object, "mu"),
    "Refit the model"
  )
})
