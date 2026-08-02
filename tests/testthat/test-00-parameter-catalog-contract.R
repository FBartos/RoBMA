context("BayesTools parameter catalog contract")


.test_parameter_catalog <- function() {

  quantities <- data.frame(
    quantity_id = c(
      "BayesTools::mu_x",
      "BayesTools::mu_f_1",
      "BayesTools::mu_f_2",
      "BayesTools::log_tau_x",
      "BayesTools::PET",
      "BayesTools::random_sd_study_intercept"
    ),
    canonical_name    = c(
      "mu_x", "mu_f[1]", "mu_f[2]", "log_tau_x", "PET",
      "random_sd_study_intercept"
    ),
    provider          = "BayesTools",
    namespace         = c("mu", "mu", "mu", "log_tau", "model", "mu"),
    role              = c(
      "fixed_coefficient", "fixed_coefficient", "fixed_coefficient",
      "fixed_coefficient", "parameter", "random_sd"
    ),
    formula_parameter = c("mu", "mu", "mu", "log_tau", "", "mu"),
    term              = c("x", "f", "f", "x", "", "intercept"),
    component         = c(
      "mu_x", "mu_f[1]", "mu_f[2]", "log_tau_x", "", "intercept"
    ),
    display_label     = c(
      "mu_x", "mu_f[1]", "mu_f[2]", "log_tau_x", "PET",
      "sd(study,intercept)"
    ),
    fitted_scale      = "fitted_original",
    display_scale     = "original",
    status            = c(
      "sampled", "sampled", "sampled", "sampled", "sampled", "derived"
    ),
    fixed_value       = NA_real_,
    internal          = FALSE,
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
  quantities[["extraction_key"]] <- I(lapply(
    quantities[["canonical_name"]],
    function(name) list(type = "coordinate", dependencies = name)
  ))
  aliases <- data.frame(
    alias       = c("f", "f", "PET"),
    quantity_id = c(
      "BayesTools::mu_f_1", "BayesTools::mu_f_2", "BayesTools::PET"
    ),
    namespace   = c("mu", "mu", "model"),
    component   = c("mu_f[1]", "mu_f[2]", ""),
    stringsAsFactors = FALSE
  )
  out <- list(
    schema_version = 1L,
    quantities     = quantities,
    aliases        = aliases
  )
  class(out) <- c("BayesTools_parameter_catalog", "list")
  return(out)
}


.test_formula_name_map <- function(parameter) {

  out <- data.frame(
    encoded_name      = c("mu_x", "mu_f", "log_tau_x"),
    jags_name         = c("mu_x", "mu_f", "log_tau_x"),
    kind              = "fixed",
    formula_parameter = c("mu", "mu", "log_tau"),
    term              = c("x", "f", "x"),
    role              = "coefficient",
    stringsAsFactors  = FALSE
  )
  out <- out[out[["formula_parameter"]] == parameter, , drop = FALSE]
  class(out) <- c("BayesTools_formula_name_map", "data.frame")
  attr(out, "schema_version") <- 1L
  return(out)
}


test_that("fitted parameter discovery is metadata-only and component-aware", {

  scale <- data.frame(x = c(0, 1))
  attr(scale, "parameter") <- "log_tau"
  attr(scale, "source")    <- "tau"
  attr(scale, "aliases")   <- "tau"
  data <- list(
    outcome  = data.frame(yi = c(0, 1), sei = c(1, 1)),
    location = data.frame(x = c(0, 1)),
    scale    = scale
  )
  attr(data, "mods")   <- FALSE
  attr(data, "random") <- TRUE
  attr(data, "scale")  <- TRUE
  object <- list(
    data   = data,
    fit    = structure(list(sentinel = TRUE), class = "BayesTools_fit"),
    priors = list(outcome = list(
      bias = BayesTools::prior_PET("normal", list(mean = 0, sd = 1))
    ))
  )
  checked <- NULL
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(fit, requires) {
      checked <<- unique(c(checked, requires))
      invisible(TRUE)
    },
    parameter_catalog = function(object, ...) .test_parameter_catalog(),
    JAGS_formula_name_map = function(fit, parameter) {
      .test_formula_name_map(parameter)
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_bundle = function(...) {
      stop("posterior discovery is forbidden")
    },
    .package = "RoBMA"
  )

  catalog <- .brma_parameter_catalog(object)
  mods    <- .brma_parameter_select_entry(object, "x", component = "mods")
  scale   <- .brma_parameter_select_entry(object, "x", component = "scale")
  factor  <- .brma_parameter_select_entry(object, "f", component = "mods")
  pet     <- .brma_parameter_select_entry(object, "PET", component = "bias")

  expect_setequal(
    checked,
    c(
      "name_encoding", "formula_name_map", "formula_design",
      "parameter_registry", "parameter_catalog"
    )
  )
  expect_true(any(catalog[["component"]] == "random"))
  expect_identical(mods[["parameter"]], "mu_x")
  expect_identical(scale[["parameter"]], "log_tau_x")
  expect_identical(factor[["parameter"]], "mu_f")
  expect_identical(
    factor[["selection"]][["quantities"]][["provider"]],
    "RoBMA"
  )
  expect_identical(
    factor[["selection"]][["quantities"]][["role"]],
    "formula_coefficient_group"
  )
  expect_identical(
    pet[["selection"]][["quantities"]][["provider"]],
    "BayesTools"
  )
  expect_s3_class(mods[["selection"]], "BayesTools_parameter_selection")
  expect_error(
    .brma_parameter_select_entry(object, "x"),
    class = "BayesTools_parameter_ambiguous"
  )
  expect_error(
    .brma_parameter_select_entry(object, "missing", component = "mods"),
    class = "BayesTools_parameter_not_found"
  )

  materialized <- NULL
  testthat::local_mocked_bindings(
    JAGS_materialize_draws = function(fit, parameters) {
      materialized <<- list(fit = fit, parameters = parameters)
      "draws"
    },
    .package = "BayesTools"
  )
  class(object) <- "brma"
  expect_identical(
    BayesTools::parameter_draws(object, factor[["selection"]]),
    "draws"
  )
  expect_identical(materialized[["fit"]], object[["fit"]])
  expect_identical(materialized[["parameters"]], c("mu_f[1]", "mu_f[2]"))
})
