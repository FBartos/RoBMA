context("BayesTools parameter catalog contract")


.test_parameter_catalog <- function() {

  quantity <- function(canonical_name, namespace, role,
                       formula_parameter = "", term = "", component = "",
                       display_label = canonical_name, status = "sampled",
                       fixed_value = NA_real_, extraction_key = NULL,
                       owner_type = "",
                       owner_name = "", quantity_name = "",
                       arguments = character(), source_type = "none") {

    if (is.null(extraction_key)) {
      extraction_key <- list(
        type         = "coordinate",
        dependencies = canonical_name
      )
    }
    BayesTools:::.bt_parameter_catalog_quantity(
      canonical_name    = canonical_name,
      namespace         = namespace,
      role              = role,
      formula_parameter = formula_parameter,
      term              = term,
      component         = component,
      display_label     = display_label,
      display_scale     = "original",
      status            = status,
      fixed_value       = fixed_value,
      owner_type        = owner_type,
      owner_name        = owner_name,
      quantity          = quantity_name,
      arguments         = arguments,
      source_type       = source_type,
      extraction_key    = extraction_key
    )
  }
  quantities <- do.call(rbind, list(
    quantity("mu_x", "mu", "fixed_coefficient", "mu", "x", "mu_x"),
    quantity(
      "mu_f[1]", "mu", "fixed_coefficient", "mu", "f", "A",
      "(mu) f[A]"
    ),
    quantity(
      "mu_f[2]", "mu", "fixed_coefficient", "mu", "f", "B",
      "(mu) f[B]"
    ),
    quantity(
      "mu_f[dif: C]", "mu", "fixed_coefficient", "mu", "f", "C",
      "(mu) f[C]", "structural", 0,
      list(type = "factor_level", dependencies = character(), weights = numeric())
    ),
    quantity(
      "log_tau_x", "log_tau", "fixed_coefficient", "log_tau", "x",
      "log_tau_x"
    ),
    quantity("PET", "model", "parameter"),
    quantity(
      "random_sd_study_intercept", "mu", "random_sd", "mu", "intercept",
      "intercept", "(mu) study: sd(intercept)",
      owner_type    = "random_block",
      owner_name    = "study",
      quantity_name = "sd",
      arguments     = "intercept",
      source_type   = "identity"
    )
  ))
  out <- BayesTools:::.bt_parameter_catalog_new(
    quantities = quantities,
    aliases    = BayesTools:::.bt_parameter_catalog_aliases(quantities)
  )
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
  random  <- .brma_parameter_select_entry(
    object,
    "study: sd(intercept)",
    component = "random"
  )

  expect_setequal(
    checked,
    c(
      "name_encoding", "formula_name_map", "formula_design",
      "parameter_map"
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
  native_factor <- factor[["selection"]][["quantities"]][["extraction_key"]][[1L]]
  expect_identical(native_factor[["dependencies"]], c("mu_f[1]", "mu_f[2]"))
  factor_hypothesis <- .hypothesis_brma_select_parameter(
    object     = object,
    hypothesis = "f[B] > f[C]",
    component  = "auto"
  )
  expect_identical(factor_hypothesis[["parameter"]], "mu_f")
  expect_setequal(
    unique(factor_hypothesis[["resolution"]][["occurrences"]][["level"]]),
    c("B", "C")
  )
  expect_true(all(
    factor_hypothesis[["resolution"]][["occurrences"]][["quantity_id"]] %in%
      factor[["member_quantity_ids"]][[1L]]
  ))
  resolved_catalog <- .brma_parameter_catalog_metadata(object)[["catalog"]]
  level_b <- BayesTools::parameter_catalog_resolve(
    resolved_catalog,
    alias     = "f",
    component = "B"
  )
  level_c <- BayesTools::parameter_catalog_resolve(
    resolved_catalog,
    alias     = "f",
    component = "C"
  )
  expect_identical(level_b[["quantities"]][["provider"]], "BayesTools")
  expect_identical(level_b[["quantities"]][["status"]], "sampled")
  expect_identical(level_c[["quantities"]][["provider"]], "BayesTools")
  expect_identical(level_c[["quantities"]][["status"]], "structural")
  expect_identical(
    pet[["selection"]][["quantities"]][["provider"]],
    "BayesTools"
  )
  expect_identical(
    random[["selection"]][["quantities"]][["provider"]],
    "BayesTools"
  )
  expect_identical(random[["status"]], "sampled")
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
