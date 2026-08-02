context("BayesTools parameter catalog contract")


.test_parameter_catalog <- function() {

  quantities <- data.frame(
    quantity_id = c(
      "BayesTools::mu_x",
      "BayesTools::log_tau_x",
      "BayesTools::random_sd_study_intercept"
    ),
    canonical_name    = c("mu_x", "log_tau_x", "random_sd_study_intercept"),
    provider          = "BayesTools",
    namespace         = c("mu", "log_tau", "mu"),
    role              = c("fixed_coefficient", "fixed_coefficient", "random_sd"),
    formula_parameter = c("mu", "log_tau", "mu"),
    term              = c("x", "x", "intercept"),
    component         = c("mu_x", "log_tau_x", "intercept"),
    display_label     = c("mu_x", "log_tau_x", "sd(study,intercept)"),
    fitted_scale      = "fitted_original",
    display_scale     = "original",
    status            = c("sampled", "sampled", "derived"),
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
    alias       = character(),
    quantity_id = character(),
    namespace   = character(),
    component   = character(),
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
    priors = list(outcome = list())
  )
  checked <- NULL
  testthat::local_mocked_bindings(
    JAGS_validate_fit_contract = function(fit, requires) {
      checked <<- requires
      invisible(TRUE)
    },
    parameter_catalog = function(object, ...) .test_parameter_catalog(),
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
  expect_s3_class(mods[["selection"]], "BayesTools_parameter_selection")
  expect_error(
    .brma_parameter_select_entry(object, "x"),
    class = "BayesTools_parameter_ambiguous"
  )
  expect_error(
    .brma_parameter_select_entry(object, "missing", component = "mods"),
    class = "BayesTools_parameter_not_found"
  )
})
