context("BayesTools hypothesis AST contract")


test_that("fitted hypotheses resolve and rewrite the parsed AST", {

  quantities <- data.frame(
    quantity_id       = "BayesTools::mu_x",
    canonical_name    = "mu_x",
    provider          = "BayesTools",
    namespace         = "mu",
    role              = "fixed_coefficient",
    formula_parameter = "mu",
    term              = "x",
    component         = "mods",
    display_label     = "mu_x",
    fitted_scale      = "fitted_original",
    display_scale     = "original",
    status            = "sampled",
    fixed_value       = NA_real_,
    internal          = FALSE,
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
  quantities[["extraction_key"]] <- I(list(list(
    type         = "registry",
    dependencies = "mu_x"
  )))
  aliases <- data.frame(
    alias       = c("x", "mu_x"),
    quantity_id = "BayesTools::mu_x",
    namespace   = "mods",
    component   = "mods",
    stringsAsFactors = FALSE
  )
  catalog <- list(
    schema_version = 1L,
    quantities     = quantities,
    aliases        = aliases
  )
  class(catalog) <- c("BayesTools_parameter_catalog", "list")
  entries <- data.frame(
    quantity_id       = "BayesTools::mu_x",
    parameter         = "mu_x",
    component         = "mods",
    term              = "x",
    source            = "mods",
    formula_parameter = "mu",
    role              = "fixed_coefficient",
    status            = "sampled",
    fixed_value       = NA_real_,
    stringsAsFactors  = FALSE,
    check.names       = FALSE
  )
  entries[["aliases"]] <- I(list(c("x", "mu_x")))
  ast <- BayesTools::hypothesis_parse("x > abs(x)")
  testthat::local_mocked_bindings(
    .brma_parameter_catalog_metadata = function(object) {
      list(catalog = catalog, entries = entries)
    },
    .brma_parameter_catalog = function(...) {
      stop("legacy catalog resolution is forbidden")
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    hypothesis_parse = function(...) stop("hypothesis was parsed twice"),
    .package = "BayesTools"
  )

  selected <- .hypothesis_brma_select_parameter(
    object     = list(),
    hypothesis = ast,
    component  = "mods"
  )
  rewritten <- .hypothesis_brma_rewrite(
    hypothesis = ast,
    aliases    = selected[["aliases"]],
    parameter  = selected[["parameter"]]
  )

  expect_identical(selected[["parameter"]], "mu_x")
  expect_s3_class(selected[["resolution"]], "BayesTools_hypothesis_resolution")
  expect_s3_class(rewritten, "BayesTools_hypothesis_ast")
  expect_identical(
    BayesTools::hypothesis_render(rewritten),
    "mu_x > abs(mu_x)"
  )
})


test_that("RoBMA no longer owns hypothesis expression parsing", {

  obsolete <- c(
    ".hypothesis_brma_symbols_normalized",
    ".hypothesis_brma_replace_symbol",
    ".hypothesis_brma_level_ref"
  )
  present <- vapply(
    obsolete,
    exists,
    logical(1),
    envir    = asNamespace("RoBMA"),
    inherits = FALSE
  )
  expect_false(any(present))
})
