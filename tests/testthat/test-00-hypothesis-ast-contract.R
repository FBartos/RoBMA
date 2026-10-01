context("BayesTools hypothesis AST contract")


test_that("fitted hypotheses resolve and rewrite the structured AST", {

  quantities <- BayesTools:::.bt_parameter_catalog_quantity(
    canonical_name    = "mu_x",
    namespace         = "mu",
    role              = "fixed_coefficient",
    formula_parameter = "mu",
    term              = "x",
    component         = "mods",
    label_parts       = catalog_label_parts("mu_x", "x", "mu"),
    display_scale     = "original",
    status            = "sampled",
    extraction_key    = list(
      type         = "coordinate",
      dependencies = "mu_x"
    )
  )
  catalog <- BayesTools:::.bt_parameter_catalog_new(
    quantities = quantities,
    aliases    = BayesTools:::.bt_parameter_catalog_aliases(quantities)
  )
  entries <- data.frame(
    quantity_id       = quantities[["quantity_id"]],
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
  entries[["member_quantity_ids"]] <- I(list(character()))
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


test_that("marginal-means aliases never silently resolve collisions", {

  object <- list(term_map = data.frame(
    term      = c("intercept", "mu"),
    parameter = c("mu_intercept", "mu_mu"),
    label     = c("intercept", "mu"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  hypothesis <- BayesTools::hypothesis_parse("mu > 0")

  expect_error(
    .hypothesis_marginal_means_select_parameter(object, hypothesis, NULL),
    "alias 'mu' is ambiguous.*Specify 'parameter'",
    class = "RoBMA_hypothesis_ambiguous"
  )
  # Statements on several marginal-means parameters are ambiguous too; both
  # have the classes of the ambiguous fitted-model references: a statement
  # problem ('parameter' selects one), not an unavailable test.
  several <- BayesTools::hypothesis_parse(c("intercept > 0", "mu_mu > 0"))
  for (statement in list(hypothesis, several)) {
    error <- tryCatch(
      .hypothesis_marginal_means_select_parameter(object, statement, NULL),
      error = identity
    )
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_ambiguous", "RoBMA_hypothesis_statement",
        "error", "condition")
    )
  }

  intercept <- .hypothesis_marginal_means_select_parameter(
    object,
    hypothesis,
    parameter = "mu"
  )
  moderator <- .hypothesis_marginal_means_select_parameter(
    object,
    hypothesis,
    parameter = "mu_mu"
  )
  expect_identical(intercept[["parameter"]], "mu_intercept")
  expect_identical(moderator[["parameter"]], "mu_mu")

  # References that resolve to no marginal-means parameter are statement
  # errors with the classes with which the fitted-model path refuses them: a
  # statement without parameter symbols (also with 'parameter') has those of
  # BayesTools_hypothesis_no_parameters, an unknown name those of
  # BayesTools_parameter_not_found.
  for (parameter in list(NULL, "mu_mu")) {
    error <- tryCatch(
      .hypothesis_marginal_means_select_parameter(
        object, BayesTools::hypothesis_parse("1 > 0"), parameter
      ),
      error = identity
    )
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_statement", "BayesTools_hypothesis_no_parameters",
        "BayesTools_parameter_resolution_error", "error", "condition")
    )
    expect_identical(
      conditionMessage(error),
      "Hypothesis must reference a marginal-means parameter."
    )
  }
  error <- tryCatch(
    .hypothesis_marginal_means_select_parameter(
      object, BayesTools::hypothesis_parse("foo > 0"), NULL
    ),
    error = identity
  )
  expect_identical(
    class(error),
    c("RoBMA_hypothesis_statement", "BayesTools_parameter_not_found",
      "BayesTools_parameter_resolution_error", "error", "condition")
  )
  expect_match(
    conditionMessage(error),
    "Could not infer a marginal-means parameter from the hypothesis.",
    fixed = TRUE
  )
})


test_that("hypothesis aliases are rewritten independently by statement", {

  rewritten <- .hypothesis_brma_rewrite(
    hypothesis = BayesTools::hypothesis_parse(c("effect > 0", "mu < 0")),
    aliases    = list(effect = "mu", mu = "mu"),
    parameter  = "mu"
  )

  expect_identical(
    BayesTools::hypothesis_render(rewritten),
    c("mu > 0", "mu < 0")
  )
})


test_that("point-null references must be direct", {

  expect_error(
    .hypothesis_brma_point_refs(
      BayesTools::hypothesis_parse("2 * mu = 0"),
      parameter = "mu"
    ),
    "direct parameter or level reference"
  )
  expect_equal(
    .hypothesis_brma_point_refs(
      BayesTools::hypothesis_parse("mu = 0"),
      parameter = "mu"
    )[["value"]],
    0
  )
  precise <- BayesTools::hypothesis_parse("mu = 0.701406683025")
  expect_identical(
    .hypothesis_brma_point_refs(precise, parameter = "mu")[["value"]],
    precise[["statements"]][[1L]][["left"]][["value"]]
  )
  expect_equal(
    nrow(.hypothesis_brma_point_refs(
      BayesTools::hypothesis_parse("2 * mu > 0"),
      parameter = "mu"
    )),
    0L
  )
})


test_that("only linear combinations of factor levels bypass the direct guard", {

  plan_of <- function(statement) {
    ast <- BayesTools::hypothesis_parse(statement)
    .hypothesis_plan_new(
      statement = ast, hypothesis = ast, parameter = "mu", label = "mu",
      component = "mods"
    )
  }
  # Scalar parameters need a direct point reference; factor levels are
  # linear targets whose combinations BayesTools compiles
  # (.hypothesis_plan_linear()).
  # The statement has to be restated: a statement error, not an unavailable
  # test.
  refusal <- .hypothesis_plan_direct_refusal(plan_of("2 * mu = 0"))
  expect_identical(refusal[["class"]], "RoBMA_hypothesis_statement")
  expect_match(refusal[["reason"]], "direct parameter or level reference", fixed = TRUE)
  expect_null(.hypothesis_plan_direct_refusal(plan_of("mu = 0")))
  expect_null(.hypothesis_plan_direct_refusal(plan_of("mu > 0")))
})

