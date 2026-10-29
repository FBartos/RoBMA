context("BayesTools hypothesis AST contract")


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
    class = "BayesTools_hypothesis_ambiguous"
  )
  # Statements resolved together that name several marginal-means parameters
  # are ambiguous too (hypothesis() resolves each statement on its own, so one
  # statement is the case that reaches this refusal); both have the classes of
  # the ambiguous fitted-model references: a statement problem, not an
  # unavailable test. Only the ambiguous alias is resolved by 'parameter'.
  several <- BayesTools::hypothesis_parse(c("intercept > 0", "mu_mu > 0"))
  expect_error(
    .hypothesis_marginal_means_select_parameter(object, several, NULL),
    "multiple marginal-means parameters (mu_intercept, mu_mu), but one statement must refer to the marginal means of one parameter",
    fixed = TRUE,
    class = "BayesTools_hypothesis_ambiguous"
  )
  for (statement in list(hypothesis, several)) {
    error <- tryCatch(
      .hypothesis_marginal_means_select_parameter(object, statement, NULL),
      error = identity
    )
    expect_identical(
      class(error),
      c("BayesTools_hypothesis_ambiguous", "BayesTools_hypothesis_statement",
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
      c("BayesTools_hypothesis_statement", "BayesTools_hypothesis_no_parameters",
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
    c("BayesTools_hypothesis_statement", "BayesTools_parameter_not_found",
      "BayesTools_parameter_resolution_error", "error", "condition")
  )
  expect_match(
    conditionMessage(error),
    "Could not infer a marginal-means parameter from the hypothesis.",
    fixed = TRUE
  )
})


