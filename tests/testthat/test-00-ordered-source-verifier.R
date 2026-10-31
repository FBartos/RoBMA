.ordered_verifier_tool_path <- function() {

  path <- testthat::test_path("..", "..", "tools", "ordered-source-verification.R")
  if (!file.exists(path) && file.exists(testthat::test_path("..", "..", ".git"))) {
    stop("The ordered verifier helper is missing from the developer checkout.", call. = FALSE)
  }
  skip_if_not(file.exists(path), "Developer ordered verifier requires a source checkout.")
  path
}


test_that("ordered verifier requires actual declared and exercised fixture roles", {

  skip_on_cran()
  env <- environment()
  sys.source(.ordered_verifier_tool_path(), env)
  definitions <- as.list(parse(testthat::test_path("test-00-iwmde-ordered-targets.R")))
  for (node in definitions) {
    if (is.call(node) && identical(node[[1L]], as.name("<-")) && is.symbol(node[[2L]]) &&
        as.character(node[[2L]]) %in% c(".ordered_target_fixture", ".ordered_narrow_gaussian_fixture")) eval(node, env)
  }
  ns <- asNamespace("RoBMA")
  normal <- .ordered_narrow_gaussian_fixture()
  selected <- .ordered_narrow_gaussian_fixture(selected = TRUE)
  gamma <- .ordered_narrow_gaussian_fixture(total = BayesTools::prior("gamma", list(2, 1)))
  role <- function(name, fixture) .ordered_source_fixture_role(name,
    fixture[["context"]][["object"]], fixture[["context"]],
    BayesTools::JAGS_ordered_parameter_spec(fixture[["context"]][["object"]][["fit"]]), ns)
  expect_silent(role("normal", normal))
  expect_silent(role("gamma", gamma))
  expect_silent(role("selected", selected))
  expect_error(role("selected", normal), "declared step-selection")
  expect_error(role("gamma", normal), "declared gamma total")
  expect_error(role("shared", selected), "multiple live declared owners")
  expect_silent(.ordered_source_block_role("selected", selected[["context"]], selected[["target"]],
    selected[["state"]], selected[["block"]], ns))
  expect_error(.ordered_source_block_role("selected", normal[["context"]], normal[["target"]],
    normal[["state"]], normal[["block"]], ns), "actual row-state selection")
})

test_that("ordered verifier counts completed checks rather than labels", {

  skip_on_cran()
  sys.source(.ordered_verifier_tool_path(), environment())
  counts <- c(source = 5L, grid = 4L, proof = 12L, physical_block = 12L, supplier = 2L)
  expect_silent(.ordered_source_verification_counts(counts))
  for (name in names(counts)) {
    incomplete <- counts
    incomplete[[name]] <- incomplete[[name]] - 1L
    expect_error(.ordered_source_verification_counts(incomplete), "coverage differs")
  }
})


test_that("ordered verifier comparisons use own-value errors and bounded finite envelopes", {

  source(.ordered_verifier_tool_path(), local = environment())
  check <- function(log_producer = log(1.01), log_reference = log(.99),
      producer_error = 1 / 101, reference_error = 1 / 99,
      floating_error = 64 * .Machine$double.eps) {
    .ordered_source_comparison(log_producer, log_reference, producer_error,
      reference_error, floating_error, "unit", 3L, "numerator")
  }
  for (values in list(c(1.01, .99, 1 / 101, 1 / 99),
      c(.99, 1.01, 1 / 99, 1 / 101), c(1, 1, 0, 0))) {
    result <- check(log(values[[1L]]), log(values[[2L]]), values[[3L]], values[[4L]])
    ratio <- values[[1L]] / values[[2L]]
    expected <- values[[3L]] * ratio + values[[4L]] +
      64 * .Machine$double.eps * max(1, ratio)
    expect_equal(result$allowance, expected, tolerance = 1e-15)
    expect_lte(result$difference, result$allowance)
  }
  expect_identical(check(0, 0, .025, 0, 0)$allowance, .025)
  prefix <- "Ordered verification for fixture 'unit', row 3, numerator: "
  for (value in list(Inf, NA_real_, NaN, numeric(), c(.01, .02))) {
    expect_error(check(producer_error = value), paste0(prefix,
      "comparison values and error estimates must be finite scalars."), fixed = TRUE)
  }
  expect_error(check(producer_error = -.01),
    paste0(prefix, "error estimates must be nonnegative."), fixed = TRUE)
  for (value in c(Inf, NA_real_, NaN)) {
    expect_error(check(log_producer = value), paste0(prefix,
      "comparison values and error estimates must be finite scalars."), fixed = TRUE)
  }
  expect_error(check(log_producer = 1000, log_reference = 0), paste0(prefix,
    "discrepancy and comparison envelope must be finite."), fixed = TRUE)
  expect_error(check(0, 0, .03, 0, 0),
    paste0(prefix, "comparison envelope exceeds the 0.025 accuracy ceiling."), fixed = TRUE)
  expect_error(check(log_producer = log(1.1), producer_error = 0,
    reference_error = 0, floating_error = 0),
    paste0(prefix, "discrepancy exceeds its comparison envelope."), fixed = TRUE)
  valid <- list(value = .99, abs.error = .01)
  expect_identical(.ordered_source_reference_integral(valid, "unit", 3L, "reference"), valid)
  for (reference in c(lapply(c(0, -1, Inf, NA_real_), function(value)
      list(value = value, abs.error = .01)),
      lapply(c(-.01, Inf, NA_real_), function(error) list(value = 1, abs.error = error)))) {
    expect_error(.ordered_source_reference_integral(reference, "unit", 3L, "reference"),
      paste0("Ordered verification for fixture 'unit', row 3, reference: reference integral ",
        "must be finite and positive with a finite nonnegative error estimate."), fixed = TRUE)
  }
})
