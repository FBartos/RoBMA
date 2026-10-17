test_that("ordered verifier requires actual declared and exercised fixture roles", {

  skip_on_cran()
  env <- environment()
  sys.source(testthat::test_path("..", "..", "tools", "ordered-source-verification.R"), env)
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
  sys.source(testthat::test_path("..", "..", "tools", "ordered-source-verification.R"), environment())
  counts <- c(source = 5L, grid = 4L, proof = 12L, physical_block = 12L, supplier = 2L)
  expect_silent(.ordered_source_verification_counts(counts))
  for (name in names(counts)) {
    incomplete <- counts
    incomplete[[name]] <- incomplete[[name]] - 1L
    expect_error(.ordered_source_verification_counts(incomplete), "coverage differs")
  }
})
