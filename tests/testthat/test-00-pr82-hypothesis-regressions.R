test_that("incomplete coefficient transforms are reported as unsupported", {

  for (sources in list(c(a = "identity"), character(), c(a = NA_character_))) {
    transform <- structure(list(target_scale = "original",
      matrix = matrix(c(1, 1), 1L, dimnames = list("mu", c("a", "b"))),
      source_transforms = sources, output_transforms = c(mu = "identity")),
      class = "BayesTools_formula_coefficient_transform")
    route <- .hypothesis_brma_formula_transform_route(
      list(transform = transform, target = "mu", target_i = 1L)
    )
    expect_identical(route[["type"]], "unsupported")
    expect_match(route[["reason"]], "not supported by hypothesis()", fixed = TRUE)
  }
})

test_that("cross-level hypotheses keep all rows and build one shared context", {

  calls <- 0L
  testthat::local_mocked_bindings(
    hypothesis_level_contrast = function(...) list(
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
  expect_identical(rownames(result), c("mu_group", "mu_group.1"))
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
