test_that("bias terms extend an empty fixed design with finite assignments", {

  object <- bPET(yi ~ 0, sei = sei,
                  data = data.frame(yi = c(.1, .2, .3), sei = c(.1, .2, .3)),
                  measure = "SMD", only_priors = TRUE, silent = TRUE)
  # Exercise an empty design, as supplied by the model-matrix fallback.
  empty <- .fitted_formula_design(object, "mu")
  empty[["model_matrix"]] <- matrix(numeric(), 3L, 0L)
  empty[["column_names"]] <- empty[["raw_column_names"]] <- character()
  empty[["assign"]] <- integer()
  empty[["rank"]] <- 0L
  empty[["aliased"]] <- logical()
  empty[["model_terms"]] <- character()
  testthat::local_mocked_bindings(
    .fitted_formula_design = function(...) empty, .package = "RoBMA"
  )
  design <- .get_model_matrix(object)
  expect_identical(attr(design, "assign"), 1L)
  expect_equal(as.numeric(design), c(.1, .2, .3), tolerance = 0)
})

