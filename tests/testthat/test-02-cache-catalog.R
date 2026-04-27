context("Cached fit catalog")

source(testthat::test_path("common-functions.R"))


test_that("fit catalog is internally consistent", {

  catalog <- fit_catalog()

  expect_equal(anyDuplicated(catalog[["name"]]), 0L,
               info = "fit catalog names must be unique")
  expect_true(all(file.exists(testthat::test_path(catalog[["source_file"]]))),
              info = "all catalog source files must exist")
  expect_true(all(catalog[["has_loo"]]),
              info = "all cached test fits should include pre-computed LOO")
  expect_true(all(!catalog[["has_waic"]]),
              info = "WAIC is not pre-computed for cached test fits")
  expect_true(all(!catalog[catalog[["class"]] %in% c("BMA.norm", "BMA.glmm", "RoBMA"), "has_marglik"]),
              info = "product-space model averaging fits must not cache marginal likelihoods")
  expect_true(all(catalog[!catalog[["class"]] %in% c("BMA.norm", "BMA.glmm", "RoBMA"), "has_marglik"]),
              info = "single-model cached fits should include marginal likelihoods")
})

test_that("cached fits are catalogued and valid", {

  cached_names <- list_fits(validate = FALSE)
  if (length(cached_names) == 0) {
    skip("No cached fits available.")
  }

  catalog_names <- fit_catalog()[["name"]]
  uncatalogued  <- setdiff(cached_names, catalog_names)

  expect_equal(uncatalogued, character(),
               info = "all cached fits must be listed in fit_catalog()")

  problems <- lapply(cached_names, validate_cached_fit, deep = TRUE)
  names(problems) <- cached_names
  problems <- problems[lengths(problems) > 0]
  problem_text <- if (length(problems) == 0) {
    ""
  } else {
    paste(
      paste0(names(problems), ": ", vapply(problems, paste, "", collapse = "; ")),
      collapse = "\n"
    )
  }

  expect_true(
    length(problems) == 0,
    info = paste(
      "cached fits must match catalog requirements; rerun test-01 after model definition changes",
      problem_text,
      sep = "\n"
    )
  )
})
