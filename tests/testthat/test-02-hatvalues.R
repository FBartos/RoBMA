context("Hatvalues")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(hatvalue_cases(), function(case) {
  test_that_case("Hatvalues match metafor", case, {
    expect_hatvalues_match_metafor(case)
  })
})

test_that("Hatvalues reject unsupported single-model families", {

  model_names <- c("nielweise2008_glmm", "dat.lehmann2018-3PSM")
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    expect_error(hatvalues(fits[[name]]), "available", info = name)
  }
})

test_that("Hatvalues for BMA.norm fits are internally consistent", {

  model_names <- c("dat.lehmann2018_BMA.norm", "dat.lehmann2018_BMA.norm_mods")
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    fit <- fits[[name]]
    expect_hatvalues_vector(hatvalues(fit), nobs(fit), info = name)
  }
})

test_that("Hatvalues use study labels as names", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit <- fits[[name]]
  expect_equal(names(hatvalues(fit)), RoBMA:::.diagnostic_study_labels(fit))
})

test_that("Hatvalues reject unsupported model-averaging families", {

  model_names <- c("bcg_BMA.glmm", "dat.lehmann2018_RoBMA")
  skip_if_missing_fits(model_names)

  expect_error(hatvalues(fits[["bcg_BMA.glmm"]]), "normal outcome models")
  expect_error(hatvalues(fits[["dat.lehmann2018_RoBMA"]]), "selection models")
})
