context("DFBETAS")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(dfbetas_metafor_cases(), function(case) {
  test_that_case("DFBETAS match metafor", case, {
    expect_dfbetas_match_metafor(case)
  })
})

test_that("DFBETAS for location-scale models expose location and scale terms", {

  name <- "bangertdrowns2004_location-scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  expect_dfbetas_table(dfbetas(fit_brma), nobs(fit_brma), info = name)
  expect_dfbetas_table(dfbetas(fit_brma, type = "scale"), nobs(fit_brma),
                       info = paste(name, "scale"))
})

test_that("DFBETAS outputs use study-label row names", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  labels   <- RoBMA:::.diagnostic_study_labels(fit_brma)

  expect_equal(rownames(dfbetas(fit_brma)), labels)
  expect_equal(
    rownames(dfbetas(fit_brma, return_loo_estimates = TRUE)),
    labels
  )
})

test_that("DFBETAS for selection models expose bias terms", {

  model_names <- c("dat.lehmann2018-3PSM", "dat.lehmann2018-3PSM_neg", "dat.lehmann2018-3PSMreg")
  skip_if_missing_fits(model_names)

  fit_pos <- fits[["dat.lehmann2018-3PSM"]]
  fit_neg <- fits[["dat.lehmann2018-3PSM_neg"]]
  dfb_pos <- dfbetas(fit_pos)
  dfb_neg <- dfbetas(fit_neg)

  expect_equal(dfb_pos[-4, 1], -dfb_neg[-4, 1], tolerance = 0.01,
               info = "positive and negative selection DFBETAS flip")

  for (name in c("dat.lehmann2018-3PSM", "dat.lehmann2018-3PSMreg")) {
    fit_brma <- fits[[name]]
    expect_dfbetas_table(dfbetas(fit_brma), nobs(fit_brma), info = name)

    bias_dfbetas <- dfbetas(fit_brma, type = "bias")
    expect_dfbetas_table(bias_dfbetas, nobs(fit_brma), info = paste(name, "bias"))
    expect_true(any(grepl("^omega", colnames(bias_dfbetas))), info = name)
  }
})

test_that("DFBETAS for PET and PEESE expose publication-bias terms", {

  bias_cases <- data.frame(
    name = c(
      "dat.lehmann2018-PET",
      "dat.lehmann2018-PETreg",
      "dat.lehmann2018-PET_neg",
      "dat.lehmann2018-PEESE",
      "dat.lehmann2018-PEESEreg",
      "dat.lehmann2018-PEESE_neg"
    ),
    column = c("PET", "PET", "PET", "PEESE", "PEESE", "PEESE"),
    stringsAsFactors = FALSE
  )
  skip_if_missing_fits(bias_cases[["name"]])

  for (i in seq_len(nrow(bias_cases))) {
    name          <- bias_cases[["name"]][[i]]
    expected_col  <- bias_cases[["column"]][[i]]
    fit_brma      <- fits[[name]]
    bias_dfbetas  <- dfbetas(fit_brma, type = "bias")

    expect_dfbetas_table(bias_dfbetas, nobs(fit_brma), info = paste(name, "bias"))
    expect_true(expected_col %in% colnames(bias_dfbetas), info = name)
  }
})

test_that("DFBETAS for GLMM fits are finite", {

  model_names <- c("nielweise2008_glmm", "bcg_glmm_reg")
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    fit_brma <- fits[[name]]
    expect_dfbetas_table(dfbetas(fit_brma), nobs(fit_brma), info = name)
  }
})

test_that("DFBETAS for model-averaging fits are internally consistent", {

  cases <- data.frame(
    name = c(
      "dat.lehmann2018_BMA.norm",
      "dat.lehmann2018_BMA.norm_mods",
      "dat.lehmann2018_BMA.norm_scale",
      "bcg_BMA.glmm",
      "nielweise2008_BMA.glmm",
      "dat.lehmann2018_RoBMA_mods",
      "dat.lehmann2018_RoBMA_3lvl_mods_scale"
    ),
    type = c(NA, NA, "scale", NA, NA, NA, NA),
    min_cols = c(1, 2, 2, 1, 1, 1, 1),
    stringsAsFactors = FALSE
  )
  skip_if_missing_fits(c(cases[["name"]], "dat.lehmann2018_RoBMA"))

  for (i in seq_len(nrow(cases))) {
    name     <- cases[["name"]][[i]]
    type     <- cases[["type"]][[i]]
    fit_brma <- fits[[name]]
    dfb      <- if (is.na(type)) dfbetas(fit_brma) else dfbetas(fit_brma, type = type)

    expect_dfbetas_table(dfb, nobs(fit_brma), min_cols = cases[["min_cols"]][[i]],
                         info = name)
  }

  bias_dfbetas <- dfbetas(fits[["dat.lehmann2018_RoBMA"]], type = "bias")
  expect_dfbetas_table(bias_dfbetas, nobs(fits[["dat.lehmann2018_RoBMA"]]),
                       min_cols = 3, info = "RoBMA bias")
  expect_true(any(grepl("^omega", colnames(bias_dfbetas))))
  expect_true("PET" %in% colnames(bias_dfbetas))
  expect_true("PEESE" %in% colnames(bias_dfbetas))
})
