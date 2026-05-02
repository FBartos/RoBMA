context("Residuals")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

REFERENCE_DIR <<- testthat::test_path("..", "results", "residuals")

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(residual_metafor_cases(), function(case) {
  test_that_case("Residuals match metafor", case, {
    expect_residuals_match_metafor(case)
  })
})

test_that("Residuals for equivalent interaction parameterizations match", {

  model_names <- c("bcg_meta-regression3", "bcg_meta-regression3b")
  skip_if_missing_fits(model_names)

  fit_brma1 <- fits[["bcg_meta-regression3"]]
  fit_brma2 <- fits[["bcg_meta-regression3b"]]

  expect_equal(residuals(fit_brma1), residuals(fit_brma2), tolerance = 0.10,
               info = "outcome residuals match")
  expect_equal(
    residuals(fit_brma1, type = "pearson"),
    residuals(fit_brma2, type = "pearson"),
    tolerance = 0.10,
    info      = "pearson residuals match"
  )
  expect_equal(rstandard(fit_brma1), rstandard(fit_brma2), tolerance = 0.10,
               info = "rstandard residuals match")
})

test_that("Extended rstudent residuals for interaction parameterizations match", {

  skip_if_not_full_diagnostics(
    "Interaction parameterization LOO-PIT checks are redundant with core residual coverage."
  )
  model_names <- c("bcg_meta-regression4", "bcg_meta-regression4b")
  skip_if_missing_fits(model_names)

  brma_rstudent1 <- suppressWarnings(rstudent(fits[["bcg_meta-regression4"]]))
  brma_rstudent2 <- suppressWarnings(rstudent(fits[["bcg_meta-regression4b"]]))

  expect_true(cor(brma_rstudent1[, "z"], brma_rstudent2[, "z"]) > 0.90,
              info = "rstudent residuals match")
})

test_that("Extended rstudent residuals align with metafor oracles", {

  skip_if_not_full_diagnostics(
    "Default residual coverage keeps one representative LOO-PIT check per family."
  )

  extended <- list(
    list(
      name      = "bcg_meta-regression4",
      resid_cor = 0.90,
      se_cor    = NULL,
      z_cor     = 0.90
    ),
    list(
      name      = "konstantopoulos2011_3lvl2",
      resid_cor = 0.60,
      se_cor    = 0.80,
      z_cor     = 0.60
    ),
    list(
      name      = "dat.lehmann2018-PETreg",
      equal     = TRUE
    )
  )
  skip_if_missing_fits(vapply(extended, `[[`, "", "name"))

  for (case in extended) {
    name <- case[["name"]]
    if (isTRUE(case[["equal"]])) {
      .expect_residual_table_equal(
        name,
        suppressWarnings(rstudent(fits[[name]])),
        rstudent(info[[name]][["metafor"]]),
        tolerance   = 0.05,
        z_tolerance = 0.10,
        label       = "rstudent"
      )
    } else {
      .expect_residual_rstudent_rank(
        name,
        suppressWarnings(rstudent(fits[[name]])),
        rstudent(info[[name]][["metafor"]]),
        resid_cor = case[["resid_cor"]],
        se_cor    = case[["se_cor"]],
        z_cor     = case[["z_cor"]]
      )
    }
  }
})

test_that("Extended rstudent residuals without direct metafor oracle are calibrated", {

  skip_if_not_full_diagnostics(
    "These checks have weaker or family-specific oracles and duplicate default branches."
  )

  scale_name <- "bangertdrowns2004_location-scale"
  sel_name   <- "dat.lehmann2018-3PSMreg"
  glmm_name  <- "bcg_glmm_reg"
  skip_if_missing_fits(c(scale_name, sel_name, glmm_name))

  scale_rstudent <- rstudent(fits[[scale_name]])
  scale_standard <- rstandard(info[[scale_name]][["metafor"]])
  expect_equal(mean(scale_rstudent$z), 0, 0.05)
  expect_equal(sd(scale_rstudent$z), 1, 0.05)
  expect_true(cor(scale_rstudent$z, scale_standard$z, method = "spearman") > 0.8)

  sel_resid    <- residuals(fits[[sel_name]])
  sel_rstudent <- suppressWarnings(rstudent(fits[[sel_name]]))
  expect_true(sel_rstudent$z[4] > 4)
  expect_true(all(abs(sel_rstudent$z[-4]) < 4))
  expect_true(cor(sel_resid, sel_rstudent$z) > 0.9)

  glmm_resid    <- residuals(fits[[glmm_name]])
  glmm_rstudent <- suppressWarnings(rstudent(fits[[glmm_name]], type = "estimate"))
  expect_true(cor(glmm_resid, glmm_rstudent$z) > 0.9)
})

test_that("Residuals for BMA.norm fits are internally consistent", {

  name <- "dat.lehmann2018_BMA.norm"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  n        <- nobs(fit_brma)

  resid_marginal <- residuals(fit_brma, conditioning_depth = "marginal")
  resid_estimate <- residuals(fit_brma, conditioning_depth = "estimate")

  expect_residual_vector(resid_marginal, n)
  expect_residual_vector(resid_estimate, n)
  expect_error(residuals(fit_brma, conditioning_depth = "cluster"),
               "is only available for multilevel models")
  expect_true(cor(resid_marginal, resid_estimate) > 0.9)
  expect_true(sd(resid_marginal) > sd(resid_estimate))

  pearson_marginal <- residuals(
    fit_brma,
    type               = "pearson",
    conditioning_depth = "marginal"
  )
  pearson_estimate <- residuals(
    fit_brma,
    type               = "pearson",
    conditioning_depth = "estimate"
  )

  expect_residual_vector(pearson_marginal, n)
  expect_residual_vector(pearson_estimate, n)
  expect_true(cor(pearson_marginal, pearson_estimate) > 0.9)

  brma_rstandard <- rstandard(fit_brma)
  expect_residual_table(brma_rstandard, n)
  expect_equal(residuals(fit_brma, type = "rstandard"), brma_rstandard[["z"]],
               tolerance = 1e-12)

  brma_rstudent <- suppressWarnings(rstudent(fit_brma))
  expect_residual_table(brma_rstudent, n)
  expect_equal(suppressWarnings(residuals(fit_brma, type = "rstudent")),
               brma_rstudent[["z"]], tolerance = 1e-12)
  expect_error(
    residuals(fit_brma, type = "rstudent", conditioning_depth = "marginal"),
    "do not set 'conditioning_depth'"
  )
})

test_that("Residuals for BMA.glmm fits are internally consistent", {

  model_names <- c("bcg_BMA.glmm", "bcg_BMA.glmm_3lvl_location_scale")
  skip_if_missing_fits(model_names)

  fit_brma <- fits[["bcg_BMA.glmm"]]
  n        <- nobs(fit_brma)

  expect_residual_vector(residuals(fit_brma), n)
  expect_error(residuals(fit_brma, type = "pearson"), "normal outcome models")
  expect_error(rstandard(fit_brma), "normal outcome models")
  expect_residual_table(suppressWarnings(rstudent(fit_brma)), n)

  fit_3lvl <- fits[["bcg_BMA.glmm_3lvl_location_scale"]]
  n_3lvl   <- nrow(fit_3lvl[["data"]][["outcome"]])

  expect_null(fit_3lvl[["marglik"]])
  expect_residual_vector(residuals(fit_3lvl, conditioning_depth = "cluster"),
                         n_3lvl)
  expect_error(rstudent(fit_3lvl, unit = "cluster"), "Cluster-unit rstudent residuals")
})

test_that("Residuals for RoBMA fits are internally consistent", {

  model_names <- c("dat.lehmann2018_RoBMA", "dat.lehmann2018_RoBMA_3lvl_mods_scale")
  skip_if_missing_fits(model_names)

  fit_brma <- fits[["dat.lehmann2018_RoBMA"]]
  n        <- nobs(fit_brma)

  expect_residual_vector(residuals(fit_brma), n)
  expect_error(residuals(fit_brma, type = "pearson"), "selection models")
  expect_error(rstandard(fit_brma), "selection models")

  brma_rstudent <- suppressWarnings(rstudent(fit_brma))
  expect_residual_table(brma_rstudent, n)
  expect_equal(suppressWarnings(residuals(fit_brma, type = "rstudent")),
               brma_rstudent[["z"]], tolerance = 1e-12)

  fit_3lvl <- fits[["dat.lehmann2018_RoBMA_3lvl_mods_scale"]]
  n_3lvl   <- nrow(fit_3lvl[["data"]][["outcome"]])

  expect_null(fit_3lvl[["marglik"]])
  expect_residual_vector(residuals(fit_3lvl, conditioning_depth = "cluster"),
                         n_3lvl)
})
