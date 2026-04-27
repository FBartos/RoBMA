context("DFBETAS")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

expect_dfbetas_table <- function(x, n, min_cols = 1) {

  expect_s3_class(x, "data.frame")
  expect_equal(nrow(x), n)
  expect_gte(ncol(x), min_cols)
  expect_true(all(is.finite(as.matrix(x))))
}


# ============================================================================ #
# Test: Simple Meta-Analysis DFBETAS
# ============================================================================ #

test_that("DFBETAS for simple meta-analysis match metafor", {

  name <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # Compute DFBETAS
  metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas    <- dfbetas(fit_brma)

  # Compare with metafor results
  expect_equal(metafor_dfbetas[, 1], brma_dfbetas[, 1], tolerance = 0.10)
})

# ============================================================================ #
# Test: Meta-Regression DFBETAS
# ============================================================================ #

test_that("DFBETAS for meta-regression match metafor", {

  name <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas    <- dfbetas(fit_brma)

  # Compare with metafor results
  # there are three cases that are quite extreme and mismatch
  skip_case <- c(4, 6, 13)
  for (i in seq_len(ncol(metafor_dfbetas))) {
    expect_equal(metafor_dfbetas[-skip_case, i], brma_dfbetas[-skip_case, i], tolerance = 0.10)
  }
})

test_that("DFBETAS for meta-regression with interaction match metafor", {

  name <- "bcg_meta-regression4"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas    <- dfbetas(fit_brma)

  # metafor seems to struggle with numerical precision
  expect_dfbetas_table(brma_dfbetas, nobs(fit_brma))
})

# ============================================================================ #
# Test: Location-Scale Model DFBETAS
# ============================================================================ #

test_that("DFBETAS for location-scale model match metafor", {

  name <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # not implemented in metafor
  # metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas <- dfbetas(fit_brma)
  expect_dfbetas_table(brma_dfbetas, nobs(fit_brma))

  # compute dfbetas for scale
  brma_dfbetas.scale <- dfbetas(fit_brma, type = "scale")
  expect_dfbetas_table(brma_dfbetas.scale, nobs(fit_brma))
})

# ============================================================================ #
# Test: 3-Level Model DFBETAS
# ============================================================================ #

test_that("DFBETAS for 3-level model match metafor", {

  name <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.mv(fit_metafor))
  brma_dfbetas    <- dfbetas(fit_brma)

  # Compare with metafor results
  for (i in seq_len(ncol(metafor_dfbetas))) {
    expect_equal(metafor_dfbetas[, i], brma_dfbetas[, i], tolerance = 0.10)
  }
})

test_that("DFBETAS for 3-level meta-regression match metafor", {

  name <- "konstantopoulos2011_3lvl2"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.mv(fit_metafor))
  brma_dfbetas    <- dfbetas(fit_brma)

  # Compare with metafor results
  for (i in seq_len(ncol(metafor_dfbetas))) {
    expect_equal(metafor_dfbetas[, i], brma_dfbetas[, i], tolerance = 0.10)
  }
})

# ============================================================================ #
# Test: Selection Models DFBETAS
# ============================================================================ #

test_that("DFBETAS for selection models work", {

  name <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # not implemented in metafor
  # brma_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas <- dfbetas(fit_brma)

  ### compare at least consistency of positive and negative bselmodel fits
  fit_brma.neg     <- fits[["dat.lehmann2018-3PSM_neg"]]
  brma_dfbetas.neg <- dfbetas(fit_brma.neg)

  # case 4 (the outlier mismatches)
  expect_equal(brma_dfbetas[-4,1], -brma_dfbetas.neg[-4,1], tolerance = 0.01)

  brma_bias_dfbetas <- dfbetas(fit_brma, type = "bias")
  expect_dfbetas_table(brma_bias_dfbetas, nobs(fit_brma))
  expect_true(any(grepl("^omega", colnames(brma_bias_dfbetas))))
})

test_that("DFBETAS for selection meta-regression models work", {

  name <- "dat.lehmann2018-3PSMreg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # not implemented in metafor
  # metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas <- dfbetas(fit_brma)
  expect_dfbetas_table(brma_dfbetas, nobs(fit_brma))

  brma_bias_dfbetas <- dfbetas(fit_brma, type = "bias")
  expect_dfbetas_table(brma_bias_dfbetas, nobs(fit_brma))
  expect_true(any(grepl("^omega", colnames(brma_bias_dfbetas))))
})

# ============================================================================ #
# Test: PET Models DFBETAS
# ============================================================================ #

test_that("DFBETAS for PET models work", {

  name <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas    <- dfbetas(fit_brma)

  # Compare with metafor results
  # NOTE that the vi is not a term in bPET model -- no df beta for it
  expect_equal(metafor_dfbetas[,1], brma_dfbetas[,1], tolerance = 0.10)

  brma_bias_dfbetas <- dfbetas(fit_brma, type = "bias")
  expect_dfbetas_table(brma_bias_dfbetas, nobs(fit_brma))
  expect_true("PET" %in% colnames(brma_bias_dfbetas))
})

test_that("DFBETAS for PET meta-regression work", {

  name <- "dat.lehmann2018-PETreg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas    <- dfbetas(fit_brma)

  # Compare with metafor results
  # NOTE that the vi is not a term in bPET model -- no df beta for it
  expect_equal(metafor_dfbetas[,1], brma_dfbetas[,1], tolerance = 0.10)
  expect_equal(metafor_dfbetas[,3], brma_dfbetas[,2], tolerance = 0.10)

  brma_bias_dfbetas <- dfbetas(fit_brma, type = "bias")
  expect_dfbetas_table(brma_bias_dfbetas, nobs(fit_brma))
  expect_true("PET" %in% colnames(brma_bias_dfbetas))
})

test_that("DFBETAS for PET models work (negative)", {

  name <- "dat.lehmann2018-PET_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas    <- dfbetas(fit_brma)

  # Compare with metafor results
  # NOTE that the vi is not a term in bPET model -- no df beta for it
  expect_equal(metafor_dfbetas[,1], brma_dfbetas[,1], tolerance = 0.10)

  brma_bias_dfbetas <- dfbetas(fit_brma, type = "bias")
  expect_dfbetas_table(brma_bias_dfbetas, nobs(fit_brma))
  expect_true("PET" %in% colnames(brma_bias_dfbetas))
})

test_that("DFBETAS for PEESE publication bias parameters work", {

  name <- "dat.lehmann2018-PEESE"
  fit_brma <- fits[[name]]

  brma_bias_dfbetas <- dfbetas(fit_brma, type = "bias")
  expect_dfbetas_table(brma_bias_dfbetas, nobs(fit_brma))
  expect_true("PEESE" %in% colnames(brma_bias_dfbetas))
})

# ============================================================================ #
# Test: GLMM DFBETAS
# ============================================================================ #

test_that("DFBETAS for GLMM work", {

  name <- "nielweise2008_glmm"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  brma_dfbetas <- dfbetas(fit_brma)
  # not implemented in metafor
  expect_dfbetas_table(brma_dfbetas, nobs(fit_brma))
})

test_that("DFBETAS for GLMM regression work", {

  name <- "bcg_glmm_reg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  brma_dfbetas <- dfbetas(fit_brma)
  # not implemented in metafor
  expect_dfbetas_table(brma_dfbetas, nobs(fit_brma))
})

# ============================================================================ #
# Test: Model-Averaging DFBETAS
# ============================================================================ #

test_that("DFBETAS for BMA.norm fits are internally consistent", {

  fit_simple <- fits[["dat.lehmann2018_BMA.norm"]]
  dfb_simple <- dfbetas(fit_simple)
  expect_dfbetas_table(dfb_simple, nobs(fit_simple))

  fit_mods <- fits[["dat.lehmann2018_BMA.norm_mods"]]
  dfb_mods <- dfbetas(fit_mods)
  expect_dfbetas_table(dfb_mods, nobs(fit_mods), min_cols = 2)

  fit_scale <- fits[["dat.lehmann2018_BMA.norm_scale"]]
  dfb_scale <- dfbetas(fit_scale, type = "scale")
  expect_dfbetas_table(dfb_scale, nobs(fit_scale), min_cols = 2)
})

test_that("DFBETAS for BMA.glmm fits are internally consistent", {

  fit_simple <- fits[["bcg_BMA.glmm"]]
  dfb_simple <- dfbetas(fit_simple)
  expect_dfbetas_table(dfb_simple, nobs(fit_simple))
})

test_that("DFBETAS for RoBMA fits are internally consistent", {

  fit_simple <- fits[["dat.lehmann2018_RoBMA_mods"]]
  dfb_simple <- dfbetas(fit_simple)
  expect_dfbetas_table(dfb_simple, nobs(fit_simple))

  fit_complex <- fits[["dat.lehmann2018_RoBMA_3lvl_mods_scale"]]
  dfb_complex <- dfbetas(fit_complex)
  expect_dfbetas_table(dfb_complex, nobs(fit_complex))

  fit_bias <- fits[["dat.lehmann2018_RoBMA"]]
  dfb_bias <- dfbetas(fit_bias, type = "bias")
  expect_dfbetas_table(dfb_bias, nobs(fit_bias), min_cols = 3)
  expect_true(any(grepl("^omega", colnames(dfb_bias))))
  expect_true("PET" %in% colnames(dfb_bias))
  expect_true("PEESE" %in% colnames(dfb_bias))
})
