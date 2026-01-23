context("DFBETAS")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

# list & load all fits
skip_if_no_fits()
skip_if_not_installed("metafor")
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


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
  expect_s3_class(brma_dfbetas, "data.frame")

  # compute dfbetas for scale
  brma_dfbetas.scale <- dfbetas(fit_brma, type = "scale")
  expect_s3_class(brma_dfbetas.scale, "data.frame")
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
  expect_equal(metafor_dfbetas[,1], brma_dfbetas[,1], tolerance = 0.10)

  # TODO:
  # Figure out how to deal with bias coefficients (add type arugment)
})

# ============================================================================ #
# Test: PET Models DFBETAS (negative)
# ============================================================================ #

test_that("DFBETAS for PET models work (negative)", {

  name <- "dat.lehmann2018-PET_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_dfbetas <- as.data.frame(metafor::dfbetas.rma.uni(fit_metafor))
  brma_dfbetas    <- dfbetas(fit_brma)

  # Compare with metafor results
  expect_equal(metafor_dfbetas[,1], brma_dfbetas[,1], tolerance = 0.10)

  # TODO:
  # Figure out how to deal with bias coefficients (add type argument)
})



# ============================================================================ #
# Test: GLMM DFBETAS
# ============================================================================ #

test_that("DFBETAS for GLMM work", {

  name <- "nielweise2008_glmm"
  fit_brma <- fits[[name]]

  brma_dfbetas <- dfbetas(fit_brma)
  expect_s3_class(brma_dfbetas, "data.frame")
})
