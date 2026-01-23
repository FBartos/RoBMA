context("Hatvalues")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "hatvalues")

# list & load all fits
skip_if_no_fits()
skip_if_not_installed("metafor")
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


# ============================================================================ #
# Test: Simple Meta-Analysis Hatvalues
# ============================================================================ #

test_that("Hatvalues for simple meta-analysis match metafor", {

    name <- "bcg_meta-analysis"
    fit_metafor <- info[[name]][["metafor"]]
    fit_brma    <- fits[[name]]

    metafor_hat <- hatvalues(fit_metafor)
    brma_hat    <- hatvalues(fit_brma)

    expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor")
})


# ============================================================================ #
# Test: Meta-Regression Hatvalues
# ============================================================================ #

test_that("Hatvalues for meta-regression match metafor", {

    name <- "bcg_meta-regression"
    fit_metafor <- info[[name]][["metafor"]]
    fit_brma    <- fits[[name]]

    metafor_hat <- hatvalues(fit_metafor)
    brma_hat    <- hatvalues(fit_brma)

    expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for meta-regression")
})


# ============================================================================ #
# Test: 3-Level Model Hatvalues
# ============================================================================ #

test_that("Hatvalues for 3-level model match metafor", {

    name <- "konstantopoulos2011_3lvl"
    fit_metafor <- info[[name]][["metafor"]]
    fit_brma    <- fits[[name]]

    metafor_hat <- hatvalues(fit_metafor)
    brma_hat    <- hatvalues(fit_brma)

    expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for 3-level model")
})


# ============================================================================ #
# Test: PET Models Hatvalues
# ============================================================================ #

test_that("Hatvalues for PET models work", {

    name <- "dat.lehmann2018-PET"
    fit_metafor <- info[[name]][["metafor"]]
    fit_brma    <- fits[[name]]

    metafor_hat <- hatvalues(fit_metafor)
    brma_hat    <- hatvalues(fit_brma)

    expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for PET model")
})

# ============================================================================ #
# Test: PET Models Hatvalues (negative)
# ============================================================================ #

test_that("Hatvalues for PET models work  (negative)", {

  name <- "dat.lehmann2018-PET_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_hat <- hatvalues(fit_metafor)
  brma_hat    <- hatvalues(fit_brma)

  expect_equal(brma_hat, as.vector(metafor_hat), tolerance = 0.05, info = "brma hatvalues should match metafor for PET model")
})


# ============================================================================ #
# Test: glmm/bselmod Hatvalues are not available
# ============================================================================ #

test_that("Hatvalues should error for other models", {

    name <- "nielweise2008_glmm"
    fit_brma <- fits[[name]]

    expect_error(hatvalues(fit_brma), "available")

    name <- "dat.lehmann2018-3PSM"
    fit_brma <- fits[[name]]

    expect_error(hatvalues(fit_brma), "available")
})
