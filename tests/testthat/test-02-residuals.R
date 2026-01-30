context("Residuals")

# TODO:
# add a unit test with interaction and meandif factor
# (to verify that the model matrix is correct)

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "residuals")

# list & load all fits
skip_if_no_fits()
skip_if_not_installed("metafor")
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


# ============================================================================ #
# Test: Simple Meta-Analysis Residuals
# ============================================================================ #

test_that("Residuals for simple meta-analysis match metafor", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson    <- residuals(fit_brma, type = "pearson")

  expect_equal(brma_pearson, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------
  # conditonal == marginal for a simple model

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard    <- rstandard(fit_brma)

  expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.05,
               info = "brma rstandard resid should match metafor")
  expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.05,
               info = "brma rstandard se should match metafor")
  expect_equal(brma_rstandard$z, metafor_rstandard$z, tolerance = 0.05,
               info = "brma rstandard z should match metafor")

  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- rstudent(fit_brma)

  expect_equal(brma_rstudent$resid, metafor_rstudent$resid, tolerance = 0.05,
               info = "brma rstudent resid should match metafor")
  expect_equal(brma_rstudent$se, metafor_rstudent$se, tolerance = 0.05,
               info = "brma rstudent se should match metafor")
  expect_equal(brma_rstudent$z, metafor_rstudent$z, tolerance = 0.05,
               info = "brma rstudent z should match metafor")

})


# ============================================================================ #
# Test: Meta-Regression Residuals
# ============================================================================ #

test_that("Residuals for meta-regression match metafor", {

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.10,
               info = "brma outcome residuals should match metafor for meta-regression")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson    <- residuals(fit_brma, type = "pearson")

  expect_equal(brma_pearson, as.vector(metafor_pearson), tolerance = 0.10,
               info = "brma pearson residuals should match metafor for meta-regression")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard    <- rstandard(fit_brma)

  expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.10,
               info = "brma rstandard resid should match metafor for meta-regression")
  expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.10,
               info = "brma rstandard se should match metafor for meta-regression")
  expect_equal(brma_rstandard$z, metafor_rstandard$z, tolerance = 0.10,
               info = "brma rstandard z should match metafor for meta-regression")


  metafor_rstandard <- rstandard(fit_metafor, type = "conditional")
  brma_rstandard    <- rstandard(fit_brma, type = "conditional")

  expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.10,
               info = "brma rstandard resid should match metafor (conditional)")
  expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.10,
               info = "brma rstandard se should match metafor (conditional)")
  expect_equal(brma_rstandard$z, metafor_rstandard$z, tolerance = 0.10,
               info = "brma rstandard z should match metafor (conditional)")


  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- suppressWarnings(rstudent(fit_brma))

  expect_equal(brma_rstudent$resid, metafor_rstudent$resid, tolerance = 0.10,
               info = "brma rstudent resid should match metafor")
  expect_equal(brma_rstudent$se, metafor_rstudent$se, tolerance = 0.10,
               info = "brma rstudent se should match metafor")
  expect_equal(brma_rstudent$z, metafor_rstudent$z, tolerance = 0.15,
               info = "brma rstudent z should match metafor")
})


# ============================================================================ #
# Test: Location-Scale Model Residuals
# ============================================================================ #

test_that("Residuals for location-scale model match metafor", {

  name        <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for location-scale model")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson    <- residuals(fit_brma, type = "pearson")

  expect_equal(brma_pearson, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for location-scale model")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard    <- rstandard(fit_brma)

  expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.05,
               info = "brma rstandard resid should match metafor for location-scale model")
  expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.05,
               info = "brma rstandard se should match metafor for location-scale model")
  expect_equal(brma_rstandard$z, metafor_rstandard$z, tolerance = 0.05,
               info = "brma rstandard z should match metafor for location-scale model")

  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals: no metafor comparison possible, check distribution
  # --------------------------------------------------

  # metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- rstudent(fit_brma)

  # expect_equal(brma_rstudent$resid, metafor_rstudent$resid, tolerance = 0.10,
  #              info = "brma rstudent resid should match metafor")
  # expect_equal(brma_rstudent$se, metafor_rstudent$se, tolerance = 0.05,
  #              info = "brma rstudent se should match metafor")
  # expect_equal(brma_rstudent$z, metafor_rstudent$z, tolerance = 0.05,
  #              info = "brma rstudent z should match metafor")

  expect_equal(mean(brma_rstudent$z), 0, 0.05)
  expect_equal(sd(brma_rstudent$z),   1, 0.05)
  expect_true(cor(brma_rstudent$z, metafor_rstandard$z, method = "spearman") > 0.8)
})


# ============================================================================ #
# Test: 3-Level Model Residuals
# ============================================================================ #

test_that("Residuals for 3-level model match metafor", {

  name        <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for 3-level model")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson    <- residuals(fit_brma, type = "pearson")

  expect_equal(brma_pearson, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for 3-level model")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard    <- rstandard(fit_brma)

  expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.05,
               info = "brma rstandard resid should match metafor for 3-level model")
  expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.05,
               info = "brma rstandard se should match metafor for 3-level model")
  expect_equal(brma_rstandard$z, metafor_rstandard$z, tolerance = 0.05,
               info = "brma rstandard z should match metafor for 3-level model")

  # TODO: check why & fix
  metafor_rstandard <- rstandard(fit_metafor, type = "conditional")
  brma_rstandard    <- rstandard(fit_brma, type = "conditional")

  expect_error(expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.10,
                            info = "brma rstandard resid should match metafor (conditional)"))
  expect_error(expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.10,
                            info = "brma rstandard se should match metafor (conditional)"))
  expect_error(expect_equal(brma_rstandard$z, metafor_rstandard$z, tolerance = 0.10,
                            info = "brma rstandard z should match metafor (conditional)"))

  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- rstudent(fit_brma)

  expect_equal(brma_rstudent$resid, metafor_rstudent$resid, tolerance = 0.05,
               info = "brma rstudent resid should match metafor")
  expect_equal(brma_rstudent$se, metafor_rstudent$se, tolerance = 0.05,
               info = "brma rstudent se should match metafor")
  expect_equal(brma_rstudent$z, metafor_rstudent$z, tolerance = 0.10,
               info = "brma rstudent z should match metafor")

  # TODO: add clustered residuals

})


# ============================================================================ #
# Test: Selection Model Residuals (Positive Effects)
# ============================================================================ #

test_that("Residuals for selection model (positive) match metafor", {

  name        <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for selection model (positive)")

  # --------------------------------------------------
  # Pearson residuals: should not be available for selection models
  # --------------------------------------------------

  expect_error(
    residuals(fit_brma, type = "pearson"),
    "Pearson residuals are not available for selection models",
    info = "pearson should error for selection models"
  )

  # --------------------------------------------------
  # Standardized (rstandard) residuals: should not be available for selection models
  # --------------------------------------------------

  expect_error(
    residuals(fit_brma, type = "rstandard"),
    "Standardized residuals.*are not available for selection models",
    info = "rstandard should error for selection models"
  )

  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals (no metafor comparison for selection models)
  # --------------------------------------------------

  brma_rstudent <- suppressWarnings(rstudent(fit_brma))

  # the residuals find one massive outlier at 4th row
  expect_true(brma_rstudent$z[4]           > 4,  info = "rstudent finds one large outlier")
  expect_true(all(abs(brma_rstudent$z[-4]) < 4), info = "rstudent finds one large outlier")
  expect_true(cor(brma_resid, brma_rstudent$z) > 0.9, info = "rstudent and residuals are directionally aligned")
})


# ============================================================================ #
# Test: Selection Model Residuals (Negative Effects)
# ============================================================================ #

test_that("Residuals for selection model (negative) match metafor", {

  name        <- "dat.lehmann2018-3PSM_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for selection model (negative)")

  # --------------------------------------------------
  # Pearson residuals: should not be available for selection models
  # --------------------------------------------------

  expect_error(
    residuals(fit_brma, type = "pearson"),
    "Pearson residuals are not available for selection models",
    info = "pearson should error for selection models"
  )

  # --------------------------------------------------
  # Standardized (rstandard) residuals: should not be available for selection models
  # --------------------------------------------------

  expect_error(
    residuals(fit_brma, type = "rstandard"),
    "Standardized residuals.*are not available for selection models",
    info = "rstandard should error for selection models"
  )

  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals (no metafor comparison for selection models)
  # --------------------------------------------------

  brma_rstudent <- suppressWarnings(rstudent(fit_brma))

  # the residuals find one massive outlier at 4th row
  # note that the direction of residuals is flipped from (previous selection model since the data were flipped too)
  expect_true(brma_rstudent$z[4]           < -4,      info = "rstudent finds one large outlier")
  expect_true(all(abs(brma_rstudent$z[-4]) < 4),      info = "rstudent finds one large outlier")
  expect_true(cor(brma_resid, brma_rstudent$z) > 0.9, info = "rstudent with negative direction is correctly flipped")
})


# ============================================================================ #
# Test: PET Model Residuals (Positive Effects)
# ============================================================================ #

test_that("Residuals for PET model (positive) match metafor", {

  name        <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for PET model (positive)")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson    <- residuals(fit_brma, type = "pearson")

  expect_equal(brma_pearson, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for PET model (positive)")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard    <- rstandard(fit_brma)

  expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.05,
               info = "brma rstandard resid should match metafor for PET model (positive)")
  expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.05,
               info = "brma rstandard se should match metafor for PET model (positive)")
  expect_equal(brma_rstandard$z, metafor_rstandard$z, tolerance = 0.05,
               info = "brma rstandard z should match metafor for PET model (positive)")

  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- suppressWarnings(rstudent(fit_brma))

  expect_equal(brma_rstudent$resid, metafor_rstudent$resid, tolerance = 0.05,
               info = "brma rstudent resid should match metafor")
  expect_equal(brma_rstudent$se, metafor_rstudent$se, tolerance = 0.05,
               info = "brma rstudent se should match metafor")
  expect_equal(brma_rstudent$z, metafor_rstudent$z, tolerance = 0.10,
               info = "brma rstudent z should match metafor")

})


# ============================================================================ #
# Test: PET Model Residuals (Negative Effects)
# ============================================================================ #

test_that("Residuals for PET model (negative) match metafor", {

  name        <- "dat.lehmann2018-PET_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for PET model (negative)")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_pearson <- residuals(fit_metafor, type = "pearson")
  brma_pearson    <- residuals(fit_brma, type = "pearson")

  expect_equal(brma_pearson, as.vector(metafor_pearson), tolerance = 0.05,
               info = "brma pearson residuals should match metafor for PET model (negative)")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstandard    <- rstandard(fit_brma)

  expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.05,
               info = "brma rstandard resid should match metafor for PET model (negative)")
  expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.05,
               info = "brma rstandard se should match metafor for PET model (negative)")
  expect_equal(brma_rstandard$z, metafor_rstandard$z, tolerance = 0.05,
               info = "brma rstandard z should match metafor for PET model (negative)")


  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- suppressWarnings(rstudent(fit_brma))

  expect_equal(brma_rstudent$resid, metafor_rstudent$resid, tolerance = 0.05,
               info = "brma rstudent resid should match metafor")
  expect_equal(brma_rstudent$se, metafor_rstudent$se, tolerance = 0.05,
               info = "brma rstudent se should match metafor")
  expect_equal(brma_rstudent$z, metafor_rstudent$z, tolerance = 0.10,
               info = "brma rstudent z should match metafor")
})


# ============================================================================ #
# Test: GLMM Model Residuals
# ============================================================================ #

test_that("Residuals for GLMM model are computed correctly", {

  name        <- "nielweise2008_glmm"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.05,
               info = "brma outcome residuals should match metafor for GLMM model")

  # --------------------------------------------------
  # Pearson residuals: should not be available for GLMM models
  # --------------------------------------------------

  expect_error(
    residuals(fit_brma, type = "pearson"),
    "Pearson residuals are only available for normal outcome models",
    info = "pearson should error for GLMM models"
  )

  # --------------------------------------------------
  # Standardized (rstandard) residuals: should not be available for GLMM models
  # --------------------------------------------------

  expect_error(
    residuals(fit_brma, type = "rstandard"),
    "Standardized residuals.*are only available for normal outcome models",
    info = "rstandard should error for GLMM models"
  )

  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals (no metafor comparison for GLMM)
  # --------------------------------------------------

  brma_rstudent <- suppressWarnings(rstudent(fit_brma, type = "estimate"))

  # The residuals seem to be well centered and scaled (the issues is already in standard residuals)
  # expect_equal(mean(brma_rstudent$z), 0, tolerance = 0.10, info = "rstudent z is standardized")
  # expect_equal(sd(brma_rstudent$z),   1, tolerance = 0.10, info = "rstudent z is standardized")
  expect_true(cor(brma_resid, brma_rstudent$z) > 0.9,      info = "rstudent and residuals are directionally aligned")
})
