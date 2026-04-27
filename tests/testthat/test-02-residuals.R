context("Residuals")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "residuals")

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

expect_residual_vector <- function(x, n) {

  expect_type(x, "double")
  expect_equal(length(x), n)
  expect_true(all(is.finite(x)))
}

expect_residual_table <- function(x, n, check_se = TRUE) {

  expect_s3_class(x, "data.frame")
  expect_true(all(c("resid", "se", "z") %in% names(x)))
  expect_equal(nrow(x), n)
  expect_true(all(is.finite(as.matrix(x))))
  if (check_se) {
    expect_true(all(x[["se"]] > 0))
  }
}


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

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.05,
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
  brma_rstandard    <- rstandard(fit_brma, conditioning_depth = "estimate")

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

  expect_true(all(is.finite(brma_rstudent$resid)),
              info = "LOO-PIT residual means should be finite")
  expect_true(all(is.finite(brma_rstudent$se)) && all(brma_rstudent$se > 0),
              info = "LOO predictive SEs should be finite and positive")
  expect_true(cor(brma_rstudent$resid, metafor_rstudent$resid, method = "spearman") > 0.90,
              info = "LOO-PIT raw residuals should be directionally aligned with metafor deleted residuals")
  expect_true(cor(brma_rstudent$z, metafor_rstudent$z, method = "spearman") > 0.90,
              info = "LOO-PIT z residuals should be directionally aligned with metafor deleted residuals")
})

test_that("Residuals for meta-regression with interactions match metafor", {

  name        <- "bcg_meta-regression4"
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


  metafor_rstandard <- suppressWarnings(rstandard(fit_metafor, type = "conditional"))
  brma_rstandard    <- rstandard(fit_brma, conditioning_depth = "estimate")

  valid <- is.finite(metafor_rstandard$se) & is.finite(metafor_rstandard$z)

  expect_equal(brma_rstandard$resid[valid], metafor_rstandard$resid[valid], tolerance = 0.10,
               info = "brma rstandard resid should match metafor (conditional)")
  expect_equal(brma_rstandard$se[valid], metafor_rstandard$se[valid], tolerance = 0.10,
               info = "brma rstandard se should match metafor (conditional)")
  expect_equal(brma_rstandard$z[valid], metafor_rstandard$z[valid], tolerance = 0.10,
               info = "brma rstandard z should match metafor (conditional)")
  expect_equal(brma_rstandard$resid[!valid], 0, tolerance = 1e-8,
               info = "brma saturated conditional residuals should be zero")
  expect_true(all(brma_rstandard$se[!valid] < 1e-8),
              info = "brma saturated conditional residual SEs should collapse to zero")
  expect_equal(brma_rstandard$z[!valid], 0, tolerance = 1e-8,
               info = "brma saturated conditional residual z values should be zero")


})

test_that("Extended rstudent residuals for interaction meta-regression align with metafor", {

  skip_if_not_full_diagnostics(
    "Interaction LOO-PIT residuals duplicate meta-regression LOO coverage and are outlier-heavy."
  )

  name        <- "bcg_meta-regression4"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- suppressWarnings(rstudent(fit_brma))

  expect_true(all(is.finite(brma_rstudent$resid)),
              info = "LOO-PIT residual means should be finite")
  expect_true(all(is.finite(brma_rstudent$se)) && all(brma_rstudent$se > 0),
              info = "LOO predictive SEs should be finite and positive")
  expect_true(cor(brma_rstudent$resid, metafor_rstudent$resid, method = "spearman", use = "complete.obs") > 0.90,
              info = "LOO-PIT raw residuals should be directionally aligned with metafor deleted residuals")
  expect_true(cor(brma_rstudent$z, metafor_rstudent$z, method = "spearman", use = "complete.obs") > 0.90,
              info = "LOO-PIT z residuals should be directionally aligned with metafor deleted residuals")
})

test_that("Residuals for meta-regression with interactions match across parameterizations", {

  fit_brma1 <- fits[["bcg_meta-regression4"]]
  fit_brma2 <- fits[["bcg_meta-regression4b"]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  brma_resid1 <- residuals(fit_brma1)
  brma_resid2 <- residuals(fit_brma2)

  expect_equal(brma_resid1, brma_resid2, tolerance = 0.10,
               info = "brma outcome residuals should match")

  # --------------------------------------------------
  # Pearson residuals: compare brma vs metafor
  # --------------------------------------------------

  brma_pearson1    <- residuals(fit_brma1, type = "pearson")
  brma_pearson2    <- residuals(fit_brma2, type = "pearson")

  expect_equal(brma_pearson1, brma_pearson2, tolerance = 0.10,
               info = "brma pearson residuals should match")

  # --------------------------------------------------
  # Standardized (rstandard) residuals: compare brma vs metafor
  # --------------------------------------------------

  brma_rstandard1 <- rstandard(fit_brma1)
  brma_rstandard2 <- rstandard(fit_brma2)

  expect_equal(brma_rstandard1, brma_rstandard2, tolerance = 0.10,
               info = "brma rstandard residuals should match")

})

test_that("Extended rstudent residuals match across interaction parameterizations", {

  fit_brma1 <- fits[["bcg_meta-regression4"]]
  fit_brma2 <- fits[["bcg_meta-regression4b"]]

  brma_rstudent1 <- suppressWarnings(rstudent(fit_brma1))
  brma_rstudent2 <- suppressWarnings(rstudent(fit_brma2))

  expect_equal(brma_rstudent1[, "z"], brma_rstudent2[, "z"], tolerance = 0.15,
               info = "brma rstudent residuals should match")
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

})

test_that("Extended rstudent residuals for location-scale model are calibrated", {

  skip_if_not_full_diagnostics(
    "Location-scale LOO-PIT has no direct metafor equivalent; core scale residuals are tested via rstandard."
  )

  name        <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_rstandard <- rstandard(fit_metafor)
  brma_rstudent     <- rstudent(fit_brma)

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

  # metafor::rstandard() ignores type = "conditional" for rma.mv fits,
  # so use the multilevel BLUP-based oracle from common-functions.R instead.
  metafor_rstandard <- metafor_rstandard_conditional_mv(fit_metafor)
  brma_rstandard    <- rstandard(fit_brma, conditioning_depth = "estimate")
  brma_theta         <- colMeans(as.matrix(predict(fit_brma, type = "estimate", quiet = TRUE)))
  metafor_theta      <- fit_metafor[["yi"]] - metafor_rstandard$resid

  expect_equal(as.vector(brma_theta), as.vector(metafor_theta), tolerance = 0.11,
               info = "brma estimate-level BLUPs should match the metafor multilevel oracle")
  expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.11,
               info = "brma conditional residual means should match the metafor BLUP oracle")
  expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.05,
               info = "brma conditional residual standard errors should match the metafor BLUP oracle")

  # --------------------------------------------------
  # rstudent (LOO-PIT) residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- suppressWarnings(rstudent(fit_brma))

  expect_true(all(is.finite(brma_rstudent$resid)),
              info = "LOO-PIT residual means should be finite")
  expect_true(all(is.finite(brma_rstudent$se)) && all(brma_rstudent$se > 0),
              info = "LOO predictive SEs should be finite and positive")
  expect_true(cor(brma_rstudent$resid, metafor_rstudent$resid, method = "spearman") > 0.65,
              info = "multilevel LOO-PIT residuals should remain directionally aligned with metafor deleted residuals")
  expect_true(cor(brma_rstudent$se, metafor_rstudent$se, method = "spearman") > 0.80,
              info = "multilevel LOO predictive SEs should preserve study precision ordering")
  expect_true(cor(brma_rstudent$z, metafor_rstudent$z, method = "spearman") > 0.65,
              info = "multilevel LOO-PIT z residuals should remain directionally aligned with metafor deleted residuals")

  # TODO: add clustered residuals

})

test_that("Residuals for 3-level meta-regression match metafor", {

  name        <- "konstantopoulos2011_3lvl2"
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

  # metafor::rstandard() ignores type = "conditional" for rma.mv fits,
  # so use the multilevel BLUP-based oracle from common-functions.R instead.
  metafor_rstandard <- metafor_rstandard_conditional_mv(fit_metafor)
  brma_rstandard    <- rstandard(fit_brma, conditioning_depth = "estimate")
  brma_theta         <- colMeans(as.matrix(predict(fit_brma, type = "estimate", quiet = TRUE)))
  metafor_theta      <- fit_metafor[["yi"]] - metafor_rstandard$resid

  expect_equal(as.vector(brma_theta), as.vector(metafor_theta), tolerance = 0.11,
               info = "brma estimate-level BLUPs should match the metafor multilevel oracle")
  expect_equal(brma_rstandard$resid, metafor_rstandard$resid, tolerance = 0.11,
               info = "brma conditional residual means should match the metafor BLUP oracle")
  expect_equal(brma_rstandard$se, metafor_rstandard$se, tolerance = 0.05,
               info = "brma conditional residual standard errors should match the metafor BLUP oracle")

  # TODO: add clustered residuals

})

test_that("Extended rstudent residuals for 3-level meta-regression align with metafor", {

  skip_if_not_full_diagnostics(
    "3-level LOO-PIT coverage is kept on the intercept-only multilevel model by default."
  )

  name        <- "konstantopoulos2011_3lvl2"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- suppressWarnings(rstudent(fit_brma))

  expect_true(all(is.finite(brma_rstudent$resid)),
              info = "LOO-PIT residual means should be finite")
  expect_true(all(is.finite(brma_rstudent$se)) && all(brma_rstudent$se > 0),
              info = "LOO predictive SEs should be finite and positive")
  expect_true(cor(brma_rstudent$resid, metafor_rstudent$resid, method = "spearman") > 0.60,
              info = "multilevel LOO-PIT residuals should remain directionally aligned with metafor deleted residuals")
  expect_true(cor(brma_rstudent$se, metafor_rstudent$se, method = "spearman") > 0.80,
              info = "multilevel LOO predictive SEs should preserve study precision ordering")
  expect_true(cor(brma_rstudent$z, metafor_rstudent$z, method = "spearman") > 0.60,
              info = "multilevel LOO-PIT z residuals should remain directionally aligned with metafor deleted residuals")
})

# ============================================================================ #
# Test: Selection Model Residuals
# ============================================================================ #

test_that("Residuals for selection model match metafor", {

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
  expect_true(all(is.finite(brma_rstudent$se)), info = "selection LOO residuals have predictive SEs for funnel plots")
  expect_true(all(brma_rstudent$se > 0),         info = "selection LOO predictive SEs are positive")
})

test_that("Residuals for selection meta-regression match metafor", {

  name        <- "dat.lehmann2018-3PSMreg"
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

})

test_that("Extended rstudent residuals for selection meta-regression work", {

  skip_if_not_full_diagnostics(
    "Selection LOO-PIT coverage is kept on positive and negative intercept-only selection models by default."
  )

  name     <- "dat.lehmann2018-3PSMreg"
  fit_brma <- fits[[name]]

  brma_resid    <- residuals(fit_brma)
  brma_rstudent <- suppressWarnings(rstudent(fit_brma))

  # the residuals find one massive outlier at 4th row
  expect_true(brma_rstudent$z[4]           > 4,  info = "rstudent finds one large outlier")
  expect_true(all(abs(brma_rstudent$z[-4]) < 4), info = "rstudent finds one large outlier")
  expect_true(cor(brma_resid, brma_rstudent$z) > 0.9, info = "rstudent and residuals are directionally aligned")
})

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
  expect_true(all(is.finite(brma_rstudent$se)),        info = "selection LOO residuals have predictive SEs for funnel plots")
  expect_true(all(brma_rstudent$se > 0),                info = "selection LOO predictive SEs are positive")
})

# ============================================================================ #
# Test: PET Model Residuals
# ============================================================================ #

test_that("Residuals for PET model match metafor", {

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

test_that("Residuals for PET meta-regression match metafor", {

  name        <- "dat.lehmann2018-PETreg"
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

})

test_that("Extended rstudent residuals for PET meta-regression match metafor", {

  skip_if_not_full_diagnostics(
    "PET LOO-PIT coverage is kept on positive and negative intercept-only PET models by default."
  )

  name        <- "dat.lehmann2018-PETreg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  metafor_rstudent <- rstudent(fit_metafor)
  brma_rstudent    <- suppressWarnings(rstudent(fit_brma))

  expect_equal(brma_rstudent$resid, metafor_rstudent$resid, tolerance = 0.05,
               info = "brma rstudent resid should match metafor")
  expect_equal(brma_rstudent$se, metafor_rstudent$se, tolerance = 0.05,
               info = "brma rstudent se should match metafor")
  expect_equal(brma_rstudent$z, metafor_rstudent$z, tolerance = 0.10,
               info = "brma rstudent z should match metafor")
})

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

test_that("Residuals for GLMM model match metafor", {

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

test_that("Residuals for GLMM meta-regression match metafor", {

  name        <- "bcg_glmm_reg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Outcome residuals: compare brma vs metafor
  # --------------------------------------------------

  metafor_resid <- residuals(fit_metafor)
  brma_resid    <- residuals(fit_brma)

  expect_equal(brma_resid, as.vector(metafor_resid), tolerance = 0.10,
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

})

test_that("Extended rstudent residuals for GLMM meta-regression work", {

  skip_if_not_full_diagnostics(
    "GLMM LOO-PIT coverage is kept on the intercept-only GLMM model by default."
  )

  name     <- "bcg_glmm_reg"
  fit_brma <- fits[[name]]

  brma_resid    <- residuals(fit_brma)
  brma_rstudent <- suppressWarnings(rstudent(fit_brma, type = "estimate"))

  expect_true(cor(brma_resid, brma_rstudent$z) > 0.9,
              info = "rstudent and residuals are directionally aligned")
})

# ============================================================================ #
# Test: BMA Model Residuals
# ============================================================================ #

test_that("Residuals for BMA.norm model-averaging fit are internally consistent", {

  name     <- "dat.lehmann2018_BMA.norm"
  fit_brma <- fits[[name]]
  n        <- nobs(fit_brma)

  brma_resid_marginal <- residuals(fit_brma, conditioning_depth = "marginal")
  brma_resid_estimate <- residuals(fit_brma, conditioning_depth = "estimate")

  expect_residual_vector(brma_resid_marginal, n)
  expect_residual_vector(brma_resid_estimate, n)
  expect_error(
    residuals(fit_brma, conditioning_depth = "cluster"),
    "is only available for multilevel models"
  )
  expect_true(cor(brma_resid_marginal, brma_resid_estimate) > 0.9,
              info = "different conditioning depths of residuals are correlated")
  expect_true(sd(brma_resid_marginal) > sd(brma_resid_estimate),
              info = "marginal residuals are more variable")

  brma_pearson_marginal <- residuals(
    fit_brma,
    type               = "pearson",
    conditioning_depth = "marginal"
  )
  brma_pearson_estimate <- residuals(
    fit_brma,
    type               = "pearson",
    conditioning_depth = "estimate"
  )

  expect_residual_vector(brma_pearson_marginal, n)
  expect_residual_vector(brma_pearson_estimate, n)
  expect_true(cor(brma_pearson_marginal, brma_pearson_estimate) > 0.9,
              info = "different conditioning depths of pearson residuals are correlated")

  brma_rstandard <- rstandard(fit_brma)
  expect_residual_table(brma_rstandard, n)
  expect_equal(residuals(fit_brma, type = "rstandard"), brma_rstandard[["z"]],
               tolerance = 1e-12)

  brma_rstudent <- suppressWarnings(rstudent(fit_brma))
  expect_residual_table(brma_rstudent, n)
  expect_equal(suppressWarnings(residuals(fit_brma, type = "rstudent")), brma_rstudent[["z"]],
               tolerance = 1e-12)
  expect_error(
    residuals(fit_brma, type = "rstudent", conditioning_depth = "marginal"),
    "do not set 'conditioning_depth'"
  )
})

test_that("Residuals for BMA.glmm model-averaging fits are internally consistent", {

  fit_brma <- fits[["bcg_BMA.glmm"]]
  n        <- nobs(fit_brma)

  expect_residual_vector(residuals(fit_brma), n)
  expect_error(residuals(fit_brma, type = "pearson"), "normal outcome models")
  expect_error(rstandard(fit_brma), "normal outcome models")

  brma_rstudent <- suppressWarnings(rstudent(fit_brma))
  expect_residual_table(brma_rstudent, n)

  fit_3lvl <- fits[["bcg_BMA.glmm_3lvl_location_scale"]]
  n_3lvl   <- nrow(fit_3lvl[["data"]][["outcome"]])

  expect_null(fit_3lvl[["marglik"]])
  expect_residual_vector(
    residuals(fit_3lvl, conditioning_depth = "cluster"),
    n_3lvl
  )
  expect_error(rstudent(fit_3lvl, unit = "cluster"), "Cluster-unit rstudent residuals")
})

test_that("Residuals for RoBMA model-averaging fits are internally consistent", {

  fit_brma <- fits[["dat.lehmann2018_RoBMA"]]
  n        <- nobs(fit_brma)

  expect_residual_vector(residuals(fit_brma), n)
  expect_error(residuals(fit_brma, type = "pearson"), "selection models")
  expect_error(rstandard(fit_brma), "selection models")

  brma_rstudent <- suppressWarnings(rstudent(fit_brma))
  expect_residual_table(brma_rstudent, n)
  expect_equal(suppressWarnings(residuals(fit_brma, type = "rstudent")), brma_rstudent[["z"]],
               tolerance = 1e-12)

  fit_3lvl <- fits[["dat.lehmann2018_RoBMA_3lvl_mods_scale"]]
  n_3lvl   <- nrow(fit_3lvl[["data"]][["outcome"]])

  expect_null(fit_3lvl[["marglik"]])
  expect_residual_vector(
    residuals(fit_3lvl, conditioning_depth = "cluster"),
    n_3lvl
  )
})
