context("Summary Heterogeneity")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "summary_heterogeneity")

# list & load all fits
skip_if_no_fits()
skip_if_not_installed("metafor")
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


# ============================================================================ #
# Test: Simple Meta-Analysis Heterogeneity
# ============================================================================ #

test_that("Heterogeneity for simple meta-analysis matches metafor", {

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Compare tau, tau2, I2, H2 against metafor
  # --------------------------------------------------

  brma_het <- summary_heterogeneity(fit_brma)

  # tau
  expect_equal(brma_het$estimates["tau", "Mean"],
               sqrt(fit_metafor$tau2), tolerance = 0.05,
               info = "brma tau should match metafor")

  # tau2
  expect_equal(brma_het$estimates["tau2", "Mean"],
               fit_metafor$tau2, tolerance = 0.05,
               info = "brma tau2 should match metafor")

  # I2
  expect_equal(brma_het$estimates["I2", "Mean"],
               fit_metafor$I2, tolerance = 0.05,
               info = "brma I2 should match metafor")

  # H2
  expect_equal(brma_het$estimates["H2", "Mean"],
               fit_metafor$H2, tolerance = 0.10,
               info = "brma H2 should match metafor")
})


# ============================================================================ #
# Test: 3-Level Model Heterogeneity (Partitioned I2)
# ============================================================================ #

test_that("Heterogeneity for 3-level model matches metafor", {

  name        <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  brma_het <- summary_heterogeneity(fit_brma)

  # --------------------------------------------------
  # Compare variance components against metafor
  # --------------------------------------------------

  # Total tau should match sqrt(sum(sigma2))
  expect_equal(brma_het$estimates["tau", "Mean"],
               sqrt(sum(fit_metafor$sigma2)), tolerance = 0.05,
               info = "brma total tau should match metafor")

  # tau [within] should match sqrt(sigma2[1]) (estimate-level)
  expect_equal(brma_het$estimates["tau [within]", "Mean"],
               sqrt(fit_metafor$tau2 * fit_metafor$rho), tolerance = 0.05,
               info = "brma tau within should match metafor sigma2[1]")

  # tau [between] should match sqrt(sigma2[2]) (study-level)
  expect_equal(brma_het$estimates["tau [between]", "Mean"],
               sqrt(fit_metafor$tau2 * (1-fit_metafor$rho)), tolerance = 0.05,
               info = "brma tau between should match metafor sigma2[2]")

  # --------------------------------------------------
  # Compare partitioned I2 using metafor formula
  # --------------------------------------------------

  # Compute metafor I2 for multilevel models
  # (following https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate)
  W         <- diag(1/fit_metafor$vi)
  X         <- model.matrix(fit_metafor)
  P         <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typical_v <- (fit_metafor$k - fit_metafor$p) / sum(diag(P))
  total_var <- fit_metafor$tau2

  I2_total   <- 100 * total_var / (total_var + typical_v)
  I2_within  <- 100 * (total_var * fit_metafor$rho) / (total_var + typical_v)
  I2_between <- 100 * (total_var * (1-fit_metafor$rho)) / (total_var + typical_v)

  expect_equal(brma_het$estimates["I2", "Mean"],
               I2_total, tolerance = 0.05,
               info = "brma total I2 should match metafor formula")

  expect_equal(brma_het$estimates["I2 [within]", "Mean"],
               I2_within, tolerance = 0.10,
               info = "brma I2 within should match metafor formula")

  expect_equal(brma_het$estimates["I2 [between]", "Mean"],
               I2_between, tolerance = 0.10,
               info = "brma I2 between should match metafor formula")
})


# ============================================================================ #
# Test: Meta-Regression Heterogeneity
# ============================================================================ #

test_that("Heterogeneity for meta-regression matches metafor", {

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  brma_het <- summary_heterogeneity(fit_brma)

  # tau (should match residual heterogeneity from metafor)
  expect_equal(brma_het$estimates["tau", "Mean"],
               sqrt(fit_metafor$tau2), tolerance = 0.05,
               info = "brma tau should match metafor for meta-regression")

  # I2 and H2 should also be computed
  expect_true("I2" %in% rownames(brma_het$estimates),
              info = "I2 should be present for meta-regression")
  expect_true("H2" %in% rownames(brma_het$estimates),
              info = "H2 should be present for meta-regression")
})


# ============================================================================ #
# Test: Location-Scale Model Heterogeneity
# ============================================================================ #

test_that("Heterogeneity for location-scale model is computed", {

  name     <- "bangertdrowns2004_location-scale"
  fit_brma <- fits[[name]]

  # For location-scale models, heterogeneity varies by observation
  # summary_heterogeneity should aggregate across the scale model matrix
  brma_het <- summary_heterogeneity(fit_brma)

  expect_true("tau" %in% rownames(brma_het$estimates),
              info = "summary_heterogeneity should return tau for location-scale models")
  expect_true("I2" %in% rownames(brma_het$estimates),
              info = "summary_heterogeneity should return I2 for location-scale models")
  expect_true("H2" %in% rownames(brma_het$estimates),
              info = "summary_heterogeneity should return H2 for location-scale models")
})


# ============================================================================ #
# Test: Output Structure
# ============================================================================ #

test_that("Output structure is correct", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  brma_het <- summary_heterogeneity(fit_brma)

  # check class

  expect_s3_class(brma_het, "summary.brma_heterogeneity")

  # check components
  expect_true("estimates" %in% names(brma_het),
              info = "output should contain estimates")
  expect_s3_class(brma_het$estimates, "BayesTools_table")

  # check parameter names for non-multilevel
  expect_equal(rownames(brma_het$estimates),
               c("tau", "tau2", "I2", "H2"),
               info = "non-multilevel should have tau, tau2, I2, H2")
})


test_that("Multilevel output structure is correct", {

  name     <- "konstantopoulos2011_3lvl"
  fit_brma <- fits[[name]]

  brma_het <- summary_heterogeneity(fit_brma)

  # check parameter names for multilevel
  expected_params <- c("tau", "tau [within]", "tau [between]",
                       "tau2", "tau2 [within]", "tau2 [between]",
                       "I2", "I2 [within]", "I2 [between]", "H2")
  expect_equal(rownames(brma_het$estimates),
               expected_params,
               info = "multilevel should have partitioned heterogeneity")
})
