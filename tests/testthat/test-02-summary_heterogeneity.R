context("Summary Heterogeneity")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "summary_heterogeneity")

# list cached fits lazily
skip_if_no_fits()
skip_if_not_installed("metafor")
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

summary_heterogeneity_columns <- c("Mean", "Median", "0.025", "0.975")
model_averaged_classes        <- c("BMA.norm", "BMA.glmm", "RoBMA")

expect_summary_heterogeneity_structure <- function(heterogeneity, expected_rows, name) {

  expect_true(inherits(heterogeneity, "summary_heterogeneity.brma"),
              info = paste0("summary_heterogeneity class for '", name, "'"))
  expect_equal(sort(rownames(heterogeneity$estimates)), sort(expected_rows),
               info = paste0("summary_heterogeneity rows for '", name, "'"))

  estimates <- heterogeneity$estimates[
    expected_rows, summary_heterogeneity_columns, drop = FALSE
  ]
  values <- as.matrix(estimates)

  expect_true(all(is.finite(values)),
              info = paste0("summary_heterogeneity finite estimates for '", name, "'"))
  expect_true(all(values >= 0),
              info = paste0("summary_heterogeneity non-negative estimates for '", name, "'"))

  i2_rows <- grep("^I2", expected_rows, value = TRUE)
  if (length(i2_rows) > 0) {
    i2_values <- as.matrix(heterogeneity$estimates[
      i2_rows, summary_heterogeneity_columns, drop = FALSE
    ])
    expect_true(all(i2_values >= 0 & i2_values <= 100),
                info = paste0("summary_heterogeneity I2 bounds for '", name, "'"))
  }

  expect_true(all(heterogeneity$estimates["H2", summary_heterogeneity_columns] >= 1),
              info = paste0("summary_heterogeneity H2 bounds for '", name, "'"))
}


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
  expect_equal(brma_het$estimates["I2", "Median"],
               fit_metafor$I2, tolerance = 0.05,
               info = "brma I2 should match metafor")

  # H2
  expect_equal(brma_het$estimates["H2", "Median"],
               fit_metafor$H2, tolerance = 0.05,
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
               sqrt(fit_metafor$tau2), tolerance = 0.05,
               info = "brma total tau should match metafor")

  # tau [within] should match sqrt(sigma2[1]) (estimate-level)
  expect_equal(brma_het$estimates["tau [within]", "Mean"],
               sqrt(fit_metafor$tau2 * (1-fit_metafor$rho)), tolerance = 0.05,
               info = "brma tau within should match metafor sigma2[1]")

  # tau [between] should match sqrt(sigma2[2]) (study-level)
  expect_equal(brma_het$estimates["tau [between]", "Mean"],
               sqrt(fit_metafor$tau2 * fit_metafor$rho), tolerance = 0.05,
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
  I2_within  <- 100 * (total_var * (1-fit_metafor$rho)) / (total_var + typical_v)
  I2_between <- 100 * (total_var * fit_metafor$rho) / (total_var + typical_v)

  expect_equal(brma_het$estimates["I2", "Median"],
               I2_total, tolerance = 0.05,
               info = "brma total I2 should match metafor formula")

  expect_equal(brma_het$estimates["I2 [within]", "Median"],
               I2_within, tolerance = 0.10,
               info = "brma I2 within should match metafor formula")

  expect_equal(brma_het$estimates["I2 [between]", "Median"],
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
               sqrt(fit_metafor$tau2), tolerance = 0.10,
               info = "brma tau should match metafor for meta-regression")

  # tau2
  expect_equal(brma_het$estimates["tau2", "Mean"],
               fit_metafor$tau2, tolerance = 0.10,
               info = "brma tau2 should match metafor for meta-regression")

  # I2
  expect_equal(brma_het$estimates["I2", "Median"],
               fit_metafor$I2, tolerance = 0.10,
               info = "brma I2 should match metafor for meta-regression")

  # H2
  expect_equal(brma_het$estimates["H2", "Median"],
               fit_metafor$H2, tolerance = 0.20,
               info = "brma H2 should match metafor for meta-regression")
})

# ============================================================================ #
# Test: Location-Scale Model Heterogeneity
# ============================================================================ #

test_that("Heterogeneity for location-scale model matches metafor", {

  name        <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # For location-scale models, heterogeneity varies by observation
  # summary_heterogeneity should aggregate across the scale model matrix
  brma_het <- summary_heterogeneity(fit_brma)

  # I2 and H2 - compare against metafor's values
  # Note: metafor's I2 and H2 for location-scale models are computed differently
  # (they use the average tau across observations)
  expect_equal(brma_het$estimates["I2", "Median"],
               mean(fit_metafor$I2), tolerance = 0.05,
               info = "brma I2 should match metafor for location-scale models")

  expect_equal(brma_het$estimates["H2", "Median"],
               mean(fit_metafor$H2), tolerance = 0.20,
               info = "brma H2 should match metafor for location-scale models")
})

# ============================================================================ #
# Test: bPET Model Heterogeneity
# ============================================================================ #

test_that("Heterogeneity for bPET model matches metafor", {

  name        <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  brma_het <- summary_heterogeneity(fit_brma)

  # tau (should match residual heterogeneity from metafor)
  expect_equal(brma_het$estimates["tau", "Mean"],
               sqrt(fit_metafor$tau2), tolerance = 0.05,
               info = "brma tau should match metafor for bPET model")

  # tau2
  expect_equal(brma_het$estimates["tau2", "Mean"],
               fit_metafor$tau2, tolerance = 0.05,
               info = "brma tau2 should match metafor for bPET model")

  # I2
  expect_equal(brma_het$estimates["I2", "Median"],
               fit_metafor$I2, tolerance = 0.05,
               info = "brma I2 should match metafor for bPET model")

  # H2
  expect_equal(brma_het$estimates["H2", "Median"],
               fit_metafor$H2, tolerance = 0.05,
               info = "brma H2 should match metafor for bPET model")
})

test_that("Heterogeneity for bPET model (negative effects) matches metafor", {

  name        <- "dat.lehmann2018-PET_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  brma_het <- summary_heterogeneity(fit_brma)

  # tau (should match residual heterogeneity from metafor)
  expect_equal(brma_het$estimates["tau", "Mean"],
               sqrt(fit_metafor$tau2), tolerance = 0.05,
               info = "brma tau should match metafor for bPET model (negative)")

  # tau2
  expect_equal(brma_het$estimates["tau2", "Mean"],
               fit_metafor$tau2, tolerance = 0.05,
               info = "brma tau2 should match metafor for bPET model (negative)")

  # I2
  expect_equal(brma_het$estimates["I2", "Median"],
               fit_metafor$I2, tolerance = 0.05,
               info = "brma I2 should match metafor for bPET model (negative)")

  # H2
  expect_equal(brma_het$estimates["H2", "Median"],
               fit_metafor$H2, tolerance = 0.05,
               info = "brma H2 should match metafor for bPET model (negative)")
})

# ============================================================================ #
# Test: bselmodel Heterogeneity
# ============================================================================ #

test_that("Heterogeneity for bselmodel matches metafor", {

  name        <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  brma_het <- summary_heterogeneity(fit_brma)

  # tau (should match heterogeneity from metafor selmodel)
  expect_equal(brma_het$estimates["tau", "Mean"],
               sqrt(fit_metafor$tau2), tolerance = 0.05,
               info = "brma tau should match metafor for bselmodel")

  # tau2
  expect_equal(brma_het$estimates["tau2", "Mean"],
               fit_metafor$tau2, tolerance = 0.05,
               info = "brma tau2 should match metafor for bselmodel")

  # I2 / H2 are not reported for selection models
})

test_that("Heterogeneity for bselmodel (negative effects) matches metafor", {

  name        <- "dat.lehmann2018-3PSM_neg"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  brma_het <- summary_heterogeneity(fit_brma)

  # tau (should match heterogeneity from metafor selmodel)
  expect_equal(brma_het$estimates["tau", "Mean"],
               sqrt(fit_metafor$tau2), tolerance = 0.05,
               info = "brma tau should match metafor for bselmodel (negative)")

  # tau2
  expect_equal(brma_het$estimates["tau2", "Mean"],
               fit_metafor$tau2, tolerance = 0.05,
               info = "brma tau2 should match metafor for bselmodel (negative)")

  # I2 / H2 are not reported for selection models
})

# ============================================================================ #
# Test: GLMM Heterogeneity
# ============================================================================ #

test_that("Heterogeneity for GLMM model matches metafor", {

  name        <- "nielweise2008_glmm"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  brma_het <- summary_heterogeneity(fit_brma)

  # tau (should match heterogeneity from metafor selmodel)
  expect_equal(brma_het$estimates["tau", "Mean"],
               sqrt(fit_metafor$tau2), tolerance = 0.10,
               info = "brma tau should match metafor for glmm")

  # tau2
  expect_equal(brma_het$estimates["tau2", "Mean"],
               fit_metafor$tau2, tolerance = 0.10,
               info = "brma tau2 should match metafor for glmm")

  # I2
  expect_equal(brma_het$estimates["I2", "Median"],
               fit_metafor$I2, tolerance = 0.20, # the actual value is quite close
               info = "brma I2 should match metafor for glmm")

  # H2
  expect_equal(brma_het$estimates["H2", "Median"],
               fit_metafor$H2, tolerance = 0.10,
               info = "brma H2 should match metafor for glmm")
})

# ============================================================================ #
# Test: Model-Averaged 2-Level Heterogeneity
# ============================================================================ #

test_that("Heterogeneity for model-averaged 2-level models is summarized", {

  model_averaged_names <- setdiff(
    catalog_fits(class = model_averaged_classes),
    catalog_fits(class = model_averaged_classes, feature = "multilevel")
  )
  skip_if_missing_fits(model_averaged_names)

  expected_rows <- c("tau", "tau2", "I2", "H2")

  for (name in model_averaged_names) {
    heterogeneity <- summary_heterogeneity(fits[[name]])
    expect_summary_heterogeneity_structure(heterogeneity, expected_rows, name)
  }
})

# ============================================================================ #
# Test: Model-Averaged 3-Level Heterogeneity
# ============================================================================ #

test_that("Heterogeneity for model-averaged 3-level models is partitioned", {

  model_averaged_names <- catalog_fits(
    class   = model_averaged_classes,
    feature = "multilevel"
  )
  skip_if_missing_fits(model_averaged_names)

  expected_rows <- c(
    "tau",
    "tau [within]",
    "tau [between]",
    "tau2",
    "tau2 [within]",
    "tau2 [between]",
    "I2",
    "I2 [within]",
    "I2 [between]",
    "H2"
  )

  for (name in model_averaged_names) {
    heterogeneity <- summary_heterogeneity(fits[[name]])
    expect_summary_heterogeneity_structure(heterogeneity, expected_rows, name)

    expect_equal(
      heterogeneity$estimates["tau2", "Mean"],
      heterogeneity$estimates["tau2 [within]", "Mean"] +
        heterogeneity$estimates["tau2 [between]", "Mean"],
      tolerance = 0.01,
      info      = paste0("tau2 partitions sum to total for '", name, "'")
    )
    expect_equal(
      heterogeneity$estimates["I2", "Mean"],
      heterogeneity$estimates["I2 [within]", "Mean"] +
        heterogeneity$estimates["I2 [between]", "Mean"],
      tolerance = 0.01,
      info      = paste0("I2 partitions sum to total for '", name, "'")
    )
  }
})

test_that("summary_heterogeneity print keeps leading blank line", {

  skip_if_missing_fits("bcg_meta-analysis")

  output <- capture.output(print(summary_heterogeneity(fits[["bcg_meta-analysis"]])))

  expect_equal(output[[1]], "")
})
