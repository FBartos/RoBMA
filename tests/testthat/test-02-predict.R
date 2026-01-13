context("Predict and wrappers")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "predict")

# list & load all fits
skip_if_no_fits()
fits <- lapply(list_fits(), load_fit)
info <- lapply(list_fits(), load_info)
names(fits) <- list_fits()
names(info) <- list_fits()


# ============================================================================ #
# Test: Simple Meta-Analysis Predictions
# ============================================================================ #

test_that("Predictions for simple meta-analysis match metafor", {

  skip_if_not_installed("metafor")

  name        <- "bcg_meta-analysis"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Pooled effect: compare brma vs metafor
  # --------------------------------------------------

  # metafor: pooled effect is in beta
  metafor_mu <- fit_metafor$beta[[1]]

  # brma: using predict directly
  brma_pred_terms <- predict(fit_brma, type = "terms", newdata = TRUE)
  brma_mu_predict <- brma_pred_terms$summary["mu", "Mean"]

  # brma: using pooled_effect wrapper
  brma_pooled <- pooled_effect(fit_brma)
  brma_mu_wrapper <- brma_pooled$summary["mu", "Mean"]

  # comparisons
  expect_equal(brma_mu_predict, brma_mu_wrapper,
               info = "predict(type='terms') should match pooled_effect()")
  expect_equal(brma_mu_predict, metafor_mu, tolerance = 0.05,
               info = "brma pooled effect should match metafor beta")

  # --------------------------------------------------
  # Pooled heterogeneity: compare brma vs metafor
  # --------------------------------------------------

  # metafor: tau is sqrt(tau2)
  metafor_tau <- sqrt(fit_metafor$tau2)

  # brma: using predict directly
  brma_pred_scale <- predict(fit_brma, type = "terms.scale", newdata = TRUE)
  brma_tau_predict <- brma_pred_scale$summary["tau", "Mean"]

  # brma: using pooled_heterogeneity wrapper
  brma_pooled_het <- pooled_heterogeneity(fit_brma)
  brma_tau_wrapper <- brma_pooled_het$summary["tau", "Mean"]

  # comparisons
  expect_equal(brma_tau_predict, brma_tau_wrapper,
               info = "predict(type='terms.scale') should match pooled_heterogeneity()")
  expect_equal(brma_tau_predict, metafor_tau, tolerance = 0.05,
               info = "brma pooled tau should match metafor tau")

  # --------------------------------------------------
  # BLUPs (true study effects): compare brma vs metafor
  # --------------------------------------------------

  # metafor: blup returns list with pred element
  metafor_blup <- metafor::blup(fit_metafor)
  metafor_theta <- metafor_blup$pred

  # brma: using predict directly
  brma_pred_effect <- predict(fit_brma, type = "effect")
  brma_theta_predict <- brma_pred_effect$summary[, "Mean"]

  # brma: using blup wrapper
  brma_blup <- blup(fit_brma)
  brma_theta_blup <- brma_blup$summary[, "Mean"]

  # brma: using true_effects wrapper (alias)
  brma_true_eff <- true_effects(fit_brma)
  brma_theta_true <- brma_true_eff$summary[, "Mean"]

  # comparisons: wrappers should be identical
  expect_equal(brma_theta_predict, brma_theta_blup,
               info = "predict(type='effect') should match blup()")
  expect_equal(brma_theta_blup, brma_theta_true,
               info = "blup() should match true_effects()")

  # comparison with metafor (allow MCMC tolerance)
  expect_equal(brma_theta_blup, metafor_theta, tolerance = 0.05,
               info = "brma BLUPs should match metafor BLUPs")

})


# ============================================================================ #
# Test: Meta-Regression Predictions
# ============================================================================ #

test_that("Predictions for meta-regression match metafor", {

  skip_if_not_installed("metafor")

  name        <- "bcg_meta-regression"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Study-specific predictions: compare brma vs metafor
  # --------------------------------------------------

  # metafor: predict returns fitted values for each study
  metafor_pred <- predict(fit_metafor)
  metafor_fitted <- metafor_pred$pred

  # brma: predict with type = "terms" and newdata = NULL gives per-study mu
  brma_pred_terms <- predict(fit_brma, type = "terms")
  brma_fitted <- brma_pred_terms$summary[, "Mean"]

  # comparison (moderator regression can have larger deviations)
  expect_equal(brma_fitted, metafor_fitted, tolerance = 0.05,
               info = "brma per-study predictions should match metafor fitted values")

  # --------------------------------------------------
  # Pooled effect: average across moderator levels
  # --------------------------------------------------

  # metafor: to get pooled effect, we average across model matrix
  # X %*% beta gives fitted, mean of that is pooled
  metafor_pooled_mu <- mean(metafor_fitted)

  # brma: using pooled_effect wrapper (newdata = TRUE averages across model matrix)
  brma_pooled <- pooled_effect(fit_brma)
  brma_pooled_mu <- brma_pooled$summary["mu", "Mean"]

  expect_equal(brma_pooled_mu, metafor_pooled_mu, tolerance = 0.05,
               info = "brma pooled effect should match mean of metafor fitted values")

  # --------------------------------------------------
  # BLUPs for meta-regression
  # --------------------------------------------------

  # metafor: blup for meta-regression
  metafor_blup <- metafor::blup(fit_metafor)
  metafor_theta <- metafor_blup$pred

  # brma: using blup wrapper
  brma_blup <- blup(fit_brma)
  brma_theta <- brma_blup$summary[, "Mean"]

  expect_equal(brma_theta, metafor_theta, tolerance = 0.05,
               info = "brma BLUPs should match metafor BLUPs for meta-regression")

})


# ============================================================================ #
# Test: Location-Scale Model Predictions
# ============================================================================ #

test_that("Predictions for location-scale model match metafor", {

  skip_if_not_installed("metafor")

  name        <- "bangertdrowns2004_location-scale"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Study-specific location predictions: compare brma vs metafor

  # --------------------------------------------------

  # metafor: predict returns fitted values for each study
  metafor_pred <- predict(fit_metafor)
  metafor_fitted <- metafor_pred$pred

  # brma: predict with type = "terms" and newdata = NULL gives per-study mu
  brma_pred_terms <- predict(fit_brma, type = "terms")
  brma_fitted <- brma_pred_terms$summary[, "Mean"]

  expect_equal(brma_fitted, metafor_fitted, tolerance = 0.05,
               info = "brma per-study location predictions should match metafor")

  # --------------------------------------------------
  # Study-specific tau values: compare brma vs metafor
  # --------------------------------------------------

  # metafor: tau^2_i = exp(Z %*% alpha), so tau_i = sqrt(exp(Z %*% alpha))
  metafor_tau2_studies <- as.vector(exp(fit_metafor$Z %*% fit_metafor$alpha))
  metafor_tau_studies  <- sqrt(metafor_tau2_studies)

  # brma: predict with type = "terms.scale" and newdata = NULL
  brma_pred_scale  <- predict(fit_brma, type = "terms.scale")
  brma_tau_studies <- brma_pred_scale$summary[, "Mean"]

  expect_equal(brma_tau_studies, metafor_tau_studies, tolerance = 0.05,
               info = "brma per-study tau should match metafor")

  # --------------------------------------------------
  # Pooled heterogeneity: compare brma vs metafor
  # --------------------------------------------------

  # metafor: pooled tau is mean of study-specific tau values
  metafor_tau_pooled <- mean(metafor_tau_studies)

  # brma: pooled_heterogeneity averages tau across scale model matrix
  brma_pooled_het <- pooled_heterogeneity(fit_brma)
  brma_tau_pooled <- brma_pooled_het$summary["tau", "Mean"]

  expect_equal(brma_tau_pooled, metafor_tau_pooled, tolerance = 0.05,
               info = "brma pooled tau should match mean of metafor study-specific tau")

})


# ============================================================================ #
# Test: 3-Level Model Predictions
# ============================================================================ #

test_that("Predictions for 3-level model match metafor", {

  skip_if_not_installed("metafor")

  name        <- "konstantopoulos2011_3lvl"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Pooled effect
  # --------------------------------------------------

  metafor_mu <- fit_metafor$beta[[1]]

  brma_pooled <- pooled_effect(fit_brma)
  brma_mu <- brma_pooled$summary["mu", "Mean"]

  expect_equal(brma_mu, metafor_mu, tolerance = 0.05,
               info = "brma pooled effect should match metafor for 3-level model")

  # --------------------------------------------------
  # Total heterogeneity (tau = sqrt(tau_within^2 + tau_between^2))
  # --------------------------------------------------

  # metafor: tau2 is total variance
  metafor_tau <- sqrt(fit_metafor$tau2)

  brma_pooled_het <- pooled_heterogeneity(fit_brma)
  brma_tau <- brma_pooled_het$summary["tau", "Mean"]

  expect_equal(brma_tau, metafor_tau, tolerance = 0.05,
               info = "brma total tau should match metafor for 3-level model")

  # --------------------------------------------------
  # BLUPs for 3-level model
  # --------------------------------------------------

  # The metafor function does not allow computing BLUPs for multivariate models
  # metafor_blup <- metafor::blup(fit_metafor)
  # metafor_theta <- metafor_blup$pred

  # brma: using blup wrapper
  brma_blup <- blup(fit_brma)
  brma_theta <- brma_blup$summary[, "Mean"]

  # Verify BLUPs are computed for all observations
  n_obs <- nrow(fit_brma$data$outcome)
  expect_equal(length(brma_theta), n_obs,
               info = "BLUPs should be computed for all observations")

  # Verify BLUPs show shrinkage (should be closer to pooled mean than raw yi)
  yi <- fit_brma$data$outcome$yi
  shrinkage_blup <- mean(abs(brma_theta - brma_mu))
  shrinkage_yi   <- mean(abs(yi - brma_mu))
  expect_true(shrinkage_blup < shrinkage_yi,
              info = "BLUPs should show shrinkage toward pooled estimate")

})


# ============================================================================ #
# Test: GLMM Model Predictions
# ============================================================================ #

test_that("Predictions for GLMM model match metafor", {

  skip_if_not_installed("metafor")

  name        <- "bcg_glmm"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Pooled effect: compare brma vs metafor
  # --------------------------------------------------

  metafor_mu <- fit_metafor$beta[[1]]

  brma_pooled <- pooled_effect(fit_brma)
  brma_mu <- brma_pooled$summary["mu", "Mean"]

  expect_equal(brma_mu, metafor_mu, tolerance = 0.05,
               info = "brma pooled effect should match metafor for GLMM model")

  # --------------------------------------------------
  # Pooled heterogeneity: compare brma vs metafor
  # --------------------------------------------------

  metafor_tau <- sqrt(fit_metafor$tau2)

  brma_pooled_het <- pooled_heterogeneity(fit_brma)
  brma_tau <- brma_pooled_het$summary["tau", "Mean"]

  expect_equal(brma_tau, metafor_tau, tolerance = 0.05,
               info = "brma pooled tau should match metafor for GLMM model")

  # --------------------------------------------------
  # BLUPs for GLMM model
  # --------------------------------------------------

  # The metafor function does not allow computing BLUPs for multivariate models
  # metafor_blup <- metafor::blup(fit_metafor)
  # metafor_theta <- metafor_blup$pred

  brma_blup <- blup(fit_brma)
  brma_theta <- brma_blup$summary[, "Mean"]

  # Verify BLUPs are computed for all observations
  n_obs <- nrow(fit_brma$data$outcome)
  expect_equal(length(brma_theta), n_obs,
               info = "BLUPs should be computed for all observations")

})


# ============================================================================ #
# Test: Selection Model Predictions
# ============================================================================ #

test_that("Predictions for selection model match metafor", {

  skip_if_not_installed("metafor")

  name        <- "dat.lehmann2018-3PSM"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Pooled effect: compare brma vs metafor
  # --------------------------------------------------

  metafor_mu <- fit_metafor$beta[[1]]

  brma_pooled <- pooled_effect(fit_brma)
  brma_mu <- brma_pooled$summary["mu", "Mean"]

  expect_equal(brma_mu, metafor_mu, tolerance = 0.05,
               info = "brma pooled effect should match metafor for selection model")

  # --------------------------------------------------
  # Pooled heterogeneity: compare brma vs metafor
  # --------------------------------------------------

  metafor_tau <- sqrt(fit_metafor$tau2)

  brma_pooled_het <- pooled_heterogeneity(fit_brma)
  brma_tau <- brma_pooled_het$summary["tau", "Mean"]

  expect_equal(brma_tau, metafor_tau, tolerance = 0.05,
               info = "brma pooled tau should match metafor for selection model")


  # --------------------------------------------------
  # BLUPs for selmodel model
  # --------------------------------------------------

  # The metafor function does not allow computing BLUPs for selmodel
  # metafor_blup <- metafor::blup(fit_metafor)
  # metafor_theta <- metafor_blup$pred

  brma_blup <- blup(fit_brma)
  brma_theta <- brma_blup$summary[, "Mean"]

  # Verify BLUPs are computed for all observations
  n_obs <- nrow(fit_brma$data$outcome)
  expect_equal(length(brma_theta), n_obs,
               info = "BLUPs should be computed for all observations")

  # Verify BLUPs show shrinkage (should be closer to pooled mean than raw yi)
  yi <- fit_brma$data$outcome$yi
  shrinkage_blup <- mean(abs(brma_theta - brma_mu))
  shrinkage_yi   <- mean(abs(yi - brma_mu))
  expect_true(shrinkage_blup < shrinkage_yi,
              info = "BLUPs should show shrinkage toward pooled estimate")

})


# ============================================================================ #
# Test: PET Model Predictions
# ============================================================================ #

test_that("Predictions for PET model match metafor", {

  skip_if_not_installed("metafor")

  name        <- "dat.lehmann2018-PET"
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]

  # --------------------------------------------------
  # Pooled effect at sei = 0: compare brma vs metafor
  # --------------------------------------------------

  # For PET model, pooled effect is the intercept (prediction at sei = 0)
  # metafor: predict at sei = 0 (newmods = 0 for the sei moderator)
  metafor_pred_sei0 <- predict(fit_metafor, newmods = 0)
  metafor_mu <- metafor_pred_sei0$pred

  # brma: pooled_effect automatically adjusts for sei = 0
  brma_pooled <- pooled_effect(fit_brma)
  brma_mu <- brma_pooled$summary["mu", "Mean"]

  expect_equal(brma_mu, metafor_mu, tolerance = 0.05,
               info = "brma pooled effect (at sei=0) should match metafor for PET model")

  # --------------------------------------------------
  # Pooled heterogeneity: compare brma vs metafor
  # --------------------------------------------------

  metafor_tau <- sqrt(fit_metafor$tau2)

  brma_pooled_het <- pooled_heterogeneity(fit_brma)
  brma_tau <- brma_pooled_het$summary["tau", "Mean"]

  expect_equal(brma_tau, metafor_tau, tolerance = 0.05,
               info = "brma pooled tau should match metafor for PET model")

  # --------------------------------------------------
  # Study-specific predictions: compare brma vs metafor
  # --------------------------------------------------

  # metafor: predict returns fitted values for each study
  metafor_pred <- predict(fit_metafor, newmods = rep(0, 81))
  metafor_fitted <- metafor_pred$pred

  # brma: predict with type = "terms" and newdata = NULL gives per-study mu
  brma_pred_terms <- predict(fit_brma, type = "terms")
  brma_fitted <- brma_pred_terms$summary[, "Mean"]

  expect_equal(brma_fitted, metafor_fitted, tolerance = 0.05,
               info = "brma per-study predictions should match metafor for PET model")

  # --------------------------------------------------
  # Study-specific predictions: compare brma vs metafor (unadjusted for bias)
  # --------------------------------------------------

  # metafor: predict returns fitted values for each study
  metafor_pred <- predict(fit_metafor)
  metafor_fitted <- metafor_pred$pred

  # brma: predict with type = "terms" and newdata = NULL gives per-study mu
  brma_pred_terms <- predict(fit_brma, type = "terms", bias_adjusted = FALSE)
  brma_fitted <- brma_pred_terms$summary[, "Mean"]

  expect_equal(brma_fitted, metafor_fitted, tolerance = 0.06,
               info = "brma per-study predictions should match metafor for PET model")

  # metafor: blup
  metafor_blup <- metafor::blup(fit_metafor)
  metafor_theta <- metafor_blup$pred

  # brma: using blup wrapper
  brma_blup <- blup(fit_brma, bias_adjusted = FALSE)
  brma_theta <- brma_blup$summary[, "Mean"]

  # brma: using blup wrapper
  brma_blup_adjusted <- blup(fit_brma)
  brma_theta_adjusted <- brma_blup_adjusted$summary[, "Mean"]

  expect_equal(brma_theta, metafor_theta, tolerance = 0.05,
               info = "brma BLUPs should match metafor BLUPs for meta-regression")
  expect_true(mean(abs(brma_theta - brma_theta_adjusted)) > 0.10,
               info = "adjusted and unadjusted BLUPs should not match")

})


# ============================================================================ #
# Test: Wrapper Function Interface
# ============================================================================ #

test_that("Wrapper functions have correct interface", {

  name     <- "bcg_meta-analysis"
  fit_brma <- fits[[name]]

  # --------------------------------------------------
  # Test as_samples = TRUE returns matrix
  # --------------------------------------------------

  samples_effect <- pooled_effect(fit_brma, as_samples = TRUE)
  expect_true(is.matrix(samples_effect),
              info = "pooled_effect with as_samples = TRUE should return matrix")
  expect_equal(ncol(samples_effect), 1,
               info = "pooled_effect should return 1 column for aggregated estimate")

  samples_het <- pooled_heterogeneity(fit_brma, as_samples = TRUE)
  expect_true(is.matrix(samples_het),
              info = "pooled_heterogeneity with as_samples = TRUE should return matrix")

  samples_blup <- blup(fit_brma, as_samples = TRUE)
  expect_true(is.matrix(samples_blup),
              info = "blup with as_samples = TRUE should return matrix")
  n_studies <- nrow(fit_brma$data$outcome)
  expect_equal(ncol(samples_blup), n_studies,
               info = "blup should return one column per study")

  # --------------------------------------------------
  # Test probs argument
  # --------------------------------------------------

  pred_90ci <- pooled_effect(fit_brma, probs = c(0.05, 0.95))
  ci_names <- colnames(pred_90ci$summary)
  # BayesTools may format as "5%" or "0.05" - check for either format
  has_lower <- any(grepl("5%|0\\.05", ci_names))
  has_upper <- any(grepl("95%|0\\.95", ci_names))
  expect_true(has_lower && has_upper,
              info = "probs argument should control CI quantiles")

  # --------------------------------------------------
  # Test true_effects is alias for blup
  # --------------------------------------------------

  blup_result <- blup(fit_brma)
  true_eff_result <- true_effects(fit_brma)
  expect_equal(blup_result$summary, true_eff_result$summary,
               info = "true_effects should be identical to blup")

})


# ============================================================================ #
# Test: Quiet Operation of Wrappers
# ============================================================================ #

test_that("Wrappers suppress aggregation messages", {

  name     <- "bcg_meta-regression"
  fit_brma <- fits[[name]]

  # Wrappers should not produce messages
  expect_silent(pooled_effect(fit_brma))
  expect_silent(pooled_heterogeneity(fit_brma))
  expect_silent(blup(fit_brma))
  expect_silent(true_effects(fit_brma))

})
