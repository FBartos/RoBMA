context("Predict and wrappers")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

REFERENCE_DIR <<- testthat::test_path("..", "results", "predict")

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(prediction_metafor_cases(), function(case) {
  test_that_case("Predictions match metafor", case, {
    expect_prediction_matches_metafor(case)
  })
})

for_each_case(prediction_newdata_metafor_cases(), function(case) {
  test_that_case("Newdata predictions match metafor", case, {
    expect_newdata_prediction_matches_metafor(case)
  })
})

test_that("Predictions for equivalent interaction parameterizations match", {

  model_names <- c("bcg_meta-regression4", "bcg_meta-regression4b")
  skip_if_missing_fits(model_names)

  fit_brma1 <- fits[["bcg_meta-regression4"]]
  fit_brma2 <- fits[["bcg_meta-regression4b"]]

  expect_equal(
    .sample_means(predict(fit_brma1, type = "terms")),
    .sample_means(predict(fit_brma2, type = "terms")),
    tolerance = 0.15,
    info      = "per-study predictions match"
  )
  expect_equal(
    .sample_mean(pooled_effect(fit_brma1), "mu"),
    .sample_mean(pooled_effect(fit_brma2), "mu"),
    tolerance = 0.05,
    info      = "pooled effects match"
  )
  expect_equal(
    .sample_means(blup(fit_brma1)),
    .sample_means(blup(fit_brma2)),
    tolerance = 0.05,
    info      = "BLUPs match"
  )
})

test_that("Wrapper functions have correct interface", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  samples_effect <- pooled_effect(fit_brma)
  expect_brma_samples_matrix(samples_effect, 1, "pooled_effect")

  samples_het <- pooled_heterogeneity(fit_brma)
  expect_s3_class(samples_het, "brma_samples")
  expect_true(is.matrix(samples_het), info = "pooled_heterogeneity matrix")

  samples_blup <- blup(fit_brma)
  expect_brma_samples_matrix(samples_blup, nrow(fit_brma$data$outcome), "blup")

  pred_90ci <- pooled_effect(fit_brma, probs = c(0.05, 0.95))
  ci_names  <- colnames(summary(pred_90ci))
  expect_true(any(grepl("5%|0\\.05", ci_names)) &&
                any(grepl("95%|0\\.95", ci_names)),
              info = "probs argument controls CI quantiles")

  expect_equal(summary(blup(fit_brma)), summary(true_effects(fit_brma)),
               info = "true_effects are identical to blup")
})

test_that("Conditional pooled wrappers condition RoBMA draws", {

  product_names <- c("dat.lehmann2018_RoBMA", "dat.lehmann2018_RoBMA_mods")
  skip_if_missing_fits(c(product_names, "bcg_meta-analysis"))

  for (name in product_names) {
    fit_brma <- fits[[name]]

    effect_rows <- .conditional_parameter_rows(
      object     = fit_brma,
      parameters = .conditional_effect_parameters(fit_brma)
    )
    pooled_effect_averaged <- pooled_effect(fit_brma)
    pooled_effect_cond     <- pooled_effect(fit_brma, conditional = TRUE)
    pooled_effect_predict  <- predict(
      fit_brma,
      newdata       = TRUE,
      type          = "terms",
      bias_adjusted = TRUE,
      conditional   = TRUE,
      quiet         = TRUE
    )

    expect_equal(
      unname(as.matrix(pooled_effect_cond)),
      unname(as.matrix(pooled_effect_averaged)[effect_rows, , drop = FALSE]),
      info = paste(name, "pooled_effect conditional rows")
    )
    expect_equal(
      unname(as.matrix(pooled_effect_cond)),
      unname(as.matrix(pooled_effect_predict)),
      info = paste(name, "pooled_effect matches conditional predict")
    )
    expect_match(attr(pooled_effect_cond, "title"), "Conditional")
  }

  fit_brma       <- fits[["dat.lehmann2018_RoBMA"]]
  summary_brma   <- summary(fit_brma, conditional = TRUE)
  summary_effect <- summary(pooled_effect(fit_brma, conditional = TRUE))
  expect_equal(
    unname(summary_effect["mu", "Mean"]),
    unname(summary_brma[["estimates_conditional"]]["mu", "Mean"]),
    info = "conditional pooled_effect matches conditional summary mu"
  )

  heterogeneity_rows <- .conditional_parameter_rows(
    object     = fit_brma,
    parameters = .conditional_heterogeneity_parameters(fit_brma)
  )
  pooled_het_averaged <- pooled_heterogeneity(fit_brma)
  pooled_het_cond     <- pooled_heterogeneity(fit_brma, conditional = TRUE)
  pooled_het_predict  <- predict(
    fit_brma,
    newdata     = TRUE,
    type        = "terms.scale",
    conditional = TRUE,
    quiet       = TRUE
  )

  expect_equal(
    unname(as.matrix(pooled_het_cond)),
    unname(as.matrix(pooled_het_averaged)[heterogeneity_rows, , drop = FALSE]),
    info = "pooled_heterogeneity conditional rows"
  )
  expect_equal(
    unname(as.matrix(pooled_het_cond)),
    unname(as.matrix(pooled_het_predict)),
    info = "pooled_heterogeneity matches conditional predict"
  )
  expect_equal(
    unname(summary(pooled_het_cond)["tau", "Mean"]),
    unname(summary_brma[["estimates_conditional"]]["tau", "Mean"]),
    info = "conditional pooled_heterogeneity matches conditional summary tau"
  )
  expect_match(attr(pooled_het_cond, "title"), "Conditional")

  expect_error(pooled_effect(fits[["bcg_meta-analysis"]], conditional = TRUE),
               "RoBMA objects")
  expect_error(pooled_heterogeneity(fits[["bcg_meta-analysis"]], conditional = TRUE),
               "RoBMA objects")
})

test_that("Model-averaged predictions cover BMA.norm, BMA.glmm, and RoBMA", {

  product_names <- c(
    "dat.lehmann2018_BMA.norm",
    "bcg_BMA.glmm",
    "dat.lehmann2018_RoBMA"
  )
  skip_if_missing_fits(product_names)

  for (name in product_names) {
    fit_brma  <- fits[[name]]
    n_studies <- nobs(fit_brma)

    terms  <- predict(fit_brma, type = "terms")
    scale  <- predict(fit_brma, type = "terms.scale")
    effect <- predict(fit_brma, type = "effect")
    set.seed(1)
    response <- predict(fit_brma, type = "response")

    expect_brma_samples_matrix(terms, n_studies, paste(name, "terms"))
    expect_brma_samples_matrix(scale, n_studies, paste(name, "terms.scale"))
    expect_brma_samples_matrix(effect, n_studies, paste(name, "effect"))
    expect_brma_samples_matrix(response, n_studies, paste(name, "response"))

    pooled <- pooled_effect(fit_brma)
    pooled_predict <- predict(
      fit_brma,
      newdata       = TRUE,
      type          = "terms",
      bias_adjusted = TRUE,
      quiet         = TRUE
    )
    expect_brma_samples_matrix(pooled, 1, paste(name, "pooled_effect"))
    expect_equal(unname(as.matrix(pooled)), unname(as.matrix(pooled_predict)))

    pooled_het <- pooled_heterogeneity(fit_brma)
    pooled_het_predict <- predict(
      fit_brma,
      newdata = TRUE,
      type    = "terms.scale",
      quiet   = TRUE
    )
    expect_brma_samples_matrix(pooled_het, 1, paste(name, "pooled_heterogeneity"))
    expect_equal(unname(as.matrix(pooled_het)), unname(as.matrix(pooled_het_predict)))

    blup_samples <- blup(fit_brma)
    true_samples <- true_effects(fit_brma)
    expect_brma_samples_matrix(blup_samples, n_studies, paste(name, "blup"))
    expect_brma_samples_matrix(true_samples, n_studies, paste(name, "true_effects"))
    expect_equal(unname(as.matrix(blup_samples)), unname(as.matrix(true_samples)))

    ranef_samples <- ranef(fit_brma)
    expect_brma_samples_matrix(ranef_samples, n_studies, paste(name, "ranef"))
    expect_equal(
      unname(as.matrix(ranef_samples)),
      unname(as.matrix(blup_samples) - as.matrix(terms))
    )
  }
})

test_that("Model-averaged multilevel ranef decomposes cluster and estimate effects", {

  product_names <- c(
    "bcg_BMA.glmm_3lvl_location_scale",
    "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  )
  skip_if_missing_fits(product_names)

  for (name in product_names) {
    fit_brma  <- fits[[name]]
    n_studies <- nobs(fit_brma)
    out       <- ranef(fit_brma)

    expect_type(out, "list")
    expect_equal(names(out), c("cluster", "estimate"), info = name)
    expect_brma_samples_matrix(out[["cluster"]], n_studies, paste(name, "cluster ranef"))
    expect_brma_samples_matrix(out[["estimate"]], n_studies, paste(name, "estimate ranef"))
  }
})

test_that("Wrappers suppress aggregation messages", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  expect_silent(pooled_effect(fit_brma))
  expect_silent(pooled_heterogeneity(fit_brma))
  expect_silent(blup(fit_brma))
  expect_silent(true_effects(fit_brma))
})
