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

test_that("Newdata prediction preserves duplicate rows and rejects novel factor levels", {

  skip_if_missing_fits(c("bcg_meta-regression", "bcg_meta-regression2"))

  fit_cont <- fits[["bcg_meta-regression"]]
  newdata_cont <- data.frame(
    ablat = c(10, 10, 70),
    year  = c(1950, 1950, 1980)
  )
  pred_cont <- predict(
    fit_cont,
    newdata = newdata_cont,
    type    = "terms",
    quiet   = TRUE
  )
  pred_single <- predict(
    fit_cont,
    newdata = newdata_cont[3, , drop = FALSE],
    type    = "terms",
    quiet   = TRUE
  )

  expect_equal(unname(pred_cont[, 1]), unname(pred_cont[, 2]),
               info = "duplicated newdata rows produce identical predictions")
  expect_equal(unname(pred_cont[, 3]), unname(pred_single[, 1]),
               info = "newdata row order is preserved")

  fit_factor <- fits[["bcg_meta-regression2"]]
  old_levels <- levels(fit_factor[["data"]][["mods"]][["alloc"]])
  newdata_factor <- data.frame(
    alloc = factor("not-in-training", levels = c(old_levels, "not-in-training"))
  )

  expect_error(
    predict(fit_factor, newdata = newdata_factor, type = "terms", quiet = TRUE),
    info = "novel factor levels are rejected"
  )
})

test_that("Predictions can use an already thinned posterior", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  fit_brma          <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
  selected_rows     <- unique(round(seq(
    from       = 1,
    to         = nrow(posterior_samples),
    length.out = min(1000L, nrow(posterior_samples))
  )))

  full_prediction <- as.matrix(predict(
    fit_brma,
    type  = "terms",
    quiet = TRUE
  ))
  thinned_prediction <- as.matrix(predict(
    fit_brma,
    type               = "terms",
    quiet              = TRUE,
    .posterior_samples = posterior_samples[selected_rows, , drop = FALSE]
  ))

  expect_equal(
    thinned_prediction,
    full_prediction[selected_rows, , drop = FALSE],
    tolerance = 1e-12,
    info      = "formula-based predictions use the caller-supplied posterior rows"
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


test_that("predict preserves the released positional type argument", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  set.seed(481)
  positional <- predict(fit_brma, NULL, "response")
  set.seed(481)
  named <- predict(fit_brma, newdata = NULL, type = "response")

  expect_identical(names(formals(predict.brma))[3L], "type")
  expect_equal(positional, named, tolerance = 0)
})


test_that("pooled_effect aggregates fitted-design location draws directly", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  terms <- predict(
    fit_brma,
    type          = "terms",
    bias_adjusted = TRUE,
    quiet         = TRUE
  )
  pooled   <- pooled_effect(fit_brma)
  expected <- matrix(rowMeans(as.matrix(terms)), ncol = 1L)
  expect_equal(unname(as.matrix(pooled)), unname(expected), tolerance = 1e-12)

  testthat::local_mocked_bindings(
    predict.brma = function(...) stop("pooled_effect used predict.brma"),
    .package = "RoBMA"
  )
  expect_silent(pooled_effect(fit_brma))
})

test_that("fitted returns in-sample posterior means", {

  model_names <- c(
    "bcg_meta-regression",
    "bangertdrowns2004_location-scale",
    "konstantopoulos2011_3lvl",
    "dat.lehmann2018_RoBMA"
  )
  skip_if_missing_fits(model_names)

  fit_reg <- fits[["bcg_meta-regression"]]
  n_reg   <- nobs(fit_reg)

  fitted_reg <- fitted(fit_reg)
  terms_reg  <- predict(fit_reg, type = "terms", quiet = TRUE)
  expect_finite_vector(fitted_reg, n_reg, "fitted vector")
  expect_equal(
    unname(fitted_reg),
    unname(.sample_means(terms_reg)),
    tolerance = 1e-12,
    info      = "default fitted values wrap marginal predictions"
  )
  expect_equal(
    unname(fitted_reg),
    unname(stats::fitted(info[["bcg_meta-regression"]][["metafor"]])),
    tolerance = 0.05,
    info      = "default fitted values match metafor"
  )
  expect_equal(
    unname(residuals(fit_reg)),
    unname(.outcome_data_yi(fit_reg) - fitted_reg),
    tolerance = 1e-12,
    info      = "default residuals use default fitted values"
  )

  fit_3lvl <- fits[["konstantopoulos2011_3lvl"]]
  expect_equal(
    unname(fitted(fit_3lvl, conditioning_depth = "cluster")),
    unname(.sample_means(predict(fit_3lvl, type = "cluster", quiet = TRUE))),
    tolerance = 1e-12,
    info      = "cluster fitted values wrap cluster predictions"
  )
  expect_equal(
    unname(fitted(fit_3lvl, conditioning_depth = "estimate")),
    unname(.sample_means(blup(fit_3lvl))),
    tolerance = 1e-12,
    info      = "estimate fitted values wrap BLUPs"
  )

  fit_scale <- fits[["bangertdrowns2004_location-scale"]]
  fitted_scale <- fitted(fit_scale, component = "scale")
  expect_equal(
    unname(fitted_scale),
    unname(.sample_means(predict(fit_scale, type = "terms.scale", quiet = TRUE))),
    tolerance = 1e-12,
    info      = "scale fitted values wrap scale predictions"
  )
  fitted_all <- fitted(fit_scale, component = "all")
  expect_equal(names(fitted_all), c("location", "scale"),
               info = "component all returns location and scale")
  expect_equal(unname(fitted_all[["location"]]), unname(fitted(fit_scale)),
               tolerance = 1e-12)
  expect_equal(unname(fitted(fit_scale, component = "mods")),
               unname(fitted(fit_scale)),
               tolerance = 1e-12,
               info = "component mods aliases location fitted values")
  expect_equal(unname(fitted_all[["scale"]]), unname(fitted_scale),
               tolerance = 1e-12)

  fit_robma <- fits[["dat.lehmann2018_RoBMA"]]
  expect_equal(
    unname(fitted(fit_robma, bias_adjusted = TRUE, conditional = TRUE)),
    unname(.sample_means(predict(
      fit_robma,
      type          = "terms",
      bias_adjusted = TRUE,
      conditional   = TRUE,
      quiet         = TRUE
    ))),
    tolerance = 1e-12,
    info      = "conditional fitted values wrap conditional predictions"
  )

  expect_error(fitted(fit_reg, unit = "cluster"), "multilevel models")
  expect_error(
    fitted(fit_scale, component = "scale", conditioning_depth = "estimate"),
    "component = 'scale'"
  )
  expect_error(
    fitted(fit_scale, component = "scale", transform = "EXP"),
    "location fitted values"
  )
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
    set.seed(247)
    pooled_effect_averaged <- pooled_effect(fit_brma)
    set.seed(247)
    pooled_effect_cond     <- expect_silent(
      pooled_effect(fit_brma, conditional = TRUE)
    )
    pooled_effect_terms <- expect_silent(
      predict(
        fit_brma,
        type          = "terms",
        bias_adjusted = TRUE,
        conditional   = TRUE,
        quiet         = TRUE
      )
    )
    pooled_effect_expected <- matrix(
      rowMeans(as.matrix(pooled_effect_terms)),
      ncol = 1L
    )

    expect_equal(
      unname(as.matrix(pooled_effect_cond)),
      unname(as.matrix(pooled_effect_averaged)[effect_rows, , drop = FALSE]),
      info = paste(name, "pooled_effect conditional rows")
    )
    expect_equal(
      unname(attr(pooled_effect_cond, "prediction_samples")),
      unname(attr(pooled_effect_averaged, "prediction_samples")[
        effect_rows,
        ,
        drop = FALSE
      ]),
      info = paste(name, "pooled prediction conditional rows")
    )
    expect_equal(
      unname(as.matrix(pooled_effect_cond)),
      unname(pooled_effect_expected),
      info = paste(name, "pooled_effect matches fitted-design row mean")
    )
    expect_match(attr(pooled_effect_cond, "title"), "Conditional")
  }

  fit_brma       <- fits[["dat.lehmann2018_RoBMA"]]
  summary_brma   <- summary(fit_brma, conditional = TRUE)
  summary_effect <- summary(suppressMessages(pooled_effect(fit_brma, conditional = TRUE)))
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
  pooled_het_cond     <- expect_silent(
    pooled_heterogeneity(fit_brma, conditional = TRUE)
  )
  expect_silent(predict(
    fit_brma,
    newdata     = NULL,
    type        = "terms.scale",
    conditional = TRUE,
    quiet       = TRUE
  ))
  expect_message(
    predict(
      fit_brma,
      newdata     = NULL,
      type        = "terms.scale",
      conditional = TRUE,
      quiet       = FALSE
    ),
    "flattened"
  )

  expect_equal(
    unname(as.matrix(pooled_het_cond)),
    unname(as.matrix(pooled_het_averaged)[heterogeneity_rows, , drop = FALSE]),
    info = "pooled_heterogeneity conditional rows"
  )
  expect_equal(
    unname(as.matrix(pooled_het_cond)),
    unname(.pooled_heterogeneity_total_samples(
      object            = fit_brma,
      posterior_samples = .get_posterior_samples(fit_brma[["fit"]])[
        heterogeneity_rows,
        ,
        drop = FALSE
      ]
    )),
    info = "pooled_heterogeneity matches conditional average scale design"
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
    "nielweise2008_BMA.glmm",
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
    pooled_terms <- predict(
      fit_brma,
      type          = "terms",
      bias_adjusted = TRUE,
      quiet         = TRUE
    )
    pooled_expected <- matrix(rowMeans(as.matrix(pooled_terms)), ncol = 1L)
    expect_brma_samples_matrix(pooled, 1, paste(name, "pooled_effect"))
    expect_equal(unname(as.matrix(pooled)), unname(pooled_expected))

    pooled_het <- pooled_heterogeneity(fit_brma)
    pooled_het_expected <- matrix(
      sqrt(rowMeans(as.matrix(scale)^2)),
      ncol = 1L
    )
    expect_brma_samples_matrix(pooled_het, 1, paste(name, "pooled_heterogeneity"))
    expect_equal(unname(as.matrix(pooled_het)), unname(pooled_het_expected))

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
    expect_equal(
      unname(as.matrix(ranef(fit_brma, component = "total"))),
      unname(as.matrix(ranef_samples))
    )
    ranef_list <- ranef(fit_brma, simplify = FALSE)
    expect_type(ranef_list, "list")
    expect_equal(names(ranef_list), "estimate")
    expect_equal(
      unname(as.matrix(ranef_list[["estimate"]])),
      unname(as.matrix(ranef_samples))
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
    expect_equal(
      unname(as.matrix(ranef(fit_brma, component = "cluster"))),
      unname(as.matrix(out[["cluster"]])),
      info = paste(name, "component cluster")
    )
    expect_equal(
      unname(as.matrix(ranef(fit_brma, component = "total"))),
      unname(as.matrix(out[["cluster"]]) + as.matrix(out[["estimate"]])),
      info = paste(name, "component total")
    )
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
