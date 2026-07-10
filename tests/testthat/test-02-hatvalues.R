context("Hatvalues")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(hatvalue_cases(), function(case) {
  test_that_case("Hatvalues match metafor", case, {
    expect_hatvalues_match_metafor(case)
  })
})

test_that("Hatvalues reject unsupported single-model families", {

  model_names <- c("nielweise2008_glmm", "dat.lehmann2018-3PSM")
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    expect_error(hatvalues(fits[[name]]), "available", info = name)
  }
})

test_that("Hatvalues for BMA.norm fits are internally consistent", {

  model_names <- c("dat.lehmann2018_BMA.norm", "dat.lehmann2018_BMA.norm_mods")
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    fit <- fits[[name]]
    expect_hatvalues_vector(hatvalues(fit), nobs(fit), info = name)
  }
})

test_that("Hatvalues support brma.mv known-V marginal GLS targets", {

  model_names <- c(
    "brma.mv_block_mvn",
    "brma.mv_block_mvn_random_scale",
    "brma.mv_block_mvn_known_R"
  )
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    fit <- fits[[name]]
    expected <- colMeans(.compute_hat_matrix_samples(
      object             = fit,
      conditioning_depth = "marginal",
      return_full_H      = FALSE,
      return_se          = FALSE
    )[["H_diag"]])

    actual <- hatvalues(fit)
    expect_hatvalues_vector(actual, nobs(fit), info = name)
    expect_equal(unname(actual), unname(expected), tolerance = 1e-12,
                 info = name)
  }
})

test_that("Hatvalues chunk known-V covariance without changing results", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit          <- fits[[name]]
  expected     <- hatvalues(fit)
  K            <- nobs(fit)
  one_draw_mem <- .known_v_covariance_peak_bytes(1L, K)
  old_options  <- options(
    RoBMA.known_v_covariance_max_bytes = one_draw_mem
  )
  on.exit(options(old_options), add = TRUE)

  actual   <- hatvalues(fit)
  metadata <- attr(actual, "known_v_diagnostic")
  attr(actual, "known_v_diagnostic") <- NULL

  expect_equal(unname(actual), unname(expected), tolerance = 1e-12)
  expect_true(metadata[["n_chunks"]] > 1L)
  expect_equal(metadata[["n_used_samples"]], metadata[["n_posterior_samples"]])
})

test_that("Hatvalues expose deterministic known-V draw thinning", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit <- fits[[name]]

  actual <- expect_warning(
    hatvalues(fit, max_samples = 10),
    "deterministically thinned"
  )
  metadata <- attr(actual, "known_v_diagnostic")

  expect_hatvalues_vector(actual, nobs(fit), info = name)
  expect_equal(metadata[["n_used_samples"]], 10L)
  expect_true(metadata[["thinned"]])
})

test_that("Hatvalues use study labels as names", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit <- fits[[name]]
  expect_equal(names(hatvalues(fit)), RoBMA:::.diagnostic_study_labels(fit))
})

test_that("Hatvalues reject unsupported model-averaging families", {

  model_names <- c("bcg_BMA.glmm", "dat.lehmann2018_RoBMA")
  skip_if_missing_fits(model_names)

  expect_error(hatvalues(fits[["bcg_BMA.glmm"]]), "normal outcome models")
  expect_error(hatvalues(fits[["dat.lehmann2018_RoBMA"]]), "selection models")
})
