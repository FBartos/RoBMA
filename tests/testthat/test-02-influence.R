context("Influence diagnostics")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

test_that("PSIS fitted-value influence helpers match direct calculations", {

  fit_samples <- matrix(
    c(
      0, 1, 2,
      1, 2, 3,
      2, 3, 4,
      3, 4, 5
    ),
    nrow  = 4,
    byrow = TRUE
  )
  colnames(fit_samples) <- paste0("study", 1:3)

  weights <- matrix(
    c(
      .4, .3, .2, .1,
      .1, .2, .3, .4,
      .25, .25, .25, .25
    ),
    nrow  = 4,
    byrow = FALSE
  )

  full_fit <- colMeans(fit_samples)
  loo_fit  <- crossprod(weights, fit_samples)
  loo_m2   <- crossprod(weights, fit_samples^2)
  loo_var  <- loo_m2 - loo_fit^2

  expected_dffits <- (full_fit - diag(loo_fit)) / sqrt(diag(loo_var))
  expected_cook_delta <- sweep(loo_fit, 2, full_fit, "-")
  expected_cook_delta <- -expected_cook_delta
  expected_cook <- rowSums(
    (expected_cook_delta %*% .symmetric_ginv(stats::cov(fit_samples))) *
      expected_cook_delta
  ) / 2

  summary <- .psis_fit_influence_summary(fit_samples, weights)
  expect_equal(summary[["full_fit"]], full_fit)
  expect_equal(summary[["loo_fit"]], loo_fit)
  expect_equal(summary[["loo_var"]], pmax(loo_var, 0))
  expect_equal(.dffits_internal(fit_samples, weights), expected_dffits)
  expect_equal(unname(.cooks.distance_internal(fit_samples, weights, P = 2)),
               expected_cook)
})

test_that("PSIS tau deletion helper aggregates remaining scale rows", {

  tau_samples <- matrix(
    c(
      1, 3, 5,
      2, 4, 6,
      3, 5, 7,
      4, 6, 8
    ),
    nrow  = 4,
    byrow = TRUE
  )
  weights <- matrix(
    c(
      .4, .3, .2, .1,
      .1, .2, .3, .4,
      .25, .25, .25, .25
    ),
    nrow  = 4,
    byrow = FALSE
  )

  deleted_tau <- cbind(
    rowMeans(tau_samples[, -1, drop = FALSE]),
    rowMeans(tau_samples[, -2, drop = FALSE]),
    rowMeans(tau_samples[, -3, drop = FALSE])
  )
  expected <- colSums(weights * deleted_tau)

  expect_equal(.influence_tau_del_from_samples(tau_samples, weights), expected)
})

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(influence_metafor_cases(), function(case) {
  test_that_case("Influence diagnostics match metafor", case, {
    expect_influence_matches_metafor(case)
  })
})

test_that("Influence tau.del for location-scale models uses aggregate deleted tau", {

  name <- "bangertdrowns2004_location-scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  tau_result <- .evaluate.brma.tau(
    fit           = fit_brma[["fit"]],
    scale_data    = fit_brma[["data"]][["scale"]],
    scale_formula = .create_fit_formula_list(data = fit_brma[["data"]], "scale"),
    scale_priors  = fit_brma[["priors"]][["scale"]],
    is_scale      = TRUE,
    is_multilevel = FALSE,
    K             = nrow(fit_brma[["data"]][["outcome"]])
  )
  tau_samples <- tau_result[["tau_total"]]
  weights     <- loo_weights(fit_brma)
  deleted_tau <- (rowSums(tau_samples) - tau_samples) / (ncol(tau_samples) - 1L)
  expected    <- colSums(weights * deleted_tau)

  inf_brma <- suppressWarnings(influence(fit_brma))
  scale_loo_coef_1 <- suppressWarnings(
    dfbetas(fit_brma, type = "scale", return_loo_estimates = TRUE)[[1]]
  )

  expect_equal(unname(inf_brma[["inf"]][["tau.del"]]), unname(expected),
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unname(inf_brma[["inf"]][["tau.del"]]),
                                unname(scale_loo_coef_1))))
})

test_that("Influence stats for selection model are available", {

  name <- "dat.lehmann2018-3PSM"
  skip_if_missing_fits(name)

  inf_brma <- suppressWarnings(influence(fits[[name]]))
  expect_true(!is.null(inf_brma))
})

test_that("Influence diagnostics and scalar outputs use study labels", {

  name <- "bcg_meta-analysis"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  labels   <- RoBMA:::.diagnostic_study_labels(fit_brma)
  inf_brma <- suppressWarnings(influence(fit_brma))

  expect_equal(rownames(inf_brma[["inf"]]), labels)
  expect_equal(rownames(inf_brma[["dfbs"]]), labels)
  expect_equal(names(covratio(fit_brma)), labels)
  expect_equal(names(dffits(fit_brma)), labels)
  expect_equal(names(cooks.distance(fit_brma)), labels)
})

test_that("DFFITS rejects unsupported model families", {

  model_names <- c(
    "bcg_glmm",
    "bcg_BMA.glmm",
    "dat.lehmann2018-3PSM",
    "dat.lehmann2018_RoBMA"
  )
  skip_if_missing_fits(model_names)

  expect_error(dffits(fits[["bcg_glmm"]]), "only available for normal outcome models")
  expect_error(dffits(fits[["bcg_BMA.glmm"]]), "only available for normal outcome models")
  expect_error(dffits(fits[["dat.lehmann2018-3PSM"]]), "not available for selection models")
  expect_error(dffits(fits[["dat.lehmann2018_RoBMA"]]), "not available for selection models")
})

test_that("Influence stats for model-averaging fits are internally consistent", {

  cases <- data.frame(
    name = c("dat.lehmann2018_BMA.norm", "bcg_BMA.glmm", "dat.lehmann2018_RoBMA"),
    unsupported_cook = c(NA, "normal outcome models", "selection models"),
    stringsAsFactors = FALSE
  )
  cases[["inf_cols"]] <- I(list(
    c("rstudent", "dffits", "cook.d", "cov.r", "tau.del", "hat"),
    c("rstudent", "cov.r", "tau.del"),
    c("rstudent", "cov.r", "tau.del")
  ))
  skip_if_missing_fits(cases[["name"]])

  for (i in seq_len(nrow(cases))) {
    name     <- cases[["name"]][[i]]
    fit_brma <- fits[[name]]
    inf_brma <- suppressWarnings(influence(fit_brma))

    expect_influence_object(inf_brma, nobs(fit_brma), cases[["inf_cols"]][[i]],
                            info = name)

    unsupported_cook <- cases[["unsupported_cook"]][[i]]
    if (!is.na(unsupported_cook)) {
      expect_false(any(c("dffits", "cook.d", "hat") %in% names(inf_brma[["inf"]])))
      expect_error(cooks.distance(fit_brma), unsupported_cook)
    }
  }
})
