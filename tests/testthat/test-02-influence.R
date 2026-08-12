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
  # Meta-analytic Cook's distance is the unscaled squared Mahalanobis distance.
  expected_cook <- rowSums(
    (expected_cook_delta %*% .symmetric_ginv(stats::cov(fit_samples))) *
      expected_cook_delta
  )

  expect_equal(.dffits_internal(fit_samples, weights), expected_dffits)
  expect_equal(unname(.cooks.distance_internal(fit_samples, weights)),
               expected_cook)
})

test_that("influence diagnostics retain tiny non-zero posterior variation", {

  fit_samples <- matrix(
    c(
      0, 1, 3,
      1, 2, 1,
      2, 4, 5,
      4, 3, 2,
      7, 6, 8
    ),
    nrow  = 5,
    byrow = TRUE,
    dimnames = list(NULL, paste0("study", 1:3))
  )
  weights <- cbind(
    c(.40, .25, .15, .10, .10),
    c(.10, .20, .30, .25, .15),
    rep(.20, 5)
  )

  scales      <- c(2^-500, 2^-20, 2^20)
  offsets     <- c(0, 2^20, -2^40)
  transformed <- sweep(fit_samples, 2, scales, "*")
  transformed <- sweep(transformed, 2, offsets, "+")

  expect_equal(
    .dffits_internal(transformed, weights),
    .dffits_internal(fit_samples, weights),
    tolerance = 1e-9
  )
  expect_equal(
    unname(.cooks.distance_internal(transformed, weights)),
    unname(.cooks.distance_internal(fit_samples, weights)),
    tolerance = 1e-9
  )

  dfbetas_base <- .dfbetas_internal(fit_samples[, 1:2], weights)
  dfbetas_tiny <- .dfbetas_internal(transformed[, 1:2], weights)
  expect_equal(
    dfbetas_tiny[["values"]],
    dfbetas_base[["values"]],
    tolerance = 1e-9
  )

  with_constant <- cbind(transformed[, 1:2], fixed = 1)
  dfbetas_fixed <- .dfbetas_internal(with_constant, weights)
  expect_true(all(is.nan(dfbetas_fixed[["values"]][, "fixed"])))
  expect_identical(dfbetas_fixed[["undefined"]][, "fixed"], rep(TRUE, 3))
})

test_that("COVRATIO is invariant to parameter location and scale", {

  samples <- cbind(
    beta_1 = c(0, 1, 2, 4, 7),
    beta_2 = c(3, 1, 5, 2, 8)
  )
  weights <- cbind(
    c(.40, .25, .15, .10, .10),
    c(.10, .20, .30, .25, .15),
    rep(.20, 5)
  )
  transformed <- sweep(samples, 2, c(2^-500, 2^-20), "*")
  transformed <- sweep(transformed, 2, c(0, 2^20), "+")

  base <- .covratio_internal(samples, weights)
  tiny <- .covratio_internal(transformed, weights)

  expect_equal(tiny[["values"]], base[["values"]], tolerance = 1e-9)
  expect_identical(tiny[["excluded"]], c(beta_1 = FALSE, beta_2 = FALSE))

  fixed <- .covratio_internal(cbind(transformed, fixed = 1), weights)
  expect_identical(fixed[["excluded"]][["fixed"]], TRUE)
  expect_equal(fixed[["values"]], tiny[["values"]], tolerance = 1e-12)
})

test_that("symmetric inverse uses a relative numerical-rank threshold", {

  x <- diag(c(1, 1e-12))
  expect_equal(.symmetric_ginv(x), diag(c(1, 1e12)), tolerance = 1e-12)
})

test_that("PSIS tau deletion helper uses within-draw RMS", {

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
    sqrt(rowMeans(tau_samples[, -1, drop = FALSE]^2)),
    sqrt(rowMeans(tau_samples[, -2, drop = FALSE]^2)),
    sqrt(rowMeans(tau_samples[, -3, drop = FALSE]^2))
  )
  expected <- colSums(weights * deleted_tau)
  arithmetic <- cbind(
    rowMeans(tau_samples[, -1, drop = FALSE]),
    rowMeans(tau_samples[, -2, drop = FALSE]),
    rowMeans(tau_samples[, -3, drop = FALSE])
  )

  expect_equal(.influence_tau_del_from_samples(tau_samples, weights), expected)
  expect_false(isTRUE(all.equal(expected, colSums(weights * arithmetic))))
})

test_that("influence mv branch uses fixed-location availability guard", {

  model <- structure(
    list(data = structure(list(), outcome_type = "norm"), priors = list()),
    class = c("brma.mv", "brma.norm", "brma")
  )

  testthat::local_mocked_bindings(
    .outcome_type = function(model) {
      "norm"
    },
    .is_weightfunction = function(model) {
      TRUE
    },
    .package = "RoBMA"
  )

  expect_error(
    influence(model),
    "selection models"
  )
})

test_that("Standalone DFBETAS and COVRATIO check Pareto k diagnostics", {

  model <- structure(
    list(
      data = list(
        outcome = data.frame(
          yi   = c(0, 1),
          sei  = c(1, 1),
          slab = c("a", "b")
        )
      ),
      fit = list()
    ),
    class = "brma"
  )
  attr(model[["data"]], "outcome_type") <- "norm"

  psis_weights <- matrix(
    c(
      0.2, 0.5, 0.3,
      0.5, 0.3, 0.2
    ),
    nrow = 3
  )
  samples <- matrix(c(0, 1, 2), ncol = 1)
  colnames(samples) <- "(mu) intercept"
  fit_samples <- matrix(
    c(
      0, 1,
      1, 2,
      2, 3
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b"))
  )
  checked       <- 0L
  context_calls <- 0L

  testthat::local_mocked_bindings(
    .diagnostic_check_loo = function(model, context = NULL, unit = "estimate") {

      checked <<- checked + 1L
      invisible(NULL)
    },
    .diagnostic_psis_context = function(model, context = NULL) {

      context_calls <<- context_calls + 1L
      list(psis_weights = psis_weights)
    },
    .diagnostic_psis_weights = function(model, weights = NULL) {

      if (!is.null(weights)) {
        return(weights)
      }
      psis_weights
    },
    .diagnostic_location_parameter_samples = function(model,
                                                      standardized_coefficients = FALSE,
                                                      transform_factors = TRUE) {

      samples
    },
    .influence_fit_samples = function(model) {

      fit_samples
    },
    .get_model_matrix = function(model) {

      cbind(intercept = c(1, 1), slope = c(0, 1))
    },
    .get_estimate_labels = function(model) {

      c("a", "b")
    },
    .is_scale = function(model) {

      FALSE
    },
    .is_bias = function(model) {

      FALSE
    },
    .package = "RoBMA"
  )

  expect_s3_class(dfbetas(model), "dfbetas.brma")
  expect_length(covratio(model), 2L)
  expect_length(dffits(model), 2L)
  expect_length(cooks.distance(model), 2L)
  expect_identical(checked, 4L)
  expect_identical(context_calls, 4L)

  dfbetas(model, .weights = psis_weights)
  covratio(model, .weights = psis_weights)
  expect_identical(checked, 4L)
  expect_identical(context_calls, 4L)
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

test_that("Influence tau.del for location-scale models uses deleted-row RMS", {

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
  deleted_tau <- sqrt(
    (rowSums(tau_samples^2) - tau_samples^2) / (ncol(tau_samples) - 1L)
  )
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


test_that("Influence rejects GLMMs without a discrete PIT convention", {

  model_names <- c("bcg_glmm", "bcg_BMA.glmm")
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    expect_error(
      influence(fits[[name]]),
      "discrete PIT convention",
      info = name
    )
  }
})


test_that("Influence stats for model-averaging fits are internally consistent", {

  cases <- data.frame(
    name = c("dat.lehmann2018_BMA.norm", "dat.lehmann2018_RoBMA"),
    unsupported_cook = c(NA, "selection models"),
    stringsAsFactors = FALSE
  )
  cases[["inf_cols"]] <- I(list(
    c("rstudent", "dffits", "cook.d", "cov.r", "tau.del", "hat"),
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

test_that("Influence for brma.mv returns implemented estimate-unit components", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  inf_brma <- suppressWarnings(influence(fit_brma))
  dffits_mv <- dffits(fit_brma)
  cooks_mv  <- cooks.distance(fit_brma)

  expect_s3_class(inf_brma, "infl.brma")
  expect_named(inf_brma[["inf"]], c("rstudent", "dffits", "cook.d", "cov.r", "hat"))
  expect_equal(nrow(inf_brma[["inf"]]), nobs(fit_brma))
  expect_true(is.data.frame(inf_brma[["dfbs"]]))
  expect_equal(nrow(inf_brma[["dfbs"]]), nobs(fit_brma))
  expect_equal(ncol(inf_brma[["dfbs"]]), 1)
  expect_equal(
    unname(inf_brma[["inf"]][["dffits"]]),
    unname(as.numeric(dffits_mv)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(inf_brma[["inf"]][["cook.d"]]),
    unname(as.numeric(cooks_mv)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(inf_brma[["inf"]][["cov.r"]]),
    unname(as.numeric(covratio(fit_brma))),
    tolerance = 1e-12
  )
  expect_equal(
    unname(inf_brma[["inf"]][["hat"]]),
    unname(hatvalues(fit_brma)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(inf_brma[["dfbs"]])),
    unname(as.matrix(dfbetas(fit_brma))),
    tolerance = 1e-12
  )
  expect_null(attr(inf_brma, "note"))
  expect_equal(
    attr(dffits_mv, "brma_mv_target")[["reported_target"]],
    "fixed_location_fitted_value"
  )
  expect_equal(
    attr(cooks_mv, "brma_mv_target")[["reported_target"]],
    "fixed_location_fitted_value"
  )
})

test_that("Influence diagnostics support random-formula brma.mv fixed coefficients", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  dfbs     <- dfbetas(fit_brma)
  covr     <- covratio(fit_brma)
  dff      <- dffits(fit_brma)
  cook     <- cooks.distance(fit_brma)

  expect_s3_class(dfbs, "dfbetas.brma")
  expect_equal(nrow(dfbs), nobs(fit_brma))
  expect_equal(colnames(dfbs), "(mu) intercept")
  expect_named(covr, .diagnostic_study_labels(fit_brma))
  expect_named(dff, .diagnostic_study_labels(fit_brma))
  expect_named(cook, .diagnostic_study_labels(fit_brma))
  expect_true(all(is.finite(covr)))
  expect_true(all(is.finite(dff)))
  expect_true(all(is.finite(cook)))
  expect_null(attr(dfbs, "note"))
  expect_null(attr(covr, "note"))
  expect_equal(attr(dff, "brma_mv_target")[["reported_target"]],
               "fixed_location_fitted_value")
  expect_equal(attr(cook, "brma_mv_target")[["reported_target"]],
               "fixed_location_fitted_value")
})

test_that("v14 brma.mv influence diagnostics return finite estimate-unit output", {

  skip_if_not_certification("This case exercises the high-draw v14 fixtures.")
  mv_names <- c(
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_assink2016_nested",
    "brma.mv_v14_ishak2007_har",
    "brma.mv_v14_begg1989_study_treatment"
  )
  skip_if_missing_fits(mv_names)

  for (name in mv_names) {
    fit_brma <- fits[[name]]
    inf_brma <- suppressWarnings(influence(fit_brma))
    dfbs     <- .with_expected_pareto_warnings_muffled(stats::dfbetas(fit_brma))
    dff      <- suppressWarnings(dffits(fit_brma))
    cook     <- suppressWarnings(cooks.distance(fit_brma))
    covr     <- .with_expected_pareto_warnings_muffled(covratio(fit_brma))
    inf_mat  <- as.matrix(inf_brma[["inf"]])
    dfbs_mat <- as.matrix(inf_brma[["dfbs"]])

    expect_s3_class(inf_brma, "infl.brma")
    expect_equal(nrow(inf_brma[["inf"]]), nobs(fit_brma), info = name)
    expect_equal(nrow(inf_brma[["dfbs"]]), nobs(fit_brma), info = name)
    expect_equal(nrow(dfbs), nobs(fit_brma), info = name)
    expect_equal(length(dff), nobs(fit_brma), info = name)
    expect_equal(length(cook), nobs(fit_brma), info = name)
    expect_equal(length(covr), nobs(fit_brma), info = name)
    expect_true(all(is.finite(inf_mat[, setdiff(colnames(inf_mat), "cov.r"), drop = FALSE])),
                info = name)
    expect_true(all(is.finite(dfbs_mat)), info = name)
    expect_false(any(grepl("rho|var_frac|var_ratio|sd\\(", colnames(dfbs_mat))),
                 info = name)
    expect_true(all(is.finite(dff)), info = name)
    expect_true(all(is.finite(cook)), info = name)
    expect_true(all(is.finite(as.matrix(dfbs))), info = name)
    expect_true(all(is.finite(inf_mat[, "cov.r"])), info = name)
    expect_true(all(is.finite(covr)), info = name)
    expect_null(attr(covr, "note"), info = name)
    expect_null(attr(inf_brma, "note"), info = name)
  }
})
