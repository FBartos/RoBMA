context("Model fitting for BMA.glmm")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_installed("metadat")
skip_refit_if_cached("BMA.glmm")


test_that("BMA.glmm fits binomial model (OR)", {
  data(dat.bcg, package = "metadat")

  fit <- BMA.glmm(
    ai = tpos, bi = tneg, ci = cpos, di = cneg,
    data = dat.bcg, measure = "OR",
    seed = 1, silent = TRUE
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit("bcg_BMA.glmm", fit)

  expect_s3_class(fit, "BMA.glmm")
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$mu))
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$tau))
  expect_null(fit$priors$outcome$bias)
  expect_no_error(summary(fit))
})


test_that("BMA.glmm handles custom priors", {
  data(dat.bcg, package = "metadat")

  fit <- BMA.glmm(
    ai = tpos, bi = tneg, ci = cpos, di = cneg,
    data = dat.bcg, measure = "OR",
    prior_effect = prior("normal", list(mean = 0, sd = 1)),
    prior_effect_null = prior("spike", list(location = 0)),
    prior_heterogeneity = prior("normal", list(mean = 0, sd = 0.5), truncation = list(lower = 0)),
    prior_heterogeneity_null = NULL,  # no null hypothesis for heterogeneity
    seed = 1, silent = TRUE
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit("bcg_BMA.glmm_custom", fit)

  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$mu))
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$tau))
  expect_equal(length(fit$priors$outcome$tau), 1)
})

test_that("BMA.glmm handles 3lvl location-scale meta-regression", {
  data(dat.bcg, package = "metadat")

  fit <- BMA.glmm(
    ai = tpos, bi = tneg, ci = cpos, di = cneg,
    mods = ~ year, scale = ~ year, cluster = alloc,
    data = dat.bcg, measure = "OR",
    prior_effect = prior("normal", list(mean = 0, sd = 1)),
    prior_effect_null = prior("spike", list(location = 0)),
    prior_heterogeneity = prior("normal", list(mean = 0, sd = 0.5), truncation = list(lower = 0)),
    prior_heterogeneity_null = NULL,  # no null hypothesis for heterogeneity
    seed = 1, silent = TRUE
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit("bcg_BMA.glmm_3lvl_location_scale", fit, info = list(mods = "year", scale = "year"))

  expect_equal(length(fit$priors$outcome$rho), 1)
  expect_equal(length(fit$priors$mods), 2)
  expect_equal(length(fit$priors$scale), 2)
})
