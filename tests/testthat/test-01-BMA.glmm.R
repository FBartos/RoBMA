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
  save_fit("bcg_BMA.glmm", fit)

  # check that the model fits and has the expected class
  expect_s3_class(fit, "BMA.glmm")

  # check that mixture priors are present
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$mu))
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$tau))

  # check that no bias prior is present
  expect_null(fit$priors$outcome$bias)

  # check that summary works
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
  save_fit("bcg_BMA.glmm_custom", fit)

  # check that custom priors are applied
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$mu))
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$tau))

  # tau mixture should only have alternative (no null)
  expect_equal(length(fit$priors$outcome$tau), 1)
})

# TODO: add single group models
# TODO: add single group models
