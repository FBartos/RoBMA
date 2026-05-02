context("Model fitting for BMA.norm")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("BMA.norm")


test_that("BMA.norm handles default model", {
  data(dat.lehmann2018, package = "metadat")
  fit <- BMA.norm(
    yi = yi, vi = vi,
    data = dat.lehmann2018, measure = "SMD",
    seed = 1, silent = TRUE
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018_BMA.norm", fit)

  expect_s3_class(fit, "BMA.norm")
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$mu))
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$tau))
  expect_null(fit$priors$outcome$bias)
  expect_no_error(summary(fit))
})

test_that("BMA.norm handles custom priors", {
  data(dat.lehmann2018, package = "metadat")
  fit <- BMA.norm(
    yi = yi, vi = vi,
    data = dat.lehmann2018, measure = "SMD",
    prior_effect = prior("normal", list(mean = 0, sd = 0.5)),
    prior_effect_null = prior("spike", list(location = 0)),
    prior_heterogeneity = prior("normal", list(mean = 0, sd = 0.25), truncation = list(lower = 0)),
    prior_heterogeneity_null = NULL,  # no null hypothesis for heterogeneity
    seed = 1, silent = TRUE
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018_BMA.norm_custom", fit)

  expect_s3_class(fit, "BMA.norm")
  expect_equal(length(fit$priors$outcome$mu), 2)
  expect_equal(length(fit$priors$outcome$tau), 1)
  expect_equal(fit$priors$outcome$mu[[2]]$parameters, list(mean = 0, sd = 0.5))
  expect_equal(fit$priors$outcome$tau[[1]]$parameters, list(mean = 0, sd = 0.25))
})

test_that("BMA.norm handles meta-regression", {
  data(dat.lehmann2018, package = "metadat")
  fit <- BMA.norm(
    yi = yi, vi = vi, mods = ~ Preregistered,
    data = dat.lehmann2018, measure = "SMD",
    seed = 1, silent = TRUE
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018_BMA.norm_mods", fit, info = list(mods = c("Preregistered")))

  expect_s3_class(fit, "BMA.norm")
  expect_equal(length(fit$priors$mods), 2)
})

test_that("BMA.norm handles scale-regression", {
  data(dat.lehmann2018, package = "metadat")
  fit <- suppressWarnings(BMA.norm(
    yi = yi, vi = vi, scale = ~ Preregistered,
    data = dat.lehmann2018, measure = "SMD",
    seed = 1, silent = TRUE
  )) # suppress warning about removing spike at heterogeneity
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018_BMA.norm_scale", fit, info = list(scale = c("Preregistered")))

  expect_s3_class(fit, "BMA.norm")
  expect_equal(length(fit$priors$scale), 2)
})
