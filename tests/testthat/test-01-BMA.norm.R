context("Model fitting for BMA.norm")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("BMA.norm")


test_that("BMA.norm handles default model", {
  data(Bem2011, package = "RoBMA")

  fit <- BMA.norm(
    yi = d, sei = se,
    data = Bem2011, measure = "SMD",
    seed = 1, silent = TRUE
  )
  save_fit("bem2011_BMA.norm", fit)

  ## simple quick consistency checks
  # check that the model fits and has expected class
  expect_s3_class(fit, "BMA.norm")

  # check that mixture priors are present
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$mu))
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$tau))

  # check that no bias prior is present
  expect_null(fit$priors$outcome$bias)

  # check that summary works
  expect_no_error(summary(fit))
})

test_that("BMA.norm handles custom priors", {
  data(Bem2011, package = "RoBMA")

  fit <- BMA.norm(
    yi = d, sei = se,
    data = Bem2011, measure = "SMD",
    prior_effect = prior("normal", list(mean = 0, sd = 0.5)),
    prior_effect_null = prior("spike", list(location = 0)),
    prior_heterogeneity = prior("normal", list(mean = 0, sd = 0.25), truncation = list(lower = 0)),
    prior_heterogeneity_null = NULL,  # no null hypothesis for heterogeneity
    seed = 1, silent = TRUE
  )
  save_fit("bem2011_BMA.norm_custom", fit)

  ## simple quick consistency checks
  # check that the model fits and has expected class
  expect_s3_class(fit, "BMA.norm")

  # check that custom priors are applied
  expect_equal(length(fit$priors$outcome$mu), 2)
  expect_equal(length(fit$priors$outcome$tau), 1)
  expect_equal(fit$priors$outcome$mu[[2]]$parameters, list(mean = 0, sd = 0.5))
  expect_equal(fit$priors$outcome$tau[[1]]$parameters, list(mean = 0, sd = 0.25))
})

test_that("BMA.norm handles meta-regression", {
  data(Bem2011, package = "RoBMA")

  fit <- BMA.norm(
    yi = d, sei = se, mods = ~ se,
    data = Bem2011, measure = "SMD",
    seed = 1, silent = TRUE
  )
  save_fit("bem2011_BMA.norm_mods", fit)

  ## simple quick consistency checks
  # check that the model fits and has expected class
  expect_s3_class(fit, "BMA.norm")

  # check that 2 mods parameters are present
  expect_equal(length(fit$priors$mods), 2)
})

test_that("BMA.norm handles scale-regression", {
  data(Bem2011, package = "RoBMA")
  # TODO: fix error
  fit <- suppressWarnings(BMA.norm(
    yi = d, sei = se, scale = ~ se,
    data = Bem2011, measure = "SMD",
    seed = 1, silent = TRUE
  )) # suppress warning about removing spike at heterogeneity
  save_fit("bem2011_BMA.norm_scale", fit)

  ## simple quick consistency checks
  # check that the model fits and has expected class
  expect_s3_class(fit, "BMA.norm")

  # check that 2 mods parameters are present
  expect_equal(length(fit$priors$scale), 2)
})
