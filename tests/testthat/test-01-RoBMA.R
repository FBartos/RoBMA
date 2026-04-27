context("Model fitting for RoBMA")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("RoBMA")


test_that("RoBMA handles default model", {
  data(dat.lehmann2018, package = "metadat")
  fit <- RoBMA(
    yi = yi, vi = vi,
    data = dat.lehmann2018, measure = "SMD",
    seed = 1, silent = TRUE
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018_RoBMA", fit)

  ## simple quick consistency checks
  # check that the model fits and has expected class
  expect_s3_class(fit, "RoBMA")

  # check that mixture priors are present
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$mu))
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$tau))
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$bias))

  # check that summary works
  expect_no_error(summary(fit))
})

test_that("RoBMA handles custom priors", {
  data(dat.lehmann2018, package = "metadat")
  fit <- RoBMA(
    yi = yi, vi = vi,
    data = dat.lehmann2018, measure = "SMD",
    prior_effect = prior("normal", list(mean = 0, sd = 0.5)),
    prior_effect_null = prior("spike", list(location = 0)),
    prior_heterogeneity = prior("normal", list(mean = 0, sd = 0.25), truncation = list(lower = 0)),
    prior_heterogeneity_null = NULL,  # no null hypothesis for heterogeneity
    seed = 1, silent = TRUE
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018_RoBMA_custom", fit)

  ## simple quick consistency checks
  # check that the model fits and has expected class
  expect_s3_class(fit, "RoBMA")

  # check that custom priors are applied
  expect_equal(length(fit$priors$outcome$mu), 2)
  expect_equal(length(fit$priors$outcome$tau), 1)
  expect_equal(fit$priors$outcome$mu[[2]]$parameters, list(mean = 0, sd = 0.5))
  expect_equal(fit$priors$outcome$tau[[1]]$parameters, list(mean = 0, sd = 0.25))
})

test_that("RoBMA handles meta-regression", {
  data(dat.lehmann2018, package = "metadat")
  fit <- suppressWarnings(RoBMA(
    yi = yi, vi = vi, mods = ~ Preregistered,
    data = dat.lehmann2018, measure = "SMD",
    seed = 1, silent = TRUE
  ))
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018_RoBMA_mods", fit, info = list(mods = c("Preregistered")))

  ## simple quick consistency checks
  # check that the model fits and has expected class
  expect_s3_class(fit, "RoBMA")

  # check that 2 mods parameters are present
  expect_equal(length(fit$priors$mods), 2)
})

test_that("RoBMA handles meta-regression with interaction", {
  data(dat.lehmann2018, package = "metadat")
  fit <- RoBMA(
    yi = yi, vi = vi, mods = ~ Preregistered * Gender,
    data = dat.lehmann2018, measure = "SMD",
    seed = 1, silent = TRUE
  )
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018_RoBMA_mods2", fit, info = list(mods = c("Preregistered", "Gender", "Preregistered:Gender")))

  ## simple quick consistency checks
  # check that the model fits and has expected class
  expect_s3_class(fit, "RoBMA")

  # check that 2 mods parameters are present
  expect_equal(length(fit$priors$mods), 4)
})

test_that("RoBMA handles multilevel location-scale meta-regression", {
  data(dat.lehmann2018, package = "metadat")
  fit <- suppressWarnings(RoBMA(
    yi = yi, vi = vi, mods = ~ Preregistered, scale = ~ Preregistered, cluster = Full_Citation,
    data = dat.lehmann2018, measure = "SMD",
    seed = 1, silent = TRUE, adapt = 1000
  ))
  fit <- suppressWarnings(add_loo(fit))
  save_fit("dat.lehmann2018_RoBMA_3lvl_mods_scale", fit, info = list(mods = c("Preregistered"), scale = c("Preregistered")))

  ## simple quick consistency checks
  # check that the model fits and has expected class
  expect_s3_class(fit, "RoBMA")

  # check that 2 mods parameters are present
  expect_equal(length(fit$priors$mods), 2)
  expect_equal(length(fit$priors$scale), 2)
  expect_equal(length(fit$priors$outcome$rho), 1)
})
