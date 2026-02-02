context("Model fitting for BMA.norm")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("BMA.norm")


test_that("BMA.norm fits model without bias adjustment", {
  data(Bem2011, package = "RoBMA")

  fit <- BMA.norm(
    yi = d, sei = se,
    data = Bem2011, measure = "SMD",
    seed = 1, silent = TRUE
  )
  save_fit("bem2011_BMA.norm", fit)

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

  # check that custom priors are applied
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$mu))
  expect_true(BayesTools::is.prior.mixture(fit$priors$outcome$tau))

  # tau mixture should only have alternative (no null)
  expect_equal(length(fit$priors$outcome$tau), 1)
})
