context("Cached zero-ordinate marginal-means fit")

source(testthat::test_path("common-functions.R"))


test_that("zero-ordinate marginal means fit is prepared for cached post-fit tests", {

  skip_on_cran()
  skip_if_fit_not_active("marginal_means_zero_ordinate")
  skip_refit_if_cached("marginal_means_zero_ordinate")
  set.seed(3)
  k    <- 30L
  data <- data.frame(x = stats::rnorm(k), sei = stats::runif(k, 0.1, 0.3))
  data[["yi"]] <- stats::rnorm(k, 0.3, data[["sei"]])
  # The effect prior is truncated at 0: the intercept (and the mean at the
  # mean of 'x') has a zero prior density at -0.5, the means at -1 and +1 SD
  # of 'x' do not.
  fit <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
    prior_effect = prior("normal", list(0, 1), truncation = list(0, Inf)),
    chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
    silent = TRUE
  ))
  save_fit("marginal_means_zero_ordinate", fit)
  expect_true(is_cached_fit_valid("marginal_means_zero_ordinate"))
})
