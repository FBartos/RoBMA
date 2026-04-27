context("Prior input handling for brma.glmm")

skip_on_cran()

test_data_bin <- data.frame(
  ai = c(4L, 6L, 8L),
  bi = c(16L, 14L, 12L),
  ci = c(3L, 5L, 7L),
  di = c(17L, 15L, 13L)
)

test_data_pois <- data.frame(
  x1i = c(4L, 6L, 8L),
  x2i = c(3L, 5L, 7L),
  t1i = c(20, 25, 30),
  t2i = c(21, 24, 31)
)


test_that("Binomial GLMM baserate priors are assigned", {

  result_default <- brma.glmm(
    ai = ai, bi = bi, ci = ci, di = di,
    data = test_data_bin, measure = "OR",
    prior_baserate = NULL, only_priors = TRUE
  )[["priors"]]

  expect_equal(result_default$outcome$pi$distribution, "beta")
  expect_equal(result_default$outcome$pi$parameters$alpha, 1)
  expect_equal(result_default$outcome$pi$parameters$beta,  1)

  custom_prior <- BayesTools::prior("beta", parameters = list(alpha = 2, beta = 3))
  result_custom <- brma.glmm(
    ai = ai, bi = bi, ci = ci, di = di,
    data = test_data_bin, measure = "OR",
    prior_baserate = custom_prior, only_priors = TRUE
  )[["priors"]]

  expect_equal(result_custom$outcome$pi$parameters$alpha, 2)
  expect_equal(result_custom$outcome$pi$parameters$beta,  3)
})


test_that("Poisson GLMM lograte priors are assigned", {

  result_default <- brma.glmm(
    x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
    data = test_data_pois, measure = "IRR",
    prior_lograte = NULL, only_priors = TRUE
  )[["priors"]]

  expect_equal(result_default$outcome$phi$distribution, "normal")

  custom_prior <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 2))
  result_custom <- brma.glmm(
    x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
    data = test_data_pois, measure = "IRR",
    prior_lograte = custom_prior, only_priors = TRUE
  )[["priors"]]

  expect_equal(result_custom$outcome$phi$parameters$sd, 2)
})
