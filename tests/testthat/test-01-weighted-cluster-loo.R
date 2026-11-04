source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("weighted-cluster-loo")

test_that("the weighted Gaussian cluster LOO fixture is fitted", {

  data <- data.frame(
    yi = c(.1, -.1, .3, .5, .2, .4, -.2, .1, .2, .4, .3, .6),
    sei = rep(c(.2, .3, .25), 4L),
    weights = c(.5, 1, 2, 1.5, .75, 1, 2, 1, .5, 1, 1.25, .8),
    study = rep(1:4, each = 3L)
  )
  object <- brma(yi = yi, sei = sei, weights = weights, cluster = study, data = data,
    measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("point", list(.4)),
    prior_heterogeneity_allocation = BayesTools::prior("point", list(.5)),
    chains = 2, adapt = 500, burnin = 500, sample = 3000,
    seed = 724, parallel = FALSE, silent = TRUE)
  object <- add_loo(object)
  object <- add_waic(object)
  save_fit("weighted_gaussian_cluster_loo", object, info = list(data = data))
  expect_true(is_cached_fit_valid("weighted_gaussian_cluster_loo"))
})
