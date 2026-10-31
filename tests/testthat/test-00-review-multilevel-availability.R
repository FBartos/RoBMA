test_that("multilevel sufficient-statistic underflow is unavailable rather than a zero effect", {

  expect_error(.evaluate.brma.multilevel_blup.norm(matrix(0, 1L, 1L),
    matrix(0, 1L, 1L), matrix(2^511, 1L, 1L), 2^-563, 2^1022, 1L),
    "Multilevel Gaussian conditional arithmetic is numerically unavailable.", fixed = TRUE)
  zero <- .evaluate.brma.multilevel_blup.norm(matrix(0, 1L, 1L),
    matrix(0, 1L, 1L), matrix(2^511, 1L, 1L), 0, 2^1022, 1L)
  expect_equal(zero$cluster, matrix(0, 1L, 1L), tolerance = 0)
  expect_equal(zero$estimate, matrix(0, 1L, 1L), tolerance = 0)
})


test_that("multilevel total conditional means use the stable weighted expression", {

  value <- .evaluate.brma.multilevel_blup.norm(matrix(0, 1L, 2L),
    matrix(1, 1L, 2L), matrix(1, 1L, 2L), c(1, -2^56), c(0, 1), c(1L, 1L))
  # The first fitted response has exactly zero sampling noise, so its total
  # latent conditional mean is y1 even when its two contributions nearly cancel.
  expect_equal(value$total[1L, 1L], 1, tolerance = 0)
  expect_equal(value$total[1L, 2L], .5 * (-2^56) + .5 * value$cluster[1L, 2L], tolerance = 0)
  positive <- .evaluate.brma.multilevel_blup.norm(matrix(0, 1L, 2L),
    matrix(1, 1L, 2L), matrix(1, 1L, 2L), c(1, -2^56), c(1e-20, 1), c(1L, 1L))
  cluster_mean <- (1 - 2^55) / 2.5
  expected <- 1 + 1e-20 * cluster_mean
  expect_lt(abs(positive$total[1L, 1L] - expected), 1e-14)
  object <- brma(yi = c(1, -2^56), vi = c(1e-20, 1), cluster = c("a", "a"),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    prior_effect = BayesTools::prior("point", list(location = 0)),
    prior_heterogeneity = BayesTools::prior("point", list(location = 1 / sqrt(.5))),
    prior_heterogeneity_allocation = BayesTools::prior("point", list(location = .5)))
  # This genuine constructor admits positive sampling variance. Zero sampling
  # noise is an internal conditional boundary control above, not public admission.
  public <- predict(object, type = "blup", quiet = TRUE,
    .posterior_samples = cbind(mu = 0, tau = 1 / sqrt(.5), rho = .5))
  expect_lt(abs(public[1L, 1L] - expected), 1e-14)
})
