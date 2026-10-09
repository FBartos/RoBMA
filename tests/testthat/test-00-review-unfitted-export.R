context("Data-frame export of unfitted constructor objects")

test_that("as.data.frame refuses objects without a fitted model", {

  dat <- data.frame(yi = seq(-.1, .2, length.out = 6), sei = .1)
  data_only <- brma(
    yi = yi, sei = sei, data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_data = TRUE, silent = TRUE
  )
  priors_only <- brma(
    yi = yi, sei = sei, data = dat, measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
    only_priors = TRUE, silent = TRUE
  )

  for (object in list(data_only, priors_only)) {
    expect_error(as.data.frame(object), class = "RoBMA_unfitted_object")
  }
  # summary() keeps returning the data or the priors
  expect_identical(summary(data_only), data_only[["data"]])
  expect_identical(summary(priors_only), priors_only[["priors"]])
})

test_that("as.data.frame exports a fitted object", {

  object <- review_factor_replay(
    c("A", "B"), "treatment", independent = TRUE, distribution = "normal"
  )
  object[["summary"]] <- .object_summary(object)

  exported <- as.data.frame(object, include_mcmc_diagnostics = FALSE)

  expect_s3_class(exported, "data.frame")
  expect_gt(nrow(exported), 0L)
  expect_true(all(c("component", "parameter") %in% names(exported)))
})
