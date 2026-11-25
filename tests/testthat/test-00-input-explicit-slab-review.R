context("Multivariate model averaging retains explicitly supplied slabs")

test_that("explicit multivariate mixture priors do not require an unused UISD", {

  dat <- data.frame(
    yi = c(-0.2, 0.1, 0.4, 0.5),
    study = factor(c("a", "a", "b", "b"))
  )
  slab <- BayesTools::prior(
    "normal", list(mean = 0, sd = 0.3),
    truncation = list(0, Inf), prior_weights = 3
  )
  arguments <- list(
    yi = dat$yi, vi = rep(0.04, 4), random = ~ 1 | study,
    data = dat, measure = "GEN", only_priors = TRUE,
    prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    prior_heterogeneity = slab
  )

  for (constructor in list(BMA.mv, RoBMA.mv)) {
    constructor_arguments <- arguments
    if (identical(constructor, RoBMA.mv)) {
      constructor_arguments[["prior_bias"]] <- FALSE
      constructor_arguments[["prior_bias_null"]] <- BayesTools::prior_none()
    }
    object <- do.call(constructor, constructor_arguments)
    allocation <- object[["priors"]][["random"]][["allocation"]][[1L]]

    expect_equal(allocation[["sd"]][["parameters"]], slab[["parameters"]])
    expect_equal(allocation[["sd"]][["truncation"]], slab[["truncation"]])
    expect_equal(mean(allocation[["inclusion"]][[1L]]), 0.75)

    with_uisd <- do.call(constructor, c(
      constructor_arguments, list(prior_unit_information_sd = 2)
    ))
    expect_equal(object[["priors"]], with_uisd[["priors"]])
  }
})

test_that("the fallback slab remains available for fixed-exclusion components", {

  dat <- data.frame(yi = c(0.1, 0.2), study = factor(c("a", "a")))
  arguments <- list(
    yi = dat$yi, vi = c(0.04, 0.09), random = ~ 1 | study,
    data = dat, measure = "GEN", only_priors = TRUE,
    prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    prior_heterogeneity = FALSE
  )

  expect_error(do.call(BMA.mv, arguments), "unit information sd", fixed = TRUE)
  object <- do.call(BMA.mv, c(arguments, list(prior_unit_information_sd = 1)))
  allocation <- object[["priors"]][["random"]][["allocation"]][[1L]]
  expect_true(BayesTools::is.prior.simple(allocation[["sd"]]))
  expect_equal(mean(allocation[["inclusion"]][[1L]]), 0)
})
