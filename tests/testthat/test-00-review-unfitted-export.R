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
    # the family class of the requests that an object cannot serve
    expect_error(as.data.frame(object), class = "RoBMA_unsupported_request")
  }
  # summary() keeps returning the data, and shows a prior-only object by its title
  expect_identical(summary(data_only), data_only[["data"]])
  expect_s3_class(summary(priors_only), "summary.brma")
  expect_identical(summary(priors_only)[["priors"]], priors_only[["priors"]])
})

test_that("summary of an only_priors object prints the title line of the release", {

  # 4.0.0 printed exactly the model title and no table for a prior-only model,
  # for example "Bayesian Random-Effects Model (k = 13)" for 13 estimates.
  dat <- data.frame(yi = seq(-.1, .2, length.out = 6), sei = .1, x = 1:6, study = rep(1:2, 3))
  prior_only <- function(fun = brma, ...) fun(
    yi = yi, sei = sei, data = dat, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE, ...)

  simple <- summary(prior_only())
  expect_s3_class(simple, "summary.brma")
  expect_identical(simple[["name"]], "Bayesian Random-Effects Model (k = 6)")
  expect_identical(capture.output(print(simple)),
                   c("", "Bayesian Random-Effects Model (k = 6)", ""))
  expect_identical(simple[["priors"]], prior_only()[["priors"]])

  titles <- c(
    "Bayesian Mixed-Effect Model (k = 6)" = summary(prior_only(mods = ~ x))[["name"]],
    "Bayesian Location-Scale Model (k = 6)" = summary(prior_only(scale = ~ x))[["name"]],
    "Bayesian Model-Averaged Random-Effects Model (k = 6)" = summary(prior_only(BMA))[["name"]],
    "Bayesian Multilevel Random-Effects Model (k = 6, clusters = 2)" =
      summary(brma(yi = yi, sei = sei, cluster = study, data = dat, measure = "GEN",
        prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE))[["name"]]
  )
  expect_identical(unname(titles), names(titles))
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
