context("Prior input handling for BMA product-space constructors")

skip_on_cran()

.get_is_null <- function(prior) {
  components <- attr(prior, "components")
  if (is.null(components)) return(NULL)
  return(components == "null")
}

.get_prior_weights <- function(prior) {
  return(vapply(prior, function(x) x[["prior_weights"]], numeric(1)))
}

bma_test_data <- data.frame(
  effect    = c(0.10, 0.25, 0.15, 0.30, 0.05),
  std_err   = sqrt(c(0.04, 0.06, 0.05, 0.08, 0.03)),
  n         = c(50L, 75L, 60L, 40L, 90L),
  mod_cont  = c(1.5, 2.3, 1.8, 3.1, 0.9),
  scale_var = c(0.5, 1.0, 0.8, 1.2, 0.6),
  stringsAsFactors = FALSE
)

bma_glmm_data <- data.frame(
  ai = c(4L, 7L, 3L, 9L, 6L),
  bi = c(46L, 68L, 57L, 31L, 84L),
  ci = c(8L, 6L, 5L, 7L, 9L),
  di = c(42L, 69L, 55L, 33L, 81L),
  stringsAsFactors = FALSE
)


test_that("BMA alias and staged normal constructor return product-space objects", {

  object_alias <- BMA(
    yi = effect, sei = std_err, data = bma_test_data,
    measure = "SMD", only_priors = TRUE
  )
  object_norm <- BMA.norm(
    yi = effect, sei = std_err, data = bma_test_data,
    measure = "SMD", only_priors = TRUE
  )

  expect_s3_class(object_alias, "BMA.norm")
  expect_s3_class(object_alias, "RoBMA")
  expect_s3_class(object_alias, "brma")
  expect_identical(class(object_alias), class(object_norm))
  expect_null(object_alias[["priors"]][["outcome"]][["bias"]])
  expect_true(BayesTools::is.prior.mixture(object_alias[["priors"]][["outcome"]][["mu"]]))
  expect_true(BayesTools::is.prior.mixture(object_alias[["priors"]][["outcome"]][["tau"]]))

  object_data <- BMA(
    yi = effect, sei = std_err, data = bma_test_data,
    measure = "SMD", only_data = TRUE
  )

  expect_s3_class(object_data, "BMA.norm")
  expect_true("data" %in% names(object_data))
  expect_null(object_data[["priors"]])
  expect_null(object_data[["fit"]])
})


test_that("BMA.glmm staged constructor returns product-space objects without bias priors", {

  object_data <- BMA.glmm(
    ai = ai, bi = bi, ci = ci, di = di,
    data = bma_glmm_data, measure = "OR",
    only_data = TRUE
  )
  object_priors <- BMA.glmm(
    ai = ai, bi = bi, ci = ci, di = di,
    data = bma_glmm_data, measure = "OR",
    only_priors = TRUE
  )

  expect_s3_class(object_data, "BMA.glmm")
  expect_s3_class(object_data, "RoBMA")
  expect_s3_class(object_data, "brma.glmm")
  expect_s3_class(object_priors, "BMA.glmm")
  expect_null(object_priors[["priors"]][["outcome"]][["bias"]])
  expect_true(BayesTools::is.prior.mixture(object_priors[["priors"]][["outcome"]][["mu"]]))
  expect_true(BayesTools::is.prior.mixture(object_priors[["priors"]][["outcome"]][["tau"]]))
})


test_that("BMA null prior NULL and FALSE omit null components", {

  result_null <- BMA(
    yi = effect, sei = std_err, data = bma_test_data,
    prior_effect_null = NULL,
    prior_heterogeneity_null = NULL,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]
  result_false <- BMA(
    yi = effect, sei = std_err, data = bma_test_data,
    prior_effect_null = FALSE,
    prior_heterogeneity_null = FALSE,
    measure = "SMD", only_priors = TRUE
  )[["priors"]]

  expect_equal(length(result_null[["outcome"]][["mu"]]), 1)
  expect_equal(length(result_null[["outcome"]][["tau"]]), 1)
  expect_equal(sum(.get_is_null(result_null[["outcome"]][["mu"]])), 0)
  expect_equal(sum(.get_is_null(result_null[["outcome"]][["tau"]])), 0)

  expect_equal(length(result_false[["outcome"]][["mu"]]), 1)
  expect_equal(length(result_false[["outcome"]][["tau"]]), 1)
  expect_equal(sum(.get_is_null(result_false[["outcome"]][["mu"]])), 0)
  expect_equal(sum(.get_is_null(result_false[["outcome"]][["tau"]])), 0)
})


test_that("BMA keeps custom prior weights in mixture constructors", {

  effect_alt <- list(
    BayesTools::prior(
      "normal", parameters = list(mean = 0, sd = 0.25),
      prior_weights = 3
    ),
    BayesTools::prior(
      "normal", parameters = list(mean = 0, sd = 0.75),
      prior_weights = 2
    )
  )
  effect_null <- BayesTools::prior(
    "spike", parameters = list(location = 0),
    prior_weights = 5
  )

  result <- BMA(
    yi = effect, sei = std_err, data = bma_test_data,
    prior_effect = effect_alt,
    prior_effect_null = effect_null,
    measure = "SMD", only_priors = TRUE
  )[["priors"]][["outcome"]][["mu"]]

  expect_equal(length(result), 3)
  expect_equal(.get_prior_weights(result), c(5, 3, 2))
  expect_equal(.get_is_null(result), c(TRUE, FALSE, FALSE))
})


test_that("BMA formula null prior settings apply to moderator and scale terms", {

  result <- suppressWarnings(BMA(
    yi = effect, sei = std_err,
    mods = ~ mod_cont, scale = ~ scale_var,
    data = bma_test_data, measure = "SMD",
    prior_mods_null  = list(mod_cont = NULL),
    prior_scale_null = list(scale_var = FALSE),
    only_priors = TRUE
  )[["priors"]])

  expect_true(BayesTools::is.prior.mixture(result[["mods"]][["mod_cont"]]))
  expect_true(BayesTools::is.prior.mixture(result[["scale"]][["scale_var"]]))
  expect_equal(sum(.get_is_null(result[["mods"]][["mod_cont"]])), 0)
  expect_equal(sum(.get_is_null(result[["scale"]][["scale_var"]])), 0)
})


test_that("BMA constructors reject bias and model arguments passed through dots", {

  expect_error(
    BMA(
      yi = effect, sei = std_err, data = bma_test_data,
      measure = "SMD", model_type = "PSMA", only_priors = TRUE
    ),
    "Unused argument.*model_type"
  )

  expect_error(
    BMA(
      yi = effect, sei = std_err, data = bma_test_data,
      measure = "SMD", prior_bias = BayesTools::prior_none(),
      only_priors = TRUE
    ),
    "Unused argument.*prior_bias"
  )

  expect_error(
    BMA.glmm(
      ai = ai, bi = bi, ci = ci, di = di,
      data = bma_glmm_data, measure = "OR",
      prior_bias = BayesTools::prior_none(),
      only_priors = TRUE
    ),
    "Unused argument.*prior_bias"
  )
})
