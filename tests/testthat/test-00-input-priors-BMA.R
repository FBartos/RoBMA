context("Prior input handling for BMA product-space constructors")

test_that("sole ordered term components retain their total-owned mixture and bind", {

  dat <- data.frame(yi = c(0.1, 0.2, 0.3, 0.4), sei = rep(0.2, 4), f = ordered(c("a", "b", "c", "a")))
  normal <- BayesTools::prior("normal", list(0, 1), prior_weights = 3)
  zero <- BayesTools::prior("point", list(0), prior_weights = 1)
  priors <- list(
    ordinary = prior_ordered(total = normal, allocation = c(0.25, 0.75), id = "retained"),
    mixed = prior_ordered(total = BayesTools::prior_mixture(list(zero, normal), is_null = c(TRUE, FALSE)),
      allocation = c(0.25, 0.75), id = "retained")
  )
  for (constructor in list(BMA, BMA.norm, RoBMA)) for (prior in priors) {
    expect_silent(object <- constructor(
      yi = yi, sei = sei, mods = ~ f, data = dat, measure = "SMD",
      prior_mods = list(f = prior), prior_mods_null = list(f = NULL), only_priors = TRUE, silent = TRUE
    ))
    expect_s3_class(object[["priors"]][["mods"]][["f"]], "prior.ordered")
    expect_false(BayesTools::is.prior.mixture(object[["priors"]][["mods"]][["f"]]))
    args <- .create_jags_formula_args(object[["data"]], object[["priors"]])
    expect_silent(compiled <- BayesTools::JAGS_formula(
      formula = args$formula_list$mu, parameter = "mu", data = args$formula_data_list$mu,
      prior_list = args$formula_prior_list$mu, formula_scale = args$formula_scale_list$mu
    ))
    bound <- compiled[["prior_list"]][["mu_f"]]
    for (field in c("total", "allocation", "contrast", "id", "prior_weights")) {
      expect_identical(bound[[field]], prior[[field]])
    }
    expect_identical(attr(bound[["total"]], "components"), attr(prior[["total"]], "components"))
    expect_identical(attr(bound[["total"]], "prior_weights"), attr(prior[["total"]], "prior_weights"))
    expect_error(constructor(
      yi = yi, sei = sei, mods = ~ f, data = dat, measure = "SMD",
      prior_mods = list(f = prior), only_priors = TRUE, silent = TRUE
    ), "prior_ordered(total = ...)", class = "RoBMA_prior_error", fixed = TRUE)
  }
})

test_that("a factor term whose null and alternative contrasts differ is refused with its remedies", {

  dat <- data.frame(yi = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), sei = rep(0.2, 6),
    g = factor(rep(c("a", "b", "c"), 2)))
  factor_prior <- function(contrast) BayesTools::prior_factor(
    if (contrast %in% c("meandif", "orthonormal")) "mnormal" else "normal",
    list(0, 1), contrast = contrast)
  spike <- BayesTools::prior("spike", list(0))
  build <- function(..., constructor = BMA.norm) constructor(
    yi = yi, sei = sei, data = dat, measure = "SMD", only_priors = TRUE, silent = TRUE, ...)

  # An explicit global null takes the default contrast ("meandif" here), which
  # BayesTools refused for the other contrasts without naming the term.
  for (contrast in c("treatment", "independent", "orthonormal")) {
    for (constructor in list(BMA.norm, RoBMA)) {
      condition <- tryCatch(build(mods = ~ g, prior_mods = list(g = factor_prior(contrast)),
        prior_mods_null = spike, constructor = constructor), error = identity)
      expect_s3_class(condition, "RoBMA_prior_contrast_conflict")
      expect_s3_class(condition, "RoBMA_prior_error")
      message <- conditionMessage(condition)
      expect_match(message, paste0("The factor term 'g' uses the '", contrast,
        "' contrast in 'prior_mods' but the 'meandif' contrast in 'prior_mods_null'"), fixed = TRUE)
      expect_match(message, paste0("prior_mods_null = list(g = prior_factor(\"spike\", list(0), contrast = \"",
        contrast, "\"))"), fixed = TRUE)
      expect_match(message, "'set_contrast_factor_predictors'", fixed = TRUE)
      expect_match(message, "leave 'prior_mods_null' unspecified", fixed = TRUE)
    }
  }
  # The same conflict of a scale term names the scale arguments.
  condition <- tryCatch(build(scale = ~ g, prior_scale = list(g = factor_prior("independent")),
    prior_scale_null = spike, prior_heterogeneity_null = NULL), error = identity)
  expect_s3_class(condition, "RoBMA_prior_contrast_conflict")
  expect_match(conditionMessage(condition), "the 'independent' contrast in 'prior_scale' but the 'meandif' contrast in 'prior_scale_null'", fixed = TRUE)
  expect_match(conditionMessage(condition), "leave 'prior_scale_null' unspecified", fixed = TRUE)

  # A null given for the term itself is not built from the default contrast, so
  # setting the default contrast is no remedy there; leaving the term out is.
  for (case in list(
    list(alt = "independent", null = "meandif", prior_mods = list(g = factor_prior("independent"))),
    list(alt = "meandif", null = "independent", prior_mods = NULL)
  )) {
    args <- list(mods = ~ g, prior_mods_null = list(g = BayesTools::prior_factor(
      "spike", list(0), contrast = case$null)))
    args$prior_mods <- case$prior_mods
    condition <- tryCatch(do.call(build, args), error = identity)
    expect_s3_class(condition, "RoBMA_prior_contrast_conflict")
    message <- conditionMessage(condition)
    expect_match(message, paste0("uses the '", case$alt, "' contrast in 'prior_mods' but the '",
      case$null, "' contrast in 'prior_mods_null'"), fixed = TRUE)
    expect_match(message, paste0("contrast = \"", case$alt, "\"))"), fixed = TRUE)
    expect_match(message, "leave the term out of 'prior_mods_null'", fixed = TRUE)
    expect_false(grepl("set_contrast_factor_predictors", message, fixed = TRUE))
  }

  # A null with the contrast of the alternative sets the conflict aside, as does
  # setting the default contrast or leaving the null unspecified, and a
  # default-contrast alternative needs none of them.
  independent <- list(g = factor_prior("independent"))
  matching_null <- list(g = BayesTools::prior_factor("spike", list(0), contrast = "independent"))
  for (object in list(
    build(mods = ~ g, prior_mods = independent, prior_mods_null = matching_null),
    build(mods = ~ g, prior_mods = independent, prior_mods_null = spike,
          set_contrast_factor_predictors = "independent"),
    build(mods = ~ g, prior_mods = independent),
    build(mods = ~ g, prior_mods = list(g = factor_prior("meandif")), prior_mods_null = spike)
  )) {
    expect_s3_class(object, "BMA.norm")
  }
})

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
  expect_equal(attr(object_data[["data"]], "set_contrast_factor_predictors"), "meandif")
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


test_that("default independent formula warning is emitted once per ensemble formula", {

  dat <- data.frame(yi = c(.1, .2, .3, .4), sei = rep(.2, 4),
                    group = factor(rep(c("a", "b"), 2)))
  messages <- character()
  withCallingHandlers(BMA(yi = yi, sei = sei, mods = ~ group, data = dat,
    set_contrast_factor_predictors = "independent", measure = "SMD", only_priors = TRUE),
    warning = function(warning) {
      messages <<- c(messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    })
  expect_equal(sum(grepl("The 'mods' formula is overspecified", messages, fixed = TRUE)), 1L)
  messages <- character()
  withCallingHandlers(BMA(yi = yi, sei = sei, scale = ~ 0 + group, data = dat,
    set_contrast_factor_predictors = "independent", measure = "SMD", only_priors = TRUE),
    warning = function(warning) {
      messages <<- c(messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    })
  scale_message <- messages[grepl("The 'scale' formula is overspecified", messages, fixed = TRUE)]
  expect_length(scale_message, 1L)
  expect_match(scale_message, "'prior_scale'", fixed = TRUE)
  expect_false(any(grepl("scale = ~ 0", scale_message, fixed = TRUE)))
})
