context("Response formulas preserve structural zero intercepts")

test_that("response-only formulas preserve omitted intercepts in model constructors", {

  dat <- data.frame(yi = c(-0.2, 0.1, 0.4, 0.5), sei = rep(0.2, 4))
  constructors <- list(brma.norm, brma.mv, BMA.norm, BMA.mv)

  for (constructor in constructors) {
    for (formula in list(yi ~ 0, yi ~ -1)) {
      arguments <- list(
        sei = dat$sei, data = dat, measure = "GEN",
        prior_unit_information_sd = 1, only_priors = TRUE
      )
      from_response <- do.call(constructor, c(list(yi = formula), arguments))
      from_mods <- do.call(constructor, c(
        list(yi = dat$yi, mods = formula[-2L]), arguments
      ))

      expect_equal(from_response[["priors"]], from_mods[["priors"]])
      expect_null(from_response[["priors"]][["outcome"]][["mu"]])
      expect_equal(
        attr(stats::terms(attr(from_response[["data"]][["mods"]], "formula")),
             "intercept"),
        0L
      )
      expect_equal(
        .create_model_syntax(from_response[["data"]], from_response[["priors"]]),
        .create_model_syntax(from_mods[["data"]], from_mods[["priors"]])
      )
    }
  }
})

test_that("response-only zero intercepts survive random-formula compilation and row selection", {

  dat <- data.frame(
    yi = c(-0.2, 0.1, NA_real_, 0.5, 0.3),
    study = factor(c("a", "a", "b", "b", "c"))
  )
  expect_warning(
    object <- brma.mv(
      yi = yi ~ 0, vi = rep(0.04, 5), random = ~ 1 | study,
      data = dat, subset = 1:4, measure = "GEN",
      prior_unit_information_sd = 1, only_priors = TRUE
    ),
    "1 observation(s) removed due to missing values.", fixed = TRUE
  )

  expect_equal(nrow(object[["data"]][["mods"]]), 3L)
  expect_true(BayesTools::is.prior.point(object[["priors"]][["location"]][["intercept"]]))
  expect_equal(mean(object[["priors"]][["location"]][["intercept"]]), 0)
  design <- .fitted_formula_design(object, "mu", required = TRUE)
  expect_equal(mean(design[["prior_list"]][["mu_intercept"]]), 0)
})

test_that("an explicit intercept prior cannot override a structural zero intercept", {

  dat <- data.frame(yi = c(.1, .3, .2), x = c(1, 2, 4))
  slope <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  nonzero <- BayesTools::prior("normal", list(mean = 1, sd = 1))
  construct <- function(constructor, formula, intercept) {
    prior_mods <- list(x = slope)
    if (!is.null(intercept)) prior_mods <- c(list(intercept = intercept), prior_mods)
    constructor(
      yi = yi, sei = rep(.2, 3), data = dat, mods = formula, measure = "SMD",
      prior_mods = prior_mods, only_priors = TRUE
    )
  }

  for (constructor in list(brma, RoBMA)) {
    # the structural zero of '~ 0 + x' is not replaceable by a nonzero prior
    expect_error(
      construct(constructor, ~ 0 + x, nonzero),
      class = "RoBMA_prior_intercept_conflict"
    )
    # the family class of the prior-input errors
    expect_error(
      construct(constructor, ~ 0 + x, nonzero),
      class = "RoBMA_prior_error"
    )
    # a point prior at zero repeats the structural zero
    control <- construct(constructor, ~ 0 + x, NULL)
    for (intercept in list(BayesTools::prior("spike", list(0)),
                           BayesTools::prior("point", list(0)))) {
      accepted <- construct(constructor, ~ 0 + x, intercept)
      expect_identical(accepted[["priors"]][["mods"]][["intercept"]],
                       control[["priors"]][["mods"]][["intercept"]])
      expect_equal(mean(accepted[["priors"]][["mods"]][["intercept"]]), 0)
    }
    # a formula with an intercept keeps accepting an explicit intercept prior
    with_intercept <- construct(constructor, ~ x, nonzero)
    expect_equal(mean(with_intercept[["priors"]][["mods"]][["intercept"]]), 1)
  }
  expect_error(
    construct(brma, ~ 0 + x, BayesTools::prior("point", list(0.5))),
    class = "RoBMA_prior_intercept_conflict"
  )
})
