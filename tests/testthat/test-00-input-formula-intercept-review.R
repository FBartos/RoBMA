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
