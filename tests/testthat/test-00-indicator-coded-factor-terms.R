context("Indicator-coded factor terms under mean-difference contrasts")

# An interaction without one of its lower-order terms (`g:x` in `~ g + g:x` or
# `~ g / x`) codes `g` by level indicators and has one coefficient per level.
# BayesTools defines mean-difference and orthonormal priors on contrast
# coefficients only, so the model-averaging constructors, whose default factor
# contrast is "meandif", stop for such terms. The documented remedies are an
# independent prior for that term, the full interaction `~ g * x`, or
# treatment contrasts. Only the prior construction is checked (no fitting).

.indicator_coded_data <- function() {

  data.frame(
    yi  = c(-0.2, 0.1, 0.4, 0.5, 0.3, 0.0, 0.2, 0.6, -0.1),
    sei = rep(0.2, 9),
    g   = factor(rep(c("a", "b", "c"), 3)),
    x   = c(0.1, 0.5, -0.3, 1.2, 0.4, -0.8, 0.9, 0.0, -0.5)
  )
}

.indicator_coded_arguments <- function(dat) {

  list(
    yi = dat[["yi"]], sei = dat[["sei"]], data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE
  )
}

.indicator_coded_message <- function(contrast) {

  paste0(
    "The '", contrast, "' prior of the factor term 'g:x' is unavailable: the ",
    "formula has no term 'x', so 'g:x' codes 'g' by level indicators and has ",
    "one coefficient per level instead of '", contrast, "' contrast coefficients."
  )
}

.all_components <- function(prior, predicate) {

  all(c(predicate(prior), vapply(prior, predicate, logical(1))))
}


test_that("model-averaging constructors stop for indicator-coded factor terms by default", {

  dat    <- .indicator_coded_data()
  common <- .indicator_coded_arguments(dat)

  for (constructor in c("BMA.norm", "RoBMA", "BMA.mv", "RoBMA.mv")) {
    for (mods in list(~ g + g:x, ~ g / x)) {
      expect_error(
        do.call(constructor, c(common, list(mods = mods))),
        .indicator_coded_message("meandif"),
        fixed = TRUE
      )
    }
  }

  expect_error(
    do.call("BMA.norm", c(common, list(
      mods = ~ g + g:x, set_contrast_factor_predictors = "orthonormal"
    ))),
    .indicator_coded_message("orthonormal"),
    fixed = TRUE
  )
  # scale formulas use the same default contrast
  expect_error(
    suppressWarnings(do.call("BMA.norm", c(common, list(
      mods = ~ x, scale = ~ g + g:x
    )))),
    .indicator_coded_message("meandif"),
    fixed = TRUE
  )

  counts <- data.frame(
    ai = c(4, 6, 3, 62, 33, 180, 8, 505, 29),
    bi = c(119, 300, 228, 13536, 5036, 1361, 2537, 87886, 7470),
    ci = c(11, 29, 11, 248, 47, 372, 10, 499, 45),
    di = c(128, 274, 209, 12619, 5761, 1079, 619, 87892, 7232),
    g  = dat[["g"]],
    x  = dat[["x"]]
  )
  expect_error(
    BMA.glmm(
      ai = ai, bi = bi, ci = ci, di = di, data = counts, measure = "OR",
      mods = ~ g + g:x, only_priors = TRUE
    ),
    .indicator_coded_message("meandif"),
    fixed = TRUE
  )

  # single-model constructors default to treatment contrasts and are unaffected
  expect_no_error(do.call("brma", c(common, list(mods = ~ g + g:x))))
})


test_that("documented remedies for indicator-coded factor terms construct their priors", {

  dat    <- .indicator_coded_data()
  common <- .indicator_coded_arguments(dat)

  for (constructor in c("BMA.norm", "RoBMA")) {
    for (mods in list(~ g + g:x, ~ g / x)) {

      independent <- do.call(constructor, c(common, list(
        mods       = mods,
        prior_mods = list("g:x" = prior_factor(
          "normal", list(0, 0.5), contrast = "independent"
        ))
      )))
      priors <- independent[["priors"]][["mods"]]
      expect_true(.all_components(priors[["g:x"]], BayesTools::is.prior.independent))
      expect_true(.all_components(priors[["g"]], BayesTools::is.prior.meandif))

      treatment <- do.call(constructor, c(common, list(
        mods = mods, set_contrast_factor_predictors = "treatment"
      )))
      priors <- treatment[["priors"]][["mods"]]
      expect_true(.all_components(priors[["g:x"]], BayesTools::is.prior.treatment))
      expect_true(.all_components(priors[["g"]], BayesTools::is.prior.treatment))
    }

    full <- do.call(constructor, c(common, list(mods = ~ g * x)))
    priors <- full[["priors"]][["mods"]]
    expect_setequal(names(priors), c("intercept", "g", "x", "g:x"))
    expect_true(.all_components(priors[["g:x"]], BayesTools::is.prior.meandif))
  }
})
