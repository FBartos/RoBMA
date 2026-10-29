context("BayesTools formula coefficient density contract")

test_that("conditional-normal ordinates retain exact structural classification", {

  priors <- list(
    a = BayesTools::prior("normal", list(0, 1)),
    b = BayesTools::prior("normal", list(0, 1)),
    s = BayesTools::prior("cauchy", list(0, 1), list(0, 5))
  )
  attr(priors$b, "multiply_by") <- "s"
  density <- BayesTools:::.prior_linear_combination_density(
    priors, c(a = 1, b = -1), n_grid = 512
  )
  ordinate <- .iwmde_prior_ordinate_classifications(density, 0)[[1L]]
  expect_identical(ordinate$method, "conditional_normal_mixture")
  expect_true(ordinate$exact)
  expect_true(ordinate$eligible)
  expect_false(ordinate$provenance$integration$exact)
  expect_identical(.iwmde_validate_prior_ordinate(ordinate, 0), ordinate)
  expect_length(.iwmde_ordinate_prior_warnings("mu_intercept", list(ordinate)), 0L)
})


test_that("IWMDE classifies prior ordinates of every route by eligibility", {

  density <- BayesTools:::.prior_linear_combination_density(
    list(
      a = BayesTools::prior("cauchy", list(0, 1)),
      b = BayesTools::prior("t", list(0, 1, 3))
    ),
    c(a = 1, b = 1),
    n_grid = 256
  )
  # A route outside a fixed method vocabulary is classified by the BayesTools
  # exactness rule.
  ordinate <- .iwmde_prior_ordinate_classifications(density, 0.3)[[1L]]
  expect_identical(ordinate$method, "convolution")
  expect_true(ordinate$eligible)
  expect_true(is.na(ordinate$condition))
  expect_identical(.iwmde_validate_prior_ordinate(ordinate, 0.3), ordinate)
  expect_length(.iwmde_ordinate_prior_warnings("mu_intercept", list(ordinate)), 0L)

  # A regular ordinate that the rule refuses warns before estimation.
  testthat::local_mocked_bindings(
    prior_ordinate_status = function(prior_density, values, labels = NULL) {
      out <- eligible_ordinate_status(values)
      out[["eligible"]]  <- FALSE
      out[["condition"]] <- "BayesTools_inexact_ordinate"
      out[["reason"]]    <- "Inexact."
      out
    },
    .package = "BayesTools"
  )
  inexact <- .iwmde_prior_ordinate_classifications(density, 0.3)[[1L]]
  expect_identical(inexact$behavior, "regular")
  expect_false(inexact$eligible)
  expect_identical(inexact$condition, "BayesTools_inexact_ordinate")
  expect_match(
    .iwmde_ordinate_prior_warnings("mu_intercept", list(inexact)),
    "could not be classified from deterministic provenance",
    fixed = TRUE
  )
})


test_that("transformed coefficient plots use exact structural weights", {

  # Synthetic compiler/consumer rows, not posterior-inference evidence.
  compiled <- BayesTools::JAGS_formula(
    formula = ~ 1 + x, parameter = "mu", data = data.frame(x = c(2, 4, 6)),
    prior_list = list(
      intercept = BayesTools::prior("normal", list(0, 1)),
      x = BayesTools::prior("normal", list(0, 1))
    ), formula_scale = TRUE
  )
  fit <- coda::mcmc.list(coda::mcmc(cbind(
    mu_intercept = c(-.2, 0, .2), mu_x = c(-.1, 0, .1)
  )))
  attr(fit, "prior_list") <- compiled[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = compiled[["formula_design"]])
  attr(fit, "formula_scale") <- list(mu = compiled[["formula_scale"]])
  fit <- as_bayestools_fit(fit)

  entry <- list(
    component         = "mods",
    role              = "fixed_coefficient",
    formula_parameter = "mu"
  )

  observed <- .plot_brma_formula_parameter_spec(
    object                    = list(fit = fit),
    parameter                 = "mu_x",
    parameter_entry           = entry,
    standardized_coefficients = FALSE
  )
  standardized <- .plot_brma_formula_parameter_spec(
    object                    = list(fit = fit),
    parameter                 = "mu_x",
    parameter_entry           = entry,
    standardized_coefficients = TRUE
  )

  expect_identical(observed[["type"]], "linear")
  expect_identical(observed[["weights"]], c(mu_x = 0.5))
  expect_null(standardized)
})


test_that("ordinary parameters retain legacy density-scale alignment", {

  entry <- list(
    component         = "mods",
    role              = "fixed_coefficient",
    formula_parameter = NA_character_
  )

  observed <- .plot_brma_formula_parameter_spec(
    object                    = list(),
    parameter                 = "mu",
    parameter_entry           = entry,
    standardized_coefficients = FALSE
  )

  expect_null(observed)
})


