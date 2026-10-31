context("Formula declaration ownership through selected draws")

test_that("selected formula draw carriers preserve owned declarations", {
  slope <- BayesTools::prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "sigma"
  draws <- cbind(mu_x = c(.5, 1.5))
  attr(draws, "prior_list") <- list(sigma = BayesTools::prior("point", list(3)))
  carrier <- BayesTools::JAGS_formula_draws(draws, ~ x, "mu", data.frame(x = c(-2, 1)),
    list(intercept = BayesTools::prior("point", list(2)), x = slope))
  selected <- .posterior_formula_fit(carrier, as.matrix(carrier)[2, , drop = FALSE])
  expect_identical(attr(selected, "prior_list"), attr(carrier, "prior_list"))
  expect_equal(unname(BayesTools::JAGS_evaluate_formula(selected, parameter = "mu")), matrix(2 + 4.5 * c(-2, 1), 2, 1))
  expect_null(.point_prior_value(BayesTools::prior("point", list(location = expression(theta)))))
  expect_identical(.point_prior_value(BayesTools::prior("point", list(3))), 3)
})

test_that("actual prior-only model formula owners retain their full declarations", {
  data <- data.frame(x = c(10, 20, 30), yi = c(.1, .2, .3), sei = rep(.5, 3))
  object <- brma(yi = yi, sei = sei, mods = ~ x, data = data, measure = "GEN",
    prior_effect = BayesTools::prior("point", list(2)),
    prior_mods = list(x = BayesTools::prior("normal", list(0, 1))),
    prior_heterogeneity = BayesTools::prior("point", list(.4)), only_priors = TRUE)
  formula_fit <- .object_formula_fit(object)
  expect_true(all(c("mu_intercept", "mu_x") %in% names(attr(formula_fit, "prior_list"))))
  posterior <- .posterior_formula_fit(formula_fit, cbind(mu_x = c(.5, 1.5)))
  expect_identical(attr(posterior, "prior_list"), attr(formula_fit, "prior_list"))
  expect_true(all(is.finite(BayesTools::JAGS_evaluate_formula(posterior, parameter = "mu"))))
})
