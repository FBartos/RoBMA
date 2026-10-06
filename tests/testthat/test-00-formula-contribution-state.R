context("Compiled coefficient and contribution state adapters")

test_that("actual brma coefficients and static-method routes use the compiled owner", {

  skip_on_cran()
  slope <- BayesTools::prior("point", list(2))
  attr(slope, "multiply_by") <- "tau"
  data <- data.frame(x = c(10, 20, 30), yi = c(5, 5, 5), sei = rep(1, 3))
  fit <- expect_warning(brma(yi = yi, sei = sei, mods = ~ x, data = data, measure = "GEN",
    prior_effect = BayesTools::prior("point", list(5)), prior_mods = list(x = slope),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1), truncation = list(lower = 0)),
    chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 281, silent = TRUE),
    "Only one chain was run", fixed = TRUE)
  draws <- as.matrix(coda::as.mcmc(fit$fit))
  original <- .brma_as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  expect_equal(as.numeric(original$mu_intercept), 5 - 4 * draws[, "tau"], tolerance = 1e-12)
  expect_equal(as.numeric(original$mu_x), rep(.2, nrow(draws)), tolerance = 1e-14)
  selected <- list(parameter = "mu_intercept", component = "mods",
    entry = list(formula_parameter = "mu", component = "mods", role = "fixed_coefficient", status = "structural", fixed_value = 5))
  route <- .brma_formula_coefficient_route(fit, selected)
  expect_identical(route$type, "state_dependent")
  expect_null(route$cause)
  expect_true("tau" %in% names(route$state_sources))
  expect_error(.plot_brma_formula_parameter_spec(fit, "mu_intercept", selected$entry, FALSE),
    class = "RoBMA_density_method_original_scale")
  info <- .plot_prior_formula_info(fit, list(source = "mods", formula_parameter = "mu"))
  expect_true("tau" %in% names(info$prior_list))
  law <- BayesTools::JAGS_formula_prior_density(fit$fit, "mu", "mu_intercept")
  expect_equal(BayesTools::prior_density_ordinate(law, 1)$log_density,
    stats::dnorm(1, log = TRUE) + log(2) - log(4), tolerance = 1e-12)
})

test_that("formula contribution weights fold constants and refuse stochastic products", {

  compiled <- BayesTools::JAGS_formula(~ x, "mu", data.frame(x = c(10, 20, 30)),
    list(intercept = BayesTools::prior("normal", list(0, 1)), x = BayesTools::prior("normal", list(0, 1))),
    formula_scale = TRUE)
  for(multiplier in list(-3, "sigma")){
    priors <- compiled$prior_list
    attr(priors$mu_x, "multiply_by") <- multiplier
    priors$sigma <- BayesTools::prior("normal", list(0, 1))
    context <- BayesTools:::.prior_density_context(priors, c("mu_intercept", "mu_x"))
    sample <- with_draw_metadata(c(1, 2), linear_weights = c(mu_intercept = 1, mu_x = 1),
      linear_weight_space = "formula_contribution", prior_context = context)
    if(is.numeric(multiplier)){
      expect_identical(.brma_fitted_linear_weights(sample, c(mu_intercept = 1, mu_x = 1)),
        c(mu_intercept = 1, mu_x = -3))
    }else{
      expect_error(.brma_fitted_linear_weights(sample, c(mu_intercept = 1, mu_x = 1)),
        class = "RoBMA_density_method_original_scale")
    }
  }
  unavailable <- data.frame(column = "mu_intercept", measure = "prior_density",
    reason = "unsupported_contribution_measure")
  sample <- with_draw_metadata(c(1, 2), measure_unavailable = unavailable)
  expect_error(.brma_formula_measure_check(sample, "prior_density"),
    class = "BayesTools_formula_prior_density_unavailable")
  expect_silent(.brma_formula_measure_check(sample, "atoms"))
})
