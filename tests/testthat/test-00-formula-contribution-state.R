source(testthat::test_path("common-functions.R"))

context("Compiled coefficient and contribution state adapters")

test_that("actual brma coefficients and static-method routes use the compiled owner", {

  skip_on_cran()
  fit <- load_fit("lane_contribution_point")
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

test_that("public contribution combinations fold constants once for both ordinate methods", {

  skip_on_cran()
  data <- data.frame(g = factor(rep(c("A", "B"), 3)),
    yi = c(-.2, .3, 0, .4, .1, .2), sei = rep(.5, 6))
  for(multiplier in list(2, 1, "tau")){
    slope <- BayesTools::prior_factor("normal", list(0, 1), contrast = "treatment")
    attr(slope, "multiply_by") <- multiplier
    key <- if (is.numeric(multiplier)) paste0("lane_contribution_constant", multiplier) else "lane_contribution_tau"
    fit <- load_fit(key)
    means <- marginal_means(fit, bf = FALSE)
    plan <- .hypothesis_plan_marginal_means_build(means, "g[B] - g[A] = 0.1")
    spec <- plan$targets[[1L]]$spec
    if(is.numeric(multiplier)){
      expect_equal(as.numeric(spec$weights), multiplier, tolerance = 1e-14)
      for(method in c("qCMDE", "IWMDE")){
        result <- hypothesis(means, "g[B] - g[A] = 0.1", density_method = method, seed = 19, logBF = TRUE)
        expect_true(is.finite(as.numeric(result$BF)))
        expect_true(length(attr(result, "raw_log_BF")) == 1L)
      }
    }else{
      expect_identical(spec$type, "unsupported_formula_transform")
      for(method in c("qCMDE", "IWMDE")){
        expect_error(hypothesis(means, "g[B] - g[A] = 0.1", density_method = method),
          class = "RoBMA_density_method_original_scale")
      }
      expect_true(is.finite(as.numeric(hypothesis(means,
        "g[B] - g[A] = 0.1", density_method = "KDE")$BF)))
    }
  }
})

test_that("scalar discovery retains declared refusal rows and usable sibling quantities", {

  skip_on_cran()
  object <- load_fit("lane_contribution_discovery")
  transform <- BayesTools::JAGS_formula_coefficient_transform(object$fit, "mu")
  expect_identical(transform$prior_recipes$mu_x$reason, "state_dependent_map")
  values <- BayesTools::transform_scale_samples(object$fit)
  expect_true(all(is.finite(values[, "mu_x"])))
  quantities <- hypothesis_quantities(object)
  refused <- quantities[quantities$parameter == "mu_x", ]
  expect_true(nrow(refused) > 0L)
  expect_false(any(refused$point_test | refused$direction_test))
  expect_true(all(nzchar(refused$reason)))
  sibling <- quantities[quantities$parameter == "tau", ]
  expect_true(any(sibling$direction_test))
  expect_error(hypothesis(object, "mu_x = 0", density_method = "KDE"), class = "BayesTools_formula_measure_unavailable")
})

test_that("scalar and exp-affine plans catch only declared measure refusal conditions", {

  condition <- tryCatch(BayesTools:::.bt_formula_density_stop("Declared law unavailable.", reason = "numerical_scale_unavailable"), error = identity)
  plan <- list(parameter = "mu_x", selected = list(entry = list()), standardized = FALSE)
  for(exp in c(FALSE, TRUE)){
    testthat::with_mocked_bindings({
      result <- if(exp) .hypothesis_plan_exp_affine(plan, list(), list(), list()) else .hypothesis_plan_scalar(plan, list(), list())
      expect_identical(result$refusal, .hypothesis_refusal_from_condition(condition))
      expect_identical(result$kind, if(exp) "exp_affine" else "scalar")
    }, .hypothesis_plan_scalar_available = function(...) stop(condition),
       .hypothesis_plan_exp_affine_available = function(...) stop(condition))
    for(class in c("BayesTools_refit_required", "BayesTools_formula_transform_unavailable", "simpleError")){
      error <- structure(list(message = "Strict owner or programming failure.", call = NULL), class = c(class, "error", "condition"))
      testthat::with_mocked_bindings({
        expect_error(if(exp) .hypothesis_plan_exp_affine(plan, list(), list(), list()) else .hypothesis_plan_scalar(plan, list(), list()), class = class)
      }, .hypothesis_plan_scalar_available = function(...) stop(error),
         .hypothesis_plan_exp_affine_available = function(...) stop(error))
    }
  }
})
