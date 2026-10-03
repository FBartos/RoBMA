.ordinate_control_test_object <- function(n_rows) {

  object <- brma(yi = c(.1, .2), sei = c(.2, .3), measure = "GEN",
    prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("point", list(.2)),
    only_priors = TRUE, silent = TRUE)
  args <- .create_jags_formula_args(object$data, object$priors)
  compiled <- lapply(names(args$formula_list), function(parameter) {
    BayesTools::JAGS_formula(formula = args$formula_list[[parameter]], parameter = parameter,
      data = args$formula_data_list[[parameter]], prior_list = args$formula_prior_list[[parameter]],
      formula_scale = args$formula_scale_list[[parameter]], prior_random = args$formula_random_prior_list[[parameter]],
      random_effects_compile = args$formula_random_effects_compile_list[[parameter]])
  })
  names(compiled) <- names(args$formula_list)
  object$formula_design <- lapply(compiled, `[[`, "formula_design")
  samples <- cbind(mu = seq(-.1, .1, length.out = n_rows), tau = rep(.2, n_rows))
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = n_rows),
    class = c("runjags", "list"))
  attr(fit, "formula_design") <- object$formula_design
  priors <- list()
  for (formula in compiled) priors <- c(priors, formula$prior_list)
  attr(fit, "prior_list") <- c(priors, .create_fit_priors(object$data, object$priors))
  object$fit <- as_bayestools_fit(fit)
  object
}

test_that("unavailable ordinates retain the actual normalized request", {

  object <- .ordinate_control_test_object(10L)
  error <- tryCatch(hypothesis(object, "mu = 0", density_method = "qCMDE",
    density_control = list(samples = 70, target_relative_mcse = .01)), error = identity)
  expect_s3_class(error, "RoBMA_density_ordinate_error")
  expect_identical(conditionMessage(error), paste0(
    "qCMDE posterior ordinate is unavailable for 'effect = 0': fewer than 20 continuous active samples. ",
    "Inspect 'density_diagnostics()' for the computation stage."))
  actual <- density_diagnostics(error)
  expect_equal(actual$requested_samples, 70)
  expect_equal(actual$target_relative_mcse, .01)
  expect_identical(actual$status, "unavailable")
  expect_identical(actual$numerical_status, "not_computed")
  expect_true(is.na(actual$ordinate))
  expect_false(actual$bf_grade_met)
  raw <- .iwmde_estimate(.iwmde_context(object), "mu", "qCMDE",
    list(samples = 70, target_relative_mcse = .01, display_grid = "ordinate"),
    outputs = "ordinate", values = 0)
  expect_equal(raw$plan$control$samples, actual$requested_samples)
  expect_equal(raw$plan$control$target_relative_mcse, actual$target_relative_mcse)

  diagnostic <- raw$diagnostics$ordinate
  diagnostic$diagnostics$requested_samples <- 99
  diagnostic$diagnostics$target_relative_mcse <- .2
  plan <- raw$plan
  authoritative <- .iwmde_estimate_ordinate_failures(plan, diagnostic)
  expect_equal(authoritative$requested_samples, 70)
  expect_equal(authoritative$target_relative_mcse, .01)
  plan$control$samples <- NULL
  plan$control$target_relative_mcse <- NULL
  historical <- .iwmde_estimate_ordinate_failures(plan, diagnostic)
  expect_equal(historical$requested_samples, 99)
  expect_equal(historical$target_relative_mcse, .2)
  diagnostic$diagnostics$requested_samples <- NULL
  diagnostic$diagnostics$target_relative_mcse <- NULL
  missing <- .iwmde_estimate_ordinate_failures(plan, diagnostic)
  expect_true(is.na(missing$requested_samples))
  expect_true(is.na(missing$target_relative_mcse))
})

test_that("successful ordinate requests retain their sampling controls", {

  result <- hypothesis(.ordinate_control_test_object(80L), "mu = 0",
    density_method = "qCMDE",
    density_control = list(samples = 70, target_relative_mcse = .01))
  actual <- density_diagnostics(result)
  expect_equal(actual$requested_samples, 70)
  expect_equal(actual$target_relative_mcse, .01)
  expect_identical(actual$status, "ok")
  expect_true(actual$bf_grade_met)
})
