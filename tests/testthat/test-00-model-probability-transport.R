test_that("replacement draws retain declared measure refusals", {
  draws <- structure(seq_len(5) / 5, class = c("mixed_posteriors.simple", "numeric"))
  unavailable <- data.frame(column = "theta", measure = "support",
    reason = "The declared support is unavailable.", cause = "unsupported_contribution_measure")
  BayesTools::posterior_metadata(draws, "measure_unavailable") <- unavailable
  replaced <- .brma_draws_with_values(draws, as.numeric(draws))
  expect_identical(BayesTools::posterior_metadata(replaced, "measure_unavailable"), unavailable)
})

test_that("model probability diagnostics survive RoBMA refusal transport", {
  draws <- structure(seq_len(5) / 5, class = c("mixed_posteriors.simple", "numeric"))
  attr(draws, "parameter") <- "theta"
  diagnostics <- list(model_indices = 1:2, log_prior_probabilities = log(c(.5, .5)),
    log_posterior_probabilities = c(0, -1000), stage = "posterior")
  unavailable <- data.frame(column = "theta", measure = "atoms", reason = "Positive continuous mass is unavailable.",
    cause = "numerical_model_probability_unavailable")
  unavailable$diagnostics <- list(diagnostics)
  BayesTools::posterior_metadata(draws, "measure_unavailable") <- unavailable
  replaced <- .brma_draws_with_values(draws, as.numeric(draws))
  expect_identical(BayesTools::posterior_metadata(replaced, "measure_unavailable")$diagnostics[[1L]], diagnostics)
  condition <- tryCatch(BayesTools::posterior_atoms_free(draws), error = identity)
  refusal <- .hypothesis_refusal_from_condition(condition)
  expect_identical(refusal$fields$diagnostics, diagnostics)
  expect_identical(refusal$reason, condition$message)
  returned <- tryCatch(.hypothesis_stop(refusal), error = identity)
  expect_identical(returned$diagnostics, diagnostics)
  malformed <- unavailable
  malformed$diagnostics[[1L]]$model_indices <- c(1L, 1L)
  expect_error(BayesTools::posterior_metadata(draws, "measure_unavailable") <- malformed, "invalid")
})

test_that("linear coordinate readers distinguish current model and primitive contexts", {
  compiled <- BayesTools::JAGS_formula(~ x, "mu", data.frame(x = c(10, 20, 30)),
    list(intercept = BayesTools::prior("normal", list(0, 1)),
      x = BayesTools::prior("normal", list(0, 1))), formula_scale = TRUE)
  columns <- c("mu_intercept", "mu_x")
  scales <- list(mu = compiled$formula_scale)
  primitive <- BayesTools:::.prior_density_context(compiled$prior_list, columns, scales)
  model <- BayesTools:::.prior_density_model_mixture_context(
    lapply(compiled$prior_list, function(prior) list(prior, prior)), columns)
  model$transforms <- primitive$transforms
  conditional <- BayesTools:::.prior_density_conditional_context(
    compiled$prior_list, columns, "mu_x", formula_scale = scales)
  draws <- function(context, space = "coefficient") with_draw_metadata(c(1, 2),
    linear_weights = c(mu_intercept = 1, mu_x = 1), linear_weight_space = space,
    prior_context = context)
  expect_identical(.brma_linear_weight_space(draws(primitive)), "coefficient")
  expect_identical(.brma_linear_weight_space(draws(model)), "coefficient")
  expect_identical(.brma_linear_weight_space(draws(conditional)), "coefficient")
  old <- model
  old$schema_version <- 1L
  expect_error(draws(old), "Recompute or refit")
  corrupt <- model
  corrupt$model_log_weights[[1L]] <- -10
  expect_error(draws(corrupt), "ownership")
  corrupt <- model
  corrupt$model_probability_declaration <- NULL
  expect_error(draws(corrupt), "ownership")
  old <- primitive
  old$schema_version <- 2L
  expect_error(.brma_linear_weight_space(draws(old)), class = "RoBMA_refit_required")
  for(context in list(primitive, model, conditional)){
    context$transforms$mu$descriptor$schema_version <- 2L
    expect_error(.brma_linear_weight_space(draws(context)), class = "RoBMA_refit_required")
  }
  expect_error(.brma_linear_weight_space(draws(model, "formula_contribution")),
    class = "RoBMA_refit_required")
})
