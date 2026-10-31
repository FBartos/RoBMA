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
  # the refusal of the draws carries these diagnostics, which the engine
  # transports (BayesTools test-hypothesis-BF-model.R)
  condition <- tryCatch(BayesTools::posterior_atoms_free(draws), error = identity)
  expect_identical(condition$diagnostics, diagnostics)
  malformed <- unavailable
  malformed$diagnostics[[1L]]$model_indices <- c(1L, 1L)
  expect_error(BayesTools::posterior_metadata(draws, "measure_unavailable") <- malformed, "invalid")
})

test_that("linear coordinate readers distinguish current model and primitive contexts", {
  # The reader of BayesTools words the refit request with the build and the
  # condition class RoBMA passes (as the RoBMA entry points do).
  linear_weight_space <- function(draws) {
    BayesTools::posterior_linear_weight_space(
      draws, build = "RoBMA/BayesTools", refit_class = "RoBMA_refit_required"
    )
  }
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
  expect_identical(linear_weight_space(draws(primitive)), "coefficient")
  expect_identical(linear_weight_space(draws(model)), "coefficient")
  expect_identical(linear_weight_space(draws(conditional)), "coefficient")
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
  expect_error(linear_weight_space(draws(old)), class = "RoBMA_refit_required")
  for(context in list(primitive, model, conditional)){
    for(version in c(2L, 3L)){
      context$transforms$mu$descriptor$schema_version <- version
      expect_error(linear_weight_space(draws(context)), class = "RoBMA_refit_required")
    }
  }
  expect_error(linear_weight_space(draws(model, "formula_contribution")),
    class = "RoBMA_refit_required")
})

test_that("current coefficient transforms support identity-scale public hypotheses", {
  set.seed(1L)
  object <- review_factor_replay(c("a", "b", "c"), "treatment",
    distribution = "normal", independent = TRUE, n = 24L)
  descriptor <- BayesTools::JAGS_formula_coefficient_transform(object$fit, "mu", "original")
  expect_identical(descriptor$schema_version, 4L)
  expect_identical(unname(descriptor$matrix), diag(3L))
  original <- RoBMA::hypothesis(object, "intercept = 0", component = "mods",
    standardized_coefficients = FALSE, density_method = "KDE", seed = 1L)
  standardized <- RoBMA::hypothesis(object, "intercept = 0", component = "mods",
    standardized_coefficients = TRUE, density_method = "KDE", seed = 1L)
  expect_s3_class(original, "BayesTools_hypothesis_BF")
  expect_identical(nrow(original), 1L)
  expect_true(is.finite(original$BF) && original$BF > 0)
  expect_identical(original$Null, "intercept = 0")
  expect_equal(original$BF, standardized$BF, tolerance = 1e-12)
})