source(testthat::test_path("common-functions.R"))

.pr82_random_mixture_oracle <- function(object, input, samples) {

  design <- BayesTools::JAGS_formula_design(object[["fit"]], "mu")
  allocation <- design[["random_allocations"]][[1L]]
  source <- allocation[["source"]][["name"]]
  weights <- BayesTools::JAGS_indexed_parameter_matrix(samples, allocation[["weight_name"]])
  study <- input[["data"]][["study"]]
  study_covariance <- outer(study, study, `==`) * 1
  K <- nrow(input[["V"]])
  X <- design[["model_matrix"]]
  covariance <- function(row, value = samples[row, source]) {
    study_gate <- samples[row, allocation[["inclusion"]][["study"]][["indicator_name"]]]
    observation_gate <- samples[row, allocation[["inclusion"]][["observation"]][["indicator_name"]]]
    input[["V"]] + value^2 * (weights[row, 1L] * study_gate * study_covariance +
      weights[row, 2L] * observation_gate * diag(K))
  }
  log_likelihood <- function(row, value = samples[row, source]) {
    M <- covariance(row, value)
    mean <- as.vector(X %*% samples[row, c("mu_intercept", "mu_x")])
    factor <- chol(M)
    z <- forwardsolve(t(factor), input[["data"]][["yi"]] - mean)
    -.5 * K * log(2 * pi) - sum(log(diag(factor))) - sum(z^2) / 2
  }
  list(source = source, covariance = covariance, log_likelihood = log_likelihood,
       design = X, allocation = allocation)
}

test_that("random-allocation prior branches preserve Gaussian likelihood and prior ratios", {

  skip_if_missing_fits("BMA.mv_random_components")
  object <- load_fit("BMA.mv_random_components")
  input <- load_info("BMA.mv_random_components")
  context <- .iwmde_context(object)
  samples <- context[["posterior_samples"]]
  oracle <- .pr82_random_mixture_oracle(object, input, samples)
  expect_identical(oracle[["allocation"]][["scale"]], "total_variance")
  source <- oracle[["source"]]
  indicator <- .iwmde_indicator_name(source)
  indices <- as.integer(samples[, indicator])
  rows <- as.integer(vapply(sort(unique(indices)), function(index) which(indices == index)[[1L]], integer(1L)))
  expect_setequal(indices[rows], c(1L, 2L))
  states <- .iwmde_row_states(context, rows, source, list(type = "primitive"))
  replacement <- list(type = "scalar", name = source)
  source_prior <- context[["flat_prior_list"]][[source]]
  for (i in seq_along(rows)) {
    row <- rows[[i]]
    state <- states[[i]]
    expect_equal(state[["baseline_log_lik"]], oracle[["log_likelihood"]](row), tolerance = 1e-12)
    localized <- .iwmde_likelihood_posterior_samples(context, samples[row, , drop = FALSE], state[["active_setup"]])
    expect_identical(localized[, indicator], samples[row, indicator])
    gates <- vapply(oracle[["allocation"]][["inclusion"]], `[[`, character(1L), "indicator_name")
    expect_identical(localized[, gates, drop = FALSE], samples[row, gates, drop = FALSE])
    prior_sd <- source_prior[[indices[[row]]]][["parameters"]][["sd"]]
    for (value in c(.2, .4)) {
      actual <- .iwmde_log_q_replacement(context, state, source, value, replacement) - state[["baseline_log_q"]]
      expected <- oracle[["log_likelihood"]](row, value) - oracle[["log_likelihood"]](row) +
        stats::dnorm(value, 0, prior_sd, log = TRUE) -
        stats::dnorm(samples[row, source], 0, prior_sd, log = TRUE)
      expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)
    }
  }
})

test_that("public qCMDE random-mixture hypotheses match conditional Normal ordinates", {

  skip_if_missing_fits("BMA.mv_random_components")
  object <- load_fit("BMA.mv_random_components")
  input <- load_info("BMA.mv_random_components")
  context <- .iwmde_context(object)
  samples <- context[["posterior_samples"]]
  oracle <- .pr82_random_mixture_oracle(object, input, samples)
  active <- which(.iwmde_parameter_active_rows(context, "mu_intercept"))
  rows <- .nested_srs_rows(active, 20L)
  prior <- .iwmde_focal_prior(context, "mu_intercept", samples[rows[[1L]], ])
  prior_mean <- prior[["parameters"]][["mean"]]
  prior_sd <- prior[["parameters"]][["sd"]]
  x <- oracle[["design"]][, 1L]
  expected <- vapply(rows, function(row) {
    precision_x <- solve(oracle[["covariance"]](row), x)
    variance <- 1 / (1 / prior_sd^2 + sum(x * precision_x))
    residual <- input[["data"]][["yi"]] - oracle[["design"]][, 2L] * samples[row, "mu_x"]
    mean <- variance * (prior_mean / prior_sd^2 + sum(precision_x * residual))
    stats::dnorm(.1, mean, sqrt(variance))
  }, numeric(1L))
  result <- hypothesis(object, "mu = 0.1", conditional = TRUE,
    standardized_coefficients = TRUE, density_method = "qCMDE",
    density_control = list(samples = 20L, n_points = 20L,
      normalization_points = 300L, normalization_prob = 1 - 1e-8), columns = "all")
  diagnostic <- density_diagnostics(result)
  expect_identical(diagnostic[["status"]], "ok")
  expect_true(diagnostic[["bf_grade_met"]])
  expect_identical(diagnostic[["achieved_row_budget"]], 20L)
  # A fixed row sample leaves only deterministic grid-normalization error in
  # this comparison with exact conditional Normal densities.
  expect_lt(abs(result[["posterior"]] / mean(expected) - 1), 1e-6)
  # The public prior height uses BayesTools' shared linear-density interpolation
  # with a 1e-4 relative refinement target, unlike the analytic reference here.
  expect_lt(abs(as.numeric(result[["prior"]]) / stats::dnorm(.1, prior_mean, prior_sd) - 1), 1e-4)
  expect_lt(abs(attr(result, "raw_BF") / (stats::dnorm(.1, prior_mean, prior_sd) / mean(expected)) - 1), 1e-4)
})
