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

# The exact conditional Normal law of the intercept in every row: the Gaussian
# likelihood in the intercept times its Normal prior.
.pr82_random_mixture_conditionals <- function(oracle, input, samples, rows,
                                              prior_mean, prior_sd) {

  x <- oracle[["design"]][, 1L]
  out <- t(vapply(rows, function(row) {
    precision_x <- solve(oracle[["covariance"]](row), x)
    variance <- 1 / (1 / prior_sd^2 + sum(x * precision_x))
    residual <- input[["data"]][["yi"]] - oracle[["design"]][, 2L] * samples[row, "mu_x"]
    c(variance * (prior_mean / prior_sd^2 + sum(precision_x * residual)), sqrt(variance))
  }, numeric(2L)))
  colnames(out) <- c("mean", "sd")
  out
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

test_that("qCMDE normalization classifies random-mixture rows as exact Normal laws", {

  skip_if_missing_fits("BMA.mv_random_components")
  object <- load_fit("BMA.mv_random_components")
  input <- load_info("BMA.mv_random_components")
  context <- .iwmde_context(object)
  samples <- context[["posterior_samples"]]
  oracle <- .pr82_random_mixture_oracle(object, input, samples)
  plan <- .iwmde_plan(context, "mu_intercept", "qCMDE",
    list(samples = 20L, n_points = 20L), outputs = "ordinate", values = .1)
  states <- plan[["rows"]][["row_states"]]
  rows <- plan[["rows"]][["estimator_rows"]]
  log_q <- .iwmde_log_q_grid(context, "mu_intercept", .1, states,
    plan[["replacement"]])
  laws <- .iwmde_qcmde_row_laws(context, states, plan[["replacement"]],
    attr(log_q, "gaussian_kernel", exact = TRUE))
  prior <- .iwmde_focal_prior(context, "mu_intercept", samples[rows[[1L]], ])
  conditionals <- .pr82_random_mixture_conditionals(oracle, input, samples, rows,
    prior[["parameters"]][["mean"]], prior[["parameters"]][["sd"]])

  # The package's Gaussian kernel with the Normal prior reproduces each row's
  # conditional law computed from the full covariance.
  expect_identical(laws[["kind"]], rep("exact", length(rows)))
  expect_lt(max(abs(laws[["mean"]] - conditionals[, "mean"]) /
    conditionals[, "sd"]), 1e-8)
  expect_lt(max(abs(laws[["sd"]] / conditionals[, "sd"] - 1)), 1e-8)

  # The qCMDE range is the union of the rows' central intervals: the X7
  # reference for these rows at 1 - 1e-8 is [-0.56928, 0.77954].
  intervals <- .iwmde_qcmde_law_intervals(laws, 1 - 1e-8)
  union <- range(intervals[["intervals"]])
  expect_equal(union, c(-0.56928, 0.77954), tolerance = 1e-4)
  expect_equal(union, range(c(
    conditionals[, "mean"] - stats::qnorm(1 - 5e-9) * conditionals[, "sd"],
    conditionals[, "mean"] + stats::qnorm(1 - 5e-9) * conditionals[, "sd"]
  )), tolerance = 1e-8)
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
  conditionals <- .pr82_random_mixture_conditionals(oracle, input, samples, rows,
    prior_mean, prior_sd)
  expected <- stats::dnorm(.1, conditionals[, "mean"], conditionals[, "sd"])
  ordinate <- function(probability) {
    hypothesis(object, "mu = 0.1", conditional = TRUE,
      standardized_coefficients = TRUE, density_method = "qCMDE",
      density_control = list(samples = 20L, n_points = 20L,
        normalization_points = 300L, normalization_prob = probability),
      columns = "all")
  }
  result <- ordinate(1 - 1e-8)
  diagnostic <- density_diagnostics(result)
  expect_identical(diagnostic[["status"]], "ok")
  expect_true(diagnostic[["bf_grade_met"]])
  expect_identical(diagnostic[["achieved_row_budget"]], 20L)
  # Error model. A fixed row sample leaves only deterministic normalization
  # error in this comparison with exact conditional Normal densities: each
  # row's tail truncation t_j, exact for these Normal rows, overstates its
  # density by the factor 1 / (1 - t_j), at most t / (1 - t) for the largest
  # t, and the trapezoid discretization, which the change between the nested
  # grids measures (it overstates the error of the selected, finer grid).
  # With 'normalization_prob' 1 - 1e-8 every row loses at most 1e-8.
  error <- result[["posterior"]] / mean(expected) - 1
  expect_identical(diagnostic[["normalization_truncation_status"]], "exact")
  expect_lte(diagnostic[["normalization_truncation"]],
             (1 - (1 - 1e-8)) * (1 + 1e-6))
  expect_lte(error, diagnostic[["truncation_ordinate_bound"]] +
    diagnostic[["ordinate_relative_change"]])
  expect_gte(error, -diagnostic[["ordinate_relative_change"]])
  expect_lt(abs(error), 1e-6)

  # The coverage is what 'normalization_prob' sets: a lower probability leaves
  # more truncation, moves the ordinate, and stays within its reported bound.
  loose <- ordinate(.99)
  loose_diagnostic <- density_diagnostics(loose)
  loose_error <- loose[["posterior"]] / mean(expected) - 1
  expect_gt(loose_diagnostic[["normalization_truncation"]], 1e-4)
  expect_gt(loose_error, 1e-5)
  expect_lte(loose_error, loose_diagnostic[["truncation_ordinate_bound"]] +
    loose_diagnostic[["ordinate_relative_change"]])
  # The public prior height uses BayesTools' shared linear-density interpolation
  # with a 1e-4 relative refinement target, unlike the analytic reference here.
  expect_lt(abs(as.numeric(result[["prior"]]) / stats::dnorm(.1, prior_mean, prior_sd) - 1), 1e-4)
  expect_lt(abs(attr(result, "raw_BF") / (stats::dnorm(.1, prior_mean, prior_sd) / mean(expected)) - 1), 1e-4)
})
