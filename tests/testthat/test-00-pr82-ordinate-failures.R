.pr82_failure_plan <- function(value = 0) {

  list(target = list(parameter = "mu", metadata = list(parameter = "mu")),
    outputs = list(requested_values = value), grids = list(requested_values = value),
    density_method = "qCMDE", method = "q_grid_cmde",
    control = list(samples = 40L, n_points = 20L))
}

.pr82_failure_diagnostic <- function(value = 0, ordinate = 0, log_ordinate = -1000) {

  list(status = "ok", diagnostics = list(bf_included = TRUE,
    bf_value = value, bf_evaluation_value = value, bf_ordinate = ordinate,
    bf_log_ordinate = log_ordinate, bf_ordinate_relative_change = .3,
    bf_mcse = 0, bf_relative_mcse = NA_real_, bf_finite_terms = 40L,
    n_evaluated_rows = 40L, estimator = "q_grid_cmde"))
}

test_that("density aggregation retains actual logs before ordinary-scale loss", {

  for (log_value in c(-10, -1000, -Inf)) {
    out <- .iwmde_density_aggregate(matrix(log_value, 1L, 40L), 1, 40L)
    expect_equal(out[["log_y"]], log_value, tolerance = 1e-13)
    expect_equal(out[["y"]], exp(log_value), tolerance = 1e-13)
  }
  error <- tryCatch(.iwmde_density_aggregate(matrix(1000, 1L, 40L), 1, 40L,
    evaluation_values = 0), error = identity)
  expect_s3_class(error, "iwmde_ordinate_numerical_error")
  expect_equal(error[["log_ordinate"]], 1000, tolerance = 1e-13)
  expect_identical(error[["ordinate"]], Inf)
  intermediate <- tryCatch(.iwmde_density_aggregate(matrix(1000, 1L, 40L), exp(-500), 40L,
    evaluation_values = 0), error = identity)
  expect_equal(intermediate[["log_ordinate"]], 500, tolerance = 1e-13)
  expect_identical(.iwmde_ordinate_numerical_status(intermediate[["ordinate"]],
    intermediate[["log_ordinate"]]), "arithmetic_failure")
})

test_that("representability requires log evidence and does not grade tiny values", {

  expect_identical(.iwmde_ordinate_numerical_status(0, -1000), "underflow")
  expect_identical(.iwmde_ordinate_numerical_status(Inf, 1000), "overflow")
  expect_identical(.iwmde_ordinate_numerical_status(0, -10), "arithmetic_failure")
  expect_identical(.iwmde_ordinate_numerical_status(Inf, 10), "arithmetic_failure")
  expect_identical(.iwmde_ordinate_numerical_status(0, -Inf), "zero_without_log_evidence")
  expect_identical(.iwmde_ordinate_numerical_status(0), "zero_without_log_evidence")
  expect_identical(.iwmde_ordinate_numerical_status(NA_real_, computed = FALSE), "not_computed")
  expect_identical(.iwmde_ordinate_numerical_status(NA_real_), "nonfinite")
  expect_identical(.iwmde_ordinate_numerical_status(-1, 0), "invalid_density")
  expect_identical(.iwmde_ordinate_numerical_status(1e-300, log(1e-300)), "finite")
})

test_that("mixed requested points retain only the genuinely failed records", {

  accepted <- BayesTools::posterior_ordinate_attribute(0, 1, "q_grid_cmde", "qCMDE")
  diagnostic <- list(status = "ok", parameter = "mu",
    diagnostics = list(estimator = "q_grid_cmde"),
    iwmde = list(x = c(0, 1), y = c(1, 0), log_y = c(0, -1000)))
  records <- .iwmde_estimate_ordinate_failures(.pr82_failure_plan(c(0, 1, 2)), diagnostic, accepted)
  expect_identical(records[["requested_value"]], c(1, 2))
  expect_identical(records[["numerical_status"]], c("underflow", "not_computed"))
  expect_identical(records[["ordinate"]], c(0, NA_real_))
  expect_identical(records[["log_ordinate"]], c(-1000, NA_real_))
  expect_false(any(records[["bf_grade_met"]]))
})

test_that("display transformations preserve computed log-ordinate evidence", {

  source <- list(value = .5, evaluation_value = .5, ordinate = .3,
    diagnostics = list(log_ordinate = log(.3)), iwmde_provenance = list())
  actual <- .hypothesis_brma_transform_iwmde_ordinate(source, list(type = "tanh"))
  expect_equal(actual[["diagnostics"]][["log_ordinate"]],
    log(.3) - log(1 - tanh(.5)^2), tolerance = 1e-14)
  expect_equal(exp(actual[["diagnostics"]][["log_ordinate"]]), actual[["ordinate"]], tolerance = 1e-14)
  failure <- .iwmde_estimate_ordinate_failures(.pr82_failure_plan(.5), .pr82_failure_diagnostic(.5))
  transformed <- .hypothesis_brma_transform_iwmde_failures(failure, list(type = "tanh"))
  expect_equal(transformed[["requested_value"]], tanh(.5))
  expect_equal(transformed[["evaluation_value"]], tanh(.5))
  expect_equal(transformed[["log_ordinate"]], -1000 - log(1 - tanh(.5)^2), tolerance = 1e-14)
  expect_identical(transformed[["numerical_status"]], "underflow")
})

test_that("untyped programming errors are not converted to ordinate failures", {

  testthat::local_mocked_bindings(.iwmde_plan = function(...) stop("untyped programming error"),
    .package = "RoBMA")
  error <- tryCatch(.iwmde_estimate_fixed_ordinate(list(), "mu", "qCMDE",
    list(samples = 20L, n_points = 20L), 0), error = identity)
  expect_identical(class(error), c("simpleError", "error", "condition"))
  expect_null(error[["density_diagnostics"]])
})

test_that("display-Jacobian range loss retains transformed failure evidence", {

  for (direction in c(-1, 1)) {
    source_value <- 10^(direction * 300)
    source_ordinate <- 10^(direction * 100)
    scale <- 10^(-direction * 300)
    source <- BayesTools::posterior_ordinate_attribute(
      source_value, source_ordinate, "q_grid_cmde", "qCMDE",
      evaluation_value = source_value,
      diagnostics = list(estimator = "q_grid_cmde", log_ordinate = log(source_ordinate),
        relative_mcse = .1, precision_target_met = TRUE, ordinate_relative_change = 0))
    error <- tryCatch(.hypothesis_brma_transform_iwmde_ordinate(source,
      list(type = "affine", offset = 0, scale = scale)), error = identity)
    expect_s3_class(error, "RoBMA_density_ordinate_error")
    record <- density_diagnostics(error)
    expect_equal(record[["requested_value"]], 1)
    expect_equal(record[["evaluation_value"]], 1)
    expect_equal(record[["log_ordinate"]], log(source_ordinate) - log(scale), tolerance = 1e-14)
    expect_identical(record[["numerical_status"]], if (direction < 0) "underflow" else "overflow")
    expect_identical(record[["relative_mcse"]], .1)
    expect_false(record[["bf_grade_met"]])
    expect_false(record[["target_met"]])
  }
})

test_that("failed point records survive without a BayesTools ordinate", {

  diagnostic <- .pr82_failure_diagnostic(value = 3)
  result <- .iwmde_estimate_result(.pr82_failure_plan(3), ordinate_diagnostic = diagnostic)
  expect_null(result[["posterior_ordinate"]])
  expect_null(result[["rejected_posterior_ordinate"]])
  failure <- result[["ordinate_failures"]]
  expect_identical(failure[["requested_value"]], 3)
  expect_identical(failure[["ordinate"]], 0)
  expect_identical(failure[["log_ordinate"]], -1000)
  expect_identical(failure[["numerical_status"]], "underflow")
  expect_identical(failure[["status"]], "unavailable")
  expect_false(failure[["bf_grade_met"]])
  expect_false(failure[["precision_target_met"]])
  expect_true(is.na(failure[["ess"]]))
  expect_true(is.na(failure[["ordinate_relative_change"]]))
  expect_false(grepl("increas|more.*samples", failure[["failure_reason"]]))
  expect_identical(failure[["evaluated_rows"]], 40L)
  expect_match(failure[["failure_reason"]], "computed log ordinate -1000", fixed = TRUE)
  expect_error(BayesTools::posterior_ordinate_attribute(3, 0, "q_grid_cmde", "qCMDE"),
    "Posterior ordinates must be finite and positive.", fixed = TRUE)
  error <- tryCatch(.iwmde_stop_ordinate_unavailable("Point computation unavailable", result), error = identity)
  expect_identical(density_diagnostics(error), failure)
})

test_that("marginal means retain target-specific unavailable reasons", {

  estimate <- .iwmde_estimate_result(.pr82_failure_plan(3),
    ordinate_diagnostic = .pr82_failure_diagnostic(3))
  object <- list(inference = list(conditional = list(mu = list(a = 1:3))))
  attached <- .marginal_means_attach_iwmde_ordinate_type(object, "conditional",
    list(a = list(parameter = "mu", level = "a")), list(a = estimate))
  posterior <- attached[["inference"]][["conditional"]][["mu"]][["a"]]
  expect_null(attr(posterior, "posterior_ordinate", exact = TRUE))
  failure <- attr(posterior, "ordinate_failures", exact = TRUE)
  expect_identical(failure[["numerical_status"]], "underflow")
  bf <- .marginal_means_iwmde_bf_scalar(posterior, 3, NULL, "qCMDE", "Unavailable.")
  expect_true(is.na(bf))
  expect_match(attr(bf, "warnings"), "underflowed to zero", fixed = TRUE)
  other <- .marginal_means_iwmde_bf_scalar(posterior, 4, NULL, "qCMDE", "Other target unavailable.")
  expect_identical(attr(other, "warnings"), "Other target unavailable.")
})

test_that("multi-level marginal means keep each failure beside valid Bayes factors", {

  underflow <- .iwmde_estimate_ordinate_failures(.pr82_failure_plan(3), .pr82_failure_diagnostic(3))
  unknown <- .iwmde_estimate_ordinate_failures(.pr82_failure_plan(3),
    .pr82_failure_diagnostic(3, log_ordinate = NA_real_))
  bad <- list(underflow = structure(1:3, ordinate_failures = underflow),
              unknown = structure(1:3, ordinate_failures = unknown))
  result <- .marginal_means_iwmde_bf(bad, 3, provenance = list(), density_method = "qCMDE")
  expect_true(all(vapply(result, is.na, logical(1L))))
  expect_match(attr(result[["underflow"]], "warnings"), "underflowed to zero", fixed = TRUE)
  expect_match(attr(result[["unknown"]], "warnings"), "without a finite log estimate", fixed = TRUE)
  provenance <- .iwmde_provenance_request("qCMDE", "q_grid_cmde", value = 3,
    metadata = list(parameter = "mu", level = "good"), attribute = "ordinate")
  ordinate <- BayesTools::posterior_ordinate_attribute(3, 1, "q_grid_cmde", "qCMDE",
    diagnostics = list(estimator = "q_grid_cmde", ordinate_relative_change = 0),
    iwmde_provenance = provenance)
  posterior <- c(list(good = structure(1:3, posterior_ordinate = ordinate)), bad)
  calls <- 0L
  testthat::local_mocked_bindings(Savage_Dickey_BF = function(posterior, ...) {
    calls <<- calls + 1L
    expect_identical(names(posterior), "good")
    list(good = 2)
  }, .package = "BayesTools")
  mixed <- .marginal_means_iwmde_bf(posterior, 3,
    provenance = list(good = provenance), density_method = "qCMDE")
  expect_identical(calls, 1L)
  expect_equal(as.numeric(mixed[["good"]]), 2)
  expect_match(attr(mixed[["underflow"]], "warnings"), "underflowed to zero", fixed = TRUE)
  expect_match(attr(mixed[["unknown"]], "warnings"), "without a finite log estimate", fixed = TRUE)
})

test_that("public marginal-means hypotheses distinguish unavailable point computations", {

  sample <- structure(seq(-1, 1, length.out = 40L),
    class = c("marginal_posterior.simple", "numeric"), linear_weights = c(mu = 1),
    posterior_atoms = BayesTools::posterior_atom_attribute())
  levels <- structure(list(A = sample), class = c("marginal_posterior.factor", "marginal_posterior", "list"),
    parameter = "mu_alloc")
  object <- structure(list(
    inference = structure(list(averaged = list(mu_alloc = levels), conditional = list(mu_alloc = levels)),
      class = "marginal_inference"),
    term_map = data.frame(term = "alloc", parameter = "mu_alloc", label = "alloc"),
    density_method = "qCMDE", model_averaged = FALSE,
    source_object = structure(list(fit = list(TRUE)), class = "brma")
  ), class = "marginal_means.brma")
  estimate <- .iwmde_estimate_result(.pr82_failure_plan(3),
    ordinate_diagnostic = .pr82_failure_diagnostic(3))
  testthat::local_mocked_bindings(
    .check_iwmde_available = function(...) invisible(TRUE),
    .iwmde_check_point_ordinate_supported = function(...) invisible(TRUE),
    .iwmde_context = function(...) list(),
    .iwmde_request_provenance = function(...) list(),
    .iwmde_estimate = function(...) estimate,
    .package = "RoBMA"
  )
  error <- tryCatch(hypothesis(object, "alloc[A] = 3", density_method = "qCMDE",
    density_control = list(samples = 40L, n_points = 20L)), error = identity)
  expect_s3_class(error, "RoBMA_density_ordinate_error")
  expect_match(conditionMessage(error), "is unavailable", fixed = TRUE)
  expect_match(conditionMessage(error), "underflowed to zero", fixed = TRUE)
  expect_false(grepl("was rejected by diagnostics", conditionMessage(error), fixed = TRUE))
  expect_identical(density_diagnostics(error)[["numerical_status"]], "underflow")
})

test_that("typed construction errors preserve their cause and requested points", {

  original <- tryCatch(.iwmde_stop_construction_failure("q_grid_cmde", "mu", 1:2,
    "joint-density evaluation", "joint log density was undefined"), error = identity)
  enriched <- tryCatch(.iwmde_enrich_ordinate_error(original, .pr82_failure_plan(c(0, 1)),
    list(), "mu", "qCMDE", list(), c(0, 1)), error = identity)
  expect_s3_class(enriched, "iwmde_construction_error")
  expect_s3_class(enriched, "RoBMA_density_ordinate_error")
  expect_identical(conditionMessage(enriched), conditionMessage(original))
  expect_identical(enriched[["stage"]], original[["stage"]])
  records <- density_diagnostics(enriched)
  expect_identical(records[["numerical_status"]], rep("not_computed", 2L))
  expect_true(all(is.na(records[["log_ordinate"]])))
  expect_match(records[["failure_reason"]][[1L]], "joint-density evaluation", fixed = TRUE)
})

test_that("legacy diagnostic tables remain readable without invented log evidence", {

  current <- .iwmde_estimate_ordinate_failures(.pr82_failure_plan(), .pr82_failure_diagnostic())
  legacy <- current[setdiff(names(current), c("ordinate", "log_ordinate", "numerical_status", "failure_reason"))]
  upgraded <- .density_diagnostics_validate(legacy)
  expect_identical(upgraded[["requested_value"]], current[["requested_value"]])
  expect_true(is.na(upgraded[["ordinate"]]))
  expect_true(is.na(upgraded[["log_ordinate"]]))
  expect_true(is.na(upgraded[["numerical_status"]]))
})

test_that("declared prior boundaries do not invent posterior numerical classifications", {

  plan <- .pr82_failure_plan()
  plan[["prior_ordinates"]] <- .iwmde_prior_ordinate_classifications(
    BayesTools::prior("gamma", list(shape = .5, rate = 1)), 0
  )
  records <- .iwmde_estimate_ordinate_failures(plan,
    list(status = "unsupported", reason = "No posterior point computation was obtained"))
  expect_identical(records[["numerical_status"]], "not_computed")
  expect_true(is.na(records[["log_ordinate"]]))
  expect_match(records[["warnings"]], "target prior density is singular", fixed = TRUE)
})

test_that("transformed hypotheses preserve early typed errors on the displayed scale", {

  record <- .iwmde_estimate_ordinate_failures(.pr82_failure_plan(2), .pr82_failure_diagnostic(2))
  original <- structure(list(message = "Original typed computation failure.", call = NULL,
    stage = "density aggregation", density_diagnostics = record),
    class = c("RoBMA_density_ordinate_error", "iwmde_ordinate_numerical_error", "error", "condition"))
  testthat::local_mocked_bindings(.iwmde_estimate = function(...) stop(original), .package = "RoBMA")
  error <- tryCatch(.hypothesis_brma_attach_iwmde_scalar(
    posterior = c(1, 4), raw_posterior = c(1, 2), context = list(), estimate_cache = NULL,
    parameter = "tau", parameter_label = "tau2", value = 4, conditional = FALSE,
    n_points = 20L, samples = 20L, target_relative_mcse = .05,
    normalization_points = 20L, normalization_prob = .999, density_method = "qCMDE",
    parameter_spec = list(type = "primitive"), display_transform = list(type = "square")
  ), error = identity)
  expect_identical(class(error), class(original))
  expect_identical(conditionMessage(error), conditionMessage(original))
  expect_identical(error[["stage"]], original[["stage"]])
  displayed <- density_diagnostics(error)
  expect_identical(displayed[["requested_value"]], 4)
  expect_identical(displayed[["evaluation_value"]], 4)
  expect_equal(displayed[["log_ordinate"]], -1000 - log(4), tolerance = 1e-14)
})
