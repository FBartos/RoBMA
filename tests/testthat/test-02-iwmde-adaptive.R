context("IWMDE fixed-budget ordinates and public diagnostics")

source(testthat::test_path("common-functions.R"))


.iwmde_fixed_test_diagnostics <- function(row_budget, relative_mcse,
                                          ess = row_budget,
                                          max_weight_share = 1 / row_budget,
                                          sampling_relative_mcse = 0) {

  return(list(
    bf_relative_mcse             = relative_mcse,
    bf_sampling_relative_mcse    = sampling_relative_mcse,
    bf_finite_terms              = row_budget,
    bf_ess                       = ess,
    bf_max_weight_share          = max_weight_share,
    bf_ordinate_relative_change  = 0,
    max_quadrature_relative_change = 0,
    estimator                    = "q_grid_cmde"
  ))
}


test_that("density and ordinate controls have separate row policies", {

  density_qcmde <- .density_control_normalize("qCMDE", purpose = "density")
  density_iwmde <- .density_control_normalize("IWMDE", purpose = "density")
  ordinate      <- .density_control_normalize("qCMDE", purpose = "ordinate")

  expect_equal(density_qcmde[["max_samples"]], 500L)
  expect_equal(density_iwmde[["max_samples"]], 1000L)
  expect_equal(ordinate[["samples"]], 500L)
  expect_equal(ordinate[["target_relative_mcse"]], .05)
  expect_equal(.iwmde_bf_warning_relative_mcse(), .05)

  expect_error(
    .density_control_normalize(
      "qCMDE",
      list(samples = 19),
      purpose = "ordinate"
    ),
    "integer at least 20 or Inf",
    fixed = TRUE
  )
  expect_error(
    .density_control_normalize(
      "qCMDE",
      list(initial_samples = 100),
      purpose = "ordinate"
    ),
    "unrecognized setting",
    fixed = TRUE
  )
})


test_that("fixed-budget ordinates evaluate one contribution-independent plan", {

  n_plan_calls     <- 0L
  n_estimate_calls <- 0L
  testthat::local_mocked_bindings(
    .iwmde_plan = function(context, parameter, density_method,
                           density_control, outputs, values,
                           parameter_spec, metadata) {

      n_plan_calls <<- n_plan_calls + 1L
      list(
        rows = list(
          continuous_rows = seq_len(100L),
          n_candidate_rows = density_control[["samples"]]
        )
      )
    },
    .iwmde_estimate_from_plan = function(context, plan, cache = NULL) {

      n_estimate_calls <<- n_estimate_calls + 1L
      row_budget <- plan[["rows"]][["n_candidate_rows"]]
      diagnostics <- .iwmde_fixed_test_diagnostics(
        row_budget    = row_budget,
        relative_mcse = .04
      )
      list(
        diagnostics = list(
          ordinate = list(
            status      = "ok",
            diagnostics = diagnostics
          )
        ),
        posterior_ordinate = list(
          value       = 0,
          ordinate    = 1,
          diagnostics = diagnostics
        )
      )
    },
    .package = "RoBMA"
  )

  estimate <- .iwmde_estimate_fixed_ordinate(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      samples              = 80L,
      target_relative_mcse = .05
    ),
    values          = 0
  )
  design <- estimate[["sampling_design"]]

  expect_equal(n_plan_calls, 1L)
  expect_equal(n_estimate_calls, 1L)
  expect_true(design[["fixed_budget"]])
  expect_equal(design[["requested_samples"]], 80L)
  expect_equal(design[["achieved_row_budget"]], 80L)
  expect_true(design[["precision_target_met"]])
  expect_true(design[["sampling_target_met"]])
  expect_true(design[["target_met"]])
  expect_false(design[["all_rows_used"]])
})


test_that("fixed-budget diagnostics never select ordinate availability", {

  n_estimate_calls <- 0L
  sampling_error   <- 0
  testthat::local_mocked_bindings(
    .iwmde_plan = function(context, parameter, density_method,
                           density_control, outputs, values,
                           parameter_spec, metadata) {

      list(
        rows = list(
          continuous_rows = seq_len(100L),
          n_candidate_rows = density_control[["samples"]]
        )
      )
    },
    .iwmde_estimate_from_plan = function(context, plan, cache = NULL) {

      n_estimate_calls <<- n_estimate_calls + 1L
      row_budget <- plan[["rows"]][["n_candidate_rows"]]
      diagnostics <- .iwmde_fixed_test_diagnostics(
        row_budget             = row_budget,
        relative_mcse          = .01,
        sampling_relative_mcse = sampling_error
      )
      list(
        diagnostics = list(
          ordinate = list(
            status = "ok",
            diagnostics = diagnostics
          )
        ),
        posterior_ordinate = list(
          value       = 0,
          ordinate    = if (sampling_error == 0) 1 else 2,
          diagnostics = diagnostics
        )
      )
    },
    .package = "RoBMA"
  )

  estimate_precise <- .iwmde_estimate_fixed_ordinate(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      samples              = 20L,
      target_relative_mcse = .05
    ),
    values          = 0
  )
  sampling_error <- .10
  estimate_imprecise <- .iwmde_estimate_fixed_ordinate(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      samples              = 20L,
      target_relative_mcse = .05
    ),
    values          = 0
  )

  expect_equal(n_estimate_calls, 2L)
  expect_true(is.list(estimate_precise[["posterior_ordinate"]]))
  expect_true(is.list(estimate_imprecise[["posterior_ordinate"]]))
  expect_null(estimate_precise[["rejected_posterior_ordinate"]])
  expect_null(estimate_imprecise[["rejected_posterior_ordinate"]])
  expect_equal(estimate_precise[["posterior_ordinate"]][["ordinate"]], 1)
  expect_equal(estimate_imprecise[["posterior_ordinate"]][["ordinate"]], 2)
  expect_true(estimate_precise[["sampling_design"]][["sampling_target_met"]])
  expect_false(estimate_imprecise[["sampling_design"]][["sampling_target_met"]])
  expect_match(
    paste(.iwmde_posterior_ordinate_warnings(
      estimate_imprecise[["posterior_ordinate"]]
    ), collapse = " "),
    "via the 'density_control' argument",
    fixed = TRUE
  )
})


test_that("fixed-budget census imprecision warns to obtain more draws", {

  testthat::local_mocked_bindings(
    .iwmde_plan = function(context, parameter, density_method,
                           density_control, outputs, values,
                           parameter_spec, metadata) {
      list(
        rows = list(
          continuous_rows = seq_len(20L),
          n_candidate_rows = 20L
        )
      )
    },
    .iwmde_estimate_from_plan = function(context, plan, cache = NULL) {
      diagnostics <- .iwmde_fixed_test_diagnostics(
        row_budget    = 20L,
        relative_mcse = .10
      )
      list(
        diagnostics = list(ordinate = list(
          status = "ok",
          diagnostics = diagnostics
        )),
        posterior_ordinate = list(diagnostics = diagnostics)
      )
    },
    .package = "RoBMA"
  )

  estimate <- .iwmde_estimate_fixed_ordinate(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      samples              = Inf,
      target_relative_mcse = .05
    ),
    values          = 0
  )

  expect_true(is.list(estimate[["posterior_ordinate"]]))
  expect_match(
    paste(.iwmde_posterior_ordinate_warnings(
      estimate[["posterior_ordinate"]]
    ), collapse = " "),
    "increasing the number of posterior draws",
    fixed = TRUE
  )
})


test_that("self-weighting thinning preserves rare-row mass in expectation", {

  population <- c(rep(1, 9L), 1001)
  population_rows <- seq_along(population)
  samples <- utils::combn(population_rows, 4L, simplify = FALSE)
  thinned <- vapply(samples, function(rows) {
    .iwmde_density_aggregate(
      log_terms                = matrix(log(population[rows]), nrow = 1L),
      active_mass              = 1,
      denominator              = length(rows),
      contribution_rows        = rows,
      sampling_population_rows = population_rows,
      chain_id                 = rep(1L, length(rows)),
      expected_chain_ids       = 1L
    )[["y"]][[1L]]
  }, numeric(1))
  census <- .iwmde_density_aggregate(
    log_terms                = matrix(log(population), nrow = 1L),
    active_mass              = 1,
    denominator              = length(population),
    contribution_rows        = population_rows,
    sampling_population_rows = population_rows,
    chain_id                 = rep(1L, length(population)),
    expected_chain_ids       = 1L
  )[["y"]][[1L]]

  expect_equal(mean(thinned), census, tolerance = 1e-12)
  expect_equal(census, mean(population), tolerance = 1e-12)
})


test_that("density_diagnostics exposes fixed-sample numerical diagnostics", {

  diagnostics <- list(
    evaluation_value                  = 0,
    relative_mcse                     = .03,
    finite_terms                      = 500L,
    ess                               = 250,
    max_weight_share                  = .02,
    active_mass                       = 1,
    n_candidate_rows                  = 500L,
    n_evaluated_rows                  = 500L,
    sampling_relative_mcse            = .02,
    sampling_fraction                 = .5,
    mcmc_uncertainty_scope            = "full_conditioned_rows",
    mcmc_uncertainty_status           = "available",
    mcmc_uncertainty_reason           = NULL,
    sampling_uncertainty_type         = "finite_population_srswor",
    normalization_relative_error      = .001,
    ordinate_relative_change          = .002,
    max_quadrature_relative_change    = .003,
    target_relative_mcse              = .05,
    requested_samples                 = 500L,
    achieved_row_budget               = 500L,
    eligible_rows                     = 1000L,
    all_rows_used                     = FALSE,
    target_met                        = TRUE,
    precision_target_met              = TRUE,
    sampling_target_met               = TRUE,
    bf_grade_met                      = TRUE,
    n_weight_fallbacks                = 1L,
    weight_fallback_reasons           = c(singular_covariance = 1L),
    estimator                         = "iwmde",
    weight_method                     = "chen_marginal_normal"
  )
  ordinate <- BayesTools::posterior_ordinate_attribute(
    value          = 0,
    ordinate       = .4,
    method         = "iwmde",
    density_method = "IWMDE",
    diagnostics    = diagnostics,
    parameter      = "mu"
  )
  ordinate[["iwmde_provenance"]] <- list(
    request_key        = "request",
    schema_version     = .iwmde_schema_version(),
    algorithm_version  = .iwmde_algorithm_version(),
    source_fingerprint = list(draws = "hash")
  )
  posterior <- stats::setNames(1:10, paste0("draw", 1:10))
  attr(posterior, "posterior_ordinate") <- ordinate
  table <- structure(
    data.frame(BF = 1, row.names = "mu = 0"),
    class = c("BayesTools_hypothesis_BF", "data.frame")
  )
  table <- .hypothesis_brma_append_iwmde_warnings(table, posterior)
  out   <- density_diagnostics(table)

  expect_identical(class(out), c("RoBMA_density_diagnostics", "data.frame"))
  expect_named(out, c(
    "schema_version", "algorithm_version", "source_fingerprint", "estimator",
    "density_method", "parameter", "level", "requested_value",
    "evaluation_value", "requested_samples", "achieved_row_budget",
    "eligible_rows",
    "evaluated_rows", "finite_terms",
    "active_mass", "relative_mcse", "sampling_relative_mcse",
    "sampling_fraction", "mcmc_uncertainty_scope",
    "sampling_uncertainty_type", "ess", "max_weight_share",
    "normalization_relative_error", "stability_metric",
    "stability_relative_error", "ordinate_relative_change",
    "quadrature_relative_change", "target_relative_mcse",
    "stability_warning_threshold", "stability_rejection_threshold",
    "quadrature_warning_threshold", "quadrature_rejection_threshold",
    "warning_relative_mcse", "warning_min_finite_terms",
    "warning_min_ess", "warning_max_weight_share", "all_rows_used",
    "target_met",
    "precision_target_met", "sampling_target_met", "bf_grade_met",
    "n_weight_fallbacks",
    "weight_fallback_reasons", "status", "warnings"
  ))
  expect_equal(nrow(out), 1L)
  expect_equal(out[["estimator"]], "iwmde")
  expect_equal(out[["requested_samples"]], 500)
  expect_equal(out[["achieved_row_budget"]], 500L)
  expect_equal(out[["relative_mcse"]], .03)
  expect_equal(out[["sampling_relative_mcse"]], .02)
  expect_equal(out[["sampling_fraction"]], .5)
  expect_equal(
    out[["mcmc_uncertainty_scope"]],
    "full_conditioned_rows"
  )
  expect_true(out[["sampling_target_met"]])
  expect_equal(out[["stability_metric"]], "normalization_relative_error")
  expect_equal(out[["stability_warning_threshold"]], .05)
  expect_equal(out[["stability_rejection_threshold"]], .10)
  expect_equal(out[["quadrature_warning_threshold"]], .025)
  expect_equal(out[["quadrature_rejection_threshold"]], .05)
  expect_equal(out[["warning_relative_mcse"]], .05)
  expect_true(out[["precision_target_met"]])
  expect_true(out[["bf_grade_met"]])
  expect_equal(out[["n_weight_fallbacks"]], 1L)
  expect_match(out[["weight_fallback_reasons"]], "singular_covariance=1")
  expect_equal(out[["status"]], "ok")
})


test_that("rejected ordinate errors retain public diagnostics", {

  diagnostics <- .iwmde_fixed_test_diagnostics(
    row_budget       = 40L,
    relative_mcse    = .30,
    ess              = 10,
    max_weight_share = .60
  )
  diagnostics[["bf_ordinate_relative_change"]] <- .06
  ordinate <- BayesTools::posterior_ordinate_attribute(
    value          = 0,
    ordinate       = .4,
    method         = "q_grid_cmde",
    density_method = "qCMDE",
    diagnostics    = diagnostics,
    parameter      = "mu"
  )
  estimate <- list(
    rejected_posterior_ordinate = ordinate,
    diagnostics = list(ordinate = list(status = "ok"))
  )
  error <- tryCatch(
    .iwmde_stop_ordinate_unavailable("ordinate rejected", estimate),
    error = identity
  )

  expect_s3_class(error, "RoBMA_density_ordinate_error")
  expect_match(conditionMessage(error), "ordinate rejected")
  expect_false("estimate" %in% names(error))
  out <- density_diagnostics(error)
  expect_s3_class(out, "RoBMA_density_diagnostics")
  expect_equal(nrow(out), 1L)
  expect_equal(out[["status"]], "rejected")
})


test_that("density_diagnostics uses explicit methods and validates its schema", {

  diagnostics <- .iwmde_empty_public_density_diagnostics()
  generic_object <- structure(list(), density_diagnostics = diagnostics)
  malformed <- structure(
    data.frame(BF = 1),
    class = c("BayesTools_hypothesis_BF", "data.frame"),
    density_diagnostics = data.frame(status = "ok")
  )

  expect_error(
    density_diagnostics(generic_object),
    "No density diagnostics method"
  )
  expect_error(
    density_diagnostics(malformed),
    "public schema"
  )
  expect_warning(
    density_diagnostics(
      structure(
        data.frame(BF = 1),
        class = c("BayesTools_hypothesis_BF", "data.frame"),
        density_diagnostics = diagnostics
      ),
      obsolete = TRUE
    ),
    "Unused argument.*'obsolete'"
  )
})


test_that("point-BF policy enforces adaptive quadrature stability", {

  diagnostics <- list(
    estimator                       = "q_grid_cmde",
    relative_mcse                   = .01,
    sampling_relative_mcse          = .01,
    finite_terms                    = 500L,
    ess                             = 250,
    max_weight_share                = .01,
    ordinate_relative_change        = .001,
    max_quadrature_relative_change  = .03
  )

  warning <- .iwmde_diagnostics_bf_warning(diagnostics)
  expect_match(warning, "quadrature sensitivity")
  expect_null(.iwmde_diagnostics_bf_failure_reason(diagnostics))

  diagnostics[["max_quadrature_relative_change"]] <- .06
  expect_match(
    .iwmde_diagnostics_bf_failure_reason(diagnostics),
    "quadrature sensitivity"
  )
})
