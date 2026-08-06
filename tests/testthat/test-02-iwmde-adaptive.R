context("IWMDE adaptive ordinates and public diagnostics")

source(testthat::test_path("common-functions.R"))


.iwmde_adaptive_test_diagnostics <- function(row_budget, relative_mcse,
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
  expect_equal(ordinate[["max_samples"]], Inf)
  expect_equal(ordinate[["initial_samples"]], 500L)
  expect_equal(ordinate[["target_relative_mcse"]], .05)
  expect_equal(.iwmde_bf_warning_relative_mcse(), .05)
  expect_equal(.iwmde_bf_max_relative_mcse(), .25)
  expect_equal(.iwmde_bf_min_ess(), 20)
  expect_equal(.iwmde_bf_max_weight_share(), .50)

  expect_error(
    .density_control_normalize(
      "qCMDE",
      list(initial_samples = 100, max_samples = 50),
      purpose = "ordinate"
    ),
    "cannot exceed"
  )
})


test_that("adaptive ordinate evaluation records cap exhaustion", {

  skip_on_cran()
  skip_if_missing_fits("bcg_meta-analysis")

  estimate <- .iwmde_estimate(
    context         = .iwmde_context(load_fit("bcg_meta-analysis")),
    parameter       = "mu",
    density_method  = "IWMDE",
    density_control = list(
      n_points             = 20,
      initial_samples      = 20,
      max_samples          = 40,
      target_relative_mcse = 1e-6,
      normalization_points = 20,
      display_grid         = "ordinate"
    ),
    outputs         = "ordinate",
    values          = 0,
    parameter_spec  = list(type = "primitive"),
    metadata        = list(parameter = "mu"),
    cache           = .iwmde_estimate_cache()
  )

  adaptation <- estimate[["adaptation"]]
  expect_true(adaptation[["adaptive"]])
  expect_equal(adaptation[["initial_row_budget"]], 20L)
  expect_equal(adaptation[["achieved_row_budget"]], 40L)
  expect_true(adaptation[["hard_cap_reached"]])
  expect_false(adaptation[["all_rows_used"]])
  expect_false(adaptation[["target_met"]])
  expect_equal(adaptation[["n_steps"]], 2L)
  expect_equal(adaptation[["history"]][["requested_row_budget"]], c(20L, 40L))
})


test_that("adaptive qCMDE stops when the precision target is met", {

  testthat::local_mocked_bindings(
    .iwmde_plan = function(context, parameter, density_method,
                           density_control, outputs, values,
                           parameter_spec, metadata, row_budget) {

      list(
        row_budget = row_budget,
        rows = list(
          continuous_rows = seq_len(100L),
          n_candidate_rows = row_budget
        )
      )
    },
    .iwmde_estimate_from_plan = function(context, plan, cache = NULL) {

      relative_mcse <- if (plan[["row_budget"]] < 80L) .20 else .04
      list(
        diagnostics = list(
          ordinate = list(
            status      = "ok",
            diagnostics = .iwmde_adaptive_test_diagnostics(
              row_budget   = plan[["row_budget"]],
              relative_mcse = relative_mcse
            )
          )
        ),
        posterior_ordinate = list(diagnostics = list())
      )
    },
    .package = "RoBMA"
  )

  estimate <- .iwmde_estimate_adaptive_ordinate(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      initial_samples      = 20L,
      max_samples          = 100L,
      target_relative_mcse = .05
    ),
    values          = 0
  )
  adaptation <- estimate[["adaptation"]]

  expect_equal(
    adaptation[["history"]][["requested_row_budget"]],
    c(20L, 80L)
  )
  expect_equal(adaptation[["achieved_row_budget"]], 80L)
  expect_equal(adaptation[["n_steps"]], 2L)
  expect_true(adaptation[["target_met"]])
  expect_false(adaptation[["hard_cap_reached"]])
  expect_false(adaptation[["all_rows_used"]])
})


test_that("adaptive qCMDE continues until budget-rescuable BF gates pass", {

  testthat::local_mocked_bindings(
    .iwmde_plan = function(context, parameter, density_method,
                           density_control, outputs, values,
                           parameter_spec, metadata, row_budget) {

      list(
        row_budget = row_budget,
        rows = list(
          continuous_rows = seq_len(100L),
          n_candidate_rows = row_budget
        )
      )
    },
    .iwmde_estimate_from_plan = function(context, plan, cache = NULL) {

      passing <- plan[["row_budget"]] >= 40L
      diagnostics <- .iwmde_adaptive_test_diagnostics(
        row_budget       = plan[["row_budget"]],
        relative_mcse    = .20,
        ess              = if (passing) 30 else 10,
        max_weight_share = if (passing) .10 else .60
      )
      diagnostics[["bf_finite_terms"]] <- if (passing) 40L else 10L

      list(
        diagnostics = list(
          ordinate = list(status = "ok", diagnostics = diagnostics)
        ),
        posterior_ordinate = list(diagnostics = list())
      )
    },
    .package = "RoBMA"
  )

  estimate <- .iwmde_estimate_adaptive_ordinate(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      initial_samples      = 20L,
      max_samples          = 100L,
      target_relative_mcse = .50
    ),
    values          = 0
  )
  adaptation <- estimate[["adaptation"]]

  expect_equal(
    adaptation[["history"]][["requested_row_budget"]],
    c(20L, 40L)
  )
  expect_true(all(adaptation[["history"]][["precision_target_met"]]))
  expect_equal(adaptation[["history"]][["bf_grade_met"]], c(FALSE, TRUE))
  expect_true(adaptation[["precision_target_met"]])
  expect_true(adaptation[["bf_grade_met"]])
  expect_true(adaptation[["target_met"]])
})


test_that("adaptive ordinates do not stop on selected-row MCSE alone", {

  testthat::local_mocked_bindings(
    .iwmde_plan = function(context, parameter, density_method,
                           density_control, outputs, values,
                           parameter_spec, metadata, row_budget) {
      list(
        row_budget = row_budget,
        rows = list(
          continuous_rows = seq_len(100L),
          n_candidate_rows = row_budget
        )
      )
    },
    .iwmde_estimate_from_plan = function(context, plan, cache = NULL) {
      sampling_error <- if (plan[["row_budget"]] < 80L) .20 else .04
      list(
        diagnostics = list(ordinate = list(
          status = "ok",
          diagnostics = .iwmde_adaptive_test_diagnostics(
            row_budget = plan[["row_budget"]],
            relative_mcse = .01,
            sampling_relative_mcse = sampling_error
          )
        )),
        posterior_ordinate = list(diagnostics = list())
      )
    },
    .package = "RoBMA"
  )

  estimate <- .iwmde_estimate_adaptive_ordinate(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      initial_samples      = 20L,
      max_samples          = 100L,
      target_relative_mcse = .05
    ),
    values          = 0
  )

  expect_equal(
    estimate[["adaptation"]][["history"]][["requested_row_budget"]],
    c(20L, 80L)
  )
  expect_true(estimate[["adaptation"]][["sampling_target_met"]])
  expect_true(estimate[["adaptation"]][["target_met"]])
})


test_that("adaptive ordinates recover per-chain batch coverage", {

  testthat::local_mocked_bindings(
    .iwmde_plan = function(context, parameter, density_method,
                           density_control, outputs, values,
                           parameter_spec, metadata, row_budget) {
      list(
        row_budget = row_budget,
        rows = list(
          continuous_rows = seq_len(100L),
          n_candidate_rows = row_budget
        )
      )
    },
    .iwmde_estimate_from_plan = function(context, plan, cache = NULL) {
      row_budget <- plan[["row_budget"]]
      chain_one_rows <- if (row_budget < 40L) 3L else 4L
      contributions <- matrix(1, nrow = 1L, ncol = row_budget)
      attr(contributions, "chain_id") <- c(
        rep(1L, chain_one_rows),
        rep(2L, row_budget - chain_one_rows)
      )
      attr(contributions, "expected_chain_ids") <- c(1L, 2L)
      attr(contributions, "target") <- 1
      mcse <- .iwmde_batch_mcse(contributions)

      list(
        diagnostics = list(ordinate = list(
          status = "ok",
          diagnostics = .iwmde_adaptive_test_diagnostics(
            row_budget    = row_budget,
            relative_mcse = mcse[["relative_mcse"]][[1L]],
            ess           = mcse[["ess"]][[1L]]
          )
        )),
        posterior_ordinate = list(diagnostics = list())
      )
    },
    .package = "RoBMA"
  )

  estimate <- .iwmde_estimate_adaptive_ordinate(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      initial_samples      = 20L,
      max_samples          = 100L,
      target_relative_mcse = .05
    ),
    values          = 0
  )

  expect_equal(
    estimate[["adaptation"]][["history"]][["requested_row_budget"]],
    c(20L, 40L)
  )
  expect_equal(
    estimate[["adaptation"]][["history"]][["precision_target_met"]],
    c(FALSE, TRUE)
  )
  expect_equal(
    estimate[["adaptation"]][["history"]][["bf_grade_met"]],
    c(FALSE, TRUE)
  )
  expect_true(estimate[["adaptation"]][["target_met"]])
})


test_that("adaptive mixed ordinates require the active-row census", {

  conditioned_rows <- seq_len(80L)
  conditioned_chain_id <- rep(1:2, each = 40L)
  active_rows <- c(seq(1L, 39L, by = 2L), seq(41L, 79L, by = 2L))
  scopes <- character()

  testthat::local_mocked_bindings(
    .iwmde_plan = function(context, parameter, density_method,
                           density_control, outputs, values,
                           parameter_spec, metadata, row_budget) {
      list(
        row_budget = row_budget,
        rows = list(
          continuous_rows = active_rows,
          n_candidate_rows = row_budget
        )
      )
    },
    .iwmde_estimate_from_plan = function(context, plan, cache = NULL) {
      row_budget <- plan[["row_budget"]]
      selected_rows <- active_rows[seq_len(row_budget)]
      density <- .iwmde_density_aggregate(
        log_terms               = matrix(log(2), nrow = 1L, ncol = row_budget),
        active_mass             = .5,
        denominator             = row_budget,
        contribution_rows       = selected_rows,
        sampling_population_rows = active_rows,
        chain_id                = conditioned_chain_id[selected_rows],
        expected_chain_ids      = 1:2,
        conditioned_rows        = conditioned_rows,
        conditioned_chain_id    = conditioned_chain_id
      )
      mcse <- .iwmde_batch_mcse(density[["mcmc_contributions"]])
      scopes <<- c(scopes, mcse[["uncertainty_scope"]])

      list(
        diagnostics = list(ordinate = list(
          status = "ok",
          diagnostics = .iwmde_adaptive_test_diagnostics(
            row_budget    = row_budget,
            relative_mcse = mcse[["relative_mcse"]][[1L]],
            ess           = mcse[["ess"]][[1L]],
            sampling_relative_mcse =
              density[["sampling_relative_mcse"]][[1L]]
          )
        )),
        posterior_ordinate = list(diagnostics = list())
      )
    },
    .package = "RoBMA"
  )

  estimate <- .iwmde_estimate_adaptive_ordinate(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      initial_samples      = 20L,
      max_samples          = 40L,
      target_relative_mcse = .05
    ),
    values          = 0
  )

  expect_equal(
    estimate[["adaptation"]][["history"]][["requested_row_budget"]],
    c(20L, 40L)
  )
  expect_equal(
    estimate[["adaptation"]][["history"]][["precision_target_met"]],
    c(FALSE, TRUE)
  )
  expect_equal(
    scopes,
    c("unavailable_incomplete_active_census", "full_conditioned_rows")
  )
  expect_true(estimate[["adaptation"]][["target_met"]])
  expect_true(estimate[["adaptation"]][["all_rows_used"]])
})


test_that("density_diagnostics exposes compact BF-grade diagnostics", {

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
    achieved_row_budget               = 500L,
    eligible_rows                     = 1000L,
    hard_cap                          = Inf,
    hard_cap_reached                  = FALSE,
    all_rows_used                     = FALSE,
    n_steps                           = 1L,
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
    "evaluation_value", "achieved_row_budget", "eligible_rows",
    "evaluated_rows", "finite_terms",
    "active_mass", "relative_mcse", "sampling_relative_mcse",
    "sampling_fraction", "mcmc_uncertainty_scope",
    "sampling_uncertainty_type", "ess", "max_weight_share",
    "normalization_relative_error", "stability_metric",
    "stability_relative_error", "ordinate_relative_change",
    "quadrature_relative_change", "target_relative_mcse",
    "stability_warning_threshold", "stability_rejection_threshold",
    "quadrature_warning_threshold", "quadrature_rejection_threshold",
    "warning_relative_mcse", "rejection_relative_mcse",
    "warning_min_finite_terms", "rejection_min_finite_terms",
    "warning_min_ess", "rejection_min_ess", "warning_max_weight_share",
    "rejection_max_weight_share", "hard_cap", "hard_cap_reached",
    "all_rows_used", "adaptation_steps", "target_met",
    "precision_target_met", "sampling_target_met", "bf_grade_met",
    "n_weight_fallbacks",
    "weight_fallback_reasons", "status", "warnings"
  ))
  expect_equal(nrow(out), 1L)
  expect_equal(out[["estimator"]], "iwmde")
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
  expect_equal(out[["rejection_relative_mcse"]], .25)
  expect_true(out[["precision_target_met"]])
  expect_true(out[["bf_grade_met"]])
  expect_equal(out[["n_weight_fallbacks"]], 1L)
  expect_match(out[["weight_fallback_reasons"]], "singular_covariance=1")
  expect_equal(out[["status"]], "ok")
})


test_that("rejected ordinate errors retain public diagnostics", {

  diagnostics <- .iwmde_adaptive_test_diagnostics(
    row_budget       = 40L,
    relative_mcse    = .30,
    ess              = 10,
    max_weight_share = .60
  )
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
