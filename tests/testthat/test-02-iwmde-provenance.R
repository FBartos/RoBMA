context("IWMDE provenance")

source(testthat::test_path("common-functions.R"))


.iwmde_provenance_test_context <- function(samples) {

  return(list(
    posterior_samples = samples,
    object = list(
      fit        = structure(list(), class = "mock_fit"),
      likelihood = list(family = "normal"),
      data       = list(measure = "GEN")
    ),
    data           = list(yi = c(0.1, 0.2)),
    priors         = list(mu = "normal"),
    selection_spec = list(active = FALSE)
  ))
}


test_that("source fingerprint covers every posterior draw", {

  samples <- matrix(
    seq_len(36),
    nrow     = 12,
    ncol     = 3,
    dimnames = list(NULL, c("mu", "tau", "eta[1]"))
  )
  changed_samples       <- samples
  changed_samples[6, 2] <- changed_samples[6, 2] + 1

  fingerprint <- .iwmde_source_fingerprint(
    .iwmde_provenance_test_context(samples)
  )
  changed_fingerprint <- .iwmde_source_fingerprint(
    .iwmde_provenance_test_context(changed_samples)
  )

  expect_false(identical(fingerprint, changed_fingerprint))
  expect_match(fingerprint[["posterior_values"]], "^iwmde_posterior_draws\\|")
})


test_that("IWMDE contexts retain one source fingerprint across sibling uses", {

  calls <- 0L
  testthat::local_mocked_bindings(
    .get_posterior_samples = function(fit, posterior_samples = NULL) {
      matrix(
        seq_len(12),
        nrow     = 4,
        dimnames = list(NULL, c("mu", "tau", "theta[1]"))
      )
    },
    .iwmde_selection_spec = function(data, priors) list(active = FALSE),
    .iwmde_formula_inputs = function(data, priors) list(),
    .iwmde_compute_source_fingerprint = function(context) {
      calls <<- calls + 1L
      list(posterior_values = "computed-once")
    },
    .package = "RoBMA"
  )
  object <- list(
    fit        = structure(list(), prior_list = list()),
    data       = list(measure = "GEN"),
    priors     = list(),
    likelihood = list(family = "normal")
  )

  context <- .iwmde_context(object)
  expect_equal(calls, 1L)
  expect_equal(
    .iwmde_source_fingerprint(context)[["posterior_values"]],
    "computed-once"
  )
  context <- .iwmde_context_ensure_caches(context)
  expect_equal(calls, 1L)
})


test_that("request keys carry schema and algorithm versions", {

  request <- .iwmde_provenance_request(
    density_method = "qCMDE",
    method         = "q_grid_cmde",
    target_key     = "primitive|mu"
  )

  expect_equal(request[["schema_version"]], "6")
  expect_equal(request[["algorithm_version"]], "14")

  testthat::local_mocked_bindings(
    .iwmde_algorithm_version = function() "changed",
    .package                 = "RoBMA"
  )
  changed_request <- .iwmde_provenance_request(
    density_method = "qCMDE",
    method         = "q_grid_cmde",
    target_key     = "primitive|mu"
  )

  expect_equal(changed_request[["algorithm_version"]], "changed")
  expect_false(identical(request[["request_key"]], changed_request[["request_key"]]))
})


test_that("ordinate request keys include prior classification provenance", {

  zero <- .iwmde_prior_ordinate_classifications(
    BayesTools::prior("beta", parameters = list(alpha = 2, beta = 2)),
    0
  )
  regular <- .iwmde_prior_ordinate_classifications(
    BayesTools::prior("normal", parameters = list(mean = 0, sd = 1)),
    0
  )
  zero_request <- .iwmde_provenance_request(
    density_method  = "qCMDE",
    method          = "q_grid_cmde",
    value           = 0,
    attribute       = "ordinate",
    target_key      = "primitive|mu",
    prior_ordinates = zero
  )
  regular_request <- .iwmde_provenance_request(
    density_method  = "qCMDE",
    method          = "q_grid_cmde",
    value           = 0,
    attribute       = "ordinate",
    target_key      = "primitive|mu",
    prior_ordinates = regular
  )

  expect_false(identical(
    zero_request[["request_key"]],
    regular_request[["request_key"]]
  ))
})


test_that("semantic request keys do not depend on achieved plan budgets", {

  first <- .iwmde_provenance_request(
    density_method    = "qCMDE",
    method            = "q_grid_cmde",
    target_key        = "primitive|mu",
    plan_key          = "budget-500",
    density_control   = list(samples = Inf),
    source_fingerprint = list(posterior_values = "draws")
  )
  final <- .iwmde_provenance_request(
    density_method    = "qCMDE",
    method            = "q_grid_cmde",
    target_key        = "primitive|mu",
    plan_key          = "budget-2000",
    density_control   = list(samples = Inf),
    source_fingerprint = list(posterior_values = "draws")
  )

  expect_equal(first[["request_key"]], final[["request_key"]])
  expect_false(identical(first[["plan_key"]], final[["plan_key"]]))
})


test_that("semantic request provenance does not construct posterior-row plans", {

  testthat::local_mocked_bindings(
    .iwmde_context_ensure_caches = identity,
    .iwmde_parameter_spec = function(context, parameter,
                                     parameter_spec = NULL) {
      list(type = "primitive", status = "ok")
    },
    .iwmde_source_fingerprint = function(context) {
      list(posterior_values = "draws")
    },
    .iwmde_plan = function(...) {
      stop("posterior-row plan must not be constructed")
    },
    .package = "RoBMA"
  )

  request <- .iwmde_request_provenance(
    context         = list(),
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = list(
      samples      = 20L,
      display_grid = "ordinate"
    ),
    attribute       = "ordinate",
    value           = 0,
    parameter_spec  = list(type = "primitive")
  )

  expect_match(request[["request_key"]], "^iwmde_request\\|")
  expect_equal(request[["target_key"]], "primitive|mu")
  expect_equal(request[["source_fingerprint"]][["posterior_values"]], "draws")
  expect_false("plan_key" %in% names(request))
})


test_that("internal schemas are constructed and versioned", {

  state <- .iwmde_new_row_state(list(baseline_log_q = 0))
  expect_s3_class(state, "iwmde_row_state")
  expect_equal(state[["schema_version"]], .iwmde_schema_version())
  expect_invisible(.iwmde_validate_row_states(list(state)))

  expect_error(
    .iwmde_validate_row_states(list(list(baseline_log_q = 0))),
    "invalid baseline"
  )
  expect_error(
    .iwmde_new_row_state(list(baseline_log_q = "invalid")),
    "invalid baseline"
  )
})


test_that("plan, density, and diagnostic schemas reject malformed fields", {

  plan_fields <- list(
    target             = list(target_key = "primitive|mu"),
    outputs            = list(need_density = TRUE, need_ordinate = FALSE),
    control            = list(n_points = 20L),
    row_budget         = 20L,
    method             = "q_grid_cmde",
    density_method     = "qCMDE",
    source_fingerprint = list(posterior_values = "draws"),
    prior_ordinates    = list(),
    ordinate_warnings  = character(),
    status             = "ok",
    rows               = list(),
    support            = list(),
    replacement        = list(),
    grids              = list()
  )
  expect_s3_class(.iwmde_new_plan(plan_fields), "iwmde_plan")
  invalid_plan             <- plan_fields
  invalid_plan[["status"]] <- "invalid"
  expect_error(.iwmde_new_plan(invalid_plan), "unknown status")
  missing_plan           <- plan_fields
  missing_plan[["rows"]] <- NULL
  expect_error(.iwmde_new_plan(missing_plan), "missing required")
  missing_plan <- plan_fields
  missing_plan[["prior_ordinates"]] <- NULL
  expect_error(.iwmde_new_plan(missing_plan), "missing required")

  density_fields <- list(
    x                 = c(0, 1),
    y                 = c(.4, .3),
    finite_terms      = c(100L, 100L),
    max_log_ratio     = c(0, 0),
    ess               = c(80, 75),
    max_weight_share  = c(.02, .03),
    mcse              = c(.01, .01),
    relative_mcse     = c(.025, .033),
    active_branch_mcse = c(.01, .01),
    active_branch_relative_mcse = c(.025, .033),
    active_mass_mcse = 0,
    active_mass_relative_mcse = 0,
    active_mass_component_mcse = c(0, 0),
    mixture_mcse_type = "selected_continuous_rows_batch_means",
    sampling_mcse          = c(.005, .005),
    sampling_relative_mcse = c(.0125, .0165),
    sampling_fraction      = .5,
    sampling_uncertainty_type = "finite_population_srswor",
    mcmc_uncertainty_scope = "selected_continuous_rows_only",
    mcmc_uncertainty_status = "available",
    mcmc_uncertainty_reason = NULL,
    n_candidate_rows  = 100L,
    n_evaluated_rows  = 100L,
    normalization_points = 50L,
    normalization_range  = c(-3, 3),
    normalization_relative_error = 0,
    normalization_scale = "identity",
    normalization_mass_ratio = 1,
    max_normalizer_relative_change = 0,
    max_quadrature_relative_change = 0,
    median_normalizer_relative_change = 0,
    normalization_refined_points = 50L,
    normalization_refined_range  = c(-3, 3),
    integral_mcse          = .01,
    integral_relative_mcse = .01,
    batch_size              = 10L,
    n_batches               = 10L,
    estimator         = "q_grid_cmde",
    weight_method     = "uniform",
    log_normalizer    = rep(0, 100L),
    pilot_log_normalizer = rep(0, 100L),
    normalization_initial_points = 40L,
    normalization_initial_range  = c(-2, 2),
    pilot_normalization_integral = 1,
    final_normalization_integral = 1,
    pilot_y                       = c(.4, .3),
    validation_y                  = c(.4, .3),
    ordinate_relative_change      = c(0, 0),
    ordinate_log_change           = c(0, 0),
    pilot_ordinate_relative_change = c(0, 0),
    pilot_ordinate_log_change      = c(0, 0),
    p95_normalizer_relative_change = 0,
    n_rescued_normalizer           = 0L,
    n_initial_dropped_normalizer   = 0L,
    n_refinement_steps              = 0L
  )
  expect_s3_class(
    .iwmde_new_density_result(density_fields),
    "iwmde_density_result"
  )
  invalid_density        <- density_fields
  invalid_density[["y"]] <- "invalid"
  expect_error(
    .iwmde_new_density_result(invalid_density),
    "non-numeric grid"
  )
  missing_density                 <- density_fields
  missing_density[["finite_terms"]] <- NULL
  expect_error(
    .iwmde_new_density_result(missing_density),
    "missing required"
  )
  unordered_density        <- density_fields
  unordered_density[["x"]] <- rev(unordered_density[["x"]])
  expect_error(
    .iwmde_new_density_result(unordered_density),
    "non-numeric grid"
  )
  duplicate_density        <- density_fields
  duplicate_density[["x"]] <- c(0, 0)
  expect_error(
    .iwmde_new_density_result(duplicate_density),
    "non-numeric grid"
  )
  short_density        <- density_fields
  short_density[["x"]] <- numeric()
  short_density[["y"]] <- numeric()
  expect_error(
    .iwmde_new_density_result(short_density),
    "non-numeric grid"
  )
  negative_density        <- density_fields
  negative_density[["y"]] <- c(.4, -.3)
  expect_error(
    .iwmde_new_density_result(negative_density),
    "non-numeric grid"
  )
  invalid_ess        <- density_fields
  invalid_ess[["ess"]] <- c(101, 75)
  expect_error(
    .iwmde_new_density_result(invalid_ess),
    "effective sample sizes"
  )
  invalid_mixture_type <- density_fields
  invalid_mixture_type[["mixture_mcse_type"]] <- "independent_components"
  expect_error(
    .iwmde_new_density_result(invalid_mixture_type),
    "invalid uncertainty metadata"
  )
  invalid_counts                       <- density_fields
  invalid_counts[["n_evaluated_rows"]] <- 101L
  expect_error(
    .iwmde_new_density_result(invalid_counts),
    "inconsistent row counts"
  )
  invalid_qcmde_vector              <- density_fields
  invalid_qcmde_vector[["pilot_y"]] <- .4
  expect_error(
    .iwmde_new_density_result(invalid_qcmde_vector),
    "inconsistent diagnostic vectors"
  )
  invalid_qcmde_count                         <- density_fields
  invalid_qcmde_count[["n_refinement_steps"]] <- -1L
  expect_error(
    .iwmde_new_density_result(invalid_qcmde_count),
    "invalid count"
  )
  invalid_qcmde_integral <- density_fields
  invalid_qcmde_integral[["final_normalization_integral"]] <- "invalid"
  expect_error(
    .iwmde_new_density_result(invalid_qcmde_integral),
    "invalid normalization metadata"
  )
  invalid_qcmde_range <- density_fields
  invalid_qcmde_range[["normalization_range"]] <- "invalid"
  expect_error(
    .iwmde_new_density_result(invalid_qcmde_range),
    "invalid normalization range"
  )
  invalid_qcmde_y <- density_fields
  invalid_qcmde_y[["pilot_y"]] <- c(Inf, .3)
  expect_error(
    .iwmde_new_density_result(invalid_qcmde_y),
    "invalid comparison densities"
  )
  invalid_qcmde_normalizer <- density_fields
  invalid_qcmde_normalizer[["log_normalizer"]] <- rep(NA_real_, 100L)
  expect_error(
    .iwmde_new_density_result(invalid_qcmde_normalizer),
    "non-finite final normalizer"
  )

  iwmde_fields <- density_fields
  iwmde_fields[["estimator"]]      <- "iwmde"
  iwmde_fields[["log_normalizer"]] <- numeric()
  qcmde_only <- c(
    "pilot_log_normalizer", "normalization_initial_points",
    "normalization_initial_range", "pilot_normalization_integral",
    "final_normalization_integral", "pilot_y", "validation_y",
    "ordinate_relative_change", "ordinate_log_change",
    "pilot_ordinate_relative_change", "pilot_ordinate_log_change",
    "p95_normalizer_relative_change", "n_rescued_normalizer",
    "n_initial_dropped_normalizer", "n_refinement_steps"
  )
  iwmde_fields[qcmde_only] <- NULL
  iwmde_fields[["support_grid_normalization_integral"]] <- 1
  iwmde_fields[["weight_partitions"]]       <- list()
  iwmde_fields[["n_weight_fallbacks"]]      <- 0L
  iwmde_fields[["n_weight_fallback_rows"]]  <- 0L
  iwmde_fields[["weight_fallback_from"]]    <- character()
  iwmde_fields[["weight_fallback_reasons"]] <- structure(
    integer(),
    names = character()
  )
  expect_s3_class(
    .iwmde_new_density_result(iwmde_fields),
    "iwmde_density_result"
  )
  invalid_fallback <- iwmde_fields
  invalid_fallback[["n_weight_fallbacks"]] <- 1L
  expect_error(
    .iwmde_new_density_result(invalid_fallback),
    "inconsistent fallback metadata"
  )

  diagnostic_fields <- list(
    parameter   = "mu",
    status      = "ok",
    target_key  = "primitive|mu",
    iwmde       = list(estimator = "q_grid_cmde"),
    diagnostics = list(estimator = "q_grid_cmde")
  )
  expect_s3_class(
    .iwmde_new_diagnostic(diagnostic_fields),
    "iwmde_parameter_diagnostic"
  )
  invalid_diagnostic             <- diagnostic_fields
  invalid_diagnostic[["status"]] <- "invalid"
  expect_error(
    .iwmde_new_diagnostic(invalid_diagnostic),
    "unknown status"
  )
  missing_diagnostic                 <- diagnostic_fields
  missing_diagnostic[["diagnostics"]] <- NULL
  expect_error(
    .iwmde_new_diagnostic(missing_diagnostic),
    "missing required"
  )
  expect_s3_class(
    .iwmde_unsupported("mu", "unsupported target"),
    "iwmde_parameter_diagnostic"
  )
  expect_s3_class(
    .iwmde_point_only_diagnostic(
      parameter = "mu",
      samples   = rep(0, 2),
      component = list(point_masses = c("0" = 1))
    ),
    "iwmde_parameter_diagnostic"
  )
})


test_that("stored plan specs exclude density construction state", {

  parameter_spec <- list(
    type                  = "linear",
    weights               = c(mu = 1),
    conditional           = "mu",
    conditional_rule      = "AND",
    condition_key         = "condition",
    status                = "ok",
    prior_density_context = list(large = TRUE),
    marginal_samples      = structure(1:3, prior_density = list(raw = TRUE)),
    source                = "marginal_means"
  )
  stored <- .iwmde_plan_parameter_spec(parameter_spec)

  expect_named(stored, c(
    "type", "weights", "conditional", "conditional_rule", "condition_key",
    "status"
  ))
  expect_null(stored[["prior_density_context"]])
  expect_null(stored[["marginal_samples"]])
})


test_that("IWMDE planning propagates row-state programming errors", {

  testthat::local_mocked_bindings(
    .iwmde_row_states = function(...) {
      stop("row-state programming defect", call. = FALSE)
    },
    .package = "RoBMA"
  )

  expect_error(
    .iwmde_plan_baseline_contract(
      context          = list(),
      plan             = list(
        target         = list(parameter = "mu"),
        execution_spec = list(type = "primitive")
      ),
      candidate_rows   = 1:20,
      candidate_values = rep(0, 20)
    ),
    "row-state programming defect"
  )
})


test_that("plan keys carry schema and algorithm versions", {

  prior_ordinates <- .iwmde_prior_ordinate_classifications(
    BayesTools::prior("normal", parameters = list(mean = 0, sd = 1)),
    0
  )
  plan <- list(
    target             = list(target_key = "primitive|mu"),
    outputs            = list(need_density = TRUE),
    control            = list(n_points = 20L),
    method             = "q_grid_cmde",
    density_method     = "qCMDE",
    status             = "ok",
    prior_ordinates    = prior_ordinates,
    source_fingerprint = list(posterior_values = "draws")
  )
  payload <- .iwmde_plan_key_payload(plan)

  expect_equal(payload[["schema_version"]], "6")
  expect_equal(payload[["algorithm_version"]], "14")
  expect_identical(
    payload[["prior_ordinates"]],
    .iwmde_compact_nulls(prior_ordinates)
  )

  testthat::local_mocked_bindings(
    .iwmde_algorithm_version = function() "changed",
    .package                 = "RoBMA"
  )
  changed_payload <- .iwmde_plan_key_payload(plan)

  expect_equal(changed_payload[["algorithm_version"]], "changed")
  expect_false(identical(
    .iwmde_hash("iwmde_plan", payload),
    .iwmde_hash("iwmde_plan", changed_payload)
  ))
})
