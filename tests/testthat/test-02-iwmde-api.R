# ============================================================================ #
# test-02-iwmde-api.R
# ============================================================================ #

context("IWMDE public API")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


test_that("density methods are case-insensitive", {

  expect_equal(.density_method_normalize("kde"), "KDE")
  expect_equal(.density_method_normalize("QcMdE"), "qCMDE")
  expect_equal(.density_method_normalize("iwmde"), "IWMDE")
  expect_equal(.density_method_normalize("NORMAL", allow_normal = TRUE), "normal")
  expect_equal(.density_method_normalize_precomputed("qcmde"), "qCMDE")
  expect_equal(.density_method_normalize_precomputed("Iwmde"), "IWMDE")
})

test_that("IWMDE histogram uses global active-sample mass", {

  values <- c(-2, -1, -.5, 0, .5, 1, 2)
  hist_data <- .iwmde_histogram(values, xlim = c(-1, 1), mass = .8)
  hist_mass <- sum(hist_data[["density"]] * diff(hist_data[["breaks"]]))

  expect_equal(hist_mass, .8 * 5 / 7, tolerance = 1e-12)

  empty <- .iwmde_histogram(values, xlim = c(3, 4), mass = .8)
  expect_equal(sum(empty[["density"]] * diff(empty[["breaks"]])), 0)
  expect_true(all(empty[["density"]] == 0))

  all_missing <- .iwmde_histogram(c(NA_real_, NaN), xlim = c(0, 1), mass = .8)
  expect_equal(sum(all_missing[["density"]] * diff(all_missing[["breaks"]])), 0)
  expect_true(all(all_missing[["density"]] == 0))
})


test_that("plot qCMDE/IWMDE matrix samples use exact parameter columns", {

  sample <- matrix(
    c(1, 2, 3, 4),
    nrow     = 2,
    dimnames = list(NULL, c("theta_extra", "theta"))
  )
  samples <- list(theta = sample)

  expect_equal(
    .plot_brma_plotted_samples(samples, "theta", "theta"),
    c(3, 4)
  )

  sample_missing <- sample[, 1L, drop = FALSE]
  colnames(sample_missing) <- "theta_extra"
  samples <- list(theta = sample_missing)

  expect_null(.plot_brma_plotted_samples(samples, "theta", "theta"))
})


test_that("fixed nonzero tau priors fill missing posterior tau", {

  priors <- list(
    outcome = list(
      tau = BayesTools::prior("spike", list(0.05))
    )
  )
  samples <- matrix(
    seq_len(4) / 10,
    ncol     = 1L,
    dimnames = list(NULL, "mu")
  )

  expect_equal(.fixed_tau_prior_value(priors), 0.05)
  expect_false(.has_fixed_zero_tau_prior(priors))

  expect_error(
    .evaluate.brma.tau(
      fit               = NULL,
      scale_data        = NULL,
      scale_formula     = NULL,
      scale_priors      = NULL,
      is_scale          = FALSE,
      is_multilevel     = FALSE,
      K                 = 3L,
      posterior_samples = samples
    ),
    "Missing posterior tau columns"
  )

  tau_result <- .evaluate.brma.tau(
    fit               = NULL,
    scale_data        = NULL,
    scale_formula     = NULL,
    scale_priors      = NULL,
    is_scale          = FALSE,
    is_multilevel     = FALSE,
    K                 = 3L,
    posterior_samples = samples,
    allow_missing_tau = .fixed_tau_prior_value(priors)
  )
  expect_equal(tau_result[["tau_total"]], matrix(0.05, nrow = 4L, ncol = 3L))

  active_setup <- list(
    priors            = priors,
    is_PET            = FALSE,
    is_PEESE          = FALSE,
    is_weightfunction = FALSE
  )
  localized <- .iwmde_likelihood_posterior_samples(
    context      = list(indicator_names = character(), selection_spec = NULL),
    samples      = samples,
    active_setup = active_setup
  )
  expect_equal(localized[, "tau"], rep(0.05, nrow(samples)))

  data <- list(
    outcome = data.frame(yi = seq_len(3)),
    measure = "GEN"
  )
  attr(data, "outcome_type") <- "norm"
  row_parameters <- .iwmde_row_parameters(
    context      = list(data = data),
    row          = samples[1L, ],
    active_setup = c(active_setup, list(fit_priors = list()))
  )
  expect_equal(row_parameters[["tau"]], 0.05)
})


test_that("factor qCMDE/IWMDE densities attach only when all plotted columns resolve", {

  expect_true(.plot_brma_factor_density_complete(
    density_columns  = c("b", "a"),
    expected_columns = c("a", "b")
  ))
  expect_false(.plot_brma_factor_density_complete(
    density_columns  = "a",
    expected_columns = c("a", "b")
  ))
  expect_false(.plot_brma_factor_density_complete(
    density_columns  = character(),
    expected_columns = "a"
  ))
})


test_that("IWMDE diagnostics compute for representative cached fits", {

  fit_names <- c(
    "bcg_meta-analysis",
    "bcg_meta-regression",
    "bangertdrowns2004_location-scale",
    "konstantopoulos2011_3lvl",
    "dat.lehmann2018-PET",
    "dat.lehmann2018-PEESE",
    "dat.lehmann2018-3PSM",
    "bcg_glmm",
    "nielweise2008_glmm",
    "bcg_BMA.glmm"
  )
  .skip_if_missing_raw_fits(fit_names)

  cases <- list(
    "bcg_meta-analysis"                 = c("mu", "tau"),
    "bcg_meta-regression"               = c("mu_intercept", "mu_ablat", "tau"),
    "bangertdrowns2004_location-scale"  = c("mu_intercept", "log_tau_intercept"),
    "konstantopoulos2011_3lvl"          = c("mu", "tau", "rho", "gamma[1]"),
    "dat.lehmann2018-PET"               = c("mu", "tau", "PET"),
    "dat.lehmann2018-PEESE"             = c("mu", "tau", "PEESE"),
    "dat.lehmann2018-3PSM"              = c("mu", "tau"),
    "bcg_glmm"                          = c("mu", "tau", "pi[1]", "theta[1]"),
    "nielweise2008_glmm"                = c("mu", "tau", "phi[1]", "theta[1]"),
    "bcg_BMA.glmm"                      = c("mu", "tau", "pi[1]", "theta[1]")
  )

  for (name in names(cases)) {
    out <- plot_iwmde_diagnostics(
      object      = load_fit(name, validate = FALSE),
      parameters  = cases[[name]],
      n_points    = 20,
      max_samples = 20,
      plot        = FALSE,
      as_data     = TRUE
    )

    .expect_iwmde_ok(out, cases[[name]])
  }
})

test_that("IWMDE diagnostics can include requested density support values", {

  .skip_if_missing_raw_fits("bcg_meta-analysis")

  fit     <- load_fit("bcg_meta-analysis", validate = FALSE)
  context <- .iwmde_context(fit)
  base <- .iwmde_parameter_diagnostic(
    context              = context,
    parameter            = "mu",
    n_points             = 20,
    max_samples          = 20,
    normalization_points = 20,
    normalization_prob   = .999,
    diagnostic_cache     = .iwmde_diagnostic_cache()
  )
  include_value <- base[["xlim"]][2] + .05 * diff(base[["xlim"]])
  out <- .iwmde_parameter_diagnostic(
    context              = context,
    parameter            = "mu",
    n_points             = 20,
    max_samples          = 20,
    normalization_points = 20,
    normalization_prob   = .999,
    include_values       = include_value,
    diagnostic_cache     = .iwmde_diagnostic_cache()
  )

  expect_equal(out[["status"]], "ok")
  expect_true(include_value >= min(out[["iwmde"]][["x"]]))
  expect_true(include_value <= max(out[["iwmde"]][["x"]]))
  expect_true(any(abs(out[["iwmde"]][["x"]] - include_value) <=
                    sqrt(.Machine$double.eps) * max(1, abs(include_value))))
  expect_true(out[["diagnostics"]][["bf_included"]])
  expect_equal(out[["diagnostics"]][["bf_value"]], include_value)
  expect_true(is.finite(out[["diagnostics"]][["bf_ordinate"]]))
  expect_true(is.finite(out[["diagnostics"]][["bf_error_percent"]]))
})

test_that("IWMDE BF diagnostics retain boundary evaluation point", {

  evaluation_value <- 1e-6
  boundary_meanlog <- -.2
  boundary_sdlog   <- .5
  expected_ordinate <- stats::dlnorm(
    evaluation_value,
    meanlog = boundary_meanlog,
    sdlog   = boundary_sdlog
  )
  density <- list(
    x                = 0,
    evaluation_x     = evaluation_value,
    y                = expected_ordinate,
    mcse             = .1,
    relative_mcse    = .05,
    finite_terms     = 40L,
    ess              = 25,
    max_weight_share = .2,
    max_log_ratio    = 1,
    pilot_y          = 2,
    ordinate_relative_change = 0,
    ordinate_log_change      = 0
  )

  bf_diagnostics <- .iwmde_density_bf_diagnostics(density, 0)
  diagnostic <- list(
    status      = "ok",
    diagnostics = c(
      bf_diagnostics,
      list(
        estimator              = "q_grid_cmde",
        weight_method          = "conditional_grid",
        active_mass            = 1,
        normalization_integral = 1,
        normalization_final_integral = 1,
        normalization_mass_ratio = 1,
        max_ordinate_relative_change = 0,
        max_normalizer_relative_change = 0,
        normalization_range    = c(evaluation_value / 2, evaluation_value * 2)
      )
    )
  )
  ordinate <- .iwmde_posterior_ordinate_attribute(diagnostic, "qCMDE")

  expect_true(bf_diagnostics[["bf_included"]])
  expect_equal(bf_diagnostics[["bf_value"]], 0)
  expect_equal(bf_diagnostics[["bf_evaluation_value"]], evaluation_value)
  expect_equal(ordinate[["value"]], 0)
  expect_equal(ordinate[["evaluation_value"]], evaluation_value)
  expect_equal(ordinate[["ordinate"]], expected_ordinate)
  expect_equal(
    ordinate[["iwmde_provenance"]][["result"]][["evaluation_value"]],
    evaluation_value
  )
  expect_true(.iwmde_posterior_ordinate_supports_bf(ordinate))
})


test_that("qCMDE/IWMDE posterior attributes carry RoBMA provenance", {

  density_control <- list(
    n_points             = 20,
    max_samples          = 30,
    normalization_points = 40,
    normalization_prob   = .95,
    display_grid         = "adaptive"
  )
  metadata <- list(parameter = "mu", conditional_rule = "AND")
  diagnostic <- list(
    status       = "ok",
    target_key   = "primitive|mu",
    active_mass  = 1,
    point_masses = data.frame(x = numeric(), mass = numeric()),
    iwmde        = list(
      x                       = c(-1, 0, 1),
      y                       = c(.1, .4, .1),
      estimator               = "q_grid_cmde",
      normalization_integral  = 1,
      normalization_mass_ratio = 1,
      ordinate_relative_change = c(0, 0, 0)
    ),
    diagnostics = list(
      bf_included         = TRUE,
      bf_value            = 0,
      bf_evaluation_value = 1e-6,
      bf_ordinate         = .4,
      bf_mcse             = .01,
      bf_relative_mcse    = .01,
      bf_error_percent    = 1,
      bf_finite_terms     = 40,
      bf_ess              = 30,
      bf_max_weight_share = .2,
      bf_max_log_ratio    = 0,
      active_mass         = 1,
      normalization_integral = 1,
      bf_ordinate_relative_change = 0,
      max_ordinate_relative_change = 0,
      max_normalizer_relative_change = 0,
      normalization_range = c(-1, 1),
      estimator           = "q_grid_cmde",
      weight_method       = "mock"
    )
  )

  density_attr <- .iwmde_posterior_density_attribute(
    diagnostic      = diagnostic,
    density_method  = "qCMDE",
    metadata        = metadata,
    density_control = density_control
  )
  ordinate_attr <- .iwmde_posterior_ordinate_attribute(
    diagnostic      = diagnostic,
    density_method  = "qCMDE",
    metadata        = metadata,
    density_control = density_control
  )
  expected_provenance <- .iwmde_provenance_request(
    density_method   = "qCMDE",
    method           = "q_grid_cmde",
    metadata         = metadata,
    density_control  = density_control,
    value            = 0,
    attribute        = "ordinate",
    target_key       = "primitive|mu"
  )
  expected_density_provenance <- .iwmde_provenance_request(
    density_method   = "qCMDE",
    method           = "q_grid_cmde",
    metadata         = metadata,
    density_control  = density_control,
    attribute        = "density",
    target_key       = "primitive|mu"
  )

  expect_true("iwmde_provenance" %in% names(density_attr))
  expect_true("iwmde_provenance" %in% names(ordinate_attr))
  expect_true(.iwmde_posterior_density_matches_request(
    posterior_density = density_attr,
    provenance        = expected_density_provenance
  ))
  no_provenance_density <- density_attr
  no_provenance_density[["iwmde_provenance"]] <- NULL
  expect_false(.iwmde_posterior_density_matches_request(
    posterior_density = no_provenance_density,
    provenance        = expected_density_provenance
  ))
  expect_false(.iwmde_posterior_density_matches_request(
    posterior_density = density_attr,
    provenance        = .iwmde_provenance_request(
      density_method   = "qCMDE",
      method           = "q_grid_cmde",
      metadata         = metadata,
      density_control  = modifyList(density_control, list(n_points = 25L)),
      attribute        = "density",
      target_key       = "primitive|mu"
    )
  ))
  expect_false(.iwmde_posterior_density_matches_request(
    posterior_density = density_attr,
    provenance        = .iwmde_provenance_request(
      density_method   = "qCMDE",
      method           = "q_grid_cmde",
      metadata         = metadata,
      density_control  = density_control,
      attribute        = "density",
      target_key       = "primitive|tau"
    )
  ))

  provenance <- ordinate_attr[["iwmde_provenance"]]
  expect_equal(provenance[["schema_version"]], "1")
  expect_equal(provenance[["provenance_level"]], "diagnostic_adapter")
  expect_equal(provenance[["density_method"]], "qCMDE")
  expect_equal(provenance[["internal_method"]], "q_grid_cmde")
  expect_equal(provenance[["target"]][["parameter"]], "mu")
  expect_equal(provenance[["target"]][["conditional_rule"]], "AND")
  expect_equal(provenance[["requested_value"]], .iwmde_key_number(0))
  expect_equal(
    provenance[["result"]][["evaluation_value"]],
    1e-6
  )
  expect_equal(
    provenance[["result"]][["evaluation_value_key"]],
    .iwmde_key_number(1e-6)
  )
  expect_true(isTRUE(provenance[["result"]][["bf_grade"]]))
  expect_true(.iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = ordinate_attr,
    value              = 0,
    provenance         = expected_provenance
  ))
  expect_false(.iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = ordinate_attr,
    value              = 0,
    provenance         = .iwmde_provenance_request(
      density_method   = "IWMDE",
      method           = "iwmde",
      metadata         = metadata,
      density_control  = density_control,
      value            = 0,
      attribute        = "ordinate",
      target_key       = "primitive|mu"
    )
  ))
})


test_that("iwmde_estimate returns plan-backed attributes and caches by provenance", {

  density_control <- list(
    n_points             = 20,
    max_samples          = 20,
    normalization_points = 20,
    normalization_prob   = .95,
    display_grid         = "adaptive"
  )
  context <- list(
    posterior_samples = matrix(
      stats::rnorm(60),
      ncol     = 1L,
      dimnames = list(NULL, "mu")
    ),
    object = list(
      fit        = list(),
      likelihood = list(family = "normal"),
      data       = list(measure = "yi")
    )
  )
  mock_diagnostic <- function(parameter, target_key) {

    list(
      parameter    = parameter,
      status       = "ok",
      target_key   = target_key,
      active_mass  = 1,
      point_masses = data.frame(x = numeric(), mass = numeric()),
      iwmde        = list(
        x                         = c(-1, 0, 1),
        y                         = c(.1, .4, .1),
        estimator                 = "q_grid_cmde",
        normalization_integral    = 1,
        normalization_mass_ratio  = 1,
        ordinate_relative_change  = c(0, 0, 0)
      ),
      diagnostics = list(
        bf_included         = TRUE,
        bf_value            = 0,
        bf_evaluation_value = 0,
        bf_ordinate         = .4,
        bf_mcse             = .01,
        bf_relative_mcse    = .01,
        bf_error_percent    = 1,
        bf_finite_terms     = 40,
        bf_ess              = 30,
        bf_max_weight_share = .2,
        bf_max_log_ratio    = 0,
        active_mass         = 1,
        normalization_integral = 1,
        bf_ordinate_relative_change = 0,
        max_ordinate_relative_change = 0,
        max_normalizer_relative_change = 0,
        normalization_range = c(-1, 1),
        estimator           = "q_grid_cmde",
        weight_method       = "mock"
      )
    )
  }
  density_calls  <- 0L
  ordinate_calls <- 0L

  testthat::local_mocked_bindings(
    .iwmde_execute_plan_diagnostic = function(context, plan, output,
                                              execution_cache = NULL,
                                              diagnostic_cache = NULL) {

      if (identical(output, "density")) {
        density_calls <<- density_calls + 1L
      } else {
        ordinate_calls <<- ordinate_calls + 1L
      }
      mock_diagnostic(
        parameter  = plan[["target"]][["parameter"]],
        target_key = plan[["target"]][["target_key"]]
      )
    },
    .package = "RoBMA"
  )

  cache <- .iwmde_estimate_cache()
  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = density_control,
    outputs         = c("density", "ordinate"),
    values          = 0,
    parameter_spec  = list(type = "primitive"),
    metadata        = list(parameter = "mu"),
    cache           = cache
  )
  estimate_again <- .iwmde_estimate(
    context         = context,
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = density_control,
    outputs         = c("density", "ordinate"),
    values          = 0,
    parameter_spec  = list(type = "primitive"),
    metadata        = list(parameter = "mu"),
    cache           = cache
  )

  expect_s3_class(estimate, "iwmde_estimate")
  expect_s3_class(estimate[["plan"]], "iwmde_plan")
  expect_equal(estimate_again[["plan"]][["plan_key"]], estimate[["plan"]][["plan_key"]])
  expect_equal(density_calls, 1L)
  expect_equal(ordinate_calls, 1L)
  expect_true("iwmde_provenance" %in% names(estimate[["posterior_density"]]))
  expect_true("iwmde_provenance" %in% names(estimate[["posterior_ordinate"]]))

  provenance <- estimate[["posterior_ordinate"]][["iwmde_provenance"]]
  expect_equal(provenance[["provenance_level"]], "iwmde_plan")
  expect_equal(provenance[["plan_key"]], estimate[["plan"]][["plan_key"]])
  expect_true(is.list(provenance[["source_fingerprint"]]))
  expect_true(.iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = estimate[["posterior_ordinate"]],
    value              = 0,
    provenance         = .iwmde_estimate_request_provenance(
      plan      = estimate[["plan"]],
      attribute = "ordinate",
      value     = 0
    )
  ))

  context_changed <- context
  context_changed[["posterior_samples"]][1, 1] <- 99
  changed_plan <- .iwmde_plan(
    context         = context_changed,
    parameter       = "mu",
    density_method  = "qCMDE",
    density_control = density_control,
    outputs         = "ordinate",
    values          = 0,
    parameter_spec  = list(type = "primitive"),
    metadata        = list(parameter = "mu")
  )
  expect_false(.iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = estimate[["posterior_ordinate"]],
    value              = 0,
    provenance         = .iwmde_estimate_request_provenance(
      plan      = changed_plan,
      attribute = "ordinate",
      value     = 0
    )
  ))
})


test_that("qCMDE ordinate and IWMDE mass thresholds warn before failing", {

  ordinate <- function(estimator, normalization_integral,
                       ordinate_relative_change = NA_real_) {

    BayesTools::posterior_ordinate_attribute(
      value          = 0,
      ordinate       = .4,
      method         = estimator,
      density_method = if (identical(estimator, "q_grid_cmde")) {
        "qCMDE"
      } else {
        "IWMDE"
      },
      diagnostics    = list(
        evaluation_value       = 0,
        relative_mcse          = .1,
        finite_terms           = 60,
        ess                    = 30,
        max_weight_share       = .2,
        active_mass            = 1,
        normalization_integral = normalization_integral,
        ordinate_relative_change = ordinate_relative_change,
        max_normalizer_relative_change = ordinate_relative_change,
        normalization_range    = c(-1, 1),
        estimator              = estimator,
        weight_method          = "test"
      )
    )
  }

  qcmde_ok <- ordinate("q_grid_cmde", 1, .02)
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_ok))
  expect_length(.iwmde_posterior_ordinate_warnings(qcmde_ok), 0L)

  qcmde_warn <- ordinate("q_grid_cmde", 1, .03)
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_warn))
  expect_true(any(grepl(
    "qCMDE.*ordinate.*3%.*2.5%.*5%",
    .iwmde_posterior_ordinate_warnings(qcmde_warn)
  )))

  qcmde_fail <- ordinate("q_grid_cmde", 1, .06)
  expect_false(.iwmde_posterior_ordinate_supports_bf(qcmde_fail))
  expect_match(
    .iwmde_posterior_ordinate_failure_reasons(qcmde_fail),
    "qCMDE.*ordinate.*6%.*5%"
  )

  iwmde_warn <- ordinate("iwmde", .94)
  expect_true(.iwmde_posterior_ordinate_supports_bf(iwmde_warn))
  expect_true(any(grepl(
    "IWMDE.*6%.*5%.*10%",
    .iwmde_posterior_ordinate_warnings(iwmde_warn)
  )))

  iwmde_fail <- ordinate("iwmde", .89)
  expect_false(.iwmde_posterior_ordinate_supports_bf(iwmde_fail))
  expect_match(
    .iwmde_posterior_ordinate_failure_reasons(iwmde_fail),
    "IWMDE.*11%.*10%"
  )

  bf <- .iwmde_bf_append_warning(1, iwmde_warn)
  expect_true(any(grepl("IWMDE.*6%", attr(bf, "warnings"))))

  posterior <- structure(
    stats::rnorm(20),
    posterior_ordinate = qcmde_warn
  )
  table <- data.frame(BF = "1", row.names = "mu")
  table <- .hypothesis_brma_append_iwmde_warnings(
    table     = table,
    posterior = posterior
  )
  expect_true(any(grepl("qCMDE.*ordinate.*3%", attr(table, "warnings"))))

  multi_ordinate <- list(
    status    = "ok",
    ordinates = list(
      list(
        value       = 0,
        ordinate    = 1,
        diagnostics = list(warning = "warning for zero"),
        parameter   = "theta"
      ),
      list(
        value       = 1,
        ordinate    = 1,
        diagnostics = list(warning = "warning for one"),
        parameter   = "theta"
      )
    )
  )
  posterior <- structure(
    stats::rnorm(20),
    posterior_ordinate = multi_ordinate
  )
  table <- data.frame(
    Alternative = c("theta != 0", "theta != 1"),
    Null        = c("theta = 0", "theta = 1"),
    BF          = c("1", "1"),
    row.names   = c("theta", "theta1"),
    check.names = FALSE
  )
  attr(table, "parsed") <- list(
    list(
      left  = list(type = "point", expr = "theta", value = 0, label = "theta = 0"),
      right = list(type = "not_point", expr = "theta", value = 0, label = "theta != 0")
    ),
    list(
      left  = list(type = "point", expr = "theta", value = 1, label = "theta = 1"),
      right = list(type = "not_point", expr = "theta", value = 1, label = "theta != 1")
    )
  )
  table <- .hypothesis_brma_append_iwmde_warnings(
    table     = table,
    posterior = posterior
  )
  expect_equal(
    attr(table, "warnings"),
    c(theta = "warning for zero", theta1 = "warning for one")
  )
  compact_table <- data.frame(
    BF          = c("1", "1"),
    row.names   = c("theta", "theta1"),
    check.names = FALSE
  )
  attr(compact_table, "parsed") <- attr(table, "parsed", exact = TRUE)
  compact_table <- .hypothesis_brma_append_iwmde_warnings(
    table     = compact_table,
    posterior = posterior
  )
  expect_equal(
    attr(compact_table, "warnings"),
    c(theta = "warning for zero", theta1 = "warning for one")
  )

  factor_posterior <- list(
    a = structure(
      stats::rnorm(20),
      posterior_ordinate = list(
        value       = 0,
        ordinate    = 1,
        diagnostics = list(warning = "warning for a"),
        parameter   = "mu",
        level       = "a"
      )
    ),
    b = structure(
      stats::rnorm(20),
      posterior_ordinate = list(
        value       = 0,
        ordinate    = 1,
        diagnostics = list(warning = "warning for b"),
        parameter   = "mu",
        level       = "b"
      )
    )
  )
  table <- data.frame(
    Alternative = c("mu != 0", "mu != 0"),
    Null        = c("mu = 0", "mu = 0"),
    BF          = c("1", "1"),
    row.names   = c("mu[a]", "mu[b]"),
    check.names = FALSE
  )
  attr(table, "parsed") <- list(
    list(
      left  = list(type = "point", expr = "mu", value = 0, label = "mu = 0"),
      right = list(type = "not_point", expr = "mu", value = 0, label = "mu != 0")
    )
  )
  table <- .hypothesis_brma_append_iwmde_warnings(
    table     = table,
    posterior = factor_posterior
  )
  expect_equal(
    attr(table, "warnings"),
    c("mu[a]" = "warning for a", "mu[b]" = "warning for b")
  )

  raw_diagnostic <- list(
    status      = "ok",
    diagnostics = list(
      bf_ordinate            = .4,
      bf_relative_mcse       = .1,
      bf_finite_terms        = 60,
      bf_ess                 = 30,
      bf_max_weight_share    = .2,
      active_mass            = 1,
      normalization_integral = 1,
      bf_ordinate_relative_change = .06,
      max_normalizer_relative_change = .06,
      normalization_range    = c(-1, 1),
      estimator              = "q_grid_cmde"
    )
  )
  expect_match(
    .hypothesis_brma_diagnostic_reason(raw_diagnostic),
    "qCMDE.*ordinate.*6%.*5%"
  )
})


test_that("qCMDE and IWMDE BF warnings cover Monte Carlo reliability", {

  ordinate <- BayesTools::posterior_ordinate_attribute(
    value          = 0,
    ordinate       = .4,
    method         = "q_grid_cmde",
    density_method = "qCMDE",
    diagnostics    = list(
      evaluation_value       = 0,
      relative_mcse          = .30,
      finite_terms           = 40,
      ess                    = 10,
      max_weight_share       = .60,
      active_mass            = 1,
      normalization_integral = 1,
      ordinate_relative_change = 0,
      max_normalizer_relative_change = 0,
      normalization_range    = c(-1, 1),
      estimator              = "q_grid_cmde",
      weight_method          = "test"
    )
  )
  warnings <- .iwmde_posterior_ordinate_warnings(ordinate)

  expect_true(.iwmde_posterior_ordinate_supports_bf(ordinate))
  expect_true(any(grepl("relative MCSE.*30%.*25%.*100%", warnings)))
  expect_true(any(grepl("uses only\\s+40.*finite importance terms.*50.*20", warnings)))
  expect_true(any(grepl("effective sample size.*10.*20.*4", warnings)))
  expect_true(any(grepl("largest importance weight.*60%.*50%.*80%", warnings)))
})


test_that("qCMDE and IWMDE row-loss thresholds warn before failing", {

  ordinate <- function(estimator, n_candidate_rows, n_normalized_rows) {

    row_drop_fraction <- .iwmde_row_drop_fraction(
      n_candidate_rows  = n_candidate_rows,
      n_normalized_rows = n_normalized_rows
    )

    BayesTools::posterior_ordinate_attribute(
      value          = 0,
      ordinate       = .4,
      method         = estimator,
      density_method = if (identical(estimator, "q_grid_cmde")) {
        "qCMDE"
      } else {
        "IWMDE"
      },
      diagnostics    = list(
        evaluation_value       = 0,
        relative_mcse          = .1,
        finite_terms           = n_normalized_rows,
        ess                    = 30,
        max_weight_share       = .2,
        active_mass            = 1,
        normalization_integral = 1,
        ordinate_relative_change = 0,
        max_normalizer_relative_change = 0,
        normalization_range    = c(-1, 1),
        n_candidate_rows       = n_candidate_rows,
        n_normalized_rows      = n_normalized_rows,
        row_drop_fraction      = row_drop_fraction,
        estimator              = estimator,
        weight_method          = "test"
      )
    )
  }

  qcmde_warn <- ordinate("q_grid_cmde", 100, 97)
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_warn))
  expect_true(any(grepl(
    "qCMDE.*dropped.*3%.*2.5%.*5%",
    .iwmde_posterior_ordinate_warnings(qcmde_warn)
  )))

  qcmde_fail <- ordinate("q_grid_cmde", 100, 94)
  expect_false(.iwmde_posterior_ordinate_supports_bf(qcmde_fail))
  expect_match(
    .iwmde_posterior_ordinate_failure_reasons(qcmde_fail),
    "qCMDE.*dropped.*6%.*5%"
  )

  iwmde_warn <- ordinate("iwmde", 100, 94)
  expect_true(.iwmde_posterior_ordinate_supports_bf(iwmde_warn))
  expect_true(any(grepl(
    "IWMDE.*dropped.*6%.*5%.*10%",
    .iwmde_posterior_ordinate_warnings(iwmde_warn)
  )))

  iwmde_fail <- ordinate("iwmde", 100, 89)
  expect_false(.iwmde_posterior_ordinate_supports_bf(iwmde_fail))
  expect_match(
    .iwmde_posterior_ordinate_failure_reasons(iwmde_fail),
    "IWMDE.*dropped.*11%.*10%"
  )
})


test_that("IWMDE diagnostics warn on unused dots", {

  .skip_if_missing_raw_fits("bcg_meta-analysis")

  expect_warning(
    plot_iwmde_diagnostics(
      object      = load_fit("bcg_meta-analysis", validate = FALSE),
      parameters  = "mu",
      n_points    = 20,
      max_samples = 20,
      method      = "iwmde",
      plot        = FALSE
    ),
    "Unused argument.*method"
  )
})

test_that("old IWMDE variant labels are rejected", {

  .skip_if_missing_raw_fits(c("bcg_meta-analysis", "bcg_meta-regression2"))

  expect_error(
    plot_iwmde_diagnostics(
      object         = load_fit("bcg_meta-analysis", validate = FALSE),
      parameters     = "mu",
      n_points       = 20,
      max_samples    = 20,
      density_method = "IWMDE-Chen",
      plot           = FALSE
    ),
    "must be one of"
  )
  expect_error(
    plot(
      load_fit("bcg_meta-analysis", validate = FALSE),
      parameter      = "mu",
      density_method = "IWMDE-CMDE"
    ),
    "must be one of"
  )
  expect_error(
    marginal_means(
      load_fit("bcg_meta-regression2", validate = FALSE),
      n_samples      = 1000,
      density_method = "IWMDE-Chen"
    ),
    "must be one of"
  )
})

test_that("unused dots warn in cleaned density interfaces", {

  .skip_if_missing_raw_fits(c("bcg_meta-analysis", "bcg_meta-regression2"))

  expect_warning(
    plot(
      load_fit("bcg_meta-analysis", validate = FALSE),
      parameter       = "mu",
      density_method  = "qCMDE",
      density_control = list(n_points = 20, max_samples = 20),
      iwmde_n_points  = 20
    ),
    "Unused argument.*iwmde_n_points"
  )
  expect_warning(
    marginal_means(
      load_fit("bcg_meta-regression2", validate = FALSE),
      n_samples          = 1000,
      density_method     = "qCMDE",
      density_control    = list(n_points = 20, max_samples = 20),
      iwmde_max_samples  = 20
    ),
    "Unused argument.*iwmde_max_samples"
  )
})

test_that("density_control validates public density settings", {

  valid <- .density_control_normalize(
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 20,
      max_samples          = 30,
      normalization_points = 40,
      normalization_prob   = .95,
      display_grid         = "uniform"
    )
  )

  expect_equal(valid[["n_points"]], 20)
  expect_equal(valid[["max_samples"]], 30)
  expect_equal(valid[["normalization_points"]], 40)
  expect_equal(valid[["normalization_prob"]], .95)
  expect_equal(valid[["display_grid"]], "uniform")
  expect_error(
    .density_control_normalize("qCMDE", list(unknown = 1)),
    "unrecognized"
  )
  expect_error(
    .density_control_normalize("qCMDE", list(20)),
    "fully named"
  )
  expect_error(
    .density_control_normalize(
      "qCMDE",
      structure(list(20, 30), names = c("n_points", "n_points"))
    ),
    "duplicate"
  )
  expect_error(
    .density_control_normalize("KDE", list(n_points = 20)),
    "only used"
  )
  valid_iwmde <- .density_control_normalize(
    density_method  = "IWMDE",
    density_control = list(
      normalization_points = 40,
      normalization_prob   = .95
    )
  )
  expect_equal(valid_iwmde[["normalization_points"]], 40)
  expect_equal(valid_iwmde[["normalization_prob"]], .95)
  mm_iwmde <- .hypothesis_marginal_means_density_control(
    object = list(iwmde_settings = list(
      n_points             = 20,
      max_samples          = 30,
      normalization_points = 40,
      normalization_prob   = .95,
      display_grid         = "uniform"
    )),
    density_method  = "IWMDE",
    density_control = NULL
  )
  expect_equal(mm_iwmde[["normalization_points"]], 40)
  expect_equal(mm_iwmde[["normalization_prob"]], .95)
  expect_error(
    .density_control_normalize("qCMDE", list(normalization_prob = 0)),
    "higher than 0"
  )
  expect_error(
    .density_control_normalize("qCMDE", list(normalization_prob = NA_real_)),
    "cannot contain"
  )
})

test_that("IWMDE row thinning is deterministic and equally spaced", {

  rows <- seq(2L, 80L, by = 2L)
  set.seed(123)
  old_seed <- .Random.seed

  out <- .iwmde_select_active_rows(rows = rows, max_samples = 10)

  expect_equal(out, rows[unique(round(seq(1, length(rows), length.out = 10)))])
  expect_identical(.Random.seed, old_seed)
  expect_equal(.iwmde_select_active_rows(rows = rows, max_samples = 10), out)
})


test_that("IWMDE plan freezes finite baseline row contract", {

  context <- list(
    posterior_samples = matrix(
      seq_len(25),
      ncol     = 1L,
      dimnames = list(NULL, "mu")
    ),
    object = list(
      fit        = list(),
      likelihood = list(family = "normal"),
      data       = list(measure = "GEN")
    ),
    data            = list(),
    priors          = list(),
    flat_prior_list = list()
  )
  drop_positions <- 1L

  testthat::local_mocked_bindings(
    .iwmde_row_states = function(context, rows, parameter = NULL,
                                 parameter_spec = NULL) {

      lapply(seq_along(rows), function(i) {
        list(baseline_log_q = if (i %in% drop_positions) -Inf else 0)
      })
    },
    .package = "RoBMA"
  )

  make_plan <- function() {
    .iwmde_plan(
      context         = context,
      parameter       = "mu",
      density_method  = "qCMDE",
      density_control = list(
        n_points             = 20,
        max_samples          = 21,
        normalization_points = 20
      ),
      outputs         = "ordinate",
      values          = 0,
      parameter_spec  = list(type = "primitive"),
      metadata        = NULL
    )
  }

  plan <- make_plan()
  expect_equal(plan[["status"]], "ok")
  expect_equal(plan[["rows"]][["n_denominator_rows"]], 21L)
  expect_equal(plan[["rows"]][["n_candidate_rows"]], 21L)
  expect_equal(plan[["rows"]][["n_estimator_rows"]], 20L)
  expect_equal(plan[["rows"]][["n_dropped_log_q"]], 1L)
  expect_equal(length(plan[["rows"]][["row_states"]]), 20L)
  expect_equal(
    plan[["rows"]][["estimator_rows"]],
    plan[["rows"]][["candidate_rows"]][-1L]
  )

  execution <- .iwmde_plan_row_execution(
    context = context,
    plan    = plan
  )
  expect_equal(execution[["active_rows"]], plan[["rows"]][["estimator_rows"]])
  expect_equal(execution[["n_candidate_rows"]], 21L)
  expect_equal(execution[["n_dropped_log_q"]], 1L)

  drop_positions <- 2L
  changed_plan <- make_plan()
  expect_false(identical(changed_plan[["plan_key"]], plan[["plan_key"]]))

  drop_positions <- 1:2
  unsupported_plan <- make_plan()
  expect_equal(unsupported_plan[["status"]], "unsupported")
  expect_match(
    unsupported_plan[["reason"]],
    "fewer than 20 finite baseline log-q"
  )
  expect_equal(unsupported_plan[["rows"]][["n_denominator_rows"]], 21L)
  expect_equal(unsupported_plan[["rows"]][["n_estimator_rows"]], 19L)
})


test_that("IWMDE condition rows keep one-row shape and ignore missing values", {

  context <- list(
    posterior_samples = data.frame(mu = 0),
    flat_prior_list   = list()
  )
  rows <- .iwmde_parameter_condition_rows(
    context = context,
    parameter_spec = list(
      conditional      = c("missing_a", NA, "missing_b"),
      conditional_rule = "OR"
    )
  )

  expect_length(rows, 1L)
  expect_false(rows)
})

test_that("IWMDE marginal-mean specs preserve child condition metadata", {

  condition_key <- paste(c("OR", 2L, "mu_alloc", "mu_intercept"), collapse = "\r")
  condition_event <- structure(
    list(
      conditional      = c("mu_intercept", "mu_alloc"),
      conditional_rule = "OR",
      condition_key    = condition_key
    ),
    class = "BayesTools_condition_event"
  )
  prior_density_context <- structure(
    list(condition_key = condition_key),
    class = "prior_density_context"
  )
  samples <- structure(
    stats::rnorm(20),
    linear_weights              = c(mu_intercept = 1, mu_alloc = 1),
    conditional                 = "stale_parent_condition",
    conditional_rule            = "AND",
    effective_conditional       = c("mu_intercept", "mu_alloc"),
    effective_conditional_rule  = "OR",
    condition_key               = condition_key,
    condition_event             = condition_event,
    resolved_condition_event    = condition_event,
    prior_density_context       = prior_density_context
  )
  marginal_means_object <- list(
    inference        = list(
      conditional = list(mu_alloc = list(alternate = samples))
    ),
    parameters       = "mu_alloc",
    term_map         = data.frame(
      term             = "alloc",
      parameter        = "mu_alloc",
      label            = "alloc",
      stringsAsFactors = FALSE
    ),
    conditional_rule = "AND"
  )

  specs <- .iwmde_marginal_means_specs(
    marginal_means_object = marginal_means_object,
    parameter             = "alloc",
    type                  = "conditional",
    levels                = "alternate"
  )
  spec <- specs[["alloc: alternate"]]

  expect_equal(spec[["conditional"]], c("mu_intercept", "mu_alloc"))
  expect_equal(spec[["conditional_rule"]], "OR")
  expect_equal(spec[["condition_key"]], condition_key)
  expect_identical(spec[["condition_event"]], condition_event)
  expect_identical(spec[["resolved_condition_event"]], condition_event)
  expect_identical(spec[["prior_density_context"]], prior_density_context)
  expect_equal(.iwmde_parameter_condition_key(spec), condition_key)
  expect_equal(
    .iwmde_posterior_metadata(samples, "mu_alloc", "alternate")[["condition_key"]],
    condition_key
  )

  stale_samples <- samples
  attr(stale_samples, "effective_conditional") <- NULL
  attr(stale_samples, "effective_conditional_rule") <- NULL
  attr(stale_samples, "condition_key") <- NULL
  marginal_means_object[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]] <- stale_samples

  expect_error(
    .iwmde_marginal_means_specs(
      marginal_means_object = marginal_means_object,
      parameter             = "alloc",
      type                  = "conditional",
      levels                = "alternate"
    ),
    "metadata is incomplete"
  )
})

test_that("IWMDE diagnostics report raw support-grid mass", {

  .skip_if_missing_raw_fits("bcg_meta-analysis")

  out <- plot_iwmde_diagnostics(
    object         = load_fit("bcg_meta-analysis", validate = FALSE),
    parameters     = c("mu", "tau"),
    n_points       = 20,
    max_samples    = 50,
    density_method = "IWMDE",
    plot           = FALSE,
    as_data        = TRUE
  )

  .expect_iwmde_ok(
    out            = out,
    parameters     = c("mu", "tau"),
    estimator      = "iwmde",
    mass_tolerance = .25
  )
  for (parameter in names(out)) {
    diagnostic <- out[[parameter]]
    expect_true(startsWith(diagnostic[["diagnostics"]][["weight_method"]], "chen_"))
  }
})

test_that("IWMDE restricted weights integrate without post-hoc scaling", {

  .skip_if_missing_raw_fits("bcg_meta-analysis")

  diagnostic <- .iwmde_parameter_diagnostic(
    context              = .iwmde_context(
      load_fit("bcg_meta-analysis", validate = FALSE)
    ),
    parameter            = "negative_tau",
    parameter_spec       = list(type = "linear", weights = c(tau = -1)),
    n_points             = 100,
    max_samples          = 100,
    normalization_points = 200,
    normalization_prob   = .999,
    method               = "iwmde",
    display_grid_method  = "adaptive",
    diagnostic_cache     = .iwmde_diagnostic_cache()
  )

  expect_equal(diagnostic[["status"]], "ok")
  expect_equal(diagnostic[["support"]], c(-Inf, 0))
  expect_equal(diagnostic[["diagnostics"]][["weight_method"]], "chen_gamma")
  expect_true(is.finite(diagnostic[["iwmde"]][["normalization_integral"]]))
  expect_gt(diagnostic[["iwmde"]][["normalization_mass_ratio"]], 0)
  expect_equal(diagnostic[["diagnostics"]][["normalization_scale"]], "support_grid")
  expect_equal(
    diagnostic[["diagnostics"]][["normalization_integral"]],
    diagnostic[["active_mass"]],
    tolerance = .10
  )
  expect_equal(
    diagnostic[["diagnostics"]][["integral"]],
    diagnostic[["active_mass"]],
    tolerance = .05
  )
})

test_that("IWMDE uses moment-matched recommended weights", {

  .skip_if_missing_raw_fits(c("bcg_meta-analysis", "konstantopoulos2011_3lvl"))

  meta_out <- plot_iwmde_diagnostics(
    object         = load_fit("bcg_meta-analysis", validate = FALSE),
    parameters     = c("mu", "tau"),
    n_points       = 20,
    max_samples    = 50,
    density_method = "IWMDE",
    plot           = FALSE,
    as_data        = TRUE
  )
  rho_out <- plot_iwmde_diagnostics(
    object         = load_fit("konstantopoulos2011_3lvl", validate = FALSE),
    parameters     = "rho",
    n_points       = 20,
    max_samples    = 30,
    density_method = "IWMDE",
    plot           = FALSE,
    as_data        = TRUE
  )

  expect_equal(meta_out[["mu"]][["diagnostics"]][["weight_method"]],
               "chen_conditional_normal")
  expect_equal(meta_out[["tau"]][["diagnostics"]][["weight_method"]],
               "chen_gamma")
  expect_equal(rho_out[["rho"]][["diagnostics"]][["weight_method"]],
               "chen_logit_conditional_normal")
})

test_that("restricted Chen fallback weights are proper densities", {

  gamma_weight <- .iwmde_chen_gamma_log_weight_single(
    active_values = seq(.001, 3, length.out = 500),
    weight_values = c(.25, .30, .35, .40, .45, .50),
    support       = c(0, Inf)
  )
  gamma_mass <- .iwmde_trapz(
    seq(.001, 3, length.out = 500),
    exp(gamma_weight[["log_weight"]])
  )

  beta_grid <- seq(.001, .999, length.out = 500)
  beta_weight <- .iwmde_chen_beta_log_weight(
    active_values = beta_grid,
    weight_values = c(.25, .30, .35, .40, .45, .50),
    support       = c(0, 1)
  )
  beta_mass <- .iwmde_trapz(beta_grid, exp(beta_weight[["log_weight"]]))

  expect_equal(gamma_weight[["method"]], "chen_gamma")
  expect_equal(beta_weight[["method"]], "chen_beta")
  expect_equal(gamma_mass, 1, tolerance = .01)
  expect_equal(beta_mass, 1, tolerance = .01)
})

test_that("IWMDE Chen conditional-normal weights error on failed conditioning", {

  context <- list(
    posterior_samples = data.frame(
      mu  = c(-0.2, -0.1, 0.1, 0.2),
      tau = c(0.3, 0.4, 0.5, 0.6)
    ),
    indicator_names = character(),
    selection_spec  = NULL
  )

  expect_error(
    .iwmde_chen_conditional_normal_log_weight(
      context        = context,
      parameter      = "mu",
      parameter_spec = list(type = "primitive"),
      active_rows    = 1:4,
      active_values  = context[["posterior_samples"]][["mu"]],
      weight_rows    = 1:4,
      weight_values  = c(NA_real_, NA_real_, 0.1, NA_real_)
    ),
    "fewer than three finite conditioning rows"
  )

  expect_error(
    .iwmde_chen_logit_conditional_normal_log_weight(
      context        = context,
      parameter      = "mu",
      parameter_spec = list(type = "primitive"),
      active_rows    = 1:4,
      active_values  = c(0.2, 0.3, 0.4, 0.5),
      weight_rows    = 1:4,
      weight_values  = c(NA_real_, NA_real_, 0.4, NA_real_),
      support        = c(0, 1)
    ),
    "fewer than three finite conditioning rows"
  )
})

test_that("IWMDE diagnostics compute for estimated marginal means", {

  .skip_if_missing_raw_fits("bcg_meta-regression2")

  fit <- load_fit("bcg_meta-regression2", validate = FALSE)
  mm  <- marginal_means(fit, n_samples = 1000)
  out <- plot_iwmde_marginal_means_diagnostics(
    object                = fit,
    parameter             = "alloc",
    marginal_means_object = mm,
    n_points              = 20,
    max_samples           = 20,
    plot                  = FALSE,
    as_data               = TRUE
  )

  .expect_iwmde_ok(out, names(out))
  expect_named(out, c("alloc: alternate", "alloc: random", "alloc: systematic"))
  expect_s3_class(mm[["source_object"]], "brma")
})

test_that("IWMDE diagnostics compute for conditional estimated marginal means", {

  .skip_if_missing_raw_fits("dat.lehmann2018_RoBMA_mods")

  fit <- load_fit("dat.lehmann2018_RoBMA_mods", validate = FALSE)
  mm  <- marginal_means(fit, n_samples = 1000)
  out <- plot_iwmde_marginal_means_diagnostics(
    object                = fit,
    parameter             = "Preregistered",
    type                  = "conditional",
    levels                = "Not Pre-Registered",
    marginal_means_object = mm,
    n_points              = 20,
    max_samples           = 20,
    plot                  = FALSE,
    as_data               = TRUE
  )

  .expect_iwmde_ok(out, names(out))
  expect_named(out, "Preregistered: Not Pre-Registered")
  expect_true(out[["Preregistered: Not Pre-Registered"]][["diagnostics"]][["n_total"]] <
                length(mm[["inference"]][["averaged"]][["mu_Preregistered"]][["Not Pre-Registered"]]))
})

test_that("IWMDE reuses identical marginal-mean linear targets", {

  .skip_if_missing_raw_fits("bcg_meta-regression4")

  fit <- load_fit("bcg_meta-regression4", validate = FALSE)
  mm  <- marginal_means(fit, n_samples = 1000)
  out <- plot_iwmde_marginal_means_diagnostics(
    object                = fit,
    parameter             = c("intercept", "alloc", "year_before1969", "alloc:year_before1969"),
    levels                = c("intercept", "alternate", "FALSE", "alternate, FALSE"),
    marginal_means_object = mm,
    n_points              = 20,
    max_samples           = 20,
    plot                  = FALSE,
    as_data               = TRUE
  )

  duplicate_targets <- c(
    "intercept: intercept",
    "alloc: alternate",
    "year_before1969: FALSE",
    "alloc:year_before1969: alternate, FALSE"
  )
  .expect_iwmde_ok(out[duplicate_targets], duplicate_targets)

  baseline <- out[[duplicate_targets[1L]]]
  for (target in duplicate_targets[-1L]) {
    expect_equal(out[[target]][["target_key"]], baseline[["target_key"]])
    expect_equal(out[[target]][["active_rows"]], baseline[["active_rows"]])
    expect_equal(out[[target]][["iwmde"]][["x"]], baseline[["iwmde"]][["x"]])
    expect_equal(out[[target]][["iwmde"]][["y"]], baseline[["iwmde"]][["y"]])
    expect_equal(out[[target]][["parameter"]], target)
  }
})

test_that("IWMDE marginal means use attached source objects", {

  .skip_if_missing_raw_fits(c("bcg_meta-regression2", "bcg_meta-regression2b"))

  fit    <- load_fit("bcg_meta-regression2", validate = FALSE)
  fit_b  <- load_fit("bcg_meta-regression2b", validate = FALSE)
  mm     <- marginal_means(fit, n_samples = 1000)
  mm_old <- mm
  mm_old[["source_object"]] <- NULL

  out <- plot_iwmde_marginal_means_diagnostics(
    object                = fit_b,
    parameter             = "alloc",
    marginal_means_object = mm,
    n_points              = 20,
    max_samples           = 20,
    plot                  = FALSE,
    as_data               = TRUE
  )
  .expect_iwmde_ok(out, names(out))

  expect_error(
    plot_iwmde_marginal_means_diagnostics(
      object                = fit,
      marginal_means_object = mm_old,
      plot                  = FALSE,
      as_data               = TRUE
    ),
    "source fitted brma object"
  )
  expect_error(
    plot_iwmde_marginal_means_diagnostics(
      object                = fit,
      marginal_means_object = marginal_means(fit, n_samples = 1000, transform = "EXP"),
      plot                  = FALSE,
      as_data               = TRUE
    ),
    "linear predictor scale"
  )
})

test_that("IWMDE diagnostics handle RoBMA mixture branches", {

  .skip_if_missing_raw_fits("dat.lehmann2018_RoBMA")

  out <- plot_iwmde_diagnostics(
    object      = load_fit("dat.lehmann2018_RoBMA", validate = FALSE),
    parameters  = c("mu", "tau", "PET", "PEESE", "omega[1]", "eta[1]"),
    n_points    = 20,
    max_samples = 20,
    plot        = FALSE,
    as_data     = TRUE
  )

  .expect_iwmde_ok(out[c("mu", "tau", "PET", "PEESE")], c("mu", "tau", "PET", "PEESE"))
  expect_equal(out[["omega[1]"]][["status"]], "unsupported")
  expect_equal(out[["eta[1]"]][["status"]], "unsupported")
  expect_gt(nrow(out[["PET"]][["point_masses"]]), 0)
  expect_gt(nrow(out[["PEESE"]][["point_masses"]]), 0)
  expect_equal(
    sum(out[["PET"]][["histogram"]][["density"]] * diff(out[["PET"]][["histogram"]][["breaks"]])),
    out[["PET"]][["active_mass"]],
    tolerance = 1e-8
  )
})

test_that("IWMDE diagnostics draw base plot matrix", {

  .skip_if_missing_raw_fits("dat.lehmann2018_RoBMA")

  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = file, width = 8, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_iwmde_diagnostics(
    object      = load_fit("dat.lehmann2018_RoBMA", validate = FALSE),
    parameters  = c("mu", "tau", "PET", "PEESE"),
    n_points    = 20,
    max_samples = 20
  ))
})
