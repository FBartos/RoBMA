context("Hypothesis Bayes factors")

source(testthat::test_path("common-functions.R"))


.mock_random_non_known_v_brma_mv <- function() {

  data <- structure(list(), random = TRUE)
  structure(
    list(
      fit    = list(dummy = TRUE),
      data   = data,
      priors = list()
    ),
    class = c("brma.mv", "brma")
  )
}


.hypothesis_expect_bridge_ready <- function(object) {

  mcse <- .hypothesis_bridge_mcse(object)
  expect_true(is.finite(logml(object)))
  expect_true(is.finite(mcse))
  expect_lt(mcse, 0.05)
}


.hypothesis_bridge_mcse <- function(object) {

  repetitions <- object[["marglik"]][["repetitions"]]
  included    <- repetitions[["success"]] & repetitions[["finite"]]
  if (!any(included)) {
    return(NA_real_)
  }

  return(max(repetitions[["mcse"]][included]))
}


.hypothesis_bridge_bf01 <- function(null, full) {

  exp(logml(null) - logml(full))
}


.hypothesis_expect_bridge_agreement <- function(result, null, full,
                                                hard_tolerance          = 0.05,
                                                uncertainty_multiplier = 3,
                                                index                   = 1L) {

  raw_bf           <- as.numeric(attr(result, "raw_BF", exact = TRUE))[index]
  bf_error_percent <- as.numeric(result[["BF_error"]])[index]
  expect_length(raw_bf, 1L)
  expect_length(bf_error_percent, 1L)

  log_difference     <- abs(
    log(raw_bf) -
      (logml(full) - logml(null))
  )
  estimator_log_mcse <- sqrt(log1p((bf_error_percent / 100)^2))
  bridge_log_mcse    <- sqrt(
    .hypothesis_bridge_mcse(null)^2 +
      .hypothesis_bridge_mcse(full)^2
  )
  combined_log_mcse  <- sqrt(estimator_log_mcse^2 + bridge_log_mcse^2)

  expect_true(
    log_difference <= hard_tolerance,
    info = "estimator and bridge log Bayes factors must meet the hard tolerance"
  )
  expect_true(
    log_difference <= uncertainty_multiplier * combined_log_mcse,
    info = "estimator and bridge log Bayes factors must agree within uncertainty"
  )
}


test_that("qCMDE factor point guards use display aliases", {

  expect_equal(
    .hypothesis_brma_alias_label(
      aliases   = list(mu_alloc = "mu_alloc", alloc = "mu_alloc"),
      parameter = "mu_alloc"
    ),
    "alloc"
  )
  expect_error(
    .hypothesis_brma_attach_iwmde_scalar(
      posterior                = list(random = 0),
      raw_posterior            = list(random = 0),
      context                  = NULL,
      estimate_cache           = .iwmde_estimate_cache(),
      parameter                = "mu_alloc",
      parameter_label          = "alloc",
      value                    = 0,
      conditional              = NULL,
      n_points                 = 20,
      max_samples              = 10,
      normalization_points     = 10,
      normalization_prob       = .99,
      density_method           = "qCMDE"
    ),
    "alloc\\[level\\] = 0"
  )
  expect_error(
    .hypothesis_marginal_means_attach_iwmde(
      object = list(
        density_method = "qCMDE",
        source_object  = structure(
          list(fit = list()),
          class = c("RoBMA", "brma")
        ),
        inference      = list(
          conditional = list(mu_alloc = list(random = 0))
        )
      ),
      parameter       = "mu_alloc",
      parameter_label = "alloc",
      hypothesis      = "mu_alloc = 0",
      density_method  = "qCMDE",
      density_control = NULL
    ),
    "alloc\\[level\\] = 0"
  )
})


test_that("qCMDE/IWMDE hypotheses guard non-known-V random-formula objects upfront", {

  object <- .mock_random_non_known_v_brma_mv()

  expect_error(
    hypothesis.brma(
      object,
      "mu = 0",
      density_method  = "qCMDE",
      density_control = list(n_points = 20, max_samples = 20)
    ),
    "qCMDE/IWMDE hypothesis\\(\\).*random-formula"
  )
})


test_that("GLMM IWMDE point Bayes factors fail certification upfront", {

  object <- structure(list(), class = c("brma.glmm", "brma"))

  expect_error(
    .iwmde_check_point_ordinate_supported(object, "IWMDE"),
    "do not meet the bridge-sampling certification tolerance"
  )
  expect_invisible(
    .iwmde_check_point_ordinate_supported(object, "qCMDE")
  )
})


test_that("qCMDE point attachment drops stale same-value ordinates", {

  diagnostics <- list(
    estimator                 = "q_grid_cmde",
    relative_mcse             = .01,
    finite_terms              = 100,
    ess                       = 100,
    max_weight_share          = .10,
    row_drop_fraction         = 0,
    active_mass               = 1,
    final_normalization_integral = 1,
    support_grid_normalization_integral = 1,
    ordinate_relative_change  = 0
  )
  stale <- list(
    value            = 0,
    ordinate         = 1,
    method           = "q_grid_cmde",
    density_method   = "qCMDE",
    diagnostics      = diagnostics,
    iwmde_provenance = list(request_key = "stale")
  )
  fresh <- stale
  fresh[["iwmde_provenance"]] <- list(request_key = "fresh")
  posterior <- stats::rnorm(50)
  attr(posterior, "posterior_ordinate") <- stale
  raw_posterior <- as.numeric(posterior)
  prior_density <- BayesTools:::.prior_linear_density_point(0)
  attr(raw_posterior, "prior_density") <- prior_density
  captured_parameter_spec <- NULL

  testthat::local_mocked_bindings(
    .iwmde_estimate = function(..., parameter_spec) {
      captured_parameter_spec <<- parameter_spec
      list(
        diagnostics        = list(ordinate = list(status = "ok")),
        posterior_ordinate = fresh
      )
    },
    .package = "RoBMA"
  )

  out <- .hypothesis_brma_attach_iwmde_scalar(
    posterior            = posterior,
    raw_posterior        = raw_posterior,
    context              = list(),
    estimate_cache       = .iwmde_estimate_cache(),
    parameter            = "mu",
    parameter_label      = "mu",
    value                = 0,
    conditional          = NULL,
    n_points             = 20,
    max_samples          = 20,
    normalization_points = 20,
    normalization_prob   = .99,
    density_method       = "qCMDE"
  )

  entries <- .iwmde_posterior_ordinate_entries(
    attr(out, "posterior_ordinate", exact = TRUE)
  )
  expect_length(entries, 1L)
  expect_equal(entries[[1L]][["iwmde_provenance"]][["request_key"]], "fresh")
  expect_identical(captured_parameter_spec[["prior_density"]], prior_density)
})


test_that("hypothesis alias rewriting preserves trailing level-reference backticks", {

  rewritten <- .hypothesis_brma_rewrite(
    hypothesis = "alloc > alloc[random]",
    aliases    = list(alloc = "mu_alloc", mu_alloc = "mu_alloc"),
    parameter  = "mu_alloc"
  )

  expect_equal(rewritten, "mu_alloc > mu_alloc[random]")
})


test_that("hypothesis quantities do not advertise normal point tests", {

  testthat::local_mocked_bindings(
    .brma_parameter_catalog = function(object) {

      data.frame(
        alias      = c("mu", "tau", "omega"),
        parameter  = c("mu", "tau", "omega"),
        component  = c("mods", "scale", "bias"),
        term       = c("mu", "tau", "omega"),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    },
    .package = "RoBMA"
  )

  quantities <- hypothesis_quantities(structure(list(), class = "brma"))

  expect_false(any(quantities[["component"]] == "bias"))
  expect_equal(unique(quantities[["point_test_methods"]]), "KDE, qCMDE, IWMDE")
  expect_false(any(grepl("\\bnormal\\b", quantities[["point_test_methods"]])))
})


test_that("hypothesis resolves coefficient aliases on both coefficient scales", {

  skip_on_cran()
  skip_if_missing_fits(c(
    "konstantopoulos2011_3lvl",
    "konstantopoulos2011_3lvl2"
  ))

  fit_full <- load_fit("konstantopoulos2011_3lvl2")
  fit_null <- load_fit("konstantopoulos2011_3lvl")

  .hypothesis_expect_bridge_ready(fit_full)
  .hypothesis_expect_bridge_ready(fit_null)
  .expect_bridge_nesting(fit_null, fit_full, "vi")

  bf_bridge <- .hypothesis_bridge_bf01(fit_null, fit_full)
  bf_alias <- hypothesis(
    fit_full,
    "vi = 0",
    standardized_coefficients = FALSE,
    density_method            = "KDE"
  )
  bf_internal <- hypothesis(
    fit_full,
    "mu_vi = 0",
    standardized_coefficients = FALSE,
    density_method            = "KDE"
  )
  bf_standardized <- hypothesis(
    fit_full,
    "vi = 0",
    standardized_coefficients = TRUE,
    density_method            = "KDE"
  )

  expect_s3_class(bf_alias, "BayesTools_table")
  expect_s3_class(bf_alias, "BayesTools_hypothesis_BF")
  expect_equal(colnames(bf_alias), c("Alternative", "Null", "BF", "BF_error"))
  expect_equal(attr(bf_alias[["BF_error"]], "name"), "error%(BF)")
  expect_equal(attr(bf_alias, "raw_BF"), attr(bf_internal, "raw_BF"),
               tolerance = 1e-12)
  expect_equal(attr(bf_alias, "raw_BF"), attr(bf_standardized, "raw_BF"),
               tolerance = 1e-12)
  expect_equal(log(attr(bf_alias, "raw_BF")), log(1 / bf_bridge),
               tolerance = 0.25)
})


test_that("hypothesis component disambiguates shared location-scale terms", {

  skip_on_cran()
  skip_if_missing_fits("dat.lehmann2018_RoBMA_3lvl_mods_scale")

  fit <- load_fit("dat.lehmann2018_RoBMA_3lvl_mods_scale")
  quantities <- hypothesis_quantities(fit)

  expect_true(all(c("point_test", "direction_test", "reason") %in% names(quantities)))
  expect_equal(unique(quantities[["point_test_methods"]]), "KDE, qCMDE, IWMDE")
  expect_false(any(grepl("\\bnormal\\b", quantities[["point_test_methods"]])))
  expect_true(any(
    quantities[["alias"]] == "Preregistered" &
      quantities[["parameter"]] == "mu_Preregistered" &
      quantities[["component"]] == "mods"
  ))
  expect_true(any(
    quantities[["alias"]] == "Preregistered" &
      quantities[["parameter"]] == "log_tau_Preregistered" &
      quantities[["component"]] == "scale"
  ))
  expect_true(all(quantities[["point_test"]][quantities[["component"]] == "scale"]))
  expect_false(any(nzchar(quantities[["reason"]][quantities[["component"]] == "scale"])))

  expect_error(
    hypothesis(
      fit,
      "Preregistered = 0",
      density_method = "KDE",
      n_samples      = 1000
    ),
    "multiple model parameters"
  )

  bf_mods <- hypothesis(
    fit,
    "Preregistered > 0",
    component      = "mods",
    density_method = "KDE",
    n_samples      = 1000
  )
  bf_location <- hypothesis(
    fit,
    "Preregistered > 0",
    component      = "location",
    density_method = "KDE",
    n_samples      = 1000
  )
  bf_scale <- hypothesis(
    fit,
    "Preregistered > 0",
    component      = "scale",
    density_method = "KDE",
    n_samples      = 1000
  )

  expect_s3_class(bf_mods, "BayesTools_hypothesis_BF")
  expect_s3_class(bf_location, "BayesTools_hypothesis_BF")
  expect_s3_class(bf_scale, "BayesTools_hypothesis_BF")
})


test_that("hypothesis does not advertise or test publication-bias parameters", {

  skip_on_cran()
  skip_if_missing_fits("dat.lehmann2018-3PSM")

  fit <- load_fit("dat.lehmann2018-3PSM")
  quantities <- hypothesis_quantities(fit)

  expect_false(any(quantities[["component"]] == "bias"))
  expect_error(
    hypothesis(
      fit,
      "omega = 1",
      component      = "bias",
      density_method = "KDE",
      n_samples      = 1000
    ),
    "publication-bias parameters are not supported"
  )
})


test_that("hypothesis default output is compact and attaches table warnings", {

  prior     <- data.frame(theta = c(-2, -1, 0, 1))
  posterior <- data.frame(theta = c(-3, -2, -1, 0))

  out <- hypothesis(
    posterior,
    prior      = prior,
    hypothesis = "theta > 0 vs theta <= 0"
  )

  expect_s3_class(out, "BayesTools_table")
  expect_equal(colnames(out), c("Alternative", "Null", "BF", "BF_error"))
  expect_false("warning" %in% colnames(out))
  expect_match(attr(out, "warnings"), "Posterior region mass is zero")
})


test_that("qCMDE point-null ordinates agree with bridge oracle and report error", {

  skip_if_not_certification("This bridge comparison uses a 5,000-row density budget.")
  skip_on_cran()
  skip_if_missing_fits(c(
    "konstantopoulos2011_3lvl",
    "konstantopoulos2011_3lvl2"
  ))

  fit_full <- load_fit("konstantopoulos2011_3lvl2")
  fit_null <- load_fit("konstantopoulos2011_3lvl")

  .hypothesis_expect_bridge_ready(fit_full)
  .hypothesis_expect_bridge_ready(fit_null)
  .expect_bridge_nesting(
    fit_null,
    fit_full,
    "vi"
  )

  bf_qcmde <- hypothesis(
    fit_full,
    c("vi = 0", "vi = 0.01"),
    standardized_coefficients = FALSE,
    columns                   = "all",
    density_method            = "qCMDE",
    density_control           = list(
      n_points             = 60,
      max_samples          = 5000,
      initial_samples      = 5000,
      normalization_points = 100
    ),
    n_samples                 = 1000
  )
  bf_iwmde <- hypothesis(
    fit_full,
    "vi = 0",
    standardized_coefficients = FALSE,
    columns                   = "all",
    density_method            = "IWMDE",
    density_control           = list(
      n_points        = 60,
      max_samples     = 5000,
      initial_samples = 5000
    ),
    n_samples                 = 1000
  )
  bf_point_region <- hypothesis(
    fit_full,
    "vi = 0 vs vi > 0",
    standardized_coefficients = FALSE,
    columns                   = "all",
    density_method            = "qCMDE",
    density_control           = list(
      n_points             = 40,
      max_samples          = 120,
      normalization_points = 50
    ),
    n_samples                 = 1000,
    seed                      = 1
  )

  expect_true(all(bf_qcmde[["method"]] == "Savage-Dickey (precomputed)"))
  expect_true(all(is.finite(bf_qcmde[["BF_error"]])))
  diagnostics <- density_diagnostics(bf_qcmde)
  expect_s3_class(diagnostics, "RoBMA_density_diagnostics")
  expect_equal(nrow(diagnostics), 2L)
  expect_true(all(diagnostics[["achieved_row_budget"]] == 5000L))
  expect_true(all(diagnostics[["relative_mcse"]] < .25))
  .hypothesis_expect_bridge_agreement(
    bf_qcmde,
    fit_null,
    fit_full,
    index = 1L
  )
  expect_equal(bf_iwmde[["method"]], "Savage-Dickey (precomputed)")
  expect_true(is.finite(bf_iwmde[["BF_error"]]))
  .hypothesis_expect_bridge_agreement(bf_iwmde, fit_null, fit_full)
  expect_equal(bf_point_region[["method"]], "transitive Savage-Dickey")
  expect_true(is.finite(bf_point_region[["BF_error"]]))
})


test_that("factor level hypotheses use joint priors and child ordinates", {

  skip_if_not_certification("This case validates fitted qCMDE ordinates.")
  skip_on_cran()
  skip_if_missing_fits("bcg_meta-regression2")

  fit <- load_fit("bcg_meta-regression2")

  level_bf <- hypothesis(
    fit,
    "mu_alloc[alternate] > mu_alloc[random]",
    columns   = "all",
    seed      = 11,
    n_samples = 1000
  )
  point_bf <- hypothesis(
    fit,
    c("alloc[random] = 0", "alloc[random] = 0.01"),
    columns         = "all",
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 40,
      max_samples          = 120,
      normalization_points = 50
    ),
    n_samples       = 1000
  )

  expect_equal(level_bf[["method"]], "prior-posterior odds")
  expect_true(is.finite(attr(level_bf, "raw_BF")))
  expect_true(is.finite(level_bf[["BF_error"]]))
  expect_true(all(point_bf[["method"]] == "Savage-Dickey (precomputed)"))
  expect_true(all(is.finite(point_bf[["BF_error"]])))
  expect_error(
    hypothesis(
      fit,
      "alloc = 0",
      density_method  = "qCMDE",
      density_control = list(
        n_points             = 40,
        max_samples          = 120,
        normalization_points = 50
      ),
      n_samples       = 1000
    ),
    "alloc\\[level\\] = 0"
  )
  expect_error(
    hypothesis(
      fit,
      "alloc[alternate] = 0",
      density_method  = "qCMDE",
      density_control = list(
        n_points             = 40,
        max_samples          = 120,
        normalization_points = 50
      ),
      n_samples       = 1000
    ),
    "linear weights are all zero"
  )
})


test_that("qCMDE point-null ordinates support boundary nulls", {

  skip_if_not_certification("This case certifies boundary density ordinates.")
  skip_on_cran()
  skip_if_missing_fits("bcg_meta-analysis")

  fit <- load_fit("bcg_meta-analysis")
  bf_tau_zero <- hypothesis(
    fit,
    "tau = 0",
    columns         = "all",
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 30,
      max_samples          = 240,
      normalization_points = 60
    ),
    n_samples       = 10000
  )
  bf_tau_zero_iwmde <- hypothesis(
    fit,
    "tau = 0",
    columns         = "all",
    density_method  = "IWMDE",
    density_control = list(
      n_points             = 30,
      max_samples          = 240,
      normalization_points = 80
    ),
    n_samples       = 10000
  )

  expect_equal(bf_tau_zero[["method"]], "Savage-Dickey (precomputed)")
  expect_true(is.finite(attr(bf_tau_zero, "raw_BF")))
  expect_equal(bf_tau_zero_iwmde[["method"]], "Savage-Dickey (precomputed)")
  expect_true(is.finite(attr(bf_tau_zero_iwmde, "raw_BF")))

  context <- .iwmde_context(fit)
  tau_prior_sd <- fit[["priors"]][["outcome"]][["tau"]][["parameters"]][["sd"]]
  oracle_cases <- list(
    qCMDE = list(
      bf      = bf_tau_zero,
      control = list(
        n_points             = 30,
        max_samples          = 240,
        normalization_points = 60
      )
    ),
    IWMDE = list(
      bf      = bf_tau_zero_iwmde,
      control = list(
        n_points             = 30,
        max_samples          = 240,
        normalization_points = 80
      )
    )
  )

  for (method in names(oracle_cases)) {
    case <- oracle_cases[[method]]
    estimate <- .iwmde_estimate(
      context         = context,
      parameter       = "tau",
      density_method  = method,
      density_control = c(
        case[["control"]],
        list(display_grid = "ordinate")
      ),
      outputs         = "ordinate",
      values          = 0,
      parameter_spec  = list(
        type             = "primitive",
        conditional      = NULL,
        conditional_rule = "AND"
      ),
      metadata        = list(parameter = "tau"),
      cache           = .iwmde_estimate_cache()
    )
    ordinate <- estimate[["posterior_ordinate"]]

    expect_equal(ordinate[["value"]], 0)
    expect_equal(ordinate[["evaluation_value"]], 0)
    expect_true(.iwmde_posterior_ordinate_supports_bf(ordinate))

    prior_ordinate <- 2 * stats::dnorm(
      ordinate[["evaluation_value"]],
      mean = 0,
      sd   = tau_prior_sd
    )
    expected_bf <- prior_ordinate / ordinate[["ordinate"]]
    expect_equal(
      log(as.numeric(attr(case[["bf"]], "raw_BF"))),
      log(expected_bf),
      tolerance = 1e-4
    )
  }
})


test_that("marginal means hypothesis wrapper resolves aliases and guards qCMDE", {

  skip_if_not_certification("This case computes fitted marginal qCMDE ordinates.")
  skip_on_cran()
  skip_if_missing_fits("bcg_meta-regression2")

  fit <- load_fit("bcg_meta-regression2")
  mm_qcmde <- marginal_means(
    fit,
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 40,
      max_samples          = 240,
      normalization_points = 50
    ),
    n_samples       = 1000
  )

  expect_error(
    hypothesis(
      mm_qcmde,
      "alloc[alternate] > alloc[random]",
      columns = "all",
      seed    = 11
    ),
    "different conditional posterior subsets"
  )
  point_bf <- hypothesis(
    mm_qcmde,
    "alloc[random] = 0",
    columns        = "all",
    density_method = "qCMDE"
  )

  expect_equal(point_bf[["method"]], "Savage-Dickey (precomputed)")
  expect_true(is.finite(point_bf[["BF_error"]]))
  expect_error(
    hypothesis(
      mm_qcmde,
      "alloc = 0",
      density_method = "qCMDE"
    ),
    "alloc\\[level\\] = 0"
  )
  point_bf_default <- hypothesis(
    mm_qcmde,
    "alloc[random] = 0",
    density_method = "qCMDE",
    columns        = "all"
  )
  expect_equal(point_bf_default[["method"]], "Savage-Dickey (precomputed)")

  mm_kde <- marginal_means(fit, n_samples = 1000)
  point_bf_kde <- hypothesis(
    mm_kde,
    "alloc[alternate] = 0",
    columns        = "all",
    density_method = "qCMDE",
    density_control = list(
      n_points             = 40,
      max_samples          = 240,
      normalization_points = 50
    )
  )
  expect_equal(point_bf_kde[["method"]], "Savage-Dickey (precomputed)")
})


test_that("marginal means qCMDE hypotheses compute missing ordinates on demand", {

  condition_key <- "OR\r2\rmu_alloc\rmu_intercept"
  sample <- structure(
    stats::rnorm(40),
    linear_weights             = c(mu_intercept = 1, mu_alloc = 1),
    prior_density              = BayesTools::prior(
      "normal",
      parameters = list(mean = 0, sd = 1)
    ),
    effective_conditional      = c("mu_intercept", "mu_alloc"),
    effective_conditional_rule = "OR",
    condition_key              = condition_key
  )
  inference <- structure(
    list(
      averaged    = list(mu_alloc = list(alternate = sample)),
      conditional = list(mu_alloc = list(alternate = sample)),
      inference   = list()
    ),
    class = c("marginal_inference", "list")
  )
  object <- list(
    inference        = inference,
    parameters       = "mu_alloc",
    term_map         = data.frame(
      term             = "alloc",
      parameter        = "mu_alloc",
      label            = "alloc",
      stringsAsFactors = FALSE
    ),
    conditional_rule = "OR",
    density_method   = "KDE",
    source_object    = structure(list(fit = list()), class = c("RoBMA", "brma"))
  )
  class(object) <- "marginal_means.brma"

  captured <- NULL
  minimal_context <- function(object) {

    list(
      source            = object,
      posterior_samples = cbind(
        mu_intercept = seq(-1, 1, length.out = 40),
        mu_alloc     = seq(.5, -.5, length.out = 40)
      ),
      object = list(
        fit        = list(),
        likelihood = list(family = "normal"),
        data       = list(measure = "yi")
      )
    )
  }
  testthat::local_mocked_bindings(
    .iwmde_context = minimal_context,
    .iwmde_row_states = function(context, rows, parameter = NULL,
                                 parameter_spec = NULL) {

      lapply(rows, function(row) {
        .iwmde_new_row_state(list(baseline_log_q = 0))
      })
    },
    .iwmde_execute_plan_diagnostic = function(
        context, plan, output, execution_cache = NULL,
        diagnostic_cache = NULL) {

      values <- plan[["grids"]][["requested_values"]]
      list(
        status      = "ok",
        target_key  = plan[["target"]][["target_key"]],
        diagnostics = list(
          bf_included         = TRUE,
          bf_value            = values[[1L]],
          bf_evaluation_value = values[[1L]],
          bf_ordinate         = 2,
          bf_mcse             = .01,
          bf_relative_mcse    = .01,
          bf_error_percent    = 1,
          bf_finite_terms     = 60,
          bf_ess              = 30,
          bf_max_weight_share = .2,
          bf_max_log_ratio    = 0,
          active_mass         = 1,
          final_normalization_integral = 1,
          support_grid_normalization_integral = 1,
          bf_ordinate_relative_change = 0,
          max_ordinate_relative_change = 0,
          max_normalizer_relative_change = 0,
          normalization_range    = c(-1, 1),
          estimator              = plan[["method"]],
          weight_method          = "mock"
        )
      )
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    hypothesis_BF = function(posterior, hypothesis, parameter,
                             density_method, ...) {

      captured <<- list(
        posterior      = posterior,
        hypothesis     = hypothesis,
        parameter      = parameter,
        density_method = density_method
      )
      return("ok")
    },
    .package = "BayesTools"
  )

  out <- hypothesis(
    object,
    "alloc[alternate] = 0",
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 20,
      max_samples          = 20,
      normalization_points = 20
    )
  )

  ordinate <- attr(
    captured[["posterior"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_ordinate",
    exact = TRUE
  )

  expect_equal(as.character(out), "ok")
  expect_equal(captured[["density_method"]], "precomputed")
  expect_true(BayesTools::posterior_ordinate_has_value(ordinate, 0))
  expect_equal(ordinate[["condition_key"]], condition_key)
})


test_that("marginal means qCMDE hypotheses reuse only compatible ordinates", {

  condition_key <- "OR\r2\rmu_alloc\rmu_intercept"
  density_control <- .density_control_normalize(
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 20,
      max_samples          = 20,
      normalization_points = 20
    )
  )
  sample <- structure(
    stats::rnorm(40),
    linear_weights             = c(mu_intercept = 1, mu_alloc = 1),
    prior_density              = BayesTools::prior(
      "normal",
      parameters = list(mean = 0, sd = 1)
    ),
    effective_conditional      = c("mu_intercept", "mu_alloc"),
    effective_conditional_rule = "OR",
    condition_key              = condition_key
  )
  make_object <- function(sample) {

    inference <- structure(
      list(
        averaged    = list(mu_alloc = list(alternate = sample)),
        conditional = list(mu_alloc = list(alternate = sample)),
        inference   = list()
      ),
      class = c("marginal_inference", "list")
    )
    object <- list(
      inference        = inference,
      parameters       = "mu_alloc",
      term_map         = data.frame(
        term             = "alloc",
        parameter        = "mu_alloc",
        label            = "alloc",
        stringsAsFactors = FALSE
      ),
      conditional_rule = "OR",
      density_method   = "KDE",
      source_object    = structure(list(fit = list()), class = c("RoBMA", "brma"))
    )
    class(object) <- "marginal_means.brma"

    return(object)
  }
  make_diagnostic <- function(method, ordinate) {

    list(
      status      = "ok",
      diagnostics = list(
        bf_included         = TRUE,
        bf_value            = 0,
        bf_evaluation_value = 0,
        bf_ordinate         = ordinate,
        bf_mcse             = .01,
        bf_relative_mcse    = .01,
        bf_error_percent    = 1,
        bf_finite_terms     = 60,
        bf_ess              = 30,
        bf_max_weight_share = .2,
        bf_max_log_ratio    = 0,
        active_mass         = 1,
        final_normalization_integral = 1,
        support_grid_normalization_integral = 1,
        bf_ordinate_relative_change = 0,
        max_ordinate_relative_change = 0,
        max_normalizer_relative_change = 0,
        normalization_range = c(-1, 1),
        estimator           = method,
        weight_method       = "mock"
      )
    )
  }
  metadata <- .iwmde_posterior_metadata(
    samples   = sample,
    parameter = "mu_alloc",
    level     = "alternate"
  )
  compatible_object <- make_object(sample)
  minimal_context <- function(object) {

    list(
      source            = object,
      posterior_samples = cbind(
        mu_intercept = seq(-1, 1, length.out = 40),
        mu_alloc     = seq(.5, -.5, length.out = 40)
      ),
      object = list(
        fit        = list(),
        likelihood = list(family = "normal"),
        data       = list(measure = "yi")
      )
    )
  }
  row_states_mock <- function(context, rows, parameter = NULL,
                              parameter_spec = NULL) {

    lapply(rows, function(row) {
      .iwmde_new_row_state(list(baseline_log_q = 0))
    })
  }
  testthat::local_mocked_bindings(
    .iwmde_row_states = row_states_mock,
    .package = "RoBMA"
  )
  compatible_spec <- .iwmde_marginal_means_specs(
    marginal_means_object = compatible_object,
    parameter             = "mu_alloc",
    type                  = "conditional",
    levels                = "alternate"
  )[[1L]]
  compatible_plan <- .iwmde_plan(
    context         = minimal_context(compatible_object),
    parameter       = compatible_spec[["label"]],
    density_method  = "qCMDE",
    density_control = density_control,
    outputs         = "ordinate",
    values          = 0,
    parameter_spec  = compatible_spec,
    metadata        = metadata
  )
  compatible_diagnostic <- make_diagnostic("q_grid_cmde", 2)
  compatible_diagnostic[["target_key"]] <- compatible_plan[["target"]][["target_key"]]
  compatible_diagnostic[["plan"]]       <- compatible_plan
  compatible <- .iwmde_posterior_ordinate_attribute(
    diagnostic      = compatible_diagnostic,
    density_method  = "qCMDE",
    metadata        = metadata,
    density_control = density_control
  )
  incompatible <- .iwmde_posterior_ordinate_attribute(
    diagnostic      = make_diagnostic("iwmde", 3),
    density_method  = "IWMDE",
    metadata        = metadata,
    density_control = density_control
  )

  calls <- 0L
  testthat::local_mocked_bindings(
    .iwmde_context = minimal_context,
    .iwmde_execute_plan_diagnostic = function(
        context, plan, output, execution_cache = NULL,
        diagnostic_cache = NULL) {

      calls <<- calls + 1L
      diagnostic <- make_diagnostic(plan[["method"]], 4)
      diagnostic[["target_key"]] <- plan[["target"]][["target_key"]]
      diagnostic
    },
    .package = "RoBMA"
  )

  compatible_sample <- sample
  attr(compatible_sample, "posterior_ordinate") <- compatible
  reused <- .hypothesis_marginal_means_attach_iwmde(
    object          = make_object(compatible_sample),
    parameter       = "mu_alloc",
    parameter_label = "alloc",
    hypothesis      = "mu_alloc[alternate] = 0",
    density_method  = "qCMDE",
    density_control = density_control
  )
  reused_ordinate <- attr(
    reused[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_ordinate",
    exact = TRUE
  )

  expect_equal(calls, 0L)
  expect_equal(reused_ordinate[["ordinate"]], 2)

  incompatible_sample <- sample
  attr(incompatible_sample, "posterior_ordinate") <- incompatible
  recomputed <- .hypothesis_marginal_means_attach_iwmde(
    object          = make_object(incompatible_sample),
    parameter       = "mu_alloc",
    parameter_label = "alloc",
    hypothesis      = "mu_alloc[alternate] = 0",
    density_method  = "qCMDE",
    density_control = density_control
  )
  recomputed_ordinate <- attr(
    recomputed[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_ordinate",
    exact = TRUE
  )

  expect_equal(calls, 1L)
  expect_equal(recomputed_ordinate[["ordinate"]], 4)
  expect_equal(recomputed_ordinate[["density_method"]], "qCMDE")
})


test_that("marginal means hypothesis BFs use alternative-conditioned marginals", {

  make_parameter <- function(label) {

    level <- structure(
      stats::rnorm(20),
      class = c("marginal_posterior.simple", "numeric"),
      marker = label
    )
    out <- structure(
      list(alternate = level),
      class     = c("marginal_posterior.factor", "list"),
      parameter = "mu_alloc"
    )

    return(out)
  }

  object <- list(
    inference = structure(
      list(
        averaged    = list(mu_alloc = make_parameter("averaged")),
        conditional = list(mu_alloc = make_parameter("conditional")),
        inference   = list()
      ),
      class = c("marginal_inference", "list")
    ),
    term_map = data.frame(
      term             = "alloc",
      parameter        = "mu_alloc",
      label            = "alloc",
      check.names      = FALSE,
      stringsAsFactors = FALSE
    ),
    density_method = "KDE"
  )
  class(object) <- "marginal_means.brma"

  captured <- list()
  testthat::local_mocked_bindings(
    hypothesis_BF = function(posterior, hypothesis, parameter, density_method,
                             ...) {

      captured[[length(captured) + 1L]] <<- list(
        posterior      = posterior,
        hypothesis     = hypothesis,
        parameter      = parameter,
        density_method = density_method
      )

      return("ok")
    },
    .package = "BayesTools"
  )

  expect_equal(
    hypothesis(object, "alloc[alternate] > 0"),
    "ok"
  )
  expect_equal(
    hypothesis(object, "alloc[alternate] = 0"),
    "ok"
  )
  expect_warning(
    expect_equal(
      hypothesis(object, "alloc[alternate] > 0", type = "conditional"),
      "ok"
    ),
    "Unused argument.*'type'"
  )

  for (call in captured) {
    expect_identical(
      call[["posterior"]][["conditional"]],
      object[["inference"]][["conditional"]]
    )
    expect_equal(call[["parameter"]], "mu_alloc")
    expect_equal(call[["density_method"]], "KDE")
  }
  expect_equal(
    vapply(captured, `[[`, character(1), "hypothesis"),
    c(
      "mu_alloc[alternate] > 0",
      "mu_alloc[alternate] = 0",
      "mu_alloc[alternate] > 0"
    )
  )
})


test_that("marginal means hypothesis density methods use public names", {

  expect_false("type" %in% names(formals(hypothesis.marginal_means.brma)))
  expect_null(formals(hypothesis.marginal_means.brma)[["density_method"]])
})


test_that("marginal means hypothesis inherits and can override stored density method", {

  posterior <- list(
    mu_alloc = structure(
      list(alternate = stats::rnorm(20)),
      class = c("marginal_posterior.factor", "list")
    )
  )
  object <- list(
    inference = structure(
      list(
        averaged    = posterior,
        conditional = posterior,
        inference   = list()
      ),
      class = c("marginal_inference", "list")
    ),
    term_map = data.frame(
      term             = "alloc",
      parameter        = "mu_alloc",
      label            = "alloc",
      stringsAsFactors = FALSE
    ),
    density_method = "qCMDE"
  )
  class(object) <- "marginal_means.brma"

  captured <- character()
  testthat::local_mocked_bindings(
    hypothesis_BF = function(..., density_method) {
      captured <<- c(captured, density_method)
      return("ok")
    },
    .package = "BayesTools"
  )

  expect_equal(hypothesis(object, "alloc[alternate] > 0"), "ok")
  expect_equal(
    hypothesis(
      object,
      "alloc[alternate] > 0",
      density_method = "KDE"
    ),
    "ok"
  )
  expect_equal(captured, c("precomputed", "KDE"))
})


test_that("marginal means hypothesis rejects normal density method", {

  skip_on_cran()
  skip_if_missing_fits("bcg_meta-regression2")

  fit <- load_fit("bcg_meta-regression2")
  mm  <- marginal_means(fit, n_samples = 1000)

  expect_error(
    hypothesis(
      mm,
      "alloc[random] = 0",
      columns        = "all",
      density_method = "normal"
    ),
    "density_method.*must be one of"
  )
})


test_that("marginal means hypothesis validates controls without a point null", {

  skip_on_cran()
  skip_if_missing_fits("bcg_meta-regression2")

  fit <- load_fit("bcg_meta-regression2")
  mm  <- marginal_means(fit, n_samples = 1000)

  expect_error(
    hypothesis(
      mm,
      "alloc[random] > 0",
      density_method  = "qCMDE",
      density_control = list(unknown = 1)
    ),
    "unrecognized setting"
  )
})


test_that("factor-level PET moderator point null agrees with bridge oracle", {

  skip_if_not_certification("This bridge comparison uses 2,000-row density budgets.")
  skip_on_cran()
  skip_if_missing_fits(c(
    "dat.lehmann2018-PET",
    "dat.lehmann2018-PETreg"
  ))

  fit_full <- load_fit("dat.lehmann2018-PETreg")
  fit_null <- load_fit("dat.lehmann2018-PET")

  .hypothesis_expect_bridge_ready(fit_full)
  .hypothesis_expect_bridge_ready(fit_null)
  .expect_bridge_nesting(
    fit_null,
    fit_full,
    "Preregistered"
  )

  bf_bridge <- .hypothesis_bridge_bf01(fit_null, fit_full)
  bf_kde <- hypothesis(
    fit_full,
    "Preregistered[Pre-Registered] = 0",
    density_method = "KDE",
    n_samples      = 2000
  )
  bf_qcmde <- hypothesis(
    fit_full,
    "Preregistered[Pre-Registered] = 0",
    columns         = "all",
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 60,
      max_samples          = 2000,
      initial_samples      = 2000,
      normalization_points = 100
    ),
    n_samples       = 1000
  )
  bf_iwmde <- hypothesis(
    fit_full,
    "Preregistered[Pre-Registered] = 0",
    columns         = "all",
    density_method  = "IWMDE",
    density_control = list(
      n_points        = 60,
      max_samples     = 2000,
      initial_samples = 2000
    ),
    n_samples       = 1000
  )

  expect_equal(log(attr(bf_kde, "raw_BF")), log(1 / bf_bridge),
               tolerance = 0.35)
  .hypothesis_expect_bridge_agreement(bf_qcmde, fit_null, fit_full)
  .hypothesis_expect_bridge_agreement(bf_iwmde, fit_null, fit_full)
  expect_equal(bf_qcmde[["method"]], "Savage-Dickey (precomputed)")
  expect_true(is.finite(bf_qcmde[["BF_error"]]))
  expect_equal(bf_iwmde[["method"]], "Savage-Dickey (precomputed)")
  expect_true(is.finite(bf_iwmde[["BF_error"]]))
})


test_that("qCMDE requires direct point-null expressions", {

  skip_on_cran()
  skip_if_missing_fits("konstantopoulos2011_3lvl2")

  fit <- load_fit("konstantopoulos2011_3lvl2")

  expect_error(
    hypothesis(
      fit,
      "vi + 0 = 0",
      density_method  = "qCMDE",
      density_control = list(
        n_points             = 40,
        max_samples          = 120,
        normalization_points = 50
      ),
      n_samples       = 1000
    ),
    "direct parameter or level point hypotheses"
  )
})
