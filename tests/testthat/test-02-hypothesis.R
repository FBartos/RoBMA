context("Hypothesis Bayes factors")

source(testthat::test_path("common-functions.R"))


.hypothesis_expect_bridge_ready <- function(object) {

  expect_true(is.finite(logml(object)))
  expect_true(is.finite(object[["marglik"]][["mcse_logml"]]))
  expect_lt(object[["marglik"]][["mcse_logml"]], 0.05)
}


.hypothesis_bridge_bf01 <- function(null, full) {

  exp(logml(null) - logml(full))
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
      diagnostic_cache         = NULL,
      parameter                = "mu_alloc",
      parameter_label          = "alloc",
      value                    = 0,
      conditional              = NULL,
      max_samples              = 10,
      normalization_points     = 10,
      normalization_prob       = .99,
      method                   = "q_grid_cmde",
      density_method           = "qCMDE",
      standardized_coefficients = FALSE
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

  skip_on_cran()
  skip_if_missing_fits(c(
    "konstantopoulos2011_3lvl",
    "konstantopoulos2011_3lvl2"
  ))

  fit_full <- load_fit("konstantopoulos2011_3lvl2")
  fit_null <- load_fit("konstantopoulos2011_3lvl")

  .hypothesis_expect_bridge_ready(fit_full)
  .hypothesis_expect_bridge_ready(fit_null)

  bf_bridge <- .hypothesis_bridge_bf01(fit_null, fit_full)
  bf_qcmde <- hypothesis(
    fit_full,
    c("vi = 0", "vi = 0.01"),
    standardized_coefficients = FALSE,
    columns                   = "all",
    density_method            = "qCMDE",
    density_control           = list(
      n_points             = 40,
      max_samples          = 120,
      normalization_points = 50
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
      n_points    = 40,
      max_samples = 120
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
  expect_equal(log(attr(bf_qcmde, "raw_BF")[[1L]]), log(1 / bf_bridge),
               tolerance = 0.25)
  expect_equal(bf_iwmde[["method"]], "Savage-Dickey (precomputed)")
  expect_true(is.finite(bf_iwmde[["BF_error"]]))
  expect_equal(log(attr(bf_iwmde, "raw_BF")), log(1 / bf_bridge),
               tolerance = 0.25)
  expect_equal(bf_point_region[["method"]], "transitive Savage-Dickey")
  expect_true(is.finite(bf_point_region[["BF_error"]]))
})


test_that("factor level hypotheses use joint priors and child ordinates", {

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
      max_samples          = 80,
      normalization_points = 40
    ),
    n_samples       = 1000
  )

  expect_equal(bf_tau_zero[["method"]], "Savage-Dickey (precomputed)")
  expect_true(is.finite(attr(bf_tau_zero, "raw_BF")))
})


test_that("marginal means hypothesis wrapper resolves aliases and guards qCMDE", {

  skip_on_cran()
  skip_if_missing_fits("bcg_meta-regression2")

  fit <- load_fit("bcg_meta-regression2")
  mm_qcmde <- marginal_means(
    fit,
    density_method  = "qCMDE",
    density_control = list(
      n_points             = 40,
      max_samples          = 120,
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
    type           = "conditional",
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
    type           = "conditional",
    columns        = "all",
    density_method = "qCMDE",
    density_control = list(
      n_points             = 40,
      max_samples          = 120,
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
  class(object) <- c("marginal_means.brma", "marginal_means")

  captured <- NULL
  testthat::local_mocked_bindings(
    .iwmde_context = function(object) list(source = object),
    .iwmde_parameter_ordinate_diagnostic = function(
        context, parameter, values, max_samples, normalization_points,
        normalization_prob, method, parameter_spec, diagnostic_cache) {

      list(
        status      = "ok",
        diagnostics = list(
          bf_included         = TRUE,
          bf_value            = values[[1L]],
          bf_evaluation_value = values[[1L]],
          bf_ordinate         = 2,
          bf_mcse             = .01,
          bf_relative_mcse    = .01,
          bf_error_percent    = 1,
          bf_finite_terms     = 20,
          bf_ess              = 20,
          bf_max_weight_share = .2,
          bf_max_log_ratio    = 0,
          active_mass         = 1,
          normalization_integral = 1,
          normalization_range    = c(-1, 1),
          estimator              = method,
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

  expect_equal(out, "ok")
  expect_equal(captured[["density_method"]], "precomputed")
  expect_true(BayesTools::posterior_ordinate_has_value(ordinate, 0))
  expect_equal(ordinate[["condition_key"]], condition_key)
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
  class(object) <- c("marginal_means.brma", "marginal_means")

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
    hypothesis(object, "alloc[alternate] > 0", type = "averaged"),
    "ok"
  )
  expect_equal(
    hypothesis(object, "alloc[alternate] = 0", type = "averaged"),
    "ok"
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
    c("mu_alloc[alternate] > 0", "mu_alloc[alternate] = 0")
  )
})


test_that("marginal means hypothesis density methods use public names", {

  expect_equal(
    eval(formals(hypothesis.marginal_means.brma)[["density_method"]]),
    c("KDE", "qCMDE", "IWMDE")
  )
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
      type           = "conditional",
      columns        = "all",
      density_method = "normal"
    ),
    "density_method.*must be one of"
  )
})


test_that("factor-level PET moderator point null agrees with bridge oracle", {

  skip_on_cran()
  skip_if_missing_fits(c(
    "dat.lehmann2018-PET",
    "dat.lehmann2018-PETreg"
  ))

  fit_full <- load_fit("dat.lehmann2018-PETreg")
  fit_null <- load_fit("dat.lehmann2018-PET")

  .hypothesis_expect_bridge_ready(fit_full)
  .hypothesis_expect_bridge_ready(fit_null)

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
      n_points             = 40,
      max_samples          = 120,
      normalization_points = 50
    ),
    n_samples       = 1000
  )

  expect_equal(log(attr(bf_kde, "raw_BF")), log(1 / bf_bridge),
               tolerance = 0.35)
  expect_equal(log(attr(bf_qcmde, "raw_BF")), log(1 / bf_bridge),
               tolerance = 0.35)
  expect_equal(bf_qcmde[["method"]], "Savage-Dickey (precomputed)")
  expect_true(is.finite(bf_qcmde[["BF_error"]]))
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
