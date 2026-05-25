context("Estimated marginal means")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "marginal_means")

.marginal_means_test_object <- function() {

  p <- stats::ppoints(200)

  levels <- list(
    alternate  = stats::qnorm(p, mean = -0.60, sd = 0.20),
    random     = stats::qnorm(p, mean = -0.20, sd = 0.25),
    systematic = stats::qnorm(p, mean =  0.20, sd = 0.30)
  )
  class(levels) <- c("marginal_posterior.factor", "marginal_posterior")

  samples   <- list(mu_alloc = levels)
  inference <- list(averaged = samples, conditional = samples, inference = list())
  class(inference) <- c("marginal_inference", "list")

  emm <- list(
    inference        = inference,
    parameters       = "mu_alloc",
    term_map         = data.frame(
      term             = "alloc",
      parameter        = "mu_alloc",
      label            = "alloc",
      stringsAsFactors = FALSE
    ),
    input_measure    = "GEN",
    effect_transform = .effect_output_setup_measure(input_measure = "GEN"),
    model_averaged   = FALSE,
    bf               = FALSE
  )
  class(emm) <- c("marginal_means.brma", "marginal_means")

  return(emm)
}


test_that("marginal_means base plot keeps separate level colors", {

  emm  <- .marginal_means_test_object()
  levels <- length(emm[["inference"]][["averaged"]][["mu_alloc"]])
  dots   <- .set_dots_plot(n_levels = levels)

  expect_gte(length(unique(dots[["col"]])), levels)
  .with_temp_plot_device(expect_silent(plot(emm, parameter = "alloc", legend = FALSE)))
})

test_that("plot.marginal_means warns on unused dots", {

  emm <- .marginal_means_test_object()

  expect_warning(
    plot(emm, parameter = "alloc", iwmde_n_points = 20),
    "Unused argument.*iwmde_n_points"
  )
})


# list cached fits lazily
skip_if_no_fits()
fit_names <- unique(c(
  "bcg_meta-analysis",
  marginal_means_cases()[["name"]],
  marginal_means_interaction_plot_cases()[["name"]]
))
fits <- lazy_fits(fit_names, validate = FALSE)


marginal_means_expected_stats <- function(emm, type = "averaged",
                                          effect_transform = emm[["effect_transform"]],
                                          probs = c(.025, .50, .975)) {

  samples <- .transform_marginal_samples_effect(
    samples          = emm[["inference"]][[type]],
    effect_transform = effect_transform
  )

  stats <- unlist(lapply(emm[["parameters"]], function(parameter) {

    lapply(samples[[parameter]], function(draws) {

      draws <- as.numeric(draws)
      c(
        Mean = mean(draws),
        SD   = stats::sd(draws),
        stats::quantile(draws, probs = probs, names = FALSE)
      )
    })
  }), recursive = FALSE)

  stats <- do.call(rbind, stats)
  colnames(stats) <- c("Mean", "SD", as.character(probs))

  return(stats)
}


test_that("marginal_means stores BayesTools marginal inference", {

  mm <- marginal_means(fits[["bcg_meta-regression2"]], n_samples = 1000)

  expect_s3_class(mm, "marginal_means.brma")
  expect_s3_class(mm[["inference"]], "marginal_inference")
  expect_equal(mm[["parameters"]], c("mu_intercept", "mu_alloc"))
  expect_equal(mm[["term_map"]][["term"]], c("intercept", "alloc"))
})


test_that("marginal_means attaches qCMDE densities and refreshes BFs", {

  mm <- marginal_means(
    fits[["bcg_meta-regression2"]],
    n_samples         = 1000,
    bf                = TRUE,
    density_method    = "qCMDE",
    density_control   = list(n_points = 20, max_samples = 20)
  )

  averaged_density <- attr(
    mm[["inference"]][["averaged"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  )
  conditional_density <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  )

  expect_equal(averaged_density[["status"]], "ok")
  expect_equal(conditional_density[["status"]], "ok")
  expect_equal(averaged_density[["density_method"]], "qCMDE")
  expect_true(all(is.finite(averaged_density[["x"]])))
  expect_true(all(is.finite(averaged_density[["y"]])))
  expect_true(all(averaged_density[["y"]] >= 0))
  expect_equal(
    attr(
      mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
      "posterior_ordinate"
    )[["value"]],
    mm[["null_hypothesis"]]
  )
  expect_false(is.null(attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "prior_density"
  )))

  refreshed <- .marginal_means_refresh_iwmde_bf(
    inference            = mm[["inference"]],
    parameters           = mm[["parameters"]],
    null_hypothesis      = mm[["null_hypothesis"]],
    normal_approximation = mm[["normal_approximation"]]
  )
  expected <- BayesTools::Savage_Dickey_BF(
    posterior            = mm[["inference"]][["conditional"]][["mu_alloc"]],
    null_hypothesis      = mm[["null_hypothesis"]],
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(refreshed[["inference"]][["mu_alloc"]], expected)
  expect_true(is.finite(attr(expected[["alternate"]], "BF_error_percent")))
})


test_that("marginal_means IWMDE ordinates do not expand plot densities", {

  mm <- marginal_means(
    fits[["bcg_meta-regression2"]],
    null_hypothesis  = 5,
    n_samples        = 1000,
    bf               = TRUE,
    density_method   = "qCMDE",
    density_control  = list(n_points = 20, max_samples = 20)
  )

  conditional_density <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  )
  conditional_ordinate <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_ordinate"
  )
  ordinate_diagnostic <- mm[["iwmde_ordinate_diagnostics"]][["conditional"]][["alloc: alternate"]]

  expect_true(max(conditional_density[["x"]]) < mm[["null_hypothesis"]])
  expect_null(conditional_ordinate)
  expect_equal(ordinate_diagnostic[["iwmde"]][["x"]], mm[["null_hypothesis"]])
  expect_true(max(ordinate_diagnostic[["xlim"]]) < mm[["null_hypothesis"]])
  expect_true(max(ordinate_diagnostic[["diagnostics"]][["normalization_range"]]) >=
                mm[["null_hypothesis"]])
  expect_false(.marginal_means_has_bf_posterior_ordinate(
    mm[["inference"]][["conditional"]][["mu_alloc"]]
  ))
  expect_true(all(!is.finite(unlist(mm[["inference"]][["inference"]][["mu_alloc"]]))))
})


test_that("marginal_means precomputes qCMDE ordinates when BFs are hidden", {

  mm <- marginal_means(
    fits[["bcg_meta-regression2"]],
    n_samples         = 1000,
    bf                = FALSE,
    density_method    = "qCMDE",
    density_control   = list(n_points = 20, max_samples = 20)
  )

  posterior_ordinate <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_ordinate"
  )
  expected <- BayesTools::Savage_Dickey_BF(
    posterior            = mm[["inference"]][["conditional"]][["mu_alloc"]],
    null_hypothesis      = mm[["null_hypothesis"]],
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )
  hidden_summary <- summary(mm)
  bf_summary     <- summary(mm, bf = TRUE)

  expect_false(mm[["bf"]])
  expect_equal(posterior_ordinate[["value"]], mm[["null_hypothesis"]])
  expect_true(.marginal_means_has_bf_posterior_ordinate(
    mm[["inference"]][["conditional"]][["mu_alloc"]]
  ))
  expect_equal(mm[["inference"]][["inference"]][["mu_alloc"]], expected)
  expect_false("inclusion_BF" %in% attr(hidden_summary, "type"))
  expect_true("inclusion_BF" %in% attr(bf_summary, "type"))
  expect_true(is.finite(attr(expected[["alternate"]], "BF_error_percent")))
})


test_that("marginal_means refreshes BFs from BF-grade IWMDE densities", {

  mm <- marginal_means(
    fits[["bcg_meta-regression2"]],
    n_samples         = 1000,
    bf                = TRUE,
    density_method    = "IWMDE",
    density_control   = list(n_points = 80, max_samples = 200)
  )

  conditional_density <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  )
  refreshed <- .marginal_means_refresh_iwmde_bf(
    inference            = mm[["inference"]],
    parameters           = mm[["parameters"]],
    null_hypothesis      = mm[["null_hypothesis"]],
    normal_approximation = mm[["normal_approximation"]]
  )

  expect_equal(conditional_density[["density_method"]], "IWMDE")
  expect_equal(conditional_density[["diagnostics"]][["estimator"]], "iwmde")
  expect_equal(
    attr(
      mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
      "posterior_ordinate"
    )[["method"]],
    "iwmde"
  )
  expect_true(.marginal_means_has_bf_posterior_ordinate(
    mm[["inference"]][["conditional"]][["mu_alloc"]]
  ))

  expected <- BayesTools::Savage_Dickey_BF(
    posterior            = mm[["inference"]][["conditional"]][["mu_alloc"]],
    null_hypothesis      = mm[["null_hypothesis"]],
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(refreshed[["inference"]][["mu_alloc"]], expected)
  expect_true(is.finite(attr(expected[["alternate"]], "BF_error_percent")))
})


test_that("marginal_means does not mix qCMDE and KDE BF refresh within a parameter", {

  samples <- list(
    level_a = structure(
      rnorm(20),
      posterior_density = list(diagnostics = list(estimator = "q_grid_cmde"))
    ),
    level_b = rnorm(20)
  )

  expect_false(.marginal_means_has_bf_posterior_ordinate(samples))
})


test_that("marginal_means rejects IWMDE BFs with poor ordinate diagnostics", {

  posterior_ordinate <- list(
    value       = 0,
    ordinate    = .4,
    diagnostics = list(
      estimator           = "iwmde",
      relative_mcse       = 2,
      finite_terms        = 20,
      ess                 = 10,
      max_weight_share    = .2
    )
  )

  expect_false(.iwmde_posterior_ordinate_supports_bf(posterior_ordinate))
  posterior_ordinate[["diagnostics"]][["relative_mcse"]] <- 1
  posterior_ordinate[["diagnostics"]][["ess"]] <- 10
  expect_false(.iwmde_posterior_ordinate_supports_bf(posterior_ordinate))
  posterior_ordinate[["diagnostics"]][["relative_mcse"]] <- .1
  posterior_ordinate[["diagnostics"]][["ess"]] <- 2
  expect_false(.iwmde_posterior_ordinate_supports_bf(posterior_ordinate))
  posterior_ordinate[["diagnostics"]][["ess"]] <- 10
  posterior_ordinate[["diagnostics"]][["max_weight_share"]] <- .95
  expect_false(.iwmde_posterior_ordinate_supports_bf(posterior_ordinate))

  qcmde_ordinate <- list(
    value       = 0,
    ordinate    = .4,
    diagnostics = list(
      estimator              = "q_grid_cmde",
      relative_mcse          = .1,
      finite_terms           = 20,
      ess                    = 10,
      max_weight_share       = .2,
      active_mass            = .5,
      normalization_integral = .5,
      normalization_range    = c(-1, 1)
    )
  )
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))
  qcmde_ordinate[["diagnostics"]][["normalization_integral"]] <- NA_real_
  expect_false(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))
  qcmde_ordinate[["diagnostics"]][["normalization_integral"]] <- .5
  qcmde_ordinate[["diagnostics"]][["active_mass"]] <- NA_real_
  expect_false(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))
  qcmde_ordinate[["diagnostics"]][["active_mass"]] <- 0
  expect_false(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))
  qcmde_ordinate[["diagnostics"]][["active_mass"]] <- .5
  qcmde_ordinate[["diagnostics"]][["normalization_range"]] <- c(.25, 1)
  expect_false(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))

  diagnostic <- list(
    status      = "ok",
    diagnostics = list(
      bf_included         = TRUE,
      bf_value            = 0,
      bf_ordinate         = .4,
      bf_mcse             = .8,
      bf_relative_mcse    = 2,
      bf_error_percent    = 200,
      bf_finite_terms     = 20,
      bf_ess              = 10,
      bf_max_weight_share = .2,
      bf_max_log_ratio    = 1,
      estimator           = "iwmde",
      weight_method       = "test"
    )
  )
  expect_null(.iwmde_posterior_ordinate_attribute(diagnostic, "IWMDE"))
})


test_that("marginal_means hides BFs for non-RoBMA fits by default", {

  emm_brma <- marginal_means(fits[["bcg_meta-regression2"]], n_samples = 1000)
  emm_brma_bf <- marginal_means(
    fits[["bcg_meta-regression2"]],
    n_samples = 1000,
    bf        = TRUE
  )
  emm_robma <- marginal_means(
    fits[["dat.lehmann2018_RoBMA_mods"]],
    n_samples = 1000
  )

  expect_false("inclusion_BF" %in% attr(summary(emm_brma), "type"))
  expect_false(any(grepl(
    "Savage-Dickey",
    attr(summary(emm_brma), "warnings"),
    fixed = TRUE
  )))
  expect_true("inclusion_BF" %in% attr(summary(emm_brma_bf), "type"))
  expect_true("inclusion_BF" %in% attr(summary(emm_brma, bf = TRUE), "type"))
  expect_true("inclusion_BF" %in% attr(summary(emm_robma), "type"))
})


test_that("marginal_means plot errors show formula terms", {

  emm <- marginal_means(fits[["bcg_meta-regression2"]], n_samples = 1000)

  expect_error(
    plot(emm),
    "^The 'parameter' argument must be specified\\. Available terms are: 'intercept', 'alloc'\\.$"
  )
  expect_error(
    plot(emm, parameter = "missing"),
    "^Unknown marginal means parameter 'missing'\\. Available terms are: 'intercept', 'alloc'\\.$"
  )
})


test_that("marginal_means plot labels effect axis and term legend", {

  emm  <- marginal_means(fits[["bcg_meta-regression2"]], n_samples = 1000)
  plot <- plot(emm, parameter = "alloc", plot_type = "ggplot")
  colour_scale   <- plot[["scales"]][["get_scales"]]("colour")
  linetype_scale <- plot[["scales"]][["get_scales"]]("linetype")

  expect_identical(
    plot[["scales"]][["get_scales"]]("x")[["name"]],
    "Effect Size"
  )
  expect_identical(colour_scale[["name"]], "alloc")
  expect_identical(linetype_scale[["name"]], "alloc")

  plot_custom <- plot(
    emm,
    parameter = "alloc",
    plot_type = "ggplot",
    xlab      = "Custom Label"
  )
  expect_identical(
    plot_custom[["scales"]][["get_scales"]]("x")[["name"]],
    "Custom Label"
  )

  plot_custom_legend <- plot(
    emm,
    parameter    = "alloc",
    plot_type    = "ggplot",
    legend_title = "Custom Legend"
  )
  expect_identical(
    plot_custom_legend[["scales"]][["get_scales"]]("colour")[["name"]],
    "Custom Legend"
  )
})


test_that("marginal_means plot forwards stored posterior density", {

  emm <- .marginal_means_test_object()
  density <- list(
    status       = "ok",
    x            = seq(-1, 1, length.out = 5),
    y            = rep(.5, 5),
    method       = "q_grid_cmde",
    diagnostics  = list(estimator = "q_grid_cmde"),
    point_masses = data.frame(x = numeric(), mass = numeric())
  )
  attr(
    emm[["inference"]][["averaged"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  ) <- density

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_marginal = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(emm, parameter = "alloc", plot_type = "ggplot")

  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu_alloc")
  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
  expect_identical(
    attr(captured[["samples"]][["mu_alloc"]][["alternate"]], "posterior_density"),
    density
  )
})

test_that("marginal_means summaries transform effect-size scale", {

  emm <- marginal_means(
    fits[["bcg_meta-regression2"]],
    n_samples = 1000,
    transform = "EXP"
  )

  summary_rr <- summary(emm)
  expected   <- marginal_means_expected_stats(emm)
  actual     <- as.data.frame(summary_rr)[, colnames(expected)]

  expect_equal(
    unname(as.matrix(actual)),
    unname(expected),
    tolerance = sqrt(.Machine$double.eps)
  )
  expect_equal(attr(summary_rr, "title"), "Marginal Means (risk ratio):")
  expect_match(attr(summary_rr, "footnotes"), "risk ratio")
})


test_that("marginal_means summaries convert effect-size measures", {

  emm <- marginal_means(
    fits[["dat.lehmann2018_BMA.norm_mods"]],
    n_samples = 1000
  )
  effect_transform <- .effect_output_setup_measure(
    input_measure  = emm[["input_measure"]],
    output_measure = "COR"
  )

  summary_cor <- summary(emm, output_measure = "COR")
  expected    <- marginal_means_expected_stats(
    emm              = emm,
    effect_transform = effect_transform
  )
  actual      <- as.data.frame(summary_cor)[, colnames(expected)]

  expect_equal(
    unname(as.matrix(actual)),
    unname(expected),
    tolerance = sqrt(.Machine$double.eps)
  )
  expect_equal(
    attr(summary_cor, "title"),
    "Model-Averaged Marginal Means (correlation):"
  )
  expect_match(attr(summary_cor, "footnotes"), "correlation")
})


test_that("marginal_means summaries match reference tables", {

  cases <- marginal_means_cases()
  for (i in seq_len(nrow(cases))) {

    name <- cases[["name"]][i]
    fit  <- fits[[name]]
    emm  <- marginal_means(fit, n_samples = 1000)

    if (inherits(fit, "RoBMA")) {
      test_reference_table(
        summary(emm),
        paste0("summary_marginal-", name, "-averaged.txt"),
        paste0("Averaged marginal means reference table mismatch for '", name, "'")
      )

      test_reference_table(
        summary(emm, type = "conditional"),
        paste0("summary_marginal-", name, "-conditional.txt"),
        paste0("Conditional marginal means reference table mismatch for '", name, "'")
      )
    } else {
      test_reference_table(
        summary(emm),
        paste0("summary_marginal-", name, ".txt"),
        paste0("Marginal means reference table mismatch for '", name, "'")
      )
    }
  }
})


test_that("marginal_means core plot snapshots are stable", {

  for_each_case(marginal_means_cases(), function(case) {

    name      <- case_name(case)
    parameter <- case_value(case, "parameter")
    fit       <- fits[[name]]
    emm       <- marginal_means(fit, n_samples = 1000)

    expect_brma_plot_snapshot(
      paste0(name, "-marginal_means-baseplot-no-prior"),
      function() plot(emm, parameter = parameter, prior = FALSE)
    )

    expect_brma_plot_snapshot(
      paste0(name, "-marginal_means-ggplot-no-prior"),
      plot(emm, parameter = parameter, prior = FALSE, plot_type = "ggplot")
    )
  }, tier = visual_test_tier())
})

test_that("marginal_means prior plot gallery snapshots are stable", {

  skip_if_not_full_visuals("Marginal-means prior overlays are gallery coverage.")

  for_each_case(marginal_means_cases(), function(case) {

    name      <- case_name(case)
    parameter <- case_value(case, "parameter")
    fit       <- fits[[name]]
    emm       <- marginal_means(fit, n_samples = 1000)

    expect_brma_plot_snapshot(
      paste0(name, "-marginal_means-baseplot-prior"),
      function() plot(emm, parameter = parameter, prior = TRUE)
    )

    expect_brma_plot_snapshot(
      paste0(name, "-marginal_means-ggplot-prior"),
      plot(emm, parameter = parameter, prior = TRUE, plot_type = "ggplot")
    )
  }, tier = visual_test_tier())
})


test_that("marginal_means interaction plots render moderator type combinations", {

  skip_if_missing_fits(unique(marginal_means_interaction_plot_cases()[["name"]]))

  for_each_case(marginal_means_interaction_plot_cases(), function(case) {

    name      <- case_name(case)
    parameter <- case_value(case, "parameter")
    fit       <- fits[[name]]
    emm       <- marginal_means(fit, n_samples = 1000)

    expect_true(parameter %in% emm[["term_map"]][["term"]])

    for (show_prior in c(FALSE, TRUE)) {
      expect_s3_class(
        plot(
          emm,
          parameter = parameter,
          prior     = show_prior,
          plot_type = "ggplot"
        ),
        "ggplot"
      )

      .with_temp_plot_device(expect_silent(plot(
        emm,
        parameter = parameter,
        prior     = show_prior
      )))
    }
  }, tier = test_tier())
})


test_that("marginal_means requires moderators", {

  expect_error(
    marginal_means(fits[["bcg_meta-analysis"]], n_samples = 1000),
    "moderators"
  )
})


test_that("marginal_means conditional type is RoBMA-only", {

  emm <- marginal_means(fits[["bcg_meta-regression2"]], n_samples = 1000)

  expect_error(
    summary(emm, type = "conditional"),
    "RoBMA marginal means"
  )
  expect_error(
    plot(emm, parameter = "alloc", type = "conditional", plot_type = "ggplot"),
    "RoBMA marginal means"
  )

  emm_robma <- marginal_means(
    fits[["dat.lehmann2018_RoBMA_mods"]],
    n_samples = 1000
  )
  expect_s3_class(
    summary(emm_robma, type = "conditional"),
    "summary.marginal_means.brma"
  )
})
