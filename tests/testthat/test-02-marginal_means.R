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
