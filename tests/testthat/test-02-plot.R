context("Prior and posterior plots")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))
source(testthat::test_path("helper-iwmde.R"))

test_that("plot.brma clears stale posterior density before qCMDE attach", {

  samples <- list(mu = stats::rnorm(20))
  attr(samples[["mu"]], "posterior_density") <- list(
    x      = seq(-1, 1, length.out = 5),
    y      = rep(.5, 5),
    method = "q_grid_cmde"
  )
  attr(samples[["mu"]], "posterior_densities") <- list(mu = list(
    x      = seq(-1, 1, length.out = 5),
    y      = rep(.5, 5),
    method = "q_grid_cmde"
  ))

  samples <- .plot_brma_clear_posterior_density(
    samples          = samples,
    sample_parameter = "mu"
  )

  expect_null(attr(samples[["mu"]], "posterior_density", exact = TRUE))
  expect_null(attr(samples[["mu"]], "posterior_densities", exact = TRUE))
})


test_that("plot.brma qCMDE supports marginalized random SDs", {

  fit_name <- "brma.mv_block_mvn_random"
  if (!fit_name %in% list_fits(validate = FALSE, active_only = TRUE)) {
    skip(paste0("Raw cached fit unavailable for the active profile: ", fit_name))
  }
  fit <- try(load_fit(fit_name, validate = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    skip(paste0("Raw cached fit unavailable: ", fit_name))
  }

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {

      captured <<- list(samples = samples, parameter = parameter)
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fit,
    parameter       = "mu",
    plot_type       = "ggplot",
    density_method  = "qCMDE",
    density_control = list(n_points = 20, samples = 20)
  )
  posterior_density <- attr(
    captured[["samples"]][[captured[["parameter"]]]],
    "posterior_density",
    exact = TRUE
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(posterior_density[["method"]], "q_grid_cmde")
  expect_true(all(is.finite(posterior_density[["y"]])))
})


# list cached fits lazily
skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)


factor_plot_density_methods <- function(captured) {

  plot_data_samples_factor <- get(
    ".plot_data_samples.factor",
    envir    = asNamespace("BayesTools"),
    inherits = FALSE
  )
  plot_data <- plot_data_samples_factor(
    samples                  = captured[["samples"]],
    parameter                = captured[["parameter"]],
    n_points                 = 20,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    density_method           = "precomputed"
  )
  density_names <- grep("^density[0-9]+$", names(plot_data), value = TRUE)

  return(vapply(plot_data[density_names], function(x) {
    method <- attr(x, "posterior_density_method", exact = TRUE)
    if (is.null(method)) {
      return(NA_character_)
    }
    return(method)
  }, character(1)))
}


expect_factor_precomputed_densities <- function(captured, density_method) {

  expected_method <- if (identical(density_method, "qCMDE")) {
    "q_grid_cmde"
  } else {
    "iwmde"
  }
  methods <- factor_plot_density_methods(captured)

  expect_gt(length(methods), 0L)
  expect_false(anyNA(methods))
  expect_true(all(methods == expected_method))
}


test_that("plot.brma uses parameter x-axis labels by default", {

  temp_fit <- fits[["bcg_meta-analysis"]]

  x_label <- function(plot) {

    return(plot$scales$get_scales("x")$name)
  }

  expect_identical(
    x_label(plot(temp_fit, "mu", plot_type = "ggplot")),
    "Effect Size"
  )
  expect_identical(
    x_label(plot(temp_fit, "tau", plot_type = "ggplot")),
    "Heterogeneity"
  )
  expect_identical(
    x_label(plot(temp_fit, "mu", plot_type = "ggplot", xlab = "Custom Label")),
    "Custom Label"
  )
})


test_that("plot.brma transforms coefficients and heterogeneity intercepts", {

  skip_if_missing_fits(c(
    "bcg_meta-regression",
    "bangertdrowns2004_location-scale"
  ))

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {

      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fits[["bcg_meta-regression"]],
    parameter = "ablat",
    component = "mods",
    transform = "EXP",
    plot_type = "ggplot"
  )
  transformation <- captured[["dots"]][["transformation"]]

  expect_s3_class(out, "mock_plot")
  expect_identical(captured[["parameter"]], "mu_ablat")
  expect_equal(transformation[["fun"]](c(0, log(2))), c(1, 2))
  expect_equal(transformation[["jac"]](c(0, log(2))), c(1, 2))
  expect_identical(
    captured[["dots"]][["par_name"]],
    "Effect Size: ablat (risk ratio)"
  )

  out <- plot(
    fits[["bangertdrowns2004_location-scale"]],
    parameter = "ni100",
    component = "scale",
    transform = "EXP",
    plot_type = "ggplot"
  )
  transformation <- captured[["dots"]][["transformation"]]

  expect_s3_class(out, "mock_plot")
  expect_identical(captured[["parameter"]], "log_tau_ni100")
  expect_equal(transformation[["fun"]](c(0, log(2))), c(1, 2))
  expect_equal(transformation[["jac"]](c(0, log(2))), c(1, 2))
  expect_true(captured[["dots"]][["transformation_settings"]])
  expect_identical(
    captured[["dots"]][["par_name"]],
    "Heterogeneity: ni100 (multiplicative scale)"
  )

  out <- plot(
    fits[["bangertdrowns2004_location-scale"]],
    parameter                 = "intercept",
    component                 = "scale",
    standardized_coefficients = TRUE,
    transform                 = "LOG",
    plot_type                 = "ggplot"
  )
  transformation <- captured[["dots"]][["transformation"]]

  expect_s3_class(out, "mock_plot")
  expect_identical(captured[["parameter"]], "log_tau_intercept")
  expect_equal(transformation[["fun"]](c(.5, 1, 2)), log(c(.5, 1, 2)))
  expect_equal(transformation[["jac"]](c(.5, 1, 2)), c(2, 1, .5))
  expect_identical(
    captured[["dots"]][["par_name"]],
    "Heterogeneity (log scale)"
  )
})


test_that("plot.brma limits parameter-specific transformations", {

  skip_if_missing_fits(c("bcg_meta-analysis", "bcg_meta-regression"))

  expect_error(
    plot(
      fits[["bcg_meta-regression"]],
      parameter      = "ablat",
      component      = "mods",
      output_measure = "OR"
    ),
    "output_measure"
  )
  expect_error(
    plot(
      fits[["bcg_meta-regression"]],
      parameter = "ablat",
      component = "mods",
      transform = "LOG"
    ),
    "positive heterogeneity intercepts"
  )
  expect_error(
    plot(fits[["bcg_meta-analysis"]], parameter = "tau", transform = "EXP"),
    "meta-regression coefficients"
  )
  expect_error(
    plot(fits[["bcg_meta-analysis"]], parameter = "mu", transform = "LOG"),
    "positive heterogeneity intercepts"
  )
  expect_error(
    plot(fits[["bcg_meta-analysis"]], parameter = "tau", transform = "LOG"),
    "positive heterogeneity intercepts"
  )
})


test_that("plot.brma component disambiguates shared location-scale terms", {

  skip_if_missing_fits("dat.lehmann2018_RoBMA_3lvl_mods_scale")

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  fit <- fits[["dat.lehmann2018_RoBMA_3lvl_mods_scale"]]

  expect_error(
    plot(fit, parameter = "Preregistered", plot_type = "ggplot"),
    "ambiguous"
  )

  out <- plot(
    fit,
    parameter = "Preregistered",
    component = "scale",
    plot_type = "ggplot"
  )
  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "log_tau_Preregistered")

  out <- plot(
    fit,
    parameter = "Preregistered",
    component = "mods",
    plot_type = "ggplot"
  )
  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu_Preregistered")

  out <- plot(
    fit,
    parameter = "Preregistered",
    component = "location",
    plot_type = "ggplot"
  )
  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu_Preregistered")
})


test_that("lines.brma forwards posterior overlays", {

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- lines(
    fits[["bcg_meta-analysis"]],
    parameter = "mu",
    plot_type = "ggplot",
    col = "green"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu")
  expect_true(captured[["dots"]][["add"]])
  expect_false(captured[["dots"]][["prior"]])
  expect_equal(captured[["dots"]][["plot_type"]], "ggplot")
  expect_equal(captured[["dots"]][["col"]], "green")

  expect_error(
    lines(fits[["bcg_meta-analysis"]], parameter = "mu", prior = TRUE),
    "posterior densities only"
  )
})

test_that("plot.brma forwards secondary probability-axis controls", {

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fits[["bcg_meta-analysis"]],
    parameter = "mu",
    plot_type = "ggplot",
    ylim      = c(0, 12.5),
    ylim2     = c(0, 1),
    ylab2     = "Probability mass"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["dots"]][["ylim"]], c(0, 12.5))
  expect_equal(captured[["dots"]][["ylim2"]], c(0, 1))
  expect_identical(captured[["dots"]][["ylab2"]], "Probability mass")
})


test_that("plot.brma uses KDE by default", {

  captured <- NULL
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fits[["bcg_meta-analysis"]],
    parameter = "mu",
    plot_type = "ggplot"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["dots"]][["density_method"]], "KDE")
  expect_null(attr(captured[["samples"]][["mu"]], "posterior_density"))
})


test_that("plot.brma forwards attached qCMDE posterior density", {

  captured <- NULL
  .local_mock_iwmde_estimate_success()
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fits[["bcg_meta-analysis"]],
    parameter          = "mu",
    plot_type          = "ggplot",
    density_method     = "qCMDE",
    density_control    = list(n_points = 20, samples = 20)
  )

  posterior_density <- attr(captured[["samples"]][["mu"]], "posterior_density")

  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu")
  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
  expect_equal(posterior_density[["status"]], "ok")
  expect_equal(posterior_density[["density_method"]], "qCMDE")
  expect_true(all(is.finite(posterior_density[["x"]])))
  expect_true(all(is.finite(posterior_density[["y"]])))
})


test_that("plot.brma forwards attached IWMDE posterior density", {

  captured <- NULL
  .local_mock_iwmde_estimate_success()
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fits[["bcg_meta-analysis"]],
    parameter          = "mu",
    plot_type          = "ggplot",
    density_method     = "IWMDE",
    density_control    = list(n_points = 20, samples = 50)
  )

  posterior_density <- attr(captured[["samples"]][["mu"]], "posterior_density")

  expect_s3_class(out, "mock_plot")
  expect_equal(captured[["parameter"]], "mu")
  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
  expect_equal(posterior_density[["status"]], "ok")
  expect_equal(posterior_density[["density_method"]], "IWMDE")
  expect_equal(posterior_density[["diagnostics"]][["estimator"]], "iwmde")
})


test_that("plot.brma fails closed when an explicit estimator is rejected", {

  testthat::local_mocked_bindings(
    .iwmde_estimate = function(...) {

      return(list(
        diagnostics = list(density = list(
          status = "unsupported",
          reason = "sentinel diagnostic rejection"
        )),
        posterior_density = NULL
      ))
    },
    .package = "RoBMA"
  )

  expect_error(
    plot(
      fits[["bcg_meta-analysis"]],
      parameter       = "mu",
      plot_type       = "ggplot",
      density_method  = "qCMDE",
      density_control = list(n_points = 20, samples = 20)
    ),
    "density was unavailable: sentinel diagnostic rejection",
    fixed = TRUE
  )
})


test_that("plot.brma forwards qCMDE density on the fitted coefficient scale", {

  captured <- NULL
  .local_mock_iwmde_estimate_success()
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  plot(
    fits[["bcg_meta-regression"]],
    parameter_mods            = "year",
    standardized_coefficients = TRUE,
    plot_type                 = "ggplot",
    density_method            = "qCMDE",
    density_control           = list(n_points = 20, samples = 20)
  )

  plotted_samples   <- captured[["samples"]][[captured[["parameter"]]]]
  posterior_density <- attr(plotted_samples, "posterior_density")

  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
  expect_null(posterior_density[["diagnostics"]][["plot_scale_transform"]])
  expect_lte(
    max(abs(posterior_density[["x"]])),
    max(abs(as.numeric(plotted_samples))) + .05
  )
})


test_that("plot.brma uses exact original-scale coefficient targets", {

  captured                <- NULL
  captured_parameter_spec <- NULL
  estimate                <- .mock_iwmde_estimate_success()
  testthat::local_mocked_bindings(
    .iwmde_estimate = function(..., parameter_spec) {
      captured_parameter_spec <<- parameter_spec
      estimate(..., parameter_spec = parameter_spec)
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    plot_posterior = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    fits[["bcg_meta-regression"]],
    parameter_mods   = "year",
    plot_type        = "ggplot",
    density_method   = "qCMDE",
    density_control  = list(n_points = 20, samples = 20)
  )

  plotted_samples   <- captured[["samples"]][[captured[["parameter"]]]]
  posterior_density <- attr(plotted_samples, "posterior_density")

  expect_s3_class(out, "mock_plot")
  expect_identical(captured_parameter_spec[["type"]], "linear")
  expect_equal(
    range(posterior_density[["x"]]),
    range(as.numeric(plotted_samples)),
    tolerance = 1e-12
  )
  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
})


test_that("plot.brma forwards attached qCMDE/IWMDE densities for factor terms", {

  skip_if_missing_fits("bcg_meta-regression2")

  for (case in list(
    list(method = "qCMDE", samples = 20L, standardized = FALSE),
    list(method = "IWMDE", samples = 50L, standardized = FALSE),
    list(method = "qCMDE", samples = 20L, standardized = TRUE),
    list(method = "IWMDE", samples = 50L, standardized = TRUE)
  )) {
    captured <- NULL
    .local_mock_iwmde_estimate_success()
    testthat::local_mocked_bindings(
      plot_posterior = function(samples, parameter, ...) {
        captured <<- list(samples = samples, parameter = parameter, dots = list(...))
        return(structure(list(), class = "mock_plot"))
      },
      .package = "BayesTools"
    )

    expect_warning(
      out <- plot(
        fits[["bcg_meta-regression2"]],
        parameter_mods            = "alloc",
        plot_type                 = "ggplot",
        density_method            = case[["method"]],
        standardized_coefficients = case[["standardized"]],
        density_control           = list(
          n_points = 20,
          samples  = case[["samples"]]
        )
      ),
      NA
    )

    posterior_densities <- attr(
      captured[["samples"]][[captured[["parameter"]]]],
      "posterior_densities"
    )

    expect_s3_class(out, "mock_plot")
    expect_equal(captured[["dots"]][["density_method"]], "precomputed")
    expect_false("alternate" %in% names(posterior_densities))
    expect_true(all(vapply(posterior_densities, function(density) {
      identical(density[["status"]], "ok") &&
        identical(density[["density_method"]], case[["method"]])
    }, logical(1))))
    expect_factor_precomputed_densities(captured, case[["method"]])
  }
})


test_that("plot.brma forwards qCMDE/IWMDE densities for single-column factor terms", {

  skip_if_missing_fits("dat.lehmann2018_RoBMA_mods")

  for (case in list(
    list(method = "qCMDE", samples = 20L, conditional = FALSE),
    list(method = "IWMDE", samples = 50L, conditional = FALSE),
    list(method = "qCMDE", samples = 20L, conditional = TRUE),
    list(method = "IWMDE", samples = 50L, conditional = TRUE)
  )) {
    captured <- NULL
    .local_mock_iwmde_estimate_success()
    testthat::local_mocked_bindings(
      plot_posterior = function(samples, parameter, ...) {
        captured <<- list(samples = samples, parameter = parameter, dots = list(...))
        return(structure(list(), class = "mock_plot"))
      },
      .package = "BayesTools"
    )

    expect_warning(
      out <- plot(
        fits[["dat.lehmann2018_RoBMA_mods"]],
        parameter_mods  = "Preregistered",
        conditional     = case[["conditional"]],
        plot_type       = "ggplot",
        density_method  = case[["method"]],
        density_control = list(
          n_points = 20,
          samples  = case[["samples"]]
        )
      ),
      NA
    )

    sample              <- captured[["samples"]][[captured[["parameter"]]]]
    posterior_densities <- attr(sample, "posterior_densities")

    expect_s3_class(out, "mock_plot")
    expect_equal(captured[["dots"]][["density_method"]], "precomputed")
    expect_true(length(posterior_densities) > 0L)
    expect_true(all(vapply(posterior_densities, function(density) {
      identical(density[["status"]], "ok") &&
        identical(density[["density_method"]], case[["method"]])
    }, logical(1))))
    if (isTRUE(case[["conditional"]])) {
      condition_keys <- unique(vapply(posterior_densities, function(density) {
        key <- density[["condition_key"]]
        if (is.null(key)) {
          return(NA_character_)
        }
        return(as.character(key))
      }, character(1)))
      expect_true(attr(sample, "condition_key", exact = TRUE) %in% condition_keys)
    }
    expect_factor_precomputed_densities(captured, case[["method"]])
  }
})


test_that("plot.brma forwards fitted-scale densities under BayesTools aliases", {

  skip_if_missing_fits(c(
    "bcg_meta-regression3b",
    "bcg_meta-regression4",
    "dat.lehmann2018_RoBMA_mods2"
  ))

  cases <- data.frame(
    name = c(
      "bcg_meta-regression3b",
      "bcg_meta-regression4",
      "dat.lehmann2018_RoBMA_mods2"
    ),
    parameter_mods = c(
      "alloc:year",
      "alloc:year_before1969",
      "Preregistered:Gender"
    ),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(cases))) {
    captured <- NULL
    .local_mock_iwmde_estimate_success()
    testthat::local_mocked_bindings(
      plot_posterior = function(samples, parameter, ...) {
        captured <<- list(samples = samples, parameter = parameter, dots = list(...))
        return(structure(list(), class = "mock_plot"))
      },
      .package = "BayesTools"
    )

    expect_warning(
      out <- plot(
        fits[[cases[["name"]][[i]]]],
        parameter_mods           = cases[["parameter_mods"]][[i]],
        standardized_coefficients = TRUE,
        plot_type                = "ggplot",
        density_method           = "qCMDE",
        density_control          = list(n_points = 20, samples = 20)
      ),
      NA
    )

    expect_s3_class(out, "mock_plot")
    expect_equal(captured[["dots"]][["density_method"]], "precomputed")
    expect_factor_precomputed_densities(captured, "qCMDE")
  }
})


test_that("plot.brma drops qCMDE density when sample scales differ", {

  posterior_density <- list(
    x            = seq(-1, 1, length.out = 20),
    y            = stats::dnorm(seq(-1, 1, length.out = 20)),
    diagnostics  = list(),
    point_masses = data.frame(x = numeric(), mass = numeric())
  )

  out <- .plot_brma_align_iwmde_density(
    posterior_density = posterior_density,
    raw_samples       = seq(-1, 1, length.out = 20),
    plotted_samples   = seq(-1, 1, length.out = 20)^2
  )

  expect_null(out)
})


test_that("plot.brma forwards attached qCMDE density for PET and PEESE parameters", {

  .skip_if_missing_raw_fits(c("dat.lehmann2018-PET", "dat.lehmann2018-PEESE"))

  for (case in list(
    list(fit = "dat.lehmann2018-PET", parameter = "PET"),
    list(fit = "dat.lehmann2018-PEESE", parameter = "PEESE")
  )) {
    captured <- NULL
    .local_mock_iwmde_estimate_success()
    testthat::local_mocked_bindings(
      plot_posterior = function(samples, parameter, ...) {
        captured <<- list(samples = samples, parameter = parameter, dots = list(...))
        return(structure(list(), class = "mock_plot"))
      },
      .package = "BayesTools"
    )

    out <- plot(
      load_fit(case[["fit"]], validate = FALSE),
      parameter         = case[["parameter"]],
      plot_type         = "ggplot",
      density_method    = "qCMDE",
      density_control   = list(n_points = 20, samples = 20)
    )

    density_source <- if (!is.null(captured[["samples"]][[case[["parameter"]]]])) {
      case[["parameter"]]
    } else {
      "bias"
    }
    posterior_density <- attr(captured[["samples"]][[density_source]], "posterior_density")

    expect_s3_class(out, "mock_plot")
    expect_equal(captured[["dots"]][["density_method"]], "precomputed")
    expect_equal(posterior_density[["status"]], "ok")
    expect_equal(posterior_density[["density_method"]], "qCMDE")
  }
})


test_that("Prior and posterior distributions for brma.norm models", {

  ### simple meta-analysis ----
  name <- "bcg_meta-analysis"
  temp_fit <- fits[[name]]

  ### effect size
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_no_prior"), function() plot(temp_fit, "mu", prior = FALSE))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_no_prior"), plot(temp_fit, "mu", prior = FALSE, plot_type = "ggplot"))
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_prior"), function() plot(temp_fit, "mu", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_prior"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot"))

  # change range
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_range"), function() plot(temp_fit, "mu", prior = TRUE, xlim = c(-1, 1), ylim = c(0, 5)))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_range"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot", xlim = c(-1, 1), ylim = c(0, 5)))

  # change aesthetics
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_aesthetics"), function() plot(temp_fit, "mu", prior = TRUE, lwd = 3, lty = 3, col = "blue", dots_prior = list(lwd = 3, lty = 1, col = "red")))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_aesthetics"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", dots_prior = list(lwd = 3, lty = 1, col = "red")))

  ### heterogeneity
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot_tau_pp_prior"), function() plot(temp_fit, "tau", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-tau_ggplot_pp_prior"), plot(temp_fit, "tau", prior = TRUE, plot_type = "ggplot"))

  ### meta-regression (continuous)
  name <- "bcg_meta-regression"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_orig"), function() plot(temp_fit, "mu", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_std"), function() plot(temp_fit, "mu", prior = TRUE, standardized_coefficients = TRUE))

  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_orig"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_std"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression (categorical: dummy)
  name <- "bcg_meta-regression2"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression (categorical: meandif)
  name <- "bcg_meta-regression2b"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression with interactions (categorical: dummy)
  name <- "bcg_meta-regression3"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-reg0_baseplot_orig"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg0_baseplot_std"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg1_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg1_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg2_baseplot_orig"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg2_baseplot_std"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg3_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg3_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE, standardized_coefficients = TRUE))

  ### meta-regression with interactions (categorical: meandif)
  name <- "bcg_meta-regression3b"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-reg0_baseplot_orig"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg0_baseplot_std"), function() plot(temp_fit, parameter_mods = "intercept", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg1_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg1_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg2_baseplot_orig"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg2_baseplot_std"), function() plot(temp_fit, parameter_mods = "year", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg3_baseplot_orig"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-reg3_baseplot_std"), function() plot(temp_fit, parameter_mods = "alloc:year", prior = TRUE, standardized_coefficients = TRUE))

  ### location-scale model
  name <- "bangertdrowns2004_location-scale"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-scale0_baseplot_orig"), function() plot(temp_fit, parameter_scale = "intercept", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-scale0_baseplot_std"), function() plot(temp_fit, parameter_scale = "intercept", prior = TRUE, standardized_coefficients = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-scale1_baseplot_orig"), function() plot(temp_fit, parameter_scale = "ni100", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-scale1_baseplot_std"), function() plot(temp_fit, parameter_scale = "ni100", prior = TRUE, standardized_coefficients = TRUE))

  ### between-study heterogeneity and multilevel models
  name <- "konstantopoulos2011_3lvl"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-rho_baseplot"), function() plot(temp_fit, parameter = "rho", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot"), function() plot(temp_fit, parameter = "tau", prior = TRUE))
})

test_that("Transformed effect-size plots render", {

  name     <- "bcg_meta-analysis"
  temp_fit <- fits[[name]]

  .with_temp_plot_device(
    expect_silent(plot(temp_fit, "mu", transform = "EXP", plot_type = "base"))
  )
  expect_true(.is_ggplot(
    plot(temp_fit, "mu", transform = "EXP", plot_type = "ggplot")
  ))
})

test_that("Prior and posterior plots transform effect-size axis", {

  skip_if_not_full_visuals("Effect-size transform snapshots are visual-gallery coverage.")

  name     <- "bcg_meta-analysis"
  temp_fit <- fits[[name]]

  expect_vdiffr_snapshot(paste0(name, "-mu_transform_exp_comparison"), function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]], mar = oldpar[["mar"]]))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    plot(temp_fit, "mu", plot_type = "base", main = "log RR")
    plot(temp_fit, "mu", transform = "EXP", plot_type = "base", main = "RR")
  })
})

test_that("Prior and posterior distributions for bPET / bPEESE objects", {

  skip_if_not_full_visuals("PET/PEESE prior-posterior snapshots are visual-gallery coverage.")

  ### PET model
  name <- "dat.lehmann2018-PET"
  temp_fit <- fits[[name]]
  set.seed(1)

  # no prior
  expect_vdiffr_snapshot("baseplot_pp_PETPEESE_no_prior", function() plot_pet_peese(temp_fit, prior = FALSE))
  expect_vdiffr_snapshot("ggplot_pp_PETPEESE_no_prior", plot_pet_peese(temp_fit, prior = FALSE, plot_type = "ggplot"))

  # change range
  expect_vdiffr_snapshot("baseplot_pp_PETPEESE_range", function() plot_pet_peese(temp_fit, prior = TRUE, xlim = c(0, 0.5), ylim = c(-3, 3)))
  expect_vdiffr_snapshot("ggplot_pp_PETPEESE_range", plot_pet_peese(temp_fit, prior = TRUE, plot_type = "ggplot", xlim = c(0, 0.5), ylim = c(-3, 3)))

  # change aesthetics
  expect_vdiffr_snapshot("baseplot_pp_PETPEESE_aesthetics", function() plot_pet_peese(temp_fit, prior = TRUE, lwd = 3, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                          dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))
  expect_vdiffr_snapshot("ggplot_pp_PETPEESE_aesthetics", plot_pet_peese(temp_fit, prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                             dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))

})

test_that("Prior and posterior distributions for bselmodel objects", {

  skip_if_not_full_visuals("Selection prior-posterior snapshots are visual-gallery coverage.")


  ### weight function
  name <- "dat.lehmann2018-3PSM"
  temp_fit <- fits[[name]]
  set.seed(1)

  # no prior
  expect_vdiffr_snapshot("baseplot_pp_weightfunction_no_prior", function() plot_weightfunction(temp_fit, prior = FALSE))
  expect_vdiffr_snapshot("ggplot_pp_weightfunction_no_prior", plot_weightfunction(temp_fit, prior = FALSE, plot_type = "ggplot"))

  # change range
  expect_vdiffr_snapshot("baseplot_pp_weightfunction_range", function() plot_weightfunction(temp_fit, prior = TRUE, rescale_p_values = FALSE))
  expect_vdiffr_snapshot("ggplot_pp_weightfunction_range", plot_weightfunction(temp_fit, prior = TRUE, plot_type = "ggplot", rescale_p_values = FALSE))

  # change aesthetics
  expect_vdiffr_snapshot("baseplot_pp_weightfunction_aesthetics", function() plot_weightfunction(temp_fit, prior = TRUE, lwd = 3, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                                      dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))
  expect_vdiffr_snapshot("ggplot_pp_weightfunction_aesthetics", plot_weightfunction(temp_fit, prior = TRUE, plot_type = "ggplot", lwd = 2, lty = 3, col = "blue", col.fill = scales::alpha("blue", 0.20),
                                                                                         dots_prior = list(lwd = 3, lty = 1, col = "red", col.fill = scales::alpha("red", 0.20))))

})

test_that("Weightfunction plot supports observed p-value rug", {

  name <- "dat.lehmann2018-3PSM"
  skip_if_missing_fits(name)

  .with_temp_plot_device(
    expect_silent(plot_weightfunction(fits[[name]], show_data = TRUE))
  )
  expect_true(.is_ggplot(
    plot_weightfunction(
      fits[[name]],
      show_data = TRUE,
      plot_type = "ggplot",
      dots_data = list(color = "red", alpha = .5, linewidth = .4, rug_side = "top")
    )
  ))
  expect_silent(
    plot_weightfunction(
      fits[[name]],
      show_data        = TRUE,
      rescale_p_values = FALSE,
      dots_data        = list(col = "blue", lwd = 1, side = "bottom")
    )
  )

  rug_plot <- plot_weightfunction(
    fits[[name]],
    show_data        = TRUE,
    plot_type        = "ggplot",
    rescale_p_values = FALSE
  )
  rug_data <- ggplot_geom_layer_data(rug_plot, "GeomRug")

  expect_equal(
    rug_data[["p"]],
    .weightfunction_observed_p_values(fits[[name]]),
    tolerance = sqrt(.Machine$double.eps)
  )
})

test_that("Prior and posterior distributions for BMA.norm objects", {

  skip_if_not_full_visuals("BMA prior-posterior snapshots are visual-gallery coverage.")

  name <- "dat.lehmann2018_BMA.norm_mods"
  temp_fit <- fits[[name]]
  set.seed(1)

  # effect
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_no_prior"), function() plot(temp_fit, "mu", prior = FALSE))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_no_prior"), plot(temp_fit, "mu", prior = FALSE, plot_type = "ggplot"))
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_prior"), function() plot(temp_fit, "mu", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_prior"), plot(temp_fit, "mu", prior = TRUE, plot_type = "ggplot"))

  # moderation
  expect_vdiffr_snapshot(paste0(name, "-mods_ggplot_pp_no_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = FALSE, plot_type = "ggplot"))
  expect_vdiffr_snapshot(paste0(name, "-mods_ggplot_pp_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = TRUE, plot_type = "ggplot"))

  # heterogeneity
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot_pp_no_prior"), function() plot(temp_fit, "tau", prior = FALSE))
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot_pp_prior"), function() plot(temp_fit, "tau", prior = TRUE))

})

test_that("Prior and posterior distributions for RoBMA objects", {

  skip_if_not_full_visuals("RoBMA prior-posterior snapshots are visual-gallery coverage.")

  name <- "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  temp_fit <- fits[[name]]
  set.seed(1)

  # effect
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_no_prior"), function() plot(temp_fit, "mu", prior = FALSE))
  expect_vdiffr_snapshot(paste0(name, "-mu_baseplot_pp_prior"), function() plot(temp_fit, "mu", prior = TRUE))

  # moderation
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_no_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = FALSE, plot_type = "ggplot"))
  expect_vdiffr_snapshot(paste0(name, "-mu_ggplot_pp_prior"), plot(temp_fit, parameter_mods = "Preregistered", prior = TRUE, plot_type = "ggplot"))

  # heterogeneity
  expect_vdiffr_snapshot(paste0(name, "-tau_baseplot_pp_no_prior"), function() plot(temp_fit, "tau", prior = TRUE))
  expect_vdiffr_snapshot(paste0(name, "-rho_baseplot_pp_no_prior"), function() plot(temp_fit, "rho", prior = TRUE))

})
