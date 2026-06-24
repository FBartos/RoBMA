context("Estimated marginal means")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-visuals.R"))
REFERENCE_DIR <<- testthat::test_path("..", "results", "marginal_means")

.expect_marginal_means_bf_values <- function(actual, expected) {

  expect_equal(as.numeric(unlist(actual)), as.numeric(unlist(expected)))
  for (name in names(expected)) {
    expect_equal(
      attr(actual[[name]], "BF_error_percent", exact = TRUE),
      attr(expected[[name]], "BF_error_percent", exact = TRUE)
    )
    expect_equal(
      attr(actual[[name]], "posterior_density_source", exact = TRUE),
      attr(expected[[name]], "posterior_density_source", exact = TRUE)
    )
  }
}


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


.marginal_means_iwmde_test_object <- function() {

  emm <- .marginal_means_test_object()
  for (type in c("averaged", "conditional")) {
    for (level in names(emm[["inference"]][[type]][["mu_alloc"]])) {
      sample <- emm[["inference"]][[type]][["mu_alloc"]][[level]]
      attr(sample, "linear_weights") <- c(mu_intercept = 1, mu_alloc = 1)
      if (identical(type, "conditional")) {
        attr(sample, "effective_conditional") <- c("mu_intercept", "mu_alloc")
        attr(sample, "effective_conditional_rule") <- "OR"
        attr(sample, "condition_key") <- "OR\r2\rmu_alloc\rmu_intercept"
      }
      emm[["inference"]][[type]][["mu_alloc"]][[level]] <- sample
    }
  }

  return(emm)
}


.marginal_means_iwmde_density_diagnostics <- function(estimator = "q_grid_cmde") {

  list(
    estimator                         = estimator,
    max_relative_mcse                 = .01,
    min_finite_terms                  = 100L,
    min_ess                           = 100,
    max_weight_share                  = .10,
    row_drop_fraction                 = 0,
    active_mass                       = 1,
    normalization_integral            = 1,
    normalization_mass_ratio          = 1,
    max_ordinate_relative_change      = 0,
    max_normalizer_relative_change    = 0,
    max_quadrature_relative_change    = 0
  )
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


test_that("marginal_means BF ordinate targets stay conditional", {

  emm <- .marginal_means_iwmde_test_object()

  density_specs <- .marginal_means_iwmde_specs_by_type(
    marginal_means_object = emm,
    parameter             = "alloc",
    type                  = "averaged",
    levels                = "alternate",
    targeted              = TRUE
  )
  ordinate_specs <- .marginal_means_iwmde_ordinate_specs(
    marginal_means_object = emm,
    parameter             = "alloc",
    levels                = "alternate",
    targeted              = TRUE
  )

  expect_equal(names(density_specs), "averaged")
  expect_equal(names(density_specs[["averaged"]]), "alloc: alternate")
  expect_equal(density_specs[["averaged"]][["alloc: alternate"]][["marginal_type"]],
               "averaged")
  expect_equal(names(ordinate_specs), "alloc: alternate")
  expect_equal(ordinate_specs[["alloc: alternate"]][["marginal_type"]],
               "conditional")
})


test_that("marginal_means attaches BF ordinates for averaged density targets", {

  emm <- .marginal_means_iwmde_test_object()
  estimate_calls <- list()
  refreshed <- NULL
  testthat::local_mocked_bindings(
    .iwmde_context = function(object) list(source = object),
    .iwmde_estimate = function(context, parameter, density_method,
                               density_control, outputs, values,
                               parameter_spec, metadata, cache) {

      estimate_calls[[length(estimate_calls) + 1L]] <<- list(
        outputs       = outputs,
        marginal_type = parameter_spec[["marginal_type"]]
      )
      diagnostic <- list(status = "ok")
      if (identical(outputs, "density")) {
        return(list(
          diagnostics       = list(density = diagnostic),
          posterior_density = list(
            density_method = density_method,
            marginal_type  = parameter_spec[["marginal_type"]],
            diagnostics    = .marginal_means_iwmde_density_diagnostics(
              .density_method_iwmde_estimator(density_method)
            )
          )
        ))
      }

      return(list(
        diagnostics        = list(ordinate = diagnostic),
        posterior_ordinate = list(
          density_method = density_method,
          value          = values,
          marginal_type  = parameter_spec[["marginal_type"]]
        )
      ))
    },
    .marginal_means_refresh_iwmde_bf = function(inference, parameters,
                                                null_hypothesis,
                                                density_method = "KDE",
                                                object = NULL) {

      refreshed <<- inference
      return(inference)
    },
    .package = "RoBMA"
  )

  out <- .marginal_means_attach_iwmde(
    object                = structure(list(), class = "brma"),
    marginal_means_object = emm,
    n_points              = 20,
    max_samples           = 20,
    normalization_points  = 20,
    normalization_prob    = .99,
    density_method        = "qCMDE",
    display_grid          = "support",
    null_hypothesis       = 0,
    parameter             = "alloc",
    type                  = "averaged",
    levels                = "alternate",
    targeted              = TRUE
  )

  expect_equal(
    vapply(estimate_calls, `[[`, character(1), "outputs"),
    c("density", "ordinate")
  )
  expect_equal(
    vapply(estimate_calls, `[[`, character(1), "marginal_type"),
    c("averaged", "conditional")
  )
  expect_equal(
    attr(
      out[["inference"]][["averaged"]][["mu_alloc"]][["alternate"]],
      "posterior_density"
    )[["marginal_type"]],
    "averaged"
  )
  expect_null(attr(
    out[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  ))
  expect_equal(
    attr(
      out[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
      "posterior_ordinate"
    )[["marginal_type"]],
    "conditional"
  )
  expect_equal(
    attr(
      refreshed[["conditional"]][["mu_alloc"]][["alternate"]],
      "posterior_ordinate"
    )[["marginal_type"]],
    "conditional"
  )
})


test_that("marginal_means plot defaults to KDE despite stored posterior density", {

  emm <- .marginal_means_test_object()
  density <- list(
    status         = "ok",
    x              = seq(-1, 1, length.out = 5),
    y              = rep(.5, 5),
    density_method = "qCMDE",
    method         = "q_grid_cmde",
    diagnostics    = .marginal_means_iwmde_density_diagnostics("q_grid_cmde"),
    point_masses   = data.frame(x = numeric(), mass = numeric())
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
  expect_equal(captured[["dots"]][["density_method"]], "KDE")
  expect_identical(
    attr(captured[["samples"]][["mu_alloc"]][["alternate"]], "posterior_density"),
    density
  )
})


test_that("marginal_means plot computes missing explicit qCMDE densities", {

  emm <- .marginal_means_test_object()
  emm[["source_object"]] <- structure(list(fit = TRUE), class = "brma")
  levels <- names(emm[["inference"]][["averaged"]][["mu_alloc"]])
  provenance <- stats::setNames(lapply(levels, function(level) {
    list(
      request_key    = paste0("qCMDE-", level),
      density_method = "qCMDE",
      method         = "q_grid_cmde"
    )
  }), levels)
  density <- list(
    status         = "ok",
    x              = seq(-1, 1, length.out = 5),
    y              = rep(.5, 5),
    density_method = "qCMDE",
    method         = "q_grid_cmde",
    diagnostics    = .marginal_means_iwmde_density_diagnostics("q_grid_cmde"),
    point_masses   = data.frame(x = numeric(), mass = numeric())
  )
  attr(
    emm[["inference"]][["averaged"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  ) <- density
  density_for <- function(level) {
    density[["iwmde_provenance"]] <- provenance[[level]]
    return(density)
  }

  attached <- NULL
  captured <- NULL
  testthat::local_mocked_bindings(
    .iwmde_marginal_means_source_object = function(marginal_means_object) {
      marginal_means_object[["source_object"]]
    },
    .iwmde_context = function(object) {
      list(object = object)
    },
    .marginal_means_density_request_provenance = function(...) {
      provenance
    },
    .marginal_means_attach_iwmde = function(object, marginal_means_object,
                                            n_points, max_samples,
                                            normalization_points,
                                            normalization_prob, density_method,
                                            display_grid, null_hypothesis,
                                            parameter, type, levels, targeted,
                                            include_ordinates = TRUE) {

      attached <<- list(
        density_method    = density_method,
        parameter         = parameter,
        type              = type,
        levels            = levels,
        include_ordinates = include_ordinates
      )
      for (level in levels) {
        attr(
          marginal_means_object[["inference"]][[type]][["mu_alloc"]][[level]],
          "posterior_density"
        ) <- density_for(level)
      }

      return(marginal_means_object)
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    plot_marginal = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- expect_message(
    plot(
      emm,
      parameter       = "alloc",
      plot_type       = "ggplot",
      density_method  = "qCMDE"
    ),
    "Computing qCMDE density for 3 marginal-means levels"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(attached[["density_method"]], "qCMDE")
  expect_equal(attached[["parameter"]], "alloc")
  expect_equal(attached[["type"]], "averaged")
  expect_equal(attached[["levels"]], levels)
  expect_false(attached[["include_ordinates"]])
  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
  expect_true(all(vapply(
    levels,
    function(level) {
      posterior_density <- attr(
        captured[["samples"]][["mu_alloc"]][[level]],
        "posterior_density",
        exact = TRUE
      )
      identical(
        posterior_density[["iwmde_provenance"]][["request_key"]],
        provenance[[level]][["request_key"]]
      )
    },
    logical(1)
  )))
})


test_that("marginal_means plot density coverage requires provenance", {

  emm <- .marginal_means_test_object()
  samples <- emm[["inference"]][["averaged"]][["mu_alloc"]]
  levels <- names(samples)
  provenance <- stats::setNames(lapply(levels, function(level) {
    list(
      request_key    = paste0("qCMDE-", level),
      density_method = "qCMDE",
      method         = "q_grid_cmde"
    )
  }), levels)
  qcmde_density <- list(
    x              = seq(-1, 1, length.out = 5),
    y              = rep(.5, 5),
    density_method = "qCMDE",
    method         = "q_grid_cmde",
    diagnostics    = .marginal_means_iwmde_density_diagnostics("q_grid_cmde")
  )
  for (level in levels) {
    density <- qcmde_density
    density[["iwmde_provenance"]] <- provenance[[level]]
    attr(samples[[level]], "posterior_density") <- density
  }

  expect_equal(
    .marginal_means_missing_posterior_density_levels(samples, provenance),
    character()
  )
  attr(samples[["random"]], "posterior_density") <- qcmde_density
  expect_equal(
    .marginal_means_missing_posterior_density_levels(samples, provenance),
    "random"
  )
  stale_density <- qcmde_density
  stale_density[["iwmde_provenance"]] <- list(
    request_key    = "stale",
    density_method = "qCMDE",
    method         = "q_grid_cmde"
  )
  attr(samples[["systematic"]], "posterior_density") <- stale_density
  expect_equal(
    .marginal_means_missing_posterior_density_levels(samples, provenance),
    c("random", "systematic")
  )
})


test_that("marginal_means BF refresh requires ordinate provenance", {

  diagnostics <- list(
    estimator                    = "q_grid_cmde",
    ordinate                     = 1,
    relative_mcse                = .01,
    finite_terms                 = 100,
    ess                          = 100,
    max_weight_share             = .10,
    row_drop_fraction            = 0,
    active_mass                  = 1,
    normalization_integral       = 1,
    normalization_mass_ratio     = 1,
    ordinate_relative_change     = 0,
    max_normalizer_relative_change = 0
  )
  provenance <- list(
    request_key    = "expected",
    density_method = "qCMDE",
    method         = "q_grid_cmde"
  )
  stale_provenance <- provenance
  stale_provenance[["request_key"]] <- "stale"
  ordinate_for <- function(provenance) {
    list(
      value            = 0,
      ordinate         = 1,
      evaluation_value = 0,
      method           = "q_grid_cmde",
      density_method   = "qCMDE",
      diagnostics      = diagnostics,
      iwmde_provenance = provenance
    )
  }
  posterior <- list(
    alternate  = stats::rnorm(50),
    random     = stats::rnorm(50),
    systematic = stats::rnorm(50)
  )
  attr(posterior[["alternate"]], "posterior_ordinate") <- ordinate_for(provenance)
  attr(posterior[["random"]], "posterior_ordinate")    <- ordinate_for(stale_provenance)

  testthat::local_mocked_bindings(
    Savage_Dickey_BF = function(posterior, null_hypothesis,
                                normal_approximation, silent,
                                density_method) {

      out <- as.list(seq_along(posterior))
      names(out) <- names(posterior)
      return(out)
    },
    .package = "BayesTools"
  )

  bf <- .marginal_means_iwmde_bf(
    posterior       = posterior,
    null_hypothesis = 0,
    provenance      = list(
      alternate  = provenance,
      random     = provenance,
      systematic = provenance
    )
  )

  expect_equal(as.numeric(bf[["alternate"]]), 1)
  expect_true(is.na(bf[["random"]]))
  expect_true(is.na(bf[["systematic"]]))
  expect_match(
    attr(bf[["random"]], "warnings"),
    "posterior ordinate was unavailable"
  )
})


test_that("marginal_means plot does not reuse qCMDE density for explicit IWMDE", {

  emm <- .marginal_means_test_object()
  emm[["source_object"]] <- structure(list(fit = TRUE), class = "brma")
  levels <- names(emm[["inference"]][["averaged"]][["mu_alloc"]])
  provenance <- stats::setNames(lapply(levels, function(level) {
    list(
      request_key    = paste0("IWMDE-", level),
      density_method = "IWMDE",
      method         = "iwmde"
    )
  }), levels)
  qcmde_density <- list(
    status         = "ok",
    x              = seq(-1, 1, length.out = 5),
    y              = rep(.5, 5),
    density_method = "qCMDE",
    method         = "q_grid_cmde",
    diagnostics    = .marginal_means_iwmde_density_diagnostics("q_grid_cmde"),
    point_masses   = data.frame(x = numeric(), mass = numeric())
  )
  iwmde_density <- qcmde_density
  iwmde_density[["density_method"]] <- "IWMDE"
  iwmde_density[["method"]]         <- "iwmde"
  iwmde_density[["diagnostics"]]    <-
    .marginal_means_iwmde_density_diagnostics("iwmde")
  for (level in levels) {
    attr(
      emm[["inference"]][["averaged"]][["mu_alloc"]][[level]],
      "posterior_density"
    ) <- qcmde_density
  }

  attached <- NULL
  captured <- NULL
  testthat::local_mocked_bindings(
    .iwmde_marginal_means_source_object = function(marginal_means_object) {
      marginal_means_object[["source_object"]]
    },
    .iwmde_context = function(object) {
      list(object = object)
    },
    .marginal_means_density_request_provenance = function(...) {
      provenance
    },
    .marginal_means_attach_iwmde = function(object, marginal_means_object,
                                            n_points, max_samples,
                                            normalization_points,
                                            normalization_prob, density_method,
                                            display_grid, null_hypothesis,
                                            parameter, type, levels, targeted,
                                            include_ordinates = TRUE) {

      attached <<- list(density_method = density_method, levels = levels)
      for (level in levels) {
        density <- iwmde_density
        density[["iwmde_provenance"]] <- provenance[[level]]
        attr(
          marginal_means_object[["inference"]][[type]][["mu_alloc"]][[level]],
          "posterior_density"
        ) <- density
      }

      return(marginal_means_object)
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    plot_marginal = function(samples, parameter, ...) {
      captured <<- list(samples = samples, parameter = parameter, dots = list(...))
      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- expect_message(
    plot(
      emm,
      parameter       = "alloc",
      plot_type       = "ggplot",
      density_method  = "IWMDE",
      density_control = list(n_points = 20, max_samples = 20)
    ),
    "because 'density_control' was supplied"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(attached[["density_method"]], "IWMDE")
  expect_equal(
    attached[["levels"]],
    levels
  )
  expect_equal(captured[["dots"]][["density_method"]], "precomputed")
  expect_true(all(vapply(
    captured[["samples"]][["mu_alloc"]],
    function(sample) {
      identical(
        attr(sample, "posterior_density", exact = TRUE)[["density_method"]],
        "IWMDE"
      )
    },
    logical(1)
  )))
})


test_that("marginal_means plot errors when explicit IWMDE density is unavailable", {

  emm <- .marginal_means_test_object()
  emm[["source_object"]] <- structure(list(fit = TRUE), class = "brma")
  levels <- names(emm[["inference"]][["averaged"]][["mu_alloc"]])
  provenance <- stats::setNames(lapply(levels, function(level) {
    list(
      request_key    = paste0("IWMDE-", level),
      density_method = "IWMDE",
      method         = "iwmde"
    )
  }), levels)
  testthat::local_mocked_bindings(
    .iwmde_marginal_means_source_object = function(marginal_means_object) {
      marginal_means_object[["source_object"]]
    },
    .iwmde_context = function(object) {
      list(object = object)
    },
    .marginal_means_density_request_provenance = function(...) {
      provenance
    },
    .marginal_means_attach_iwmde = function(object, marginal_means_object,
                                            n_points, max_samples,
                                            normalization_points,
                                            normalization_prob, density_method,
                                            display_grid, null_hypothesis,
                                            parameter, type, levels, targeted,
                                            include_ordinates = TRUE) {

      return(marginal_means_object)
    },
    .package = "RoBMA"
  )

  expect_message(
    expect_error(
      plot(
        emm,
        parameter       = "alloc",
        plot_type       = "ggplot",
        density_method  = "IWMDE",
        density_control = list(n_points = 20, max_samples = 20)
      ),
      "IWMDE density was unavailable"
    ),
    "Computing IWMDE density"
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
  expect_false("normal_approximation" %in% names(formals(marginal_means.brma)))
  expect_equal(mm[["density_method"]], "KDE")
  expect_identical(mm[["source_object"]], fits[["bcg_meta-regression2"]])
})


test_that("marginal_means hides normal approximation density method", {

  captured <- NULL
  testthat::local_mocked_bindings(
    as_marginal_inference = function(model, marginal_parameters, parameters,
                                     conditional_list, conditional_rule,
                                     formula, null_hypothesis,
                                     normal_approximation, n_samples, silent,
                                     force_plots) {

      captured <<- list(
        normal_approximation = normal_approximation,
        parameters           = parameters
      )
      samples <- stats::setNames(
        lapply(parameters, function(parameter) {
          structure(
            list(level = stats::rnorm(20)),
            class = c("marginal_posterior.factor", "marginal_posterior")
          )
        }),
        parameters
      )
      inference <- list(
        averaged    = samples,
        conditional = samples,
        inference   = stats::setNames(vector("list", length(parameters)), parameters)
      )
      class(inference) <- c("marginal_inference", "list")

      return(inference)
    },
    .package = "BayesTools"
  )

  expect_error(
    marginal_means(
      fits[["bcg_meta-regression2"]],
      density_method = "normal"
    ),
    "must be one of"
  )

  expect_null(captured)
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
    density_method       = mm[["density_method"]]
  )
  expected <- BayesTools::Savage_Dickey_BF(
    posterior            = mm[["inference"]][["conditional"]][["mu_alloc"]],
    null_hypothesis      = mm[["null_hypothesis"]],
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )
  hypothesis_bf <- BayesTools::hypothesis_BF(
    posterior      = mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    hypothesis     = paste0("mu_alloc = ", mm[["null_hypothesis"]]),
    parameter      = "mu_alloc",
    columns        = "all",
    density_method = "precomputed"
  )
  level_bf <- BayesTools::hypothesis_BF(
    posterior  = mm[["inference"]][["averaged"]][["mu_alloc"]],
    hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
    parameter  = "mu_alloc",
    columns    = "all",
    seed       = 11
  )

  .expect_marginal_means_bf_values(
    actual   = refreshed[["inference"]][["mu_alloc"]],
    expected = expected
  )
  expect_true(is.finite(attr(expected[["alternate"]], "BF_error_percent")))
  expect_equal(attr(hypothesis_bf, "raw_BF"), as.numeric(expected[["alternate"]]),
               tolerance = 1e-12)
  expect_equal(as.numeric(hypothesis_bf[["BF_error"]]),
               attr(expected[["alternate"]], "BF_error_percent"),
               tolerance = 1e-12)
  expect_equal(hypothesis_bf[["method"]], "Savage-Dickey (precomputed)")
  expect_equal(level_bf[["method"]], "prior-posterior odds")
  expect_true(is.finite(attr(level_bf, "raw_BF")))
})


test_that("marginal_means restricts qCMDE precomputation targets", {

  mm <- marginal_means(
    fits[["bcg_meta-regression2"]],
    n_samples       = 1000,
    bf              = TRUE,
    density_method  = "qCMDE",
    parameter       = "alloc",
    type            = "conditional",
    levels          = "alternate",
    density_control = list(n_points = 20, max_samples = 20)
  )

  averaged_density <- attr(
    mm[["inference"]][["averaged"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  )
  conditional_density <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  )
  skipped_density <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["random"]],
    "posterior_density"
  )
  skipped_ordinate <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["random"]],
    "posterior_ordinate"
  )

  expect_null(averaged_density)
  expect_equal(conditional_density[["density_method"]], "qCMDE")
  expect_null(skipped_density)
  expect_null(skipped_ordinate)
  expect_equal(mm[["iwmde_settings"]][["parameter"]], "alloc")
  expect_equal(mm[["iwmde_settings"]][["type"]], "conditional")
  expect_equal(mm[["iwmde_settings"]][["levels"]], "alternate")
  expect_true(is.finite(as.numeric(
    mm[["inference"]][["inference"]][["mu_alloc"]][["alternate"]]
  )))
  expect_true(is.na(mm[["inference"]][["inference"]][["mu_alloc"]][["random"]]))
  expect_true(is.na(mm[["inference"]][["inference"]][["mu_alloc"]][["systematic"]]))
  expect_false(is.null(attr(
    mm[["inference"]][["inference"]][["mu_alloc"]][["random"]],
    "warnings"
  )))
})


test_that("marginal_means computes BF ordinates when density target is averaged", {

  mm <- marginal_means(
    fits[["bcg_meta-regression2"]],
    n_samples       = 1000,
    bf              = TRUE,
    density_method  = "qCMDE",
    parameter       = "alloc",
    type            = "averaged",
    levels          = "alternate",
    density_control = list(n_points = 20, max_samples = 20)
  )

  averaged_density <- attr(
    mm[["inference"]][["averaged"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  )
  conditional_density <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_density"
  )
  conditional_ordinate <- attr(
    mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    "posterior_ordinate"
  )

  expect_equal(averaged_density[["density_method"]], "qCMDE")
  expect_null(conditional_density)
  expect_equal(conditional_ordinate[["density_method"]], "qCMDE")
  expect_equal(mm[["iwmde_settings"]][["type"]], "averaged")
  expect_true(.marginal_means_has_bf_posterior_ordinate(
    mm[["inference"]][["conditional"]][["mu_alloc"]]
  ))
  expect_true(is.finite(as.numeric(
    mm[["inference"]][["inference"]][["mu_alloc"]][["alternate"]]
  )))
  expect_true(is.na(mm[["inference"]][["inference"]][["mu_alloc"]][["random"]]))
  expect_true(is.na(mm[["inference"]][["inference"]][["mu_alloc"]][["systematic"]]))
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
  expect_true(max(ordinate_diagnostic[["diagnostics"]][["normalization_range"]]) <
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
  .expect_marginal_means_bf_values(
    actual   = mm[["inference"]][["inference"]][["mu_alloc"]],
    expected = expected
  )
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
    density_method       = mm[["density_method"]]
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

  posterior <- mm[["inference"]][["conditional"]][["mu_alloc"]]
  valid <- vapply(posterior, function(sample) {
    .iwmde_posterior_ordinate_supports_bf(
      attr(sample, "posterior_ordinate", exact = TRUE)
    )
  }, logical(1))
  bf_posterior <- .marginal_means_bf_posterior(posterior[valid])
  class(bf_posterior) <- class(posterior)
  expected <- BayesTools::Savage_Dickey_BF(
    posterior            = bf_posterior,
    null_hypothesis      = mm[["null_hypothesis"]],
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )
  hypothesis_bf <- BayesTools::hypothesis_BF(
    posterior      = mm[["inference"]][["conditional"]][["mu_alloc"]][["alternate"]],
    hypothesis     = paste0("mu_alloc = ", mm[["null_hypothesis"]]),
    parameter      = "mu_alloc",
    columns        = "all",
    density_method = "precomputed"
  )

  expect_true(all(valid[c("alternate", "random")]))
  expect_false(valid[["systematic"]])
  expect_equal(
    as.numeric(refreshed[["inference"]][["mu_alloc"]][["alternate"]]),
    as.numeric(expected[["alternate"]])
  )
  expect_equal(
    as.numeric(refreshed[["inference"]][["mu_alloc"]][["random"]]),
    as.numeric(expected[["random"]])
  )
  expect_true(is.na(refreshed[["inference"]][["mu_alloc"]][["systematic"]]))
  expect_match(
    attr(refreshed[["inference"]][["mu_alloc"]][["systematic"]], "warnings"),
    "posterior ordinate was unavailable"
  )
  expect_true(is.finite(attr(expected[["alternate"]], "BF_error_percent")))
  expect_equal(attr(hypothesis_bf, "raw_BF"), as.numeric(expected[["alternate"]]),
               tolerance = 1e-12)
  expect_equal(as.numeric(hypothesis_bf[["BF_error"]]),
               attr(expected[["alternate"]], "BF_error_percent"),
               tolerance = 1e-12)
  expect_equal(hypothesis_bf[["method"]], "Savage-Dickey (precomputed)")
})


test_that("marginal_means does not mix qCMDE and KDE BF refresh within a parameter", {

  samples <- list(
    level_a = structure(
      rnorm(20),
      posterior_density = list(
        diagnostics = .marginal_means_iwmde_density_diagnostics("q_grid_cmde")
      )
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
  posterior_ordinate[["diagnostics"]][["max_weight_share"]] <- .2
  posterior_ordinate[["diagnostics"]][["active_mass"]] <- 1
  posterior_ordinate[["diagnostics"]][["normalization_integral"]] <- 1
  posterior_ordinate[["diagnostics"]][["normalization_mass_ratio"]] <- 1
  expect_true(.iwmde_posterior_ordinate_supports_bf(posterior_ordinate))
  posterior_ordinate[["diagnostics"]][["normalization_integral"]] <- .70
  posterior_ordinate[["diagnostics"]][["normalization_mass_ratio"]] <- 1 / .70
  expect_false(.iwmde_posterior_ordinate_supports_bf(posterior_ordinate))
  posterior_ordinate[["diagnostics"]][["normalization_integral"]] <- .92
  posterior_ordinate[["diagnostics"]][["normalization_mass_ratio"]] <- 1 / .92
  expect_true(.iwmde_posterior_ordinate_supports_bf(posterior_ordinate))
  posterior_ordinate[["diagnostics"]][["normalization_integral"]] <- NA_real_
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
      ordinate_relative_change = 0,
      max_normalizer_relative_change = 0,
      normalization_range    = c(-1, 1)
    )
  )
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))
  qcmde_ordinate[["diagnostics"]][["ordinate_relative_change"]] <- .30
  expect_false(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))
  qcmde_ordinate[["diagnostics"]][["ordinate_relative_change"]] <- NA_real_
  expect_false(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))
  qcmde_ordinate[["diagnostics"]][["ordinate_relative_change"]] <- 0
  qcmde_ordinate[["diagnostics"]][["active_mass"]] <- NA_real_
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))
  qcmde_ordinate[["diagnostics"]][["active_mass"]] <- 0
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))
  qcmde_ordinate[["diagnostics"]][["active_mass"]] <- .5
  qcmde_ordinate[["diagnostics"]][["normalization_range"]] <- c(.25, 1)
  expect_true(.iwmde_posterior_ordinate_supports_bf(qcmde_ordinate))

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
