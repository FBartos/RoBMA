context("Conditional posterior plots")

mock_plot_fit <- function(class = c("BMA.norm", "RoBMA", "brma")) {

  data <- structure(
    list(),
    measure      = "SMD",
    outcome_type = "norm"
  )

  fit <- structure(list(), class = "BayesTools_fit")

  structure(
    list(
      data   = data,
      priors = list(outcome = list()),
      fit    = fit
    ),
    class = class
  )
}

test_that("plot.brma forwards conditional parameter for RoBMA objects", {

  calls <- list()

  testthat::local_mocked_bindings(
    .brma_parameter_select = function(object, parameter, ...) parameter,
    .brma_parameter_select_entry = function(object, parameter, ...) {
      list(parameter = parameter, component = "mods")
    },
    .package = "RoBMA"
  )

  testthat::local_mocked_bindings(
    as_mixed_posteriors = function(model, parameters, conditional = NULL,
                                   transform_scaled = FALSE, ...) {
      calls[["as_mixed_posteriors"]] <<- list(
        parameters       = parameters,
        conditional      = conditional,
        transform_scaled = transform_scaled
      )

      samples <- list()
      samples[[parameters]] <- structure(c(0, 1), class = "mixed_posteriors")

      return(samples)
    },
    plot_posterior = function(samples, parameter, ...) {
      calls[["plot_posterior"]] <<- list(
        samples   = samples,
        parameter = parameter
      )

      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    mock_plot_fit(),
    parameter   = "mu",
    conditional = TRUE,
    plot_type   = "ggplot"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(calls[["as_mixed_posteriors"]][["parameters"]], "mu")
  expect_equal(calls[["as_mixed_posteriors"]][["conditional"]], "mu")
  expect_true(calls[["as_mixed_posteriors"]][["transform_scaled"]])
  expect_equal(calls[["plot_posterior"]][["parameter"]], "mu")
})

test_that("plot.brma rejects conditional parameter for non-RoBMA objects", {

  expect_error(
    plot(
      mock_plot_fit(class = c("brma.norm", "brma")),
      parameter   = "mu",
      conditional = TRUE
    ),
    "RoBMA objects"
  )
})
