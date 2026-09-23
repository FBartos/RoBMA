context("Publication-bias mixture plots")

mock_bias_mixture_fit <- function() {

  bias <- BayesTools::prior_mixture(
    prior_list = list(
      BayesTools::prior_none(),
      BayesTools::prior_weightfunction(
        side    = "two-sided",
        steps   = c(0.05),
        weights = BayesTools::wf_cumulative(c(1, 1))
      ),
      BayesTools::prior_PET("normal", list(mean = 0, sd = 1)),
      BayesTools::prior_PEESE("normal", list(mean = 0, sd = 1))
    ),
    is_null = c(TRUE, FALSE, FALSE, FALSE)
  )

  data <- list(outcome = data.frame(
    yi  = c(0.1, 0.2),
    sei = c(0.1, 0.2)
  ))
  attr(data, "measure")          <- "SMD"
  attr(data, "outcome_type")     <- "norm"
  attr(data, "effect_direction") <- "positive"

  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "prior_list") <- list(
    mu   = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    bias = bias
  )

  structure(
    list(
      data   = data,
      priors = list(outcome = list(bias = bias)),
      fit    = fit
    ),
    class = c("RoBMA", "brma")
  )
}

test_that("publication-bias plots request mixed bias posterior", {

  calls <- list()

  testthat::local_mocked_bindings(
    as_mixed_posteriors = function(model, parameters, ...) {
      calls[["as_mixed_posteriors"]] <<- list(parameters = parameters)

      return(list(bias = structure(
        matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "omega[0,0.05]")),
        class = "mixed_posteriors"
      )))
    },
    plot_posterior = function(samples, parameter, ...) {
      calls[["plot_posterior"]] <<- list(parameter = parameter)

      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot_weightfunction(
    mock_bias_mixture_fit(),
    show_data = FALSE,
    plot_type = "ggplot"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(calls[["as_mixed_posteriors"]][["parameters"]], "bias")
  expect_equal(calls[["plot_posterior"]][["parameter"]], "weightfunction")
})

test_that("PET-PEESE plot requests location and mixed bias posterior", {

  calls <- list()

  testthat::local_mocked_bindings(
    .brma_parameter_select = function(object, parameter, ...) parameter,
    .package = "RoBMA"
  )

  testthat::local_mocked_bindings(
    as_mixed_posteriors = function(model, parameters, ...) {
      calls[["as_mixed_posteriors"]] <<- list(parameters = parameters)

      return(list(
        mu   = structure(c(0, 0), class = "mixed_posteriors"),
        bias = structure(
          matrix(0, nrow = 2, ncol = 2, dimnames = list(NULL, c("PET", "PEESE"))),
          class = "mixed_posteriors"
        )
      ))
    },
    plot_posterior = function(samples, parameter, ...) {
      calls[["plot_posterior"]] <<- list(parameter = parameter)

      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot_pet_peese(
    mock_bias_mixture_fit(),
    show_data = FALSE,
    plot_type = "ggplot"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(calls[["as_mixed_posteriors"]][["parameters"]], c("mu", "bias"))
  expect_equal(calls[["plot_posterior"]][["parameter"]], "PETPEESE")
})

test_that("PET-PEESE plot rejects moderator-dependent curves", {

  fit <- mock_bias_mixture_fit()
  attr(fit[["data"]], "mods") <- TRUE

  expect_error(
    plot_pet_peese(fit, show_data = FALSE),
    "not available for models with moderators",
    fixed = TRUE
  )
})

test_that("generic posterior plot requests mixed bias posterior", {

  calls <- list()

  testthat::local_mocked_bindings(
    .brma_parameter_select = function(object, parameter, ...) {
      if (identical(parameter, "weightfunction")) "omega" else parameter
    },
    .brma_parameter_select_entry = function(object, parameter, ...) {
      list(parameter = parameter, component = "bias")
    },
    .package = "RoBMA"
  )

  testthat::local_mocked_bindings(
    as_mixed_posteriors = function(model, parameters, conditional = NULL, ...) {
      calls[["as_mixed_posteriors"]] <<- list(
        parameters  = parameters,
        conditional = conditional
      )

      return(list(bias = structure(
        matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "omega[0,0.05]")),
        class = "mixed_posteriors"
      )))
    },
    plot_posterior = function(samples, parameter, ...) {
      calls[["plot_posterior"]] <<- list(parameter = parameter)

      return(structure(list(), class = "mock_plot"))
    },
    .package = "BayesTools"
  )

  out <- plot(
    mock_bias_mixture_fit(),
    parameter = "weightfunction",
    plot_type = "ggplot"
  )

  expect_s3_class(out, "mock_plot")
  expect_equal(calls[["as_mixed_posteriors"]][["parameters"]], "bias")
  expect_null(calls[["as_mixed_posteriors"]][["conditional"]])
  expect_equal(calls[["plot_posterior"]][["parameter"]], "omega")
})
