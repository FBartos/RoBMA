test_that("bias component plots and switching use the categorical monitor", {
  dat <- data.frame(yi = c(.2, .4, .1, .3), sei = rep(.1, 4))
  setup <- review_compile_object(RoBMA(
    yi = yi, sei = sei, data = dat, measure = "SMD",
    prior_bias = list(prior_PET("normal", list(0, 1)), prior_PEESE("normal", list(0, 1))),
    only_priors = TRUE
  ))
  samples <- cbind(mu = .1, mu_indicator = 2, tau = .2, tau_indicator = 2,
                   PET = 0, PEESE = 0, bias_indicator = c(1, 2, 2, 3, 1, 1))
  fit <- review_attach_replay(setup, samples)
  frequency <- plot_diagnostic_trace(fit, parameter = "bias", statistic = "component",
                                     batch_size = 3, plot_type = "ggplot")
  expect_length(frequency, 3L)
  y <- function(plot) ggplot2::ggplot_build(plot)$data[[1L]]$y
  expect_equal(lapply(unname(frequency), y),
               list(c(1/3, 2/3), c(2/3, 0), c(0, 1/3)))
  inclusion <- plot_diagnostic(fit, parameter = "bias", type = "trace",
                               statistic = "inclusion", batch_size = 3, plot_type = "ggplot")
  switching <- plot_diagnostic_trace(fit, parameter = "bias", statistic = "switching",
                                     batch_size = 3, plot_type = "ggplot")
  expect_equal(y(inclusion), c(2/3, 1/3))
  expect_equal(y(switching), c(.5, .5))
})

test_that("random SD traces select the gate of the requested component", {
  fit <- review_bma_gate_replay()$object
  trace <- function(parameter) plot_diagnostic_trace(
    fit, parameter = parameter, component = "random", statistic = "switching",
    batch_size = 2, plot_type = "ggplot"
  )
  expect_equal(ggplot2::ggplot_build(trace("study: tau"))$data[[1L]]$y, c(1, 1))
  expect_equal(ggplot2::ggplot_build(trace("other: tau"))$data[[1L]]$y, c(0, 0))
  # The aggregate SD has a continuous prior and depends on both block gates.
  expect_error(trace("tau_total"), class = "RoBMA_indicator_unavailable")
})
