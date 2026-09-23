test_that("correlation plot transformations declare support without weakening input validation", {

  transformation <- .effect_plot_transformation(.effect_output_setup_measure(
    input_measure = "ZCOR", output_measure = "COR"
  ))
  expect_identical(transformation$output_support, c(-1, 1))
  expect_identical(transformation$inv(c(-1, 0, 1)), c(-Inf, 0, Inf))
  expect_error(transformation$inv(1 + .Machine$double.eps), "within", fixed = TRUE)
  expect_error(.cor_to_z(NA_real_), "finite", fixed = TRUE)

  plotted <- stats::density(
    BayesTools::prior("normal", list(0, .4)),
    x_seq = c(-1.5, -1, -.5, 0, .5, 1, 1.5),
    transformation = transformation, transformation_settings = TRUE
  )
  expect_equal(plotted$x, c(-.5, 0, .5))
  expect_equal(plotted$y, stats::dnorm(atanh(plotted$x), sd = .4) / (1 - plotted$x^2))
  expect_identical(attr(plotted, "x_range"), c(-1.5, 1.5))

  log_plot <- list(
    active = TRUE, transform = "LOG", output_measure = "COR",
    transformation = .log_plot_transformation()
  )
  expect_null(.effect_plot_transformation(log_plot)$output_support)
})
