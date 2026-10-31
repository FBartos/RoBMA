source(testthat::test_path("common-functions.R"))

test_that("public correlation posterior plots retain wider ranges and bounded layers", {

  skip_if_not_installed("ggplot2")
  skip_if_missing_fits("dat.lehmann2018_RoBMA")
  fit <- load_fit("dat.lehmann2018_RoBMA")

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  graphics::par(xaxs = "i")
  expect_silent(plot(
    fit, parameter = "mu", output_measure = "COR", prior = TRUE,
    xlim = c(-1.5, 1.5), density_method = "KDE"
  ))
  expect_equal(graphics::par("usr")[1:2], c(-1.5, 1.5))

  plot <- plot(fit, parameter = "mu", output_measure = "COR", prior = TRUE,
               xlim = c(-1.5, 1.5), density_method = "KDE", plot_type = "ggplot")
  expect_equal(plot$scales$get_scales("x")$limits, c(-1.5, 1.5))
  built <- ggplot2::ggplot_build(plot)
  coordinates <- unlist(lapply(built$data, function(layer) layer$x))
  expect_true(length(coordinates) > 0L)
  expect_true(all(is.finite(coordinates)))
  expect_true(all(coordinates >= -1 & coordinates <= 1))
})
