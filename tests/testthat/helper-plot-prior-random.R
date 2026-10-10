# Helpers of the tests of plot_prior() and print_prior() for the random-effect
# selectors of fitted objects (test-00-plot-prior-random.R and
# test-02-plot-prior-random.R).

# The layers of a built ggplot by geom: the curve and the point-mass arrows.
.random_prior_layers <- function(plot) {

  built <- ggplot2::ggplot_build(plot)[["data"]]
  geoms <- vapply(plot[["layers"]], function(layer) {
    class(layer[["geom"]])[[1L]]
  }, character(1))

  list(
    curve  = built[[which(geoms == "GeomLine")[[1L]]]],
    arrows = built[geoms == "GeomSegment"]
  )
}

.random_prior_curve <- function(plot) {

  curve <- .random_prior_layers(plot)[["curve"]]
  curve[is.finite(curve[["y"]]) & curve[["y"]] > 0, c("x", "y")]
}

# The plotted curve equals the expected density at its plotted values, to a
# relative error at every point.
.expect_random_prior_curve <- function(curve, expected, tolerance, info = NULL,
                                       minimum = 100L) {

  expect_gt(nrow(curve), minimum)
  error <- abs(curve[["y"]] - expected) / pmax(abs(expected), 1e-300)
  expect_true(
    max(error) < tolerance,
    info = paste0(info, " maximum relative error ", signif(max(error), 3L))
  )
}

.random_prior_device <- function(expr) {

  file <- tempfile(fileext = ".png")
  grDevices::png(filename = file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)

  force(expr)
}
