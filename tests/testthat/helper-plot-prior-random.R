# Helpers of the tests of plot_prior() and print_prior() for the random-effect
# selectors of fitted objects (test-00-plot-prior-random.R and
# test-02-random-prior-plots.R).

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
  curve[, c("x", "y")]
}

# The plotted curve equals the expected density at its plotted values, to a
# relative error at every point.
.expect_random_prior_curve <- function(curve, expected, tolerance, info = NULL,
                                       minimum = 100L, support = NULL) {

  expect_gt(nrow(curve), minimum)
  expect_gt(length(expected), 0L)
  expect_length(expected, nrow(curve))
  if (!nrow(curve) || length(expected) != nrow(curve)) {
    return(invisible(FALSE))
  }

  # A renderer closes a density at a finite support edge with a duplicate x
  # row. Its outside-facing row is zero; the other row retains the declared
  # boundary density. Identify these rows from support and x, never observed y.
  if (!is.null(support)) {
    expect_length(support, 2L)
    expect_true(all(!is.na(support)) && support[[1L]] < support[[2L]])
    for (side in seq_len(2L)) {
      boundary <- which(curve[["x"]] == support[[side]])
      expect_lte(length(boundary), 2L)
      if (length(boundary) == 2L) {
        closure <- if (side == 1L) boundary[[1L]] else boundary[[2L]]
        expected[[closure]] <- 0
      }
    }
  }
  singular <- is.infinite(expected) & expected > 0
  expect_true(all(is.finite(expected) | singular))
  if (any(singular)) {
    expect_true(!is.null(support) && all(curve[["x"]][singular] %in% support))
    expect_true(all(is.infinite(curve[["y"]][singular]) & curve[["y"]][singular] > 0))
  }
  regular <- !singular
  expect_true(any(regular))
  error <- abs(curve[["y"]][regular] - expected[regular]) / pmax(abs(expected[regular]), 1e-300)
  expect_true(
    max(error) < tolerance,
    info = paste0(info, " maximum relative error ", signif(max(error), 3L))
  )
}
