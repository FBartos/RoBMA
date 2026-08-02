if (!exists("skip_if_no_vdiffr_snapshots", mode = "function")) {
  source(testthat::test_path("common-functions.R"))
}
if (!exists(".write_canonical_svg", mode = "function")) {
  source(testthat::test_path("helper-visual-writer.R"))
}

expect_vdiffr_snapshot <- function(title, fig, ...) {

  skip_if_no_vdiffr_snapshots()
  testthat::skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title,
    fig,
    ...,
    writer = .write_canonical_svg
  )
}

expect_brma_plot_snapshot <- function(name, plot) {

  expect_vdiffr_snapshot(name, plot)
}

.with_temp_plot_device <- function(expr) {

  file <- tempfile(fileext = ".png")
  grDevices::png(filename = file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)

  force(expr)
}

.is_ggplot <- function(x) {

  inherits(x, "ggplot")
}

ggplot_geom_layer_data <- function(plot, geom_class, index = 1L) {

  layer_index <- which(vapply(
    plot[["layers"]],
    function(layer) inherits(layer[["geom"]], geom_class),
    logical(1)
  ))

  if (length(layer_index) < index) {
    stop("Missing ggplot layer for geom class '", geom_class, "'.", call. = FALSE)
  }

  return(plot[["layers"]][[layer_index[[index]]]][["data"]])
}

expect_side_by_side_plot <- function(name, metafor, brma) {

  expect_vdiffr_snapshot(name, function() {
    old_par <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(old_par), add = TRUE)
    metafor()
    brma()
  })
}
