if (!exists("skip_if_no_vdiffr_snapshots", mode = "function")) {
  source(testthat::test_path("common-functions.R"))
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

.write_canonical_svg <- function(plot, file, title = "") {

  vdiffr::write_svg(plot, file, title)
  .canonicalize_dense_diamond_polygons(file)
}


# Metafor 5.0.0 increased each straight diamond edge from two endpoints to
# 10,000 collinear vertices. Collapse only that exact redundant encoding so
# snapshots remain sensitive to the corners and every other SVG element.
.canonicalize_dense_diamond_polygons <- function(file,
                                                  vertices_per_edge = 10000L) {

  svg      <- readLines(file, warn = FALSE)
  polygons <- grep("^<polygon points='", svg)

  for (i in polygons) {
    parts <- regmatches(
      svg[i],
      regexec("^<polygon points='([^']*)'(.*)$", svg[i])
    )[[1L]]
    if (length(parts) != 3L) {
      next
    }

    points <- strsplit(trimws(parts[2L]), " +")[[1L]]
    if (length(points) != 4L * vertices_per_edge) {
      next
    }

    point_parts <- strsplit(points, ",", fixed = TRUE)
    if (any(lengths(point_parts) != 2L)) {
      next
    }
    coordinates <- suppressWarnings(matrix(
      as.numeric(unlist(point_parts, use.names = FALSE)),
      ncol  = 2L,
      byrow = TRUE
    ))
    if (any(!is.finite(coordinates))) {
      next
    }

    corner_index <- seq.int(
      from       = 1L,
      by         = vertices_per_edge,
      length.out = 4L
    )
    corners      <- coordinates[corner_index, , drop = FALSE]
    fractions    <- seq(0, 1, length.out = vertices_per_edge)
    straight     <- vapply(seq_len(4L), function(edge) {

      edge_index <- seq.int(
        from       = corner_index[[edge]],
        length.out = vertices_per_edge
      )
      next_edge   <- if (edge == 4L) 1L else edge + 1L
      expected    <- outer(1 - fractions, corners[edge, ]) +
        outer(fractions, corners[next_edge, ])
      error_scale <- max(1, abs(expected), abs(coordinates[edge_index, ]))

      # svglite stores coordinates to 0.01 SVG unit. Permit only that
      # serialization loss plus floating-point error in this check.
      all(abs(coordinates[edge_index, ] - expected) <=
            .01 + 16 * .Machine$double.eps * error_scale)
    }, logical(1))
    if (!all(straight)) {
      next
    }

    corner_points <- points[corner_index]
    svg[i] <- paste0(
      "<polygon points='", paste(corner_points, collapse = " "), " '", parts[3L]
    )
  }

  writeLines(svg, file, useBytes = TRUE)
  return(invisible(file))
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
