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
