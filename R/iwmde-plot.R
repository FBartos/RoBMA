# ============================================================================ #
# IWMDE Plotting and Plot Data Helpers
# ============================================================================ #

.iwmde_kde <- function(values, xlim, n_points, mass = 1) {

  values <- values[is.finite(values)]
  if (length(unique(values)) < 2L) {
    return(NULL)
  }

  density <- stats::density(values, n = n_points, from = xlim[1], to = xlim[2])
  return(list(x = density[["x"]], y = density[["y"]] * mass))
}


.iwmde_histogram <- function(values, xlim, mass = 1) {

  values <- values[is.finite(values)]
  values <- values[values >= xlim[1] & values <= xlim[2]]
  if (length(values) == 0L) {
    values <- mean(xlim)
  }
  breaks <- seq(xlim[1], xlim[2], length.out = 18)

  hist_data <- graphics::hist(values, breaks = breaks, plot = FALSE)
  return(list(
    mids    = hist_data[["mids"]],
    density = hist_data[["density"]] * mass,
    breaks  = hist_data[["breaks"]]
  ))
}


.iwmde_trapz <- function(x, y) {

  if (length(x) < 2L || length(y) < 2L) {
    return(NA_real_)
  }

  return(sum(.iwmde_trapz_weights(x) * y))
}


.iwmde_trapz_columns <- function(x, y) {

  y <- as.matrix(y)
  if (nrow(y) != length(x) || length(x) < 2L) {
    return(rep(NA_real_, ncol(y)))
  }

  return(as.numeric(crossprod(.iwmde_trapz_weights(x), y)))
}


.iwmde_trapz_weights <- function(x) {

  n <- length(x)
  if (n < 2L) {
    return(rep(0, n))
  }

  dx      <- diff(x)
  weights <- numeric(n)
  weights[-n] <- weights[-n] + dx / 2
  weights[-1L] <- weights[-1L] + dx / 2

  return(weights)
}


.plot_iwmde_diagnostics_base <- function(x, ...) {

  dots <- list(...)
  n    <- length(x)
  side <- ceiling(sqrt(n))

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(side, side),
    mar   = c(2.5, 2.5, 2.2, .6),
    oma   = c(0, 0, 0, 0)
  )

  for (i in seq_along(x)) {
    .plot_iwmde_parameter_base(x[[i]], show_legend = i == 1L, dots = dots)
  }
  if (side * side > n) {
    for (i in seq_len(side * side - n)) {
      graphics::plot.new()
    }
  }

  return(invisible(NULL))
}


.plot_iwmde_parameter_base <- function(x, show_legend, dots) {

  parameter <- x[["parameter"]]
  if (!identical(x[["status"]], "ok")) {
    graphics::plot.new()
    graphics::title(main = parameter, cex.main = .85)
    graphics::text(.5, .55, x[["status"]], cex = .9)
    graphics::text(.5, .42, x[["reason"]], cex = .65)
    return(invisible(NULL))
  }

  ymax <- max(
    x[["histogram"]][["density"]],
    if (!is.null(x[["kde"]])) x[["kde"]][["y"]] else 0,
    x[["iwmde"]][["y"]],
    na.rm = TRUE
  )
  if (!is.finite(ymax) || ymax <= 0) {
    ymax <- 1
  }

  graphics::plot(
    NA_real_,
    NA_real_,
    type     = "n",
    xlim     = x[["xlim"]],
    ylim     = c(0, ymax * 1.08),
    main     = parameter,
    xlab     = "",
    ylab     = "",
    cex.main = .85
  )
  graphics::rect(
    xleft   = utils::head(x[["histogram"]][["breaks"]], -1L),
    ybottom = 0,
    xright  = utils::tail(x[["histogram"]][["breaks"]], -1L),
    ytop    = x[["histogram"]][["density"]],
    col     = if (is.null(dots[["hist_col"]])) "#D9D9D9" else dots[["hist_col"]],
    border  = if (is.null(dots[["hist_border"]])) "white" else dots[["hist_border"]]
  )

  if (!is.null(x[["kde"]])) {
    graphics::lines(
      x[["kde"]][["x"]],
      x[["kde"]][["y"]],
      col = if (is.null(dots[["kde_col"]])) "#2C7FB8" else dots[["kde_col"]],
      lwd = if (is.null(dots[["kde_lwd"]])) 1.5 else dots[["kde_lwd"]]
    )
  }
  graphics::lines(
    x[["iwmde"]][["x"]],
    x[["iwmde"]][["y"]],
    col = if (is.null(dots[["iwmde_col"]])) "#D7191C" else dots[["iwmde_col"]],
    lwd = if (is.null(dots[["iwmde_lwd"]])) 1.5 else dots[["iwmde_lwd"]]
  )

  points <- x[["point_masses"]]
  if (!is.null(points) && nrow(points) > 0L) {
    graphics::segments(
      x0  = points[["x"]],
      y0  = 0,
      x1  = points[["x"]],
      y1  = ymax * pmin(1, points[["mass"]] / max(points[["mass"]])),
      col = "#D7191C",
      lwd = 1.3,
      lty = 3
    )
  }

  if (show_legend) {
    estimator_label <- if (identical(
      x[["diagnostics"]][["estimator"]],
      "q_grid_cmde"
    )) {
      "q-grid CMDE"
    } else {
      "IWMDE"
    }
    graphics::legend(
      "topright",
      legend = c("Histogram", "KDE", estimator_label),
      col    = c("#D9D9D9", "#2C7FB8", "#D7191C"),
      lty    = c(NA, 1, 1),
      lwd    = c(NA, 1.5, 1.5),
      pch    = c(15, NA, NA),
      pt.cex = 1.2,
      bty    = "n",
      cex    = .65
    )
  }

  return(invisible(NULL))
}
