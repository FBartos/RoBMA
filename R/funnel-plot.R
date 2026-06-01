# ============================================================================ #
# funnel-plot.R
# ============================================================================ #

# .funnel_plot_base
# ---------------------------------------------------------------------------- #
#
# Create funnel plot using base R graphics.
#
# @param data list of funnel plot data from .funnel_data_* functions
# @param dots list of graphical parameters
#
# @return NULL invisibly
#
# ---------------------------------------------------------------------------- #
.funnel_plot_base <- function(data, dots) {
  # extract graphical parameters
  pch   <- dots[["pch"]]
  col   <- dots[["col"]]
  bg    <- dots[["bg"]]
  cex   <- dots[["cex"]]
  back  <- dots[["back"]]
  shade <- dots[["shade"]]
  lty   <- dots[["lty"]]
  col_line <- dots[["col.line"]]
  col_refline <- dots[["col.refline"]]
  main <- dots[["main"]]
  las  <- dots[["las"]]

  # extract data components
  df_points <- data$points
  df_funnel <- data$funnel
  df_funnel_edge1 <- data$funnel_edge1
  df_funnel_edge2 <- data$funnel_edge2
  df_background <- data$background
  x_range  <- data$x_range
  se_range <- data$y_range
  xlab <- data$xlab
  ylab <- data$ylab
  # extract data components
  df_points <- data$points
  df_funnel <- data$funnel
  df_funnel_edge1 <- data$funnel_edge1
  df_funnel_edge2 <- data$funnel_edge2
  df_background <- data$background

  # use limits from data object or dots for PLOTTING LIMITS
  if (!is.null(dots[["xlim"]])) xlim_plot <- dots[["xlim"]] else xlim_plot <- data$xlim
  if (!is.null(dots[["ylim"]])) ylim_plot <- dots[["ylim"]] else ylim_plot <- data$ylim

  # Generate pretty ticks based on the limits
  x_ticks <- pretty(xlim_plot)
  y_ticks <- pretty(ylim_plot)

  xlab <- data$xlab
  ylab <- data$ylab
  df_refline <- data$refline

  # set up the plot area with reversed y-axis
  graphics::plot(
    NA, NA,
    xlim = range(xlim_plot),
    ylim = rev(range(ylim_plot)),
    xlab = xlab,
    ylab = ylab,
    main = main,
    type = "n",
    axes = FALSE,
    las  = las
  )
  graphics::axis(1, at = x_ticks, las = las)
  graphics::axis(2, at = rev(y_ticks), las = las)

  # draw background polygon (if not suppressed)
  if (!is.na(back)) {
    graphics::polygon(df_background$x, df_background$y, col = back, border = NA)
  }

  # draw funnel polygon (if not suppressed)
  if (!is.na(shade)) {
    graphics::polygon(df_funnel$x, df_funnel$y, col = shade, border = NA)
  }

  # vertical reference line (now potentially curved)
  graphics::lines(df_refline$x, df_refline$y, lty = lty, col = col_refline)

  # funnel edge lines
  graphics::lines(df_funnel_edge1$x, df_funnel_edge1$y, lty = lty, col = col_line)
  graphics::lines(df_funnel_edge2$x, df_funnel_edge2$y, lty = lty, col = col_line)

  # plot points
  graphics::points(df_points$x, df_points$y, pch = pch, col = col, bg = bg, cex = cex)

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .funnel_plot_ggplot
# ---------------------------------------------------------------------------- #
#
# Create funnel plot using ggplot2.
#
# @param data list of funnel plot data from .funnel_data_* functions
# @param dots list of graphical parameters
#
# @return ggplot object
#
# ---------------------------------------------------------------------------- #
.funnel_plot_ggplot <- function(data, dots) {
  # extract graphical parameters
  pch   <- dots[["pch"]]
  col   <- dots[["col"]]
  bg    <- dots[["bg"]]
  size  <- dots[["size"]]
  back  <- dots[["back"]]
  shade <- dots[["shade"]]
  lty   <- dots[["lty"]]
  col_line <- dots[["col.line"]]
  col_refline <- dots[["col.refline"]]
  main <- dots[["main"]]

  # extract data components
  df_points <- data$points
  df_funnel <- data$funnel
  df_funnel_edge1 <- data$funnel_edge1
  df_funnel_edge2 <- data$funnel_edge2
  df_background <- data$background
  x_range   <- data$x_range
  se_range  <- data$y_range
  xlab      <- data$xlab
  ylab      <- data$ylab
  # extract data components
  df_points <- data$points
  df_funnel <- data$funnel
  df_funnel_edge1 <- data$funnel_edge1
  df_funnel_edge2 <- data$funnel_edge2
  df_background   <- data$background

  # use limits from data object or dots
  if (!is.null(dots[["xlim"]])) xlim_plot <- dots[["xlim"]] else xlim_plot <- data$xlim
  if (!is.null(dots[["ylim"]])) ylim_plot <- dots[["ylim"]] else ylim_plot <- data$ylim

  # Generate pretty ticks
  x_ticks <- pretty(xlim_plot)
  y_ticks <- pretty(ylim_plot)

  xlab <- data$xlab
  ylab <- data$ylab
  df_refline <- data$refline

  out <- ggplot2::ggplot()

  # add background polygon (if not suppressed)
  if (!is.na(back)) {
    out <- out +
      ggplot2::geom_polygon(
        mapping = ggplot2::aes(
          x = df_background$x,
          y = df_background$y
        ),
        fill = back
      )
  }

  # add funnel polygon (if not suppressed)
  if (!is.na(shade)) {
    out <- out +
      ggplot2::geom_polygon(
        mapping = ggplot2::aes(
          x = df_funnel$x,
          y = df_funnel$y
        ),
        fill = shade
      )
  }

  out <- .add_funnel_path_or_line(
    out    = out,
    data   = df_refline,
    lty    = lty,
    colour = col_refline
  )

  # add funnel edge lines
  out <- .add_funnel_path_or_line(
    out    = out,
    data   = df_funnel_edge1,
    lty    = lty,
    colour = col_line
  )
  out <- .add_funnel_path_or_line(
    out    = out,
    data   = df_funnel_edge2,
    lty    = lty,
    colour = col_line
  )

  # add points
  out <- out +
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = df_points$x,
        y = df_points$y
      ),
      colour = col,
      fill = bg,
      shape = pch,
      size = size
    )

  # set axis scales and labels
  out <- out +
    ggplot2::scale_x_continuous(
      breaks = x_ticks,
      limits = range(xlim_plot),
      name   = xlab
    ) +
    ggplot2::scale_y_reverse(
      breaks = rev(y_ticks),
      limits = rev(range(ylim_plot)),
      name   = ylab
    )

  # add title if specified
  if (!is.null(main)) {
    out <- out + ggplot2::ggtitle(main)
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .add_funnel_path_or_line
# ---------------------------------------------------------------------------- #
#
# Add a contour line to a ggplot object.
#
# PET/PEESE contours can be non-monotone in x because they are parameterized by
# standard error. ggplot2::geom_line() sorts by x, which can scramble such
# contours, so those cases must use geom_path().
#
# ---------------------------------------------------------------------------- #
.add_funnel_path_or_line <- function(out, data, lty, colour) {

  x_finite <- data[["x"]][is.finite(data[["x"]])]
  if (length(x_finite) < 3L) {
    is_monotone <- TRUE
  } else {
    x_diff      <- diff(x_finite)
    x_diff      <- x_diff[x_diff != 0]
    is_monotone <- length(x_diff) < 2L || all(x_diff > 0) || all(x_diff < 0)
  }

  if (is_monotone) {
    return(out +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = data[["x"]],
          y = data[["y"]]
        ),
        linetype = lty,
        colour   = colour
      )
    )
  }

  return(out +
    ggplot2::geom_path(
      data    = data,
      mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      linetype = lty,
      colour   = colour
    )
  )
}
