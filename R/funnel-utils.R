# ============================================================================ #
# funnel-utils.R
# ============================================================================ #

# .clip_line_x
# ---------------------------------------------------------------------------- #
#
# Clip a line (x, y) to x-limits, interpolating exact intersections.
# Assumes y is monotonic (funnel lines usually are monotonic with SE).
#
# @param x numeric vector of x coordinates
# @param y numeric vector of y coordinates
# @param xlim numeric vector of length 2 (min, max)
#
# @return data.frame with xe and ye
#
# ---------------------------------------------------------------------------- #
.clip_line_x <- function(x, y, xlim) {

  if (length(x) < 2) return(data.frame(x=x, y=y))

  xmin <- min(xlim)
  xmax <- max(xlim)

  x_out <- numeric(0)
  y_out <- numeric(0)

  for (i in 1:(length(x) - 1)) {
    x1 <- x[i]
    y1 <- y[i]
    x2 <- x[i+1]
    y2 <- y[i+1]

    # strictly inside
    if (x1 >= xmin && x1 <= xmax) {
      x_out <- c(x_out, x1)
      y_out <- c(y_out, y1)
    }

    if (x1 != x2) {
      t_cross <- numeric(0)
      x_cross <- numeric(0)
      y_cross <- numeric(0)
      for (boundary in c(xmin, xmax)) {
        t <- (boundary - x1) / (x2 - x1)
        if (t > 0 && t < 1) {
          t_cross <- c(t_cross, t)
          x_cross <- c(x_cross, boundary)
          y_cross <- c(y_cross, y1 + (y2 - y1) * t)
        }
      }
      if (length(t_cross) > 0L) {
        cross_order <- order(t_cross)
        x_out <- c(x_out, x_cross[cross_order])
        y_out <- c(y_out, y_cross[cross_order])
      }
    }
  }

  # handle last point
  xn <- x[length(x)]
  yn <- y[length(y)]
  if (xn >= xmin && xn <= xmax) {
    x_out <- c(x_out, xn)
    y_out <- c(y_out, yn)
  }

  return(data.frame(x = x_out, y = y_out))
}

# .normalize_funnel_max_samples
# ---------------------------------------------------------------------------- #
#
# Validate max posterior samples used by model-averaged funnel contours.
#
# ---------------------------------------------------------------------------- #
.normalize_funnel_max_samples <- function(max_samples) {

  return(.normalize_max_samples(max_samples, "max_samples"))
}


# ---------------------------------------------------------------------------- #
# .set_dots_funnel
# ---------------------------------------------------------------------------- #
#
# Process dots arguments for funnel plot with sensible defaults.
#
# Sets defaults for all graphical parameters that can be customized in the
# funnel plot. Users can override any of these by passing them in ...
#
# @param ... additional graphical arguments
#
# @return list of graphical parameters with defaults applied
#
# ---------------------------------------------------------------------------- #
.set_dots_funnel <- function(dots) {

  # point styling (use metafor-style argument names)
  dots <- .plot_point_style_defaults(dots)

  # funnel region styling
  if (is.null(dots[["back"]])) dots[["back"]] <- "grey" # background fill

  if (is.null(dots[["shade"]])) dots[["shade"]] <- "white" # funnel fill

  # line styling
  if (is.null(dots[["lty"]])) dots[["lty"]] <- "dotted" # line type
  if (is.null(dots[["col.line"]])) dots[["col.line"]] <- "black" # funnel edge color

  # reference line (position set in data functions based on mode)
  # if user provides refline, it will be passed through to override
  if (is.null(dots[["col.refline"]])) dots[["col.refline"]] <- "black" # reference line color

  # axis labels (now set separately in data functions)
  # main title (NULL = no title by default)
  if (is.null(dots[["main"]])) dots[["main"]] <- NULL
  if (is.null(dots[["las"]]))  dots[["las"]]  <- 1

  # axis limits (NULL = auto)
  # xlim and ylim are left as NULL if not provided (auto-computed)

  # data return flag
  if (is.null(dots[["as_data"]])) dots[["as_data"]] <- FALSE

  return(dots)
}
