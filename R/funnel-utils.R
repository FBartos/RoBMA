# ============================================================================ #
# funnel-utils.R
# ============================================================================ #

# .clip_line_x
# ---------------------------------------------------------------------------- #
#
# Clip a line (x, y) to x-limits, interpolating exact intersections.
# Disconnected visible runs are separated by NA coordinates.
#
# @param x numeric vector of x coordinates
# @param y numeric vector of y coordinates
# @param xlim numeric vector of length 2 (min, max)
#
# @return data.frame with x and y, including NA breaks between visible runs
#
# ---------------------------------------------------------------------------- #
.clip_line_x <- function(x, y, xlim) {

  xmin <- min(xlim)
  xmax <- max(xlim)

  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Funnel contour coordinates must be finite.", call. = FALSE)
  }
  if (length(x) < 2L) {
    keep <- x >= xmin & x <= xmax
    return(data.frame(x = x[keep], y = y[keep]))
  }

  x_out <- y_out <- numeric()
  # Each segment has one visible parameter interval. Only join adjacent
  # intervals when the common original vertex remains visible.
  previous_visible <- FALSE
  for (i in seq_len(length(x) - 1L)) {
    dx <- x[i + 1L] - x[i]
    bounds <- if (dx == 0) {
      if (x[i] >= xmin && x[i] <= xmax) c(0, 1) else c(1, 0)
    } else {
      crossing <- sort((c(xmin, xmax) - x[i]) / dx)
      c(max(0, crossing[1L]), min(1, crossing[2L]))
    }
    if (bounds[1L] > bounds[2L]) {
      previous_visible <- FALSE
      next
    }
    current_x <- x[i] + bounds * dx
    current_y <- y[i] + bounds * (y[i + 1L] - y[i])
    connected <- previous_visible && bounds[1L] == 0
    if (length(x_out) && !connected) {
      x_out <- c(x_out, NA_real_)
      y_out <- c(y_out, NA_real_)
    }
    keep <- if (connected) 2L else 1:2
    x_out <- c(x_out, current_x[keep])
    y_out <- c(y_out, current_y[keep])
    previous_visible <- bounds[2L] == 1
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

  if ("as_data" %in% names(dots)) {
    BayesTools::check_bool(dots[["as_data"]], "as_data")
  }
  .check_plot_numeric(dots[["xlim"]], "xlim", check_length = 2)
  .check_plot_numeric(dots[["ylim"]], "ylim", check_length = 2, lower = 0)
  .check_plot_numeric(dots[["refline"]], "refline", check_length = 1)
  if (!is.null(dots[["ylim"]]) && dots[["ylim"]][1L] == dots[["ylim"]][2L]) {
    stop("'ylim' must contain two distinct standard-error limits.", call. = FALSE)
  }

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
