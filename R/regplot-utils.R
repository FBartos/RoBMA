# ---------------------------------------------------------------------------- #
# .set_dots_regplot
# ---------------------------------------------------------------------------- #
#
# Process dots arguments for regression plot with sensible defaults.
#
# @param ... additional graphical arguments
#
# @return list of graphical parameters with defaults applied
#
# ---------------------------------------------------------------------------- #
.set_dots_regplot <- function(...) {

  dots <- list(...)

  # point styling
  dots <- .plot_point_style_defaults(dots)

  # line styling
  if (is.null(dots[["lcol"]]))     dots[["lcol"]]     <- "black"
  if (is.null(dots[["lwd"]]))      dots[["lwd"]]      <- 2
  .check_plot_positive_scalar(dots[["lwd"]], "lwd")

  # band styling
  if (is.null(dots[["shade"]]))    dots[["shade"]]    <- TRUE
  if (is.null(dots[["col.ci"]]))   dots[["col.ci"]]   <- "gray70"
  if (is.null(dots[["col.pi"]]))   dots[["col.pi"]]   <- "gray85"
  if (is.null(dots[["col.si"]]))   dots[["col.si"]]   <- "gray92"
  if (is.null(dots[["alpha.ci"]])) dots[["alpha.ci"]] <- 0.4
  if (is.null(dots[["alpha.pi"]])) dots[["alpha.pi"]] <- 0.2
  if (is.null(dots[["alpha.si"]])) dots[["alpha.si"]] <- 0.15
  BayesTools::check_bool(dots[["shade"]], "shade")
  .check_plot_numeric(dots[["alpha.ci"]], "alpha.ci", check_length = 1, lower = 0, upper = 1, allow_null = FALSE)
  .check_plot_numeric(dots[["alpha.pi"]], "alpha.pi", check_length = 1, lower = 0, upper = 1, allow_null = FALSE)
  .check_plot_numeric(dots[["alpha.si"]], "alpha.si", check_length = 1, lower = 0, upper = 1, allow_null = FALSE)

  # categorical moderator jitter
  if (is.null(dots[["jitter"]]))      dots[["jitter"]]      <- 0.2
  if (is.null(dots[["box.width"]]))   dots[["box.width"]]   <- 0.5
  if (is.null(dots[["dodge.width"]])) dots[["dodge.width"]] <- 0.75
  .check_plot_numeric(dots[["jitter"]], "jitter", check_length = 1, lower = 0, allow_null = FALSE)
  .check_plot_positive_scalar(dots[["box.width"]], "box.width")
  .check_plot_numeric(dots[["dodge.width"]], "dodge.width", check_length = 1, lower = 0, allow_null = FALSE)

  # title (NULL = no title by default)
  if (is.null(dots[["main"]]))     dots[["main"]]     <- NULL
  if (is.null(dots[["las"]]))      dots[["las"]]      <- 1
  .check_plot_label(dots[["main"]], "main")
  BayesTools::check_int(dots[["las"]], "las", lower = 0)
  if (!dots[["las"]] %in% 0:3) {
    stop("'las' must be one of 0, 1, 2, or 3.", call. = FALSE)
  }
  .check_plot_positive_scalar(dots[["cex"]], "cex")
  .check_plot_positive_scalar(dots[["size"]], "size")

  return(dots)
}
