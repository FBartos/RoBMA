# ============================================================================ #
# brma.funnel.R
# ============================================================================ #
#
# This file contains the funnel plot functions for brma objects.
#
# Two modes of dispatch:
# 1) No mods AND no scale: Shows outcome residuals against sampling distribution
# 2) All other cases: Shows standardized residuals against standard normal
#
# ============================================================================ #


### funnel plot functions ----
#' @export
funnel <- function(x, ...) UseMethod("funnel")


#' @title Funnel Plot for brma Object
#'
#' @description \code{funnel.brma} creates a funnel plot for a fitted brma object.
#' For intercept-only models without scale regression, this displays observed
#' effect sizes against the expected sampling distribution. For models with
#' moderators or scale regression, this displays standardized residuals against
#' the standard normal distribution.
#'
#' @param x a fitted brma object
#' @param residual whether to use standardized residuals mode. Defaults to
#' not specified, which means the function automatically determines the mode:
#' \itemize{
#'   \item \code{FALSE} (auto): For intercept-only models without scale
#'     regression, displays observed effect sizes minus the pooled estimate
#'     against the sampling distribution funnel.
#'   \item \code{TRUE} (auto or explicit): For models with moderators or scale
#'     regression, or when explicitly set to \code{TRUE}, displays standardized
#'     residuals against the standard normal funnel.
#' }
#' @param type the type of residuals to use when in residual mode.
#' Options are:
#' \itemize{
#'   \item \code{"LOO-PIT"} (alias: \code{"rstudent"}; default): Leave-one-out probability integral transform
#'     residuals (most robust, works for all model types).
#'   \item \code{"rstandard"}: Internally standardized residuals using
#'     the hat matrix. Only available for normal outcome models.
#'   \item \code{"outcome"}: Raw outcome residuals (observed - fitted).
#' }
#' Only used when funnel is in residual mode.
#' @param sampling_heterogeneity whether heterogeneity should be incorporated
#' into the sampling distribution funnel. Defaults to \code{TRUE}. Only used
#' when \code{residual = FALSE} (outcome mode).
#' @param sampling_bias whether publication bias should be incorporated into the
#' sampling distribution funnel. Defaults to \code{TRUE}. Only used when
#' \code{residual = FALSE} (outcome mode). When \code{TRUE} and the model
#' includes selection models (weightfunction), uses weighted normal quantiles.
#' When \code{TRUE} and the model includes PET/PEESE, incorporates the expected
#' skew from these regression adjustments.
#' @param plot_type whether to use a base plot \code{"base"} or ggplot2
#' \code{"ggplot"} for plotting. Defaults to \code{"base"}.
#' @param ... additional graphical arguments to customize the plot appearance:
#' \describe{
#'   \item{xlim, ylim}{numeric vectors of length 2 specifying axis limits}
#'   \item{xlab, ylab}{character strings for axis labels}
#'   \item{main}{character string for plot title (default: no title)}
#'   \item{pch}{point symbol (default: 21, filled circle). Use standard R pch values.}
#'   \item{col}{point border color (default: "black")}
#'   \item{bg}{point fill/background color (default: "black")}
#'   \item{cex}{point size multiplier for base graphics (default: 1)}
#'   \item{size}{point size for ggplot2 (default: 2)}
#'   \item{back}{background region color (default: "grey"). Set to \code{NA} to suppress.}
#'   \item{shade}{funnel region fill color (default: "white"). Set to \code{NA} to suppress.}
#'   \item{lty}{line type for funnel edges and center line (default: "dotted")}
#'   \item{col.line}{color for funnel edge lines (default: "black")}
#'   \item{refline}{position of vertical reference line (default: 0)}
#'   \item{col.refline}{color of vertical reference line (default: "black")}
#'   \item{as_data}{if \code{TRUE}, returns plot data instead of creating plot}
#' }
#'
#' @details
#' The funnel plot has two modes depending on the model complexity:
#'
#' \strong{Outcome mode} (intercept-only models without scale regression):
#' Displays residuals (observed - pooled estimate) on the x-axis and standard
#' errors on the y-axis. The funnel region represents the 95% prediction region
#' based on the sampling distribution, optionally incorporating heterogeneity
#' and publication bias.
#'
#'
#' \strong{Residual mode} (models with moderators or scale regression):
#' Displays residuals on the x-axis and standard errors on the
#' y-axis. The funnel region represents the 95% prediction region
#' based on the standard error. Under a correctly specified model,
#' residuals should fall within this region approximately 95% of the time.
#'
#' The \code{type} argument controls how residuals are computed in residual
#' mode. See \code{\link{residuals.brma}} for details on each type.
#'
#' For GLMM models, observed effect sizes are computed from the raw frequency
#' data using formulas equivalent to \code{metafor::escalc}.
#'
#' @return \code{funnel.brma} returns \code{NULL} invisibly if \code{plot_type = "base"}
#' or a ggplot object if \code{plot_type = "ggplot"}.
#'
#' @examples \dontrun{
#' # Simple meta-analysis: outcome mode (yi - mu vs sampling distribution)
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#' funnel(fit)
#'
#' # Force residual mode even for simple model
#' funnel(fit, residual = TRUE)
#'
#' # Meta-regression: automatically uses residual mode
#' fit_reg <- brma(yi ~ mod, sei = sei, data = dat)
#' funnel(fit_reg) # shows standardized residuals
#'
#' # Different residual types
#' funnel(fit_reg, type = "LOO-PIT")
#'
#' # Customize appearance
#' funnel(fit, pch = 19, col = "blue", bg = "lightblue")
#' funnel(fit, back = "lightgrey", shade = "white", lty = "dashed")
#'
#' # using ggplot2
#' funnel(fit, plot_type = "ggplot")
#' }
#'
#' @seealso [residuals.brma()], [rstandard.brma()], [predict.brma()]
#' @export
#' @rdname funnel
funnel.brma <- function(x, residual, type = "LOO-PIT",
                        sampling_heterogeneity = TRUE, sampling_bias = TRUE,
                        plot_type = "base", ...) {
  # input validation
  BayesTools::check_char(type, "type", allow_values = c("outcome", "rstandard", "LOO-PIT", "rstudent"))
  BayesTools::check_bool(sampling_heterogeneity, "sampling_heterogeneity")
  BayesTools::check_bool(sampling_bias, "sampling_bias")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))

  # set up graphical arguments with defaults
  dots <- .set_dots_funnel(...)

  # get model characteristics
  is_mods  <- .is_mods(x)
  is_scale <- .is_scale(x)

  # determine mode: residual mode if mods/scale present, or if explicitly requested

  if (missing(residual)) {
    is_residual <- is_mods || is_scale
  } else {
    BayesTools::check_bool(residual, "residual")
    is_residual <- residual
  }

  # generate funnel data based on mode
  if (is_residual) {
    # sampling heterogeneity/bias ignored in residual mode - explicitly set to FALSE
    sampling_heterogeneity <- FALSE
    sampling_bias          <- FALSE
    funnel_data <- .funnel_data_residual(x = x, type = type, dots = dots)
  } else {
    funnel_data <- .funnel_data_outcome(
      x                      = x,
      sampling_heterogeneity = sampling_heterogeneity,
      sampling_bias          = sampling_bias,
      dots                   = dots
    )
  }

  # allow data return for programmatic access
  if (isTRUE(dots[["as_data"]])) {
    return(funnel_data)
  }

  # create plot
  if (plot_type == "ggplot") {
    return(.funnel_plot_ggplot(funnel_data, dots))
  } else {
    .funnel_plot_base(funnel_data, dots)
    return(invisible(NULL))
  }
}


# ---------------------------------------------------------------------------- #
# .funnel_data_outcome
# ---------------------------------------------------------------------------- #
#
# Generate funnel plot data for outcome mode (no mods AND no scale).
#
# Shows observed residuals (yi - pooled estimate) against the sampling
# distribution funnel.
#
# @param x                      brma object
# @param sampling_heterogeneity logical; incorporate tau into funnel
# @param sampling_bias          logical; incorporate bias into funnel
# @param dots                   list of graphical parameters
#
# @return list with funnel plot data components
#
# ---------------------------------------------------------------------------- #
.funnel_data_outcome <- function(x, sampling_heterogeneity, sampling_bias, dots) {
  # get model characteristics
  is_multilevel     <- .is_multilevel(x)
  is_weightfunction <- .is_weightfunction(x)
  is_PET            <- .is_PET(x)
  is_PEESE          <- .is_PEESE(x)
  effect_direction  <- .effect_direction(x)

  # get observed effect sizes (yi) - these go on the x-axis
  yi <- .outcome_data_yi(x)

  # get pooled effect estimate - this is the funnel center
  mu      <- pooled_effect(x)
  mu_mean <- summary(mu)["mu", "Mean"]

  # get standard errors
  se <- .outcome_data_sei(x)
  K  <- length(se)

  # get tau for incorporating heterogeneity into sampling distribution
  if (sampling_heterogeneity) {
    tau <- .get_funnel_tau(x, is_multilevel)
  } else {
    tau <- 0
  }

  # compute funnel region
  se_range    <- pretty(c(0, max(se)))
  se_sequence <- seq(0, max(se_range), length.out = 21)

  # compute confidence interval offsets from center (based on sampling distribution)
  ci_offsets <- .get_funnel_quantiles(
    x                 = x,
    se_sequence       = se_sequence,
    mu                = mu_mean,
    tau               = tau,
    sampling_bias     = sampling_bias,
    is_weightfunction = is_weightfunction,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction
  )
  # funnel bounds are absolute coordinates (mu was passed inside)
  ci_left  <- ci_offsets$lower
  ci_right <- ci_offsets$upper

  # create data frames for plotting
  x_range <- pretty(c(ci_left, ci_right, yi))

  # apply custom axis limits if provided
  if (!is.null(dots[["xlim"]])) {
    x_range <- pretty(dots[["xlim"]])
  }
  if (!is.null(dots[["ylim"]])) {
    se_range <- pretty(dots[["ylim"]])
  }

  # set axis labels - using "Observed Effect Size" for outcome mode
  xlab <- if (!is.null(dots[["xlab"]])) dots[["xlab"]] else "Observed Effect Size"
  ylab <- if (!is.null(dots[["ylab"]])) dots[["ylab"]] else "Standard Error"

  # set data for reference line
  if (!is.null(dots[["refline"]])) {
    df_refline <- data.frame(
      x = rep(dots[["refline"]], 2),
      y = range(se_range)
    )
  } else {
    # use model center
    df_refline <- NULL
  }

  # determine plot limits
  if (!is.null(dots[["xlim"]])) xlim <- dots[["xlim"]] else xlim <- range(x_range)
  if (!is.null(dots[["ylim"]])) ylim <- dots[["ylim"]] else ylim <- range(se_range)

  # recompute funnel region STRICTLY within ylim to avoid artifacts
  # generate sequence within ylim
  se_sequence_clipped <- seq(min(ylim), max(ylim), length.out = 100)

  # recompute contours for clipped sequence
  ci_offsets_clipped <- .get_funnel_quantiles(
    x                 = x,
    se_sequence       = se_sequence_clipped,
    mu                = mu_mean,
    tau               = tau,
    sampling_bias     = sampling_bias,
    is_weightfunction = is_weightfunction,
    is_PET            = is_PET,
    is_PEESE          = is_PEESE,
    effect_direction  = effect_direction
  )
  ci_left  <- ci_offsets_clipped$lower
  ci_right <- ci_offsets_clipped$upper

  # Clamp CIs to xlim to prevent drawing outside plot area (for POLYGON)
  ci_left_clipped  <- pmax(ci_left, min(xlim))
  ci_right_clipped <- pmin(ci_right, max(xlim))

  # Compute reference line for plotting (if not user-specified)
  if (is.null(df_refline)) {
      df_refline <- data.frame(
        x = ci_offsets_clipped$mid,
        y = se_sequence_clipped
      )
      # Clip refline to xlim just in case
      df_refline <- .clip_line_x(df_refline$x, df_refline$y, xlim)
  }

  # create output data structures
  # Points: observed data
  df_points <- data.frame(
    x = yi,
    y = se
  )
  # Funnel: use CLIPPED data for filled polygon to respect limits
  df_funnel <- data.frame(
    x = c(rev(ci_left_clipped), ci_right_clipped),
    y = c(rev(se_sequence_clipped), se_sequence_clipped)
  )
  # Edges: use CLIPPED lines using exact intersection
  edge1 <- .clip_line_x(ci_left, se_sequence_clipped, xlim)
  edge2 <- .clip_line_x(ci_right, se_sequence_clipped, xlim)

  df_funnel_edge1 <- data.frame(
    x = edge1$x,
    y = edge1$y
  )
  df_funnel_edge2 <- data.frame(
    x = edge2$x,
    y = edge2$y
  )
  df_background <- data.frame(
    x = c(min(xlim), max(xlim), max(xlim), min(xlim)),
    y = c(min(ylim), min(ylim), max(ylim), max(ylim))
  )

  return(list(
    points       = df_points,
    funnel       = df_funnel,
    funnel_edge1 = df_funnel_edge1,
    funnel_edge2 = df_funnel_edge2,
    background   = df_background,
    x_range      = x_range, # Keep original ticks for axis drawing if needed
    y_range      = se_range,
    xlim         = xlim,
    ylim         = ylim,
    xlab         = xlab,
    ylab         = ylab,
    refline      = df_refline
  ))
}


# ---------------------------------------------------------------------------- #
# .funnel_data_residual
# ---------------------------------------------------------------------------- #
#
# Generate funnel plot data for residual mode (mods OR scale).
#
# @param x    brma object
# @param type character; type of residuals
# @param dots list of graphical parameters
#
# @return list with funnel plot data components
#
# ---------------------------------------------------------------------------- #
.funnel_data_residual <- function(x, type, dots) {

  # get residuals based on type
  # rstandard and rstudent return data.frame with resid, se, z columns
  if (type == "rstandard") {
    res_obj <- rstandard.brma(x)
    res <- res_obj$resid
    se  <- res_obj$se
  } else if (type == "LOO-PIT" || type == "rstudent") {
    res_obj <- rstudent.brma(x)
    res <- res_obj$resid
    se  <- res_obj$se
  } else if (type == "outcome") {
    # raw outcome residuals
    res <- residuals.brma(x, type = "outcome")
    se  <- .outcome_data_sei(x)
  }
  K <- length(se)

  # compute funnel region - NO tau incorporation for residual funnel
  # for residuals, the expected distribution is N(0, se^2)
  # so bounds are +/- 1.96 * se
  # se_range determines the axis ticks
  se_range    <- pretty(c(0, max(se)))

  # determine plot limits
  # For xlim: pretty range of residuals + quantiles
  # For ylim: defaults to se_range range
  x_init_range <- pretty(c(stats::qnorm(0.025) * max(se_range), stats::qnorm(0.975) * max(se_range), res))

  if (!is.null(dots[["xlim"]])) xlim <- dots[["xlim"]] else xlim <- range(x_init_range)
  if (!is.null(dots[["ylim"]])) ylim <- dots[["ylim"]] else ylim <- range(se_range)

  # generate sequence strictly within ylim for clean polygons
  se_sequence_clipped <- seq(min(ylim), max(ylim), length.out = 100)

  # funnel bounds are quantiles of N(0, se^2) -> qnorm(p) * se
  ci_left  <- stats::qnorm(0.025) * se_sequence_clipped
  ci_right <- stats::qnorm(0.975) * se_sequence_clipped

  # Clamp CIs to xlim to prevent drawing outside plot area (for POLYGON)
  ci_left_clipped  <- pmax(ci_left, min(xlim))
  ci_right_clipped <- pmin(ci_right, max(xlim))

  # set axis labels
  xlab <- if (!is.null(dots[["xlab"]])) dots[["xlab"]] else "Residual"
  ylab <- if (!is.null(dots[["ylab"]])) dots[["ylab"]] else "Standard Error"

  # set data for reference line (residuals centered at 0 or user specified)
  if (!is.null(dots[["refline"]])) {
    refline_x <- dots[["refline"]]
  } else {
    refline_x <- 0
  }
  df_refline <- data.frame(
    x = rep(refline_x, 2),
    y = range(se_range)
  )

  # create output data structures
  df_points <- data.frame(
    x = res,
    y = se
  )
  # Funnel: use CLIPPED data for filled polygon
  df_funnel <- data.frame(
    x = c(rev(ci_left_clipped), ci_right_clipped),
    y = c(rev(se_sequence_clipped), se_sequence_clipped)
  )
  # Edges: use CLIPPED lines using exact intersection
  edge1 <- .clip_line_x(ci_left, se_sequence_clipped, xlim)
  edge2 <- .clip_line_x(ci_right, se_sequence_clipped, xlim)

  df_funnel_edge1 <- data.frame(
    x = edge1$x,
    y = edge1$y
  )
  df_funnel_edge2 <- data.frame(
    x = edge2$x,
    y = edge2$y
  )
  df_background <- data.frame(
    x = c(min(xlim), max(xlim), max(xlim), min(xlim)),
    y = c(min(ylim), min(ylim), max(ylim), max(ylim))
  )

  return(list(
    points       = df_points,
    funnel       = df_funnel,
    funnel_edge1 = df_funnel_edge1,
    funnel_edge2 = df_funnel_edge2,
    background   = df_background,
    x_range      = x_init_range, # ticks
    y_range      = se_range,     # ticks
    xlim         = xlim,
    ylim         = ylim,
    xlab         = xlab,
    ylab         = ylab,
    refline      = df_refline
  ))
}


# ---------------------------------------------------------------------------- #
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
    axes = FALSE
  )
   graphics::axis(1, at = x_ticks)
   graphics::axis(2, at = rev(y_ticks), las = 1)

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

  # add vertical reference line (now potentially curved)
  out <- out +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = df_refline$x,
        y = df_refline$y
      ),
      linetype = lty,
      colour = col_refline
    )

  # add funnel edge lines
  out <- out +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = df_funnel_edge1$x,
        y = df_funnel_edge1$y
      ),
      linetype = lty,
      colour = col_line
    ) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = df_funnel_edge2$x,
        y = df_funnel_edge2$y
      ),
      linetype = lty,
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

    # check intersection with xmin
    if ((x1 < xmin && x2 > xmin) || (x1 > xmin && x2 < xmin)) {
      y_int <- y1 + (y2 - y1) * (xmin - x1) / (x2 - x1)
      x_out <- c(x_out, xmin)
      y_out <- c(y_out, y_int)
    }

    # check intersection with xmax
    if ((x1 < xmax && x2 > xmax) || (x1 > xmax && x2 < xmax)) {
      y_int <- y1 + (y2 - y1) * (xmax - x1) / (x2 - x1)
      x_out <- c(x_out, xmax)
      y_out <- c(y_out, y_int)
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


# ---------------------------------------------------------------------------- #
# .get_funnel_tau
# ---------------------------------------------------------------------------- #
#
# Extract mean tau (heterogeneity) estimate from brma object for funnel plot.
#
# For multilevel models, returns total tau (combining within and between components).
#
# @param object       brma object
# @param is_multilevel logical; whether the model is multilevel
#
# @return numeric scalar representing the mean tau estimate
#
# ---------------------------------------------------------------------------- #
.get_funnel_tau <- function(object, is_multilevel) {
  # use pooled_heterogeneity to get mean tau
  tau_samples <- pooled_heterogeneity(object)
  tau_summary <- summary(tau_samples)
  return(tau_summary["tau", "Mean"])
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles
# ---------------------------------------------------------------------------- #
#
# Compute quantiles for funnel plot contours based on the sampling distribution.
#
# When sampling_bias = TRUE:
# - For selection models (weightfunction): uses weighted normal quantiles
# - For PET/PEESE: incorporates expected skew from regression adjustment
#
# When sampling_bias = FALSE: uses standard normal quantiles.
#
# @param x                 brma object
# @param se_sequence       numeric vector of SE values for funnel contours
# @param mu                numeric scalar of pooled effect estimate
# @param tau               numeric scalar of heterogeneity estimate
# @param sampling_bias     logical; whether to incorporate bias into sampling dist
# @param is_weightfunction logical; whether model has selection model
# @param is_PET            logical; whether model has PET adjustment
# @param is_PEESE          logical; whether model has PEESE adjustment
# @param effect_direction  character; "positive" or "negative"
#
# @return list with 'lower' and 'upper' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles <- function(x, se_sequence, mu, tau, sampling_bias,
                                  is_weightfunction, is_PET, is_PEESE,
                                  effect_direction) {
  n_se   <- length(se_sequence)
  sd_seq <- sqrt(se_sequence^2 + tau^2)

  # default: standard normal quantiles centered at mu
  if (!sampling_bias || (!is_weightfunction && !is_PET && !is_PEESE)) {
    return(list(
      lower = stats::qnorm(0.025, mean = mu, sd = sd_seq),
      upper = stats::qnorm(0.975, mean = mu, sd = sd_seq),
      mid   = rep(mu, n_se)
    ))
  } else if (is_weightfunction) {
    # use weighted normal quantiles for selection models
    return(.get_funnel_quantiles_weighted(x, se_sequence, mu, sd_seq, effect_direction))
  } else if (is_PET || is_PEESE) {
    # for PET/PEESE, incorporate regression-based skew
    return(.get_funnel_quantiles_PETPEESE(x, se_sequence, mu, sd_seq, is_PET, is_PEESE, effect_direction))
  }
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles_weighted
# ---------------------------------------------------------------------------- #
#
# Compute weighted normal quantiles for funnel contours when model has
# selection (weightfunction) adjustment.
#
# Uses posterior mean omega and computes weighted quantiles via .qwnorm_fast.ss.
#
# @param x                brma object
# @param se_sequence      numeric vector of SE values
# @param mu               numeric scalar of pooled effect
# @param sd_seq           numeric vector of total SD (sqrt(se^2 + tau^2))
# @param effect_direction character; "positive" or "negative"
#
# @return list with 'lower' and 'upper' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles_weighted <- function(x, se_sequence, mu, sd_seq, effect_direction) {

  n_se <- length(se_sequence)

  # extract posterior mean omega from fitted model
  posterior_samples <- suppressWarnings(coda::as.mcmc(x[["fit"]]))
  # use more robust regex for extracting omega columns (e.g. omega[1], omega[2])
  omega_cols <- grep("^omega\\[", colnames(posterior_samples))

  omega_mean   <- colMeans(posterior_samples[, omega_cols, drop = FALSE])
  omega_matrix <- matrix(omega_mean, nrow = n_se, ncol = length(omega_mean), byrow = TRUE)

  # Extract publication bias priors to determine steps (p-value cutoffs)
  priors <- x[["priors"]]
  priors_bias <- priors[["outcome"]][["bias"]]
  if (!BayesTools::is.prior.mixture(priors_bias)) {
    priors_bias <- list(priors_bias)
  }

  # create the weightfunction mapping for effect size thresholds
  # (logic mirrored from .create_fit_data in fit.R)
  steps <- BayesTools::weightfunctions_mapping(
    priors_bias[sapply(priors_bias, BayesTools::is.prior.weightfunction)],
    cuts_only = TRUE,
    one_sided = TRUE
  )
  # remove 0 and 1 boundaries if present and reversed
  steps <- rev(steps)[c(-1, -length(steps))]

  lower <- numeric(n_se)
  upper <- numeric(n_se)

  # Adjust mu based on effect direction for calculation
  # Weighted normal calculation assumes positive significant effects
  # If effect_direction is negative, the "significant" effects are negative values
  if (effect_direction == "negative") {
    mu_calc <- -mu
  } else {
    mu_calc <- mu
  }

  # compute weighted quantiles for each SE level
  for (i in seq_len(n_se)) {

    # Calculate critical values for the current SE using the steps
    crit_x <- stats::qnorm(steps, lower.tail = FALSE) * se_sequence[i]

    # compute quantiles centered at mu_calc
    lower[i] <- .qwnorm_fast.ss(
      p      = 0.025,
      mean   = mu_calc,
      sd     = sd_seq[i],
      omega  = omega_matrix[i, , drop = FALSE],
      crit_x = crit_x
    )
    upper[i] <- .qwnorm_fast.ss(
      p      = 0.975,
      mean   = mu_calc,
      sd     = sd_seq[i],
      omega  = omega_matrix[i, , drop = FALSE],
      crit_x = crit_x
    )
  }

  # flip back if effect direction is negative
  if (effect_direction == "negative") {
    temp  <- -upper
    upper <- -lower
    lower <- temp
  }

  return(list(lower = lower, upper = upper, mid = rep(mu, n_se)))
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles_PETPEESE
# ---------------------------------------------------------------------------- #
#
# Compute quantiles incorporating PET/PEESE regression-based skew.
#
# PET adds a linear term in SE, PEESE adds a quadratic term in SE^2.
# This shifts the expected sampling distribution away from center.
#
# @param x                brma object
# @param se_sequence      numeric vector of SE values
# @param mu               numeric scalar of pooled effect
# @param sd_seq           numeric vector of total SD (sqrt(se^2 + tau^2))
# @param is_PET           logical; whether model has PET
# @param is_PEESE         logical; whether model has PEESE
# @param effect_direction character; "positive" or "negative"
#
# @return list with 'lower' and 'upper' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles_PETPEESE <- function(x, se_sequence, mu, sd_seq, is_PET, is_PEESE, effect_direction) {

  n_se <- length(se_sequence)

  # extract posterior mean PET/PEESE coefficients
  posterior_samples <- suppressWarnings(coda::as.mcmc(x[["fit"]]))

  # compute bias shift at each SE level
  bias_shift <- rep(0, n_se)
  direction <- ifelse(effect_direction == "negative", -1, 1)

  if (is_PET && "PET" %in% colnames(posterior_samples)) {
    PET_mean   <- mean(posterior_samples[, "PET"])
    bias_shift <- bias_shift + direction * PET_mean * se_sequence
  }

  if (is_PEESE && "PEESE" %in% colnames(posterior_samples)) {
    PEESE_mean <- mean(posterior_samples[, "PEESE"])
    bias_shift <- bias_shift + direction * PEESE_mean * se_sequence^2
  }

  # the funnel is centered at mu, plus bias shifts
  mid   <- mu + bias_shift
  lower <- stats::qnorm(0.025, mean = mid, sd = sd_seq)
  upper <- stats::qnorm(0.975, mean = mid, sd = sd_seq)

  return(list(lower = lower, upper = upper, mid = mid))
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
.set_dots_funnel <- function(...) {

  dots <- list(...)

  # point styling (use metafor-style argument names)
  if (is.null(dots[["pch"]])) dots[["pch"]] <- 21 # open circle with fill
  if (is.null(dots[["col"]])) dots[["col"]] <- "black" # point border color
  if (is.null(dots[["bg"]])) dots[["bg"]] <- "black" # point fill color
  if (is.null(dots[["cex"]])) dots[["cex"]] <- 1 # point size (base)
  if (is.null(dots[["size"]])) dots[["size"]] <- 2 # point size (ggplot)

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

  # axis limits (NULL = auto)
  # xlim and ylim are left as NULL if not provided (auto-computed)

  # data return flag
  if (is.null(dots[["as_data"]])) dots[["as_data"]] <- FALSE

  return(dots)
}
