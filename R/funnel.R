# ============================================================================ #
# brma.funnel.R
# ============================================================================ #
#
# This file contains the funnel plot functions for brma objects.
#
# Two modes of dispatch:
# 1) No mods AND no scale: Shows observed effect sizes against the fitted
#    sampling distribution
# 2) All other cases: Shows residuals against a standard-error funnel
#
# ============================================================================ #


### funnel plot functions ----
#' @export
funnel <- function(x, ...) UseMethod("funnel")


#' @title Funnel Plot for brma Object
#'
#' @description \code{funnel.brma} creates a funnel plot for a fitted brma object.
#' For intercept-only models without scale regression, the default outcome mode
#' displays observed effect sizes against the fitted sampling distribution. For
#' models with location or scale moderators, the default residual mode displays
#' residuals against a standard-error funnel.
#'
#' @param x a fitted brma object
#' @param residual whether to use residual mode. Defaults to not specified,
#' which means the function automatically determines the mode:
#' \itemize{
#'   \item \code{FALSE} (auto): For intercept-only models without scale
#'     regression, displays observed effect sizes against the fitted sampling
#'     distribution funnel.
#'   \item \code{TRUE} (auto or explicit): For models with moderators or scale
#'     regression, or when explicitly set to \code{TRUE}, displays residuals on
#'     the x-axis, using \code{type} to determine how those residuals are
#'     computed.
#' }
#' @param type the type of residuals to use when in residual mode.
#' Options are:
#' \itemize{
#'   \item \code{"LOO-PIT"} (alias: \code{"rstudent"}; default): Leave-one-out
#'     probability integral transform residuals returned by
#'     \code{\link{rstudent.brma}}.
#'   \item \code{"rstandard"}: Internally standardized residuals using
#'     \code{\link{rstandard.brma}}. Only available for normal outcome models.
#'   \item \code{"outcome"}: Raw outcome residuals from
#'     \code{\link{residuals.brma}} with \code{type = "outcome"}.
#' }
#' Only used when funnel is in residual mode.
#' @param unit output unit for residual mode. Only \code{"estimate"} is
#' implemented in this pass.
#' @param conditioning_depth residual conditioning depth for residual mode.
#' Options are \code{"marginal"}, \code{"cluster"}, and \code{"estimate"}.
#' @param sampling_heterogeneity whether heterogeneity should be incorporated
#' into the sampling distribution funnel. Defaults to \code{TRUE}. Only used
#' in outcome mode and ignored in residual mode.
#' @param sampling_bias whether publication bias should be incorporated into the
#' sampling distribution funnel. Defaults to \code{TRUE}. Only used when
#' \code{residual = FALSE} or when automatic mode selects outcome mode. Ignored
#' in residual mode. When \code{TRUE} and the model
#' includes selection models (weightfunction), uses selected-normal quantiles.
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
#'   \item{bg}{point fill/background color (default: "#A6A6A6")}
#'   \item{cex}{point size multiplier for base graphics (default: 1)}
#'   \item{size}{point size for ggplot2 (default: 2)}
#'   \item{back}{background region color (default: "grey"). Set to \code{NA} to suppress.}
#'   \item{shade}{funnel region fill color (default: "white"). Set to \code{NA} to suppress.}
#'   \item{lty}{line type for funnel edges and center line (default: "dotted")}
#'   \item{col.line}{color for funnel edge lines (default: "black")}
#'   \item{refline}{numeric override for the reference line. By default,
#'   residual mode uses 0, while outcome mode uses the center of the fitted
#'   sampling distribution, which may be curved when PET/PEESE bias adjustment
#'   is incorporated.}
#'   \item{col.refline}{color of vertical reference line (default: "black")}
#'   \item{as_data}{if \code{TRUE}, returns plot data instead of creating a plot}
#' }
#'
#' @details
#' The funnel plot has two modes. If \code{residual} is not specified, the mode
#' is chosen automatically from the fitted model: intercept-only models without
#' scale regression use outcome mode, whereas models with location or scale
#' moderators use residual mode.
#'
#' \strong{Outcome mode} (intercept-only models without scale regression):
#' Displays observed effect sizes on the x-axis and standard errors on the
#' y-axis. The reference line follows the center of the fitted sampling
#' distribution. When \code{sampling_bias = FALSE}, this center is the pooled
#' effect; when PET/PEESE bias adjustment is incorporated, the center line can
#' vary with the standard error. The funnel region represents the central 95\%
#' region of the sampling distribution, optionally incorporating heterogeneity
#' and publication bias.
#'
#'
#' \strong{Residual mode} (models with moderators or scale regression):
#' Displays residuals on the x-axis and standard errors on the
#' y-axis. The funnel region represents the central 95\% region of
#' \eqn{N(0, \mathrm{SE}^2)}. With \code{type = "LOO-PIT"}, the plotted
#' residuals and standard errors are the raw-scale LOO predictive companions
#' returned by \code{\link{rstudent.brma}}; the PIT-normalized \code{z} values
#' are used by \code{\link{qqnorm.brma}} and influence diagnostics. With
#' \code{type = "rstandard"}, the plotted values are internally standardized
#' residual companions from \code{\link{rstandard.brma}}. With
#' \code{type = "outcome"}, these are raw outcome residuals. Under a correctly
#' specified model, most points should fall within this region.
#'
#' The \code{type} argument controls how residuals are computed in residual
#' mode. See \code{\link{residuals.brma}} for details on each type.
#' The \code{sampling_heterogeneity} and \code{sampling_bias} arguments are
#' ignored in residual mode.
#'
#' For GLMM models, observed effect sizes are computed from the raw frequency
#' data using formulas equivalent to \code{metafor::escalc}. Residual-mode
#' GLMM funnels use approximate effect-size-scale residual/PIT companions, not
#' exact PIT diagnostics for the raw count likelihood.
#'
#' @return If \code{as_data = TRUE}, \code{funnel.brma} returns a list with the
#' data used for plotting, including the plotted points, funnel polygons,
#' plotting limits, labels, and reference line. Otherwise, it returns
#' \code{NULL} invisibly if \code{plot_type = "base"} or a ggplot object if
#' \code{plot_type = "ggplot"}.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE) &&
#'     requireNamespace("metafor", quietly = TRUE)) {
#'   data(dat.bcg, package = "metadat")
#'   dat <- metafor::escalc(
#'     measure = "RR",
#'     ai      = tpos,
#'     bi      = tneg,
#'     ci      = cpos,
#'     di      = cneg,
#'     data    = dat.bcg
#'   )
#'
#'   fit <- brma(yi = yi, vi = vi, data = dat, measure = "RR")
#'   funnel(fit)
#'   funnel(fit, pch = 19, col = "blue", bg = "lightblue")
#'
#'   fit_reg <- brma(
#'     yi      = yi,
#'     vi      = vi,
#'     mods    = ~ ablat + year,
#'     data    = dat,
#'     measure = "RR"
#'   )
#'   fit_reg <- add_loo(fit_reg)
#'   funnel(fit_reg)
#'   funnel(fit_reg, type = "outcome")
#'
#'   funnel_data <- funnel(fit, as_data = TRUE)
#'   funnel(fit, plot_type = "ggplot")
#' }
#' }
#'
#' @seealso [residuals.brma()], [rstandard.brma()], [rstudent.brma()],
#'   [predict.brma()]
#' @aliases funnel
#' @export
#' @rdname funnel
funnel.brma <- function(x, residual, type = "LOO-PIT",
                        unit = "estimate", conditioning_depth = "marginal",
                        sampling_heterogeneity = TRUE, sampling_bias = TRUE,
                        plot_type = "base", ...) {
  # input validation
  conditioning_depth_specified <- !missing(conditioning_depth)
  dots                         <- list(...)
  .check_legacy_level_arg(dots, "funnel()")

  BayesTools::check_bool(sampling_heterogeneity, "sampling_heterogeneity")
  BayesTools::check_bool(sampling_bias, "sampling_bias")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))

  # set up graphical arguments with defaults
  dots <- .set_dots_funnel(dots)

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
    BayesTools::check_char(type, "type", allow_values = c("outcome", "rstandard", "LOO-PIT", "rstudent"))
    unit               <- .normalize_unit(unit)
    conditioning_depth <- .normalize_conditioning_depth(conditioning_depth)
    .check_unit_conditioning_depth(
      object             = x,
      unit               = unit,
      conditioning_depth = conditioning_depth,
      caller             = "funnel()"
    )

    if (unit == "cluster") {
      .check_cluster_unit_deferred("funnel()")
    }

    if ((type == "LOO-PIT" || type == "rstudent") && conditioning_depth_specified) {
      stop(
        "LOO-PIT residuals use the estimate-unit LOO target; ",
        "do not set 'conditioning_depth'.",
        call. = FALSE
      )
    }

    # sampling heterogeneity/bias ignored in residual mode - explicitly set to FALSE
    sampling_heterogeneity <- FALSE
    sampling_bias          <- FALSE
    funnel_data <- .funnel_data_residual(
      x                  = x,
      type               = type,
      unit               = unit,
      conditioning_depth = conditioning_depth,
      dots               = dots
    )
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
# Shows observed effect sizes against the fitted sampling distribution funnel.
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
    tau <- .get_funnel_tau(x)
  } else {
    tau <- 0
  }

  # compute funnel region
  se_range    <- pretty(c(0, max(se)))
  se_sequence <- seq(0, max(se_range), length.out = 21)

  # compute confidence interval offsets from center (based on sampling distribution)
  ci_offsets <- .get_funnel_quantiles(
    x                      = x,
    se_sequence            = se_sequence,
    mu                     = mu_mean,
    tau                    = tau,
    sampling_heterogeneity = sampling_heterogeneity,
    sampling_bias          = sampling_bias,
    is_weightfunction      = is_weightfunction,
    is_PET                 = is_PET,
    is_PEESE               = is_PEESE,
    effect_direction       = effect_direction
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
    x                      = x,
    se_sequence            = se_sequence_clipped,
    mu                     = mu_mean,
    tau                    = tau,
    sampling_heterogeneity = sampling_heterogeneity,
    sampling_bias          = sampling_bias,
    is_weightfunction      = is_weightfunction,
    is_PET                 = is_PET,
    is_PEESE               = is_PEESE,
    effect_direction       = effect_direction
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
# @param x                  brma object
# @param type               character; type of residuals
# @param unit               character; output unit
# @param conditioning_depth character; conditioning depth
# @param dots               list of graphical parameters
#
# @return list with funnel plot data components
#
# ---------------------------------------------------------------------------- #
.funnel_data_residual <- function(x, type, unit, conditioning_depth, dots) {

  # get residuals based on type
  # rstandard and rstudent return data.frame with resid, se, z columns
  if (type == "rstandard") {
    res_obj <- rstandard.brma(
      model              = x,
      unit               = unit,
      conditioning_depth = conditioning_depth
    )
    res <- res_obj$resid
    se  <- res_obj$se
  } else if (type == "LOO-PIT" || type == "rstudent") {
    res_obj <- rstudent.brma(x, unit = unit)
    res <- res_obj$resid
    se  <- res_obj$se
  } else if (type == "outcome") {
    # raw outcome residuals
    res_obj <- residuals.brma(
      object             = x,
      type               = "outcome",
      unit               = unit,
      conditioning_depth = conditioning_depth
    )
    if (is.data.frame(res_obj)) {
      res <- res_obj$resid
      se  <- res_obj$se
    } else {
      res <- res_obj
      se  <- .outcome_data_sei(x)
    }
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
    x                  = res,
    y                  = se,
    unit               = rep(unit, K),
    conditioning_depth = rep(conditioning_depth, K)
  )
  if (exists("res_obj") && is.data.frame(res_obj) && "cluster" %in% names(res_obj)) {
    df_points[["cluster"]]     <- res_obj[["cluster"]]
    df_points[["n_estimates"]] <- res_obj[["n_estimates"]]
  }
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
# @return numeric scalar representing the mean tau estimate
#
# ---------------------------------------------------------------------------- #
.get_funnel_tau <- function(object) {
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
# When sampling_bias = TRUE, posterior rows dispatch to their active bias model
# and the contour is obtained by inverting the model-averaged CDF.
#
# When sampling_bias = FALSE: uses standard normal quantiles.
#
# @param x                      brma object
# @param se_sequence            numeric vector of SE values for funnel contours
# @param mu                     numeric scalar of pooled effect estimate
# @param tau                    numeric scalar of heterogeneity estimate
# @param sampling_heterogeneity logical; whether to incorporate tau
# @param sampling_bias          logical; whether to incorporate bias into sampling dist
# @param is_weightfunction      logical; whether model has selection model
# @param is_PET                 logical; whether model has PET adjustment
# @param is_PEESE               logical; whether model has PEESE adjustment
# @param effect_direction       character; "positive" or "negative"
#
# @return list with 'lower' and 'upper' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles <- function(x, se_sequence, mu, tau,
                                  sampling_heterogeneity, sampling_bias,
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
  }

  if (!BayesTools::is.prior.mixture(x[["priors"]][["outcome"]][["bias"]])) {
    if (is_PET || is_PEESE) {
      return(.get_funnel_quantiles_PETPEESE_plugin(
        x                = x,
        se_sequence      = se_sequence,
        mu               = mu,
        sd_seq           = sd_seq,
        is_PET           = is_PET,
        is_PEESE         = is_PEESE,
        effect_direction = effect_direction
      ))
    }
  }

  return(.get_funnel_quantiles_model_averaged(
    x                      = x,
    se_sequence            = se_sequence,
    sampling_heterogeneity = sampling_heterogeneity,
    effect_direction       = effect_direction
  ))
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles_PETPEESE_plugin
# ---------------------------------------------------------------------------- #
#
# Fast plug-in normal contours for single PET/PEESE priors.
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles_PETPEESE_plugin <- function(x, se_sequence, mu, sd_seq,
                                                  is_PET, is_PEESE,
                                                  effect_direction) {

  n_se              <- length(se_sequence)
  posterior_samples <- .get_posterior_samples(x[["fit"]])
  bias_shift        <- rep(0, n_se)
  direction         <- ifelse(effect_direction == "negative", -1, 1)

  if (is_PET && "PET" %in% colnames(posterior_samples)) {
    PET_mean   <- mean(posterior_samples[, "PET"])
    bias_shift <- bias_shift + direction * PET_mean * se_sequence
  }

  if (is_PEESE && "PEESE" %in% colnames(posterior_samples)) {
    PEESE_mean <- mean(posterior_samples[, "PEESE"])
    bias_shift <- bias_shift + direction * PEESE_mean * se_sequence^2
  }

  mid   <- mu + bias_shift
  lower <- stats::qnorm(0.025, mean = mid, sd = sd_seq)
  upper <- stats::qnorm(0.975, mean = mid, sd = sd_seq)

  return(list(lower = lower, upper = upper, mid = mid))
}


# ---------------------------------------------------------------------------- #
# .get_funnel_quantiles_model_averaged
# ---------------------------------------------------------------------------- #
#
# Compute model-averaged funnel contours for publication-bias models.
#
# Each posterior row dispatches to its active bias model:
# - no-bias rows use a normal CDF
# - PET/PEESE rows use a normal CDF with the row-specific bias offset
# - selection-model rows use the selected-normal CDF
#
# Quantiles are obtained from the averaged posterior-row CDF.
#
# @return list with 'lower', 'upper', and 'mid' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles_model_averaged <- function(x, se_sequence,
                                                 sampling_heterogeneity,
                                                 effect_direction) {

  setup <- .funnel_model_averaged_setup(
    x                      = x,
    sampling_heterogeneity = sampling_heterogeneity
  )

  lower <- vapply(
    se_sequence,
    .funnel_model_averaged_quantile,
    numeric(1),
    p                = 0.025,
    setup            = setup,
    effect_direction = effect_direction
  )
  upper <- vapply(
    se_sequence,
    .funnel_model_averaged_quantile,
    numeric(1),
    p                = 0.975,
    setup            = setup,
    effect_direction = effect_direction
  )
  mid <- vapply(
    se_sequence,
    .funnel_model_averaged_quantile,
    numeric(1),
    p                = 0.5,
    setup            = setup,
    effect_direction = effect_direction
  )

  return(list(lower = lower, upper = upper, mid = mid))
}


# ---------------------------------------------------------------------------- #
# .funnel_model_averaged_setup
# ---------------------------------------------------------------------------- #
#
# Prepare posterior-row quantities for model-averaged funnel contours.
#
# ---------------------------------------------------------------------------- #
.funnel_model_averaged_setup <- function(x, sampling_heterogeneity) {

  posterior_samples <- .get_posterior_samples(x[["fit"]])
  S                 <- nrow(posterior_samples)
  priors_bias       <- x[["priors"]][["outcome"]][["bias"]]

  if (!BayesTools::is.prior.mixture(priors_bias)) {
    priors_bias <- list(priors_bias)
  }

  branch_is_PET   <- vapply(priors_bias, BayesTools::is.prior.PET,   logical(1))
  branch_is_PEESE <- vapply(priors_bias, BayesTools::is.prior.PEESE, logical(1))
  bias_indicator  <- .extract_bias_indicator(x, posterior_samples = posterior_samples)

  if (any(is.na(bias_indicator)) ||
      any(bias_indicator < 1L | bias_indicator > length(priors_bias))) {
    stop("Invalid 'bias_indicator' values in posterior samples.", call. = FALSE)
  }

  mu_samples <- .funnel_mu_samples(x, posterior_samples)

  if (sampling_heterogeneity) {
    tau_samples <- .funnel_tau_samples(x, posterior_samples)
  } else {
    tau_samples <- rep(0, S)
  }

  PET_samples   <- .funnel_posterior_column(posterior_samples, "PET",   S)
  PEESE_samples <- .funnel_posterior_column(posterior_samples, "PEESE", S)

  PET_samples[!branch_is_PET[bias_indicator]]     <- 0
  PEESE_samples[!branch_is_PEESE[bias_indicator]] <- 0

  selection <- .selection_context(
    object            = x,
    posterior_samples = posterior_samples
  )
  use_normal <- if (is.null(selection)) {
    rep(TRUE, S)
  } else {
    selection[["use_normal"]]
  }

  return(list(
    mu                    = mu_samples,
    tau                   = tau_samples,
    PET                   = PET_samples,
    PEESE                 = PEESE_samples,
    bias_indicator        = bias_indicator,
    is_weightfunction     = !use_normal,
    selection             = selection
  ))
}


# ---------------------------------------------------------------------------- #
# .funnel_mu_samples
# ---------------------------------------------------------------------------- #
#
# Extract pooled location samples without PET/PEESE bias offsets.
#
# ---------------------------------------------------------------------------- #
.funnel_mu_samples <- function(x, posterior_samples) {

  if (!.is_mods(x)) {
    return(as.numeric(posterior_samples[, "mu"]))
  }

  mu_samples <- predict.brma(
    object        = x,
    newdata       = TRUE,
    type          = "terms",
    bias_adjusted = TRUE,
    quiet         = TRUE
  )

  return(as.numeric(as.matrix(mu_samples)[, 1]))
}


# ---------------------------------------------------------------------------- #
# .funnel_tau_samples
# ---------------------------------------------------------------------------- #
#
# Extract total heterogeneity samples for outcome-mode funnel contours.
#
# ---------------------------------------------------------------------------- #
.funnel_tau_samples <- function(x, posterior_samples) {

  tau_result <- .evaluate.brma.tau(
    fit               = x[["fit"]],
    scale_data        = x[["data"]][["scale"]],
    scale_formula     = if (.is_scale(x)) .create_fit_formula_list(data = x[["data"]], "scale") else NULL,
    scale_priors      = x[["priors"]][["scale"]],
    is_scale          = .is_scale(x),
    is_multilevel     = .is_multilevel(x),
    K                 = nrow(x[["data"]][["outcome"]]),
    posterior_samples = posterior_samples
  )

  return(rowMeans(tau_result[["tau_total"]]))
}


# ---------------------------------------------------------------------------- #
# .funnel_posterior_column
# ---------------------------------------------------------------------------- #
#
# Extract an optional posterior column, returning zeros when absent.
#
# ---------------------------------------------------------------------------- #
.funnel_posterior_column <- function(posterior_samples, column, S) {

  if (column %in% colnames(posterior_samples)) {
    return(as.numeric(posterior_samples[, column]))
  }

  return(rep(0, S))
}


# ---------------------------------------------------------------------------- #
# .funnel_model_averaged_quantile
# ---------------------------------------------------------------------------- #
#
# Invert the model-averaged posterior-row CDF at one SE value.
#
# ---------------------------------------------------------------------------- #
.funnel_model_averaged_quantile <- function(se, p, setup, effect_direction) {

  total_sd <- sqrt(se^2 + setup[["tau"]]^2)
  location <- .funnel_row_location(se, setup, effect_direction)
  eps_sd   <- sqrt(.Machine$double.eps)

  if (all(total_sd < eps_sd)) {
    return(unname(stats::quantile(location, probs = p, names = FALSE, type = 8)))
  }

  spread <- pmax(total_sd, eps_sd)
  lower  <- min(location - 10 * spread, na.rm = TRUE)
  upper  <- max(location + 10 * spread, na.rm = TRUE)

  if (!is.finite(lower) || !is.finite(upper)) {
    return(NA_real_)
  }
  if (lower >= upper) {
    lower <- lower - 1
    upper <- upper + 1
  }

  obj_fun <- function(q) {
    .funnel_model_averaged_cdf(
      q                = q,
      se               = se,
      setup            = setup,
      effect_direction = effect_direction
    ) - p
  }

  lower_value <- obj_fun(lower)
  upper_value <- obj_fun(upper)
  step        <- max(spread, na.rm = TRUE)
  if (!is.finite(step) || step <= 0) {
    step <- max(1, abs(location), na.rm = TRUE)
  }

  for (i in seq_len(25)) {
    if (lower_value <= 0 && upper_value >= 0) {
      break
    }
    if (lower_value > 0) {
      lower       <- lower - step
      lower_value <- obj_fun(lower)
    }
    if (upper_value < 0) {
      upper       <- upper + step
      upper_value <- obj_fun(upper)
    }
    step <- step * 2
  }

  if (lower_value > 0 || upper_value < 0) {
    return(.funnel_grid_quantile(p, lower, upper, se, setup, effect_direction))
  }

  out <- tryCatch(
    stats::uniroot(obj_fun, interval = c(lower, upper), tol = 1e-6)[["root"]],
    error = function(e) NA_real_
  )

  if (is.na(out)) {
    out <- .funnel_grid_quantile(p, lower, upper, se, setup, effect_direction)
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .funnel_model_averaged_cdf
# ---------------------------------------------------------------------------- #
#
# Average CDF values across posterior rows at one x-coordinate and SE.
#
# ---------------------------------------------------------------------------- #
.funnel_model_averaged_cdf <- function(q, se, setup, effect_direction) {

  S          <- length(setup[["mu"]])
  total_sd   <- sqrt(se^2 + setup[["tau"]]^2)
  location   <- .funnel_row_location(se, setup, effect_direction)
  eps_sd     <- sqrt(.Machine$double.eps)
  cdf_values <- rep(NA_real_, S)

  zero_sd <- total_sd < eps_sd
  if (any(zero_sd)) {
    cdf_values[zero_sd] <- as.numeric(q >= location[zero_sd])
  }

  normal_rows <- !setup[["is_weightfunction"]] & !zero_sd
  if (any(normal_rows)) {
    cdf_values[normal_rows] <- stats::pnorm(
      q,
      mean = location[normal_rows],
      sd   = total_sd[normal_rows]
    )
  }

  weighted_rows <- setup[["is_weightfunction"]] & !zero_sd
  if (any(weighted_rows)) {
    rows <- which(weighted_rows)
    cdf_values[rows] <- .funnel_selected_cdf(
      q                = q,
      rows             = rows,
      se               = se,
      total_sd         = total_sd,
      setup            = setup,
      effect_direction = effect_direction
    )
  }

  cdf_values <- pmin(pmax(cdf_values, 0), 1)
  return(mean(cdf_values))
}


# ---------------------------------------------------------------------------- #
# .funnel_row_location
# ---------------------------------------------------------------------------- #
#
# Row-specific expected observed location on the original effect-size scale.
#
# ---------------------------------------------------------------------------- #
.funnel_row_location <- function(se, setup, effect_direction) {

  direction <- ifelse(effect_direction == "negative", -1, 1)

  return(
    setup[["mu"]] +
      direction * setup[["PET"]] * se +
      direction * setup[["PEESE"]] * se^2
  )
}


# ---------------------------------------------------------------------------- #
# .funnel_selected_cdf
# ---------------------------------------------------------------------------- #
#
# Selected-normal CDF for selection-model posterior rows.
#
# ---------------------------------------------------------------------------- #
.funnel_selected_cdf <- function(q, rows, se, total_sd, setup, effect_direction) {

  selection <- setup[["selection"]]

  if (is.null(selection)) {
    return(stats::pnorm(q, mean = setup[["mu"]][rows], sd = total_sd[rows]))
  }

  .selection_require_step_evaluable(selection, ".funnel_selected_cdf()")

  if (se <= 0) {
    return(.funnel_selected_cdf_zero_se(
      q        = q,
      rows     = rows,
      total_sd = total_sd,
      setup    = setup
    ))
  }

  local_context <- .selection_context_subset_rows(selection, rows)
  cdf           <- .selection_step_cdf_matrix(
    q                 = q,
    mean              = matrix(setup[["mu"]][rows], ncol = 1),
    sd                = matrix(total_sd[rows], ncol = 1),
    sei               = se,
    selection_context = local_context,
    lower.tail        = TRUE
  )

  return(as.numeric(cdf[, 1]))
}


# ---------------------------------------------------------------------------- #
# .funnel_selected_cdf_zero_se
# ---------------------------------------------------------------------------- #
#
# Limit of the selected-normal CDF as the standard error approaches zero.
#
# Selection bins are defined on signed z = sign * y / se. At se -> 0, finite
# z-bin intervals collapse to zero mass. Only the two tail intervals remain:
# the upper z tail maps to signed effects above 0, and the lower z tail maps to
# signed effects below 0. This explicit limit avoids evaluating the native
# selected-normal kernel at the singular se = 0 boundary.
#
# ---------------------------------------------------------------------------- #
.funnel_selected_cdf_zero_se <- function(q, rows, total_sd, setup) {

  selection <- setup[["selection"]]
  sign      <- selection[["sign"]]
  omega     <- selection[["omega"]][rows, , drop = FALSE]
  mean      <- sign * setup[["mu"]][rows]
  sd        <- total_sd[rows]
  q_signed  <- sign * q
  S         <- length(rows)

  denom <- rep(0, S)
  mass  <- rep(0, S)

  for (b in seq_len(selection[["n_bins"]])) {
    z_lower <- selection[["z_lower"]][b]
    z_upper <- selection[["z_upper"]][b]

    if (is.infinite(z_lower) && z_lower < 0 &&
        is.infinite(z_upper) && z_upper > 0) {
      lower <- -Inf
      upper <-  Inf
    } else if (is.infinite(z_upper) && z_upper > 0) {
      lower <- 0
      upper <- Inf
    } else if (is.infinite(z_lower) && z_lower < 0) {
      lower <- -Inf
      upper <- 0
    } else {
      next
    }

    bin_mass <- .selection_interval_prob_vec(lower, upper, mean, sd)
    denom    <- denom + omega[, b] * bin_mass

    if (sign == 1L) {
      selected_mass <- .selection_interval_prob_vec(lower, min(upper, q_signed), mean, sd)
    } else {
      selected_mass <- .selection_interval_prob_vec(max(lower, q_signed), upper, mean, sd)
    }
    mass <- mass + omega[, b] * selected_mass
  }

  out <- mass / denom
  out <- pmin(pmax(out, 0), 1)
  return(out)
}


# ---------------------------------------------------------------------------- #
# .funnel_grid_quantile
# ---------------------------------------------------------------------------- #
#
# Fallback inversion for discontinuous CDFs from degenerate posterior rows.
#
# ---------------------------------------------------------------------------- #
.funnel_grid_quantile <- function(p, lower, upper, se, setup, effect_direction) {

  grid <- seq(lower, upper, length.out = 1000)
  cdf  <- vapply(
    grid,
    .funnel_model_averaged_cdf,
    numeric(1),
    se               = se,
    setup            = setup,
    effect_direction = effect_direction
  )
  index <- which(cdf >= p)[1]

  if (is.na(index)) {
    return(grid[length(grid)])
  }

  return(grid[index])
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
