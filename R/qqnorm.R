# ============================================================================ #
# brma.qqnorm.R
# ============================================================================ #
#
# Normal QQ plot functions for brma objects.
#
# The QQ plot displays sorted standardized residuals against theoretical
# standard normal quantiles. Under a correctly specified model, residuals
# should fall approximately along the y = x diagonal reference line.
#
# Supports pseudo confidence envelope via simulation (iid N(0,1) draws).
#
# ============================================================================ #


#' @title Normal QQ Plot for brma Object
#'
#' @description \code{qqnorm.brma} creates a normal QQ plot for a fitted brma
#' object. The plot displays sorted standardized residuals on the vertical axis
#' against the theoretical quantiles of a standard normal distribution on the
#' horizontal axis. Under a correctly specified model, the points should fall
#' approximately along the diagonal reference line.
#'
#' @param y a fitted brma object.
#' @param type the type of standardized residuals to use. Options are:
#' \itemize{
#'   \item \code{"rstudent"} (alias: \code{"LOO-PIT"}; default): Leave-one-out
#'     probability integral transform residuals. Works for all model types
#'     including GLMMs and selection models. This is the recommended type for
#'     Bayesian models as it properly accounts for estimation uncertainty and
#'     leverage. Note: requires that loo has been computed (see
#'     \code{\link{add_loo}}).
#'   \item \code{"rstandard"}: Internally standardized residuals using the hat
#'     matrix. Only available for normal outcome models without selection
#'     (weightfunction) bias adjustment.
#' }
#' @param envelope logical indicating whether to add a pseudo confidence
#' envelope to the plot. The envelope is created by simulating sets of
#' standard normal quantiles, providing a reference for the expected variability
#' under a correctly specified model. Defaults to \code{TRUE}.
#' @param level numeric value between 0 and 100 specifying the confidence level
#' for the envelope bounds. Defaults to \code{95}.
#' @param bonferroni logical indicating whether to apply Bonferroni correction
#' to the envelope, so it can be regarded as a simultaneous confidence region
#' for all \eqn{k} residuals. Defaults to \code{FALSE}.
#' @param reps integer specifying the number of simulation replications for
#' the envelope. Defaults to \code{1000}.
#' @param smooth logical indicating whether to smooth the envelope bounds using
#' Friedman's SuperSmoother (\code{\link[stats]{supsmu}}). Defaults to
#' \code{TRUE}.
#' @param xlim x-axis limits. If not specified, limits are computed from data.
#' @param ylim y-axis limits. If not specified, limits are computed from data.
#' @param xlab title for the x-axis. If not specified, defaults to
#' \code{"Theoretical Quantiles"}.
#' @param ylab title for the y-axis. If not specified, defaults to
#' \code{"Sample Quantiles"}.
#' @param plot_type whether to use a base plot \code{"base"} or ggplot2
#' \code{"ggplot"} for plotting. Defaults to \code{"base"}.
#' @param ... additional graphical arguments to customize the plot appearance:
#' \describe{
#'   \item{pch}{point symbol (default: 21, filled circle)}
#'   \item{col}{point border color (default: "black")}
#'   \item{bg}{point fill/background color (default: "black")}
#'   \item{cex}{point size multiplier for base graphics (default: 1)}
#'   \item{size}{point size for ggplot2 (default: 2)}
#'   \item{shade}{fill color for the envelope region (default: "grey90").
#'     Set to \code{NA} to suppress shading.}
#'   \item{col.envelope}{color for envelope boundary lines
#'     (default: "grey50")}
#'   \item{col.refline}{color of the diagonal reference line
#'     (default: "black")}
#'   \item{main}{character string for plot title (default: no title)}
#'   \item{as_data}{if \code{TRUE}, returns plot data instead of creating
#'     the plot}
#' }
#'
#' @details
#' The normal QQ plot is a diagnostic tool for assessing whether the
#' standardized residuals follow a standard normal distribution, as expected
#' under a correctly specified model. Each point represents one study,
#' plotted at coordinates \eqn{(\Phi^{-1}(p_i), z_{(i)})} where
#' \eqn{p_i} are the plotting positions and \eqn{z_{(i)}} are the sorted
#' standardized residuals.
#'
#' The default residual type is \code{"rstudent"} (LOO-PIT), which differs
#' from \code{metafor::qqnorm} (which defaults to \code{"rstandard"}).
#' LOO-PIT residuals are preferred because they are available for all model
#' types (including GLMMs and selection models) and properly account for
#' estimation uncertainty via leave-one-out cross-validation.
#'
#' The pseudo confidence envelope is constructed by simulation: for each
#' replication, \eqn{k} values are drawn from a standard normal distribution
#' and sorted. Pointwise quantile bounds are computed across replications.
#' When \code{smooth = TRUE}, bounds are smoothed using Friedman's
#' SuperSmoother. When \code{bonferroni = TRUE}, the confidence level is
#' adjusted for simultaneous inference across all \eqn{k} comparisons.
#'
#' @return For \code{plot_type = "base"}, returns an invisible list with
#' components \code{x} (theoretical quantiles) and \code{y} (sorted
#' residuals). For \code{plot_type = "ggplot"}, returns a ggplot object.
#'
#' If \code{as_data = TRUE}, returns a list with all computed plot data
#' including: \code{x}, \code{y}, \code{points}, \code{refline},
#' \code{envelope}, \code{xlim}, \code{ylim}, \code{xlab}, and \code{ylab}.
#'
#'
#' @examples \dontrun{
#' # Simple meta-analysis
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#' fit <- add_loo(fit)
#'
#' # Default QQ plot (rstudent/LOO-PIT residuals)
#' qqnorm(fit)
#'
#' # Using rstandard residuals (normal models only)
#' qqnorm(fit, type = "rstandard")
#'
#' # Without confidence envelope
#' qqnorm(fit, envelope = FALSE)
#'
#' # Using ggplot2
#' qqnorm(fit, plot_type = "ggplot")
#'
#' # Custom appearance
#' qqnorm(fit, pch = 19, col = "blue", shade = "lightblue")
#' }
#'
#' @seealso [rstudent.brma()], [rstandard.brma()], [residuals.brma()],
#' [funnel.brma()], [radial.brma()]
#' @exportS3Method
qqnorm.brma <- function(y, type = "rstudent", envelope = TRUE, level = 95,
                         bonferroni = FALSE, reps = 1000, smooth = TRUE,
                         xlim, ylim, xlab, ylab, plot_type = "base", ...) {

  # input validation
  BayesTools::check_char(type, "type", allow_values = c("rstandard", "rstudent", "LOO-PIT"))
  BayesTools::check_bool(envelope, "envelope")
  BayesTools::check_real(level, "level", lower = 0, upper = 100)
  BayesTools::check_bool(bonferroni, "bonferroni")
  BayesTools::check_int(reps, "reps", lower = 1)
  BayesTools::check_bool(smooth, "smooth")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))

  # set up graphical arguments with defaults
  dots <- .set_dots_qqnorm(...)

  # generate QQ plot data
  qq_data <- .qqnorm_data(
    x          = y,
    type       = type,
    envelope   = envelope,
    level      = level,
    bonferroni = bonferroni,
    reps       = reps,
    smooth     = smooth,
    xlim       = if (missing(xlim)) NULL else xlim,
    ylim       = if (missing(ylim)) NULL else ylim,
    xlab       = if (missing(xlab)) NULL else xlab,
    ylab       = if (missing(ylab)) NULL else ylab,
    dots       = dots
  )

  # allow data return for programmatic access
  if (isTRUE(dots[["as_data"]])) {
    return(qq_data)
  }

  # create plot
  if (plot_type == "ggplot") {
    return(.qqnorm_plot_ggplot(qq_data, dots))
  } else {
    .qqnorm_plot_base(qq_data, dots)
    return(invisible(qq_data[c("x", "y")]))
  }
}


# ---------------------------------------------------------------------------- #
# .qqnorm_data
# ---------------------------------------------------------------------------- #
#
# Generate data for normal QQ plot.
#
# Computes theoretical quantiles, sorted standardized residuals,
# reference line, and optional pseudo confidence envelope.
#
# @param x          brma object
# @param type       residual type ("rstudent", "LOO-PIT", or "rstandard")
# @param envelope   logical; draw envelope
# @param level      confidence level (0-100)
# @param bonferroni logical; Bonferroni correction
# @param reps       number of simulation replications
# @param smooth     logical; smooth envelope bounds
# @param xlim       x-axis limits (NULL for auto)
# @param ylim       y-axis limits (NULL for auto)
# @param xlab       x-axis label (NULL for default)
# @param ylab       y-axis label (NULL for default)
# @param dots       list of graphical parameters
#
# @return list with QQ plot data components
#
# ---------------------------------------------------------------------------- #
.qqnorm_data <- function(x, type, envelope, level, bonferroni, reps, smooth,
                          xlim, ylim, xlab, ylab, dots) {

  # get standardized residuals
  if (type == "rstandard") {
    res_obj <- rstandard.brma(x)
  } else {
    # type == "rstudent" or "LOO-PIT"
    res_obj <- rstudent.brma(x)
  }
  z <- res_obj$z
  K <- length(z)

  # theoretical standard normal quantiles
  qq_x <- stats::qnorm(stats::ppoints(K))

  # sorted standardized residuals
  sort_order <- order(z)
  qq_y       <- z[sort_order]

  # reference line: y = x (diagonal through origin)
  df_refline <- data.frame(
    x = range(qq_x),
    y = range(qq_x)
  )

  # compute envelope (if requested)
  env_data <- NULL
  if (envelope) {
    env_data <- .qqnorm_envelope(
      K          = K,
      level      = level,
      bonferroni = bonferroni,
      reps       = reps,
      smooth     = smooth,
      qq_x       = qq_x
    )
  }

  # axis limits
  if (is.null(xlim)) {
    xlim <- range(qq_x) * 1.1
  }
  if (is.null(ylim)) {
    all_y <- qq_y
    if (!is.null(env_data)) {
      all_y <- c(all_y, env_data$lower, env_data$upper)
    }
    ylim <- range(all_y) * 1.1
  }

  # axis labels
  if (is.null(xlab)) xlab <- "Theoretical Quantiles"
  if (is.null(ylab)) ylab <- "Sample Quantiles"

  # points data frame
  df_points <- data.frame(x = qq_x, y = qq_y)

  return(list(
    # core QQ data (for invisible return, matching metafor)
    x        = qq_x,
    y        = qq_y,
    # plot elements
    points   = df_points,
    refline  = df_refline,
    envelope = env_data,
    # limits and labels
    xlim     = xlim,
    ylim     = ylim,
    xlab     = xlab,
    ylab     = ylab,
    # metadata
    type     = type,
    level    = level,
    K        = K
  ))
}


# ---------------------------------------------------------------------------- #
# .qqnorm_envelope
# ---------------------------------------------------------------------------- #
#
# Simulate pseudo confidence envelope for the QQ plot.
#
# Under a correctly specified model, standardized residuals (especially
# LOO-PIT) should be approximately iid N(0,1). The envelope is constructed
# by simulating sorted standard normal samples and computing pointwise
# quantile bounds.
#
# @param K          number of observations
# @param level      confidence level (0-100)
# @param bonferroni logical; Bonferroni correction
# @param reps       number of simulation replications
# @param smooth     logical; smooth bounds with supsmu
# @param qq_x       theoretical quantiles
#
# @return data.frame with columns x, lower, upper
#
# ---------------------------------------------------------------------------- #
.qqnorm_envelope <- function(K, level, bonferroni, reps, smooth, qq_x) {

  # compute alpha (with optional Bonferroni correction)
  alpha <- 1 - level / 100
  if (bonferroni) {
    alpha <- alpha / K
  }

  # simulate sorted standard normal samples
  sim_matrix <- matrix(stats::rnorm(reps * K), nrow = reps, ncol = K)
  sim_matrix <- t(apply(sim_matrix, 1, sort))

  # pointwise quantile bounds
  lower <- apply(sim_matrix, 2, stats::quantile, probs = alpha / 2)
  upper <- apply(sim_matrix, 2, stats::quantile, probs = 1 - alpha / 2)

  # smooth bounds (skip if too few points)
  if (smooth && K >= 4) {
    lower <- stats::supsmu(qq_x, lower)$y
    upper <- stats::supsmu(qq_x, upper)$y
  }

  data.frame(
    x     = qq_x,
    lower = lower,
    upper = upper
  )
}


# ---------------------------------------------------------------------------- #
# .qqnorm_plot_base
# ---------------------------------------------------------------------------- #
#
# Create QQ normal plot using base R graphics.
#
# @param data list of QQ plot data from .qqnorm_data()
# @param dots list of graphical parameters
#
# @return NULL invisibly
#
# ---------------------------------------------------------------------------- #
.qqnorm_plot_base <- function(data, dots) {

  # extract graphical parameters
  pch  <- dots[["pch"]]
  col  <- dots[["col"]]
  bg   <- dots[["bg"]]
  cex  <- dots[["cex"]]
  main <- dots[["main"]]

  # set up empty plot
  graphics::plot(
    NA, NA,
    xlim = data$xlim,
    ylim = data$ylim,
    xlab = data$xlab,
    ylab = data$ylab,
    main = if (!is.null(main)) main else "",
    bty  = "l"
  )

  # draw envelope (if present)
  if (!is.null(data$envelope)) {
    env <- data$envelope

    # shaded region
    if (!is.na(dots[["shade"]])) {
      graphics::polygon(
        x      = c(env$x, rev(env$x)),
        y      = c(env$lower, rev(env$upper)),
        col    = dots[["shade"]],
        border = NA
      )
    }

    # envelope boundary lines
    graphics::lines(env$x, env$lower, lty = "dotted", col = dots[["col.envelope"]])
    graphics::lines(env$x, env$upper, lty = "dotted", col = dots[["col.envelope"]])
  }

  # reference line y = x
  graphics::abline(a = 0, b = 1, lty = "dashed", col = dots[["col.refline"]])

  # points (drawn last, on top)
  graphics::points(
    data$points$x, data$points$y,
    pch = pch, col = col, bg = bg, cex = cex
  )

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .qqnorm_plot_ggplot
# ---------------------------------------------------------------------------- #
#
# Create QQ normal plot using ggplot2.
#
# @param data list of QQ plot data from .qqnorm_data()
# @param dots list of graphical parameters
#
# @return ggplot object
#
# ---------------------------------------------------------------------------- #
.qqnorm_plot_ggplot <- function(data, dots) {

  # extract graphical parameters
  pch  <- dots[["pch"]]
  col  <- dots[["col"]]
  bg   <- dots[["bg"]]
  size <- dots[["size"]]
  main <- dots[["main"]]

  out <- ggplot2::ggplot()

  # envelope (if present)
  if (!is.null(data$envelope)) {
    env <- data$envelope

    # shaded region
    if (!is.na(dots[["shade"]])) {
      out <- out +
        ggplot2::geom_ribbon(
          data    = env,
          mapping = ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
          fill    = dots[["shade"]]
        )
    }

    # envelope boundary lines
    out <- out +
      ggplot2::geom_line(
        data     = env,
        mapping  = ggplot2::aes(x = .data$x, y = .data$lower),
        linetype = "dotted",
        colour   = dots[["col.envelope"]]
      ) +
      ggplot2::geom_line(
        data     = env,
        mapping  = ggplot2::aes(x = .data$x, y = .data$upper),
        linetype = "dotted",
        colour   = dots[["col.envelope"]]
      )
  }

  # reference line y = x
  out <- out +
    ggplot2::geom_abline(
      intercept = 0,
      slope     = 1,
      linetype  = "dashed",
      colour    = dots[["col.refline"]]
    )

  # points
  out <- out +
    ggplot2::geom_point(
      data    = data$points,
      mapping = ggplot2::aes(x = .data$x, y = .data$y),
      shape   = pch,
      colour  = col,
      fill    = bg,
      size    = size
    )

  # scales and coord
  out <- out +
    ggplot2::scale_x_continuous(name = data$xlab) +
    ggplot2::scale_y_continuous(name = data$ylab) +
    ggplot2::coord_cartesian(xlim = data$xlim, ylim = data$ylim) +
    ggplot2::theme_bw()

  # title
  if (!is.null(main)) {
    out <- out + ggplot2::ggtitle(main)
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .set_dots_qqnorm
# ---------------------------------------------------------------------------- #
#
# Process dots arguments for QQ normal plot with sensible defaults.
#
# @param ... additional graphical arguments
#
# @return list of graphical parameters with defaults applied
#
# ---------------------------------------------------------------------------- #
.set_dots_qqnorm <- function(...) {

  dots <- list(...)

  # point styling
  if (is.null(dots[["pch"]]))  dots[["pch"]]  <- 21
  if (is.null(dots[["col"]]))  dots[["col"]]  <- "black"
  if (is.null(dots[["bg"]]))   dots[["bg"]]   <- "black"
  if (is.null(dots[["cex"]]))  dots[["cex"]]  <- 1
  if (is.null(dots[["size"]])) dots[["size"]] <- 2

  # envelope styling
  if (is.null(dots[["shade"]]))        dots[["shade"]]        <- "grey90"
  if (is.null(dots[["col.envelope"]])) dots[["col.envelope"]] <- "grey50"

  # reference line
  if (is.null(dots[["col.refline"]])) dots[["col.refline"]] <- "black"

  # title
  if (is.null(dots[["main"]])) dots[["main"]] <- NULL

  # data return flag
  if (is.null(dots[["as_data"]])) dots[["as_data"]] <- FALSE

  return(dots)
}
