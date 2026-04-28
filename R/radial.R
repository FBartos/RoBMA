# ============================================================================ #
# brma.radial.R
# ============================================================================ #
#
# Radial (Galbraith) plot functions for brma objects.
#
# The radial plot displays effect sizes on a transformed scale where:
# - x-axis: precision (1/sqrt(vi + tau^2))
# - z-axis: standardized effect (yi/sqrt(vi + tau^2))
#
# A line from the origin through any point has slope equal to the
# observed effect size. An arc on the right side maps z-values back
# to the effect size scale.
#
# Only supported for intercept-only models (no moderators, no scale).
#
# ============================================================================ #


### radial plot functions ----
#' @export
radial <- function(x, ...) UseMethod("radial")

#' @export
galbraith <- function(x, ...) UseMethod("galbraith")


#' @title Radial (Galbraith) Plot for brma Object
#'
#' @description \code{radial.brma} creates a radial (Galbraith) plot for a
#' fitted brma object. The plot displays each study as a point with
#' precision on the x-axis and the standardized effect size on the z-axis.
#' A line from the origin through any point has slope equal to the observed
#' effect size, and an arc on the right side maps standardized values back
#' to the effect size scale.
#'
#' @param x a fitted brma object. Must be an intercept-only model
#' (no moderators or scale regression).
#' @param center logical indicating whether to center the plot at the
#' pooled estimate. When \code{TRUE}, the pooled estimate is shifted to
#' zero on the z-axis. Defaults to \code{FALSE}.
#' @param xlim x-axis limits. If not specified, limits are computed from data.
#' @param zlim z-axis limits. If not specified, limits are computed from data.
#' @param xlab title for the x-axis. If not specified, a default label is used.
#' @param zlab title for the z-axis. If not specified, a default label is used.
#' @param atz numeric vector of positions for tick marks on the z-axis (left
#' vertical axis). If not specified, positions are chosen automatically.
#' @param aty numeric vector of effect size values to mark on the y-axis arc
#' (right side). If not specified, values are chosen automatically.
#' @param steps integer specifying the number of tick marks on the y-axis arc
#' when \code{aty} is not specified. Defaults to \code{7}.
#' @param level numeric value between 0 and 100 specifying the confidence level
#' for the z-axis band. Defaults to \code{95}.
#' @param digits integer specifying the number of decimal places for y-axis
#' arc labels. Defaults to \code{2}.
#' @param transf optional transformation function applied to the y-axis arc
#' labels (e.g., \code{transf = exp} when effect sizes are on a log scale).
#' @param targs optional list of additional arguments passed to the
#' \code{transf} function.
#' @param plot_type whether to use a base plot \code{"base"} or ggplot2
#' \code{"ggplot"} for plotting. Defaults to \code{"base"}.
#' @param ... additional graphical arguments to customize the plot appearance:
#' \describe{
#'   \item{pch}{point symbol (default: 21, filled circle)}
#'   \item{col}{point border color (default: "black")}
#'   \item{bg}{point fill/background color (default: "black")}
#'   \item{cex}{point size multiplier for base graphics (default: 1)}
#'   \item{size}{point size for ggplot2 (default: 2)}
#'   \item{back}{background color for the confidence band (default: "grey90").
#'     Set to \code{NA} to suppress.}
#'   \item{arc.res}{integer specifying the number of line segments used to
#'     draw arcs (default: 100)}
#'   \item{main}{character string for plot title (default: no title)}
#'   \item{as_data}{if \code{TRUE}, returns plot data instead of creating
#'     the plot}
#' }
#'
#' @details
#' The radial (Galbraith) plot transforms each study's effect size and
#' standard error into a point in precision-standardized space:
#'
#' \itemize{
#'   \item x-axis: precision = \eqn{1/\sqrt{v_i + \hat{\tau}^2}}
#'   \item z-axis: standardized effect = \eqn{y_i/\sqrt{v_i + \hat{\tau}^2}}
#' }
#'
#' Under the random-effects model, studies consistent with the pooled effect
#' should fall within a horizontal confidence band centered at
#' \eqn{z = \pm z_{\alpha/2}}. The arc on the right side allows reading
#' individual effect sizes by projecting from the origin through a point
#' to the arc.
#'
#' This function requires an intercept-only model; radial plots are not
#' meaningful for meta-regression or location-scale models where the
#' pooled effect varies across studies.
#'
#' @return \code{radial.brma} returns \code{NULL} invisibly if
#' \code{plot_type = "base"} or a ggplot object if \code{plot_type = "ggplot"}.
#'
#' If \code{as_data = TRUE}, returns a list with the computed plot data
#' including: \code{points}, \code{refline}, \code{band}, \code{arc},
#' \code{arc_ticks}, \code{ci_arc}, \code{ci_ticks}, \code{xlim}, and
#' \code{zlim}.
#'
#' @references
#' Galbraith, R. F. (1988). Graphical display of estimates having differing
#' standard errors. \emph{Technometrics}, 30(3), 271-281.
#'
#' Galbraith, R. F. (1994). Some applications of radial plots. \emph{Journal
#' of the American Statistical Association}, 89(428), 1232-1242.
#'
#' @examples \dontrun{
#' # Simple meta-analysis
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#' radial(fit)
#'
#' # Centered at pooled estimate
#' radial(fit, center = TRUE)
#'
#' # Using ggplot2
#' radial(fit, plot_type = "ggplot")
#'
#' # Custom appearance
#' radial(fit, pch = 19, col = "blue")
#'
#' # galbraith() is an alias
#' galbraith(fit)
#' }
#'
#' @seealso [funnel.brma()], [pooled_effect.brma()], [pooled_heterogeneity.brma()]
#' @export
#' @rdname radial
radial.brma <- function(x, center = FALSE, xlim, zlim, xlab, zlab,
                        atz, aty, steps = 7, level = 95, digits = 2,
                        transf, targs, plot_type = "base", ...) {

  # input validation
  BayesTools::check_bool(center, "center")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))

  # check model type - error on models with moderators or scale
  if (.is_mods(x)) {
    stop("Radial plots cannot be drawn for models that contain moderators.", call. = FALSE)
  }
  if (.is_scale(x)) {
    stop("Radial plots cannot be drawn for models with scale regression.", call. = FALSE)
  }

  # set up graphical arguments with defaults
  dots <- .set_dots_radial(...)

  # generate radial plot data
  radial_data <- .radial_data(
    x       = x,
    center  = center,
    xlim    = if (missing(xlim)) NULL else xlim,
    zlim    = if (missing(zlim)) NULL else zlim,
    xlab    = if (missing(xlab)) NULL else xlab,
    zlab    = if (missing(zlab)) NULL else zlab,
    atz     = if (missing(atz)) NULL else atz,
    aty     = if (missing(aty)) NULL else aty,
    steps   = steps,
    level   = level,
    digits  = digits,
    transf  = if (missing(transf)) NULL else transf,
    targs   = if (missing(targs)) NULL else targs,
    dots    = dots
  )

  # allow data return for programmatic access
  if (isTRUE(dots[["as_data"]])) {
    return(radial_data)
  }

  # create plot
  if (plot_type == "ggplot") {
    return(.radial_plot_ggplot(radial_data, dots))
  } else {
    .radial_plot_base(radial_data, dots)
    return(invisible(NULL))
  }
}


#' @export
#' @rdname radial
galbraith.brma <- function(x, ...) {

  radial.brma(x, ...)
}


# ---------------------------------------------------------------------------- #
# .radial_data
# ---------------------------------------------------------------------------- #
#
# Generate data for radial (Galbraith) plot.
#
# Computes precision (x-axis), standardized effects (z-axis),
# confidence band, reference line, and arc data.
#
# Arc coordinates are computed using trigonometry (for ggplot and as_data).
# The base R plotting function recomputes arcs using the aspect ratio
# for correct circular appearance (matching metafor::radial).
#
# @param x       brma object
# @param center  logical; center at pooled estimate
# @param xlim    x-axis limits (NULL for auto)
# @param zlim    z-axis limits (NULL for auto)
# @param xlab    x-axis label (NULL for default)
# @param zlab    z-axis label (NULL for default)
# @param atz     z-axis tick positions (NULL for auto)
# @param aty     arc tick values (NULL for auto)
# @param steps   number of arc ticks when aty is NULL
# @param level   confidence level (0-100)
# @param digits  decimal places for arc labels
# @param transf  transformation for arc labels (NULL for none)
# @param targs   additional args for transf (NULL for none)
# @param dots    list of graphical parameters
#
# @return list with radial plot data components
#
# ---------------------------------------------------------------------------- #
.radial_data <- function(x, center, xlim, zlim, xlab, zlab,
                         atz, aty, steps, level, digits, transf, targs, dots) {

  # get observed effect sizes and variances
  yi <- .outcome_data_yi(x)
  vi <- .outcome_data_vi(x)
  K  <- length(yi)

  # confidence level
  alpha  <- 1 - level / 100
  z_crit <- stats::qnorm(1 - alpha / 2)
  probs  <- c(alpha / 2, 1 - alpha / 2)

  # get pooled effect estimate
  mu_samples <- pooled_effect(x, probs = probs)
  mu_summary <- summary(mu_samples)
  beta       <- mu_summary["mu", "Mean"]
  ci.lb      <- mu_summary["mu", as.character(probs[1])]
  ci.ub      <- mu_summary["mu", as.character(probs[2])]

  # get heterogeneity estimate
  tau_samples <- pooled_heterogeneity(x)
  tau_summary <- summary(tau_samples)
  tau         <- tau_summary["tau", "Mean"]
  tau2        <- tau^2

  # compute precision and standardized values
  wi <- vi + tau2
  xi <- 1 / sqrt(wi)
  zi <- yi / sqrt(wi)

  # save uncentered yi for arc label computation
  yi.c <- yi

  # arc range tracking (on original scale before centering)
  if (is.null(aty)) {
    atyis <- range(yi)
  } else {
    atyis <- range(aty)
    aty.c <- aty
  }

  # center if requested
  if (center) {
    yi         <- yi - beta
    zi         <- yi / sqrt(wi)
    beta_plot  <- 0
    ci.lb_plot <- ci.lb - beta
    ci.ub_plot <- ci.ub - beta
    atyis      <- atyis - beta
    if (!is.null(aty)) aty <- aty - beta
  } else {
    beta_plot  <- beta
    ci.lb_plot <- ci.lb
    ci.ub_plot <- ci.ub
  }

  # ---- x-axis limits (matching metafor: 30% wider) ----
  if (is.null(xlim)) {
    xlim <- c(0, 1.3 * max(xi))
  }
  xaxismax <- xlim[2]

  # arc positions (outside xlim, matching metafor)
  ci.xpos <- xlim[2] + 0.12 * diff(xlim)
  ya.xpos <- xlim[2] + 0.14 * diff(xlim)

  arc_res <- dots[["arc.res"]]

  # ---- axis labels ----
  if (is.null(xlab)) {
    xlab <- expression(x[i] == 1 / sqrt(v[i] + tau^2))
  }

  # zlab: fraction expression for base R mtext, simpler for ggplot
  zlab_auto <- is.null(zlab)
  if (zlab_auto) {
    if (center) {
      zlab <- expression(z[i] == frac(y[i] - hat(mu), sqrt(v[i] + tau^2)))
    } else {
      zlab <- expression(z[i] == frac(y[i], sqrt(v[i] + tau^2)))
    }
  }

  # ---- arc tick values (effect size scale) ----
  if (is.null(aty)) {
    aty   <- pretty(atyis, n = steps)
    aty   <- unique(aty)
    aty.c <- if (center) aty + beta else aty
  }
  # else: aty (centered if needed) and aty.c (original) already set

  # arc tick labels (optionally transformed)
  if (!is.null(transf)) {
    if (!is.null(targs)) {
      label_values <- do.call(transf, c(list(aty.c), targs))
    } else {
      label_values <- transf(aty.c)
    }
    aty_label_text <- formatC(label_values, digits = digits, format = "f")
  } else {
    aty_label_text <- formatC(aty.c, digits = digits, format = "f")
  }

  # ---- z-axis limits (matching metafor) ----
  if (is.null(zlim)) {
    zlim <- c(
      min(-5, 1.1 * min(zi),
          1.1 * ci.lb_plot * ci.xpos,
          1.1 * min(aty) * ya.xpos,
          1.1 * min(yi) * ya.xpos,
          -1.1 * z_crit + xaxismax * beta_plot),
      max(5, 1.1 * max(zi),
          1.1 * ci.ub_plot * ci.xpos,
          1.1 * max(aty) * ya.xpos,
          1.1 * max(yi) * ya.xpos,
          1.1 * z_crit + xaxismax * beta_plot)
    )
  }

  # z-axis tick positions
  if (is.null(atz)) {
    atz <- pretty(zlim)
  }

  # ---- plot data ----
  df_points <- data.frame(x = xi, z = zi)

  # reference line (within axis range only)
  df_refline <- data.frame(
    x = c(0, xaxismax),
    z = c(0, xaxismax * beta_plot)
  )

  # confidence band (parallelogram, within axis range)
  df_band <- data.frame(
    x = c(0, xaxismax, xaxismax, 0),
    z = c(z_crit, z_crit + xaxismax * beta_plot,
          -z_crit + xaxismax * beta_plot, -z_crit)
  )
  df_band_upper <- data.frame(
    x = c(0, xaxismax),
    z = c(z_crit, z_crit + xaxismax * beta_plot)
  )
  df_band_lower <- data.frame(
    x = c(0, xaxismax),
    z = c(-z_crit, -z_crit + xaxismax * beta_plot)
  )

  # ---- arc coordinates (trigonometric, for ggplot and as_data) ----
  # outer arc
  arc_slopes <- seq(min(aty), max(aty), length.out = arc_res)
  arc_theta  <- atan(arc_slopes)
  df_arc <- data.frame(
    x = ya.xpos * cos(arc_theta),
    z = ya.xpos * sin(arc_theta)
  )

  # outer arc tick marks
  tick_theta <- atan(aty)
  tick_len   <- 0.015 * diff(xlim)
  df_arc_ticks <- data.frame(
    x     = ya.xpos * cos(tick_theta),
    z     = ya.xpos * sin(tick_theta),
    x_end = (ya.xpos + tick_len) * cos(tick_theta),
    z_end = (ya.xpos + tick_len) * sin(tick_theta),
    value = aty,
    label = aty_label_text
  )

  # CI arc
  ci_values    <- c(ci.lb_plot, beta_plot, ci.ub_plot)
  ci_arc_slopes <- seq(ci.lb_plot, ci.ub_plot, length.out = ceiling(arc_res / 4))
  ci_arc_theta  <- atan(ci_arc_slopes)
  df_ci_arc <- data.frame(
    x = ci.xpos * cos(ci_arc_theta),
    z = ci.xpos * sin(ci_arc_theta)
  )

  # CI arc tick marks
  ci_tick_theta <- atan(ci_values)
  ci_tick_len   <- 0.007 * diff(xlim)
  df_ci_ticks <- data.frame(
    x     = (ci.xpos - ci_tick_len) * cos(ci_tick_theta),
    z     = (ci.xpos - ci_tick_len) * sin(ci_tick_theta),
    x_end = (ci.xpos + ci_tick_len) * cos(ci_tick_theta),
    z_end = (ci.xpos + ci_tick_len) * sin(ci_tick_theta),
    value = ci_values,
    type  = c("ci.lb", "beta", "ci.ub")
  )

  return(list(
    # plot elements
    points     = df_points,
    refline    = df_refline,
    band       = df_band,
    band_upper = df_band_upper,
    band_lower = df_band_lower,
    # arcs (trigonometric, for ggplot/as_data)
    arc        = df_arc,
    arc_ticks  = df_arc_ticks,
    ci_arc     = df_ci_arc,
    ci_ticks   = df_ci_ticks,
    # limits and labels
    xlim       = xlim,
    zlim       = zlim,
    xlab       = xlab,
    zlab       = zlab,
    zlab_auto  = zlab_auto,
    atz        = atz,
    # raw values for base R arc computation (aspect-ratio-aware)
    aty        = aty,
    aty_labels = aty_label_text,
    ci_values  = ci_values,
    beta_plot  = beta_plot,
    z_crit     = z_crit,
    xaxismax   = xaxismax,
    ci.xpos    = ci.xpos,
    ya.xpos    = ya.xpos,
    arc_res    = arc_res,
    level      = level,
    center     = center
  ))
}


# ---------------------------------------------------------------------------- #
# .radial_plot_base
# ---------------------------------------------------------------------------- #
#
# Create radial plot using base R graphics, matching metafor::radial.
#
# The arc is drawn outside the plot axis range using xpd = TRUE.
# Arc coordinates are computed using the plot's aspect ratio so the
# arc appears circular on screen (matching metafor's approach).
#
# @param data list of radial plot data from .radial_data()
# @param dots list of graphical parameters
#
# @return NULL invisibly
#
# ---------------------------------------------------------------------------- #
.radial_plot_base <- function(data, dots) {

  # extract graphical parameters
  pch  <- dots[["pch"]]
  col  <- dots[["col"]]
  bg   <- dots[["bg"]]
  cex  <- dots[["cex"]]
  back <- dots[["back"]]
  main <- dots[["main"]]
  las  <- dots[["las"]]

  # extract data components
  df_points  <- data$points
  df_refline <- data$refline
  df_band    <- data$band
  df_band_up <- data$band_upper
  df_band_lo <- data$band_lower
  xlim       <- data$xlim
  zlim       <- data$zlim
  xlab       <- data$xlab
  zlab       <- data$zlab
  zlab_auto  <- data$zlab_auto
  atz        <- data$atz
  xaxismax   <- data$xaxismax

  # raw arc parameters (recomputed with aspect ratio below)
  aty        <- data$aty
  aty_labels <- data$aty_labels
  ci_values  <- data$ci_values
  ci.xpos    <- data$ci.xpos
  ya.xpos    <- data$ya.xpos
  arc_res    <- data$arc_res

  # ---- save and restore graphical parameters ----
  par_mar <- graphics::par("mar")
  par_pty <- graphics::par("pty")
  par_xpd <- graphics::par("xpd")

  # wider margins: +4 left (for fraction zlab), +6 right (for arc labels)
  par_mar_adj <- par_mar + c(0, 4, 0, 6)
  par_mar_adj[par_mar_adj < 1] <- 1
  graphics::par(mar = par_mar_adj, pty = "s")

  on.exit(graphics::par(mar = par_mar, pty = par_pty, xpd = par_xpd), add = TRUE)

  # ---- empty plot: no box, no auto axes, exact limits ----
  graphics::plot(
    NA, NA,
    xlim = xlim,
    ylim = zlim,
    bty  = "n",
    xaxt = "n",
    yaxt = "n",
    xlab = xlab,
    ylab = "",
    xaxs = "i",
    yaxs = "i",
    las  = las
  )

  # title
  if (!is.null(main)) {
    graphics::title(main = main)
  }

  # ---- confidence band (parallelogram) ----
  if (!is.na(back)) {
    graphics::polygon(
      x      = df_band$x,
      y      = df_band$z,
      col    = back,
      border = NA
    )
  }

  # reference line
  graphics::segments(
    df_refline$x[1], df_refline$z[1],
    df_refline$x[2], df_refline$z[2],
    lty = "solid"
  )

  # band edges
  graphics::segments(
    df_band_up$x[1], df_band_up$z[1],
    df_band_up$x[2], df_band_up$z[2],
    lty = "dotted"
  )
  graphics::segments(
    df_band_lo$x[1], df_band_lo$z[1],
    df_band_lo$x[2], df_band_lo$z[2],
    lty = "dotted"
  )

  # ---- axes (manual) ----
  graphics::axis(side = 1, las = las)

  graphics::axis(side = 2, at = atz, labels = atz, las = las)

  # z-axis label via mtext (fraction expression renders properly)
  if (zlab_auto) {
    graphics::mtext(
      zlab, side = 2, line = par_mar_adj[2] - 1,
      at = 0, adj = 0, las = 1
    )
  } else {
    graphics::mtext(
      zlab, side = 2, line = par_mar_adj[2] - 4,
      at = 0
    )
  }

  # ---- arcs (outside plot region, using aspect ratio) ----
  graphics::par(xpd = TRUE)

  par_usr <- graphics::par("usr")
  asp_rat <- (par_usr[4] - par_usr[3]) / (par_usr[2] - par_usr[1])

  # helper: compute (x, z) on circle of radius `len` at slope `s`
  # accounting for aspect ratio (matching metafor)
  .arc_x <- function(len, s) sqrt(len^2 / (1 + (s / asp_rat)^2))
  .arc_z <- function(len, s) .arc_x(len, s) * s

  # ---- outer arc line ----
  arc_slopes <- seq(min(aty), max(aty), length.out = arc_res)
  arc_xi     <- .arc_x(ya.xpos, arc_slopes)
  arc_zi     <- .arc_z(ya.xpos, arc_slopes)

  valid <- arc_zi > zlim[1] & arc_zi < zlim[2]
  if (any(valid)) {
    graphics::lines(arc_xi[valid], arc_zi[valid])
  }

  # ---- outer arc tick marks ----
  tick_len  <- 0.015 * diff(xlim)
  tick_xi_l <- .arc_x(ya.xpos, aty)
  tick_zi_l <- .arc_z(ya.xpos, aty)
  tick_xi_u <- .arc_x(ya.xpos + tick_len, aty)
  tick_zi_u <- .arc_z(ya.xpos + tick_len, aty)

  valid <- tick_zi_l > zlim[1] & tick_zi_u > zlim[1] &
           tick_zi_l < zlim[2] & tick_zi_u < zlim[2]
  if (any(valid)) {
    graphics::segments(
      tick_xi_l[valid], tick_zi_l[valid],
      tick_xi_u[valid], tick_zi_u[valid]
    )
  }

  # ---- outer arc labels ----
  label_len <- ya.xpos + 0.02 * diff(xlim)
  label_xi  <- .arc_x(label_len, aty)
  label_zi  <- .arc_z(label_len, aty)

  valid <- label_zi > zlim[1] & label_zi < zlim[2]
  if (any(valid)) {
    graphics::text(
      label_xi[valid], label_zi[valid],
      aty_labels[valid],
      pos = 4, cex = 0.85, offset = 0.25
    )
  }

  # ---- CI arc line ----
  ci_arc_slopes <- seq(ci_values[1], ci_values[3],
                       length.out = ceiling(arc_res / 4))
  ci_arc_xi     <- .arc_x(ci.xpos, ci_arc_slopes)
  ci_arc_zi     <- .arc_z(ci.xpos, ci_arc_slopes)

  valid <- ci_arc_zi > zlim[1] & ci_arc_zi < zlim[2]
  if (any(valid)) {
    graphics::lines(ci_arc_xi[valid], ci_arc_zi[valid])
  }

  # ---- CI arc tick marks ----
  ci_tick_len  <- 0.007 * diff(xlim)
  ci_tick_xi_l <- .arc_x(ci.xpos - ci_tick_len, ci_values)
  ci_tick_zi_l <- .arc_z(ci.xpos - ci_tick_len, ci_values)
  ci_tick_xi_u <- .arc_x(ci.xpos + ci_tick_len, ci_values)
  ci_tick_zi_u <- .arc_z(ci.xpos + ci_tick_len, ci_values)

  valid <- ci_tick_zi_l > zlim[1] & ci_tick_zi_u > zlim[1] &
           ci_tick_zi_l < zlim[2] & ci_tick_zi_u < zlim[2]
  if (any(valid)) {
    graphics::segments(
      ci_tick_xi_l[valid], ci_tick_zi_l[valid],
      ci_tick_xi_u[valid], ci_tick_zi_u[valid]
    )
  }

  # ---- restore xpd and draw points last (on top) ----
  graphics::par(xpd = par_xpd)

  graphics::points(
    df_points$x, df_points$z,
    pch = pch, col = col, bg = bg, cex = cex
  )

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .radial_plot_ggplot
# ---------------------------------------------------------------------------- #
#
# Create radial plot using ggplot2.
#
# Uses pre-computed trigonometric arc coordinates from .radial_data().
# Arcs are drawn outside the panel via coord_cartesian(clip = "off").
#
# @param data list of radial plot data from .radial_data()
# @param dots list of graphical parameters
#
# @return ggplot object
#
# ---------------------------------------------------------------------------- #
.radial_plot_ggplot <- function(data, dots) {

  # extract graphical parameters
  pch  <- dots[["pch"]]
  col  <- dots[["col"]]
  bg   <- dots[["bg"]]
  size <- dots[["size"]]
  back <- dots[["back"]]
  main <- dots[["main"]]

  # extract data components
  df_points    <- data$points
  df_refline   <- data$refline
  df_band      <- data$band
  df_band_up   <- data$band_upper
  df_band_lo   <- data$band_lower
  xlim         <- data$xlim
  zlim         <- data$zlim
  xlab         <- data$xlab
  zlab         <- data$zlab
  atz          <- data$atz

  # raw arc parameters (recomputed with aspect ratio, matching base R)
  aty        <- data$aty
  aty_labels <- data$aty_labels
  ci_values  <- data$ci_values
  ci.xpos    <- data$ci.xpos
  ya.xpos    <- data$ya.xpos
  arc_res    <- data$arc_res

  # ---- compute aspect-ratio-corrected arc coordinates ----
  # matching the base R approach: square panel with data aspect ratio
  asp_rat <- diff(zlim) / diff(xlim)

  # helper: compute (x, z) on circle of radius `len` at slope `s`
  # accounting for aspect ratio (matching metafor / base R)
  .arc_x <- function(len, s) sqrt(len^2 / (1 + (s / asp_rat)^2))
  .arc_z <- function(len, s) .arc_x(len, s) * s

  # ---- outer arc line ----
  arc_slopes <- seq(min(aty), max(aty), length.out = arc_res)
  df_arc <- data.frame(
    x = .arc_x(ya.xpos, arc_slopes),
    z = .arc_z(ya.xpos, arc_slopes)
  )
  valid <- df_arc$z > zlim[1] & df_arc$z < zlim[2]
  df_arc <- df_arc[valid, ]

  # ---- outer arc tick marks ----
  tick_len <- 0.015 * diff(xlim)
  df_arc_ticks <- data.frame(
    x     = .arc_x(ya.xpos, aty),
    z     = .arc_z(ya.xpos, aty),
    x_end = .arc_x(ya.xpos + tick_len, aty),
    z_end = .arc_z(ya.xpos + tick_len, aty),
    value = aty,
    label = aty_labels
  )
  valid <- df_arc_ticks$z > zlim[1] & df_arc_ticks$z_end > zlim[1] &
           df_arc_ticks$z < zlim[2] & df_arc_ticks$z_end < zlim[2]
  df_arc_ticks <- df_arc_ticks[valid, ]

  # ---- outer arc labels ----
  label_len <- ya.xpos + 0.02 * diff(xlim)
  df_arc_labels <- data.frame(
    x     = .arc_x(label_len, aty),
    z     = .arc_z(label_len, aty),
    label = aty_labels
  )
  valid <- df_arc_labels$z > zlim[1] & df_arc_labels$z < zlim[2]
  df_arc_labels <- df_arc_labels[valid, ]

  # ---- CI arc line ----
  ci_arc_slopes <- seq(ci_values[1], ci_values[3],
                       length.out = ceiling(arc_res / 4))
  df_ci_arc <- data.frame(
    x = .arc_x(ci.xpos, ci_arc_slopes),
    z = .arc_z(ci.xpos, ci_arc_slopes)
  )
  valid <- df_ci_arc$z > zlim[1] & df_ci_arc$z < zlim[2]
  df_ci_arc <- df_ci_arc[valid, ]

  # ---- CI arc tick marks ----
  ci_tick_len <- 0.007 * diff(xlim)
  df_ci_ticks <- data.frame(
    x     = .arc_x(ci.xpos - ci_tick_len, ci_values),
    z     = .arc_z(ci.xpos - ci_tick_len, ci_values),
    x_end = .arc_x(ci.xpos + ci_tick_len, ci_values),
    z_end = .arc_z(ci.xpos + ci_tick_len, ci_values),
    value = ci_values,
    type  = c("ci.lb", "beta", "ci.ub")
  )
  valid <- df_ci_ticks$z > zlim[1] & df_ci_ticks$z_end > zlim[1] &
           df_ci_ticks$z < zlim[2] & df_ci_ticks$z_end < zlim[2]
  df_ci_ticks <- df_ci_ticks[valid, ]

  # ---- build plot ----
  out <- ggplot2::ggplot()

  # confidence band (parallelogram)
  if (!is.na(back)) {
    out <- out +
      ggplot2::geom_polygon(
        data    = df_band,
        mapping = ggplot2::aes(x = .data$x, y = .data$z),
        fill    = back
      )
  }

  # confidence band edges
  out <- out +
    ggplot2::geom_line(
      data     = df_band_up,
      mapping  = ggplot2::aes(x = .data$x, y = .data$z),
      linetype = "dotted",
      colour   = "grey50"
    ) +
    ggplot2::geom_line(
      data     = df_band_lo,
      mapping  = ggplot2::aes(x = .data$x, y = .data$z),
      linetype = "dotted",
      colour   = "grey50"
    )

  # reference line
  out <- out +
    ggplot2::geom_line(
      data    = df_refline,
      mapping = ggplot2::aes(x = .data$x, y = .data$z),
      colour  = "black"
    )

  # points
  out <- out +
    ggplot2::geom_point(
      data    = df_points,
      mapping = ggplot2::aes(x = .data$x, y = .data$z),
      shape   = pch,
      colour  = col,
      fill    = bg,
      size    = size
    )

  # outer arc (drawn outside panel via clip = "off")
  out <- out +
    ggplot2::geom_path(
      data    = df_arc,
      mapping = ggplot2::aes(x = .data$x, y = .data$z),
      colour  = "black"
    )

  # arc tick marks
  if (nrow(df_arc_ticks) > 0) {
    out <- out +
      ggplot2::geom_segment(
        data    = df_arc_ticks,
        mapping = ggplot2::aes(
          x    = .data$x,
          y    = .data$z,
          xend = .data$x_end,
          yend = .data$z_end
        ),
        colour = "black"
      )
  }

  # arc tick labels
  if (nrow(df_arc_labels) > 0) {
    out <- out +
      ggplot2::geom_text(
        data    = df_arc_labels,
        mapping = ggplot2::aes(x = .data$x, y = .data$z, label = .data$label),
        hjust   = 0,
        size    = 3
      )
  }

  # CI arc
  out <- out +
    ggplot2::geom_path(
      data      = df_ci_arc,
      mapping   = ggplot2::aes(x = .data$x, y = .data$z),
      colour    = "black",
      linewidth = 0.8
    )

  # CI arc ticks
  if (nrow(df_ci_ticks) > 0) {
    out <- out +
      ggplot2::geom_segment(
        data      = df_ci_ticks,
        mapping   = ggplot2::aes(
          x    = .data$x,
          y    = .data$z,
          xend = .data$x_end,
          yend = .data$z_end
        ),
        colour    = "black",
        linewidth = ifelse(df_ci_ticks$type == "beta", 0.8, 0.5)
      )
  }

  # scales
  out <- out +
    ggplot2::scale_x_continuous(name = xlab) +
    ggplot2::scale_y_continuous(breaks = atz, name = zlab)

  # coord: square panel matching base R pty="s"; arcs outside via clip = "off"
  out <- out +
    ggplot2::coord_cartesian(
      xlim = xlim,
      ylim = zlim,
      clip = "off"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio  = 1,
      panel.border  = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(),
      axis.line.y.left   = ggplot2::element_line(),
      plot.margin   = ggplot2::margin(t = 10, r = 60, b = 10, l = 10)
    )

  # title
  if (!is.null(main)) {
    out <- out + ggplot2::ggtitle(main)
  }

  return(out)
}


# ---------------------------------------------------------------------------- #
# .set_dots_radial
# ---------------------------------------------------------------------------- #
#
# Process dots arguments for radial plot with sensible defaults.
#
# @param ... additional graphical arguments
#
# @return list of graphical parameters with defaults applied
#
# ---------------------------------------------------------------------------- #
.set_dots_radial <- function(...) {

  dots <- list(...)

  # point styling
  if (is.null(dots[["pch"]]))  dots[["pch"]]  <- 21
  if (is.null(dots[["col"]]))  dots[["col"]]  <- "black"
  if (is.null(dots[["bg"]]))   dots[["bg"]]   <- "black"
  if (is.null(dots[["cex"]]))  dots[["cex"]]  <- 1
  if (is.null(dots[["size"]])) dots[["size"]]  <- 2

  # confidence band
  if (is.null(dots[["back"]])) dots[["back"]] <- "grey90"

  # arc resolution
  if (is.null(dots[["arc.res"]])) dots[["arc.res"]] <- 100

  # title
  if (is.null(dots[["main"]])) dots[["main"]] <- NULL
  if (is.null(dots[["las"]]))  dots[["las"]]  <- 1

  # data return flag
  if (is.null(dots[["as_data"]])) dots[["as_data"]] <- FALSE

  return(dots)
}
