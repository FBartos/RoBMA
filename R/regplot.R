# ============================================================================ #
# brma.regplot.R
# ============================================================================ #
#
# This file contains the regression plot (bubble plot) functions for brma
# objects with moderators.
#
# Bubble plots display observed effect sizes against a moderator variable,
# with point sizes proportional to study precision (inverse variance).
# For continuous moderators, a regression line with CI/PI bands is shown.
# For categorical moderators, grouped points with optional jitter are shown.
#
# ============================================================================ #


### regplot functions ----
#' @export
regplot <- function(x, ...) UseMethod("regplot")


#' @title Regression Plot (Bubble Plot) for brma Object
#'
#' @description \code{regplot.brma} creates a regression plot (also known as
#' bubble plot) for a fitted brma object with moderators. The plot displays
#' observed effect sizes against a moderator variable, with point sizes
#' proportional to study precision.
#'
#' @param x a fitted brma object with moderators
#' @param mod index or name of the moderator variable to plot on the x-axis.
#' If not specified and only one moderator is present, that moderator is used.
#' If multiple moderators are present, this argument is required.
#' @param pred logical; whether to show the prediction line. Defaults to \code{TRUE}.
#' @param ci logical; whether to show confidence interval bands. Defaults to \code{TRUE}.
#' @param pi logical; whether to show prediction interval bands. Defaults to \code{FALSE}.
#' @param level numeric; confidence/prediction interval level in percent.
#' Defaults to \code{95}.
#' @param at numeric vector; for continuous moderators, values at which to
#' evaluate the prediction. If not specified, uses a sequence across the
#' observed range.
#' @param digits integer; number of decimal places for labels. Defaults to \code{2}.
#' @param transf function; transformation to apply to the y-axis (effect sizes).
#' Defaults to \code{NULL} (no transformation).
#' @param atransf function; transformation to apply to axis labels only.
#' Defaults to \code{NULL} (same as \code{transf}).
#' @param targs list; additional arguments for transformation functions.
#' @param refline numeric; position of horizontal reference line.
#' Defaults to \code{NULL} (no reference line).
#' @param psize numeric vector or \code{NULL}; point sizes for each study.
#' If \code{NULL} (default), sizes are computed based on inverse sampling
#' variance.
#' @param plim numeric vector of length 2; range for point sizes.
#' Defaults to \code{c(0.5, 3)}.
#' @param by character; name of a moderator variable to use for separate
#' lines/colors. Defaults to \code{NULL}.
#' @param legend logical; whether to show legend when \code{by} is specified.
#' Defaults to \code{TRUE}.
#' @param xlab character; x-axis label. Defaults to the moderator name.
#' @param ylab character; y-axis label. Defaults to "Observed Effect Size".
#' @param xlim numeric vector of length 2; x-axis limits. Defaults to data range.
#' @param ylim numeric vector of length 2; y-axis limits. Defaults to data range.
#' @param plot_type character; whether to use base R graphics (\code{"base"})
#' or ggplot2 (\code{"ggplot"}). Defaults to \code{"base"}.
#' @param as_data logical; if \code{TRUE}, returns plot data instead of
#' creating plot. Defaults to \code{FALSE}.
#' @param ... additional graphical arguments:
#' \describe{
#'   \item{main}{character string for plot title}
#'   \item{pch}{point symbol (default: 21)}
#'   \item{col}{point border color (default: "black")}
#'   \item{bg}{point fill color (default: "gray")}
#'   \item{lcol}{line color (default: "black")}
#'   \item{lwd}{line width (default: 2)}
#'   \item{shade}{CI/PI band shading (default: TRUE)}
#'   \item{col.ci}{CI band color (default: "gray70")}
#'   \item{col.pi}{PI band color (default: "gray85")}
#'   \item{alpha.ci}{CI band transparency (default: 0.4)}
#'   \item{alpha.pi}{PI band transparency (default: 0.2)}
#'   \item{jitter}{jitter amount for categorical moderators (default: 0.2)}
#' }
#'
#' @details
#' The regression plot (bubble plot) is a standard visualization for
#' meta-regression results. It displays:
#' \itemize{
#'   \item Observed effect sizes (y-axis) against moderator values (x-axis)
#'   \item Point sizes proportional to study precision (inverse variance)
#'   \item Prediction line showing the estimated regression relationship
#'   \item Confidence bands showing uncertainty in the mean prediction
#'   \item Optional prediction bands showing expected range of true effects
#' }
#'
#' For continuous moderators, predictions are computed across the observed
#' range of the moderator. For categorical moderators (factors), predictions
#' are computed at each factor level with optional jittering of points.
#'
#' The \code{by} argument allows displaying separate regression lines for
#' different levels of a second moderator, useful for visualizing interactions.
#'
#' @return \code{regplot.brma} returns \code{NULL} invisibly if
#' \code{plot_type = "base"} or a ggplot object if \code{plot_type = "ggplot"}.
#' If \code{as_data = TRUE}, returns a list with plot data components.
#'
#' @examples \dontrun{
#' # fit a meta-regression model
#' fit <- brma(yi ~ ablat, sei = sei, data = dat.bcg)
#'
#' # basic bubble plot
#' regplot(fit)
#'
#' # with prediction intervals
#' regplot(fit, pi = TRUE)
#'
#' # using ggplot2
#' regplot(fit, plot_type = "ggplot")
#'
#' # customize appearance
#' regplot(fit, col = "blue", lcol = "darkblue", pch = 19)
#'
#' # categorical moderator
#' fit_cat <- brma(yi ~ alloc, sei = sei, data = dat.bcg)
#' regplot(fit_cat)
#' }
#'
#' @seealso [funnel.brma()], [predict.brma()]
#' @export
#' @rdname regplot
regplot.brma <- function(x, mod = NULL, pred = TRUE, ci = TRUE, pi = FALSE,
                         level = 95, at = NULL, digits = 2,
                         transf = NULL, atransf = NULL, targs = NULL,
                         refline = NULL, psize = NULL, plim = c(0.5, 3),
                         by = NULL, legend = TRUE,
                         xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
                         plot_type = "base", as_data = FALSE, ...) {

  # input validation
  BayesTools::check_bool(pred, "pred")
  BayesTools::check_bool(ci, "ci")
  BayesTools::check_bool(pi, "pi")
  BayesTools::check_real(level, "level", lower = 50, upper = 99.9)
  BayesTools::check_int(digits, "digits", lower = 0)
  BayesTools::check_bool(legend, "legend")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(as_data, "as_data")

  # check that model has moderators
  if (!.is_mods(x)) {
    stop("regplot requires a model with moderators. ",
         "Use funnel() for intercept-only models.", call. = FALSE)
  }

  # set up graphical arguments with defaults
  dots <- .set_dots_regplot(...)

  # identify the moderator to plot
  mod_info <- .regplot_get_moderator(x, mod)
  mod_name <- mod_info$name
  mod_type <- mod_info$type
  mod_data <- mod_info$data

  # identify grouping variable (if specified)
  by_info <- NULL
  if (!is.null(by)) {
    by_info <- .regplot_get_moderator(x, by)
    if (by_info$name == mod_name) {
      stop("'by' variable must be different from 'mod' variable.", call. = FALSE)
    }
  }

  # generate plot data
  regplot_data <- .regplot_data(
    x        = x,
    mod_name = mod_name,
    mod_type = mod_type,
    mod_data = mod_data,
    by_info  = by_info,
    pred     = pred,
    ci       = ci,
    pi       = pi,
    level    = level,
    at       = at,
    psize    = psize,
    plim     = plim,
    transf   = transf,
    xlim     = xlim,
    ylim     = ylim,
    xlab     = xlab,
    ylab     = ylab,
    refline  = refline,
    dots     = dots
  )

  # allow data return for programmatic access
  if (isTRUE(as_data)) {
    return(regplot_data)
  }

  # create plot
  if (plot_type == "ggplot") {
    return(.regplot_plot_ggplot(regplot_data, dots, legend = legend))
  } else {
    .regplot_plot_base(regplot_data, dots, legend = legend)
    return(invisible(NULL))
  }
}


# ---------------------------------------------------------------------------- #
# .regplot_get_moderator
# ---------------------------------------------------------------------------- #
#
# Identify and validate the moderator variable for plotting.
#
# @param x     brma object
# @param mod   index or name of moderator (NULL for auto-detect)
#
# @return list with name, type, and data for the moderator
#
# ---------------------------------------------------------------------------- #
.regplot_get_moderator <- function(x, mod) {

  mods_data <- x[["data"]][["mods"]]

  if (is.null(mods_data)) {
    stop("No moderator data found in the model.", call. = FALSE)
  }

  # get list of moderator variable names (excluding formula attribute)
  mod_names <- names(mods_data)

  if (is.null(mod)) {
    # auto-detect: use single moderator or error if multiple
    if (length(mod_names) == 1) {
      mod_name <- mod_names[1]
    } else {
      stop("Multiple moderators present. Please specify 'mod' argument. ",
           "Available moderators: ", paste(mod_names, collapse = ", "),
           call. = FALSE)
    }
  } else if (is.numeric(mod)) {
    # index-based selection
    if (mod < 1 || mod > length(mod_names)) {
      stop("Invalid moderator index. Must be between 1 and ",
           length(mod_names), ".", call. = FALSE)
    }
    mod_name <- mod_names[mod]
  } else if (is.character(mod)) {
    # name-based selection
    if (!mod %in% mod_names) {
      stop("Moderator '", mod, "' not found. ",
           "Available moderators: ", paste(mod_names, collapse = ", "),
           call. = FALSE)
    }
    mod_name <- mod
  } else {
    stop("'mod' must be NULL, a character name, or a numeric index.",
         call. = FALSE)
  }

  # extract moderator data and determine type
  mod_values <- mods_data[[mod_name]]

  if (is.factor(mod_values) || is.character(mod_values)) {
    mod_type <- "categorical"
    if (is.character(mod_values)) {
      mod_values <- factor(mod_values)
    }
  } else {
    mod_type <- "continuous"
  }

  return(list(
    name = mod_name,
    type = mod_type,
    data = mod_values
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_data
# ---------------------------------------------------------------------------- #
#
# Generate data for regression plot.
#
# @param x        brma object
# @param mod_name name of moderator variable
# @param mod_type "continuous" or "categorical"
# @param mod_data moderator data values
# @param by_info  list with grouping variable info (or NULL)
# @param pred     logical; show prediction line
# @param ci       logical; show CI bands
# @param pi       logical; show PI bands
# @param level    confidence level (0-100)
# @param at       values at which to evaluate predictions
# @param psize    point sizes (or NULL for auto)
# @param plim     point size range
# @param transf   transformation function
# @param xlim     x-axis limits
# @param ylim     y-axis limits
# @param xlab     x-axis label
# @param ylab     y-axis label
# @param refline  reference line position
# @param dots     graphical parameters
#
# @return list with plot data components
#
# ---------------------------------------------------------------------------- #
.regplot_data <- function(x, mod_name, mod_type, mod_data, by_info,
                          pred, ci, pi, level, at, psize, plim, transf,
                          xlim, ylim, xlab, ylab, refline, dots) {

  # get observed data
  yi  <- .outcome_data_yi(x)
  sei <- .outcome_data_sei(x)
  vi  <- sei^2
  K   <- length(yi)

  # convert level to probability for quantiles
  alpha <- (100 - level) / 100
  probs <- c(alpha / 2, 1 - alpha / 2)

  # compute point sizes based on inverse variance if not provided
  if (is.null(psize)) {
    weights   <- 1 / vi
    weights_n <- (weights - min(weights)) / (max(weights) - min(weights) + 1e-10)
    psize     <- plim[1] + weights_n * (plim[2] - plim[1])
  } else if (length(psize) == 1) {
    psize <- rep(psize, K)
  }

  # apply transformation to observed effect sizes
  yi_plot <- yi
  if (!is.null(transf)) {
    yi_plot <- transf(yi)
  }

  # determine x-axis values for prediction
  if (mod_type == "continuous") {
    if (is.null(at)) {
      x_range  <- range(mod_data)
      at_pred  <- seq(x_range[1], x_range[2], length.out = 101)
    } else {
      at_pred <- at
    }
  } else {
    # categorical: use factor levels
    at_pred <- levels(mod_data)
  }

  # create newdata for predictions
  # need to set other moderators to their means/reference levels
  mods_data    <- x[["data"]][["mods"]]
  n_pred       <- length(at_pred)
  newdata_list <- list()

  for (nm in names(mods_data)) {
    if (nm == mod_name) {
      if (mod_type == "continuous") {
        newdata_list[[nm]] <- at_pred
      } else {
        newdata_list[[nm]] <- factor(at_pred, levels = levels(mod_data))
      }
    } else {
      # set other predictors to mean/reference
      other_vals <- mods_data[[nm]]
      if (is.factor(other_vals) || is.character(other_vals)) {
        # use first level (reference)
        newdata_list[[nm]] <- factor(rep(levels(factor(other_vals))[1], n_pred),
                                     levels = levels(factor(other_vals)))
      } else {
        # use mean
        newdata_list[[nm]] <- rep(mean(other_vals), n_pred)
      }
    }
  }
  newdata <- as.data.frame(newdata_list)

  # get predictions for CI (fixed effects only)
  pred_samples <- predict.brma(
    object  = x,
    newdata = newdata,
    type    = "terms",
    quiet   = TRUE
  )

  # compute summary statistics
  pred_mean  <- colMeans(pred_samples)
  pred_lower <- apply(pred_samples, 2, quantile, probs = probs[1])
  pred_upper <- apply(pred_samples, 2, quantile, probs = probs[2])

  # get predictions for PI (includes heterogeneity)
  if (pi) {
    pi_samples <- predict.brma(
      object  = x,
      newdata = newdata,
      type    = "estimate",
      quiet   = TRUE
    )
    pi_lower <- apply(pi_samples, 2, quantile, probs = probs[1])
    pi_upper <- apply(pi_samples, 2, quantile, probs = probs[2])
  } else {
    pi_lower <- NULL
    pi_upper <- NULL
  }

  # apply transformation to predictions
  if (!is.null(transf)) {
    pred_mean  <- transf(pred_mean)
    pred_lower <- transf(pred_lower)
    pred_upper <- transf(pred_upper)
    if (pi) {
      pi_lower <- transf(pi_lower)
      pi_upper <- transf(pi_upper)
    }
  }

  # determine plot limits
  if (is.null(xlim)) {
    if (mod_type == "continuous") {
      xlim <- range(pretty(range(mod_data)))
    } else {
      xlim <- c(0.5, length(levels(mod_data)) + 0.5)
    }
  }

  if (is.null(ylim)) {
    all_y <- c(yi_plot, pred_mean)
    if (ci) all_y <- c(all_y, pred_lower, pred_upper)
    if (pi) all_y <- c(all_y, pi_lower, pi_upper)
    ylim <- range(pretty(range(all_y, na.rm = TRUE)))
  }

  # set axis labels
  if (is.null(xlab)) xlab <- mod_name
  if (is.null(ylab)) ylab <- "Observed Effect Size"

  # create output data structures
  df_points <- data.frame(
    x    = if (mod_type == "continuous") mod_data else as.numeric(mod_data),
    y    = yi_plot,
    size = psize,
    stringsAsFactors = FALSE
  )

  # for categorical, add factor info
  if (mod_type == "categorical") {
    df_points$level <- mod_data
  }

  # prediction line data
  df_pred <- NULL
  if (pred) {
    df_pred <- data.frame(
      x = if (mod_type == "continuous") at_pred else seq_along(at_pred),
      y = pred_mean
    )
  }

  # CI band data
  df_ci <- NULL
  if (ci) {
    if (mod_type == "continuous") {
      df_ci <- data.frame(
        x     = c(at_pred, rev(at_pred)),
        y     = c(pred_lower, rev(pred_upper)),
        lower = pred_lower,
        upper = pred_upper,
        xpred = at_pred
      )
    } else {
      # for categorical, store as points with error bars
      df_ci <- data.frame(
        x     = seq_along(at_pred),
        y     = pred_mean,
        lower = pred_lower,
        upper = pred_upper,
        level = at_pred
      )
    }
  }

  # PI band data
  df_pi <- NULL
  if (pi) {
    if (mod_type == "continuous") {
      df_pi <- data.frame(
        x     = c(at_pred, rev(at_pred)),
        y     = c(pi_lower, rev(pi_upper)),
        lower = pi_lower,
        upper = pi_upper,
        xpred = at_pred
      )
    } else {
      df_pi <- data.frame(
        x     = seq_along(at_pred),
        y     = pred_mean,
        lower = pi_lower,
        upper = pi_upper,
        level = at_pred
      )
    }
  }

  # reference line data
  df_refline <- NULL
  if (!is.null(refline)) {
    df_refline <- data.frame(
      y = refline
    )
  }

  return(list(
    points   = df_points,
    pred     = df_pred,
    ci       = df_ci,
    pi       = df_pi,
    refline  = df_refline,
    xlim     = xlim,
    ylim     = ylim,
    xlab     = xlab,
    ylab     = ylab,
    mod_type = mod_type,
    mod_name = mod_name,
    levels   = if (mod_type == "categorical") levels(mod_data) else NULL
  ))
}


# ---------------------------------------------------------------------------- #
# .regplot_plot_base
# ---------------------------------------------------------------------------- #
#
# Create regression plot using base R graphics.
#
# @param data   list of plot data from .regplot_data
# @param dots   list of graphical parameters
# @param legend logical; show legend
#
# @return NULL invisibly
#
# ---------------------------------------------------------------------------- #
.regplot_plot_base <- function(data, dots, legend = TRUE) {

  # extract graphical parameters
  pch       <- dots[["pch"]]
  col       <- dots[["col"]]
  bg        <- dots[["bg"]]
  lcol      <- dots[["lcol"]]
  lwd       <- dots[["lwd"]]
  col_ci    <- dots[["col.ci"]]
  col_pi    <- dots[["col.pi"]]
  alpha_ci  <- dots[["alpha.ci"]]
  alpha_pi  <- dots[["alpha.pi"]]
  shade     <- dots[["shade"]]
  main      <- dots[["main"]]
  jitter_am <- dots[["jitter"]]

  # extract data components
  df_points  <- data$points
  df_pred    <- data$pred
  df_ci      <- data$ci
  df_pi      <- data$pi
  df_refline <- data$refline
  mod_type   <- data$mod_type

  # set up the plot area
  if (mod_type == "categorical") {
    # categorical: custom x-axis
    graphics::plot(
      NA, NA,
      xlim = data$xlim,
      ylim = data$ylim,
      xlab = data$xlab,
      ylab = data$ylab,
      main = main,
      type = "n",
      xaxt = "n"
    )
    graphics::axis(1, at = seq_along(data$levels), labels = data$levels)
  } else {
    graphics::plot(
      NA, NA,
      xlim = data$xlim,
      ylim = data$ylim,
      xlab = data$xlab,
      ylab = data$ylab,
      main = main,
      type = "n"
    )
  }

  # draw reference line (if specified)
  if (!is.null(df_refline)) {
    graphics::abline(h = df_refline$y, lty = "dashed", col = "gray50")
  }

  # draw PI band (if present and shading enabled)
  if (!is.null(df_pi) && shade) {
    if (mod_type == "continuous") {
      col_pi_alpha <- grDevices::adjustcolor(col_pi, alpha.f = alpha_pi)
      graphics::polygon(df_pi$x, df_pi$y, col = col_pi_alpha, border = NA)
    } else {
      # categorical: draw error bars for PI
      for (i in seq_len(nrow(df_pi))) {
        graphics::segments(
          x0  = df_pi$x[i],
          y0  = df_pi$lower[i],
          x1  = df_pi$x[i],
          y1  = df_pi$upper[i],
          col = grDevices::adjustcolor(col_pi, alpha.f = alpha_pi + 0.3),
          lwd = 6
        )
      }
    }
  }

  # draw CI band (if present and shading enabled)
  if (!is.null(df_ci) && shade) {
    if (mod_type == "continuous") {
      col_ci_alpha <- grDevices::adjustcolor(col_ci, alpha.f = alpha_ci)
      graphics::polygon(df_ci$x, df_ci$y, col = col_ci_alpha, border = NA)
    } else {
      # categorical: draw error bars for CI
      for (i in seq_len(nrow(df_ci))) {
        graphics::segments(
          x0  = df_ci$x[i],
          y0  = df_ci$lower[i],
          x1  = df_ci$x[i],
          y1  = df_ci$upper[i],
          col = grDevices::adjustcolor(col_ci, alpha.f = alpha_ci + 0.3),
          lwd = 3
        )
      }
    }
  }

  # draw prediction line
  if (!is.null(df_pred)) {
    if (mod_type == "continuous") {
      graphics::lines(df_pred$x, df_pred$y, col = lcol, lwd = lwd)
    } else {
      # categorical: draw points at predicted means
      graphics::points(df_pred$x, df_pred$y, pch = 18, col = lcol, cex = 1.5)
    }
  }

  # draw observed points (bubbles)
  if (mod_type == "categorical") {
    # add jitter for categorical
    set.seed(42)  # for reproducibility
    x_jittered <- df_points$x + stats::runif(nrow(df_points), -jitter_am, jitter_am)
    graphics::points(x_jittered, df_points$y, pch = pch, col = col, bg = bg,
                     cex = df_points$size)
  } else {
    graphics::points(df_points$x, df_points$y, pch = pch, col = col, bg = bg,
                     cex = df_points$size)
  }

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .regplot_plot_ggplot
# ---------------------------------------------------------------------------- #
#
# Create regression plot using ggplot2.
#
# @param data   list of plot data from .regplot_data
# @param dots   list of graphical parameters
# @param legend logical; show legend
#
# @return ggplot object
#
# ---------------------------------------------------------------------------- #
.regplot_plot_ggplot <- function(data, dots, legend = TRUE) {

  # extract graphical parameters
  pch       <- dots[["pch"]]
  col       <- dots[["col"]]
  bg        <- dots[["bg"]]
  lcol      <- dots[["lcol"]]
  lwd       <- dots[["lwd"]]
  col_ci    <- dots[["col.ci"]]
  col_pi    <- dots[["col.pi"]]
  alpha_ci  <- dots[["alpha.ci"]]
  alpha_pi  <- dots[["alpha.pi"]]
  shade     <- dots[["shade"]]
  main      <- dots[["main"]]
  jitter_am <- dots[["jitter"]]
  size      <- dots[["size"]]

  # extract data components
  df_points  <- data$points
  df_pred    <- data$pred
  df_ci      <- data$ci
  df_pi      <- data$pi
  df_refline <- data$refline
  mod_type   <- data$mod_type

  out <- ggplot2::ggplot()

  # add reference line (if specified)
  if (!is.null(df_refline)) {
    out <- out +
      ggplot2::geom_hline(
        yintercept = df_refline$y,
        linetype   = "dashed",
        colour     = "gray50"
      )
  }

  # add PI band (if present and shading enabled)
  if (!is.null(df_pi) && shade) {
    if (mod_type == "continuous") {
      out <- out +
        ggplot2::geom_polygon(
          mapping = ggplot2::aes(x = df_pi$x, y = df_pi$y),
          fill    = col_pi,
          alpha   = alpha_pi
        )
    } else {
      # categorical: draw error bars
      out <- out +
        ggplot2::geom_errorbar(
          data    = df_pi,
          mapping = ggplot2::aes(x = x, ymin = lower, ymax = upper),
          width   = 0.2,
          colour  = col_pi,
          alpha   = alpha_pi + 0.3,
          linewidth = 2
        )
    }
  }

  # add CI band (if present and shading enabled)
  if (!is.null(df_ci) && shade) {
    if (mod_type == "continuous") {
      out <- out +
        ggplot2::geom_polygon(
          mapping = ggplot2::aes(x = df_ci$x, y = df_ci$y),
          fill    = col_ci,
          alpha   = alpha_ci
        )
    } else {
      # categorical: draw error bars
      out <- out +
        ggplot2::geom_errorbar(
          data    = df_ci,
          mapping = ggplot2::aes(x = x, ymin = lower, ymax = upper),
          width   = 0.1,
          colour  = col_ci,
          alpha   = alpha_ci + 0.3,
          linewidth = 1
        )
    }
  }

  # add prediction line
  if (!is.null(df_pred)) {
    if (mod_type == "continuous") {
      out <- out +
        ggplot2::geom_line(
          mapping   = ggplot2::aes(x = df_pred$x, y = df_pred$y),
          colour    = lcol,
          linewidth = lwd / 2
        )
    } else {
      # categorical: draw predicted means as diamonds
      out <- out +
        ggplot2::geom_point(
          mapping = ggplot2::aes(x = df_pred$x, y = df_pred$y),
          shape   = 18,
          colour  = lcol,
          size    = 4
        )
    }
  }

  # add observed points (bubbles)
  if (mod_type == "categorical") {
    # add jitter for categorical
    out <- out +
      ggplot2::geom_point(
        data    = df_points,
        mapping = ggplot2::aes(x = x, y = y, size = size),
        shape   = pch,
        colour  = col,
        fill    = bg,
        position = ggplot2::position_jitter(width = jitter_am, height = 0, seed = 42)
      )
  } else {
    out <- out +
      ggplot2::geom_point(
        data    = df_points,
        mapping = ggplot2::aes(x = x, y = y, size = size),
        shape   = pch,
        colour  = col,
        fill    = bg
      )
  }

  # configure size scale
  out <- out +
    ggplot2::scale_size_identity()

  # set axis scales and labels
  if (mod_type == "categorical") {
    out <- out +
      ggplot2::scale_x_continuous(
        breaks = seq_along(data$levels),
        labels = data$levels,
        limits = data$xlim,
        name   = data$xlab
      )
  } else {
    out <- out +
      ggplot2::scale_x_continuous(
        limits = data$xlim,
        name   = data$xlab
      )
  }

  out <- out +
    ggplot2::scale_y_continuous(
      limits = data$ylim,
      name   = data$ylab
    )

  # add title if specified
  if (!is.null(main)) {
    out <- out + ggplot2::ggtitle(main)
  }

  # remove legend for size (since we're using scale_size_identity)
  out <- out + ggplot2::guides(size = "none")

  return(out)
}


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
  if (is.null(dots[["pch"]]))      dots[["pch"]]      <- 21
  if (is.null(dots[["col"]]))      dots[["col"]]      <- "black"
  if (is.null(dots[["bg"]]))       dots[["bg"]]       <- "gray"
  if (is.null(dots[["size"]]))     dots[["size"]]     <- 2

  # line styling
  if (is.null(dots[["lcol"]]))     dots[["lcol"]]     <- "black"
  if (is.null(dots[["lwd"]]))      dots[["lwd"]]      <- 2

  # band styling
  if (is.null(dots[["shade"]]))    dots[["shade"]]    <- TRUE
  if (is.null(dots[["col.ci"]]))   dots[["col.ci"]]   <- "gray70"
  if (is.null(dots[["col.pi"]]))   dots[["col.pi"]]   <- "gray85"
  if (is.null(dots[["alpha.ci"]])) dots[["alpha.ci"]] <- 0.4
  if (is.null(dots[["alpha.pi"]])) dots[["alpha.pi"]] <- 0.2

  # categorical moderator jitter
  if (is.null(dots[["jitter"]]))   dots[["jitter"]]   <- 0.2

  # title (NULL = no title by default)
  if (is.null(dots[["main"]]))     dots[["main"]]     <- NULL

  return(dots)
}
