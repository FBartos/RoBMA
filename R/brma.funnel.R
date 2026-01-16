# ============================================================================ #
# brma.funnel.R
# ============================================================================ #
#
# This file contains the funnel plot functions for brma objects.
#
# ============================================================================ #


### funnel plot functions ----
#' @export
funnel <- function(x, ...) UseMethod("funnel")


#' @title Funnel Plot for brma Object
#'
#' @description \code{funnel.brma} creates a funnel plot for a fitted brma object,
#' displaying residuals against standard errors with the expected sampling
#' distribution funnel.
#'
#' @param x a fitted brma object
#' @param sampling_heterogeneity whether heterogeneity should be incorporated
#' into the sampling distribution funnel. Defaults to \code{TRUE}.
#' @param sampling_bias whether publication bias should be incorporated into the
#' sampling distribution funnel. Defaults to \code{TRUE}. When \code{TRUE} and
#' the model includes selection models (weightfunction), uses weighted normal
#' quantiles for the funnel contours. When \code{TRUE} and the model includes
#' PET/PEESE, incorporates the expected skew from these regression adjustments.
#' @param bias_adjusted whether residuals should be computed from bias-adjusted
#' fitted values. Defaults to \code{FALSE}, which means residuals represent
#' deviation from the raw (biased) effect estimate. When \code{TRUE},
#' residuals are computed from bias-corrected fitted values.
#' @param plot_type whether to use a base plot \code{"base"} or ggplot2
#' \code{"ggplot"} for plotting. Defaults to \code{"base"}.
#' @param ... additional graphical arguments to customize the plot appearance:
#' \describe{
#'   \item{xlim, ylim}{numeric vectors of length 2 specifying axis limits}
#'   \item{xlab, ylab}{character strings for axis labels (defaults: "Residual", "Standard Error")}
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
#' The funnel plot displays residuals (observed effect minus fitted value)
#' on the x-axis and standard errors on the y-axis. The funnel region represents
#' the 95% prediction region for residuals based on the sampling distribution.
#'
#' For meta-regression models, residuals account for moderator effects.
#'
#' When \code{sampling_heterogeneity = TRUE}, the funnel boundaries incorporate
#' between-study heterogeneity:
#' \deqn{CI = \pm 1.96 \times \sqrt{se^2 + \tau^2}}
#' where \eqn{\tau} is the estimated between-study heterogeneity.
#'
#' When \code{sampling_bias = TRUE} and the model includes publication bias
#' adjustment (selection models or PET/PEESE), the funnel contours reflect
#' the expected distortion in the sampling distribution due to bias. For
#' selection models, weighted normal quantiles are used. For PET/PEESE,
#' the regression-based skew is incorporated.
#'
#' The \code{bias_adjusted} parameter controls how residuals (the points) are
#' computed: when \code{TRUE}, residuals show deviation from the bias-corrected
#' pooled effect; when \code{FALSE}, they show deviation from raw predictions.
#'
#' For GLMM models, observed effect sizes are computed from the raw frequency
#' data using formulas equivalent to \code{metafor::escalc}.
#'
#' @return \code{funnel.brma} returns \code{NULL} invisibly if \code{plot_type = "base"}
#' or a ggplot object if \code{plot_type = "ggplot"}.
#'
#' @examples \dontrun{
#' # fit a brma model
#' fit <- brma(yi ~ 1, sei = sei, data = dat)
#'
#' # create funnel plot
#' funnel(fit)
#'
#' # without heterogeneity in sampling distribution
#' funnel(fit, sampling_heterogeneity = FALSE)
#'
#' # without bias in sampling distribution
#' funnel(fit, sampling_bias = FALSE)
#'
#' # residuals from bias-adjusted predictions
#' funnel(fit, bias_adjusted = TRUE)
#'
#' # customize point appearance
#' funnel(fit, pch = 19, col = "blue", bg = "lightblue", cex = 1.5)
#'
#' # customize funnel region styling
#' funnel(fit, back = "lightgrey", shade = "white", lty = "dashed")
#'
#' # custom axis labels and title
#' funnel(fit, xlab = "Effect Size Residual", main = "Funnel Plot")
#'
#' # custom axis limits
#' funnel(fit, xlim = c(-2, 2), ylim = c(0, 1))
#'
#' # using ggplot2
#' funnel(fit, plot_type = "ggplot")
#' funnel(fit, plot_type = "ggplot", pch = 19, size = 3, col = "red")
#' }
#'
#' @seealso [residuals.brma()], [predict.brma()]
#' @export
#' @rdname funnel
funnel.brma <- function(x, sampling_heterogeneity = TRUE, sampling_bias = TRUE,
                        bias_adjusted = FALSE, plot_type = "base", ...) {

  # input validation
  BayesTools::check_bool(sampling_heterogeneity, "sampling_heterogeneity")
  BayesTools::check_bool(sampling_bias, "sampling_bias")
  BayesTools::check_bool(bias_adjusted, "bias_adjusted")
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))

  # set up graphical arguments with defaults
  dots <- .set_dots_funnel(...)

  # extract graphical parameters
  pch         <- dots[["pch"]]
  col         <- dots[["col"]]
  bg          <- dots[["bg"]]
  cex         <- dots[["cex"]]
  size        <- dots[["size"]]
  back        <- dots[["back"]]
  shade       <- dots[["shade"]]
  lty         <- dots[["lty"]]
  col_line    <- dots[["col.line"]]
  refline     <- dots[["refline"]]
  col_refline <- dots[["col.refline"]]
  xlab        <- dots[["xlab"]]
  ylab        <- dots[["ylab"]]
  main        <- dots[["main"]]
  xlim        <- dots[["xlim"]]
  ylim        <- dots[["ylim"]]
  as_data     <- dots[["as_data"]]

  # get model characteristics
  outcome_type      <- .outcome_type(x)
  is_multilevel     <- .is_multilevel(x)
  is_weightfunction <- .is_weightfunction(x)
  is_PET            <- .is_PET(x)
  is_PEESE          <- .is_PEESE(x)
  effect_direction  <- .effect_direction(x)

  # get residuals (yi - fitted) as mean of posterior samples
  res <- residuals.brma(x, type = "outcome", bias_adjusted = bias_adjusted)

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

  # compute confidence intervals for the funnel based on sampling distribution
  ci_bounds <- .get_funnel_quantiles(
    x                = x,
    se_sequence      = se_sequence,
    tau              = tau,
    sampling_bias    = sampling_bias,
    is_weightfunction = is_weightfunction,
    is_PET           = is_PET,
    is_PEESE         = is_PEESE,
    effect_direction = effect_direction
  )
  ci_left  <- ci_bounds$lower
  ci_right <- ci_bounds$upper

  # create data frames for plotting
  x_range <- pretty(c(ci_left, ci_right))

  # apply custom axis limits if provided
  if (!is.null(xlim)) {
    x_range <- pretty(xlim)
  }
  if (!is.null(ylim)) {
    se_range <- pretty(ylim)
  }

  df_points <- data.frame(
    x = res,
    y = se
  )
  df_funnel <- data.frame(
    x = c(rev(ci_left), ci_right),
    y = c(rev(se_sequence), se_sequence)
  )
  df_funnel_edge1 <- data.frame(
    x = ci_left,
    y = se_sequence
  )
  df_funnel_edge2 <- data.frame(
    x = ci_right,
    y = se_sequence
  )
  df_background <- data.frame(
    x = c(min(x_range), max(x_range), max(x_range), min(x_range)),
    y = c(min(se_range), min(se_range), max(se_range), max(se_range))
  )

  # allow data return for programmatic access
  if (isTRUE(as_data)) {
    return(list(
      points       = df_points,
      funnel       = df_funnel,
      funnel_edge1 = df_funnel_edge1,
      funnel_edge2 = df_funnel_edge2,
      background   = df_background,
      x_range      = x_range,
      y_range      = se_range
    ))
  }

  # create plot
  if (plot_type == "ggplot") {

    out <- ggplot2::ggplot()

    # add background polygon (if not suppressed)
    if (!is.na(back)) {
      out <- out +
        ggplot2::geom_polygon(
          mapping = ggplot2::aes(
            x = df_background$x,
            y = df_background$y),
          fill = back
        )
    }

    # add funnel polygon (if not suppressed)
    if (!is.na(shade)) {
      out <- out +
        ggplot2::geom_polygon(
          mapping = ggplot2::aes(
            x = df_funnel$x,
            y = df_funnel$y),
          fill = shade
        )
    }

    # add vertical reference line
    out <- out +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = c(refline, refline),
          y = range(se_range)),
        linetype = lty,
        colour   = col_refline
      )

    # add funnel edge lines
    out <- out +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = df_funnel_edge1$x,
          y = df_funnel_edge1$y),
        linetype = lty,
        colour   = col_line
      ) +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = df_funnel_edge2$x,
          y = df_funnel_edge2$y),
        linetype = lty,
        colour   = col_line
      )

    # add points
    out <- out +
      ggplot2::geom_point(
        mapping = ggplot2::aes(
          x = df_points$x,
          y = df_points$y),
        colour = col,
        fill   = bg,
        shape  = pch,
        size   = size
      )

    # set axis scales and labels
    out <- out +
      ggplot2::scale_x_continuous(
        breaks = x_range,
        limits = range(x_range),
        name   = xlab
      ) +
      ggplot2::scale_y_reverse(
        breaks = rev(se_range),
        limits = rev(range(se_range)),
        name   = ylab
      )

    # add title if specified
    if (!is.null(main)) {
      out <- out + ggplot2::ggtitle(main)
    }

    return(out)

  } else if (plot_type == "base") {

    # set up the plot area with reversed y-axis
    graphics::plot(
      NA, NA,
      xlim = range(x_range),
      ylim = rev(range(se_range)),
      xlab = xlab,
      ylab = ylab,
      main = main,
      type = "n",
      axes = FALSE
    )
    graphics::axis(1, at = x_range)
    graphics::axis(2, at = rev(se_range), las = 1)

    # draw background polygon (if not suppressed)
    if (!is.na(back)) {
      graphics::polygon(df_background$x, df_background$y, col = back, border = NA)
    }

    # draw funnel polygon (if not suppressed)
    if (!is.na(shade)) {
      graphics::polygon(df_funnel$x, df_funnel$y, col = shade, border = NA)
    }

    # vertical reference line
    graphics::lines(c(refline, refline), range(se_range), lty = lty, col = col_refline)

    # funnel edge lines
    graphics::lines(df_funnel_edge1$x, df_funnel_edge1$y, lty = lty, col = col_line)
    graphics::lines(df_funnel_edge2$x, df_funnel_edge2$y, lty = lty, col = col_line)

    # plot points
    graphics::points(df_points$x, df_points$y, pch = pch, col = col, bg = bg, cex = cex)

    return(invisible(NULL))
  }
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
  tau_summary <- pooled_heterogeneity(object)
  return(tau_summary$summary["tau", "Mean"])
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
.get_funnel_quantiles <- function(x, se_sequence, tau, sampling_bias,
                                   is_weightfunction, is_PET, is_PEESE,
                                   effect_direction) {

  n_se   <- length(se_sequence)
  sd_seq <- sqrt(se_sequence^2 + tau^2)

  # default: standard normal quantiles centered at 0
  if (!sampling_bias || (!is_weightfunction && !is_PET && !is_PEESE)) {
    return(list(
      lower = stats::qnorm(0.025, mean = 0, sd = sd_seq),
      upper = stats::qnorm(0.975, mean = 0, sd = sd_seq)
    ))
  }

  # sampling_bias = TRUE with bias model present
  if (is_weightfunction) {
    # use weighted normal quantiles for selection models
    return(.get_funnel_quantiles_weighted(x, se_sequence, sd_seq, effect_direction))

  } else if (is_PET || is_PEESE) {
    # for PET/PEESE, incorporate regression-based skew
    return(.get_funnel_quantiles_PETPEESE(x, se_sequence, sd_seq, is_PET, is_PEESE, effect_direction))
  }

  # fallback to normal quantiles
  return(list(
    lower = stats::qnorm(0.025, mean = 0, sd = sd_seq),
    upper = stats::qnorm(0.975, mean = 0, sd = sd_seq)
  ))
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
# @param sd_seq           numeric vector of total SD (sqrt(se^2 + tau^2))
# @param effect_direction character; "positive" or "negative"
#
# @return list with 'lower' and 'upper' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles_weighted <- function(x, se_sequence, sd_seq, effect_direction) {

  n_se <- length(se_sequence)

  # extract posterior mean omega from fitted model
  posterior_samples <- suppressWarnings(coda::as.mcmc(x[["fit"]]))
  omega_cols        <- grep("omega", colnames(posterior_samples))
  omega_mean        <- colMeans(posterior_samples[, omega_cols, drop = FALSE])
  omega_matrix      <- matrix(omega_mean, nrow = n_se, ncol = length(omega_mean), byrow = TRUE)

  # get critical values from fit data
  # for funnel plot, we use the standard critical values at typical alpha levels
  fit_data <- x[["fit_data"]]
  if (!is.null(fit_data) && !is.null(fit_data$crit_yi)) {
    # use first row as template (same p-value cutoffs for all observations)
    crit_yi <- fit_data$crit_yi[1, , drop = FALSE]
  } else {
    # fallback: use standard two-sided critical values
    # typical selection model uses one-sided p-values: 0.025, 0.05, 0.5, 1
    crit_yi <- matrix(c(stats::qnorm(0.975), stats::qnorm(0.95), 0), nrow = 1)
  }

  lower <- numeric(n_se)
  upper <- numeric(n_se)

  # compute weighted quantiles for each SE level
  # note: weighted distribution is computed in "positive" effect space
  for (i in seq_len(n_se)) {
    # for residuals centered at 0, compute quantiles
    lower[i] <- .qwnorm_fast.ss(
      p      = 0.025,
      mean   = 0,
      sd     = sd_seq[i],
      omega  = omega_matrix[i, , drop = FALSE],
      crit_x = crit_yi * se_sequence[i]  # scale critical values by SE
    )
    upper[i] <- .qwnorm_fast.ss(
      p      = 0.975,
      mean   = 0,
      sd     = sd_seq[i],
      omega  = omega_matrix[i, , drop = FALSE],
      crit_x = crit_yi * se_sequence[i]
    )
  }

  # flip back if effect direction is negative
  if (effect_direction == "negative") {
    temp  <- -upper
    upper <- -lower
    lower <- temp
  }

  return(list(lower = lower, upper = upper))
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
# @param sd_seq           numeric vector of total SD (sqrt(se^2 + tau^2))
# @param is_PET           logical; whether model has PET
# @param is_PEESE         logical; whether model has PEESE
# @param effect_direction character; "positive" or "negative"
#
# @return list with 'lower' and 'upper' quantile vectors
#
# ---------------------------------------------------------------------------- #
.get_funnel_quantiles_PETPEESE <- function(x, se_sequence, sd_seq, is_PET, is_PEESE, effect_direction) {

  n_se <- length(se_sequence)

  # extract posterior mean PET/PEESE coefficients
  posterior_samples <- suppressWarnings(coda::as.mcmc(x[["fit"]]))

  # compute bias shift at each SE level
  bias_shift <- rep(0, n_se)
  direction  <- ifelse(effect_direction == "negative", -1, 1)

  if (is_PET && "PET" %in% colnames(posterior_samples)) {
    PET_mean   <- mean(posterior_samples[, "PET"])
    bias_shift <- bias_shift + direction * PET_mean * se_sequence
  }

  if (is_PEESE && "PEESE" %in% colnames(posterior_samples)) {
    PEESE_mean <- mean(posterior_samples[, "PEESE"])
    bias_shift <- bias_shift + direction * PEESE_mean * se_sequence^2
  }

  # the funnel is centered at residual = 0, but bias shifts the expected
  # sampling distribution, so we offset the quantiles by the bias amount
  lower <- stats::qnorm(0.025, mean = bias_shift, sd = sd_seq)
  upper <- stats::qnorm(0.975, mean = bias_shift, sd = sd_seq)

  return(list(lower = lower, upper = upper))
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
  if (is.null(dots[["pch"]]))  dots[["pch"]]  <- 21       # open circle with fill
  if (is.null(dots[["col"]]))  dots[["col"]]  <- "black"  # point border color
  if (is.null(dots[["bg"]]))   dots[["bg"]]   <- "black"  # point fill color
  if (is.null(dots[["cex"]]))  dots[["cex"]]  <- 1        # point size (base)
  if (is.null(dots[["size"]])) dots[["size"]] <- 2        # point size (ggplot)

  # funnel region styling
  if (is.null(dots[["back"]]))  dots[["back"]]  <- "grey"   # background fill

  if (is.null(dots[["shade"]])) dots[["shade"]] <- "white"  # funnel fill

  # line styling
  if (is.null(dots[["lty"]]))      dots[["lty"]]      <- "dotted"  # line type
  if (is.null(dots[["col.line"]])) dots[["col.line"]] <- "black"   # funnel edge color

  # reference line
  if (is.null(dots[["refline"]]))     dots[["refline"]]     <- 0        # vertical line position
  if (is.null(dots[["col.refline"]])) dots[["col.refline"]] <- "black"  # reference line color

  # axis labels
  if (is.null(dots[["xlab"]])) dots[["xlab"]] <- "Residual"
  if (is.null(dots[["ylab"]])) dots[["ylab"]] <- "Standard Error"
  if (is.null(dots[["main"]])) dots[["main"]] <- NULL  # no title by default

  # axis limits (NULL = auto)
  # xlim and ylim are left as NULL if not provided (auto-computed)

  # data return flag
  if (is.null(dots[["as_data"]])) dots[["as_data"]] <- FALSE

  return(dots)
}