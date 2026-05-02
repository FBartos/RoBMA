# ============================================================================ #
# brma.zcurve.R
# ============================================================================ #
#
# Z-curve diagnostics for brma class objects.
#
# Z-curve analysis provides estimates of Expected Discovery Rate (EDR) and
# related metrics for assessing replicability of meta-analytic findings.
# The implementation follows the framework from Bartos & Schimmack (2022) and
# extends it to Bayesian meta-analytic models.
#
# Design principles:
# - as_zcurve() transforms brma to zcurve_brma, enabling summary/plot methods
# - Extrapolation mode removes publication bias to estimate "true" power
# - Fitted mode includes bias adjustments for observed distribution
# - All computations are sample-based to propagate posterior uncertainty
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# as_zcurve generic and method
# ---------------------------------------------------------------------------- #

#' @export
as_zcurve <- function(object, ...) UseMethod("as_zcurve")


#' @title Transform brma Object to Z-Curve
#'
#' @description Transforms an estimated brma model into a z-curve object
#' that can be summarized and plotted to assess replicability.
#'
#' @param object a brma object.
#' @param significance_level z-value threshold for significance. Defaults
#' to \code{qnorm(0.975)} (two-sided alpha = 0.05).
#' @param max_samples maximum number of posterior samples for estimation.
#' Defaults to 1000.
#' @param ... additional arguments (currently unused).
#'
#' @details
#' Z-curve analysis estimates the Expected Discovery Rate (EDR), which
#' represents the average power of statistically significant studies. This
#' provides insight into the replicability of findings in a literature.
#'
#' The implementation uses extrapolation mode by default, which removes
#' publication bias adjustments (PET/PEESE/selection weights) to estimate
#' what the true power distribution would be without selective reporting.
#'
#' The resulting object retains all original brma properties while adding
#' z-curve results, enabling both standard meta-analytic summaries and
#' z-curve diagnostics on the same object.
#'
#' @return The input object with added class \code{"zcurve_brma"} and a new
#' \code{zcurve} list component containing:
#' \describe{
#'   \item{estimates}{posterior samples for EDR and weights}
#'   \item{data}{observed z-statistics and significance counts}
#' }
#'
#' @seealso [summary.zcurve_brma()], [plot.zcurve_brma()], [hist.zcurve_brma()]
#'
#' @export
as_zcurve.brma <- function(object, significance_level = stats::qnorm(0.975), max_samples = 1000, ...) {

  BayesTools::check_real(significance_level, "significance_level", lower = 0)
  BayesTools::check_int(max_samples, "max_samples", lower = 10)

  # compute the main estimates (extrapolate = TRUE: unbiased power)
  # This corresponds to "Expected Discovery Rate" if studies were exact replications
  # but without publication bias
  out_estimates <- .zcurve_fun.brma(
    object      = object,
    z_threshold = significance_level,
    max_samples = max_samples,
    extrapolate = TRUE
  )

  # compute data summaries (observed z-stats)
  yi  <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)
  z   <- yi / sei

  # store z-curve results
  object[["zcurve"]] <- list(
    estimates = list(
      EDR     = out_estimates[["EDR"]],
      weights = out_estimates[["weights"]]
    ),
    data = list(
      significance_level = significance_level,
      z                  = z,
      N_significant      = sum(abs(z) > significance_level),
      N_observed         = length(z)
    )
  )

  # add class
  class(object) <- c("zcurve_brma", class(object))

  return(object)
}


# ---------------------------------------------------------------------------- #
# summary.zcurve_brma
# ---------------------------------------------------------------------------- #

#' @title Summarize Z-Curve Results
#'
#' @description Creates summary tables for z-curve estimates including
#' EDR, Soric FDR, and estimated missing studies.
#'
#' @param object a zcurve_brma object.
#' @param probs quantiles of the posterior distribution to display.
#' Defaults to \code{c(.025, .975)}.
#' @param ... additional arguments (currently unused).
#'
#' @details
#' The summary includes:
#' \describe{
#'   \item{EDR}{Expected Discovery Rate - average power of significant studies}
#'   \item{Soric FDR}{Expected proportion of false discoveries among significant
#'     results, computed from EDR following Soric (1989)}
#'   \item{Missing N}{Estimated number of studies suppressed by publication bias,
#'     computed from the selection model weights}
#' }
#'
#' The footer reports the Observed Discovery Rate (ODR) with 95\% CI for
#' comparison with the model-estimated EDR.
#'
#' @return An object of class \code{"summary.zcurve_brma"} containing the
#' estimates table.
#'
#' @seealso [as_zcurve.brma()], [plot.zcurve_brma()]
#'
#' @export
summary.zcurve_brma <- function(object, probs = c(.025, .975), ...) {

  # get proportion of significant results
  obs_proportion <- stats::prop.test(
    object$zcurve$data[["N_significant"]],
    object$zcurve$data[["N_observed"]],
    conf.level = 0.95
  )

  info_text <- sprintf(
    "Estimated using %1$i estimates, %2$i significant (ODR = %3$.2f, 95%% CI [%4$.2f, %5$.2f]).",
    object$zcurve$data[["N_observed"]],
    object$zcurve$data[["N_significant"]],
    obs_proportion$estimate,
    obs_proportion$conf.int[1],
    obs_proportion$conf.int[2]
  )

  # compute estimates
  sig_level <- stats::pnorm(object$zcurve$data[["significance_level"]], lower.tail = FALSE) * 2

  estimates <- cbind.data.frame(
    "EDR"       = object$zcurve$estimates[["EDR"]],
    "Soric FDR" = .get_Soric_FDR(object$zcurve$estimates[["EDR"]], sig_level),
    "Missing N" = (object$zcurve$estimates[["weights"]] - 1) * object$zcurve$data[["N_observed"]]
  )

  estimates_table <- BayesTools::ensemble_estimates_table(
    samples    = estimates,
    parameters = names(estimates),
    probs      = probs,
    title      = "Z-Curve Estimates:",
    footnotes  = info_text
  )

  output <- list(
    estimates = estimates_table
  )

  class(output) <- "summary.zcurve_brma"
  return(output)
}


#' @title Print Z-Curve Summary
#'
#' @param x a summary.zcurve_brma object.
#' @param ... additional arguments (currently unused).
#'
#' @return Invisibly returns the summary object.
#'
#' @export
print.summary.zcurve_brma <- function(x, ...) {

  cat("\n")
  print(x[["estimates"]])
  cat("\n")

  return(invisible(x))
}


#' @export
print.zcurve_brma <- function(x, ...) {
  print(summary(x, ...))
}


# ---------------------------------------------------------------------------- #
# plot.zcurve_brma
# ---------------------------------------------------------------------------- #

#' @title Plot Z-Curve Results
#'
#' @description Plots a z-curve visualization showing the histogram of
#' observed z-statistics overlaid with model-implied densities.
#'
#' @param x a zcurve_brma object.
#' @param plot_type graphics system: \code{"base"} or \code{"ggplot"}.
#' Defaults to \code{"base"}.
#' @param probs quantiles for credible intervals. Defaults to \code{c(.025, .975)}.
#' @param max_samples maximum posterior samples for density estimation.
#' Defaults to 500.
#' @param plot_fit whether to show fitted density (with bias adjustments).
#' Defaults to \code{TRUE}.
#' @param plot_extrapolation whether to show extrapolated density (bias removed).
#' Defaults to \code{TRUE}.
#' @param plot_ci whether to show credible interval bands. Defaults to \code{TRUE}.
#' @param plot_thresholds whether to show significance threshold lines.
#' Defaults to \code{TRUE}.
#' @param from,to z-value range for plotting. Defaults to \code{-6} and \code{6}.
#' @param by.hist bin width for histogram. Defaults to 0.5.
#' @param length.out.hist number of histogram bins (alternative to \code{by.hist}).
#' @param by.lines step size for density lines. Defaults to 0.05.
#' @param length.out.lines number of density points (alternative to \code{by.lines}).
#' @param dots_hist graphical parameters for histogram (list).
#' @param dots_fit graphical parameters for fit lines (list).
#' @param dots_extrapolation graphical parameters for extrapolation lines (list).
#' @param dots_thresholds graphical parameters for threshold lines (list).
#' @param ... additional graphical parameters passed to components.
#'
#' @details
#' The plot displays two density curves:
#' \describe{
#'   \item{Fit (black)}{Model-implied density including publication bias adjustments.
#'     This represents the expected distribution of z-statistics given the estimated
#'     selection process.}
#'   \item{Extrapolation (blue)}{Bias-corrected curve representing the hypothetical
#'     distribution without selective reporting, scaled by the expected suppressed
#'     studies under selection models. This curve is not normalized to integrate
#'     to one when selection implies missing studies.}
#' }
#'
#' @return \code{NULL} invisibly for base graphics, or a ggplot2 object.
#'
#' @seealso [hist.zcurve_brma()], [lines.zcurve_brma()], [summary.zcurve_brma()]
#'
#' @export
plot.zcurve_brma <- function(x, plot_type = "base",
                             probs = c(.025, .975), max_samples = 500,
                             plot_fit = TRUE, plot_extrapolation = TRUE,
                             plot_ci = TRUE, plot_thresholds = TRUE,
                             from = -6, to = 6,
                             by.hist = 0.5, length.out.hist = NULL,
                             by.lines = 0.05, length.out.lines = NULL,
                             dots_hist = NULL, dots_fit = NULL,
                             dots_extrapolation = NULL, dots_thresholds = NULL, ...) {

  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(plot_fit, "plot_fit")
  BayesTools::check_bool(plot_extrapolation, "plot_extrapolation")
  BayesTools::check_bool(plot_ci, "plot_ci")
  BayesTools::check_bool(plot_thresholds, "plot_thresholds")
  BayesTools::check_real(probs, "probs", lower = 0, upper = 1, check_length = 2)

  dots <- list(...)

  # get line values first so we can set ylim
  ymax <- 0

  # 1. Fit lines (Fitted model - with bias)
  if (plot_fit) {
    dots_fit <- if(!is.null(dots_fit)) dots_fit else list()
    lines_fit <- do.call(lines.zcurve_brma, c(
      list(x = x, plot_type = plot_type, probs = probs, max_samples = max_samples,
           extrapolate = FALSE, plot_ci = plot_ci, from = from, to = to,
           by = by.lines, length.out = length.out.lines, as_data = TRUE),
      dots_fit, dots
    ))
    ymax <- max(c(ymax, lines_fit$y, if (plot_ci) lines_fit$y_uCI))
  }

  # 2. Extrapolation lines (Unbiased model - "True" z-curve)
  if (plot_extrapolation) {
    # prepare extrapolation dots with default blue color
    dots_extrapolation <- if(!is.null(dots_extrapolation)) dots_extrapolation else list()
    lines_extrapolation <- do.call(lines.zcurve_brma, c(
      list(x = x, plot_type = plot_type, probs = probs, max_samples = max_samples,
           extrapolate = TRUE, plot_ci = plot_ci, from = from, to = to,
           by = by.lines, length.out = length.out.lines, as_data = TRUE,
           col = "blue"), # default color if not in dots
      dots_extrapolation, dots
    ))
    ymax <- max(c(ymax, lines_extrapolation$y, if (plot_ci) lines_extrapolation$y_uCI))
  }

  # 3. Create histogram base
  plot <- hist.zcurve_brma(
    x               = x,
    plot_type       = plot_type,
    from            = from,
    to              = to,
    by              = by.hist,
    length.out      = length.out.hist,
    plot_thresholds = plot_thresholds,
    dots_thresholds = dots_thresholds,
    ylim            = if (ymax != 0) c(0, ymax * 1.05) else NULL,
    dots_hist       = dots_hist,
    dots_all        = dots
  )

  # 4. Add fit lines
  if (plot_fit) {
    dots_fit_lines <- .get_dots_lines_zcurve(c(dots, dots_fit), plot_type = plot_type)

    if (plot_type == "base") {
      if (plot_ci) {
        graphics::polygon(
          c(lines_fit$x, rev(lines_fit$x)),
          c(lines_fit$y_lCI, rev(lines_fit$y_uCI)),
          border = NA,
          col    = scales::alpha(dots_fit_lines$col, dots_fit_lines$alpha)
        )
      }
      graphics::lines(
        lines_fit$x, lines_fit$y,
        lwd = dots_fit_lines$lwd,
        col = dots_fit_lines$col,
        lty = dots_fit_lines$lty
      )
    } else if (plot_type == "ggplot") {
      if (plot_ci) {
        plot <- plot + ggplot2::geom_ribbon(
          ggplot2::aes(
            x    = lines_fit$x,
            ymin = lines_fit$y_lCI,
            ymax = lines_fit$y_uCI
          ),
          fill  = dots_fit_lines$color,
          alpha = dots_fit_lines$alpha
        )
      }
      plot <- plot + ggplot2::geom_line(
        ggplot2::aes(
          x = lines_fit$x,
          y = lines_fit$y
        ),
        color     = dots_fit_lines$color,
        linewidth = dots_fit_lines$linewidth,
        linetype  = dots_fit_lines$linetype
      )
    }
  }

  # 5. Add extrapolation lines
  if (plot_extrapolation) {
     dots_ext_lines <- .get_dots_lines_zcurve(c(dots, dots_extrapolation), plot_type = plot_type, col = "blue")

    if (plot_type == "base") {
      if (plot_ci) {
        graphics::polygon(
          c(lines_extrapolation$x, rev(lines_extrapolation$x)),
          c(lines_extrapolation$y_lCI, rev(lines_extrapolation$y_uCI)),
          border = NA,
          col    = scales::alpha(dots_ext_lines$col, dots_ext_lines$alpha)
        )
      }
      graphics::lines(
        lines_extrapolation$x, lines_extrapolation$y,
        lwd = dots_ext_lines$lwd,
        col = dots_ext_lines$col,
        lty = dots_ext_lines$lty
      )
    } else if (plot_type == "ggplot") {
      if (plot_ci) {
        plot <- plot + ggplot2::geom_ribbon(
          ggplot2::aes(
            x    = lines_extrapolation$x,
            ymin = lines_extrapolation$y_lCI,
            ymax = lines_extrapolation$y_uCI
          ),
          fill  = dots_ext_lines$color,
          alpha = dots_ext_lines$alpha
        )
      }
      plot <- plot + ggplot2::geom_line(
        ggplot2::aes(
          x = lines_extrapolation$x,
          y = lines_extrapolation$y
        ),
        color     = dots_ext_lines$color,
        linewidth = dots_ext_lines$linewidth,
        linetype  = dots_ext_lines$linetype
      )
    }
  }

  # return
  if (plot_type == "base") {
    return(invisible())
  } else if (plot_type == "ggplot") {
    return(plot)
  }
}


# ---------------------------------------------------------------------------- #
# hist.zcurve_brma
# ---------------------------------------------------------------------------- #

#' @title Histogram of Z-Statistics
#'
#' @description Plots a histogram of observed z-values from the meta-analysis.
#'
#' @param x a zcurve_brma object.
#' @param plot_type graphics system: \code{"base"} or \code{"ggplot"}.
#' Defaults to \code{"base"}.
#' @param from,to z-value range for plotting. Defaults to \code{-6} and \code{6}.
#' @param by bin width. Defaults to 0.5.
#' @param length.out number of bins (alternative to \code{by}).
#' @param add whether to add to existing plot. Defaults to \code{FALSE}.
#' @param plot_thresholds whether to show significance threshold lines.
#' Defaults to \code{TRUE}.
#' @param dots_thresholds graphical parameters for threshold lines (list).
#' @param dots_hist graphical parameters for histogram (list).
#' @param dots_all graphical parameters passed to all components (list).
#' @param ... additional graphical parameters.
#'
#' @details
#' Z-statistics are computed as effect size divided by standard error
#' (\code{yi / sei}). Histogram bins are adjusted to align with significance
#' thresholds when selection model priors are present.
#'
#' @return \code{NULL} invisibly for base graphics, or a ggplot2 object.
#'
#' @seealso [plot.zcurve_brma()], [lines.zcurve_brma()]
#'
#' @export
hist.zcurve_brma <- function(x, plot_type = "base",
                             from = -6, to = 6, by = 0.5, length.out = NULL,
                             add = FALSE, plot_thresholds = TRUE,
                             dots_thresholds = NULL, dots_hist = NULL, dots_all = NULL, ...) {

  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_real(from, "from")
  BayesTools::check_real(to, "to")
  BayesTools::check_real(by, "by", allow_NULL = TRUE)
  BayesTools::check_bool(add, "add")
  BayesTools::check_bool(plot_thresholds, "plot_thresholds")

  dots <- list(...)
  if (!is.null(dots_all)) {
    dots <- c(dots, dots_all[!names(dots_all) %in% names(dots)])
  }

  z <- x$zcurve$data[["z"]]

  # create bin sequence using z-curve helper (handles thresholds)
  z_sequence <- .zcurve_bins(priors = x[["priors"]], from = from, to = to, by = by, length.out = length.out, type = "hist")
  z_in_range <- z >= min(z_sequence) & z <= max(z_sequence)

  # message about out-of-range values
  if (sum(!z_in_range) > 0) {
    message(sprintf("%1$i z-statistics are out of the plotting range", sum(!z_in_range)))
  }

  # create histogram data
  z_hist <- graphics::hist(z[z_in_range], breaks = z_sequence, plot = FALSE)
  z_hist$density <- z_hist$density * mean(z_in_range)

  df_hist <- data.frame(
    x       = z_hist$mids,
    density = z_hist$density,
    breaks  = diff(z_hist$breaks)
  )

  # return data if requested via dots
  if (isTRUE(dots[["as_data"]])) {
    return(df_hist)
  }

  # merge dots_hist
  if (!is.null(dots_hist)) {
    dots <- c(dots_hist, dots[!names(dots) %in% names(dots_hist)])
  }
  dots_hist_params <- .get_dots_hist_zcurve(dots, plot_type = plot_type, max_density = max(z_hist$density))

  # create the plot
  if (plot_type == "ggplot") {
    out <- ggplot2::ggplot() +
      ggplot2::geom_col(
        ggplot2::aes(
          x = df_hist$x,
          y = df_hist$density
        ),
        fill  = dots_hist_params$fill,
        color = dots_hist_params$color,
        alpha = dots_hist_params$alpha,
        width = df_hist$breaks
      ) +
      ggplot2::labs(x = dots_hist_params$xlab, y = dots_hist_params$ylab) +
      ggplot2::ggtitle(dots_hist_params$main)
  } else {
    if (add) {
      graphics::rect(
        xleft   = df_hist$x - df_hist$breaks / 2,
        xright  = df_hist$x + df_hist$breaks / 2,
        ybottom = 0,
        ytop    = df_hist$density,
        border  = dots_hist_params$border,
        col     = dots_hist_params$col
      )
    } else {
      graphics::plot(
        z_hist,
        freq   = FALSE,
        border = dots_hist_params$border,
        col    = dots_hist_params$col,
        xlab   = dots_hist_params$xlab,
        ylab   = dots_hist_params$ylab,
        main   = dots_hist_params$main,
        ylim   = dots_hist_params$ylim,
        xaxt   = dots_hist_params$xaxt,
        yaxt   = dots_hist_params$yaxt,
        las    = dots_hist_params$las
      )
    }
  }

  # add threshold lines
  if (plot_thresholds) {
    thresholds <- .zcurve_threshold(x[["priors"]])
    if (length(thresholds) > 0) {
      dots_thresholds <- if (!is.null(dots_thresholds)) dots_thresholds else list()
      dots_thresholds_params <- .get_dots_thresholds_zcurve(dots_thresholds, plot_type = plot_type)

      if (plot_type == "base") {
        graphics::abline(
          v   = thresholds,
          col = dots_thresholds_params$col,
          lty = dots_thresholds_params$lty,
          lwd = dots_thresholds_params$lwd
        )
      } else if (plot_type == "ggplot") {
        out <- out + ggplot2::geom_vline(
          xintercept = thresholds,
          color      = dots_thresholds_params$color,
          linetype   = dots_thresholds_params$linetype,
          linewidth  = dots_thresholds_params$linewidth
        )
      }
    }
  }

  # return
  if (plot_type == "base") {
    return(invisible())
  } else if (plot_type == "ggplot") {
    return(out)
  }
}


# ---------------------------------------------------------------------------- #
# lines.zcurve_brma
# ---------------------------------------------------------------------------- #

#' @title Add Z-Curve Density Lines
#'
#' @description Adds model-implied density lines to an existing z-curve plot.
#'
#' @param x a zcurve_brma object.
#' @param plot_type graphics system: \code{"base"} or \code{"ggplot"}.
#' Defaults to \code{"base"}.
#' @param probs quantiles for credible intervals. Defaults to \code{c(.025, .975)}.
#' @param max_samples maximum posterior samples for density. Defaults to 500.
#' @param plot_ci whether to show credible interval bands. Defaults to \code{TRUE}.
#' @param extrapolate whether to remove bias adjustments. Defaults to \code{FALSE}.
#' @param from,to z-value range for density. Defaults to \code{-6} and \code{6}.
#' @param by step size for density points. Defaults to 0.05.
#' @param length.out number of density points (alternative to \code{by}).
#' @param col line color. Defaults to \code{"black"}.
#' @param as_data whether to return data instead of plotting. Defaults to \code{FALSE}.
#' @param ... additional graphical parameters.
#'
#' @details
#' When \code{extrapolate = FALSE}, the density includes all bias adjustments
#' (PET/PEESE regression, selection weights) representing the fitted model.
#' When \code{extrapolate = TRUE}, bias adjustments are removed to show the
#' hypothetical distribution without publication bias. Under selection models,
#' the curve is scaled by inverse selection probability and need not integrate
#' to one.
#'
#' @return \code{NULL} invisibly for base graphics, ggplot2 layers for ggplot,
#' or a data frame with columns \code{x}, \code{y}, \code{y_lCI}, \code{y_uCI}
#' if \code{as_data = TRUE}.
#'
#' @seealso [plot.zcurve_brma()], [hist.zcurve_brma()]
#'
#' @export
lines.zcurve_brma <- function(x, plot_type = "base",
                              probs = c(.025, .975), max_samples = 500,
                              plot_ci = TRUE, extrapolate = FALSE,
                              from = -6, to = 6, by = 0.05, length.out = NULL,
                              col = "black", as_data = FALSE, ...) {

  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_real(probs, "probs", lower = 0, upper = 1, check_length = 2)
  BayesTools::check_int(max_samples, "max_samples", lower = 10)
  BayesTools::check_bool(plot_ci, "plot_ci")
  BayesTools::check_bool(extrapolate, "extrapolate")
  BayesTools::check_real(from, "from")
  BayesTools::check_real(to, "to")
  BayesTools::check_real(by, "by", allow_NULL = TRUE)


  dots <- list(...)

  # prepare z sequence
  z_sequence <- .zcurve_bins(priors = x[["priors"]], from = from, to = to, by = by, length.out = length.out, type = "dens")

  # get densities from model
  z_density <- .zcurve_fun.brma(
     object      = x,
     z_sequence  = z_sequence,
     max_samples = max_samples,
     extrapolate = extrapolate
  )

  df_density <- data.frame(
    x     = z_sequence,
    y     = colMeans(z_density),
    y_lCI = apply(z_density, 2, stats::quantile, probs = probs[1]),
    y_uCI = apply(z_density, 2, stats::quantile, probs = probs[2])
  )

  # return data if requested
  if (as_data) {
    return(df_density)
  }

  # get line-specific parameters
  dots_lines <- .get_dots_lines_zcurve(dots, plot_type = plot_type, col = col)

  if (plot_type == "ggplot") {
    out <- list()
    if (plot_ci) {
      out[[1]] <- ggplot2::geom_ribbon(
        ggplot2::aes(
          x    = df_density$x,
          ymin = df_density$y_lCI,
          ymax = df_density$y_uCI
        ),
        fill  = dots_lines$color,
        alpha = dots_lines$alpha
      )
    }
    out[[length(out) + 1]] <- ggplot2::geom_line(
      ggplot2::aes(
        x = df_density$x,
        y = df_density$y
      ),
      color     = dots_lines$color,
      linewidth = dots_lines$linewidth,
      linetype  = dots_lines$linetype
    )
  } else {
    if (plot_ci) {
      graphics::polygon(
        c(df_density$x, rev(df_density$x)),
        c(df_density$y_lCI, rev(df_density$y_uCI)),
        border = NA,
        col    = scales::alpha(dots_lines$col, dots_lines$alpha)
      )
    }
    graphics::lines(
      df_density$x, df_density$y,
      lwd = dots_lines$lwd,
      col = dots_lines$col,
      lty = dots_lines$lty
    )
  }

  # return
  if (plot_type == "base") {
    return(invisible())
  } else if (plot_type == "ggplot") {
    return(out)
  }
}


# ============================================================================ #
# Internal Helper Functions
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# .zcurve_fun.brma
# ---------------------------------------------------------------------------- #
#
# Core computation function for z-curve estimates and densities.
#
# This function operates in two modes:
# 1. EDR mode (z_threshold set): Computes Expected Discovery Rate
# 2. Density mode (z_sequence set): Computes density over z-values
#
# The extrapolate parameter controls bias adjustment:
# - extrapolate=TRUE:  Removes publication bias (PET/PEESE/weights) to estimate
#                      "true" power distribution
# - extrapolate=FALSE: Includes all bias adjustments for fitted distribution
#
# @param object      brma object with fit and priors
# @param z_threshold z-value threshold for significance (EDR mode)
# @param z_sequence  vector of z-values for density evaluation (density mode)
# @param max_samples maximum posterior samples for computation
# @param extrapolate whether to remove bias adjustments
#
# @return EDR mode:     list with EDR (S-vector) and weights (S-vector)
#         Density mode: S x length(z_sequence) matrix of densities
#
# ---------------------------------------------------------------------------- #
.zcurve_fun.brma <- function(object, z_threshold = NULL, z_sequence = NULL,
                             max_samples = 1000, extrapolate = FALSE) {

  ### extract model info
  is_multilevel     <- .is_multilevel(object)
  is_weightfunction <- .is_weightfunction(object)
  outcome_type      <- .outcome_type(object)
  effect_direction  <- .effect_direction(object)
  priors            <- object[["priors"]]
  data              <- object[["data"]]

  ### get outcome data
  outcome_data <- object[["data"]][["outcome"]]
  yi           <- outcome_data[["yi"]]
  sei          <- outcome_data[["sei"]]
  K            <- length(yi)

  ### determine mu type
  # "cluster" for multilevel (mu + gamma), "terms" for marginal/single-level
  mu_type <- if (is_multilevel) "cluster" else "terms"

  ### 1. Predict mu (location) samples
  # predict.brma handles PET/PEESE via bias_adjusted argument
  # extrapolate=TRUE -> bias_adjusted=TRUE (unbiased predictions)
  mu_samples <- predict.brma(
    object        = object,
    type          = mu_type,
    bias_adjusted = extrapolate,
    quiet         = TRUE
  )

  ### 2. Get tau (heterogeneity) samples
  # specific handling to get tau_within
  # .evaluate.brma.tau handles ML splitting logic
  tau_result <- .evaluate.brma.tau(
    fit           = object[["fit"]],
    scale_data    = object[["data"]][["scale"]],
    scale_formula = if (.is_scale(object)) .create_fit_formula_list(data = object[["data"]], "scale") else NULL,
    scale_priors  = object[["priors"]][["scale"]],
    is_scale      = .is_scale(object),
    is_multilevel = is_multilevel,
    K             = K
  )
  tau_within <- tau_result[["tau_within"]]

  ### 3. Subsample for efficiency
  S <- nrow(mu_samples)
  if (S > max_samples) {
    selected_ind <- round(seq(from = 1, to = S, length.out = max_samples))
    mu_samples   <- mu_samples[selected_ind, , drop = FALSE]
    tau_within   <- tau_within[selected_ind, , drop = FALSE]
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))[selected_ind, ]
    S <- max_samples
  } else {
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
  }

  ### 4. Prepare for Computation
  selection <- .zcurve_selection_context(
    object            = object,
    data              = data,
    priors            = priors,
    posterior_samples = posterior_samples,
    is_weightfunction = is_weightfunction
  )


  ### 5. Dispatch: EDR (Threshold) vs Density (Sequence)

  if (!is.null(z_threshold)) {
    return(.zcurve_threshold_vectorized(
      z_threshold      = z_threshold,
      mu_samples       = mu_samples,
      tau_within       = tau_within,
      sei              = sei,
      selection        = selection,
      extrapolate      = extrapolate,
      effect_direction = effect_direction
    ))
  }

  if (!is.null(z_sequence)) {
    return(.zcurve_density_vectorized(
      z_sequence       = z_sequence,
      mu_samples       = mu_samples,
      tau_within       = tau_within,
      sei              = sei,
      selection        = selection,
      extrapolate      = extrapolate,
      effect_direction = effect_direction
    ))
  }

  return(NULL)
}


# ---------------------------------------------------------------------------- #
# .zcurve_threshold_vectorized
# ---------------------------------------------------------------------------- #
#
# Vectorized EDR computation for normal and selected-normal rows.
#
# ---------------------------------------------------------------------------- #
.zcurve_threshold_vectorized <- function(z_threshold, mu_samples, tau_within,
                                         sei, selection, extrapolate,
                                         effect_direction) {

  S         <- nrow(mu_samples)
  K         <- ncol(mu_samples)
  sei_mat   <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd  <- sqrt(tau_within^2 + sei_mat^2)
  q_upper   <- z_threshold * sei
  q_lower   <- -z_threshold * sei

  thresholds <- stats::pnorm(
    matrix(q_upper, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu_samples,
    sd         = total_sd,
    lower.tail = FALSE
  ) + stats::pnorm(
    matrix(q_lower, nrow = S, ncol = K, byrow = TRUE),
    mean       = mu_samples,
    sd         = total_sd,
    lower.tail = TRUE
  )
  weights       <- .zcurve_inverse_selection_weights(
    mean      = mu_samples,
    sd        = total_sd,
    sei       = sei,
    selection = selection
  )
  weighted_rows <- if (is.null(selection)) integer(0) else which(!selection[["use_normal"]])
  if (length(weighted_rows) > 0) {
    selection_weight <- .selection_context_subset_rows(selection, weighted_rows)
    mean_weight      <- mu_samples[weighted_rows, , drop = FALSE]
    sd_weight        <- total_sd[weighted_rows, , drop = FALSE]

    if (!extrapolate) {
      prob_upper <- .selection_step_cdf_matrix(
        q                 = q_upper,
        mean              = mean_weight,
        sd                = sd_weight,
        sei               = sei,
        selection_context = selection_weight,
        lower.tail        = FALSE
      )
      prob_lower <- .selection_step_cdf_matrix(
        q                 = q_lower,
        mean              = mean_weight,
        sd                = sd_weight,
        sei               = sei,
        selection_context = selection_weight,
        lower.tail        = TRUE
      )

      thresholds[weighted_rows, ] <- prob_upper + prob_lower
    }
  }

  if (!is.null(selection) && extrapolate) {
    EDR <- rowSums(thresholds * weights) / rowSums(weights)
  } else {
    EDR <- rowMeans(thresholds)
  }

  return(list(
    EDR     = EDR,
    weights = rowMeans(weights)
  ))
}


# ---------------------------------------------------------------------------- #
# .zcurve_density_vectorized
# ---------------------------------------------------------------------------- #
#
# Vectorized z-density computation for normal and selected-normal rows.
#
# ---------------------------------------------------------------------------- #
.has_native_zcurve_density <- function(selection = FALSE) {

  symbols <- "RoBMA_zcurve_normal_density_matrix"
  if (selection) {
    symbols <- c(symbols, "RoBMA_selnorm_zcurve_density_matrix")
  }

  return(all(vapply(symbols, is.loaded, logical(1), PACKAGE = "RoBMA")))
}

.zcurve_normal_density_matrix <- function(z_sequence, mean, sd, sei) {

  if (!.has_native_zcurve_density(selection = FALSE)) {
    stop("The native z-curve density kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_zcurve_normal_density_matrix",
    .native_numeric_vector(z_sequence),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    PACKAGE = "RoBMA"
  ))
}

.zcurve_selnorm_density_matrix <- function(z_sequence, mean, sd, sei,
                                           selection_context, extrapolate) {

  .selection_require_step_evaluable(selection_context, ".zcurve_density_vectorized()")

  if (!.has_native_zcurve_density(selection = TRUE)) {
    stop("The native selected-normal z-curve density kernel is not loaded.", call. = FALSE)
  }

  return(.Call(
    "RoBMA_selnorm_zcurve_density_matrix",
    .native_numeric_vector(z_sequence),
    .native_numeric_matrix(mean),
    .native_numeric_matrix(sd),
    .native_numeric_vector(sei),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["alpha"]]),
    .native_integer_vector(selection_context[["phack_kind"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    .native_numeric_vector(selection_context[["z_lower"]]),
    .native_numeric_vector(selection_context[["z_upper"]]),
    .native_integer_vector(selection_context[["sign"]]),
    .native_integer_vector(selection_context[["phack_q"]]),
    .native_numeric_vector(selection_context[["phack_z_source"]]),
    .native_numeric_vector(selection_context[["phack_z_dest"]]),
    .native_numeric_vector(selection_context[["segments"]][["bounds"]]),
    .native_integer_vector(selection_context[["segments"]][["step_bin"]]),
    .native_integer_vector(selection_context[["segments"]][["phack_region"]]),
    as.logical(extrapolate),
    PACKAGE = "RoBMA"
  ))
}

.zcurve_density_vectorized <- function(z_sequence, mu_samples, tau_within,
                                       sei, selection, extrapolate,
                                       effect_direction) {

  S           <- nrow(mu_samples)
  K           <- ncol(mu_samples)
  sei_mat     <- matrix(sei, nrow = S, ncol = K, byrow = TRUE)
  total_sd    <- sqrt(tau_within^2 + sei_mat^2)

  if (is.null(selection)) {
    return(.zcurve_normal_density_matrix(
      z_sequence = z_sequence,
      mean       = mu_samples,
      sd         = total_sd,
      sei        = sei
    ))
  }

  density <- .zcurve_selnorm_density_matrix(
    z_sequence        = z_sequence,
    mean              = mu_samples,
    sd                = total_sd,
    sei               = sei,
    selection_context = selection,
    extrapolate       = extrapolate
  )

  return(density)
}


# ---------------------------------------------------------------------------- #
# .zcurve_inverse_selection_weights
# ---------------------------------------------------------------------------- #
#
# Expected attempted studies represented by each observed study under the
# selection model. Normal/no-bias branches have weight one.
#
# ---------------------------------------------------------------------------- #
.zcurve_inverse_selection_weights <- function(mean, sd, sei, selection) {

  weights <- matrix(1, nrow = nrow(mean), ncol = ncol(mean))

  if (is.null(selection)) {
    return(weights)
  }

  weighted_rows <- which(!selection[["use_normal"]])
  if (length(weighted_rows) == 0L) {
    return(weights)
  }

  selection_weight <- .selection_context_subset_rows(selection, weighted_rows)
  log_norm <- .selection_step_log_norm_matrix(
    mean              = mean[weighted_rows, , drop = FALSE],
    sd                = sd[weighted_rows, , drop = FALSE],
    sei               = sei,
    selection_context = selection_weight
  )

  weights[weighted_rows, ] <- exp(-log_norm)
  return(weights)
}


# ---------------------------------------------------------------------------- #
# .zcurve_weighted_rows / .zcurve_constant_rows
# ---------------------------------------------------------------------------- #
#
# Row selectors for fitted and extrapolated selection-model z-curves.
#
# ---------------------------------------------------------------------------- #
.zcurve_weighted_rows <- function(selection, extrapolate) {

  if (is.null(selection) || extrapolate) {
    return(integer(0))
  }

  return(which(!selection[["use_normal"]]))
}

.zcurve_constant_rows <- function(selection, extrapolate) {

  return(integer(0))
}


# ---------------------------------------------------------------------------- #
# .zcurve_selection_context
# ---------------------------------------------------------------------------- #
#
# Prepare posterior-row selection metadata for branch-aware z-curve evaluation.
#
# ---------------------------------------------------------------------------- #
.zcurve_selection_context <- function(object, data, priors, posterior_samples,
                                      is_weightfunction) {

  if (!is_weightfunction) {
    return(NULL)
  }

  return(.selection_context(
    object            = object,
    posterior_samples = posterior_samples
  ))
}


# ---------------------------------------------------------------------------- #
# .zcurve_selection_args
# ---------------------------------------------------------------------------- #
#
# Extract active omega and cutpoints for one posterior row and estimate.
#
# ---------------------------------------------------------------------------- #
.zcurve_selection_args <- function(selection, row, estimate, n = 1L) {

  .selection_require_step_evaluable(selection, ".zcurve_selection_args()")

  omega <- selection[["omega"]][row, , drop = FALSE]

  if (n > 1L) {
    omega <- matrix(as.numeric(omega), nrow = n, ncol = ncol(omega), byrow = TRUE)
  }

  return(list(
    omega   = omega,
    crit_yi = stats::qnorm(rev(selection[["p_cuts"]])[-c(1L, length(selection[["p_cuts"]]))],
                           lower.tail = FALSE) *
      selection[["sei"]][estimate]
  ))
}


# ============================================================================ #
# Graphical Helper Functions
# ============================================================================ #
#
# These functions extract and set default graphical parameters for z-curve
# plotting components. They handle the translation between base graphics
# and ggplot2 parameter naming conventions.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .get_dots_hist_zcurve
# ---------------------------------------------------------------------------- #
#
# Extracts histogram graphical parameters with defaults.
#
# @param dots      list of user-supplied graphical parameters
# @param plot_type "base" or "ggplot"
# @param max_density maximum density value for setting ylim
#
# @return list of graphical parameters appropriate for plot_type
#
# ---------------------------------------------------------------------------- #
.get_dots_hist_zcurve <- function(dots, plot_type = "base", max_density = NULL) {
  if (plot_type == "base") {
    if (!is.null(dots[["ylim"]]) && !is.null(max_density)) {
      ylim <- range(dots[["ylim"]], max_density)
    } else if (!is.null(max_density)) {
      ylim <- c(0, max_density)
    } else {
      ylim <- NULL
    }
    hist_dots <- list(
      border = if (!is.null(dots[["border"]])) dots[["border"]] else "gray60",
      col    = if (!is.null(dots[["col"]]))    dots[["col"]]    else NA,
      xlab   = if (!is.null(dots[["xlab"]]))   dots[["xlab"]]   else "Z-Statistic",
      ylab   = if (!is.null(dots[["ylab"]]))   dots[["ylab"]]   else "Density",
      main   = if (!is.null(dots[["main"]]))   dots[["main"]]   else "",
      ylim   = ylim,
      xaxt   = if (!is.null(dots[["xaxt"]]))   dots[["xaxt"]]   else "s",
      yaxt   = if (!is.null(dots[["yaxt"]]))   dots[["yaxt"]]   else "s",
      las    = if (!is.null(dots[["las"]]))    dots[["las"]]    else 1
    )
  } else if (plot_type == "ggplot") {
    hist_dots <- list(
      color = if (!is.null(dots[["color"]])) dots[["color"]] else if (!is.null(dots[["col"]])) dots[["col"]] else "gray60",
      fill  = if (!is.null(dots[["fill"]]))  dots[["fill"]]  else NA,
      alpha = if (!is.null(dots[["alpha"]])) dots[["alpha"]] else 1,
      xlab  = if (!is.null(dots[["xlab"]]))  dots[["xlab"]]  else "Z-Statistic",
      ylab  = if (!is.null(dots[["ylab"]]))  dots[["ylab"]]  else "Density",
      main  = if (!is.null(dots[["main"]]))  dots[["main"]]  else ""
    )
  }
  return(hist_dots)
}


# ---------------------------------------------------------------------------- #
# .get_dots_lines_zcurve
# ---------------------------------------------------------------------------- #
#
# Extracts density line graphical parameters with defaults.
#
# @param dots      list of user-supplied graphical parameters
# @param plot_type "base" or "ggplot"
# @param col       default line color
# @param alpha     default alpha for CI ribbons
#
# @return list of graphical parameters appropriate for plot_type
#
# ---------------------------------------------------------------------------- #
.get_dots_lines_zcurve <- function(dots, plot_type = "base", col = "black", alpha = 0.40) {
  if (plot_type == "base") {
    line_dots <- list(
      lwd   = if (!is.null(dots[["lwd"]]))   dots[["lwd"]]   else 2,
      col   = if (!is.null(dots[["col"]]))   dots[["col"]]   else col,
      lty   = if (!is.null(dots[["lty"]]))   dots[["lty"]]   else 1,
      alpha = if (!is.null(dots[["alpha"]])) dots[["alpha"]] else alpha
    )
  } else if (plot_type == "ggplot") {
    line_dots <- list(
      linewidth = if (!is.null(dots[["linewidth"]])) dots[["linewidth"]] else if (!is.null(dots[["lwd"]])) dots[["lwd"]] else 1,
      color     = if (!is.null(dots[["color"]]))     dots[["color"]]     else if (!is.null(dots[["col"]])) dots[["col"]] else col,
      linetype  = if (!is.null(dots[["linetype"]])) dots[["linetype"]]  else if (!is.null(dots[["lty"]])) dots[["lty"]] else 1,
      alpha     = if (!is.null(dots[["alpha"]]))    dots[["alpha"]]     else alpha
    )
  }
  return(line_dots)
}


# ---------------------------------------------------------------------------- #
# .get_dots_thresholds_zcurve
# ---------------------------------------------------------------------------- #
#
# Extracts threshold line graphical parameters with defaults.
#
# @param dots      list of user-supplied graphical parameters
# @param plot_type "base" or "ggplot"
#
# @return list of graphical parameters appropriate for plot_type
#
# ---------------------------------------------------------------------------- #
.get_dots_thresholds_zcurve <- function(dots, plot_type = "base") {
  if (plot_type == "base") {
    threshold_dots <- list(
      col = if (!is.null(dots[["col"]])) dots[["col"]] else "red",
      lty = if (!is.null(dots[["lty"]])) dots[["lty"]] else 3,
      lwd = if (!is.null(dots[["lwd"]])) dots[["lwd"]] else 1
    )
  } else if (plot_type == "ggplot") {
    threshold_dots <- list(
      color     = if (!is.null(dots[["color"]])) dots[["color"]] else if (!is.null(dots[["col"]])) dots[["col"]] else "red",
      linetype  = if (!is.null(dots[["linetype"]])) dots[["linetype"]] else if (!is.null(dots[["lty"]])) dots[["lty"]] else "dashed",
      linewidth = if (!is.null(dots[["linewidth"]])) dots[["linewidth"]] else if (!is.null(dots[["lwd"]])) dots[["lwd"]] else 0.5
    )
  }
  return(threshold_dots)
}


# ---------------------------------------------------------------------------- #
# .get_Soric_FDR
# ---------------------------------------------------------------------------- #
#
# Computes Soric's (1989) False Discovery Rate estimate from EDR.
#
# @param EDR       Expected Discovery Rate (vector or scalar)
# @param sig_level two-sided significance level (e.g., 0.05)
#
# @return Soric FDR estimate
#
# ---------------------------------------------------------------------------- #
.get_Soric_FDR <- function(EDR, sig_level) {

  ((1 / EDR) - 1) * (sig_level / (1 - sig_level))
}


# ---------------------------------------------------------------------------- #
# .zcurve_bins
# ---------------------------------------------------------------------------- #
#
# Creates bin sequence for z-curve histograms and densities.
#
# For histograms, bins are shifted to align with significance thresholds
# when selection model priors are present. For densities, threshold
# values are added to the sequence to ensure accurate representation.
#
# @param priors     list of priors (used to extract thresholds)
# @param from       lower bound of z-value range
# @param to         upper bound of z-value range
# @param by         step size (NULL to use length.out)
# @param length.out number of bins (NULL to use by)
# @param type       "hist" or "dens"
#
# @return numeric vector of bin boundaries or density evaluation points
#
# ---------------------------------------------------------------------------- #
.zcurve_bins      <- function(priors, from, to, by, length.out, type = "hist"){

  if(is.null(length.out)){
    bins <- seq(from = from, to = to, by = by)
  }else{
    bins <- seq(from = from, to = to, length.out = length.out)
  }

  # obtain tresholds from the specified priors
  z_bounds <- .zcurve_threshold(priors)

  # return simple binning in case of no weightfunctions
  if(length(z_bounds) == 0){
    return(bins)
  }

  # retain bounds within the plotting range
  z_bounds <- z_bounds[z_bounds > from & z_bounds < to]

  if(type == "hist"){
    # for histogram, shift the specified bin boundaries to the closest z-treshold
    for(i in seq_along(z_bounds)){

      # get index of the first larger sequence
      i_larger <- which(bins > z_bounds[i])[1]

      # if there is none skip
      if(is.na(i_larger))
        next

      # get index of the closer one from below
      i_lower  <- i_larger - 1

      # replace the closer sequence with the boundary
      if(bins[i_larger] - z_bounds[i] < z_bounds[i] - bins[i_lower]){
        bins[i_larger] <- z_bounds[i]
      }else{
        bins[i_lower]  <- z_bounds[i]
      }
    }
  }else if(type == "dens"){
    # for density, extend the specified support at the threshold
    bins <- sort(unique(c(bins, z_bounds)))
  }

  return(bins)
}
.zcurve_bias_priors <- function(priors){

  if(is.null(priors)){
    return(list())
  }

  if(!is.null(priors[["outcome"]])){
    return(.selection_bias_priors(priors))
  }

  priors_bias <- priors[["bias"]]
  if(is.null(priors_bias)){
    return(list())
  }
  if(BayesTools::is.prior.mixture(priors_bias)){
    out <- as.list(priors_bias)
    class(out) <- "list"
    return(out)
  }

  return(list(priors_bias))
}
.zcurve_threshold <- function(priors){

  priors_bias <- .zcurve_bias_priors(priors)

  # return simple binning in case of no bias
  if(length(priors_bias) == 0L){
    return()
  }

  p_bounds <- .selection_kernel_breaks(priors_bias)

  # return simple binning in case of no selection kernels
  if(is.null(p_bounds)){
    return()
  }

  # obtain tresholds from the specified priors
  z_bounds <- stats::qnorm(rev(p_bounds), 0, 1, lower.tail = FALSE)
  z_bounds <- z_bounds[!is.infinite(z_bounds)]

  return(z_bounds)
}
