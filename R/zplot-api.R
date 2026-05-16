#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# as_zplot generic and method
# ---------------------------------------------------------------------------- #

#' @export
as_zplot <- function(object, ...) UseMethod("as_zplot")


#' @title Transform brma Object to Zplot
#'
#' @description Transforms an estimated brma model into a zplot object
#' that can be summarized and plotted to assess replicability.
#'
#' @param object a normal-outcome \code{brma} object.
#' @param significance_level z-value threshold for significance. Defaults
#' to \code{qnorm(0.975)} (two-sided alpha = 0.05).
#' @param max_samples maximum number of posterior samples for estimation.
#' Defaults to 10000. Use \code{Inf} to use all posterior samples.
#' @param ... additional arguments (currently unused).
#'
#' @details
#' Zplot analysis estimates the Expected Discovery Rate (EDR), the posterior
#' mean probability that an exact replication would be statistically significant
#' at the supplied threshold. This provides insight into the replicability of
#' findings in a literature.
#'
#' The implementation uses extrapolation mode by default, which removes
#' selection-model weights when evaluating the model-implied z-value density.
#' PET/PEESE regression terms remain part of the fitted location model.
#' Zplot diagnostics are available only for normal outcome models. GLMM
#' objects are rejected because their raw likelihood is on a count scale while
#' zplot diagnostics require observed effect-size z-statistics with standard
#' errors.
#'
#' The resulting object retains all original brma properties while adding
#' zplot results, enabling both standard meta-analytic summaries and
#' zplot diagnostics on the same object.
#'
#' @return The input object with added class \code{"zplot_brma"} and a new
#' \code{zplot} list component containing:
#' \describe{
#'   \item{estimates}{a list with posterior samples for \code{EDR} and
#'     missing-study \code{weights}}
#'   \item{data}{a list with \code{significance_level}, observed
#'     \code{z}-statistics, \code{N_significant}, and \code{N_observed}}
#' }
#'
#' @seealso [summary.zplot_brma()], [plot.zplot_brma()], [hist.zplot_brma()]
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   zfit <- as_zplot(fit)
#'   summary(zfit)
#'   plot(zfit)
#' }
#' }
#'
#' @aliases as_zplot
#' @export
as_zplot.brma <- function(object, significance_level = stats::qnorm(0.975), max_samples = 10000, ...) {

  BayesTools::check_real(significance_level, "significance_level", lower = 0)
  max_samples <- .normalize_max_samples(max_samples, "max_samples")

  if (.outcome_type(object) != "norm") {
    stop(
      "as_zplot() is only available for normal outcome models; ",
      "GLMM objects are not supported.",
      call. = FALSE
    )
  }

  # compute the main estimates (extrapolate = TRUE: unbiased power)
  # This corresponds to "Expected Discovery Rate" if studies were exact replications
  # but without publication bias
  out_estimates <- .zplot_fun.brma(
    object      = object,
    z_threshold = significance_level,
    max_samples = max_samples,
    extrapolate = TRUE
  )

  # compute data summaries (observed z-stats)
  yi  <- .outcome_data_yi(object)
  sei <- .outcome_data_sei(object)
  z   <- yi / sei

  # store zplot results
  object[["zplot"]] <- list(
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
  class(object) <- c("zplot_brma", class(object))

  return(object)
}


# ---------------------------------------------------------------------------- #
# zplot generic and methods
# ---------------------------------------------------------------------------- #

#' @export
zplot <- function(object, ...) UseMethod("zplot")


#' @title Plot Zplot Diagnostics Directly
#'
#' @description Convenience wrapper for creating and plotting zplot diagnostics
#' from a fitted \code{brma} object.
#'
#' @param object a normal-outcome \code{brma} object, or a
#' \code{zplot_brma} object.
#' @param significance_level z-value threshold for significance. Defaults
#' to \code{qnorm(0.975)} (two-sided alpha = 0.05).
#' @param summary_max_samples maximum number of posterior samples used for the
#' EDR and missing-study summaries stored in the generated zplot object.
#' This is separate from the plot-density \code{max_samples} argument accepted
#' by \code{plot.zplot_brma()} through \code{...}. Defaults to 10000. Use
#' \code{Inf} to use all posterior samples.
#' @param ... arguments passed to \code{\link[=plot.zplot_brma]{plot.zplot_brma()}}.
#'
#' @details When \code{object} already inherits from \code{zplot_brma},
#' \code{zplot()} dispatches directly to \code{plot.zplot_brma()} without
#' recomputing stored summaries.
#'
#' @return \code{NULL} invisibly for base graphics, or a ggplot2 object.
#'
#' @seealso [as_zplot.brma()], [plot.zplot_brma()], [summary.zplot_brma()]
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   zplot(fit)
#' }
#' }
#'
#' @aliases zplot
#' @export
zplot.brma <- function(object, significance_level = stats::qnorm(0.975),
                       summary_max_samples = 10000, ...) {

  zplot_object <- as_zplot(
    object             = object,
    significance_level = significance_level,
    max_samples        = summary_max_samples
  )

  return(plot(zplot_object, ...))
}


#' @rdname zplot.brma
#' @export
zplot.zplot_brma <- function(object, ...) {

  return(plot(object, ...))
}


# ---------------------------------------------------------------------------- #
# summary.zplot_brma
# ---------------------------------------------------------------------------- #

#' @title Summarize Zplot Results
#'
#' @description Creates summary tables for zplot estimates including
#' EDR, Soric FDR, and estimated missing studies.
#'
#' @param object a zplot_brma object.
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
#' @return An object of class \code{"summary.zplot_brma"} containing the
#' estimates table.
#'
#' @seealso [as_zplot.brma()], [plot.zplot_brma()]
#'
#' @export
summary.zplot_brma <- function(object, probs = c(.025, .975), ...) {

  # get proportion of significant results
  obs_proportion <- stats::prop.test(
    object$zplot$data[["N_significant"]],
    object$zplot$data[["N_observed"]],
    conf.level = 0.95
  )

  info_text <- sprintf(
    "Estimated using %1$i estimates, %2$i significant (ODR = %3$.2f, 95%% CI [%4$.2f, %5$.2f]).",
    object$zplot$data[["N_observed"]],
    object$zplot$data[["N_significant"]],
    obs_proportion$estimate,
    obs_proportion$conf.int[1],
    obs_proportion$conf.int[2]
  )

  # compute estimates
  sig_level <- stats::pnorm(object$zplot$data[["significance_level"]], lower.tail = FALSE) * 2

  estimates <- cbind.data.frame(
    "EDR"       = object$zplot$estimates[["EDR"]],
    "Soric FDR" = .get_Soric_FDR(object$zplot$estimates[["EDR"]], sig_level),
    "Missing N" = (object$zplot$estimates[["weights"]] - 1) * object$zplot$data[["N_observed"]]
  )

  estimates_table <- BayesTools::ensemble_estimates_table(
    samples    = estimates,
    parameters = names(estimates),
    probs      = probs,
    title      = "Zplot Estimates:",
    footnotes  = info_text
  )

  output <- list(
    estimates = estimates_table
  )

  class(output) <- "summary.zplot_brma"
  return(output)
}


#' @title Print Zplot Summary
#'
#' @param x a summary.zplot_brma object.
#' @param ... additional arguments (currently unused).
#'
#' @return Invisibly returns the summary object.
#'
#' @export
print.summary.zplot_brma <- function(x, ...) {

  cat("\n")
  print(x[["estimates"]])
  cat("\n")

  return(invisible(x))
}


#' @export
print.zplot_brma <- function(x, ...) {
  print(summary(x, ...))
}


# ---------------------------------------------------------------------------- #
# plot.zplot_brma
# ---------------------------------------------------------------------------- #

#' @title Plot Zplot Results
#'
#' @description Plots a zplot visualization showing the histogram of
#' observed z-statistics overlaid with model-implied densities.
#'
#' @param x a zplot_brma object.
#' @param plot_type graphics system: \code{"base"} or \code{"ggplot"}.
#' Defaults to \code{"base"}.
#' @param probs quantiles for credible intervals. Defaults to \code{c(.025, .975)}.
#' @param max_samples maximum posterior samples for density estimation.
#' Defaults to 10000. Use \code{Inf} to use all posterior samples.
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
#' @seealso [hist.zplot_brma()], [lines.zplot_brma()], [summary.zplot_brma()]
#'
#' @export
plot.zplot_brma <- function(x, plot_type = "base",
                             probs = c(.025, .975), max_samples = 10000,
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
  max_samples <- .normalize_max_samples(max_samples, "max_samples")

  dots <- list(...)

  # get line values first so we can set ylim
  ymax <- 0

  # 1. Fit lines (Fitted model - with bias)
  if (plot_fit) {
    dots_fit <- if(!is.null(dots_fit)) dots_fit else list()
    lines_fit <- do.call(lines.zplot_brma, c(
      list(x = x, plot_type = plot_type, probs = probs, max_samples = max_samples,
           extrapolate = FALSE, plot_ci = plot_ci, from = from, to = to,
           by = by.lines, length.out = length.out.lines, as_data = TRUE),
      dots_fit, dots
    ))
    ymax <- max(c(ymax, lines_fit$y, if (plot_ci) lines_fit$y_uCI))
  }

  # 2. Extrapolation lines (unbiased zplot density)
  if (plot_extrapolation) {
    # prepare extrapolation dots with default blue color
    dots_extrapolation <- if(!is.null(dots_extrapolation)) dots_extrapolation else list()
    lines_extrapolation <- do.call(lines.zplot_brma, c(
      list(x = x, plot_type = plot_type, probs = probs, max_samples = max_samples,
           extrapolate = TRUE, plot_ci = plot_ci, from = from, to = to,
           by = by.lines, length.out = length.out.lines, as_data = TRUE,
           col = "blue"), # default color if not in dots
      dots_extrapolation, dots
    ))
    ymax <- max(c(ymax, lines_extrapolation$y, if (plot_ci) lines_extrapolation$y_uCI))
  }

  # 3. Create histogram base
  plot <- hist.zplot_brma(
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
    dots_fit_lines <- .get_dots_lines_zplot(c(dots, dots_fit), plot_type = plot_type)

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
     dots_ext_lines <- .get_dots_lines_zplot(c(dots, dots_extrapolation), plot_type = plot_type, col = "blue")

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
# hist.zplot_brma
# ---------------------------------------------------------------------------- #

#' @title Histogram of Z-Statistics
#'
#' @description Plots a histogram of observed z-values from the meta-analysis.
#'
#' @param x a zplot_brma object.
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
#' @seealso [plot.zplot_brma()], [lines.zplot_brma()]
#'
#' @export
hist.zplot_brma <- function(x, plot_type = "base",
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

  z <- x$zplot$data[["z"]]

  # create bin sequence using zplot helper (handles thresholds)
  z_sequence <- .zplot_bins(priors = x[["priors"]], from = from, to = to, by = by, length.out = length.out, type = "hist")
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
  dots_hist_params <- .get_dots_hist_zplot(dots, plot_type = plot_type, max_density = max(z_hist$density))

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
    thresholds <- .zplot_threshold(x[["priors"]])
    if (length(thresholds) > 0) {
      dots_thresholds <- if (!is.null(dots_thresholds)) dots_thresholds else list()
      dots_thresholds_params <- .get_dots_thresholds_zplot(dots_thresholds, plot_type = plot_type)

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
# lines.zplot_brma
# ---------------------------------------------------------------------------- #

#' @title Add Zplot Density Lines
#'
#' @description Adds model-implied density lines to an existing zplot.
#'
#' @param x a zplot_brma object.
#' @param plot_type graphics system: \code{"base"} or \code{"ggplot"}.
#' Defaults to \code{"base"}.
#' @param probs quantiles for credible intervals. Defaults to \code{c(.025, .975)}.
#' @param max_samples maximum posterior samples for density. Defaults to 10000.
#' Use \code{Inf} to use all posterior samples.
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
#' @seealso [plot.zplot_brma()], [hist.zplot_brma()]
#'
#' @export
lines.zplot_brma <- function(x, plot_type = "base",
                              probs = c(.025, .975), max_samples = 10000,
                              plot_ci = TRUE, extrapolate = FALSE,
                              from = -6, to = 6, by = 0.05, length.out = NULL,
                              col = "black", as_data = FALSE, ...) {

  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_real(probs, "probs", lower = 0, upper = 1, check_length = 2)
  max_samples <- .normalize_max_samples(max_samples, "max_samples")
  BayesTools::check_bool(plot_ci, "plot_ci")
  BayesTools::check_bool(extrapolate, "extrapolate")
  BayesTools::check_real(from, "from")
  BayesTools::check_real(to, "to")
  BayesTools::check_real(by, "by", allow_NULL = TRUE)


  dots <- list(...)

  # prepare z sequence
  z_sequence <- .zplot_bins(priors = x[["priors"]], from = from, to = to, by = by, length.out = length.out, type = "dens")

  # get densities from model
  z_density <- .zplot_fun.brma(
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
  dots_lines <- .get_dots_lines_zplot(dots, plot_type = plot_type, col = col)

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
