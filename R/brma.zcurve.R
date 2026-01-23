# ============================================================================ #
# brma.zcurve.R
# ============================================================================ #
#
# This file implements Z-curve diagnostics for brma class objects.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# as_zcurve generic and method
# ---------------------------------------------------------------------------- #

#' @export
as_zcurve <- function(object, ...) UseMethod("as_zcurve")


#' @title Transform brma object into z-curve
#'
#' @description
#' This function transforms the estimated brma model into a
#' z-curve object that can be further summarized and plotted.
#'
#' @param object A brma object
#' @param significance_level Significance level used for computation of z-curve
#' estimates. Defaults to \code{qnorm(0.975)} (two-sided alpha = 0.05).
#' @param max_samples Maximum number of samples from the posterior distribution
#' that will be used for estimating z-curve estimates. Defaults to 1000.
#' @param ... additional arguments (currently unused).
#'
#' @return \code{as_zcurve} returns the input object with class "zcurve_brma"
#' and a new \code{"zcurve"} list component containing the results.
#'
#' @seealso \code{\link{summary.zcurve_brma}}, \code{\link{plot.zcurve_brma}}
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
#' @description Creates summary tables for a zcurve_brma object.
#'
#' @param object a zcurve_brma object.
#' @param probs quantiles of the posterior distribution to display.
#' Defaults to \code{c(.025, .975)}.
#' @param ... additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.zcurve_brma"}.
#'
#' @seealso \code{\link{as_zcurve.brma}}, \code{\link{plot.zcurve_brma}}
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
#' @description Plots a z-curve visualization for a zcurve_brma object,
#' showing the histogram of z-statistics with the model-implied density.
#'
#' @param x a zcurve_brma object.
#' @param plot_type whether to use base graphics (\code{"base"}) or
#' ggplot2 (\code{"ggplot"}). Defaults to \code{"base"}.
#' @param probs quantiles for credible intervals. Defaults to \code{c(.025, .975)}.
#' @param max_samples maximum number of posterior samples for plotting.
#' Defaults to 500.
#' @param plot_fit whether to include the model fit (biases included). Defaults to \code{TRUE}.
#' @param plot_extrapolation whether to include the extrapolated model fit (no bias). Defaults to \code{TRUE}.
#' @param plot_CI whether to include credible intervals. Defaults to \code{TRUE}.
#' @param plot_thresholds whether to display significance thresholds.
#' Defaults to \code{TRUE}.
#' @param from lower bound of z-value range for plotting. Defaults to -6.
#' @param to upper bound of z-value range for plotting. Defaults to 6.
#' @param by.hist bin width for the histogram. Defaults to 0.5.
#' @param length.out.hist number of bins for histogram (alternative to by.hist).
#' @param by.lines step size for plotting density lines. Defaults to 0.05.
#' @param length.out.lines number of points for density lines (alternative to by.lines).
#' @param dots_hist list of graphical parameters for the histogram.
#' @param dots_fit list of graphical parameters for fit lines.
#' @param dots_extrapolation list of graphical parameters for extrapolation lines.
#' @param dots_thresholds list of graphical parameters for threshold lines.
#' @param ... additional graphical parameters.
#'
#' @return Returns \code{NULL} invisibly for base graphics, or a ggplot2
#' object if \code{plot_type = "ggplot"}.
#'
#' @seealso \code{\link{hist.zcurve_brma}}, \code{\link{lines.zcurve_brma}}
#'
#' @export
plot.zcurve_brma <- function(x, plot_type = "base",
                             probs = c(.025, .975), max_samples = 500,
                             plot_fit = TRUE, plot_extrapolation = TRUE,
                             plot_CI = TRUE, plot_thresholds = TRUE,
                             from = -6, to = 6,
                             by.hist = 0.5, length.out.hist = NULL,
                             by.lines = 0.05, length.out.lines = NULL,
                             dots_hist = NULL, dots_fit = NULL,
                             dots_extrapolation = NULL, dots_thresholds = NULL, ...) {

  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(plot_fit, "plot_fit")
  BayesTools::check_bool(plot_extrapolation, "plot_extrapolation")
  BayesTools::check_bool(plot_CI, "plot_CI")
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
           extrapolate = FALSE, plot_CI = plot_CI, from = from, to = to,
           by = by.lines, length.out = length.out.lines, as_data = TRUE),
      dots_fit, dots
    ))
    ymax <- max(c(ymax, lines_fit$y, if (plot_CI) lines_fit$y_uCI))
  }

  # 2. Extrapolation lines (Unbiased model - "True" z-curve)
  if (plot_extrapolation) {
    # prepare extrapolation dots with default blue color
    dots_extrapolation <- if(!is.null(dots_extrapolation)) dots_extrapolation else list()
    lines_extrapolation <- do.call(lines.zcurve_brma, c(
      list(x = x, plot_type = plot_type, probs = probs, max_samples = max_samples,
           extrapolate = TRUE, plot_CI = plot_CI, from = from, to = to,
           by = by.lines, length.out = length.out.lines, as_data = TRUE,
           col = "blue"), # default color if not in dots
      dots_extrapolation, dots
    ))
    ymax <- max(c(ymax, lines_extrapolation$y, if (plot_CI) lines_extrapolation$y_uCI))
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
      if (plot_CI) {
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
      if (plot_CI) {
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
      if (plot_CI) {
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
      if (plot_CI) {
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
#' @description Plots a histogram of observed z-values for a zcurve_brma object.
#'
#' @param x a zcurve_brma object.
#' @param plot_type whether to use base graphics (\code{"base"}) or
#' ggplot2 (\code{"ggplot"}). Defaults to \code{"base"}.
#' @param from lower bound of z-value range. Defaults to -6.
#' @param to upper bound of z-value range. Defaults to 6.
#' @param by bin width. Defaults to 0.5.
#' @param length.out number of bins.
#' @param add logical; if TRUE, adds to existing plot. Defaults to FALSE.
#' @param plot_thresholds whether to display significance thresholds.
#' Defaults to TRUE.
#' @param dots_thresholds list of graphical parameters for threshold lines.
#' @param dots_hist list of graphical parameters for histogram.
#' @param dots_all list of all graphical parameters.
#' @param ... additional graphical parameters.
#'
#' @return Returns \code{NULL} invisibly for base graphics, or a ggplot2
#' object if \code{plot_type = "ggplot"}.
#'
#' @seealso \code{\link{as_zcurve.brma}}, \code{\link{plot.zcurve_brma}}
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
        las    = 1,
        border = dots_hist_params$border,
        col    = dots_hist_params$col,
        xlab   = dots_hist_params$xlab,
        ylab   = dots_hist_params$ylab,
        main   = dots_hist_params$main,
        ylim   = dots_hist_params$ylim,
        xaxt   = dots_hist_params$xaxt,
        yaxt   = dots_hist_params$yaxt
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
#' @description Adds posterior predictive density lines to a z-curve plot.
#'
#' @param x a zcurve_brma object.
#' @param plot_type whether to use base graphics (\code{"base"}) or
#' ggplot2 (\code{"ggplot"}). Defaults to \code{"base"}.
#' @param probs quantiles for credible intervals. Defaults to \code{c(.025, .975)}.
#' @param max_samples maximum number of posterior samples. Defaults to 500.
#' @param plot_CI whether to include credible intervals. Defaults to TRUE.
#' @param extrapolate whether to extrapolate (unbiased) or show fit (biased). Defaults to FALSE.
#' @param from lower bound of z-value range. Defaults to -6.
#' @param to upper bound of z-value range. Defaults to 6.
#' @param by step size for density evaluation. Defaults to 0.05.
#' @param length.out number of points.
#' @param col line color. Defaults to "black".
#' @param as_data if TRUE, returns data frame instead of plotting.
#' @param ... additional graphical parameters.
#'
#' @return Returns \code{NULL} invisibly for base graphics, a list of ggplot2
#' layers if \code{plot_type = "ggplot"}, or a data frame if \code{as_data = TRUE}.
#'
#' @seealso \code{\link{as_zcurve.brma}}, \code{\link{plot.zcurve_brma}}
#'
#' @export
lines.zcurve_brma <- function(x, plot_type = "base",
                              probs = c(.025, .975), max_samples = 500,
                              plot_CI = TRUE, extrapolate = FALSE,
                              from = -6, to = 6, by = 0.05, length.out = NULL,
                              col = "black", as_data = FALSE, ...) {

  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_real(probs, "probs", lower = 0, upper = 1, check_length = 2)
  BayesTools::check_int(max_samples, "max_samples", lower = 10)
  BayesTools::check_bool(plot_CI, "plot_CI")
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
    if (plot_CI) {
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
    if (plot_CI) {
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
# @param object      brma object
# @param z_threshold significance threshold (z-value) for EDR computation
# @param z_sequence  vector of z-values for density computation
# @param max_samples maximum number of posterior samples to use
# @param extrapolate logical; if TRUE, computes unbiased estimates (removes PET/PEESE/Weights)
#                    if FALSE, computes fitted estimates (includes biases)
#
# @return if z_threshold is set: list(EDR, weights)
#         if z_sequence is set: S x length(z_sequence) matrix of densities
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
  # "study" for multilevel (mu + gamma), "terms" for marginal/single-level
  mu_type <- if (is_multilevel) "study" else "terms"

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
  # Need fit_data for selection models
  fit_data <- if (is_weightfunction) .create_fit_data(data = data, priors = priors) else NULL
  omega_samples <- if (is_weightfunction) posterior_samples[, grep("omega", colnames(posterior_samples)), drop = FALSE] else NULL


  ### 5. Dispatch: EDR (Threshold) vs Density (Sequence)

  if (!is.null(z_threshold)) {

    # --- EDR Computation ---

    outcome_thresholds <- rep(NA, S)
    outcome_weights    <- rep(NA, S)

    for (j in seq_len(S)) {
      temp_thresholds <- rep(NA, K)
      temp_weights    <- rep(NA, K)

      for (i in seq_len(K)) {

        # Total SD for prediction: sqrt(tau_within^2 + se^2)
        total_sd <- sqrt(tau_within[j, i]^2 + sei[i]^2)

        # Logic:
        # 1. Standard calculation using pnorm/dnorm
        # 2. If weightfunction AND extrapolate=FALSE: use weighted version (.dwnorm_fast)
        #    Else: use standard unweighted

        use_weighted <- is_weightfunction && !extrapolate

        if (!use_weighted) {
          # Standard Normal Probability > Threshold (two-tailed)
          # P(|Z| > z_crit) = P(Z > z_c) + P(Z < -z_c)
          prob <-
            stats::pnorm(z_threshold * sei[i], mu_samples[j, i], total_sd, lower.tail = FALSE) +
            stats::pnorm(-z_threshold * sei[i], mu_samples[j, i], total_sd, lower.tail = TRUE)

          temp_thresholds[i] <- prob
          temp_weights[i]    <- 1

        } else {
          # Weighted Normal Probability (Selection Models)
          # We need to account for the weight function which modifies the density

          # The density is f_w(y) = f(y) * w(y) / C
          # Probability in region is Integral(f_w(y))
          # .pwnorm_fast.ss computes cumulative probability for weighted normal

          # Adjust for negative effect direction (flip mu) if needed for wnorm helpers
          mu_j_i <- mu_samples[j, i]
          if (effect_direction == "negative") {
            mu_j_i <- -mu_j_i
          }

          # Calculate Prob(|Y| > z_crit * se)
          prob_upper <- .pwnorm_fast.ss(
            q          = z_threshold * sei[i],
            mean       = mu_j_i,
            sd         = total_sd,
            omega      = omega_samples[j, , drop = FALSE],
            crit_x     = fit_data$crit_yi[, i],
            lower.tail = FALSE
          )

          prob_lower <- .pwnorm_fast.ss(
            q          = -z_threshold * sei[i],
            mean       = mu_j_i,
            sd         = total_sd,
            omega      = omega_samples[j, , drop = FALSE],
            crit_x     = fit_data$crit_yi[, i],
            lower.tail = TRUE
          )

          temp_thresholds[i] <- prob_upper + prob_lower

          # Weight for normalization (if needed for output but EDR is just prob)
          # RoBMA returns 1/constant as weight?

          # Calculate normalizing constant C to return as weight
          # (used for "Missing N" calculation)
          temp_consts <- .dwnorm_fast.ss(
            x               = 0, # dummy
            mean            = mu_j_i,
            sd              = total_sd,
            omega           = omega_samples[j, , drop = FALSE],
            crit_x          = fit_data$crit_yi[, i],
            attach_constant = TRUE
          )
          temp_weights[i] <- 1 / attr(temp_consts, "constant")

        }
      }

      if (is_weightfunction && !extrapolate) {
         # for selection models, RoBMA weights the EDR computation by the inverse probability of publication?
         # "outcome_thresholds[j] <- stats::weighted.mean(temp_thresholds, temp_weights)" in zcurve.R line 1089
         outcome_thresholds[j] <- stats::weighted.mean(temp_thresholds, temp_weights)
      } else {
         mean(temp_thresholds)
      }

      outcome_weights[j] <- mean(temp_weights)
    }

    return(list(
      EDR     = outcome_thresholds,
      weights = outcome_weights
    ))


  } else if (!is.null(z_sequence)) {

    # --- Density Computation ---

    outcome_densities <- matrix(NA, nrow = S, ncol = length(z_sequence))

    for (j in seq_len(S)) {
      temp_densities_mat <- matrix(NA, nrow = K, ncol = length(z_sequence))

      for (i in seq_len(K)) {

        total_sd <- sqrt(tau_within[j, i]^2 + sei[i]^2)
        use_weighted <- is_weightfunction && !extrapolate

        if (!use_weighted) {

          # Unweighted density: dnorm(z*se, mu, sd) * se
          # Jacobian * se because we transform from y to z scale: z = y/se -> dy = se * dz
          temp_densities_mat[i, ] <- stats::dnorm(
            z_sequence * sei[i],
            mean = mu_samples[j, i],
            sd   = total_sd
          ) * sei[i]

          # Special case: Extrapolate with Selection Model
          # If we are extrapolating a selection model, we need to adjust the density magnitude
          if (is_weightfunction && extrapolate) {
             # Calculate constant C
             mu_j_i <- mu_samples[j, i]
             if (effect_direction == "negative") {
              mu_j_i_flipped <- -mu_j_i
             } else {
              mu_j_i_flipped <- mu_j_i  
             }

             temp_consts <- .dwnorm_fast.ss(
                x               = 0,
                mean            = mu_j_i_flipped,
                sd              = total_sd,
                omega           = omega_samples[j, , drop = FALSE],
                crit_x          = fit_data$crit_yi[, i],
                attach_constant = TRUE
             )
             constant <- attr(temp_consts, "constant")

             temp_densities_mat[i, ] <- temp_densities_mat[i, ] / constant
          }


        } else {
          # Weighted density (Fitted Model)
          # Use .dwnorm_fast.ss vectorization over x (z_sequence)

          mu_j_i <- mu_samples[j, i]
          y_seq  <- z_sequence * sei[i]
          if (effect_direction == "negative") {
            mu_j_i_flipped <- -mu_j_i
            y_seq_flipped  <- -y_seq
          } else {
            mu_j_i_flipped <- mu_j_i
            y_seq_flipped  <- y_seq
          }
          
          # .dwnorm_fast.ss can handle vector x
          dens_y <- .dwnorm_fast.ss(
            x      = y_seq_flipped,
            mean   = mu_j_i_flipped,
            sd     = total_sd,
            omega  = matrix(omega_samples[j, ], nrow = length(y_seq), ncol = ncol(omega_samples), byrow = TRUE),
            crit_x = fit_data$crit_yi[, i]
          )

          # transform density back: f_Z(z) = f_Y(z*se) * se
          temp_densities_mat[i, ] <- dens_y * sei[i]

        }
      }

      outcome_densities[j, ] <- colMeans(temp_densities_mat)
    }

    return(outcome_densities)
  }
}


# ---------------------------------------------------------------------------- #
# Graphical helper functions (reused verbatim from previous brma.zcurve.R)
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
      yaxt   = if (!is.null(dots[["yaxt"]]))   dots[["yaxt"]]   else "s"
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

.get_Soric_FDR <- function(EDR, sig_level) {
  ((1 / EDR) - 1) * (sig_level / (1 - sig_level))
}

# Logic from zcurve.R for binning
.zcurve_bins <- function(priors, from, to, by, length.out, type = "hist") {
  if (is.null(length.out)) {
    bins <- seq(from = from, to = to, by = by)
  } else {
    bins <- seq(from = from, to = to, length.out = length.out)
  }

  priors_bias <- priors[["bias"]]
  if (is.null(priors_bias)) return(bins)

  # Check if .zcurve_threshold exists
  if (exists(".zcurve_threshold", mode = "function")) {
    z_bounds <- .zcurve_threshold(priors)
    if (length(z_bounds) == 0) return(bins)

    z_bounds <- z_bounds[z_bounds > from & z_bounds < to]

    if (type == "hist") {
      # shift bins to match thresholds
      for (i in seq_along(z_bounds)) {
        i_larger <- which(bins > z_bounds[i])[1]
        if (is.na(i_larger)) next
        i_lower <- i_larger - 1
        if (bins[i_larger] - z_bounds[i] < z_bounds[i] - bins[i_lower]) {
          bins[i_larger] <- z_bounds[i]
        } else {
          bins[i_lower] <- z_bounds[i]
        }
      }
    } else if (type == "dens") {
      bins <- sort(unique(c(bins, z_bounds)))
    }
  }

  return(bins)
}
