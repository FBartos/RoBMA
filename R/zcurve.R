#' @title Transform RoBMA object into z-curve
#'
#' @description
#' This function transforms the estimated RoBMA model into a
#' z-curve object that can be further summarized and plotted.
#' Only available for normal-normal models estimated using the spike-and-slab
#' algorithm (i.e., \code{algorithm = "ss"}). See
#' \insertCite{bartos2025zcurve;textual}{RoBMA} and
#' \href{../doc/ZCurveDiagnostics.html}{\code{vignette("ZCurveDiagnostics", package = "RoBMA")}}
#' for more detail.
#'
#'
#' @param x A RoBMA object
#' @param significance_level Significance level used for computation of z-curve
#' estimates.
#' @param max_samples Maximum number of samples from the posterior distribution
#' that will be used for estimating z-curve estimates.
#'
#' @return \code{as_zcurve} returns a list of tables of class 'zcurve_RoBMA'.
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n,
#'              study_names = Anderson2010$labels, algorithm = "ss")
#'
#' zcurve_fit <- as_zcurve(fit)
#' summary(zcurve_fit)
#' }
#'
#' @export
as_zcurve <- function(x, significance_level = stats::qnorm(0.975), max_samples = 1000){

  if(x[["add_info"]][["algorithm"]] != "ss")
    stop("Predictions can only be computed for spike and slab models.")
  if(inherits(x, "BiBMA") || inherits(x, "BiBMA.reg"))
    stop("The true effects can only be computed for normal-normal (NoBMA / RoBMA) models.")
  BayesTools::check_real(significance_level, "significance_level", lower = 0)
  BayesTools::check_int(max_samples, "max_samples", lower = 10)

  # compute the main estimates and CIs
  # plotting densities are generated only if figures are requested
  out             <- .zcurve_fun(x, z_threshold = significance_level, max_samples = max_samples, conditional = FALSE)
  out_conditional <- try(.zcurve_fun(x, z_threshold = significance_level, max_samples = max_samples, conditional = TRUE))


  # compute data summaries
  outcome_data <- .get_outcome_data(x)

  # store new z-curve objects
  x[["coefficients"]] <- c("EDR" = mean(out[["EDR"]]))
  x[["zcurve"]]       <- list(
    estimates = list(
      EDR     = out[["EDR"]],
      weights = out[["weights"]]
    ),
    estimates_conditional = list(
      EDR     = out_conditional[["EDR"]],
      weights = out_conditional[["weights"]]
    ),
    data = list(
      significance_level = significance_level,
      z                  = outcome_data[["y"]] / outcome_data[["se"]],
      N_significant      = sum(abs(outcome_data[["y"]] / outcome_data[["se"]]) > significance_level),
      N_observed         = nrow(outcome_data)
    )
  )

  # add class
  class(x) <- c("zcurve_RoBMA", class(x))

  return(x)
}


#' @title Summarize fitted zcurve_RoBMA object
#'
#' @description \code{summary.zcurve_RoBMA} creates summary tables for a
#' zcurve_RoBMA object.
#'
#' @inheritParams summary.RoBMA
#'
#' @return \code{summary.RoBMA} returns a list of tables of class 'BayesTools_table'.
#'
#' @seealso [as_zcurve()]
#' @export
summary.zcurve_RoBMA <- function(object, conditional = FALSE, probs = c(.025, .975), ...){

  # compute summary for the z-curve object
  # most functions are based on the zcurve package

  # get proportion of significant results
  obs_proportion <- stats::prop.test(object$zcurve$data[["N_significant"]], object$zcurve$data[["N_observed"]], conf.level = 0.95)
  info_text      <- sprintf("Estimated using %1$i estimates, %2$i significant (ODR = %3$.2f, 95 CI [%4$.2f, %5$.2f]).",
                            object$zcurve$data[["N_observed"]],
                            object$zcurve$data[["N_significant"]],
                            obs_proportion$estimate,
                            obs_proportion$conf.int[1],
                            obs_proportion$conf.int[2]
  )

  estimates <- cbind.data.frame(
    "EDR"           = object$zcurve$estimates[["EDR"]],
    "Soric FDR"     = .get_Soric_FDR(object$zcurve$estimates[["EDR"]], stats::pnorm(object$zcurve$data[["significance_level"]], lower.tail = FALSE) * 2),
    "Missing N"     = (object$zcurve$estimates[["weights"]] - 1) * object$zcurve$data[["N_observed"]]
  )

  estimates <- BayesTools::ensemble_estimates_table(
    samples    = estimates,
    parameters = names(estimates),
    probs      = probs,
    title      = "Model-averaged estimates:",
    footnotes  = info_text,
    warnings   = .collect_errors_and_warnings(object)
  )

  if(conditional){

    estimates_conditional <- cbind.data.frame(
      "EDR"           = object$zcurve$estimates_conditional[["EDR"]],
      "Soric FDR"     = .get_Soric_FDR(object$zcurve$estimates_conditional[["EDR"]], stats::pnorm(object$zcurve$data[["significance_level"]], lower.tail = FALSE) * 2),
      "Missing N"     = (object$zcurve$estimates_conditional[["weights"]] - 1) * object$zcurve$data[["N_observed"]]
    )

    estimates_conditional <- BayesTools::ensemble_estimates_table(
      samples    = estimates_conditional,
      parameters = names(estimates_conditional),
      probs      = probs,
      title      = "Conditional estimates:",
      footnotes  = info_text,
      warnings   = .collect_errors_and_warnings(object)
    )

  }

  # create the output object
  output <- list(
    call       = object[["call"]],
    title      = "Z-curve",
    estimates  = estimates
  )

  if(conditional){
    output$estimates_conditional <- estimates_conditional
  }

  class(output) <- "summary.zcurve_RoBMA"
  return(output)
}

#' @title Create Z-Curve Meta-Analytic Plot
#'
#' @description
#' Plots a fitted \code{zcurve_RoBMA} object, visualizing the z-curve, model fit, extrapolation, and credible intervals.
#'
#'
#' @param x A zcurve_RoBMA object to be plotted.
#' @param max_samples Maximum number of posterior samples to use for plotting credible intervals. Defaults to 500.
#' @param plot_fit Should the model fit be included in the plot? Defaults to TRUE.
#' @param plot_extrapolation Should model extrapolation be included in the plot? Defaults to TRUE.
#' @param plot_CI Should credible intervals be included in the plot? Defaults to TRUE.
#' @param plot_thresholds Should significance thresholds be displayed in the plot? Defaults to TRUE.
#' @param from Lower bound of the z-value range for plotting. Defaults to -6.
#' @param to Upper bound of the z-value range for plotting. Defaults to 6.
#' @param by.hist Bin width for the histogram of observed z-values. Defaults to 0.5.
#' @param length.out.hist Number of bins for the histogram. If NULL, determined by by.hist. Defaults to NULL.
#' @param by.lines Step size for plotting model fit and extrapolation lines. Defaults to 0.05.
#' @param length.out.lines Number of points for plotting lines. If NULL, determined by by.lines. Defaults to NULL.
#' @param ... Additional arguments passed to the underlying plotting functions.
#' @inheritParams as_zcurve
#' @inheritParams plot.RoBMA
#' @inheritParams summary.RoBMA
#'
#' @return
#' Returns \code{NULL} if \code{plot_type = "base"}, or a \code{ggplot2} object if \code{plot_type = "ggplot2"}.
#'
#' @seealso [as_zcurve()], [lines.zcurve_RoBMA()], [hist.zcurve_RoBMA()]
#'
#' @examples \dontrun{
#' # using the example data from Anderson et al. 2010 and fitting the default model
#' # (note that the model can take a while to fit)
#' fit <- RoBMA(r = Anderson2010$r, n = Anderson2010$n,
#'              study_names = Anderson2010$labels, algorithm = "ss")
#'
#' zcurve_fit <- as_zcurve(fit)
#' plot(zcurve_fit)
#' }
#'
#' @export
plot.zcurve_RoBMA <- function(x, conditional = FALSE, plot_type = "base",
                              probs = c(.025, .975), max_samples = 500,
                              plot_fit = TRUE, plot_extrapolation = TRUE, plot_CI = TRUE, plot_thresholds = TRUE,
                              from = -6, to = 6, by.hist = 0.5, length.out.hist = NULL, by.lines = 0.05, length.out.lines = NULL, ...){

  # plots the z-curve object
  # most functions are based on the zcurve package
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_real(from, "from", allow_NULL = TRUE)
  BayesTools::check_real(to, "to", allow_NULL = TRUE)
  BayesTools::check_real(by.hist, "by.hist", allow_NULL = TRUE)
  BayesTools::check_real(length.out.hist, "length.out.hist", allow_NULL = TRUE)
  BayesTools::check_real(by.lines, "by.lines", allow_NULL = TRUE)
  BayesTools::check_real(length.out.lines, "length.out.lines", allow_NULL = TRUE)
  BayesTools::check_int(max_samples, "max_samples", lower = 10)
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_bool(plot_extrapolation, "plot_extrapolation")
  BayesTools::check_bool(plot_CI, "plot_CI")
  BayesTools::check_bool(plot_thresholds, "plot_thresholds")
  BayesTools::check_real(probs, "probs", lower = 0, upper = 1, check_length = 2)

  # construct the plot
  # get the line values first so we can set-up ylim of the histogram
  ymax <- 0
  if(plot_fit){
    lines_fit <- lines.zcurve_RoBMA(x, conditional = conditional, plot_type = plot_type,
                                    probs = probs, max_samples = max_samples,
                                    extrapolate = FALSE, plot_CI = plot_CI,
                                    from = from, to = to, by = by.lines, length.out = length.out.lines, as_data = TRUE)
    ymax <- max(c(ymax, lines_fit$y, if(plot_CI) lines_fit$y_uCI))
  }
  if(plot_extrapolation){
    lines_extrapolation <- lines.zcurve_RoBMA(x, conditional = conditional, plot_type = plot_type,
                                              probs = probs, max_samples = max_samples,
                                              extrapolate = TRUE, plot_CI = plot_CI,
                                              from = from, to = to, by = by.lines, length.out = length.out.lines, as_data = TRUE)
    ymax <- max(c(ymax, lines_extrapolation$y, if(plot_CI) lines_extrapolation$y_uCI))
  }


  plot <- hist.zcurve_RoBMA(x, plot_type = plot_type, from = from, to = to, by = by.hist, length.out = length.out.hist,
                            ylim = if(ymax != 0) c(0, ymax) else NULL, plot_thresholds = plot_thresholds)

  # add plot lines
  if(plot_fit){
    if(plot_type == "base"){
      if(plot_CI){
        graphics::polygon(
          c(lines_fit$x,     rev(lines_fit$x)),
          c(lines_fit$y_lCI, rev(lines_fit$y_uCI)),
          border = NA, col = scales::alpha("black", 0.40))
      }
      graphics::lines(lines_fit$x, lines_fit$y, lwd = 2, col = "black", lty = 1)
    }else if(plot_type == "ggplot"){
      if(plot_CI){
        plot <- plot + ggplot2::geom_ribbon(
          ggplot2::aes(
            x    = lines_fit$x,
            ymin = lines_fit$y_lCI,
            ymax = lines_fit$y_uCI),
          fill = scales::alpha("black", 0.40)
        )
      }
      plot <- plot + ggplot2::geom_line(
        ggplot2::aes(
          x = lines_fit$x,
          y = lines_fit$y),
        color     = "black",
        linewidth = 1,
        linetype  = 1
      )
    }
  }

  # add extrapolation lines
  if(plot_extrapolation){
    if(plot_type == "base"){
      if(plot_CI){
        graphics::polygon(
          c(lines_extrapolation$x,     rev(lines_extrapolation$x)),
          c(lines_extrapolation$y_lCI, rev(lines_extrapolation$y_uCI)),
          border = NA, col = scales::alpha("blue", 0.40))
      }
      graphics::lines(lines_extrapolation$x, lines_extrapolation$y, lwd = 2, col = "blue", lty = 1)
    }else if(plot_type == "ggplot"){
      if(plot_CI){
        plot <- plot + ggplot2::geom_ribbon(
          ggplot2::aes(
            x    = lines_extrapolation$x,
            ymin = lines_extrapolation$y_lCI,
            ymax = lines_extrapolation$y_uCI),
          fill = scales::alpha("blue", 0.40)
        )
      }
      plot <- plot + ggplot2::geom_line(
        ggplot2::aes(
          x = lines_extrapolation$x,
          y = lines_extrapolation$y),
        color     = "blue",
        linewidth = 1,
        linetype  = 1
      )
    }
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}

#' @title Create Histogram of Z-Statistics
#'
#' @description
#' Plots a histogram of observed z-values for a fitted \code{zcurve_RoBMA} object, with options to customize the plotting range, bin width, and display of significance thresholds.
#'
#' @param x A \code{zcurve_RoBMA} object containing the fitted model.
#' @param by Numeric value specifying the bin width for the histogram. Defaults to 0.5.
#' @param length.out Optional integer specifying the number of bins. If NULL, determined by \code{by}. Defaults to NULL.
#' @param plot_thresholds Logical; should significance thresholds be displayed on the plot? Defaults to TRUE.
#' @param ... Additional arguments passed to the underlying plotting functions.
#' @inheritParams as_zcurve
#' @inheritParams plot.RoBMA
#' @inheritParams summary.RoBMA
#' @inheritParams plot.zcurve_RoBMA
#'
#' @return
#' Returns \code{NULL} if \code{plot_type = "base"}, or a \code{ggplot2} object if \code{plot_type = "ggplot2"}.
#'
#' @seealso [as_zcurve()], [plot.zcurve_RoBMA()], [hist.zcurve_RoBMA()]
#'
#' @export
hist.zcurve_RoBMA  <- function(x, plot_type = "base", from = -6, to = 6, by = 0.5, length.out = NULL, plot_thresholds = TRUE, ...){

  # most functions are based on the zcurve package
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_real(from, "from", allow_NULL = TRUE)
  BayesTools::check_real(to, "to", allow_NULL = TRUE)
  BayesTools::check_real(by, "by", allow_NULL = TRUE)
  BayesTools::check_real(length.out, "length.out", allow_NULL = TRUE)
  BayesTools::check_bool(plot_thresholds, "plot_thresholds")

  # extract the data
  dots <- list(...)
  z    <- x$zcurve$data[["z"]]

  # specify plotting range
  z_sequence <- .zcurve_bins(priors = x[["priors"]], from = from, to = to, by = by, length.out = length.out, type = "hist")
  z_in_range <- z >= min(z_sequence) & z <= max(z_sequence)

  # compute the number of z-statistics outside of the plotting range
  if(sum(!z_in_range) > 0){
    message(sprintf("%1$i z-statistics are out of the plotting range", sum(!z_in_range)))
  }

  # create z-curve histogram data
  z_hist <- graphics::hist(z[z_in_range], breaks = z_sequence, plot = FALSE)

  # adjust histogram for values outside of plotting range
  z_hist$density <- z_hist$density * mean(z_in_range)

  # create plotting objects
  df_hist <- data.frame(
    x       = z_hist$mids,
    density = z_hist$density,
    breaks  = diff(z_hist$breaks)
  )

  # return the plotting object if requested
  if(isTRUE(dots[["as_data"]])){
    return(df_hist)
  }

  # create the plot otherwise
  if(!is.null(dots[["ylim"]])){
    ylim <- range(dots[["ylim"]], max(z_hist$density))
  }else{
    ylim <- c(0, max(z_hist$density))
  }
  if(!is.null(dots[["xlab"]])){
    xlab <- dots[["xlab"]]
  }else{
    xlab <- "Z-Statistic"
  }

  if(plot_type == "ggplot"){
    out <- ggplot2::ggplot() +
      ggplot2::geom_col(
        ggplot2::aes(
          x = df_hist$x,
          y = df_hist$density),
        fill  = NA,
        color = "black",
        width = df_hist$breaks
      ) +
      ggplot2::labs(x = xlab, y = "Density")
  } else {
    graphics::plot(z_hist, freq = FALSE, las = 1, density = 0, angle = 0,
                   border = "black", xlab = xlab, main = "", ylim = ylim)
  }

  if(plot_thresholds){
    tresholds <- .zcurve_threshold(x[["priors"]])
    if(length(tresholds) > 0){
      if(plot_type == "base"){
        graphics::abline(v = tresholds, col = "red", lty = 3)
      }else if(plot_type == "ggplot"){
        out <- out + ggplot2::geom_vline(xintercept = tresholds, color = "red", linetype = "dashed")
      }
    }
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(out)
  }
}

#' @title Add Lines With Posterior Predictive Distribution of Z-Statistics
#'
#' @description
#' Adds lines to a plot of a fitted zcurve_RoBMA object. This function is typically used to overlay additional information or model fits on an existing plot.
#'
#' @param extrapolate Logical indicating whether to extrapolate values beyond the observed data range.
#' @param by Numeric value specifying the increment for the sequence.
#' @param length.out Optional integer specifying the desired length of the output sequence.
#' @param col Color of the plotted line.
#' @inheritParams as_zcurve
#' @inheritParams plot.RoBMA
#' @inheritParams summary.RoBMA
#' @inheritParams plot.zcurve_RoBMA
#'
#' @seealso [as_zcurve()], [plot.zcurve_RoBMA()], [hist.zcurve_RoBMA()]
#'
#' @export
lines.zcurve_RoBMA <- function(x, conditional = FALSE, plot_type = "base",
                               probs = c(.025, .975), max_samples = 500,
                               extrapolate = FALSE, plot_CI = TRUE,
                               from = -6, to = 6, by = 0.05, length.out = NULL, col = if(extrapolate) "blue" else "black", ...){

  # most functions are based on the zcurve package
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_real(from, "from", allow_NULL = TRUE)
  BayesTools::check_real(to, "to", allow_NULL = TRUE)
  BayesTools::check_real(by, "by", allow_NULL = TRUE)
  BayesTools::check_real(length.out, "length.out", allow_NULL = TRUE)
  BayesTools::check_int(max_samples, "max_samples", lower = 10)
  BayesTools::check_bool(conditional, "conditional")
  BayesTools::check_bool(extrapolate, "extrapolate")
  BayesTools::check_bool(plot_CI, "plot_CI")
  BayesTools::check_real(probs, "probs", lower = 0, upper = 1, check_length = 2)

  # extract the data
  dots <- list(...)
  z    <- x$zcurve$data[["z"]]

  # specify plotting range
  z_sequence <- .zcurve_bins(priors = x[["priors"]], from = from, to = to, by = by, length.out = length.out, type = "dens")

  # create z-curve histogram data
  z_density <- .zcurve_fun(x, z_sequence = z_sequence, max_samples = max_samples, conditional = conditional, extrapolate = extrapolate)

  # create plotting objects
  df_density <- data.frame(
    x     = z_sequence,
    y     = colMeans(z_density),
    y_lCI = apply(z_density, 2, stats::quantile, probs = probs[1]),
    y_uCI = apply(z_density, 2, stats::quantile, probs = probs[2])
  )

  # return the plotting object if requested
  if(isTRUE(dots[["as_data"]])){
    return(df_density)
  }

  # create the plot otherwise
  if(!is.null(dots[["lwd"]])){
    lwd <- dots[["lwd"]]
  }else{
    lwd <- 1
  }
  if(!is.null(dots[["lty"]])){
    lty <- dots[["lty"]]
  }else{
    lty <- 1
  }
  if(!is.null(dots[["alpha"]])){
    alpha <- dots[["alpha"]]
  }else{
    alpha <- 0.4
  }

  if(plot_type == "ggplot"){
    out <- list()
    if(plot_CI){
      out[[1]] <- ggplot2::geom_ribbon(
        ggplot2::aes(
          x    = df_density$x,
          ymin = df_density$y_lCI,
          ymax = df_density$y_uCI),
        fill = scales::alpha(col, alpha)
      )
    }
    out[[length(out) + 1]] <- ggplot2::geom_line(
      ggplot2::aes(
        x = df_density$x,
        y = df_density$y),
      color     = col,
      linewidth = lwd,
      linetype  = lty
    )
  }else{
    if(plot_CI){
      graphics::polygon(
        c(df_density$x,     rev(df_density$x)),
        c(df_density$y_lCI, rev(df_density$y_uCI)),
        border = NA, col = scales::alpha(col, alpha))
    }
    graphics::lines(df_density$x, df_density$y, lwd = lwd, col = col, lty = lty)
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(out)
  }
}



.zcurve_bins      <- function(priors, from, to, by, length.out, type = "hist"){

  if(is.null(length.out)){
    bins <- seq(from = from, to = to, by = by)
  }else{
    bins <- seq(from = from, to = to, length.out = length.out)
  }

  priors_bias <- priors[["bias"]]

  # return simple binning in case of no bias
  if(is.null(priors_bias)){
    return(bins)
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
.zcurve_threshold <- function(priors){

  priors_bias <- priors[["bias"]]

  # return simple binning in case of no bias
  if(is.null(priors_bias)){
    return()
  }

  priors_weightfunctions <- priors_bias[sapply(priors_bias, is.prior.weightfunction)]

  # return simple binning in case of no weightfunctions
  if(length(priors_weightfunctions) == 0){
    return()
  }

  # obtain tresholds from the specified priors
  p_bounds <- BayesTools::weightfunctions_mapping(priors_weightfunctions, cuts_only = TRUE)
  z_bounds <- stats::qnorm(rev(p_bounds), 0, 1, lower.tail = FALSE)
  z_bounds <- z_bounds[!is.infinite(z_bounds)]

  return(z_bounds)
}

#' @title Prints summary object for zcurve_RoBMA method
#'
#' @param x a summary of a zcurve_RoBMA object
#' @param ... additional arguments
#'
#'
#' @return \code{summary.zcurve_RoBMA} invisibly returns the print statement.
#'
#' @seealso [as_zcurve()]
#' @export
print.summary.zcurve_RoBMA <- function(x, ...){

  cat("Call:\nas_zcurve: ")
  print(x[["call"]])

  cat("\n")
  cat(x[["title"]])

  for(type in c("estimates", "estimates_conditional")){
    if(!is.null(x[[type]])){
      cat("\n")
      print(x[[type]])
    }
  }
}


# imported from zcurve
.get_Soric_FDR     <- function(EDR, sig_level){
  ((1/EDR) - 1)*(sig_level/(1-sig_level))
}



#' @title Prints a fitted zcurve_RoBMA object
#'
#' @param x a fitted zcurve_RoBMA object.
#' @param ... additional arguments.
#'
#'
#' @return \code{print.zcurve_RoBMA} invisibly returns the print statement.
#'
#' @seealso [as_zcurve()]
#' @export
print.zcurve_RoBMA <- function(x, ...){
  cat("Call:\nas_zcurve: ")
  print(x$call)
  cat("\nEstimates:\n")
  print(stats::coef(x))
}


.zcurve_fun <- function(x, z_threshold = NULL, z_sequence = NULL, max_samples = 1000, conditional = FALSE, extrapolate = FALSE){

  # tree modes of the function
  # z_threshold is specified: compute the EDR (automatically sets extrapolate to TRUE to incorporate weights)
  # extrapolate = FALSE: produces z-density of the fitted model as is == fit assessment
  #                      - requires incorporating publication bias indices (i.e., PET, PEESE, weighted likelihoods)
  # extrapolate = TRUE: produces z-density of the expected unbiased distribution
  #                     - assumes that selection and SE dependent bias would not come into play
  #                     - wherever weighted likelihoods are involved, they need to be inverse-weighted to extrapolate to missing estimates
  if(!is.null(z_threshold)){
    extrapolate <- TRUE
  }

  # get the model fitting scale
  if (is.BiBMA(x)) {
    model_scale <- "logOR"
  } else {
    model_scale <- x$add_info[["effect_measure"]]
  }

  # extract posterior samples (and obtain conditional indicator)
  posterior_samples <- suppressWarnings(coda::as.mcmc(x[["model"]][["fit"]]))
  priors            <- x[["priors"]]

  # subset the samples to speed up the computation
  if(conditional){
    # if conditional output is to be provided, condition first and then subset
    # otherwise there might be no conditional samples left

    # select the indicator
    if(.is_regression(x)){
      mu_indicator <- posterior_samples[,"mu_intercept_indicator"]
      mu_is_null   <- attr(x[["model"]]$priors$terms[["intercept"]], "components") == "null"
    }else{
      mu_indicator <- posterior_samples[,"mu_indicator"]
      mu_is_null   <- attr(x[["model"]]$priors$mu, "components") == "null"
    }
    mu_indicator <- mu_indicator %in% which(!mu_is_null)

    selected_samples_ind <- unique(round(seq(from = 1, to = sum(mu_indicator), length.out = max_samples)))
    selected_samples_ind <- seq_len(nrow(posterior_samples))[mu_indicator][selected_samples_ind]
    n_samples            <- length(selected_samples_ind)
    posterior_samples    <- posterior_samples[selected_samples_ind,,drop=FALSE]

  }else{

    selected_samples_ind <- unique(round(seq(from = 1, to = nrow(posterior_samples), length.out = max_samples)))
    n_samples            <- length(selected_samples_ind)
    posterior_samples    <- posterior_samples[selected_samples_ind,,drop=FALSE]
  }

  # dispatch between meta-regression / meta-analysis input
  if(.is_regression(x)){
    newdata.predictors <- do.call(cbind.data.frame, x$data[["predictors"]])
  }
  newdata.outcome <- .get_outcome_data(x)

  # obtain the (study-specific) mu estimate
  # meta-regression and meta-analysis separately
  if(.is_regression(x)){

    mu_samples  <- t(BayesTools::JAGS_evaluate_formula(
      fit         = x$model$fit,
      formula     = x$formula,
      parameter   = "mu",
      data        = newdata.predictors,
      prior_list  = attr(x$model$fit, "prior_list")
    ))[selected_samples_ind,,drop=FALSE]

  }else{

    mu_samples   <- matrix(posterior_samples[,"mu"], ncol = nrow(newdata.outcome), nrow = n_samples)

  }

  # transform the samples to the model fitting scale (the data are at the model fitting scale)
  mu_samples  <- .scale(mu_samples,  x$add_info[["output_scale"]], model_scale)

  # add conditional warnings
  if(conditional){
    if(sum(mu_indicator) < 100)
      warning(gettextf("There is only a very small number of posterior samples (%1s) assuming presence of the effect. The resulting estimates are not reliable.",
                       sum(mu_indicator)), call. = FALSE, immediate. = TRUE)
    if(sum(mu_indicator) <= 2)
      stop("Less or equal to 2 posterior samples assuming presence of the effects. The estimates could not be computed.")
  }

  # adjust effect samples and data to match original observed data direction
  # (to undo fitting adjustment for one-sided weightfunction)
  if(!is.null(x$add_info[["effect_direction"]]) && x$add_info[["effect_direction"]] == "negative"){
    newdata.outcome[,"y"] <- -newdata.outcome[,"y"]
    mu_samples            <- -mu_samples
    if(!is.null(z_sequence)){
      z_sequence          <- -z_sequence
    }
  }

  # extract the list of priors
  priors_bias  <- priors[["bias"]]
  if(!extrapolate){
    # add PET/PEESE adjustment
    if(any(sapply(priors_bias, is.prior.PET))){
      PET_samples <- posterior_samples[,"PET"]
      # PET is scale invariant (no-scaling needed)
    }else{
      PET_samples <- rep(0, n_samples)
    }
    if(any(sapply(priors_bias, is.prior.PEESE))){
      PEESE_samples <- posterior_samples[,"PEESE"]
      # PEESE scales with the inverse
      PEESE_samples <- .scale(PEESE_samples, model_scale, x$add_info[["output_scale"]])
    }else{
      PEESE_samples <- rep(0, n_samples)
    }

    for(i in seq_len(ncol(mu_samples))){
      mu_samples[,i] <- mu_samples[,i] + PET_samples * newdata.outcome[i,"se"] + PEESE_samples * newdata.outcome[i,"se"]^2
    }
  }

  if(inherits(x, "NoBMA.reg") || inherits(x, "BiBMA.reg") || (length(priors[["bias"]]) == 1 && is.prior.none(priors[["bias"]][[1]]))){

    # predicting responses without selection models does not require incorporating the between-study random-effects
    # (the marginalized and non-marginalized parameterization are equivalent)
    tau_samples <- posterior_samples[,"tau"]
    tau_samples <- .scale(tau_samples, x$add_info[["output_scale"]], model_scale)

    # dispatch between EDR computation / densities
    if(!is.null(z_threshold)){

      # compute the proportion of estimates larger than threshold under the population parameters
      outcome_thresholds <- rep(NA, nrow(mu_samples))
      outcome_weights    <- rep(NA, nrow(mu_samples))

      for(j in 1:nrow(mu_samples)){
        # create containers for temporal samples from the posterior distribution
        temp_thresholds <- rep(NA, ncol(mu_samples))
        temp_weights    <- rep(NA, ncol(mu_samples))

        # compute the densities and threshold for each observation
        for(i in seq_len(ncol(mu_samples))){
          temp_thresholds[i] <-
            stats::pnorm(z_threshold * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = FALSE) +
            stats::pnorm(-z_threshold * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = TRUE)
          temp_weights[i]    <- 1
        }

        # store the results
        outcome_thresholds[j] <- mean(temp_thresholds)
        outcome_weights[j]    <- mean(temp_weights)
      }

    }else{

      # compute the density ad specified support
      # outcome_lower     <- rep(NA, nrow(mu_samples))
      # outcome_higher    <- rep(NA, nrow(mu_samples))
      outcome_densities <- matrix(NA, nrow = nrow(mu_samples), ncol = length(z_sequence))

      for(j in 1:nrow(mu_samples)){
        # create containers for temporal samples from the posterior distribution
        # temp_lower     <- rep(NA, ncol(mu_samples))
        # temp_higher    <- rep(NA, ncol(mu_samples))
        temp_densities <- matrix(NA, nrow = ncol(mu_samples), ncol = length(z_sequence))

        # compute the densities and threshold for each observation
        for(i in seq_len(ncol(mu_samples))){
          # temp_lower[i]      <- stats::pnorm(z_sequence[1]                  * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = TRUE)
          # temp_higher[i]     <- stats::pnorm(z_sequence[length(z_sequence)] * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = FALSE)
          temp_densities[i,] <- stats::dnorm(z_sequence * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2)) * newdata.outcome[i,"se"]
          # the density needs to be transformed due to the support change
        }

        # store the results
        # outcome_lower[j]      <- mean(temp_lower)
        # outcome_higher[j]     <- mean(temp_higher)
        outcome_densities[j,] <- colMeans(temp_densities)
      }

    }


  }else{

    # required for study ids / crit_x values in selection models
    fit_data <- .fit_data_ss(
      data             = newdata.outcome,
      priors           = priors,
      effect_direction = x$add_info[["effect_direction"]],
      prior_scale      = x$add_info[["prior_scale"]],
      weighted         = FALSE,
      weighted_type    = FALSE,
      multivariate     = .is_multivariate(x)
    )

    # predicting response requires incorporating the between-study random effects if selection models are present
    # (we use approximate selection likelihood which samples the true study effects instead of marginalizing them)
    if(.is_multivariate(x)){

      tau_samples <- posterior_samples[,"tau"]
      rho_samples <- posterior_samples[,"rho"]
      # deal with computer precision errors from JAGS
      rho_samples[rho_samples>1] <- 1
      rho_samples[rho_samples<0] <- 0
      # tau_within  = tau * sqrt(rho)
      # tau_between = tau * sqrt(1-rho)
      tau_within_samples  <- tau_samples * sqrt(rho_samples)
      tau_between_samples <- tau_samples * sqrt(1-rho_samples)
      gamma_samples       <- posterior_samples[,grep("gamma", colnames(posterior_samples)),drop = FALSE]

      tau_between_samples <- .scale(tau_between_samples, x$add_info[["output_scale"]], model_scale)
      tau_within_samples  <- .scale(tau_within_samples,  x$add_info[["output_scale"]], model_scale)

      # incorporate within study heterogeneity into the predictor
      # either estimated for prediction on the same data or integrated over for new data
      for(i in seq_len(nrow(newdata.outcome))){
        mu_samples[,i] <- mu_samples[,i] + gamma_samples[,fit_data$study_ids[i]] * tau_within_samples
      }

      # tau_between samples work as tau for the final sampling step
      tau_samples <- tau_between_samples

    }else{

      tau_samples  <- posterior_samples[,"tau"]
      tau_samples  <- .scale(tau_samples,  x$add_info[["output_scale"]], model_scale)

    }

    # selection models are sampled separately for increased efficiency
    bias_indicator           <- posterior_samples[,"bias_indicator"]
    weightfunction_indicator <- bias_indicator %in% which(sapply(priors[["bias"]], is.prior.weightfunction))


    # dispatch between EDR computation / densities
    if(!is.null(z_threshold)){

      # compute the proportion of estimates larger than threshold under the population parameters
      outcome_thresholds <- rep(NA, nrow(mu_samples))
      outcome_weights    <- rep(NA, nrow(mu_samples))

      for(j in 1:nrow(mu_samples)){
        # create containers for temporal samples from the posterior distribution
        temp_thresholds <- rep(NA, ncol(mu_samples))
        temp_weights    <- rep(NA, ncol(mu_samples))

        for(i in seq_len(ncol(mu_samples))){

          # sample normal models/PET/PEESE
          if(!weightfunction_indicator[j]){
            temp_thresholds[i] <-
              stats::pnorm(z_threshold * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = FALSE) +
              stats::pnorm(-z_threshold * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = TRUE)
            temp_weights[i]    <- 1
          }

          # sample selection models
          if(weightfunction_indicator[j]){
            temp_consts <- .dwnorm_fast.ss(
              x      = 0,
              mean   = mu_samples[j,i],
              sd     = sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2),
              omega  = posterior_samples[j, grep("omega", colnames(posterior_samples)),drop = FALSE],
              crit_x = fit_data$crit_y[, i],
              attach_constant = TRUE
            )
            temp_thresholds[i] <-
              stats::pnorm(z_threshold * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = FALSE) +
              stats::pnorm(-z_threshold * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = TRUE)
            temp_weights[i]    <- 1/attr(temp_consts, "constant")
          }
        }

        # store the results
        outcome_thresholds[j] <- stats::weighted.mean(temp_thresholds, temp_weights)
        outcome_weights[j]    <- mean(temp_weights)
      }

    }else{

      # compute the density add specified support
      # outcome_lower     <- rep(NA, nrow(mu_samples))
      # outcome_higher    <- rep(NA, nrow(mu_samples))
      outcome_densities <- matrix(NA, nrow = nrow(mu_samples), ncol = length(z_sequence))

      for(j in 1:nrow(mu_samples)){
        # create containers for temporal samples from the posterior distribution
        # temp_lower     <- rep(NA, ncol(mu_samples))
        # temp_higher    <- rep(NA, ncol(mu_samples))
        temp_densities <- matrix(NA, nrow = ncol(mu_samples), ncol = length(z_sequence))

        for(i in seq_len(ncol(mu_samples))){

          # sample normal models/PET/PEESE
          if(!weightfunction_indicator[j]){
            # temp_lower[i]      <- stats::pnorm(z_sequence[1]                  * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = TRUE)
            # temp_higher[i]     <- stats::pnorm(z_sequence[length(z_sequence)] * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2), lower.tail = FALSE)
            temp_densities[i,] <- stats::dnorm(z_sequence * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2)) * newdata.outcome[i,"se"]
          }

          # sample selection models
          if(weightfunction_indicator[j]){
            # # if including probability > & < than plotting range, those would also need to be re-standardized
            # temp_lower[i]  <- .pwnorm_fast.ss(
            #   q          = z_sequence[1] * newdata.outcome[i,"se"],
            #   mean       = mu_samples[j,i],
            #   sd         = sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2),
            #   omega      = posterior_samples[j, grep("omega", colnames(posterior_samples)),drop = FALSE],
            #   crit_x     = fit_data$crit_y[, i, drop=FALSE],
            #   lower.tail = TRUE)
            # temp_higher[i] <- .pwnorm_fast.ss(
            #   q          = z_sequence[length(z_sequence)] * newdata.outcome[i,"se"],
            #   mean       = mu_samples[j,i],
            #   sd         = sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2),
            #   omega      = posterior_samples[j, grep("omega", colnames(posterior_samples)),drop = FALSE],
            #   crit_x     = fit_data$crit_y[, i, drop=FALSE],
            #   lower.tail = FALSE)


            # the density needs to be transformed due to the support change
            # (+ extrapolation by removing standardizing constant in case of selection models)
            if(extrapolate){
              temp_consts <- .dwnorm_fast.ss(
                x      = 0,
                mean   = mu_samples[j,i],
                sd     = sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2),
                omega  = posterior_samples[j, grep("omega", colnames(posterior_samples)),drop = FALSE],
                crit_x = fit_data$crit_y[, i],
                attach_constant = TRUE
              )
              temp_densities[i,] <- stats::dnorm(z_sequence * newdata.outcome[i,"se"], mu_samples[j,i], sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2)) * newdata.outcome[i,"se"] / attr(temp_consts, "constant")
            }else{
              temp_out <- .dwnorm_fast.ss(
                x      = z_sequence * newdata.outcome[i,"se"],
                mean   = mu_samples[j,i],
                sd     = sqrt(tau_samples[j]^2 + newdata.outcome[i,"se"]^2),
                omega  = matrix(posterior_samples[j, grep("omega", colnames(posterior_samples)),drop = FALSE],
                                nrow = length(z_sequence), ncol = length(grep("omega", colnames(posterior_samples))), byrow = TRUE),
                crit_x = fit_data$crit_y[, i]
              )
              temp_densities[i,] <- temp_out * newdata.outcome[i,"se"]
            }

          }
        }

        # store the results
        # outcome_lower[j]      <- mean(temp_lower)
        # outcome_higher[j]     <- mean(temp_higher)
        outcome_densities[j,] <- colMeans(temp_densities)
      }

    }

  }

  if(!is.null(z_threshold)){
    return(list(
      EDR     = outcome_thresholds,
      weights = outcome_weights
    ))
  }else{
    return(outcome_densities)
  }
}
