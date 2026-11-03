.get_dots_hist_zplot <- function(dots, plot_type = "base", max_density = NULL) {

  .check_plot_limits(dots[["xlim"]], "xlim")
  ylim <- dots[["ylim"]]
  if (is.null(ylim) && !is.null(max_density) &&
      (plot_type == "base" || !is.null(dots[[".zplot_auto_ymax"]]))) {
    ylim <- c(0, max(c(max_density, dots[[".zplot_auto_ymax"]])))
    # an empty plotting range has no density to scale the axis to
    if (ylim[2] == 0) {
      ylim[2] <- 1
    }
  }
  if (plot_type == "base") {
    hist_dots <- list(
      border = if (!is.null(dots[["border"]])) dots[["border"]] else "gray60",
      col    = if (!is.null(dots[["col"]]))    dots[["col"]]    else NA,
      xlab   = if (!is.null(dots[["xlab"]]))   dots[["xlab"]]   else "Z-Statistic",
      ylab   = if (!is.null(dots[["ylab"]]))   dots[["ylab"]]   else "Density",
      main   = if (!is.null(dots[["main"]]))   dots[["main"]]   else "",
      xlim   = dots[["xlim"]],
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
      main  = if (!is.null(dots[["main"]]))  dots[["main"]]  else "",
      xlim  = dots[["xlim"]],
      ylim  = ylim
    )
  }
  return(hist_dots)
}


# ---------------------------------------------------------------------------- #
# .get_dots_lines_zplot
# ---------------------------------------------------------------------------- #
#
# Extracts density line graphical parameters with defaults.
#
# @param dots      list of user-supplied graphical parameters
# @param plot_type "base" or "ggplot"
# @param col       default line color
# @param alpha     default alpha for CI ribbons
# @param defaults  include default values, or normalize supplied styles only
#
# @return list of graphical parameters appropriate for plot_type
#
# ---------------------------------------------------------------------------- #
.get_dots_lines_zplot <- function(dots, plot_type = "base", col = "black",
                                   alpha = 0.40, defaults = TRUE) {
  if (plot_type == "base") {
    line_dots <- list(
      lwd   = if (!is.null(dots[["lwd"]]))   dots[["lwd"]]   else if (defaults) 2,
      col   = if (!is.null(dots[["col"]]))   dots[["col"]]   else if (defaults) col,
      lty   = if (!is.null(dots[["lty"]]))   dots[["lty"]]   else if (defaults) 1,
      alpha = if (!is.null(dots[["alpha"]])) dots[["alpha"]] else if (defaults) alpha
    )
  } else if (plot_type == "ggplot") {
    line_dots <- list(
      linewidth = if (!is.null(dots[["linewidth"]])) dots[["linewidth"]] else if (!is.null(dots[["lwd"]])) dots[["lwd"]] else if (defaults) 1,
      color     = if (!is.null(dots[["color"]]))     dots[["color"]]     else if (!is.null(dots[["col"]])) dots[["col"]] else if (defaults) col,
      linetype  = if (!is.null(dots[["linetype"]])) dots[["linetype"]]  else if (!is.null(dots[["lty"]])) dots[["lty"]] else if (defaults) 1,
      alpha     = if (!is.null(dots[["alpha"]]))    dots[["alpha"]]     else if (defaults) alpha
    )
  }
  if (!defaults) {
    line_dots <- line_dots[!vapply(line_dots, is.null, logical(1L))]
  }
  return(line_dots)
}


# ---------------------------------------------------------------------------- #
# .get_dots_thresholds_zplot
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
.get_dots_thresholds_zplot <- function(dots, plot_type = "base") {
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
# .zplot_check_grid
# ---------------------------------------------------------------------------- #
#
# Validates the z-value range and grid controls shared by the zplot methods.
#
# @param names argument names of the step and length controls in messages
#
# ---------------------------------------------------------------------------- #
.zplot_check_grid <- function(from, to, by, length.out,
                              names = c("by", "length.out")) {

  BayesTools::check_real(from, "from", allow_NA = FALSE)
  BayesTools::check_real(to, "to", allow_NA = FALSE)
  if (!is.finite(from) || !is.finite(to) || from >= to) {
    stop("'from' and 'to' must be finite with 'from' smaller than 'to'.", call. = FALSE)
  }
  BayesTools::check_int(length.out, names[2], lower = 1,
                        allow_NULL = TRUE, allow_NA = FALSE)
  if (is.null(length.out)) {
    BayesTools::check_real(by, names[1], allow_NA = FALSE)
    if (!is.finite(by) || by <= 0) {
      stop("'", names[1], "' must be a positive finite number.", call. = FALSE)
    }
  }

  return(invisible(TRUE))
}


# ---------------------------------------------------------------------------- #
# .zplot_bins
# ---------------------------------------------------------------------------- #
#
# Creates bin sequence for zplot histograms and densities.
#
# For histograms, bins are shifted to align with significance thresholds
# when selection model priors are present. For densities, threshold
# values are added to the sequence to ensure accurate representation.
#
# @param from       lower bound of z-value range
# @param to         upper bound of z-value range
# @param by         step size (NULL to use length.out)
# @param length.out number of bins (NULL to use by)
# @param type       "hist" or "dens"
# @param thresholds signed z-statistic thresholds of the selection kernel (see
#                   `.zplot_threshold()`), empty for none
#
# @return numeric vector of bin boundaries or density evaluation points
#
# ---------------------------------------------------------------------------- #
.zplot_bins      <- function(from, to, by, length.out, type = "hist",
                             thresholds){

  if(is.null(length.out)){
    bins <- seq(from = from, to = to, by = by)
    # keep the requested upper limit when 'by' does not divide the range,
    # including a step wider than the range (a single boundary)
    if (length(bins) == 1L || to - bins[length(bins)] > 1e-10 * by) {
      bins <- c(bins, to)
    } else {
      bins[length(bins)] <- to
    }
  }else{
    # histogram 'length.out' counts bins, which have one more boundary
    length.out <- as.integer(round(length.out))
    bins <- seq(from = from, to = to,
                length.out = if (type == "hist") length.out + 1L else length.out)
  }

  # return simple binning in case of no weightfunctions
  if(length(thresholds) == 0){
    return(bins)
  }

  # retain bounds within the plotting range
  z_bounds <- thresholds[thresholds > from & thresholds < to]

  if(type == "hist"){
    # for histogram, shift the specified bin boundaries to the closest z-treshold
    original_bins <- bins
    claimed <- c(1L, length(bins))
    for(i in seq_along(z_bounds)){

      closest <- which.min(abs(original_bins - z_bounds[i]))
      # Preserve endpoints and thresholds already assigned to this edge.
      if(closest %in% claimed){
        bins <- c(bins, z_bounds[i])
      }else{
        bins[closest] <- z_bounds[i]
        claimed <- c(claimed, closest)
      }
    }
    bins <- sort(unique(bins))
  }else if(type == "dens"){
    # for density, extend the specified support at the threshold
    bins <- sort(unique(c(bins, z_bounds)))
  }

  return(bins)
}
.zplot_bias_priors <- function(priors){

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
.zplot_threshold <- function(priors, effect_direction = "positive"){

  priors_bias <- .zplot_bias_priors(priors)

  # return simple binning in case of no bias
  if(length(priors_bias) == 0L){
    return(numeric())
  }

  p_bounds <- .selection_kernel_breaks(priors_bias)

  # return simple binning in case of no selection kernels
  if(is.null(p_bounds)){
    return(numeric())
  }

  # obtain tresholds from the specified priors
  z_bounds <- stats::qnorm(rev(p_bounds), 0, 1, lower.tail = FALSE)
  sign <- if (identical(effect_direction, "negative")) -1 else 1
  z_bounds <- sort(unique(sign * z_bounds[is.finite(z_bounds)]))

  return(z_bounds)
}
