.get_dots_hist_zplot <- function(dots, plot_type = "base", max_density = NULL) {
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
# .get_dots_lines_zplot
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
.get_dots_lines_zplot <- function(dots, plot_type = "base", col = "black", alpha = 0.40) {
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
# .zplot_bins
# ---------------------------------------------------------------------------- #
#
# Creates bin sequence for zplot histograms and densities.
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
.zplot_bins      <- function(priors, from, to, by, length.out, type = "hist"){

  if(is.null(length.out)){
    bins <- seq(from = from, to = to, by = by)
  }else{
    bins <- seq(from = from, to = to, length.out = length.out)
  }

  # obtain tresholds from the specified priors
  z_bounds <- .zplot_threshold(priors)

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
.zplot_threshold <- function(priors){

  priors_bias <- .zplot_bias_priors(priors)

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
