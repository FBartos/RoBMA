
### basic plotting functions ----
#' @title Plots brma Object
#'
#' @description \code{plot.brma} visualizes posterior
#' (and prior) distribution a brma object.
#'
#' @param x a fitted RoBMA object
#' @param parameter a parameter to be plotted. Defaults to
#' \code{"mu"} (for the effect size). The additional options
#' are \code{"tau"} (for the heterogeneity),
#' \code{"weightfunction"} (for the estimated weightfunction),
#' or \code{"PET-PEESE"} (for the PET-PEESE regression).
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' }
#'
#'
#' @return \code{plot.brma} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBMA()]
#' @export
plot.brma  <- function(
    x, parameter, parameter_mods, parameter_scale,
    prior = FALSE, plot_type = "base", dots_prior = NULL, ...) {


  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")

  ### select and validate the parameter to be plotted
  parameter <- .check_and_select_plot_parameter(
    parameter       = parameter,
    parameter_mods  = parameter_mods,
    parameter_scale = parameter_scale,
    object          = x
  )

  ### obtain posterior samples in the plotting format
  samples <- BayesTools::as_mixed_posteriors(
    model      = x[["fit"]],
    parameters = parameter
  )

  ### set up plotting arguments
  n_levels   <- .get_samples_n_levels(samples, parameter)
  dots       <- .set_dots_plot(..., n_levels = n_levels)
  dots_prior <- .set_dots_prior(dots_prior)

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- parameter
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$dots_prior               <- dots_prior
  args$individual               <- TRUE
  args$show_figures             <- NULL

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


#' @title Plots Weight Function of brma Object
#'
#' @description \code{plot.brma} visualizes posterior
#' (and prior) weight function of a brma object.
#'
#' @param x a fitted RoBMA object
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' }
#'
#'
#' @return \code{plot.brma} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBMA()]
#' @export
plot_weightfunction <- function(x, ...)  UseMethod("plot_weightfunction")

#' @export
#' @rdname plot_weightfunction
plot_weightfunction.brma  <- function(
    x, rescale_p_values = TRUE, show_data = TRUE,
    prior = FALSE, plot_type = "base", dots_prior = NULL, ...) {

  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(rescale_p_values, "rescale_p_values")
  BayesTools::check_bool(show_data, "show_data")

  if (!.is_weightfunction(x)) {
    stop("'plot_weightfunction' is available only for models with a weightfunction component.", call. = FALSE)
  }

  ### obtain posterior samples in the plotting format
  samples <- BayesTools::as_mixed_posteriors(
    model      = x[["fit"]],
    parameters = "omega"
  )

  ### set up plotting arguments
  dots       <- .set_dots_plot(...)
  dots_prior <- .set_dots_prior(dots_prior)

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- "weightfunction"
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$dots_prior               <- dots_prior
  args$individual               <- FALSE
  args$rescale_x                <- rescale_p_values

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}


#' @title Plots Weight Function of brma Object
#'
#' @description \code{plot.brma} visualizes posterior
#' (and prior) weight function of a brma object.
#'
#' @param x a fitted RoBMA object
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' }
#'
#'
#' @return \code{plot.brma} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object object of class 'ggplot2' if \code{plot_type = "ggplot2"}.
#'
#' @seealso [RoBMA()]
#' @export
plot_PETPEESE <- function(x, ...)  UseMethod("plot_PETPEESE")

#' @export
#' @rdname plot_PETPEESE
plot_PETPEESE.brma  <- function(
    x, show_data = TRUE,
    prior = FALSE, plot_type = "base", dots_prior = NULL, ...) {


  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(show_data, "show_data")

  if (!(.is_PET(x) || .is_PEESE(x))) {
    stop("'plot_PETPEESE' is available only for models with a PET or PEESE component.", call. = FALSE)
  }

  ### obtain posterior samples in the plotting format
  samples <- BayesTools::as_mixed_posteriors(
    model      = x[["fit"]],
    parameters = c("mu", if (.is_PET(x)) "PET", if (.is_PEESE(x)) "PEESE")
  )

  ### set up plotting arguments
  dots       <- .set_dots_plot(...)
  dots_prior <- .set_dots_prior(dots_prior)

  # set plotting range (make sure all effect sizes are within the plotting range)
  data_outcome <- x[["data"]][["outcome"]]
  if (is.null(dots[["ylim"]])) {
    dots[["ylim"]] <- range(pretty(c(0, range(data_outcome$yi))))
  }

  # prepare the argument call
  args                          <- dots
  args$samples                  <- samples
  args$parameter                <- "PETPEESE"
  args$plot_type                <- plot_type
  args$prior                    <- prior
  args$n_points                 <- 1000
  args$n_samples                <- 10000
  args$force_samples            <- FALSE
  args$dots_prior               <- dots_prior
  args$individual               <- FALSE
  args$effect_direction         <- .effect_direction(x)

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # include data
  if (show_data) {
    if (plot_type == "ggplot") {
      plot <- plot + ggplot2::geom_point(
        ggplot2::aes(
          x  = data_outcome$sei,
          y  = data_outcome$yi),
        size  = if(!is.null(dots[["cex"]])) dots[["cex"]] else 2,
        shape = if(!is.null(dots[["pch"]])) dots[["pch"]] else 16
      )
    } else if(plot_type == "base") {
      graphics::points(data_outcome$sei, data_outcome$yi, cex = if(!is.null(dots[["cex"]])) dots[["cex"]] else 1, pch = if(!is.null(dots[["pch"]])) dots[["pch"]] else 16)
    }
  }

  # return the plots
  if (plot_type == "base") {
    return(invisible(plot))
  } else if (plot_type == "ggplot") {
    return(plot)
  }
}



.clean_PET_PEESE_samples <- function(samples, parameter) {
  for (i in seq_along(attr(samples[["PET"]], "prior_list"))) {
    if (is.prior.PET(attr(samples[["PET"]], "prior_list")[[i]])) {
      class(attr(samples[["PET"]], "prior_list")[[i]])   <- class(attr(samples[["PET"]], "prior_list")[[i]])[!class(attr(samples[["PET"]], "prior_list")[[i]]) %in% "prior.PET"]
    } else if (!is.prior.point(attr(samples[["PET"]], "prior_list")[[i]])){
      attr(samples[["PET"]], "prior_list")[[i]] <- prior("point", parameter = list("location" = 0), prior_weights = attr(samples[["PET"]], "prior_list")[[i]][["prior_weights"]])
    }
  }
  for (i in seq_along(attr(samples[["PEESE"]], "prior_list"))) {
    if (is.prior.PEESE(attr(samples[["PEESE"]], "prior_list")[[i]])) {
      class(attr(samples[["PEESE"]], "prior_list")[[i]])   <- class(attr(samples[["PEESE"]], "prior_list")[[i]])[!class(attr(samples[["PEESE"]], "prior_list")[[i]]) %in% "prior.PEESE"]
    } else if (!is.prior.point(attr(samples[["PEESE"]], "prior_list")[[i]])) {
      attr(samples[["PEESE"]], "prior_list")[[i]] <- prior("point", parameter = list("location" = 0), prior_weights = attr(samples[["PEESE"]], "prior_list")[[i]][["prior_weights"]])
    }
  }
  return(samples)
}




### basic MCMC diagnostics functions ----

#' @export
diagnostic_plots <- function(x, ...) UseMethod("diagnostic_plots")

#' @export
#' @rdname diagnostic_plots
diagnostic_plots.brma <- function(x, parameter, parameter_mods, parameter_scale, type, plot_type = "base", lags = 30, ...){

  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_char(type, "type", allow_values = c("autocorrelation", "trace", "density"))

  ### select and validate the parameter to be plotted
  parameter <- .check_and_select_plot_parameter(
    parameter       = parameter,
    parameter_mods  = parameter_mods,
    parameter_scale = parameter_scale,
    object          = x
  )

  # a message with info about multiple plots
  if(plot_type == "base" && parameter == "omega") {
    message("Multiple plots will be produced. See '?layout' for help with setting multiple plots.")
  }

  dots  <- .set_dots_diagnostics(..., type = type, chains = x[["fit_control"]][["chains"]])

  # get the parameter name
  args                   <- dots
  args$fit               <- x[["fit"]]
  args$parameter         <- parameter
  args$type              <- type
  args$plot_type         <- plot_type
  args$lags              <- lags
  args$transform_factors <- TRUE
  args$short_name        <- FALSE
  args$parameter_names   <- FALSE
  args$formula_prefix    <- FALSE

  plots <- do.call(BayesTools::JAGS_diagnostics, args)

  # return the plots
  if(plot_type == "base"){
    return(invisible(plots))
  }else if(plot_type == "ggplot"){
    return(plots)
  }
}

#' @export
#' @rdname diagnostic_plots
diagnostic_plots_autocorrelation <- function(x, ...) UseMethod("diagnostic_plots_autocorrelation")

#' @export
#' @rdname diagnostic_plots
diagnostic_plots_autocorrelation.brma <- function(x, parameter = NULL, plot_type = "base", lags = 30, ...) {
  diagnostic_plots(x = x, parameter = parameter, type = "autocorrelation", plot_type = plot_type, lags = lags, ...)
}

#' @export
#' @rdname diagnostic_plots
diagnostic_plots_trace <- function(x, ...) UseMethod("diagnostic_plots_trace")

#' @export
#' @rdname diagnostic_plots
diagnostic_plots_trace.brma           <- function(x, parameter = NULL, plot_type = "base", ...) {
  diagnostic_plots(x = x, parameter = parameter, type = "trace", plot_type = plot_type, ...)
}

#' @export
#' @rdname diagnostic_plots
diagnostic_plots_density <- function(x, ...) UseMethod("diagnostic_plots_density")

#' @export
#' @rdname diagnostic_plots
diagnostic_plots_density.brma         <- function(x, parameter = NULL, plot_type = "base", ...) {
  diagnostic_plots(x = x, parameter = parameter, type = "density", plot_type = plot_type, ...)
}


### Common plotting helper functions ----
# Internal helper function for selecting and validating a plot parameter
#
# Processes parameter selection for plotting based on three mutually exclusive arguments:
# parameter, parameter_mods, parameter_scale. Only one should be specified.
#
# @param parameter Base parameter name ("mu", "tau", "rho")
# @param parameter_mods Moderator parameter for location regression
# @param parameter_scale Moderator parameter for scale regression
# @param object A brma object
#
# @return A character string with the validated parameter name (with appropriate prefix).
#   Attributes:
#   - parameter_samples: The JAGS samples name
#   - type: One of "base", "mods", "scale" indicating the parameter source
.check_and_select_plot_parameter <- function(parameter, parameter_mods, parameter_scale, object) {

  # Check for model characteristics
  is_mods           <- .is_mods(object)
  is_scale          <- .is_scale(object)
  is_multilevel     <- .is_multilevel(object)
  is_PET            <- .is_PET(object)
  is_PEESE          <- .is_PEESE(object)
  is_weightfunction <- .is_weightfunction(object)

  # Determine which arguments were provided
  has_parameter       <- !missing(parameter)       && !is.null(parameter)
  has_parameter_mods  <- !missing(parameter_mods)  && !is.null(parameter_mods)
  has_parameter_scale <- !missing(parameter_scale) && !is.null(parameter_scale)

  # Count how many were specified

  n_specified <- sum(c(has_parameter, has_parameter_mods, has_parameter_scale))

  # Error if more than one specified
  if (n_specified > 1) {
    stop("Only one of 'parameter', 'parameter_mods', or 'parameter_scale' can be specified.", call. = FALSE)
  }

  # Default behavior if none specified
  if (n_specified == 0) {
    if (is_mods) {
      # Default to intercept for meta-regression models
      return("mu_intercept")
    } else {
      # Default to mu for simple models
      return("mu")
    }
  }

  # Process based on which argument was specified
  if (has_parameter) {

    # Validate base parameters
    BayesTools::check_char(parameter, "parameter")

    # Handle mu parameter
    if (parameter == "mu") {
      if (is_mods) {
        # mu becomes mu_intercept when mods are present
        return("mu_intercept")
      } else {
        return("mu")
      }
    }

    # Handle tau parameter
    if (parameter == "tau") {
      if (is_scale) {
        # tau becomes log_tau_intercept when scale is present
        return("log_tau_intercept")
      } else {
        return("tau")
      }
    }

    # Handle rho parameter (only for multilevel models)
    if (parameter == "rho") {
      if (!is_multilevel) {
        stop("The 'rho' parameter is only available for multilevel models.", call. = FALSE)
      }
      return("rho")
    }

    ### Handle publication bias parameters
    if (parameter == "PET") {
      if (!is_PET) {
        stop("The 'PET' parameter is only available for PET models.", call. = FALSE)
      }
      return("PET")
    }

    if (parameter == "PEESE") {
      if (!is_PEESE) {
        stop("The 'PEESE' parameter is only available for PEESE models.", call. = FALSE)
      }
      return("PEESE")
    }

    if (parameter == "omega" || parameter == "weightfunction") {
      if (!is_weightfunction) {
        stop("The 'omega'/'weightfunction' parameter is only available for selection models.", call. = FALSE)
      }
      return("omega")
    }

    # Unknown base parameter
    # Build list of valid parameters dynamically
    valid_params <- c("'mu'", "'tau'")
    if (is_multilevel)     valid_params <- c(valid_params, "'rho'")
    if (is_PET)            valid_params <- c(valid_params, "'PET'")
    if (is_PEESE)          valid_params <- c(valid_params, "'PEESE'")
    if (is_weightfunction) valid_params <- c(valid_params, "'omega'/'weightfunction'")

    stop(sprintf(
      "Unknown parameter '%s'. Valid base parameters are: %s.",
      parameter,
      paste(valid_params, collapse = ", ")
    ), call. = FALSE)

  } else if (has_parameter_mods) {

    # Validate that mods are present
    if (!is_mods) {
      stop("The 'parameter_mods' argument can only be used for models with moderators.", call. = FALSE)
    }

    BayesTools::check_char(parameter_mods, "parameter_mods")

    # Get available terms from the mods formula
    mods_formula    <- attr(object[["data"]][["mods"]], "formula")
    available_terms <- c("intercept", attr(stats::terms(mods_formula), "term.labels"))

    # Validate the specified term exists
    if (!parameter_mods %in% available_terms) {
      stop(sprintf(
        "The specified 'parameter_mods' term '%s' is not available. Available terms are: %s.",
        parameter_mods,
        paste0("'", available_terms, "'", collapse = ", ")
      ), call. = FALSE)
    }

    # Prefix with "mu_"
    return(paste0("mu_", parameter_mods))

  } else if (has_parameter_scale) {

    # Validate that scale is present
    if (!is_scale) {
      stop("The 'parameter_scale' argument can only be used for location-scale models.", call. = FALSE)
    }

    BayesTools::check_char(parameter_scale, "parameter_scale")

    # Get available terms from the scale formula
    scale_formula <- attr(object[["data"]][["scale"]], "formula")
    available_terms <- c("intercept", attr(stats::terms(scale_formula), "term.labels"))

    # Validate the specified term exists
    if (!parameter_scale %in% available_terms) {
      stop(sprintf(
        "The specified 'parameter_scale' term '%s' is not available. Available terms are: %s.",
        parameter_scale,
        paste0("'", available_terms, "'", collapse = ", ")
      ), call. = FALSE)
    }

    # Prefix with "log_tau_"
    return(paste0("log_tau_", parameter_scale))
  }
}

.set_dots_plot        <- function(..., n_levels = 1) {

  dots <- list(...)
  if (is.null(dots[["col"]]) & n_levels == 1) {
    dots[["col"]]      <- "black"
  }else if (is.null(dots[["col"]]) & n_levels > 1) {
    dots[["col"]]      <- grDevices::palette.colors(n = n_levels + 1, palette = "Okabe-Ito")[-1]
  }
  if (is.null(dots[["col.fill"]])) {
    dots[["col.fill"]] <- "#4D4D4D4C" # scales::alpha("grey30", .30)
  }

  return(dots)
}
.set_dots_prior       <- function(dots_prior) {

  if (is.null(dots_prior)) {
    dots_prior <- list()
  }

  if (is.null(dots_prior[["col"]])) {
    dots_prior[["col"]]      <- "grey60"
  }
  if (is.null(dots_prior[["lty"]])) {
    dots_prior[["lty"]]      <- 1
  }
  if (is.null(dots_prior[["col.fill"]])) {
    dots_prior[["col.fill"]] <- "#B3B3B34C" # scales::alpha("grey70", .30)
  }

  return(dots_prior)
}
.set_dots_diagnostics <- function(..., type, chains) {

  dots <- list(...)
  if (is.null(dots[["col"]])) {
    dots[["col"]] <- if(type == "autocorrelation") "black" else rev(scales::viridis_pal()(chains))
  }

  return(dots)
}
.get_samples_n_levels <- function(samples, parameter) {

  if (inherits(samples[[parameter]], "mixed_posteriors.factor")) {
    if (attr(samples[[parameter]], "orthonormal") || attr(samples[[parameter]], "meandif")) {
      n_levels <- length(attr(samples[[parameter]], "level_names"))
    } else if (attr(samples[[parameter]], "treatment")) {
      n_levels <- length(attr(samples[[parameter]], "level_names")) - 1
    } else {
      n_levels <- 1
    }
  } else {
    n_levels <- 1
  }

  return(n_levels)
}
