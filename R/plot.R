
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
    prior = FALSE, standardized_coefficients = FALSE, plot_type = "base", dots_prior = NULL, ...) {


  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(standardized_coefficients, "standardized_coefficients")

  ### select and validate the parameter to be plotted
  parameter <- .check_and_select_plot_parameter(
    parameter       = parameter,
    parameter_mods  = parameter_mods,
    parameter_scale = parameter_scale,
    object          = x
  )

  ### obtain posterior samples in the plotting format
  samples <- BayesTools::as_mixed_posteriors(
    model            = x[["fit"]],
    parameters       = parameter,
    transform_scaled = !standardized_coefficients
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


#' @title Plot Prior Distributions
#'
#' @description \code{plot_prior} visualizes prior distributions stored in
#' \code{brma}, \code{BMA}, and \code{RoBMA} objects.
#' This is especially useful for objects created with \code{only_priors = TRUE}.
#'
#' @param x a \code{brma}, \code{BMA}, or \code{RoBMA} object, or a prior
#' distribution object.
#' @param parameter character. Base parameter to plot. Defaults to \code{"mu"}.
#' Common options are \code{"mu"}, \code{"tau"}, \code{"rho"}, \code{"PET"},
#' \code{"PEESE"}, \code{"omega"}, and \code{"bias"}. Moderator and scale
#' terms can also be selected by name when unambiguous.
#' @param parameter_mods character. Moderator term to plot.
#' Use \code{"intercept"} for the adjusted effect in meta-regression models.
#' @param parameter_scale character. Scale-regression term to plot.
#' Use \code{"intercept"} for the heterogeneity intercept in location-scale models.
#' @param standardized_coefficients whether to plot moderator and scale-regression
#' priors on the standardized predictor scale. Defaults to \code{TRUE}, which
#' shows the priors as specified. Set to \code{FALSE} to transform them to the
#' original predictor scale when continuous predictors were standardized.
#' @param plot_type whether to use a base plot \code{"base"} or ggplot2
#' \code{"ggplot"} for plotting. Defaults to \code{"base"}.
#' @param ... additional arguments passed to the prior plotting method.
#'
#' @return \code{plot_prior} returns either \code{NULL} invisibly if
#' \code{plot_type = "base"} or a \code{ggplot2} object if
#' \code{plot_type = "ggplot"}. If multiple parameters are requested with
#' \code{plot_type = "ggplot"}, a named list of plots is returned.
#'
#' @examples \dontrun{
#' priors <- BMA(yi = yi, sei = sei, data = dat, measure = "SMD", only_priors = TRUE)
#' plot_prior(priors, parameter = "mu")
#' plot_prior(priors, parameter = "tau")
#' plot_prior(priors, parameter_mods = "x", standardized_coefficients = FALSE)
#' }
#'
#' @seealso [BMA()] [RoBMA()] [brma()] [prior()]
#' @export
plot_prior <- function(x, ...) UseMethod("plot_prior")

#' @export
#' @rdname plot_prior
plot_prior.prior <- function(x, plot_type = "base", ...) {

  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))

  plot <- plot(x, plot_type = plot_type, ...)

  if (plot_type == "base") {
    return(invisible(plot))
  } else if (plot_type == "ggplot") {
    return(plot)
  }
}

#' @export
#' @rdname plot_prior
plot_prior.brma <- function(
    x, parameter = "mu", parameter_mods, parameter_scale,
    standardized_coefficients = TRUE, plot_type = "base", ...) {

  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(standardized_coefficients, "standardized_coefficients")

  if (is.null(x[["priors"]])) {
    stop("The object does not contain prior distributions.", call. = FALSE)
  }

  has_parameter       <- !missing(parameter)       && !is.null(parameter)
  has_parameter_mods  <- !missing(parameter_mods)  && !is.null(parameter_mods)
  has_parameter_scale <- !missing(parameter_scale) && !is.null(parameter_scale)

  n_specified <- sum(c(has_parameter, has_parameter_mods, has_parameter_scale))
  if (n_specified > 1) {
    stop("Only one of 'parameter', 'parameter_mods', or 'parameter_scale' can be specified.", call. = FALSE)
  }

  if (has_parameter && length(parameter) > 1) {
    if (!is.character(parameter) || anyNA(parameter) || any(!nzchar(parameter))) {
      stop("The 'parameter' argument must be a character vector without missing or empty values.", call. = FALSE)
    }

    plots <- lapply(parameter, function(parameter_i) {
      plot_prior(
        x                         = x,
        parameter                 = parameter_i,
        standardized_coefficients = standardized_coefficients,
        plot_type                 = plot_type,
        ...
      )
    })
    names(plots) <- parameter

    if (plot_type == "base") {
      return(invisible(plots))
    } else if (plot_type == "ggplot") {
      return(plots)
    }
  }

  selected <- .select_plot_prior_parameter(
    object          = x,
    parameter       = if (has_parameter) parameter else NULL,
    parameter_mods  = if (has_parameter_mods) parameter_mods else NULL,
    parameter_scale = if (has_parameter_scale) parameter_scale else NULL
  )

  dots <- list(...)
  if (is.null(dots[["par_name"]])) {
    dots[["par_name"]] <- selected[["label"]]
  }

  if (!standardized_coefficients) {
    prior_plot <- .plot_prior_unstandardized(
      object    = x,
      selected  = selected,
      plot_type = plot_type,
      dots      = dots
    )
    if (!is.null(prior_plot)) {
      if (plot_type == "base") {
        return(invisible(prior_plot))
      } else if (plot_type == "ggplot") {
        return(prior_plot)
      }
    }
  }

  if (.use_plot_prior_list_dispatch(selected[["prior"]])) {
    plot <- do.call(
      BayesTools::plot_prior_list,
      c(
        list(prior_list = list(selected[["prior"]]), plot_type = plot_type),
        dots
      )
    )
  } else {
    plot <- do.call(
      "plot",
      c(
        list(x = selected[["prior"]], plot_type = plot_type),
        dots
      )
    )
  }

  if (plot_type == "base") {
    return(invisible(plot))
  } else if (plot_type == "ggplot") {
    return(plot)
  }
}


#' @title Print Prior Distributions
#'
#' @description \code{print_prior} prints prior distributions stored in
#' \code{brma}, \code{BMA}, and \code{RoBMA} objects.
#' This is a user-facing helper for inspecting priors without extracting them
#' from the internal \code{$priors} list.
#'
#' @param x a \code{brma}, \code{BMA}, or \code{RoBMA} object, or a prior
#' distribution object.
#' @param parameter character. Base parameter to print. Defaults to \code{"mu"}.
#' Common options are \code{"mu"}, \code{"tau"}, \code{"rho"}, \code{"PET"},
#' \code{"PEESE"}, \code{"omega"}, and \code{"bias"}. Moderator and scale
#' terms can also be selected by name when unambiguous.
#' @param parameter_mods character. Moderator term to print.
#' Use \code{"intercept"} for the adjusted effect in meta-regression models.
#' @param parameter_scale character. Scale-regression term to print.
#' Use \code{"intercept"} for the heterogeneity intercept in location-scale models.
#' @param ... additional arguments passed to the prior printing method.
#'
#' @return \code{print_prior} invisibly returns the selected prior distribution.
#' If multiple parameters are requested, a named list of prior distributions is
#' returned invisibly.
#'
#' @examples \dontrun{
#' priors <- BMA(yi = yi, sei = sei, data = dat, measure = "SMD", only_priors = TRUE)
#' print_prior(priors, parameter = "mu")
#' print_prior(priors, parameter = "tau")
#' print_prior(priors, parameter_mods = "x")
#' }
#'
#' @seealso [plot_prior()] [BMA()] [RoBMA()] [brma()] [prior()]
#' @export
print_prior <- function(x, ...) UseMethod("print_prior")

#' @export
#' @rdname print_prior
print_prior.prior <- function(x, ...) {

  .print_prior_object(x, ...)

  return(invisible(x))
}

#' @export
#' @rdname print_prior
print_prior.brma <- function(
    x, parameter = "mu", parameter_mods, parameter_scale, ...) {

  if (is.null(x[["priors"]])) {
    stop("The object does not contain prior distributions.", call. = FALSE)
  }

  has_parameter       <- !missing(parameter)       && !is.null(parameter)
  has_parameter_mods  <- !missing(parameter_mods)  && !is.null(parameter_mods)
  has_parameter_scale <- !missing(parameter_scale) && !is.null(parameter_scale)

  n_specified <- sum(c(has_parameter, has_parameter_mods, has_parameter_scale))
  if (n_specified > 1) {
    stop("Only one of 'parameter', 'parameter_mods', or 'parameter_scale' can be specified.", call. = FALSE)
  }

  if (has_parameter && length(parameter) > 1) {
    if (!is.character(parameter) || anyNA(parameter) || any(!nzchar(parameter))) {
      stop("The 'parameter' argument must be a character vector without missing or empty values.", call. = FALSE)
    }

    selected <- lapply(parameter, function(parameter_i) {
      .select_plot_prior_parameter(
        object           = x,
        parameter        = parameter_i,
        parameter_mods   = NULL,
        parameter_scale  = NULL,
        allow_mixed_bias = TRUE
      )
    })
    names(selected) <- parameter

    priors <- lapply(selected, `[[`, "prior")
    for (i in seq_along(selected)) {
      .print_prior_object(priors[[i]], label = selected[[i]][["label"]], ...)
    }

    return(invisible(priors))
  }

  selected <- .select_plot_prior_parameter(
    object           = x,
    parameter        = if (has_parameter) parameter else NULL,
    parameter_mods   = if (has_parameter_mods) parameter_mods else NULL,
    parameter_scale  = if (has_parameter_scale) parameter_scale else NULL,
    allow_mixed_bias = TRUE
  )

  .print_prior_object(selected[["prior"]], label = selected[["label"]], ...)

  return(invisible(selected[["prior"]]))
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


#' @title Plots PET-PEESE Fit of brma Object
#'
#' @description \code{plot.brma} visualizes posterior
#' (and prior) PET-PEESE fit of a brma object.
#'
#' @param x a fitted RoBMA object
#' @param show_data whether the observed effect sizes should be added to
#' figure. Defaults to \code{TRUE}.
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

    return(BayesTools::JAGS_parameter_names(
      parameters        = parameter_mods,
      formula_parameter = "mu"
    ))

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

    return(BayesTools::JAGS_parameter_names(
      parameters        = parameter_scale,
      formula_parameter = "log_tau"
    ))
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

.select_plot_prior_parameter <- function(
    object, parameter, parameter_mods, parameter_scale, allow_mixed_bias = FALSE) {

  n_specified <- sum(!vapply(
    list(parameter, parameter_mods, parameter_scale),
    is.null,
    logical(1)
  ))

  if (n_specified == 0) {
    parameter <- "mu"
  }

  priors <- object[["priors"]]

  if (!is.null(parameter_mods)) {
    .check_plot_prior_name(parameter_mods, "parameter_mods")
    return(.select_plot_prior_term(
      prior_list = priors[["mods"]],
      term       = parameter_mods,
      argument   = "parameter_mods",
      prefix     = "mu",
      source     = "mods"
    ))
  }

  if (!is.null(parameter_scale)) {
    .check_plot_prior_name(parameter_scale, "parameter_scale")
    return(.select_plot_prior_term(
      prior_list = priors[["scale"]],
      term       = parameter_scale,
      argument   = "parameter_scale",
      prefix     = "log_tau",
      source     = "scale"
    ))
  }

  .check_plot_prior_name(parameter, "parameter")

  parameter <- switch(
    parameter,
    "effect"        = "mu",
    "heterogeneity" = "tau",
    "weightfunction" = "omega",
    parameter
  )

  if (parameter == "mu" && is.null(priors[["outcome"]][["mu"]]) && !is.null(priors[["mods"]][["intercept"]])) {
    return(list(
      prior  = priors[["mods"]][["intercept"]],
      label  = "mu_intercept",
      source = "mods",
      term   = "intercept"
    ))
  }

  if (parameter == "tau" && is.null(priors[["outcome"]][["tau"]]) && !is.null(priors[["scale"]][["intercept"]])) {
    return(list(
      prior  = priors[["scale"]][["intercept"]],
      label  = "log_tau_intercept",
      source = "scale",
      term   = "intercept"
    ))
  }

  if (parameter %in% c("omega", "PET", "PEESE", "bias") && !is.null(priors[["outcome"]][["bias"]])) {
    prior <- .select_plot_prior_bias(
      prior            = priors[["outcome"]][["bias"]],
      parameter        = parameter,
      allow_mixed_bias = allow_mixed_bias
    )
    return(list(prior = prior, label = parameter, source = "bias", term = parameter))
  }

  if (parameter %in% names(priors[["outcome"]]) && !is.null(priors[["outcome"]][[parameter]])) {
    return(list(prior = priors[["outcome"]][[parameter]], label = parameter, source = "outcome", term = parameter))
  }

  in_mods  <- !is.null(priors[["mods"]])  && parameter %in% names(priors[["mods"]])
  in_scale <- !is.null(priors[["scale"]]) && parameter %in% names(priors[["scale"]])

  if (in_mods && in_scale) {
    stop(sprintf(
      "The term '%s' is available in both moderator and scale priors. Use 'parameter_mods' or 'parameter_scale'.",
      parameter
    ), call. = FALSE)
  }

  if (in_mods) {
    return(list(
      prior  = priors[["mods"]][[parameter]],
      label  = paste0("mu_", parameter),
      source = "mods",
      term   = parameter
    ))
  }

  if (in_scale) {
    return(list(
      prior  = priors[["scale"]][[parameter]],
      label  = paste0("log_tau_", parameter),
      source = "scale",
      term   = parameter
    ))
  }

  stop(sprintf(
    "Unknown prior parameter '%s'. Available parameters are: %s.",
    parameter,
    .collapse_plot_prior_names(priors)
  ), call. = FALSE)
}

.select_plot_prior_term <- function(prior_list, term, argument, prefix, source) {

  if (is.null(prior_list)) {
    stop(sprintf("The '%s' argument can only be used when the object contains the corresponding priors.", argument), call. = FALSE)
  }

  if (!term %in% names(prior_list)) {
    stop(sprintf(
      "The specified '%s' term '%s' is not available. Available terms are: %s.",
      argument,
      term,
      paste0("'", names(prior_list), "'", collapse = ", ")
    ), call. = FALSE)
  }

  return(list(
    prior  = prior_list[[term]],
    label  = paste0(prefix, "_", term),
    source = source,
    term   = term
  ))
}

.select_plot_prior_bias <- function(prior, parameter, allow_mixed_bias = FALSE) {

  if (!BayesTools::is.prior.mixture(prior)) {
    if (parameter == "bias" ||
        (parameter == "omega" && BayesTools::is.prior.weightfunction(prior)) ||
        (parameter == "PET" && BayesTools::is.prior.PET(prior)) ||
        (parameter == "PEESE" && BayesTools::is.prior.PEESE(prior))) {
      return(prior)
    }

    stop(sprintf("The publication-bias prior does not contain a '%s' component.", parameter), call. = FALSE)
  }

  has_weightfunction <- any(vapply(prior, BayesTools::is.prior.weightfunction, logical(1)))
  has_petpeese       <- any(vapply(prior, function(x) BayesTools::is.prior.PET(x) || BayesTools::is.prior.PEESE(x), logical(1)))

  if (parameter == "bias") {
    if (has_weightfunction && has_petpeese) {
      if (!allow_mixed_bias) {
        stop(
          "The publication-bias prior mixes weight-function and PET-PEESE components. ",
          "Use parameter = 'omega', 'PET', or 'PEESE' to plot one component type.",
          call. = FALSE
        )
      }
    }
    return(prior)
  }

  keep <- switch(
    parameter,
    "omega" = vapply(prior, BayesTools::is.prior.weightfunction, logical(1)),
    "PET"   = vapply(prior, BayesTools::is.prior.PET, logical(1)),
    "PEESE" = vapply(prior, BayesTools::is.prior.PEESE, logical(1))
  )

  if (!any(keep)) {
    stop(sprintf("The publication-bias prior does not contain a '%s' component.", parameter), call. = FALSE)
  }

  selected <- unclass(prior)[keep]
  if (length(selected) == 1) {
    return(selected[[1]])
  }

  class(selected) <- c("prior", "prior.mixture")
  attr(selected, "components")    <- rep("alternative", length(selected))
  attr(selected, "prior_weights") <- rep(1, length(selected))

  return(selected)
}

.print_prior_object <- function(x, label = NULL, ...) {

  dots <- list(...)
  silent <- isTRUE(dots[["silent"]])
  dots[["silent"]] <- TRUE

  output <- do.call(print, c(list(x = x), dots))

  if (!silent) {
    if (!is.null(label)) {
      cat(label, ": ", sep = "")
    }
    cat(output, "\n", sep = "")
  }

  return(invisible(x))
}

.plot_prior_unstandardized <- function(object, selected, plot_type, dots) {

  if (!(selected[["source"]] %in% c("mods", "scale"))) {
    return(NULL)
  }

  if (!isTRUE(.standardize_continuous_predictors(object))) {
    return(NULL)
  }

  formula_info <- .plot_prior_formula_info(
    object = object,
    source = selected[["source"]]
  )

  if (is.null(formula_info) || is.null(formula_info[["formula_scale"]])) {
    return(NULL)
  }

  parameter <- gsub(":", "__xXx__", selected[["label"]], fixed = TRUE)
  if (!parameter %in% formula_info[["column_names"]]) {
    return(NULL)
  }

  n_points <- if (!is.null(dots[["n_points"]])) dots[["n_points"]] else 1000
  BayesTools::check_int(n_points, "n_points", lower = 2)

  densities <- BayesTools:::.generate_transformed_prior_densities(
    prior_list    = formula_info[["prior_list"]],
    column_names  = formula_info[["column_names"]],
    formula_scale = formula_info[["formula_scale"]],
    n_grid        = max(n_points, 16)
  )

  if (!parameter %in% names(densities)) {
    return(NULL)
  }

  density <- densities[[parameter]]
  if (.plot_prior_identity_transform(density, parameter)) {
    return(NULL)
  }

  plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    x                         = density,
    n_points                  = n_points,
    x_range                   = dots[["xlim"]],
    transformation            = dots[["transformation"]],
    transformation_arguments  = dots[["transformation_arguments"]],
    transformation_settings   = if (!is.null(dots[["transformation_settings"]])) dots[["transformation_settings"]] else FALSE
  )

  par_name <- dots[["par_name"]]
  dots[c(
    "par_name", "n_points", "n_samples", "force_samples", "x_range_quant",
    "individual", "show_figures", "rescale_x", "transformation",
    "transformation_arguments", "transformation_settings"
  )] <- NULL

  plot <- do.call(
    BayesTools:::.plot_prior_list.both,
    c(
      list(
        plot_data = plot_data,
        plot_type = plot_type,
        par_name  = par_name
      ),
      dots
    )
  )

  return(plot)
}

.plot_prior_formula_info <- function(object, source) {

  parameter <- switch(
    source,
    "mods"  = "mu",
    "scale" = "log_tau"
  )

  formula_output <- BayesTools::JAGS_formula(
    formula       = .create_fit_formula_list(data = object[["data"]], parameter = source),
    parameter     = parameter,
    data          = .create_fit_formula_data_list(data = object[["data"]], parameter = source),
    prior_list    = .create_fit_formula_prior_list(priors = object[["priors"]], parameter = source),
    formula_scale = .data_standardize_continuous_predictors(object[["data"]])
  )

  if (is.null(formula_output[["formula_scale"]])) {
    return(NULL)
  }

  formula_scale <- list(formula_output[["formula_scale"]])
  names(formula_scale) <- parameter

  return(list(
    prior_list    = formula_output[["prior_list"]],
    column_names  = names(formula_output[["prior_list"]]),
    formula_scale = formula_scale
  ))
}

.plot_prior_identity_transform <- function(density, parameter) {

  weights <- attr(density, "weights")

  if (is.null(weights) || length(weights) != 1) {
    return(FALSE)
  }

  if (!identical(names(weights), parameter)) {
    return(FALSE)
  }

  return(isTRUE(all.equal(unname(weights), 1)))
}

.use_plot_prior_list_dispatch <- function(prior) {

  return(
    BayesTools::is.prior.factor(prior) &&
      !BayesTools::is.prior.mixture(prior) &&
      (BayesTools::is.prior.meandif(prior) || BayesTools::is.prior.orthonormal(prior))
  )
}

.check_plot_prior_name <- function(x, argument) {

  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop(sprintf("The '%s' argument must be a single non-empty character string.", argument), call. = FALSE)
  }

  return(invisible(TRUE))
}

.collapse_plot_prior_names <- function(priors) {

  names_available <- names(priors[["outcome"]])

  if (!is.null(priors[["outcome"]][["bias"]])) {
    names_available <- unique(c(names_available, "omega", "PET", "PEESE"))
  }

  if (!is.null(priors[["mods"]])) {
    names_available <- c(names_available, paste0("parameter_mods = '", names(priors[["mods"]]), "'"))
  }

  if (!is.null(priors[["scale"]])) {
    names_available <- c(names_available, paste0("parameter_scale = '", names(priors[["scale"]]), "'"))
  }

  return(paste0("'", names_available, "'", collapse = ", "))
}
