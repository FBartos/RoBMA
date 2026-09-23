### basic MCMC diagnostics functions ----

#' @title Plot MCMC Diagnostics
#'
#' @description \code{plot_diagnostic} creates visual MCMC diagnostics for a
#' fitted brma object. Convenience wrappers are available for trace, density,
#' and autocorrelation plots.
#'
#' @param x a fitted brma object
#' @param parameter base parameter to plot. Defaults to \code{NULL}, which uses
#'   \code{"mu"} or the meta-regression intercept. Valid values include
#'   \code{"mu"}, \code{"tau"}, \code{"rho"}, \code{"PET"}, \code{"PEESE"},
#'   and \code{"omega"} or \code{"weightfunction"} when present.
#' @param parameter_mods legacy moderator selector. Prefer \code{parameter}
#'   with \code{component = "mods"}.
#' @param parameter_scale legacy scale-regression selector. Prefer
#'   \code{parameter} with \code{component = "scale"}.
#' @param component parameter component. Defaults to \code{"auto"}, which
#'   infers the component when possible. Use \code{"mods"} (alias
#'   \code{"location"}), \code{"scale"}, \code{"random"}, or \code{"bias"} to
#'   disambiguate terms used in multiple model components. The random component
#'   selects semantic standard deviation, correlation, and allocation parameters
#'   from `brma.mv()` random formulas.
#' @param type diagnostic plot type. Convenience wrappers set a type-specific
#'   default but still forward this argument to \code{plot_diagnostic.brma()}.
#' @param plot_type whether to use a base plot \code{"base"} or ggplot2
#'   \code{"ggplot"} for plotting. Defaults to \code{"base"}.
#' @param lags number of lags for autocorrelation plots. Defaults to 30.
#' @param ... additional graphical arguments passed through RoBMA's diagnostic
#'   setup to \code{BayesTools::JAGS_diagnostics()}.
#'
#' @return \code{plot_diagnostic} returns the object returned by
#'   \code{BayesTools::JAGS_diagnostics()}, invisibly for base graphics.
#'
#' @seealso [summary.brma()]
#'
#' @export
plot_diagnostic <- function(x, ...) UseMethod("plot_diagnostic")

#' @export
#' @rdname plot_diagnostic
plot_diagnostic.brma <- function(
    x, parameter = NULL, parameter_mods = NULL, parameter_scale = NULL,
    type, plot_type = "base", lags = 30, component = "auto", ...) {

  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_char(type, "type", allow_values = c("autocorrelation", "trace", "density"))

  ### select and validate the parameter to be plotted
  parameter <- .check_and_select_plot_parameter(
    parameter        = parameter,
    parameter_mods   = parameter_mods,
    parameter_scale  = parameter_scale,
    component        = component,
    object           = x
  )
  parameter_entry <- .brma_parameter_select_entry(x, parameter)
  is_random       <- identical(parameter_entry[["component"]], "random")

  # a message with info about multiple plots
  if(plot_type == "base" && parameter == "omega") {
    message("Multiple plots will be produced. See '?layout' for help with setting multiple plots.")
  }

  dots  <- .set_dots_diagnostics(..., type = type, chains = x[["fit_control"]][["chains"]])

  # get the parameter name
  if (is_random) {
    random_diagnostic <- .brma_random_parameter_diagnostic_fit(
      object    = x,
      parameter = parameter
    )
    diagnostic_fit       <- random_diagnostic[["fit"]]
    diagnostic_parameter <- random_diagnostic[["parameter"]]
    if (is.null(dots[["main"]])) {
      dots[["main"]] <- random_diagnostic[["label"]]
    }
  } else {
    diagnostic_fit       <- x[["fit"]]
    diagnostic_parameter <- parameter
  }

  args                   <- dots
  args$fit               <- diagnostic_fit
  args$parameter         <- diagnostic_parameter
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
#' @rdname plot_diagnostic
plot_diagnostic_autocorrelation <- function(x, ...) UseMethod("plot_diagnostic_autocorrelation")

#' @export
#' @rdname plot_diagnostic
plot_diagnostic_autocorrelation.brma <- function(
    x, parameter = NULL, parameter_mods = NULL, parameter_scale = NULL,
    type = "autocorrelation", plot_type = "base", lags = 30,
    component = "auto", ...) {

  plot_diagnostic(
    x               = x,
    parameter       = parameter,
    parameter_mods   = parameter_mods,
    parameter_scale  = parameter_scale,
    component        = component,
    type             = type,
    plot_type        = plot_type,
    lags             = lags,
    ...
  )
}

#' @export
#' @rdname plot_diagnostic
plot_diagnostic_trace <- function(x, ...) UseMethod("plot_diagnostic_trace")

#' @export
#' @rdname plot_diagnostic
plot_diagnostic_trace.brma           <- function(
    x, parameter = NULL, parameter_mods = NULL, parameter_scale = NULL,
    type = "trace", plot_type = "base", lags = 30,
    component = "auto", ...) {

  plot_diagnostic(
    x               = x,
    parameter       = parameter,
    parameter_mods   = parameter_mods,
    parameter_scale  = parameter_scale,
    component        = component,
    type             = type,
    plot_type        = plot_type,
    lags             = lags,
    ...
  )
}

#' @export
#' @rdname plot_diagnostic
plot_diagnostic_density <- function(x, ...) UseMethod("plot_diagnostic_density")

#' @export
#' @rdname plot_diagnostic
plot_diagnostic_density.brma         <- function(
    x, parameter = NULL, parameter_mods = NULL, parameter_scale = NULL,
    type = "density", plot_type = "base", lags = 30,
    component = "auto", ...) {

  plot_diagnostic(
    x               = x,
    parameter       = parameter,
    parameter_mods   = parameter_mods,
    parameter_scale  = parameter_scale,
    component        = component,
    type             = type,
    plot_type        = plot_type,
    lags             = lags,
    ...
  )
}


### Common plotting helper functions ----
# Internal helper function for selecting and validating a plot parameter
#
# Processes parameter selection for plotting based on a parameter plus component
# namespace. The legacy parameter_mods/parameter_scale arguments are still
# supported.
#
# @param parameter Base parameter name ("mu", "tau", "rho")
# @param parameter_mods Moderator parameter for location regression
# @param parameter_scale Moderator parameter for scale regression
# @param object A brma object
#
# @return A character string with the validated parameter name (with appropriate prefix).
#   Attributes:
#   - parameter_samples: The JAGS samples name
