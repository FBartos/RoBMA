
### basic plotting functions ----
#' @title Plots brma Object
#'
#' @description \code{plot.brma} visualizes posterior
#' (and prior) distribution a brma object.
#'
#' @param x a fitted \code{brma}, \code{BMA}, or \code{RoBMA} object.
#' @param parameter a parameter to be plotted. Defaults to \code{"mu"} for
#' the effect size, or to the meta-regression intercept when moderators are
#' present. Additional options are \code{"tau"}, \code{"rho"} for multilevel
#' models, \code{"PET"}, \code{"PEESE"}, and \code{"omega"} or
#' \code{"weightfunction"} for selection models. Use \code{plot_pet_peese()}
#' for PET/PEESE regression plots.
#' @param parameter_mods character. Moderator term to plot. Use
#' \code{"intercept"} for the adjusted effect in meta-regression models.
#' @param parameter_scale character. Scale-regression term to plot. Use
#' \code{"intercept"} for the heterogeneity intercept in location-scale models.
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting. Defaults to
#' \code{"base"}.
#' @param prior whether prior distribution should be added to
#' figure. Defaults to \code{FALSE}.
#' @param standardized_coefficients whether to plot moderator and
#' scale-regression coefficients on the standardized predictor scale. Defaults
#' to \code{FALSE}.
#' @param conditional whether to plot the conditional posterior distribution
#' for RoBMA product-space objects. Defaults to \code{FALSE}.
#' @param dots_prior list of additional graphical arguments
#' to be passed to the plotting function of the prior
#' distribution. Supported arguments are \code{lwd},
#' \code{lty}, \code{col}, and \code{col.fill}, to adjust
#' the line thickness, line type, line color, and fill color
#' of the prior distribution respectively.
#' @inheritParams predict.brma
#' @param ... list of additional graphical arguments
#' to be passed to the plotting function. Supported arguments
#' are \code{lwd}, \code{lty}, \code{col}, \code{col.fill},
#' \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}
#' to adjust the line thickness, line type, line color, fill color,
#' x-label, y-label, title, x-axis range, and y-axis range
#' respectively.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   plot(fit, parameter = "mu")
#'   plot(fit, parameter = "tau", prior = TRUE)
#'   plot(fit, parameter = "PET")
#' }
#' }
#'
#'
#' @return \code{plot.brma} returns either \code{NULL} if \code{plot_type = "base"}
#' or a \code{ggplot2} object if \code{plot_type = "ggplot"}.
#'
#' @seealso [RoBMA()]
#' @export
plot.brma  <- function(
    x, parameter, parameter_mods, parameter_scale,
    prior = FALSE, standardized_coefficients = FALSE,
    conditional = FALSE,
    output_measure = NULL, transform = NULL,
    plot_type = "base", dots_prior = NULL, ...) {


  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(standardized_coefficients, "standardized_coefficients")
  BayesTools::check_bool(conditional, "conditional")
  if (conditional && !.is_RoBMA(x)) {
    stop("'conditional' plots are available only for RoBMA objects.", call. = FALSE)
  }

  ### select and validate the parameter to be plotted
  parameter <- .check_and_select_plot_parameter(
    parameter       = parameter,
    parameter_mods  = parameter_mods,
    parameter_scale = parameter_scale,
    object          = x
  )
  effect_transform <- .effect_output_setup(
    object         = x,
    output_measure = output_measure,
    transform      = transform
  )
  if (.effect_output_requested(effect_transform) &&
      !.is_effect_location_parameter(parameter)) {
    stop(
      "'output_measure' and 'transform' are only available for effect-size ",
      "location parameters ('mu' or the meta-regression intercept).",
      call. = FALSE
    )
  }

  ### obtain posterior samples in the plotting format
  sample_parameter <- .as_mixed_posteriors_parameters(x, parameter)
  samples <- BayesTools::as_mixed_posteriors(
    model            = x[["fit"]],
    parameters       = sample_parameter,
    conditional      = if (conditional) parameter else NULL,
    transform_scaled = !standardized_coefficients
  )

  ### set up plotting arguments
  n_levels   <- .get_samples_n_levels(samples, parameter)
  dots       <- .set_dots_plot(..., n_levels = n_levels)
  dots_prior <- .set_dots_prior(dots_prior)
  if (is.null(dots[["par_name"]])) {
    dots[["par_name"]] <- .plot_parameter_label(parameter, effect_transform)
  }

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
  if (.effect_output_requested(effect_transform)) {
    args$transformation           <- .effect_plot_transformation(effect_transform)
    args$transformation_arguments <- NULL
    args$transformation_settings  <- FALSE
  }

  # suppress messages about transformations
  plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}

