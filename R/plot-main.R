
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
#' @param parameter_mods legacy moderator selector. Prefer \code{parameter}
#' with \code{component = "mods"}. Use \code{"intercept"} for the
#' adjusted effect in meta-regression models.
#' @param parameter_scale legacy scale-regression selector. Prefer
#' \code{parameter} with \code{component = "scale"}. Use
#' \code{"intercept"} for the heterogeneity intercept in location-scale models.
#' @param component parameter component. Defaults to \code{"auto"}, which
#' infers the component when possible. Use \code{"mods"} (alias
#' \code{"location"}), \code{"scale"}, or \code{"bias"} to disambiguate
#' terms used in multiple model components.
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
#' @param density_method posterior density method. \code{"KDE"} uses the
#' standard BayesTools kernel density estimate. \code{"qCMDE"} attaches RoBMA
#' row-normalized q-grid conditional densities. \code{"IWMDE"} attaches
#' Chen-style moment-matched IWMDE densities. Matching is case-insensitive.
#' @param density_control named list of density-estimation settings. Supported
#' entries are \code{n_points}, \code{max_samples}, \code{display_grid},
#' \code{normalization_points}, and \code{normalization_prob}. The
#' normalization entries are only used with \code{density_method = "qCMDE"}.
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
    x, parameter = NULL, parameter_mods = NULL, parameter_scale = NULL,
    prior = FALSE, standardized_coefficients = FALSE,
    conditional = FALSE,
    output_measure = NULL, transform = NULL,
    plot_type = "base", dots_prior = NULL,
    density_method = c("KDE", "qCMDE", "IWMDE"),
    density_control = NULL, component = "auto", ...) {


  ### check user input
  BayesTools::check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  BayesTools::check_bool(prior, "prior")
  BayesTools::check_bool(standardized_coefficients, "standardized_coefficients")
  BayesTools::check_bool(conditional, "conditional")
  dots_raw <- list(...)
  .warn_unused_dots(
    dots    = dots_raw,
    allowed = .plot_dots_allowed(),
    caller  = "plot.brma()"
  )
  dots_raw <- .keep_allowed_dots(dots_raw, .plot_dots_allowed())
  density_method <- .density_method_normalize(density_method)
  if (.density_method_uses_precomputed(density_method) ||
      !is.null(density_control)) {
    density_control <- .density_control_normalize(
      density_method  = density_method,
      density_control = density_control
    )
  }
  if (conditional && !.is_RoBMA(x)) {
    stop("'conditional' plots are available only for RoBMA objects.", call. = FALSE)
  }

  ### select and validate the parameter to be plotted
  parameter <- .check_and_select_plot_parameter(
    parameter        = parameter,
    parameter_mods   = parameter_mods,
    parameter_scale  = parameter_scale,
    component        = component,
    object           = x
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
  density_sample_parameter <- .plot_brma_density_sample_parameter(
    samples          = samples,
    parameter        = parameter,
    sample_parameter = sample_parameter
  )
  if (.density_method_uses_precomputed(density_method)) {
    samples <- .plot_brma_attach_iwmde(
      object                  = x,
      samples                 = samples,
      parameter               = parameter,
      sample_parameter        = density_sample_parameter,
      conditional             = if (conditional) parameter else NULL,
      n_points                = density_control[["n_points"]],
      max_samples             = density_control[["max_samples"]],
      normalization_points    = density_control[["normalization_points"]],
      normalization_prob      = density_control[["normalization_prob"]],
      density_method          = density_method,
      display_grid            = density_control[["display_grid"]]
    )
  }

  ### set up plotting arguments
  n_levels   <- .get_samples_n_levels(samples, parameter)
  dots       <- do.call(.set_dots_plot, c(dots_raw, list(n_levels = n_levels)))
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
  args$density_method           <- if (
    .plot_brma_has_posterior_density(samples, density_sample_parameter)
  ) {
    "precomputed"
  } else {
    if (.density_method_uses_precomputed(density_method)) {
      warning(density_method, " density was not available; using KDE.",
              call. = FALSE)
    }
    "KDE"
  }
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


.plot_brma_attach_iwmde <- function(object, samples, parameter, sample_parameter,
                                    conditional,
                                    n_points, max_samples,
                                    normalization_points,
                                    normalization_prob, density_method,
                                    display_grid) {

  plotted_samples <- .plot_brma_plotted_samples(
    samples          = samples,
    sample_parameter = sample_parameter,
    parameter        = parameter
  )
  if (is.null(plotted_samples)) {
    return(samples)
  }
  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }

  context <- .iwmde_context(object)
  method  <- .density_method_iwmde_estimator(density_method)
  diagnostic <- .iwmde_parameter_diagnostic(
    context              = context,
    parameter            = parameter,
    n_points             = n_points,
    max_samples          = max_samples,
    normalization_points = normalization_points,
    normalization_prob   = normalization_prob,
    method               = method,
    display_grid_method  = display_grid,
    parameter_spec       = list(
      type             = "primitive",
      conditional      = conditional,
      conditional_rule = "AND"
    ),
    diagnostic_cache     = .iwmde_diagnostic_cache()
  )
  attr(samples, "iwmde_diagnostics") <- list(parameter = diagnostic)

  if (identical(diagnostic[["status"]], "ok")) {
    posterior_density <- .iwmde_posterior_density_attribute(
      diagnostic     = diagnostic,
      density_method = density_method
    )
    posterior_density <- .plot_brma_align_iwmde_density(
      posterior_density = posterior_density,
      raw_samples       = diagnostic[["samples"]],
      plotted_samples   = plotted_samples
    )
    if (!is.null(posterior_density)) {
      attr(samples[[sample_parameter]], "posterior_density") <- posterior_density
    }
  }

  return(samples)
}


.plot_brma_density_sample_parameter <- function(samples, parameter,
                                                sample_parameter) {

  if (parameter %in% c("PET", "PEESE", "omega") &&
      "bias" %in% sample_parameter &&
      !is.null(samples[["bias"]])) {
    return("bias")
  }
  if (!is.null(samples[[parameter]])) {
    return(parameter)
  }

  return(sample_parameter[[1L]])
}


.plot_brma_plotted_samples <- function(samples, sample_parameter, parameter) {

  sample <- samples[[sample_parameter]]
  if (is.null(sample) || is.list(sample)) {
    return(NULL)
  }
  if (is.matrix(sample)) {
    column <- grep(parameter, colnames(sample), fixed = TRUE)
    if (length(column) != 1L) {
      return(NULL)
    }
    return(as.numeric(sample[, column]))
  }
  if (!is.numeric(sample)) {
    return(NULL)
  }

  return(as.numeric(sample))
}


.plot_brma_has_posterior_density <- function(samples, sample_parameter) {

  if (is.null(samples[[sample_parameter]])) {
    return(FALSE)
  }

  return(!is.null(attr(samples[[sample_parameter]], "posterior_density")))
}


.plot_brma_align_iwmde_density <- function(posterior_density, raw_samples,
                                           plotted_samples) {

  transform <- .plot_brma_affine_sample_transform(raw_samples, plotted_samples)
  if (is.null(posterior_density) || is.null(transform)) {
    return(NULL)
  }
  if (isTRUE(all.equal(transform[["intercept"]], 0)) &&
      isTRUE(all.equal(transform[["slope"]], 1))) {
    return(posterior_density)
  }

  posterior_density[["x"]] <- transform[["intercept"]] +
    transform[["slope"]] * posterior_density[["x"]]
  posterior_density[["y"]] <- posterior_density[["y"]] / abs(transform[["slope"]])
  order_x <- order(posterior_density[["x"]])
  posterior_density[["x"]] <- posterior_density[["x"]][order_x]
  posterior_density[["y"]] <- posterior_density[["y"]][order_x]
  if (!is.null(posterior_density[["point_masses"]]) &&
      nrow(posterior_density[["point_masses"]]) > 0L) {
    posterior_density[["point_masses"]][["x"]] <- transform[["intercept"]] +
      transform[["slope"]] * posterior_density[["point_masses"]][["x"]]
  }
  posterior_density[["diagnostics"]][["plot_scale_transform"]] <- transform

  return(posterior_density)
}


.plot_brma_affine_sample_transform <- function(raw_samples, plotted_samples,
                                               tolerance = 1e-7) {

  raw_samples     <- as.numeric(raw_samples)
  plotted_samples <- as.numeric(plotted_samples)
  if (length(raw_samples) != length(plotted_samples) ||
      length(raw_samples) < 2L) {
    return(NULL)
  }

  keep <- is.finite(raw_samples) & is.finite(plotted_samples)
  raw_samples     <- raw_samples[keep]
  plotted_samples <- plotted_samples[keep]
  if (length(raw_samples) < 2L ||
      diff(range(raw_samples)) <= sqrt(.Machine$double.eps)) {
    return(NULL)
  }

  raw_center     <- raw_samples - mean(raw_samples)
  plotted_center <- plotted_samples - mean(plotted_samples)
  denominator    <- sum(raw_center^2)
  if (!is.finite(denominator) || denominator <= 0) {
    return(NULL)
  }

  slope     <- sum(raw_center * plotted_center) / denominator
  intercept <- mean(plotted_samples) - slope * mean(raw_samples)
  if (!is.finite(slope) || !is.finite(intercept) ||
      abs(slope) <= sqrt(.Machine$double.eps)) {
    return(NULL)
  }

  fitted <- intercept + slope * raw_samples
  scale  <- max(1, diff(range(plotted_samples)), max(abs(plotted_samples)))
  if (max(abs(fitted - plotted_samples)) > tolerance * scale) {
    return(NULL)
  }

  return(list(intercept = intercept, slope = slope))
}

