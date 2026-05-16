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
#' \code{"PEESE"}, \code{"omega"}, and \code{"bias"}, with aliases
#' \code{"effect"} = \code{"mu"}, \code{"heterogeneity"} = \code{"tau"},
#' and \code{"weightfunction"} = \code{"omega"}. \code{"bias"} plots only
#' non-mixed or homogeneous bias priors; for mixed weightfunction and PET/PEESE
#' mixtures use \code{"omega"}, \code{"PET"}, or \code{"PEESE"}. Moderator and
#' scale terms can also be selected by name when unambiguous. A character vector
#' requests multiple base parameters.
#' @param parameter_mods character. Moderator term to plot.
#' Use \code{"intercept"} for the adjusted effect in meta-regression models.
#' @param parameter_scale character. Scale-regression term to plot.
#' Use \code{"intercept"} for the heterogeneity intercept in location-scale models.
#' @param standardized_coefficients whether to plot moderator and scale-regression
#' priors on the standardized predictor scale. Defaults to \code{TRUE}, which
#' shows the priors as specified. Set to \code{FALSE} to transform them to the
#' original predictor scale when continuous predictors were standardized.
#' @inheritParams predict.brma
#' @param plot_type whether to use a base plot \code{"base"} or ggplot2
#' \code{"ggplot"} for plotting. Defaults to \code{"base"}.
#' @param ... additional arguments passed to the prior plotting method.
#'
#' @details \code{output_measure} and \code{transform} transform the prior
#' plotting scale only for effect-size location priors (\code{"mu"} or the
#' meta-regression intercept).
#'
#' @return \code{plot_prior} returns either \code{NULL} invisibly if
#' \code{plot_type = "base"} or a \code{ggplot2} object if
#' \code{plot_type = "ggplot"}. If multiple parameters are requested, a named
#' list is returned, invisibly for base plots.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   priors <- BMA(
#'     yi           = yi,
#'     vi           = vi,
#'     mods         = ~ Preregistered,
#'     data         = dat.lehmann2018,
#'     measure      = "SMD",
#'     only_priors  = TRUE
#'   )
#'
#'   plot_prior(priors, parameter = "mu")
#'   plot_prior(priors, parameter = "tau")
#'   plot_prior(priors, parameter_mods = "Preregistered")
#' }
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
    standardized_coefficients = TRUE,
    output_measure = NULL, transform = NULL,
    plot_type = "base", ...) {

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
        output_measure            = output_measure,
        transform                 = transform,
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
  effect_transform <- .effect_output_setup(
    object         = x,
    output_measure = output_measure,
    transform      = transform
  )
  if (.effect_output_requested(effect_transform) &&
      !.is_effect_location_parameter(selected[["label"]])) {
    stop(
      "'output_measure' and 'transform' are only available for effect-size ",
      "location priors ('mu' or the meta-regression intercept).",
      call. = FALSE
    )
  }

  dots <- list(...)
  if (is.null(dots[["par_name"]])) {
    dots[["par_name"]] <- selected[["label"]]
  }
  if (.effect_output_requested(effect_transform)) {
    dots[["transformation"]]           <- .effect_plot_transformation(effect_transform)
    dots[["transformation_arguments"]] <- NULL
    dots[["transformation_settings"]]  <- FALSE
    if (.effect_output_active(effect_transform)) {
      dots[["par_name"]] <- effect_transform[["label"]]
    }
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

  if (BayesTools::is.prior.mixture(selected[["prior"]])) {
    plot <- do.call(
      BayesTools::plot_prior_list,
      c(
        list(
          prior_list = .prior_mixture_plot_list(selected[["prior"]]),
          plot_type  = plot_type
        ),
        dots
      )
    )
  } else if (.use_plot_prior_list_dispatch(selected[["prior"]])) {
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

#' @export
plot.only_priors.brma <- function(x, ...) {

  plot_prior.brma(x = x, ...)
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
#' @param parameter character. Base parameter to print. If omitted, all stored
#' prior distributions are printed.
#' Common options are \code{"mu"}, \code{"tau"}, \code{"rho"}, \code{"PET"},
#' \code{"PEESE"}, \code{"omega"}, and \code{"bias"}, with aliases
#' \code{"effect"} = \code{"mu"}, \code{"heterogeneity"} = \code{"tau"},
#' and \code{"weightfunction"} = \code{"omega"}. GLMM outcome priors
#' \code{"pi"} and \code{"phi"} are available when present. Moderator and
#' scale terms can also be selected by name when unambiguous. A character vector
#' requests multiple base parameters.
#' @param parameter_mods character. Moderator term to print.
#' Use \code{"intercept"} for the adjusted effect in meta-regression models.
#' @param parameter_scale character. Scale-regression term to print.
#' Use \code{"intercept"} for the heterogeneity intercept in location-scale models.
#' @param ... additional arguments passed to the prior printing method. Use
#' \code{silent = TRUE} for programmatic inspection without console output.
#'
#' @return \code{print_prior} invisibly returns the selected prior distribution.
#' If multiple parameters are requested, a named list of prior distributions is
#' returned invisibly.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   priors <- BMA(
#'     yi          = yi,
#'     vi          = vi,
#'     mods        = ~ Preregistered,
#'     data        = dat.lehmann2018,
#'     measure     = "SMD",
#'     only_priors = TRUE
#'   )
#'
#'   print_prior(priors)
#'   print_prior(priors, parameter = "mu")
#'   print_prior(priors, parameter = "tau")
#'   print_prior(priors, parameter_mods = "Preregistered")
#' }
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
    x, parameter, parameter_mods, parameter_scale, ...) {

  if (is.null(x[["priors"]])) {
    stop("The object does not contain prior distributions.", call. = FALSE)
  }

  parameter_missing   <- missing(parameter)
  has_parameter       <- !parameter_missing        && !is.null(parameter)
  has_parameter_mods  <- !missing(parameter_mods)  && !is.null(parameter_mods)
  has_parameter_scale <- !missing(parameter_scale) && !is.null(parameter_scale)

  n_specified <- sum(c(has_parameter, has_parameter_mods, has_parameter_scale))
  if (n_specified > 1) {
    stop("Only one of 'parameter', 'parameter_mods', or 'parameter_scale' can be specified.", call. = FALSE)
  }

  if (parameter_missing && !has_parameter_mods && !has_parameter_scale) {
    selected <- .select_print_prior_all(x)
    priors   <- lapply(selected, `[[`, "prior")
    names(priors) <- vapply(selected, `[[`, character(1), "label")

    for (i in seq_along(selected)) {
      .print_prior_object(priors[[i]], label = names(priors)[i], ...)
    }

    return(invisible(priors))
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

#' @export
print.only_priors.brma <- function(x, ...) {

  print_prior.brma(x = x, ...)
}

