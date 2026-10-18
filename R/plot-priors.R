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
#' requests multiple base parameters. A quantity of the fitted model that is no
#' model parameter (e.g. an inclusion indicator or a weight-function
#' coordinate such as \code{"omega[1]"}) stops with the classes
#' \code{RoBMA_hypothesis_target} and \code{RoBMA_hypothesis_unavailable}, as
#' in [plot.brma()].
#' @param parameter_mods legacy moderator selector. Prefer \code{parameter}
#' with \code{component = "mods"}. Use \code{"intercept"} for the
#' adjusted effect in meta-regression models.
#' @param parameter_scale legacy scale-regression selector. Prefer
#' \code{parameter} with \code{component = "scale"}. Use
#' \code{"intercept"} for the heterogeneity intercept in location-scale models.
#' @param component parameter component. Defaults to \code{"auto"}, which
#' infers the component when possible. Use \code{"mods"} (alias
#' \code{"location"}), \code{"scale"}, \code{"bias"}, or \code{"random"}
#' (the random-effect quantities of fitted objects) to disambiguate
#' terms used in multiple prior components. A \code{parameter_mods} or
#' \code{parameter_scale} selection with another \code{component} stops with
#' an error of class \code{RoBMA_component_mismatch}.
#' @param standardized_coefficients whether to plot moderator and scale-regression
#' priors on the standardized predictor scale. Defaults to \code{TRUE}, which
#' shows the priors as specified. Set to \code{FALSE} to transform them to the
#' original predictor scale when continuous predictors were standardized (see
#' Details for terms without an original-scale prior).
#' @inheritParams predict.brma
#' @param plot_type whether to use a base plot \code{"base"} or ggplot2
#' \code{"ggplot"} for plotting. Defaults to \code{"base"}.
#' @param ... additional arguments passed to the prior plotting method. For
#' mixed continuous and point priors, \code{ylim} controls the density axis,
#' \code{ylim2} controls the probability-mass axis, and \code{ylab2} controls
#' its label.
#'
#' @details \code{output_measure} and \code{transform} transform the prior
#' plotting scale only for effect-size location priors (\code{"mu"} or the
#' meta-regression intercept). For base plots, the probability mapping created
#' by \code{plot_prior()} is reused by subsequent posterior
#' \code{lines.brma()} calls, including overlays from other fitted objects.
#'
#' The random-effect quantities of fitted \code{brma.mv}, \code{BMA.mv}, and
#' \code{RoBMA.mv} objects are selected by their names in \code{parameter},
#' such as \code{"tau"}, \code{"tau_total"}, \code{"tau2"},
#' \code{"tau2_prop(study)"}, or \code{"study: tau"} (\code{"tau"} stops as
#' ambiguous when several blocks have an SD). A quantity that is one declared
#' prior (the SD of a block, the total SD of a variance allocation, or a
#' correlation with its own prior) shows that prior. Any other quantity
#' (a variance, a share of the total variance, the SD of a component of an
#' allocation, and the quantities of gated allocations) shows the prior law it
#' follows from the declared priors, the density that \code{plot()} overlays
#' with \code{prior = TRUE}, drawn without a posterior; point masses (for
#' example of inclusion gates) are drawn as arrows. For these laws,
#' \code{xlim}, \code{n_points}, \code{par_name}, and graphical arguments
#' apply as for other priors, while \code{x_seq} and \code{x_range_quant}
#' stop. A quantity without a prior density stops with the classes
#' \code{RoBMA_hypothesis_target} and \code{RoBMA_hypothesis_unavailable}.
#' \code{standardized_coefficients} selects the fitted (default) or the
#' original predictor scale of these laws. Objects created with
#' \code{only_priors = TRUE} have no fitted random-effect quantities; print
#' their random-effect priors with
#' \code{print_prior(x, parameter = "random")}.
#'
#' A term with several fitted coordinates (a factor with three or more levels
#' or with independent contrasts, and the interactions of such a factor) has
#' an original-scale prior only when the standardization of continuous
#' predictors leaves its coordinates unchanged. With
#' \code{standardized_coefficients = FALSE}, the prior of such a term whose
#' coordinates the standardization changes (a factor or factor-by-factor
#' interaction next to its interaction with a standardized predictor, and that
#' interaction term) stops with an error of class
#' \code{RoBMA_density_method_original_scale} and its parent
#' \code{RoBMA_density_method_unavailable}; use
#' \code{standardized_coefficients = TRUE}. The other terms (including factors
#' and factor-by-factor interactions without an interaction with a
#' standardized predictor) are plotted as before.
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
#'   plot_prior(priors, parameter = "Preregistered", component = "mods")
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
    x, parameter = "mu", parameter_mods = NULL, parameter_scale = NULL,
    standardized_coefficients = TRUE,
    output_measure = NULL, transform = NULL,
    plot_type = "base", component = "auto", ...) {

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
        component                 = component,
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
    object           = x,
    parameter        = if (has_parameter) parameter else NULL,
    parameter_mods   = if (has_parameter_mods) parameter_mods else NULL,
    parameter_scale  = if (has_parameter_scale) parameter_scale else NULL,
    component        = component,
    standardized_coefficients = standardized_coefficients
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

  if (inherits(selected[["prior"]], "prior_linear_density")) {
    plot <- .plot_prior_density_law(
      law       = selected[["prior"]],
      plot_type = plot_type,
      dots      = dots
    )
  } else if (BayesTools::is.prior.spike_and_slab(selected[["prior"]])) {
    plot <- do.call(
      "plot",
      c(list(x = selected[["prior"]], plot_type = plot_type), dots)
    )
  } else if (BayesTools::is.prior.mixture(selected[["prior"]])) {
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
#' scale terms can also be selected by name when unambiguous. For
#' \code{brma.mv()} objects with random-effect formulas, \code{"random"} prints
#' the structured random-effect prior specification, and the random-effect
#' quantities of fitted objects (such as \code{"tau"}, \code{"tau_total"},
#' \code{"tau2"}, \code{"tau2_prop(study)"}, or \code{"study: tau"}) print
#' their declared prior, or, for a quantity that follows from several declared
#' priors, one line giving the derivation (for example
#' \code{tau2 = tau^2} or \code{study: tau = sigma_total * sqrt(w[1])}, in the
#' notation of the \code{"random"} specification) followed by the declared
#' priors it refers to. A quantity without a prior density stops as in
#' [plot_prior()]. A character vector requests multiple base parameters. A
#' quantity of the fitted model that is no model parameter (e.g. an inclusion
#' indicator or a weight-function coordinate such as \code{"omega[1]"}) stops
#' with the classes \code{RoBMA_hypothesis_target} and
#' \code{RoBMA_hypothesis_unavailable}, as in [plot.brma()].
#' @param parameter_mods legacy moderator selector. Prefer \code{parameter}
#' with \code{component = "mods"}. Use \code{"intercept"} for the
#' adjusted effect in meta-regression models.
#' @param parameter_scale legacy scale-regression selector. Prefer
#' \code{parameter} with \code{component = "scale"}. Use
#' \code{"intercept"} for the heterogeneity intercept in location-scale models.
#' @param component parameter component. Defaults to \code{"auto"}, which
#' infers the component when possible. Use \code{"mods"} (alias
#' \code{"location"}), \code{"scale"}, \code{"bias"}, or \code{"random"}
#' (the random-effect quantities of fitted objects) to disambiguate
#' terms used in multiple prior components. A \code{parameter_mods} or
#' \code{parameter_scale} selection with another \code{component} stops with
#' an error of class \code{RoBMA_component_mismatch}.
#' @param ... additional arguments passed to the prior printing method. Use
#' \code{silent = TRUE} for programmatic inspection without console output.
#'
#' @return \code{print_prior} invisibly returns the selected prior distribution
#' (for a random-effect quantity that is not one declared prior, its prior
#' density, a BayesTools \code{prior_linear_density}).
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
#'   print_prior(priors, parameter = "Preregistered", component = "mods")
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
    x, parameter = NULL, parameter_mods = NULL, parameter_scale = NULL,
    component = "auto", ...) {

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

  if (parameter_missing && !has_parameter_mods && !has_parameter_scale &&
      identical(.parameter_component_normalize(component), "auto")) {
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
      .select_print_prior_parameter(
        object           = x,
        parameter        = parameter_i,
        parameter_mods   = NULL,
        parameter_scale  = NULL,
        component        = component,
        allow_mixed_bias = TRUE
      )
    })
    names(selected) <- parameter

    priors <- lapply(selected, `[[`, "prior")
    for (i in seq_along(selected)) {
      .print_prior_selected(selected[[i]], ...)
    }

    return(invisible(priors))
  }

  selected <- .select_print_prior_parameter(
    object           = x,
    parameter        = if (has_parameter) parameter else NULL,
    parameter_mods   = if (has_parameter_mods) parameter_mods else NULL,
    parameter_scale  = if (has_parameter_scale) parameter_scale else NULL,
    component        = component,
    allow_mixed_bias = TRUE
  )

  .print_prior_selected(selected, ...)

  return(invisible(selected[["prior"]]))
}

#' @export
print.only_priors.brma <- function(x, ...) {

  print_prior.brma(x = x, ...)
}
