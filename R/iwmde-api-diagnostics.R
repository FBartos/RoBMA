# ============================================================================ #
# IWMDE Developer Diagnostic API
# ============================================================================ #

#' Plot IWMDE posterior-density diagnostics
#'
#' @description Compare sample KDE, histogram, and RoBMA-side IWMDE curves for
#' raw posterior parameters in a fitted model. This is a diagnostic helper for
#' checking likelihood-aware posterior-density estimates before routing them
#' through the regular plotting interface. The \code{"qCMDE"} curve uses a
#' grid-normalized conditional weight on an internal identity/log/logit scale.
#'
#' @details Weightfunction coordinates named by the active selection backend
#' and cumulative \code{eta} coordinates are reported as unsupported until a
#' joint replacement map is implemented. For spike-mixture parameters, KDE and
#' histogram layers use the continuous active samples scaled by their posterior
#' mass; point masses are drawn separately.
#'
#' @section Stability:
#' This is a developer-facing diagnostic API for validating the IWMDE feature.
#' The plotting contract and tuning arguments may be simplified after review.
#'
#' @param object a fitted \code{brma}, \code{BMA}, or \code{RoBMA} object.
#' @param parameters character vector of raw posterior parameter columns. If
#' \code{NULL}, all supported non-indicator scalar columns are attempted.
#' @param n_points number of grid points per density curve.
#' @param max_samples maximum posterior rows used by the density estimator per
#' parameter.
#' @param normalization_points number of initial internal-scale grid points used
#' for qCMDE row-wise conditional normalization and IWMDE support-grid
#' diagnostics. qCMDE internally refines and extends this grid for the final
#' row normalizers. Defaults to \code{max(50, n_points)}.
#' @param normalization_prob posterior central probability used to construct the
#' estimator-owned hidden normalization grid. Display ranges and requested point
#' ordinates do not force this grid; qCMDE refines and extends it internally for
#' validation.
#' @param density_method density estimator. \code{"qCMDE"} uses
#' row-normalized q-grid conditional densities. \code{"IWMDE"} uses
#' Chen-style moment-matched weights: conditional normal weights for
#' unconstrained targets, Gamma weights for one-sided targets, and
#' logit-conditional-normal weights for bounded targets. Matching is
#' case-insensitive.
#' @param display_grid support-point placement for the plotted IWMDE curve.
#' \code{"adaptive"} combines a uniform backbone, posterior quantiles, and
#' pilot-density height/curvature points. \code{"uniform"} uses equally spaced
#' points on the displayed parameter scale.
#' @param plot whether to draw the square diagnostic plot matrix.
#' @param as_data whether to return only the computed plot data.
#' @param ... additional graphical arguments passed to the base plotting path.
#'
#' @return A list with one entry per requested parameter. Returned visibly when
#' \code{as_data = TRUE} or \code{plot = FALSE}, and invisibly after plotting.
#'
#' @noRd
plot_iwmde_diagnostics <- function(object, parameters = NULL, n_points = 100,
                                   max_samples = 500,
                                   normalization_points = NULL,
                                   normalization_prob = .999,
                                   density_method = c("qCMDE", "IWMDE"),
                                   display_grid = c("adaptive", "uniform"),
                                   plot = TRUE, as_data = FALSE, ...) {

  dots_raw <- list(...)
  .warn_unused_dots(
    dots    = dots_raw,
    allowed = .iwmde_diagnostics_dots_allowed(),
    caller  = "plot_iwmde_diagnostics()"
  )
  dots_raw <- .keep_allowed_dots(dots_raw, .iwmde_diagnostics_dots_allowed())
  if (!inherits(object, "brma")) {
    stop("'object' must be a fitted brma, BMA, or RoBMA object.", call. = FALSE)
  }
  BayesTools::check_int(n_points, "n_points", lower = 20)
  BayesTools::check_int(max_samples, "max_samples", lower = 20)
  BayesTools::check_int(
    normalization_points, "normalization_points",
    lower       = 20,
    allow_NULL  = TRUE,
    check_length = 1
  )
  BayesTools::check_real(
    normalization_prob, "normalization_prob",
    lower        = 0,
    upper        = 1,
    check_length = 1
  )
  BayesTools::check_bool(plot, "plot")
  BayesTools::check_bool(as_data, "as_data")
  if (!is.null(parameters)) {
    BayesTools::check_char(parameters, "parameters", check_length = 0, allow_NA = FALSE)
  }
  density_method <- .density_method_normalize_precomputed(density_method)
  display_grid <- .iwmde_normalize_display_grid(display_grid)
  .check_iwmde_available(object, "plot_iwmde_diagnostics()")

  context <- .iwmde_context(object)
  if (is.null(parameters)) {
    parameters <- .iwmde_default_parameters(context)
  }
  if (length(parameters) == 0L) {
    stop("No posterior parameters are available for IWMDE diagnostics.", call. = FALSE)
  }
  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }
  estimate_cache <- .iwmde_estimate_cache()
  density_control <- list(
    n_points             = n_points,
    max_samples          = max_samples,
    normalization_points = normalization_points,
    normalization_prob   = normalization_prob,
    display_grid         = display_grid
  )

  out <- lapply(parameters, function(parameter) {
    estimate <- .iwmde_estimate(
      context         = context,
      parameter       = parameter,
      density_method  = density_method,
      density_control = density_control,
      outputs         = "density",
      cache           = estimate_cache
    )

    estimate[["diagnostics"]][["density"]]
  })
  names(out) <- parameters
  class(out) <- c("iwmde_diagnostics", "list")

  if (as_data || !plot) {
    return(out)
  }

  do.call(.plot_iwmde_diagnostics_base, c(list(x = out), dots_raw))
  return(invisible(out))
}


#' Plot IWMDE diagnostics for estimated marginal means
#'
#' @description Compare sample KDE, histogram, and requested density-estimator
#' curves for estimated marginal means. Each marginal mean is treated as the
#' linear function of the fitted primitive formula coefficients stored by
#' \code{marginal_means()}. Diagnostics are computed on the fitted linear
#' predictor scale.
#'
#' @section Stability:
#' This is a developer-facing diagnostic API for validating the IWMDE feature.
#' The plotting contract and tuning arguments may be simplified after review.
#'
#' @param object a fitted \code{brma}, \code{BMA}, or \code{RoBMA} object.
#' @param parameter optional moderator term(s) to plot. Defaults to all
#' available marginal-mean terms.
#' @param type marginal-mean type. For RoBMA product-space objects, use
#' \code{"averaged"} or \code{"conditional"}.
#' @param levels optional marginal-mean level names to keep.
#' @param marginal_means_object optional precomputed \code{marginal_means.brma}
#' object. If \code{NULL}, it is computed from \code{object}.
#' @param n_marginal_samples number of samples/grid points passed to
#' \code{marginal_means()} when \code{marginal_means_object = NULL}.
#' @inheritParams plot_iwmde_diagnostics
#'
#' @return A list with one entry per marginal-mean level. Returned visibly when
#' \code{as_data = TRUE} or \code{plot = FALSE}, and invisibly after plotting.
#'
#' @noRd
plot_iwmde_marginal_means_diagnostics <- function(object, parameter = NULL,
                                                  type = NULL, levels = NULL,
                                                  marginal_means_object = NULL,
                                                  n_marginal_samples = 10000,
                                                  n_points = 100,
                                                  max_samples = 500,
                                                  normalization_points = NULL,
                                                  normalization_prob = .999,
                                                  density_method = c("qCMDE", "IWMDE"),
                                                  display_grid = c("adaptive", "uniform"),
                                                  plot = TRUE, as_data = FALSE, ...) {

  dots_raw <- list(...)
  .warn_unused_dots(
    dots    = dots_raw,
    allowed = .iwmde_diagnostics_dots_allowed(),
    caller  = "plot_iwmde_marginal_means_diagnostics()"
  )
  dots_raw <- .keep_allowed_dots(dots_raw, .iwmde_diagnostics_dots_allowed())
  if (!inherits(object, "brma")) {
    stop("'object' must be a fitted brma, BMA, or RoBMA object.", call. = FALSE)
  }
  BayesTools::check_int(n_marginal_samples, "n_marginal_samples", lower = 2)
  BayesTools::check_int(n_points, "n_points", lower = 20)
  BayesTools::check_int(max_samples, "max_samples", lower = 20)
  BayesTools::check_int(
    normalization_points, "normalization_points",
    lower        = 20,
    allow_NULL   = TRUE,
    check_length = 1
  )
  BayesTools::check_real(
    normalization_prob, "normalization_prob",
    lower        = 0,
    upper        = 1,
    check_length = 1
  )
  BayesTools::check_bool(plot, "plot")
  BayesTools::check_bool(as_data, "as_data")
  if (!is.null(parameter)) {
    BayesTools::check_char(parameter, "parameter", check_length = 0, allow_NA = FALSE)
  }
  if (!is.null(levels)) {
    BayesTools::check_char(levels, "levels", check_length = 0, allow_NA = FALSE)
  }
  density_method <- .density_method_normalize_precomputed(density_method)
  display_grid <- .iwmde_normalize_display_grid(display_grid)

  if (is.null(marginal_means_object)) {
    marginal_means_object <- marginal_means(
      object    = object,
      n_samples = n_marginal_samples
    )
  } else {
    object <- .iwmde_marginal_means_source_object(marginal_means_object)
  }
  .check_iwmde_available(object, "plot_iwmde_marginal_means_diagnostics()")
  if (!inherits(marginal_means_object, "marginal_means.brma")) {
    stop("'marginal_means_object' must inherit from 'marginal_means.brma'.",
         call. = FALSE)
  }
  if (.effect_output_active(marginal_means_object[["effect_transform"]])) {
    stop("IWMDE marginal-mean diagnostics operate on the fitted linear ",
         "predictor scale. Recompute 'marginal_means_object' without ",
         "'output_measure' or 'transform'.", call. = FALSE)
  }
  type <- .marginal_means_type(marginal_means_object, type)
  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }
  estimate_cache <- .iwmde_estimate_cache()
  density_control <- list(
    n_points             = n_points,
    max_samples          = max_samples,
    normalization_points = normalization_points,
    normalization_prob   = normalization_prob,
    display_grid         = display_grid
  )

  context <- .iwmde_context(object)
  specs   <- .iwmde_marginal_means_specs(
    marginal_means_object = marginal_means_object,
    parameter             = parameter,
    type                  = type,
    levels                = levels
  )
  if (length(specs) == 0L) {
    stop("No marginal means are available for IWMDE diagnostics.", call. = FALSE)
  }

  out <- lapply(specs, function(spec) {
    estimate <- .iwmde_estimate(
      context         = context,
      parameter       = spec[["label"]],
      density_method  = density_method,
      density_control = density_control,
      outputs         = "density",
      parameter_spec  = spec,
      cache           = estimate_cache
    )

    estimate[["diagnostics"]][["density"]]
  })
  names(out) <- names(specs)
  class(out) <- c("iwmde_diagnostics", "list")

  if (as_data || !plot) {
    return(out)
  }

  do.call(.plot_iwmde_diagnostics_base, c(list(x = out), dots_raw))
  return(invisible(out))
}


.iwmde_marginal_means_source_object <- function(marginal_means_object) {

  if (!inherits(marginal_means_object, "marginal_means.brma")) {
    stop("'marginal_means_object' must inherit from 'marginal_means.brma'.",
         call. = FALSE)
  }

  source_object <- marginal_means_object[["source_object"]]
  if (is.null(source_object) ||
      !inherits(source_object, "brma") ||
      is.null(source_object[["fit"]])) {
    stop("'marginal_means_object' does not contain the source fitted brma object.",
         call. = FALSE)
  }

  return(source_object)
}
