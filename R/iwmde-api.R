# ============================================================================ #
# IWMDE Public API and Diagnostic Orchestration
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


.iwmde_context <- function(object) {

  posterior_samples <- as.matrix(.get_posterior_samples(object[["fit"]]))
  if (is.null(colnames(posterior_samples))) {
    stop("Posterior samples must have column names.", call. = FALSE)
  }

  data   <- object[["data"]]
  priors <- object[["priors"]]

  context <- list(
    object            = object,
    data              = data,
    priors            = priors,
    posterior_samples = posterior_samples,
    flat_prior_list   = attr(object[["fit"]], "prior_list"),
    selection_spec     = .iwmde_selection_spec(data, priors),
    formula_fit        = object[["fit"]],
    formula_inputs     = .iwmde_formula_inputs(data, priors),
    indicator_names    = .iwmde_indicator_names(posterior_samples),
    active_cache       = new.env(parent = emptyenv()),
    focal_prior_cache  = new.env(parent = emptyenv()),
    support_cache      = new.env(parent = emptyenv()),
    likelihood_cache   = new.env(parent = emptyenv()),
    row_cache          = new.env(parent = emptyenv()),
    predictor_cache    = new.env(parent = emptyenv())
  )

  class(context) <- "iwmde_context"
  return(context)
}


.check_iwmde_available <- function(object, caller) {

  if (.is_random(object) && !.is_data_known_v(object[["data"]])) {
    .check_random_formula_postfit_deferred(object, caller)
  }

  return(invisible(TRUE))
}


.iwmde_context_unavailable_reason <- function(context) {

  object <- context[["object"]]
  data   <- context[["data"]]
  if (!is.null(object)) {
    data <- object[["data"]]
  }

  if (is.null(data) || !.is_data_random(data)) {
    return(NULL)
  }
  if (.is_data_known_v(data)) {
    return(NULL)
  }

  "qCMDE/IWMDE is not implemented for brma.mv() random-formula models yet."
}


.iwmde_context_ensure_caches <- function(context) {

  cache_names <- c(
    "active_cache",
    "focal_prior_cache",
    "support_cache",
    "likelihood_cache",
    "row_cache",
    "predictor_cache"
  )
  for (cache_name in cache_names) {
    if (!is.environment(context[[cache_name]])) {
      context[[cache_name]] <- new.env(parent = emptyenv())
    }
  }
  if (is.null(context[["indicator_names"]])) {
    context[["indicator_names"]] <- character()
  }
  if (is.null(context[["priors"]])) {
    context[["priors"]] <- list()
  }
  if (is.null(context[["flat_prior_list"]])) {
    context[["flat_prior_list"]] <- list()
  }
  if (is.null(context[["formula_inputs"]])) {
    context[["formula_inputs"]] <- list()
  }

  return(context)
}


.iwmde_indicator_names <- function(posterior_samples) {

  indicator_names <- grep("(^|_)indicator$", colnames(posterior_samples), value = TRUE)
  indicator_names <- c(
    indicator_names,
    intersect("bias_indicator", colnames(posterior_samples))
  )
  indicator_names <- unique(indicator_names)

  return(indicator_names)
}


.iwmde_default_parameters <- function(context) {

  parameters <- colnames(context[["posterior_samples"]])
  parameters <- parameters[!grepl("(^|_)indicator$", parameters)]
  parameters <- parameters[parameters != "bias_indicator"]

  keep <- !vapply(parameters, function(parameter) {
    .iwmde_parameter_is_weightfunction_coordinate(
      parameter = parameter,
      context   = context
    )
  }, logical(1))

  return(parameters[keep])
}


.iwmde_normalize_method <- function(method) {

  return(match.arg(method, c("q_grid_cmde", "iwmde")))
}


.density_method_normalize <- function(density_method, allow_normal = FALSE) {

  allowed <- c("KDE", "qCMDE", "IWMDE")
  if (isTRUE(allow_normal)) {
    allowed <- c("KDE", "normal", "qCMDE", "IWMDE")
  }

  return(.density_method_match(density_method, allowed))
}


.density_method_normalize_precomputed <- function(density_method) {

  density_method <- .density_method_match(density_method, c("qCMDE", "IWMDE"))

  return(density_method)
}


.density_method_match <- function(density_method, allowed) {

  if (is.null(density_method)) {
    density_method <- allowed[[1L]]
  }

  if (is.character(density_method) &&
      length(density_method) > 0L &&
      !anyNA(density_method)) {
    density_method_lower <- tolower(density_method)
    allowed_lower        <- tolower(allowed)

    if (length(density_method) > 1L) {
      if (identical(density_method_lower, allowed_lower)) {
        density_method <- allowed
      }
    } else {
      index <- pmatch(density_method_lower, allowed_lower, duplicates.ok = TRUE)
      if (!is.na(index)) {
        density_method <- allowed[[index]]
      }
    }
  }

  return(BayesTools::posterior_density_method_match(
    method  = density_method,
    allowed = allowed,
    name    = "density_method"
  ))
}


.iwmde_posterior_metadata <- function(samples, parameter, level = NULL) {

  metadata <- c(
    list(
      parameter = parameter,
      level     = level
    ),
    .iwmde_sample_condition_metadata(samples)
  )
  metadata <- metadata[!vapply(metadata, is.null, logical(1))]

  return(metadata)
}


.iwmde_sample_condition_metadata <- function(samples,
                                            include_prior_density_context = FALSE,
                                            require_child_condition = FALSE) {

  conditional <- .iwmde_first_nonempty_attr(
    samples,
    c("effective_conditional", "conditional")
  )
  conditional_rule <- .iwmde_first_nonempty_attr(
    samples,
    c("effective_conditional_rule", "conditional_rule")
  )
  condition_key <- attr(samples, "condition_key", exact = TRUE)

  if (isTRUE(require_child_condition) &&
      (is.null(conditional) || is.null(conditional_rule) || is.null(condition_key))) {
    stop(
      "Conditional marginal-means qCMDE/IWMDE metadata is incomplete; recompute the marginal_means object with current BayesTools.",
      call. = FALSE
    )
  }

  metadata <- list(
    conditional              = conditional,
    conditional_rule         = conditional_rule,
    condition_key            = condition_key,
    condition_event          = attr(samples, "condition_event", exact = TRUE),
    resolved_condition_event = attr(samples, "resolved_condition_event", exact = TRUE)
  )
  if (isTRUE(include_prior_density_context)) {
    metadata[["prior_density_context"]] <- attr(
      samples,
      "prior_density_context",
      exact = TRUE
    )
  }
  metadata <- metadata[!vapply(metadata, is.null, logical(1))]

  return(metadata)
}


.iwmde_first_nonempty_attr <- function(x, names) {

  for (name in names) {
    value <- attr(x, name, exact = TRUE)
    if (!is.null(value) && length(value) > 0L) {
      return(value)
    }
  }

  return(NULL)
}


.density_method_uses_precomputed <- function(density_method, allow_normal = FALSE) {

  density_method <- .density_method_normalize(
    density_method = density_method,
    allow_normal   = allow_normal
  )

  return(BayesTools::posterior_density_method_uses_precomputed(density_method))
}


.density_method_iwmde_estimator <- function(density_method) {

  density_method <- .density_method_normalize(density_method)
  if (identical(density_method, "KDE")) {
    return(NULL)
  }
  if (identical(density_method, "qCMDE")) {
    return("q_grid_cmde")
  }

  return("iwmde")
}


.density_control_normalize <- function(density_method, density_control = NULL,
                                       allow_normal = FALSE) {

  density_method <- .density_method_normalize(
    density_method = density_method,
    allow_normal   = allow_normal
  )
  allowed_names  <- c(
    "n_points",
    "max_samples",
    "normalization_points",
    "normalization_prob",
    "display_grid"
  )
  defaults <- list(
    n_points             = 100L,
    max_samples          = 500L,
    normalization_points = NULL,
    normalization_prob   = .999,
    display_grid         = "adaptive"
  )

  if (is.null(density_control)) {
    return(defaults)
  }
  if (!is.list(density_control)) {
    stop("'density_control' must be a named list.", call. = FALSE)
  }
  if (length(density_control) == 0L) {
    return(defaults)
  }

  control_names <- names(density_control)
  if (is.null(control_names) || any(!nzchar(control_names))) {
    stop("'density_control' must be a fully named list.", call. = FALSE)
  }
  duplicated_names <- unique(control_names[duplicated(control_names)])
  if (length(duplicated_names) > 0L) {
    stop("'density_control' contains duplicate setting(s): ",
         paste0("'", duplicated_names, "'", collapse = ", "), ".",
         call. = FALSE)
  }
  unknown_names <- setdiff(control_names, allowed_names)
  if (length(unknown_names) > 0L) {
    stop("'density_control' contains unrecognized setting(s): ",
         paste0("'", unknown_names, "'", collapse = ", "), ".",
         call. = FALSE)
  }
  if (density_method %in% c("KDE", "normal")) {
    stop("'density_control' is only used when 'density_method' is ",
         "'qCMDE' or 'IWMDE'.", call. = FALSE)
  }

  for (name in control_names) {
    defaults[[name]] <- density_control[[name]]
  }

  BayesTools::check_int(defaults[["n_points"]], "density_control$n_points", lower = 20)
  BayesTools::check_int(defaults[["max_samples"]], "density_control$max_samples", lower = 20)
  BayesTools::check_int(
    defaults[["normalization_points"]], "density_control$normalization_points",
    lower        = 20,
    allow_NULL   = TRUE,
    check_length = 1
  )
  BayesTools::check_real(
    defaults[["normalization_prob"]], "density_control$normalization_prob",
    lower        = 0,
    upper        = 1,
    check_length = 1,
    allow_NA     = FALSE
  )
  if (defaults[["normalization_prob"]] <= 0) {
    stop("'density_control$normalization_prob' must be higher than 0.",
         call. = FALSE)
  }
  defaults[["display_grid"]] <- .iwmde_normalize_display_grid(defaults[["display_grid"]])

  return(defaults)
}


.iwmde_posterior_density_attribute <- function(diagnostic, density_method,
                                               metadata = NULL,
                                               density_control = NULL) {

  if (!identical(diagnostic[["status"]], "ok") ||
      is.null(diagnostic[["iwmde"]])) {
    return(NULL)
  }

  density <- diagnostic[["iwmde"]]
  diagnostics <- diagnostic[["diagnostics"]]
  density_failure <- .iwmde_diagnostics_density_failure_reason(diagnostics)
  if (!is.null(density_failure)) {
    return(NULL)
  }
  density_warning <- .iwmde_diagnostics_density_warning(diagnostics)
  if (length(density_warning) > 0L) {
    diagnostics[["warning"]] <- density_warning
  }
  provenance <- .iwmde_result_provenance(
    diagnostic      = diagnostic,
    density_method  = density_method,
    metadata        = metadata,
    density_control = density_control,
    attribute       = "density"
  )
  args <- c(
    list(
      x              = density[["x"]],
      y              = density[["y"]],
      method         = density[["estimator"]],
      density_method = .density_method_normalize(density_method),
      diagnostics    = diagnostics,
      point_masses   = diagnostic[["point_masses"]],
      iwmde_provenance = provenance
    ),
    metadata
  )
  out <- do.call(
    BayesTools::posterior_density_attribute,
    args
  )

  return(out)
}


.iwmde_posterior_ordinate_attribute <- function(diagnostic, density_method,
                                                metadata = NULL,
                                                density_control = NULL) {

  if (!identical(diagnostic[["status"]], "ok") ||
      is.null(diagnostic[["diagnostics"]])) {
    return(NULL)
  }

  diagnostics <- diagnostic[["diagnostics"]]
  if (!isTRUE(diagnostics[["bf_included"]]) ||
      !is.finite(diagnostics[["bf_ordinate"]]) ||
      diagnostics[["bf_ordinate"]] <= 0) {
    return(NULL)
  }

  ordinate_diagnostics <- list(
    evaluation_value                  = diagnostics[["bf_evaluation_value"]],
    mcse                              = diagnostics[["bf_mcse"]],
    relative_mcse                     = diagnostics[["bf_relative_mcse"]],
    BF_error_percent                  = diagnostics[["bf_error_percent"]],
    finite_terms                      = diagnostics[["bf_finite_terms"]],
    ess                               = diagnostics[["bf_ess"]],
    max_weight_share                  = diagnostics[["bf_max_weight_share"]],
    max_log_ratio                     = diagnostics[["bf_max_log_ratio"]],
    active_mass                       = diagnostics[["active_mass"]],
    n_candidate_rows                  = diagnostics[["n_candidate_rows"]],
    n_denominator_rows                = diagnostics[["n_denominator_rows"]],
    n_estimator_rows                  = diagnostics[["n_estimator_rows"]],
    n_evaluated_rows                  = diagnostics[["n_evaluated_rows"]],
    n_normalized_rows                 = diagnostics[["n_normalized_rows"]],
    n_dropped_rows                    = diagnostics[["n_dropped_rows"]],
    row_drop_fraction                 = diagnostics[["row_drop_fraction"]],
    n_dropped_log_q                   = diagnostics[["n_dropped_log_q"]],
    n_dropped_weight                  = diagnostics[["n_dropped_weight"]],
    weight_partitions                 = diagnostics[["weight_partitions"]],
    n_dropped_normalizer              = diagnostics[["n_dropped_normalizer"]],
    n_validation_dropped_normalizer   = diagnostics[["n_validation_dropped_normalizer"]],
    pilot_ordinate                    = diagnostics[["bf_pilot_ordinate"]],
    validation_ordinate               = diagnostics[["bf_validation_ordinate"]],
    ordinate_relative_change          = diagnostics[["bf_ordinate_relative_change"]],
    ordinate_log_change               = diagnostics[["bf_ordinate_log_change"]],
    pilot_ordinate_relative_change    =
      diagnostics[["bf_pilot_ordinate_relative_change"]],
    pilot_ordinate_log_change         =
      diagnostics[["bf_pilot_ordinate_log_change"]],
    normalization_integral            = diagnostics[["normalization_integral"]],
    normalization_final_integral      = diagnostics[["normalization_final_integral"]],
    normalization_mass_ratio          = diagnostics[["normalization_mass_ratio"]],
    normalization_range               = diagnostics[["normalization_range"]],
    normalization_initial_points      = diagnostics[["normalization_initial_points"]],
    normalization_initial_range       = diagnostics[["normalization_initial_range"]],
    max_normalizer_relative_change    = diagnostics[["max_normalizer_relative_change"]],
    p95_normalizer_relative_change    = diagnostics[["p95_normalizer_relative_change"]],
    median_normalizer_relative_change = diagnostics[["median_normalizer_relative_change"]],
    normalization_refined_points      = diagnostics[["normalization_refined_points"]],
    normalization_refined_range       = diagnostics[["normalization_refined_range"]],
    n_refinement_steps                = diagnostics[["n_refinement_steps"]],
    estimator                         = diagnostics[["estimator"]],
    weight_method                     = diagnostics[["weight_method"]]
  )
  ordinate_warning <- .iwmde_diagnostics_bf_warning(ordinate_diagnostics)
  if (length(ordinate_warning) > 0L) {
    ordinate_diagnostics[["warning"]] <- ordinate_warning
  }
  provenance <- .iwmde_result_provenance(
    diagnostic       = diagnostic,
    density_method   = density_method,
    metadata         = metadata,
    density_control  = density_control,
    value            = diagnostics[["bf_value"]],
    evaluation_value = diagnostics[["bf_evaluation_value"]],
    attribute        = "ordinate"
  )

  args <- c(
    list(
      value            = diagnostics[["bf_value"]],
      ordinate         = diagnostics[["bf_ordinate"]],
      method           = diagnostics[["estimator"]],
      density_method   = .density_method_normalize(density_method),
      diagnostics      = ordinate_diagnostics,
      evaluation_value = diagnostics[["bf_evaluation_value"]],
      iwmde_provenance = provenance
    ),
    metadata
  )
  out <- do.call(
    BayesTools::posterior_ordinate_attribute,
    args
  )

  if (!.iwmde_posterior_ordinate_supports_bf(out)) {
    return(NULL)
  }

  return(out)
}


.iwmde_posterior_ordinate_supports_bf <- function(posterior_ordinate) {

  BayesTools::posterior_ordinate_supports_bf(
    ordinate  = posterior_ordinate,
    validator = .iwmde_posterior_ordinate_bf_validator
  )
}


.iwmde_posterior_ordinate_bf_validator <- function(posterior_ordinate) {

  return(is.null(
    .iwmde_posterior_ordinate_bf_failure_reason(posterior_ordinate)
  ))
}


.iwmde_posterior_ordinate_bf_failure_reason <- function(posterior_ordinate) {

  diagnostics <- posterior_ordinate[["diagnostics"]]
  if (is.null(diagnostics) || !is.list(diagnostics)) {
    return("missing qCMDE/IWMDE BF diagnostics")
  }

  ordinate <- .iwmde_ordinate_scalar(posterior_ordinate, "ordinate")
  if (!is.finite(ordinate)) {
    ordinate <- .iwmde_ordinate_scalar(posterior_ordinate, "y")
  }
  if (!is.finite(ordinate) || ordinate <= 0) {
    return("posterior ordinate is zero or non-finite")
  }

  value <- .iwmde_ordinate_scalar(posterior_ordinate, "evaluation_value")
  if (!is.finite(value)) {
    value <- .iwmde_ordinate_scalar(posterior_ordinate, "value")
  }
  if (is.finite(value)) {
    diagnostics[["evaluation_value"]] <- value
  }

  return(.iwmde_diagnostics_bf_failure_reason(diagnostics))
}


.iwmde_diagnostics_bf_failure_reason <- function(diagnostics) {

  estimator <- diagnostics[["estimator"]]
  if (length(estimator) != 1L) {
    return("missing qCMDE/IWMDE estimator diagnostics")
  }
  estimator <- as.character(estimator)
  if (!estimator %in% c("q_grid_cmde", "iwmde")) {
    return("unknown qCMDE/IWMDE estimator")
  }

  ordinate <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("ordinate", "bf_ordinate")
  )
  if (!is.na(ordinate) && (!is.finite(ordinate) || ordinate <= 0)) {
    return("posterior ordinate is zero or non-finite")
  }

  relative_mcse <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("relative_mcse", "bf_relative_mcse")
  )
  if (!is.finite(relative_mcse) || relative_mcse < 0 ||
      relative_mcse >= .iwmde_bf_max_relative_mcse()) {
    return(paste0(
      "relative MCSE is ",
      .iwmde_percent(relative_mcse),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_bf_max_relative_mcse()),
      ")"
    ))
  }

  finite_terms <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("finite_terms", "bf_finite_terms")
  )
  if (!is.finite(finite_terms) || finite_terms < 20) {
    return(paste0(
      "only ", .iwmde_count(finite_terms),
      " finite importance terms are available (minimum 20)"
    ))
  }

  ess <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("ess", "bf_ess")
  )
  if (!is.finite(ess) || ess < .iwmde_bf_min_ess()) {
    return(paste0(
      "effective sample size is ", .iwmde_count(ess),
      " (minimum ", .iwmde_count(.iwmde_bf_min_ess()), ")"
    ))
  }

  max_weight_share <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_weight_share", "bf_max_weight_share")
  )
  if (!is.finite(max_weight_share) ||
      max_weight_share >= .iwmde_bf_max_weight_share()) {
    return(paste0(
      "largest importance weight contributes ",
      .iwmde_percent(max_weight_share),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_bf_max_weight_share()),
      ")"
    ))
  }

  row_loss_reason <- .iwmde_diagnostics_row_loss_failure_reason(
    diagnostics = diagnostics,
    estimator   = estimator
  )
  if (!is.null(row_loss_reason)) {
    return(row_loss_reason)
  }

  if (identical(estimator, "q_grid_cmde")) {
    mass_reason <- .iwmde_diagnostics_mass_failure_reason(
      diagnostics = diagnostics,
      estimator   = estimator
    )
    if (!is.null(mass_reason)) {
      return(mass_reason)
    }
  } else if (identical(estimator, "iwmde")) {
    mass_reason <- .iwmde_diagnostics_mass_failure_reason(
      diagnostics = diagnostics,
      estimator   = estimator
    )
    if (!is.null(mass_reason)) {
      return(mass_reason)
    }
  }

  return(NULL)
}


.iwmde_diagnostics_density_failure_reason <- function(diagnostics) {

  if (is.null(diagnostics) || !is.list(diagnostics)) {
    return("missing qCMDE/IWMDE density diagnostics")
  }
  density_diagnostics <- .iwmde_density_diagnostics_as_ordinate(diagnostics)

  estimator <- diagnostics[["estimator"]]
  if (length(estimator) != 1L) {
    return("missing qCMDE/IWMDE estimator diagnostics")
  }
  estimator <- as.character(estimator)
  if (!estimator %in% c("q_grid_cmde", "iwmde")) {
    return("unknown qCMDE/IWMDE estimator")
  }

  reason <- .iwmde_diagnostics_density_sample_failure_reason(density_diagnostics)
  if (!is.null(reason)) {
    return(reason)
  }

  row_loss_reason <- .iwmde_diagnostics_row_loss_failure_reason(
    diagnostics = density_diagnostics,
    estimator   = estimator
  )
  if (!is.null(row_loss_reason)) {
    return(row_loss_reason)
  }

  mass_reason <- .iwmde_diagnostics_mass_failure_reason(
    diagnostics = density_diagnostics,
    estimator   = estimator
  )
  if (!is.null(mass_reason)) {
    return(mass_reason)
  }

  return(.iwmde_diagnostics_quadrature_failure_reason(density_diagnostics))
}


.iwmde_diagnostics_mass_failure_reason <- function(diagnostics, estimator) {

  normalization_error <- .iwmde_diagnostics_normalization_relative_error(
    diagnostics = diagnostics,
    estimator   = estimator
  )
  if (!is.finite(normalization_error)) {
    return("normalization diagnostics are unavailable")
  }

  fail_tolerance <- .iwmde_bf_mass_fail_tolerance(estimator)
  if (normalization_error > fail_tolerance) {
    return(paste0(
      .iwmde_estimator_label(estimator),
      .iwmde_diagnostics_normalization_error_phrase(estimator),
      .iwmde_percent(normalization_error),
      " (maximum allowed ",
      .iwmde_percent(fail_tolerance),
      ")"
    ))
  }

  return(NULL)
}


.iwmde_diagnostics_density_warning <- function(diagnostics) {

  if (is.null(diagnostics) || !is.list(diagnostics)) {
    return(character())
  }

  density_diagnostics <- .iwmde_density_diagnostics_as_ordinate(diagnostics)
  estimator <- diagnostics[["estimator"]]
  if (length(estimator) != 1L) {
    return(character())
  }
  estimator <- as.character(estimator)
  if (!estimator %in% c("q_grid_cmde", "iwmde")) {
    return(character())
  }

  warnings <- .iwmde_diagnostics_density_sample_warning(density_diagnostics)
  normalization_warning <- .iwmde_diagnostics_density_mass_warning(
    diagnostics = density_diagnostics,
    estimator   = estimator
  )
  quadrature_warning <- .iwmde_diagnostics_quadrature_warning(density_diagnostics)
  warnings <- c(warnings, normalization_warning, quadrature_warning)

  return(unique(warnings[nzchar(warnings)]))
}


.iwmde_diagnostics_density_sample_failure_reason <- function(diagnostics) {

  estimator_rows <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("n_estimator_rows", "n_normalized_rows", "n_evaluated_rows")
  )
  n_active_state_keys <- .iwmde_diagnostic_scalar(
    diagnostics,
    "n_active_state_keys"
  )
  if (is.finite(estimator_rows) &&
      is.finite(n_active_state_keys) &&
      n_active_state_keys > 1 &&
      estimator_rows < .iwmde_density_min_estimator_rows()) {
    return(paste0(
      "model-averaged density uses only ",
      .iwmde_count(estimator_rows),
      " estimator rows across ",
      .iwmde_count(n_active_state_keys),
      " active states (minimum ",
      .iwmde_count(.iwmde_density_min_estimator_rows()),
      ")"
    ))
  }

  relative_mcse <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_relative_mcse", "relative_mcse", "bf_relative_mcse")
  )
  if (!is.finite(relative_mcse) || relative_mcse < 0 ||
      relative_mcse >= .iwmde_density_max_relative_mcse()) {
    return(paste0(
      "density relative MCSE is ",
      .iwmde_percent(relative_mcse),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")"
    ))
  }

  ess <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("min_ess", "ess", "bf_ess")
  )
  min_ess <- .iwmde_density_min_ess(estimator_rows)
  if (!is.finite(ess) || ess < min_ess) {
    return(paste0(
      "density effective sample size is ", .iwmde_count(ess),
      " (minimum ", .iwmde_count(min_ess), ")"
    ))
  }

  max_weight_share <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_weight_share", "bf_max_weight_share")
  )
  if (!is.finite(max_weight_share) ||
      max_weight_share >= .iwmde_density_max_weight_share()) {
    return(paste0(
      "largest density importance weight contributes ",
      .iwmde_percent(max_weight_share),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_density_max_weight_share()),
      ")"
    ))
  }

  return(NULL)
}


.iwmde_diagnostics_density_sample_warning <- function(diagnostics) {

  warnings <- character()
  estimator_rows <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("n_estimator_rows", "n_normalized_rows", "n_evaluated_rows")
  )
  if (is.finite(estimator_rows) &&
      estimator_rows < .iwmde_density_warning_min_estimator_rows()) {
    warnings <- c(warnings, paste0(
      "Density uses only ",
      .iwmde_count(estimator_rows),
      " estimator rows",
      " (warning threshold ",
      .iwmde_count(.iwmde_density_warning_min_estimator_rows()),
      "; recommended minimum ",
      .iwmde_count(.iwmde_density_min_estimator_rows()),
      ")."
    ))
  }

  relative_mcse <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_relative_mcse", "relative_mcse", "bf_relative_mcse")
  )
  if (is.finite(relative_mcse) &&
      relative_mcse >= .iwmde_density_warning_relative_mcse() &&
      relative_mcse < .iwmde_density_max_relative_mcse()) {
    warnings <- c(warnings, paste0(
      "Density relative MCSE is ",
      .iwmde_percent(relative_mcse),
      " (warning threshold ",
      .iwmde_percent(.iwmde_density_warning_relative_mcse()),
      "; rejection threshold ",
      .iwmde_percent(.iwmde_density_max_relative_mcse()),
      ")."
    ))
  }

  ess <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("min_ess", "ess", "bf_ess")
  )
  min_ess <- .iwmde_density_min_ess(estimator_rows)
  warning_min_ess <- .iwmde_density_warning_min_ess(estimator_rows)
  if (is.finite(ess) &&
      ess >= min_ess &&
      ess < warning_min_ess) {
    warnings <- c(warnings, paste0(
      "Density effective sample size is ",
      .iwmde_count(ess),
      " (warning threshold ",
      .iwmde_count(warning_min_ess),
      "; rejection threshold ",
      .iwmde_count(min_ess),
      ")."
    ))
  }

  max_weight_share <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_weight_share", "bf_max_weight_share")
  )
  if (is.finite(max_weight_share) &&
      max_weight_share >= .iwmde_density_warning_weight_share() &&
      max_weight_share < .iwmde_density_max_weight_share()) {
    warnings <- c(warnings, paste0(
      "Largest density importance weight contributes ",
      .iwmde_percent(max_weight_share),
      " (warning threshold ",
      .iwmde_percent(.iwmde_density_warning_weight_share()),
      "; rejection threshold ",
      .iwmde_percent(.iwmde_density_max_weight_share()),
      ")."
    ))
  }

  return(warnings)
}


.iwmde_diagnostics_density_mass_warning <- function(diagnostics, estimator) {

  warnings <- character()
  row_loss <- .iwmde_diagnostics_row_loss_fraction(diagnostics)
  if (is.finite(row_loss)) {
    warning_tolerance <- .iwmde_bf_mass_warning_tolerance(estimator)
    fail_tolerance    <- .iwmde_bf_mass_fail_tolerance(estimator)
    if (row_loss > warning_tolerance &&
        row_loss <= fail_tolerance) {
      warnings <- c(warnings, paste0(
        .iwmde_estimator_label(estimator),
        " dropped ",
        .iwmde_percent(row_loss),
        " of target rows during density estimation",
        " (warning threshold ",
        .iwmde_percent(warning_tolerance),
        "; rejection threshold ",
        .iwmde_percent(fail_tolerance),
        ")."
      ))
    }
  }

  normalization_error <- .iwmde_diagnostics_normalization_relative_error(
    diagnostics = diagnostics,
    estimator   = estimator
  )
  warning_tolerance <- .iwmde_bf_mass_warning_tolerance(estimator)
  fail_tolerance    <- .iwmde_bf_mass_fail_tolerance(estimator)
  if (is.finite(normalization_error) &&
      normalization_error > warning_tolerance &&
      normalization_error <= fail_tolerance) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      .iwmde_diagnostics_normalization_error_phrase(estimator),
      .iwmde_percent(normalization_error),
      " (warning threshold ",
      .iwmde_percent(warning_tolerance),
      "; rejection threshold ",
      .iwmde_percent(fail_tolerance),
      ")."
    ))
  }

  return(warnings)
}


.iwmde_diagnostics_bf_warning <- function(diagnostics) {

  estimator <- diagnostics[["estimator"]]
  if (length(estimator) != 1L) {
    return(character())
  }
  estimator <- as.character(estimator)
  if (!estimator %in% c("q_grid_cmde", "iwmde")) {
    return(character())
  }

  warnings <- character()

  relative_mcse <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("relative_mcse", "bf_relative_mcse")
  )
  if (is.finite(relative_mcse) &&
      relative_mcse >= .iwmde_bf_warning_relative_mcse() &&
      relative_mcse < .iwmde_bf_max_relative_mcse()) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      " relative MCSE is ",
      .iwmde_percent(relative_mcse),
      " (warning threshold ",
      .iwmde_percent(.iwmde_bf_warning_relative_mcse()),
      "; BF rejection threshold ",
      .iwmde_percent(.iwmde_bf_max_relative_mcse()),
      ")."
    ))
  }

  finite_terms <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("finite_terms", "bf_finite_terms")
  )
  if (is.finite(finite_terms) &&
      finite_terms >= 20 &&
      finite_terms < .iwmde_bf_warning_min_finite_terms()) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      " uses only ",
      .iwmde_count(finite_terms),
      " finite importance terms",
      " (warning threshold ",
      .iwmde_count(.iwmde_bf_warning_min_finite_terms()),
      "; BF rejection threshold 20)."
    ))
  }

  ess <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("ess", "bf_ess")
  )
  if (is.finite(ess) &&
      ess >= .iwmde_bf_min_ess() &&
      ess < .iwmde_bf_warning_min_ess()) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      " effective sample size is ",
      .iwmde_count(ess),
      " (warning threshold ",
      .iwmde_count(.iwmde_bf_warning_min_ess()),
      "; BF rejection threshold ",
      .iwmde_count(.iwmde_bf_min_ess()),
      ")."
    ))
  }

  max_weight_share <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_weight_share", "bf_max_weight_share")
  )
  if (is.finite(max_weight_share) &&
      max_weight_share >= .iwmde_bf_warning_weight_share() &&
      max_weight_share < .iwmde_bf_max_weight_share()) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      " largest importance weight contributes ",
      .iwmde_percent(max_weight_share),
      " (warning threshold ",
      .iwmde_percent(.iwmde_bf_warning_weight_share()),
      "; BF rejection threshold ",
      .iwmde_percent(.iwmde_bf_max_weight_share()),
      ")."
    ))
  }

  row_loss <- .iwmde_diagnostics_row_loss_fraction(diagnostics)
  if (is.finite(row_loss)) {
    warning_tolerance <- .iwmde_bf_mass_warning_tolerance(estimator)
    fail_tolerance    <- .iwmde_bf_mass_fail_tolerance(estimator)
    if (row_loss > warning_tolerance &&
        row_loss <= fail_tolerance) {
      warnings <- c(warnings, paste0(
        .iwmde_estimator_label(estimator),
        " dropped ",
        .iwmde_percent(row_loss),
        " of target rows during density estimation",
        " (warning threshold ",
        .iwmde_percent(warning_tolerance),
        "; BF rejection threshold ",
        .iwmde_percent(fail_tolerance),
        ")."
      ))
    }
  }

  normalization_error <- .iwmde_diagnostics_normalization_relative_error(
    diagnostics = diagnostics,
    estimator   = estimator
  )
  warning_tolerance <- .iwmde_bf_mass_warning_tolerance(estimator)
  fail_tolerance    <- .iwmde_bf_mass_fail_tolerance(estimator)
  if (is.finite(normalization_error) &&
      normalization_error > warning_tolerance &&
      normalization_error <= fail_tolerance) {
    warnings <- c(warnings, paste0(
      .iwmde_estimator_label(estimator),
      .iwmde_diagnostics_normalization_error_phrase(estimator),
      .iwmde_percent(normalization_error),
      " (warning threshold ",
      .iwmde_percent(warning_tolerance),
      "; BF rejection threshold ",
      .iwmde_percent(fail_tolerance),
      ")."
    ))
  }

  return(warnings)
}


.iwmde_density_diagnostics_as_ordinate <- function(diagnostics) {

  out <- diagnostics
  out[["relative_mcse"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_relative_mcse", "relative_mcse")
  )
  out[["finite_terms"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("min_finite_terms", "finite_terms")
  )
  out[["ess"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("min_ess", "ess")
  )
  out[["max_weight_share"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_weight_share", "bf_max_weight_share")
  )
  out[["ordinate_relative_change"]] <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_ordinate_relative_change", "ordinate_relative_change")
  )

  return(out)
}


.iwmde_diagnostics_quadrature_failure_reason <- function(diagnostics) {

  quadrature_change <- .iwmde_diagnostics_quadrature_relative_change(diagnostics)
  if (!is.finite(quadrature_change)) {
    return(NULL)
  }
  if (quadrature_change > .iwmde_quadrature_fail_tolerance()) {
    return(paste0(
      "adaptive quadrature sensitivity is ",
      .iwmde_percent(quadrature_change),
      " (maximum allowed ",
      .iwmde_percent(.iwmde_quadrature_fail_tolerance()),
      ")"
    ))
  }

  return(NULL)
}


.iwmde_diagnostics_quadrature_warning <- function(diagnostics) {

  quadrature_change <- .iwmde_diagnostics_quadrature_relative_change(diagnostics)
  if (!is.finite(quadrature_change) ||
      quadrature_change <= .iwmde_quadrature_warning_tolerance() ||
      quadrature_change > .iwmde_quadrature_fail_tolerance()) {
    return(character())
  }

  return(paste0(
    "Adaptive quadrature sensitivity is ",
    .iwmde_percent(quadrature_change),
    " (warning threshold ",
    .iwmde_percent(.iwmde_quadrature_warning_tolerance()),
    "; rejection threshold ",
    .iwmde_percent(.iwmde_quadrature_fail_tolerance()),
    ")."
  ))
}


.iwmde_diagnostics_quadrature_relative_change <- function(diagnostics) {

  quadrature_change <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("max_quadrature_relative_change", "quadrature_relative_change")
  )
  if (!is.finite(quadrature_change) || quadrature_change < 0) {
    return(NA_real_)
  }

  return(quadrature_change)
}


.iwmde_diagnostics_row_loss_failure_reason <- function(diagnostics, estimator) {

  row_loss <- .iwmde_diagnostics_row_loss_fraction(diagnostics)
  if (!is.finite(row_loss)) {
    return(NULL)
  }

  fail_tolerance <- .iwmde_bf_mass_fail_tolerance(estimator)
  if (row_loss > fail_tolerance) {
    return(paste0(
      .iwmde_estimator_label(estimator),
      " dropped ",
      .iwmde_percent(row_loss),
      " of target rows during density estimation",
      " (maximum allowed ",
      .iwmde_percent(fail_tolerance),
      ")"
    ))
  }

  return(NULL)
}


.iwmde_diagnostics_row_loss_fraction <- function(diagnostics) {

  row_loss <- .iwmde_diagnostic_scalar(diagnostics, "row_drop_fraction")
  if (is.finite(row_loss) && row_loss >= 0) {
    return(min(1, row_loss))
  }

  n_candidate_rows <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("n_candidate_rows", "n_active")
  )
  n_normalized_rows <- .iwmde_diagnostic_scalar(
    diagnostics,
    "n_normalized_rows"
  )
  if (!is.finite(n_candidate_rows) || n_candidate_rows <= 0 ||
      !is.finite(n_normalized_rows) || n_normalized_rows < 0) {
    return(NA_real_)
  }

  return(.iwmde_row_drop_fraction(
    n_candidate_rows  = n_candidate_rows,
    n_normalized_rows = n_normalized_rows
  ))
}


.iwmde_diagnostics_normalization_relative_error <- function(diagnostics,
                                                            estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return(.iwmde_diagnostics_qcmde_ordinate_relative_change(diagnostics))
  }

  return(.iwmde_diagnostics_mass_relative_error(diagnostics))
}


.iwmde_diagnostics_qcmde_ordinate_relative_change <- function(diagnostics) {

  ordinate_change <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    c("ordinate_relative_change", "bf_ordinate_relative_change")
  )
  if (!is.finite(ordinate_change) || ordinate_change < 0) {
    return(NA_real_)
  }

  return(ordinate_change)
}


.iwmde_diagnostics_normalization_error_phrase <- function(estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return(" posterior ordinate changes by ")
  }

  return(" normalization mass differs from active posterior mass by ")
}


.iwmde_diagnostics_mass_relative_error <- function(diagnostics) {

  normalization_integral <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    "normalization_integral"
  )
  active_mass <- .iwmde_diagnostic_scalar_any(
    diagnostics,
    "active_mass"
  )
  if (!is.finite(normalization_integral) || !is.finite(active_mass) ||
      active_mass <= 0) {
    return(NA_real_)
  }

  return(abs(normalization_integral / active_mass - 1))
}


.iwmde_ordinate_scalar <- function(posterior_ordinate, name) {

  value <- posterior_ordinate[[name]]
  if (is.null(value)) {
    return(NA_real_)
  }

  value <- suppressWarnings(as.numeric(value))[1]
  if (length(value) == 0L) {
    return(NA_real_)
  }

  return(value)
}


.iwmde_bf_max_relative_mcse <- function() {

  return(1)
}


.iwmde_bf_warning_relative_mcse <- function() {

  return(.25)
}


.iwmde_bf_min_ess <- function() {

  return(4)
}


.iwmde_bf_warning_min_ess <- function() {

  return(20)
}


.iwmde_bf_max_weight_share <- function() {

  return(.80)
}


.iwmde_bf_warning_weight_share <- function() {

  return(.50)
}


.iwmde_bf_warning_min_finite_terms <- function() {

  return(50)
}


.iwmde_bf_mass_warning_tolerance <- function(estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return(.025)
  }
  if (identical(estimator, "iwmde")) {
    return(.05)
  }

  return(Inf)
}


.iwmde_bf_mass_fail_tolerance <- function(estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return(.05)
  }
  if (identical(estimator, "iwmde")) {
    return(.10)
  }

  return(0)
}


.iwmde_density_min_estimator_rows <- function() {

  return(300)
}


.iwmde_density_warning_min_estimator_rows <- function() {

  return(500)
}


.iwmde_density_max_relative_mcse <- function() {

  return(.25)
}


.iwmde_density_warning_relative_mcse <- function() {

  return(.10)
}


.iwmde_density_min_ess <- function(estimator_rows = NA_real_) {

  estimator_rows <- as.numeric(estimator_rows)[1]
  if (is.finite(estimator_rows) && estimator_rows > 0) {
    return(min(50, max(4, .25 * estimator_rows)))
  }

  return(50)
}


.iwmde_density_warning_min_ess <- function(estimator_rows = NA_real_) {

  estimator_rows <- as.numeric(estimator_rows)[1]
  if (is.finite(estimator_rows) && estimator_rows > 0) {
    return(min(100, max(20, .50 * estimator_rows)))
  }

  return(100)
}


.iwmde_density_max_weight_share <- function() {

  return(.25)
}


.iwmde_density_warning_weight_share <- function() {

  return(.10)
}


.iwmde_quadrature_warning_tolerance <- function() {

  return(.025)
}


.iwmde_quadrature_fail_tolerance <- function() {

  return(.05)
}


.iwmde_diagnostic_scalar <- function(diagnostics, name) {

  value <- diagnostics[[name]]
  if (is.null(value)) {
    return(NA_real_)
  }
  if (is.list(value)) {
    value <- unlist(value, use.names = FALSE)
  }

  value <- suppressWarnings(as.numeric(value))[1]
  if (length(value) == 0L) {
    return(NA_real_)
  }

  return(value)
}


.iwmde_diagnostic_scalar_any <- function(diagnostics, names) {

  for (name in names) {
    value <- .iwmde_diagnostic_scalar(diagnostics, name)
    if (!is.na(value)) {
      return(value)
    }
  }

  return(NA_real_)
}


.iwmde_estimator_label <- function(estimator) {

  if (identical(estimator, "q_grid_cmde")) {
    return("qCMDE")
  }
  if (identical(estimator, "iwmde")) {
    return("IWMDE")
  }

  return(as.character(estimator)[1])
}


.iwmde_percent <- function(value) {

  if (!is.finite(value)) {
    return("NA")
  }

  paste0(formatC(100 * value, digits = 3, format = "fg"), "%")
}


.iwmde_count <- function(value) {

  if (!is.finite(value)) {
    return("NA")
  }

  formatC(value, digits = 3, format = "fg")
}


.iwmde_posterior_ordinate_warnings <- function(posterior_ordinate) {

  entries <- .iwmde_posterior_ordinate_entries(posterior_ordinate)
  warnings <- unlist(lapply(entries, function(entry) {
    diagnostics <- entry[["diagnostics"]]
    if (is.null(diagnostics) || !is.list(diagnostics)) {
      return(character())
    }
    explicit <- diagnostics[["warning"]]
    if (!is.null(explicit)) {
      explicit <- as.character(explicit)
      explicit <- explicit[nzchar(explicit)]
      if (length(explicit) > 0L) {
        return(explicit)
      }
    }
    .iwmde_diagnostics_bf_warning(diagnostics)
  }), use.names = FALSE)

  unique(warnings[nzchar(warnings)])
}


.iwmde_posterior_ordinate_failure_reasons <- function(posterior_ordinate) {

  entries <- .iwmde_posterior_ordinate_entries(posterior_ordinate)
  reasons <- unlist(lapply(entries, function(entry) {
    reason <- .iwmde_posterior_ordinate_bf_failure_reason(entry)
    if (is.null(reason)) {
      return(character())
    }
    reason
  }), use.names = FALSE)

  unique(reasons[nzchar(reasons)])
}


.iwmde_posterior_ordinate_entries <- function(posterior_ordinate) {

  if (is.null(posterior_ordinate)) {
    return(list())
  }
  if (is.list(posterior_ordinate) &&
      !is.data.frame(posterior_ordinate[["ordinates"]]) &&
      is.list(posterior_ordinate[["ordinates"]]) &&
      is.null(posterior_ordinate[["ordinates"]][["x"]]) &&
      is.null(posterior_ordinate[["ordinates"]][["value"]]) &&
      is.null(posterior_ordinate[["ordinates"]][["null_hypothesis"]])) {
    return(posterior_ordinate[["ordinates"]])
  }

  return(list(posterior_ordinate))
}


.iwmde_bf_append_warning <- function(bf, posterior_ordinate) {

  warnings <- .iwmde_posterior_ordinate_warnings(posterior_ordinate)
  if (length(warnings) == 0L) {
    return(bf)
  }

  attr(bf, "warnings") <- unique(c(
    attr(bf, "warnings", exact = TRUE),
    warnings
  ))

  return(bf)
}


.iwmde_density_bf_diagnostics <- function(density, values) {

  out <- list(
    bf_value            = NA_real_,
    bf_evaluation_value = NA_real_,
    bf_included         = FALSE,
    bf_grid_index       = NA_integer_,
    bf_ordinate         = NA_real_,
    bf_pilot_ordinate   = NA_real_,
    bf_validation_ordinate = NA_real_,
    bf_ordinate_relative_change = NA_real_,
    bf_ordinate_log_change      = NA_real_,
    bf_pilot_ordinate_relative_change = NA_real_,
    bf_pilot_ordinate_log_change      = NA_real_,
    bf_mcse             = NA_real_,
    bf_relative_mcse    = NA_real_,
    bf_error_percent    = NA_real_,
    bf_finite_terms     = NA_integer_,
    bf_ess              = NA_real_,
    bf_max_weight_share = NA_real_,
    bf_max_log_ratio    = NA_real_
  )

  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(out)
  }
  out[["bf_value"]] <- values[1]

  x <- as.numeric(density[["x"]])
  if (length(x) == 0L) {
    return(out)
  }
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(out[["bf_value"]]))
  index <- which(abs(x - out[["bf_value"]]) <= tolerance)
  if (length(index) == 0L) {
    return(out)
  }
  index <- index[1]

  out[["bf_included"]]         <- TRUE
  out[["bf_grid_index"]]       <- index
  out[["bf_evaluation_value"]] <- .iwmde_density_evaluation_value(density, index)
  out[["bf_ordinate"]]         <- .iwmde_density_index_value(density, "y", index)
  out[["bf_pilot_ordinate"]]   <- .iwmde_density_index_value(density, "pilot_y", index)
  out[["bf_validation_ordinate"]] <-
    .iwmde_density_index_value(density, "validation_y", index)
  out[["bf_ordinate_relative_change"]] <-
    .iwmde_density_index_value(density, "ordinate_relative_change", index)
  out[["bf_ordinate_log_change"]] <-
    .iwmde_density_index_value(density, "ordinate_log_change", index)
  out[["bf_pilot_ordinate_relative_change"]] <-
    .iwmde_density_index_value(density, "pilot_ordinate_relative_change", index)
  out[["bf_pilot_ordinate_log_change"]] <-
    .iwmde_density_index_value(density, "pilot_ordinate_log_change", index)
  out[["bf_mcse"]]             <- .iwmde_density_index_value(density, "mcse", index)
  out[["bf_relative_mcse"]]    <- .iwmde_density_index_value(density, "relative_mcse", index)
  out[["bf_finite_terms"]]     <- as.integer(.iwmde_density_index_value(density, "finite_terms", index))
  out[["bf_ess"]]              <- .iwmde_density_index_value(density, "ess", index)
  out[["bf_max_weight_share"]] <- .iwmde_density_index_value(density, "max_weight_share", index)
  out[["bf_max_log_ratio"]]    <- .iwmde_density_index_value(density, "max_log_ratio", index)

  if (is.finite(out[["bf_relative_mcse"]])) {
    out[["bf_error_percent"]] <- 100 * out[["bf_relative_mcse"]]
  }

  return(out)
}


.iwmde_density_evaluation_value <- function(density, index) {

  value <- .iwmde_density_index_value(density, "evaluation_x", index)
  if (is.finite(value)) {
    return(value)
  }

  return(.iwmde_density_index_value(density, "x", index))
}


.iwmde_density_index_value <- function(density, name, index) {

  value <- density[[name]]
  if (is.null(value) || length(value) < index) {
    return(NA_real_)
  }

  value <- suppressWarnings(as.numeric(value[index]))[1]
  if (length(value) == 0L) {
    return(NA_real_)
  }

  return(value)
}


.iwmde_parameter_ordinate_diagnostic <- function(context, parameter, values,
                                                 max_samples,
                                                 normalization_points,
                                                 normalization_prob,
                                                 method = c("q_grid_cmde", "iwmde"),
                                                 parameter_spec = NULL,
                                                 diagnostic_cache = NULL) {

  method <- .iwmde_normalize_method(method)
  plan <- .iwmde_plan(
    context         = context,
    parameter       = parameter,
    density_method  = if (identical(method, "q_grid_cmde")) "qCMDE" else "IWMDE",
    density_control = list(
      n_points             = max(20L, length(values)),
      max_samples          = max_samples,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      display_grid         = "ordinate"
    ),
    outputs         = "ordinate",
    values          = values,
    parameter_spec  = parameter_spec,
    metadata        = NULL
  )

  return(.iwmde_execute_plan_diagnostic(
    context          = context,
    plan             = plan,
    output           = "ordinate",
    execution_cache  = new.env(parent = emptyenv()),
    diagnostic_cache = diagnostic_cache
  ))
}


.iwmde_ordinate_interior_values <- function(values, support, xlim) {

  out   <- values
  width <- diff(xlim)
  if (!is.finite(width) || width <= 0) {
    width <- 1
  }
  eps <- sqrt(.Machine$double.eps) * max(1, width)

  if (is.finite(support[1])) {
    out[out <= support[1]] <- support[1] + eps
  }
  if (is.finite(support[2])) {
    out[out >= support[2]] <- support[2] - eps
  }

  return(out)
}


.iwmde_normalize_display_grid <- function(display_grid) {

  return(match.arg(display_grid, c("adaptive", "uniform")))
}


.iwmde_parameter_diagnostic <- function(context, parameter, n_points,
                                        max_samples, normalization_points,
                                        normalization_prob,
                                        method = c("q_grid_cmde", "iwmde"),
                                        display_grid_method = "adaptive",
                                        include_values = NULL,
                                        parameter_spec = NULL,
                                        diagnostic_cache = NULL) {

  method <- .iwmde_normalize_method(method)
  plan <- .iwmde_plan(
    context         = context,
    parameter       = parameter,
    density_method  = if (identical(method, "q_grid_cmde")) "qCMDE" else "IWMDE",
    density_control = list(
      n_points             = n_points,
      max_samples          = max_samples,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      display_grid         = display_grid_method
    ),
    outputs         = "density",
    values          = include_values,
    parameter_spec  = parameter_spec,
    metadata        = NULL
  )

  return(.iwmde_execute_plan_diagnostic(
    context          = context,
    plan             = plan,
    output           = "density",
    execution_cache  = new.env(parent = emptyenv()),
    diagnostic_cache = diagnostic_cache
  ))
}


.iwmde_diagnostic_cache <- function() {

  cache <- new.env(parent = emptyenv())

  return(cache)
}


.iwmde_target_key <- function(parameter, parameter_spec) {

  condition_key <- .iwmde_parameter_condition_key(parameter_spec)

  if (identical(parameter_spec[["type"]], "linear")) {
    weights <- parameter_spec[["weights"]]
    weights <- weights[order(names(weights))]
    parts   <- paste0(
      names(weights),
      "=",
      .iwmde_key_number(weights)
    )

    return(paste(c("linear", parts, condition_key), collapse = "|"))
  }

  if (identical(parameter_spec[["type"]], "primitive")) {
    return(paste(c("primitive", parameter, condition_key), collapse = "|"))
  }

  if (identical(parameter_spec[["status"]], "unsupported")) {
    return(paste0("unsupported|", parameter))
  }

  return(paste0("primitive|", parameter))
}


.iwmde_parameter_condition_key <- function(parameter_spec) {

  condition_key <- parameter_spec[["condition_key"]]
  if (!is.null(condition_key) && length(condition_key) > 0L) {
    condition_key <- as.character(condition_key[[1L]])
    if (!is.na(condition_key) && nzchar(condition_key)) {
      return(condition_key)
    }
  }

  conditional <- parameter_spec[["conditional"]]
  if (is.null(conditional) || length(conditional) == 0L) {
    return(NULL)
  }

  conditional <- sort(unique(as.character(conditional)))
  rule        <- parameter_spec[["conditional_rule"]]
  if (is.null(rule) || length(rule) == 0L) {
    rule <- "OR"
  }

  return(paste0("conditional=", rule, ":", paste(conditional, collapse = ",")))
}


.iwmde_key_number <- function(x) {

  x <- zapsmall(as.numeric(x), digits = 14)
  x[abs(x) < sqrt(.Machine$double.eps)] <- 0

  return(formatC(x, digits = 15, format = "fg", flag = "#"))
}


.iwmde_cache_has <- function(cache, key) {

  return(!is.null(cache) && exists(key, envir = cache, inherits = FALSE))
}


.iwmde_cache_get <- function(cache, key) {

  return(get(key, envir = cache, inherits = FALSE))
}


.iwmde_cache_set <- function(cache, key, value) {

  if (!is.null(cache)) {
    assign(key, value, envir = cache)
  }

  return(invisible(value))
}


.iwmde_relabel_diagnostic <- function(diagnostic, parameter) {

  diagnostic[["parameter"]] <- parameter

  return(diagnostic)
}


.iwmde_select_active_rows <- function(rows, max_samples, context = NULL) {

  if (length(rows) <= max_samples) {
    return(rows)
  }

  selected <- NULL
  if (!is.null(context) && length(context[["indicator_names"]]) > 0L) {
    samples <- context[["posterior_samples"]]
    group <- vapply(rows, function(row) {
      .iwmde_active_key(context, samples[row, ])
    }, character(1))
    selected <- .thin_sample_rows_by_group(group, max_samples)
  }
  if (is.null(selected)) {
    selected <- .thin_sample_rows(length(rows), max_samples)
  }

  return(rows[selected])
}


.iwmde_marginal_means_specs <- function(marginal_means_object, parameter,
                                        type, levels) {

  if (is.null(parameter)) {
    selected <- lapply(marginal_means_object[["parameters"]], function(parameter) {
      .marginal_means_select_parameter(marginal_means_object, parameter)
    })
  } else {
    selected <- lapply(parameter, function(parameter) {
      .marginal_means_select_parameter(marginal_means_object, parameter)
    })
  }

  samples <- marginal_means_object[["inference"]][[type]]
  specs   <- list()

  for (selected_parameter in selected) {
    parameter_name <- selected_parameter[["parameter"]]
    if (!parameter_name %in% names(samples)) {
      next
    }

    level_samples <- samples[[parameter_name]]
    if (!is.list(level_samples)) {
      level_samples <- list(parameter_name = level_samples)
      names(level_samples) <- parameter_name
    }
    keep_levels <- names(level_samples)
    if (!is.null(levels)) {
      keep_levels <- intersect(keep_levels, levels)
    }

    for (level in keep_levels) {
      level_sample       <- level_samples[[level]]
      condition_metadata <- .iwmde_sample_condition_metadata(
        samples                       = level_sample,
        include_prior_density_context = TRUE,
        require_child_condition       = identical(type, "conditional")
      )
      label <- paste0(selected_parameter[["label"]], ": ", level)
      key   <- label
      specs[[key]] <- c(
        list(
          type      = "linear",
          label     = label,
          parameter = parameter_name,
          level     = level,
          weights   = attr(level_sample, "linear_weights")
        ),
        condition_metadata,
        list(
          source           = "marginal_means",
          selected         = selected_parameter,
          marginal_type    = type,
          marginal_samples = level_sample
        )
      )
    }
  }

  return(specs)
}
