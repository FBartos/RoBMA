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
#' @details Weightfunction \code{omega}, \code{log_omega}, and cumulative
#' \code{eta} coordinates are reported as unsupported until a joint replacement
#' map is implemented. For spike-mixture parameters, KDE and histogram layers
#' use the continuous active samples scaled by their posterior mass; point
#' masses are drawn separately.
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
#' @param normalization_points number of internal-scale grid points used for
#' row-wise conditional normalization. Defaults to \code{max(50, n_points)}.
#' @param normalization_prob posterior central probability used to construct
#' the normalization grid before forcing it to include the displayed range.
#' @param density_method density estimator. \code{"qCMDE"} uses
#' row-normalized q-grid conditional densities. \code{"IWMDE"} uses
#' Chen-style moment-matched weights: conditional normal weights for
#' unconstrained targets and power/exponential weights for bounded or
#' one-sided targets.
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
plot_iwmde_diagnostics <- function(object, parameters = NULL, n_points = 150,
                                   max_samples = 400,
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
  method <- .density_method_iwmde_estimator(density_method)
  display_grid <- .iwmde_normalize_display_grid(display_grid)

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
  diagnostic_cache <- .iwmde_diagnostic_cache()

  out <- lapply(parameters, function(parameter) {
    .iwmde_parameter_diagnostic(
      context     = context,
      parameter   = parameter,
      n_points    = n_points,
      max_samples = max_samples,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      method               = method,
      display_grid_method  = display_grid,
      diagnostic_cache     = diagnostic_cache
    )
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
                                                  n_points = 150,
                                                  max_samples = 400,
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
  method <- .density_method_iwmde_estimator(density_method)
  display_grid <- .iwmde_normalize_display_grid(display_grid)

  if (is.null(marginal_means_object)) {
    marginal_means_object <- marginal_means(
      object    = object,
      n_samples = n_marginal_samples
    )
  } else {
    .iwmde_validate_marginal_means_object(object, marginal_means_object)
  }
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
  diagnostic_cache <- .iwmde_diagnostic_cache()

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
    .iwmde_parameter_diagnostic(
      context              = context,
      parameter            = spec[["label"]],
      n_points             = n_points,
      max_samples          = max_samples,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob,
      method               = method,
      parameter_spec       = spec,
      display_grid_method  = display_grid,
      diagnostic_cache     = diagnostic_cache
    )
  })
  names(out) <- names(specs)
  class(out) <- c("iwmde_diagnostics", "list")

  if (as_data || !plot) {
    return(out)
  }

  do.call(.plot_iwmde_diagnostics_base, c(list(x = out), dots_raw))
  return(invisible(out))
}


.iwmde_validate_marginal_means_object <- function(object,
                                                  marginal_means_object) {

  if (!inherits(marginal_means_object, "marginal_means.brma")) {
    return(invisible(FALSE))
  }

  signature <- marginal_means_object[["iwmde_signature"]]
  if (is.null(signature)) {
    stop("'marginal_means_object' does not contain an IWMDE fit signature. ",
         "Omit it or recompute it with the current RoBMA version.",
         call. = FALSE)
  }
  if (!identical(signature, .iwmde_fit_signature(object))) {
    stop("'marginal_means_object' was not computed from 'object'.",
         call. = FALSE)
  }

  return(invisible(TRUE))
}


.iwmde_fit_signature <- function(object) {

  samples <- as.matrix(.get_posterior_samples(object[["fit"]]))
  means   <- vapply(seq_len(ncol(samples)), function(i) {
    .iwmde_signature_stat(samples[, i], mean)
  }, numeric(1))
  sds <- vapply(seq_len(ncol(samples)), function(i) {
    .iwmde_signature_stat(samples[, i], stats::sd)
  }, numeric(1))

  return(list(
    class             = class(object),
    posterior_dim     = dim(samples),
    posterior_columns = colnames(samples),
    posterior_means   = means,
    posterior_sds     = sds,
    outcome_type      = .data_outcome_type(object[["data"]]),
    is_mods           = .is_data_mods(object[["data"]]),
    is_scale          = .is_data_scale(object[["data"]]),
    is_multilevel     = .is_data_multilevel(object[["data"]]),
    is_weights        = .is_data_weights(object[["data"]])
  ))
}


.iwmde_signature_stat <- function(x, fun) {

  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }

  return(round(fun(x), 12))
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

  # These are constrained/derived publication-weight coordinates. IWMDE needs a
  # coherent replacement map that recomputes their companion coordinates.
  unsupported_prefix <- c("omega", "log_omega", "eta")
  keep <- !vapply(parameters, function(parameter) {
    any(parameter == unsupported_prefix) ||
      any(startsWith(parameter, paste0(unsupported_prefix, "[")))
  }, logical(1))

  return(parameters[keep])
}


.iwmde_normalize_method <- function(method) {

  return(match.arg(method, c("q_grid_cmde", "iwmde")))
}


.density_method_normalize <- function(density_method) {

  return(match.arg(density_method, c("KDE", "qCMDE", "IWMDE")))
}


.density_method_normalize_precomputed <- function(density_method) {

  density_method <- match.arg(density_method, c("qCMDE", "IWMDE"))

  return(density_method)
}


.density_method_uses_precomputed <- function(density_method) {

  density_method <- .density_method_normalize(density_method)

  return(!identical(density_method, "KDE"))
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


.density_control_normalize <- function(density_method, density_control = NULL) {

  density_method <- .density_method_normalize(density_method)
  allowed_names  <- c(
    "n_points",
    "max_samples",
    "normalization_points",
    "normalization_prob",
    "display_grid"
  )
  qcmde_only <- c("normalization_points", "normalization_prob")
  defaults <- list(
    n_points             = 150L,
    max_samples          = 400L,
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
  if (identical(density_method, "KDE")) {
    stop("'density_control' is only used when 'density_method' is ",
         "'qCMDE' or 'IWMDE'.", call. = FALSE)
  }

  unavailable_names <- intersect(control_names, qcmde_only)
  if (identical(density_method, "IWMDE") && length(unavailable_names) > 0L) {
    stop("'density_control' setting(s) ",
         paste0("'", unavailable_names, "'", collapse = ", "),
         " are only available when 'density_method = \"qCMDE\"'.",
         call. = FALSE)
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


.iwmde_posterior_density_attribute <- function(diagnostic, density_method) {

  if (!identical(diagnostic[["status"]], "ok") ||
      is.null(diagnostic[["iwmde"]])) {
    return(NULL)
  }

  density <- diagnostic[["iwmde"]]
  out <- list(
    status         = "ok",
    x              = density[["x"]],
    y              = density[["y"]],
    method         = density[["estimator"]],
    density_method = .density_method_normalize(density_method),
    diagnostics    = diagnostic[["diagnostics"]],
    point_masses   = diagnostic[["point_masses"]]
  )
  class(out) <- c("RoBMA_posterior_density", "list")

  return(out)
}


.iwmde_posterior_ordinate_attribute <- function(diagnostic, density_method) {

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

  out <- list(
    status         = "ok",
    value          = diagnostics[["bf_value"]],
    ordinate       = diagnostics[["bf_ordinate"]],
    method         = diagnostics[["estimator"]],
    density_method = .density_method_normalize(density_method),
    diagnostics    = list(
      mcse                   = diagnostics[["bf_mcse"]],
      relative_mcse          = diagnostics[["bf_relative_mcse"]],
      BF_error_percent       = diagnostics[["bf_error_percent"]],
      finite_terms           = diagnostics[["bf_finite_terms"]],
      ess                    = diagnostics[["bf_ess"]],
      max_weight_share       = diagnostics[["bf_max_weight_share"]],
      max_log_ratio          = diagnostics[["bf_max_log_ratio"]],
      active_mass            = diagnostics[["active_mass"]],
      normalization_integral = diagnostics[["normalization_integral"]],
      normalization_range    = diagnostics[["normalization_range"]],
      estimator              = diagnostics[["estimator"]],
      weight_method          = diagnostics[["weight_method"]]
    )
  )
  class(out) <- c("RoBMA_posterior_ordinate", "list")

  if (!.iwmde_posterior_ordinate_supports_bf(out)) {
    return(NULL)
  }

  return(out)
}


.iwmde_posterior_ordinate_supports_bf <- function(posterior_ordinate) {

  if (is.null(posterior_ordinate) || !is.list(posterior_ordinate)) {
    return(FALSE)
  }

  diagnostics <- posterior_ordinate[["diagnostics"]]
  if (is.null(diagnostics) || !is.list(diagnostics)) {
    return(FALSE)
  }

  estimator <- diagnostics[["estimator"]]
  if (length(estimator) != 1L) {
    return(FALSE)
  }
  estimator <- as.character(estimator)
  if (!estimator %in% c("q_grid_cmde", "iwmde")) {
    return(FALSE)
  }

  ordinate         <- .iwmde_ordinate_scalar(posterior_ordinate, "ordinate")
  relative_mcse    <- .iwmde_diagnostic_scalar(diagnostics, "relative_mcse")
  finite_terms     <- .iwmde_diagnostic_scalar(diagnostics, "finite_terms")
  ess              <- .iwmde_diagnostic_scalar(diagnostics, "ess")
  max_weight_share <- .iwmde_diagnostic_scalar(diagnostics, "max_weight_share")
  if (!is.finite(ordinate) || ordinate <= 0 ||
      !is.finite(relative_mcse) || relative_mcse < 0 ||
      relative_mcse >= .iwmde_bf_max_relative_mcse() ||
      !is.finite(finite_terms) || finite_terms < 20 ||
      !is.finite(ess) || ess < .iwmde_bf_min_ess() ||
      !is.finite(max_weight_share) ||
      max_weight_share >= .iwmde_bf_max_weight_share()) {
    return(FALSE)
  }

  if (identical(estimator, "q_grid_cmde")) {
    normalization_integral <- .iwmde_diagnostic_scalar(diagnostics, "normalization_integral")
    active_mass            <- .iwmde_diagnostic_scalar(diagnostics, "active_mass")
    if (!is.finite(normalization_integral) || !is.finite(active_mass) ||
        active_mass <= 0 ||
        abs(normalization_integral / active_mass - 1) >
          .iwmde_bf_mass_tolerance()) {
      return(FALSE)
    }
    normalization_range <- suppressWarnings(as.numeric(diagnostics[["normalization_range"]]))
    if (length(normalization_range) != 2L ||
        !all(is.finite(normalization_range))) {
      return(FALSE)
    }
    value <- .iwmde_ordinate_scalar(posterior_ordinate, "value")
    tolerance <- sqrt(.Machine$double.eps) * max(1, abs(value))
    if (!is.finite(value) ||
        value < min(normalization_range) - tolerance ||
        value > max(normalization_range) + tolerance) {
      return(FALSE)
    }
  }

  return(TRUE)
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


.iwmde_bf_min_ess <- function() {

  return(4)
}


.iwmde_bf_max_weight_share <- function() {

  return(.80)
}


.iwmde_bf_mass_tolerance <- function() {

  return(.25)
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


.iwmde_density_bf_diagnostics <- function(density, values) {

  out <- list(
    bf_value            = NA_real_,
    bf_included         = FALSE,
    bf_grid_index       = NA_integer_,
    bf_ordinate         = NA_real_,
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
  out[["bf_ordinate"]]         <- .iwmde_density_index_value(density, "y", index)
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

  method         <- .iwmde_normalize_method(method)
  values         <- as.numeric(values)
  values         <- values[is.finite(values)]
  parameter_spec <- .iwmde_parameter_spec(context, parameter, parameter_spec)
  target_key     <- .iwmde_target_key(parameter, parameter_spec)
  cache_key      <- .iwmde_diagnostic_cache_key(
    target_key            = paste0("ordinate|", target_key),
    n_points              = length(values),
    max_samples           = max_samples,
    normalization_points  = normalization_points,
    normalization_prob    = normalization_prob,
    method                = method,
    display_grid_method   = "ordinate",
    include_values        = values
  )
  if (.iwmde_cache_has(diagnostic_cache, cache_key)) {
    return(.iwmde_relabel_diagnostic(
      diagnostic = .iwmde_cache_get(diagnostic_cache, cache_key),
      parameter  = parameter
    ))
  }

  if (length(values) == 0L) {
    return(.iwmde_unsupported(parameter, "no finite ordinate values"))
  }
  if (identical(parameter_spec[["status"]], "unsupported")) {
    out <- .iwmde_unsupported(parameter, parameter_spec[["reason"]])
    .iwmde_cache_set(diagnostic_cache, cache_key, out)
    return(out)
  }

  posterior_values <- .iwmde_parameter_values(context, parameter, parameter_spec)
  finite           <- is.finite(posterior_values)
  rows             <- .iwmde_parameter_condition_rows(context, parameter_spec)
  finite_rows      <- finite & rows
  if (!any(finite_rows)) {
    return(.iwmde_unsupported(parameter, "posterior samples are not finite"))
  }

  component <- .iwmde_parameter_components(context, parameter, parameter_spec)
  component <- .iwmde_restrict_parameter_component(component, finite_rows)
  if (!any(component[["active"]])) {
    return(.iwmde_unsupported(
      parameter,
      "parameter has no continuous active posterior component"
    ))
  }

  continuous_rows <- which(component[["active"]] & finite_rows)
  if (length(continuous_rows) < 20L) {
    return(.iwmde_unsupported(parameter, "fewer than 20 continuous active samples"))
  }

  active_rows <- continuous_rows
  if (length(active_rows) > max_samples) {
    active_rows <- .iwmde_select_active_rows(
      rows        = active_rows,
      max_samples = max_samples
    )
  }

  active_values <- posterior_values[active_rows]
  if (stats::sd(active_values) <= sqrt(.Machine$double.eps)) {
    return(.iwmde_unsupported(parameter, "active samples have zero variance"))
  }

  support <- .iwmde_parameter_support(context, parameter, active_rows, parameter_spec)
  values  <- values[values > support[1] & values < support[2]]
  if (length(values) == 0L) {
    return(.iwmde_unsupported(parameter, "ordinate values are outside support"))
  }

  transform <- .iwmde_parameter_transform(support)
  xlim      <- .iwmde_plot_range(posterior_values[finite_rows], support)
  if (!all(is.finite(xlim)) || xlim[1] >= xlim[2]) {
    return(.iwmde_unsupported(parameter, "could not construct a finite plotting range"))
  }

  active_mass       <- mean(component[["active"]][finite_rows])
  continuous_values <- posterior_values[continuous_rows]
  row_states        <- .iwmde_row_states(context, active_rows, parameter, parameter_spec)
  baseline_log_q    <- vapply(row_states, function(state) {
    state[["baseline_log_q"]]
  }, numeric(1))
  keep_rows         <- is.finite(baseline_log_q)
  n_dropped_log_q   <- sum(!keep_rows)
  active_rows       <- active_rows[keep_rows]
  active_values     <- active_values[keep_rows]
  row_states        <- row_states[keep_rows]
  if (length(active_rows) < 20L) {
    return(.iwmde_unsupported(parameter, "fewer than 20 finite baseline log-q values"))
  }

  replacement <- .iwmde_replacement_spec(context, parameter, parameter_spec)
  if (identical(method, "q_grid_cmde")) {
    normalization_grid <- .iwmde_normalization_grid(
      values               = continuous_values,
      display_grid         = c(xlim, values),
      support              = support,
      transform            = transform,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob
    )
    if (is.null(normalization_grid)) {
      return(.iwmde_unsupported(parameter, "could not construct a normalization grid"))
    }
    density <- .iwmde_density_grid(
      context            = context,
      parameter          = parameter,
      display_grid       = values,
      normalization_grid = normalization_grid,
      transform          = transform,
      row_states         = row_states,
      active_mass        = active_mass,
      replacement        = replacement
    )
  } else {
    density <- .iwmde_density_iwmde(
      context        = context,
      parameter      = parameter,
      parameter_spec = parameter_spec,
      display_grid   = values,
      row_states     = row_states,
      active_rows    = active_rows,
      active_values  = active_values,
      weight_rows    = active_rows,
      weight_values  = active_values,
      support        = support,
      active_mass    = active_mass,
      replacement    = replacement
    )
  }
  if (density[["n_normalized_rows"]] < 20L) {
    return(.iwmde_unsupported(
      parameter,
      "fewer than 20 rows had finite conditional normalizers"
    ))
  }
  if (min(density[["finite_terms"]]) == 0L) {
    return(.iwmde_unsupported(
      parameter,
      "at least one IWMDE ordinate had no finite importance terms"
    ))
  }

  bf_diagnostics <- .iwmde_density_bf_diagnostics(density, values)
  out <- list(
    parameter    = parameter,
    status       = "ok",
    samples      = posterior_values[finite_rows],
    target_key   = target_key,
    active_rows  = active_rows,
    active_mass  = active_mass,
    point_masses = component[["point_masses"]],
    support      = support,
    xlim         = xlim,
    iwmde        = density,
    diagnostics  = list(
      integral               = NA_real_,
      plot_integral          = NA_real_,
      point_mass_total       = sum(component[["point_masses"]][["mass"]]),
      plot_total_mass        = NA_real_,
      display_mass_fraction  = NA_real_,
      normalization_integral = density[["normalization_integral"]],
      normalization_points   = density[["normalization_points"]],
      normalization_range    = density[["normalization_range"]],
      normalization_scale    = density[["normalization_scale"]],
      active_mass            = active_mass,
      n_active               = length(active_rows),
      n_total                = sum(finite_rows),
      n_dropped_log_q        = n_dropped_log_q,
      max_log_ratio          = density[["max_log_ratio"]],
      min_finite_terms       = min(density[["finite_terms"]]),
      n_normalized_rows      = density[["n_normalized_rows"]],
      min_ess                = min(density[["ess"]]),
      max_weight_share       = max(density[["max_weight_share"]]),
      max_mcse               = .iwmde_max_or_na(density[["mcse"]]),
      max_relative_mcse      = .iwmde_max_or_na(density[["relative_mcse"]]),
      estimator              = density[["estimator"]],
      weight_method          = density[["weight_method"]],
      display_grid           = "ordinate",
      bf_value               = bf_diagnostics[["bf_value"]],
      bf_included            = bf_diagnostics[["bf_included"]],
      bf_grid_index          = bf_diagnostics[["bf_grid_index"]],
      bf_ordinate            = bf_diagnostics[["bf_ordinate"]],
      bf_mcse                = bf_diagnostics[["bf_mcse"]],
      bf_relative_mcse       = bf_diagnostics[["bf_relative_mcse"]],
      bf_error_percent       = bf_diagnostics[["bf_error_percent"]],
      bf_finite_terms        = bf_diagnostics[["bf_finite_terms"]],
      bf_ess                 = bf_diagnostics[["bf_ess"]],
      bf_max_weight_share    = bf_diagnostics[["bf_max_weight_share"]],
      bf_max_log_ratio       = bf_diagnostics[["bf_max_log_ratio"]]
    )
  )
  class(out) <- c("iwmde_parameter_diagnostic", "list")

  .iwmde_cache_set(diagnostic_cache, cache_key, out)
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

  samples <- context[["posterior_samples"]]
  method <- .iwmde_normalize_method(method)
  display_grid_method <- .iwmde_normalize_display_grid(display_grid_method)
  include_values <- as.numeric(include_values)
  include_values <- include_values[is.finite(include_values)]
  parameter_spec <- .iwmde_parameter_spec(context, parameter, parameter_spec)
  target_key     <- .iwmde_target_key(parameter, parameter_spec)
  cache_key      <- .iwmde_diagnostic_cache_key(
    target_key            = target_key,
    n_points              = n_points,
    max_samples           = max_samples,
    normalization_points  = normalization_points,
    normalization_prob    = normalization_prob,
    method                = method,
    display_grid_method   = display_grid_method,
    include_values        = include_values
  )
  if (.iwmde_cache_has(diagnostic_cache, cache_key)) {
    return(.iwmde_relabel_diagnostic(
      diagnostic = .iwmde_cache_get(diagnostic_cache, cache_key),
      parameter  = parameter
    ))
  }

  if (identical(parameter_spec[["status"]], "unsupported")) {
    out <- .iwmde_unsupported(parameter, parameter_spec[["reason"]])
    .iwmde_cache_set(diagnostic_cache, cache_key, out)
    return(out)
  }

  values <- .iwmde_parameter_values(context, parameter, parameter_spec)
  finite <- is.finite(values)
  rows   <- .iwmde_parameter_condition_rows(context, parameter_spec)
  finite_rows <- finite & rows
  if (!any(finite_rows)) {
    return(.iwmde_unsupported(parameter, "posterior samples are not finite"))
  }

  component <- .iwmde_parameter_components(context, parameter, parameter_spec)
  component <- .iwmde_restrict_parameter_component(component, finite_rows)
  if (!any(component[["active"]])) {
    return(.iwmde_point_only_diagnostic(
      parameter = parameter,
      samples   = values[finite_rows],
      component = component
    ))
  }

  continuous_rows <- which(component[["active"]] & finite_rows)
  if (length(continuous_rows) < 20L) {
    return(.iwmde_unsupported(parameter, "fewer than 20 continuous active samples"))
  }

  active_rows <- continuous_rows
  if (length(active_rows) > max_samples) {
    active_rows <- .iwmde_select_active_rows(
      rows        = active_rows,
      max_samples = max_samples
    )
  }

  active_values <- values[active_rows]
  if (stats::sd(active_values) <= sqrt(.Machine$double.eps)) {
    return(.iwmde_unsupported(parameter, "active samples have zero variance"))
  }

  support   <- .iwmde_parameter_support(context, parameter, active_rows, parameter_spec)
  transform <- .iwmde_parameter_transform(support)
  xlim      <- .iwmde_plot_range(values[finite_rows], support)
  xlim      <- .iwmde_include_plot_values(xlim, include_values, support)
  if (!all(is.finite(xlim)) || xlim[1] >= xlim[2]) {
    return(.iwmde_unsupported(parameter, "could not construct a finite plotting range"))
  }

  active_mass       <- mean(component[["active"]][finite_rows])
  continuous_values <- values[continuous_rows]
  display_grid      <- .iwmde_display_grid(
    xlim        = xlim,
    n_points    = n_points,
    transform   = transform,
    values      = c(continuous_values, include_values),
    grid_method = display_grid_method
  )
  display_grid <- .iwmde_include_display_values(
    grid    = display_grid,
    values  = include_values,
    xlim    = xlim,
    support = support
  )

  row_states        <- .iwmde_row_states(context, active_rows, parameter, parameter_spec)
  baseline_log_q    <- vapply(row_states, function(state) {
    state[["baseline_log_q"]]
  }, numeric(1))
  keep_rows         <- is.finite(baseline_log_q)
  n_dropped_log_q   <- sum(!keep_rows)
  active_rows      <- active_rows[keep_rows]
  active_values    <- active_values[keep_rows]
  row_states       <- row_states[keep_rows]
  if (length(active_rows) < 20L) {
    return(.iwmde_unsupported(parameter, "fewer than 20 finite baseline log-q values"))
  }
  replacement      <- .iwmde_replacement_spec(context, parameter, parameter_spec)
  if (identical(method, "q_grid_cmde")) {
    normalization_grid <- .iwmde_normalization_grid(
      values               = continuous_values,
      display_grid         = display_grid,
      support              = support,
      transform            = transform,
      normalization_points = normalization_points,
      normalization_prob   = normalization_prob
    )
    if (is.null(normalization_grid)) {
      return(.iwmde_unsupported(parameter, "could not construct a normalization grid"))
    }
    density <- .iwmde_density_grid(
      context            = context,
      parameter          = parameter,
      display_grid       = display_grid,
      normalization_grid = normalization_grid,
      transform          = transform,
      row_states         = row_states,
      active_mass        = active_mass,
      replacement        = replacement
    )
  } else {
    density <- .iwmde_density_iwmde(
      context        = context,
      parameter      = parameter,
      parameter_spec = parameter_spec,
      display_grid   = display_grid,
      row_states     = row_states,
      active_rows    = active_rows,
      active_values  = active_values,
      weight_rows    = active_rows,
      weight_values  = active_values,
      support        = support,
      active_mass    = active_mass,
      replacement    = replacement
    )
  }
  if (density[["n_normalized_rows"]] < 20L) {
    return(.iwmde_unsupported(
      parameter,
      "fewer than 20 rows had finite conditional normalizers"
    ))
  }
  if (min(density[["finite_terms"]]) == 0L) {
    return(.iwmde_unsupported(
      parameter,
      "at least one IWMDE grid point had no finite importance terms"
    ))
  }

  kde              <- .iwmde_kde(continuous_values, xlim, n_points, mass = active_mass)
  hist_data        <- .iwmde_histogram(continuous_values, xlim, mass = active_mass)
  plot_integral    <- .iwmde_trapz(density[["x"]], density[["y"]])
  point_mass_total <- sum(component[["point_masses"]][["mass"]])
  max_mcse          <- .iwmde_max_or_na(density[["mcse"]])
  max_relative_mcse <- .iwmde_max_or_na(density[["relative_mcse"]])
  bf_diagnostics    <- .iwmde_density_bf_diagnostics(density, include_values)

  out <- list(
    parameter    = parameter,
    status       = "ok",
    samples      = values[finite_rows],
    target_key   = target_key,
    active_rows  = active_rows,
    active_mass  = active_mass,
    point_masses = component[["point_masses"]],
    support      = support,
    xlim         = xlim,
    histogram    = hist_data,
    kde          = kde,
    iwmde        = density,
    diagnostics  = list(
      integral                    = plot_integral,
      plot_integral               = plot_integral,
      point_mass_total            = point_mass_total,
      plot_total_mass             = plot_integral + point_mass_total,
      display_mass_fraction       = plot_integral / active_mass,
      normalization_integral      = density[["normalization_integral"]],
      normalization_points        = density[["normalization_points"]],
      normalization_range         = density[["normalization_range"]],
      normalization_scale         = transform[["type"]],
      n_active                    = length(active_rows),
      n_total                     = sum(finite_rows),
      n_dropped_log_q             = n_dropped_log_q,
      max_log_ratio               = density[["max_log_ratio"]],
      min_finite_terms            = min(density[["finite_terms"]]),
      n_normalized_rows           = density[["n_normalized_rows"]],
      min_ess                     = min(density[["ess"]]),
      max_weight_share            = max(density[["max_weight_share"]]),
      max_mcse                    = max_mcse,
      max_relative_mcse           = max_relative_mcse,
      plot_integral_mcse          = density[["integral_mcse"]],
      plot_integral_relative_mcse = density[["integral_relative_mcse"]],
      batch_size                  = density[["batch_size"]],
      n_batches                   = density[["n_batches"]],
      estimator                   = density[["estimator"]],
      weight_method               = density[["weight_method"]],
      display_grid                = display_grid_method,
      bf_value                    = bf_diagnostics[["bf_value"]],
      bf_included                 = bf_diagnostics[["bf_included"]],
      bf_grid_index               = bf_diagnostics[["bf_grid_index"]],
      bf_ordinate                 = bf_diagnostics[["bf_ordinate"]],
      bf_mcse                     = bf_diagnostics[["bf_mcse"]],
      bf_relative_mcse            = bf_diagnostics[["bf_relative_mcse"]],
      bf_error_percent            = bf_diagnostics[["bf_error_percent"]],
      bf_finite_terms             = bf_diagnostics[["bf_finite_terms"]],
      bf_ess                      = bf_diagnostics[["bf_ess"]],
      bf_max_weight_share         = bf_diagnostics[["bf_max_weight_share"]],
      bf_max_log_ratio            = bf_diagnostics[["bf_max_log_ratio"]]
    )
  )
  class(out) <- c("iwmde_parameter_diagnostic", "list")

  .iwmde_cache_set(diagnostic_cache, cache_key, out)
  return(out)
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


.iwmde_diagnostic_cache_key <- function(target_key, n_points, max_samples,
                                        normalization_points,
                                        normalization_prob, method,
                                        display_grid_method,
                                        include_values = NULL) {

  include_values <- sort(unique(as.numeric(include_values)))
  include_values <- include_values[is.finite(include_values)]
  include_key    <- if (length(include_values) > 0L) {
    paste(.iwmde_key_number(include_values), collapse = ",")
  } else {
    ""
  }

  return(paste(
    target_key,
    paste0("n_points=", n_points),
    paste0("max_samples=", max_samples),
    paste0("normalization_points=", normalization_points),
    paste0("normalization_prob=", .iwmde_key_number(normalization_prob)),
    paste0("method=", method),
    paste0("display_grid=", display_grid_method),
    paste0("include_values=", include_key),
    sep = "\r"
  ))
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


.iwmde_select_active_rows <- function(rows, max_samples) {

  if (length(rows) <= max_samples) {
    return(rows)
  }

  selected <- .thin_sample_rows(length(rows), max_samples)

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
      label <- paste0(selected_parameter[["label"]], ": ", level)
      key   <- label
      specs[[key]] <- list(
        type            = "linear",
        label           = label,
        parameter       = parameter_name,
        level           = level,
        weights         = attr(level_samples[[level]], "linear_weights"),
        conditional     = attr(level_samples[[level]], "effective_conditional"),
        conditional_rule = marginal_means_object[["conditional_rule"]],
        source          = "marginal_means",
        selected        = selected_parameter,
        marginal_type   = type,
        marginal_samples = level_samples[[level]]
      )
    }
  }

  return(specs)
}


