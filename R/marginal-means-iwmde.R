# ============================================================================ #
# Marginal-Means qCMDE/IWMDE Helpers
# ============================================================================ #

.marginal_means_precompute_type <- function(type, model_averaged) {

  if (is.null(type)) {
    return(c("averaged", "conditional"))
  }

  BayesTools::check_char(type, "type", check_length = 0, allow_NA = FALSE)
  type <- match.arg(type, c("averaged", "conditional"), several.ok = TRUE)
  if (!isTRUE(model_averaged) && "conditional" %in% type) {
    stop("The 'type' argument is available only for RoBMA marginal means.",
         call. = FALSE)
  }

  return(unique(type))
}


.marginal_means_precompute_levels <- function(levels) {

  if (is.null(levels)) {
    return(NULL)
  }

  BayesTools::check_char(levels, "levels", check_length = 0, allow_NA = FALSE)

  return(unique(levels))
}


# Extract the fitted source object stored by marginal_means().
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


.marginal_means_attach_iwmde <- function(object, marginal_means_object,
                                         n_points, max_samples,
                                         normalization_points,
                                         normalization_prob, density_method,
                                         display_grid, null_hypothesis,
                                         parameter, type, levels, targeted,
                                         include_ordinates = TRUE,
                                         ordinate_control = NULL) {

  .check_iwmde_available(object, "qCMDE/IWMDE marginal_means()")
  if (isTRUE(include_ordinates)) {
    .iwmde_check_point_ordinate_supported(object, density_method)
  }

  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }

  context              <- .iwmde_context(object)
  estimate_cache       <- .iwmde_estimate_cache()
  diagnostics          <- list()
  ordinate_diagnostics <- list()
  density_control_list <- list(
    n_points             = n_points,
    max_samples          = max_samples,
    normalization_points = normalization_points,
    normalization_prob   = normalization_prob,
    display_grid         = display_grid
  )
  if (is.null(ordinate_control)) {
    ordinate_control_list <- density_control_list
  } else {
    ordinate_control_list <- ordinate_control
    if (is.null(ordinate_control_list[["normalization_points"]])) {
      ordinate_control_list[["normalization_points"]] <- max(
        50L,
        ordinate_control_list[["n_points"]]
      )
    }
  }
  ordinate_control_list[["display_grid"]] <- "ordinate"
  include_values       <- NULL
  specs_by_type        <- .marginal_means_iwmde_specs_by_type(
    marginal_means_object = marginal_means_object,
    parameter             = parameter,
    type                  = type,
    levels                = levels,
    targeted              = targeted
  )
  ordinate_specs <- if (isTRUE(include_ordinates)) {
    .marginal_means_iwmde_ordinate_specs(
      marginal_means_object = marginal_means_object,
      parameter             = parameter,
      levels                = levels,
      targeted              = targeted
    )
  } else {
    list()
  }

  for (current_type in names(specs_by_type)) {
    specs <- specs_by_type[[current_type]]

    estimates <- lapply(specs, function(spec) {
      .iwmde_estimate(
        context         = context,
        parameter       = spec[["label"]],
        density_method  = density_method,
        density_control = density_control_list,
        outputs         = "density",
        values          = include_values,
        parameter_spec  = spec,
        metadata        = .marginal_means_iwmde_metadata(
          marginal_means_object = marginal_means_object,
          type                  = current_type,
          spec                  = spec
        ),
        cache           = estimate_cache
      )
    })
    names(estimates) <- names(specs)
    diagnostics[[current_type]] <- lapply(estimates, function(estimate) {
      estimate[["diagnostics"]][["density"]]
    })

    marginal_means_object <- .marginal_means_attach_iwmde_type(
      marginal_means_object = marginal_means_object,
      type                  = current_type,
      specs                 = specs,
      estimates             = estimates
    )
  }

  if (length(ordinate_specs) > 0L) {
    ordinate_estimates <- lapply(ordinate_specs, function(spec) {
      .iwmde_estimate(
        context         = context,
        parameter       = spec[["label"]],
        density_method  = density_method,
        density_control = ordinate_control_list,
        outputs         = "ordinate",
        values          = null_hypothesis,
        parameter_spec  = spec,
        metadata        = .marginal_means_iwmde_metadata(
          marginal_means_object = marginal_means_object,
          type                  = "conditional",
          spec                  = spec
        ),
        cache           = estimate_cache
      )
    })
    names(ordinate_estimates) <- names(ordinate_specs)
    ordinate_diagnostics[["conditional"]] <- lapply(ordinate_estimates, function(estimate) {
      estimate[["diagnostics"]][["ordinate"]]
    })
    marginal_means_object <- .marginal_means_attach_iwmde_ordinate_type(
      marginal_means_object = marginal_means_object,
      type                  = "conditional",
      specs                 = ordinate_specs,
      estimates             = ordinate_estimates
    )
  }

  marginal_means_object[["density_diagnostics"]] <- diagnostics
  marginal_means_object[["ordinate_diagnostics"]] <- ordinate_diagnostics
  marginal_means_object[["density_settings"]] <- list(
    n_points             = n_points,
    max_samples          = max_samples,
    normalization_points = normalization_points,
    normalization_prob   = normalization_prob,
    density_method       = density_method,
    method               = .density_method_iwmde_estimator(density_method),
    display_grid         = display_grid,
    include_values       = include_values,
    ordinate_values      = null_hypothesis,
    parameter            = parameter,
    type                 = type,
    levels               = levels,
    include_ordinates    = isTRUE(include_ordinates)
  )
  marginal_means_object[["ordinate_settings"]] <- ordinate_control_list
  if (isTRUE(include_ordinates)) {
    marginal_means_object[["inference"]] <- .marginal_means_refresh_iwmde_bf(
      inference            = marginal_means_object[["inference"]],
      parameters           = marginal_means_object[["parameters"]],
      null_hypothesis      = marginal_means_object[["null_hypothesis"]],
      density_method       = marginal_means_object[["density_method"]],
      object               = marginal_means_object
    )
  }

  return(marginal_means_object)
}


.marginal_means_iwmde_specs_by_type <- function(marginal_means_object,
                                                parameter, type, levels,
                                                targeted) {

  specs_by_type <- lapply(type, function(current_type) {
    .iwmde_marginal_means_specs(
      marginal_means_object = marginal_means_object,
      parameter             = parameter,
      type                  = current_type,
      levels                = levels
    )
  })
  names(specs_by_type) <- type

  empty <- vapply(specs_by_type, function(specs) {
    length(specs) == 0L
  }, logical(1))
  if (any(empty) && isTRUE(targeted)) {
    stop(
      "No qCMDE/IWMDE marginal-means targets matched 'parameter', ",
      "'type', and 'levels'.",
      call. = FALSE
    )
  }

  return(specs_by_type[!empty])
}


.marginal_means_iwmde_ordinate_specs <- function(marginal_means_object,
                                                 parameter, levels,
                                                 targeted) {

  specs_by_type <- .marginal_means_iwmde_specs_by_type(
    marginal_means_object = marginal_means_object,
    parameter             = parameter,
    type                  = "conditional",
    levels                = levels,
    targeted              = targeted
  )

  return(specs_by_type[["conditional"]])
}


.marginal_means_attach_iwmde_type <- function(marginal_means_object, type,
                                              specs, estimates) {

  for (name in names(specs)) {
    estimate <- estimates[[name]]
    diagnostic <- estimate[["diagnostics"]][["density"]]
    if (!identical(diagnostic[["status"]], "ok") ||
        is.null(estimate[["posterior_density"]])) {
      next
    }

    parameter <- specs[[name]][["parameter"]]
    level     <- specs[[name]][["level"]]
    samples   <- marginal_means_object[["inference"]][[type]][[parameter]]
    if (is.null(samples)) {
      next
    }

    if (is.list(samples) && !is.null(samples[[level]])) {
      attr(samples[[level]], "posterior_density") <- estimate[["posterior_density"]]
    } else {
      attr(samples, "posterior_density") <- estimate[["posterior_density"]]
    }
    marginal_means_object[["inference"]][[type]][[parameter]] <- samples
  }

  return(marginal_means_object)
}


.marginal_means_attach_iwmde_ordinate_type <- function(marginal_means_object, type,
                                                       specs, estimates) {

  for (name in names(specs)) {
    estimate <- estimates[[name]]
    diagnostic <- estimate[["diagnostics"]][["ordinate"]]
    if (!identical(diagnostic[["status"]], "ok")) {
      next
    }

    parameter <- specs[[name]][["parameter"]]
    level     <- specs[[name]][["level"]]
    samples   <- marginal_means_object[["inference"]][[type]][[parameter]]
    if (is.null(samples)) {
      next
    }

    if (is.list(samples) && !is.null(samples[[level]])) {
      attr(samples[[level]], "posterior_ordinate") <- estimate[["posterior_ordinate"]]
    } else {
      attr(samples, "posterior_ordinate") <- estimate[["posterior_ordinate"]]
    }
    marginal_means_object[["inference"]][[type]][[parameter]] <- samples
  }

  return(marginal_means_object)
}


.marginal_means_iwmde_metadata <- function(marginal_means_object, type, spec) {

  parameter <- spec[["parameter"]]
  level     <- spec[["level"]]
  samples   <- marginal_means_object[["inference"]][[type]][[parameter]]
  if (is.null(samples)) {
    return(list(parameter = parameter, level = level))
  }

  if (is.list(samples) && !is.null(samples[[level]])) {
    return(.iwmde_posterior_metadata(
      samples   = samples[[level]],
      parameter = parameter,
      level     = level
    ))
  }

  return(.iwmde_posterior_metadata(
    samples   = samples,
    parameter = parameter,
    level     = level
  ))
}


.marginal_means_iwmde_settings_control <- function(object, density_method,
                                                  display_grid) {

  settings <- if (identical(display_grid, "ordinate") &&
                  is.list(object[["ordinate_settings"]])) {
    object[["ordinate_settings"]]
  } else {
    object[["density_settings"]]
  }
  control <- settings[intersect(
    c(
      "n_points",
      "max_samples",
      "initial_samples",
      "target_relative_mcse",
      "normalization_points",
      "normalization_prob",
      "display_grid"
    ),
    names(settings)
  )]
  control[["display_grid"]] <- display_grid
  control <- .iwmde_density_control_resolve(
    density_method  = density_method,
    density_control = control,
    purpose         = if (identical(display_grid, "ordinate")) {
      "ordinate"
    } else {
      "density"
    }
  )

  .marginal_means_density_control_effective(control)
}


.marginal_means_ordinate_request_provenance <- function(object, parameter,
                                                        null_hypothesis,
                                                        density_method) {

  if (is.null(object)) {
    stop(
      "Cannot validate a stored qCMDE/IWMDE ordinate without its ",
      "marginal-means object.",
      call. = FALSE
    )
  }

  source_object <- .iwmde_marginal_means_source_object(object)
  context       <- .iwmde_context(source_object)
  specs <- .iwmde_marginal_means_specs(
    marginal_means_object = object,
    parameter             = parameter,
    type                  = "conditional",
    levels                = NULL
  )
  density_control <- .marginal_means_iwmde_settings_control(
    object         = object,
    density_method = density_method,
    display_grid   = "ordinate"
  )
  out <- list()

  for (name in names(specs)) {
    spec <- specs[[name]]
    metadata <- .marginal_means_iwmde_metadata(
      marginal_means_object = object,
      type                  = "conditional",
      spec                  = spec
    )
    out[[spec[["level"]]]] <- .iwmde_request_provenance(
      context         = context,
      parameter       = spec[["label"]],
      density_method  = density_method,
      density_control = density_control,
      attribute       = "ordinate",
      value           = null_hypothesis,
      parameter_spec  = spec,
      metadata        = metadata
    )
  }

  return(out)
}


.marginal_means_refresh_iwmde_bf <- function(inference, parameters,
                                             null_hypothesis,
                                             density_method = "KDE",
                                             object = NULL) {

  density_method <- .density_method_normalize(
    density_method = density_method
  )
  if (!.density_method_uses_precomputed(density_method) ||
      is.null(inference[["conditional"]]) ||
      is.null(inference[["inference"]])) {
    return(inference)
  }

  for (parameter in parameters) {
    posterior <- inference[["conditional"]][[parameter]]
    if (is.null(posterior)) {
      next
    }
    provenance <- .marginal_means_ordinate_request_provenance(
      object          = object,
      parameter       = parameter,
      null_hypothesis = null_hypothesis,
      density_method  = density_method
    )
    bf <- .marginal_means_iwmde_bf(
      posterior       = posterior,
      null_hypothesis = null_hypothesis,
      provenance      = provenance,
      density_method  = density_method
    )

    inference[["inference"]][[parameter]] <- bf
  }

  return(inference)
}


.marginal_means_iwmde_bf <- function(posterior, null_hypothesis,
                                     provenance = NULL,
                                     density_method = NULL,
                                     warning = .marginal_means_iwmde_bf_warning()) {

  if (!is.list(posterior)) {
    return(.marginal_means_iwmde_bf_scalar(
      posterior       = posterior,
      null_hypothesis = null_hypothesis,
      provenance      = provenance,
      density_method  = density_method,
      warning         = warning
    ))
  }

  out <- vector("list", length(posterior))
  names(out) <- names(posterior)
  valid <- vapply(names(posterior), function(level) {
    .iwmde_posterior_ordinate_matches_request(
      posterior_ordinate = attr(posterior[[level]], "posterior_ordinate",
                                exact = TRUE),
      value              = null_hypothesis,
      provenance         = provenance[[level]]
    )
  }, logical(1))
  if (!any(valid)) {
    for (i in seq_along(out)) {
      out[[i]] <- .marginal_means_unavailable_bf_scalar(
        .marginal_means_iwmde_bf_warning(
          attr(posterior[[i]], "posterior_ordinate", exact = TRUE),
          default = warning
        )
      )
    }
    return(out)
  }

  bf_posterior <- .marginal_means_bf_posterior(posterior[valid])
  class(bf_posterior) <- class(posterior)
  bf <- BayesTools::Savage_Dickey_BF(
    posterior            = bf_posterior,
    null_hypothesis      = null_hypothesis,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  for (i in seq_along(out)) {
    if (isTRUE(valid[[i]])) {
      out[[i]] <- .iwmde_bf_append_warning(
        bf                 = bf[[names(out)[[i]]]],
        posterior_ordinate = attr(posterior[[i]], "posterior_ordinate",
                                  exact = TRUE)
      )
    } else {
      out[[i]] <- .marginal_means_unavailable_bf_scalar(
        .marginal_means_iwmde_bf_warning(
          attr(posterior[[i]], "posterior_ordinate", exact = TRUE),
          default = warning
        )
      )
    }
  }

  return(out)
}


.marginal_means_iwmde_bf_scalar <- function(posterior, null_hypothesis,
                                            provenance, density_method,
                                            warning) {

  if (!.iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = attr(posterior, "posterior_ordinate", exact = TRUE),
    value              = null_hypothesis,
    provenance         = provenance
  )) {
    return(.marginal_means_unavailable_bf_scalar(
      .marginal_means_iwmde_bf_warning(
        attr(posterior, "posterior_ordinate", exact = TRUE),
        default = warning
      )
    ))
  }

  posterior <- .marginal_means_bf_posterior(posterior)
  bf <- BayesTools::Savage_Dickey_BF(
    posterior            = posterior,
    null_hypothesis      = null_hypothesis,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  .iwmde_bf_append_warning(
    bf                 = bf,
    posterior_ordinate = attr(posterior, "posterior_ordinate", exact = TRUE)
  )
}


.marginal_means_iwmde_bf_warning <- function(posterior_ordinate = NULL,
                                             default = NULL) {

  if (is.null(default)) {
    default <- paste0(
      "Precomputed posterior ordinate was unavailable or failed diagnostics; ",
      "Savage-Dickey Bayes factor is not reported."
    )
  }

  reasons <- .iwmde_posterior_ordinate_failure_reasons(posterior_ordinate)
  if (length(reasons) == 0L) {
    return(default)
  }

  paste0(
    "Precomputed posterior ordinate failed diagnostics (",
    paste(reasons, collapse = "; "),
    "); ",
    "Savage-Dickey Bayes factor is not reported."
  )
}


.marginal_means_unavailable_bf_scalar <- function(warning) {

  out <- NA_real_
  attr(out, "warnings") <- warning
  return(out)
}


.marginal_means_posterior_density_matches <- function(samples, provenance) {

  density <- attr(samples, "posterior_density", exact = TRUE)
  .iwmde_posterior_density_matches_request(
    posterior_density = density,
    provenance        = provenance
  )
}


.marginal_means_missing_posterior_density_levels <- function(samples,
                                                             provenance,
                                                             parameter = NULL) {

  if (!is.list(samples)) {
    level <- .marginal_means_posterior_density_levels(samples, parameter)
    if (.marginal_means_posterior_density_matches(
          samples    = samples,
          provenance = provenance[[level]]
        )) {
      return(character())
    }

    return(level)
  }

  levels <- names(samples)
  missing <- !vapply(levels, function(level) {
    .marginal_means_posterior_density_matches(
      samples    = samples[[level]],
      provenance = provenance[[level]]
    )
  }, logical(1))

  return(levels[missing])
}


.marginal_means_posterior_density_levels <- function(samples, parameter = NULL) {

  if (is.list(samples)) {
    return(names(samples))
  }

  level <- names(samples)
  if (is.null(level) || length(level) == 0L || !nzchar(level[[1L]])) {
    if (is.null(parameter)) {
      return("")
    }
    return(parameter)
  }

  return(level[[1L]])
}


.marginal_means_density_control_effective <- function(density_control) {

  if (is.null(density_control[["normalization_points"]])) {
    density_control[["normalization_points"]] <- max(
      50L,
      density_control[["n_points"]]
    )
  }

  return(density_control)
}


.marginal_means_density_request_provenance <- function(x, selected, type,
                                                       density_method,
                                                       density_control,
                                                       context) {

  specs_by_type <- .marginal_means_iwmde_specs_by_type(
    marginal_means_object = x,
    parameter             = selected[["term"]],
    type                  = type,
    levels                = NULL,
    targeted              = TRUE
  )
  specs <- specs_by_type[[type]]
  out   <- list()

  for (name in names(specs)) {
    spec <- specs[[name]]
    metadata <- .marginal_means_iwmde_metadata(
      marginal_means_object = x,
      type                  = type,
      spec                  = spec
    )
    out[[spec[["level"]]]] <- .iwmde_request_provenance(
      context         = context,
      parameter       = spec[["label"]],
      density_method  = density_method,
      density_control = density_control,
      attribute       = "density",
      parameter_spec  = spec,
      metadata        = metadata
    )
  }

  return(out)
}


.marginal_means_plot_iwmde_message <- function(density_method, n_levels,
                                               force) {

  level_label <- if (n_levels == 1L) "level" else "levels"
  reason <- if (isTRUE(force)) {
    "because 'density_control' was supplied"
  } else {
    "because stored densities are missing"
  }

  message(
    "Computing ", density_method, " density for ", n_levels,
    " marginal-means ", level_label, " ", reason, ". ",
    "Precompute with marginal_means(..., density_method = \"",
    density_method, "\") to reuse densities across plots."
  )

  return(invisible(NULL))
}


.marginal_means_plot_prepare_iwmde <- function(x, selected, type,
                                               density_method,
                                               density_control,
                                               force = FALSE) {

  parameter <- selected[["parameter"]]
  samples   <- x[["inference"]][[type]][[parameter]]
  density_control <- .marginal_means_density_control_effective(density_control)
  source_object   <- .iwmde_marginal_means_source_object(x)
  context         <- .iwmde_context(source_object)
  provenance <- .marginal_means_density_request_provenance(
    x               = x,
    selected        = selected,
    type            = type,
    density_method  = density_method,
    density_control = density_control,
    context         = context
  )
  missing_levels <- if (isTRUE(force)) {
    .marginal_means_posterior_density_levels(samples, parameter)
  } else {
    .marginal_means_missing_posterior_density_levels(
      samples    = samples,
      provenance = provenance,
      parameter  = parameter
    )
  }
  if (length(missing_levels) > 0L) {
    .marginal_means_plot_iwmde_message(
      density_method = density_method,
      n_levels       = length(missing_levels),
      force          = force
    )
    x <- .marginal_means_attach_iwmde(
      object                = source_object,
      marginal_means_object = x,
      n_points              = density_control[["n_points"]],
      max_samples           = density_control[["max_samples"]],
      normalization_points  = density_control[["normalization_points"]],
      normalization_prob    = density_control[["normalization_prob"]],
      density_method        = density_method,
      display_grid          = density_control[["display_grid"]],
      null_hypothesis       = x[["null_hypothesis"]],
      parameter             = selected[["term"]],
      type                  = type,
      levels                = missing_levels,
      targeted              = TRUE,
      include_ordinates     = FALSE
    )
  }

  missing_levels <- .marginal_means_missing_posterior_density_levels(
    samples    = x[["inference"]][[type]][[parameter]],
    provenance = provenance,
    parameter  = parameter
  )
  if (length(missing_levels) > 0L) {
    stop(
      density_method, " density was unavailable for marginal-means level(s): ",
      paste0("'", missing_levels, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }

  return(x)
}


.marginal_means_bf_posterior <- function(samples) {

  if (!is.list(samples)) {
    if (!.iwmde_posterior_ordinate_supports_bf(attr(samples, "posterior_ordinate"))) {
      attr(samples, "posterior_ordinate") <- NULL
      attr(samples, "posterior_density") <- NULL
    }
    return(samples)
  }

  for (i in seq_along(samples)) {
    if (!.iwmde_posterior_ordinate_supports_bf(attr(samples[[i]], "posterior_ordinate"))) {
      attr(samples[[i]], "posterior_ordinate") <- NULL
      attr(samples[[i]], "posterior_density") <- NULL
    }
  }

  return(samples)
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
          weights   = attr(level_sample, "linear_weights"),
          prior_density = attr(level_sample, "prior_density", exact = TRUE)
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
