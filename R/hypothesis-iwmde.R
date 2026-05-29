.hypothesis_brma_attach_iwmde <- function(object, posterior, parameter,
                                          hypothesis, conditional,
                                          n_points, max_samples,
                                          normalization_points,
                                          normalization_prob, density_method,
                                          standardized_coefficients,
                                          n_samples) {

  point_refs <- .hypothesis_brma_point_refs(hypothesis, parameter)
  if (nrow(point_refs) == 0L) {
    return(posterior)
  }

  sample_parameter <- .as_mixed_posteriors_parameters(object, parameter)
  raw_samples <- BayesTools::as_mixed_posteriors(
    model            = object[["fit"]],
    parameters       = sample_parameter,
    conditional      = conditional,
    conditional_rule = "AND",
    transform_scaled = FALSE,
    n_prior_samples  = n_samples
  )
  raw_sample_parameter <- .plot_brma_density_sample_parameter(
    samples          = raw_samples,
    parameter        = parameter,
    sample_parameter = sample_parameter
  )
  raw_posterior <- BayesTools::marginal_posterior(
    samples       = raw_samples,
    parameter     = raw_sample_parameter,
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = n_samples
  )

  if (is.null(normalization_points)) {
    normalization_points <- max(50L, n_points)
  }

  context          <- .iwmde_context(object)
  diagnostic_cache <- .iwmde_diagnostic_cache()
  method           <- .density_method_iwmde_estimator(density_method)

  for (i in seq_len(nrow(point_refs))) {
    ref <- point_refs[i, , drop = FALSE]
    if (is.na(ref[["level"]])) {
      posterior <- .hypothesis_brma_attach_iwmde_scalar(
        posterior                = posterior,
        raw_posterior            = raw_posterior,
        context                  = context,
        diagnostic_cache         = diagnostic_cache,
        parameter                = parameter,
        value                    = ref[["value"]],
        conditional              = conditional,
        max_samples              = max_samples,
        normalization_points     = normalization_points,
        normalization_prob       = normalization_prob,
        method                   = method,
        density_method           = density_method,
        standardized_coefficients = standardized_coefficients
      )
    } else {
      posterior <- .hypothesis_brma_attach_iwmde_level(
        posterior            = posterior,
        raw_posterior        = raw_posterior,
        context              = context,
        diagnostic_cache     = diagnostic_cache,
        parameter            = parameter,
        level                = ref[["level"]],
        value                = ref[["value"]],
        conditional          = conditional,
        max_samples          = max_samples,
        normalization_points = normalization_points,
        normalization_prob   = normalization_prob,
        method               = method,
        density_method       = density_method
      )
    }
  }

  return(posterior)
}


.hypothesis_brma_attach_iwmde_scalar <- function(
    posterior, raw_posterior, context, diagnostic_cache, parameter, value,
    conditional, max_samples, normalization_points, normalization_prob, method,
    density_method, standardized_coefficients) {

  raw_values <- as.numeric(raw_posterior)
  display_values <- as.numeric(posterior)
  transform <- .plot_brma_affine_sample_transform(raw_values, display_values)
  if (is.null(transform)) {
    if (isTRUE(standardized_coefficients)) {
      transform <- list(intercept = 0, slope = 1)
    } else {
      stop("Could not align qCMDE/IWMDE ordinate to the displayed coefficient scale.",
           call. = FALSE)
    }
  }

  raw_value <- (value - transform[["intercept"]]) / transform[["slope"]]
  diagnostic <- .iwmde_parameter_ordinate_diagnostic(
    context              = context,
    parameter            = parameter,
    values               = raw_value,
    max_samples          = max_samples,
    normalization_points = normalization_points,
    normalization_prob   = normalization_prob,
    method               = method,
    parameter_spec       = list(
      type             = "primitive",
      conditional      = conditional,
      conditional_rule = "AND"
    ),
    diagnostic_cache     = diagnostic_cache
  )

  ordinate <- .iwmde_posterior_ordinate_attribute(
    diagnostic     = diagnostic,
    density_method = density_method
  )
  if (is.null(ordinate)) {
    stop(
      "Precomputed ", density_method, " posterior ordinate is unavailable for '",
      parameter, " = ", value, "': ",
      .hypothesis_brma_diagnostic_reason(diagnostic),
      call. = FALSE
    )
  }

  ordinate <- .hypothesis_brma_transform_ordinate(
    ordinate  = ordinate,
    transform = transform,
    value     = value
  )
  attr(posterior, "posterior_ordinate") <- .hypothesis_brma_append_ordinate(
    existing = attr(posterior, "posterior_ordinate", exact = TRUE),
    ordinate = ordinate
  )

  return(posterior)
}


.hypothesis_brma_attach_iwmde_level <- function(
    posterior, raw_posterior, context, diagnostic_cache, parameter, level,
    value, conditional, max_samples, normalization_points, normalization_prob,
    method, density_method) {

  if (!is.list(posterior) || !level %in% names(posterior)) {
    stop("Hypothesis references unknown level '", level,
         "' for parameter '", parameter, "'.", call. = FALSE)
  }
  if (!is.list(raw_posterior) || !level %in% names(raw_posterior)) {
    stop("Raw posterior for level '", level, "' is unavailable.",
         call. = FALSE)
  }

  raw_values     <- as.numeric(raw_posterior[[level]])
  display_values <- as.numeric(posterior[[level]])
  transform <- .plot_brma_affine_sample_transform(raw_values, display_values)
  if (is.null(transform)) {
    transform <- list(intercept = 0, slope = 1)
  }
  raw_value <- (value - transform[["intercept"]]) / transform[["slope"]]

  weights <- attr(raw_posterior[[level]], "linear_weights", exact = TRUE)
  diagnostic <- .iwmde_parameter_ordinate_diagnostic(
    context              = context,
    parameter            = paste0(parameter, "[", level, "]"),
    values               = raw_value,
    max_samples          = max_samples,
    normalization_points = normalization_points,
    normalization_prob   = normalization_prob,
    method               = method,
    parameter_spec       = list(
      type             = "linear",
      weights          = weights,
      conditional      = if (is.null(conditional)) {
        attr(raw_posterior[[level]], "effective_conditional", exact = TRUE)
      } else {
        conditional
      },
      conditional_rule = "AND"
    ),
    diagnostic_cache     = diagnostic_cache
  )

  ordinate <- .iwmde_posterior_ordinate_attribute(
    diagnostic     = diagnostic,
    density_method = density_method
  )
  if (is.null(ordinate)) {
    stop(
      "Precomputed ", density_method, " posterior ordinate is unavailable for '",
      parameter, "[", level, "] = ", value, "': ",
      .hypothesis_brma_diagnostic_reason(diagnostic),
      call. = FALSE
    )
  }

  ordinate <- .hypothesis_brma_transform_ordinate(
    ordinate  = ordinate,
    transform = transform,
    value     = value
  )
  attr(posterior[[level]], "posterior_ordinate") <- .hypothesis_brma_append_ordinate(
    existing = attr(posterior[[level]], "posterior_ordinate", exact = TRUE),
    ordinate = ordinate
  )

  return(posterior)
}


.hypothesis_brma_append_ordinate <- function(existing, ordinate) {

  if (is.null(existing)) {
    return(ordinate)
  }

  children <- c(
    .hypothesis_brma_ordinate_children(existing),
    list(ordinate)
  )
  out <- list(
    ordinates = children
  )
  class(out) <- c("RoBMA_posterior_ordinates", "list")

  return(out)
}


.hypothesis_brma_ordinate_children <- function(ordinate) {

  children <- ordinate[["ordinates"]]
  if (is.list(children) &&
      !is.data.frame(children) &&
      is.null(children[["x"]]) &&
      is.null(children[["value"]]) &&
      is.null(children[["null_hypothesis"]])) {
    return(children)
  }

  return(list(ordinate))
}


.hypothesis_brma_transform_ordinate <- function(ordinate, transform, value) {

  slope <- abs(transform[["slope"]])
  ordinate[["value"]]    <- value
  ordinate[["ordinate"]] <- ordinate[["ordinate"]] / slope

  diagnostics <- ordinate[["diagnostics"]]
  if (!is.null(diagnostics[["mcse"]])) {
    diagnostics[["mcse"]] <- diagnostics[["mcse"]] / slope
  }
  if (!is.null(diagnostics[["normalization_range"]])) {
    diagnostics[["normalization_range"]] <- transform[["intercept"]] +
      transform[["slope"]] * diagnostics[["normalization_range"]]
    diagnostics[["normalization_range"]] <- sort(diagnostics[["normalization_range"]])
  }
  diagnostics[["plot_scale_transform"]] <- transform
  ordinate[["diagnostics"]] <- diagnostics

  return(ordinate)
}


.hypothesis_brma_diagnostic_reason <- function(diagnostic) {

  reason <- diagnostic[["reason"]]
  if (is.null(reason) || !nzchar(reason)) {
    reason <- "failed BF-grade diagnostics"
  }

  return(reason)
}

