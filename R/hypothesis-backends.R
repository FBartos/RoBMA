# ============================================================================ #
# hypothesis-backends.R
# ============================================================================ #
#
# The qCMDE and IWMDE density backends of the BayesTools hypothesis engine
# (BayesTools::hypothesis_BF_model()). A backend attaches likelihood-aware
# posterior ordinates to the draws of a tested quantity; BayesTools owns the
# bookkeeping (which stored ordinates are kept, the evaluation of the Bayes
# factors, the lifting of warnings and diagnostics to the rows of the table),
# RoBMA owns the estimators (R/iwmde*.R), their acceptance rules and thresholds
# (R/iwmde-control.R), their refusals, and their typed failures
# ("RoBMA_density_ordinate_error"). The providers that describe the models to
# the engine are in hypothesis-provider.R.
#
# ============================================================================ #


# The backends of RoBMA, in the order in which they are registered and
# advertised by hypothesis_quantities().
.hypothesis_density_backend_names <- function() {

  c("qCMDE", "IWMDE")
}


.hypothesis_density_backend <- function(name) {

  BayesTools::hypothesis_density_backend(
    name             = name,
    type             = "ordinate",
    available        = function(object, plan, control) {
      capability <- .iwmde_capability(object = object, density_method = name)
      if (capability[["available"]]) {
        return(NULL)
      }
      list(
        reason = capability[["reason"]],
        type   = "method",
        class  = capability[["class"]]
      )
    },
    control          = function(control) {
      .hypothesis_density_control_fill(.iwmde_density_control_resolve(
        density_method  = name,
        density_control = control,
        purpose         = "ordinate"
      ))
    },
    attach           = function(posterior, request, control, cache) {
      .hypothesis_density_attach(posterior, request, control, cache)
    },
    reuse            = function(entry, request, control, posterior = NULL) {
      .hypothesis_density_reuse(entry, request, control, posterior)
    },
    warnings         = function(posterior, request) {
      .iwmde_collect_posterior_ordinate_warning_records(posterior)
    },
    diagnostics      = function(posterior, request) {
      .iwmde_collect_public_density_diagnostics(posterior)
    },
    bind_diagnostics = .iwmde_bind_density_diagnostics,
    rescale_failure  = function(condition, requested, evaluated, label,
                                density_method) {
      .hypothesis_density_rescale_failure(
        condition, requested, evaluated, label, density_method
      )
    },
    refusal_classes  = function(cause) {
      c(paste0("RoBMA_density_method_", cause), .iwmde_unavailable_class())
    },
    clear_cache      = .iwmde_qcmde_cache_clear,
    policy           = c(
      .hypothesis_density_policy(name),
      list(instance = .hypothesis_density_backends_instance())
    ),
    package          = "RoBMA"
  )
}


# The thresholds of the acceptance of an ordinate that the backend applies
# (R/iwmde-control.R), as the backend declares them to the registry.
.hypothesis_density_policy <- function(density_method) {

  estimator <- .density_method_iwmde_estimator(density_method)

  list(
    ordinate_error_warning     = .iwmde_bf_mass_warning_tolerance(estimator),
    ordinate_error_rejection   = .iwmde_bf_mass_fail_tolerance(estimator),
    relative_mcse_warning      = .iwmde_bf_warning_relative_mcse(),
    min_ess_warning            = .iwmde_bf_warning_min_ess(),
    max_weight_share_warning   = .iwmde_bf_warning_weight_share(),
    min_finite_terms_warning   = .iwmde_bf_warning_min_finite_terms(),
    quadrature_change_warning  = .iwmde_quadrature_warning_tolerance(),
    quadrature_change_rejection = .iwmde_quadrature_fail_tolerance()
  )
}


# The identity of this namespace instance of RoBMA, declared with the backends
# it registers (as the "instance" of their policy) so that the registry shows
# which instance provides them. It is made once per namespace instance and uses
# no random numbers.
.hypothesis_density_backends_instance <- function() {

  instance <- RoBMA.private[["hypothesis_backends_instance"]]
  if (is.null(instance)) {
    instance <- paste(
      format(Sys.time(), "%Y%m%d%H%M%OS6"), Sys.getpid(),
      proc.time()[["elapsed"]], sep = "-"
    )
    RoBMA.private[["hypothesis_backends_instance"]] <- instance
  }

  instance
}


# Registers the backends of RoBMA (replacing the ones of a previous load) and
# drops the plans of the engine, which were built with the previous code.
.hypothesis_density_backends_register <- function() {

  for (name in .hypothesis_density_backend_names()) {
    BayesTools::register_hypothesis_density_backend(
      .hypothesis_density_backend(name),
      overwrite = TRUE
    )
  }
  BayesTools::hypothesis_cache_clear()

  invisible(NULL)
}


# Registers the backends again when the registry of BayesTools lacks them or
# holds the ones of another instance of the RoBMA namespace: reloading
# BayesTools alone empties its registry, while RoBMA registers them when it
# loads. The providers call it, so that every call of the engine finds them.
# Nothing happens (and no plan is dropped) while the registry holds the
# backends of this instance. Returns whether they were registered.
.hypothesis_density_backends_ensure <- function() {

  registered <- BayesTools::hypothesis_density_backends()
  instance   <- .hypothesis_density_backends_instance()
  current    <- vapply(.hypothesis_density_backend_names(), function(name) {
    row <- match(name, registered[["name"]])
    !is.na(row) && identical(registered[["package"]][[row]], "RoBMA") &&
      identical(registered[["policy"]][[row]][["instance"]], instance)
  }, logical(1))
  if (all(current)) {
    return(invisible(FALSE))
  }
  .hypothesis_density_backends_register()

  invisible(TRUE)
}


.hypothesis_density_backends_unregister <- function() {

  for (name in .hypothesis_density_backend_names()) {
    BayesTools::unregister_hypothesis_density_backend(name)
  }

  invisible(NULL)
}


# The grid of the normalization of a point ordinate defaults to the points of
# the density, at least 50.
.hypothesis_density_control_fill <- function(density_control) {

  if (is.null(density_control[["normalization_points"]])) {
    density_control[["normalization_points"]] <- max(
      50L,
      density_control[["n_points"]]
    )
  }

  density_control
}


# ---------------------------------------------------------------------------- #
# Attachment of ordinates
# ---------------------------------------------------------------------------- #

# The ordinates of the requested values of one tested quantity: a scalar
# posterior (also the linear target of a combination of levels, or a random
# effect), a level of the marginal posterior of a factor term, or a level of
# the stored marginal means. The context of the fitted model and the cache of
# the estimates live in the cache of the group of statements.
.hypothesis_density_attach <- function(posterior, request, control, cache) {

  density_method <- request[["density_method"]]
  control        <- .hypothesis_density_control_fill(control)
  .iwmde_check_point_ordinate_supported(request[["object"]], density_method)
  if (is.null(cache[["context"]])) {
    cache[["context"]]        <- .iwmde_context(
      request[["object"]], control[["integration_control"]]
    )
    cache[["estimate_cache"]] <- .iwmde_estimate_cache()
  }

  if (identical(request[["kind"]], "scalar")) {
    return(.hypothesis_density_attach_scalar(
      posterior, request, control, cache[["context"]], cache[["estimate_cache"]]
    ))
  }
  if (identical(request[["route"]], "marginal_mean")) {
    return(.hypothesis_density_attach_marginal_level(
      posterior, request, control, cache[["context"]], cache[["estimate_cache"]]
    ))
  }

  .hypothesis_density_attach_level(
    posterior, request, control, cache[["context"]], cache[["estimate_cache"]]
  )
}


# The control with which the estimator computes a point ordinate of a fitted
# object: the setting of the display grid is immaterial for point ordinates and
# marked as such.
.hypothesis_density_ordinate_control <- function(control) {

  list(
    n_points             = control[["n_points"]],
    samples              = control[["samples"]],
    target_relative_mcse = control[["target_relative_mcse"]],
    normalization_points = control[["normalization_points"]],
    normalization_prob   = control[["normalization_prob"]],
    integration_control  = control[["integration_control"]],
    display_grid         = "ordinate"
  )
}


.hypothesis_density_attach_scalar <- function(posterior, request, control,
                                              context, estimate_cache) {

  parameter         <- request[["parameter"]]
  parameter_label   <- request[["label"]]
  value             <- request[["values"]]
  conditional       <- request[["conditional"]]
  density_method    <- request[["density_method"]]
  parameter_spec    <- request[["spec"]]
  display_transform <- request[["display_transform"]]

  if (is.list(posterior)) {
    stop(
      "qCMDE/IWMDE point hypotheses for factor parameters must specify ",
      "a level, e.g. '", parameter_label, "[level] = ", value, "'.",
      call. = FALSE
    )
  }

  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "unsupported_formula_transform")) {
    .iwmde_stop_unavailable(.iwmde_unavailable(parameter_spec[["reason"]], "original_scale"))
  }
  identity_transform <- is.null(display_transform) || identical(display_transform[["type"]], "identity")
  if (is.null(parameter_spec)) {
    parameter_spec <- list(
      type          = "primitive",
      prior_density = if (identity_transform)
        BayesTools::posterior_metadata(posterior, "prior_density")
    )
  } else if (is.null(parameter_spec[["prior_density"]]) && identity_transform) {
    parameter_spec[["prior_density"]] <- BayesTools::posterior_metadata(
      posterior,
      "prior_density"
    )
  }
  condition <- .iwmde_sample_condition_metadata(posterior)
  if(length(conditional) > 0L && length(condition[["conditional"]]) > 0L &&
     !setequal(conditional, condition[["conditional"]])){
    stop("The scalar ordinate condition does not match its posterior declaration. Recompute the posterior for the requested condition.", call. = FALSE)
  }
  for(field in c("conditional", "conditional_rule", "condition_key", "condition_event", "resolved_condition_event")) parameter_spec[[field]] <- NULL
  parameter_spec[["conditional"]] <- conditional
  complete_condition <- length(condition[["conditional"]]) > 0L &&
    length(condition[["conditional_rule"]]) > 0L && length(condition[["condition_key"]]) > 0L
  if(length(conditional) > 0L && complete_condition){
    for(field in names(condition)) parameter_spec[[field]] <- condition[[field]]
  }else if(length(conditional) > 0L) parameter_spec[["conditional_rule"]] <- "AND"
  source_value <- if (is.null(display_transform)) {
    value
  } else {
    BayesTools::parameter_transform_inverse(value, display_transform)
  }
  source_keys <- vapply(source_value, .iwmde_key_number, character(1))

  estimate <- tryCatch(.iwmde_estimate(
    context         = context,
    parameter       = parameter,
    density_method  = density_method,
    density_control = .hypothesis_density_ordinate_control(control),
    outputs        = "ordinate",
    values         = source_value,
    parameter_spec = parameter_spec,
    metadata       = .iwmde_posterior_metadata(
      samples   = posterior,
      parameter = parameter
    ),
    cache          = estimate_cache
  ), RoBMA_density_ordinate_error = function(error) {
    if (!is.null(display_transform)) {
      error[["density_diagnostics"]] <- .hypothesis_brma_transform_iwmde_failures(
        error[["density_diagnostics"]], display_transform,
        source_values = source_value, requested_values = value
      )
    }
    stop(error)
  })
  diagnostic <- estimate[["diagnostics"]][["ordinate"]]
  marginal_parameter <- attr(posterior, "parameter", exact = TRUE)
  values <- .iwmde_sorted_ordinate_values(source_value)
  for (requested_source_value in values) {
    requested_index <- match(
      .iwmde_key_number(requested_source_value),
      source_keys
    )
    requested_value <- value[[requested_index]]
    ordinate <- .iwmde_posterior_ordinate_keep_values(
      posterior_ordinate = estimate[["posterior_ordinate"]],
      values             = requested_source_value
    )
    if (is.null(ordinate)) {
      failure_estimate <- estimate
      failure_value <- requested_source_value
      if (!is.null(display_transform)) {
        failure_estimate[["ordinate_failures"]] <- .hypothesis_brma_transform_iwmde_failures(
          estimate[["ordinate_failures"]], display_transform,
          source_values = source_value, requested_values = value
        )
        failure_value <- requested_value
      }
      .iwmde_stop_ordinate_unavailable(
        message = .hypothesis_brma_iwmde_ordinate_failure_message(
          density_method = density_method,
          target         = paste0(parameter_label, " = ", requested_value),
          diagnostic     = diagnostic,
          reason         = .hypothesis_brma_estimate_ordinate_reason(
            estimate = failure_estimate,
            value    = failure_value
          ),
          numerical_status = .iwmde_ordinate_failure_record(
            failure_estimate, failure_value)[["numerical_status"]]
        ),
        estimate = failure_estimate
      )
    }
    if (!is.null(display_transform)) {
      ordinate <- .hypothesis_brma_transform_iwmde_ordinate(
        ordinate,
        display_transform
      )
      ordinate[["value"]] <- requested_value
    }
    if (is.character(marginal_parameter) &&
        length(marginal_parameter) == 1L &&
        !is.na(marginal_parameter) && nzchar(marginal_parameter)) {
      ordinate[["parameter"]] <- marginal_parameter
    }

    existing <- .iwmde_posterior_ordinate_drop_value(
      posterior_ordinate = BayesTools::posterior_metadata(
        posterior,
        "posterior_ordinate"
      ),
      value              = requested_value
    )
    BayesTools::posterior_metadata(posterior, "posterior_ordinate") <-
      BayesTools::posterior_ordinate_append(
        existing = existing,
        ordinate = ordinate
      )
  }

  return(posterior)
}


# One level of the marginal posterior of a factor term of a fitted object.
.hypothesis_density_attach_level <- function(posterior, request, control,
                                             context, estimate_cache) {

  parameter      <- request[["parameter"]]
  level          <- request[["level"]]
  value          <- request[["values"]]
  conditional    <- request[["conditional"]]
  density_method <- request[["density_method"]]
  parameter_spec <- request[["spec"]]

  if (!is.list(posterior) || !level %in% names(posterior)) {
    stop("Hypothesis references unknown level '", level,
         "' for parameter '", parameter, "'.", call. = FALSE)
  }
  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "unsupported_formula_transform")) {
    .iwmde_stop_unavailable(.iwmde_unavailable(parameter_spec[["reason"]], "original_scale"))
  }

  weights <- if (is.null(parameter_spec)) {
    BayesTools::posterior_metadata(posterior[[level]], "linear_weights")
  } else {
    parameter_spec[["weights"]]
  }
  linear_weights <- .iwmde_linear_weights(weights)
  if(is.null(parameter_spec) && !is.null(linear_weights)){
    linear_weights <- .iwmde_fitted_linear_weights(posterior[[level]], linear_weights)
    weights <- linear_weights
  }
  if (is.null(linear_weights) || length(linear_weights) == 0L) {
    reason <- if (is.null(linear_weights)) {
      "linear weights are unavailable"
    } else {
      "linear weights are all zero"
    }
    stop(
      "Precomputed ", density_method, " posterior ordinate is unavailable for '",
      parameter, "[", level, "] = ", value, "': ", reason,
      call. = FALSE
    )
  }

  if (is.null(parameter_spec)) {
    parameter_spec <- list(
      type          = "linear",
      weights       = weights,
      prior_density = BayesTools::posterior_metadata(
        posterior[[level]],
        "prior_density"
      )
    )
  }
  parameter_spec[["conditional"]] <- if (is.null(conditional)) {
    BayesTools::posterior_metadata(
      posterior[[level]],
      "condition"
    )[["effective_conditional"]]
  } else {
    conditional
  }
  parameter_spec[["conditional_rule"]] <- "AND"

  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = paste0(parameter, "[", level, "]"),
    density_method  = density_method,
    density_control = .hypothesis_density_ordinate_control(control),
    outputs        = "ordinate",
    values         = value,
    parameter_spec = parameter_spec,
    metadata       = .iwmde_posterior_metadata(
      samples   = posterior[[level]],
      parameter = paste0(parameter, "[", level, "]"),
      level     = level
    ),
    cache          = estimate_cache
  )
  diagnostic <- estimate[["diagnostics"]][["ordinate"]]

  ordinate <- estimate[["posterior_ordinate"]]
  if (is.null(ordinate)) {
    .iwmde_stop_ordinate_unavailable(
      message = .hypothesis_brma_iwmde_ordinate_failure_message(
        density_method = density_method,
        target         = paste0(parameter, "[", level, "] = ", value),
        diagnostic     = diagnostic,
        reason         = .hypothesis_brma_estimate_ordinate_reason(estimate, value),
        numerical_status = .iwmde_ordinate_failure_record(estimate, value)[["numerical_status"]]
      ),
      estimate = estimate
    )
  }

  existing <- .iwmde_posterior_ordinate_drop_value(
    posterior_ordinate = BayesTools::posterior_metadata(
      posterior[[level]],
      "posterior_ordinate"
    ),
    value              = value
  )
  BayesTools::posterior_metadata(
    posterior[[level]],
    "posterior_ordinate"
  ) <- BayesTools::posterior_ordinate_append(
    existing = existing,
    ordinate = ordinate
  )

  return(posterior)
}


# The marginal-means object of a request, with the marginal means of the
# request in their current state (the stored ordinates a request keeps).
.hypothesis_density_marginal_object <- function(request, posterior) {

  object <- request[["model"]]
  if (is.list(posterior) && !is.null(names(posterior))) {
    object[["inference"]][[request[["inference_type"]]]][[request[["parameter"]]]] <-
      posterior
  }

  object
}


# One level of the stored marginal means of a marginal-means object: the
# ordinate is computed from its source fitted model (the object of the request).
.hypothesis_density_attach_marginal_level <- function(posterior, request,
                                                      control, context,
                                                      estimate_cache) {

  parameter      <- request[["parameter"]]
  level          <- request[["level"]]
  value          <- request[["values"]]
  density_method <- request[["density_method"]]
  label          <- paste0(parameter, "[", level, "]")
  sample         <- posterior[[level]]
  existing       <- BayesTools::posterior_metadata(sample, "posterior_ordinate")

  specs <- .iwmde_marginal_means_specs(
    marginal_means_object = .hypothesis_density_marginal_object(request, posterior),
    parameter             = parameter,
    type                  = request[["inference_type"]],
    levels                = level
  )
  if (length(specs) != 1L) {
    stop(
      "Could not construct a qCMDE/IWMDE target for '", label, "'.",
      call. = FALSE
    )
  }
  spec <- specs[[1L]]
  metadata <- .iwmde_posterior_metadata(
    samples   = sample,
    parameter = label,
    level     = spec[["level"]]
  )
  existing <- .iwmde_posterior_ordinate_drop_value(existing, value)

  estimate <- .iwmde_estimate(
    context         = context,
    parameter       = spec[["label"]],
    density_method  = density_method,
    density_control = control,
    outputs         = "ordinate",
    values          = value,
    parameter_spec  = spec,
    metadata        = metadata,
    cache           = estimate_cache
  )
  diagnostic <- estimate[["diagnostics"]][["ordinate"]]
  ordinate   <- estimate[["posterior_ordinate"]]
  if (is.null(ordinate)) {
    .iwmde_stop_ordinate_unavailable(
      message = .hypothesis_brma_iwmde_ordinate_failure_message(
        density_method = density_method,
        target         = paste0(label, " = ", value),
        diagnostic     = diagnostic,
        reason         = .hypothesis_brma_estimate_ordinate_reason(estimate, value),
        numerical_status = .iwmde_ordinate_failure_record(estimate, value)[["numerical_status"]]
      ),
      estimate = estimate
    )
  }

  BayesTools::posterior_metadata(
    sample,
    "posterior_ordinate"
  ) <- BayesTools::posterior_ordinate_append(
    existing = existing,
    ordinate = ordinate
  )
  posterior[[level]] <- sample

  return(posterior)
}


# Whether an ordinate stored on the stored marginal means of a marginal-means
# object is the ordinate that the request would compute (the same method, control,
# source of the fit, target, and conditioning), so that it is not computed
# again. The ordinates of fitted objects are always computed again.
.hypothesis_density_reuse <- function(entry, request, control, posterior) {

  if (!identical(request[["route"]], "marginal_mean") || is.null(posterior) ||
      !identical(request[["kind"]], "level")) {
    return(FALSE)
  }
  level  <- request[["level"]]
  sample <- if (is.list(posterior) && !is.null(names(posterior))) {
    posterior[[level]]
  } else {
    posterior
  }
  if (is.null(sample)) {
    return(FALSE)
  }
  label <- paste0(request[["parameter"]], "[", level, "]")
  specs <- .iwmde_marginal_means_specs(
    marginal_means_object = .hypothesis_density_marginal_object(request, posterior),
    parameter             = request[["parameter"]],
    type                  = request[["inference_type"]],
    levels                = level
  )
  if (length(specs) != 1L) {
    return(FALSE)
  }
  spec     <- specs[[1L]]
  value    <- suppressWarnings(as.numeric(entry[["value"]]))[1L]
  context  <- .iwmde_context(request[["object"]], control[["integration_control"]])
  metadata <- .iwmde_posterior_metadata(
    samples   = sample,
    parameter = label,
    level     = spec[["level"]]
  )
  expected_provenance <- .iwmde_request_provenance(
    context         = context,
    parameter       = spec[["label"]],
    density_method  = request[["density_method"]],
    density_control = control,
    attribute       = "ordinate",
    value           = value,
    parameter_spec  = spec,
    metadata        = metadata
  )

  .iwmde_posterior_ordinate_matches_request(
    posterior_ordinate = entry,
    value              = value,
    provenance         = expected_provenance
  )
}


# A variance tested through its standard deviation: the failure of the
# ordinates of the standard deviations is restated for the variances that were
# requested (their values, the scales, and the message). A condition of another
# kind is not restated.
.hypothesis_density_rescale_failure <- function(condition, requested, evaluated,
                                                label, density_method) {

  if (!inherits(condition, "RoBMA_density_ordinate_error")) {
    return(NULL)
  }
  diagnostics <- condition[["density_diagnostics"]]
  if (!is.data.frame(diagnostics) || nrow(diagnostics) == 0L) return(NULL)
  index <- match(diagnostics[["requested_value"]], evaluated)
  matched <- !is.na(index)
  if (!any(matched)) return(NULL)
  diagnostics[["requested_value"]][matched] <- requested[index[matched]]
  if (is.null(diagnostics[["requested_scale"]])) diagnostics[["requested_scale"]] <- NA_character_
  if (is.null(diagnostics[["evaluation_scale"]])) diagnostics[["evaluation_scale"]] <- NA_character_
  diagnostics[["requested_scale"]][matched] <- "variance"
  diagnostics[["evaluation_scale"]][matched] <- "sd"
  condition[["density_diagnostics"]] <- diagnostics
  first <- which(matched)[[1L]]
  condition[["message"]] <- .hypothesis_brma_iwmde_ordinate_failure_message(
    density_method = density_method,
    target = paste0(label, " = ", diagnostics[["requested_value"]][[first]]),
    diagnostic = list(status = if (diagnostics[["status"]][[first]] == "rejected") "ok" else "unavailable"),
    reason = diagnostics[["failure_reason"]][[first]],
    numerical_status = diagnostics[["numerical_status"]][[first]]
  )

  condition
}


# ---------------------------------------------------------------------------- #
# Display transforms and failure messages of the ordinates
# ---------------------------------------------------------------------------- #

.hypothesis_brma_transform_iwmde_ordinate <- function(
    ordinate, display_transform) {

  source_evaluation_value <- ordinate[["evaluation_value"]]
  jacobian <- BayesTools::parameter_transform_jacobian(
    source_evaluation_value,
    display_transform
  )
  if (length(jacobian) != 1L || !is.finite(jacobian) || jacobian <= 0) {
    stop("Unsupported qCMDE/IWMDE display transform.", call. = FALSE)
  }
  ordinate[["value"]] <- BayesTools::parameter_transform_forward(
    ordinate[["value"]],
    display_transform
  )
  ordinate[["evaluation_value"]] <- BayesTools::parameter_transform_forward(
    ordinate[["evaluation_value"]],
    display_transform
  )
  ordinate[["ordinate"]] <- ordinate[["ordinate"]] / jacobian
  diagnostics <- ordinate[["diagnostics"]]
  if (!is.null(diagnostics[["log_ordinate"]])) {
    diagnostics[["log_ordinate"]] <- diagnostics[["log_ordinate"]] - log(jacobian)
  }
  value_fields <- intersect(
    c("evaluation_value"),
    names(diagnostics)
  )
  for (field in value_fields) {
    diagnostics[[field]] <- BayesTools::parameter_transform_forward(
      diagnostics[[field]],
      display_transform
    )
  }
  density_fields <- intersect(
    c(
      "mcse", "active_branch_mcse", "active_mass_component_mcse", "sampling_mcse",
      "pilot_ordinate", "validation_ordinate"
    ),
    names(diagnostics)
  )
  for (field in density_fields) {
    diagnostics[[field]] <- diagnostics[[field]] / jacobian
  }
  ordinate[["diagnostics"]] <- diagnostics
  provenance <- ordinate[["iwmde_provenance"]]
  if (is.list(provenance)) {
    for (field in intersect(c("value", "evaluation_value"), names(provenance))) {
      provenance[[field]] <- BayesTools::parameter_transform_forward(
        provenance[[field]],
        display_transform
      )
    }
    ordinate[["iwmde_provenance"]] <- provenance
  }

  if (!is.finite(ordinate[["ordinate"]]) || ordinate[["ordinate"]] <= 0) {
    entry <- unclass(ordinate)
    entry[["computed"]] <- TRUE
    entry[["diagnostics"]][["evaluation_value"]] <- ordinate[["evaluation_value"]]
    failure <- .iwmde_ordinate_failure_row(entry, preserve_accuracy = TRUE)
    .iwmde_stop_ordinate_unavailable(
      paste0("Transformed posterior ordinate is unavailable: ", failure[["failure_reason"]]),
      list(ordinate_failures = failure)
    )
  }

  ordinate
}


.hypothesis_brma_transform_iwmde_failures <- function(records, transform,
                                                       source_values = NULL,
                                                       requested_values = NULL) {

  if (is.null(records) || !nrow(records)) return(records)
  sources <- records[["requested_value"]]
  records[["requested_value"]] <- BayesTools::parameter_transform_forward(sources, transform)
  if (!is.null(source_values) && !is.null(requested_values)) {
    index <- match(sources, source_values)
    matched <- !is.na(index)
    records[["requested_value"]][matched] <- requested_values[index[matched]]
  }
  for (i in seq_len(nrow(records))) {
    value <- records[["evaluation_value"]][[i]]
    if (!is.finite(value)) next
    jacobian <- BayesTools::parameter_transform_jacobian(value, transform)
    if (length(jacobian) != 1L || !is.finite(jacobian) || jacobian <= 0) {
      stop("Unsupported qCMDE/IWMDE display transform.", call. = FALSE)
    }
    records[["evaluation_value"]][[i]] <- BayesTools::parameter_transform_forward(value, transform)
    records[["ordinate"]][[i]] <- records[["ordinate"]][[i]] / jacobian
    records[["log_ordinate"]][[i]] <- records[["log_ordinate"]][[i]] - log(jacobian)
    previous <- records[["numerical_status"]][[i]]
    current <- .iwmde_ordinate_numerical_status(records[["ordinate"]][[i]], records[["log_ordinate"]][[i]])
    records[["numerical_status"]][[i]] <- current
    if (!identical(current, "finite")) records[["status"]][[i]] <- "unavailable"
    if (!identical(previous, "finite") || !identical(current, "finite")) {
      records[["failure_reason"]][[i]] <- .iwmde_ordinate_numerical_reason(current,
        records[["log_ordinate"]][[i]], records[["failure_reason"]][[i]])
    }
  }
  records
}


.hypothesis_brma_estimate_ordinate_reason <- function(estimate, value) {

  failure <- .iwmde_ordinate_failure_record(estimate, value)
  if (!is.null(failure)) return(failure[["failure_reason"]])
  rejected <- .iwmde_posterior_ordinate_keep_values(
    posterior_ordinate = estimate[["rejected_posterior_ordinate"]],
    values             = value
  )
  if (!is.null(rejected)) {
    reason <- .iwmde_posterior_ordinate_bf_failure_reason(rejected)
    if (!is.null(reason) && nzchar(reason)) {
      return(reason)
    }
  }

  return(.hypothesis_brma_diagnostic_reason(
    estimate[["diagnostics"]][["ordinate"]]
  ))
}


.hypothesis_brma_iwmde_ordinate_failure_message <- function(
    density_method, target, diagnostic, reason, numerical_status = NULL) {

  reason <- sub("[.]$", "", reason)
  # A reason that already tells the reader to inspect the diagnostics keeps that
  # sentence; the remedy then adds only what the reason does not say.
  inspected <- grepl("inspect", reason, ignore.case = TRUE)
  remedy <- if (!is.null(numerical_status) && !identical(numerical_status, "finite")) {
    paste0(" ", if (inspected) "Use" else "Inspect 'density_diagnostics()'; use",
           " a region hypothesis when a point ordinate is not representable.")
  } else if (inspected) {
    ""
  } else {
    " Inspect 'density_diagnostics()' for the failed criterion and its computation settings."
  }
  if (identical(diagnostic[["status"]], "ok") &&
      (is.null(numerical_status) || identical(numerical_status, "finite"))) {
    return(paste0(
      density_method, " posterior ordinate for '", target,
      "' was rejected by diagnostics: ", reason, ".", remedy
    ))
  }

  return(paste0(
    density_method, " posterior ordinate is unavailable for '", target,
    "': ", reason, ".", remedy
  ))
}


.hypothesis_brma_diagnostic_reason <- function(diagnostic) {

  reason <- diagnostic[["reason"]]
  if (is.null(reason) || !nzchar(reason)) {
    diagnostic_reason <- .iwmde_diagnostics_bf_failure_reason(
      diagnostic[["diagnostics"]]
    )
    if (!is.null(diagnostic_reason) && nzchar(diagnostic_reason)) {
      reason <- diagnostic_reason
    } else {
      reason <- "failed qCMDE/IWMDE numerical diagnostics"
    }
  }

  return(reason)
}


# ---------------------------------------------------------------------------- #
# Warnings of the stored ordinates
# ---------------------------------------------------------------------------- #

# One record per warning of the ordinates stored on 'posterior' (recursively
# over a list of level posteriors): the warning text, the parameter, level,
# label, and value of the ordinate. The engine puts each on the rows of the
# statements it belongs to.
.iwmde_collect_posterior_ordinate_warning_records <- function(posterior) {

  records <- .iwmde_posterior_ordinate_warning_records(
    BayesTools::posterior_metadata(posterior, "posterior_ordinate")
  )
  if (is.list(posterior)) {
    child_records <- lapply(posterior, function(sample) {
      .iwmde_collect_posterior_ordinate_warning_records(sample)
    })
    child_records <- child_records[vapply(child_records, nrow, integer(1)) > 0L]
    if (length(child_records) > 0L) {
      records <- rbind(records, do.call(rbind, child_records))
    }
  }

  records <- records[!duplicated(records), , drop = FALSE]
  rownames(records) <- NULL

  return(records)
}


.iwmde_posterior_ordinate_warning_records <- function(posterior_ordinate) {

  entries <- .iwmde_posterior_ordinate_entries(posterior_ordinate)
  if (length(entries) == 0L) {
    return(.iwmde_empty_posterior_ordinate_warning_records())
  }

  records <- lapply(entries, function(entry) {

    warnings <- .iwmde_posterior_ordinate_warnings(entry)
    if (length(warnings) == 0L) {
      return(.iwmde_empty_posterior_ordinate_warning_records())
    }

    data.frame(
      warning   = warnings,
      parameter = .hypothesis_brma_ordinate_character(entry, "parameter"),
      level     = .hypothesis_brma_ordinate_character(entry, "level"),
      label     = .hypothesis_brma_ordinate_label(entry),
      value     = .hypothesis_brma_ordinate_value(entry),
      stringsAsFactors = FALSE
    )
  })
  records <- records[vapply(records, nrow, integer(1)) > 0L]
  if (length(records) == 0L) {
    return(.iwmde_empty_posterior_ordinate_warning_records())
  }

  do.call(rbind, records)
}


.iwmde_empty_posterior_ordinate_warning_records <- function() {

  data.frame(
    warning   = character(),
    parameter = character(),
    level     = character(),
    label     = character(),
    value     = numeric(),
    stringsAsFactors = FALSE
  )
}


.hypothesis_brma_ordinate_character <- function(entry, name) {

  value <- entry[[name]]
  if (is.null(value) || length(value) == 0L) {
    return(NA_character_)
  }

  value <- as.character(value[[1L]])
  if (is.na(value) || !nzchar(value)) {
    return(NA_character_)
  }

  return(value)
}


.hypothesis_brma_ordinate_label <- function(entry) {

  parameter <- .hypothesis_brma_ordinate_character(entry, "parameter")
  level     <- .hypothesis_brma_ordinate_character(entry, "level")
  if (!is.na(parameter) && nzchar(parameter) &&
      !is.na(level) && nzchar(level)) {
    return(paste0(parameter, "[", level, "]"))
  }
  if (!is.na(parameter) && nzchar(parameter)) {
    return(parameter)
  }

  return(NA_character_)
}


.hypothesis_brma_ordinate_value <- function(entry) {

  for (name in c("value", "x", "null_hypothesis")) {
    value <- suppressWarnings(as.numeric(entry[[name]]))
    value <- value[is.finite(value)]
    if (length(value) > 0L) {
      return(value[[1L]])
    }
  }

  return(NA_real_)
}
