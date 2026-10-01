# ============================================================================ #
# IWMDE Context and Availability
# ============================================================================ #

.iwmde_context <- function(object, integration_control = NULL) {

  posterior_samples <- as.matrix(.get_posterior_samples(object[["fit"]]))
  if (is.null(colnames(posterior_samples))) {
    stop("Posterior samples must have column names.", call. = FALSE)
  }
  chain_id <- .iwmde_chain_id(
    fit       = object[["fit"]],
    n_samples = nrow(posterior_samples)
  )

  data   <- object[["data"]]
  priors <- object[["priors"]]
  flat_prior_list <- attr(object[["fit"]], "prior_list")
  internal_priors <- BayesTools::JAGS_formula_internal_coordinate_priors(
    object[["fit"]]
  )
  collisions <- intersect(names(flat_prior_list), names(internal_priors))
  if (length(collisions) > 0L) {
    stop(
      "Internal formula-coordinate priors collide with fitted priors: ",
      paste(collisions, collapse = ", "), ".",
      call. = FALSE
    )
  }
  flat_prior_list <- c(flat_prior_list, internal_priors)

  context <- list(
    object            = object,
    data              = data,
    priors            = priors,
    posterior_samples = posterior_samples,
    chain_id          = chain_id,
    flat_prior_list   = flat_prior_list,
    selection_spec     = .iwmde_selection_spec(data, priors),
    formula_fit        = object[["fit"]],
    formula_inputs     = .iwmde_formula_inputs(data, priors),
    indicator_names    = .iwmde_indicator_names(
      posterior_samples,
      prior_list = flat_prior_list
    ),
    active_cache       = new.env(parent = emptyenv()),
    focal_prior_cache  = new.env(parent = emptyenv()),
    support_cache      = new.env(parent = emptyenv()),
    likelihood_cache   = new.env(parent = emptyenv()),
    prior_cache        = new.env(parent = emptyenv()),
    row_cache          = new.env(parent = emptyenv()),
    predictor_cache    = new.env(parent = emptyenv()),
    evaluator_cache    = new.env(parent = emptyenv())
  )

  class(context) <- "iwmde_context"
  context <- .iwmde_context_ensure_caches(context)
  return(.iwmde_context_with_integration_control(context, integration_control))
}


.iwmde_context_with_integration_control <- function(context, control) {

  if (is.null(control)) {
    return(context)
  }
  control <- .check_selection_likelihood_control(
    control, argument = "density_control$integration_control"
  )
  if (!.is_data_joint_selection(context[["data"]])) {
    stop(
      "'density_control$integration_control' is unavailable for this model: ",
      "post-fit integration controls require a bound Gaussian selection model.",
      call. = FALSE
    )
  }

  plan <- .data_selection_execution_plan(context[["data"]])
  execution_plan <- .selection_joint_execution_plan_with_control(plan, control)
  if (identical(execution_plan, plan)) {
    return(context)
  }

  attr(context[["data"]], "selection_execution_plan") <- execution_plan
  return(.iwmde_context_ensure_caches(context, reset = TRUE))
}


.iwmde_chain_id <- function(fit, n_samples) {

  if (is.null(fit) || length(fit) == 0L) {
    return(rep(1L, n_samples))
  }

  geometry <- BayesTools::JAGS_draw_geometry(fit)
  chain_lengths <- geometry[["chains"]][["iterations"]]
  if (geometry[["total_draws"]] != n_samples) {
    stop(
      "Posterior chain metadata does not match the materialized draws.",
      call. = FALSE
    )
  }

  return(rep(geometry[["chain_order"]], chain_lengths))
}


.iwmde_capability <- function(object = NULL, data = NULL,
                              density_method = NULL) {

  if (!is.null(object) && !is.null(object[["data"]])) {
    data <- object[["data"]]
  }
  if (!is.null(data) && .is_data_random(data) &&
      !.is_data_known_v(data)) {
    return(.iwmde_unavailable(
      reason = paste0(
        "qCMDE/IWMDE is not implemented for brma.mv() random-formula ",
        "models without known V yet."
      ),
      type   = "random_unknown_v"
    ))
  }
  scale_capability <- .iwmde_scale_capability(data)
  if (!is.null(scale_capability)) {
    return(scale_capability)
  }
  if (!is.null(density_method)) {
    density_method <- .density_method_normalize(density_method)
    is_glmm        <- inherits(object, "brma.glmm") ||
      (!is.null(data) &&
        isTRUE(.data_outcome_type(data) %in% c("bin", "pois")))
    if (is_glmm && identical(density_method, "IWMDE")) {
      return(.iwmde_unavailable(
        reason = paste0(
          "IWMDE density estimation is unavailable for binomial and ",
          "Poisson GLMMs. Use density_method = 'qCMDE'."
        ),
        type   = "glmm"
      ))
    }
  }

  return(list(available = TRUE, reason = ""))
}


# A qCMDE/IWMDE refusal: its reason and the classes of the condition it
# stops with (.iwmde_stop_unavailable()). The class
# "RoBMA_density_method_<type>" names the cause and the parent class
# "RoBMA_density_method_unavailable" (.iwmde_unavailable_class()) the
# family. The capability refusals of the fitted model are
# "random_unknown_v" (brma.mv() random-formula models without known V),
# "scale_components" (component-specific scale formulas), and "glmm" (IWMDE
# for binomial and Poisson GLMMs); hypothesis() refuses these causes with
# its own classes followed by these classes. plot() also refuses
# "conditional_random" (conditional random-effect plots), "random_target"
# (a random-effect quantity without a supported scalar random-component
# coordinate), and "original_scale" (an original-scale coefficient or factor
# cell that is no linear combination of its fitted coordinates, e.g. an
# exp(affine) scale intercept); hypothesis() refuses the same causes with
# its own classes followed by these classes as well. plot_prior() refuses
# "original_scale" for a term with several fitted coordinates whose
# original-scale coordinates the standardization changes
# (.plot_prior_check_coordinates_unchanged()).
.iwmde_unavailable <- function(reason, type) {

  list(
    available = FALSE,
    reason    = reason,
    class     = c(
      paste0("RoBMA_density_method_", type),
      .iwmde_unavailable_class()
    )
  )
}


.iwmde_unavailable_class <- function() {

  "RoBMA_density_method_unavailable"
}


# Stops with a capability refusal of .iwmde_capability(); 'caller' prefixes
# its reason.
.iwmde_stop_unavailable <- function(capability, caller = NULL) {

  message <- capability[["reason"]]
  if (!is.null(caller)) {
    message <- paste0(caller, ": ", message)
  }

  stop(structure(
    class = c(capability[["class"]], "error", "condition"),
    list(message = message, call = NULL)
  ))
}


# qCMDE/IWMDE evaluate the heterogeneity of one scale formula ('log_tau',
# .iwmde_formula_inputs()); component-specific scale formulas (a named
# 'scale' list of brma.mv(), one formula and 'log_tau_<component>' parameter
# per random component or block) are unavailable. NULL when the scale
# formula of 'data' is supported, otherwise the capability refusal.
.iwmde_scale_capability <- function(data) {

  if (is.null(data) || !.is_data_scale(data) ||
      !inherits(data[["scale"]], "RoBMA_scale_components")) {
    return(NULL)
  }
  formulas <- if (length(.data_scale_components(data)) > 1L) {
    "several scale formulas"
  } else {
    "a component-specific scale formula"
  }

  .iwmde_unavailable(
    reason = paste0(
      "qCMDE/IWMDE density estimation is unavailable for models with ",
      formulas, ". Use density_method = 'KDE'."
    ),
    type   = "scale_components"
  )
}


.check_iwmde_available <- function(object, caller) {

  capability <- .iwmde_capability(object = object)
  if (!capability[["available"]]) {
    .iwmde_stop_unavailable(capability, caller = caller)
  }

  return(invisible(TRUE))
}


.iwmde_context_unavailable_reason <- function(context) {

  object <- context[["object"]]
  data   <- context[["data"]]
  if (!is.null(object)) {
    data <- object[["data"]]
  }

  capability <- .iwmde_capability(object = object, data = data)
  if (capability[["available"]]) {
    return(NULL)
  }

  return(capability[["reason"]])
}


.iwmde_context_ensure_caches <- function(context, reset = FALSE) {

  cache_names <- c(
    "active_cache",
    "focal_prior_cache",
    "support_cache",
    "likelihood_cache",
    "prior_cache",
    "row_cache",
    "predictor_cache",
    "evaluator_cache"
  )
  for (cache_name in cache_names) {
    if (reset || !is.environment(context[[cache_name]])) {
      context[[cache_name]] <- new.env(parent = emptyenv())
    }
  }
  if (is.null(context[["indicator_names"]])) {
    context[["indicator_names"]] <- character()
  }
  if (is.null(context[["chain_id"]])) {
    context[["chain_id"]] <- rep(1L, nrow(context[["posterior_samples"]]))
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
  if (reset || is.null(context[["source_fingerprint"]])) {
    context[["source_fingerprint"]] <-
      .iwmde_compute_source_fingerprint(context)
  }

  return(context)
}


# The prepared BayesTools evaluator of the deterministic nodes 'nodes' of the
# fitted model (JAGS_deterministic_evaluator()), resolved once per context and
# node set. NULL without a fit or without nodes.
.iwmde_deterministic_evaluator <- function(context, nodes) {

  fit   <- context[["object"]][["fit"]]
  nodes <- unique(nodes)
  if (is.null(fit) || length(nodes) == 0L) {
    return(NULL)
  }

  cache <- context[["evaluator_cache"]]
  key   <- paste(nodes, collapse = "|")
  if (is.environment(cache) && !is.null(cache[[key]])) {
    return(cache[[key]])
  }
  evaluator <- BayesTools::JAGS_deterministic_evaluator(fit, nodes = nodes)
  if (is.environment(cache)) {
    cache[[key]] <- evaluator
  }

  return(evaluator)
}


.iwmde_indicator_names <- function(posterior_samples, prior_list = NULL) {

  indicator_names <- grep("(^|_)indicator$", colnames(posterior_samples), value = TRUE)
  indicator_names <- c(
    indicator_names,
    intersect("bias_indicator", colnames(posterior_samples))
  )
  indicator_names <- unique(indicator_names)
  if (!is.null(prior_list)) {
    random_gate_names <- unique(unname(
      .random_allocation_inclusion_indicators(prior_list)
    ))
    random_gate_names <- random_gate_names[nzchar(random_gate_names)]
    indicator_names <- setdiff(indicator_names, random_gate_names)
  }

  return(indicator_names)
}
