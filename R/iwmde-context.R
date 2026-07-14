# ============================================================================ #
# IWMDE Context and Availability
# ============================================================================ #

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
  return(.iwmde_context_ensure_caches(context))
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
  if (is.null(context[["source_fingerprint"]])) {
    context[["source_fingerprint"]] <-
      .iwmde_compute_source_fingerprint(context)
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
