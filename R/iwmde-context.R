# ============================================================================ #
# IWMDE Context and Availability
# ============================================================================ #

.iwmde_context <- function(object) {

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
  lkj_pair_priors <- .iwmde_lkj_pair_prior_list(
    fit               = object[["fit"]],
    posterior_columns = colnames(posterior_samples)
  )
  collisions <- intersect(names(flat_prior_list), names(lkj_pair_priors))
  if (length(collisions) > 0L) {
    stop(
      "LKJ primitive prior names collide with fitted scalar priors: ",
      paste(collisions, collapse = ", "), ".",
      call. = FALSE
    )
  }
  flat_prior_list <- c(flat_prior_list, lkj_pair_priors)

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
    indicator_names    = .iwmde_indicator_names(posterior_samples),
    active_cache       = new.env(parent = emptyenv()),
    focal_prior_cache  = new.env(parent = emptyenv()),
    support_cache      = new.env(parent = emptyenv()),
    likelihood_cache   = new.env(parent = emptyenv()),
    prior_cache        = new.env(parent = emptyenv()),
    row_cache          = new.env(parent = emptyenv()),
    predictor_cache    = new.env(parent = emptyenv())
  )

  class(context) <- "iwmde_context"
  return(.iwmde_context_ensure_caches(context))
}


.iwmde_lkj_pair_prior_list <- function(fit, posterior_columns) {

  formula_design <- attr(fit, "formula_design", exact = TRUE)
  if (!is.list(formula_design)) {
    return(list())
  }
  terms <- unlist(lapply(formula_design, `[[`, "random_effects"), recursive = FALSE)
  out <- list()
  for (term in terms) {
    correlation <- term[["correlation"]]
    if (!is.list(correlation) || !identical(correlation[["type"]], "lkj") ||
        !identical(as.integer(term[["n_columns"]]), 2L)) {
      next
    }
    primitive <- correlation[["primitive_names"]]
    eta       <- correlation[["eta"]]
    if (!is.character(primitive) || length(primitive) != 1L ||
        is.na(primitive) || !nzchar(primitive) ||
        !primitive %in% posterior_columns ||
        !is.numeric(eta) || length(eta) != 1L || !is.finite(eta) || eta <= 0) {
      next
    }
    out[[primitive]] <- BayesTools::prior(
      "beta",
      parameters = list(alpha = eta, beta = eta)
    )
  }

  out
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
    return(list(
      available = FALSE,
      reason    = paste0(
        "qCMDE/IWMDE is not implemented for brma.mv() random-formula ",
        "models without known V yet."
      )
    ))
  }
  if (!is.null(density_method)) {
    density_method <- .density_method_normalize(density_method)
    is_glmm        <- inherits(object, "brma.glmm") ||
      (!is.null(data) &&
        isTRUE(.data_outcome_type(data) %in% c("bin", "pois")))
    if (is_glmm && identical(density_method, "IWMDE")) {
      return(list(
        available = FALSE,
        reason    = paste0(
          "IWMDE density estimation is unavailable for binomial and ",
          "Poisson GLMMs. Use density_method = 'qCMDE'."
        )
      ))
    }
  }

  return(list(available = TRUE, reason = ""))
}


.check_iwmde_available <- function(object, caller) {

  capability <- .iwmde_capability(object = object)
  if (!capability[["available"]]) {
    stop(caller, ": ", capability[["reason"]], call. = FALSE)
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


.iwmde_context_ensure_caches <- function(context) {

  cache_names <- c(
    "active_cache",
    "focal_prior_cache",
    "support_cache",
    "likelihood_cache",
    "prior_cache",
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
