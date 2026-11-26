# ============================================================================ #
# sampler-joint-location-scale.R
# ============================================================================ #
#
# Heterogeneity scale interweaving in the location block sampler
# (src/samplers/LocationScale.cc): when the block owns the standardized latents
# of the sampled random-effect term, it can also own the heterogeneity sources
# (the total scale a, the weight sources and their gates' targets), give them
# their ordinary samplers unchanged and follow them with an exact scale-group
# move of (a, w, Z) that leaves every child density unchanged. R declares the
# structure from the bound priors and the formula design; the native factory
# verifies it on the compiled graph and declines otherwise.
#
# ============================================================================ #


# The scale marker of a location declaration with latents, or NULL when the
# heterogeneity has no supported layout: the allocation layout of the `.mv`
# interfaces (total SD, inclusion indicator, Dirichlet weights held as gamma
# sources) or the tau and rho of clustered fits.
.sampler_location_scale <- function(prior_list, design, latents) {

  if (length(latents) == 0L) return(NULL)
  bias <- 0
  if (BayesTools::is.prior.mixture(prior_list[["bias"]])) {
    bias <- tryCatch(BayesTools::JAGS_mixture_sources(prior_list, "bias")[["gate"]], error = function(e) NULL)
    if (!is.character(bias) || length(bias) != 1L) return(NULL)
  }
  sampled <- if (is.null(design)) list() else
    Filter(function(term) identical(term[["compile_mode"]], "sampled"), design[["random_effects"]])
  layout <- if (length(sampled) == 1L && identical(latents, .sampler_location_latents(design))) {
    .sampler_location_scale_allocation(prior_list, sampled[[1L]])
  } else if (identical(latents, .sampler_location_cluster_latents(prior_list))) {
    .sampler_location_scale_clustered(prior_list)
  }
  if (is.null(layout)) return(NULL)
  list(name = "RoBMA_location_scale", tag = "robma_location_scale",
       arguments = list(layout[["a_source"]], layout[["a"]], layout[["sigma"]], layout[["weights"]],
                        layout[["share"]], layout[["indicators"]], layout[["levels"]], bias,
                        c(layout[["L"]], layout[["s"]], length(latents))))
}


# The allocation layout of the sampled term: one allocated heterogeneity
# component whose total SD is a simple prior, the term's share of it from the
# component's Dirichlet split (gamma sources) or the whole component (L = 1),
# and the component's inclusion indicator.
.sampler_location_scale_allocation <- function(prior_list, term) {

  component <- term[["component_label"]]
  sigma     <- term[["sd_parameter_names"]]
  if (!is.character(component) || length(component) != 1L || !is.character(sigma) || length(sigma) != 1L) {
    return(NULL)
  }
  has <- function(name, value) vapply(prior_list, function(prior) {
    identical(unname(attr(prior, name, exact = TRUE)), value)
  }, logical(1))
  # the allocated total SD of exactly this component (its terms are named by
  # component; a component with one term names the term itself)
  total <- names(prior_list)[vapply(prior_list, function(prior) {
    isTRUE(attr(prior, "random_allocation_sd", exact = TRUE)) &&
      identical(names(attr(prior, "random_allocation_terms", exact = TRUE)), component)
  }, logical(1))]
  if (length(total) != 1L || BayesTools::is.prior.mixture(prior_list[[total]])) return(NULL)
  allocation <- attr(prior_list[[total]], "random_allocation", exact = TRUE)
  a_source <- .sampler_location_scale_source(prior_list, total)
  if (is.null(a_source) || !identical(a_source, total)) return(NULL)

  split <- names(prior_list)[vapply(prior_list, function(prior) {
    parent <- attr(prior, "random_allocation_parent", exact = TRUE)
    identical(parent[["allocation"]], allocation) && identical(parent[["component"]], component)
  }, logical(1))]
  if (length(split) > 1L) return(NULL)
  if (length(split) == 1L) {
    terms   <- unname(attr(prior_list[[split]], "random_allocation_terms", exact = TRUE))
    sources <- BayesTools::JAGS_sampler_sources(prior_list, split)
    s       <- match(term[["block_name"]], terms)
    if (!inherits(prior_list[[split]], "prior.simplex") || is.na(s) || nrow(sources) != length(terms) ||
        !all(sources[["status"]] == "eligible")) return(NULL)
    weights <- paste0(split, "[", seq_along(terms), "]")
    L <- length(terms)
  } else {
    if (!identical(unname(attr(prior_list[[total]], "random_allocation_terms", exact = TRUE)), term[["block_name"]])) {
      return(NULL)
    }
    weights <- 0
    s <- 1L
    L <- 1L
  }

  inclusion <- names(prior_list)[has("random_allocation", allocation) & has("random_component", component)]
  indicator <- if (length(inclusion) == 1L) attr(prior_list[[inclusion]], "random_allocation_indicator", exact = TRUE)
  if (length(inclusion) > 1L || (length(inclusion) == 1L && (!is.character(indicator) || length(indicator) != 1L))) {
    return(NULL)
  }
  list(a_source = total, a = total, sigma = sigma, weights = weights, share = 0,
       indicators = if (is.null(indicator)) 0 else indicator, levels = if (is.null(indicator)) 0 else 1,
       L = L, s = s)
}


# The clustered layout: tau the total SD and rho the between-cluster share of
# its variance (tau_between = tau sqrt(rho) scales the cluster effects), each
# with one stochastic source behind its gate (or itself without a mixture).
.sampler_location_scale_clustered <- function(prior_list) {

  if (is.null(prior_list[["tau"]]) || is.null(prior_list[["rho"]])) return(NULL)
  gates <- character()
  levels <- numeric()
  source <- function(name) {
    prior <- prior_list[[name]]
    if (!BayesTools::is.prior.mixture(prior)) return(.sampler_location_scale_source(prior_list, name))
    mixture    <- BayesTools::JAGS_mixture_sources(prior_list, name)
    stochastic <- which(vapply(mixture[["branches"]], function(branch) isTRUE(branch[["stochastic"]]), logical(1)))
    if (length(stochastic) != 1L || !identical(mixture[["branches"]][[stochastic]][["family"]], "scalar") ||
        !is.character(mixture[["gate"]])) return(NULL)
    gates  <<- c(gates, mixture[["gate"]])
    levels <<- c(levels, stochastic)
    mixture[["branches"]][[stochastic]][["source"]]
  }
  tau <- source("tau")
  rho <- source("rho")
  if (is.null(tau) || is.null(rho)) return(NULL)
  list(a_source = tau, a = "tau", sigma = "tau_between", weights = "rho", share = rho,
       indicators = if (length(gates)) gates else 0, levels = if (length(levels)) levels else 0, L = 2L, s = 1L)
}


# The one eligible stochastic source of a simple prior, or NULL.
.sampler_location_scale_source <- function(prior_list, name) {

  sources <- tryCatch(BayesTools::JAGS_sampler_sources(prior_list, name), error = function(e) NULL)
  if (is.null(sources) || nrow(sources) != 1L || !identical(sources[["status"]], "eligible")) return(NULL)
  sources[["node"]]
}


# The scale owners' counters and factory decisions in this process (read and
# cleared in one call).
.location_scale_report_read <- function() {

  .Call("RoBMA_location_scale_report", TRUE, PACKAGE = "RoBMA")
}


# The scale ownership of a fit's retained runs: the distinct factory decisions
# (layout, L, K, status, reason, invariance check, embedded samplers) and the
# counters summed over runs and processes (the `_max` entries are maxima);
# NULL when no report covers it.
.location_scale_record <- function(values) {

  reports <- list()
  for (run in values) {
    for (process in run) {
      if (!is.null(process[["location_scale"]])) reports[[length(reports) + 1L]] <- process[["location_scale"]]
    }
  }
  if (length(reports) == 0L) return(NULL)
  decisions <- do.call(rbind, lapply(reports, function(report) as.data.frame(report[["decisions"]])))
  counts    <- do.call(rbind, lapply(reports, `[[`, "counters"))
  counters  <- colSums(counts)
  maxima    <- endsWith(colnames(counts), "_max")
  counters[maxima] <- apply(counts[, maxima, drop = FALSE], 2L, max)
  list(decisions = unique(decisions), counters = counters)
}
