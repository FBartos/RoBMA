# ============================================================================ #
# sampler-joint-location.R
# ============================================================================ #
#
# The location block sampler (RoBMA::LocationBlock,
# src/samplers/LocationBlock.cc) of model-averaged `.mv` fits owns the location
# gates, their normal slabs, the ungated normal location coefficients and the
# standardized latents of the sampled random-effect term. R declares the
# candidates from the bound priors and formula design; the native factory
# verifies the compiled graph and declines (the ordinary samplers keep the
# nodes) when it does not qualify.
#
# ============================================================================ #


# The automatic policy routes multilevel and random-effect fits of a Gaussian
# outcome (the `.mv` interfaces and clustered fits) to the location block; on
# request every fit of a Gaussian outcome qualifies. Independent fits reach it
# under the automatic policy only when no other joint sampler takes them
# (.sampler_joint_declared).
.sampler_joint_location_data <- function(data, policy) {

  if (.data_outcome_type(data) != "norm") return(FALSE)
  random <- .is_data_multilevel(data) || .is_data_random(data)
  return(random || policy == "location")
}


# The markers of an eligible fit or NULL. `probabilities` are the requested
# exact probability readouts: a location gate among them keeps the
# enumerating sampler, so the whole location block is declined. With
# `moderator`, a gated term other than the intercept is required.
.sampler_joint_location <- function(prior_list, formula_design, data, probabilities = character(),
                                    moderator = FALSE) {

  design <- formula_design[["mu"]]
  if (is.null(design)) return(NULL)
  map    <- design[["name_map"]]
  fixed  <- map[map[["kind"]] == "fixed" & map[["role"]] == "coefficient", , drop = FALSE]
  groups <- lapply(fixed[["jags_name"]], .sampler_location_group, prior_list = prior_list)
  groups <- groups[!vapply(groups, is.null, logical(1))]
  gates  <- unlist(lapply(groups, `[[`, "gate"))
  if (length(gates) == 0L) return(NULL)
  if (moderator) {
    intercept <- fixed[["jags_name"]][fixed[["term"]] == "intercept"]
    gated <- vapply(groups, function(group) !is.null(group[["gate"]]) && !group[["physical"]] %in% intercept, logical(1))
    if (!any(gated)) return(NULL)
  }
  if (any(gates %in% sub("\\[1\\]$", "", probabilities))) return(NULL)

  latents <- .sampler_location_latents(design)
  if (length(latents) == 0L) latents <- .sampler_location_cluster_latents(prior_list)
  anchors <- vapply(groups, function(group) {
    if (is.null(group[["gate"]])) group[["slab"]] else group[["gate"]]
  }, character(1))
  markers <- c(
    lapply(seq_along(groups), function(i) {
      group <- groups[[i]]
      list(name = paste0("RoBMA_location_group_", i), tag = "robma_location_group",
           arguments = list(anchors[[i]], group[["slab"]], group[["label"]], group[["physical"]], i))
    }),
    lapply(seq_along(latents), function(i) {
      list(name = paste0("RoBMA_location_latent_", i), tag = "robma_location_latent",
           arguments = list(latents[[i]], i))
    }),
    list(list(name = "RoBMA_location_setup", tag = "robma_location_setup",
              arguments = list(anchors[[1L]], length(groups), length(latents)))),
    Filter(Negate(is.null), list(.sampler_location_scale(prior_list, design, latents)))
  )
  intercept <- fixed[["jags_name"]][fixed[["term"]] == "intercept"]
  slot      <- which(vapply(groups, `[[`, character(1), "physical") %in% intercept)
  bias      <- .sampler_location_bias(prior_list, slot, probabilities)
  if (is.null(bias)) return(list(markers = markers))
  return(list(markers = c(markers, list(bias[["marker"]])), data = bias[["data"]]))
}


# The bias label of the location block (src/samplers/LocationBlockBias.h): the
# bias indicator, each label's family (0 none, 1 step selection, 2 PET,
# 3 PEESE) and the PET and PEESE slope sources, declared when the fit has a
# publication-bias mixture of these families with positive prior weights, a
# location group for the intercept and no exact probability readout of the
# indicator; NULL otherwise. The native factory verifies the compiled graph and
# keeps the ordinary samplers of the bias nodes when it declines.
.sampler_location_bias <- function(prior_list, intercept_slot, probabilities) {

  if (length(intercept_slot) != 1L || !inherits(prior_list[["bias"]], "prior.bias_mixture") ||
      length(prior_list[["bias"]]) < 2L) return(NULL)
  mixture  <- BayesTools::JAGS_mixture_sources(prior_list, "bias")
  if (mixture[["gate"]] %in% sub("\\[1\\]$", "", probabilities)) return(NULL)
  branches <- mixture[["branches"]]
  codes    <- c(none = 0, selection = 1, PET = 2, PEESE = 3)[vapply(branches, `[[`, character(1), "family")]
  if (anyNA(codes) || !any(codes == 1) || sum(codes == 2) > 1L || sum(codes == 3) > 1L ||
      !all(is.finite(mixture[["probabilities"]]) & mixture[["probabilities"]] > 0)) return(NULL)
  zero   <- "RoBMA_location_zero"
  source <- function(code) {
    branch <- branches[codes == code]
    if (length(branch) == 0L) return(zero)
    branch <- branch[[1L]]
    if (!isTRUE(branch[["stochastic"]]) || is.na(branch[["source"]])) return(NA_character_)
    branch[["source"]]
  }
  sources <- c(source(2), source(3))
  if (anyNA(sources)) return(NULL)
  data <- list(RoBMA_location_bias_codes = unname(codes))
  if (zero %in% sources) data[[zero]] <- 0
  list(
    marker = list(name = "RoBMA_location_bias", tag = "robma_location_bias",
                  arguments = list(mixture[["gate"]], "RoBMA_location_bias_codes", sources[[1L]], sources[[2L]],
                                   intercept_slot)),
    data   = data
  )
}


# One location term: a two-label mixture of an exact zero and a normal slab
# (gated), or a normal prior (ungated scalar); NULL for any other prior.
.sampler_location_group <- function(name, prior_list) {

  prior <- prior_list[[name]]
  if (is.null(prior) || !is.null(attr(prior, "multiply_by", exact = TRUE))) return(NULL)
  if (!BayesTools::is.prior.mixture(prior)) {
    if (!.sampler_prior_scalar(prior) || !.sampler_location_normal(prior)) return(NULL)
    return(list(gate = NULL, slab = name, label = 0L, physical = name))
  }
  if (BayesTools::is.prior.spike_and_slab(prior) || length(prior) != 2L) return(NULL)

  mixture  <- BayesTools::JAGS_mixture_sources(prior_list, name)
  branches <- mixture[["branches"]]
  family   <- vapply(branches, `[[`, character(1), "family")
  spike    <- which(family == "point")
  slab     <- which(family %in% c("scalar", "vector", "scalar_array"))
  if (length(spike) != 1L || length(slab) != 1L ||
      !all(is.finite(mixture[["probabilities"]]) & mixture[["probabilities"]] > 0)) return(NULL)
  location <- branches[[spike]][["location"]]
  if (!is.numeric(location) || !length(location) || !all(location == 0) ||
      !isTRUE(branches[[slab]][["stochastic"]]) || !.sampler_location_normal(branches[[slab]][["prior"]])) {
    return(NULL)
  }
  return(list(gate = mixture[["gate"]], slab = branches[[slab]][["source"]], label = slab, physical = name))
}


# A normal law with finite fixed parameters and no truncation: a scalar
# normal, or the normal or multivariate normal contrasts of a factor.
.sampler_location_normal <- function(prior) {

  if (!is.null(attr(prior, "multiply_by", exact = TRUE)) || BayesTools::is.prior.R2(prior) ||
      BayesTools::is.prior.regularized(prior) || BayesTools::is.prior.regularized_term(prior) ||
      BayesTools::is.prior.ordered(prior) || !isTRUE(prior[["distribution"]] %in% c("normal", "mnormal"))) {
    return(FALSE)
  }
  parameters <- prior[["parameters"]]
  truncation <- prior[["truncation"]]
  return(is.numeric(parameters[["mean"]]) && length(parameters[["mean"]]) == 1L &&
           is.finite(parameters[["mean"]]) && is.numeric(parameters[["sd"]]) &&
           length(parameters[["sd"]]) == 1L && isTRUE(parameters[["sd"]] > 0) &&
           is.finite(1 / parameters[["sd"]]^2) &&
           isTRUE(truncation[["lower"]] == -Inf) && isTRUE(truncation[["upper"]] == Inf))
}


# The standardized latents of the one sampled (non-centered) random-effect term
# with one column, or none: a child whose mean depends on two latents would
# make the native factory decline them anyway.
.sampler_location_latents <- function(design) {

  sampled <- Filter(function(term) identical(term[["compile_mode"]], "sampled"), design[["random_effects"]])
  if (length(sampled) != 1L) return(character())
  term <- sampled[[1L]]
  if (!identical(term[["parameterization_resolved"]], "noncentered") || !isTRUE(term[["n_columns"]] == 1L) ||
      !isTRUE(term[["n_groups"]] >= 1L) || !is.character(term[["parameter_stem"]])) return(character())
  return(paste0(term[["parameter_stem"]], "_xRE_Zx[", seq_len(term[["n_groups"]]), ",1]"))
}


# The standardized cluster effects gamma[k] ~ N(0, 1) of a clustered RoBMA() or
# BMA() fit (scaled by the between-cluster SD in each estimate's mean), one
# latent per cluster; none when the fit has no sampled cluster effects. The
# native factory verifies their laws and that each mean depends on one of them.
.sampler_location_cluster_latents <- function(prior_list) {

  if (is.null(prior_list[["gamma"]])) return(character())
  sources <- BayesTools::JAGS_sampler_sources(prior_list, "gamma")
  if (nrow(sources) == 0L || !all(sources[["status"]] == "eligible")) return(character())
  return(sources[["node"]])
}


# The per-fit counters of the location block samplers in this process
# (sweeps, the exact correction, per component type the proposals, stage
# passes and acceptances, child evaluations), read and cleared in one call.
.location_block_report_read <- function() {

  .Call("RoBMA_location_block_report", TRUE, PACKAGE = "RoBMA")
}


# The evidence record of the location gates in this process (one row per
# declaration slot and chain, after adaptation), read and cleared in one call.
.location_block_evidence_read <- function() {

  .Call("RoBMA_location_block_evidence", TRUE, PACKAGE = "RoBMA")
}


# The evidence record of the location block's bias label in this process (per
# chain the label moves' counts and the swap's stage-1 log ratios by
# direction, and the sweeps ending in each label), read and cleared in one call.
.location_bias_evidence_read <- function() {

  .Call("RoBMA_location_bias_evidence", TRUE, PACKAGE = "RoBMA")
}


# The counters of a fit's retained runs, summed over runs and processes (the
# 'max_' and '_max' entries are maxima), or NULL when no report covers them.
.location_block_record <- function(values) {

  reports <- list()
  for (run in values) {
    for (process in run) {
      if (!is.null(process[["location"]])) reports[[length(reports) + 1L]] <- process[["location"]]
    }
  }
  if (length(reports) == 0L) return(NULL)
  counts <- do.call(rbind, reports)
  maxima <- startsWith(colnames(counts), "max_") | endsWith(colnames(counts), "_max")
  record <- colSums(counts)
  record[maxima] <- apply(counts[, maxima, drop = FALSE], 2L, max)
  return(record)
}


# The runtime report of a fit (BayesTools::JAGS_fit(runtime_report)): the
# corrected selection sampler's counters for selection fits and the location
# block's counters when the fit declares it; NULL when neither applies.
.sampler_report_callback <- function(selection, location) {

  if (!selection && !location) return(NULL)
  callback <- function(context) {
    report <- if (selection) utils::getFromNamespace(".selection_sampler_report_read", "RoBMA")() else list()
    if (location) {
      report[["location"]] <- utils::getFromNamespace(".location_block_report_read", "RoBMA")()
      report[["location_evidence"]] <- utils::getFromNamespace(".location_block_evidence_read", "RoBMA")()
      report[["location_bias"]] <- utils::getFromNamespace(".location_bias_evidence_read", "RoBMA")()
      report[["location_scale"]] <- utils::getFromNamespace(".location_scale_report_read", "RoBMA")()
    }
    report
  }
  environment(callback) <- list2env(list(selection = selection, location = location), parent = baseenv())
  return(callback)
}
