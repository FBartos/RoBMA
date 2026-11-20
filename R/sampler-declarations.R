# ============================================================================ #
# sampler-declarations.R
# ============================================================================ #
#
# The sampler declarations RoBMA passes to BayesTools::JAGS_fit()
# ('sampler_gates', 'sampler_blocks'). They choose which sampler of the
# BayesTools JAGS module updates which node; the model, the priors, the
# likelihood and the posterior are unchanged. RoBMA only decides what to
# declare: BayesTools writes the markers, expands each prior to its nodes and
# declines blocks of fewer than two nodes.
#
# ============================================================================ #


# The declarations of a fresh fit, from the data, priors and fit priors of
# .fit(): gates for the branches of a publication-bias mixture, and one factor
# slice block for selection models with estimated heterogeneity.
.fit_sampler_declarations <- function(data, priors, fit_priors) {

  return(list(
    sampler_gates  = .sampler_gates_declared(fit_priors),
    sampler_blocks = .sampler_blocks_declared(data, priors)
  ))
}


# A fitting error raised by a BayesTools sampler or by RoBMA's corrected
# selection sampler keeps its message and gains a hint: the option that returns
# the parameters to the stock JAGS samplers, and a request to report the model
# (the density could not be computed, or a declaration is wrong).
.sampler_error_hint <- function(message) {

  # Errors of either block sampler (RoBMA's corrected block sampler and
  # BayesTools' factor slice sampler, which takes the same block under "slice"
  # and stops the same way): only undeclaring the blocks helps. A member on a
  # bound of its support is avoided by a start strictly inside it.
  bound <- if (grepl("current value lies on a bound of its support", message, fixed = TRUE)) {
    "A start strictly inside the support of the parameter avoids this error; otherwise fit"
  } else {
    "Fit"
  }
  if (grepl("Corrected slice (block)", message, fixed = TRUE) ||
      grepl("the support bounds of this node are too wide", message, fixed = TRUE)) {
    return(paste0(
      message, "\nThe error was raised by the sampler of a parameter block that RoBMA uses by default. ",
      bound, " the model again after 'RoBMA.options(sampler.blocks = FALSE)' to sample these parameters with ",
      "the ordinary samplers, and please report the model to the RoBMA developers."
    ))
  }
  if (grepl("Factor slice sampler:", message, fixed = TRUE)) {
    return(paste0(
      message, "\nThe error was raised by a sampler of the BayesTools JAGS module that RoBMA uses by default. ",
      bound, " the model again after 'RoBMA.options(sampler.blocks = FALSE)' to sample these parameters with ",
      "the ordinary JAGS samplers, and please report the model to the RoBMA developers."
    ))
  }

  if (grepl("Corrected slice", message, fixed = TRUE) ||
      grepl("Surrogate slice transition failed", message, fixed = TRUE)) {
    return(paste0(
      message, "
The error was raised by the corrected selection sampler that RoBMA uses by default. ",
      "Fit the model again after 'RoBMA.options(selection.sampler = \"slice\")' to sample these parameters ",
      "with the ordinary samplers, and please report the model to the RoBMA developers."
    ))
  }

  if (!grepl("gated node changes the likelihood", message, fixed = TRUE)) {
    return(message)
  }

  return(paste0(
    message, "\nThe error was raised by a sampler of the BayesTools JAGS module that RoBMA uses by default. ",
    "Fit the model again after 'RoBMA.options(sampler.gates = FALSE)' to sample these parameters with the ordinary JAGS ",
    "samplers, and please report the model to the RoBMA developers."
  ))
}


# Whether the samplers of the JAGS modules may be declared at all:
# 'selection.sampler = "jags.default"' keeps only the samplers JAGS assigns,
# whatever 'sampler.gates' and 'sampler.blocks' say.
.sampler_module_samplers <- function() {

  return(!isTRUE(RoBMA.get_option("selection.sampler") %in% "jags.default"))
}


# Gates: a bias-mixture prior, with or without a weight function branch.
.sampler_gates_declared <- function(fit_priors) {

  return(.sampler_module_samplers() && isTRUE(RoBMA.get_option("sampler.gates")) &&
           inherits(fit_priors[["bias"]], "prior.bias_mixture"))
}


# Blocks: NULL unless a selection likelihood is present; otherwise a function
# of the combined prior list (formula priors under their final names and with
# their random-effect attributes), which BayesTools::JAGS_fit() calls once.
# Under 'selection.sampler = "slice.coarse"' RoBMA's corrected block sampler
# takes a declared block with a coarse-capable selection child; BayesTools'
# factor slice sampler takes the others.
.sampler_blocks_declared <- function(data, priors) {

  if (!.sampler_module_samplers() || !isTRUE(RoBMA.get_option("sampler.blocks")) ||
      !.is_priors_weightfunction(priors) || .selection_all_sources_conditioned(data)) {
    return(NULL)
  }

  scale_parameters <- vapply(
    .data_scale_component_specs(data), function(spec) spec[["parameter"]], character(1)
  )

  return(.sampler_block_function(scale_parameters))
}


# A closure over the names of the scale formulas only (it keeps no data).
.sampler_block_function <- function(scale_parameters) {

  force(scale_parameters)
  return(function(prior_list) .sampler_block_members(prior_list, scale_parameters))
}


# The block of a combined prior list: list(<prior names>) or NULL.
.sampler_block_members <- function(prior_list, scale_parameters = character()) {

  if (!.sampler_heterogeneity_estimated(prior_list, scale_parameters)) {
    return(NULL)
  }

  # RoBMA's unit-level latent parameters are never members, whatever their prior
  unit_latents <- c("theta", "gamma", "sampling_z", "pi", "phi")
  eligible     <- vapply(prior_list, .sampler_block_eligible, logical(1))
  members      <- setdiff(names(prior_list)[eligible], unit_latents)

  if (.sampler_is_ensemble(prior_list)) {
    # Model-averaged fits: the plain allocation SD and the Dirichlet allocation
    # weights of a random-effect variance allocation only.
    members <- members[vapply(prior_list[members], function(prior) {
      .sampler_is_allocation_sd(prior) || .sampler_is_allocation_weights(prior)
    }, logical(1))]
  }

  if (length(members) < 2L) {
    return(NULL)
  }

  return(list(members))
}


# Priors that are never declared, whatever their role: mixtures (bias
# mixtures included), spike-and-slab, ordered (totals included), point and
# discrete priors, and inclusion priors of random-effect components.
.sampler_prior_guarded <- function(prior) {

  return(BayesTools::is.prior.mixture(prior) ||
           BayesTools::is.prior.spike_and_slab(prior) ||
           BayesTools::is.prior.ordered(prior) ||
           BayesTools::is.prior.point(prior) ||
           BayesTools::is.prior.discrete(prior) ||
           !is.null(attr(prior, "random_allocation_inclusion", exact = TRUE)))
}

.sampler_prior_scalar <- function(prior) {

  return(BayesTools::is.prior.simple(prior) &&
           !BayesTools::is.prior.factor(prior) &&
           !BayesTools::is.prior.vector(prior))
}

# Random-effect metadata BayesTools attaches to formula priors.
.sampler_prior_random_metadata <- function(prior) {

  attributes <- c("random_sd", "random_correlation", "random_allocation_sd",
                  "random_allocation", "random_allocation_inclusion", "random_summary")
  return(any(vapply(attributes, function(attribute) {
    !is.null(attr(prior, attribute, exact = TRUE))
  }, logical(1))))
}

.sampler_is_allocation_sd <- function(prior) {

  return(!.sampler_prior_guarded(prior) && .sampler_prior_scalar(prior) &&
           isTRUE(attr(prior, "random_allocation_sd", exact = TRUE)))
}

.sampler_is_allocation_weights <- function(prior) {

  return(!.sampler_prior_guarded(prior) && BayesTools::is.prior.vector(prior) &&
           !is.null(attr(prior, "random_allocation", exact = TRUE)) &&
           !isTRUE(attr(prior, "random_allocation_sd", exact = TRUE)) &&
           identical(prior[["distribution"]], "dirichlet"))
}

# A direct random-effect SD or a scalar structured correlation (e.g. of cs() or
# hcs()) with a simple continuous scalar prior, outside any allocation.
.sampler_is_random_scalar <- function(prior) {

  return(!.sampler_prior_guarded(prior) && .sampler_prior_scalar(prior) &&
           (isTRUE(attr(prior, "random_sd", exact = TRUE)) ||
              isTRUE(attr(prior, "random_correlation", exact = TRUE))) &&
           is.null(attr(prior, "random_allocation", exact = TRUE)) &&
           is.null(attr(prior, "random_summary", exact = TRUE)))
}

# Priors that may join a block: simple continuous scalar priors and treatment
# or independent factor priors without random-effect metadata, a single weight
# function, the plain allocation SD and the Dirichlet allocation weights, and
# direct random-effect SDs and scalar structured correlations with simple
# continuous scalar priors. Vector priors (mean-difference and orthonormal
# factors, mnormal, mt), other random-effect priors (factor or vector SDs; LKJ
# correlations are not prior-list entries) and the guarded priors are not.
# The ensemble filter of .sampler_block_members() keeps only the allocation
# parameters of model-averaged fits.
.sampler_block_eligible <- function(prior) {

  if (.sampler_prior_guarded(prior)) {
    return(FALSE)
  }
  if (.sampler_is_allocation_sd(prior) || .sampler_is_allocation_weights(prior)) {
    return(TRUE)
  }
  if (.sampler_prior_random_metadata(prior)) {
    return(.sampler_is_random_scalar(prior))
  }

  return(.sampler_prior_scalar(prior) ||
           (BayesTools::is.prior.factor(prior) && !BayesTools::is.prior.vector(prior) &&
              (BayesTools::is.prior.treatment(prior) || BayesTools::is.prior.independent(prior))) ||
           BayesTools::is.prior.weightfunction(prior))
}


# A model-averaged fit: a mixture, a spike-and-slab or an inclusion prior.
.sampler_is_ensemble <- function(prior_list) {

  return(any(vapply(prior_list, function(prior) {
    BayesTools::is.prior.mixture(prior) || BayesTools::is.prior.spike_and_slab(prior) ||
      !is.null(attr(prior, "random_allocation_inclusion", exact = TRUE))
  }, logical(1))))
}


# A heterogeneity parameter is estimated: a non-point prior of 'tau', 'rho', a
# scale formula (location-scale models), a random-effect SD or an allocation SD.
.sampler_heterogeneity_estimated <- function(prior_list, scale_parameters = character()) {

  heterogeneity <- vapply(names(prior_list), function(name) {
    prior <- prior_list[[name]]
    name %in% c("tau", "rho") ||
      isTRUE(attr(prior, "parameter", exact = TRUE) %in% scale_parameters) ||
      isTRUE(attr(prior, "random_sd", exact = TRUE)) ||
      isTRUE(attr(prior, "random_allocation_sd", exact = TRUE))
  }, logical(1))
  point <- vapply(prior_list, BayesTools::is.prior.point, logical(1))

  return(any(heterogeneity & !point))
}
