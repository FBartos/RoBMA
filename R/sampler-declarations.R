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

  if (grepl("Corrected slice", message, fixed = TRUE) ||
      grepl("Surrogate slice transition failed", message, fixed = TRUE)) {
    return(paste0(
      message, "
The error was raised by the corrected selection sampler that RoBMA uses by default. ",
      "Fit the model again after 'RoBMA.options(selection.sampler = \"default\")' to sample these parameters ",
      "with the ordinary samplers, and please report the model to the RoBMA developers."
    ))
  }

  option <- if (grepl("Factor slice sampler:", message, fixed = TRUE)) {
    "sampler.blocks"
  } else if (grepl("gated node changes the likelihood", message, fixed = TRUE)) {
    "sampler.gates"
  } else {
    return(message)
  }

  return(paste0(
    message, "\nThe error was raised by a sampler of the BayesTools JAGS module that RoBMA uses by default. ",
    "Fit the model again after 'RoBMA.options(", option, " = FALSE)' to sample these parameters with the ordinary JAGS ",
    "samplers, and please report the model to the RoBMA developers."
  ))
}


# Gates: a bias-mixture prior, with or without a weight function branch.
.sampler_gates_declared <- function(fit_priors) {

  return(isTRUE(RoBMA.get_option("sampler.gates")) &&
           inherits(fit_priors[["bias"]], "prior.bias_mixture"))
}


# Blocks: NULL unless a selection likelihood is present; otherwise a function
# of the combined prior list (formula priors under their final names and with
# their random-effect attributes), which BayesTools::JAGS_fit() calls once.
# Model-averaged fits declare their allocation block unless the corrected
# selection sampler can take the allocation parameters (.sampler_ensemble_block()).
.sampler_blocks_declared <- function(data, priors) {

  if (!isTRUE(RoBMA.get_option("sampler.blocks")) ||
      !.is_priors_weightfunction(priors) || .selection_all_sources_conditioned(data)) {
    return(NULL)
  }

  scale_parameters <- vapply(
    .data_scale_component_specs(data), function(spec) spec[["parameter"]], character(1)
  )

  return(.sampler_block_function(scale_parameters, ensemble_block = .sampler_ensemble_block(data)))
}


# Whether a model-averaged fit declares its allocation block: yes unless the
# corrected selection sampler is on ('selection.sampler = "coarse_corrected"')
# and can take the allocation parameters, which needs a dependent selection
# block (dense, rank-one or factor; singletons never engage it) and, unless one
# of them is dense, a selection branch with the product rule (cluster and
# factor blocks under best-weight rules only have the exact density as their
# surrogate, so every corrected sampler would fall back).
.sampler_ensemble_block <- function(data) {

  if (!identical(RoBMA.get_option("selection.sampler"), "coarse_corrected") ||
      !.is_data_joint_selection(data)) {
    return(TRUE)
  }
  methods <- .data_selection_execution_plan(data)[["block_methods"]]
  if (!any(methods %in% c("dense", "rank_one", "factor"))) {
    return(TRUE)
  }
  model    <- .data_selection_model(data)
  branches <- model[["branches"]][model[["active_branches"]]]
  best     <- length(branches) > 0L && all(vapply(branches, function(branch) {
    identical(branch[["weight_rule"]], "best")
  }, logical(1L)))

  return(best && !any(methods == "dense"))
}


# A closure over the names of the scale formulas and the ensemble switch only
# (it keeps no data).
.sampler_block_function <- function(scale_parameters, ensemble_block) {

  force(scale_parameters)
  force(ensemble_block)
  return(function(prior_list) .sampler_block_members(prior_list, scale_parameters, ensemble_block))
}


# The block of a combined prior list: list(<prior names>) or NULL. With
# 'ensemble_block = FALSE' a model-averaged fit declares no block.
.sampler_block_members <- function(prior_list, scale_parameters = character(), ensemble_block = TRUE) {

  if (!.sampler_heterogeneity_estimated(prior_list, scale_parameters)) {
    return(NULL)
  }

  # RoBMA's unit-level latent parameters are never members, whatever their prior
  unit_latents <- c("theta", "gamma", "sampling_z", "pi", "phi")
  eligible     <- vapply(prior_list, .sampler_block_eligible, logical(1))
  members      <- setdiff(names(prior_list)[eligible], unit_latents)

  if (.sampler_is_ensemble(prior_list)) {
    if (!ensemble_block) {
      return(NULL)
    }
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

# Priors that may join a block: simple continuous scalar priors and treatment
# or independent factor priors without random-effect metadata, a single weight
# function, the plain allocation SD and the Dirichlet allocation weights.
# Vector priors (mean-difference and orthonormal factors, mnormal, mt), other
# random-effect priors (SDs, correlations) and the guarded priors are not.
.sampler_block_eligible <- function(prior) {

  if (.sampler_prior_guarded(prior)) {
    return(FALSE)
  }
  if (.sampler_is_allocation_sd(prior) || .sampler_is_allocation_weights(prior)) {
    return(TRUE)
  }
  if (.sampler_prior_random_metadata(prior)) {
    return(FALSE)
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
