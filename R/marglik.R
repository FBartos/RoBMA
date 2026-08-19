# ============================================================================ #
# marglik.R
# ============================================================================ #
#
# This file implements marginal likelihood computation for brma class objects
# via bridge sampling. The marginal likelihood is used for Bayesian model
# comparison via Bayes factors.
#
# The implementation reuses the evaluate and pdf helper functions from
# evaluate.R and pdf.R by converting single-sample parameters
# to 1-row matrices.
#
# ============================================================================ #


### RoBMA 4.0.0


# ---------------------------------------------------------------------------- #
# add_marglik generic and method
# ---------------------------------------------------------------------------- #

#' @export
add_marglik <- function(object, ...) UseMethod("add_marglik")


#' @title Add Marginal Likelihood to brma Objects
#'
#' @description Compute the marginal likelihood of a brma model using
#' bridge sampling and store the result in the object.
#'
#' @param object a brma model object.
#' @param parallel whether bridge-density evaluations run in parallel. The
#' default, \code{NULL}, inherits the setting used to fit \code{object}.
#' @param cores number of bridge-sampling worker processes. The default,
#' \code{NULL}, inherits the fitted core count and is capped by
#' \code{RoBMA.get_option("max_cores")}. It is used only when
#' \code{parallel = TRUE}.
#' @param repetitions number of independent bridge-sampling repetitions.
#' @param method bridge transformation; either \code{"normal"} or
#' \code{"warp3"}.
#' @param maxiter maximum number of bridge iterations per repetition.
#' @param silent whether bridge-sampling progress is suppressed.
#' @param ... reserved for future use.
#'
#' @details
#' The marginal likelihood is computed using the \code{bridgesampling} package
#' via the \code{BayesTools::JAGS_bridgesampling} wrapper. The result is stored
#' in the object and can be extracted using \code{\link{bridge_sampler.brma}}.
#'
#' Product-space model-averaging objects (\code{BMA.norm}, \code{BMA.glmm},
#' and \code{RoBMA}) do not expose a bridge-sampling marginal likelihood;
#' use predictive comparison methods such as \code{\link{loo.brma}} instead.
#' For a single model whose bridge parameters are all fixed by point priors,
#' the log marginal likelihood is evaluated exactly at those fixed values;
#' no bridge-sampling repetitions are required.
#' For \code{brma.mv()} known-\code{V} objects, bridge sampling evaluates the
#' joint likelihood corresponding to the fitted known-\code{V} backend, not the
#' conditional estimate-wise target used by LOO/WAIC diagnostics. Sampled
#' Gaussian location random-effect blocks are integrated exactly during bridge
#' evaluation as their draw-specific \eqn{ZGZ'} covariance. Their SD,
#' allocation, and correlation parameters remain bridge coordinates; only the
#' standardized latent effects are removed. This is a reparameterization of the
#' same marginal-likelihood target, not a likelihood approximation. Fitted
#' estimate-level marginalized blocks remain in the diagonal row variance.
#' For Gaussian models fitted with the specialized \code{cluster} interface,
#' the standardized cluster effects are likewise integrated exactly. The bridge
#' retains total heterogeneity and its allocation while evaluating the implied
#' diagonal-plus-cluster-rank-one covariance, including row-specific scale
#' regression and likelihood weights.
#' Selection likelihoods retain the fitted joint latent parameterization because
#' they are not Gaussian in the shared random effects.
#'
#' @return The brma object with the marginal likelihood result stored in
#' \code{object[["marglik"]]}.
#'
#' @seealso \code{\link{bridge_sampler.brma}}, \code{\link{logml.brma}},
#' \code{\link{bf.brma}}, \code{\link{post_prob.brma}}
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   fit <- add_marglik(fit)
#'
#'   bridge <- bridge_sampler(fit)
#'   print(bridge)
#'
#'   logml(fit)
#' }
#' }
#'
#' @aliases add_marglik
#' @export
add_marglik.brma <- function(object, parallel = NULL, cores = NULL,
                             repetitions = 1L,
                             method = c("normal", "warp3"),
                             maxiter = 10000L, silent = TRUE, ...) {

  .check_marglik_available(object, "add_marglik()")
  method <- match.arg(method)
  parallel_control <- .marglik_parallel_control(
    object   = object,
    parallel = parallel,
    cores    = cores
  )
  marglik <- .marglik(
    object      = object,
    cores       = parallel_control[["cores"]],
    repetitions = repetitions,
    method      = method,
    maxiter     = maxiter,
    silent      = silent
  )
  if (inherits(marglik, "error")) {
    stop(conditionMessage(marglik), call. = FALSE)
  }
  marglik <- .brma_mv_attach_marglik_target_metadata(marglik, object)
  object[["marglik"]] <- marglik
  return(object)
}


# ---------------------------------------------------------------------------- #
# .marglik internal function
# ---------------------------------------------------------------------------- #

# Compute marginal likelihood for brma objects using bridge sampling.
# This is an internal function called by add_marglik.brma.
#
# @param object a brma model object.
#
# @return A `BayesTools_marglik` object.
#
# @keywords internal
.marglik <- function(object, cores = 1L, repetitions = 1L,
                     method = c("normal", "warp3"), maxiter = 10000L,
                     silent = TRUE) {

  data   <- object[["data"]]
  priors <- object[["priors"]]
  fit    <- object[["fit"]]
  method <- match.arg(method)

  .check_marglik_available(object, ".marglik()")

  ### create model base
  fit_priors       <- .create_fit_priors(data = data, priors = priors)
  fit_data         <- .create_fit_data(data = data, priors = priors)
  fit_data         <- .marglik_add_selection_bridge_data(
    fit_data         = fit_data,
    priors           = priors,
    effect_direction = .data_effect_direction(data)
  )
  bridge_setup <- .marglik_fixed_zero_random_setup(
    object     = object,
    fit        = fit,
    fit_priors = fit_priors
  )
  fit        <- bridge_setup[["fit"]]
  fit_priors <- bridge_setup[["fit_priors"]]
  cluster_marginalization <- .marglik_cluster_effects_setup(
    data       = data,
    priors     = priors,
    fit_priors = fit_priors
  )
  fit_priors <- cluster_marginalization[["fit_priors"]]
  sampling_latent_setup <- .marglik_sampling_latent_setup(
    data       = data,
    priors     = priors,
    fit_priors = fit_priors
  )
  fit_priors <- sampling_latent_setup[["fit_priors"]]
  bridge_random_marginalization <- .marglik_bridge_random_marginalization(
    object            = object,
    fit               = fit,
    fixed_zero_random = bridge_setup[["fixed_zero_random"]],
    sampling_latent_marginalized = sampling_latent_setup[["marginalized"]]
  )
  if (sampling_latent_setup[["marginalized"]] &&
      is.null(bridge_random_marginalization[["dependency_blocks"]])) {
    bridge_random_marginalization[["dependency_blocks"]] <-
      .marglik_random_dependency_blocks(
        model_data                  = data,
        formula_design              = NULL,
        blocks                     = character(),
        sampling_latent_marginalized = TRUE
      )
  }
  fit_data <- .marglik_add_random_covariance_bridge_data(
    fit_data          = fit_data,
    model_data        = data,
    marginalizing     = length(bridge_random_marginalization[["blocks"]]) > 0L ||
      sampling_latent_setup[["marginalized"]],
    sampling_latent_marginalized = sampling_latent_setup[["marginalized"]],
    dependency_blocks = bridge_random_marginalization[["dependency_blocks"]]
  )
  bridge_sd_source_spec <- .marglik_bridge_sd_source_spec(
    add_parameters = .marglik_formula_source_parameters(data),
    fit            = fit,
    K              = nrow(data[["outcome"]])
  )
  known_v_setup <- .marglik_known_v_setup(data)
  marginalized_variance_plan <- if (bridge_setup[["fixed_zero_random"]]) {
    NULL
  } else {
    .marglik_marginalized_variance_plan(data)
  }
  bridge_context_mode <- if (length(bridge_random_marginalization[["blocks"]]) > 0L) {
    "marginal"
  } else if (!bridge_setup[["fixed_zero_random"]] &&
             .marglik_needs_bridge_context(data)) {
    "nodes"
  } else {
    FALSE
  }
  covariance_plan_cache <- new.env(parent = emptyenv())

  ### compute marginal likelihood
  marglik <- BayesTools::JAGS_bridgesampling(
    fit                = fit,
    log_posterior      = .log_posterior,
    data               = fit_data,
    prior_list         = fit_priors,
    formula_random_effects_marginalize_list = .optional_jags_list(
      bridge_random_marginalization[["request"]]
    ),
    add_parameters                      = .optional_jags_character(bridge_sd_source_spec[["parameters"]]),
    add_bounds                          = bridge_sd_source_spec[["bounds"]],
    bridge_context                      = bridge_context_mode,
    bridge_context_node_names           = .marglik_variance_plan_node_names(
      marginalized_variance_plan
    ),
    repetitions                         = repetitions,
    method                              = method,
    maxiter                             = maxiter,
    silent                              = silent,
    cores                               = cores,
    packages                            = .marglik_bridge_packages(fit, cores),
    # additional arguments passed to .log_posterior via ...
    is_mods                  = .is_data_mods(data),
    is_scale                 = .is_data_scale(data),
    is_random                = .is_data_random(data),
    is_multilevel            = .is_data_multilevel(data),
    is_weights               = .is_data_weights(data),
    is_known_v               = .is_data_known_v(data),
    model_data               = data,
    known_V                  = known_v_setup[["known_V"]],
    known_v_backend          = known_v_setup[["backend"]],
    marginalized_variance_plan = marginalized_variance_plan,
    is_PET                   = .is_priors_PET(priors),
    is_PEESE                 = .is_priors_PEESE(priors),
    is_weightfunction        = .is_priors_weightfunction(priors),
    fixed_tau                = .fixed_tau_prior_value(priors),
    fixed_rho                = bridge_setup[["fixed_rho"]],
    fixed_zero_random        = bridge_setup[["fixed_zero_random"]],
    cluster_effects_marginalized = cluster_marginalization[["marginalized"]],
    sampling_latent_marginalized = sampling_latent_setup[["marginalized"]],
    covariance_plan_cache    = covariance_plan_cache,
    effect_direction         = .data_effect_direction(data),
    outcome_type             = .data_outcome_type(data)
  )

  marglik[["diagnostics"]][["random_effect_marginalization"]] <-
    bridge_random_marginalization[["diagnostics"]]
  marglik[["diagnostics"]][["sampling_latent_marginalization"]] <-
    sampling_latent_setup[["diagnostics"]]
  marglik[["diagnostics"]][["cluster_effects_marginalization"]] <-
    cluster_marginalization[["diagnostics"]]

  return(marglik)
}


.marglik_known_v_setup <- function(data) {

  if (!.is_data_known_v(data)) {
    return(list(known_V = NULL, backend = NULL))
  }

  known_V <- .data_known_v_data(data)
  list(
    known_V = known_V,
    backend = .known_v_effective_backend(known_V)
  )
}


.marglik_marginalized_variance_plan <- function(data) {

  if (!.is_data_known_v(data) || !.is_data_random(data)) {
    return(NULL)
  }

  terms <- .data_marginalized_random_effects(data)
  K     <- nrow(data[["outcome"]])
  plans <- vector("list", length(terms))
  for (term_i in seq_along(terms)) {
    term      <- terms[[term_i]]
    parameter <- term[["sd_parameter_names"]]
    if (length(parameter) != 1L || is.na(parameter) || !nzchar(parameter)) {
      return(NULL)
    }
    multiplier <- .marginalized_random_effect_row_multiplier(term, K = K)
    plans[[term_i]] <- list(
      parameter  = parameter,
      multiplier = if (is.null(multiplier)) rep(1, K) else multiplier
    )
  }

  list(K = K, terms = plans)
}


.marglik_variance_plan_node_names <- function(plan) {

  if (is.null(plan)) {
    return(NULL)
  }
  unique(vapply(plan[["terms"]], `[[`, character(1), "parameter"))
}


.marglik_cluster_effects_setup <- function(data, priors, fit_priors) {

  eligible <- .is_data_multilevel(data) &&
    .data_outcome_type(data) == "norm" &&
    !.is_priors_weightfunction(priors) &&
    "gamma" %in% names(fit_priors)
  if (!eligible) {
    reason <- if (!.is_data_multilevel(data)) {
      "model has no cluster effect from the specialized multilevel interface"
    } else if (.data_outcome_type(data) != "norm") {
      "cluster effects enter a non-Gaussian likelihood"
    } else if (.is_priors_weightfunction(priors)) {
      "selection normalization is non-Gaussian in the cluster effect"
    } else {
      "cluster effect is already structurally absent"
    }
    return(list(
      fit_priors  = fit_priors,
      marginalized = FALSE,
      diagnostics = list(
        requested = FALSE,
        exact     = TRUE,
        reason    = reason
      )
    ))
  }

  list(
    fit_priors  = fit_priors[setdiff(names(fit_priors), "gamma")],
    marginalized = TRUE,
    diagnostics = list(
      requested = TRUE,
      included  = "gamma",
      exact     = TRUE,
      target    = paste(
        "standard-normal cluster effects integrated as diagonal plus",
        "cluster-block rank-one covariance"
      )
    )
  )
}


.marglik_parallel_control <- function(object, parallel = NULL, cores = NULL) {

  fit_control <- object[["fit_control"]]
  if (is.null(parallel)) {
    parallel <- fit_control[["parallel"]]
    if (is.null(parallel)) {
      parallel <- FALSE
    }
  }
  BayesTools::check_bool(parallel, "parallel", allow_NA = FALSE)

  if (is.null(cores)) {
    cores <- fit_control[["cores"]]
    if (is.null(cores)) {
      cores <- RoBMA.get_option("max_cores")
    }
  }
  BayesTools::check_int(cores, "cores", lower = 1, allow_NA = FALSE)

  max_cores <- as.integer(RoBMA.get_option("max_cores"))
  cores <- if (isTRUE(parallel)) {
    min(as.integer(cores), max_cores)
  } else {
    1L
  }

  list(parallel = isTRUE(parallel), cores = cores)
}


.marglik_bridge_packages <- function(fit, cores) {

  if (cores <= 1L) {
    return(NULL)
  }

  packages <- c("BayesTools", "RoBMA", attr(fit, "required_packages"))
  unique(packages[!is.na(packages) & nzchar(packages)])
}


# Remove bridge coordinates that are provably irrelevant because their SD
# ancestor is point-fixed at zero. The fitted object is copied before its
# formula-design replay metadata are changed.
.marglik_fixed_zero_random_setup <- function(object, fit, fit_priors) {

  fixed_tau <- .fixed_tau_prior_value(object[["priors"]])
  fixed_rho <- .fixed_rho_prior_value(object[["priors"]])

  if (.is_data_multilevel(object[["data"]]) &&
      !is.null(fixed_tau) && identical(fixed_tau, 0)) {
    fit_priors <- fit_priors[setdiff(names(fit_priors), c("gamma", "rho"))]
    if (is.null(fixed_rho)) {
      fixed_rho <- 0
    }
  }

  if (!.is_data_random(object[["data"]])) {
    return(list(
      fit               = fit,
      fit_priors        = fit_priors,
      fixed_rho         = fixed_rho,
      fixed_zero_random = FALSE
    ))
  }

  formula_design <- attr(fit, "formula_design", exact = TRUE)
  mu_design      <- formula_design[["mu"]]
  if (!.marglik_formula_random_effects_fixed_zero(
    design = mu_design,
    data   = object[["data"]]
  )) {
    return(list(
      fit               = fit,
      fit_priors        = fit_priors,
      fixed_rho         = fixed_rho,
      fixed_zero_random = FALSE
    ))
  }

  # This path is rare and exact: a fitted sampled random block whose complete
  # scale is point-fixed at zero. Rebuild only that fixed design so BayesTools
  # does not need irrelevant latent coordinates. Ordinary bridge setup consumes
  # the authoritative fitted design directly and never rebuilds formulas.
  fit_formula_args <- .create_jags_formula_args(
    data   = object[["data"]],
    priors = object[["priors"]]
  )
  fixed_formula <- mu_design[["formula"]]
  attr(fixed_formula, "random_terms")      <- NULL
  attr(fixed_formula, "random_components") <- NULL
  environment(fixed_formula) <- environment(
    fit_formula_args[["formula_list"]][["mu"]]
  )

  fixed_design <- BayesTools::JAGS_formula(
    formula       = fixed_formula,
    parameter     = "mu",
    data          = fit_formula_args[["formula_data_list"]][["mu"]],
    prior_list    = fit_formula_args[["formula_prior_list"]][["mu"]],
    formula_scale = fit_formula_args[["formula_scale_list"]][["mu"]]
  )[["formula_design"]]

  formula_design[["mu"]] <- fixed_design
  attr(fit, "formula_design") <- formula_design

  return(list(
    fit               = fit,
    fit_priors        = fit_priors,
    fixed_rho         = fixed_rho,
    fixed_zero_random = TRUE
  ))
}


.marglik_formula_random_effects_fixed_zero <- function(design, data) {

  if (is.null(design) || length(design[["random_effects"]]) == 0L) {
    return(FALSE)
  }

  prior_list <- design[["prior_list"]]
  K          <- nrow(data[["outcome"]])
  all(vapply(
    design[["random_effects"]],
    .marglik_random_effect_fixed_zero,
    logical(1),
    data       = data,
    prior_list = prior_list,
    K          = K
  ))
}


.marglik_random_effect_fixed_zero <- function(term, data, prior_list, K) {

  binding <- term[["sd_binding"]]
  if (is.null(binding)) {
    return(FALSE)
  }

  source <- binding[["source"]]
  if (!is.null(source)) {
    values <- .random_sd_source_fixed_values(
      source     = source,
      data       = data,
      prior_list = prior_list,
      K          = K
    )
    return(!is.null(values) && all(values == 0))
  }

  sources <- binding[["sources_by_column"]]
  if (length(sources) > 0L) {
    fixed_zero <- vapply(sources, function(column_source) {
      values <- .random_sd_source_fixed_values(
        source     = column_source,
        data       = data,
        prior_list = prior_list,
        K          = K
      )
      !is.null(values) && all(values == 0)
    }, logical(1))
    return(all(fixed_zero))
  }

  parameters <- term[["sd_parameter_names"]]
  parameters <- parameters[!is.na(parameters) & nzchar(parameters)]
  if (length(parameters) == 0L) {
    return(FALSE)
  }

  all(vapply(parameters, function(parameter) {
    values <- .prior_fixed_values(prior_list[[parameter]])
    !is.null(values) && all(values == 0)
  }, logical(1)))
}


.marglik_needs_bridge_context <- function(data) {

  .is_data_known_v(data) &&
    .is_data_random(data) &&
    .data_has_marginalized_random_effects(data)
}


.marglik_formula_source_parameters <- function(data) {

  if (!.is_data_random(data) || !.is_data_scale(data)) {
    return(character())
  }

  unique(.data_scale_formula_sources(data))
}


.marglik_sampling_latent_setup <- function(data, priors, fit_priors) {

  eligible <- .is_data_known_v_backend(data, "latent") &&
    .data_outcome_type(data) == "norm" &&
    !.is_priors_weightfunction(priors) &&
    .data_known_v_rank(data) > 0L
  if (!eligible) {
    return(list(
      fit_priors = fit_priors,
      marginalized = FALSE,
      diagnostics = list(
        requested = FALSE,
        included = character(),
        exact = TRUE,
        reason = "not an eligible Gaussian latent known-V likelihood"
      )
    ))
  }

  sampling_prior <- fit_priors[["sampling_z"]]
  valid_standard_normal <- inherits(sampling_prior, "prior.simple") &&
    identical(sampling_prior[["distribution"]], "normal") &&
    identical(as.numeric(sampling_prior[["parameters"]][["mean"]]), 0) &&
    identical(as.numeric(sampling_prior[["parameters"]][["sd"]]), 1) &&
    identical(as.numeric(sampling_prior[["truncation"]][["lower"]]), -Inf) &&
    identical(as.numeric(sampling_prior[["truncation"]][["upper"]]), Inf) &&
    identical(as.numeric(sampling_prior[["prior_weights"]]), 1)
  if (!isTRUE(valid_standard_normal)) {
    stop(
      "Known-V sampling latent marginalization requires its fitted normalized standard-normal prior.",
      call. = FALSE
    )
  }

  fit_priors[["sampling_z"]] <- NULL
  list(
    fit_priors = fit_priors,
    marginalized = TRUE,
    diagnostics = list(
      requested = TRUE,
      included = paste0("sampling_z[", seq_len(.data_known_v_rank(data)), "]"),
      exact = TRUE,
      target = "Gaussian sampling factors integrated into V"
    )
  )
}


.marglik_bridge_random_marginalization <- function(object, fit,
                                                     fixed_zero_random,
                                                     sampling_latent_marginalized = FALSE) {

  data   <- object[["data"]]
  priors <- object[["priors"]]
  empty <- function(reason) {
    list(
      blocks = list(),
      diagnostics = list(
        requested = FALSE,
        included = character(),
        skipped = data.frame(
          block_name = character(),
          reason = character(),
          stringsAsFactors = FALSE
        ),
        reason = reason
      )
    )
  }

  if (!.is_data_random(data) || !.is_data_known_v(data) ||
      .data_outcome_type(data) != "norm") {
    return(empty("not an eligible Gaussian known-V random-formula model"))
  }
  if (.is_priors_weightfunction(priors)) {
    return(empty("selection likelihood is not Gaussian in the shared random effects"))
  }
  if (isTRUE(fixed_zero_random)) {
    return(empty("all random-effect contributions are point-fixed at zero"))
  }

  formula_design <- attr(fit, "formula_design", exact = TRUE)
  mu_design      <- formula_design[["mu"]]
  if (is.null(mu_design)) {
    return(empty("fitted mu formula design is unavailable"))
  }
  all_terms      <- mu_design[["random_effects"]]
  sampled_blocks <- .formula_design_sampled_random_effect_blocks(mu_design)
  if (length(all_terms) == 0L || length(sampled_blocks) == 0L) {
    return(empty("no fitted sampled Gaussian random-effect blocks"))
  }

  block_names <- vapply(all_terms, .random_effect_term_block_name, character(1))
  compile_modes <- vapply(
    all_terms,
    .random_effect_term_compile_mode,
    character(1)
  )
  eligible <- compile_modes == "sampled"
  sampled_blocks <- block_names[eligible]
  skipped <- data.frame(
    block_name = block_names[compile_modes != "sampled"],
    reason = rep(
      "already marginalized by the fitted likelihood",
      sum(compile_modes != "sampled")
    ),
    stringsAsFactors = FALSE
  )
  if (length(sampled_blocks) == 0L) {
    out <- empty("no sampled random-effect block has an exact covariance contract")
    out[["diagnostics"]][["skipped"]] <- skipped
    return(out)
  }
  dependency_blocks <- .marglik_random_dependency_blocks(
    model_data    = data,
    formula_design = mu_design,
    blocks        = sampled_blocks,
    sampling_latent_marginalized = sampling_latent_marginalized
  )

  list(
    blocks = list(mu = sampled_blocks),
    request = list(
      mu = list(
        blocks = sampled_blocks,
        row_blocks = dependency_blocks,
        factor_state = TRUE
      )
    ),
    dependency_blocks = dependency_blocks,
    diagnostics = list(
      requested = TRUE,
      included = sampled_blocks,
      skipped = skipped,
      exact = TRUE,
      target = "Gaussian random effects integrated as ZGZ'"
    )
  )
}


.marglik_add_random_covariance_bridge_data <- function(
    fit_data, model_data, marginalizing,
    sampling_latent_marginalized = FALSE, dependency_blocks = NULL) {

  if (!isTRUE(marginalizing)) {
    return(fit_data)
  }
  if (!.is_data_known_v(model_data)) {
    stop("Random-covariance bridge data require known-V model metadata.",
         call. = FALSE)
  }

  known_V <- .data_known_v_data(model_data)
  backend <- .data_known_v_effective_backend(model_data)
  covariance <- if (backend == "latent" && !sampling_latent_marginalized) {
    diag(.known_v_residual_variance(known_V), nrow = .known_v_nrow(known_V))
  } else {
    .marglik_known_v_covariance_matrix(known_V)
  }
  fit_data[["marglik_sampling_covariance"]] <- covariance
  fit_data[["marglik_dependency_blocks"]] <- dependency_blocks

  fit_data
}


.marglik_random_dependency_blocks <- function(model_data, formula_design,
                                               blocks,
                                               sampling_latent_marginalized = FALSE) {

  known_V <- .data_known_v_data(model_data)
  K       <- .known_v_nrow(known_V)
  backend <- .data_known_v_effective_backend(model_data)
  adjacency <- if (backend == "latent" && !sampling_latent_marginalized) {
    diag(TRUE, nrow = K, ncol = K)
  } else {
    .marglik_known_v_covariance_matrix(known_V) != 0
  }
  terms <- if (is.null(formula_design)) NULL else
    formula_design[["random_effects"]]
  term_names <- vapply(terms, .random_effect_term_block_name, character(1))
  terms <- terms[term_names %in% blocks]

  for (term in terms) {
    group_map <- as.integer(term[["group_map"]])
    if (length(group_map) != K || anyNA(group_map) || any(group_map < 1L)) {
      return(list(seq_len(K)))
    }
    if (.random_effect_term_has_known_group_covariance(term)) {
      group_covariance <- term[["group_covariance"]]
      kernel <- group_covariance[["kernel"]]
      if (!is.matrix(kernel) || any(group_map > nrow(kernel))) {
        return(list(seq_len(K)))
      }
      adjacency <- adjacency | kernel[group_map, group_map, drop = FALSE] != 0
    } else {
      adjacency <- adjacency | outer(group_map, group_map, "==")
    }
  }
  diag(adjacency) <- TRUE

  .known_v_block_indices(adjacency * 1)
}


.marglik_known_v_covariance_matrix <- function(known_V) {

  K <- .known_v_nrow(known_V)
  covariance <- matrix(0, nrow = K, ncol = K)
  for (block in .known_v_blocks(known_V)) {
    index <- block[["index"]]
    covariance[index, index] <- block[["covariance"]]
  }

  covariance
}


.check_marglik_available <- function(object, caller) {

  if (inherits(object, "RoBMA")) {
    stop(
      "Marginal likelihood is not available for product-space ",
      "model-averaging objects (BMA.norm, BMA.glmm, RoBMA).",
      call. = FALSE
    )
  }
  if (.is_random(object) &&
      !(inherits(object, "brma.mv") && .is_data_known_v(object[["data"]]))) {
    .check_random_formula_postfit_deferred(object, caller)
  }

  # Public constructors reject p-hacking/composed selection kernels earlier.
  if (.marglik_has_composed_bias(object[["priors"]][["outcome"]][["bias"]])) {
    stop(
      "Marginal likelihood is not available for combined ",
      "prior_bias(selection, phacking) models yet.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.marglik_has_composed_bias <- function(prior) {

  if (is.null(prior)) {
    return(FALSE)
  }
  if (BayesTools::is.prior.mixture(prior)) {
    return(any(vapply(as.list(prior), .marglik_has_composed_bias, logical(1))))
  }
  if (!.is_prior_bias_kernel(prior)) {
    return(FALSE)
  }

  return(.prior_has_phacking(prior))
}


.marglik_add_selection_bridge_data <- function(fit_data, priors, effect_direction) {

  if (!.is_priors_weightfunction(priors) || is.null(fit_data[["yi"]])) {
    return(fit_data)
  }

  selection_spec <- .selection_spec(
    priors           = priors,
    yi               = fit_data[["yi"]],
    sei              = fit_data[["sei"]],
    effect_direction = effect_direction,
    signed_data      = TRUE
  )

  if (is.null(selection_spec)) {
    return(fit_data)
  }

  fit_data[["sel_kernel_mode"]]           <- selection_spec[["kernel_mode"]]
  fit_data[["sel_segment_bounds"]]        <- .selection_jags_bounds(selection_spec[["segments"]][["bounds"]])
  fit_data[["sel_segment_step_bin"]]      <- selection_spec[["segments"]][["step_bin"]]
  fit_data[["sel_segment_phack_region"]]  <- selection_spec[["segments"]][["phack_region"]]
  fit_data[["sel_phack_q"]]               <- selection_spec[["phack_q"]]
  fit_data[["sel_phack_z_source"]]        <- selection_spec[["phack_z_source"]]
  fit_data[["sel_phack_z_dest"]]          <- selection_spec[["phack_z_dest"]]

  return(fit_data)
}


.marglik_bridge_sd_source_spec <- function(add_parameters, fit, K) {

  posterior_names <- colnames(suppressWarnings(as.matrix(coda::as.mcmc(fit))))
  parameters      <- character()
  lb              <- numeric()
  ub              <- numeric()

  add_parameter <- function(parameter, lower, upper) {

    if (parameter %in% parameters) {
      return(invisible(NULL))
    }

    parameters <<- c(parameters, parameter)
    lb <<- c(lb, stats::setNames(lower, parameter))
    ub <<- c(ub, stats::setNames(upper, parameter))

    invisible(NULL)
  }

  for (parameter in add_parameters) {
    if (parameter %in% posterior_names) {
      add_parameter(parameter, 0, Inf)
      next
    }

    indexed_parameter <- paste0(parameter, "[", seq_len(K), "]")
    if (all(indexed_parameter %in% posterior_names)) {
      for (indexed in indexed_parameter) {
        add_parameter(indexed, 0, Inf)
      }
      next
    }

    add_parameter(parameter, 0, Inf)
  }

  if (length(parameters) == 0L) {
    return(list(
      parameters = character(),
      bounds     = NULL
    ))
  }

  return(list(
    parameters = parameters,
    bounds     = list(
      lb = lb[parameters],
      ub = ub[parameters]
    )
  ))
}


# ---------------------------------------------------------------------------- #
# .log_posterior
# ---------------------------------------------------------------------------- #
#
# Compute the log-likelihood of the data given a single posterior sample.
# This function is called by BayesTools::JAGS_bridgesampling for each
# posterior sample during bridge sampling.
#
# The function converts single-sample parameters to 1-row matrices to reuse
# the existing evaluate and pdf helper functions from evaluate.R and pdf.R.
# (in contrast to .log_lik.brma, this function does not add likelihood contributions from
#  estimate-level effects / baserate/ lograte for GLMMs since those are automatically 
#  handled via JAGS_bridgesampling function -- i.e., all parameters generated directly 
#  by BayesTools are automatically evaluated in the likelihood calculation, the only 
#  part missing is the final p(data | parameters) term) 
#
# @param parameters       named list of parameter values (preprocessed by BayesTools)
#                         - mu: scalar (no mods) or vector of length K (with mods)
#                         - tau: scalar (no scale regression)
#                         - log_tau: vector of length K (scale regression)
#                         - rho: scalar (if multilevel)
#                         - gamma: vector of length n_clusters (if multilevel)
#                         - PET: scalar (if PET model)
#                         - PEESE: scalar (if PEESE model)
#                         - omega: vector of weights (if weightfunction)
#                         - pi: vector of baserates (if binomial)
#                         - phi: vector of log-rates (if Poisson)
#                         - theta: vector of random effects (if GLMM)
# @param data             list containing fit_data (from .create_fit_data)
# @param is_mods          logical; whether model has moderators
# @param is_scale         logical; whether model has scale regression
# @param is_multilevel    logical; whether model is multilevel
# @param is_weights       logical; whether model uses weights (currently unused)
# @param is_PET           logical; whether model includes PET
# @param is_PEESE         logical; whether model includes PEESE
# @param is_weightfunction logical; whether model includes selection weights
# @param effect_direction character; "positive" or "negative"
# @param outcome_type     character; "norm", "bin", or "pois"
#
# @return scalar; sum of log-likelihoods across all observations
#
# ---------------------------------------------------------------------------- #
.log_posterior <- function(
    parameters, data,
    is_mods, is_scale, is_multilevel, is_weights,
    is_known_v, is_PET, is_PEESE, is_weightfunction, effect_direction,
    outcome_type, is_random = FALSE, model_data = NULL,
    known_V = NULL, known_v_backend = NULL,
    marginalized_variance_plan = NULL,
    bridge_context = NULL, fixed_tau = NULL, fixed_rho = NULL,
    fixed_zero_random = FALSE, cluster_effects_marginalized = FALSE,
    sampling_latent_marginalized = FALSE,
    covariance_plan_cache = NULL) {

  ### extract number of observations
  K <- data[["K"]]

  if (is_known_v) {
    if (is.null(known_V)) {
      if (is.null(model_data) || !.is_data_known_v(model_data)) {
        stop("Known-V bridge likelihood requires current model metadata.",
             call. = FALSE)
      }
      known_V <- .data_known_v_data(model_data)
    }
    if (is.null(known_v_backend)) {
      known_v_backend <- .known_v_effective_backend(known_V)
    }
  }
  if (is_known_v && is_weights) {
    stop(
      "Known-V bridge likelihoods do not support likelihood weights.",
      call. = FALSE
    )
  }

  ### compute mu samples as 1 x K matrix
  # BayesTools evaluates the formula and returns:
  # - scalar mu (no mods) -> replicate to K columns

  # - vector mu of length K (with mods) -> use directly
  mu_samples <- .marglik_get_mu_samples(
    parameters       = parameters,
    is_mods          = is_mods,
    is_PET           = is_PET,
    is_PEESE         = is_PEESE,
    effect_direction = effect_direction,
    sei              = if (outcome_type == "norm") data[["sei"]] else NULL,
    K                = K
  )

  ### compute tau samples as 1 x K matrix
  # BayesTools returns:
  # - scalar tau (no scale regression) -> replicate to K columns
  # - vector log_tau of length K (scale regression) -> exponentiate, use directly
  tau_result <- if (is_random) {
    .marglik_get_zero_tau_samples(K)
  } else {
    .marglik_get_tau_samples(
      parameters      = parameters,
      is_scale        = is_scale,
      is_multilevel   = is_multilevel,
      K               = K,
      fixed_tau       = fixed_tau,
      fixed_rho       = fixed_rho
    )
  }
  tau_within_samples  <- tau_result[["tau_within"]]
  tau_between_samples <- tau_result[["tau_between"]]

  ### add cluster-level (gamma) contribution for multilevel models
  if (is_multilevel && !isTRUE(cluster_effects_marginalized)) {
    cluster_contribution <- .marglik_get_cluster_effects(
      parameters  = parameters,
      tau_between = tau_between_samples,
      cluster     = data[["cluster"]]
    )
    mu_samples <- mu_samples + cluster_contribution
  }

  if (is_known_v && !isTRUE(sampling_latent_marginalized)) {
    sampling_dependency <- .marglik_get_sampling_dependency(
      parameters       = parameters,
      known_V          = known_V,
      effect_direction = effect_direction,
      K                = K
    )
    mu_samples <- mu_samples + sampling_dependency
  }

  ### dispatch to appropriate log-likelihood computation based on outcome type
  if (outcome_type == "norm") {

    # for selection models, PET, and PEESE, the outcome is in "positive" space
    # (yi flipped for negative effect direction in .create_fit_data)
    # so we need to flip mu_samples to match
    if (effect_direction == "negative") {
      mu_samples <- -mu_samples
    }

    if (is_known_v) {

      log_lik <- .marglik_known_v_norm_log_lik(
        parameters               = parameters,
        data                     = data,
        model_data               = model_data,
        known_V                  = known_V,
        known_v_backend          = known_v_backend,
        bridge_context           = bridge_context,
        fixed_zero_random        = fixed_zero_random,
        marginalized_variance_plan = marginalized_variance_plan,
        sampling_latent_marginalized = sampling_latent_marginalized,
        covariance_plan_cache    = covariance_plan_cache,
        mu_samples               = mu_samples,
        tau_within_samples       = tau_within_samples,
        is_random                = is_random,
        is_weightfunction        = is_weightfunction,
        effect_direction         = effect_direction,
        K                        = K
      )

    } else if (is_weightfunction) {

      selection_context <- .marglik_selection_context(parameters, data)
      log_lik <- .outcome_pdf.selnorm(
        yi                = data[["yi"]],
        mu_samples        = mu_samples,
        tau_within        = tau_within_samples,
        sei               = data[["sei"]],
        selection_sei     = data[["sei"]],
        selection_context = selection_context
      )

    } else if (isTRUE(cluster_effects_marginalized)) {

      log_lik <- .marglik_cluster_norm_log_lik(
        yi          = data[["yi"]],
        mu          = as.numeric(mu_samples),
        tau_within  = as.numeric(tau_within_samples),
        tau_between = as.numeric(tau_between_samples),
        sei         = data[["sei"]],
        cluster     = data[["cluster"]],
        weights     = if (is_weights) data[["weight"]] else NULL
      )

    } else {

      # standard normal likelihood
      log_lik <- .outcome_pdf.norm(
        yi         = data[["yi"]],
        mu_samples = mu_samples,
        tau_within = tau_within_samples,
        sei        = data[["sei"]]
      )

    }

  } else if (outcome_type == "bin") {

    # add GLMM random effects (theta * tau_within) for binomial models
    theta_contribution <- .marglik_get_theta_samples(
      parameters = parameters,
      tau_within = tau_within_samples,
      K          = K
    )
    mu_samples <- mu_samples + theta_contribution

    # get baserate samples
    logit_baserate <- .marglik_get_baserate_samples(
      parameters = parameters,
      K          = K
    )

    log_lik <- .outcome_pdf.binom_conditional(
      ai             = data[["ai"]],
      ci             = data[["ci"]],
      n1i            = data[["n1i"]],
      n2i            = data[["n2i"]],
      mu_samples     = mu_samples,
      logit_baserate = logit_baserate
    )

  } else if (outcome_type == "pois") {

    # add GLMM random effects (theta * tau_within) for Poisson models
    theta_contribution <- .marglik_get_theta_samples(
      parameters = parameters,
      tau_within = tau_within_samples,
      K          = K
    )
    mu_samples <- mu_samples + theta_contribution

    # get log-rate samples
    log_phi <- .marglik_get_lograte_samples(
      parameters = parameters,
      K          = K
    )

    log_lik <- .outcome_pdf.pois_conditional(
      x1i        = data[["x1i"]],
      x2i        = data[["x2i"]],
      t1i        = data[["t1i"]],
      t2i        = data[["t2i"]],
      mu_samples = mu_samples,
      log_phi    = log_phi
    )

  }

  if (is_weights && !isTRUE(cluster_effects_marginalized)) {
    log_lik <- .apply_log_lik_weights(log_lik, data[["weight"]])
  }

  ### return sum of log-likelihoods (scalar)
  return(sum(log_lik))
}


# ---------------------------------------------------------------------------- #
# Helper functions for .log_posterior
# ---------------------------------------------------------------------------- #
#
# These functions convert single-sample parameters from BayesTools to 1-row
# matrices compatible with the evaluate and pdf helper functions.
# ---------------------------------------------------------------------------- #


#' @keywords internal
.marglik_get_mu_samples <- function(parameters, is_mods, is_PET, is_PEESE,
                                    effect_direction, sei, K) {

  # BayesTools returns:
  # - scalar mu (no mods)
  #   => replicate to 1 x K matrix
  # - vector mu of length K (with mods, formula already evaluated)
  #   => vector of length K -> 1 x K matrix
  mu_samples <- matrix(parameters[["mu"]], nrow = 1, ncol = K)

  # add PET adjustment (PET * sei)
  # direction multiplier: +1 for positive, -1 for negative
  if (is_PET) {
    direction <- ifelse(effect_direction == "negative", -1, 1)
    mu_samples <- mu_samples + direction * parameters[["PET"]] * sei
  }

  # add PEESE adjustment (PEESE * sei^2)
  if (is_PEESE) {
    direction <- ifelse(effect_direction == "negative", -1, 1)
    mu_samples <- mu_samples + direction * parameters[["PEESE"]] * sei^2
  }

  return(mu_samples)
}


#' @keywords internal
.marglik_get_tau_samples <- function(parameters, is_scale, is_multilevel, K,
                                     fixed_tau = NULL, fixed_rho = NULL) {

  # BayesTools returns:
  # - scalar tau (no scale regression)
  # - vector log_tau of length K (scale regression, need to exponentiate)
  if (is_scale) {
    # log_tau vector -> exponentiate and convert to 1 x K matrix
    tau_samples <- matrix(exp(parameters[["log_tau"]]), nrow = 1, ncol = K)
  } else {
    if (!is.null(fixed_tau) && is.null(parameters[["tau"]])) {
      tau_samples <- matrix(fixed_tau, nrow = 1, ncol = K)
    } else {
      # scalar tau -> replicate to 1 x K matrix
      tau_samples <- matrix(parameters[["tau"]], nrow = 1, ncol = K)
    }
  }

  rho <- if (is_multilevel && !is.null(fixed_rho)) {
    fixed_rho
  } else if (is_multilevel) {
    parameters[["rho"]]
  } else {
    NULL
  }

  return(.heterogeneity_components(
    tau_total     = tau_samples,
    rho           = rho,
    is_multilevel = is_multilevel,
    context       = "Marginal-likelihood heterogeneity"
  ))
}


.marglik_get_zero_tau_samples <- function(K) {

  tau_samples <- matrix(0, nrow = 1L, ncol = K)

  return(list(
    tau_total   = tau_samples,
    tau_within  = tau_samples,
    tau_between = tau_samples
  ))
}


#' @keywords internal
.marglik_get_cluster_effects <- function(parameters, tau_between, cluster) {

  K <- ncol(tau_between)

  if (is.null(parameters[["gamma"]]) && all(tau_between == 0)) {
    return(matrix(0, nrow = 1L, ncol = K))
  }

  # extract gamma samples (vector of length n_clusters)
  # BayesTools returns gamma as a vector
  gamma <- parameters[["gamma"]]

  # gamma[cluster] maps each observation to its cluster-level effect
  # multiply by tau_between to get the contribution
  cluster_contribution <- matrix(gamma[cluster] * tau_between, nrow = 1, ncol = K)

  return(cluster_contribution)
}


# Restore exact IEEE infinities from selection bounds serialized for JAGS.
.marglik_restore_jags_selection_bounds <- function(x) {

  if (is.null(x)) {
    return(NULL)
  }

  negative_sentinel <- !is.na(x) & x == -1e300
  positive_sentinel <- !is.na(x) & x ==  1e300
  x[negative_sentinel] <- -Inf
  x[positive_sentinel] <-  Inf

  return(x)
}


.marglik_selection_context <- function(parameters, data) {

  z_lower        <- .marglik_restore_jags_selection_bounds(data[["sel_z_lower"]])
  z_upper        <- .marglik_restore_jags_selection_bounds(data[["sel_z_upper"]])
  segment_bounds <- .marglik_restore_jags_selection_bounds(
    data[["sel_segment_bounds"]]
  )
  n_bins         <- length(z_lower)
  phack_z_source <- if (!is.null(data[["phack_z_source"]])) {
    data[["phack_z_source"]]
  } else if (!is.null(data[["sel_phack_z_source"]])) {
    data[["sel_phack_z_source"]]
  } else {
    c(0, 0)
  }
  phack_z_dest <- if (!is.null(data[["phack_z_dest"]])) {
    data[["phack_z_dest"]]
  } else if (!is.null(data[["sel_phack_z_dest"]])) {
    data[["sel_phack_z_dest"]]
  } else {
    c(0, 0)
  }
  omega <- if (!is.null(parameters[["omega"]])) {
    matrix(parameters[["omega"]], nrow = 1)
  } else if (!is.null(data[["sel_omega"]])) {
    matrix(data[["sel_omega"]], nrow = 1)
  } else {
    matrix(1, nrow = 1, ncol = n_bins)
  }

  alpha <- if (!is.null(parameters[["alpha"]])) parameters[["alpha"]] else {
    if (!is.null(data[["sel_phack_alpha"]])) data[["sel_phack_alpha"]] else 0
  }
  phack_kind <- if (!is.null(parameters[["phack_kind"]])) {
    as.integer(parameters[["phack_kind"]])
  } else {
    if (!is.null(data[["sel_phack_kind"]])) {
      as.integer(data[["sel_phack_kind"]])
    } else if (!is.null(data[["sel_phack_q"]]) &&
               data[["sel_kernel_mode"]] %in% c(SELKERNEL_PHACK_POWER, SELKERNEL_STEP_PHACK_POWER)) {
      as.integer(data[["sel_phack_q"]])
    } else {
      0L
    }
  }

  selection_context <- list(
    kernel_mode    = data[["sel_kernel_mode"]],
    z_lower        = z_lower,
    z_upper        = z_upper,
    obs_bin        = data[["sel_obs_bin"]],
    sign           = data[["sel_sign"]],
    n_bins         = n_bins,
    has_step       = data[["sel_kernel_mode"]] %in% c(SELKERNEL_STEP, SELKERNEL_STEP_PHACK_POWER),
    has_phack      = data[["sel_kernel_mode"]] %in% c(SELKERNEL_PHACK_POWER, SELKERNEL_STEP_PHACK_POWER),
    phack_q        = if (phack_kind > 0L) phack_kind else 1L,
    phack_z_source = phack_z_source,
    phack_z_dest   = phack_z_dest,
    segments       = list(
      bounds       = segment_bounds,
      step_bin     = data[["sel_segment_step_bin"]],
      phack_region = data[["sel_segment_phack_region"]]
    ),
    omega          = omega,
    alpha          = alpha,
    phack_kind     = phack_kind
  )

  return(.selection_reset_native_cache(selection_context))
}


.marglik_get_sampling_dependency <- function(parameters, effect_direction, K,
                                             model_data = NULL, known_V = NULL) {

  if (is.null(known_V)) {
    if (is.null(model_data) || !.is_data_known_v(model_data)) {
      stop("Known-V latent dependency requires current model metadata.",
           call. = FALSE)
    }
    known_V <- .data_known_v_data(model_data)
  }
  sampling_rank <- .known_v_rank(known_V)
  if (sampling_rank == 0L) {
    return(matrix(0, nrow = 1, ncol = K))
  }
  latent_blocks <- .known_v_backend_blocks(known_V, "latent")

  sampling_z <- parameters[["sampling_z"]]
  if (is.null(sampling_z)) {
    sampling_z_names <- paste0("sampling_z[", seq_len(sampling_rank), "]")
    if (!all(sampling_z_names %in% names(parameters))) {
      stop("Missing known-V latent sampling factors.", call. = FALSE)
    }
    sampling_z <- unlist(parameters[sampling_z_names], use.names = FALSE)
  }

  sampling_dependency <- numeric(K)
  for (block in latent_blocks) {
    if (block[["rank"]] == 0L) {
      next
    }
    index <- block[["index"]]
    sampling_dependency[index] <- as.vector(
      block[["B"]] %*% sampling_z[block[["z_start"]]:block[["z_end"]]]
    )
  }
  sampling_dependency <- matrix(
    sampling_dependency,
    nrow = 1,
    ncol = K
  )

  if (effect_direction == "negative") {
    sampling_dependency <- -sampling_dependency
  }

  return(sampling_dependency)
}


.marglik_known_v_norm_log_lik <- function(parameters, data, model_data,
                                          known_V, known_v_backend,
                                          bridge_context,
                                          fixed_zero_random,
                                          marginalized_variance_plan,
                                          sampling_latent_marginalized,
                                          covariance_plan_cache,
                                          mu_samples, tau_within_samples,
                                          is_random, is_weightfunction,
                                          effect_direction, K) {

  extra_variance <- .marglik_known_v_extra_variance(
    parameters         = parameters,
    model_data         = model_data,
    bridge_context     = bridge_context,
    fixed_zero_random  = fixed_zero_random,
    marginalized_variance_plan = marginalized_variance_plan,
    tau_within_samples = tau_within_samples,
    is_random          = is_random,
    K                  = K
  )
  tau_within <- sqrt(extra_variance)
  marginal_random_covariance <- .marglik_bridge_random_covariance(
    bridge_context = bridge_context,
    K              = K,
    validation_cache = covariance_plan_cache
  )

  if (isTRUE(sampling_latent_marginalized) ||
      !is.null(marginal_random_covariance)) {
    if (is_weightfunction) {
      stop(
        "Shared Gaussian random-effect marginalization is not available for selection likelihoods.",
        call. = FALSE
      )
    }
    return(.marglik_known_v_integrated_gaussian_log_lik(
      data                       = data,
      model_data                 = model_data,
      mu_samples                 = mu_samples,
      extra_variance             = extra_variance,
      marginal_random_covariance = marginal_random_covariance,
      covariance_plan_cache      = covariance_plan_cache,
      effect_direction           = effect_direction,
      K                          = K
    ))
  }

  if (known_v_backend %in% c("latent", "diagonal")) {
    if (is_weightfunction) {
      selection_context <- .marglik_selection_context(parameters, data)
      return(.outcome_pdf.selnorm(
        yi                = data[["yi"]],
        mu_samples        = mu_samples,
        tau_within        = tau_within,
        sei               = sqrt(data[["sampling_var"]]),
        selection_sei     = data[["sei"]],
        selection_context = selection_context
      ))
    }

    return(.outcome_pdf.norm(
      yi         = data[["yi"]],
      mu_samples = mu_samples,
      tau_within = tau_within,
      sei        = sqrt(data[["sampling_var"]])
    ))
  }

  if (is_weightfunction) {
    stop(
      "Selection-model bridge likelihoods for known-V data require ",
      "known_v_parameterization = 'latent'.",
      call. = FALSE
    )
  }
  if (known_v_backend == "whitened") {
    return(.marglik_known_v_whitened_log_lik(
      data           = data,
      known_V        = known_V,
      mu_samples     = mu_samples,
      extra_variance = extra_variance,
      K              = K
    ))
  }
  if (known_v_backend == "block_mvn") {
    return(.marglik_known_v_block_mvn_log_lik_sum(
      model_data       = model_data,
      known_V          = known_V,
      mu_samples       = mu_samples,
      extra_variance   = extra_variance,
      effect_direction = effect_direction,
      covariance_plan_cache = covariance_plan_cache
    ))
  }

  stop("Unknown known-V parameterization: ", known_v_backend,
       call. = FALSE)
}


.marglik_cluster_norm_log_lik <- function(yi, mu, tau_within, tau_between,
                                           sei, cluster, weights = NULL) {

  K <- length(yi)
  if (length(mu) != K || length(tau_within) != K ||
      length(tau_between) != K || length(sei) != K ||
      length(cluster) != K) {
    stop("Internal error: invalid clustered bridge likelihood dimensions.",
         call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- rep(1, K)
  }
  if (length(weights) != K || anyNA(weights) || any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("Internal error: invalid clustered bridge likelihood weights.",
         call. = FALSE)
  }

  variance <- sei^2 + tau_within^2
  if (anyNA(variance) || any(!is.finite(variance)) || any(variance <= 0)) {
    return(-Inf)
  }
  cluster <- as.integer(cluster)
  if (anyNA(cluster) || any(cluster < 1L) ||
      !identical(sort(unique(cluster)), seq_len(max(cluster)))) {
    stop("Internal error: invalid clustered bridge likelihood indices.",
         call. = FALSE)
  }

  residual  <- yi - mu
  precision <- weights / variance
  adjustment_precision <- 1 + as.numeric(rowsum(
    precision * tau_between^2,
    group   = cluster,
    reorder = FALSE
  ))
  adjustment_score <- as.numeric(rowsum(
    precision * tau_between * residual,
    group   = cluster,
    reorder = FALSE
  ))

  normalization <- -0.5 * sum(weights * (log(2 * pi) + log(variance)))
  quadratic <- sum(precision * residual^2) -
    sum(adjustment_score^2 / adjustment_precision)

  normalization - 0.5 * sum(log(adjustment_precision)) - 0.5 * quadratic
}


.marglik_bridge_random_covariance <- function(bridge_context, K,
                                               validation_cache = NULL) {

  if (is.null(bridge_context) ||
      is.null(bridge_context[["marginalized_random"]]) ||
      is.null(bridge_context[["marginalized_random"]][["mu"]])) {
    return(NULL)
  }

  value <- bridge_context[["marginalized_random"]][["mu"]]
  representation <- value[["representation"]]
  if (is.null(representation) && !is.null(value[["covariance"]])) {
    representation <- "dense"
  }
  if (identical(representation, "factor")) {
    return(.marglik_bridge_random_covariance_factors(
      value            = value,
      K                = K,
      validation_cache = validation_cache
    ))
  }
  if (identical(representation, "factor_state")) {
    return(.marglik_bridge_random_covariance_states(
      value            = value,
      K                = K,
      validation_cache = validation_cache
    ))
  }
  if (!identical(representation, "dense")) {
    stop(
      "Bridge-marginalized random-effect covariance has an unknown representation.",
      call. = FALSE
    )
  }

  covariance <- value[["covariance"]]
  covariance <- .marglik_validate_random_covariance_matrix(covariance, K)

  list(
    representation = "dense",
    covariance = covariance
  )
}


.marglik_bridge_random_covariance_states <- function(
    value, K, validation_cache = NULL) {

  contract_id <- value[["contract_id"]]
  if (!is.environment(contract_id)) {
    stop(
      "Bridge-marginalized random-effect factor-state contract is missing its identity.",
      call. = FALSE
    )
  }
  process_id <- Sys.getpid()
  cached <- if (is.environment(validation_cache)) {
    validation_cache[["factor_state_validation"]]
  } else {
    NULL
  }
  if (!is.null(cached) && identical(cached[["process_id"]], process_id)) {
    if (!identical(contract_id, cached[["contract_id"]])) {
      stop(
        "Bridge-marginalized random-effect factor-state contract changed between evaluations.",
        call. = FALSE
      )
    }
    return(list(
      representation = "factor_state",
      row_blocks = cached[["row_blocks"]],
      factor_plans = cached[["factor_plans"]],
      factor_states = .marglik_validate_random_covariance_factor_states(
        states    = value[["factor_states"]],
        templates = cached[["factor_plans"]],
        K         = K
      )
    ))
  }

  factor_plans  <- value[["factor_plans"]]
  factor_states <- value[["factor_states"]]
  if (!is.list(factor_plans) || length(factor_plans) == 0L ||
      !is.list(factor_states) ||
      length(factor_states) != length(factor_plans)) {
    stop(
      "Bridge-marginalized random-effect factor-state contract is invalid.",
      call. = FALSE
    )
  }
  factors <- Map(c, factor_plans, factor_states)
  validated <- .marglik_bridge_random_covariance_factors(
    value = list(
      row_blocks = value[["row_blocks"]],
      factors = factors
    ),
    K = K,
    validation_cache = NULL
  )
  validated_plans  <- validated[["factor_plans"]]
  validated_states <- validated[["factor_states"]]
  if (is.environment(validation_cache)) {
    validation_cache[["factor_state_validation"]] <- list(
      process_id = process_id,
      contract_id = contract_id,
      row_blocks = validated[["row_blocks"]],
      factor_plans = validated_plans
    )
  }

  list(
    representation = "factor_state",
    row_blocks = validated[["row_blocks"]],
    factor_plans = validated_plans,
    factor_states = validated_states
  )
}


.marglik_validate_random_covariance_factor_states <- function(
    states, templates, K) {

  if (!is.list(states) || length(states) != length(templates)) {
    stop(
      "Bridge-marginalized random-effect covariance factor states changed between evaluations.",
      call. = FALSE
    )
  }

  out <- vector("list", length(states))
  for (factor_i in seq_along(states)) {
    state    <- states[[factor_i]]
    template <- templates[[factor_i]]
    type     <- template[["type"]]
    coefficient_structure <- template[["coefficient_structure"]]
    expected_names <- if (identical(type, "dense")) {
      "covariance"
    } else c(
      "coefficient_factor",
      if (identical(coefficient_structure, "markov")) c(
        "coefficient_scale",
        "markov_transition",
        "markov_innovation_variance"
      ) else character(),
      if (identical(type, "row_group")) "row_scale" else character()
    )
    if (!is.list(state) || !identical(names(state), expected_names)) {
      stop(
        "Bridge-marginalized random-effect covariance factor state structure changed between evaluations.",
        call. = FALSE
      )
    }
    if (identical(type, "dense")) {
      out[[factor_i]] <- list(
        covariance = .marglik_validate_random_covariance_matrix(
          state[["covariance"]],
          K
        )
      )
      next
    }
    n_columns <- ncol(template[["model_matrix"]])
    coefficient_factor <- state[["coefficient_factor"]]
    if (!is.matrix(coefficient_factor) ||
        !is.numeric(coefficient_factor) ||
        !identical(dim(coefficient_factor), c(n_columns, n_columns)) ||
        anyNA(coefficient_factor) || any(!is.finite(coefficient_factor))) {
      stop(
        "Bridge-marginalized random-effect coefficient factor is invalid.",
        call. = FALSE
      )
    }

    value <- list(coefficient_factor = coefficient_factor)
    if (identical(coefficient_structure, "markov")) {
      value <- c(
        value,
        .marglik_validate_random_covariance_markov_state(
          state,
          n_columns
        )
      )
    }
    if (identical(type, "row_group")) {
      row_scale <- state[["row_scale"]]
      if (!is.numeric(row_scale) || length(row_scale) != K ||
          anyNA(row_scale) || any(!is.finite(row_scale)) ||
          any(row_scale < 0)) {
        stop(
          "Bridge-marginalized random-effect row scale is invalid.",
          call. = FALSE
        )
      }
      value[["row_scale"]] <- as.double(row_scale)
    }
    out[[factor_i]] <- value
  }

  out
}


.marglik_bridge_random_covariance_factors <- function(
    value, K, validation_cache = NULL) {

  row_blocks <- value[["row_blocks"]]
  factors    <- value[["factors"]]
  process_id <- Sys.getpid()
  cached <- if (is.environment(validation_cache)) {
    validation_cache[["factor_validation"]]
  } else {
    NULL
  }
  if (!is.null(cached) && identical(cached[["process_id"]], process_id)) {
    if (!identical(row_blocks, cached[["source_row_blocks"]])) {
      stop(
        "Bridge-marginalized random-effect row blocks changed between evaluations.",
        call. = FALSE
      )
    }
    factors <- .marglik_validate_random_covariance_factor_values(
      factors   = factors,
      templates = cached[["factors"]],
      K         = K
    )
    return(.marglik_canonical_random_covariance_factors(
      row_blocks = cached[["row_blocks"]],
      factors    = factors
    ))
  }

  if (!is.list(row_blocks) || length(row_blocks) == 0L) {
    stop(
      "Bridge-marginalized random-effect row blocks are missing.",
      call. = FALSE
    )
  }
  row_blocks <- lapply(row_blocks, function(index) {
    if (!is.numeric(index) || length(index) == 0L || anyNA(index) ||
        any(!is.finite(index)) || any(index != as.integer(index)) ||
        any(index < 1L) || any(index > K) || anyDuplicated(index)) {
      stop(
        "Bridge-marginalized random-effect covariance block indices are invalid.",
        call. = FALSE
      )
    }
    as.integer(index)
  })
  if (!identical(sort(as.integer(unlist(row_blocks))), seq_len(K))) {
    stop(
      "Bridge-marginalized random-effect row blocks must partition every observation.",
      call. = FALSE
    )
  }
  if (!is.list(factors) || length(factors) == 0L) {
    stop(
      "Bridge-marginalized random-effect covariance factors are missing.",
      call. = FALSE
    )
  }
  factors <- lapply(
    factors,
    .marglik_validate_random_covariance_factor,
    K = K
  )
  if (is.environment(validation_cache)) {
    validation_cache[["factor_validation"]] <- list(
      process_id        = process_id,
      source_row_blocks = value[["row_blocks"]],
      row_blocks        = row_blocks,
      factors           = factors
    )
  }

  .marglik_canonical_random_covariance_factors(
    row_blocks = row_blocks,
    factors    = factors
  )
}


.marglik_canonical_random_covariance_factors <- function(row_blocks, factors) {

  list(
    representation = "factor_state",
    row_blocks      = row_blocks,
    factor_plans    = lapply(factors, .marglik_covariance_factor_plan),
    factor_states   = lapply(factors, .marglik_covariance_factor_state)
  )
}


.marglik_validate_random_covariance_factor_values <- function(
    factors, templates, K) {

  if (!is.list(factors) || length(factors) != length(templates)) {
    stop(
      "Bridge-marginalized random-effect covariance factors changed between evaluations.",
      call. = FALSE
    )
  }

  out <- vector("list", length(factors))
  for (factor_i in seq_along(factors)) {
    factor   <- factors[[factor_i]]
    template <- templates[[factor_i]]
    factor_coefficient_structure <- factor[["coefficient_structure"]]
    if (is.null(factor_coefficient_structure)) {
      factor_coefficient_structure <- "dense"
    }
    if (!is.list(factor) || !identical(factor[["type"]], template[["type"]]) ||
        !identical(factor_coefficient_structure,
                   template[["coefficient_structure"]])) {
      stop(
        "Bridge-marginalized random-effect covariance factor structure changed between evaluations.",
        call. = FALSE
      )
    }
    type <- template[["type"]]
    if (identical(type, "dense")) {
      out[[factor_i]] <- list(
        type = type,
        covariance = .marglik_validate_random_covariance_matrix(
          factor[["covariance"]],
          K
        )
      )
      next
    }

    model_matrix <- factor[["model_matrix"]]
    group_map    <- factor[["group_map"]]
    if (!is.matrix(model_matrix) || !is.numeric(model_matrix) ||
        !identical(dim(model_matrix), dim(template[["model_matrix"]]))) {
      stop(
        "Bridge-marginalized random-effect design matrix structure changed between evaluations.",
        call. = FALSE
      )
    }
    if (!identical(model_matrix, template[["model_matrix"]])) {
      stop(
        "Bridge-marginalized random-effect design matrix changed between evaluations.",
        call. = FALSE
      )
    }
    if (!is.numeric(group_map) || length(group_map) != K || anyNA(group_map) ||
        any(!is.finite(group_map)) || any(group_map != as.integer(group_map)) ||
        any(group_map < 1L)) {
      stop("Bridge-marginalized random-effect group mapping is invalid.",
           call. = FALSE)
    }
    group_map <- as.integer(group_map)
    if (!identical(group_map, template[["group_map"]])) {
      stop(
        "Bridge-marginalized random-effect group mapping changed between evaluations.",
        call. = FALSE
      )
    }
    n_columns <- ncol(template[["model_matrix"]])
    coefficient_factor <- factor[["coefficient_factor"]]
    if (!is.matrix(coefficient_factor) || !is.numeric(coefficient_factor) ||
        !identical(dim(coefficient_factor), c(n_columns, n_columns)) ||
        anyNA(coefficient_factor) || any(!is.finite(coefficient_factor))) {
      stop("Bridge-marginalized random-effect coefficient factor is invalid.",
           call. = FALSE)
    }
    coefficient_covariance <- factor[["coefficient_covariance"]]
    if (!is.null(coefficient_covariance)) {
      coefficient_covariance <- .marglik_validate_random_covariance_matrix(
        coefficient_covariance,
        n_columns
      )
      if (!isTRUE(all(coefficient_covariance ==
                      tcrossprod(coefficient_factor)))) {
        stop(
          "Bridge-marginalized random-effect coefficient covariance and factor disagree.",
          call. = FALSE
        )
      }
    }

    value <- list(
      type                   = type,
      model_matrix           = model_matrix,
      group_map              = group_map,
      coefficient_structure  = template[["coefficient_structure"]],
      coefficient_factor     = coefficient_factor
    )
    if (identical(template[["coefficient_structure"]], "markov")) {
      value <- c(
        value,
        .marglik_validate_random_covariance_markov_state(
          factor,
          n_columns
        )
      )
    }
    if (!is.null(coefficient_covariance)) {
      value[["coefficient_covariance"]] <- coefficient_covariance
    }
    if (identical(type, "row_group")) {
      row_scale <- factor[["row_scale"]]
      if (!is.numeric(row_scale) || length(row_scale) != K ||
          anyNA(row_scale) || any(!is.finite(row_scale)) ||
          any(row_scale < 0)) {
        stop("Bridge-marginalized random-effect row scale is invalid.",
             call. = FALSE)
      }
      value[["row_scale"]] <- as.double(row_scale)
    } else if (identical(type, "known_group")) {
      group_covariance <- .marglik_validate_random_covariance_matrix(
        factor[["group_covariance"]],
        nrow(template[["group_covariance"]])
      )
      if (!identical(group_covariance, template[["group_covariance"]])) {
        stop(
          "Bridge-marginalized known group covariance changed between evaluations.",
          call. = FALSE
        )
      }
      value[["group_covariance"]] <- group_covariance
    }
    out[[factor_i]] <- value
  }

  out
}


.marglik_validate_random_covariance_markov_state <- function(factor,
                                                              n_columns) {

  coefficient_scale <- factor[["coefficient_scale"]]
  transition <- factor[["markov_transition"]]
  innovation <- factor[["markov_innovation_variance"]]
  if (!is.numeric(coefficient_scale) ||
      length(coefficient_scale) != n_columns ||
      anyNA(coefficient_scale) || any(!is.finite(coefficient_scale)) ||
      any(coefficient_scale < 0) ||
      !is.numeric(transition) || length(transition) != n_columns - 1L ||
      anyNA(transition) || any(!is.finite(transition)) ||
      !is.numeric(innovation) || length(innovation) != n_columns - 1L ||
      anyNA(innovation) || any(!is.finite(innovation)) ||
      any(innovation <= 0)) {
    stop(
      "Bridge-marginalized random-effect Markov state is invalid.",
      call. = FALSE
    )
  }

  list(
    coefficient_scale = as.double(coefficient_scale),
    markov_transition = as.double(transition),
    markov_innovation_variance = as.double(innovation)
  )
}


.marglik_validate_random_covariance_factor <- function(factor, K) {

  if (!is.list(factor) || !is.character(factor[["type"]]) ||
      length(factor[["type"]]) != 1L || is.na(factor[["type"]])) {
    stop("Bridge-marginalized random-effect covariance factor is invalid.",
         call. = FALSE)
  }
  type <- factor[["type"]]
  if (identical(type, "dense")) {
    return(list(
      type = type,
      covariance = .marglik_validate_random_covariance_matrix(
        factor[["covariance"]],
        K
      )
    ))
  }
  if (!type %in% c("group", "row_group", "known_group")) {
    stop("Bridge-marginalized random-effect covariance factor type is unknown.",
         call. = FALSE)
  }
  model_matrix <- factor[["model_matrix"]]
  group_map    <- factor[["group_map"]]
  if (!is.matrix(model_matrix) || !is.numeric(model_matrix) ||
      nrow(model_matrix) != K || ncol(model_matrix) < 1L ||
      anyNA(model_matrix) || any(!is.finite(model_matrix))) {
    stop("Bridge-marginalized random-effect design matrix is invalid.",
         call. = FALSE)
  }
  if (!is.numeric(group_map) || length(group_map) != K || anyNA(group_map) ||
      any(!is.finite(group_map)) || any(group_map != as.integer(group_map)) ||
      any(group_map < 1L)) {
    stop("Bridge-marginalized random-effect group mapping is invalid.",
         call. = FALSE)
  }
  group_map <- as.integer(group_map)
  n_columns <- ncol(model_matrix)
  coefficient_structure <- factor[["coefficient_structure"]]
  if (is.null(coefficient_structure)) {
    coefficient_structure <- "dense"
  }
  if (!is.character(coefficient_structure) ||
      length(coefficient_structure) != 1L ||
      is.na(coefficient_structure) ||
      !coefficient_structure %in% c("dense", "diagonal", "markov")) {
    stop(
      "Bridge-marginalized random-effect coefficient structure is invalid.",
      call. = FALSE
    )
  }
  coefficient_factor <- factor[["coefficient_factor"]]
  if (!is.matrix(coefficient_factor) || !is.numeric(coefficient_factor) ||
      !identical(dim(coefficient_factor), c(n_columns, n_columns)) ||
      anyNA(coefficient_factor) || any(!is.finite(coefficient_factor))) {
    stop("Bridge-marginalized random-effect coefficient factor is invalid.",
         call. = FALSE)
  }
  coefficient_covariance <- factor[["coefficient_covariance"]]
  if (!is.null(coefficient_covariance)) {
    coefficient_covariance <- .marglik_validate_random_covariance_matrix(
      coefficient_covariance,
      n_columns
    )
    if (!isTRUE(all(coefficient_covariance ==
                    tcrossprod(coefficient_factor)))) {
      stop(
        "Bridge-marginalized random-effect coefficient covariance and factor disagree.",
        call. = FALSE
      )
    }
  }
  if (type %in% c("group", "row_group")) {
    out <- list(
      type = type,
      model_matrix = model_matrix,
      group_map = group_map,
      coefficient_structure = coefficient_structure,
      coefficient_factor = coefficient_factor
    )
    if (identical(coefficient_structure, "markov")) {
      markov_state <- .marglik_validate_random_covariance_markov_state(
        factor,
        n_columns
      )
      out <- c(out, markov_state)
    }
    if (!is.null(coefficient_covariance)) {
      out[["coefficient_covariance"]] <- coefficient_covariance
    }
    if (identical(type, "row_group")) {
      row_scale <- factor[["row_scale"]]
      if (!is.numeric(row_scale) || length(row_scale) != K ||
          anyNA(row_scale) || any(!is.finite(row_scale)) ||
          any(row_scale < 0)) {
        stop("Bridge-marginalized random-effect row scale is invalid.",
             call. = FALSE)
      }
      out[["row_scale"]] <- as.double(row_scale)
    }
    return(out)
  }

  group_covariance <- factor[["group_covariance"]]
  if (!is.matrix(group_covariance) || !is.numeric(group_covariance) ||
      nrow(group_covariance) != ncol(group_covariance) ||
      any(group_map > nrow(group_covariance))) {
    stop("Bridge-marginalized known group covariance is invalid.",
         call. = FALSE)
  }
  group_covariance <- .marglik_validate_random_covariance_matrix(
    group_covariance,
    nrow(group_covariance)
  )
  out <- list(
    type = type,
    model_matrix = model_matrix,
    group_map = group_map,
    group_covariance = group_covariance,
    coefficient_structure = coefficient_structure,
    coefficient_factor = coefficient_factor
  )
  if (!is.null(coefficient_covariance)) {
    out[["coefficient_covariance"]] <- coefficient_covariance
  }
  out
}


.marglik_validate_random_covariance_matrix <- function(covariance, K) {

  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      !identical(dim(covariance), c(K, K))) {
    stop(
      "Bridge-marginalized random-effect covariance has inconsistent dimensions.",
      call. = FALSE
    )
  }
  if (anyNA(covariance) || any(!is.finite(covariance))) {
    stop(
      "Bridge-marginalized random-effect covariance must be finite.",
      call. = FALSE
    )
  }
  if (K > 1L && any(covariance != t(covariance))) {
    stop(
      "Bridge-marginalized random-effect covariance must be symmetric.",
      call. = FALSE
    )
  }

  covariance
}


.marglik_known_v_integrated_gaussian_log_lik <- function(
    data, model_data, mu_samples, extra_variance,
    marginal_random_covariance, covariance_plan_cache,
    effect_direction, K) {

  sampling_covariance <- data[["marglik_sampling_covariance"]]
  if (!is.matrix(sampling_covariance) ||
      !identical(dim(sampling_covariance), c(K, K))) {
    stop(
      "Known-V shared-random bridge likelihood is missing its sampling covariance.",
      call. = FALSE
    )
  }

  yi <- model_data[["outcome"]][["yi"]]
  if (effect_direction == "negative") {
    yi <- -yi
  }

  block_indices <- data[["marglik_dependency_blocks"]]
  factor_representation <- !is.null(marginal_random_covariance) &&
    identical(
      marginal_random_covariance[["representation"]],
      "factor_state"
    )
  if (factor_representation) {
    random_indices <- marginal_random_covariance[["row_blocks"]]
    if (!is.null(block_indices) && !identical(block_indices, random_indices)) {
      stop(
        "Known-V bridge dependency blocks disagree with the random covariance contract.",
        call. = FALSE
      )
    }
    block_indices <- random_indices
  } else if (is.null(block_indices)) {
    covariance <- sampling_covariance +
      if (is.null(marginal_random_covariance)) 0 else
        marginal_random_covariance[["covariance"]]
    diag(covariance) <- diag(covariance) +
      as.numeric(extra_variance[1L, ])
    block_indices <- .known_v_block_indices(covariance)
  }
  if (!is.list(block_indices) ||
      !identical(sort(as.integer(unlist(block_indices))), seq_len(K))) {
    stop(
      "Known-V shared-random bridge dependency blocks are invalid.",
      call. = FALSE
    )
  }
  random_covariance_plans <- if (factor_representation) {
    marginal_random_covariance[["factor_plans"]]
  } else if (!is.null(marginal_random_covariance)) {
    list(list(type = "dense"))
  } else {
    list()
  }
  random_covariance_states <- if (factor_representation) {
    marginal_random_covariance[["factor_states"]]
  } else if (!is.null(marginal_random_covariance)) {
    list(list(covariance = marginal_random_covariance[["covariance"]]))
  } else {
    list()
  }

  .marglik_covariance_plan_loglik(
    cache                    = covariance_plan_cache,
    y                        = as.double(yi),
    mean                     = as.double(mu_samples[1L, ]),
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = random_covariance_plans,
    random_covariance_states = random_covariance_states,
    block_indices            = block_indices,
    extra_variance           = as.double(extra_variance[1L, ])
  )
}


.marglik_covariance_plan_loglik <- function(
    cache, y, mean, sampling_covariance, random_covariance_plans,
    random_covariance_states, block_indices, extra_variance) {

  plan <- .marglik_covariance_plan_get(
    cache                    = cache,
    y                        = y,
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = random_covariance_plans,
    block_indices            = block_indices
  )
  .Call(
    "RoBMA_known_v_covariance_plan_loglik",
    plan,
    mean,
    random_covariance_states,
    extra_variance,
    PACKAGE = "RoBMA"
  )
}


.marglik_covariance_plan_loglik_batch <- function(
    cache, y, means, sampling_covariance, random_covariance_plans,
    random_covariance_states, block_indices, extra_variances) {

  plan <- .marglik_covariance_plan_get(
    cache                    = cache,
    y                        = y,
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = random_covariance_plans,
    block_indices            = block_indices
  )
  .Call(
    "RoBMA_known_v_covariance_plan_loglik_batch",
    plan,
    t(means),
    random_covariance_states,
    t(extra_variances),
    PACKAGE = "RoBMA"
  )
}


.marglik_covariance_plan_group_iid_variance_grid_loglik <- function(
    cache, y, means, sampling_covariance, block_indices, group_variances,
    diagonal_variances) {

  plan <- .marglik_covariance_plan_get(
    cache                    = cache,
    y                        = y,
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = list(),
    block_indices            = block_indices
  )
  .Call(
    "RoBMA_known_v_covariance_plan_group_iid_variance_grid_loglik",
    plan,
    t(means),
    group_variances,
    diagonal_variances,
    PACKAGE = "RoBMA"
  )
}


.marglik_covariance_plan_location_quadratic_batch <- function(
    cache, y, means, bases, sampling_covariance, random_covariance_plans,
    random_covariance_states, block_indices, extra_variances) {

  plan <- .marglik_covariance_plan_get(
    cache                    = cache,
    y                        = y,
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = random_covariance_plans,
    block_indices            = block_indices
  )
  .Call(
    "RoBMA_known_v_covariance_plan_location_quadratic_batch",
    plan,
    t(means),
    t(bases),
    random_covariance_states,
    t(extra_variances),
    PACKAGE = "RoBMA"
  )
}


.marglik_covariance_plan_conditional_loglik_batch <- function(
    cache, y, means, sampling_covariance, random_covariance_plans,
    random_covariance_states, block_indices, extra_variances) {

  plan <- .marglik_covariance_plan_get(
    cache                    = cache,
    y                        = y,
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = random_covariance_plans,
    block_indices            = block_indices
  )
  t(.Call(
    "RoBMA_known_v_covariance_plan_conditional_loglik_batch",
    plan,
    t(means),
    random_covariance_states,
    t(extra_variances),
    PACKAGE = "RoBMA"
  ))
}


.marglik_covariance_plan_conditional_summary_batch <- function(
    cache, y, means, sampling_covariance, random_covariance_plans,
    random_covariance_states, block_indices, extra_variances) {

  plan <- .marglik_covariance_plan_get(
    cache                    = cache,
    y                        = y,
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = random_covariance_plans,
    block_indices            = block_indices
  )
  summary <- .Call(
    "RoBMA_known_v_covariance_plan_conditional_summary_batch",
    plan,
    t(means),
    random_covariance_states,
    t(extra_variances),
    PACKAGE = "RoBMA"
  )
  summary[["residual"]] <- t(summary[["residual"]])
  summary[["variance"]] <- t(summary[["variance"]])

  return(summary)
}


.marglik_covariance_plan_get <- function(
    cache, y, sampling_covariance, random_covariance_plans, block_indices) {

  if (!is.null(cache) && !is.environment(cache)) {
    stop("Known-V covariance plan cache must be an environment.",
         call. = FALSE)
  }

  process_id <- Sys.getpid()
  rebuild_plan <- is.null(cache) ||
    !identical(cache[["process_id"]], process_id) ||
    is.null(cache[["plan"]])
  if (rebuild_plan) {
    plan <- .Call(
      "RoBMA_known_v_covariance_plan_create",
      y,
      sampling_covariance,
      random_covariance_plans,
      block_indices,
      PACKAGE = "RoBMA"
    )
    if (is.environment(cache)) {
      cache[["plan"]]       <- plan
      cache[["process_id"]] <- process_id
    }
  } else {
    plan <- cache[["plan"]]
  }

  return(plan)
}


.marglik_covariance_factor_plan <- function(factor) {

  if (identical(factor[["type"]], "dense")) {
    return(list(type = "dense"))
  }

  plan <- factor[c(
    "type",
    "model_matrix",
    "group_map",
    "coefficient_structure"
  )]
  if (identical(factor[["type"]], "known_group")) {
    plan[["group_covariance"]] <- factor[["group_covariance"]]
  }
  plan
}


.marglik_covariance_factor_state <- function(factor) {

  if (identical(factor[["type"]], "dense")) {
    return(list(covariance = factor[["covariance"]]))
  }

  state <- list(
    coefficient_factor = factor[["coefficient_factor"]]
  )
  if (identical(factor[["coefficient_structure"]], "markov")) {
    state[["coefficient_scale"]] <- factor[["coefficient_scale"]]
    state[["markov_transition"]] <- factor[["markov_transition"]]
    state[["markov_innovation_variance"]] <-
      factor[["markov_innovation_variance"]]
  }
  if (identical(factor[["type"]], "row_group")) {
    state[["row_scale"]] <- factor[["row_scale"]]
  }
  state
}


.marglik_known_v_extra_variance <- function(parameters, model_data,
                                            bridge_context,
                                            tau_within_samples, is_random, K,
                                            fixed_zero_random = FALSE,
                                            marginalized_variance_plan = NULL) {

  extra_variance <- if (isTRUE(fixed_zero_random)) {
    matrix(0, nrow = 1L, ncol = K)
  } else if (is_random && !is.null(marginalized_variance_plan)) {
    .marglik_evaluate_marginalized_variance_plan(
      plan           = marginalized_variance_plan,
      bridge_context = bridge_context,
      K              = K
    )
  } else if (is_random) {
    if (is.null(model_data)) {
      stop(
        "Random-formula known-V bridge likelihood requires model metadata.",
        call. = FALSE
      )
    }
    .evaluate_marginalized_random_variance(
      data              = model_data,
      posterior_samples = .marglik_bridge_posterior_samples(
        parameters     = parameters,
        bridge_context = bridge_context
      ),
      K                 = K
    )
  } else {
    tau_within_samples^2
  }
  extra_variance <- as.matrix(extra_variance)

  if (nrow(extra_variance) != 1L || ncol(extra_variance) != K) {
    stop(
      "Known-V bridge diagonal variance contributions have inconsistent dimensions.",
      call. = FALSE
    )
  }
  if (any(!is.finite(extra_variance)) || any(extra_variance < 0)) {
    stop(
      "Known-V bridge diagonal variance contributions must be finite and non-negative.",
      call. = FALSE
    )
  }

  return(extra_variance)
}


.marglik_evaluate_marginalized_variance_plan <- function(
    plan, bridge_context, K) {

  if (!is.list(plan) || !identical(plan[["K"]], K) ||
      !is.list(plan[["terms"]])) {
    stop("Known-V bridge marginalized-variance plan is invalid.",
         call. = FALSE)
  }
  nodes <- if (inherits(bridge_context, "BayesTools_bridge_context")) {
    bridge_context[["nodes"]]
  } else {
    NULL
  }
  variance <- numeric(K)
  for (term in plan[["terms"]]) {
    parameter <- term[["parameter"]]
    if (is.null(nodes) || !parameter %in% names(nodes)) {
      stop(
        "Bridge context is missing marginalized random-effect SD node: ",
        parameter,
        call. = FALSE
      )
    }
    scale <- nodes[[parameter]]
    if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) ||
        !is.finite(scale) || scale < 0) {
      stop(
        "Bridge marginalized random-effect SD node must be finite and non-negative: ",
        parameter,
        call. = FALSE
      )
    }
    variance <- variance + scale^2 * term[["multiplier"]]
  }

  matrix(variance, nrow = 1L, ncol = K)
}


.marglik_known_v_whitened_log_lik <- function(data, mu_samples,
                                              extra_variance, K,
                                              known_V) {

  whitening_blocks <- .known_v_backend_blocks(known_V, "whitened")
  log_lik          <- numeric(K)
  variance         <- numeric(K)

  independent <- .known_v_independent_indices(known_V)
  if (length(independent) > 0L) {
    variance[independent] <- .known_v_diagonal(known_V)[independent] +
      as.numeric(extra_variance[1L, independent])
    if (any(!is.finite(variance[independent])) ||
        any(variance[independent] <= 0)) {
      stop("Known-V whitened bridge variances must be positive.",
           call. = FALSE)
    }
    log_lik[independent] <- stats::dnorm(
      x    = data[["known_v_independent_y"]],
      mean = as.numeric(mu_samples[1L, independent]),
      sd   = sqrt(variance[independent]),
      log  = TRUE
    )
  }

  for (b in seq_along(whitening_blocks)) {
    block        <- whitening_blocks[[b]]
    index        <- block[["index"]]
    whitening_mu <- as.vector(
      block[["rotation"]] %*% as.numeric(mu_samples[1L, index])
    )
    variance[index] <- block[["variance"]] +
      as.numeric(extra_variance[1L, index])
    if (any(!is.finite(variance[index])) || any(variance[index] <= 0)) {
      stop("Known-V whitened bridge variances must be positive.",
           call. = FALSE)
    }
    log_lik[index] <- stats::dnorm(
      x    = data[[paste0("whitening_y_", b)]],
      mean = whitening_mu,
      sd   = sqrt(variance[index]),
      log  = TRUE
    )
  }

  if (any(!is.finite(variance)) || any(variance <= 0)) {
    stop("Known-V whitened bridge variances must be positive.", call. = FALSE)
  }

  return(matrix(log_lik, nrow = 1L, ncol = K))
}


.marglik_known_v_block_mvn_log_lik_sum <- function(model_data, known_V,
                                                   mu_samples,
                                                   extra_variance,
                                                   effect_direction,
                                                   covariance_plan_cache) {

  yi <- model_data[["outcome"]][["yi"]]

  if (effect_direction == "negative") {
    yi <- -yi
  }

  .marglik_covariance_plan_loglik(
    cache                    = covariance_plan_cache,
    y                        = as.double(yi),
    mean                     = as.double(mu_samples[1L, ]),
    sampling_covariance      = .marglik_known_v_covariance_matrix(known_V),
    random_covariance_plans  = list(),
    random_covariance_states = list(),
    block_indices            = lapply(.known_v_blocks(known_V), `[[`, "index"),
    extra_variance           = as.double(extra_variance[1L, ])
  )
}


.marglik_mvn_log_density <- function(y, mean, covariance,
                                     context = "bridge") {

  chol_covariance <- .known_v_chol_covariance(
    covariance = covariance,
    context    = context
  )

  residual <- y - mean
  z        <- backsolve(chol_covariance, residual, transpose = TRUE)
  size     <- length(y)

  -0.5 * (
    size * log(2 * pi) +
      2 * sum(log(diag(chol_covariance))) +
      sum(z^2)
  )
}


.parameters_as_sample_matrix <- function(parameters) {

  values <- numeric()
  for (name in names(parameters)) {
    value <- parameters[[name]]
    if (length(value) == 0L) {
      next
    }
    value <- as.numeric(value)
    names(value) <- if (length(value) == 1L) {
      name
    } else {
      paste0(name, "[", seq_along(value), "]")
    }
    values <- c(values, value)
  }

  posterior_row <- attr(parameters, "posterior_samples", exact = TRUE)
  if (!is.null(posterior_row)) {
    posterior_row <- as.matrix(posterior_row)
    if (nrow(posterior_row) != 1L || is.null(colnames(posterior_row))) {
      stop("Bridge posterior row metadata have invalid dimensions.",
           call. = FALSE)
    }
    extra <- setdiff(colnames(posterior_row), names(values))
    if (length(extra) > 0L) {
      extra_values <- as.numeric(posterior_row[1L, extra, drop = TRUE])
      names(extra_values) <- extra
      values <- c(values, extra_values)
    }
  }

  out <- matrix(values, nrow = 1L)
  colnames(out) <- names(values)
  return(out)
}


.marglik_bridge_posterior_samples <- function(parameters, bridge_context = NULL) {

  posterior_samples <- .parameters_as_sample_matrix(parameters)

  if (!inherits(bridge_context, "BayesTools_bridge_context")) {
    return(posterior_samples)
  }

  nodes <- bridge_context[["nodes"]]
  if (length(nodes) == 0L) {
    return(posterior_samples)
  }
  if (!is.numeric(nodes) || is.null(names(nodes))) {
    stop("Bridge context nodes must be a named numeric vector.",
         call. = FALSE)
  }

  nodes <- nodes[!is.na(names(nodes)) & nzchar(names(nodes))]
  extra <- setdiff(names(nodes), colnames(posterior_samples))
  if (length(extra) == 0L) {
    return(posterior_samples)
  }

  extra_samples <- matrix(
    as.numeric(nodes[extra]),
    nrow     = 1L,
    dimnames = list(NULL, extra)
  )
  posterior_samples <- cbind(posterior_samples, extra_samples)

  return(posterior_samples)
}


#' @keywords internal
.marglik_get_theta_samples <- function(parameters, tau_within, K) {

  # theta is a vector of length K (estimate-level random effects for GLMM)
  theta <- parameters[["theta"]]

  # theta contribution = theta[k] * tau_within[k]
  theta_contribution <- matrix(theta * tau_within, nrow = 1, ncol = K)

  return(theta_contribution)
}


#' @keywords internal
.marglik_get_baserate_samples <- function(parameters, K) {

  # pi is a vector of length K (baserates for binomial models)
  pi <- parameters[["pi"]]

  # return logit(pi) as 1 x K matrix
  logit_baserate <- matrix(.logit(pi), nrow = 1, ncol = K)

  return(logit_baserate)
}


#' @keywords internal
.marglik_get_lograte_samples <- function(parameters, K) {

  # phi is a vector of length K (log-rates for Poisson models)
  phi <- parameters[["phi"]]

  # return log(phi) as 1 x K matrix
  # Note: phi is already on log scale in JAGS model
  log_phi <- matrix(phi, nrow = 1, ncol = K)

  return(log_phi)
}
