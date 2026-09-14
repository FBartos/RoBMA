# ============================================================================ #
# selection-likelihood.R
# ============================================================================ #
#
# Conditional Gaussian vector selection preparation. Generic numerical and
# random-covariance compilation is owned by BayesTools; this file supplies the
# meta-analytic covariance and constructor integration.
#
# ============================================================================ #


SELNORM_CLUSTER_QUADRATURE_ORDERS <-
  c(15L, 31L, 63L, 127L, 255L, 511L, 1023L)
SELNORM_FACTOR_QUADRATURE_ORDERS <- list(
  `2` = c(3L, 5L, 7L, 9L, 11L, 15L, 21L, 31L, 43L, 63L, 95L),
  `3` = c(3L, 5L, 7L, 9L, 11L, 15L, 21L, 31L, 43L, 63L),
  `4` = c(3L, 5L, 7L, 9L, 11L, 15L, 21L)
)


#' Control numerical integration of selection likelihoods
#'
#' @description
#' Creates numerical integration settings for the product and best-p-value
#' selection likelihoods used by [bselmodel()], [bselmodel.mv()], [RoBMA()], and
#' [RoBMA.mv()]. Covariance factors use deterministic sequences of
#' Gauss-Hermite rules. Certified factors of structural rank one through four
#' fall back to randomized quasi-Monte Carlo in the factor dimension when the
#' quadrature sequence does not converge. General covariance blocks with
#' product weights first try same-diagonal covariance approximations with a
#' common Gaussian factor and independent group factors. Monotone weights
#' permit covariance bounds on the normalizer. Two-sided and other nonmonotone
#' weights require a covariance-perturbation bound; agreement of auxiliary
#' normalizers alone is insufficient. The absolute bound uses the total
#' variation of the step weights. Strictly positive weights also permit a
#' relative Gaussian covariance bound, including reconstructed diagonal and
#' factor rounding, which remains useful for tiny normalizers. Their original
#' full covariance remains authoritative.
#' The covariance error bound, estimated
#' quadrature error, and a bound on Gaussian mass outside the quadrature nodes
#' share one numerical error budget. Interior quadrature error is estimated,
#' so this is not a rigorous probability interval. Unsupported or unresolved
#' cases use randomized quasi-Monte Carlo in the dense conditional-normal
#' representation.
#' Each fit or post-fit request uses a fixed integration design, so repeated
#' likelihood evaluations are deterministic. One-observation blocks are evaluated
#' analytically.
#'
#' Post-fit dense two-bin product normalizers can use a bounded interpolation
#' of the standard-normal CDF. Its deterministic error is included in the same
#' numerical error budget, separately from Monte Carlo standard error. Other
#' weights, arithmetic outside its supported range, or insufficient error
#' budgets retain the ordinary scalar primitive. An unresolved interpolation
#' attempt retries the ordinary covariance-envelope calculation before QMC.
#' Fitting and random generation retain their ordinary scalar primitive.
#'
#' The numerical defaults are package operating settings, not statistical
#' estimands or literature-mandated constants. `points_per_scramble`,
#' `max_points_per_scramble`, and `scrambles` set the fixed computational
#' design, while
#' `relative_tolerance` sets its diagnostic acceptance rule. Increasing a
#' budget changes only the numerical evaluation of the same likelihood target.
#'
#' For qCMDE/IWMDE post-fit densities and hypotheses in Gaussian selection
#' models, pass this object as
#' `density_control$integration_control`, for example
#' `density_control = list(integration_control = set_selection_likelihood_control(max_points_per_scramble = 32768))`.
#' The default `NULL` retains the fitted integration settings. An explicit
#' control rebuilds a temporary execution plan from the stored covariance
#' metadata without changing the fitted object or its posterior draws.
#' In fitting and qCMDE/IWMDE likelihood evaluation, the maximum point budget
#' controls factor QMC fallback. Post-fit full-event and partial-event
#' calculations, including z-plots and sensitivity diagnostics, also refine
#' dense QMC evaluations up to that budget. Analytic and deterministic
#' quadrature rules remain unchanged.
#'
#' @param points_per_scramble nominal number of integrand evaluations per
#'   randomized quasi-Monte Carlo scramble for certified factor
#'   fallbacks, and the number of shifted Halton points for general covariance
#'   blocks.
#'   Factor likelihood-normalizer fallbacks combine the prior Gaussian with
#'   proposals optimized from zero and selection-boundary starting points.
#'   Every finite optimized candidate is retained, including candidates that
#'   locate the same mode. The total count is divided as evenly as possible
#'   across these proposals, and importance weights use the actual allocation.
#'   Odd likelihood-normalizer factor budgets use the next even total because
#'   the fixed design uses two streams. This count
#'   is the initial nested-design budget, not an accuracy guarantee. The
#'   integration diagnostics determine whether it is adequate for each
#'   evaluated likelihood state.
#'   Selected-density projection fallbacks instead retain the prior Gaussian
#'   and one proposal optimized from zero.
#' @param max_points_per_scramble largest nested-design budget used for a
#'   certified factor fallback, or a post-fit Gaussian event
#'   calculation, when the initial `points_per_scramble` budget fails its
#'   integration diagnostic. Refinement doubles the point count up to this
#'   value while reusing the same fixed design.
#' @param scrambles number of independent randomized shifts used to estimate
#'   integration error for certified factors and general covariance
#'   blocks.
#' @param relative_tolerance largest accepted successive relative quadrature
#'   change for factor integration, or the largest nested-design relative
#'   change and relative Monte Carlo standard error for a factor
#'   fallback. General covariance blocks use the relative Monte Carlo standard
#'   error. Dense covariance-envelope quadrature additionally requires the
#'   sum of covariance-bound width, twice the successive quadrature change,
#'   the relative exterior-mass bound, and any post-fit CDF interpolation error
#'   to satisfy this tolerance. Factor quadrature also checks exterior Gaussian
#'   mass before accepting rule agreement; an unresolved rule sequence uses the
#'   existing QMC fallback. Zplot projections additionally check density error
#'   relative to the peak within each posterior draw and raw normalizer-mass
#'   diagnostics, including their bounded approximation and omission terms.
#'   These diagnostics do not certify the total integrated error of a compressed
#'   displayed curve; see [as_zplot.brma()].
#' @param seed non-negative integer used only to create fixed quasi-Monte Carlo
#'   designs. It does not alter R's random-number state.
#'
#' @return A `RoBMA_selection_likelihood_control` object.
#'
#' @export
set_selection_likelihood_control <- function(
    points_per_scramble = 512L, max_points_per_scramble = 8192L,
    scrambles = 8L,
    relative_tolerance = 0.005, seed = 1L) {

  BayesTools::check_int(
    points_per_scramble,
    "points_per_scramble",
    lower = 8L,
    check_length = 1L,
    allow_NA = FALSE
  )
  BayesTools::check_int(
    max_points_per_scramble,
    "max_points_per_scramble",
    lower = points_per_scramble,
    check_length = 1L,
    allow_NA = FALSE
  )
  BayesTools::check_int(
    scrambles,
    "scrambles",
    lower = 2L,
    check_length = 1L,
    allow_NA = FALSE
  )
  BayesTools::check_real(
    relative_tolerance,
    "relative_tolerance",
    lower = 0,
    upper = 1,
    check_length = 1L,
    allow_NA = FALSE
  )
  if (relative_tolerance == 0) {
    stop("'relative_tolerance' must be positive.", call. = FALSE)
  }
  BayesTools::check_int(
    seed,
    "seed",
    lower = 0L,
    check_length = 1L,
    allow_NA = FALSE
  )

  structure(
    list(
      points_per_scramble     = as.integer(points_per_scramble),
      max_points_per_scramble = as.integer(max_points_per_scramble),
      scrambles               = as.integer(scrambles),
      relative_tolerance      = as.numeric(relative_tolerance),
      seed                    = as.integer(seed)
    ),
    class = c("RoBMA_selection_likelihood_control", "list")
  )
}


.check_selection_likelihood_control <- function(
    control, argument = "selection_control") {

  if (!inherits(control, "RoBMA_selection_likelihood_control")) {
    stop(
      "'", argument, "' must be created by ",
      "set_selection_likelihood_control().",
      call. = FALSE
    )
  }
  do.call(set_selection_likelihood_control, unclass(control))
}


# Reuse fitted covariance metadata when a post-fit call requests another design.
.selection_joint_execution_plan_with_control <- function(plan, control) {

  if (all(vapply(names(control), function(name) {

    identical(plan[[name]], control[[name]])
  }, logical(1)))) {
    return(plan)
  }

  if (identical(plan[["statistical_target"]], "whole_sampling_error_selection")) {
    for (name in names(control)) plan[[name]] <- control[[name]]
    plan[["designs"]] <- lapply(plan[["designs"]], function(design) {
      if (length(design) == 1L) return(design)
      BayesTools::selection_qmc_design(
        dimensions = dim(design)[[3L]], points = control[["max_points_per_scramble"]],
        scrambles = control[["scrambles"]], seed = control[["seed"]]
      )
    })
    return(plan)
  }

  updated <- .selection_joint_execution_plan(
    row_blocks             = plan[["row_blocks"]],
    block_methods          = plan[["block_methods"]],
    factor_ranks           = plan[["factor_ranks"]],
    selection_control      = control,
    sampling               = plan[["sampling"]],
    sampling_factor_blocks = plan[["sampling_factor_blocks"]],
    random_covariance      = plan[["random_covariance"]]
  )
  updated[["selection_spec"]] <- plan[["selection_spec"]]
  updated
}


.data_selection_model <- function(data) {

  model <- attr(data, "selection_model", exact = TRUE)
  if (!is.null(model) &&
      (!inherits(model, "RoBMA_selection_model") ||
       !identical(model[["schema_version"]], 3L))) {
    stop("The fitted selection-model specification is invalid.", call. = FALSE)
  }
  model
}


.is_data_joint_selection <- function(data) {

  !is.null(.data_selection_model(data))
}


.selection_retains_estimate <- function(data) {

  model <- .data_selection_model(data)
  !is.null(model) && identical(model[["estimate_random_effects"]], "condition") &&
    isTRUE(model[["applicability"]][["estimate_random_effects"]])
}


.selection_integrates_estimate <- function(data) {

  model <- .data_selection_model(data)
  !is.null(model) && identical(model[["estimate_random_effects"]], "integrate") &&
    isTRUE(model[["applicability"]][["estimate_random_effects"]])
}


.selection_retains_other_random <- function(data) {

  model <- .data_selection_model(data)
  !is.null(model) && identical(model[["other_random_effects"]], "condition") &&
    isTRUE(model[["applicability"]][["other_random_effects"]])
}


.selection_retains_any_random <- function(data) {

  .selection_retains_estimate(data) || .selection_retains_other_random(data)
}


.selection_retains_sampling <- function(data) {

  model <- .data_selection_model(data)
  !is.null(model) && identical(model[["known_sampling_variance"]], "condition") &&
    isTRUE(model[["applicability"]][["known_sampling_variance"]])
}


.selection_all_sources_conditioned <- function(data) {

  if (!.selection_retains_sampling(data)) return(FALSE)
  sources <- .data_selection_model(data)[["sources"]][["random"]]
  all(vapply(sources, function(source) isTRUE(source[["retained"]]), logical(1L)))
}


.setup_uses_joint_selection_likelihood <- function(setup) {

  .is_data_joint_selection(setup[["data"]])
}


.data_selection_execution_plan <- function(data) {

  plan <- attr(data, "selection_execution_plan", exact = TRUE)
  if (.is_data_joint_selection(data) && is.null(plan)) {
    stop(
      "The fitted selection execution plan is unavailable.",
      call. = FALSE
    )
  }
  if (!is.null(plan) &&
      (!inherits(plan, "RoBMA_selection_execution_plan") ||
       !identical(plan[["schema_version"]], 5L))) {
    stop(
      "The fitted selection execution plan is invalid.",
      call. = FALSE
    )
  }
  plan
}


.prepare_selection_model_object <- function(object) {

  data   <- object[["data"]]
  priors <- .selection_bias_priors(object[["priors"]])
  models <- lapply(priors, BayesTools::selection_model_spec)
  active <- vapply(seq_along(priors), function(i) {

    !is.null(models[[i]]) && priors[[i]][["prior_weights"]] > 0
  }, logical(1))
  if (!any(active)) {
    attr(object[["data"]], "selection_binding") <- NULL
    return(object)
  }
  common <- models[[which(active)[[1L]]]]
  modes  <- c("estimate_random_effects", "other_random_effects", "known_sampling_variance")
  if (any(vapply(models[active], function(model) {

    !identical(model[modes], common[modes])
  }, logical(1)))) {
    stop(
      "Active weightfunction branches must use the same 'estimate_random_effects', ",
      "'other_random_effects', and 'known_sampling_variance' settings.",
      call. = FALSE
    )
  }

  random_sources <- list()
  if (.is_data_random(data)) {
    design <- .object_bayestools_formula_design(
      object                 = object,
      parameter              = "mu",
      source                 = "location",
      random_effects_compile = NULL
    )
    object[["formula_design"]][["mu"]] <- design
    roles <- BayesTools::random_effects_level_roles(
      design[["random_effects"]], n_rows = nrow(data[["outcome"]])
    )
    random_sources <- lapply(design[["random_effects"]], function(term) {

      name <- .random_effect_term_block_name(term)
      role <- roles[[name]]
      mode <- if (identical(role, "estimate")) {
        "estimate_random_effects"
      } else "other_random_effects"
      list(
        name     = name,
        role     = role,
        retained = identical(common[[mode]], "condition")
      )
    })
    for (term in design[["random_effects"]]) {
      if (!identical(term[["parameterization_resolved"]], "mean_centered")) next
      name <- .random_effect_term_block_name(term)
      source <- random_sources[[match(name, vapply(random_sources, `[[`, character(1L), "name"))]]
      if (!isTRUE(source[["retained"]])) {
        stop("Mean-centered parameterization is unavailable for integrated selection source '",
             name, "'.", call. = FALSE)
      }
    }
  } else if (!inherits(object, "brma.mv")) {
    random_sources <- list(list(
      name     = "estimate",
      role     = "estimate",
      retained = identical(common[["estimate_random_effects"]], "condition")
    ))
    if (.is_data_multilevel(data)) {
      random_sources[[2L]] <- list(
        name     = "cluster",
        role     = "other",
        retained = identical(common[["other_random_effects"]], "condition")
      )
    }
  }

  has_estimate <- any(vapply(random_sources, function(source) {

    identical(source[["role"]], "estimate")
  }, logical(1)))
  has_random_context <- any(vapply(random_sources, function(source) {

    identical(source[["role"]], "other")
  }, logical(1)))
  has_sampling_context <- any(data[["outcome"]][["sei"]] > 0)
  has_sampling_dependence <- FALSE
  if (.is_data_known_v(data)) {
    known_V <- .data_known_v_data(data)
    if (identical(common[["known_sampling_variance"]], "condition")) {
      known_V <- .known_v_resolve_selection_structure(known_V)
      attr(data, "known_V_data") <- known_V
    }
    has_sampling_dependence <- length(.known_v_correlated_blocks(known_V)) > 0L
  }

  binding <- attr(data, "selection_binding", exact = TRUE)
  if (is.null(binding)) {
    stop("Selection group binding data are unavailable.", call. = FALSE)
  }
  groups <- lapply(models[active], function(model) {

    .selection_bind_groups(
      model            = model,
      row_index        = binding[["row_index"]],
      input_data       = binding[["input_data"]],
      cluster          = binding[["cluster"]],
      allow_singletons = !.is_data_random(data) &&
        !has_random_context && !has_sampling_dependence
    )
  })
  best <- vapply(models[active], function(model) {

    identical(model[["weight_rule"]], "best")
  }, logical(1))
  publication_groups <- if (any(best)) groups[best] else groups[1L]
  partition <- publication_groups[[1L]][["group_index"]]
  if (any(vapply(publication_groups, function(group) {

    !identical(group[["group_index"]], partition)
  }, logical(1)))) {
    stop(
      "Active best-weight branches must use the same publication partition.",
      call. = FALSE
    )
  }

  model <- structure(list(
    schema_version          = 3L,
    estimate_random_effects = common[["estimate_random_effects"]],
    other_random_effects    = common[["other_random_effects"]],
    known_sampling_variance = common[["known_sampling_variance"]],
    branches         = models,
    active_branches  = which(active),
    groups           = publication_groups[[1L]],
    branch_groups    = groups,
    sources          = list(random = random_sources),
    applicability    = list(
      estimate_random_effects = has_estimate,
      other_random_effects    = has_random_context,
      known_sampling_variance = has_sampling_context
    )
  ), class = c("RoBMA_selection_model", "list"))
  attr(data, "selection_model")   <- model
  attr(data, "selection_binding") <- NULL
  object[["data"]] <- data
  object
}


.prepare_selection_likelihood_object <- function(
    object, selection_control) {

  selection_control <- .check_selection_likelihood_control(selection_control)

  if (!.is_priors_weightfunction(object[["priors"]])) {
    stop(
      "Selection-likelihood preparation requires a selection-kernel prior.",
      call. = FALSE
    )
  }
  model <- .data_selection_model(object[["data"]])
  if (is.null(model)) {
    return(object)
  }

  signed_yi <- object[["data"]][["outcome"]][["yi"]]
  if (identical(.data_effect_direction(object[["data"]]), "negative")) {
    signed_yi <- -signed_yi
  }
  selection_spec <- .selection_spec(
    priors           = object[["priors"]],
    yi               = signed_yi,
    sei              = object[["data"]][["outcome"]][["sei"]],
    effect_direction = .data_effect_direction(object[["data"]]),
    signed_data      = TRUE
  )
  if (is.null(selection_spec) || !identical(selection_spec[["mode"]], "step")) {
    stop(
      "Joint selection requires a step selection kernel. ",
      "P-hacking and other non-step kernels are unavailable.",
      call. = FALSE
    )
  }

  formula_design <- NULL
  if (.is_data_random(object[["data"]])) {
    formula_design <- object[["formula_design"]][["mu"]]
    if (!inherits(formula_design, "BayesTools_formula_design")) {
      stop(
        "The compiled selection random-effect metadata are unavailable.",
        call. = FALSE
      )
    }
    integrated_names <- vapply(Filter(function(source) {

      !source[["retained"]]
    }, model[["sources"]][["random"]]), `[[`, character(1), "name")
    formula_design[["random_effects"]] <- Filter(function(term) {

      .random_effect_term_block_name(term) %in% integrated_names
    }, formula_design[["random_effects"]])
    if (length(formula_design[["random_effects"]]) == 0L) {
      formula_design <- NULL
    }
  }

  if (.selection_retains_sampling(object[["data"]])) {
    return(.prepare_conditioned_sampling_likelihood(
      object, formula_design, selection_spec, selection_control
    ))
  }

  sampling <- .selection_joint_sampling_plan(object[["data"]])
  row_blocks <- .selection_joint_dependency_blocks(
    data         = object[["data"]],
    random_terms = if (is.null(formula_design)) {
      list()
    } else {
      formula_design[["random_effects"]]
    },
    sampling     = sampling
  )
  uses_best <- any(vapply(
    model[["branches"]][model[["active_branches"]]],
    function(branch) identical(branch[["weight_rule"]], "best"),
    logical(1)
  ))
  if (uses_best) {
    partition <- model[["groups"]][["group_index"]]
    for (rows in row_blocks) {
      if (length(unique(partition[rows])) > 1L) {
        cross_group <- outer(partition, partition, `!=`)
        sources <- character()
        if (!is.null(formula_design)) {
          for (term in formula_design[["random_effects"]]) {
            dependency <- BayesTools::random_effects_dependency_matrix(
              list(term), n_rows = length(partition)
            )
            if (any(dependency[cross_group])) {
              sources <- c(sources, paste0(
                "Random effects (from random): ", .random_effect_term_block_name(term)
              ))
            }
          }
        }
        if (.is_data_multilevel(object[["data"]]) &&
            !.selection_retains_other_random(object[["data"]])) {
          cluster <- object[["data"]][["outcome"]][["cluster"]]
          if (any(outer(cluster, cluster, `==`)[cross_group])) {
            sources <- c(sources, "Random effects: cluster")
          }
        }
        if (any(.selection_joint_sampling_block(
          sampling, seq_along(partition)
        )[cross_group] != 0) || length(sources) == 0L) {
          sources <- c(sources, "Sampling context (from V)")
        }
        stop(
          "'weight_rule = \"best\"' is unavailable because integrated ",
          "source(s) ", paste(sources, collapse = "; "),
          " connect publication groups ",
          paste(unique(partition[rows]), collapse = ", "),
          ".",
          call. = FALSE
        )
      }
    }
    row_blocks <- model[["groups"]][["row_blocks"]]
  }
  if (.is_data_weights(object[["data"]]) &&
      (uses_best || any(lengths(row_blocks) > 1L))) {
    stop(
      "'weights' are unavailable for jointly normalized selection vectors. ",
      "Omit 'weights'.",
      call. = FALSE
    )
  }
  sampling_factor_blocks <- if (!identical(
      sampling[["representation"]],
      "diagonal_factor"
    )) {
    NULL
  } else {
    .selection_joint_sampling_factor_blocks(sampling, row_blocks)
  }
  random_covariance <- NULL
  if (!is.null(formula_design)) {
    random_covariance <- BayesTools::JAGS_formula_random_marginal_covariance(
      formula_design = formula_design,
      row_blocks     = row_blocks,
      prefix         = "sel_joint_random",
      representation = if (!is.null(sampling_factor_blocks)) {
        "auto"
      } else "dense"
    )
  }
  block_routing <- .selection_joint_block_routing(
    data                   = object[["data"]],
    row_blocks             = row_blocks,
    sampling_factor_blocks = sampling_factor_blocks,
    random_covariance      = random_covariance
  )
  plan <- .selection_joint_execution_plan(
    row_blocks        = row_blocks,
    block_methods     = block_routing[["methods"]],
    factor_ranks      = block_routing[["ranks"]],
    selection_control = selection_control,
    sampling          = sampling,
    sampling_factor_blocks = sampling_factor_blocks,
    random_covariance = random_covariance
  )
  plan[["selection_spec"]] <- selection_spec
  attr(object[["data"]], "selection_execution_plan") <- plan
  object[["selection_control"]] <- selection_control

  object
}


.selection_joint_sampling_plan <- function(data) {

  if (.selection_retains_sampling(data)) {
    K <- nrow(data[["outcome"]])
    return(structure(list(
      representation = "diagonal_factor",
      diagonal       = numeric(K),
      loading        = matrix(0, nrow = K, ncol = 0L)
    ), class = c("RoBMA_selection_sampling_plan", "list")))
  }

  factor <- .selection_joint_sampling_factor(data)
  if (!is.null(factor)) {
    return(structure(
      c(list(representation = "diagonal_factor"), factor),
      class = c("RoBMA_selection_sampling_plan", "list")
    ))
  }

  structure(
    list(
      representation = "dense",
      covariance     = .known_v_covariance_matrix(.data_known_v_data(data))
    ),
    class = c("RoBMA_selection_sampling_plan", "list")
  )
}


.prepare_conditioned_sampling_likelihood <- function(
    object, formula_design, selection_spec, selection_control) {

  data  <- object[["data"]]
  model <- .data_selection_model(data)
  K     <- nrow(data[["outcome"]])
  if (.is_data_weights(data) && any(data[["outcome"]][["weights"]] != 1)) {
    stop("Non-unit 'weights' are unavailable with 'known_sampling_variance = \"condition\"'. ",
      "Omit 'weights' or set 'known_sampling_variance = \"integrate\"' in 'selection_model()'.",
      call. = FALSE)
  }
  full_sampling <- .selection_joint_sampling_factor(data)
  full_sampling <- if (is.null(full_sampling)) {
    list(representation = "dense",
         covariance = .known_v_covariance_matrix(.data_known_v_data(data)))
  } else c(list(representation = "diagonal_factor"), full_sampling)
  random_terms <- if (is.null(formula_design)) list() else
    formula_design[["random_effects"]]
  all_terms <- if (.is_data_random(data)) {
    object[["formula_design"]][["mu"]][["random_effects"]]
  } else list()
  retained_terms <- Filter(function(term) {
    !.random_effect_term_block_name(term) %in%
      vapply(random_terms, .random_effect_term_block_name, character(1))
  }, all_terms)
  row_blocks <- .selection_joint_dependency_blocks(data, all_terms, full_sampling)
  if (.is_data_multilevel(data) && !.is_data_random(data)) {
    adjacency <- outer(data[["outcome"]][["cluster"]], data[["outcome"]][["cluster"]], `==`)
    for (rows in row_blocks) adjacency[rows, rows] <- TRUE
    row_blocks <- .known_v_block_indices(adjacency * 1)
  }
  uses_best <- any(vapply(model[["branches"]][model[["active_branches"]]],
    function(branch) identical(branch[["weight_rule"]], "best"), logical(1)))
  if (uses_best) {
    candidate_blocks <- .selection_joint_dependency_blocks(
      data, random_terms, .selection_joint_sampling_plan(data)
    )
    partition <- model[["groups"]][["group_index"]]
    if (any(vapply(candidate_blocks, function(rows) {
      length(unique(partition[rows])) > 1L
    }, logical(1)))) {
      stop("'weight_rule = \"best\"' is unavailable because integrated random effects connect publication groups.",
           call. = FALSE)
    }
    adjacency <- diag(TRUE, K)
    for (rows in c(row_blocks, model[["groups"]][["row_blocks"]])) {
      adjacency[rows, rows] <- TRUE
    }
    row_blocks <- .known_v_block_indices(adjacency * 1)
  }
  random_covariance <- if (is.null(formula_design)) NULL else
    BayesTools::JAGS_formula_random_marginal_covariance(
      formula_design = formula_design, row_blocks = row_blocks,
      prefix = "sel_joint_random", representation = "auto"
    )
  retained_covariance <- NULL
  if (length(retained_terms)) {
    retained_design <- object[["formula_design"]][["mu"]]
    retained_design[["random_effects"]] <- retained_terms
    retained_covariance <- BayesTools::JAGS_formula_random_marginal_covariance(
      formula_design = retained_design, row_blocks = row_blocks,
      prefix = "sel_retained_random", representation = "auto"
    )
  }
  ranks <- if (is.null(random_covariance)) integer(length(row_blocks)) else
    as.integer(random_covariance[["loading_ranks"]])
  if (.is_data_multilevel(data) && !.is_data_random(data) &&
      !.selection_retains_other_random(data)) {
    ranks <- ranks + vapply(row_blocks, function(rows) {
      length(unique(data[["outcome"]][["cluster"]][rows]))
    }, integer(1))
  }
  plan <- .selection_joint_execution_plan(
    row_blocks = row_blocks,
    block_methods = rep("conditioned_sampling", length(row_blocks)),
    factor_ranks = ranks, selection_control = selection_control,
    sampling = .selection_joint_sampling_plan(data),
    sampling_factor_blocks = NULL, random_covariance = random_covariance
  )
  dimensions <- pmax(2L * lengths(row_blocks), ranks)
  keys       <- ifelse(ranks == 0L, "conditioned_analytic", paste0("conditioned_", dimensions))
  plan[["design_keys"]] <- keys
  plan[["designs"]] <- stats::setNames(lapply(unique(keys), function(key) {
    index <- match(key, keys)
    # A diagonal candidate law has an analytic normalizer. No QMC values are
    # used, but JAGS still requires a nonempty argument to the distribution.
    if (ranks[[index]] == 0L) return(array(0.5, c(1L, 1L, 1L)))
    BayesTools::selection_qmc_design(
      dimensions = dimensions[[index]],
      points = selection_control[["max_points_per_scramble"]],
      scrambles = selection_control[["scrambles"]], seed = selection_control[["seed"]]
    )
  }), unique(keys))
  plan[["statistical_target"]] <- "whole_sampling_error_selection"
  plan[["exactness"]] <- if (any(ranks > 4L)) "E2" else if (any(ranks > 1L)) {
    "EF"
  } else if (any(ranks == 1L)) "E1" else "E0"
  plan[["sampling_covariance"]] <- if (.is_data_known_v(data)) {
    .known_v_covariance_matrix(.data_known_v_data(data))
  } else diag(data[["outcome"]][["sei"]]^2, K)
  plan[["sampling_auxiliary"]] <- .selection_sampling_structure(data)
  plan[["integrated_parameter_stems"]] <- vapply(random_terms,
    `[[`, character(1), "parameter_stem")
  plan[["retained_parameter_stems"]] <- vapply(retained_terms,
    `[[`, character(1), "parameter_stem")
  plan[["retained_random_covariance"]] <- retained_covariance
  plan[["quadrature"]] <- .selection_joint_cluster_quadrature_rules(
    SELNORM_CLUSTER_QUADRATURE_ORDERS)
  plan[["factor_quadrature"]] <- .selection_joint_factor_quadrature_rules(4L)
  plan[["selection_spec"]] <- selection_spec
  attr(data, "selection_execution_plan") <- plan
  object[["data"]] <- data
  object[["selection_control"]] <- selection_control
  object
}


.selection_joint_sampling_diagonal <- function(sampling) {

  if (identical(sampling[["representation"]], "dense")) {
    return(diag(sampling[["covariance"]]))
  }
  if (!identical(sampling[["representation"]], "diagonal_factor")) {
    stop("Selection sampling-covariance metadata are invalid.",
         call. = FALSE)
  }

  sampling[["diagonal"]] + rowSums(sampling[["loading"]]^2)
}


.selection_conditioned_random_expression <- function(data, retained = FALSE) {

  plan  <- .data_selection_execution_plan(data)
  stems <- plan[[if (retained) "retained_parameter_stems" else "integrated_parameter_stems"]]
  terms <- if (length(stems)) paste0(stems, "[i]") else character()
  if (!.is_data_random(data)) {
    if (if (retained) .selection_retains_estimate(data) else .selection_integrates_estimate(data)) {
      tau <- if (.is_data_multilevel(data)) "tau_within" else "tau"
      if (.is_data_scale(data)) tau <- paste0(tau, "[i]")
      terms <- c(terms, paste0("theta[i] * ", tau))
    }
    if (.is_data_multilevel(data) && .selection_retains_other_random(data) == retained) {
      tau   <- if (.is_data_scale(data)) "tau_between[i]" else "tau_between"
      terms <- c(terms, paste0("gamma[cluster[i]] * ", tau))
    }
  }
  if (!length(terms)) return("0")
  expression <- paste(terms, collapse = " + ")
  if (identical(.data_effect_direction(data), "negative")) {
    expression <- paste0("-(", expression, ")")
  }
  expression
}


.selection_conditioned_sampling_fit_data <- function(data, priors) {

  plan         <- .data_selection_execution_plan(data)
  all_retained <- .selection_all_sources_conditioned(data)
  yi           <- data[["outcome"]][["yi"]]
  if (identical(.data_effect_direction(data), "negative")) yi <- -yi
  sei  <- data[["outcome"]][["sei"]]
  selection_data <- list()
  if (!all_retained) {
    spec <- .selection_spec(priors = priors, yi = yi, sei = sei,
      effect_direction = .data_effect_direction(data), signed_data = TRUE)
    selection_data <- spec[["jags_data"]]
    selection_data[["sel_obs_bin"]] <- NULL
  }
  out <- c(list(K = length(yi)), selection_data)
  if (.is_data_multilevel(data) && !all_retained) out[["cluster"]] <- data[["outcome"]][["cluster"]]
  if (.is_priors_PET(priors) || .is_priors_PEESE(priors)) out[["sei"]] <- sei
  if (!all_retained) {
    auxiliary   <- plan[["sampling_auxiliary"]]
    blocks      <- auxiliary[["latent_blocks"]]
    independent <- setdiff(seq_along(yi), unlist(lapply(blocks, `[[`, "index")))
    if (length(independent) && auxiliary[["rank"]] > 0L) {
      out[["known_v_independent_n"]] <- length(independent)
      out[["known_v_independent_index"]] <- independent
    }
    for (index in seq_along(blocks)) {
      out[[paste0("sampling_B_", index)]] <- blocks[[index]][["B"]]
    }
    for (key in names(plan[["designs"]])) {
      out[[.selection_joint_qmc_name(key)]] <- plan[["designs"]][[key]]
    }
    out[["sel_cond_cluster_nodes"]] <- plan[["quadrature"]][["nodes"]]
    out[["sel_cond_cluster_log_weights"]] <- plan[["quadrature"]][["log_weights"]]
    out[["sel_cond_cluster_orders"]] <- plan[["quadrature"]][["orders"]]
    out[["sel_cond_factor_nodes"]] <- plan[["factor_quadrature"]][["4"]][["nodes"]]
    out[["sel_cond_factor_log_weights"]] <- plan[["factor_quadrature"]][["4"]][["log_weights"]]
    out[["sel_cond_factor_orders"]] <- plan[["factor_quadrature"]][["4"]][["orders"]]
    out[["sel_cond_factor_rule_counts"]] <- plan[["factor_quadrature"]][["4"]][["rule_counts"]]
  }
  groups <- .data_selection_model(data)[["groups"]][["group_index"]]
  for (index in seq_along(plan[["row_blocks"]])) {
    rows       <- plan[["row_blocks"]][[index]]
    prefix     <- paste0("sel_joint_block_", index)
    covariance <- plan[["sampling_covariance"]][rows, rows, drop = FALSE]
    out[[paste0(prefix, "_y")]] <- yi[rows]
    if (!all_retained) {
      out[[paste0(prefix, "_sei")]] <- sei[rows]
      out[[paste0(prefix, "_obs_bin")]] <- spec[["obs_bin"]][rows]
    }
    out[[paste0(prefix, "_row")]] <- rows
    out[[paste0(prefix, "_sampling_lower")]] <-
      unname(covariance[lower.tri(covariance, diag = TRUE)])
    pairs               <- which(lower.tri(covariance, diag = TRUE), arr.ind = TRUE)
    retained_covariance <- plan[["retained_random_covariance"]]
    specialized_scale   <- !.is_data_random(data) && .is_data_scale(data)
    if (!is.null(retained_covariance) || (specialized_scale &&
        (.selection_retains_estimate(data) || .selection_retains_other_random(data)))) {
      out[[paste0(prefix, "_local_row_1")]] <- pairs[, 1L]
    }
    if ((!is.null(retained_covariance) && retained_covariance[["loading_ranks"]][[index]] > 0L) ||
        (specialized_scale && .is_data_multilevel(data) && .selection_retains_other_random(data))) {
      out[[paste0(prefix, "_local_row_2")]] <- pairs[, 2L]
    }
    if (!is.null(retained_covariance) || .selection_retains_estimate(data)) {
      out[[paste0(prefix, "_pair_diagonal")]] <- as.integer(pairs[, 1L] == pairs[, 2L])
    }
    if (.is_data_multilevel(data) && !.is_data_random(data) && .selection_retains_other_random(data)) {
      cluster <- data[["outcome"]][["cluster"]][rows]
      out[[paste0(prefix, "_same_cluster")]] <- as.integer(cluster[pairs[, 1L]] == cluster[pairs[, 2L]])
    }
    if (!all_retained) out[[paste0(prefix, "_group")]] <- match(groups[rows], unique(groups[rows]))
    if (.is_data_multilevel(data) && !.is_data_random(data) &&
        !.selection_retains_other_random(data)) {
      cluster <- data[["outcome"]][["cluster"]][rows]
      out[[paste0(prefix, "_cluster_loading")]] <-
        1L * outer(cluster, unique(cluster), `==`)
    }
  }
  if (!is.null(plan[["random_covariance"]])) {
    out <- c(out, plan[["random_covariance"]][["data"]])
  }
  if (!is.null(plan[["retained_random_covariance"]])) {
    out <- c(out, plan[["retained_random_covariance"]][["data"]])
  }
  out
}


.selection_conditioned_sampling_model_syntax <- function(data, selection_spec) {

  plan                <- .data_selection_execution_plan(data)
  all_retained        <- .selection_all_sources_conditioned(data)
  covariance          <- plan[["random_covariance"]]
  retained_covariance <- plan[["retained_random_covariance"]]
  syntax              <- if (is.null(covariance)) "" else covariance[["syntax"]]
  if (!is.null(retained_covariance)) syntax <- paste0(syntax, retained_covariance[["syntax"]])
  for (index in seq_along(plan[["row_blocks"]])) {
    rows     <- plan[["row_blocks"]][[index]]
    K        <- length(rows)
    rank     <- plan[["factor_ranks"]][[index]]
    prefix   <- paste0("sel_joint_block_", index)
    row      <- paste0(prefix, "_row[j]")
    diagonal <- "0"
    loading  <- character()
    if (!is.null(covariance)) {
      diagonal <- paste0(covariance[["diagonal_names"]][[index]], "[j]")
      loading  <- paste0(covariance[["loading_names"]][[index]],
        "[j,", seq_len(rank), "]")
      if (!rank) loading <- character()
    } else if (!.is_data_random(data)) {
      if (.selection_integrates_estimate(data)) {
        tau <- if (.is_data_multilevel(data)) "tau_within" else "tau"
        if (.is_data_scale(data)) tau <- paste0(tau, "[", row, "]")
        diagonal <- paste0("pow(", tau, ",2)")
      }
      if (.is_data_multilevel(data) && !.selection_retains_other_random(data)) {
        tau <- if (.is_data_scale(data)) paste0("tau_between[", row, "]") else
          "tau_between"
        loading <- paste0(tau, " * ", prefix, "_cluster_loading[j,",
          seq_len(rank), "]")
      }
    }
    syntax <- paste0(syntax, "for(j in 1:", K, "){\n",
      "  ", prefix, "_mu[j] = sel_joint_mu[", row, "]\n",
      if (!all_retained) paste0(
        "  ", prefix, "_e0[j] = sampling_dependency[", row, "] + sel_joint_retained[", row, "]\n",
        "  ", prefix, "_u0[j] = sel_joint_integrated[", row, "]\n") else "",
      if (!all_retained || K > 1L) paste0("  ", prefix, "_diagonal[j] = ", diagonal, "\n") else "",
      if (rank) paste0("  ", prefix, "_loading[j,", seq_len(rank), "] = ",
        loading, "\n", collapse = "") else "",
      "}\n")
    local_1        <- paste0(prefix, "_local_row_1[l]")
    local_2        <- paste0(prefix, "_local_row_2[l]")
    retained_parts <- paste0(prefix, "_sampling_lower[l]")
    if (!is.null(retained_covariance)) {
      retained_parts <- c(retained_parts, paste0(prefix, "_pair_diagonal[l] * ",
        retained_covariance[["diagonal_names"]][[index]], "[", local_1, "]"))
      retained_rank <- retained_covariance[["loading_ranks"]][[index]]
      if (retained_rank) retained_parts <- c(retained_parts, paste0("inprod(",
        retained_covariance[["loading_names"]][[index]], "[", local_1, ",1:", retained_rank, "],",
        retained_covariance[["loading_names"]][[index]], "[", local_2, ",1:", retained_rank, "])") )
    } else if (!.is_data_random(data)) {
      global_1 <- paste0(prefix, "_row[", local_1, "]")
      global_2 <- paste0(prefix, "_row[", local_2, "]")
      if (.selection_retains_estimate(data)) {
        tau <- if (.is_data_multilevel(data)) "tau_within" else "tau"
        if (.is_data_scale(data)) tau <- paste0(tau, "[", global_1, "]")
        retained_parts <- c(retained_parts, paste0(prefix, "_pair_diagonal[l] * pow(", tau, ",2)"))
      }
      if (.is_data_multilevel(data) && .selection_retains_other_random(data)) {
        between <- if (.is_data_scale(data)) paste0("tau_between[", global_1,
          "] * tau_between[", global_2, "]") else "pow(tau_between,2)"
        retained_parts <- c(retained_parts, paste0(prefix, "_same_cluster[l] * ", between))
      }
    }
    syntax <- paste0(syntax, "for(l in 1:", K * (K + 1L) / 2L, "){\n",
      "  ", prefix, "_retained_covariance[l] = ", paste(retained_parts, collapse = " + "), "\n}\n")
    if (all_retained) {
      # The source-support preflight has ruled out zero acceptance. Positive
      # weights then cancel, leaving the full observed Gaussian covariance.
      # Keep normalized auxiliary priors for latent reconstruction, but do not
      # make them or omega parents of the observed likelihood.
      syntax <- paste0(syntax, if (K == 1L) {
        paste0(prefix, "_y[1] ~ dnorm(", prefix, "_mu[1], 1/",
          prefix, "_retained_covariance[1])\n")
      } else {
        paste0(prefix, "_y[1:", K, "] ~ dknown_v_mnorm(",
          prefix, "_mu[1:", K, "],", prefix, "_diagonal[1:", K, "],",
          prefix, "_retained_covariance)\n")
      })
      next
    }
    loading_expression <- if (rank) paste0(prefix, "_loading[1:", K,
      ",1:", rank, "]") else "0"
    syntax <- paste0(syntax, prefix, "_y[1:", K,
      "] ~ dselnorm_sampling_conditioned(",
      prefix, "_mu[1:", K, "],", prefix, "_retained_covariance,",
      prefix, "_diagonal[1:", K, "],", loading_expression, ",", rank, ",",
      prefix, "_e0[1:", K, "],", prefix, "_u0[1:", K, "],",
      prefix, "_sei,", selection_spec[["jags_omega"]], ",",
      "sel_z_lower,sel_z_upper,", prefix, "_obs_bin,sel_sign,",
      .selection_joint_kernel_mode_expression(selection_spec), ",",
      "sel_telescope_probabilities,", selection_spec[["jags_vector_rule"]], ",",
      prefix, "_group,", .selection_joint_qmc_name(plan[["design_keys"]][[index]]), ",",
      plan[["points_per_scramble"]], ",", plan[["max_points_per_scramble"]], ",",
      plan[["scrambles"]], ",", format(plan[["relative_tolerance"]], scientific = FALSE),
      ",sel_cond_cluster_nodes,sel_cond_cluster_log_weights,sel_cond_cluster_orders,",
      "sel_cond_factor_nodes,sel_cond_factor_log_weights,sel_cond_factor_orders,",
      "sel_cond_factor_rule_counts)\n")
  }
  syntax
}


.selection_joint_sampling_block <- function(sampling, rows) {

  if (identical(sampling[["representation"]], "dense")) {
    return(sampling[["covariance"]][rows, rows, drop = FALSE])
  }
  if (!identical(sampling[["representation"]], "diagonal_factor")) {
    stop("Selection sampling-covariance metadata are invalid.",
         call. = FALSE)
  }

  covariance <- diag(
    sampling[["diagonal"]][rows],
    nrow = length(rows),
    ncol = length(rows)
  )
  loading <- sampling[["loading"]][rows, , drop = FALSE]
  if (ncol(loading) > 0L) {
    covariance <- covariance + tcrossprod(loading)
  }
  covariance
}


.selection_sampling_auxiliary <- function(fit, data, posterior_samples = NULL) {

  samples   <- .get_posterior_samples(fit, posterior_samples)
  structure <- .selection_sampling_structure(data)
  out       <- matrix(0, nrow(samples), nrow(data[["outcome"]]))
  if (!structure[["rank"]]) return(out)
  z <- .extract_indexed_parameter_samples(
    samples, "sampling_z", n_expected = structure[["rank"]]
  )
  for (block in structure[["latent_blocks"]]) {
    columns <- seq.int(block[["z_start"]], block[["z_end"]])
    out[, block[["index"]]] <- out[, block[["index"]], drop = FALSE] +
      z[, columns, drop = FALSE] %*% t(block[["B"]])
  }
  if (identical(.data_effect_direction(data), "negative")) out <- -out
  out
}


.selection_conditioned_sampling_factors <- function(setup) {

  data   <- setup[["data"]]
  plan   <- .data_selection_execution_plan(data)
  S      <- setup[["S"]]
  K      <- setup[["K"]]
  result <- setup[["selection_conditioned_factors"]]
  if (!is.null(result)) {
    if (nrow(result[["diagonal"]]) == 1L && S > 1L) {
      result[["diagonal"]] <- result[["diagonal"]][rep(1L, S), , drop = FALSE]
      result[["loadings"]] <- lapply(result[["loadings"]], function(loading) {
        loading[rep(1L, S), , , drop = FALSE]
      })
    }
    return(result)
  }
  if (!is.null(plan[["random_covariance"]])) {
    return(.selection_joint_random_factor_samples(setup))
  }
  diagonal <- matrix(0, S, K)
  if (!.is_data_random(data) && .selection_integrates_estimate(data)) {
    scale <- setup[["tau_within"]]
    if (nrow(scale) == 1L && S > 1L) scale <- scale[rep(1L, S), , drop = FALSE]
    diagonal <- diagonal + scale^2
  }
  loadings <- lapply(seq_along(plan[["row_blocks"]]), function(index) {
    rows    <- plan[["row_blocks"]][[index]]
    rank    <- plan[["factor_ranks"]][[index]]
    loading <- array(0, c(S, length(rows), rank))
    if (rank && .is_data_multilevel(data) && !.is_data_random(data) &&
        !.selection_retains_other_random(data)) {
      cluster <- data[["outcome"]][["cluster"]][rows]
      scale   <- setup[["tau_between"]]
      if (nrow(scale) == 1L && S > 1L) scale <- scale[rep(1L, S), , drop = FALSE]
      for (column in seq_len(rank)) {
        local <- which(cluster == unique(cluster)[[column]])
        loading[, local, column] <- scale[, rows[local], drop = FALSE]
      }
    }
    loading
  })
  list(diagonal = diagonal, loadings = loadings,
       ranks = as.integer(plan[["factor_ranks"]]))
}


.selection_conditioned_sampling_context <- function(setup) {

  .selection_context_from_parts(
    fit = setup[["fit"]], data = setup[["data"]], priors = setup[["priors"]],
    posterior_samples = setup[["posterior_samples"]],
    effect_direction = .data_effect_direction(setup[["data"]])
  )
}


.selection_conditioned_sampling_diagnostics <- function(result, plan) {

  mcse   <- max(result[["relative_mcse"]])
  change <- max(result[["relative_change"]])
  if (!is.finite(mcse) || !is.finite(change) ||
      max(mcse, change) > plan[["relative_tolerance"]]) {
    stop("Selection normalizer was rejected by diagnostics: relative Monte Carlo standard error was ",
      format(mcse), " and nested-design relative change was ", format(change),
      ". Increase 'max_points_per_scramble' or 'scrambles' in 'selection_control'.",
      call. = FALSE)
  }
  invisible(NULL)
}


# For C = V + sum(retained random covariances) and integrated covariance I,
# prior auxiliaries C0 and I0 give delta = solve(C + I, y - mu - C0 - I0).
# Every source j is reconstructed as source_j0 + covariance_j %*% delta.
# The joint likelihood factor is phi_(C+I)(y-mu) W(y) / A(mu + C0 + C delta).
# Backend sampling of an auxiliary does not change its selection source role.
.selection_conditioned_sampling_state <- function(setup) {

  data <- setup[["data"]]
  S    <- nrow(setup[["posterior_samples"]])
  K    <- nrow(data[["outcome"]])
  setup[["S"]] <- S
  setup[["K"]] <- K
  plan                       <- .data_selection_execution_plan(data)
  sources                    <- .data_selection_model(data)[["sources"]][["random"]]
  contributions              <- .selection_random_source_contributions(setup)
  baseline                   <- setup[["mu"]]
  retained                   <- matrix(0, S, K)
  retained_random_covariance <- array(0, c(S, K, K))
  integrated                 <- matrix(0, S, K)
  for (source in sources) {
    contribution <- contributions[[source[["name"]]]]
    if (source[["retained"]]) {
      retained                   <- retained + contribution
      retained_random_covariance <- retained_random_covariance +
        .selection_random_source_covariance(setup, source)
    } else {
      integrated <- integrated + contribution
    }
  }
  auxiliary <- .selection_sampling_auxiliary(setup[["fit"]], data,
    setup[["posterior_samples"]])
  factors          <- .selection_conditioned_sampling_factors(setup)
  context          <- .selection_conditioned_sampling_context(setup)
  covariance       <- array(0, c(S, K, K))
  e                <- correction <- matrix(0, S, K)
  block_log_lik    <- log_normalizer <- matrix(0, S, length(plan[["row_blocks"]]))
  diagnostics      <- vector("list", length(plan[["row_blocks"]]))
  groups           <- .data_selection_model(data)[["groups"]][["group_index"]]
  total_covariance <- retained_random_covariance
  for (index in seq_along(plan[["row_blocks"]])) {
    rows           <- plan[["row_blocks"]][[index]]
    loading        <- factors[["loadings"]][[index]]
    rank           <- dim(loading)[[3L]]
    V              <- plan[["sampling_covariance"]][rows, rows, drop = FALSE]
    retained_lower <- matrix(0, S, length(rows) * (length(rows) + 1L) / 2L)
    for (draw in seq_len(S)) {
      L <- matrix(loading[draw, , ], length(rows), rank)
      covariance[draw, rows, rows] <- diag(factors[["diagonal"]][draw, rows],
        length(rows)) + tcrossprod(L)
      C <- V + matrix(retained_random_covariance[draw, rows, rows], length(rows))
      retained_lower[draw, ] <- C[lower.tri(C, diag = TRUE)]
      total_covariance[draw, rows, rows] <- C +
        matrix(covariance[draw, rows, rows], length(rows))
    }
    result <- .Call("RoBMA_selnorm_sampling_conditioned_batch",
      as.numeric(data[["outcome"]][["yi"]][rows]), baseline[, rows, drop = FALSE],
      retained_lower,
      factors[["diagonal"]][, rows, drop = FALSE], matrix(loading, nrow = S),
      as.integer(rank), auxiliary[, rows, drop = FALSE] + retained[, rows, drop = FALSE],
      integrated[, rows, drop = FALSE],
      as.numeric(context[["sei"]][rows]), context[["omega"]],
      as.numeric(context[["z_lower"]]), as.numeric(context[["z_upper"]]),
      as.integer(context[["obs_bin"]][rows]), as.integer(context[["sign"]]),
      as.integer(context[["kernel_mode"]]), isTRUE(context[["telescope_probabilities"]]),
      as.integer(context[["vector_rule"]]), as.integer(match(groups[rows], unique(groups[rows]))),
      plan[["designs"]][[plan[["design_keys"]][[index]]]],
      as.integer(plan[["points_per_scramble"]]), as.integer(plan[["max_points_per_scramble"]]),
      as.integer(plan[["scrambles"]]), as.numeric(plan[["relative_tolerance"]]),
      as.numeric(plan[["quadrature"]][["nodes"]]),
      as.numeric(plan[["quadrature"]][["log_weights"]]),
      as.numeric(plan[["quadrature"]][["orders"]]),
      as.numeric(plan[["factor_quadrature"]][["4"]][["nodes"]]),
      as.numeric(plan[["factor_quadrature"]][["4"]][["log_weights"]]),
      as.numeric(plan[["factor_quadrature"]][["4"]][["orders"]]),
      as.numeric(plan[["factor_quadrature"]][["4"]][["rule_counts"]]), PACKAGE = "RoBMA")
    .selection_conditioned_sampling_diagnostics(result, plan)
    correction[, rows] <- result[["delta"]]
    for (draw in seq_len(S)) {
      e[draw, rows] <- auxiliary[draw, rows] + as.vector(V %*% result[["delta"]][draw, ])
      baseline[draw, rows] <- setup[["mu"]][draw, rows] + retained[draw, rows] +
        as.vector(matrix(retained_random_covariance[draw, rows, rows], length(rows)) %*%
          result[["delta"]][draw, ])
    }
    block_log_lik[, index] <- result[["log_lik"]]
    log_normalizer[, index] <- result[["log_normalizer"]]
    diagnostics[[index]] <- result[c("relative_mcse", "relative_change")]
  }
  list(e = e, correction = correction, integrated_covariance = covariance,
       baseline_mu = baseline, sampling_covariance = plan[["sampling_covariance"]],
       total_covariance = total_covariance,
       log_lik = rowSums(block_log_lik), block_log_lik = block_log_lik,
       log_normalizer = log_normalizer, candidate_factors = factors,
       diagnostics = diagnostics)
}


.selection_conditioned_sampling_normalizer <- function(setup, means) {

  means <- as.matrix(means)
  S     <- nrow(means)
  K     <- ncol(means)
  data  <- setup[["data"]]
  plan  <- .data_selection_execution_plan(data)
  setup[["S"]] <- S
  setup[["K"]] <- K
  if (nrow(setup[["posterior_samples"]]) == 1L && S > 1L) {
    setup[["posterior_samples"]] <- setup[["posterior_samples"]][rep(1L, S), , drop = FALSE]
  }
  factors     <- .selection_conditioned_sampling_factors(setup)
  context     <- .selection_conditioned_sampling_context(setup)
  groups      <- .data_selection_model(data)[["groups"]][["group_index"]]
  value       <- numeric(S)
  diagnostics <- vector("list", length(plan[["row_blocks"]]))
  for (index in seq_along(plan[["row_blocks"]])) {
    rows    <- plan[["row_blocks"]][[index]]
    loading <- factors[["loadings"]][[index]]
    rank    <- dim(loading)[[3L]]
    result  <- .Call("RoBMA_selnorm_conditioned_normalizer_batch",
      means[, rows, drop = FALSE], factors[["diagonal"]][, rows, drop = FALSE],
      matrix(loading, nrow = S), as.integer(rank), as.numeric(context[["sei"]][rows]),
      context[["omega"]], as.numeric(context[["z_lower"]]), as.numeric(context[["z_upper"]]),
      as.integer(context[["sign"]]), as.integer(context[["kernel_mode"]]),
      isTRUE(context[["telescope_probabilities"]]), as.integer(context[["vector_rule"]]),
      as.integer(match(groups[rows], unique(groups[rows]))),
      plan[["designs"]][[plan[["design_keys"]][[index]]]],
      as.integer(plan[["points_per_scramble"]]), as.integer(plan[["max_points_per_scramble"]]),
      as.integer(plan[["scrambles"]]), as.numeric(plan[["relative_tolerance"]]),
      as.numeric(plan[["quadrature"]][["nodes"]]),
      as.numeric(plan[["quadrature"]][["log_weights"]]),
      as.numeric(plan[["quadrature"]][["orders"]]),
      as.numeric(plan[["factor_quadrature"]][["4"]][["nodes"]]),
      as.numeric(plan[["factor_quadrature"]][["4"]][["log_weights"]]),
      as.numeric(plan[["factor_quadrature"]][["4"]][["orders"]]),
      as.numeric(plan[["factor_quadrature"]][["4"]][["rule_counts"]]), PACKAGE = "RoBMA")
    .selection_conditioned_sampling_diagnostics(result, plan)
    value <- value + result[["log_normalizer"]]
    diagnostics[[index]] <- result[c("relative_mcse", "relative_change")]
  }
  list(log_mass = value, diagnostics = diagnostics)
}


.selection_joint_sampling_factor <- function(data) {

  if (!.is_data_known_v(data)) {
    return(list(
      diagonal = data[["outcome"]][["sei"]]^2,
      loading  = matrix(
        0,
        nrow = length(data[["outcome"]][["sei"]]),
        ncol = 0L
      )
    ))
  }
  known_V <- .data_known_v_data(data)
  storage <- .known_v_storage(known_V)
  if (storage == "diagonal" ||
      (storage == "blocks" &&
       length(.known_v_correlated_blocks(known_V)) == 0L)) {
    return(list(
      diagonal = .known_v_diagonal(known_V),
      loading  = matrix(0, nrow = .known_v_nrow(known_V), ncol = 0L)
    ))
  }
  if (storage == "factor") {
    loading <- known_V[["factor_loading"]]
    blocks  <- known_V[["block_indices"]]
    .known_v_validate_dependency_blocks(blocks, .known_v_nrow(known_V))
    partition <- integer(.known_v_nrow(known_V))
    for (block in seq_along(blocks)) partition[blocks[[block]]] <- block
    for (column in seq_len(ncol(loading))) {
      if (length(unique(partition[loading[, column] != 0])) > 1L) {
        # Exact covariance cancellation can separate rows sharing auxiliary
        # columns. The existing dense plan preserves V across those blocks.
        return(NULL)
      }
    }
    return(list(
      diagonal = as.numeric(known_V[["factor_diagonal"]]),
      loading  = unname(loading)
    ))
  }
  NULL
}


.selection_joint_sampling_factor_blocks <- function(factor, row_blocks) {

  diagonal <- factor[["diagonal"]]
  loading  <- factor[["loading"]]
  K        <- length(diagonal)
  if (!is.numeric(diagonal) || anyNA(diagonal) || any(!is.finite(diagonal)) ||
      any(diagonal < 0) || !is.numeric(loading) || !is.matrix(loading) ||
      nrow(loading) != K || anyNA(loading) || any(!is.finite(loading))) {
    stop("Selection sampling-factor metadata are invalid.",
         call. = FALSE)
  }

  support <- lapply(seq_len(ncol(loading)), function(column) {
    which(loading[, column] != 0)
  })
  singleton <- which(lengths(support) == 1L)
  if (length(singleton) > 0L) {
    for (column in singleton) {
      row <- support[[column]][[1L]]
      diagonal[[row]] <- diagonal[[row]] + loading[row, column]^2
    }
  }
  dependent <- which(lengths(support) > 1L)

  lapply(row_blocks, function(rows) {
    columns <- dependent[vapply(
      support[dependent],
      function(column_rows) any(column_rows %in% rows),
      logical(1L)
    )]
    if (length(columns) > 0L && any(vapply(
      support[columns],
      function(column_rows) !all(column_rows %in% rows),
      logical(1L)
    ))) {
      stop(
        "Internal error: a certified sampling factor crosses selection blocks.",
        call. = FALSE
      )
    }
    list(
      diagonal = diagonal[rows],
      loading  = loading[rows, columns, drop = FALSE],
      rank     = length(columns)
    )
  })
}


.selection_joint_dependency_blocks <- function(
    data, random_terms = list(), sampling = NULL) {

  if (!is.list(random_terms)) {
    stop("Random-effect dependency metadata must be a list.",
         call. = FALSE)
  }
  if (is.null(sampling)) {
    sampling <- .selection_joint_sampling_plan(data)
  }
  structural_sampling <- if (identical(
      sampling[["representation"]],
      "dense"
    )) {
    sampling[["covariance"]]
  } else {
    K <- length(sampling[["diagonal"]])
    adjacency <- diag(TRUE, K)
    loading <- sampling[["loading"]]
    for (column in seq_len(ncol(loading))) {
      support <- loading[, column] != 0
      adjacency <- adjacency | outer(support, support, "&")
    }
    adjacency * 1
  }

  if (length(random_terms) > 0L) {
    block_names <- unique(vapply(
      random_terms,
      .random_effect_term_block_name,
      character(1)
    ))
    return(.random_effect_dependency_blocks(
      sampling_covariance = structural_sampling,
      formula_design      = list(random_effects = random_terms),
      blocks              = block_names
    ))
  }

  adjacency <- structural_sampling != 0

  if (.is_data_multilevel(data) && !.selection_retains_other_random(data)) {
    cluster   <- data[["outcome"]][["cluster"]]
    adjacency <- adjacency | outer(cluster, cluster, "==")
  }
  diag(adjacency) <- TRUE

  .known_v_block_indices(adjacency * 1)
}


.selection_joint_lower_pairs <- function(plan, rows) {

  pairs <- plan[["lower_pairs"]][[as.character(length(rows))]]
  if (is.null(pairs)) {
    row_1 <- unlist(lapply(seq_along(rows), function(column) {
      column:length(rows)
    }), use.names = FALSE)
    row_2 <- rep.int(seq_along(rows), rev(seq_along(rows)))
    pairs <- data.frame(row_1 = row_1, row_2 = row_2)
  }
  if (!is.data.frame(pairs) ||
      !identical(names(pairs), c("row_1", "row_2")) ||
      nrow(pairs) != length(rows) * (length(rows) + 1L) / 2L ||
      anyNA(pairs) ||
      any(unlist(pairs, use.names = FALSE) < 1L) ||
      any(unlist(pairs, use.names = FALSE) > length(rows))) {
    stop(
      "Internal error: selection covariance ordering is invalid.",
      call. = FALSE
    )
  }
  list(
    row_1 = rows[pairs[["row_1"]]],
    row_2 = rows[pairs[["row_2"]]]
  )
}


.selection_joint_uses_diagonal_indicator <- function(
    data, method, random_covariance) {

  method == "dense" &&
    (!.is_data_random(data) ||
      identical(random_covariance[["representation"]], "diagonal_factor"))
}


.selection_joint_fit_data <- function(data, priors) {

  if (.selection_retains_sampling(data)) {
    return(.selection_conditioned_sampling_fit_data(data, priors))
  }
  plan             <- .data_selection_execution_plan(data)
  yi               <- data[["outcome"]][["yi"]]
  effect_direction <- .data_effect_direction(data)
  if (identical(effect_direction, "negative")) {
    yi <- -yi
  }
  sei <- data[["outcome"]][["sei"]]
  selection_spec <- .selection_spec(
    priors           = priors,
    yi               = yi,
    sei              = sei,
    effect_direction = effect_direction,
    signed_data      = TRUE
  )
  selection_data <- selection_spec[["jags_data"]]
  selection_data[["sel_obs_bin"]] <- NULL
  fit_data <- c(list(K = length(yi)), selection_data)
  if (.is_data_multilevel(data) && .selection_retains_other_random(data)) {
    fit_data[["cluster"]] <- data[["outcome"]][["cluster"]]
  }
  if (.is_priors_PET(priors) || .is_priors_PEESE(priors)) {
    fit_data[["sei"]] <- sei
  }

  block_methods <- plan[["block_methods"]]
  if (any(block_methods %in% c("rank_one", "dense"))) {
    quadrature <- plan[["quadrature"]]
    fit_data[["sel_joint_cluster_nodes"]] <- quadrature[["nodes"]]
    fit_data[["sel_joint_cluster_log_weights"]] <-
      quadrature[["log_weights"]]
    fit_data[["sel_joint_cluster_orders"]] <- quadrature[["orders"]]
  }
  for (rank_name in names(plan[["factor_quadrature"]])) {
    quadrature <- plan[["factor_quadrature"]][[rank_name]]
    fit_data[[.selection_joint_factor_quadrature_name(
      "nodes", rank_name
    )]] <- quadrature[["nodes"]]
    fit_data[[.selection_joint_factor_quadrature_name(
      "log_weights", rank_name
    )]] <- quadrature[["log_weights"]]
    fit_data[[.selection_joint_factor_quadrature_name(
      "orders", rank_name
    )]] <- quadrature[["orders"]]
    fit_data[[.selection_joint_factor_quadrature_name(
      "rule_counts", rank_name
    )]] <- quadrature[["rule_counts"]]
  }
  for (design_key in names(plan[["designs"]])) {
    fit_data[[.selection_joint_qmc_name(design_key)]] <-
      plan[["designs"]][[design_key]]
  }

  if (length(plan[["singleton_blocks"]]) > 0L) {
    rows <- plan[["singleton_rows"]]
    fit_data[["sel_joint_singleton_n"]] <- length(rows)
    fit_data[["sel_joint_singleton_y"]] <- yi[rows]
    fit_data[["sel_joint_singleton_sei"]] <- sei[rows]
    fit_data[["sel_joint_singleton_weights"]] <- if (.is_data_weights(data)) {
      data[["outcome"]][["weights"]][rows]
    } else rep(1, length(rows))
    fit_data[["sel_joint_singleton_obs_bin"]] <-
      selection_spec[["obs_bin"]][rows]
    fit_data[["sel_joint_singleton_row"]] <- rows
    fit_data[["sel_joint_singleton_sampling_variance"]] <-
      .selection_joint_sampling_diagonal(plan[["sampling"]])[rows]
  }

  for (block_index in plan[["dependent_blocks"]]) {
    rows   <- plan[["row_blocks"]][[block_index]]
    prefix <- paste0("sel_joint_block_", block_index)
    method <- block_methods[[block_index]]
    fit_data[[paste0(prefix, "_y")]]       <- yi[rows]
    fit_data[[paste0(prefix, "_sei")]]     <- sei[rows]
    fit_data[[paste0(prefix, "_obs_bin")]] <- selection_spec[["obs_bin"]][rows]
    fit_data[[paste0(prefix, "_row")]]     <- rows
    uses_diagonal <- .selection_joint_uses_diagonal_indicator(
      data              = data,
      method            = method,
      random_covariance = plan[["random_covariance"]]
    )
    if (method %in% c("rank_one", "factor")) {
      sampling_factor <- plan[["sampling_factor_blocks"]][[block_index]]
      fit_data[[paste0(prefix, "_sampling_variance")]] <-
        sampling_factor[["diagonal"]]
      if (sampling_factor[["rank"]] > 0L) {
        fit_data[[paste0(prefix, "_sampling_loading")]] <-
          unname(sampling_factor[["loading"]])
      }
      if (.is_data_multilevel(data) && !.selection_retains_other_random(data)) {
        cluster <- data[["outcome"]][["cluster"]][rows]
        fit_data[[paste0(prefix, "_cluster_loading")]] <-
          1L * outer(cluster, unique(cluster), `==`)
      }
    }
    if (method == "dense") {
      pairs <- .selection_joint_lower_pairs(plan, rows)
      if (!is.null(plan[["random_covariance"]]) && uses_diagonal) {
        fit_data[[paste0(prefix, "_local_row_1")]] <-
          match(pairs[["row_1"]], rows)
        if (plan[["random_covariance"]][["loading_ranks"]][[block_index]] > 0L) {
          fit_data[[paste0(prefix, "_local_row_2")]] <-
            match(pairs[["row_2"]], rows)
        }
      }
    }
    if (uses_diagonal) {
      # The dense covariance syntax references the pair-diagonal indicator only
      # for a compiled random covariance or an integrated specialized estimate
      # variance. Emitting it otherwise leaves unused JAGS data behind.
      if (!is.null(plan[["random_covariance"]]) ||
          (!.is_data_random(data) && .selection_integrates_estimate(data))) {
        fit_data[[paste0(prefix, "_diagonal")]] <- as.integer(
          pairs[["row_1"]] == pairs[["row_2"]]
        )
      }
      if (!.is_data_random(data) && .is_data_scale(data)) {
        fit_data[[paste0(prefix, "_row_1")]] <- pairs[["row_1"]]
        if (.is_data_multilevel(data) && !.selection_retains_other_random(data)) {
          fit_data[[paste0(prefix, "_row_2")]] <- pairs[["row_2"]]
        }
      }
      if (.is_data_multilevel(data) && !.selection_retains_other_random(data)) {
        cluster <- data[["outcome"]][["cluster"]]
        fit_data[[paste0(prefix, "_same_cluster")]] <- as.integer(
          cluster[pairs[["row_1"]]] == cluster[pairs[["row_2"]]]
        )
      }
    }
    if (method == "dense") {
      sampling_covariance <- .selection_joint_sampling_block(
        plan[["sampling"]],
        rows
      )
      fit_data[[paste0(prefix, "_sampling_lower")]] <- unname(
        sampling_covariance[cbind(
          match(pairs[["row_1"]], rows),
          match(pairs[["row_2"]], rows)
        )]
      )
    }
  }
  if (!is.null(plan[["random_covariance"]])) {
    fit_data <- c(fit_data, plan[["random_covariance"]][["data"]])
  }

  fit_data
}


.selection_joint_block_routing <- function(
    data, row_blocks, sampling_factor_blocks, random_covariance) {

  block_sizes <- lengths(row_blocks)
  methods <- rep("dense", length(row_blocks))
  ranks   <- rep(NA_integer_, length(row_blocks))
  methods[block_sizes == 1L] <- "singleton"
  ranks[block_sizes == 1L]   <- 0L

  if (is.null(sampling_factor_blocks)) {
    return(list(methods = methods, ranks = ranks))
  }
  random_factor <- is.null(random_covariance) || identical(
    random_covariance[["representation"]],
    "diagonal_factor"
  )
  if (!random_factor) {
    return(list(methods = methods, ranks = ranks))
  }

  for (block_index in which(block_sizes > 1L)) {
    sampling <- sampling_factor_blocks[[block_index]]
    if (any(sampling[["diagonal"]] <= 0)) {
      next
    }
    random_rank <- if (is.null(random_covariance)) {
      0L
    } else {
      random_covariance[["loading_ranks"]][[block_index]]
    }
    specialized_rank <- if (.is_data_multilevel(data) &&
        !.is_data_random(data) && !.selection_retains_other_random(data)) {
      length(unique(data[["outcome"]][["cluster"]][row_blocks[[block_index]]]))
    } else 0L
    rank <- sampling[["rank"]] + random_rank + specialized_rank
    ranks[[block_index]] <- rank
    methods[[block_index]] <- if (rank == 1L) {
      "rank_one"
    } else if (rank >= 2L && rank <= 4L) {
      "factor"
    } else {
      "dense"
    }
  }
  list(methods = methods, ranks = ranks)
}


.selection_joint_execution_plan <- function(
    row_blocks, block_methods, factor_ranks, selection_control, sampling,
    sampling_factor_blocks, random_covariance) {

  block_sizes <- lengths(row_blocks)
  if (length(block_sizes) != length(block_methods) ||
      length(block_sizes) != length(factor_ranks)) {
    stop("Internal error: selection integration routing is invalid.",
         call. = FALSE)
  }
  design_keys <- rep(NA_character_, length(block_sizes))
  dense <- which(block_methods == "dense" & block_sizes > 1L)
  factor <- which(block_methods %in% c("rank_one", "factor"))
  factor_points_per_proposal <- as.integer(ceiling(
    selection_control[["points_per_scramble"]] / 2
  ))
  factor_max_points_per_proposal <- as.integer(ceiling(
    selection_control[["max_points_per_scramble"]] / 2
  ))
  design_keys[dense]  <- as.character(block_sizes[dense])
  design_keys[factor] <- paste0("factor_", factor_ranks[factor])
  unique_keys <- unique(design_keys[!is.na(design_keys)])
  designs <- stats::setNames(lapply(unique_keys, function(key) {
    block_index <- which(design_keys == key)[[1L]]
    dimensions <- if (block_methods[[block_index]] %in% c("rank_one", "factor")) {
      2L * factor_ranks[[block_index]]
    } else {
      2L * block_sizes[[block_index]]
    }
    points <- if (block_methods[[block_index]] %in% c("rank_one", "factor")) {
      factor_max_points_per_proposal
    } else {
      selection_control[["points_per_scramble"]]
    }
    BayesTools::selection_qmc_design(
      dimensions = dimensions,
      points     = points,
      scrambles  = selection_control[["scrambles"]],
      seed       = selection_control[["seed"]]
    )
  }), unique_keys)
  pair_sizes <- sort(unique(block_sizes[block_methods == "dense"]))
  lower_pairs <- stats::setNames(lapply(pair_sizes, function(block_n) {
    row_1 <- unlist(lapply(seq_len(block_n), function(column) {
      column:block_n
    }), use.names = FALSE)
    row_2 <- rep.int(seq_len(block_n), rev(seq_len(block_n)))
    data.frame(row_1 = row_1, row_2 = row_2)
  }), as.character(pair_sizes))

  exactness <- if (any(block_methods == "dense")) {
    "E2"
  } else if (any(block_methods == "factor")) {
    "EF"
  } else if (any(block_methods == "rank_one")) {
    "E1"
  } else {
    "E0"
  }
  structure(list(
    schema_version      = 5L,
    statistical_target  = "conditional_gaussian_vector_selection",
    row_blocks          = row_blocks,
    block_sizes         = as.integer(block_sizes),
    block_methods       = block_methods,
    factor_ranks        = as.integer(factor_ranks),
    singleton_blocks    = which(block_sizes == 1L),
    singleton_rows      = as.integer(vapply(
      row_blocks[block_sizes == 1L],
      `[[`,
      integer(1L),
      1L
    )),
    dependent_blocks    = which(block_sizes > 1L),
    sampling            = sampling,
    sampling_factor_blocks = sampling_factor_blocks,
    random_covariance   = random_covariance,
    design_keys         = design_keys,
    designs             = designs,
    points_per_scramble = selection_control[["points_per_scramble"]],
    max_points_per_scramble =
      selection_control[["max_points_per_scramble"]],
    factor_points_per_proposal = factor_points_per_proposal,
    factor_max_points_per_proposal = factor_max_points_per_proposal,
    scrambles            = selection_control[["scrambles"]],
    seed                 = selection_control[["seed"]],
    relative_tolerance   = selection_control[["relative_tolerance"]],
    quadrature          = if (any(block_methods %in% c("rank_one", "dense"))) {
      orders <- SELNORM_CLUSTER_QUADRATURE_ORDERS
      if (any(block_methods == "dense") && !any(block_methods == "rank_one")) {
        orders <- c(7L, orders)
      }
      .selection_joint_cluster_quadrature_rules(orders)
    } else NULL,
    factor_quadrature   = .selection_joint_factor_quadrature_rules(
      unique(factor_ranks[factor])
    ),
    lower_pairs         = lower_pairs,
    exactness           = exactness
  ), class = c("RoBMA_selection_execution_plan", "list"))
}


.selection_joint_cluster_quadrature_rules <- function(orders) {

  rules  <- lapply(orders, .gauss_hermite_nodes)
  list(
    orders = orders,
    nodes = unlist(lapply(rules, `[[`, "nodes"), use.names = FALSE),
    log_weights = unlist(
      lapply(rules, `[[`, "log_weights"),
      use.names = FALSE
    )
  )
}


.selection_joint_factor_quadrature_rules <- function(ranks) {

  rank_names <- intersect(names(SELNORM_FACTOR_QUADRATURE_ORDERS),
                          as.character(ranks))
  if (length(rank_names) == 0L) return(list())
  shared_orders <- SELNORM_FACTOR_QUADRATURE_ORDERS[["2"]]
  rule_counts   <- lengths(SELNORM_FACTOR_QUADRATURE_ORDERS)
  for (rank in names(rule_counts)) {
    if (!identical(SELNORM_FACTOR_QUADRATURE_ORDERS[[rank]],
                   utils::head(shared_orders, rule_counts[[rank]]))) {
      stop("Internal error: factor quadrature rules must share a prefix.",
           call. = FALSE)
    }
  }
  quadrature <- .selection_joint_cluster_quadrature_rules(shared_orders)
  stats::setNames(lapply(rank_names, function(rank) {
    c(quadrature, list(rule_counts = unname(rule_counts[
      as.character(seq.int(2L, as.integer(rank)))
    ])))
  }), rank_names)
}


.selection_joint_qmc_name <- function(block_n) {

  paste0("sel_joint_qmc_", block_n)
}


.selection_joint_factor_quadrature_name <- function(quantity, rank) {

  paste0("sel_joint_factor_", quantity, "_", rank)
}


.selection_joint_kernel_mode_expression <- function(selection_spec) {

  if (isTRUE(selection_spec[["jags_use_step_switch"]])) {
    return(selection_spec[["jags_kernel_mode"]])
  }
  as.character(SELKERNEL_STEP)
}


.selection_joint_signed_context <- function(setup, signed_yi) {

  context <- .selection_context_from_parts(
    fit               = setup[["fit"]],
    data              = setup[["data"]],
    priors            = setup[["priors"]],
    posterior_samples = setup[["posterior_samples"]],
    effect_direction  = setup[["effect_direction"]]
  )
  context[["sign"]]    <- 1L
  context[["yi"]]      <- signed_yi
  context[["obs_bin"]] <- .selection_obs_bin(
    yi     = signed_yi,
    sei    = setup[["selection_sei"]],
    p_cuts = context[["p_cuts"]],
    sign   = 1L
  )
  context[["jags_data"]][["sel_sign"]]    <- 1L
  context[["jags_data"]][["sel_obs_bin"]] <- context[["obs_bin"]]

  .selection_reset_native_cache(context)
}


.selection_joint_random_covariance_samples <- function(setup) {

  if (!.is_data_random(setup[["data"]])) {
    return(NULL)
  }
  plan <- .data_selection_execution_plan(setup[["data"]])
  if (is.null(plan[["random_covariance"]])) {
    return(NULL)
  }
  object <- list(
    fit    = setup[["fit"]],
    data   = setup[["data"]],
    priors = setup[["priors"]]
  )
  random_vcov <- .brma_mv_random_effects_marginal_vcov(
    object            = object,
    posterior_samples = setup[["posterior_samples"]],
    blocks            = plan[["random_covariance"]][["term_names"]]
  )
  samples <- random_vcov[["samples"]]
  if (!identical(dim(samples), c(setup[["S"]], setup[["K"]], setup[["K"]])) ||
      any(!is.finite(samples))) {
    stop(
      "Selection random-effect covariance samples are invalid.",
      call. = FALSE
    )
  }

  samples
}


.selection_joint_random_factor_samples <- function(setup, inputs = NULL) {

  if (!.is_data_random(setup[["data"]])) {
    return(NULL)
  }
  execution_plan <- .data_selection_execution_plan(setup[["data"]])
  random_covariance <- execution_plan[["random_covariance"]]
  if (is.null(random_covariance)) {
    return(NULL)
  }
  blocks <- random_covariance[["term_names"]]
  if (!is.character(blocks) || !length(blocks) || anyNA(blocks) ||
      any(!nzchar(blocks)) || anyDuplicated(blocks)) {
    stop("Selection random-effect source metadata are invalid.", call. = FALSE)
  }
  object <- list(
    fit    = setup[["fit"]],
    data   = setup[["data"]],
    priors = setup[["priors"]]
  )
  result <- setup[["selection_random_factor_samples"]]
  if (is.null(result)) {
    # The JAGS syntax representation does not determine post-fit covariance
    # eligibility. Retain its declared source names/row partition, and let
    # BayesTools compile the corresponding exact post-fit factor geometry.
    arguments <- list(
      object = object, posterior_samples = setup[["posterior_samples"]],
      blocks = blocks, row_blocks = execution_plan[["row_blocks"]]
    )
    if (!is.null(inputs)) arguments[["inputs"]] <- inputs
    factors <- do.call(.brma_mv_random_effects_marginal_factor_states, arguments)
    result <- tryCatch(
      BayesTools::random_effects_marginal_diagonal_factor(factors),
      BayesTools_random_effects_marginal_factor_unavailable = function(e) NULL
    )
    if (is.null(result)) {
      return(NULL)
    }
  }
  expected <- c(setup[["S"]], setup[["K"]])
  ranks <- if (is.list(result)) result[["ranks"]] else NULL
  if (!is.list(result) || !identical(dim(result[["diagonal"]]), expected) ||
      any(!is.finite(result[["diagonal"]])) ||
      any(result[["diagonal"]] < 0) ||
      !identical(result[["row_blocks"]], execution_plan[["row_blocks"]]) ||
      length(result[["loadings"]]) != length(execution_plan[["row_blocks"]]) ||
      length(result[["loading_supports"]]) != length(result[["loadings"]]) ||
      !is.integer(ranks) || length(ranks) != length(result[["loadings"]]) ||
      anyNA(ranks) || any(ranks < 0L)) {
    stop("Selection random-effect factors are invalid.",
         call. = FALSE)
  }
  for (block_index in seq_along(result[["loadings"]])) {
    expected_dim <- c(
      setup[["S"]],
      length(execution_plan[["row_blocks"]][[block_index]]),
      result[["ranks"]][[block_index]]
    )
    if (!identical(dim(result[["loadings"]][[block_index]]), expected_dim) ||
        any(!is.finite(result[["loadings"]][[block_index]]))) {
      stop("Selection random-effect factor loadings are invalid.",
           call. = FALSE)
    }
    support <- result[["loading_supports"]][[block_index]]
    if (!is.logical(support) || anyNA(support) ||
        !identical(dim(support), expected_dim[-1L])) {
      stop("Selection random-effect loading supports are invalid.",
           call. = FALSE)
    }
  }
  result
}


.selection_joint_singleton_variances <- function(
    setup, rows, block_indices, random_covariance_samples = NULL,
    random_factor_samples = NULL) {

  execution_plan <- .data_selection_execution_plan(setup[["data"]])
  S           <- setup[["S"]]
  variances   <- matrix(
    .selection_joint_sampling_diagonal(execution_plan[["sampling"]])[rows],
    nrow = S,
    ncol = length(rows),
    byrow = TRUE
  )

  if (!.is_data_random(setup[["data"]])) {
    if (.selection_integrates_estimate(setup[["data"]])) {
      variances <- variances + setup[["tau_within"]][, rows, drop = FALSE]^2
    }
    if (isTRUE(setup[["is_multilevel"]]) &&
        !.selection_retains_other_random(setup[["data"]])) {
      variances <- variances + setup[["tau_between"]][, rows, drop = FALSE]^2
    }
    return(variances)
  }

  if (!is.null(random_factor_samples)) {
    variances <- variances +
      random_factor_samples[["diagonal"]][, rows, drop = FALSE]
    for (index in seq_along(rows)) {
      block_index <- block_indices[[index]]
      if (random_factor_samples[["ranks"]][[block_index]] > 0L) {
        loading <- matrix(
          random_factor_samples[["loadings"]][[block_index]],
          nrow = S
        )
        variances[, index] <- variances[, index] + rowSums(loading^2)
      }
    }
    return(variances)
  }

  if (is.null(random_covariance_samples)) {
    if (is.null(execution_plan[["random_covariance"]])) {
      return(variances)
    }
    stop(
      "Selection random-effect covariance samples are missing.",
      call. = FALSE
    )
  }
  array_indices <- cbind(
    rep(seq_len(S), times = length(rows)),
    rep(rows, each = S),
    rep(rows, each = S)
  )
  variances + matrix(
    random_covariance_samples[array_indices],
    nrow = S,
    ncol = length(rows)
  )
}


.selection_joint_factor_block_samples <- function(
    setup, block_index, random_factor_samples = NULL) {

  execution_plan <- .data_selection_execution_plan(setup[["data"]])
  rows        <- execution_plan[["row_blocks"]][[block_index]]
  sampling    <- execution_plan[["sampling_factor_blocks"]][[block_index]]
  S           <- setup[["S"]]
  block_n     <- length(rows)
  residual_variance <- matrix(
    sampling[["diagonal"]],
    nrow = S,
    ncol = block_n,
    byrow = TRUE
  )
  loading_parts <- list()
  support_parts <- list()
  if (sampling[["rank"]] > 0L) {
    part <- array(0, dim = c(S, block_n, sampling[["rank"]]))
    for (factor in seq_len(sampling[["rank"]])) {
      part[, , factor] <- matrix(
        sampling[["loading"]][, factor],
        nrow = S,
        ncol = block_n,
        byrow = TRUE
      )
    }
    loading_parts[[length(loading_parts) + 1L]] <- part
    support_parts[[length(support_parts) + 1L]] <- sampling[["loading"]] != 0
  }

  if (.is_data_random(setup[["data"]])) {
    if (is.null(random_factor_samples) &&
        !is.null(execution_plan[["random_covariance"]])) {
      stop("Selection random-effect factors are missing.",
           call. = FALSE)
    }
    if (!is.null(random_factor_samples)) {
      residual_variance <- residual_variance +
        random_factor_samples[["diagonal"]][, rows, drop = FALSE]
    }
    if (!is.null(random_factor_samples) &&
        random_factor_samples[["ranks"]][[block_index]] > 0L) {
      loading_parts[[length(loading_parts) + 1L]] <-
        random_factor_samples[["loadings"]][[block_index]]
      support_parts[[length(support_parts) + 1L]] <-
        random_factor_samples[["loading_supports"]][[block_index]]
    }
  } else {
    if (.selection_integrates_estimate(setup[["data"]])) {
      residual_variance <- residual_variance +
        setup[["tau_within"]][, rows, drop = FALSE]^2
    }
    if (isTRUE(setup[["is_multilevel"]]) &&
        !.selection_retains_other_random(setup[["data"]])) {
      cluster <- setup[["data"]][["outcome"]][["cluster"]][rows]
      support <- outer(cluster, unique(cluster), `==`)
      part <- array(0, dim = c(S, block_n, ncol(support)))
      for (j in seq_len(ncol(support))) {
        part[, , j] <- sweep(
          setup[["tau_between"]][, rows, drop = FALSE], 2L, support[, j], `*`
        )
      }
      loading_parts[[length(loading_parts) + 1L]] <- part
      support_parts[[length(support_parts) + 1L]] <- support
    }
  }
  if (any(!is.finite(residual_variance)) || any(residual_variance <= 0)) {
    stop(
      "Selection factor residual variances must be finite and positive.",
      call. = FALSE
    )
  }
  loading <- if (length(loading_parts) == 0L) {
    matrix(numeric(), nrow = S, ncol = 0L)
  } else {
    part_ranks <- vapply(
      loading_parts,
      function(part) dim(part)[[3L]],
      integer(1L)
    )
    combined <- array(0, dim = c(S, block_n, sum(part_ranks)))
    start <- 1L
    for (part_index in seq_along(loading_parts)) {
      columns <- seq.int(start, length.out = part_ranks[[part_index]])
      combined[, , columns] <- loading_parts[[part_index]]
      start <- start + part_ranks[[part_index]]
    }
    matrix(combined, nrow = S)
  }
  list(
    residual_sd = sqrt(residual_variance),
    loading     = loading,
    loading_support = if (length(support_parts)) do.call(cbind, support_parts) else
      matrix(logical(), block_n, 0L),
    rank        = if (block_n == 0L) 0L else ncol(loading) %/% block_n
  )
}


.selection_joint_covariance_lower <- function(
    setup, block_index, random_covariance_samples = NULL,
    random_factor_samples = NULL) {

  execution_plan <- .data_selection_execution_plan(setup[["data"]])
  rows        <- execution_plan[["row_blocks"]][[block_index]]
  sampling    <- .selection_joint_sampling_block(
    execution_plan[["sampling"]],
    rows
  )
  block_n     <- length(rows)
  pairs       <- .selection_joint_lower_pairs(
    execution_plan,
    seq_len(block_n)
  )
  lower_n     <- length(pairs[["row_1"]])
  lower       <- matrix(
    sampling[cbind(pairs[["row_1"]], pairs[["row_2"]])],
    nrow = setup[["S"]],
    ncol = lower_n,
    byrow = TRUE
  )
  global_row_1 <- rows[pairs[["row_1"]]]
  global_row_2 <- rows[pairs[["row_2"]]]

  if (.is_data_random(setup[["data"]])) {
    if (!is.null(random_factor_samples)) {
      diagonal <- pairs[["row_1"]] == pairs[["row_2"]]
      lower[, diagonal] <- lower[, diagonal, drop = FALSE] +
        random_factor_samples[["diagonal"]][
          , global_row_1[diagonal], drop = FALSE
        ]
      loading <- random_factor_samples[["loadings"]][[block_index]]
      if (dim(loading)[[3L]] == 0L) {
        return(lower)
      }
      for (column in seq_len(lower_n)) {
        row_1_loading <- matrix(
          loading[, pairs[["row_1"]][[column]], , drop = FALSE],
          nrow = setup[["S"]]
        )
        row_2_loading <- matrix(
          loading[, pairs[["row_2"]][[column]], , drop = FALSE],
          nrow = setup[["S"]]
        )
        lower[, column] <- lower[, column] +
          rowSums(row_1_loading * row_2_loading)
      }
      return(lower)
    }
    if (is.null(random_covariance_samples)) {
      if (is.null(execution_plan[["random_covariance"]])) {
        return(lower)
      }
      stop(
        "Selection random-effect covariance samples are missing.",
        call. = FALSE
      )
    }
    for (column in seq_len(lower_n)) {
      lower[, column] <- lower[, column] + random_covariance_samples[
        , global_row_1[[column]], global_row_2[[column]]
      ]
    }
    return(lower)
  }

  diagonal <- pairs[["row_1"]] == pairs[["row_2"]]
  if (.selection_integrates_estimate(setup[["data"]])) {
    lower[, diagonal] <- lower[, diagonal, drop = FALSE] +
      setup[["tau_within"]][, global_row_1[diagonal], drop = FALSE]^2
  }
  if (isTRUE(setup[["is_multilevel"]]) &&
      !.selection_retains_other_random(setup[["data"]])) {
    cluster      <- setup[["data"]][["outcome"]][["cluster"]]
    same_cluster <- cluster[global_row_1] == cluster[global_row_2]
    lower[, same_cluster] <- lower[, same_cluster, drop = FALSE] +
      setup[["tau_between"]][
        , global_row_1[same_cluster], drop = FALSE
      ] * setup[["tau_between"]][
        , global_row_2[same_cluster], drop = FALSE
      ]
  }

  lower
}


.selection_joint_dense_loglik_block <- function(
    yi, means, covariance_lower, sei, selection_context, execution_plan,
    block_size, return_normalizer = FALSE, normalizer_grid = NULL,
    covariance_grid = NULL) {

  S <- nrow(means)
  if (!return_normalizer && !is.null(covariance_grid)) {
    result <- .selection_covariance_grid_loglik(yi, means, covariance_lower, sei,
      selection_context, execution_plan, block_size, covariance_grid)
    if (!is.null(result)) return(result)
  }
  if (!return_normalizer && !is.null(normalizer_grid)) {
    result <- .selection_normalizer_grid_loglik(yi, means, covariance_lower, sei,
      selection_context, execution_plan, block_size, normalizer_grid)
    if (!is.null(result)) return(result)
  }
  native_static <- BayesTools::selection_native_static_args(selection_context)
  # This local post-fit request never mutates the fitted/JAGS quadrature plan.
  quadrature <- execution_plan[["quadrature"]]
  if (!is.null(quadrature)) attr(quadrature, "bounded_cdf") <- TRUE
  result <- .Call(
    "RoBMA_selnorm_mnorm_step_loglik_batch",
    .native_numeric_vector(yi),
    .native_numeric_matrix(means),
    .native_numeric_matrix(covariance_lower),
    .native_numeric_vector(sei),
    .native_numeric_matrix(selection_context[["omega"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    .native_integer_vector(selection_context[["obs_bin"]]),
    native_static[["sign"]],
    native_static[["telescope_probabilities"]],
    .native_integer_vector(selection_context[["kernel_mode"]]),
    .native_numeric_vector(execution_plan[["designs"]][[
      as.character(block_size)
    ]]),
    .native_integer_vector(execution_plan[["points_per_scramble"]]),
    .native_integer_vector(execution_plan[["scrambles"]]),
    .native_numeric_vector(execution_plan[["relative_tolerance"]]),
    as.logical(return_normalizer),
    .native_integer_vector(selection_context[["vector_rule"]]),
    quadrature,
    PACKAGE = "RoBMA"
  )
  expected_names <- c("log_density", "relative_mcse", "log_normalizer",
                      "integration_diagnostics")
  if (!is.list(result) || !identical(names(result), expected_names) ||
      any(lengths(result[1:3]) != S) ||
      !identical(dim(result[["integration_diagnostics"]]), c(S, 4L))) {
    stop("The selection native kernel returned invalid output.",
         call. = FALSE)
  }
  if (!is.null(quadrature)) {
    cdf_error <- attr(result[["integration_diagnostics"]], "cdf_relative_error", exact = TRUE)
    if (!is.numeric(cdf_error) || length(cdf_error) != S || !is.null(dim(cdf_error)) ||
        any(!is.finite(cdf_error)) || any(cdf_error < 0)) {
      stop("Selection CDF integration diagnostics are unavailable.", call. = FALSE)
    }
  }
  if (return_normalizer) return(result)
  .selection_joint_dense_loglik_check_mcse(result[["relative_mcse"]], execution_plan)

  result[["log_density"]]
}


.selection_joint_singleton_loglik_matrix <- function(
    yi, means, variances, sei, selection_context) {

  sigma <- sqrt(variances)
  .selnorm_kernel_loglik_matrix(
    yi             = yi,
    mu_num         = means,
    sigma_num      = sigma,
    mu_norm        = means,
    sigma_norm     = sigma,
    sei            = sei,
    omega          = selection_context[["omega"]],
    selection_spec = selection_context,
    alpha          = selection_context[["alpha"]],
    phack_kind     = selection_context[["phack_kind"]],
    kernel_mode    = selection_context[["kernel_mode"]]
  )
}


.selection_joint_cluster_loglik_block <- function(
    yi, means, residual_sd, loading, sei, selection_context,
    execution_plan, return_normalizer = FALSE) {

  S <- nrow(means)
  native_static <- BayesTools::selection_native_static_args(selection_context)
  quadrature <- execution_plan[["quadrature"]]
  result <- .Call(
    "RoBMA_selnorm_cluster_step_loglik_batch",
    .native_numeric_vector(yi),
    .native_numeric_matrix(means),
    .native_numeric_matrix(residual_sd),
    .native_numeric_matrix(loading),
    .native_numeric_vector(sei),
    .native_numeric_matrix(selection_context[["omega"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    .native_integer_vector(selection_context[["obs_bin"]]),
    native_static[["sign"]],
    native_static[["telescope_probabilities"]],
    .native_integer_vector(selection_context[["kernel_mode"]]),
    .native_numeric_vector(quadrature[["nodes"]]),
    .native_numeric_vector(quadrature[["log_weights"]]),
    .native_numeric_vector(quadrature[["orders"]]),
    .native_numeric_vector(execution_plan[["designs"]][["factor_1"]]),
    .native_integer_vector(execution_plan[["factor_points_per_proposal"]]),
    .native_integer_vector(execution_plan[["factor_max_points_per_proposal"]]),
    .native_integer_vector(execution_plan[["scrambles"]]),
    .native_numeric_vector(execution_plan[["relative_tolerance"]]),
    as.logical(return_normalizer),
    .native_integer_vector(selection_context[["vector_rule"]]),
    PACKAGE = "RoBMA"
  )
  expected_names <- c("log_density", "relative_mcse", "relative_change",
                      "log_normalizer")
  if (!is.list(result) || !identical(names(result), expected_names) ||
      any(lengths(result) != S)) {
    stop("The cluster selection kernel returned invalid output.",
         call. = FALSE)
  }
  if (return_normalizer) return(result)
  diagnostic <- pmax(result[["relative_mcse"]], result[["relative_change"]])
  failed <- which(!is.finite(diagnostic) |
    diagnostic > execution_plan[["relative_tolerance"]])
  if (length(failed) > 0L) {
    stop(
      "Selection cluster normalizer was rejected by diagnostics: ",
      "relative Monte Carlo standard error was ",
      format(result[["relative_mcse"]][failed[[1L]]], digits = 4),
      " and nested-design relative change was ",
      format(result[["relative_change"]][failed[[1L]]], digits = 4),
      ". Increase 'max_points_per_scramble' or 'scrambles' in ",
      "'set_selection_likelihood_control()'. Pass the control as ",
      "'selection_control' when fitting, 'integration_control' in zplot(), ",
      "or 'density_control$integration_control' for supported post-fit ",
      "densities and hypotheses.",
      call. = FALSE
    )
  }
  result[["log_density"]]
}


.selection_joint_factor_loglik_block <- function(
    yi, means, residual_sd, loading, sei, selection_context,
    execution_plan, block_index, return_normalizer = FALSE) {

  S <- nrow(means)
  factor_rank <- execution_plan[["factor_ranks"]][[block_index]]
  design_key  <- execution_plan[["design_keys"]][[block_index]]
  quadrature <- execution_plan[["factor_quadrature"]][[
    as.character(factor_rank)
  ]]
  if (factor_rank < 2L || factor_rank > 4L ||
      ncol(loading) != length(yi) * factor_rank || is.na(design_key) ||
      length(quadrature[["orders"]]) < 3L) {
    stop("Internal error: selection factor inputs are inconsistent.",
         call. = FALSE)
  }
  native_static <- BayesTools::selection_native_static_args(selection_context)
  result <- .Call(
    "RoBMA_selnorm_factor_step_loglik_batch",
    .native_numeric_vector(yi),
    .native_numeric_matrix(means),
    .native_numeric_matrix(residual_sd),
    .native_numeric_matrix(loading),
    .native_numeric_vector(sei),
    .native_numeric_matrix(selection_context[["omega"]]),
    native_static[["z_lower"]],
    native_static[["z_upper"]],
    .native_integer_vector(selection_context[["obs_bin"]]),
    native_static[["sign"]],
    native_static[["telescope_probabilities"]],
    .native_integer_vector(selection_context[["kernel_mode"]]),
    .native_numeric_vector(quadrature[["nodes"]]),
    .native_numeric_vector(quadrature[["log_weights"]]),
    .native_numeric_vector(quadrature[["orders"]]),
    .native_numeric_vector(quadrature[["rule_counts"]]),
    .native_numeric_vector(execution_plan[["designs"]][[design_key]]),
    .native_integer_vector(
      execution_plan[["factor_points_per_proposal"]]
    ),
    .native_integer_vector(
      execution_plan[["factor_max_points_per_proposal"]]
    ),
    .native_integer_vector(execution_plan[["scrambles"]]),
    .native_numeric_vector(execution_plan[["relative_tolerance"]]),
    as.logical(return_normalizer),
    .native_integer_vector(selection_context[["vector_rule"]]),
    PACKAGE = "RoBMA"
  )
  expected_names <- c("log_density", "relative_mcse", "relative_change",
                      "log_normalizer")
  if (!is.list(result) || !identical(names(result), expected_names) ||
      any(lengths(result) != S)) {
    stop("The selection factor kernel returned invalid output.",
         call. = FALSE)
  }
  if (return_normalizer) return(result)
  diagnostic <- pmax(result[["relative_mcse"]], result[["relative_change"]])
  failed <- which(!is.finite(diagnostic) |
    diagnostic > execution_plan[["relative_tolerance"]])
  if (length(failed) > 0L) {
    stop(
      "Selection factor normalizer was rejected by diagnostics: ",
      "relative Monte Carlo standard error was ",
      format(result[["relative_mcse"]][failed[[1L]]], digits = 4),
      " and nested-design relative change was ",
      format(result[["relative_change"]][failed[[1L]]], digits = 4),
      ". Increase 'max_points_per_scramble' or 'scrambles' in ",
      "'set_selection_likelihood_control()'. Pass the control as ",
      "'selection_control' when fitting, 'integration_control' in zplot(), ",
      "or ",
      "'density_control$integration_control' for supported post-fit ",
      "densities and hypotheses.",
      call. = FALSE
    )
  }
  result[["log_density"]]
}


.selection_joint_block_loglik_from_setup <- function(setup) {

  if (.selection_retains_sampling(setup[["data"]])) {
    return(.selection_conditioned_sampling_state(setup)[["block_log_lik"]])
  }
  if (!.is_data_joint_selection(setup[["data"]])) {
    stop("Selection-likelihood metadata are unavailable.",
         call. = FALSE)
  }
  location <- .estimate_normal_covariance_target_location_from_setup(setup)
  selection_context <- .selection_joint_signed_context(
    setup     = setup,
    signed_yi = location[["y"]]
  )
  random_factor <- .selection_joint_random_factor_samples(setup)
  random_covariance <- if (is.null(random_factor)) {
    .selection_joint_random_covariance_samples(setup)
  } else {
    NULL
  }
  .selection_joint_block_loglik(
    setup             = setup,
    yi                = location[["y"]],
    means             = location[["means"]],
    selection_context = selection_context,
    random_covariance = random_covariance,
    random_factor     = random_factor
  )
}


# The joint bridge target and posterior block scores share the same evaluator.
.selection_joint_block_loglik <- function(
    setup, yi, means, selection_context, random_covariance, random_factor) {

  execution_plan <- .data_selection_execution_plan(setup[["data"]])
  log_lik <- matrix(
    0,
    nrow = setup[["S"]],
    ncol = length(execution_plan[["row_blocks"]])
  )

  singleton_blocks <- execution_plan[["singleton_blocks"]]
  if (length(singleton_blocks) > 0L) {
    rows <- execution_plan[["singleton_rows"]]
    singleton_context <- selection_context
    singleton_context[["obs_bin"]] <- selection_context[["obs_bin"]][rows]
    variances <- .selection_joint_singleton_variances(
      setup                     = setup,
      rows                      = rows,
      block_indices             = singleton_blocks,
      random_covariance_samples = random_covariance,
      random_factor_samples     = random_factor
    )
    log_lik[, singleton_blocks] <- .selection_joint_singleton_loglik_matrix(
      yi                = yi[rows],
      means             = means[, rows, drop = FALSE],
      variances         = variances,
      sei               = setup[["selection_sei"]][rows],
      selection_context = singleton_context
    )
    if (.is_data_weights(setup[["data"]])) {
      log_lik[, singleton_blocks] <- sweep(
        log_lik[, singleton_blocks, drop = FALSE], 2L,
        setup[["data"]][["outcome"]][["weights"]][rows], `*`
      )
    }
  }
  for (block_index in execution_plan[["dependent_blocks"]]) {
    rows <- execution_plan[["row_blocks"]][[block_index]]
    method <- execution_plan[["block_methods"]][[block_index]]
    block_context <- if (is.null(setup[["normalizer_grid"]])) selection_context else
      BayesTools::selection_context_subset_observations(selection_context, rows)
    block_context[["obs_bin"]] <- selection_context[["obs_bin"]][rows]
    if (method %in% c("rank_one", "factor")) {
      components <- .selection_joint_factor_block_samples(
        setup                 = setup,
        block_index           = block_index,
        random_factor_samples = random_factor
      )
    }
    if (method == "rank_one") {
      log_lik[, block_index] <- .selection_joint_cluster_loglik_block(
        yi                 = yi[rows],
        means              = means[, rows, drop = FALSE],
        residual_sd        = components[["residual_sd"]],
        loading            = components[["loading"]],
        sei                = setup[["selection_sei"]][rows],
        selection_context  = block_context,
        execution_plan     = execution_plan
      )
      next
    }
    if (method == "factor") {
      log_lik[, block_index] <- .selection_joint_factor_loglik_block(
        yi                = yi[rows],
        means             = means[, rows, drop = FALSE],
        residual_sd       = components[["residual_sd"]],
        loading           = components[["loading"]],
        sei               = setup[["selection_sei"]][rows],
        selection_context = block_context,
        execution_plan    = execution_plan,
        block_index       = block_index
      )
      next
    }
    log_lik[, block_index] <- .selection_joint_dense_loglik_block(
      yi                = yi[rows],
      means             = means[, rows, drop = FALSE],
      covariance_lower  = .selection_joint_covariance_lower(
        setup                     = setup,
        block_index               = block_index,
        random_covariance_samples = random_covariance,
        random_factor_samples     = random_factor
      ),
      sei               = setup[["selection_sei"]][rows],
      selection_context = block_context,
      execution_plan    = execution_plan,
      block_size        = length(rows),
      normalizer_grid = if (is.null(setup[["normalizer_grid"]])) NULL else list(
        state = setup[["normalizer_grid"]], block = block_index, rows = rows,
        sign = if (identical(setup[["effect_direction"]], "negative")) -1 else 1),
      covariance_grid = if (is.null(setup[["covariance_grid"]])) NULL else list(
        state = setup[["covariance_grid"]], block = block_index, rows = rows,
        sign = if (identical(setup[["effect_direction"]], "negative")) -1 else 1)
    )
  }

  log_lik
}


.selection_joint_loglik_from_setup <- function(setup) {

  rowSums(.selection_joint_block_loglik_from_setup(setup))
}


.selection_joint_singleton_random_variance_expression <- function(
    plan, block_index) {

  random_covariance <- plan[["random_covariance"]]
  if (is.null(random_covariance)) {
    return(NULL)
  }
  if (identical(random_covariance[["representation"]], "dense")) {
    return(paste0(random_covariance[["lower_names"]][[block_index]], "[1]"))
  }
  if (!identical(
      random_covariance[["representation"]],
      "diagonal_factor"
    )) {
    stop("Selection random-covariance metadata are invalid.",
         call. = FALSE)
  }

  terms <- paste0(
    random_covariance[["diagonal_names"]][[block_index]],
    "[1]"
  )
  rank <- random_covariance[["loading_ranks"]][[block_index]]
  if (rank > 0L) {
    loading <- paste0(
      random_covariance[["loading_names"]][[block_index]],
      "[1,1:", rank, "]"
    )
    terms <- c(terms, paste0("inprod(", loading, ",", loading, ")"))
  }
  paste(terms, collapse = " + ")
}


.selection_joint_singleton_model_syntax <- function(
    data, selection_spec, plan) {

  if (length(plan[["singleton_blocks"]]) == 0L) {
    return("")
  }

  random_assignments <- ""
  if (!is.null(plan[["random_covariance"]])) {
    random_assignments <- paste0(vapply(
      seq_along(plan[["singleton_blocks"]]),
      function(index) {
        block_index <- plan[["singleton_blocks"]][[index]]
        paste0(
          "sel_joint_singleton_random_variance[", index, "] = ",
          .selection_joint_singleton_random_variance_expression(
            plan,
            block_index
          ),
          "\n"
        )
      },
      character(1L)
    ), collapse = "")
  }

  row <- "sel_joint_singleton_row[s]"
  variance_terms <- "sel_joint_singleton_sampling_variance[s]"
  if (!is.null(plan[["random_covariance"]])) {
    variance_terms <- c(
      variance_terms,
      "sel_joint_singleton_random_variance[s]"
    )
  } else if (!.is_data_random(data)) {
    tau_within <- if (.is_data_multilevel(data)) "tau_within" else "tau"
    if (.is_data_scale(data)) {
      tau_within <- paste0(tau_within, "[", row, "]")
    }
    if (.selection_integrates_estimate(data)) {
      variance_terms <- c(variance_terms, paste0("pow(", tau_within, ",2)"))
    }
    if (.is_data_multilevel(data) && !.selection_retains_other_random(data)) {
      tau_between <- "tau_between"
      if (.is_data_scale(data)) {
        tau_between <- paste0(tau_between, "[", row, "]")
      }
      variance_terms <- c(
        variance_terms,
        paste0("pow(", tau_between, ",2)")
      )
    }
  }

  paste0(
    random_assignments,
    "for(s in 1:sel_joint_singleton_n){\n",
    "  sel_joint_singleton_variance[s] = ",
    paste(variance_terms, collapse = " + "), "\n",
    "  sel_joint_singleton_y[s] ~ dselnorm_step_switch(",
    "sel_joint_mu[", row, "],sqrt(sel_joint_singleton_variance[s]),",
    "sel_joint_singleton_sei[s],sel_joint_singleton_weights[s],",
    selection_spec[["jags_omega"]], ",",
    "sel_z_lower,sel_z_upper,sel_joint_singleton_obs_bin[s],sel_sign,",
    .selection_joint_kernel_mode_expression(selection_spec), ",",
    "sel_telescope_probabilities)\n",
    "}\n"
  )
}


.selection_joint_model_syntax <- function(data, selection_spec) {

  if (.selection_retains_sampling(data)) {
    return(.selection_conditioned_sampling_model_syntax(data, selection_spec))
  }
  plan <- .data_selection_execution_plan(data)
  random_covariance <- plan[["random_covariance"]]
  random_representation <- if (is.null(random_covariance)) {
    NULL
  } else {
    random_covariance[["representation"]]
  }
  syntax <- ""
  if (!is.null(random_covariance)) {
    syntax <- paste0(syntax, random_covariance[["syntax"]])
  }
  syntax <- paste0(
    syntax,
    .selection_joint_singleton_model_syntax(data, selection_spec, plan)
  )

  for (block_index in plan[["dependent_blocks"]]) {
    rows      <- plan[["row_blocks"]][[block_index]]
    block_n   <- length(rows)
    method    <- plan[["block_methods"]][[block_index]]
    factor_rank <- plan[["factor_ranks"]][[block_index]]
    lower_n <- if (method %in% c("rank_one", "factor")) {
      0L
    } else {
      length(.selection_joint_lower_pairs(plan, rows)[["row_1"]])
    }
    prefix    <- paste0("sel_joint_block_", block_index)
    design_key <- plan[["design_keys"]][[block_index]]
    qmc_name  <- if (!is.na(design_key)) {
      .selection_joint_qmc_name(design_key)
    } else {
      NULL
    }
    covariance_syntax <- ""
    if (method == "dense") {
      random_lower <- if (is.null(random_covariance) ||
          identical(random_representation, "diagonal_factor")) {
        NULL
      } else {
        paste0(
          random_covariance[["lower_names"]][[block_index]],
          "[l]"
        )
      }
      random_factor_covariance <- NULL
      uses_diagonal <- .selection_joint_uses_diagonal_indicator(
        data              = data,
        method            = method,
        random_covariance = random_covariance
      )
      if (!is.null(random_covariance) && uses_diagonal) {
        row_1_expression <- paste0(prefix, "_local_row_1[l]")
        row_2_expression <- paste0(prefix, "_local_row_2[l]")
        random_diagonal <- paste0(
          prefix, "_diagonal[l] * ",
          random_covariance[["diagonal_names"]][[block_index]],
          "[", row_1_expression, "]"
        )
        random_rank <- random_covariance[["loading_ranks"]][[block_index]]
        random_product <- if (random_rank == 0L) {
          NULL
        } else {
          paste0(
            "inprod(", random_covariance[["loading_names"]][[block_index]],
            "[", row_1_expression, ",1:", random_rank, "],",
            random_covariance[["loading_names"]][[block_index]],
            "[", row_2_expression, ",1:", random_rank, "])"
          )
        }
        random_factor_covariance <- paste(
          c(random_diagonal, random_product),
          collapse = " + "
        )
      }

      extra_covariance <- NULL
      if (!.is_data_random(data)) {
        within <- if (.is_data_multilevel(data)) "tau_within" else "tau"
        if (.is_data_scale(data)) {
          within <- paste0(within, "[", prefix, "_row_1[l]]")
        }
        if (.selection_integrates_estimate(data)) {
          extra_covariance <- paste0(prefix, "_diagonal[l] * pow(", within, ",2)")
        }
        if (.is_data_multilevel(data) && !.selection_retains_other_random(data)) {
          between <- if (.is_data_scale(data)) {
            paste0("tau_between[", prefix, "_row_1[l]] * tau_between[",
                   prefix, "_row_2[l]]")
          } else "pow(tau_between,2)"
          extra_covariance <- c(
            extra_covariance, paste0(prefix, "_same_cluster[l] * ", between)
          )
        }
      }
      covariance_terms <- c(
        paste0(prefix, "_sampling_lower[l]"),
        random_lower,
        random_factor_covariance,
        extra_covariance
      )
      covariance_syntax <- paste0(
        "for(l in 1:", lower_n, "){\n",
        "  ", prefix, "_covariance[l] = ",
        paste(covariance_terms, collapse = " + "), "\n",
        "}\n"
      )
    }

    factor_setup_syntax <- ""
    if (method %in% c("rank_one", "factor")) {
      residual_terms <- paste0(prefix, "_sampling_variance[j]")
      loading_expressions <- character()
      sampling_rank <- plan[["sampling_factor_blocks"]][[
        block_index
      ]][["rank"]]
      if (sampling_rank > 0L) {
        loading_expressions <- c(
          loading_expressions,
          paste0(prefix, "_sampling_loading[j,", seq_len(sampling_rank), "]")
        )
      }
      if (.is_data_random(data) && !is.null(random_covariance)) {
        residual_terms <- c(
          residual_terms,
          paste0(
            random_covariance[["diagonal_names"]][[block_index]],
            "[j]"
          )
        )
        random_rank <- random_covariance[["loading_ranks"]][[block_index]]
        if (random_rank > 0L) {
          loading_expressions <- c(
            loading_expressions,
            paste0(
              random_covariance[["loading_names"]][[block_index]],
              "[j,", seq_len(random_rank), "]"
            )
          )
        }
      } else if (!.is_data_random(data)) {
        tau_within <- if (.is_data_multilevel(data)) "tau_within" else "tau"
        if (.is_data_scale(data)) {
          tau_within <- paste0(tau_within, "[", prefix, "_row[j]]")
        }
        if (.selection_integrates_estimate(data)) {
          residual_terms <- c(
            residual_terms,
            paste0("pow(", tau_within, ",2)")
          )
        }
        if (.is_data_multilevel(data) && !.selection_retains_other_random(data)) {
          tau_between <- "tau_between"
          if (.is_data_scale(data)) {
            tau_between <- paste0(
              tau_between, "[", prefix, "_row[j]]"
            )
          }
          cluster_rank <- length(unique(data[["outcome"]][["cluster"]][rows]))
          loading_expressions <- c(
            loading_expressions,
            paste0(tau_between, " * ", prefix, "_cluster_loading[j,",
                   seq_len(cluster_rank), "]")
          )
        }
      }
      if (length(loading_expressions) != factor_rank) {
        stop(
          "Internal error: selection factor syntax rank is inconsistent.",
          call. = FALSE
        )
      }
      loading_assignments <- if (method == "rank_one") {
        paste0("  ", prefix, "_loading[j] = ", loading_expressions, "\n")
      } else {
        paste0(
          "  ", prefix, "_loading[j,", seq_len(factor_rank), "] = ",
          loading_expressions, "\n",
          collapse = ""
        )
      }
      factor_setup_syntax <- paste0(
        "for(j in 1:", block_n, "){\n",
        "  ", prefix, "_residual_sd[j] = sqrt(",
        paste(residual_terms, collapse = " + "), ")\n",
        loading_assignments,
        "}\n"
      )
    }

    density_syntax <- if (method == "rank_one") {
      paste0(
        prefix, "_y[1:", block_n, "] ~ dselnorm_cluster_step(",
        prefix, "_mu[1:", block_n, "],",
        prefix, "_residual_sd[1:", block_n, "],",
        prefix, "_loading[1:", block_n, "],",
        prefix, "_sei[1:", block_n, "],",
        selection_spec[["jags_omega"]], ",",
        "sel_z_lower,sel_z_upper,",
        prefix, "_obs_bin[1:", block_n, "],",
        "sel_sign,sel_telescope_probabilities,",
        .selection_joint_kernel_mode_expression(selection_spec), ",",
        "sel_joint_cluster_nodes,sel_joint_cluster_log_weights,",
        "sel_joint_cluster_orders,",
        qmc_name, "[1:", plan[["scrambles"]], ",1:",
        plan[["factor_max_points_per_proposal"]], ",1:2],",
        plan[["factor_points_per_proposal"]], ",",
        plan[["factor_max_points_per_proposal"]], ",",
        plan[["scrambles"]], ",",
        format(plan[["relative_tolerance"]], scientific = FALSE), ",",
        selection_spec[["jags_vector_rule"]], ")\n"
      )
    } else if (method == "factor") {
      paste0(
        prefix, "_y[1:", block_n, "] ~ dselnorm_factor_step(",
        prefix, "_mu[1:", block_n, "],",
        prefix, "_residual_sd[1:", block_n, "],",
        prefix, "_loading[1:", block_n, ",1:", factor_rank, "],",
        prefix, "_sei[1:", block_n, "],",
        selection_spec[["jags_omega"]], ",",
        "sel_z_lower,sel_z_upper,",
        prefix, "_obs_bin[1:", block_n, "],",
        "sel_sign,sel_telescope_probabilities,",
        .selection_joint_kernel_mode_expression(selection_spec), ",",
        .selection_joint_factor_quadrature_name("nodes", factor_rank), ",",
        .selection_joint_factor_quadrature_name(
          "log_weights", factor_rank
        ), ",",
        .selection_joint_factor_quadrature_name("orders", factor_rank), ",",
        .selection_joint_factor_quadrature_name("rule_counts", factor_rank), ",",
        qmc_name, "[1:", plan[["scrambles"]], ",1:",
        plan[["factor_max_points_per_proposal"]], ",1:",
        2L * factor_rank, "],",
        plan[["factor_points_per_proposal"]], ",",
        plan[["factor_max_points_per_proposal"]], ",",
        plan[["scrambles"]], ",",
        format(plan[["relative_tolerance"]], scientific = FALSE), ",",
        selection_spec[["jags_vector_rule"]], ")\n"
      )
    } else {
      paste0(
        prefix, "_y[1:", block_n, "] ~ dselnorm_mnorm_step(",
        prefix, "_mu[1:", block_n, "],",
        prefix, "_covariance[1:", lower_n, "],",
        prefix, "_sei[1:", block_n, "],",
        selection_spec[["jags_omega"]], ",",
        "sel_z_lower,sel_z_upper,",
        prefix, "_obs_bin[1:", block_n, "],",
        "sel_sign,sel_telescope_probabilities,",
        .selection_joint_kernel_mode_expression(selection_spec), ",",
        qmc_name, "[1:", plan[["scrambles"]], ",1:",
        plan[["points_per_scramble"]], ",1:", 2L * block_n, "],",
        plan[["points_per_scramble"]], ",",
        plan[["scrambles"]], ",",
        format(plan[["relative_tolerance"]], scientific = FALSE), ",",
        selection_spec[["jags_vector_rule"]], ",",
        "sel_joint_cluster_nodes,sel_joint_cluster_log_weights,",
        "sel_joint_cluster_orders)\n"
      )
    }

    syntax <- paste0(
      syntax,
      "for(j in 1:", block_n, "){\n",
      "  ", prefix, "_mu[j] = sel_joint_mu[", prefix, "_row[j]]\n",
      "}\n",
      factor_setup_syntax,
      covariance_syntax,
      density_syntax
    )
  }

  syntax
}


.selection_joint_dense_loglik_check_mcse <- function(relative_mcse, execution_plan) {

  failed <- which(
    !is.finite(relative_mcse) |
      relative_mcse > execution_plan[["relative_tolerance"]]
  )
  if (length(failed) > 0L) {
    observed <- relative_mcse[failed[[1L]]]
    stop(
      "Selection normalizer was rejected by diagnostics: relative ",
      "Monte Carlo standard error was ", format(observed, digits = 4),
      ". Increase 'points_per_scramble' in ",
      "'set_selection_likelihood_control()'. Pass the control as ",
      "'selection_control' when fitting, 'integration_control' in zplot(), ",
      "or ",
      "'density_control$integration_control' for supported post-fit ",
      "densities and hypotheses.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
