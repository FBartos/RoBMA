# ============================================================================ #
# selection-likelihood.R
# ============================================================================ #
#
# Exact finite-vector selection likelihood preparation. Generic numerical and
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


#' Control numerical integration of exact selection likelihoods
#'
#' @description
#' Creates numerical integration settings for the finite-vector product
#' selection likelihood used by [bselmodel()], [bselmodel.mv()], [RoBMA()], and
#' [RoBMA.mv()]. Covariance factors use deterministic sequences of
#' Gauss-Hermite rules. Certified factors of structural rank two through four
#' fall back to randomized quasi-Monte Carlo in the factor dimension when the
#' quadrature sequence does not converge. General covariance blocks use
#' randomized quasi-Monte Carlo in the dense conditional-normal representation.
#' Every integration design is fixed before fitting, so repeated likelihood
#' evaluations are deterministic. One-observation blocks are evaluated
#' analytically.
#'
#' The numerical defaults are package operating settings, not statistical
#' estimands or literature-mandated constants. `points_per_scramble`,
#' `max_points_per_scramble`, and `scrambles` set the fixed computational
#' design, while
#' `relative_tolerance` sets its diagnostic acceptance rule. Increasing a
#' budget changes only the numerical evaluation of the same likelihood target.
#'
#' @param points_per_scramble nominal number of integrand evaluations per
#'   randomized quasi-Monte Carlo scramble for certified higher-rank factor
#'   fallbacks, and the number of shifted Halton points for general covariance
#'   blocks.
#'   Factor evaluations are allocated equally to the prior- and mode-centered
#'   proposal components; odd values therefore use the next even total. This
#'   count is the initial nested-design budget, not an accuracy guarantee. The
#'   integration diagnostics determine whether it is adequate for each
#'   evaluated likelihood state.
#' @param max_points_per_scramble largest nested-design budget used for a
#'   certified higher-rank factor fallback when the initial
#'   `points_per_scramble` budget fails its integration diagnostic. Refinement
#'   doubles the point count up to this value while reusing the same fixed
#'   design.
#' @param scrambles number of independent randomized shifts used to estimate
#'   integration error for certified higher-rank factors and general covariance
#'   blocks.
#' @param relative_tolerance largest accepted successive relative quadrature
#'   change for factor integration, or the largest nested-design relative
#'   change and relative Monte Carlo standard error for a higher-rank factor
#'   fallback. General covariance blocks use the relative Monte Carlo standard
#'   error.
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


.check_selection_likelihood_control <- function(control) {

  if (!inherits(control, "RoBMA_selection_likelihood_control")) {
    stop(
      "'selection_control' must be created by ",
      "set_selection_likelihood_control().",
      call. = FALSE
    )
  }
  do.call(set_selection_likelihood_control, unclass(control))
}


.set_data_selection_likelihood <- function(data, likelihood, setup = NULL) {

  attr(data, "selection_likelihood")       <- likelihood
  attr(data, "exact_selection_likelihood") <- setup
  data
}


.data_selection_likelihood <- function(data) {

  likelihood <- attr(data, "selection_likelihood", exact = TRUE)
  if (is.null(likelihood)) {
    return("approximate")
  }
  likelihood
}


.is_data_exact_selection <- function(data) {

  identical(.data_selection_likelihood(data), "exact")
}


.uses_exact_selection_likelihood <- function(data, priors) {

  .is_data_exact_selection(data) && .is_priors_weightfunction(priors)
}


.setup_uses_exact_selection_likelihood <- function(setup) {

  .uses_exact_selection_likelihood(setup[["data"]], setup[["priors"]])
}


.data_exact_selection_setup <- function(data) {

  plan <- attr(data, "exact_selection_likelihood", exact = TRUE)
  if (.is_data_exact_selection(data) && is.null(plan)) {
    stop(
      "Internal error: exact selection-likelihood metadata are missing.",
      call. = FALSE
    )
  }
  if (!is.null(plan) &&
      (!inherits(plan, "RoBMA_selection_execution_plan") ||
       !identical(plan[["schema_version"]], 1L))) {
    stop(
      "Internal error: exact selection execution-plan metadata are invalid.",
      call. = FALSE
    )
  }
  plan
}


.prepare_selection_likelihood_object <- function(
    object, selection_likelihood, selection_control) {

  selection_likelihood <- match.arg(
    selection_likelihood,
    c("exact", "approximate")
  )
  selection_control <- .check_selection_likelihood_control(selection_control)

  if (!.is_priors_weightfunction(object[["priors"]])) {
    stop(
      "Selection-likelihood preparation requires a selection-kernel prior.",
      call. = FALSE
    )
  }
  if (identical(selection_likelihood, "approximate")) {
    object[["data"]] <- .set_data_selection_likelihood(
      data       = object[["data"]],
      likelihood = "approximate"
    )
    object[["selection_likelihood"]] <- list(
      type    = "approximate",
      target  = "row_selected_normal_conditional_on_random_effects",
      control = selection_control
    )
    return(object)
  }

  if (.is_data_weights(object[["data"]])) {
    stop(
      "'weights' are unavailable with 'selection_likelihood = \"exact\"' ",
      "because a powered joint selection density is not a generative ",
      "finite-vector selection model. Use 'selection_likelihood = ",
      "\"approximate\"' or omit 'weights'.",
      call. = FALSE
    )
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
      "'selection_likelihood = \"exact\"' currently requires a step ",
      "selection kernel. P-hacking and other non-step kernels are unavailable.",
      call. = FALSE
    )
  }

  formula_design <- NULL
  if (.is_data_random(object[["data"]])) {
    formula_design <- object[["formula_design"]][["mu"]]
    if (!inherits(formula_design, "BayesTools_formula_design") ||
        any(vapply(
          formula_design[["random_effects"]],
          function(term) !identical(term[["compile_mode"]], "marginalized"),
          logical(1)
        ))) {
      stop(
        "Internal error: exact selection random effects were not compiled ",
        "as marginalized.",
        call. = FALSE
      )
    }
  }

  sampling <- .selection_exact_sampling_plan(object[["data"]])
  row_blocks <- .selection_exact_dependency_blocks(
    data         = object[["data"]],
    random_terms = if (is.null(formula_design)) {
      list()
    } else {
      formula_design[["random_effects"]]
    },
    sampling     = sampling
  )
  sampling_factor_blocks <- if (!identical(
      sampling[["representation"]],
      "diagonal_factor"
    )) {
    NULL
  } else {
    .selection_exact_sampling_factor_blocks(sampling, row_blocks)
  }
  random_covariance <- NULL
  if (!is.null(formula_design)) {
    random_covariance <- BayesTools::JAGS_formula_random_marginal_covariance(
      formula_design = formula_design,
      row_blocks     = row_blocks,
      prefix         = "sel_exact_random",
      representation = if (!is.null(sampling_factor_blocks)) {
        "auto"
      } else "dense"
    )
  }
  block_routing <- .selection_exact_block_routing(
    data                   = object[["data"]],
    row_blocks             = row_blocks,
    sampling_factor_blocks = sampling_factor_blocks,
    random_covariance      = random_covariance
  )
  plan <- .selection_exact_execution_plan(
    row_blocks        = row_blocks,
    block_methods     = block_routing[["methods"]],
    factor_ranks      = block_routing[["ranks"]],
    selection_control = selection_control,
    sampling          = sampling,
    sampling_factor_blocks = sampling_factor_blocks,
    random_covariance = random_covariance
  )
  object[["data"]] <- .set_data_selection_likelihood(
    data       = object[["data"]],
    likelihood = "exact",
    setup      = plan
  )
  object[["selection_likelihood"]] <- list(
    type      = "exact",
    target    = plan[["statistical_target"]],
    exactness = plan[["exactness"]],
    control   = selection_control
  )

  object
}


.prepare_approximate_selection_known_v <- function(
    object, selection_likelihood) {

  if (!identical(selection_likelihood, "approximate") ||
      !.is_priors_weightfunction(object[["priors"]]) ||
      !.is_data_known_v(object[["data"]])) {
    return(object)
  }

  known_V   <- .data_known_v_data(object[["data"]])
  requested <- .known_v_requested_parameterization(known_V)
  if (!requested %in% c("auto", "latent")) {
    stop(
      "'selection_likelihood = \"approximate\"' requires ",
      "'known_v_parameterization = \"latent\"'.",
      call. = FALSE
    )
  }
  if (.known_v_effective_backend(known_V) %in% c("latent", "diagonal")) {
    return(object)
  }

  residual_fraction <- .known_v_requested_residual_fraction(known_V)
  new_known_V <- .known_v_prepare(
    V                                   = .known_v_as_input(known_V),
    keep_rows                           = rep(TRUE, .known_v_nrow(known_V)),
    known_v_parameterization            = "latent",
    known_v_residual_fraction           = residual_fraction,
    known_v_residual_fraction_specified = !is.null(residual_fraction),
    known_v_is_scale                    = .is_data_scale(object[["data"]]),
    warn_singular                       = FALSE
  )
  if (identical(requested, "auto")) {
    new_known_V <- .known_v_update(
      new_known_V,
      list(parameterization_requested = "auto")
    )
  }
  attr(object[["data"]], "known_V_data") <- new_known_V
  object[["data"]][["outcome"]][["sei"]] <- sqrt(
    .known_v_diagonal(new_known_V)
  )

  object
}


.selection_exact_sampling_plan <- function(data) {

  factor <- .selection_exact_sampling_factor(data)
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


.selection_exact_sampling_diagonal <- function(sampling) {

  if (identical(sampling[["representation"]], "dense")) {
    return(diag(sampling[["covariance"]]))
  }
  if (!identical(sampling[["representation"]], "diagonal_factor")) {
    stop("Exact selection sampling-covariance metadata are invalid.",
         call. = FALSE)
  }

  sampling[["diagonal"]] + rowSums(sampling[["loading"]]^2)
}


.selection_exact_sampling_block <- function(sampling, rows) {

  if (identical(sampling[["representation"]], "dense")) {
    return(sampling[["covariance"]][rows, rows, drop = FALSE])
  }
  if (!identical(sampling[["representation"]], "diagonal_factor")) {
    stop("Exact selection sampling-covariance metadata are invalid.",
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


.selection_exact_sampling_factor <- function(data) {

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
    return(list(
      diagonal = as.numeric(known_V[["factor_diagonal"]]),
      loading  = unname(known_V[["factor_loading"]])
    ))
  }
  NULL
}


.selection_exact_sampling_factor_blocks <- function(factor, row_blocks) {

  diagonal <- factor[["diagonal"]]
  loading  <- factor[["loading"]]
  K        <- length(diagonal)
  if (!is.numeric(diagonal) || anyNA(diagonal) || any(!is.finite(diagonal)) ||
      any(diagonal < 0) || !is.numeric(loading) || !is.matrix(loading) ||
      nrow(loading) != K || anyNA(loading) || any(!is.finite(loading))) {
    stop("Exact selection sampling-factor metadata are invalid.",
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
        "Internal error: a certified sampling factor crosses exact-selection blocks.",
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


.selection_exact_dependency_blocks <- function(
    data, random_terms = list(), sampling = NULL) {

  if (!is.list(random_terms)) {
    stop("Random-effect dependency metadata must be a list.",
         call. = FALSE)
  }
  if (is.null(sampling)) {
    sampling <- .selection_exact_sampling_plan(data)
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

  if (.is_data_multilevel(data)) {
    cluster   <- data[["outcome"]][["cluster"]]
    adjacency <- adjacency | outer(cluster, cluster, "==")
  }
  diag(adjacency) <- TRUE

  .known_v_block_indices(adjacency * 1)
}


.selection_exact_lower_pairs <- function(plan, rows) {

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
      "Internal error: exact selection covariance ordering is invalid.",
      call. = FALSE
    )
  }
  list(
    row_1 = rows[pairs[["row_1"]]],
    row_2 = rows[pairs[["row_2"]]]
  )
}


.selection_exact_uses_diagonal_indicator <- function(
    data, method, random_covariance) {

  method %in% c("singleton", "dense") &&
    (!.is_data_random(data) ||
      identical(random_covariance[["representation"]], "diagonal_factor"))
}


.selection_exact_fit_data <- function(data, priors) {

  plan             <- .data_exact_selection_setup(data)
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
  if (.is_priors_PET(priors) || .is_priors_PEESE(priors)) {
    fit_data[["sei"]] <- sei
  }

  block_methods <- plan[["block_methods"]]
  if (any(block_methods == "rank_one")) {
    quadrature <- .selection_exact_cluster_quadrature_rules(
      plan[["quadrature_orders"]]
    )
    fit_data[["sel_exact_cluster_nodes"]] <- quadrature[["nodes"]]
    fit_data[["sel_exact_cluster_log_weights"]] <-
      quadrature[["log_weights"]]
    fit_data[["sel_exact_cluster_orders"]] <- quadrature[["orders"]]
  }
  for (rank_name in names(plan[["factor_quadrature_orders"]])) {
    quadrature <- .selection_exact_cluster_quadrature_rules(
      plan[["factor_quadrature_orders"]][[rank_name]]
    )
    fit_data[[.selection_exact_factor_quadrature_name(
      "nodes", rank_name
    )]] <- quadrature[["nodes"]]
    fit_data[[.selection_exact_factor_quadrature_name(
      "log_weights", rank_name
    )]] <- quadrature[["log_weights"]]
    fit_data[[.selection_exact_factor_quadrature_name(
      "orders", rank_name
    )]] <- quadrature[["orders"]]
  }
  for (design_key in names(plan[["designs"]])) {
    fit_data[[.selection_exact_qmc_name(design_key)]] <-
      plan[["designs"]][[design_key]]
  }

  if (length(plan[["singleton_blocks"]]) > 0L) {
    rows <- plan[["singleton_rows"]]
    fit_data[["sel_exact_singleton_n"]] <- length(rows)
    fit_data[["sel_exact_singleton_y"]] <- yi[rows]
    fit_data[["sel_exact_singleton_sei"]] <- sei[rows]
    fit_data[["sel_exact_singleton_obs_bin"]] <-
      selection_spec[["obs_bin"]][rows]
    fit_data[["sel_exact_singleton_row"]] <- rows
    fit_data[["sel_exact_singleton_sampling_variance"]] <-
      .selection_exact_sampling_diagonal(plan[["sampling"]])[rows]
  }

  for (block_index in plan[["dependent_blocks"]]) {
    rows   <- plan[["row_blocks"]][[block_index]]
    prefix <- paste0("sel_exact_block_", block_index)
    method <- block_methods[[block_index]]
    fit_data[[paste0(prefix, "_y")]]       <- yi[rows]
    fit_data[[paste0(prefix, "_sei")]]     <- sei[rows]
    fit_data[[paste0(prefix, "_obs_bin")]] <- selection_spec[["obs_bin"]][rows]
    fit_data[[paste0(prefix, "_row")]]     <- rows
    uses_diagonal <- .selection_exact_uses_diagonal_indicator(
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
    }
    if (method == "dense") {
      pairs <- .selection_exact_lower_pairs(plan, rows)
      if (.is_data_random(data) && method == "dense" && uses_diagonal) {
        fit_data[[paste0(prefix, "_local_row_1")]] <-
          match(pairs[["row_1"]], rows)
        fit_data[[paste0(prefix, "_local_row_2")]] <-
          match(pairs[["row_2"]], rows)
      }
    }
    if (uses_diagonal) {
      fit_data[[paste0(prefix, "_diagonal")]] <- as.integer(
        pairs[["row_1"]] == pairs[["row_2"]]
      )
      if (!.is_data_random(data) && .is_data_scale(data)) {
        fit_data[[paste0(prefix, "_row_1")]] <- pairs[["row_1"]]
        if (.is_data_multilevel(data)) {
          fit_data[[paste0(prefix, "_row_2")]] <- pairs[["row_2"]]
        }
      }
    }
    if (method == "dense") {
      sampling_covariance <- .selection_exact_sampling_block(
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


.selection_exact_block_routing <- function(
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
    specialized_rank <- as.integer(
      .is_data_multilevel(data) && !.is_data_random(data)
    )
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


.selection_exact_execution_plan <- function(
    row_blocks, block_methods, factor_ranks, selection_control, sampling,
    sampling_factor_blocks, random_covariance) {

  block_sizes <- lengths(row_blocks)
  if (length(block_sizes) != length(block_methods) ||
      length(block_sizes) != length(factor_ranks)) {
    stop("Internal error: exact-selection integration routing is invalid.",
         call. = FALSE)
  }
  design_keys <- rep(NA_character_, length(block_sizes))
  dense <- which(block_methods == "dense" & block_sizes > 1L)
  factor <- which(block_methods == "factor")
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
    dimensions <- if (block_methods[[block_index]] == "factor") {
      2L * factor_ranks[[block_index]]
    } else {
      2L * block_sizes[[block_index]]
    }
    points <- if (block_methods[[block_index]] == "factor") {
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
    schema_version      = 1L,
    statistical_target  = "finite_vector_product_selection",
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
    quadrature_orders    = if (any(block_methods == "rank_one")) {
      SELNORM_CLUSTER_QUADRATURE_ORDERS
    } else integer(),
    factor_quadrature_orders = SELNORM_FACTOR_QUADRATURE_ORDERS[
      intersect(names(SELNORM_FACTOR_QUADRATURE_ORDERS),
                as.character(unique(factor_ranks[factor])))
    ],
    lower_pairs         = lower_pairs,
    exactness           = exactness
  ), class = c("RoBMA_selection_execution_plan", "list"))
}


.selection_exact_cluster_quadrature_rules <- function(orders) {

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


.selection_exact_qmc_name <- function(block_n) {

  paste0("sel_exact_qmc_", block_n)
}


.selection_exact_factor_quadrature_name <- function(quantity, rank) {

  paste0("sel_exact_factor_", quantity, "_", rank)
}


.selection_exact_kernel_mode_expression <- function(selection_spec) {

  if (isTRUE(selection_spec[["jags_use_step_switch"]])) {
    return(selection_spec[["jags_kernel_mode"]])
  }
  as.character(SELKERNEL_STEP)
}


.selection_exact_signed_context <- function(setup, signed_yi) {

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


.selection_exact_random_covariance_samples <- function(setup) {

  if (!.is_data_random(setup[["data"]])) {
    return(NULL)
  }
  object <- list(
    fit    = setup[["fit"]],
    data   = setup[["data"]],
    priors = setup[["priors"]]
  )
  random_vcov <- .brma_mv_random_effects_marginal_vcov(
    object            = object,
    posterior_samples = setup[["posterior_samples"]]
  )
  samples <- random_vcov[["samples"]]
  if (!identical(dim(samples), c(setup[["S"]], setup[["K"]], setup[["K"]])) ||
      any(!is.finite(samples))) {
    stop(
      "Exact selection random-effect covariance samples are invalid.",
      call. = FALSE
    )
  }

  samples
}


.selection_exact_random_factor_samples <- function(setup) {

  if (!.is_data_random(setup[["data"]])) {
    return(NULL)
  }
  exact_setup <- .data_exact_selection_setup(setup[["data"]])
  representation <- exact_setup[["random_covariance"]][["representation"]]
  if (!identical(representation, "diagonal_factor")) {
    return(NULL)
  }
  object <- list(
    fit    = setup[["fit"]],
    data   = setup[["data"]],
    priors = setup[["priors"]]
  )
  factors <- .brma_mv_random_effects_marginal_factor_states(
    object            = object,
    posterior_samples = setup[["posterior_samples"]],
    blocks            = NULL,
    row_blocks        = exact_setup[["row_blocks"]]
  )
  result <- BayesTools::random_effects_marginal_diagonal_factor(factors)
  expected <- c(setup[["S"]], setup[["K"]])
  if (!identical(dim(result[["diagonal"]]), expected) ||
      any(!is.finite(result[["diagonal"]])) ||
      any(result[["diagonal"]] < 0) ||
      length(result[["loadings"]]) != length(exact_setup[["row_blocks"]]) ||
      !identical(
        result[["ranks"]],
        as.integer(exact_setup[["random_covariance"]][["loading_ranks"]])
      )) {
    stop("Exact selection random-effect factors are invalid.",
         call. = FALSE)
  }
  for (block_index in seq_along(result[["loadings"]])) {
    expected_dim <- c(
      setup[["S"]],
      length(exact_setup[["row_blocks"]][[block_index]]),
      result[["ranks"]][[block_index]]
    )
    if (!identical(dim(result[["loadings"]][[block_index]]), expected_dim) ||
        any(!is.finite(result[["loadings"]][[block_index]]))) {
      stop("Exact selection random-effect factor loadings are invalid.",
           call. = FALSE)
    }
  }
  result
}


.selection_exact_singleton_variances <- function(
    setup, rows, block_indices, random_covariance_samples = NULL,
    random_factor_samples = NULL) {

  exact_setup <- .data_exact_selection_setup(setup[["data"]])
  S           <- setup[["S"]]
  variances   <- matrix(
    .selection_exact_sampling_diagonal(exact_setup[["sampling"]])[rows],
    nrow = S,
    ncol = length(rows),
    byrow = TRUE
  )

  if (!.is_data_random(setup[["data"]])) {
    variances <- variances + setup[["tau_within"]][, rows, drop = FALSE]^2
    if (isTRUE(setup[["is_multilevel"]])) {
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
    stop(
      "Exact selection random-effect covariance samples are missing.",
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


.selection_exact_factor_block_samples <- function(
    setup, block_index, random_factor_samples = NULL) {

  exact_setup <- .data_exact_selection_setup(setup[["data"]])
  rows        <- exact_setup[["row_blocks"]][[block_index]]
  sampling    <- exact_setup[["sampling_factor_blocks"]][[block_index]]
  S           <- setup[["S"]]
  block_n     <- length(rows)
  residual_variance <- matrix(
    sampling[["diagonal"]],
    nrow = S,
    ncol = block_n,
    byrow = TRUE
  )
  loading_parts <- list()
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
  }

  if (.is_data_random(setup[["data"]])) {
    if (is.null(random_factor_samples)) {
      stop("Exact selection random-effect factors are missing.",
           call. = FALSE)
    }
    residual_variance <- residual_variance +
      random_factor_samples[["diagonal"]][, rows, drop = FALSE]
    if (random_factor_samples[["ranks"]][[block_index]] > 0L) {
      loading_parts[[length(loading_parts) + 1L]] <-
        random_factor_samples[["loadings"]][[block_index]]
    }
  } else {
    residual_variance <- residual_variance +
      setup[["tau_within"]][, rows, drop = FALSE]^2
    if (isTRUE(setup[["is_multilevel"]])) {
      part <- array(
        setup[["tau_between"]][, rows, drop = FALSE],
        dim = c(S, block_n, 1L)
      )
      loading_parts[[length(loading_parts) + 1L]] <- part
    }
  }
  if (any(!is.finite(residual_variance)) || any(residual_variance <= 0)) {
    stop(
      "Exact selection factor residual variances must be finite and positive.",
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
    rank        = if (block_n == 0L) 0L else ncol(loading) %/% block_n
  )
}


.selection_exact_covariance_lower <- function(
    setup, rows, random_covariance_samples = NULL,
    random_factor_samples = NULL, block_index = NULL) {

  exact_setup <- .data_exact_selection_setup(setup[["data"]])
  sampling    <- .selection_exact_sampling_block(
    exact_setup[["sampling"]],
    rows
  )
  block_n     <- length(rows)
  pairs       <- .selection_exact_lower_pairs(
    exact_setup,
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
      if (is.null(block_index)) {
        block_index <- which(vapply(
          exact_setup[["row_blocks"]],
          identical,
          logical(1L),
          rows
        ))
      }
      if (length(block_index) != 1L) {
        stop("Internal error: exact-selection factor block is ambiguous.",
             call. = FALSE)
      }
      diagonal <- pairs[["row_1"]] == pairs[["row_2"]]
      lower[, diagonal] <- lower[, diagonal, drop = FALSE] +
        random_factor_samples[["diagonal"]][
          , global_row_1[diagonal], drop = FALSE
        ]
      loading <- random_factor_samples[["loadings"]][[block_index]]
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
      stop(
        "Exact selection random-effect covariance samples are missing.",
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
  lower[, diagonal] <- lower[, diagonal, drop = FALSE] +
    setup[["tau_within"]][, global_row_1[diagonal], drop = FALSE]^2
  if (isTRUE(setup[["is_multilevel"]])) {
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


.selection_exact_joint_loglik_block <- function(
    yi, means, covariance_lower, sei, selection_context, execution_plan,
    block_size) {

  S <- nrow(means)
  native_args <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_context,
    S              = S,
    kernel_mode    = selection_context[["kernel_mode"]]
  )
  native_static <- native_args[["static"]]
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
    .native_integer_vector(native_args[["kernel_mode"]]),
    .native_numeric_vector(execution_plan[["designs"]][[
      as.character(block_size)
    ]]),
    .native_integer_vector(execution_plan[["points_per_scramble"]]),
    .native_integer_vector(execution_plan[["scrambles"]]),
    .native_numeric_vector(execution_plan[["relative_tolerance"]]),
    PACKAGE = "RoBMA"
  )
  if (!is.list(result) ||
      !identical(names(result), c("log_density", "relative_mcse")) ||
      length(result[["log_density"]]) != S ||
      length(result[["relative_mcse"]]) != S) {
    stop("The exact selection native kernel returned invalid output.",
         call. = FALSE)
  }
  failed <- which(
    !is.finite(result[["relative_mcse"]]) |
      result[["relative_mcse"]] > execution_plan[["relative_tolerance"]]
  )
  if (length(failed) > 0L) {
    observed <- result[["relative_mcse"]][failed[[1L]]]
    stop(
      "Exact selection normalizer was rejected by diagnostics: relative ",
      "Monte Carlo standard error was ", format(observed, digits = 4),
      ". Refit with more integration points in 'selection_control', for ",
      "example by increasing 'points_per_scramble'.",
      call. = FALSE
    )
  }

  result[["log_density"]]
}


.selection_exact_singleton_loglik_matrix <- function(
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


.selection_exact_cluster_loglik_block <- function(
    yi, means, residual_sd, loading, sei, selection_context,
    execution_plan) {

  S <- nrow(means)
  native_args <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_context,
    S              = S,
    kernel_mode    = selection_context[["kernel_mode"]]
  )
  native_static <- native_args[["static"]]
  quadrature <- .selection_exact_cluster_quadrature_rules(
    execution_plan[["quadrature_orders"]]
  )
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
    .native_integer_vector(native_args[["kernel_mode"]]),
    .native_numeric_vector(quadrature[["nodes"]]),
    .native_numeric_vector(quadrature[["log_weights"]]),
    .native_numeric_vector(quadrature[["orders"]]),
    .native_numeric_vector(execution_plan[["relative_tolerance"]]),
    PACKAGE = "RoBMA"
  )
  if (!is.list(result) ||
      !identical(names(result), c("log_density", "relative_change")) ||
      length(result[["log_density"]]) != S ||
      length(result[["relative_change"]]) != S) {
    stop("The exact cluster selection kernel returned invalid output.",
         call. = FALSE)
  }
  failed <- which(
    !is.finite(result[["relative_change"]]) |
      result[["relative_change"]] > execution_plan[["relative_tolerance"]]
  )
  if (length(failed) > 0L) {
    stop(
      "Exact selection cluster normalizer was rejected by diagnostics: ",
      "successive quadrature relative change was ",
      format(result[["relative_change"]][failed[[1L]]], digits = 4), ".",
      call. = FALSE
    )
  }
  result[["log_density"]]
}


.selection_exact_factor_loglik_block <- function(
    yi, means, residual_sd, loading, sei, selection_context,
    execution_plan, block_index) {

  S <- nrow(means)
  factor_rank <- execution_plan[["factor_ranks"]][[block_index]]
  design_key  <- execution_plan[["design_keys"]][[block_index]]
  quadrature_orders <- execution_plan[["factor_quadrature_orders"]][[
    as.character(factor_rank)
  ]]
  if (factor_rank < 2L || factor_rank > 4L ||
      ncol(loading) != length(yi) * factor_rank || is.na(design_key) ||
      length(quadrature_orders) < 3L) {
    stop("Internal error: exact selection factor inputs are inconsistent.",
         call. = FALSE)
  }
  quadrature <- .selection_exact_cluster_quadrature_rules(quadrature_orders)
  native_args <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_context,
    S              = S,
    kernel_mode    = selection_context[["kernel_mode"]]
  )
  native_static <- native_args[["static"]]
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
    .native_integer_vector(native_args[["kernel_mode"]]),
    .native_numeric_vector(quadrature[["nodes"]]),
    .native_numeric_vector(quadrature[["log_weights"]]),
    .native_numeric_vector(quadrature[["orders"]]),
    .native_numeric_vector(execution_plan[["designs"]][[design_key]]),
    .native_integer_vector(
      execution_plan[["factor_points_per_proposal"]]
    ),
    .native_integer_vector(
      execution_plan[["factor_max_points_per_proposal"]]
    ),
    .native_integer_vector(execution_plan[["scrambles"]]),
    .native_numeric_vector(execution_plan[["relative_tolerance"]]),
    PACKAGE = "RoBMA"
  )
  expected_names <- c("log_density", "relative_mcse", "relative_change")
  if (!is.list(result) || !identical(names(result), expected_names) ||
      any(lengths(result) != S)) {
    stop("The exact selection factor kernel returned invalid output.",
         call. = FALSE)
  }
  diagnostic <- pmax(result[["relative_mcse"]], result[["relative_change"]])
  failed <- which(!is.finite(diagnostic) |
    diagnostic > execution_plan[["relative_tolerance"]])
  if (length(failed) > 0L) {
    stop(
      "Exact selection factor normalizer was rejected by diagnostics: ",
      "relative Monte Carlo standard error was ",
      format(result[["relative_mcse"]][failed[[1L]]], digits = 4),
      " and nested-design relative change was ",
      format(result[["relative_change"]][failed[[1L]]], digits = 4),
      ". Refit with more integration points in 'selection_control', for ",
      "example by increasing 'max_points_per_scramble' or 'scrambles'.",
      call. = FALSE
    )
  }
  result[["log_density"]]
}


.selection_exact_block_loglik_from_setup <- function(setup) {

  if (!.is_data_exact_selection(setup[["data"]])) {
    stop("Exact selection-likelihood metadata are unavailable.",
         call. = FALSE)
  }
  location <- .estimate_normal_covariance_target_location_from_setup(setup)
  exact_setup <- .data_exact_selection_setup(setup[["data"]])
  selection_context <- .selection_exact_signed_context(
    setup     = setup,
    signed_yi = location[["y"]]
  )
  random_factor <- .selection_exact_random_factor_samples(setup)
  random_covariance <- if (is.null(random_factor)) {
    .selection_exact_random_covariance_samples(setup)
  } else {
    NULL
  }
  covariance_lower <- function(rows, block_index) {

    arguments <- list(
      setup                     = setup,
      rows                      = rows,
      random_covariance_samples = random_covariance
    )
    if (!is.null(random_factor)) {
      arguments[["random_factor_samples"]] <- random_factor
      arguments[["block_index"]] <- block_index
    }
    do.call(.selection_exact_covariance_lower, arguments)
  }
  log_lik <- matrix(
    0,
    nrow = setup[["S"]],
    ncol = length(exact_setup[["row_blocks"]])
  )

  block_sizes      <- lengths(exact_setup[["row_blocks"]])
  singleton_blocks <- which(block_sizes == 1L)
  if (length(singleton_blocks) > 0L) {
    rows <- as.integer(vapply(
      exact_setup[["row_blocks"]][singleton_blocks],
      `[[`,
      integer(1L),
      1L
    ))
    singleton_context <- selection_context
    singleton_context[["obs_bin"]] <- selection_context[["obs_bin"]][rows]
    variances <- .selection_exact_singleton_variances(
      setup                     = setup,
      rows                      = rows,
      block_indices             = singleton_blocks,
      random_covariance_samples = random_covariance,
      random_factor_samples     = random_factor
    )
    log_lik[, singleton_blocks] <- .selection_exact_singleton_loglik_matrix(
      yi                = location[["y"]][rows],
      means             = location[["means"]][, rows, drop = FALSE],
      variances         = variances,
      sei               = setup[["selection_sei"]][rows],
      selection_context = singleton_context
    )
  }
  if (length(singleton_blocks) == length(block_sizes)) {
    return(log_lik)
  }

  for (block_index in which(block_sizes > 1L)) {
    rows <- exact_setup[["row_blocks"]][[block_index]]
    method <- exact_setup[["block_methods"]][[block_index]]
    block_context <- selection_context
    block_context[["obs_bin"]] <- selection_context[["obs_bin"]][rows]
    if (method %in% c("rank_one", "factor")) {
      components <- .selection_exact_factor_block_samples(
        setup                 = setup,
        block_index           = block_index,
        random_factor_samples = random_factor
      )
    }
    if (method == "rank_one") {
      log_lik[, block_index] <- .selection_exact_cluster_loglik_block(
        yi                 = location[["y"]][rows],
        means              = location[["means"]][, rows, drop = FALSE],
        residual_sd        = components[["residual_sd"]],
        loading            = components[["loading"]],
        sei                = setup[["selection_sei"]][rows],
        selection_context  = block_context,
        execution_plan = exact_setup
      )
      next
    }
    if (method == "factor") {
      log_lik[, block_index] <- .selection_exact_factor_loglik_block(
        yi                = location[["y"]][rows],
        means             = location[["means"]][, rows, drop = FALSE],
        residual_sd       = components[["residual_sd"]],
        loading           = components[["loading"]],
        sei               = setup[["selection_sei"]][rows],
        selection_context = block_context,
        execution_plan  = exact_setup,
        block_index       = block_index
      )
      next
    }
    log_lik[, block_index] <- .selection_exact_joint_loglik_block(
      yi               = location[["y"]][rows],
      means            = location[["means"]][, rows, drop = FALSE],
      covariance_lower = covariance_lower(rows, block_index),
      sei               = setup[["selection_sei"]][rows],
      selection_context = block_context,
      execution_plan  = exact_setup,
      block_size        = length(rows)
    )
  }

  log_lik
}


.selection_exact_joint_loglik_from_setup <- function(setup) {

  rowSums(.selection_exact_block_loglik_from_setup(setup))
}


.selection_exact_singleton_random_variance_expression <- function(
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
    stop("Exact selection random-covariance metadata are invalid.",
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


.selection_exact_singleton_model_syntax <- function(
    data, selection_spec, plan) {

  if (length(plan[["singleton_blocks"]]) == 0L) {
    return("")
  }

  random_assignments <- ""
  if (.is_data_random(data)) {
    random_assignments <- paste0(vapply(
      seq_along(plan[["singleton_blocks"]]),
      function(index) {
        block_index <- plan[["singleton_blocks"]][[index]]
        paste0(
          "sel_exact_singleton_random_variance[", index, "] = ",
          .selection_exact_singleton_random_variance_expression(
            plan,
            block_index
          ),
          "\n"
        )
      },
      character(1L)
    ), collapse = "")
  }

  row <- "sel_exact_singleton_row[s]"
  variance_terms <- "sel_exact_singleton_sampling_variance[s]"
  if (.is_data_random(data)) {
    variance_terms <- c(
      variance_terms,
      "sel_exact_singleton_random_variance[s]"
    )
  } else {
    tau_within <- if (.is_data_multilevel(data)) "tau_within" else "tau"
    if (.is_data_scale(data)) {
      tau_within <- paste0(tau_within, "[", row, "]")
    }
    variance_terms <- c(variance_terms, paste0("pow(", tau_within, ",2)"))
    if (.is_data_multilevel(data)) {
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
    "for(s in 1:sel_exact_singleton_n){\n",
    "  sel_exact_singleton_variance[s] = ",
    paste(variance_terms, collapse = " + "), "\n",
    "  sel_exact_singleton_y[s] ~ dselnorm_step_switch(",
    "sel_exact_mu[", row, "],sqrt(sel_exact_singleton_variance[s]),",
    "sel_exact_singleton_sei[s],1,", selection_spec[["jags_omega"]], ",",
    "sel_z_lower,sel_z_upper,sel_exact_singleton_obs_bin[s],sel_sign,",
    .selection_exact_kernel_mode_expression(selection_spec), ",",
    "sel_telescope_probabilities)\n",
    "}\n"
  )
}


.selection_exact_model_syntax <- function(data, selection_spec) {

  plan <- .data_exact_selection_setup(data)
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
    .selection_exact_singleton_model_syntax(data, selection_spec, plan)
  )

  for (block_index in plan[["dependent_blocks"]]) {
    rows      <- plan[["row_blocks"]][[block_index]]
    block_n   <- length(rows)
    method    <- plan[["block_methods"]][[block_index]]
    factor_rank <- plan[["factor_ranks"]][[block_index]]
    lower_n <- if (method %in% c("rank_one", "factor")) {
      0L
    } else {
      length(.selection_exact_lower_pairs(plan, rows)[["row_1"]])
    }
    prefix    <- paste0("sel_exact_block_", block_index)
    design_key <- plan[["design_keys"]][[block_index]]
    qmc_name  <- if (!is.na(design_key)) {
      .selection_exact_qmc_name(design_key)
    } else {
      NULL
    }
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
    uses_diagonal <- .selection_exact_uses_diagonal_indicator(
      data              = data,
      method            = method,
      random_covariance = random_covariance
    )
    if (.is_data_random(data) && uses_diagonal) {
      row_1_expression <- if (method == "singleton") {
        "1"
      } else {
        paste0(prefix, "_local_row_1[l]")
      }
      row_2_expression <- if (method == "singleton") {
        "1"
      } else {
        paste0(prefix, "_local_row_2[l]")
      }
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

    extra_covariance <- if (.is_data_random(data)) {
      "0"
    } else if (.is_data_multilevel(data)) {
      paste0(
        prefix, "_diagonal[l] * pow(tau_within[",
        prefix, "_row_1[l]],2) + tau_between[",
        prefix, "_row_1[l]] * tau_between[",
        prefix, "_row_2[l]]"
      )
    } else {
      paste0(
        prefix, "_diagonal[l] * pow(tau[",
        prefix, "_row_1[l]],2)"
      )
    }
    if (!.is_data_scale(data) && !.is_data_random(data)) {
      extra_covariance <- if (.is_data_multilevel(data)) {
        paste0(
          prefix, "_diagonal[l] * pow(tau_within,2) + ",
          "pow(tau_between,2)"
        )
      } else {
        paste0(prefix, "_diagonal[l] * pow(tau,2)")
      }
    }
    covariance_terms <- c(
      paste0(prefix, "_sampling_lower[l]"),
      random_lower,
      random_factor_covariance,
      extra_covariance
    )
    covariance_terms <- covariance_terms[!is.na(covariance_terms) &
      nzchar(covariance_terms) & covariance_terms != "0"]

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
      if (.is_data_random(data)) {
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
      } else {
        tau_within <- if (.is_data_multilevel(data)) "tau_within" else "tau"
        if (.is_data_scale(data)) {
          tau_within <- paste0(tau_within, "[", prefix, "_row[j]]")
        }
        residual_terms <- c(
          residual_terms,
          paste0("pow(", tau_within, ",2)")
        )
        if (.is_data_multilevel(data)) {
          tau_between <- "tau_between"
          if (.is_data_scale(data)) {
            tau_between <- paste0(
              tau_between, "[", prefix, "_row[j]]"
            )
          }
          loading_expressions <- c(loading_expressions, tau_between)
        }
      }
      if (length(loading_expressions) != factor_rank) {
        stop(
          "Internal error: exact-selection factor syntax rank is inconsistent.",
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
        .selection_exact_kernel_mode_expression(selection_spec), ",",
        "sel_exact_cluster_nodes,sel_exact_cluster_log_weights,",
        "sel_exact_cluster_orders,",
        format(plan[["relative_tolerance"]], scientific = FALSE), ")\n"
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
        .selection_exact_kernel_mode_expression(selection_spec), ",",
        .selection_exact_factor_quadrature_name("nodes", factor_rank), ",",
        .selection_exact_factor_quadrature_name(
          "log_weights", factor_rank
        ), ",",
        .selection_exact_factor_quadrature_name("orders", factor_rank), ",",
        qmc_name, "[1:", plan[["scrambles"]], ",1:",
        plan[["factor_max_points_per_proposal"]], ",1:",
        2L * factor_rank, "],",
        plan[["factor_points_per_proposal"]], ",",
        plan[["factor_max_points_per_proposal"]], ",",
        plan[["scrambles"]], ",",
        format(plan[["relative_tolerance"]], scientific = FALSE), ")\n"
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
        .selection_exact_kernel_mode_expression(selection_spec), ",",
        qmc_name, "[1:", plan[["scrambles"]], ",1:",
        plan[["points_per_scramble"]], ",1:", 2L * block_n, "],",
        plan[["points_per_scramble"]], ",",
        plan[["scrambles"]], ",",
        format(plan[["relative_tolerance"]], scientific = FALSE), ")\n"
      )
    }

    covariance_syntax <- if (method %in% c("rank_one", "factor")) {
      ""
    } else {
      paste0(
        "for(l in 1:", lower_n, "){\n",
        "  ", prefix, "_covariance[l] = ",
        paste(covariance_terms, collapse = " + "), "\n",
        "}\n"
      )
    }
    syntax <- paste0(
      syntax,
      "for(j in 1:", block_n, "){\n",
      "  ", prefix, "_mu[j] = sel_exact_mu[", prefix, "_row[j]]\n",
      "}\n",
      factor_setup_syntax,
      covariance_syntax,
      density_syntax
    )
  }

  syntax
}
