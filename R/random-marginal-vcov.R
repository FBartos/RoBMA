# ============================================================================ #
# random-marginal-vcov.R
# ============================================================================ #
#
# Shared adapter from brma.mv fitted objects to the BayesTools random-effect
# marginal covariance backend.
#
# ============================================================================ #


.brma_mv_random_effects_formula_design <- function(
    object, data = object[["data"]], include_known_group_covariance = TRUE) {

  formula_design <- if (.is_scale(object)) {
    .predict_known_v_formula_design_with_row_source_values(
      object = object,
      data   = data
    )
  } else {
    .fitted_formula_design(object, "mu", required = TRUE)
  }
  if (!include_known_group_covariance) {
    formula_design <- .brma_mv_remove_known_group_covariance(formula_design)
  }

  return(formula_design)
}


.brma_mv_random_effects_block_names <- function(formula_design) {

  random_terms <- formula_design[["random_effects"]]
  if (length(random_terms) == 0L) {
    stop("Random-formula metadata contains no random-effect blocks.",
         call. = FALSE)
  }
  block_names <- vapply(
    random_terms,
    .random_effect_term_block_name,
    character(1)
  )
  if (anyNA(block_names) || any(!nzchar(block_names)) ||
      anyDuplicated(block_names)) {
    stop("Random-formula block names must be non-empty and unique.",
         call. = FALSE)
  }

  return(block_names)
}


.random_effect_dependency_blocks <- function(
    sampling_covariance, formula_design, blocks) {

  if (!is.matrix(sampling_covariance) ||
      !is.numeric(sampling_covariance) ||
      nrow(sampling_covariance) != ncol(sampling_covariance) ||
      anyNA(sampling_covariance) || any(!is.finite(sampling_covariance)) ||
      any(sampling_covariance != t(sampling_covariance))) {
    stop(
      "Sampling covariance is invalid for random-effect dependency ",
      "construction.",
      call. = FALSE
    )
  }
  K         <- nrow(sampling_covariance)
  adjacency <- sampling_covariance != 0
  terms <- if (is.null(formula_design)) NULL else
    formula_design[["random_effects"]]
  term_names <- vapply(terms, .random_effect_term_block_name, character(1))
  terms <- terms[term_names %in% blocks]

  for (term in terms) {
    group_map <- as.integer(term[["group_map"]])
    if (length(group_map) != K || anyNA(group_map) || any(group_map < 1L)) {
      stop(
        "Random-effect grouping metadata are invalid for dependency ",
        "construction.",
        call. = FALSE
      )
    }
    if (.random_effect_term_has_known_group_covariance(term)) {
      kernel <- term[["group_covariance"]][["kernel"]]
      if (!is.matrix(kernel) || any(group_map > nrow(kernel))) {
        stop(
          "Known random-effect group covariance is invalid for dependency ",
          "construction.",
          call. = FALSE
        )
      }
      adjacency <- adjacency | kernel[group_map, group_map, drop = FALSE] != 0
    } else {
      adjacency <- adjacency | outer(group_map, group_map, "==")
    }
  }
  diag(adjacency) <- TRUE

  .known_v_block_indices(adjacency * 1)
}


.brma_mv_random_effects_marginal_inputs <- function(
    object, posterior_samples, data = object[["data"]],
    include_known_group_covariance = TRUE) {

  formula_design <- .brma_mv_random_effects_formula_design(
    object                         = object,
    data                           = data,
    include_known_group_covariance = include_known_group_covariance
  )
  formula_fit <- .posterior_formula_fit(
    fit               = object[["fit"]],
    posterior_samples = posterior_samples,
    formula_design    = TRUE
  )
  attr(formula_fit, "formula_design") <- list(mu = formula_design)

  location_priors <- attr(object[["fit"]], "prior_list")
  if (is.null(location_priors)) {
    location_priors <- formula_design[["prior_list"]]
  }
  if (is.null(location_priors)) {
    location_priors <- object[["priors"]][["location"]]
  }

  return(list(
    formula_fit     = formula_fit,
    formula_design  = formula_design,
    location_priors = location_priors
  ))
}


.brma_mv_random_effects_marginal_vcov <- function(
    object, posterior_samples, blocks = NULL, diagonal_only = FALSE,
    data = object[["data"]], new_levels = NULL,
    include_known_group_covariance = TRUE, inputs = NULL) {

  if (is.null(inputs)) {
    inputs <- .brma_mv_random_effects_marginal_inputs(
      object                         = object,
      posterior_samples              = posterior_samples,
      data                           = data,
      include_known_group_covariance = include_known_group_covariance
    )
  }

  return(BayesTools::random_effects_marginal_vcov(
    fit               = inputs[["formula_fit"]],
    parameter         = "mu",
    data              = data[["location"]],
    posterior_samples = posterior_samples,
    prior_list        = inputs[["location_priors"]],
    blocks            = blocks,
    new_levels        = new_levels,
    diagonal_only     = diagonal_only
  ))
}


.brma_mv_random_effects_marginal_factor_states <- function(
    object, posterior_samples, blocks, row_blocks,
    data = object[["data"]], inputs = NULL,
    include_known_group_covariance = TRUE) {

  if (is.null(inputs)) {
    inputs <- .brma_mv_random_effects_marginal_inputs(
      object                         = object,
      posterior_samples              = posterior_samples,
      data                            = data,
      include_known_group_covariance = include_known_group_covariance
    )
  }

  return(BayesTools::random_effects_marginal_factor_states(
    fit               = inputs[["formula_fit"]],
    parameter         = "mu",
    posterior_samples = posterior_samples,
    prior_list        = inputs[["location_priors"]],
    blocks            = blocks,
    row_blocks        = row_blocks
  ))
}


.brma_mv_random_effects_marginal_diagonal_by_block <- function(
    object, posterior_samples, blocks, data = object[["data"]],
    include_known_group_covariance = TRUE) {

  inputs <- .brma_mv_random_effects_marginal_inputs(
    object                         = object,
    posterior_samples              = posterior_samples,
    data                            = data,
    include_known_group_covariance = include_known_group_covariance
  )
  random_terms <- inputs[["formula_design"]][["random_effects"]]
  n_rows       <- nrow(random_terms[[1L]][["model_matrix"]])
  factors      <- .brma_mv_random_effects_marginal_factor_states(
    object            = object,
    posterior_samples = posterior_samples,
    blocks            = blocks,
    row_blocks        = list(seq_len(n_rows)),
    data              = data,
    inputs            = inputs
  )

  variance <- BayesTools::random_effects_marginal_factor_diagonal(
    factors,
    by_block = TRUE
  )
  row_names <- rownames(random_terms[[1L]][["model_matrix"]])
  if (is.null(row_names)) {
    row_names <- as.character(seq_len(n_rows))
  }
  for (block in names(variance)) {
    colnames(variance[[block]]) <- row_names
    names(dimnames(variance[[block]])) <- c("draw", "row")
  }

  return(variance)
}


.brma_mv_random_effects_marginal_factor_plan <- function(
    object, posterior_samples, blocks = NULL, row_blocks = NULL,
    data = object[["data"]], sampling_latent_marginalized = TRUE) {

  inputs <- .brma_mv_random_effects_marginal_inputs(
    object            = object,
    posterior_samples = posterior_samples,
    data              = data
  )
  available_blocks <- .brma_mv_random_effects_block_names(
    inputs[["formula_design"]]
  )
  if (is.null(blocks)) {
    blocks <- available_blocks
  }
  if (!is.character(blocks) || length(blocks) == 0L || anyNA(blocks) ||
      any(!nzchar(blocks)) || anyDuplicated(blocks) ||
      any(!blocks %in% available_blocks)) {
    stop("Requested random-formula blocks are invalid.", call. = FALSE)
  }
  if (is.null(row_blocks)) {
    row_blocks <- .random_effect_dependency_blocks(
      sampling_covariance = .known_v_dependency_covariance(
        data,
        sampling_latent_marginalized = sampling_latent_marginalized
      ),
      formula_design = inputs[["formula_design"]],
      blocks         = blocks
    )
  }

  random_factors <- .brma_mv_random_effects_marginal_factor_states(
    object            = object,
    posterior_samples = posterior_samples,
    blocks            = blocks,
    row_blocks        = row_blocks,
    data              = data,
    inputs            = inputs
  )
  if (!identical(random_factors[["row_blocks"]], row_blocks) ||
      !identical(random_factors[["metadata"]][["included_blocks"]], blocks) ||
      length(random_factors[["factor_states"]]) != nrow(posterior_samples)) {
    stop("Random-effect covariance returned invalid factor states.",
         call. = FALSE)
  }

  return(random_factors)
}


.brma_mv_remove_known_group_covariance <- function(formula_design) {

  terms <- formula_design[["random_effects"]]
  for (i in seq_along(terms)) {
    if (!.random_effect_term_has_known_group_covariance(terms[[i]])) {
      next
    }

    terms[[i]][["group_covariance"]]         <- NULL
    terms[[i]][["marginal_variance_factor"]] <- NULL
    terms[[i]][["row_multiplier"]]           <- NULL
    terms[[i]][["row_multiplier_name"]]      <- NULL
  }

  formula_design[["random_effects"]] <- terms
  return(formula_design)
}
