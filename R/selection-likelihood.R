# ============================================================================ #
# selection-likelihood.R
# ============================================================================ #
#
# Exact finite-vector selection likelihood preparation. Generic numerical and
# random-covariance compilation is owned by BayesTools; this file supplies the
# meta-analytic covariance and constructor integration.
#
# ============================================================================ #


SELNORM_CLUSTER_QUADRATURE_ORDERS <- c(31L, 63L, 127L, 255L, 511L)


#' Control numerical integration of exact selection likelihoods
#'
#' @description
#' Creates numerical integration settings for the finite-vector product
#' selection likelihood used by [bselmodel()], [bselmodel.mv()], [RoBMA()], and
#' [RoBMA.mv()]. Rank-one multilevel blocks use a deterministic sequence of
#' Gauss-Hermite rules; general covariance blocks use randomized quasi-Monte
#' Carlo. Every integration design is fixed before fitting, so repeated
#' likelihood evaluations are deterministic. One-dimensional likelihood blocks
#' are evaluated analytically.
#'
#' @param points_per_scramble number of shifted Halton points per randomized
#'   quasi-Monte Carlo scramble for general covariance blocks.
#' @param scrambles number of independent randomized shifts used to estimate
#'   integration error for general covariance blocks.
#' @param relative_tolerance largest accepted successive relative quadrature
#'   change for rank-one multilevel blocks or relative Monte Carlo standard
#'   error for general covariance blocks.
#' @param seed non-negative integer used only to create the fixed integration
#'   design for general covariance blocks. It does not alter R's random-number
#'   state.
#'
#' @return A `RoBMA_selection_likelihood_control` object.
#'
#' @export
set_selection_likelihood_control <- function(
    points_per_scramble = 512L, scrambles = 8L,
    relative_tolerance = 0.005, seed = 1L) {

  BayesTools::check_int(
    points_per_scramble,
    "points_per_scramble",
    lower = 8L,
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
      points_per_scramble = as.integer(points_per_scramble),
      scrambles            = as.integer(scrambles),
      relative_tolerance   = as.numeric(relative_tolerance),
      seed                 = as.integer(seed)
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

  setup <- attr(data, "exact_selection_likelihood", exact = TRUE)
  if (.is_data_exact_selection(data) && is.null(setup)) {
    stop(
      "Internal error: exact selection-likelihood metadata are missing.",
      call. = FALSE
    )
  }
  setup
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

  row_blocks <- .selection_exact_dependency_blocks(
    data         = object[["data"]],
    random_terms = if (is.null(formula_design)) {
      list()
    } else {
      formula_design[["random_effects"]]
    }
  )
  block_sizes       <- lengths(row_blocks)
  closed_form       <- all(block_sizes == 1L)
  cluster_reduction <- !closed_form && .selection_exact_uses_cluster_reduction(
    object[["data"]]
  )
  plan <- if (closed_form) {
    .selection_exact_closed_form_plan(block_sizes)
  } else if (cluster_reduction) {
    .selection_exact_cluster_plan(
      block_sizes        = block_sizes,
      relative_tolerance = selection_control[["relative_tolerance"]]
    )
  } else {
    BayesTools::selection_likelihood_plan(
      block_sizes         = block_sizes,
      points_per_scramble = selection_control[["points_per_scramble"]],
      scrambles            = selection_control[["scrambles"]],
      seed                 = selection_control[["seed"]],
      relative_tolerance   = selection_control[["relative_tolerance"]]
    )
  }
  random_covariance <- NULL
  if (!is.null(formula_design)) {
    random_covariance <- BayesTools::JAGS_formula_random_marginal_covariance(
      formula_design = formula_design,
      row_blocks     = row_blocks,
      prefix         = "sel_exact_random"
    )
  }
  setup <- list(
    schema_version    = 1L,
    target            = "finite_vector_product_selection",
    exactness         = plan[["exactness"]],
    row_blocks        = row_blocks,
    sampling_covariance = .selection_exact_sampling_covariance(
      object[["data"]]
    ),
    integration_plan  = plan,
    random_covariance = random_covariance
  )
  class(setup) <- c("RoBMA_exact_selection_likelihood", "list")
  object[["data"]] <- .set_data_selection_likelihood(
    data       = object[["data"]],
    likelihood = "exact",
    setup      = setup
  )
  object[["selection_likelihood"]] <- c(
    list(type = "exact"),
    setup[c("target", "exactness", "row_blocks", "integration_plan")]
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


.selection_exact_sampling_covariance <- function(data) {

  if (.is_data_known_v(data)) {
    return(.known_v_covariance_matrix(.data_known_v_data(data)))
  }
  diag(data[["outcome"]][["sei"]]^2)
}


.selection_exact_dependency_blocks <- function(data, random_terms = list()) {

  if (!is.list(random_terms)) {
    stop("Random-effect dependency metadata must be a list.",
         call. = FALSE)
  }
  if (length(random_terms) > 0L) {
    block_names <- unique(vapply(
      random_terms,
      .random_effect_term_block_name,
      character(1)
    ))
    return(.random_effect_dependency_blocks(
      sampling_covariance = .selection_exact_sampling_covariance(data),
      formula_design      = list(random_effects = random_terms),
      blocks              = block_names
    ))
  }

  sampling_covariance <- .selection_exact_sampling_covariance(data)
  adjacency <- sampling_covariance != 0

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


.selection_exact_fit_data <- function(data, priors) {

  setup            <- .data_exact_selection_setup(data)
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

  plan <- setup[["integration_plan"]]
  cluster_reduction <- .selection_exact_uses_cluster_reduction(data)
  if (cluster_reduction) {
    quadrature <- .selection_exact_cluster_quadrature_rules(
      plan[["quadrature_orders"]]
    )
    fit_data[["sel_exact_cluster_nodes"]] <- quadrature[["nodes"]]
    fit_data[["sel_exact_cluster_log_weights"]] <-
      quadrature[["log_weights"]]
    fit_data[["sel_exact_cluster_orders"]] <- quadrature[["orders"]]
  }
  dependent_sizes <- if (cluster_reduction) {
    integer()
  } else {
    sort(unique(
      lengths(setup[["row_blocks"]])[lengths(setup[["row_blocks"]]) > 1L]
    ))
  }
  for (block_n in dependent_sizes) {
    fit_data[[.selection_exact_qmc_name(block_n)]] <- plan[["designs"]][[
      as.character(block_n)
    ]]
  }
  sampling_covariance <- setup[["sampling_covariance"]]
  for (block_index in seq_along(setup[["row_blocks"]])) {
    rows  <- setup[["row_blocks"]][[block_index]]
    prefix <- paste0("sel_exact_block_", block_index)
    fit_data[[paste0(prefix, "_y")]]       <- yi[rows]
    fit_data[[paste0(prefix, "_sei")]]     <- sei[rows]
    fit_data[[paste0(prefix, "_obs_bin")]] <- selection_spec[["obs_bin"]][rows]
    fit_data[[paste0(prefix, "_row")]]     <- rows
    if (cluster_reduction && length(rows) > 1L) {
      fit_data[[paste0(prefix, "_sampling_sd")]] <- sei[rows]
    }
    if (!cluster_reduction || length(rows) == 1L) {
      pairs <- .selection_exact_lower_pairs(plan, rows)
    }
    if (!.is_data_random(data) &&
        (!cluster_reduction || length(rows) == 1L)) {
      fit_data[[paste0(prefix, "_diagonal")]] <- as.integer(
        pairs[["row_1"]] == pairs[["row_2"]]
      )
      if (.is_data_scale(data)) {
        fit_data[[paste0(prefix, "_row_1")]] <- pairs[["row_1"]]
        if (.is_data_multilevel(data)) {
          fit_data[[paste0(prefix, "_row_2")]] <- pairs[["row_2"]]
        }
      }
    }
    if (!cluster_reduction || length(rows) == 1L) {
      fit_data[[paste0(prefix, "_sampling_lower")]] <- unname(
        sampling_covariance[cbind(pairs[["row_1"]], pairs[["row_2"]])]
      )
    }
  }
  if (!is.null(setup[["random_covariance"]])) {
    fit_data <- c(fit_data, setup[["random_covariance"]][["data"]])
  }

  fit_data
}


.selection_exact_uses_cluster_reduction <- function(data) {

  setup <- attr(data, "exact_selection_likelihood", exact = TRUE)
  if (!is.null(setup)) {
    return(identical(setup[["exactness"]], "E1"))
  }
  .is_data_multilevel(data) && !.is_data_random(data) &&
    !.is_data_known_v(data)
}


.selection_exact_closed_form_plan <- function(block_sizes) {

  structure(list(
    schema_version     = 1L,
    block_sizes        = as.integer(block_sizes),
    lower_pairs        = list(
      "1" = data.frame(row_1 = 1L, row_2 = 1L)
    ),
    exactness          = "E0",
    statistical_target = "finite_vector_product_selection"
  ), class = c("RoBMA_selection_closed_form_plan", "list"))
}


.selection_exact_cluster_plan <- function(block_sizes, relative_tolerance) {

  structure(list(
    schema_version     = 1L,
    block_sizes        = as.integer(block_sizes),
    relative_tolerance = as.numeric(relative_tolerance),
    quadrature_orders  = SELNORM_CLUSTER_QUADRATURE_ORDERS,
    lower_pairs        = list(
      "1" = data.frame(row_1 = 1L, row_2 = 1L)
    ),
    exactness          = "E1",
    statistical_target = "finite_vector_product_selection"
  ), class = c("RoBMA_selection_cluster_plan", "list"))
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


.selection_exact_covariance_lower <- function(
    setup, rows, random_covariance_samples = NULL) {

  exact_setup <- .data_exact_selection_setup(setup[["data"]])
  sampling    <- exact_setup[["sampling_covariance"]][rows, rows, drop = FALSE]
  block_n     <- length(rows)
  pairs       <- .selection_exact_lower_pairs(
    exact_setup[["integration_plan"]],
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
    yi, means, covariance_lower, sei, selection_context, integration_plan,
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
    .native_numeric_vector(integration_plan[["designs"]][[
      as.character(block_size)
    ]]),
    .native_integer_vector(integration_plan[["points_per_scramble"]]),
    .native_integer_vector(integration_plan[["scrambles"]]),
    .native_numeric_vector(integration_plan[["relative_tolerance"]]),
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
      result[["relative_mcse"]] > integration_plan[["relative_tolerance"]]
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
    integration_plan) {

  S <- nrow(means)
  native_args <- BayesTools::selection_native_kernel_args(
    selection_spec = selection_context,
    S              = S,
    kernel_mode    = selection_context[["kernel_mode"]]
  )
  native_static <- native_args[["static"]]
  quadrature <- .selection_exact_cluster_quadrature_rules(
    integration_plan[["quadrature_orders"]]
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
    .native_numeric_vector(integration_plan[["relative_tolerance"]]),
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
      result[["relative_change"]] > integration_plan[["relative_tolerance"]]
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
  random_covariance <- .selection_exact_random_covariance_samples(setup)
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
    variances <- do.call(cbind, lapply(rows, function(row) {
      .selection_exact_covariance_lower(
        setup                     = setup,
        rows                      = row,
        random_covariance_samples = random_covariance
      )[, 1L]
    }))
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
    block_context <- selection_context
    block_context[["obs_bin"]] <- selection_context[["obs_bin"]][rows]
    if (.selection_exact_uses_cluster_reduction(setup[["data"]]) &&
        length(rows) > 1L) {
      residual_sd <- sqrt(
        setup[["tau_within"]][, rows, drop = FALSE]^2 +
          matrix(
            setup[["selection_sei"]][rows]^2,
            nrow = setup[["S"]],
            ncol = length(rows),
            byrow = TRUE
          )
      )
      log_lik[, block_index] <- .selection_exact_cluster_loglik_block(
        yi                 = location[["y"]][rows],
        means              = location[["means"]][, rows, drop = FALSE],
        residual_sd        = residual_sd,
        loading            = setup[["tau_between"]][, rows, drop = FALSE],
        sei                = setup[["selection_sei"]][rows],
        selection_context  = block_context,
        integration_plan = exact_setup[["integration_plan"]]
      )
      next
    }
    log_lik[, block_index] <- .selection_exact_joint_loglik_block(
      yi               = location[["y"]][rows],
      means            = location[["means"]][, rows, drop = FALSE],
      covariance_lower = .selection_exact_covariance_lower(
        setup                     = setup,
        rows                      = rows,
        random_covariance_samples = random_covariance
      ),
      sei               = setup[["selection_sei"]][rows],
      selection_context = block_context,
      integration_plan  = exact_setup[["integration_plan"]],
      block_size        = length(rows)
    )
  }

  log_lik
}


.selection_exact_joint_loglik_from_setup <- function(setup) {

  rowSums(.selection_exact_block_loglik_from_setup(setup))
}


.selection_exact_model_syntax <- function(data, selection_spec) {

  setup <- .data_exact_selection_setup(data)
  plan  <- setup[["integration_plan"]]
  cluster_reduction <- .selection_exact_uses_cluster_reduction(data)
  syntax <- ""
  if (!is.null(setup[["random_covariance"]])) {
    syntax <- paste0(syntax, setup[["random_covariance"]][["syntax"]])
  }

  for (block_index in seq_along(setup[["row_blocks"]])) {
    rows      <- setup[["row_blocks"]][[block_index]]
    block_n   <- length(rows)
    lower_n <- if (cluster_reduction && block_n > 1L) {
      0L
    } else {
      length(.selection_exact_lower_pairs(plan, rows)[["row_1"]])
    }
    prefix    <- paste0("sel_exact_block_", block_index)
    qmc_name  <- if (block_n > 1L) {
      .selection_exact_qmc_name(block_n)
    } else {
      NULL
    }
    random_lower <- if (is.null(setup[["random_covariance"]])) {
      NULL
    } else {
      paste0(
        setup[["random_covariance"]][["lower_names"]][[block_index]],
        "[l]"
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
      extra_covariance
    )
    covariance_terms <- covariance_terms[!is.na(covariance_terms) &
      nzchar(covariance_terms) & covariance_terms != "0"]

    density_syntax <- if (block_n == 1L) {
      paste0(
        prefix, "_y[1] ~ dselnorm_step_switch(",
        prefix, "_mu[1],sqrt(", prefix, "_covariance[1]),",
        prefix, "_sei[1],1,",
        selection_spec[["jags_omega"]], ",",
        "sel_z_lower,sel_z_upper,",
        prefix, "_obs_bin[1],",
        "sel_sign,",
        .selection_exact_kernel_mode_expression(selection_spec), ",",
        "sel_telescope_probabilities)\n"
      )
    } else if (cluster_reduction) {
      tau_within <- if (.is_data_scale(data)) {
        paste0("tau_within[", prefix, "_row[j]]")
      } else {
        "tau_within"
      }
      tau_between <- if (.is_data_scale(data)) {
        paste0("tau_between[", prefix, "_row[j]]")
      } else {
        "tau_between"
      }
      paste0(
        "for(j in 1:", block_n, "){\n",
        "  ", prefix, "_residual_sd[j] = sqrt(pow(",
        prefix, "_sampling_sd[j],2) + pow(", tau_within, ",2))\n",
        "  ", prefix, "_loading[j] = ", tau_between, "\n",
        "}\n",
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

    covariance_syntax <- if (cluster_reduction && block_n > 1L) {
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
      covariance_syntax,
      density_syntax
    )
  }

  syntax
}
