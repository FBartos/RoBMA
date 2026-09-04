# ============================================================================ #
# Approximate selection-likelihood diagnostics
# ============================================================================ #


#' Diagnose the Conditional Selection-Likelihood Approximation
#'
#' @description
#' Quantifies how strongly estimate-level selection would reweight the latent
#' Gaussian effects that are conditioned upon by an approximate multivariate
#' selection likelihood.
#'
#' @param object A fitted [bselmodel.mv()] object using
#'   `selection_likelihood = "approximate"`.
#' @param max_posterior_samples Maximum number of posterior draws to evaluate.
#'   Draws are deterministically thinned over posterior row order. Defaults to
#'   `256`.
#' @param latent_samples Number of fresh pre-selection Gaussian latent draws per
#'   posterior draw. Defaults to `512`.
#' @param seed Integer seed for the diagnostic simulation. The calling R
#'   session's random-number state is restored before returning.
#'
#' @details
#' For posterior draw \eqn{\theta} and dependency block \eqn{b}, the diagnostic
#' simulates fresh Gaussian effects \eqn{g} from their pre-selection
#' distribution and evaluates
#' \deqn{W_b(g; \theta) = \prod_{i \in b} A_i(g; \theta),}
#' where \eqn{A_i} is the exact univariate selection normalizer used by the
#' fitted conditional likelihood. The known sampling-covariance factors and all
#' sampled random-effect blocks are obtained from the same compiled covariance
#' metadata used by post-fit likelihood calculations. Estimate-level effects
#' marginalized during fitting remain in the conditional row variance.
#'
#' `ess_fraction` is the population importance-sampling ESS fraction estimated
#' by `(mean(W)^2 / mean(W^2))`. `total_variation` estimates the total-variation
#' distance between the simulated pre-selection latent distribution and its
#' normalized `W`-reweighted version. Values near one for `ess_fraction` and
#' near zero for `total_variation` indicate weak latent reweighting. The
#' diagnostic is descriptive: it is neither an error bound nor a comparison to
#' a fitted exact model. Because \eqn{\theta} is sampled from the approximate
#' posterior, the diagnostic cannot reveal posterior regions that the
#' approximation itself failed to explore.
#'
#' The simulation budgets control diagnostic Monte Carlo precision, not the
#' quadrature or integration rule of either fitted likelihood. Monte Carlo
#' standard errors are reported so that a budget can be assessed. The largest
#' block-level posterior-median total-variation distance controls notification:
#' values up to `0.05` are silent, values above `0.05` and below `0.10` produce
#' a message, and values of at least `0.10` produce a warning. These are
#' operational notification boundaries, not universal statistical adequacy
#' standards.
#'
#' @return A data frame with one row per connected dependency block. The
#'   `q05`, `median`, and `q95` columns summarize the posterior distribution of
#'   the ESS fraction, total-variation distance, and interquartile range of
#'   `log(W)`. The `mcse_median` and `mcse_max` columns summarize conditional
#'   Monte Carlo standard errors across posterior draws. Remaining columns
#'   identify the block and record the posterior and latent simulation budgets.
#'
#' @seealso [bselmodel.mv()], [set_selection_likelihood_control()]
#'
#' @export
selection_approximation_diagnostics <- function(
    object, max_posterior_samples = 256L, latent_samples = 512L,
    seed = 1L) {

  if (!inherits(object, "bselmodel.mv")) {
    stop(
      "'object' must be a fitted 'bselmodel.mv' object.",
      call. = FALSE
    )
  }
  if (!identical(.data_selection_likelihood(object[["data"]]), "approximate")) {
    stop(
      "'selection_approximation_diagnostics()' requires a model fitted with ",
      "'selection_likelihood = \"approximate\"'.",
      call. = FALSE
    )
  }

  max_posterior_samples <- .normalize_max_samples(
    max_posterior_samples,
    argument = "max_posterior_samples"
  )
  latent_samples <- .normalize_max_samples(
    latent_samples,
    argument = "latent_samples",
    minimum = 2L
  )
  if (is.infinite(latent_samples)) {
    stop("'latent_samples' must be a finite positive integer.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != floor(seed)) {
    stop("'seed' must be a single non-negative integer.", call. = FALSE)
  }
  seed <- as.integer(seed)

  rng_kind <- RNGkind()
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    do.call(RNGkind, as.list(rng_kind))
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)

  sample_info <- .known_v_diagnostic_posterior_samples(
    object            = object,
    max_samples       = max_posterior_samples,
    caller            = "selection_approximation_diagnostics()",
    warn              = FALSE
  )
  posterior_samples <- sample_info[["posterior_samples"]]
  setup <- .estimate_likelihood_setup.brma(
    object                  = object,
    posterior_samples       = posterior_samples,
    condition_local_effects = FALSE
  )
  selection_context <- .selection_context(object, posterior_samples)
  if (is.null(selection_context)) {
    stop(
      "Selection-approximation diagnostics are unavailable without a ",
      "selection-model prior.",
      call. = FALSE
    )
  }
  .selection_require_step_evaluable(
    selection_context,
    "selection_approximation_diagnostics()"
  )

  known_V <- .data_known_v_data(object[["data"]])
  plan <- .known_v_marginal_factor_plan(
    object            = object,
    posterior_samples = posterior_samples,
    known_V           = known_V,
    extra_variances   = NULL
  )
  residual_covariance <- plan[["sampling_covariance"]]
  off_diagonal <- row(residual_covariance) != col(residual_covariance)
  if (any(residual_covariance[off_diagonal] != 0)) {
    stop(
      "Selection-approximation diagnostics are unavailable because the ",
      "conditional sampling covariance is not diagonal.",
      call. = FALSE
    )
  }

  residual_variance <- sweep(
    plan[["extra_variances"]],
    2L,
    diag(residual_covariance),
    "+"
  )
  if (any(!is.finite(residual_variance)) || any(residual_variance <= 0)) {
    stop(
      "Selection-approximation diagnostics require finite positive ",
      "conditional row variances.",
      call. = FALSE
    )
  }

  result <- .selection_approximation_run(
    means               = setup[["mu"]],
    residual_variance   = residual_variance,
    selection_sei       = setup[["selection_sei"]],
    selection_context   = selection_context,
    factor_plans        = plan[["random_covariance_plans"]],
    factor_states       = plan[["random_covariance_states"]],
    row_blocks          = plan[["block_indices"]],
    latent_samples      = latent_samples
  )

  result[["posterior_samples"]]       <- sample_info[["n_used"]]
  result[["total_posterior_samples"]] <- sample_info[["n_total"]]
  result[["latent_samples"]]          <- latent_samples
  .selection_approximation_notify(result[["total_variation_median"]])
  result
}


.selection_approximation_notify <- function(total_variation) {

  maximum_tv <- max(total_variation)
  tv_percent <- sprintf("%.1f", 100 * maximum_tv)
  if (maximum_tv >= 0.10) {
    warning(
      "The approximate selection likelihood showed substantial ",
      "latent-distribution reweighting: the largest block-level ",
      "posterior-median total variation distance was ", tv_percent, "%. ",
      "The approximate and exact likelihoods may yield different inference. ",
      "Consider refitting with 'selection_likelihood' set to 'exact'.",
      call. = FALSE
    )
  } else if (maximum_tv > 0.05) {
    message(
      "Approximate selection-likelihood diagnostic: the largest block-level ",
      "posterior-median total variation distance was ", tv_percent, "%. ",
      "Inspect 'selection_approximation_diagnostics()' for block-level ",
      "results."
    )
  }
  invisible(NULL)
}


.selection_approximation_run <- function(
    means, residual_variance, selection_sei, selection_context,
    factor_plans, factor_states, row_blocks, latent_samples) {

  means             <- as.matrix(means)
  residual_variance <- as.matrix(residual_variance)
  S                 <- nrow(means)
  K                 <- ncol(means)
  B                 <- length(row_blocks)

  if (!identical(dim(residual_variance), c(S, K)) ||
      length(selection_sei) != K || !is.list(factor_plans) ||
      !is.list(factor_states) || length(factor_states) != S || B == 0L ||
      !identical(
        sort(as.integer(unlist(row_blocks, use.names = FALSE))),
        seq_len(K)
      )) {
    stop(
      "Selection-approximation diagnostic inputs are inconsistent.",
      call. = FALSE
    )
  }

  metric_names <- c(
    "ess_fraction",
    "ess_fraction_mcse",
    "total_variation",
    "total_variation_mcse",
    "log_weight_iqr"
  )
  metrics <- stats::setNames(
    lapply(metric_names, function(x) matrix(NA_real_, nrow = S, ncol = B)),
    metric_names
  )

  for (draw in seq_len(S)) {
    latent <- .selection_approximation_latent_draws(
      factor_plans   = factor_plans,
      factor_states  = factor_states[[draw]],
      row_blocks     = row_blocks,
      latent_samples = latent_samples,
      K              = K
    )
    conditional_means <- sweep(latent, 2L, means[draw, ], "+")
    conditional_sd <- matrix(
      sqrt(residual_variance[draw, ]),
      nrow  = latent_samples,
      ncol  = K,
      byrow = TRUE
    )
    draw_context <- BayesTools::selection_context_subset_rows(
      selection_context,
      rep.int(draw, latent_samples)
    )
    log_normalizers <- .selection_step_log_norm_matrix(
      mean              = conditional_means,
      sd                = conditional_sd,
      sei               = selection_sei,
      selection_context = draw_context
    )
    if (!identical(dim(log_normalizers), c(latent_samples, K)) ||
        any(!is.finite(log_normalizers))) {
      stop(
        "Selection-approximation normalizers are invalid.",
        call. = FALSE
      )
    }

    for (block in seq_along(row_blocks)) {
      block_metrics <- .selection_approximation_weight_metrics(
        rowSums(log_normalizers[, row_blocks[[block]], drop = FALSE])
      )
      for (metric in metric_names) {
        metrics[[metric]][draw, block] <- block_metrics[[metric]]
      }
    }
  }

  .selection_approximation_summarize(metrics, row_blocks)
}


.selection_approximation_latent_draws <- function(
    factor_plans, factor_states, row_blocks, latent_samples, K) {

  if (length(factor_plans) != length(factor_states)) {
    stop(
      "Selection-approximation covariance factor states are inconsistent.",
      call. = FALSE
    )
  }
  latent <- matrix(0, nrow = latent_samples, ncol = K)
  for (factor_index in seq_along(factor_plans)) {
    plan  <- factor_plans[[factor_index]]
    state <- factor_states[[factor_index]]
    basis <- if (identical(plan[["type"]], "dense")) {
      NULL
    } else {
      .selection_approximation_factor_basis(plan, state, K)
    }
    for (rows in row_blocks) {
      latent[, rows] <- latent[, rows, drop = FALSE] +
        .selection_approximation_factor_draws(
          plan           = plan,
          state          = state,
          rows           = rows,
          latent_samples = latent_samples,
          K              = K,
          basis          = basis
        )
    }
  }
  latent
}


.selection_approximation_factor_draws <- function(
    plan, state, rows, latent_samples, K, basis = NULL) {

  type <- plan[["type"]]
  if (identical(type, "dense")) {
    covariance <- state[["covariance"]]
    if (!is.matrix(covariance) || nrow(covariance) != K ||
        ncol(covariance) != K) {
      stop(
        "Selection-approximation dense covariance state is invalid.",
        call. = FALSE
      )
    }
    return(.selection_approximation_covariance_draws(
      covariance     = covariance[rows, rows, drop = FALSE],
      latent_samples = latent_samples
    ))
  }
  if (!type %in% c("group", "row_group", "known_group")) {
    stop(
      "Selection-approximation covariance factor type is unsupported.",
      call. = FALSE
    )
  }

  if (is.null(basis)) {
    basis <- .selection_approximation_factor_basis(plan, state, K)
  }
  group_map <- plan[["group_map"]]
  if (length(group_map) != K) {
    stop(
      "Selection-approximation covariance group mapping is invalid.",
      call. = FALSE
    )
  }
  out <- matrix(0, nrow = latent_samples, ncol = length(rows))

  if (!identical(type, "known_group")) {
    for (group in unique(group_map[rows])) {
      local_rows <- which(group_map[rows] == group)
      active <- which(colSums(basis[rows[local_rows], , drop = FALSE] != 0) > 0)
      if (length(active) == 0L) {
        next
      }
      coefficients <- matrix(
        stats::rnorm(latent_samples * length(active)),
        nrow = latent_samples,
        ncol = length(active)
      )
      out[, local_rows] <- coefficients %*%
        t(basis[rows[local_rows], active, drop = FALSE])
    }
    return(out)
  }

  active_rows <- which(rowSums(basis[rows, , drop = FALSE] != 0) > 0)
  if (length(active_rows) == 0L) {
    return(out)
  }
  active_columns <- which(
    colSums(basis[rows[active_rows], , drop = FALSE] != 0) > 0
  )
  groups <- unique(group_map[rows[active_rows]])
  group_covariance <- plan[["group_covariance"]][
    groups,
    groups,
    drop = FALSE
  ]
  group_factor <- .covariance_sampling_factor(
    .covariance_factorization(group_covariance)
  )
  if (is.null(group_factor)) {
    stop(
      "Selection-approximation known group covariance is not positive ",
      "semidefinite.",
      call. = FALSE
    )
  }
  coefficients <- array(
    0,
    dim = c(latent_samples, length(groups), length(active_columns))
  )
  for (column in seq_along(active_columns)) {
    coefficients[, , column] <- matrix(
      stats::rnorm(latent_samples * length(groups)),
      nrow = latent_samples,
      ncol = length(groups)
    ) %*% group_factor
  }
  for (group_index in seq_along(groups)) {
    local_rows <- which(group_map[rows] == groups[[group_index]])
    group_coefficients <- matrix(
      coefficients[, group_index, , drop = FALSE],
      nrow = latent_samples,
      ncol = length(active_columns)
    )
    out[, local_rows] <- group_coefficients %*%
      t(basis[rows[local_rows], active_columns, drop = FALSE])
  }
  out
}


.selection_approximation_factor_basis <- function(plan, state, K) {

  basis <- plan[["model_matrix"]] %*% state[["coefficient_factor"]]
  if (identical(plan[["type"]], "row_group")) {
    basis <- basis * state[["row_scale"]]
  }
  if (!is.matrix(basis) || nrow(basis) != K || any(!is.finite(basis))) {
    stop(
      "Selection-approximation covariance factor state is invalid.",
      call. = FALSE
    )
  }
  basis
}


.selection_approximation_covariance_draws <- function(
    covariance, latent_samples) {

  if (all(covariance == 0)) {
    return(matrix(0, nrow = latent_samples, ncol = nrow(covariance)))
  }
  factor <- .covariance_sampling_factor(
    .covariance_factorization(covariance)
  )
  if (is.null(factor)) {
    stop(
      "Selection-approximation latent covariance is not positive semidefinite.",
      call. = FALSE
    )
  }
  matrix(
    stats::rnorm(latent_samples * nrow(covariance)),
    nrow = latent_samples,
    ncol = nrow(covariance)
  ) %*% factor
}


.selection_approximation_weight_metrics <- function(log_weights) {

  if (!is.numeric(log_weights) || length(log_weights) < 2L ||
      any(!is.finite(log_weights))) {
    stop(
      "Selection-approximation log weights are invalid.",
      call. = FALSE
    )
  }
  weights <- exp(log_weights - max(log_weights))
  mean_weight <- mean(weights)
  mean_square <- mean(weights^2)
  normalized  <- weights / mean_weight

  ess_fraction <- mean_weight^2 / mean_square
  ess_influence <-
    (2 * mean_weight / mean_square) * (weights - mean_weight) -
    (mean_weight^2 / mean_square^2) * (weights^2 - mean_square)
  ess_mcse <- stats::sd(ess_influence) / sqrt(length(weights))

  absolute_difference <- abs(normalized - 1)
  total_variation     <- 0.5 * mean(absolute_difference)
  derivative <- -0.5 * mean(
    weights * sign(normalized - 1)
  ) / mean_weight^2
  tv_influence <- 0.5 *
    (absolute_difference - mean(absolute_difference)) +
    derivative * (weights - mean_weight)
  tv_mcse <- stats::sd(tv_influence) / sqrt(length(weights))

  list(
    ess_fraction        = ess_fraction,
    ess_fraction_mcse   = ess_mcse,
    total_variation     = total_variation,
    total_variation_mcse = tv_mcse,
    log_weight_iqr      = unname(stats::IQR(log_weights, type = 8))
  )
}


.selection_approximation_summarize <- function(metrics, row_blocks) {

  probabilities <- c(0.05, 0.5, 0.95)
  summarize <- function(values) {
    unname(stats::quantile(values, probabilities, names = FALSE, type = 8))
  }
  B <- length(row_blocks)
  out <- data.frame(
    block       = seq_len(B),
    rows        = vapply(
      row_blocks,
      function(rows) paste(rows, collapse = ","),
      character(1L)
    ),
    n_estimates = vapply(row_blocks, length, integer(1L)),
    stringsAsFactors = FALSE
  )
  for (metric in c("ess_fraction", "total_variation", "log_weight_iqr")) {
    values <- t(vapply(
      seq_len(B),
      function(block) summarize(metrics[[metric]][, block]),
      numeric(3L)
    ))
    out[[paste0(metric, "_q05")]]    <- values[, 1L]
    out[[paste0(metric, "_median")]] <- values[, 2L]
    out[[paste0(metric, "_q95")]]    <- values[, 3L]
  }
  for (metric in c("ess_fraction_mcse", "total_variation_mcse")) {
    out[[paste0(metric, "_median")]] <- apply(
      metrics[[metric]],
      2L,
      stats::median
    )
    out[[paste0(metric, "_max")]] <- apply(
      metrics[[metric]],
      2L,
      max
    )
  }
  out
}
