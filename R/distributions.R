# Native interface helpers ----------------------------------------------------

# Compute sqrt(x^2 + y^2) without intermediate overflow or underflow.
.root_sum_squares <- function(x, y) {

  scale    <- pmax(abs(x), abs(y))
  positive <- !is.na(scale) & is.finite(scale) & scale > 0
  out      <- scale
  x_scaled <- x / scale
  y_scaled <- y / scale
  out[positive] <- scale[positive] * sqrt(
    x_scaled[positive]^2 + y_scaled[positive]^2
  )
  zero <- !is.na(scale) & scale == 0
  out[zero] <- 0
  return(out)
}

.native_numeric_matrix <- function(x) {

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.double(x)) {
    storage.mode(x) <- "double"
  }
  return(x)
}

.native_numeric_vector <- function(x) {

  if (!is.double(x)) {
    x <- as.double(x)
  }
  return(x)
}

.native_integer_vector <- function(x) {

  if (!is.integer(x)) {
    x <- as.integer(x)
  }
  return(x)
}

# Project the complete selected event using its declared Gaussian factor.
# Zero-variance focal coordinates contribute zero here; callers mixing retained
# contexts add their unchanged Gaussian marginals analytically.
.selection_factor_projection <- function(
    mean, diagonal, loading, selection_se, selection_context, execution_plan,
    z, probability, groups, qmc) {

  mean <- .native_numeric_matrix(mean)
  S <- nrow(mean)
  K <- ncol(mean)
  rank <- ncol(loading) %/% K
  selection_context <- BayesTools::selection_context_validate(
    selection_context, n_samples = S,
    required = c("omega", "kernel_mode", "vector_rule")
  )
  cluster <- execution_plan[["quadrature"]]
  factor <- execution_plan[["factor_quadrature"]]
  result <- .Call("RoBMA_selnorm_factor_projection_batch",
    mean, .native_numeric_matrix(diagonal), .native_numeric_matrix(loading),
    as.integer(rank), as.numeric(selection_se), .native_numeric_matrix(selection_context[["omega"]]),
    as.numeric(selection_context[["z_lower"]]), as.numeric(selection_context[["z_upper"]]),
    as.integer(selection_context[["sign"]]), as.integer(selection_context[["kernel_mode"]]),
    isTRUE(selection_context[["telescope_probabilities"]]),
    as.integer(selection_context[["vector_rule"]]), as.integer(match(groups, unique(groups))),
    qmc, as.integer(execution_plan[["points_per_scramble"]]),
    as.integer(execution_plan[["max_points_per_scramble"]]),
    as.integer(execution_plan[["scrambles"]]), as.numeric(execution_plan[["relative_tolerance"]]),
    as.numeric(cluster[["nodes"]]), as.numeric(cluster[["log_weights"]]), as.numeric(cluster[["orders"]]),
    as.numeric(factor[["nodes"]]), as.numeric(factor[["log_weights"]]),
    as.numeric(factor[["orders"]]),
    as.numeric(z), if (probability) 2L else 0L, PACKAGE = "RoBMA")
  if (!is.list(result) || !identical(dim(result[["density"]]), c(S, length(z))) ||
      length(result[["relative_mcse"]]) != S || length(result[["log_density"]]) != S) {
    stop("The selected factor projection returned invalid output.", call. = FALSE)
  }
  error <- max(result[["relative_mcse"]])
  if (!is.finite(error) || error > execution_plan[["relative_tolerance"]]) {
    stop("Zplot selection projection was rejected by diagnostics: relative integration error was ",
      format(error, digits = 4), ". Increase 'max_points_per_scramble' or 'scrambles' in ",
      "'integration_control = set_selection_likelihood_control()'.", call. = FALSE)
  }
  if (any(!is.finite(result[["density"]])) || any(result[["density"]] < 0) ||
      any(!is.finite(result[["log_density"]]))) {
    stop("The selected factor projection returned invalid output.", call. = FALSE)
  }
  result
}


# Integrate the original complete selection event inside Gaussian box limits.
# A caller handling partial observations supplies their conditional Gaussian law
# and the correctly transformed event weights; unspecified rows retain infinite
# limits. Statistical conditioning/source resolution belongs to the caller.
.selection_gaussian_event_mass <- function(
    mean, covariance_lower = NULL, selection_se, selection_context, execution_plan,
    lower = NULL, upper = NULL, qmc = NULL, rank_one_loading = NULL) {

  mean <- .native_numeric_matrix(mean)
  S    <- nrow(mean)
  K    <- ncol(mean)
  selection_context <- BayesTools::selection_context_validate(
    selection_context, n_samples = S,
    required = c("omega", "kernel_mode", "vector_rule")
  )
  limits <- function(x, default, name) {

    if (is.null(x)) return(matrix(default, S, K))
    if (!is.numeric(x) || anyNA(x)) {
      stop("'", name, "' must contain numeric Gaussian event limits.", call. = FALSE)
    }
    if (is.null(dim(x)) && length(x) == K) {
      return(matrix(as.double(x), S, K, byrow = TRUE))
    }
    if (!is.matrix(x) || !identical(dim(x), c(S, K))) {
      stop("'", name, "' must have one value per observation or match 'mean'.",
           call. = FALSE)
    }
    return(.native_numeric_matrix(x))
  }
  lower <- limits(lower, -Inf, "lower")
  upper <- limits(upper, Inf, "upper")
  if (!is.null(rank_one_loading) && !is.null(covariance_lower)) {
    stop("'rank_one_loading' and 'covariance_lower' are alternative Gaussian covariance representations.",
         call. = FALSE)
  }
  points    <- execution_plan[["points_per_scramble"]]
  scrambles <- execution_plan[["scrambles"]]
  if (is.null(points) || is.null(scrambles)) {
    stop("Gaussian selection events require explicit integration settings.", call. = FALSE)
  }
  if (!is.numeric(points) || length(points) != 1L || !is.finite(points) ||
      points < 1 || points > .Machine$integer.max || points != as.integer(points) ||
      !is.numeric(scrambles) || length(scrambles) != 1L || !is.finite(scrambles) ||
      scrambles < 2 || scrambles > .Machine$integer.max || scrambles != as.integer(scrambles)) {
    stop("Gaussian selection event controls are invalid.", call. = FALSE)
  }
  if (is.null(qmc)) {
    qmc <- execution_plan[["designs"]][[as.character(K)]]
    if (is.null(qmc)) {
      qmc <- if (K == 1L || !is.null(rank_one_loading)) numeric() else BayesTools::selection_qmc_design(
        dimensions = 2L * K,
        points     = points,
        scrambles  = scrambles,
        seed       = execution_plan[["seed"]]
      )
    }
  }
  quadrature <- execution_plan[["quadrature"]]
  if (!is.null(quadrature)) attr(quadrature, "bounded_cdf") <- TRUE
  result <- .Call(
    "RoBMA_selnorm_gaussian_event_mass_batch",
    mean, if (is.null(covariance_lower)) NULL else .native_numeric_matrix(covariance_lower),
    .native_numeric_vector(selection_se),
    .native_numeric_matrix(selection_context[["omega"]]),
    .native_numeric_vector(selection_context[["z_lower"]]),
    .native_numeric_vector(selection_context[["z_upper"]]),
    .native_integer_vector(selection_context[["sign"]]),
    .native_integer_vector(selection_context[["kernel_mode"]]),
    .native_integer_vector(selection_context[["vector_rule"]]),
    lower, upper, .native_numeric_vector(qmc),
    as.integer(points), as.integer(scrambles),
    if (is.null(rank_one_loading)) NULL else .native_numeric_matrix(rank_one_loading),
    quadrature,
    .native_numeric_vector(execution_plan[["relative_tolerance"]]),
    PACKAGE = "RoBMA"
  )
  if (!is.null(quadrature)) {
    cdf_error <- attr(result[["relative_quadrature_error"]], "cdf_relative_error", exact = TRUE)
    if (!is.numeric(cdf_error) || length(cdf_error) != S || !is.null(dim(cdf_error)) ||
        any(!is.finite(cdf_error)) || any(cdf_error < 0)) {
      stop("Selection CDF integration diagnostics are unavailable.", call. = FALSE)
    }
  }
  return(result)
}
