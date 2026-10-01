# ============================================================================ #
# Partial Gaussian Selection Events
# ============================================================================ #
#
# Shared full-event weighting and Gaussian conditioning for partial-vector
# densities, deletion scores, and selected-density projections.
# ============================================================================ #
# A partial vector is evaluated inside its original publication event. The
# returned numerator still needs that event's original normalizing constant.
.selection_joint_event_numerator <- function(
    mean, covariance_lower = NULL, selection_se, selection_context, execution_plan,
    observed, values, lower = NULL, upper = NULL, qmc = NULL,
    rank_one_loading = NULL) {

  mean <- as.matrix(mean)
  S <- nrow(mean)
  K <- ncol(mean)
  if (!is.numeric(observed) || anyNA(observed) ||
      any(observed != as.integer(observed)) || anyDuplicated(observed) ||
      any(observed < 1L | observed > K)) {
    stop("Observed selection-event indices must be unique valid row indices.", call. = FALSE)
  }
  observed <- as.integer(observed)
  missing <- setdiff(seq_len(K), observed)
  shared_values <- is.null(dim(values))
  if (is.null(dim(values)) && length(values) != length(observed)) {
    stop("Observed selection-event values must match the observed indices.", call. = FALSE)
  }
  values <- if (is.null(dim(values))) {
    matrix(values, S, length(observed), byrow = TRUE)
  } else {
    as.matrix(values)
  }
  if (!identical(dim(values), c(S, length(observed))) ||
      any(!is.finite(values))) {
    stop("Observed selection-event values must match the Gaussian sample rows.", call. = FALSE)
  }
  pairs <- which(lower.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  if (!is.null(rank_one_loading)) {
    rank_one_loading <- as.matrix(rank_one_loading)
    if (!is.null(covariance_lower) ||
        !identical(dim(rank_one_loading), c(S, K)) ||
        any(!is.finite(rank_one_loading))) {
      stop("A declared rank-one selection event requires finite loadings and no dense covariance.", call. = FALSE)
    }
  } else {
    covariance_lower <- as.matrix(covariance_lower)
    if (!nrow(covariance_lower) %in% c(1L, S) ||
        ncol(covariance_lower) != nrow(pairs)) {
      stop("Selection-event covariance must use the packed lower triangle.", call. = FALSE)
    }
  }
  if (!length(observed)) {
    if (!is.null(covariance_lower) && nrow(covariance_lower) == 1L && S > 1L) {
      covariance_lower <- covariance_lower[rep.int(1L, S), , drop = FALSE]
    }
    result <- .selection_gaussian_event_mass(
      mean, covariance_lower, selection_se, selection_context, execution_plan,
      lower = lower, upper = upper, qmc = qmc,
      rank_one_loading = rank_one_loading
    )
    result[["log_numerator"]] <- result[["log_mass"]]
    return(result)
  }
  conditional_mean <- mean[, missing, drop = FALSE]
  missing_pairs <- which(lower.tri(matrix(0, length(missing), length(missing)),
                                  diag = TRUE), arr.ind = TRUE)
  conditional_lower <- matrix(0, S, nrow(missing_pairs))
  conditional_rank_one <- NULL
  gaussian_log_density <- numeric(S)
  if (!is.null(rank_one_loading)) {
    if (length(observed) != 1L || any(rank_one_loading[, observed] == 0)) {
      stop("Selected partial-vector density is unavailable because the observed rank-one law has no Lebesgue density.", call. = FALSE)
    }
    latent <- (values[, 1L] - mean[, observed]) / rank_one_loading[, observed]
    gaussian_log_density <- stats::dnorm(
      values[, 1L], mean[, observed], abs(rank_one_loading[, observed]), log = TRUE
    )
    conditional_mean <- conditional_mean + rank_one_loading[, missing, drop = FALSE] * latent
    conditional_lower <- NULL
    conditional_rank_one <- matrix(0, S, length(missing))
  } else if (nrow(covariance_lower) == 1L) {
    # A one-row covariance explicitly declares the same law for every mean.
    # Reuse its Gaussian conditioning algebra without inspecting sample values.
    covariance <- matrix(0, K, K)
    covariance[pairs] <- covariance_lower[1L, ]
    covariance[pairs[, 2:1, drop = FALSE]] <- covariance_lower[1L, ]
    factor <- tryCatch(chol(covariance[observed, observed, drop = FALSE]),
                       error = function(e) NULL)
    if (is.null(factor)) {
      stop("Observed selection-event covariance must be positive definite.", call. = FALSE)
    }
    residual <- forwardsolve(t(factor), t(values - mean[, observed, drop = FALSE]))
    gaussian_log_density <- -length(observed) * log(2 * pi) / 2 -
      sum(log(diag(factor))) - colSums(residual^2) / 2
    if (length(missing)) {
      cross <- forwardsolve(t(factor), covariance[observed, missing, drop = FALSE])
      conditional_mean <- conditional_mean + t(crossprod(cross, residual))
      conditional <- covariance[missing, missing, drop = FALSE] - crossprod(cross)
      conditional_lower <- matrix(conditional[missing_pairs], S, nrow(missing_pairs),
                                  byrow = TRUE)
    }
  } else {
    for (s in seq_len(S)) {
      covariance <- matrix(0, K, K)
      covariance[pairs] <- covariance_lower[s, ]
      covariance[pairs[, 2:1, drop = FALSE]] <- covariance_lower[s, ]
      factor <- tryCatch(chol(covariance[observed, observed, drop = FALSE]),
                         error = function(e) NULL)
      if (is.null(factor)) {
        stop("Observed selection-event covariance must be positive definite.", call. = FALSE)
      }
      residual <- forwardsolve(t(factor), values[s, ] - mean[s, observed])
      gaussian_log_density[s] <- -length(observed) * log(2 * pi) / 2 -
        sum(log(diag(factor))) - sum(residual^2) / 2
      if (length(missing)) {
        cross <- forwardsolve(t(factor), covariance[observed, missing, drop = FALSE])
        conditional_mean[s, ] <- conditional_mean[s, ] + as.vector(crossprod(cross, residual))
        conditional <- covariance[missing, missing, drop = FALSE] - crossprod(cross)
        conditional_lower[s, ] <- conditional[missing_pairs]
      }
    }
  }
  if (!length(missing)) {
    return(list(
      log_numerator = gaussian_log_density + .selection_joint_log_weight(
        values, selection_se[observed], selection_context
      ), relative_mcse = numeric(S)
    ))
  }
  conditional_context <- .selection_joint_condition_event_context(
    selection_context,
    if (shared_values) values[1L, , drop = FALSE] else values,
    selection_se[observed]
  )
  result <- .selection_gaussian_event_mass(
    conditional_mean, conditional_lower, selection_se[missing],
    conditional_context, execution_plan, lower = lower, upper = upper, qmc = qmc,
    rank_one_loading = conditional_rank_one
  )
  product <- rep_len(selection_context[["vector_rule"]], S) == 0L
  observed_weight <- numeric(S)
  if (any(product)) {
    observed_weight[product] <- .selection_joint_log_weight(
      values, selection_se[observed], selection_context
    )[product]
  }
  result[["log_numerator"]] <- gaussian_log_density + observed_weight + result[["log_mass"]]
  result
}


.selection_joint_log_weight <- function(y, sei, context) {

  y <- as.matrix(y)
  S <- nrow(y)
  rule <- context[["vector_rule"]]
  if (!length(rule) || !length(rule) %in% c(1L, S) || anyNA(rule) ||
      any(!rule %in% 0:2)) {
    stop("Selection vector-weight metadata are unavailable.", call. = FALSE)
  }
  rule <- rep_len(rule, S)
  mode <- rep_len(context[["kernel_mode"]], S)
  z <- sweep(y * context[["sign"]], 2L, sei, "/")
  out <- numeric(S)
  for (s in which(mode != SELKERNEL_NORMAL)) {
    bins <- if (rule[s] == 0L) {
      .selection_step_bin_from_z(z[s, ], context[["p_cuts"]])
    } else {
      .selection_joint_best_bin(y[s, ], sei, context, rule[s])
    }
    out[s] <- sum(log(context[["omega"]][s, bins]))
  }
  out
}


.selection_joint_best_bin <- function(y, sei, context, rule) {

  # Match the declared physical thresholds at equality, including either tail
  # of a two-sided event. Selection uses the smallest p, regardless of weight.
  threshold <- stats::qnorm(context[["p_cuts"]][-1L], lower.tail = FALSE)
  score <- if (rule == 2L) abs(y) else context[["sign"]] * y
  min(vapply(seq_along(y), function(row) {
    which(score[row] >= sei[row] * threshold)[1L]
  }, integer(1L)))
}


.selection_joint_condition_event_context <- function(context, observed, sei) {

  observed <- as.matrix(observed)
  S <- nrow(context[["omega"]])
  rule <- context[["vector_rule"]]
  if (!length(rule) || !length(rule) %in% c(1L, S) || anyNA(rule) ||
      any(!rule %in% 0:2)) {
    stop("Selection vector-weight metadata are unavailable.", call. = FALSE)
  }
  if (!nrow(observed) %in% c(1L, S) || ncol(observed) != length(sei)) {
    stop("Retained selection-event values must match their sample rows and standard errors.", call. = FALSE)
  }
  if (!ncol(observed)) return(context)
  rule <- rep_len(rule, S)
  cuts <- context[["p_cuts"]]
  bounds <- stats::qnorm(cuts, lower.tail = FALSE)
  midpoint <- vapply(seq_len(length(bounds) - 1L), function(bin) {
    .selection_segment_midpoint(bounds[[bin + 1L]], bounds[[bin]])
  }, numeric(1L))
  for (vector_rule in unique(rule[rule != 0L])) {
    rows <- which(rule == vector_rule)
    scores <- if (vector_rule == 2L) abs(midpoint) else midpoint
    candidate_bins <- .selection_step_bin_from_z(scores, cuts)
    if (nrow(observed) == 1L) {
      retained_bin <- .selection_joint_best_bin(observed[1L, ], sei, context, vector_rule)
      mapped <- pmin(candidate_bins, retained_bin)
      context[["omega"]][rows, ] <- context[["omega"]][rows, mapped, drop = FALSE]
    } else {
      for (s in rows) {
        retained_bin <- .selection_joint_best_bin(observed[s, ], sei, context, vector_rule)
        mapped <- pmin(candidate_bins, retained_bin)
        context[["omega"]][s, ] <- context[["omega"]][s, mapped]
      }
    }
  }
  context
}
