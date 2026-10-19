# ============================================================================ #
# IWMDE Weight Kernels
# ============================================================================ #

.iwmde_chen_conditional_gaussian <- function(z_fit, x_fit, x_eval) {

  if (!is.numeric(z_fit) || !is.matrix(x_fit) || !is.matrix(x_eval) ||
      length(z_fit) != nrow(x_fit) || ncol(x_fit) != ncol(x_eval)) {
    stop("Internal IWMDE conditioning inputs are inconsistent.",
         call. = FALSE)
  }
  if (length(z_fit) < 3L) {
    .iwmde_chen_conditional_stop("fewer than three conditioning rows")
  }
  if (any(!is.finite(z_fit)) || any(!is.finite(x_fit))) {
    .iwmde_chen_conditional_stop(
      "conditioning fit rows contain non-finite values"
    )
  }

  column_scale <- apply(abs(x_fit), 2L, max)
  column_scale[column_scale == 0] <- 1
  x_fit  <- sweep(x_fit, 2L, column_scale, "/")
  x_eval <- sweep(x_eval, 2L, column_scale, "/")
  center <- colMeans(x_fit)
  scale  <- apply(x_fit, 2L, stats::sd)
  keep   <- is.finite(scale) & scale > 0
  if (!any(keep)) {
    .iwmde_chen_conditional_stop("all conditioning columns have zero variance")
  }

  x_fit  <- sweep(x_fit[, keep, drop = FALSE], 2L, center[keep], "-")
  x_fit  <- sweep(x_fit, 2L, scale[keep], "/")
  x_eval <- sweep(x_eval[, keep, drop = FALSE], 2L, center[keep], "-")
  x_eval <- sweep(x_eval, 2L, scale[keep], "/")

  if (length(z_fit) <= ncol(x_fit) + 1L) {
    .iwmde_chen_conditional_stop(
      "too few rows for a positive conditional variance"
    )
  }

  focal_scale <- max(abs(z_fit))
  if (focal_scale == 0) {
    .iwmde_chen_conditional_stop("the conditional variance is not positive")
  }
  z_fit     <- z_fit / focal_scale
  covariance <- stats::cov(cbind(z_fit, x_fit))
  if (!all(is.finite(covariance))) {
    .iwmde_chen_conditional_stop(
      "the conditional covariance matrix is not finite"
    )
  }

  covariance_zx <- matrix(covariance[1L, -1L], nrow = 1L)
  covariance_x  <- covariance[-1L, -1L, drop = FALSE]
  factor <- try(chol(covariance_x), silent = TRUE)
  if (inherits(factor, "try-error")) {
    .iwmde_chen_conditional_stop(
      "the conditional covariance matrix is not positive definite"
    )
  }

  coefficients <- backsolve(
    factor,
    forwardsolve(t(factor), as.numeric(covariance_zx))
  )
  residuals <- as.numeric(z_fit - mean(z_fit) - x_fit %*% coefficients)
  variance  <- sum(residuals^2) / (length(residuals) - 1L)
  if (!is.finite(variance) || variance <= 0) {
    .iwmde_chen_conditional_stop("the conditional variance is not positive")
  }

  eval_keep <- stats::complete.cases(x_eval)
  means <- rep(NA_real_, nrow(x_eval))
  means[eval_keep] <- as.numeric(
    mean(z_fit) + x_eval[eval_keep, , drop = FALSE] %*% coefficients
  )

  return(list(
    means       = means,
    sd          = sqrt(variance),
    focal_scale = focal_scale,
    eval_keep   = eval_keep
  ))
}

.iwmde_iwmde_normalization <- function(normalization_grid,
                                       log_q_norm,
                                       baseline_log_q,
                                       log_weight, active_mass,
                                       denominator) {

  empty <- list(
    mass_ratio = NA_real_,
    integral   = NA_real_,
    points     = 0L,
    range      = c(NA_real_, NA_real_),
    scale_type = NA_character_
  )
  if (is.null(normalization_grid) ||
      is.null(log_q_norm) ||
      length(normalization_grid[["x"]]) < 2L ||
      nrow(log_q_norm) != length(normalization_grid[["x"]]) ||
      ncol(log_q_norm) == 0L ||
      length(baseline_log_q) != ncol(log_q_norm) ||
      length(log_weight) != ncol(log_q_norm) ||
      length(log_weight) == 0L) {
    return(empty)
  }

  log_terms <- sweep(
    sweep(log_q_norm, 2L, baseline_log_q, "-"),
    2L,
    log_weight,
    "+"
  )
  log_terms <- sweep(
    log_terms,
    1L,
    normalization_grid[["log_jacobian"]],
    "+"
  )
  density_terms <- .iwmde_density_aggregate(
    log_terms   = log_terms,
    active_mass = active_mass,
    denominator = denominator
  )
  integral <- .iwmde_trapz(
    normalization_grid[["z"]],
    density_terms[["y"]]
  )
  mass_ratio <- if (is.finite(integral) && integral > 0) {
    active_mass / integral
  } else {
    NA_real_
  }

  return(list(
    mass_ratio = mass_ratio,
    integral   = integral,
    points     = length(normalization_grid[["x"]]),
    range      = range(normalization_grid[["x"]]),
    scale_type = "support_grid"
  ))
}


.iwmde_chen_log_weight <- function(context, parameter, parameter_spec,
                                   active_rows, active_values, weight_rows,
                                   weight_values) {

  .iwmde_require_target_type(parameter_spec)
  if (identical(parameter_spec[["type"]], "ordered")) {
    .iwmde_ordered_stop(parameter, "working_distribution_unavailable",
      "the normalized conditional source uses its source normalizer directly")
  }
  active_supports <- .iwmde_parameter_row_supports(
    context        = context,
    parameter      = parameter,
    rows           = active_rows,
    parameter_spec = parameter_spec
  )
  weight_supports <- .iwmde_parameter_row_supports(
    context        = context,
    parameter      = parameter,
    rows           = weight_rows,
    parameter_spec = parameter_spec
  )
  for (population in c("active", "weight")) {
    supports <- if (identical(population, "active")) active_supports else weight_supports
    rows <- if (identical(population, "active")) active_rows else weight_rows
    invalid <- !vapply(seq_len(nrow(supports)), function(i) {
      .iwmde_chen_valid_support(supports[i, ])
    }, logical(1L))
    if (any(invalid)) {
      .iwmde_stop_construction_failure(
        estimator = "iwmde",
        parameter = parameter,
        rows      = rows[invalid],
        stage     = "proposal-support construction",
        detail    = "the declared conditional target support has invalid bounds"
      )
    }
  }
  # Linear and allocated-component targets have row-dependent supports and
  # a proposal for every row, so each row is standardized by its own bounds and
  # the rows are partitioned by which bounds are finite. A fit on the raw values
  # of draw-dependent supports would be a fit on a mixture of supports: no
  # proposal would integrate to one on any row's support. Elsewhere a support is
  # constant within a partition of equal bounds, and is used as it is.
  standardized <- parameter_spec[["type"]] %in% c("linear", "random_component_sd")
  if (standardized) {
    active_chart <- .iwmde_chen_standardize(active_values, active_supports)
    weight_chart <- .iwmde_chen_standardize(weight_values, weight_supports)
    active_keys  <- active_chart[["key"]]
    weight_keys  <- weight_chart[["key"]]
  } else {
    active_keys  <- .iwmde_chen_support_keys(active_supports)
    weight_keys  <- .iwmde_chen_support_keys(weight_supports)
  }
  out          <- rep(-Inf, length(active_values))
  zero_support <- stats::setNames(rep(FALSE, length(active_values)), as.character(active_rows))
  methods      <- character()
  partitions   <- list()

  for (key in unique(active_keys[!is.na(active_keys)])) {
    active_index <- which(active_keys == key)
    weight_index <- which(weight_keys == key)
    if (length(active_index) == 0L || length(weight_index) == 0L) {
      next
    }

    if (standardized) {
      row_support   <- .iwmde_chen_standard_support(key)
      row_active    <- active_chart[["value"]][active_index]
      row_weight    <- weight_chart[["value"]][weight_index]
      log_jacobian  <- active_chart[["log_jacobian"]][active_index]
    } else {
      row_support   <- active_supports[active_index[1L], ]
      row_active    <- active_values[active_index]
      row_weight    <- weight_values[weight_index]
      log_jacobian  <- 0
    }
    if (!.iwmde_chen_valid_support(row_support)) {
      next
    }
    if (!is.finite(row_support[1]) && !is.finite(row_support[2])) {
      weight <- .iwmde_chen_try_weight(
        expr          = .iwmde_chen_conditional_normal_log_weight(
          context        = context,
          parameter      = parameter,
          parameter_spec = parameter_spec,
          active_rows    = active_rows[active_index],
          active_values  = row_active,
          weight_rows    = weight_rows[weight_index],
          weight_values  = row_weight
        ),
        fallback      = .iwmde_chen_marginal_normal_log_weight(
          active_values = row_active,
          weight_values = row_weight
        ),
        fallback_from = "chen_conditional_normal"
      )
    } else if (all(is.finite(row_support))) {
      weight <- .iwmde_chen_try_weight(
        expr          = .iwmde_chen_logit_conditional_normal_log_weight(
          context        = context,
          parameter      = parameter,
          parameter_spec = parameter_spec,
          active_rows    = active_rows[active_index],
          active_values  = row_active,
          weight_rows    = weight_rows[weight_index],
          weight_values  = row_weight,
          support        = row_support
        ),
        fallback      = .iwmde_chen_beta_log_weight(
          active_values = row_active,
          weight_values = row_weight,
          support       = row_support
        ),
        fallback_from = "chen_logit_conditional_normal"
      )
    } else {
      weight <- .iwmde_chen_try_weight(
        expr = .iwmde_chen_compact_conditional_normal_log_weight(
          context, parameter, parameter_spec, active_rows[active_index],
          row_active, weight_rows[weight_index], row_weight, row_support
        ),
        fallback = .iwmde_chen_compact_marginal_normal_log_weight(
          row_active, row_weight, row_support
        ),
        fallback_from = "chen_compact_conditional_normal"
      )
    }
    weight <- .iwmde_chen_weight_fallback_fields(weight)

    out[active_index] <- weight[["log_weight"]] + log_jacobian
    certified <- weight[["zero_support"]]
    if (!is.null(certified)) {
      if (!is.logical(certified) || !is.null(dim(certified)) ||
          length(certified) != length(active_index) || anyNA(certified)) {
        stop("Internal IWMDE zero-support certificate is inconsistent.", call. = FALSE)
      }
      zero_support[active_index] <- certified
    }
    methods <- c(methods, weight[["method"]])
    partitions[[length(partitions) + 1L]] <- list(
      key             = key,
      support         = row_support,
      method          = weight[["method"]],
      fallback_from   = weight[["fallback_from"]],
      fallback_reason = weight[["fallback_reason"]],
      n_eval_rows     = length(active_index),
      n_fit_rows      = length(weight_index),
      n_finite_terms  = sum(is.finite(weight[["log_weight"]]))
    )
  }

  fallback_summary <- .iwmde_chen_fallback_summary(partitions)

  return(list(
    log_weight       = out,
    zero_support     = zero_support,
    method           = .iwmde_chen_method_label(methods, "chen"),
    fallback_from    = fallback_summary[["from"]],
    fallback_reason  = fallback_summary[["reason"]],
    fallback_count   = fallback_summary[["count"]],
    fallback_rows    = fallback_summary[["rows"]],
    fallback_reasons = fallback_summary[["reasons"]],
    partitions       = partitions
  ))
}


.iwmde_chen_try_weight <- function(expr, fallback, fallback_from) {

  tryCatch(
    .iwmde_chen_weight_fallback_fields(expr),
    iwmde_chen_weight_unavailable = function(e) {

      .iwmde_chen_weight_fallback_fields(
        weight          = fallback,
        fallback_from   = fallback_from,
        fallback_reason = gsub("[\r\n\t]+", " ", conditionMessage(e))
      )
    }
  )
}


.iwmde_chen_weight_fallback_fields <- function(
    weight, fallback_from = NULL, fallback_reason = NULL) {

  if (!is.null(fallback_from) || is.null(weight[["fallback_from"]])) {
    weight[["fallback_from"]] <- if (is.null(fallback_from)) {
      NA_character_
    } else {
      as.character(fallback_from)[[1L]]
    }
  }
  if (!is.null(fallback_reason) || is.null(weight[["fallback_reason"]])) {
    weight[["fallback_reason"]] <- if (is.null(fallback_reason)) {
      NA_character_
    } else {
      as.character(fallback_reason)[[1L]]
    }
  }

  return(weight)
}


.iwmde_chen_fallback_summary <- function(partitions) {

  if (length(partitions) == 0L) {
    return(list(
      count   = 0L,
      rows    = 0L,
      from    = character(),
      reason  = character(),
      reasons = structure(integer(), names = character())
    ))
  }

  fallback_reason <- vapply(partitions, function(partition) {
    partition[["fallback_reason"]]
  }, character(1))
  fallback_from <- vapply(partitions, function(partition) {
    partition[["fallback_from"]]
  }, character(1))
  fallback_index <- !is.na(fallback_reason)
  reason_table   <- table(fallback_reason[fallback_index])
  reason_counts  <- as.integer(reason_table)
  names(reason_counts) <- names(reason_table)

  return(list(
    count   = as.integer(sum(fallback_index)),
    rows    = as.integer(sum(vapply(
      partitions[fallback_index],
      `[[`,
      integer(1),
      "n_eval_rows"
    ))),
    from    = unique(fallback_from[fallback_index]),
    reason  = unique(fallback_reason[fallback_index]),
    reasons = reason_counts
  ))
}


.iwmde_chen_valid_support <- function(support) {

  return(
    length(support) == 2L &&
      !any(is.na(support)) &&
      support[1] < support[2]
  )
}


# Standardize values by the bounds of their own row: the unit interval when
# both are finite (log Jacobian -log(upper - lower)), the distance from the
# finite bound when one is (Jacobian 1), and the value itself when none is. A
# row with missing or empty bounds has no key and is left without a weight. The
# key is the pattern of finite bounds.
.iwmde_chen_standardize <- function(values, supports) {

  lower <- supports[, 1L]
  upper <- supports[, 2L]
  valid <- !is.na(lower) & !is.na(upper) & lower < upper
  both  <- valid & is.finite(lower) & is.finite(upper)
  low   <- valid & is.finite(lower) & !is.finite(upper)
  high  <- valid & !is.finite(lower) & is.finite(upper)
  none  <- valid & !is.finite(lower) & !is.finite(upper)
  key   <- rep(NA_character_, length(values))
  key[both] <- "both"
  key[low]  <- "lower"
  key[high] <- "upper"
  key[none] <- "none"
  value <- rep(NA_real_, length(values))
  value[both] <- (values[both] - lower[both]) / (upper[both] - lower[both])
  value[low]  <- values[low] - lower[low]
  value[high] <- upper[high] - values[high]
  value[none] <- values[none]
  log_jacobian <- rep(0, length(values))
  log_jacobian[both] <- -log(upper[both] - lower[both])

  list(key = unname(key), value = unname(value),
       log_jacobian = unname(log_jacobian))
}


# The support of the standardized values of a partition.
.iwmde_chen_standard_support <- function(key) {

  switch(key,
    both  = c(0, 1),
    lower = c(0, Inf),
    upper = c(0, Inf),
    none  = c(-Inf, Inf)
  )
}


.iwmde_chen_support_keys <- function(supports) {

  # One vectorized encoding per bound instead of one closure call per row.
  paste(
    .iwmde_key_number(supports[, 1L]),
    .iwmde_key_number(supports[, 2L]),
    sep = ","
  )
}


.iwmde_chen_method_label <- function(methods, fallback) {

  methods <- unique(methods[nzchar(methods)])
  if (length(methods) == 0L) {
    return(fallback)
  }
  if (length(methods) == 1L) {
    return(methods)
  }

  return(paste0("chen_mixed(", paste(methods, collapse = ","), ")"))
}


.iwmde_chen_support_distance <- function(values, support) {

  if (is.finite(support[[1L]])) values - support[[1L]] else support[[2L]] - values
}


# Endpoints and membership stay in r = distance / focal_scale coordinates.
# The actual rounded endpoints determine the normalizing mass, not nominal .95.
.iwmde_chen_compact_normal_log_weight <- function(distance, mean, sd,
                                                 focal_scale, method) {

  n <- length(distance)
  mean <- rep_len(mean, n)
  sd <- rep_len(sd, n)
  out <- rep(-Inf, n)
  certified <- rep(FALSE, n)
  if (!is.numeric(focal_scale) || length(focal_scale) != 1L ||
      !is.finite(focal_scale) || focal_scale <= 0) {
    return(list(log_weight = out, zero_support = certified, method = method))
  }
  valid <- is.finite(mean) & is.finite(sd) & sd > 0
  log_s0 <- rep(NA_real_, n)
  log_s0[valid] <- stats::pnorm(0, mean[valid], sd[valid], lower.tail = FALSE, log.p = TRUE)
  valid <- valid & is.finite(log_s0)
  lower <- upper <- rep(NA_real_, n)
  lower[valid] <- stats::qnorm(log_s0[valid] + log(.975), mean[valid], sd[valid], lower.tail = FALSE, log.p = TRUE)
  upper[valid] <- stats::qnorm(log_s0[valid] + log(.025), mean[valid], sd[valid], lower.tail = FALSE, log.p = TRUE)
  valid <- valid & is.finite(lower) & is.finite(upper) & lower > 0 & lower < upper
  log_mass <- rep(NA_real_, n)
  log_mass[valid] <- .iwmde_qcmde_log_normal_mass(lower[valid], upper[valid], mean[valid], sd[valid])
  valid <- valid & is.finite(log_mass) & log_mass <= 0
  r <- distance / focal_scale
  evaluated <- valid & is.finite(distance) & distance >= 0 & !is.na(r)
  certified[evaluated] <- r[evaluated] < lower[evaluated] | r[evaluated] > upper[evaluated]
  inside <- evaluated & !certified
  out[inside] <- stats::dnorm(r[inside], mean[inside], sd[inside], log = TRUE) -
    log_mass[inside] - log(focal_scale)
  list(log_weight = out, zero_support = certified, method = method)
}


.iwmde_chen_compact_marginal_normal_log_weight <- function(active_values,
                                                          weight_values,
                                                          support) {

  distance_fit <- .iwmde_chen_support_distance(weight_values, support)
  distance_eval <- .iwmde_chen_support_distance(active_values, support)
  method <- "chen_compact_marginal_normal"
  if (length(distance_fit) < 3L || any(!is.finite(distance_fit)) || any(distance_fit <= 0)) {
    return(list(log_weight = rep(-Inf, length(active_values)),
      zero_support = rep(FALSE, length(active_values)), method = method))
  }
  focal_scale <- max(abs(distance_fit))
  scaled <- distance_fit / focal_scale
  .iwmde_chen_compact_normal_log_weight(distance_eval, mean(scaled), stats::sd(scaled), focal_scale, method)
}


.iwmde_chen_compact_conditional_normal_log_weight <- function(context, parameter,
    parameter_spec, active_rows, active_values, weight_rows, weight_values, support) {

  conditioning <- .iwmde_chen_conditioning_matrix(context, parameter, parameter_spec,
    active_rows, weight_rows)
  if (!ncol(conditioning[["fit"]])) {
    return(.iwmde_chen_compact_marginal_normal_log_weight(active_values, weight_values, support))
  }
  distance_fit <- .iwmde_chen_support_distance(weight_values, support)
  if (any(!is.finite(distance_fit)) || any(distance_fit <= 0)) {
    .iwmde_chen_conditional_stop("conditioning fit rows lie outside the declared positive distance support")
  }
  gaussian <- .iwmde_chen_conditional_gaussian(distance_fit, conditioning[["fit"]], conditioning[["eval"]])
  .iwmde_chen_compact_normal_log_weight(.iwmde_chen_support_distance(active_values, support),
    gaussian[["means"]], gaussian[["sd"]], gaussian[["focal_scale"]], "chen_compact_conditional_normal")
}


.iwmde_chen_logit_conditional_normal_log_weight <- function(context,
                                                            parameter,
                                                            parameter_spec,
                                                            active_rows,
                                                            active_values,
                                                            weight_rows,
                                                            weight_values,
                                                            support) {

  active <- .iwmde_chen_logit_transform(active_values, support)
  weight <- .iwmde_chen_logit_transform(weight_values, support)
  out    <- rep(-Inf, length(active_values))

  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    active_rows    = active_rows,
    weight_rows    = weight_rows
  )
  if (ncol(conditioning[["fit"]]) == 0L) {
    return(.iwmde_chen_beta_log_weight(
      active_values = active_values,
      weight_values = weight_values,
      support       = support
    ))
  }

  x_fit  <- conditioning[["fit"]]
  x_eval <- conditioning[["eval"]]
  z_fit  <- weight[["z"]]
  gaussian <- .iwmde_chen_conditional_gaussian(
    z_fit  = z_fit,
    x_fit  = x_fit,
    x_eval = x_eval
  )
  eval_keep <- gaussian[["eval_keep"]] & is.finite(active[["z"]]) &
    is.finite(active[["log_jacobian"]]) &
    is.finite(gaussian[["means"]])
  out[eval_keep] <- stats::dnorm(
    active[["z"]][eval_keep] / gaussian[["focal_scale"]],
    mean = gaussian[["means"]][eval_keep],
    sd   = gaussian[["sd"]],
    log  = TRUE
  ) - log(gaussian[["focal_scale"]]) +
    active[["log_jacobian"]][eval_keep]

  return(list(log_weight = out, method = "chen_logit_conditional_normal"))
}


.iwmde_chen_logit_transform <- function(values, support) {

  lower <- support[1]
  upper <- support[2]
  width <- upper - lower
  prob  <- (values - lower) / width

  return(list(
    z            = stats::qlogis(prob),
    log_jacobian = -log(width) - log(prob) - log1p(-prob)
  ))
}


.iwmde_chen_beta_log_weight <- function(active_values, weight_values,
                                        support) {

  out   <- rep(-Inf, length(active_values))
  lower <- support[1]
  upper <- support[2]
  width <- upper - lower

  prob_fit <- (weight_values - lower) / width
  prob_fit <- prob_fit[
    is.finite(prob_fit) & prob_fit > 0 & prob_fit < 1
  ]
  prob_eval <- (active_values - lower) / width
  keep      <- is.finite(prob_eval) & prob_eval > 0 & prob_eval < 1
  if (length(prob_fit) < 3L || !any(keep)) {
    return(list(log_weight = out, method = "chen_beta"))
  }

  location <- mean(prob_fit)
  variance <- stats::var(prob_fit)
  maximum  <- location * (1 - location)
  if (!is.finite(location) || !is.finite(variance) ||
      location <= 0 || location >= 1 ||
      variance <= 0 ||
      variance >= maximum) {
    return(list(log_weight = out, method = "chen_beta"))
  }

  common <- maximum / variance - 1
  alpha  <- location * common
  beta   <- (1 - location) * common
  if (!is.finite(alpha) || !is.finite(beta) ||
      alpha <= 0 || beta <= 0) {
    return(list(log_weight = out, method = "chen_beta"))
  }

  out[keep] <- stats::dbeta(
    prob_eval[keep],
    shape1 = alpha,
    shape2 = beta,
    log    = TRUE
  ) - log(width)

  return(list(log_weight = out, method = "chen_beta"))
}


.iwmde_chen_conditional_stop <- function(reason) {

  condition <- structure(
    list(
      message = paste0(
        "IWMDE Chen conditional-normal weights are unavailable: ",
        reason,
        "."
      ),
      call = NULL
    ),
    class = c(
      "iwmde_chen_weight_unavailable",
      "error",
      "condition"
    )
  )
  stop(condition)
}

.iwmde_chen_conditional_normal_log_weight <- function(context, parameter,
                                                      parameter_spec,
                                                      active_rows,
                                                      active_values,
                                                      weight_rows,
                                                      weight_values) {

  conditioning <- .iwmde_chen_conditioning_matrix(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec,
    active_rows    = active_rows,
    weight_rows    = weight_rows
  )
  if (ncol(conditioning[["fit"]]) == 0L) {
    return(.iwmde_chen_marginal_normal_log_weight(
      active_values = active_values,
      weight_values = weight_values
    ))
  }

  x_fit  <- conditioning[["fit"]]
  x_eval <- conditioning[["eval"]]
  z_fit  <- weight_values
  gaussian <- .iwmde_chen_conditional_gaussian(
    z_fit  = z_fit,
    x_fit  = x_fit,
    x_eval = x_eval
  )
  out       <- rep(-Inf, length(active_values))
  eval_keep <- gaussian[["eval_keep"]] & is.finite(active_values) &
    is.finite(gaussian[["means"]])
  out[eval_keep] <- stats::dnorm(
    active_values[eval_keep] / gaussian[["focal_scale"]],
    mean = gaussian[["means"]][eval_keep],
    sd   = gaussian[["sd"]],
    log  = TRUE
  ) - log(gaussian[["focal_scale"]])

  return(list(log_weight = out, method = "chen_conditional_normal"))
}


.iwmde_chen_marginal_normal_log_weight <- function(active_values,
                                                   weight_values) {

  out    <- rep(-Inf, length(active_values))
  values <- as.numeric(weight_values)
  if (length(values) < 2L || any(!is.finite(values))) {
    return(list(log_weight = out, method = "chen_marginal_normal"))
  }

  focal_scale <- max(abs(values))
  if (focal_scale == 0) {
    return(list(log_weight = out, method = "chen_marginal_normal"))
  }
  scaled_values <- values / focal_scale
  location      <- mean(scaled_values)
  scale         <- stats::sd(scaled_values)
  finite        <- is.finite(active_values)
  if (!is.finite(location) || !is.finite(scale) ||
      scale <= 0) {
    return(list(log_weight = out, method = "chen_marginal_normal"))
  }

  out[finite] <- stats::dnorm(
    active_values[finite] / focal_scale,
    mean = location,
    sd   = scale,
    log  = TRUE
  ) - log(focal_scale)

  return(list(log_weight = out, method = "chen_marginal_normal"))
}


.iwmde_chen_conditioning_matrix <- function(context, parameter,
                                            parameter_spec, active_rows,
                                            weight_rows) {

  context <- .iwmde_context_ensure_caches(context)
  samples <- context[["posterior_samples"]]
  columns <- .iwmde_chen_conditioning_columns(
    context        = context,
    parameter      = parameter,
    parameter_spec = parameter_spec
  )
  if (length(columns) == 0L) {
    return(list(
      fit     = matrix(numeric(), nrow = length(weight_rows), ncol = 0L),
      eval    = matrix(numeric(), nrow = length(active_rows), ncol = 0L),
      columns = character()
    ))
  }

  constrained <- vapply(
    columns,
    .iwmde_chen_conditioning_column_is_constrained,
    logical(1),
    context = context
  )
  columns <- columns[order(!constrained)]

  fit_values  <- matrix(NA_real_, nrow = length(weight_rows), ncol = length(columns))
  eval_values <- matrix(NA_real_, nrow = length(active_rows), ncol = length(columns))
  colnames(fit_values)  <- columns
  colnames(eval_values) <- columns
  keep_columns          <- rep(TRUE, length(columns))

  for (i in seq_along(columns)) {
    column     <- columns[[i]]
    raw_fit    <- .iwmde_parameter_column_values(
      context   = context,
      samples   = samples[weight_rows, , drop = FALSE],
      parameter = column
    )
    raw_eval   <- .iwmde_parameter_column_values(
      context   = context,
      samples   = samples[active_rows, , drop = FALSE],
      parameter = column
    )
    raw_values <- c(raw_fit, raw_eval)
    if (all(is.finite(raw_values)) && length(unique(raw_values)) == 1L) {
      keep_columns[[i]] <- FALSE
      next
    }

    transformed <- .iwmde_chen_transform_conditioning_column(
      context     = context,
      fit_values  = raw_fit,
      eval_values = raw_eval,
      column      = column
    )
    if (any(!is.finite(transformed[["fit"]]))) {
      .iwmde_chen_conditional_stop(
        "conditioning fit columns contain non-finite transformed values"
      )
    }
    fit_values[, column]  <- transformed[["fit"]]
    eval_values[, column] <- transformed[["eval"]]
  }

  columns     <- columns[keep_columns]
  fit_values  <- fit_values[, keep_columns, drop = FALSE]
  eval_values <- eval_values[, keep_columns, drop = FALSE]

  if (length(columns) == 0L) {
    return(list(
      fit     = matrix(numeric(), nrow = length(weight_rows), ncol = 0L),
      eval    = matrix(numeric(), nrow = length(active_rows), ncol = 0L),
      columns = character()
    ))
  }

  keep <- vapply(seq_along(columns), function(i) {
    values <- fit_values[, i]
    length(unique(values)) > 1L
  }, logical(1))

  if (!any(keep)) {
    return(list(
      fit     = matrix(numeric(), nrow = length(weight_rows), ncol = 0L),
      eval    = matrix(numeric(), nrow = length(active_rows), ncol = 0L),
      columns = character()
    ))
  }

  columns     <- columns[keep]
  fit_values  <- fit_values[, keep, drop = FALSE]
  eval_values <- eval_values[, keep, drop = FALSE]

  return(list(
    fit     = fit_values,
    eval    = eval_values,
    columns = columns
  ))
}


.iwmde_chen_conditioning_columns <- function(context, parameter,
                                             parameter_spec) {

  samples <- context[["posterior_samples"]]
  columns <- colnames(samples)
  if (!is.null(parameter_spec[["target_columns"]])) {
    target_columns <- parameter_spec[["target_columns"]]
  } else if (!is.null(parameter_spec) &&
             identical(parameter_spec[["type"]], "linear")) {
    target_columns <- names(parameter_spec[["weights"]])
  } else {
    target_columns <- parameter
  }
  conditioning_exclude <- parameter_spec[["conditioning_exclude"]]
  if (is.null(conditioning_exclude)) {
    conditioning_exclude <- character()
  }

  candidates <- columns[
    !columns %in% target_columns &
      !columns %in% conditioning_exclude &
      !columns %in% context[["indicator_names"]] &
      !vapply(columns, .iwmde_parameter_is_local_latent, logical(1)) &
      vapply(columns, function(column) {
        .iwmde_chen_is_global_conditioning_column(
          context = context,
          column  = column
        )
      }, logical(1))
  ]

  return(candidates)
}


.iwmde_chen_is_global_conditioning_column <- function(context, column) {

  column %in% .iwmde_chen_global_conditioning_columns(context)
}

.iwmde_chen_global_conditioning_columns <- function(context) {

  context <- .iwmde_context_ensure_caches(context)
  cache   <- context[["prior_cache"]]
  key     <- "chen_global_conditioning_columns"
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }

  sample_columns <- colnames(context[["posterior_samples"]])
  prior_list     <- context[["flat_prior_list"]]
  columns <- sample_columns[vapply(sample_columns, function(column) {
    prior_name <- .iwmde_parameter_prior_name(context, column)
    if (is.null(prior_name)) {
      return(FALSE)
    }

    prior <- prior_list[[prior_name]]
    !BayesTools::is.prior(prior) ||
      (!BayesTools::is.prior.point(prior) &&
       !BayesTools::is.prior.none(prior))
  }, logical(1))]

  simplex_names <- names(prior_list)[vapply(
    prior_list,
    inherits,
    logical(1),
    what = "prior.simplex"
  )]
  for (parameter in simplex_names) {
    prior <- prior_list[[parameter]]
    alpha <- prior[["parameters"]][["alpha"]]
    if (!is.numeric(alpha) || length(alpha) < 2L) {
      stop(
        "Simplex conditioning metadata are invalid for '", parameter, "'.",
        call. = FALSE
      )
    }
    expected <- paste0(parameter, "[", seq_along(alpha), "]")
    if (!all(expected %in% sample_columns)) {
      .stop_refit_required(
        "Simplex conditioning coordinates are missing for '", parameter,
        "'. Refit the model with the current RoBMA/BayesTools build."
      )
    }

    columns <- setdiff(columns, expected[[length(expected)]])
  }

  selection_columns <- sample_columns[vapply(
    sample_columns,
    .iwmde_parameter_is_weightfunction_coordinate,
    logical(1),
    context         = context,
    include_private = FALSE
  )]
  columns <- unique(c(columns, selection_columns))
  columns <- columns[
    !columns %in% context[["indicator_names"]] &
      !vapply(columns, .iwmde_parameter_is_local_latent, logical(1))
  ]
  assign(key, columns, envir = cache)

  return(columns)
}


.iwmde_chen_transform_conditioning_column <- function(context,
                                                      fit_values,
                                                      eval_values,
                                                      column) {

  omega_coordinates <- if (!is.null(context[["selection_spec"]])) {
    context[["selection_spec"]][["jags_omega"]]
  } else {
    character()
  }
  if (.iwmde_parameter_matches_coordinate(column, omega_coordinates)) {
    support <- .iwmde_chen_omega_support(context)
    if (all(is.finite(support))) {
      # Declared finite omega supports include structural endpoint states in
      # product-space rows. Keep them as finite regressors; the conditional
      # Gaussian fit already centers/scales columns and removes constants.
      fit  <- as.numeric(fit_values)
      eval <- as.numeric(eval_values)
      fit[!is.finite(fit) | fit < support[[1L]] | fit > support[[2L]]] <- NA_real_
      eval[!is.finite(eval) | eval < support[[1L]] | eval > support[[2L]]] <- NA_real_
      return(list(fit = fit, eval = eval))
    }
    return(.iwmde_chen_transform_lower_bounded(
      fit_values,
      eval_values,
      support[[1L]]
    ))
  }

  support <- .iwmde_chen_conditioning_support(context, column)
  if (!is.null(support)) {
    if (all(is.finite(support))) {
      return(.iwmde_chen_transform_bounded(
        fit_values,
        eval_values,
        support
      ))
    }
    if (is.finite(support[[1L]])) {
      return(.iwmde_chen_transform_lower_bounded(
        fit_values,
        eval_values,
        support[[1L]]
      ))
    }
    if (is.finite(support[[2L]])) {
      return(.iwmde_chen_transform_upper_bounded(
        fit_values,
        eval_values,
        support[[2L]]
      ))
    }
  }

  return(list(
    fit  = as.numeric(fit_values),
    eval = as.numeric(eval_values)
  ))
}


.iwmde_chen_conditioning_support <- function(context, column) {

  prior_name <- .iwmde_parameter_prior_name(context, column)
  if (is.null(prior_name)) {
    return(NULL)
  }
  prior <- context[["flat_prior_list"]][[prior_name]]
  if (identical(prior_name, "bias") &&
      column %in% c("PET", "PEESE")) {
    predicate <- if (identical(column, "PET")) {
      BayesTools::is.prior.PET
    } else {
      BayesTools::is.prior.PEESE
    }
    branches <- if (BayesTools::is.prior.mixture(prior)) prior else list(prior)
    branches <- branches[vapply(branches, predicate, logical(1))]
    if (length(branches) == 0L) {
      return(NULL)
    }
    supports <- vapply(branches, .iwmde_prior_support, numeric(2))
    return(c(min(supports[1L, ]), max(supports[2L, ])))
  }
  if (inherits(prior, "prior.simplex")) {
    return(c(0, 1))
  }
  if (BayesTools::is.prior.mixture(prior)) {
    supports <- vapply(prior, .iwmde_prior_support, numeric(2))
    same_support <- ncol(supports) > 0L && all(vapply(
      seq_len(ncol(supports)),
      function(i) identical(supports[, i], supports[, 1L]),
      logical(1)
    ))
    if (!same_support) {
      return(NULL)
    }
    return(supports[, 1L])
  }
  if (!BayesTools::is.prior(prior)) {
    return(NULL)
  }

  .iwmde_prior_support(prior)
}


.iwmde_chen_conditioning_column_is_constrained <- function(column, context) {

  omega_coordinates <- if (!is.null(context[["selection_spec"]])) {
    context[["selection_spec"]][["jags_omega"]]
  } else {
    character()
  }
  if (.iwmde_parameter_matches_coordinate(column, omega_coordinates)) {
    return(TRUE)
  }

  support <- .iwmde_chen_conditioning_support(context, column)
  !is.null(support) && any(is.finite(support))
}


.iwmde_chen_omega_support <- function(context) {

  prior <- context[["flat_prior_list"]][["bias"]]
  if (is.null(prior)) {
    return(c(0, Inf))
  }

  branches <- if (BayesTools::is.prior.mixture(prior)) prior else list(prior)
  supports <- lapply(branches, function(branch) {
    if (BayesTools::is.prior.weightfunction(branch)) {
      return(.iwmde_prior_support(branch))
    }
    c(1, 1)
  })
  support_matrix <- do.call(cbind, supports)
  c(min(support_matrix[1L, ]), max(support_matrix[2L, ]))
}


.iwmde_chen_transform_lower_bounded <- function(fit_values, eval_values,
                                                lower) {

  fit  <- as.numeric(fit_values)
  eval <- as.numeric(eval_values)
  fit[fit < lower]   <- NA_real_
  eval[eval < lower] <- NA_real_

  return(list(
    fit  = log1p(fit - lower),
    eval = log1p(eval - lower)
  ))
}


.iwmde_chen_transform_upper_bounded <- function(fit_values, eval_values,
                                                upper) {

  fit  <- as.numeric(fit_values)
  eval <- as.numeric(eval_values)
  fit[fit > upper]   <- NA_real_
  eval[eval > upper] <- NA_real_

  return(list(
    fit  = -log1p(upper - fit),
    eval = -log1p(upper - eval)
  ))
}


.iwmde_chen_transform_bounded <- function(fit_values, eval_values, support) {

  lower <- support[[1L]]
  upper <- support[[2L]]
  width <- upper - lower
  fit   <- (as.numeric(fit_values) - lower) / width
  eval  <- (as.numeric(eval_values) - lower) / width
  fit[fit < 0 | fit > 1]     <- NA_real_
  eval[eval < 0 | eval > 1] <- NA_real_

  return(list(
    fit  = stats::qlogis(fit),
    eval = stats::qlogis(eval)
  ))
}
