# Sufficient resolution bound for the UNWEIGHTED Gaussian basis.
# Weighted selection still requires the independent model/mass/tail checks.
.zplot_context_gaussian_bound <- function(beta, order) {

  result <- numeric(length(beta))
  if (any(is.na(beta) | beta < 0)) stop("Gaussian resolution ratios must be nonnegative.", call. = FALSE)
  result[is.infinite(beta)] <- Inf
  active <- which(is.finite(beta) & beta > 0)
  if (!length(active)) return(result)
  log_beta <- log(beta[active])
  log_x <- (log(2) + lgamma(2 * order + 1) - lgamma(order + 1)) / order - log(4) - log_beta
  x <- exp(log_x)
  first <- order * log_beta + stats::pgamma(x, shape = order + .5, log.p = TRUE)
  second <- log(4) + stats::pnorm(-sqrt(2 * x), log.p = TRUE)
  largest <- pmax(first, second)
  result[active] <- exp(largest + log(exp(first - largest) + exp(second - largest)))
  result
}


# The emitted diagonal-plus-factor covariance is a numerical integration
# representation. Bound its difference from the original covariance; never
# repair the original matrix or infer a statistical source split from it.
.zplot_context_geometry <- function(covariance, packed, sei, cholesky) {

  geometry <- .Call("RoBMA_selnorm_covariance_envelope_components", packed,
    sei, FALSE, PACKAGE = "RoBMA")
  if (!isTRUE(geometry$available)) return(NULL)
  used <- colSums(geometry$loading != 0) > 0
  geometry$loading <- geometry$loading[, used, drop = FALSE]
  if (!ncol(geometry$loading) || ncol(geometry$loading) > nrow(covariance)) return(NULL)
  radius <- geometry$covariance_residual_radius
  if (!is.matrix(radius) || !identical(dim(radius), dim(covariance)) ||
      any(!is.finite(radius)) || any(radius < 0)) return(NULL)
  difference <- norm(abs(geometry$covariance_residual) + radius, type = "F")
  identity <- diag(nrow(covariance))
  gamma <- function(operations) {

    product <- operations * .Machine$double.eps / 2
    if (product >= 1) Inf else product / (1 - product)
  }
  factor <- t(cholesky)
  inverse <- forwardsolve(factor, identity)
  # These operation-count allowances retain the existing covariance bound.
  # They do not turn library transcendental evaluations into interval proofs.
  inverse_residual <- norm(identity - factor %*% inverse, type = "F") +
    gamma(nrow(covariance) + 2L) *
    norm(abs(factor) %*% abs(inverse) + identity, type = "F")
  if (!is.finite(inverse_residual) || inverse_residual >= 1) return(NULL)
  inverse_bound <- norm(inverse, type = "F") / (1 - inverse_residual)
  reconstruction <- tcrossprod(factor) - covariance
  reconstruction_bound <- norm(reconstruction, type = "F") +
    gamma(nrow(covariance) + 4L) *
    norm(tcrossprod(abs(factor)) + abs(covariance), type = "F")
  relative_reconstruction <- inverse_bound^2 * reconstruction_bound
  if (!is.finite(relative_reconstruction) || relative_reconstruction >= 1) return(NULL)
  epsilon <- inverse_bound^2 * difference / (1 - relative_reconstruction)
  if (!is.finite(epsilon) || epsilon < 0 || epsilon >= 1) return(NULL)
  list(loading = geometry$loading, residual_sd = geometry$residual_sd,
    residual_variance = geometry$residual_variance, groups = geometry$groups,
    type_index = geometry$type_index, relative_error = epsilon)
}


# Called with a validated single-row selection context and metadata-prepared
# retained Gaussian factors. NULL delegates this cell to the existing full-event
# implementation; no alternative target, partial result or false pass is returned.
.zplot_context_projection <- function(z, mean, covariance, context_factor, sei,
                                      selection, control) {

  if (isTRUE(selection$use_normal) || selection$kernel_mode != SELKERNEL_STEP ||
      selection$vector_rule != 0L) return(NULL)
  omega <- as.double(selection$omega)
  if (any(omega <= 0)) return(NULL)
  cdf_enabled <- TRUE
  k <- length(mean)
  stopifnot(k > 0L, length(sei) == k, all(is.finite(mean)),
    is.matrix(covariance), identical(dim(covariance), c(k, k)),
    is.matrix(context_factor), ncol(context_factor) == k,
    nrow(context_factor) >= 1L, all(is.finite(context_factor)),
    length(z) > 0L, all(is.finite(z)), all(is.finite(sei)), all(sei > 0))
  if (!identical(covariance, t(covariance))) return(NULL)
  cholesky <- tryCatch(chol(covariance), error = function(e) NULL)
  if (is.null(cholesky)) return(NULL)
  sei <- as.double(sei)
  if (all(omega == omega[[1L]])) {
    sd <- sqrt(diag(covariance) + colSums(context_factor^2))
    return(list(density = .zplot_normal_density_matrix(z, matrix(mean, 1L), matrix(sd, 1L), sei),
      relative_error = 0, mass_error = 0, mass_omission_error = 0, integration_error = c(absolute = 0)))
  }
  precision <- chol2inv(cholesky)
  packed <- as.double(covariance[lower.tri(covariance, diag = TRUE)])
  point_context <- all(context_factor == 0)
  if (!point_context && nrow(context_factor) != 1L) {
    return(NULL)
  }
  geometry <- .zplot_context_geometry(covariance, packed, sei, cholesky)
  if (is.null(geometry)) return(NULL)
  epsilon_actual <- geometry$relative_error
  epsilon <- epsilon_actual / (1 - epsilon_actual)
  if (!is.finite(epsilon) || epsilon < 0 || epsilon > 1 / 3) {
    return(NULL)
  }
  beta <- epsilon / (2 * (1 - epsilon))
  log_weight_ratio <- log(max(omega)) - log(min(omega))
  tilted_moment <- function(dimension) 4 * dimension * log_weight_ratio + 2 * dimension * log(2)
  determinant_bound <- function(dimension) -.5 * dimension * log1p(-epsilon)
  covariance_log_A_bound <- determinant_bound(k) + beta * tilted_moment(k)
  row_density_bound <- max(omega) / min(omega) * sei / geometry$residual_sd / sqrt(2 * pi)
  compact_density_bound <- mean(row_density_bound)
  # The actual algorithm integrates one common axis and one axis per
  # nonsingleton child. Exactly absorbed singleton children have no axis.
  integration_rank <- if (geometry$groups == 1L) 1L else
    1L + sum(tabulate(geometry$type_index, nbins = geometry$groups) > 1L)
  static <- BayesTools::selection_native_static_args(selection)
  tolerance <- control$relative_tolerance
  stopifnot(is.finite(tolerance), tolerance > 0)
  context_orders <- c(3L, 7L, SELNORM_CLUSTER_QUADRATURE_ORDERS)
  inner_orders <- c(7L, SELNORM_CLUSTER_QUADRATURE_ORDERS)
  grid_range <- range(z)
  # Gaussian basis resolution only. It does not certify the additional
  # selection factors or 1/A(b); all existing weighted-model gates still apply.
  variance <- geometry$residual_variance
  factor_beta <- geometry$loading^2 / (2 * variance)
  context_beta <- colSums(context_factor^2) / (2 * variance)
  height_ratio <- sqrt(1 + (rowSums(geometry$loading^2) + colSums(context_factor^2)) / variance)
  gaussian_resolution <- function(inner_order, outer_order) {

    inner <- height_ratio * rowSums(matrix(.zplot_context_gaussian_bound(factor_beta, inner_order), k))
    outer <- if (point_context) rep(0, k) else height_ratio * .zplot_context_gaussian_bound(context_beta, outer_order)
    c(inner = max(inner), outer = max(outer), total = max(inner + outer))
  }
  final_resolution <- gaussian_resolution(tail(inner_orders, 1L), tail(context_orders, 1L))
  if (!is.finite(final_resolution[["total"]]) || final_resolution[["total"]] > tolerance) {
    return(NULL)
  }
  # Numerator and denominator use the same original ladder but advance
  # independently. Initially N15/N7 share D31; D15 is returned by N15.
  denominator_orders <- inner_orders
  denominator_index <- min(3L, length(denominator_orders))
  context_index <- inner_index <- 2L
  compression_refinements <- 0L
  pruning_refinements <- 0L
  tables <- new.env(parent = emptyenv())
  results <- new.env(parent = emptyenv())
  normalizers <- new.env(parent = emptyenv())
  log_sum <- function(values) {

    largest <- max(values)
    if (is.infinite(largest) && largest < 0) -Inf else largest + log(sum(exp(values - largest)))
  }
  gaussian_tail <- function(nodes) {

    log_sum(c(stats::pnorm(min(nodes), log.p = TRUE),
      stats::pnorm(max(nodes), lower.tail = FALSE, log.p = TRUE)))
  }
  normal_interval <- function(lower, upper) {

    # Independent scalar normal-CDF evaluations, with the original tail choice
    # and arithmetic retained at every endpoint.
    positive <- lower >= 0
    big <- small <- numeric(length(lower))
    big[positive] <- stats::pnorm(lower[positive], lower.tail = FALSE, log.p = TRUE)
    small[positive] <- stats::pnorm(upper[positive], lower.tail = FALSE, log.p = TRUE)
    big[!positive] <- stats::pnorm(upper[!positive], log.p = TRUE)
    small[!positive] <- stats::pnorm(lower[!positive], log.p = TRUE)
    big + log(-expm1(small - big))
  }
  context_table <- function(order) {

    key <- as.character(order)
    if (exists(key, tables, inherits = FALSE)) return(get(key, tables))
    context_rule <- if (point_context) list(nodes = 0, log_weights = 0) else .gauss_hermite_nodes(order)
    means <- if (point_context) matrix(mean, 1L) else
      sweep(outer(context_rule$nodes, as.numeric(context_factor)), 2L, mean, "+")
    log_context_mass <- log_sum(context_rule$log_weights)
    # Bound the normalized compact/full-C density ratio over the complete
    # requested grid/context rectangle. Original C remains the target.
    # The means are affine in one sorted context coordinate. Their actual
    # endpoint rows and the positive-SE grid endpoints give the same extrema
    # as the former repeated full-grid/row scans.
    lower_mean <- pmin(means[1L, ], means[nrow(means), ])
    upper_mean <- pmax(means[1L, ], means[nrow(means), ])
    standardized_distance <- pmax(abs(grid_range[[1L]] * sei - upper_mean),
      abs(grid_range[[2L]] * sei - lower_mean)) /
      sqrt((1 - epsilon_actual) * diag(covariance))
    marginal <- determinant_bound(1L) + beta * standardized_distance^2
    companion <- rep(0, k)
    if (k > 1L) {
      mean_shift <- epsilon / (1 - epsilon) * standardized_distance
      inflated_mean_shift <- mean_shift / (1 - epsilon)
      moment_bound <- tilted_moment(k - 1L)
      if (any(!is.finite(inflated_mean_shift)) || any(beta + inflated_mean_shift / (2 * sqrt(moment_bound)) > 1 / 4)) return(NULL)
      companion <- determinant_bound(k - 1L) + mean_shift^2 / (2 * (1 - epsilon)) +
        beta * moment_bound + inflated_mean_shift * sqrt(moment_bound)
    }
    log_density_ratio <- covariance_log_A_bound + marginal + companion
    covariance_error <- mean(row_density_bound * expm1(log_density_ratio)) * exp(log_context_mass)
    log_density_bound <- log(max(omega)) - log(min(omega)) +
      log(mean(sei * sqrt(diag(precision)))) - .5 * log(2 * pi)
    if (!is.finite(covariance_error)) return(NULL)
    value <- list(means = means, rule = context_rule,
      interval = if (point_context) NULL else normal_interval(
        c(-Inf, context_rule$nodes), c(context_rule$nodes, Inf)),
      mass = exp(log_context_mass),
      covariance_error = covariance_error,
      log_density_bound = log_density_bound,
      outer_tail = if (point_context) 0 else exp(log_density_bound + gaussian_tail(context_rule$nodes)))
    assign(key, value, tables)
    value
  }
  inverse_rule_bound <- function(table, log_A, inner) {

    # These are the raw finite-rule A_q(b), not inferred exact A values.
    # A_q is monotone under the same structural conditions and has a known
    # positive floor including the literal inner quadrature-rule mass.
    log_floor_inverse <- -k * log(min(omega)) - integration_rank * log_sum(inner$log_weights)
    continuous <- if (point_context) -log_A[[1L]] else log_floor_inverse
    if (!point_context) {
      signed_loading <- static$sign * as.numeric(context_factor)
      increasing_weight <- all(diff(omega) <= 0)
      decreasing_weight <- all(diff(omega) >= 0)
      increasing_A <- (increasing_weight && all(signed_loading >= 0)) ||
        (decreasing_weight && all(signed_loading <= 0))
      decreasing_A <- (decreasing_weight && all(signed_loading >= 0)) ||
        (increasing_weight && all(signed_loading <= 0))
      if (increasing_A || decreasing_A) {
        interval <- table$interval
        rectangle <- if (increasing_A) c(log_floor_inverse, -log_A) else
          c(-log_A, log_floor_inverse)
        rectangle_bound <- log_sum(interval + rectangle)
        if (is.finite(rectangle_bound)) continuous <- min(continuous, rectangle_bound)
      }
    }
    max(continuous, log_sum(table$rule$log_weights - log_A))
  }
  log_abs_expm1 <- function(value) {

    answer <- numeric(length(value))
    large <- value > log(2)
    answer[large] <- value[large] + log1p(-exp(-value[large]))
    answer[!large] <- log(abs(expm1(value[!large])))
    answer
  }
  prune_context <- function(table) {

    N <- nrow(table$means)
    keep <- seq_len(N)
    log_omitted <- -Inf
    if (!point_context && pruning_allowance > 0 && N > 1L) {
      log_cap <- log(pruning_allowance) - table$log_density_bound
      candidates <- order(table$rule$log_weights, seq_len(N))
      removed <- integer()
      # Preserve at least one positive rule node. The original weights of all
      # retained nodes are unchanged; neither rule mass nor A is renormalized.
      for (index in candidates[-length(candidates)]) {
        next_mass <- log_sum(c(log_omitted, table$rule$log_weights[[index]]))
        if (next_mass > log_cap) break
        log_omitted <- next_mass
        removed <- c(removed, index)
      }
      keep <- setdiff(keep, removed)
    }
    list(keep = keep, mass = exp(log_omitted),
      retained_mass = exp(log_sum(table$rule$log_weights[keep])),
      error = if (is.infinite(log_omitted) && log_omitted < 0) 0 else
        exp(table$log_density_bound + log_omitted))
  }
  log_add <- function(first, second) {

    largest <- pmax(first, second)
    value <- largest + log1p(exp(-abs(first - second)))
    value[largest == -Inf] <- -Inf
    value
  }
  native_nodes <- function(quadrature) {

    value <- as.double(quadrature$nodes)
    attr(value, "bounded_cdf") <- cdf_enabled
    value
  }
  cdf_record <- function(value, contexts) {

    log_A <- value$compact_log_normalizers
    fields <- c("cdf_factor_log_error", "cdf_log_error", "cdf_numerator_log_error")
    bounds <- lapply(fields, function(name) attr(log_A, name, exact = TRUE))
    valid <- vapply(bounds, function(bound) is.numeric(bound) && length(bound) == contexts &&
      all(is.finite(bound)) && all(bound >= 0), logical(1L))
    if (!all(valid)) stop("Zplot CDF integration diagnostics are unavailable.", call. = FALSE)
    if (!cdf_enabled && any(unlist(bounds, use.names = FALSE) != 0)) {
      stop("Zplot CDF integration diagnostics are inconsistent.", call. = FALSE)
    }
    list(log_A = as.double(log_A),
      mass_error = bounds[[2L]], numerator_error = bounds[[3L]])
  }
  component_omission <- function(compression_error) {

    value <- attr(compression_error, "omitted_mass_error", exact = TRUE)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0) {
      stop("Zplot mixture omission diagnostics are unavailable.", call. = FALSE)
    }
    as.double(value)
  }
  normalizer_key <- function(outer_order, order) paste(outer_order, order, cdf_enabled, sep = "/")
  raw_normalizer <- function(outer_order, order) {

    key <- normalizer_key(outer_order, order)
    if (exists(key, normalizers, inherits = FALSE)) return(get(key, normalizers))
    table <- context_table(outer_order)
    if (is.null(table)) return(NULL)
    quadrature <- .gauss_hermite_nodes(order)
    # The existing native seam calculates every raw mass but collects no
    # density components. Its final positive-mass gate intentionally reports
    # unavailable. Accept only this complete, explicit mass-only result.
    value <- .Call("RoBMA_selnorm_context_star_zplot", table$means,
      packed, sei, omega,
      static$z_lower, static$z_upper, as.integer(static$sign),
      as.logical(static$telescope_probabilities), NULL,
      rep(-Inf, nrow(table$means)), native_nodes(quadrature), as.double(quadrature$log_weights),
      0, 0, PACKAGE = "RoBMA")
    log_A <- value$compact_log_normalizers
    omitted_mass <- component_omission(value$compression_error)
    if (!identical(value$available, FALSE) || !identical(value$mass, 0) ||
        !identical(value$mixture_components, 0) || !identical(as.double(value$compression_error), 0) ||
        omitted_mass != 0 ||
        length(log_A) != nrow(table$means) || any(!is.finite(log_A))) return(NULL)
    record <- cdf_record(value, nrow(table$means))
    assign(key, record, normalizers)
    record
  }
  evaluate <- function(outer_order, inner_order, allowance) {

    denominator_order <- denominator_orders[denominator_index]
    key <- paste(outer_order, inner_order, denominator_order, cdf_enabled, format(allowance, digits = 17),
      format(pruning_allowance, digits = 17), sep = "/")
    if (exists(key, results, inherits = FALSE)) return(get(key, results))
    table <- context_table(outer_order)
    if (is.null(table)) return(NULL)
    denominator <- raw_normalizer(outer_order, denominator_order)
    if (is.null(denominator)) return(NULL)
    log_D <- denominator$log_A
    inner <- .gauss_hermite_nodes(inner_order)
    denominator_rule <- .gauss_hermite_nodes(denominator_order)
    pruned <- prune_context(table)
    keep <- pruned$keep
    context_weights <- table$rule$log_weights
    context_weights[-keep] <- -Inf
    # Both numerator rules receive exactly the same context-specific D.
    # The explicit post-fit policy is carried only on this native node vector.
    out <- .Call("RoBMA_selnorm_context_star_zplot", table$means,
      packed, sei, omega,
      static$z_lower, static$z_upper, as.integer(static$sign),
      as.logical(static$telescope_probabilities), as.double(log_D),
      as.double(context_weights), native_nodes(inner), as.double(inner$log_weights),
      as.double(z), as.double(allowance), PACKAGE = "RoBMA")
    if (!isTRUE(out$available)) return(NULL)
    omitted_mass <- component_omission(out$compression_error)
    if (allowance == 0 && omitted_mass != 0) {
      stop("Zplot mixture omission diagnostics are inconsistent.", call. = FALSE)
    }
    log_A <- out$compact_log_normalizers
    stopifnot(length(log_A) == nrow(table$means), all(is.finite(log_A)),
      all(is.finite(out$density)), all(out$density >= 0),
      is.finite(out$mass), out$mass > 0, is.finite(out$compression_error))
    # Raw Aq is independent of the supplied denominator and compression.
    # Reuse N15's raw A15 as the initial denominator's coarse comparison.
    record <- cdf_record(out, nrow(table$means))
    assign(normalizer_key(outer_order, inner_order), record, normalizers)
    log_tail_A <- gaussian_tail(inner$nodes) + log(integration_rank) + k * log(max(omega))
    # Exact finite D is monotone under the structural gate. Its lower anchors
    # are Dhat*exp(-ED); monotonicity of the polynomial surrogate is not assumed.
    log_inverse <- inverse_rule_bound(table, log_D - denominator$mass_error, denominator_rule)
    numerator_mass_tail <- exp(log_tail_A + log_inverse)
    numerator_tail <- numerator_mass_tail * mean(sei / geometry$residual_sd) / sqrt(2 * pi)
    log_mass_ratio <- log_A[keep] - log_D[keep]
    mass_log_error <- record$mass_error[keep] + denominator$mass_error[keep]
    curve_log_error <- record$numerator_error[keep] + denominator$mass_error[keep]
    log_weighted_ratio <- table$rule$log_weights[keep] + log_mass_ratio
    # Shared-D finite curves are bounded by B0*(Af/D), not by B0 alone.
    # Inflate Ahatf/Dhat to bound Af/D, then apply the numerator/D log error.
    cdf_error <- compact_density_bound * exp(log_sum(log_weighted_ratio + mass_log_error +
      log_abs_expm1(curve_log_error)))
    mass_cdf_error <- exp(log_sum(log_weighted_ratio + log_abs_expm1(mass_log_error)))
    value <- list(density = out$density, mass = out$mass,
      compression_error = as.double(out$compression_error), component_omission_mass = omitted_mass,
      cdf_error = cdf_error, mass_cdf_error = mass_cdf_error,
      denominator_error = denominator$mass_error,
      tail = numerator_tail, mass_tail = numerator_mass_tail,
      log_D = log_D, table = table, keep = keep, omitted_mass = pruned$mass,
      retained_mass = pruned$retained_mass, pruning_error = pruned$error)
    assign(key, value, results)
    value
  }
  inner_table <- function(outer_order, allowance) {

    fine <- evaluate(outer_order, inner_orders[inner_index], allowance)
    coarse <- evaluate(outer_order, inner_orders[inner_index - 1L], allowance)
    if (is.null(fine) || is.null(coarse)) return(NULL)
    stopifnot(identical(fine$keep, coarse$keep), identical(fine$log_D, coarse$log_D),
      identical(fine$denominator_error, coarse$denominator_error))
    # Always compare adjacent DENOMINATOR rules. After D31 becomes D63,
    # use D31/D63 even when the numerator remains at N15.
    previous_D <- raw_normalizer(outer_order, denominator_orders[denominator_index - 1L])
    if (is.null(previous_D)) return(NULL)
    keep <- fine$keep
    log_weights <- fine$table$rule$log_weights[keep]
    log_difference <- log_abs_expm1(previous_D$log_A[keep] - fine$log_D[keep])
    adjacent_error <- previous_D$mass_error[keep] + fine$denominator_error[keep]
    # |1-Dprev/D| <= exp(Eprev+ED)*dhat + expm1(Eprev+ED).
    if (any(adjacent_error > 0)) log_difference <- log_add(
      adjacent_error + log_difference, log_abs_expm1(adjacent_error))
    weighted_difference <- exp(log_sum(log_weights + log_difference))
    denominator_rule <- .gauss_hermite_nodes(denominator_orders[denominator_index])
    log_denominator_tail <- gaussian_tail(denominator_rule$nodes) +
      log(integration_rank) + k * log(max(omega))
    weighted_denominator_tail <- exp(log_sum(log_weights + log_denominator_tail -
      fine$log_D[keep] + fine$denominator_error[keep]))
    weighted_eta <- 2 * weighted_difference + weighted_denominator_tail
    # Common D makes the numerator comparison direct; there is no coarse
    # denominator rescaling penalty. Denominator uncertainty remains separate.
    inner_change <- 2 * max(abs(fine$density - coarse$density))
    normalizer <- compact_density_bound * weighted_eta
    compressed <- 3 * fine$compression_error + 2 * coarse$compression_error
    cdf <- 3 * fine$cdf_error + 2 * coarse$cdf_error
    error <- inner_change + normalizer + compressed + cdf + fine$tail +
      fine$table$covariance_error + fine$pruning_error
    # Both native masses are now meaningful raw ratios An/D. Keep the
    # fine target-mass check AND the common-denominator coarse-rule change.
    # Additional loss from globally omitted INPUT Gaussian components only.
    # Raw An/D is unchanged; this is not a Taylor-polynomial L1 certificate.
    mass_omission_error <- 3 * fine$component_omission_mass + 2 * coarse$component_omission_mass
    mass_error <- abs(fine$mass - fine$retained_mass) +
      2 * abs(fine$mass - coarse$mass) + weighted_eta + fine$mass_tail + fine$omitted_mass +
      3 * fine$mass_cdf_error + 2 * coarse$mass_cdf_error + mass_omission_error
    list(density = fine$density, error = error, mass_error = mass_error / fine$table$mass,
      mass_omission_error = mass_omission_error / fine$table$mass,
      inner_change = inner_change, normalizer = normalizer, compressed = compressed, cdf = cdf,
      covariance = fine$table$covariance_error, tail = fine$tail, pruning = fine$pruning_error,
      outer_tail = fine$table$outer_tail, fine = fine)
  }
  reference_sd <- sqrt(diag(covariance) + colSums(context_factor^2))
  reference_at <- function(index) mean(sei * stats::dnorm(z[[index]] * sei, mean, reference_sd))
  component_height <- sei / reference_sd
  reference_indices <- integer()
  if (all(is.finite(component_height)) && any(component_height > 0)) {
    highest <- which.max(component_height)
    anchor <- mean[[highest]] / sei[[highest]]
    central <- sum(mean / reference_sd) / sum(component_height)
    if (is.finite(anchor) && is.finite(central)) {
      reference_indices <- unique(c(which.min(abs(z - anchor)), which.min(abs(z - central))))
    }
  }
  # These are actual requested-grid scores, hence a lower bound on the old
  # full-grid Gaussian maximum. They only tighten the initial allocations
  # (and the reference contribution to refinement routing); acceptance below
  # still uses the unchanged selected-density error and adjusted-peak gate.
  reference_peak <- if (length(reference_indices))
    max(vapply(reference_indices, reference_at, numeric(1L))) else NA_real_
  if (!is.finite(reference_peak) || reference_peak <= 0) {
    reference_peak <- max(vapply(seq_along(z), reference_at, numeric(1L)))
  }
  allowance <- tolerance * reference_peak / 128
  # Initial allocation only; the unchanged total error gate decides acceptance.
  # Raw A0 anchors and the continuous inverse-Aq tail rectangles remain unpruned.
  pruning_allowance <- if (point_context) 0 else tolerance * reference_peak / 32
  repeat {
    resolution <- gaussian_resolution(inner_orders[inner_index],
      if (point_context) 0L else context_orders[context_index])
    if (!is.finite(resolution[["total"]]) || resolution[["total"]] > tolerance) {
      if (inner_index < length(inner_orders) &&
          (point_context || resolution[["inner"]] >= resolution[["outer"]])) {
        inner_index <- inner_index + 1L
      } else if (!point_context && context_index < length(context_orders)) {
        context_index <- context_index + 1L
      } else if (inner_index < length(inner_orders)) {
        inner_index <- inner_index + 1L
      } else return(NULL)
      next
    }
    fine <- inner_table(if (point_context) 0L else context_orders[context_index], allowance)
    coarse <- if (point_context) NULL else inner_table(context_orders[context_index - 1L], allowance)
    if (is.null(fine) || (!point_context && is.null(coarse))) return(NULL)
    outer_change <- if (point_context) 0 else 2 * max(abs(fine$density - coarse$density))
    combine <- function(field) if (point_context) fine[[field]] else 3 * fine[[field]] + 2 * coarse[[field]]
    error <- outer_change + combine("error") + fine$outer_tail
    # Adjust the finite outer rule's peak for its inner uncertainty. This is
    # a relative-error scale, not a confidence bound on the continuous curve.
    peak_lower <- max(fine$density) - fine$error
    relative_error <- if (peak_lower > 0) error / peak_lower else Inf
    limiting_mass <- if (point_context || fine$mass_error >= coarse$mass_error) fine else coarse
    mass_error <- limiting_mass$mass_error
    mass_omission_error <- limiting_mass$mass_omission_error
    normalizer_part <- combine("normalizer")
    compression_part <- combine("compressed")
    cdf_part <- combine("cdf")
    covariance_part <- combine("covariance")
    pruning_part <- combine("pruning")
    inner_part <- combine("inner_change") + combine("tail")
    if (is.finite(relative_error) && relative_error <= tolerance && mass_error <= tolerance) {
      return(list(density = matrix(fine$density, 1L), relative_error = relative_error,
        mass_error = mass_error, mass_omission_error = mass_omission_error,
        integration_error = c(absolute = error, numerator = inner_part,
          denominator = normalizer_part, compression = compression_part,
          cdf_curve = cdf_part, covariance = covariance_part, pruning = pruning_part,
          context = outer_change + fine$outer_tail),
        quadrature = c(context = if (point_context) 0L else context_orders[context_index],
          numerator = inner_orders[inner_index], denominator = denominator_orders[denominator_index]),
        bounded_cdf = cdf_enabled))
    }
    # Refinement uses a positive observed/reference scale. The acceptance
    # gate above still requires the original uncertainty-adjusted peak.
    # A dominant unresolved inner tail must be refined before tiny ancillary
    # errors trigger pruning restoration or a covariance fallback.
    # More quadrature nodes do not reliably reduce an omission allocation.
    # Tighten the shared compression/input-pruning allowance when the new
    # contribution dominates the limiting mass check, eventually using zero.
    if (mass_error > tolerance && mass_omission_error >= mass_error / 2 && allowance > 0) {
      compression_refinements <- compression_refinements + 1L
      allowance <- if (compression_refinements < 8L) allowance / 2 else 0
      next
    }
    routing_budget <- tolerance * max(max(fine$density), reference_peak)
    priority <- c(inner = inner_part, outer = outer_change + fine$outer_tail,
      normalizer = normalizer_part, compression = compression_part, cdf = cdf_part,
      pruning = pruning_part, covariance = covariance_part)
    advanced <- FALSE
    for (action in names(sort(priority, decreasing = TRUE))) {
      if (action == "normalizer" && denominator_index < length(denominator_orders)) {
        denominator_index <- denominator_index + 1L
      } else if (action == "inner" && inner_index < length(inner_orders)) {
        inner_index <- inner_index + 1L
      } else if (action == "outer" && !point_context && context_index < length(context_orders)) {
        context_index <- context_index + 1L
      } else if (action == "cdf" && cdf_enabled && cdf_part > routing_budget / 4) {
        # Policy is part of both cache keys: no approximate anchors or curves
        # enter the exact rerun. Gaussian/context geometry remains reusable.
        cdf_enabled <- FALSE
      } else if (action == "compression" && allowance > 0 &&
          compression_part > routing_budget / 4) {
        compression_refinements <- compression_refinements + 1L
        allowance <- if (compression_refinements < 8L) allowance / 2 else 0
      } else if (action == "pruning" && pruning_allowance > 0 &&
          (pruning_part > routing_budget / 4 || fine$fine$omitted_mass > tolerance / 4)) {
        pruning_refinements <- pruning_refinements + 1L
        pruning_allowance <- if (pruning_refinements < 8L) pruning_allowance / 2 else 0
      } else if (action == "covariance" && covariance_part > routing_budget) {
        return(NULL)
      } else next
      advanced <- TRUE
      break
    }
    if (!advanced) return(NULL)
  }
}
