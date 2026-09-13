# Call-owned approximation of fixed scalar normalizer grids. Native anchors
# are actual requested points; the original priors and Gaussian law stay intact.
.selection_normalizer_grid <- function(context, queries, rows) {

  if (!.is_data_joint_selection(context[["data"]]) ||
      !.is_data_known_v(context[["data"]]) ||
      .selection_retains_sampling(context[["data"]])) return(NULL)
  plan <- .data_selection_execution_plan(context[["data"]])
  dense <- sum(plan[["block_methods"]] == "dense")
  queries <- sort(unique(queries[is.finite(queries)]))
  if (!dense || length(queries) < 16L || !length(rows)) return(NULL)
  # Reserve for native records, reported error, binding data and R overhead.
  bytes <- 4 * (length(queries) * length(rows) * dense * 48 +
    length(rows) * sum(lengths(plan[["row_blocks"]])^2) * 8)
  if (!is.finite(bytes) || bytes > .known_v_covariance_max_bytes()) return(NULL)
  out <- new.env(parent = emptyenv())
  out$queries <- queries
  out$rows <- rows
  out$cases <- new.env(parent = emptyenv())
  out$geometries <- new.env(parent = emptyenv())
  out$stats <- c(anchor_evaluations = 0, fallback_evaluations = 0,
    constant_evaluations = 0, gaussian_evaluations = 0, interpolated_points = 0,
    max_relative_error = 0, max_anchor_relative_error = 0,
    max_interpolation_relative_error = 0, max_native_mcse = 0)
  out$untracked <- FALSE
  out$geometry_limit <- .known_v_covariance_max_bytes() - bytes
  out$geometry_bytes <- 0
  out
}

.selection_normalizer_grid_geometry <- function(shared, ids) {

  ids <- sort(unique(ids))
  key <- paste(ids, collapse = ":")
  if (exists(key, shared$geometries, inherits = FALSE)) return(get(key, shared$geometries))
  # A geometry is retained only through this bounded shared cache. Cases
  # reference its immutable matrices; uncached geometries are not retained.
  estimate <- 4 * (8 * length(shared$queries) * (3 * 32 + 6) + 4096)
  if (estimate > shared$geometry_limit - shared$geometry_bytes) return(list())
  x <- shared$queries[ids]
  midpoint <- (min(x) + max(x)) / 2
  leaves <- list()
  unit <- .Machine$double.eps
  gamma <- function(n) (n * unit) / (1 - n * unit)
  normal <- function(value) is.finite(value) & abs(value) >= .Machine$double.xmin
  for (side in c(FALSE, TRUE)) {
    available <- ids[(x > midpoint) == side]
    n <- min(32L, length(available) %/% 2L)
    if (n < 8L) next
    values <- shared$queries[available]
    center <- (min(values) + max(values)) / 2
    radius <- (max(values) - min(values)) / 2
    if (!is.finite(center) || !is.finite(radius) || radius <= 0) next
    ideal <- center + radius * cos(pi * (seq_len(n) - 0.5) / n)
    nearest <- vapply(ideal, function(value) available[[which.min(abs(values - value))]], integer(1L))
    nodes <- sort(unique(c(available[[1L]], available[[length(available)]], nearest[-c(1L, n)])))
    m <- length(nodes)
    if (m < 8L) next
    query_ids <- which(shared$queries >= min(values) & shared$queries <= max(values))
    query <- shared$queries[query_ids]
    anchors <- shared$queries[nodes]
    matching <- match(query_ids, nodes)
    ordinary <- which(is.na(matching))
    L <- matrix(0, length(query), m)
    safe <- rep(TRUE, length(query))
    for (j in seq_len(m)) {
      product <- rep(1, length(ordinary))
      good <- rep(TRUE, length(ordinary))
      for (k in setdiff(seq_len(m), j)) {
        denominator <- anchors[[j]] - anchors[[k]]
        if (!normal(denominator)) { good[] <- FALSE; break }
        active <- which(good)
        if (!length(active)) break
        numerator <- query[ordinary[active]] - anchors[[k]]
        valid <- normal(numerator)
        ratio <- numerator / denominator
        valid <- valid & normal(ratio)
        updated <- product[active] * ratio
        valid <- valid & normal(updated)
        good[active[!valid]] <- FALSE
        product[active[valid]] <- updated[valid]
      }
      safe[ordinary[!good]] <- FALSE
      L[ordinary[good], j] <- product[good]
    }
    exact <- which(!is.na(matching))
    if (length(exact)) L[cbind(exact, matching[exact])] <- 1
    # Actual-double differences define the interpolation nodes. Relative
    # difference/division/product error is bounded by gamma_(4*(m-1)).
    coefficient_gamma <- gamma(4 * (m - 1))
    coefficient_error <- abs(L) * coefficient_gamma / (1 - coefficient_gamma)
    if (length(exact)) coefficient_error[exact, ] <- 0
    log_product <- rep(-Inf, length(query))
    for (i in ordinary[safe[ordinary]]) {
      distance <- abs(query[[i]] - anchors)
      if (any(!normal(distance))) { safe[[i]] <- FALSE; next }
      logarithms <- log(distance)
      # Conventional log/subtraction/summation allowances, evaluated upward
      # in log space to avoid losing small products through underflow.
      allowance <- (gamma(m) + unit / (1 - unit)) * sum(abs(logarithms)) - m * log1p(-unit)
      log_product[[i]] <- sum(logarithms) + allowance
    }
    if (!any(safe)) next
    leaves[[length(leaves) + 1L]] <- list(nodes = nodes, ids = query_ids[safe], m = m,
      cardinal = L[safe, , drop = FALSE],
      absolute_cardinal = (abs(L) + coefficient_error)[safe, , drop = FALSE],
      coefficient_error = coefficient_error[safe, , drop = FALSE],
      log_product = log_product[safe])
  }
  bytes <- as.numeric(utils::object.size(leaves))
  if (bytes > shared$geometry_limit - shared$geometry_bytes) return(list())
  assign(key, leaves, shared$geometries)
  shared$geometry_bytes <- shared$geometry_bytes + bytes
  leaves
}


.selection_normalizer_grid_loglik <- function(yi, means, covariance_lower, sei,
    selection_context, execution_plan, block_size, metadata) {

  shared <- metadata[["state"]][["shared"]]
  if (!is.environment(shared)) return(NULL)
  S <- nrow(means)
  modes <- rep(selection_context[["kernel_mode"]], length.out = S)
  rules <- rep(selection_context[["vector_rule"]], length.out = S)
  omega <- as.matrix(selection_context[["omega"]])
  if (all(modes == 0L)) return(NULL)
  if (any(modes != 1L) || any(rules != 0L) || any(!is.finite(omega)) || any(omega <= 0)) {
    shared$untracked <- TRUE
    return(NULL)
  }
  # Constant positive weights cancel from the selected density exactly.
  if (all(omega == omega[, 1L])) return(NULL)
  affine <- metadata[["state"]]
  state_index <- affine[["state_index"]]
  qid <- match(affine[["values"]], shared$queries)
  if (length(qid) != S || anyNA(qid)) {
    shared$untracked <- TRUE
    return(NULL)
  }
  sign <- metadata[["sign"]]
  block_rows <- metadata[["rows"]]
  groups <- split(seq_len(S), state_index)
  entries <- vector("list", length(groups))
  requests <- integer()
  request_group <- integer()
  request_qid <- integer()
  request_kind <- character()
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    state <- state_index[positions[[1L]]]
    index <- match(affine[["rows"]][[state]], shared$rows)
    if (is.na(index)) stop("Normalizer grid rows do not match their request.", call. = FALSE)
    key <- paste(metadata[["block"]], index, sep = ":")
    packed <- as.numeric(covariance_lower[positions[[1L]], ])
    weights <- as.numeric(omega[positions[[1L]], ])
    if (any(covariance_lower[positions, , drop = FALSE] !=
        matrix(packed, length(positions), length(packed), byrow = TRUE)) ||
        any(omega[positions, , drop = FALSE] != matrix(weights, length(positions), length(weights), byrow = TRUE))) {
      shared$untracked <- TRUE
      return(NULL)
    }
    b <- sign * as.numeric(affine[["basis"]][state, block_rows])
    binding <- list(mean = sign * as.numeric(affine[["mean"]][state, block_rows]),
      b = b, current = affine[["current"]][[state]], covariance = packed,
      omega = weights, sei = as.numeric(sei), obs_bin = selection_context[["obs_bin"]])
    if (exists(key, shared$cases, inherits = FALSE)) {
      entry <- get(key, shared$cases)
      if (!identical(entry$binding, binding)) {
        shared$untracked <- TRUE
        return(NULL)
      }
    } else {
      C <- matrix(0, block_size, block_size)
      C[lower.tri(C, diag = TRUE)] <- packed
      C[upper.tri(C)] <- t(C)[upper.tri(C)]
      root <- tryCatch(chol(C), error = function(e) NULL)
      if (is.null(root)) { shared$untracked <- TRUE; return(NULL) }
      entry <- new.env(parent = emptyenv())
      entry$binding <- binding
      entry$row <- index
      entry$d <- sqrt(sum(forwardsolve(t(root), b)^2))
      entry$constant <- all(b == 0)
      if (!entry$constant && (!is.finite(entry$d) || entry$d <= 0)) {
        shared$untracked <- TRUE
        return(NULL)
      }
      entry$analytic <- all(C[upper.tri(C)] == 0)
      entry$log_scale <- block_size * log(max(weights))
      entry$range <- -expm1(block_size * log1p((min(weights) - max(weights)) / max(weights)))
      entry$data <- matrix(NA_real_, length(shared$queries), 4L,
        dimnames = list(NULL, c("log_density", "log_A", "mcse", "eta")))
      entry$have <- logical(length(shared$queries))
      entry$seen <- logical(length(shared$queries))
      entry$reported <- numeric(length(shared$queries))
      entry$unknown <- logical(length(shared$queries))
      entry$leaves <- if (entry$constant || entry$analytic || entry$range == 0) list() else
        .selection_normalizer_grid_geometry(shared, qid[positions])
      assign(key, entry, shared$cases)
    }
    entries[[g]] <- entry
    entry$seen[unique(qid[positions])] <- TRUE
    anchors <- if (entry$constant && !any(entry$have)) qid[positions[[1L]]] else
      unique(unlist(lapply(entry$leaves, `[[`, "nodes"), use.names = FALSE))
    anchors <- anchors[!entry$have[anchors]]
    if (length(anchors)) {
      # Every new anchor is an actually requested point in this call.
      at <- match(anchors, qid[positions])
      if (anyNA(at)) stop("Normalizer anchors were not requested.", call. = FALSE)
      requests <- c(requests, positions[at])
      request_group <- c(request_group, rep(g, length(at)))
      request_qid <- c(request_qid, anchors)
      request_kind <- c(request_kind, rep(if (entry$constant) "constant" else "anchor", length(at)))
    }
  }
  fetch <- function(positions, group, ids, kind) {
    if (!length(positions)) return(invisible(NULL))
    result <- .selection_joint_dense_loglik_block(yi, means[positions, , drop = FALSE],
      covariance_lower[positions, , drop = FALSE], sei,
      BayesTools::selection_context_subset_rows(selection_context, positions),
      execution_plan, block_size, return_normalizer = TRUE)
    diagnostic <- result[["integration_diagnostics"]]
    cdf <- attr(diagnostic, "cdf_relative_error", exact = TRUE)
    if (is.null(cdf)) cdf <- numeric(length(positions))
    eta <- 2 * diagnostic[, "quadrature_change"] + diagnostic[, "covariance_width"] +
      diagnostic[, "tail_bound"] + cdf
    eta[diagnostic[, "used_covariance_envelope"] != 1] <- NA_real_
    for (g in unique(group)) {
      at <- which(group == g)
      entry <- entries[[g]]
      if (entry$analytic) eta[at] <- 0
      value <- cbind(result$log_density[at], result$log_normalizer[at], result$relative_mcse[at], eta[at])
      rows <- if (entry$constant) seq_along(entry$have) else ids[at]
      if (entry$constant) value <- matrix(value[1L, ], length(rows), 4L, byrow = TRUE)
      entry$data[rows, ] <- value
      entry$have[rows] <- TRUE
    }
    for (name in c("anchor", "fallback", "constant")) {
      field <- paste0(name, "_evaluations")
      shared$stats[[field]] <- shared$stats[[field]] + sum(kind == name)
    }
    shared$stats[["max_native_mcse"]] <- max(shared$stats[["max_native_mcse"]], result$relative_mcse)
    invisible(NULL)
  }
  fetch(requests, request_group, request_qid, request_kind)
  log_A <- rep(NA_real_, S)
  error <- rep(NA_real_, S)
  result <- rep(NA_real_, S)
  mcse <- numeric(S)
  interpolated <- logical(S)
  requests <- integer(); request_group <- integer(); request_qid <- integer()
  tolerance <- execution_plan[["relative_tolerance"]]
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    entry <- entries[[g]]
    for (leaf in entry$leaves) {
      anchors <- entry$data[leaf$nodes, , drop = FALSE]
      A <- exp(anchors[, "log_A"] - entry$log_scale)
      if (any(!is.finite(A)) || any(A < .Machine$double.xmin) ||
          any(!is.finite(anchors[, "eta"])) || any(anchors[, "eta"] < 0) ||
          any(anchors[, "eta"] > tolerance) || any(anchors[, "mcse"] != 0)) next
      ids <- intersect(qid[positions], leaf$ids)
      ids <- ids[!entry$have[ids]]
      if (!length(ids)) next
      within <- match(ids, leaf$ids)
      L <- leaf$cardinal[within, , drop = FALSE]
      absolute <- leaf$absolute_cardinal[within, , drop = FALSE]
      P <- as.numeric(L %*% A)
      unit <- .Machine$double.eps
      conversion <- expm1(unit * (abs(anchors[, "log_A"]) + abs(entry$log_scale)) - log1p(-unit))
      anchor_error <- as.numeric(absolute %*% (A * (anchors[, "eta"] + conversion * (1 + anchors[, "eta"]))))
      terms <- c(log(entry$range / 2), leaf$m * log(entry$d), -0.5 * lgamma(leaf$m + 1))
      log_remainder <- sum(terms) + leaf$log_product[within]
      log_remainder <- log_remainder + (5 * unit / (1 - 5 * unit)) *
        (sum(abs(terms)) + abs(leaf$log_product[within]))
      remainder <- exp(log_remainder) / (1 - unit)
      dot_gamma <- 2 * leaf$m * unit / (1 - 2 * leaf$m * unit)
      absolute_products <- sweep(abs(L), 2L, A, `*`)
      dot_safe <- rowSums((abs(L) > 0) &
        (!is.finite(absolute_products) | absolute_products < .Machine$double.xmin)) == 0
      rounding <- as.numeric(leaf$coefficient_error[within, , drop = FALSE] %*% A) +
        dot_gamma / (1 - dot_gamma) * rowSums(absolute_products)
      total <- remainder + anchor_error + rounding
      accept <- dot_safe & remainder >= .Machine$double.xmin & is.finite(P) & P > 0 &
        is.finite(total) & total < P & total / P <= tolerance
      if (!any(accept)) next
      used <- positions[qid[positions] %in% ids[accept]]
      which_value <- match(qid[used], ids)
      log_A[used] <- log(P[which_value]) + entry$log_scale
      error[used] <- total[which_value] / P[which_value]
      interpolated[used] <- TRUE
      shared$stats[["max_anchor_relative_error"]] <- max(shared$stats[["max_anchor_relative_error"]],
        anchor_error[accept] / P[accept])
      shared$stats[["max_interpolation_relative_error"]] <- max(shared$stats[["max_interpolation_relative_error"]],
        (remainder[accept] + rounding[accept]) / P[accept])
    }
    pending <- positions[!entry$have[qid[positions]] & !interpolated[positions]]
    pending <- pending[!duplicated(qid[pending])]
    requests <- c(requests, pending)
    request_group <- c(request_group, rep(g, length(pending)))
    request_qid <- c(request_qid, qid[pending])
  }
  fetch(requests, request_group, request_qid, rep("fallback", length(requests)))
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    entry <- entries[[g]]
    direct <- positions[!interpolated[positions]]
    result[direct] <- entry$data[qid[direct], "log_density"]
    error[direct] <- entry$data[qid[direct], "eta"]
    mcse[direct] <- entry$data[qid[direct], "mcse"]
  }
  .selection_joint_dense_loglik_check_mcse(mcse, execution_plan)
  if (any(interpolated)) {
    rows <- which(interpolated)
    gaussian_context <- BayesTools::selection_context_subset_rows(selection_context, rows)
    # This requests only the exact Gaussian component; fitted bias metadata
    # and the selected normalizer remain in the original context.
    gaussian_context[["kernel_mode"]] <- rep(0L, length(rows))
    gaussian <- .selection_joint_dense_loglik_block(yi, means[rows, , drop = FALSE],
      covariance_lower[rows, , drop = FALSE], sei, gaussian_context,
      execution_plan, block_size)
    for (bin in selection_context[["obs_bin"]]) gaussian <- gaussian + log(omega[rows, bin])
    result[rows] <- gaussian - log_A[rows]
    shared$stats[["gaussian_evaluations"]] <- shared$stats[["gaussian_evaluations"]] + length(rows)
    shared$stats[["interpolated_points"]] <- shared$stats[["interpolated_points"]] + length(rows)
  }
  if (any(is.finite(error))) shared$stats[["max_relative_error"]] <-
    max(shared$stats[["max_relative_error"]], error[is.finite(error)])
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    entry <- entries[[g]]
    # Repeated ids share one cached/computed normalizer and its error.
    at <- positions[!duplicated(qid[positions])]
    ids <- qid[at]
    values <- error[at]
    bad <- !is.finite(values) | values < 0 | values >= 1
    entry$unknown[ids[bad]] <- TRUE
    entry$reported[ids[!bad]] <- pmax(entry$reported[ids[!bad]], -log1p(-values[!bad]))
  }
  result
}

.selection_normalizer_grid_diagnostics <- function(shared) {

  if (!is.environment(shared)) return(NULL)
  total <- matrix(0, length(shared$queries), length(shared$rows))
  seen <- unknown <- matrix(FALSE, nrow(total), ncol(total))
  unique_requests <- 0
  for (name in ls(shared$cases, all.names = TRUE)) {
    entry <- get(name, shared$cases)
    total[, entry$row] <- total[, entry$row] + entry$reported
    seen[, entry$row] <- seen[, entry$row] | entry$seen
    unknown[, entry$row] <- unknown[, entry$row] | entry$unknown
    unique_requests <- unique_requests + sum(entry$seen)
  }
  error_known <- any(seen) && !shared$untracked && !any(unknown & seen)
  log_error <- if (error_known) max(total) else NA_real_
  c(list(used = shared$stats[["interpolated_points"]] > 0,
    error_scope = "finite requested ordinates and their positive quadrature normalization; native errors remain estimated, with actual-node interpolation and ordinary floating allowances",
    unique_requested_normalizers = unique_requests,
    max_log_likelihood_error = log_error,
    conditional_relative_error = if (is.finite(log_error)) expm1(2 * log_error) else NA_real_,
    unknown_error_points = sum(unknown & seen), untracked_path = shared$untracked),
    as.list(shared$stats))
}
