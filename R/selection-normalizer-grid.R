# Call-owned approximation of fixed scalar normalizer grids. Native anchors
# are actual requested points; the original priors and Gaussian law stay intact.
.selection_normalizer_grid <- function(context, queries, rows) {

  if (!.is_data_joint_selection(context[["data"]]) ||
      !.is_data_known_v(context[["data"]]) ||
      .selection_retains_sampling(context[["data"]])) return(NULL)
  plan <- .data_selection_execution_plan(context[["data"]])
  # The transport identity and its certificate do not depend on how a block
  # normalizer is evaluated, only on the accuracy each route reports, so every
  # correlated route is eligible.
  tracked <- sum(plan[["block_methods"]] %in%
    c("dense", "rank_one", "factor"))
  queries <- sort(unique(queries[is.finite(queries)]))
  if (!tracked || length(queries) < 16L || !length(rows)) return(NULL)
  # Reserve for native records, reported error, binding data and R overhead.
  bytes <- 4 * (length(queries) * length(rows) * tracked * 48 +
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


# The packed lower triangle of diag(residual_sd^2) + L L' for every draw, in
# the column-major order the grid's block reconstruction expects.
.selection_factor_covariance_lower <- function(residual_sd, loading, block_size) {

  S     <- nrow(residual_sd)
  rank  <- if (block_size == 0L) 0L else ncol(loading) %/% block_size
  index <- which(lower.tri(matrix(0, block_size, block_size), diag = TRUE),
                 arr.ind = TRUE)
  rows  <- index[, 1L]
  cols  <- index[, 2L]
  out   <- matrix(0, nrow = S, ncol = nrow(index))
  # One matrix-indexed pass per packed column set instead of a loop per packed
  # column: every element keeps the diagonal term first and then the factor
  # products in increasing factor order, so the additions are unchanged.
  diagonal <- which(rows == cols)
  if (length(diagonal)) {
    out[, diagonal] <- residual_sd[, rows[diagonal], drop = FALSE]^2
  }
  for (factor in seq_len(rank)) {
    offset <- (factor - 1L) * block_size
    out <- out + loading[, offset + rows, drop = FALSE] *
      loading[, offset + cols, drop = FALSE]
  }

  out
}


# Packed block covariances of a candidate batch, kept at one row per distinct
# state when the batch only repeats one posterior row per state. The packed
# values are a row-wise function of the factor inputs, so representing repeated
# input rows once leaves every packed value unchanged; a batch whose factor
# inputs differ inside a state keeps the full per-row matrix and therefore the
# per-group covariance comparison in .selection_normalizer_grid_loglik().
.selection_factor_covariance_rows <- function(residual_sd, loading, block_size,
                                              state_index) {

  S       <- nrow(residual_sd)
  compact <- NULL
  if (S > 1L && !is.null(state_index) && length(state_index) == S) {
    first <- match(state_index, state_index)
    same  <- !any(residual_sd != residual_sd[first, , drop = FALSE]) &&
      !any(loading != loading[first, , drop = FALSE])
    if (isTRUE(same)) compact <- first
  }
  if (is.null(compact)) {
    return(list(
      values = .selection_factor_covariance_lower(residual_sd, loading, block_size),
      rows   = NULL
    ))
  }
  distinct <- which(compact == seq_len(S))

  list(
    values = .selection_factor_covariance_lower(
      residual_sd = residual_sd[distinct, , drop = FALSE],
      loading     = loading[distinct, , drop = FALSE],
      block_size  = block_size
    ),
    rows = match(compact, distinct)
  )
}


# Run starts of consecutive draws sharing a packed covariance. Exact equality
# is transitive, so comparing neighbours marks the runs the per-draw comparison
# against the run's own covariance marked. A missing comparison keeps the
# per-draw test, whose "missing value" error the caller saw before.
.selection_grid_covariance_runs <- function(covariance_lower, S) {

  if (S <= 1L) return(seq_len(S))
  difference <- covariance_lower[-1L, , drop = FALSE] !=
    covariance_lower[-S, , drop = FALSE]
  if (!anyNA(difference)) {
    return(c(1L, which(rowSums(difference) > 0) + 1L))
  }
  starts   <- 1L
  previous <- covariance_lower[1L, ]
  for (draw in seq.int(2L, S)) {
    packed <- covariance_lower[draw, ]
    if (any(packed != previous)) {
      starts   <- c(starts, draw)
      previous <- packed
    }
  }

  starts
}


# The exact Gaussian component of each block state. Interpolated normalizers
# divide into this density, so it is evaluated from the same covariance the
# block's own route integrates. `covariance_rows` maps result rows to rows of
# `covariance_lower`; NULL is the identity mapping.
.selection_grid_gaussian_lpdf <- function(yi, means, covariance_lower, block_size,
                                          covariance_rows = NULL) {

  S <- nrow(means)
  if (!S) return(numeric(0))
  out       <- numeric(S)
  lower     <- lower.tri(matrix(0, block_size, block_size), diag = TRUE)
  upper     <- upper.tri(matrix(0, block_size, block_size))
  packed_at <- matrix(0L, block_size, block_size)
  packed_at[lower] <- seq_len(sum(lower))
  packed_at[upper] <- t(packed_at)[upper]
  # Residuals for the whole batch and the transpose it needs are built once.
  residuals <- yi - t(means)
  index     <- if (is.null(covariance_rows)) seq_len(S) else covariance_rows
  starts    <- if (is.null(covariance_rows)) {
    .selection_grid_covariance_runs(covariance_lower, S)
  } else if (S == 1L) 1L else c(1L, which(index[-1L] != index[-S]) + 1L)
  ends      <- c(starts[-1L] - 1L, S)
  # Consecutive draws that share a packed covariance are factored once and
  # whitened together. backsolve() solves each right-hand-side column with the
  # same triangular arithmetic as the single-column call, and colSums() sums the
  # squares of one column in the same order, so the values are unchanged.
  for (segment in seq_along(starts)) {
    packed <- covariance_lower[index[[starts[[segment]]]], ]
    root   <- tryCatch(chol(matrix(packed[as.integer(packed_at)],
                                   block_size, block_size)),
                       error = function(e) NULL)
    if (is.null(root)) {
      return(NULL)
    }
    run      <- starts[[segment]]:ends[[segment]]
    whitened <- backsolve(root, residuals[, run, drop = FALSE], transpose = TRUE)
    out[run] <- -0.5 * (block_size * log(2 * pi) + 2 * sum(log(diag(root))) +
      colSums(whitened^2))
  }

  out
}


.selection_normalizer_grid_loglik <- function(yi, means, covariance_lower, sei,
    selection_context, execution_plan, block_size, metadata, factor = NULL,
    covariance_rows = NULL) {

  shared <- metadata[["state"]][["shared"]]
  if (!is.environment(shared)) return(NULL)
  if (!is.null(covariance_rows) && is.null(factor)) {
    stop("A compacted normalizer-grid covariance requires a factor block.",
         call. = FALSE)
  }
  S <- nrow(means)
  # Recycling only repeats these fields, so a source no longer than the batch
  # already holds every value the recycled vector would; anything else keeps
  # the explicit expansion.
  modes <- selection_context[["kernel_mode"]]
  rules <- selection_context[["vector_rule"]]
  if (!length(modes) || length(modes) > S) modes <- rep(modes, length.out = S)
  if (!length(rules) || length(rules) > S) rules <- rep(rules, length.out = S)
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
  # Everything a state group needs that does not depend on the group's own
  # candidate rows is prepared once for the whole call: the grid row of each
  # state, its cache key, its packed covariance, weights and affine row, and
  # the group's query ids.
  first_positions <- vapply(groups, `[[`, integer(1L), 1L)
  group_states    <- state_index[first_positions]
  state_rows      <- match(affine[["rows"]], shared$rows)
  group_rows      <- state_rows[group_states]
  if (anyNA(group_rows)) stop("Normalizer grid rows do not match their request.", call. = FALSE)
  group_keys    <- paste(metadata[["block"]], group_rows, sep = ":")
  group_packed  <- covariance_lower[
    if (is.null(covariance_rows)) first_positions else covariance_rows[first_positions], ,
    drop = FALSE]
  group_weights <- omega[first_positions, , drop = FALSE]
  group_mean    <- sign * affine[["mean"]][group_states, block_rows, drop = FALSE]
  group_basis   <- sign * affine[["basis"]][group_states, block_rows, drop = FALSE]
  group_current <- affine[["current"]][group_states]
  group_qid     <- lapply(groups, function(positions) qid[positions])
  group_unique  <- lapply(group_qid, unique)
  sei_values    <- as.numeric(sei)
  obs_bin       <- selection_context[["obs_bin"]]
  requests      <- vector("list", length(groups))
  request_group <- vector("list", length(groups))
  request_qid   <- vector("list", length(groups))
  request_kind  <- vector("list", length(groups))
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    index <- group_rows[[g]]
    key <- group_keys[[g]]
    packed <- as.numeric(group_packed[g, ])
    weights <- as.numeric(group_weights[g, ])
    # A compacted covariance carries one packed row per state by construction:
    # it is only built after the factor inputs of every state were shown to
    # agree, which is what the per-row comparison below establishes otherwise.
    covariance_differs <- if (!is.null(covariance_rows)) FALSE else
      any(covariance_lower[positions, , drop = FALSE] !=
        matrix(packed, length(positions), length(packed), byrow = TRUE))
    if (covariance_differs ||
        any(omega[positions, , drop = FALSE] != matrix(weights, length(positions), length(weights), byrow = TRUE))) {
      shared$untracked <- TRUE
      return(NULL)
    }
    b <- as.numeric(group_basis[g, ])
    binding <- list(mean = as.numeric(group_mean[g, ]),
      b = b, current = group_current[[g]], covariance = packed,
      omega = weights, sei = sei_values, obs_bin = obs_bin)
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
        .selection_normalizer_grid_geometry(shared, group_qid[[g]])
      # The anchor set of a case is fixed with its geometry.
      entry$anchor_ids <- unique(unlist(lapply(entry$leaves, `[[`, "nodes"),
        use.names = FALSE))
      assign(key, entry, shared$cases)
    }
    entries[[g]] <- entry
    entry$seen[group_unique[[g]]] <- TRUE
    anchors <- if (entry$constant && !any(entry$have)) qid[positions[[1L]]] else
      entry$anchor_ids
    anchors <- anchors[!entry$have[anchors]]
    if (length(anchors)) {
      # Every new anchor is an actually requested point in this call.
      at <- match(anchors, group_qid[[g]])
      if (anyNA(at)) stop("Normalizer anchors were not requested.", call. = FALSE)
      requests[[g]] <- positions[at]
      request_group[[g]] <- rep(g, length(at))
      request_qid[[g]] <- anchors
      request_kind[[g]] <- rep(if (entry$constant) "constant" else "anchor", length(at))
    }
  }
  requests <- unlist(requests, use.names = FALSE)
  request_group <- unlist(request_group, use.names = FALSE)
  request_qid <- unlist(request_qid, use.names = FALSE)
  request_kind <- unlist(request_kind, use.names = FALSE)
  if (is.null(requests)) {
    requests <- integer(); request_group <- integer()
    request_qid <- integer(); request_kind <- character()
  }
  fetch <- function(positions, group, ids, kind) {
    if (!length(positions)) return(invisible(NULL))
    context <- BayesTools::selection_context_subset_rows(selection_context, positions)
    result <- if (is.null(factor)) {
      .selection_joint_dense_loglik_block(yi, means[positions, , drop = FALSE],
        covariance_lower[positions, , drop = FALSE], sei, context,
        execution_plan, block_size, return_normalizer = TRUE)
    } else if (identical(factor[["method"]], "rank_one")) {
      .selection_joint_cluster_loglik_block(yi, means[positions, , drop = FALSE],
        factor[["residual_sd"]][positions, , drop = FALSE],
        factor[["loading"]][positions, , drop = FALSE], sei, context,
        execution_plan, return_normalizer = TRUE)
    } else {
      .selection_joint_factor_loglik_block(yi, means[positions, , drop = FALSE],
        factor[["residual_sd"]][positions, , drop = FALSE],
        factor[["loading"]][positions, , drop = FALSE], sei, context,
        execution_plan, factor[["block_index"]], return_normalizer = TRUE)
    }
    eta <- if (is.null(factor)) {
      diagnostic <- result[["integration_diagnostics"]]
      cdf <- attr(diagnostic, "cdf_relative_error", exact = TRUE)
      if (is.null(cdf)) cdf <- numeric(length(positions))
      value <- 2 * diagnostic[, "quadrature_change"] + diagnostic[, "covariance_width"] +
        diagnostic[, "tail_bound"] + cdf
      value[diagnostic[, "used_covariance_envelope"] != 1] <- NA_real_
      value
    } else {
      # The factor routes certify their successive quadrature change; keep the
      # same doubling the dense envelope applies to its own change term.
      value <- 2 * result[["relative_change"]]
      value[!is.finite(value)] <- NA_real_
      value
    }
    # The request vectors are built group by group in increasing group order,
    # so splitting their positions visits the same groups in the same order as
    # the previous scan over unique(group), without a pass per group.
    for (at in split(seq_along(group), group)) {
      g <- group[[at[[1L]]]]
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
  requests <- vector("list", length(groups))
  request_group <- vector("list", length(groups))
  request_qid <- vector("list", length(groups))
  tolerance <- execution_plan[["relative_tolerance"]]
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    own_qid <- group_qid[[g]]
    entry <- entries[[g]]
    for (leaf in entry$leaves) {
      anchors <- entry$data[leaf$nodes, , drop = FALSE]
      A <- exp(anchors[, "log_A"] - entry$log_scale)
      if (any(!is.finite(A)) || any(A < .Machine$double.xmin) ||
          any(!is.finite(anchors[, "eta"])) || any(anchors[, "eta"] < 0) ||
          any(anchors[, "eta"] > tolerance) || any(anchors[, "mcse"] != 0)) next
      # intersect() on these two plain integer vectors, whose second argument
      # is already unique and whose first is uniquified once per group.
      unique_qid <- group_unique[[g]]
      ids <- unique_qid[match(unique_qid, leaf$ids, 0L) > 0L]
      ids <- ids[!entry$have[ids]]
      if (!length(ids)) next
      within <- match(ids, leaf$ids)
      L <- leaf$cardinal[within, , drop = FALSE]
      absolute_L <- abs(L)
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
      absolute_products <- absolute_L * rep(A, each = nrow(absolute_L))
      dot_safe <- rowSums((absolute_L > 0) &
        (!is.finite(absolute_products) | absolute_products < .Machine$double.xmin)) == 0
      rounding <- as.numeric(leaf$coefficient_error[within, , drop = FALSE] %*% A) +
        dot_gamma / (1 - dot_gamma) * rowSums(absolute_products)
      total <- remainder + anchor_error + rounding
      accept <- dot_safe & remainder >= .Machine$double.xmin & is.finite(P) & P > 0 &
        is.finite(total) & total < P & total / P <= tolerance
      if (!any(accept)) next
      at_used <- own_qid %in% ids[accept]
      used <- positions[at_used]
      which_value <- match(own_qid[at_used], ids)
      log_A[used] <- log(P[which_value]) + entry$log_scale
      error[used] <- total[which_value] / P[which_value]
      interpolated[used] <- TRUE
      shared$stats[["max_anchor_relative_error"]] <- max(shared$stats[["max_anchor_relative_error"]],
        anchor_error[accept] / P[accept])
      shared$stats[["max_interpolation_relative_error"]] <- max(shared$stats[["max_interpolation_relative_error"]],
        (remainder[accept] + rounding[accept]) / P[accept])
    }
    at_pending <- !entry$have[own_qid] & !interpolated[positions]
    at_pending[at_pending] <- !duplicated(own_qid[at_pending])
    if (any(at_pending)) {
      requests[[g]] <- positions[at_pending]
      request_group[[g]] <- rep(g, sum(at_pending))
      request_qid[[g]] <- own_qid[at_pending]
    }
  }
  requests <- unlist(requests, use.names = FALSE)
  request_group <- unlist(request_group, use.names = FALSE)
  request_qid <- unlist(request_qid, use.names = FALSE)
  if (is.null(requests)) {
    requests <- integer(); request_group <- integer(); request_qid <- integer()
  }
  fetch(requests, request_group, request_qid, rep("fallback", length(requests)))
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    own_qid <- group_qid[[g]]
    entry <- entries[[g]]
    at_direct <- !interpolated[positions]
    direct <- positions[at_direct]
    values <- entry$data[own_qid[at_direct], , drop = FALSE]
    result[direct] <- values[, "log_density"]
    error[direct] <- values[, "eta"]
    mcse[direct] <- values[, "mcse"]
  }
  .selection_joint_dense_loglik_check_mcse(mcse, execution_plan)
  if (any(interpolated)) {
    rows <- which(interpolated)
    gaussian_context <- BayesTools::selection_context_subset_rows(selection_context, rows)
    # This requests only the exact Gaussian component; fitted bias metadata
    # and the selected normalizer remain in the original context.
    gaussian_context[["kernel_mode"]] <- rep(0L, length(rows))
    gaussian <- if (is.null(factor)) {
      .selection_joint_dense_loglik_block(yi, means[rows, , drop = FALSE],
        covariance_lower[rows, , drop = FALSE], sei, gaussian_context,
        execution_plan, block_size)
    } else if (is.null(covariance_rows)) {
      .selection_grid_gaussian_lpdf(yi, means[rows, , drop = FALSE],
        covariance_lower[rows, , drop = FALSE], block_size)
    } else {
      .selection_grid_gaussian_lpdf(yi, means[rows, , drop = FALSE],
        covariance_lower, block_size, covariance_rows = covariance_rows[rows])
    }
    if (is.null(gaussian)) {
      shared$untracked <- TRUE
      return(NULL)
    }
    for (bin in selection_context[["obs_bin"]]) gaussian <- gaussian + log(omega[rows, bin])
    result[rows] <- gaussian - log_A[rows]
    shared$stats[["gaussian_evaluations"]] <- shared$stats[["gaussian_evaluations"]] + length(rows)
    shared$stats[["interpolated_points"]] <- shared$stats[["interpolated_points"]] + length(rows)
  }
  if (any(is.finite(error))) shared$stats[["max_relative_error"]] <-
    max(shared$stats[["max_relative_error"]], error[is.finite(error)])
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    own_qid <- group_qid[[g]]
    entry <- entries[[g]]
    # Repeated ids share one cached/computed normalizer and its error.
    first <- !duplicated(own_qid)
    ids <- own_qid[first]
    values <- error[positions[first]]
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
