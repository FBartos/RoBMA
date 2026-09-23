# ============================================================================ #
# Selection Normalizer Interpolation along Single-Component SD Families
# ============================================================================ #
#
# Call-owned certified interpolation of selection normalizers along the
# one-parameter families traced by single random-component SD sweeps:
#
# 1. "variance" (integrated components): Sigma(tau) = D + tau^2 m m' with D
#    the block covariance with the selected component removed and m the
#    component's indicator column restricted to the block.
# 2. "mean" (conditioned components, the default for study-level effects):
#    Sigma fixed and mu(tau) = mu_0 + (tau - tau_0) b with b the component's
#    standardized latent draws.
#
# Both families are declared by the BayesTools affine marginal-update metadata
# for single-column component-SD sweeps; the structure is never inferred from
# evaluated covariance matrices.
#
# Derivation: conditioning on the latent factor z ~ N(0,1) of the rank-one
# augmentation (family 1) or of the latent direction (family 2) gives the
# transport identity A(tau) = E_z[a(tau z)] with a(t) the mean-shift
# normalizer of the fixed base covariance along the direction (m resp. b).
# Iterated Stein's lemma on X ~ N(0, Sigma) bounds
# |a^(j)(t)| <= (R/2) d^j sqrt(j!) with d = ||Sigma^{-1/2} direction|| and R
# the range of the product weight; differentiating the transport under the
# integral adds the exact absolute moment
# E|z|^j = 2^{j/2} Gamma((j+1)/2) / sqrt(pi), so
#
#   |A^(j)(tau)| <= (R/2) d^j sqrt(j!) E|z|^j.
#
# E|z|^j / sqrt(j!) < 1 for j >= 2, so the interpolation remainder keeps the
# Chebyshev-style decay of the reviewed mean-shift grid. Native anchors are
# actual requested grid points; cardinal interpolation, rounding allowances,
# and the per-point acceptance gate mirror that implementation. One extra
# eps-scale term covers the generic builder's floating-point arithmetic
# against the exact family through a first-order selection-sensitivity bound.

# Shared-state creation and eligibility. Mirrors `.selection_normalizer_grid`.
.selection_covariance_grid <- function(context, queries, rows, row_states,
                                       replacement, parameter) {

  data <- context[["data"]]
  if (!.is_data_joint_selection(data) || !.is_data_known_v(data) ||
      .selection_retains_sampling(data) ||
      .data_outcome_type(data) != "norm" || .is_data_weights(data)) {
    return(NULL)
  }
  update <- replacement[["covariance_update"]]
  if (!inherits(update, "BayesTools_random_effects_marginal_update_plan") ||
      !identical(update[["family"]], "affine") ||
      !replacement[["type"]] %in% c("primitive", "scalar", "random_component_sd")) {
    return(NULL)
  }
  block_name <- update[["blocks"]]
  if (!is.character(block_name) || length(block_name) != 1L ||
      is.na(block_name) || !nzchar(block_name)) {
    return(NULL)
  }
  plan <- .data_selection_execution_plan(data)
  # Only the dense route. Extending the interpolation to the rank-one and
  # factor routes was implemented and measured: their kernels are cheap enough
  # after the vectorized tails and the budgeted quadrature that the anchor
  # bookkeeping costs about four times what it saves. See the follow-up
  # implementation log.
  if (!any(plan[["block_methods"]] == "dense")) {
    return(NULL)
  }
  integrated <- isTRUE(block_name %in% plan[["random_covariance"]][["term_names"]])
  loading <- if (integrated) {
    .selection_covariance_grid_loading(context, block_name)
  } else {
    list(model_matrix = NULL)
  }
  if (integrated && is.null(loading)) {
    return(NULL)
  }
  queries <- sort(unique(queries[is.finite(queries) & queries >= 0]))
  if (length(queries) < 16L || !length(rows)) {
    return(NULL)
  }
  bytes <- 4 * (length(queries) * length(rows) * 8 +
    length(rows) * sum(lengths(plan[["row_blocks"]])^2) * 8)
  if (!is.finite(bytes) || bytes > .known_v_covariance_max_bytes()) {
    return(NULL)
  }
  out <- new.env(parent = emptyenv())
  out$queries <- queries
  out$rows <- rows
  out$row_states <- stats::setNames(row_states, as.character(rows))
  out$update <- update
  out$parameter <- parameter
  out$replacement <- replacement
  out$family <- if (integrated) "variance" else "mean"
  out$loading <- loading
  out$cases <- new.env(parent = emptyenv())
  out$geometries <- new.env(parent = emptyenv())
  out$bases <- new.env(parent = emptyenv())
  out$stats <- c(anchor_evaluations = 0, fallback_evaluations = 0,
    constant_evaluations = 0, gaussian_evaluations = 0, interpolated_points = 0,
    max_relative_error = 0, max_anchor_relative_error = 0,
    max_interpolation_relative_error = 0, max_native_mcse = 0)
  out$untracked <- FALSE
  out$geometry_limit <- .known_v_covariance_max_bytes() - bytes
  out$geometry_bytes <- 0
  out
}


# Compile the selected component's model matrix from the formula design. The
# entries must be exact 0/1 indicators; per selection block at most one column
# may be nonzero for the rank-one basis certificate. Compiled once per grid
# evaluation and cached in the iwmde context.
.selection_covariance_grid_loading <- function(context, block_name) {

  cache <- context[["covariance_grid_loading"]]
  if (!is.null(cache) && identical(cache[["block_name"]], block_name)) {
    return(cache[["loading"]])
  }
  design <- tryCatch(
    .fitted_formula_design(context[["object"]], "mu", required = TRUE),
    error = function(e) NULL
  )
  loading <- NULL
  if (!is.null(design)) {
    term <- Filter(function(term) {
      identical(term[["block_name"]], block_name)
    }, design[["random_effects"]])
    if (length(term) == 1L) {
      model_matrix <- term[[1L]][["model_matrix"]]
      K <- nrow(context[["data"]][["outcome"]])
      if (is.matrix(model_matrix) && nrow(model_matrix) == K &&
          all(model_matrix %in% c(0, 1)) && all(is.finite(model_matrix))) {
        loading <- model_matrix
      }
    }
  }
  context[["covariance_grid_loading"]] <- list(
    block_name = block_name, loading = loading)
  loading
}


# The state's current family quantity on the swept SD scale, from the update
# plan's recorded source coordinate and transform.
.selection_covariance_grid_quantity <- function(shared, samples_row) {

  update <- shared[["update"]]
  source <- update[["source_parameter"]]
  # The plan records the transform's name and the transform itself; only the
  # latter can be applied to a draw.
  transform <- update[["source_transform_spec"]]
  if (!is.character(source) || length(source) != 1L || is.na(source) ||
      !source %in% names(samples_row)) {
    return(NA_real_)
  }
  value <- samples_row[[source]]
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    return(NA_real_)
  }
  if (is.null(transform)) {
    return(value)
  }
  tryCatch(
    BayesTools::parameter_transform_forward(value, transform),
    error = function(e) NA_real_
  )
}


# Fixed per-state family ingredients, identical across the grid evaluation's
# calls because they are built from the state's original samples only.
#
# "variance": the total block covariance (which embeds the component's current
# coefficient) together with that coefficient.
# "mean": the fixed total block covariance, the family origin mean, and the
# latent direction obtained from one canonical perturbed predictor build.
.selection_covariance_grid_base <- function(shared, context, active_setup,
                                            state_row) {

  key <- as.character(state_row)
  if (exists(key, shared$bases, inherits = FALSE)) {
    return(get(key, shared$bases, inherits = FALSE))
  }
  data <- context[["data"]]
  plan <- .data_selection_execution_plan(data)
  object <- list(
    fit    = context[["object"]][["fit"]],
    data   = data,
    priors = active_setup[["priors"]]
  )
  samples <- context[["posterior_samples"]][state_row, , drop = FALSE]
  random_vcov <- tryCatch(
    .brma_mv_random_effects_marginal_vcov(
      object            = object,
      posterior_samples = samples,
      blocks            = plan[["random_covariance"]][["term_names"]]
    ),
    error = function(e) NULL
  )
  out <- NULL
  if (!is.null(random_vcov)) {
    samples_vcov <- random_vcov[["samples"]]
    K <- nrow(data[["outcome"]])
    if (identical(dim(samples_vcov), c(1L, K, K)) && all(is.finite(samples_vcov))) {
      total <- samples_vcov[1L, , ] +
        .known_v_covariance_matrix(.data_known_v_data(data))
      quantity <- .selection_covariance_grid_quantity(shared, samples[1L, ])
      if (is.finite(quantity)) {
        if (identical(shared[["family"]], "variance")) {
          out <- list(covariance = total, coefficient = quantity^2)
        } else {
          state <- shared[["row_states"]][[as.character(state_row)]]
          target <- if (quantity > 0) 2 * quantity else quantity + 1
          replaced <- tryCatch(
            .iwmde_build_replacement_samples(
              context     = context,
              parameter   = shared[["parameter"]],
              values      = target,
              row_states  = list(state),
              replacement = shared[["replacement"]]
            ),
            error = function(e) NULL
          )
          direction <- NULL
          if (!is.null(replaced) && isTRUE(all(replaced[["valid"]]))) {
            mu_0 <- .selection_covariance_grid_predictor(
              context, active_setup, samples)
            mu_1 <- .selection_covariance_grid_predictor(
              context, active_setup, replaced[["samples"]])
            if (!is.null(mu_0) && !is.null(mu_1)) {
              direction <- (mu_1 - mu_0) / (target - quantity)
            }
          }
          if (!is.null(direction) && all(is.finite(direction)) &&
              any(direction != 0)) {
            out <- list(covariance = total, coefficient = quantity^2,
                        mu0 = mu_0, direction = direction)
          }
        }
      }
    }
  }
  assign(key, out, shared$bases)
  out
}


# Full-predictor means for one state's samples, including conditioned latent
# effects, matching the joint-likelihood setup's mean construction.
.selection_covariance_grid_predictor <- function(context, active_setup, samples) {

  tryCatch({
    conditioned <- .iwmde_conditioned_random_effects_from_latent(
      context           = context,
      posterior_samples = samples,
      unit              = "estimate"
    )
    setup <- .log_lik_posterior_setup(
      fit                        = context[["object"]][["fit"]],
      posterior_samples          = samples,
      data                       = context[["data"]],
      priors                     = active_setup[["priors"]],
      unit                       = "estimate",
      data_hash                  = NULL,
      conditioned_random_effects = conditioned
    )
    if (identical(setup[["effect_direction"]], "negative")) {
      return(-setup[["mu"]][1L, ])
    }
    setup[["mu"]][1L, ]
  }, error = function(e) NULL)
}


# Outward-rounded whitened norm of a direction against the fixed block
# covariance. The relative inflation covers the Cholesky and solve roundings.
.selection_covariance_grid_direction_norm <- function(S, direction) {

  root <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(root)) {
    return(NA_real_)
  }
  solved <- backsolve(root, direction, transpose = TRUE)
  k <- length(direction)
  unit <- .Machine$double.eps
  gamma <- (3 * k + 8) * unit / (1 - (3 * k + 8) * unit)
  sqrt(sum(solved^2)) * (1 + gamma)
}


# First-order selection-sensitivity bound for eps-scale departures from the
# exact family. The Price-style arithmetic mirrors
# `covariance_price_gap_bound`'s per-pair prefactors.
.selection_covariance_grid_sensitivity <- function(block_covariance, omega) {

  block_size <- nrow(block_covariance)
  bound <- (2 * block_size + 8) * .Machine$double.eps * max(abs(block_covariance))
  if (!(bound > 0)) {
    return(NA_real_)
  }
  variation <- sum(abs(diff(sort(omega))))
  pairs <- block_size * (block_size - 1) / 2
  determinant <- min(vapply(seq_len(block_size)[-1L], function(i) {
    min(vapply(seq_len(i - 1L), function(j) {
      block_covariance[i, i] * block_covariance[j, j] -
        abs(block_covariance[i, j])^2
    }, numeric(1)))
  }, numeric(1)))
  determinant <- determinant * (1 - 1e-8)
  if (!(determinant > 0)) {
    return(NA_real_)
  }
  variation * pairs * bound / (2 * pi * sqrt(determinant)) *
    max(max(omega), 1)^{max(block_size - 2L, 0L)}
}


# Per (block, posterior-row) binding for one call's state group.
.selection_covariance_grid_binding <- function(shared, rows, base, omega) {

  block_size <- length(rows)
  weights <- as.numeric(omega)
  maximum <- max(weights)
  minimum <- min(weights)
  if (!(maximum > 0) || any(!is.finite(weights)) || minimum < 0) {
    return(NULL)
  }
  # Equal weights remove selection, but the Gaussian likelihood still varies
  # with the swept mean or covariance. Let the exact Gaussian route handle it.
  if (maximum == minimum) {
    return(NULL)
  }
  if (identical(shared[["family"]], "variance")) {
    column <- shared[["loading"]][rows, , drop = FALSE]
    active <- which(colSums(column != 0) > 0)
    if (length(active) > 1L) {
      return(NULL)
    }
    m <- if (length(active) == 1L) column[, active[[1L]]] else rep(0, block_size)
    if (length(unique(m)) > 1L) {
      # A partially supported block spans several groups of the component and
      # the update basis is not rank one on it.
      return(NULL)
    }
    constant <- all(m == 0)
    D <- base[["covariance"]][rows, rows, drop = FALSE] -
      base[["coefficient"]] * (m %*% t(m))
    if (any(!is.finite(D)) || any(diag(D) <= 0)) {
      return(NULL)
    }
    d <- if (constant) 0 else .selection_covariance_grid_direction_norm(D, m)
    if (!constant && (!is.finite(d) || d <= 0)) {
      return(NULL)
    }
    sensitivity <- if (constant) 0 else
      .selection_covariance_grid_sensitivity(
        base[["covariance"]][rows, rows, drop = FALSE], weights)
    if (!is.finite(sensitivity)) {
      return(NULL)
    }
    return(list(
      constant = constant,
      d = d,
      log_scale = block_size * log(maximum),
      range = -expm1(block_size * log1p((minimum - maximum) / maximum)),
      sensitivity = sensitivity
    ))
  }

  # "mean" family: the fixed block covariance with the latent direction.
  S <- base[["covariance"]][rows, rows, drop = FALSE]
  b <- base[["direction"]][rows]
  constant <- all(b == 0)
  d <- if (constant) 0 else .selection_covariance_grid_direction_norm(S, b)
  if (!constant && (!is.finite(d) || d <= 0)) {
    return(NULL)
  }
  sensitivity <- if (constant) 0 else
    .selection_covariance_grid_sensitivity(S, weights)
  if (!is.finite(sensitivity)) {
    return(NULL)
  }
  list(
    constant = constant,
    d = d,
    log_scale = block_size * log(maximum),
    range = -expm1(block_size * log1p((minimum - maximum) / maximum)),
    sensitivity = sensitivity
  )
}


# Evaluate the dense-block log likelihood with the single-component family
# interpolation. Anchors are exact native rows of this call; every non-anchor
# point either passes its certified error gate or falls back to an exact row.
.selection_covariance_grid_loglik <- function(yi, means, covariance_lower, sei,
    selection_context, execution_plan, block_size, metadata) {

  shared <- metadata[["state"]][["shared"]]
  if (!is.environment(shared)) return(NULL)
  S <- nrow(means)
  modes <- rep(selection_context[["kernel_mode"]], length.out = S)
  rules <- rep(selection_context[["vector_rule"]], length.out = S)
  omega <- as.matrix(selection_context[["omega"]])
  if (all(modes == 0L)) return(NULL)
  if (any(modes != 1L) || any(rules != 0L) || any(!is.finite(omega)) ||
      any(omega <= 0)) {
    shared$untracked <- TRUE
    return(NULL)
  }
  state <- metadata[["state"]]
  state_index <- state[["state_index"]]
  values <- state[["values"]]
  qid <- match(values, shared$queries)
  if (length(qid) != S || anyNA(qid)) {
    shared$untracked <- TRUE
    return(NULL)
  }
  rows <- metadata[["rows"]]
  groups <- split(seq_len(S), state_index)
  entries <- vector("list", length(groups))
  requests <- integer()
  request_group <- integer()
  request_qid <- integer()
  request_kind <- character()
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    state_row <- state[["row_index"]][positions[[1L]]]
    index <- match(state_row, shared$rows)
    if (is.na(index) ||
        any(state[["row_index"]][positions] != state_row)) {
      shared$untracked <- TRUE
      return(NULL)
    }
    key <- paste(metadata[["block"]], index, sep = ":")
    weights <- as.numeric(omega[positions[[1L]], ])
    if (any(omega[positions, , drop = FALSE] !=
            matrix(weights, length(positions), length(weights), byrow = TRUE))) {
      shared$untracked <- TRUE
      return(NULL)
    }
    if (identical(shared[["family"]], "mean") &&
        any(covariance_lower[positions, , drop = FALSE] !=
            matrix(covariance_lower[positions[[1L]], ], length(positions),
                   ncol(covariance_lower)))) {
      # The mean family requires the covariance fixed across the sweep; the
      # generic builder's rows must agree exactly.
      shared$untracked <- TRUE
      return(NULL)
    }
    base <- .selection_covariance_grid_base(
      shared = shared,
      context = state[["context"]],
      active_setup = state[["active_setup"]],
      state_row = state_row
    )
    if (is.null(base)) {
      shared$untracked <- TRUE
      return(NULL)
    }
    binding <- .selection_covariance_grid_binding(
      shared = shared,
      rows = rows,
      base = base,
      omega = weights
    )
    if (is.null(binding)) {
      shared$untracked <- TRUE
      return(NULL)
    }
    if (exists(key, shared$cases, inherits = FALSE)) {
      entry <- get(key, shared$cases, inherits = FALSE)
      if (!identical(entry$binding, binding)) {
        shared$untracked <- TRUE
        return(NULL)
      }
    } else {
      entry <- new.env(parent = emptyenv())
      entry$binding <- binding
      entry$row <- index
      entry$constant <- isTRUE(binding$constant)
      entry$data <- matrix(NA_real_, length(shared$queries), 4L,
        dimnames = list(NULL, c("log_density", "log_A", "mcse", "eta")))
      entry$have <- logical(length(shared$queries))
      entry$seen <- logical(length(shared$queries))
      entry$reported <- numeric(length(shared$queries))
      entry$unknown <- logical(length(shared$queries))
      entry$leaves <- if (entry$constant || binding$range == 0) list() else
        .selection_normalizer_grid_geometry(shared, qid[positions])
      assign(key, entry, shared$cases)
    }
    entries[[g]] <- entry
    entry$seen[unique(qid[positions])] <- TRUE
    anchors <- if (entry$constant && !any(entry$have)) qid[positions[[1L]]] else
      unique(unlist(lapply(entry$leaves, `[[`, "nodes"), use.names = FALSE))
    anchors <- anchors[!entry$have[anchors]]
    if (length(anchors)) {
      at <- match(anchors, qid[positions])
      if (anyNA(at)) {
        shared$untracked <- TRUE
        return(NULL)
      }
      requests <- c(requests, positions[at])
      request_group <- c(request_group, rep(g, length(at)))
      request_qid <- c(request_qid, anchors)
      request_kind <- c(request_kind, rep(
        if (entry$constant) "constant" else "anchor", length(at)))
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
      value <- cbind(result$log_density[at], result$log_normalizer[at],
        result$relative_mcse[at], eta[at])
      rows_out <- if (entry$constant) seq_along(entry$have) else ids[at]
      if (entry$constant) value <- matrix(value[1L, ], length(rows_out), 4L, byrow = TRUE)
      entry$data[rows_out, ] <- value
      entry$have[rows_out] <- TRUE
    }
    for (name in c("anchor", "fallback", "constant")) {
      field <- paste0(name, "_evaluations")
      shared$stats[[field]] <- shared$stats[[field]] + sum(kind == name)
    }
    shared$stats[["max_native_mcse"]] <-
      max(shared$stats[["max_native_mcse"]], result$relative_mcse)
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
    binding <- entry$binding
    for (leaf in entry$leaves) {
      anchors <- entry$data[leaf$nodes, , drop = FALSE]
      A <- exp(anchors[, "log_A"] - binding$log_scale)
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
      conversion <- expm1(unit * (abs(anchors[, "log_A"]) + abs(binding$log_scale)) - log1p(-unit))
      anchor_error <- as.numeric(absolute %*% (A * (anchors[, "eta"] + conversion * (1 + anchors[, "eta"]))))
      m <- leaf$m
      # The variance family sweeps the covariance through the transport
      # A(tau) = E_z[a(tau z)], contributing the exact factor E|z|^m; the
      # mean family sweeps a deterministic shift A(tau) = a(tau) with the
      # plain mean-shift bound. E|z|^m < 1 for every m, so the factor must
      # not leak into the deterministic-shift remainder.
      log_Ez <- if (identical(shared[["family"]], "variance")) {
        (m / 2) * log(2) + lgamma((m + 1) / 2) - 0.5 * log(pi)
      } else {
        0
      }
      terms <- c(log(binding$range / 2), m * log(binding$d), log_Ez,
        -0.5 * lgamma(m + 1))
      log_remainder <- sum(terms) + leaf$log_product[within]
      log_remainder <- log_remainder + (5 * unit / (1 - 5 * unit)) *
        (sum(abs(terms)) + abs(leaf$log_product[within]))
      remainder <- exp(log_remainder)
      dot_gamma <- 2 * leaf$m * unit / (1 - 2 * leaf$m * unit)
      absolute_products <- sweep(abs(L), 2L, A, `*`)
      dot_safe <- rowSums((abs(L) > 0) &
        (!is.finite(absolute_products) | absolute_products < .Machine$double.xmin)) == 0
      rounding <- as.numeric(leaf$coefficient_error[within, , drop = FALSE] %*% A) +
        dot_gamma / (1 - dot_gamma) * rowSums(absolute_products)
      total <- remainder + anchor_error + rounding + binding$sensitivity
      accept <- dot_safe & remainder >= .Machine$double.xmin & is.finite(P) & P > 0 &
        is.finite(total) & total < P & total / P <= tolerance
      if (!any(accept)) next
      used <- positions[qid[positions] %in% ids[accept]]
      which_value <- match(qid[used], ids)
      log_A[used] <- log(P[which_value]) + binding$log_scale
      error[used] <- total[which_value] / P[which_value]
      interpolated[used] <- TRUE
      shared$stats[["max_anchor_relative_error"]] <-
        max(shared$stats[["max_anchor_relative_error"]],
            anchor_error[accept] / P[accept])
      shared$stats[["max_interpolation_relative_error"]] <- max(
        shared$stats[["max_interpolation_relative_error"]],
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
    interpolated_rows <- which(interpolated)
    gaussian_context <- BayesTools::selection_context_subset_rows(
      selection_context, interpolated_rows)
    gaussian_context[["kernel_mode"]] <- rep(0L, length(interpolated_rows))
    gaussian <- .selection_joint_dense_loglik_block(yi,
      means[interpolated_rows, , drop = FALSE],
      covariance_lower[interpolated_rows, , drop = FALSE], sei, gaussian_context,
      execution_plan, block_size)
    for (bin in selection_context[["obs_bin"]]) {
      gaussian <- gaussian + log(omega[interpolated_rows, bin])
    }
    result[interpolated_rows] <- gaussian - log_A[interpolated_rows]
    shared$stats[["gaussian_evaluations"]] <-
      shared$stats[["gaussian_evaluations"]] + length(interpolated_rows)
    shared$stats[["interpolated_points"]] <-
      shared$stats[["interpolated_points"]] + length(interpolated_rows)
  }
  if (any(is.finite(error))) shared$stats[["max_relative_error"]] <-
    max(shared$stats[["max_relative_error"]], error[is.finite(error)])
  for (g in seq_along(groups)) {
    positions <- groups[[g]]
    entry <- entries[[g]]
    at <- positions[!duplicated(qid[positions])]
    ids <- qid[at]
    values_out <- error[at]
    bad <- !is.finite(values_out) | values_out < 0 | values_out >= 1
    entry$unknown[ids[bad]] <- TRUE
    entry$reported[ids[!bad]] <- pmax(entry$reported[ids[!bad]], -log1p(-values_out[!bad]))
  }
  result
}


.selection_covariance_grid_diagnostics <- function(shared) {

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
    family = shared[["family"]],
    error_scope = paste("finite requested component-SD ordinates and their",
      "positive quadrature normalization; native errors remain estimated, with",
      "actual-node interpolation, ordinary floating allowances, and the",
      "covariance-builder accumulation allowance"),
    unique_requested_normalizers = unique_requests,
    max_log_likelihood_error = log_error,
    conditional_relative_error = if (is.finite(log_error)) expm1(2 * log_error) else NA_real_,
    unknown_error_points = sum(unknown & seen), untracked_path = shared$untracked),
    as.list(shared$stats))
}
