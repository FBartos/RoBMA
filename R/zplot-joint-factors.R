# Prepare the same integrated and retained Gaussian sources without covariance
# cubes. NULL retains the existing route for structurally unsupported targets.
.zplot_joint_factor_preparation <- function(object, posterior_samples, selection,
                                           probability, control) {

  data <- object[["data"]]
  if (!inherits(object, "brma.mv") || is.null(selection) || probability ||
      .selection_retains_sampling(data)) return(NULL)
  S_all <- nrow(posterior_samples)
  active <- which(!rep_len(selection[["use_normal"]], S_all))
  if (!length(active)) return(list(active = active))
  if (any(rep_len(selection[["vector_rule"]], S_all)[active] != 0L) ||
      any(rep_len(selection[["kernel_mode"]], S_all)[active] != SELKERNEL_STEP) || !.is_data_random(data)) return(NULL)
  model <- .data_selection_model(data)
  plan <- .selection_joint_execution_plan_with_control(.data_selection_execution_plan(data), control)
  K <- nrow(data[["outcome"]])
  blocks <- plan[["row_blocks"]]
  if (!identical(sort(as.integer(unlist(blocks, use.names = FALSE))), seq_len(K))) {
    stop("Zplot execution blocks do not partition the fitted observations.", call. = FALSE)
  }
  packed_bytes <- 8 * length(active) * max(lengths(blocks) * (lengths(blocks) + 1) / 2)
  if (packed_bytes > .known_v_covariance_max_bytes()) return(NULL)
  # Positive-definite sampling blocks guarantee positive-definite candidate
  # covariance after any integrated Gaussian random contribution. Singular
  # sampling structures keep their existing projection/factor route.
  for (rows in blocks) {
    sampling <- .selection_joint_sampling_block(plan[["sampling"]], rows)
    if (is.null(tryCatch(chol(sampling), error = function(e) NULL))) return(NULL)
  }
  sources <- model[["sources"]][["random"]]
  source_names <- vapply(sources, `[[`, character(1L), "name")
  retained_source <- vapply(sources, function(source) {
    flag <- source[["retained"]]
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) {
      stop("Zplot source retention metadata are invalid.", call. = FALSE)
    }
    flag
  }, logical(1L))
  if (anyNA(source_names) || any(!nzchar(source_names)) || anyDuplicated(source_names)) {
    stop("Zplot random-source names are invalid.", call. = FALSE)
  }
  retained_names <- source_names[retained_source]
  if (!length(retained_names)) return(NULL)
  samples <- posterior_samples[active, , drop = FALSE]
  setup <- list(fit = object[["fit"]], data = data, priors = object[["priors"]],
    posterior_samples = samples, S = nrow(samples), K = K, is_multilevel = .is_multilevel(object))
  inputs <- .brma_mv_random_effects_marginal_inputs(object, samples)
  dependency <- BayesTools::random_effects_dependency_matrix(
    random_effects = inputs[["formula_design"]][["random_effects"]], n_rows = K,
    blocks = retained_names)
  block_id <- integer(K)
  for (block in seq_along(blocks)) block_id[blocks[[block]]] <- block
  if (any(dependency & outer(block_id, block_id, `!=`))) return(NULL)
  retained_states <- .brma_mv_random_effects_marginal_factor_states(object, samples,
    blocks = retained_names, row_blocks = blocks, inputs = inputs)
  if (!identical(retained_states[["row_blocks"]], blocks) ||
      !identical(retained_states[["metadata"]][["included_blocks"]], retained_names) ||
      length(retained_states[["factor_states"]]) != nrow(samples)) {
    stop("Retained zplot factor states disagree with the fitted source metadata.", call. = FALSE)
  }
  # Avoid materializing unsupported retained loading arrays. The context
  # projection admits one compiled scalar coefficient and one group per block;
  # public reduction metadata below still makes the final eligibility decision.
  plans <- retained_states[["factor_plans"]]
  if (length(plans) != 1L || !is.matrix(plans[[1L]][["model_matrix"]]) ||
      ncol(plans[[1L]][["model_matrix"]]) != 1L ||
      any(vapply(blocks, function(rows) length(unique(plans[[1L]][["group_map"]][rows])) > 1L,
        logical(1L)))) return(NULL)
  retained <- tryCatch(BayesTools::random_effects_marginal_diagonal_factor(retained_states),
    BayesTools_random_effects_marginal_factor_unavailable = function(e) NULL)
  if (is.null(retained)) return(NULL)
  S <- nrow(samples)
  support <- retained[["diagonal_support"]]
  if (!is.logical(support) || length(support) != K || anyNA(support) ||
      !identical(dim(retained[["diagonal"]]), c(S, K)) ||
      any(!is.finite(retained[["diagonal"]])) || any(retained[["diagonal"]] < 0) ||
      any(retained[["diagonal"]][, !support, drop = FALSE] != 0) ||
      !identical(retained[["row_blocks"]], blocks) ||
      length(retained[["loadings"]]) != length(blocks) ||
      !is.integer(retained[["ranks"]]) || length(retained[["ranks"]]) != length(blocks)) {
    stop("Retained zplot factors or structural diagonal supports are invalid.", call. = FALSE)
  }
  for (block in seq_along(blocks)) {
    rows <- blocks[[block]]
    rank <- retained[["ranks"]][[block]]
    loading <- retained[["loadings"]][[block]]
    loading_support <- retained[["loading_supports"]][[block]]
    if (is.na(rank) || rank < 0L || !identical(dim(loading), c(S, length(rows), rank)) ||
        any(!is.finite(loading)) || !is.logical(loading_support) || anyNA(loading_support) ||
        !identical(dim(loading_support), c(length(rows), rank))) {
      stop("Retained zplot loading metadata are invalid.", call. = FALSE)
    }
    if (length(rows) > 1L && (any(support[rows]) || rank > 1L)) return(NULL)
  }
  integrated <- .selection_joint_random_factor_samples(setup, inputs = inputs)
  if (is.null(integrated) && !is.null(plan[["random_covariance"]])) return(NULL)
  list(active = active, setup = setup, execution_plan = plan, integrated = integrated,
    retained = retained, selection = BayesTools::selection_context_subset_rows(selection, active))
}



# Evaluate one dependency block's fitted zplot density for every posterior row
# at once. The integrated Gaussian law reaches the batched, threaded kernel
# through the block's certified factor; the retained context is one standard
# normal axis, integrated on the same Gauss-Hermite ladder the rest of the
# selection code walks. Rows leave the ladder as they converge, so a slow row
# does not force the whole block to a finer rule.
#
# Returns NULL whenever the structure or the diagnostics do not support this
# route, and the caller keeps its per-draw projection.
.zplot_joint_factor_block_batch <- function(
    z, block_mean, covariance_lower, block_sei, context_loading, context_rank,
    block_context, setup, plan, block, integrated, control, designs, threads) {

  S <- nrow(block_mean)
  k <- ncol(block_mean)
  # The batch trades more quadrature nodes for work the threaded kernel can
  # spread over rows. With one thread that trade loses to the per-draw
  # projection, which is serial but evaluates a compressed geometry.
  if (threads < 2L) {
    return(NULL)
  }
  if (k < 2L || context_rank > 1L || is.null(plan[["sampling_factor_blocks"]]) ||
      !plan[["block_methods"]][[block]] %in% c("rank_one", "factor")) {
    return(NULL)
  }
  factors <- tryCatch(
    .selection_joint_factor_block_samples(setup, block, integrated),
    error = function(e) NULL)
  if (is.null(factors) || !identical(dim(factors[["residual_sd"]]), c(S, k)) ||
      factors[["rank"]] < 1L || factors[["rank"]] > k ||
      any(!is.finite(factors[["residual_sd"]])) ||
      any(factors[["residual_sd"]] <= 0) ||
      !.selection_factor_support_is_forest(factors[["loading_support"]])) {
    return(NULL)
  }
  context <- if (context_rank == 1L) matrix(context_loading, S, k) else
    matrix(0, S, k)
  if (any(!is.finite(context))) {
    return(NULL)
  }

  # A bound on the block's selected density at any context realization: the
  # largest Gaussian ordinate of the integrated law, inflated by the weight
  # ratio. It bounds the mass this ladder leaves outside its extreme nodes.
  omega <- block_context[["omega"]]
  if (!is.matrix(omega) || any(!is.finite(omega)) || any(omega <= 0)) {
    return(NULL)
  }
  density_bound <- numeric(S)
  covariance <- matrix(0, k, k)
  lower_index <- lower.tri(covariance, diag = TRUE)
  upper_index <- upper.tri(covariance)
  for (draw in seq_len(S)) {
    covariance[lower_index] <- covariance_lower[draw, ]
    covariance[upper_index] <- t(covariance)[upper_index]
    root <- tryCatch(chol(covariance), error = function(e) NULL)
    if (is.null(root)) {
      return(NULL)
    }
    density_bound[[draw]] <- max(omega[draw, ]) / min(omega[draw, ]) *
      mean(block_sei * sqrt(diag(chol2inv(root)))) / sqrt(2 * pi)
  }
  if (any(!is.finite(density_bound))) {
    return(NULL)
  }

  # The first rules of the ladder the per-draw projection walks for this axis.
  # Beyond these the outer nodes sit far enough into the tail that the inner
  # kernel leaves its own quadrature rules for quasi-Monte Carlo, which costs
  # far more than the per-draw projection those rows fall back to.
  orders    <- c(3L, 7L, 15L)
  accepted  <- matrix(0, S, length(z))
  pending   <- seq_len(S)
  previous  <- NULL
  previous_error <- NULL
  for (order in orders) {
    if (!length(pending)) {
      break
    }
    rule  <- .gauss_hermite_nodes(order)
    nodes <- length(rule[["nodes"]])
    rows  <- length(pending)
    if (!is.finite(rows * nodes * length(z) * 8) ||
        rows * nodes * length(z) * 8 > .known_v_covariance_max_bytes()) {
      return(NULL)
    }
    row_index  <- rep.int(pending, nodes)
    node_index <- rep(seq_len(nodes), each = rows)
    projected <- tryCatch(.zplot_joint_block(
      z = z,
      mean = block_mean[row_index, , drop = FALSE] +
        context[row_index, , drop = FALSE] * rule[["nodes"]][node_index],
      covariance_lower = covariance_lower[row_index, , drop = FALSE],
      sei = block_sei,
      selection = BayesTools::selection_context_subset_rows(block_context, row_index),
      probability = FALSE, control = control, designs = designs,
      factors = list(
        residual_sd = factors[["residual_sd"]][row_index, , drop = FALSE],
        loading     = factors[["loading"]][row_index, , drop = FALSE],
        loading_support = factors[["loading_support"]]
      )), error = function(e) NULL)
    if (is.null(projected) || any(!is.finite(projected[["density"]])) ||
        any(projected[["density"]] < 0)) {
      return(NULL)
    }

    weights <- exp(rule[["log_weights"]])
    current <- matrix(0, rows, length(z))
    inner   <- numeric(rows)
    for (node in seq_len(nodes)) {
      index  <- (node - 1L) * rows + seq_len(rows)
      values <- projected[["density"]][index, , drop = FALSE]
      current <- current + weights[[node]] * values
      inner   <- inner + weights[[node]] * projected[["relative_mcse"]][index] *
        apply(values, 1L, max)
    }
    # Gaussian mass outside the extreme nodes, bounded by the density bound.
    tail  <- density_bound[pending] * 2 *
      stats::pnorm(max(rule[["nodes"]]), lower.tail = FALSE)
    peak  <- apply(current, 1L, max)
    error <- if (is.null(previous)) {
      rep(Inf, rows)
    } else {
      apply(abs(current - previous), 1L, max) + inner + previous_error
    }
    error    <- error + tail
    relative <- ifelse(peak - inner > 0, error / (peak - inner), Inf)
    settled  <- is.finite(relative) & relative <= control[["relative_tolerance"]]
    if (any(settled)) {
      accepted[pending[settled], ] <- current[settled, , drop = FALSE]
    }
    pending        <- pending[!settled]
    previous       <- current[!settled, , drop = FALSE]
    previous_error <- inner[!settled]
  }

  list(density = accepted, done = !seq_len(S) %in% pending)
}


.zplot_joint_factor_marginal <- function(object, posterior_samples, predictive,
    selection, z, probability, control) {

  prepared <- .zplot_joint_factor_preparation(object, posterior_samples, selection, probability, control)
  if (is.null(prepared)) return(NULL)
  gaussian <- .zplot_gaussian_marginal_reference(object, posterior_samples, predictive)
  sd <- sqrt(gaussian[["variance"]])
  reference <- .zplot_normal_density_matrix(z, gaussian[["mu"]], sd, gaussian[["sei"]])
  fitted <- if (identical(predictive[["mu"]], gaussian[["mu"]])) reference else
    .zplot_normal_density_matrix(z, predictive[["mu"]], sd, gaussian[["sei"]])
  active <- prepared[["active"]]
  if (length(active)) {
    setup <- prepared[["setup"]]
    plan <- prepared[["execution_plan"]]
    retained <- prepared[["retained"]]
    means <- predictive[["mu"]][active, , drop = FALSE]
    S <- nrow(means)
    K <- ncol(means)
    selected <- matrix(0, S, length(z))
    designs <- new.env(parent = emptyenv())
    threads <- .resolve_native_threads(object)
    publication <- .data_selection_model(object[["data"]])[["groups"]][["group_index"]]
    for (block in seq_along(plan[["row_blocks"]])) {
      observations <- plan[["row_blocks"]][[block]]
      k <- length(observations)
      lower <- .selection_joint_covariance_lower(setup, block,
        random_factor_samples = prepared[["integrated"]])
      block_context <- BayesTools::selection_context_subset_observations(prepared[["selection"]], observations)
      static <- BayesTools::selection_native_static_args(block_context)
      fallback_factors <- NULL
      fallback_prepared <- FALSE
      loading <- retained[["loadings"]][[block]]
      rank <- retained[["ranks"]][[block]]
      if (k == 1L) {
        variance <- retained[["diagonal"]][, observations] + rowSums(matrix(loading, S)^2)
        scalar <- .zplot_latent_mixture(z, means[, observations, drop = FALSE],
          sqrt(lower), matrix(sqrt(variance), S, 1L), predictive[["sei"]][observations],
          block_context, FALSE, control, fitted_only = TRUE)
        selected <- selected + scalar[["fitted"]] / K
        next
      }
      pairs <- .selection_joint_lower_pairs(plan, seq_len(k))
      lower_index <- pairs[["row_1"]] + (pairs[["row_2"]] - 1L) * k
      upper_index <- pairs[["row_2"]] + (pairs[["row_1"]] - 1L) * k
      block_mean <- means[, observations, drop = FALSE]
      block_sei <- predictive[["sei"]][observations]
      # A certified sampling factor lets the whole block reach the batched,
      # threaded kernel for every posterior row at once, with the retained
      # context as one outer quadrature axis, instead of rebuilding the
      # per-draw projection geometry and rule tables for each row in turn.
      batched <- .zplot_joint_factor_block_batch(
        z = z, block_mean = block_mean, covariance_lower = lower,
        block_sei = block_sei, context_loading = loading, context_rank = rank,
        block_context = block_context, setup = setup, plan = plan,
        block = block, integrated = prepared[["integrated"]],
        control = control, designs = designs, threads = threads)
      streamed_rows <- seq_len(S)
      if (!is.null(batched)) {
        done <- batched[["done"]]
        selected[done, ] <- selected[done, , drop = FALSE] +
          batched[["density"]][done, , drop = FALSE] * k / K
        streamed_rows <- which(!done)
        if (!length(streamed_rows)) {
          next
        }
      }
      # Preparation has validated every posterior row and admitted only active
      # STEP/product selection. The projection consumes these four row fields
      # and the immutable native specification; full context reconstruction is
      # needed only by its generic fallback.
      projection_context <- list(use_normal = FALSE,
        kernel_mode = SELKERNEL_STEP, vector_rule = 0L,
        omega = NULL, native_cache = block_context[["native_cache"]])
      for (draw in streamed_rows) {
        lower_draw <- lower[draw, ]
        covariance <- matrix(0, k, k)
        covariance[lower_index] <- lower_draw
        covariance[upper_index] <- lower_draw
        context_factor <- if (rank > 0L) t(matrix(loading[draw, , , drop = FALSE], k, rank)) else
          matrix(0, 1L, k)
        projection_context[["omega"]] <- block_context[["omega"]][draw, , drop = FALSE]
        projected <- .zplot_context_projection(z, as.numeric(block_mean[draw, ]),
          covariance, context_factor, block_sei, projection_context, control)
        if (is.null(projected)) {
          context <- BayesTools::selection_context_subset_rows(block_context, draw)
          assign("static", static, envir = context[["native_cache"]])
          if (!fallback_prepared) {
            if (!is.null(plan[["sampling_factor_blocks"]]) &&
                plan[["block_methods"]][[block]] %in% c("rank_one", "factor")) {
              fallback_factors <- .selection_joint_factor_block_samples(setup, block, prepared[["integrated"]])
            }
            fallback_prepared <- TRUE
          }
          local_factors <- fallback_factors
          if (!is.null(local_factors)) {
            for (name in c("residual_sd", "loading")) {
              if (!is.null(local_factors[[name]])) local_factors[[name]] <- local_factors[[name]][draw, , drop = FALSE]
            }
          }
          # Non-scalar streamed contexts have a structurally absent diagonal.
          # The generic fallback receives their complete retained covariance.
          density <- .zplot_full_event_context_draw(z, means[draw, observations], covariance,
            crossprod(context_factor), predictive[["sei"]][observations], context, FALSE, control, plan,
            factors = local_factors, publication_groups = publication[observations], designs = designs)
        } else {
          density <- .zplot_context_projection_density(projected, z, control)
        }
        if (!identical(dim(density), c(1L, length(z))) ||
            any(!is.finite(density)) || any(density < 0)) {
          stop("Zplot context projection densities are unavailable.", call. = FALSE)
        }
        selected[draw, ] <- selected[draw, ] + density[1L, ] * k / K
      }
    }
    fitted[active, ] <- selected
  }
  list(fitted = fitted, extrapolated = reference, weights = rep(1, nrow(posterior_samples)), EDR = NULL)
}
