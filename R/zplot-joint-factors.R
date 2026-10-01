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
      streamed_rows <- seq_len(S)
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
