# ============================================================================ #
# known-v-gls.R
# ============================================================================ #
#
# Shared known-V GLS covariance/projection helpers used by brma.mv diagnostics.
#
# ============================================================================ #


.known_v_extra_variance_samples <- function(object,
                                            posterior_samples = NULL,
                                            max_samples = Inf) {

  if (!inherits(object, "brma.mv") || !.is_data_known_v(object[["data"]])) {
    stop("Known-V extra variance samples require a brma.mv known-V object.",
         call. = FALSE)
  }

  known_V <- .data_known_v_data(object[["data"]])
  K       <- nrow(object[["data"]][["outcome"]])

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V covariance metadata is inconsistent with the outcome data.",
         call. = FALSE)
  }

  sample_info <- .known_v_diagnostic_posterior_samples(
    object            = object,
    posterior_samples = posterior_samples,
    max_samples       = max_samples,
    caller            = ".known_v_extra_variance_samples()",
    warn              = FALSE
  )
  posterior_samples <- sample_info[["posterior_samples"]]
  S                 <- nrow(posterior_samples)

  if (!.is_random(object)) {
    extra_variance <- .known_v_diagonal_extra_variance_samples(
      object            = object,
      posterior_samples = posterior_samples,
      K                 = K
    )
    attr(extra_variance, "known_v_diagnostic") <-
      .known_v_diagnostic_metadata(sample_info)
    return(extra_variance)
  }

  random_vcov <- .brma_mv_random_effects_marginal_vcov(
    object            = object,
    posterior_samples = posterior_samples,
    diagonal_only     = TRUE
  )
  random_variance <- .brma_mv_validate_random_marginal_variance_samples(
    random_vcov = random_vcov,
    S           = S,
    K           = K
  )
  random_variance <- unname(random_variance)

  metadata <- .known_v_diagnostic_metadata(sample_info)
  attr(random_variance, "known_v_diagnostic") <- metadata

  return(random_variance)
}


# ---------------------------------------------------------------------------- #
# .known_v_marginal_factor_plan
# ---------------------------------------------------------------------------- #
#
# Compile the exact marginal covariance representation used by known-V GLS
# diagnostics. Random-formula models use BayesTools' metadata-declared factor
# plans; diagonal heterogeneity models use the covariance plan's native
# per-observation variance input. Both representations retain the fitted
# sampling covariance without constructing draw x row x row arrays.
#
# ---------------------------------------------------------------------------- #
.known_v_marginal_factor_plan <- function(
    object, posterior_samples, known_V,
    extra_variances = NULL) {

  data <- object[["data"]]
  S    <- nrow(posterior_samples)
  K    <- nrow(data[["outcome"]])

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V covariance metadata is inconsistent with the outcome data.",
         call. = FALSE)
  }

  if (.is_data_random(data)) {
    joint_selection <- .is_data_joint_selection(data)
    sampled_blocks <- if (joint_selection) {
      NULL
    } else {
      .data_sampled_random_effect_blocks(data)
    }
    extra_variance <- if (is.null(extra_variances)) {
      if (joint_selection) {
        matrix(0, nrow = S, ncol = K)
      } else {
        .evaluate_marginalized_random_variance(
          data              = data,
          posterior_samples = posterior_samples,
          K                 = K,
          source_samples    = .known_v_marginalized_random_source_samples(
            fit               = object[["fit"]],
            data              = data,
            priors            = object[["priors"]],
            posterior_samples = posterior_samples
          ),
          sd_evaluator      = .marginalized_random_sd_evaluator(
            fit   = object[["fit"]],
            terms = .data_marginalized_random_effects(data)
          )
        )
      }
    } else {
      as.matrix(extra_variances)
    }

    if (joint_selection || length(sampled_blocks) > 0L) {
      random_factors <- .brma_mv_random_effects_marginal_factor_plan(
        object            = object,
        posterior_samples = posterior_samples,
        blocks            = sampled_blocks,
        row_blocks        = if (joint_selection) {
          .data_selection_execution_plan(data)[["row_blocks"]]
        } else {
          NULL
        },
        data              = data
      )
      dependency_blocks <- random_factors[["row_blocks"]]
      random_variance  <- BayesTools::random_effects_marginal_factor_diagonal(
        random_factors
      )
    } else {
      dependency_blocks <- .known_v_dependency_blocks(data, K)
      random_factors   <- list(
        factor_plans  = list(),
        factor_states = rep(list(list()), S)
      )
      random_variance <- matrix(0, nrow = S, ncol = K)
    }
  } else {
    dependency_blocks <- .known_v_dependency_blocks(data, K)
    random_factors   <- list(
      factor_plans  = list(),
      factor_states = rep(list(list()), S)
    )
    extra_variance <- if (is.null(extra_variances)) {
      .known_v_diagonal_extra_variance_samples(
        object            = object,
        posterior_samples = posterior_samples,
        K                 = K
      )
    } else {
      as.matrix(extra_variances)
    }
    random_variance <- matrix(0, nrow = S, ncol = K)
  }

  if (!identical(dim(extra_variance), c(S, K)) ||
      any(!is.finite(extra_variance)) || any(extra_variance < 0)) {
    stop("Known-V covariance factor plan received invalid extra variances.",
         call. = FALSE)
  }

  covariance_diagonal <- sweep(
    random_variance + extra_variance,
    2L,
    .known_v_diagonal(known_V),
    "+"
  )
  if (!identical(dim(covariance_diagonal), c(S, K)) ||
      any(!is.finite(covariance_diagonal)) ||
      any(covariance_diagonal < 0)) {
    stop("Known-V covariance factor plan returned an invalid diagonal.",
         call. = FALSE)
  }

  sampling_factors <- .known_v_sampling_factor_plan(
    known_V = known_V
  )
  latent_variances <- random_variance
  if (is.null(sampling_factors)) {
    sampling_covariance <- .known_v_covariance_matrix(known_V)
  } else {
    sampling_covariance <- sampling_factors[["sampling_covariance"]]
    latent_variances <- sweep(
      latent_variances, 2L,
      rowSums(sampling_factors[["factor_plan"]][["model_matrix"]]^2), "+"
    )
    random_factors[["factor_plans"]] <- c(
      random_factors[["factor_plans"]],
      list(sampling_factors[["factor_plan"]])
    )
    for (s in seq_len(S)) {
      random_factors[["factor_states"]][[s]] <- c(
        random_factors[["factor_states"]][[s]],
        list(sampling_factors[["factor_state"]])
      )
    }
  }

  return(list(
    sampling_covariance      = sampling_covariance,
    random_covariance_plans  = random_factors[["factor_plans"]],
    random_covariance_states = random_factors[["factor_states"]],
    block_indices            = dependency_blocks,
    extra_variances          = extra_variance,
    latent_variances         = latent_variances,
    covariance_diagonal      = covariance_diagonal
  ))
}


# Compile the exact fixed sampling factor independently of the fitting backend.
# Stored sampling factors and accepted whitening zeros retain their support
# independently of the post-fit native eigensolve.
.known_v_sampling_factor_plan <- function(known_V) {

  storage <- .known_v_storage(known_V)
  K       <- .known_v_nrow(known_V)
  sampling_covariance <- NULL
  if (identical(storage, "factor")) {
    model_matrix      <- known_V[["factor_loading"]]
    residual_variance <- known_V[["factor_diagonal"]]
  } else if (identical(.known_v_effective_backend(known_V), "latent")) {
    blocks <- .known_v_backend_blocks(known_V, "latent")
    rank   <- sum(vapply(blocks, `[[`, integer(1), "rank"))
    if (rank == 0L) {
      return(NULL)
    }

    model_matrix <- matrix(0, nrow = K, ncol = rank)
    offset       <- 0L
    for (block in blocks) {
      columns <- offset + seq_len(block[["rank"]])
      model_matrix[block[["index"]], columns] <- block[["B"]]
      offset <- offset + block[["rank"]]
    }
    residual_variance <- .known_v_residual_variance(known_V)
  } else if (identical(.known_v_effective_backend(known_V), "whitened")) {
    sampling_covariance <- .known_v_covariance_matrix(known_V)
    blocks <- .known_v_backend_blocks(known_V, "whitened")
    model_matrix <- matrix(0, nrow = K, ncol = 0L)
    seen <- integer()
    for (block in blocks) {
      index <- block[["index"]]
      variance <- block[["variance"]]
      rotation <- block[["rotation"]]
      if (!is.numeric(index) || !length(index) || anyNA(index) ||
          any(!is.finite(index)) || any(index != as.integer(index)) ||
          any(index < 1L | index > K) || anyDuplicated(index) ||
          any(index %in% seen) || !identical(as.integer(block[["size"]]), length(index)) ||
          !is.numeric(variance) || length(variance) != length(index) ||
          anyNA(variance) || any(!is.finite(variance)) || any(variance < 0) ||
          !is.matrix(rotation) || !is.numeric(rotation) ||
          !identical(dim(rotation), rep.int(length(index), 2L)) ||
          anyNA(rotation) || any(!is.finite(rotation))) {
        stop("Known-V whitening support metadata are invalid.", call. = FALSE)
      }
      seen <- c(seen, index)
      if (!any(variance == 0)) next
      positive <- which(variance > 0)
      loading <- sweep(t(rotation[positive, , drop = FALSE]),
        2L, sqrt(variance[positive]), "*")
      fixed <- matrix(0, nrow = K, ncol = length(positive))
      fixed[index, ] <- loading
      model_matrix <- cbind(model_matrix, fixed)
      sampling_covariance[index, index] <- 0
    }
    if (ncol(model_matrix) > 0L) {
      original <- .known_v_diagonal(known_V)
      reconstructed <- diag(sampling_covariance) + rowSums(model_matrix^2)
      if (any(abs(reconstructed - original) >
          8 * K * .Machine$double.eps * max(original))) {
        stop("Known-V whitening support no longer reproduces its covariance diagonal.",
          call. = FALSE)
      }
    }
    residual_variance <- diag(sampling_covariance)
  } else {
    return(NULL)
  }

  rank <- ncol(model_matrix)
  if (rank == 0L) {
    return(NULL)
  }
  if (length(residual_variance) != K || any(!is.finite(residual_variance)) ||
      any(residual_variance < 0)) {
    stop("Known-V factor residual variance metadata are invalid.",
         call. = FALSE)
  }

  return(list(
    sampling_covariance = if (is.null(sampling_covariance)) {
      diag(residual_variance, nrow = K, ncol = K)
    } else sampling_covariance,
    factor_plan = list(
      type                  = "group",
      model_matrix          = model_matrix,
      group_map             = rep(1L, K),
      coefficient_structure = "diagonal"
    ),
    factor_state = list(
      coefficient_factor = diag(1, nrow = rank, ncol = rank)
    )
  ))
}


.known_v_covariance_plan_precision_rhs_batch <- function(plan_data, rhs) {

  S <- nrow(plan_data[["extra_variances"]])
  K <- ncol(plan_data[["extra_variances"]])
  if (!is.matrix(rhs) || nrow(rhs) != K || ncol(rhs) == 0L ||
      !is.numeric(rhs) || any(!is.finite(rhs))) {
    stop("Known-V precision right-hand sides are invalid.", call. = FALSE)
  }

  zero_means <- matrix(0, nrow = S, ncol = K)
  out        <- array(NA_real_, dim = c(S, K, ncol(rhs)))
  for (column in seq_len(ncol(rhs))) {
    out[, , column] <- .marglik_covariance_plan_precision_residual_batch(
      cache                    = new.env(parent = emptyenv()),
      y                        = rhs[, column],
      means                    = zero_means,
      sampling_covariance      = plan_data[["sampling_covariance"]],
      random_covariance_plans  = plan_data[["random_covariance_plans"]],
      random_covariance_states = plan_data[["random_covariance_states"]],
      block_indices            = plan_data[["block_indices"]],
      extra_variances          = plan_data[["extra_variances"]]
    )
  }

  return(out)
}
.known_v_factor_gls_projection_batch <- function(
    object, posterior_samples, known_V, X, y,
    return_full_H = FALSE, return_se = FALSE, return_resid = FALSE,
    zero_candidates = integer()) {

  S <- nrow(posterior_samples)
  K <- length(y)
  p <- ncol(X)
  if (!is.matrix(X) || nrow(X) != K || p == 0L ||
      !is.numeric(X) || any(!is.finite(X)) ||
      !is.numeric(y) || any(!is.finite(y))) {
    stop("Known-V GLS projection inputs are invalid.", call. = FALSE)
  }

  plan_data <- .known_v_marginal_factor_plan(
    object            = object,
    posterior_samples = posterior_samples,
    known_V           = known_V
  )
  rhs <- if (return_resid) cbind(y, X) else X
  precision_rhs <- .known_v_covariance_plan_precision_rhs_batch(
    plan_data = plan_data,
    rhs       = rhs
  )
  x_offset <- if (return_resid) 1L else 0L

  H_diag_samples <- matrix(NA_real_, nrow = S, ncol = K)
  H_samples <- if (return_full_H) {
    array(NA_real_, dim = c(S, K, K))
  } else {
    NULL
  }
  residual_samples <- if (return_resid) {
    matrix(NA_real_, nrow = S, ncol = K)
  } else {
    NULL
  }
  residual_variance_samples <- if (return_se) {
    matrix(NA_real_, nrow = S, ncol = K)
  } else {
    NULL
  }

  for (s in seq_len(S)) {
    W_X <- matrix(
      precision_rhs[s, , x_offset + seq_len(p)],
      nrow = K,
      ncol = p
    )
    XtWX_inv <- .hat_solve_crossprod(
      crossprod(X, W_X),
      rank = attr(X, "rank")
    )
    X_cov    <- X %*% XtWX_inv

    H_diag_samples[s, ] <- rowSums(X_cov * W_X)
    if (return_full_H) {
      H_samples[s, , ] <- X_cov %*% t(W_X)
    }
    if (return_resid) {
      W_y <- precision_rhs[s, , 1L]
      beta_hat <- as.vector(XtWX_inv %*% crossprod(X, W_y))
      residual_samples[s, ] <- y - as.vector(X %*% beta_hat)
    }
    if (return_se) {
      # M_ii - x_i' G^-1 x_i subtracts nearly equal numbers once the model all
      # but interpolates a row. Such rows (below a millionth of M_ii) use the
      # factor form ||(e_i - H_i) F||^2 instead, which has no cancellation.
      diagonal_s <- plan_data[["covariance_diagonal"]][s, ]
      variance_s <- diagonal_s - rowSums(X_cov * X)
      tight      <- which(!(variance_s >= 1e-6 * diagonal_s))
      if (length(tight)) {
        # e_i - H_i for the tight rows i, without forming the K x K identity
        directions <- -t(X_cov[tight, , drop = FALSE] %*% t(W_X))
        unit <- cbind(tight, seq_along(tight))
        directions[unit] <- directions[unit] + 1
        variance_s[tight] <- .known_v_factor_quadratic(
          plan_data = plan_data,
          draw      = s,
          vectors   = directions
        )
      }
      residual_variance_samples[s, ] <- variance_s
    }
    zero_rows <- .known_v_isolated_projection_rows(plan_data, s, zero_candidates)
    if (return_resid) residual_samples[s, zero_rows] <- 0
    if (return_se) residual_variance_samples[s, zero_rows] <- 0
  }

  return(list(
    H_diag             = H_diag_samples,
    H                  = H_samples,
    M_diag             = plan_data[["covariance_diagonal"]],
    residual           = residual_samples,
    residual_variance  = residual_variance_samples,
    covariance_path    = "factor_plan_gls"
  ))
}


# A design certificate Q e_i = 0 also certifies a zero marginal residual row
# only when its covariance multiplier has no other partners. Inspect candidate
# covariance rows through the existing factors, without a new dense covariance.
.known_v_isolated_projection_rows <- function(plan_data, draw, candidates) {

  K <- nrow(plan_data[["sampling_covariance"]])
  if (length(candidates) == K) return(candidates) # wholly saturated design
  plans <- plan_data[["random_covariance_plans"]]
  states <- plan_data[["random_covariance_states"]][[draw]]
  candidates[vapply(candidates, function(i) {
    others <- setdiff(seq_len(K), i)
    if (any(plan_data[["sampling_covariance"]][i, others] != 0)) return(FALSE)
    for (factor in seq_along(plans)) {
      plan <- plans[[factor]]
      state <- states[[factor]]
      # Unrecognized representations retain the computed residual instead of
      # claiming a structural zero from incomplete covariance information.
      if (!plan[["type"]] %in% c("group", "known_group", "row_group") ||
          is.null(state[["coefficient_factor"]])) return(FALSE)
      design <- plan[["model_matrix"]]
      root <- state[["coefficient_factor"]]
      loading <- design[i, ]
      if (identical(plan[["type"]], "row_group")) loading <- loading * state[["row_scale"]][i]
      covariance_row <- as.vector(design %*% (root %*% crossprod(root, loading)))
      if (identical(plan[["type"]], "row_group")) covariance_row <- covariance_row * state[["row_scale"]]
      groups <- plan[["group_map"]]
      covariance_row <- covariance_row * if (identical(plan[["type"]], "known_group")) {
        plan[["group_covariance"]][groups, groups[i]]
      } else groups == groups[i]
      if (any(covariance_row[others] != 0)) return(FALSE)
    }
    TRUE
  }, logical(1))]
}


.known_v_diagnostic_posterior_samples <- function(object,
                                                  posterior_samples = NULL,
                                                  max_samples = Inf,
                                                  caller = "known-V diagnostic",
                                                  warn = TRUE) {

  max_samples       <- .normalize_max_samples(max_samples, "max_samples")
  posterior_samples <- .get_posterior_samples(object[["fit"]], posterior_samples)
  n_total           <- nrow(posterior_samples)
  selected_rows     <- .thin_sample_rows(n_total, max_samples)

  if (!is.null(selected_rows)) {
    posterior_samples <- posterior_samples[selected_rows, , drop = FALSE]
    if (warn) {
      warning(
        caller, " uses ", nrow(posterior_samples), " of ", n_total,
        " posterior draws because 'max_samples' was set; draws were ",
        "deterministically thinned across posterior row order.",
        call. = FALSE
      )
    }
  }

  return(list(
    posterior_samples = posterior_samples,
    n_total           = n_total,
    n_used            = nrow(posterior_samples),
    max_samples       = max_samples,
    thinned           = !is.null(selected_rows)
  ))
}


.known_v_covariance_max_bytes <- function(max_bytes = NULL) {

  if (is.null(max_bytes)) {
    max_bytes <- getOption(
      "RoBMA.known_v_covariance_max_bytes",
      4 * 1024^3
    )
  }

  valid <- is.numeric(max_bytes) &&
    length(max_bytes) == 1L &&
    !is.na(max_bytes) &&
    max_bytes > 0

  if (!valid) {
    stop(
      "'RoBMA.known_v_covariance_max_bytes' must be a single positive number ",
      "of bytes or Inf.",
      call. = FALSE
    )
  }

  return(as.double(max_bytes))
}


.known_v_covariance_bytes <- function(S, K) {

  return(8 * as.double(S) * as.double(K) * as.double(K))
}


.known_v_covariance_peak_bytes <- function(S, K) {

  bytes_per_matrix <- .known_v_covariance_bytes(1L, K)

  # Input and output draw arrays coexist while base covariance is added. Two
  # extra matrices conservatively cover the current block and R temporaries.
  (2 * as.double(S) + 2) * bytes_per_matrix
}


.known_v_format_bytes <- function(bytes) {

  units <- c("B", "KB", "MB", "GB", "TB")
  value <- as.double(bytes)
  unit  <- 1L
  while (is.finite(value) && value >= 1024 && unit < length(units)) {
    value <- value / 1024
    unit  <- unit + 1L
  }

  return(paste0(format(round(value, 2), trim = TRUE), " ", units[unit]))
}


.known_v_check_full_covariance_allocation <- function(S, K, max_bytes = NULL,
                                                      caller = "known-V covariance helper") {

  max_bytes <- .known_v_covariance_max_bytes(max_bytes)
  required  <- .known_v_covariance_peak_bytes(S, K)

  if (is.finite(max_bytes) && required > max_bytes) {
    stop(
      caller, " would require approximately ",
      .known_v_format_bytes(required), " at peak while constructing a full ",
      "draw x row x row covariance array, exceeding the configured budget of ",
      .known_v_format_bytes(max_bytes), ". Use a chunked diagnostic path, ",
      "reduce 'max_samples', or increase option ",
      "'RoBMA.known_v_covariance_max_bytes'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.known_v_covariance_chunk_indices <- function(S, K, max_bytes = NULL) {

  .known_v_covariance_cells_chunk_indices(S, as.double(K) * K, max_bytes)
}


.known_v_covariance_cells_chunk_indices <- function(S, cells_per_matrix, max_bytes = NULL) {

  max_bytes      <- .known_v_covariance_max_bytes(max_bytes)
  bytes_per_matrix <- 8 * as.double(cells_per_matrix)
  one_draw_bytes <- (2 * as.double(1L) + 2) * bytes_per_matrix

  if (is.finite(max_bytes) && one_draw_bytes > max_bytes) {
    stop(
      "A single known-V marginal covariance draw requires approximately ",
      .known_v_format_bytes(one_draw_bytes), ", exceeding the configured ",
      "budget of ", .known_v_format_bytes(max_bytes), ". Increase option ",
      "'RoBMA.known_v_covariance_max_bytes'.",
      call. = FALSE
    )
  }

  chunk_size <- if (is.infinite(max_bytes)) {
    S
  } else {
    max_draws         <- floor((max_bytes / bytes_per_matrix - 2) / 2)
    max(1L, min(S, max_draws))
  }
  starts <- seq.int(1L, S, by = chunk_size)

  return(lapply(starts, function(start) {
    seq.int(start, min(S, start + chunk_size - 1L))
  }))
}


.known_v_diagnostic_metadata <- function(sample_info, chunk_info = NULL) {

  metadata <- list(
    n_posterior_samples = sample_info[["n_total"]],
    n_used_samples      = sample_info[["n_used"]],
    max_samples         = sample_info[["max_samples"]],
    thinned             = sample_info[["thinned"]]
  )

  if (!is.null(chunk_info)) {
    metadata <- c(metadata, chunk_info)
  }

  return(metadata)
}


.known_v_attach_diagnostic_metadata <- function(x, metadata) {

  active <- isTRUE(metadata[["thinned"]]) ||
    (!is.null(metadata[["n_chunks"]]) && metadata[["n_chunks"]] > 1L)
  if (!active) {
    return(x)
  }

  attr(x, "known_v_diagnostic") <- metadata
  return(x)
}


.known_v_diagonal_extra_variance_samples <- function(object,
                                                      posterior_samples,
                                                      K) {

  S <- nrow(posterior_samples)

  if (.is_scale(object)) {
    tau_result <- .evaluate.brma.tau(
      fit               = object[["fit"]],
      scale_data        = object[["data"]][["scale"]],
      scale_formula     = .create_fit_formula_list(data = object[["data"]], "scale"),
      scale_priors      = object[["priors"]][["scale"]],
      is_scale          = TRUE,
      is_multilevel     = FALSE,
      K                 = K,
      posterior_samples = posterior_samples
    )
    extra_variance <- tau_result[["tau_within"]]^2
  } else if ("tau" %in% colnames(posterior_samples)) {
    extra_variance <- matrix(posterior_samples[, "tau"]^2,
                             nrow = S, ncol = K)
  } else {
    fixed_tau <- .fixed_tau_prior_value(object[["priors"]])
    if (is.null(fixed_tau)) {
      fixed_tau <- 0
    }
    extra_variance <- matrix(fixed_tau^2, nrow = S, ncol = K)
  }

  return(extra_variance)
}


# ---------------------------------------------------------------------------- #
# .known_v_factor_quadratic
# ---------------------------------------------------------------------------- #
#
# w' M w for the columns w of `vectors`, with M the marginal covariance of
# draw `draw` in the factor form M = V + diag(extra) + sum_f (B_f B_f') *
# [same group], B_f = model_matrix %*% coefficient_factor: the random-effect
# terms are sums of squares of group-aggregated loadings (a quadratic form in
# the known group covariance for known groups), so the result has no
# cancellation between a total and an explained part.
#
# ---------------------------------------------------------------------------- #
.known_v_factor_quadratic <- function(plan_data, draw, vectors) {

  value <- colSums(vectors * (plan_data[["sampling_covariance"]] %*% vectors)) +
    colSums(vectors^2 * plan_data[["extra_variances"]][draw, ])
  plans  <- plan_data[["random_covariance_plans"]]
  states <- plan_data[["random_covariance_states"]][[draw]]
  for (factor in seq_along(plans)) {
    plan  <- plans[[factor]]
    state <- states[[factor]]
    basis <- plan[["model_matrix"]] %*% state[["coefficient_factor"]]
    if (identical(plan[["type"]], "row_group")) {
      basis <- basis * state[["row_scale"]]
    }
    group_map <- as.integer(plan[["group_map"]])
    for (column in seq_len(ncol(vectors))) {
      aggregated <- rowsum(vectors[, column] * basis, group_map, reorder = TRUE)
      if (identical(plan[["type"]], "known_group")) {
        groups     <- as.integer(rownames(aggregated))
        covariance <- plan[["group_covariance"]][groups, groups, drop = FALSE]
        value[column] <- value[column] + sum(aggregated * (covariance %*% aggregated))
      } else {
        value[column] <- value[column] + sum(aggregated^2)
      }
    }
  }

  return(value)
}
