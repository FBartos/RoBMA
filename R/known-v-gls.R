# ============================================================================ #
# known-v-gls.R
# ============================================================================ #
#
# Shared known-V GLS covariance/projection helpers used by brma.mv diagnostics.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .known_v_marginal_covariance_samples
# ---------------------------------------------------------------------------- #
#
# Posterior samples of the marginal row covariance for known-V brma.mv()
# diagnostics:
#   M_s = V + diag(tau_s^2)          for non-random models
#   M_s = V + Z G_s Z'               for random-formula models
#
# ---------------------------------------------------------------------------- #
.known_v_marginal_covariance_samples <- function(object,
                                                 posterior_samples = NULL,
                                                 max_samples = Inf,
                                                 max_bytes = NULL) {

  if (!inherits(object, "brma.mv") || !.is_data_known_v(object[["data"]])) {
    stop("Known-V marginal covariance samples require a brma.mv known-V object.",
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
    caller            = ".known_v_marginal_covariance_samples()",
    warn              = FALSE
  )
  posterior_samples <- sample_info[["posterior_samples"]]
  .known_v_check_full_covariance_allocation(
    S         = nrow(posterior_samples),
    K         = K,
    max_bytes = max_bytes,
    caller    = ".known_v_marginal_covariance_samples()"
  )

  covariance_samples <- .known_v_marginal_covariance_samples_raw(
    object            = object,
    posterior_samples = posterior_samples,
    known_V           = known_V,
    K                 = K
  )

  if (dim(covariance_samples)[1L] != nrow(posterior_samples)) {
    stop("Known-V marginal covariance draw count does not match posterior samples.",
         call. = FALSE)
  }

  return(covariance_samples)
}


.known_v_marginal_covariance_samples_raw <- function(object,
                                                     posterior_samples,
                                                     known_V,
                                                     K) {

  covariance_samples <- if (.is_random(object)) {
    random_vcov <- .brma_mv_random_effects_marginal_vcov(
      object            = object,
      posterior_samples = posterior_samples
    )
    random_vcov[["samples"]]
  } else {
    .known_v_diagonal_extra_covariance_samples(
      object            = object,
      posterior_samples = posterior_samples,
      K                 = K
    )
  }

  covariance_samples <- .known_v_add_base_covariance(
    known_V            = known_V,
    covariance_samples = covariance_samples
  )

  return(covariance_samples)
}


.known_v_extra_variance_samples <- function(object,
                                            posterior_samples = NULL,
                                            max_samples = Inf,
                                            max_bytes = NULL) {

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
    sampled_blocks <- .data_sampled_random_effect_blocks(data)
    extra_variance <- if (is.null(extra_variances)) {
      .evaluate_marginalized_random_variance(
        data              = data,
        posterior_samples = posterior_samples,
        K                 = K,
        source_samples    = .known_v_marginalized_random_source_samples(
          fit               = object[["fit"]],
          data              = data,
          priors            = object[["priors"]],
          posterior_samples = posterior_samples
        )
      )
    } else {
      as.matrix(extra_variances)
    }

    if (length(sampled_blocks) > 0L) {
      random_factors <- .brma_mv_random_effects_marginal_factor_plan(
        object            = object,
        posterior_samples = posterior_samples,
        blocks            = sampled_blocks,
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

  sampling_covariance <- .marglik_known_v_covariance_matrix(known_V)
  sampling_factors    <- .known_v_latent_sampling_factor_plan(
    known_V = known_V
  )
  if (!is.null(sampling_factors)) {
    sampling_covariance <- sampling_factors[["sampling_covariance"]]
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
    covariance_diagonal      = covariance_diagonal
  ))
}


.known_v_latent_sampling_factor_plan <- function(known_V) {

  if (!identical(.known_v_effective_backend(known_V), "latent")) {
    return(NULL)
  }
  blocks <- .known_v_backend_blocks(known_V, "latent")
  K      <- .known_v_nrow(known_V)
  rank   <- sum(vapply(blocks, `[[`, integer(1), "rank"))
  if (rank == 0L) {
    return(NULL)
  }

  model_matrix <- matrix(0, nrow = K, ncol = rank)
  offset       <- 0L
  for (block in blocks) {
    columns <- offset + seq_len(block[["rank"]])
    model_matrix[block[["index"]], columns] <- block[["B"]]
    offset <- max(columns)
  }

  residual_variance <- .known_v_residual_variance(known_V)
  if (length(residual_variance) != K || any(!is.finite(residual_variance)) ||
      any(residual_variance < 0)) {
    stop("Known-V latent residual variance metadata are invalid.",
         call. = FALSE)
  }

  return(list(
    sampling_covariance = diag(residual_variance, nrow = K, ncol = K),
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
  K <- nrow(rhs)
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
    return_full_H = FALSE, return_se = FALSE, return_resid = FALSE) {

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
      residual_variance_samples[s, ] <-
        plan_data[["covariance_diagonal"]][s, ] - rowSums(X_cov * X)
    }
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

  max_bytes      <- .known_v_covariance_max_bytes(max_bytes)
  one_draw_bytes <- .known_v_covariance_peak_bytes(1L, K)

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
    bytes_per_matrix <- .known_v_covariance_bytes(1L, K)
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


.known_v_diagonal_extra_covariance_samples <- function(object,
                                                       posterior_samples,
                                                       K) {

  extra_variance <- .known_v_diagonal_extra_variance_samples(
    object            = object,
    posterior_samples = posterior_samples,
    K                 = K
  )
  S <- nrow(extra_variance)

  covariance_samples <- array(0, dim = c(S, K, K))
  for (s in seq_len(S)) {
    covariance_samples[s, , ] <- diag(extra_variance[s, ],
                                      nrow = K, ncol = K)
  }

  return(covariance_samples)
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


.known_v_add_base_covariance <- function(known_V, covariance_samples) {

  .validate_known_v(known_V)
  K <- .known_v_nrow(known_V)

  if (length(dim(covariance_samples)) != 3L ||
      dim(covariance_samples)[2L] != K ||
      dim(covariance_samples)[3L] != K) {
    stop("Known-V covariance samples must have dimensions draw x row x row.",
         call. = FALSE)
  }

  out <- covariance_samples
  for (s in seq_len(dim(out)[1L])) {
    current <- matrix(out[s, , ], nrow = K, ncol = K)
    diag(current) <- diag(current) + .known_v_diagonal(known_V)
    for (block in .known_v_correlated_blocks(known_V)) {
      index     <- block[["index"]]
      increment <- block[["covariance"]]
      diag(increment) <- 0
      current[index, index] <- current[index, index] + increment
    }
    out[s, , ] <- current
  }

  return(out)
}


