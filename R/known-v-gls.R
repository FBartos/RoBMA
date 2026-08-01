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


.known_v_marginal_variance_samples <- function(object,
                                               posterior_samples = NULL,
                                               max_samples = Inf,
                                               max_bytes = NULL) {

  extra_variance <- .known_v_extra_variance_samples(
    object            = object,
    posterior_samples = posterior_samples,
    max_samples       = max_samples,
    max_bytes         = max_bytes
  )
  known_V <- .data_known_v_data(object[["data"]])

  marginal_variance <- sweep(
    extra_variance,
    2L,
    .known_v_diagonal(known_V),
    "+"
  )
  attr(marginal_variance, "known_v_diagnostic") <-
    attr(extra_variance, "known_v_diagnostic", exact = TRUE)

  return(marginal_variance)
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


.known_v_covariance_diagonal_samples <- function(covariance_samples) {

  dims <- dim(covariance_samples)
  if (length(dims) != 3L || dims[[2L]] != dims[[3L]]) {
    stop("Known-V covariance samples must be a draw x row x row array.",
         call. = FALSE)
  }

  S <- dims[[1L]]
  K <- dims[[2L]]
  diagonal_index <- cbind(
    rep(seq_len(S), each = K),
    rep(seq_len(K), times = S),
    rep(seq_len(K), times = S)
  )

  matrix(
    covariance_samples[diagonal_index],
    nrow  = S,
    ncol  = K,
    byrow = TRUE
  )
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


.known_v_apply_marginal_covariance_chunks <- function(object,
                                                      posterior_samples,
                                                      max_bytes = NULL,
                                                      FUN) {

  known_V <- .data_known_v_data(object[["data"]])
  K       <- nrow(object[["data"]][["outcome"]])

  if (.known_v_nrow(known_V) != K) {
    stop("Known-V covariance metadata is inconsistent with the outcome data.",
         call. = FALSE)
  }

  chunks    <- .known_v_covariance_chunk_indices(
    S         = nrow(posterior_samples),
    K         = K,
    max_bytes = max_bytes
  )
  max_bytes <- .known_v_covariance_max_bytes(max_bytes)

  for (rows in chunks) {
    covariance_samples <- .known_v_marginal_covariance_samples_raw(
      object            = object,
      posterior_samples = posterior_samples[rows, , drop = FALSE],
      known_V           = known_V,
      K                 = K
    )
    FUN(covariance_samples, rows)
  }

  return(list(
    max_bytes       = max_bytes,
    chunk_size      = max(lengths(chunks)),
    n_chunks        = length(chunks),
    full_bytes      = .known_v_covariance_bytes(nrow(posterior_samples), K),
    peak_bytes      = .known_v_covariance_peak_bytes(nrow(posterior_samples), K),
    one_draw_bytes  = .known_v_covariance_peak_bytes(1L, K)
  ))
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


# ---------------------------------------------------------------------------- #
# .known_v_gls_projection
# ---------------------------------------------------------------------------- #
#
# Shared GLS projection pieces for known-V residual diagnostics.
#
# ---------------------------------------------------------------------------- #
.known_v_gls_projection <- function(X, y, covariance) {

  factorization   <- .covariance_factorization(covariance)
  covariance      <- factorization[["covariance"]]
  chol_covariance <- .covariance_cholesky(factorization)
  if (is.null(chol_covariance)) {
    stop("Known-V residual covariance is not positive definite.",
         call. = FALSE)
  }

  W        <- chol2inv(chol_covariance)
  WX       <- W %*% X
  Wy       <- as.vector(W %*% y)
  XtWX     <- crossprod(X, WX)
  XtWX_inv <- .hat_solve_crossprod(XtWX)
  H        <- X %*% XtWX_inv %*% t(WX)
  beta_hat <- as.vector(XtWX_inv %*% crossprod(X, Wy))
  residual <- y - as.vector(X %*% beta_hat)

  return(list(
    covariance        = covariance,
    covariance_factor = t(chol_covariance),
    W                 = W,
    WX                = WX,
    XtWX_inv          = XtWX_inv,
    H                 = H,
    beta_hat          = beta_hat,
    residual          = residual
  ))
}
