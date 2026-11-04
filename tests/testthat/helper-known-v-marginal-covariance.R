# Dense known-V covariance references used only by tests.

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


.known_v_certified_factor_blocks <- function(known_V) {

  if (identical(.known_v_storage(known_V), "factor")) {
    loading  <- known_V[["factor_loading"]]
    diagonal <- as.numeric(known_V[["factor_diagonal"]])
    return(lapply(known_V[["block_indices"]], function(index) {
      block_loading <- loading[index, , drop = FALSE]
      list(
        index    = index,
        diagonal = diagonal[index],
        loading  = unname(block_loading[, colSums(abs(block_loading)) > 0, drop = FALSE])
      )
    }))
  }

  factor   <- known_V[["certified_factor"]]
  diagonal <- as.numeric(factor[["diagonal"]])
  lapply(factor[["blocks"]], function(block) {
    list(
      index    = block[["index"]],
      diagonal = diagonal[block[["index"]]],
      loading  = block[["loading"]]
    )
  })
}
