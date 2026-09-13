# Internal covariance factorization policy.

.covariance_factorization <- function(covariance, strict = FALSE) {

  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      nrow(covariance) != ncol(covariance) || nrow(covariance) == 0L ||
      anyNA(covariance) || any(!is.finite(covariance))) {
    stop("Internal error: covariance must be a finite non-empty numeric square matrix.",
         call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("Internal error: invalid strict covariance policy.", call. = FALSE)
  }

  if (any(covariance != t(covariance))) {
    stop("Covariance must be symmetric.", call. = FALSE)
  }

  size             <- nrow(covariance)
  diagonal         <- diag(covariance)
  zero_diagonal    <- which(diagonal == 0)
  invalid_variance <- any(diagonal < 0) ||
    any(covariance[zero_diagonal, , drop = FALSE] != 0)
  rank_one_factor  <- if (!invalid_variance && size > 1L) .covariance_exact_rank_one_factor(covariance) else NULL
  cholesky         <- NULL
  sampling_factor  <- NULL
  if (!is.null(rank_one_factor)) {
    sampling_factor <- matrix(rank_one_factor, nrow = 1L)
  }

  # Keep the raw values-only spectrum as a diagnostic. Rank and null-space
  # calculations below use the accepted factor, not this separate eigensolve.
  values <- eigen(covariance, symmetric = TRUE, only.values = TRUE)[["values"]]

  if (!invalid_variance && is.null(sampling_factor)) {
    positive <- which(diagonal > 0)
    if (!length(positive)) {
      sampling_factor <- matrix(0, nrow = 0L, ncol = size)
    } else {
      # Apply the existing numerical PSD policy in dimensionless correlation
      # coordinates. A small marginal variance is not a numerical null axis.
      sd          <- sqrt(diagonal[positive])
      correlation <- covariance[positive, positive, drop = FALSE] / tcrossprod(sd)
      diag(correlation) <- 1
      invalid_variance <- any(!is.finite(correlation))
      if (!invalid_variance) {
        decomposition <- eigen(correlation, symmetric = TRUE)
        # The vector solver can round a negative eigenvalue to zero, so retain
        # the values-only solver for the strict sign-classification contract.
        correlation_values <- eigen(correlation, symmetric = TRUE, only.values = TRUE)[["values"]]
        spectral_values <- decomposition[["values"]]
        invalid_variance <- any(!is.finite(correlation_values)) || any(!is.finite(spectral_values))
        if (!invalid_variance) {
          operation_count <- 4 * length(positive)
          correlation_tolerance <- operation_count * .Machine$double.eps /
            (1 - operation_count * .Machine$double.eps) * max(abs(correlation_values))
          spectral_values[abs(spectral_values) <= correlation_tolerance] <- 0
          invalid_variance <- min(correlation_values) <
            (if (strict) 0 else -correlation_tolerance) || any(spectral_values < 0)
        }
        if (!invalid_variance) {
          keep <- which(spectral_values > 0)
          # A rounded Cholesky pivot can be positive for an exactly
          # dependent covariance. Establish support in correlation units
          # first; Cholesky chooses a factor only when all axes are retained.
          if (length(keep) == size) {
            cholesky <- tryCatch(chol(covariance), error = function(e) NULL)
          }
          if (!is.null(cholesky)) {
            sampling_factor <- cholesky
          } else {
            sampling_factor <- matrix(0, nrow = length(keep), ncol = size)
            sampling_factor[, positive] <- sweep(sweep(
              t(decomposition[["vectors"]][, keep, drop = FALSE]),
              1L, sqrt(spectral_values[keep]), "*"), 2L, sd, "*")
          }
        }
      }
    }
  }

  status <- if (invalid_variance) "indefinite" else if (nrow(sampling_factor) == size) {
    "positive_definite"
  } else {
    "positive_semidefinite"
  }
  if (invalid_variance) {
    decomposition   <- eigen(covariance, symmetric = TRUE)
    spectral_values <- decomposition[["values"]]
    vectors         <- decomposition[["vectors"]]
  } else if (nrow(sampling_factor) == 0L) {
    spectral_values <- rep(0, size)
    vectors         <- diag(size)
  } else {
    # One accepted basis owns whitening, pseudoinverses and null-space tests.
    # The omitted rows of the compact factor are structural zeros; do not
    # rediscover them by thresholding a second full covariance eigensolve.
    decomposition   <- svd(sampling_factor, nu = 0L, nv = size)
    spectral_values <- c(decomposition[["d"]]^2, rep(0, size - nrow(sampling_factor)))
    vectors         <- decomposition[["v"]]
  }

  structure(
    list(
      covariance           = covariance,
      cholesky             = cholesky,
      sampling_factor      = sampling_factor,
      eigenvalues          = values,
      spectral_values      = spectral_values,
      eigenvectors         = vectors,
      singular             = !identical(status, "positive_definite"),
      status               = status
    ),
    class = c("brma_covariance_factorization", "list")
  )
}


# Return u only when the stored covariance is exactly u u' in binary64.
.covariance_exact_rank_one_factor <- function(covariance) {

  if (!is.matrix(covariance) || nrow(covariance) != ncol(covariance) ||
      anyNA(covariance) || any(!is.finite(covariance)) ||
      any(covariance != t(covariance)) || any(diag(covariance) <= 0)) {
    return(NULL)
  }

  magnitude <- sqrt(diag(covariance))
  pivot     <- which.max(magnitude)
  direction <- sign(covariance[pivot, ])
  if (any(direction == 0)) {
    return(NULL)
  }
  direction[[pivot]] <- 1
  factor             <- magnitude * direction

  if (!all(tcrossprod(factor) == covariance)) {
    return(NULL)
  }

  unname(factor)
}


.covariance_is_positive_semidefinite <- function(factorization) {

  inherits(factorization, "brma_covariance_factorization") &&
    factorization[["status"]] != "indefinite"
}


.covariance_is_numerically_positive_definite <- function(factorization) {

  inherits(factorization, "brma_covariance_factorization") &&
    identical(factorization[["status"]], "positive_definite")
}


.covariance_cholesky <- function(factorization) {

  if (!.covariance_is_numerically_positive_definite(factorization)) {
    return(NULL)
  }

  factorization[["cholesky"]]
}


# Return the accepted right factor, with K rows to preserve the RNG draw count.
.covariance_sampling_factor <- function(factorization) {

  if (!.covariance_is_positive_semidefinite(factorization)) {
    return(NULL)
  }

  factor <- factorization[["sampling_factor"]]
  size   <- ncol(factor)
  if (nrow(factor) == size) {
    return(factor)
  }
  sampling_factor <- matrix(0, nrow = size, ncol = size)
  sampling_factor[seq_len(nrow(factor)), ] <- factor
  sampling_factor
}
