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

  covariance      <- (covariance + t(covariance)) / 2
  rank_one_factor <- .covariance_exact_rank_one_factor(covariance)
  decomposition   <- eigen(covariance, symmetric = TRUE)
  values          <- eigen(
    covariance,
    symmetric   = TRUE,
    only.values = TRUE
  )[["values"]]
  scale           <- max(abs(values))
  tolerance       <- if (scale == 0) {
    0
  } else {
    # Backward-error envelope for symmetrization and the eigensolve. This
    # classifies roundoff without changing the supplied covariance matrix.
    operation_count <- 4 * max(1, nrow(covariance))
    relative_error  <- operation_count * .Machine$double.eps /
      (1 - operation_count * .Machine$double.eps)
    relative_error * scale
  }

  status <- if (!is.null(rank_one_factor)) {
    if (nrow(covariance) == 1L) "positive_definite" else "positive_semidefinite"
  } else if (min(values) < if (strict) 0 else -tolerance) {
    "indefinite"
  } else if (scale == 0 || min(values) <= tolerance) {
    "positive_semidefinite"
  } else {
    "positive_definite"
  }

  # A covariance accepted as positive semidefinite can have tiny signed
  # eigensolver artifacts in its null space. Preserve the submitted covariance
  # and normalize only the spectral representation used to construct factors.
  spectral_values <- decomposition[["values"]]
  if (!identical(status, "indefinite")) {
    spectral_values[abs(spectral_values) <= tolerance] <- 0
  }

  structure(
    list(
      covariance           = covariance,
      eigenvalues          = values,
      decomposition_values = decomposition[["values"]],
      spectral_values      = spectral_values,
      eigenvectors         = decomposition[["vectors"]],
      rank_one_factor      = rank_one_factor,
      scale                = scale,
      psd_tolerance        = tolerance,
      pd_tolerance         = tolerance,
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

  tryCatch(
    chol(factorization[["covariance"]]),
    error = function(e) NULL
  )
}


# Return a right-multiplication factor F with crossprod(F) equal to covariance.
.covariance_sampling_factor <- function(factorization) {

  if (!.covariance_is_positive_semidefinite(factorization)) {
    return(NULL)
  }

  rank_one_factor <- factorization[["rank_one_factor"]]
  if (!is.null(rank_one_factor)) {
    size            <- length(rank_one_factor)
    sampling_factor <- matrix(0, nrow = size, ncol = size)
    sampling_factor[1L, ] <- rank_one_factor
    return(sampling_factor)
  }

  chol_factor <- tryCatch(
    chol(factorization[["covariance"]]),
    error = function(e) NULL
  )
  if (!is.null(chol_factor)) {
    return(chol_factor)
  }

  values  <- factorization[["spectral_values"]]
  if (any(values < 0)) {
    return(NULL)
  }
  vectors <- factorization[["eigenvectors"]]
  vectors %*% diag(sqrt(values), nrow = length(values)) %*% t(vectors)
}
