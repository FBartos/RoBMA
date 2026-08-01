# Internal covariance factorization policy.

.covariance_factorization <- function(covariance) {

  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      nrow(covariance) != ncol(covariance) || nrow(covariance) == 0L ||
      anyNA(covariance) || any(!is.finite(covariance))) {
    stop("Internal error: covariance must be a finite non-empty numeric square matrix.",
         call. = FALSE)
  }

  covariance <- (covariance + t(covariance)) / 2
  eigen      <- eigen(covariance, symmetric = TRUE)
  values     <- eigen[["values"]]
  scale      <- max(abs(values))
  tolerance  <- if (scale == 0) {
    0
  } else {
    # Backward-error envelope for symmetrization and the eigensolve. This
    # classifies roundoff without changing the supplied covariance matrix.
    operation_count <- 4 * max(1, nrow(covariance))
    relative_error  <- operation_count * .Machine$double.eps /
      (1 - operation_count * .Machine$double.eps)
    relative_error * scale
  }

  status <- if (min(values) < -tolerance) {
    "indefinite"
  } else if (scale == 0 || min(values) <= tolerance) {
    "positive_semidefinite"
  } else {
    "positive_definite"
  }

  structure(
    list(
      covariance    = covariance,
      eigenvalues   = values,
      decomposition_values = eigen[["values"]],
      eigenvectors  = eigen[["vectors"]],
      scale         = scale,
      psd_tolerance = tolerance,
      pd_tolerance  = tolerance,
      singular      = !identical(status, "positive_definite"),
      status        = status
    ),
    class = c("brma_covariance_factorization", "list")
  )
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

  chol_factor <- tryCatch(
    chol(factorization[["covariance"]]),
    error = function(e) NULL
  )
  if (!is.null(chol_factor)) {
    return(chol_factor)
  }

  values <- pmax(factorization[["decomposition_values"]], 0)
  vectors <- factorization[["eigenvectors"]]
  vectors %*% diag(sqrt(values), nrow = length(values)) %*% t(vectors)
}
