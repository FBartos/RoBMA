# Internal covariance factorization policy.

.covariance_factorization <- function(covariance) {

  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      nrow(covariance) != ncol(covariance) || nrow(covariance) == 0L ||
      anyNA(covariance) || any(!is.finite(covariance))) {
    stop("Internal error: covariance must be a finite non-empty numeric square matrix.",
         call. = FALSE)
  }

  covariance     <- (covariance + t(covariance)) / 2
  eigen          <- eigen(covariance, symmetric = TRUE)
  # Classification historically used LAPACK's values-only path. Keep that
  # spectrum separate from the eigenvectors used for reconstruction.
  values         <- eigen(covariance, symmetric = TRUE, only.values = TRUE)[["values"]]
  scale          <- max(abs(values))
  psd_tolerance  <- sqrt(.Machine$double.eps) * max(1, scale)
  pd_tolerance   <- if (scale == 0) {
    0
  } else {
    .Machine$double.eps * max(1, nrow(covariance)) * scale
  }

  status <- if (min(values) < -psd_tolerance) {
    "indefinite"
  } else if (scale == 0 || min(values) <= pd_tolerance) {
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
      psd_tolerance = psd_tolerance,
      pd_tolerance  = pd_tolerance,
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
