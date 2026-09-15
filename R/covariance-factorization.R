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


# Exact diagonal-plus-block-constant recovery for one dependency block.
#
# Returns a certified `diag(D) + U U'` representation of the supplied block
# when its correlation matrix is block constant on nested partitions, and NULL
# otherwise. Acceptance is an algebraic identity checked to working precision,
# never a numerical-rank decision: the reconstruction must reproduce every
# supplied entry within the factorization roundoff envelope.
.covariance_block_constant_factor <- function(covariance, max_levels = 3L,
                                              separation = 1e-8) {

  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      nrow(covariance) != ncol(covariance) || nrow(covariance) < 2L ||
      anyNA(covariance) || any(!is.finite(covariance)) ||
      any(covariance != t(covariance))) {
    return(NULL)
  }

  size     <- nrow(covariance)
  variance <- unname(diag(covariance))
  if (any(variance <= 0)) {
    return(NULL)
  }

  sd          <- sqrt(variance)
  correlation <- unname(covariance / tcrossprod(sd))
  diag(correlation) <- 1
  if (anyNA(correlation) || any(!is.finite(correlation))) {
    return(NULL)
  }

  # The same roundoff convention as .covariance_factorization(), in
  # dimensionless correlation units: vcalc()-style entries reach working
  # precision along different rounding paths and differ by one or two ulps.
  operation_count <- 4 * size
  tolerance       <- operation_count * .Machine$double.eps /
    (1 - operation_count * .Machine$double.eps)

  off_diagonal <- correlation[lower.tri(correlation)]
  if (any(off_diagonal < -tolerance)) {
    return(NULL)
  }

  levels <- .covariance_correlation_levels(off_diagonal, tolerance)
  if (is.null(levels)) {
    return(NULL)
  }
  levels <- sort(levels[levels > tolerance], decreasing = TRUE)
  if (length(levels) == 0L || length(levels) > max_levels) {
    return(NULL)
  }
  if (length(levels) > 1L && any(-diff(levels) < separation)) {
    return(NULL)
  }
  # rho = 1 is the exactly singular structure owned by the rank-one detector.
  if (levels[[1L]] >= 1 - tolerance) {
    return(NULL)
  }

  partitions <- .covariance_nested_level_partitions(correlation, levels, tolerance)
  if (is.null(partitions)) {
    return(NULL)
  }

  increments  <- -diff(c(levels, 0))
  if (any(increments < 0)) {
    return(NULL)
  }
  diagonal    <- (1 - levels[[1L]]) * variance
  loadings    <- list()
  supports    <- list()
  parents     <- integer()
  for (level in seq_along(levels)) {
    partition <- partitions[[level]]
    for (group in unique(partition)) {
      rows <- which(partition == group)
      if (length(rows) == 1L) {
        # A singleton group carries no dependence; absorb it exactly.
        diagonal[rows] <- diagonal[rows] + increments[[level]] * variance[rows]
        next
      }
      column <- numeric(size)
      column[rows] <- sqrt(increments[[level]]) * sd[rows]
      loadings[[length(loadings) + 1L]] <- column
      supports[[length(supports) + 1L]] <- rows
      parents[[length(parents) + 1L]]   <- level
    }
  }
  if (length(loadings) == 0L) {
    return(NULL)
  }

  merged <- .covariance_merge_equal_supports(loadings, supports)
  loading  <- merged[["loading"]]
  supports <- merged[["supports"]]

  reconstruction <- tcrossprod(loading)
  diag(reconstruction) <- diag(reconstruction) + diagonal
  scale <- max(abs(covariance))
  if (!(scale > 0)) {
    return(NULL)
  }
  residual <- max(abs(reconstruction - covariance)) / scale
  if (!is.finite(residual) || residual > 8 * size * .Machine$double.eps) {
    return(NULL)
  }
  if (any(diagonal < 0) || anyNA(diagonal) || any(!is.finite(diagonal))) {
    return(NULL)
  }

  list(
    diagonal    = diagonal,
    loading     = loading,
    rank        = ncol(loading),
    levels      = levels,
    supports    = supports,
    support     = .covariance_support_shape(supports),
    residual    = residual,
    depth       = length(levels)
  )
}


# Cluster the observed off-diagonal correlations into distinct levels.
# Returns NULL when two candidate levels are separated by less than the
# roundoff tolerance yet more than one ulp apart within a run, which would
# make the level assignment ambiguous.
.covariance_correlation_levels <- function(values, tolerance) {

  values <- sort(unique(values))
  if (length(values) == 0L) {
    return(NULL)
  }

  level_values <- numeric(0)
  members      <- values[[1L]]
  for (value in values[-1L]) {
    if (value - members[[length(members)]] <= tolerance) {
      members <- c(members, value)
    } else {
      level_values <- c(level_values, mean(members))
      members      <- value
    }
  }
  level_values <- c(level_values, mean(members))

  level_values
}


# Build the nested partition induced by each correlation level, or NULL when
# a level is not an equivalence relation or the partitions are not nested.
.covariance_nested_level_partitions <- function(correlation, levels, tolerance) {

  size       <- nrow(correlation)
  partitions <- vector("list", length(levels))
  previous   <- NULL
  for (level in seq_along(levels)) {
    related <- correlation >= levels[[level]] - tolerance
    diag(related) <- TRUE
    # Disjoint cliques: the relation must be transitive as supplied, never
    # completed by a connected-component search.
    if (!all(((related %*% related) > 0) == related)) {
      return(NULL)
    }
    membership <- apply(related, 1L, function(row) which(row)[[1L]])
    if (!is.null(previous) &&
        any(tapply(membership, previous, function(x) length(unique(x))) != 1L)) {
      return(NULL)
    }
    partitions[[level]] <- membership
    previous            <- membership
  }

  partitions
}


# Merge columns with identical support so that a single-level block yields one
# column rather than parallel columns of the same group.
.covariance_merge_equal_supports <- function(loadings, supports) {

  keys   <- vapply(supports, function(rows) paste(rows, collapse = ","), character(1))
  unique_keys <- unique(keys)
  columns <- lapply(unique_keys, function(key) {
    selected <- which(keys == key)
    if (length(selected) == 1L) {
      return(loadings[[selected]])
    }
    sqrt(Reduce(`+`, lapply(loadings[selected], function(column) column^2)))
  })

  list(
    loading  = matrix(unlist(columns, use.names = FALSE),
                      nrow = length(loadings[[1L]]), ncol = length(columns)),
    supports = supports[match(unique_keys, keys)]
  )
}


# Classify the loading supports as a chain, a tree, or neither. The nested
# quadrature rules depend on this shape.
.covariance_support_shape <- function(supports) {

  if (length(supports) <= 1L) {
    return("chain")
  }

  order_by_size <- order(lengths(supports), decreasing = TRUE)
  supports      <- supports[order_by_size]
  is_chain <- all(vapply(seq_along(supports)[-1L], function(position) {
    all(supports[[position]] %in% supports[[position - 1L]])
  }, logical(1)))
  if (is_chain) {
    return("chain")
  }

  # A forest: every pair of supports is either disjoint or strictly nested.
  pairs_ok <- TRUE
  for (i in seq_along(supports)) {
    for (j in seq_along(supports)) {
      if (i >= j) next
      shared <- intersect(supports[[i]], supports[[j]])
      if (length(shared) == 0L) next
      if (!all(supports[[j]] %in% supports[[i]])) {
        pairs_ok <- FALSE
        break
      }
    }
    if (!pairs_ok) break
  }

  if (pairs_ok) "tree" else "general"
}
