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

  # C = (1 - rho_1) I + sum_l (rho_l - rho_{l+1}) sum_{g in P_l} 1_g 1_g',
  # with rho_{L+1} = 0. Every increment is positive because the levels are
  # distinct, ordered and strictly above the roundoff floor.
  increments <- -diff(c(levels, 0))
  diagonal   <- (1 - levels[[1L]]) * variance
  loadings   <- list()
  supports   <- list()
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
    method      = "block_constant",
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


# Exact minimum-rank recovery for one dependency block.
#
# Covariance from overlapping samples, shared control arms or partially shared
# raters is exactly `diag(d) + U U'` with a small number of columns without
# being block constant on nested partitions, which is the structure
# .covariance_block_constant_factor() recovers. This searches for the smallest
# rank that reproduces every supplied entry within the same factorization
# roundoff envelope. Acceptance remains an algebraic identity: an approximate
# fit at any rank is declined and the block keeps its supplied entries and its
# general route. No approximation is ever accepted.
#
# The search stops at the Ledermann bound as well as at the kernel rank cap.
# Above that bound the off-diagonal entries no longer outnumber the free
# loadings, every positive-definite block fits exactly -- `diag(lambda_min) +
# U U'` at rank `n - 1` reproduces any of them -- and an exact fit stops being
# evidence about the supplied covariance. Recovery is a statement about that
# covariance, so it stops where the statement stops being one.
.covariance_minimum_rank_factor <- function(covariance,
                                            max_rank = SELNORM_FACTOR_MAX_RANK) {

  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      nrow(covariance) != ncol(covariance) || nrow(covariance) < 2L ||
      anyNA(covariance) || any(!is.finite(covariance)) ||
      any(covariance != t(covariance))) {
    return(NULL)
  }

  covariance <- unname(covariance)
  size       <- nrow(covariance)
  variance   <- diag(covariance)
  scale      <- max(abs(covariance))
  if (any(variance <= 0) || !(scale > 0)) {
    return(NULL)
  }
  # The same convention as .covariance_factorization(): entries that reach
  # working precision along different rounding paths differ by a few ulps.
  tolerance <- 8 * size * .Machine$double.eps

  # Every accepted representation keeps each residual variance above the
  # certificate floor, so `D + U U'` is positive definite with its smallest
  # eigenvalue at least that floor. A block at or below it -- exactly singular,
  # like a rank-one `s s'` with no residual spread, or rank deficient to
  # working precision -- has no acceptable representation at any rank, and
  # the ladder would only refine towards representations the certificate
  # declines.
  smallest <- min(eigen(covariance, symmetric = TRUE,
                        only.values = TRUE)[["values"]])
  if (!is.finite(smallest) || smallest <= tolerance * max(variance)) {
    return(NULL)
  }

  # Rank one is always in scope: for two rows it is the only representation
  # there is, and it is the one the cluster route wants.
  limit <- max(min(
    as.integer(max_rank), size - 1L,
    max(.covariance_ledermann_rank(size), 1L)
  ), 0L)
  # Ranks below the off-diagonal bound cannot hold an exact representation, and
  # a bound above the limit settles the block without any refinement at all.
  lower <- .covariance_off_diagonal_rank_bound(covariance, tolerance)
  if (lower > limit) {
    return(NULL)
  }
  first <- max(lower, 1L)
  ranks <- if (limit >= first) seq.int(first, limit) else integer(0)
  for (rank in ranks) {
    accepted <- Filter(Negate(is.null), lapply(
      .covariance_low_rank_starts(covariance, rank),
      function(loading) {
        .covariance_low_rank_accept(
          covariance, .covariance_low_rank_polish(covariance, loading),
          scale, tolerance
        )
      }
    ))
    if (length(accepted) > 0L) {
      # Several exact representations of the same rank can differ in how much
      # residual spread they leave on a row, and a route conditions on that
      # spread. Prefer the best-conditioned one rather than the first found.
      margin <- vapply(accepted, function(candidate) {
        min(candidate[["diagonal"]] / diag(covariance))
      }, numeric(1))
      return(accepted[[which.max(margin)]])
    }
  }

  NULL
}


# Smallest rank any exact `D + U U'` representation of the block could have.
# Splitting the rows into two disjoint groups makes the off-diagonal submatrix
# `V[I, J]` equal `U[I, ] U[J, ]'` exactly, with no diagonal term left in it, so
# its rank never exceeds the factor rank whatever the residual variances are.
# Two splits are taken because each is blind to a different structure: halves
# see an unstructured block, and the odd/even split sees the decaying bands that
# separate across halves. One small SVD each settles a block that no refinement
# could have recovered.
.covariance_off_diagonal_rank_bound <- function(covariance, tolerance) {

  size <- nrow(covariance)
  if (size < 4L) {
    return(0L)
  }
  scale <- max(abs(covariance))
  # A singular value moves by at most the spectral norm of the perturbation,
  # itself at most `size` times its largest entry, so a block whose
  # reconstruction the certificate would accept keeps every singular value past
  # its true rank below this threshold.
  threshold <- 16 * size * tolerance * scale
  rows      <- seq_len(size)
  splits    <- list(rows <= size %/% 2L, rows %% 2L == 1L)

  bound <- 0L
  for (split in splits) {
    block <- covariance[split, !split, drop = FALSE]
    if (!length(block)) {
      next
    }
    values <- tryCatch(svd(block, nu = 0L, nv = 0L)[["d"]],
                       error = function(e) NULL)
    if (is.null(values) || anyNA(values) || any(!is.finite(values))) {
      # An undecided split constrains nothing; leave the search unbounded.
      return(0L)
    }
    bound <- max(bound, sum(values > threshold))
  }

  as.integer(bound)
}


# Largest rank at which the off-diagonal entries of a block still constrain the
# loadings, `(2n + 1 - sqrt(8n + 1)) / 2` (Ledermann, 1937).
.covariance_ledermann_rank <- function(size) {

  as.integer(floor((2 * size + 1 - sqrt(8 * size + 1)) / 2))
}


# Certify one candidate and shape it like a recovered factor, or return NULL.
.covariance_low_rank_accept <- function(covariance, loading, scale, tolerance) {

  size     <- nrow(covariance)
  diagonal <- diag(covariance) - rowSums(loading^2)

  # A column supported by one row contributes only to that row's variance.
  # Absorbing it is exact and keeps every retained column a real dependence.
  alone <- which(colSums(loading != 0) < 2L)
  if (length(alone) > 0L) {
    diagonal <- diagonal + rowSums(loading[, alone, drop = FALSE]^2)
    loading  <- loading[, -alone, drop = FALSE]
  }
  # A residual variance at the roundoff floor is not a certifiably positive
  # one, and the representation it belongs to is singular: the conditional
  # distribution the factor routes evaluate would have no residual spread on
  # that row. Decline and let a higher rank, or the general route, own it.
  if (ncol(loading) == 0L || anyNA(diagonal) || any(!is.finite(diagonal)) ||
      any(diagonal <= tolerance * max(diag(covariance)))) {
    return(NULL)
  }

  reconstruction <- tcrossprod(loading)
  diag(reconstruction) <- diag(reconstruction) + diagonal
  residual <- max(abs(reconstruction - covariance)) / scale
  if (!is.finite(residual) || residual > tolerance) {
    return(NULL)
  }

  supports <- lapply(seq_len(ncol(loading)), function(column) {
    which(loading[, column] != 0)
  })
  list(
    method   = "minimum_rank",
    diagonal = diagonal,
    loading  = loading,
    rank     = ncol(loading),
    levels   = numeric(0),
    supports = supports,
    support  = .covariance_support_shape(supports),
    residual = residual,
    depth    = NA_integer_
  )
}


# Candidate loadings for one rank. Rank one is determined in closed form by any
# triple of distinct rows; higher ranks start from the classical principal
# factor iteration, from the plain truncated spectrum, and from fixed residual
# fractions of the supplied variances. The exact representations of a block are
# a manifold when the off-diagonal system is underdetermined, and only part of
# it keeps every residual variance positive, so several starts are tried and
# each is also offered shrunk inside the feasible region.
.covariance_low_rank_starts <- function(covariance, rank) {

  if (rank == 1L) {
    start <- .covariance_rank_one_start(covariance)
    return(if (is.null(start)) list() else list(start))
  }

  size     <- nrow(covariance)
  variance <- diag(covariance)
  truncate <- function(diagonal) {
    reduced <- covariance
    diag(reduced) <- variance - diagonal
    decomposition <- eigen(reduced, symmetric = TRUE)
    values <- pmax(decomposition[["values"]][seq_len(rank)], 0)
    if (anyNA(values) || any(!is.finite(values))) {
      return(NULL)
    }
    decomposition[["vectors"]][, seq_len(rank), drop = FALSE] *
      rep(sqrt(values), each = size)
  }

  spectral <- truncate(rep(0, size))
  if (is.null(spectral)) {
    return(list())
  }
  diagonal <- pmin(pmax(variance - rowSums(spectral^2), 0), variance)
  factored <- spectral
  for (iteration in seq_len(200L)) {
    candidate <- truncate(diagonal)
    if (is.null(candidate)) {
      break
    }
    factored <- candidate
    updated  <- pmin(pmax(variance - rowSums(factored^2), 0), variance)
    if (max(abs(updated - diagonal)) <= .Machine$double.eps * max(variance)) {
      break
    }
    diagonal <- updated
  }

  starts <- c(list(factored, spectral),
              lapply(c(.25, .6), function(fraction) {
                truncate(fraction * variance)
              }))
  starts <- Filter(Negate(is.null), starts)
  # Row-scaling a start into the interior costs nothing and changes which
  # exact representation the refinement converges to when several exist.
  interior <- lapply(starts, function(start) {
    total <- rowSums(start^2)
    scale <- sqrt(pmin(1, .5 * variance / pmax(total, .Machine$double.xmin)))
    start * scale
  })

  c(starts, interior)
}


.covariance_rank_one_start <- function(covariance) {

  size <- nrow(covariance)
  off  <- abs(covariance)
  diag(off) <- 0
  if (max(off) <= 0) {
    return(NULL)
  }

  if (size == 2L) {
    # Two rows leave one degree of freedom. Splitting the covariance in the
    # proportion of the supplied variances keeps both residual variances
    # non-negative whenever any exact representation does.
    correlation <- covariance[1L, 2L] /
      sqrt(covariance[1L, 1L] * covariance[2L, 2L])
    if (!is.finite(correlation) || abs(correlation) > 1) {
      return(NULL)
    }
    magnitude <- sqrt(abs(correlation) * diag(covariance))
    return(matrix(
      c(magnitude[[1L]], sign(correlation) * magnitude[[2L]]), 2L, 1L
    ))
  }

  # u_p^2 = V_pj V_pk / V_jk holds for any three distinct rows of an exact
  # rank-one representation. Choose the most strongly related rows for it.
  pivot  <- which.max(apply(off, 1L, max))
  others <- setdiff(seq_len(size), pivot)
  block  <- off[others, others, drop = FALSE]
  if (max(block) <= 0) {
    return(NULL)
  }
  pair <- which(block == max(block), arr.ind = TRUE)[1L, ]
  j <- others[[pair[[1L]]]]
  k <- others[[pair[[2L]]]]
  squared <- covariance[pivot, j] * covariance[pivot, k] / covariance[j, k]
  if (!is.finite(squared) || squared <= 0) {
    return(NULL)
  }

  value   <- sqrt(squared)
  loading <- covariance[pivot, ] / value
  loading[[pivot]] <- value
  matrix(loading, ncol = 1L)
}


# Refine the loadings against the off-diagonal entries alone. The diagonal is
# not a residual: it is recovered from the loadings afterwards, so the system
# has `n (n - 1) / 2` equations in `n r` unknowns and Gauss-Newton converges
# quadratically near an exact representation. Rotations of U leave the residual
# unchanged, so the normal equations are damped rather than solved exactly.
# That quadratic convergence is also what bounds the work: a start outside the
# basin stalls, and only an exact representation keeps earning its iterations.
.covariance_low_rank_polish <- function(covariance, loading,
                                        iterations = 80L, stall = 12L) {

  size  <- nrow(covariance)
  rank  <- ncol(loading)
  scale <- max(abs(covariance))
  pairs <- which(upper.tri(covariance), arr.ind = TRUE)
  rows  <- pairs[, 1L]
  cols  <- pairs[, 2L]
  target <- covariance[cbind(rows, cols)]

  residual_of <- function(value) {
    target - rowSums(value[rows, , drop = FALSE] * value[cols, , drop = FALSE])
  }
  current <- residual_of(loading)
  best    <- max(abs(current))
  if (!is.finite(best)) {
    return(loading)
  }
  damping <- 1e-9 * max(scale, 1e-300)
  stalled <- 0L

  for (iteration in seq_len(iterations)) {
    if (best <= .Machine$double.eps * scale) {
      break
    }
    # The normal equations of the Gauss-Newton system have a closed block form,
    # so neither the Jacobian nor its cross-product is ever materialized from
    # the pair list: (A'A)[(k,a),(l,b)] is U[l,a] U[k,b] off the diagonal and
    # the k-th deleted Gram entry on it, and (A'r)[k,a] is (R U)[k,a].
    residual_matrix <- matrix(0, size, size)
    residual_matrix[cbind(rows, cols)] <- current
    residual_matrix[cbind(cols, rows)] <- current
    gradient <- as.vector(residual_matrix %*% loading)
    gram     <- crossprod(loading)
    normal   <- matrix(0, size * rank, size * rank)
    for (a in seq_len(rank)) {
      for (b in seq_len(rank)) {
        value <- tcrossprod(loading[, b], loading[, a])
        diag(value) <- gram[a, b] - loading[, a] * loading[, b]
        normal[(a - 1L) * size + seq_len(size),
               (b - 1L) * size + seq_len(size)] <- value
      }
    }

    step <- tryCatch(
      solve(normal + diag(damping, nrow(normal)), gradient),
      error = function(e) NULL
    )
    if (is.null(step) || anyNA(step) || any(!is.finite(step))) {
      break
    }
    candidate <- loading + matrix(step, size, rank)
    proposed  <- residual_of(candidate)
    if (all(is.finite(proposed)) && max(abs(proposed)) < best) {
      halved  <- max(abs(proposed)) <= best / 2
      loading <- candidate
      current <- proposed
      best    <- max(abs(proposed))
      damping <- damping / 10
      # Near an exact representation the residual is squared every step, so a
      # run of steps that do not even halve it means the start was not in that
      # basin and the remaining iterations would only refine an inexact fit.
      stalled <- if (halved) 0L else stalled + 1L
    } else {
      damping <- damping * 10
      stalled <- stalled + 1L
      if (damping > scale) {
        break
      }
    }
    if (stalled >= stall) {
      break
    }
  }

  loading
}


# Cluster the observed off-diagonal correlations into distinct levels. Values
# within the roundoff tolerance of their neighbour join one level and are
# represented by its mean; the reconstruction certificate, not this clustering,
# decides whether the resulting representation is accepted.
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

  partitions <- vector("list", length(levels))
  previous   <- NULL
  for (level in seq_along(levels)) {
    related <- correlation >= levels[[level]] - tolerance
    diag(related) <- TRUE
    # Disjoint cliques: label each row by the first row it relates to, then
    # require the relation to be exactly that labelling. A relation that only
    # connects rows transitively is declined, never completed.
    membership <- max.col(related, ties.method = "first")
    if (!all(related == outer(membership, membership, `==`))) {
      return(NULL)
    }
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

  # A forest: every pair of supports is either disjoint or nested. Nested
  # partitions can only produce this shape, so "general" is a guard.
  for (outer in seq_len(length(supports) - 1L)) {
    for (inner in seq.int(outer + 1L, length(supports))) {
      if (!any(supports[[inner]] %in% supports[[outer]])) {
        next
      }
      if (!all(supports[[inner]] %in% supports[[outer]])) {
        return("general")
      }
    }
  }

  "tree"
}
