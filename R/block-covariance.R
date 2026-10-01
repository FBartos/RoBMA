# ---------------------------------------------------------------------------- #
# RoBMA_block_covariance
# ---------------------------------------------------------------------------- #
#
# A posterior-draw covariance over K outcomes that is block-diagonal by
# construction, stored one block at a time instead of as a dense S x K x K cube.
#
# The joint selection prediction parts are built from sources whose dependency
# structure the fitted model already declares: a per-row random variance, a
# cluster or known-R covariance, and the sampling covariance. Their dense cubes
# are 1.12 GiB each at S = 15000 and K = 100 while the block-diagonal content is
# about a tenth of that, and no consumer reads an element whose two rows lie in
# different blocks. This representation keeps only the blocks.
#
# Fields:
#   blocks        list of sorted integer row vectors partitioning seq_len(K)
#   values        one double array per block, dim c(R_b, n_b, n_b), where R_b is
#                 S for a covariance that varies over draws and 1 for one that
#                 does not (the sampling covariance is stored once, not
#                 replicated over the draws)
#   S, K          the dense dimensions this object stands for
#   zero          per block, whether every stored value is exactly zero
#   row_block     per row, the index of the block that holds it
#   row_position  per row, its position inside that block
#
# Results are the dense cube's results. Two operations are deliberately *not*
# done block-wise: `.block_covariance_dense()` assembles the K x K matrix for a
# per-draw Cholesky or spectral factorization, because LAPACK's blocked
# factorization of the assembled matrix and of its blocks differ in the last
# bits (measured: 8.9e-16 on a K = 100 partition), and callers that factorize
# require bit-identical results. `.block_covariance_chol_solve()` is the
# block-wise solve, for callers that do not.
#
# ---------------------------------------------------------------------------- #


.block_covariance <- function(blocks, values, S, K) {

  S <- as.integer(S)
  K <- as.integer(K)
  if (!is.list(blocks) || length(blocks) == 0L ||
      length(values) != length(blocks)) {
    stop("Block covariance metadata must pair one value array per block.",
         call. = FALSE)
  }
  .known_v_validate_dependency_blocks(blocks, K)

  row_block    <- integer(K)
  row_position <- integer(K)
  zero         <- logical(length(blocks))
  for (index in seq_along(blocks)) {
    rows  <- as.integer(blocks[[index]])
    value <- values[[index]]
    n     <- length(rows)
    if (!is.numeric(value) || !identical(length(dim(value)), 3L) ||
        !dim(value)[[1L]] %in% c(1L, S) ||
        dim(value)[[2L]] != n || dim(value)[[3L]] != n) {
      stop("Block covariance values must be draw x row x row arrays.",
           call. = FALSE)
    }
    if (any(!is.finite(value))) {
      stop("Block covariance values must be finite.", call. = FALSE)
    }
    blocks[[index]]       <- rows
    values[[index]]       <- value
    row_block[rows]       <- index
    row_position[rows]    <- seq_len(n)
    zero[[index]]         <- !any(value != 0)
  }

  structure(
    list(
      blocks       = blocks,
      values       = values,
      S            = S,
      K            = K,
      zero         = zero,
      row_block    = row_block,
      row_position = row_position
    ),
    class = "RoBMA_block_covariance"
  )
}


.is_block_covariance <- function(x) {

  inherits(x, "RoBMA_block_covariance")
}


.block_covariance_singleton_blocks <- function(K) {

  lapply(seq_len(K), function(row) as.integer(row))
}


# An all-zero covariance. Its blocks carry the partition the consumers read,
# so a singleton-only prediction stays singleton-blocked.
.block_covariance_zero <- function(S, K, blocks = NULL) {

  if (is.null(blocks)) blocks <- .block_covariance_singleton_blocks(K)
  .block_covariance(
    blocks = blocks,
    values = lapply(blocks, function(rows) {
      array(0, dim = c(1L, length(rows), length(rows)))
    }),
    S = S, K = K
  )
}


# From an S x K matrix of per-row variances. With no blocks the rows are
# independent and the partition is the singletons.
.block_covariance_from_diagonal <- function(diagonal, blocks = NULL) {

  diagonal <- as.matrix(diagonal)
  S <- nrow(diagonal)
  K <- ncol(diagonal)
  if (is.null(blocks)) blocks <- .block_covariance_singleton_blocks(K)
  .block_covariance(
    blocks = blocks,
    values = lapply(blocks, function(rows) {
      n     <- length(rows)
      value <- array(0, dim = c(S, n, n))
      for (position in seq_len(n)) {
        value[, position, position] <- diagonal[, rows[[position]]]
      }
      value
    }),
    S = S, K = K
  )
}


# From a K x K matrix that does not vary over the posterior draws. The matrix
# is stored once; nothing is replicated over S.
.block_covariance_from_matrix <- function(matrix, S, blocks = NULL) {

  K <- nrow(matrix)
  if (!is.matrix(matrix) || ncol(matrix) != K) {
    stop("A constant block covariance must come from a square matrix.",
         call. = FALSE)
  }
  if (is.null(blocks)) blocks <- .known_v_block_indices(matrix)
  .block_covariance(
    blocks = blocks,
    values = lapply(blocks, function(rows) {
      n <- length(rows)
      array(matrix[rows, rows], dim = c(1L, n, n))
    }),
    S = S, K = K
  )
}


# From a dense S x K x K cube whose declared dependency structure is `blocks`.
# The declaration is verified: every element outside the blocks must be exactly
# zero, because the consumers of this representation no longer scan for one.
.block_covariance_from_array <- function(samples, blocks) {

  dimension <- dim(samples)
  if (!identical(length(dimension), 3L) || dimension[[2L]] != dimension[[3L]]) {
    stop("A block covariance must come from a draw x row x row array.",
         call. = FALSE)
  }
  S <- dimension[[1L]]
  K <- dimension[[2L]]
  .known_v_validate_dependency_blocks(blocks, K)

  membership <- integer(K)
  for (index in seq_along(blocks)) {
    membership[as.integer(blocks[[index]])] <- index
  }
  for (column in seq_len(K)) {
    for (row in seq_len(K)) {
      if (membership[[row]] == membership[[column]]) {
        next
      }
      if (any(samples[, row, column] != 0)) {
        stop("Covariance samples are nonzero outside their declared ",
             "dependency blocks.", call. = FALSE)
      }
    }
  }

  .block_covariance(
    blocks = blocks,
    values = lapply(blocks, function(rows) {
      n <- length(rows)
      array(samples[, rows, rows], dim = c(S, n, n))
    }),
    S = S, K = K
  )
}


.block_covariance_dim <- function(x) {

  c(x[["S"]], x[["K"]], x[["K"]])
}


# The coarsest partition on which both partitions are block-diagonal: the
# connected components of the union of the two block structures. A covariance
# that is block-diagonal on either one is block-diagonal on this.
.block_covariance_join_blocks <- function(left, right, K) {

  component <- seq_len(K)
  merge <- function(component, blocks) {
    for (rows in blocks) {
      rows   <- as.integer(rows)
      target <- min(component[rows])
      component[component %in% component[rows]] <- target
    }
    component
  }
  component <- merge(component, left)
  component <- merge(component, right)
  # Repeat until stable: merging on the right can connect components the left
  # pass had kept apart, and vice versa.
  repeat {
    updated <- merge(merge(component, left), right)
    if (identical(updated, component)) break
    component <- updated
  }
  unname(lapply(split(seq_len(K), component), as.integer))
}


# Re-express the covariance on a coarser partition. Every block of `x` has to
# lie inside one block of `blocks`; the assembled block is zero between them.
.block_covariance_coarsen <- function(x, blocks) {

  if (identical(x[["blocks"]], blocks)) {
    return(x)
  }
  K <- x[["K"]]
  .known_v_validate_dependency_blocks(blocks, K)
  values <- lapply(blocks, function(rows) {
    rows   <- as.integer(rows)
    n      <- length(rows)
    source <- unique(x[["row_block"]][rows])
    draws  <- max(vapply(source, function(index) {
      dim(x[["values"]][[index]])[[1L]]
    }, integer(1L)))
    value <- array(0, dim = c(draws, n, n))
    for (index in source) {
      inside <- which(x[["row_block"]][rows] == index)
      if (!identical(sort(rows[inside]), sort(as.integer(x[["blocks"]][[index]])))) {
        stop("A block covariance can only be coarsened onto a partition that ",
             "contains its own blocks.", call. = FALSE)
      }
      position <- x[["row_position"]][rows[inside]]
      part     <- x[["values"]][[index]]
      if (dim(part)[[1L]] == draws) {
        value[, inside, inside] <- part[, position, position, drop = FALSE]
      } else {
        value[, inside, inside] <- array(
          rep(part[1L, position, position], each = draws),
          dim = c(draws, length(inside), length(inside))
        )
      }
    }
    value
  })

  .block_covariance(blocks = blocks, values = values, S = x[["S"]], K = K)
}


# The sum of two block covariances, on the coarsest partition that holds both.
.block_covariance_add <- function(x, y) {

  if (!.is_block_covariance(x) || !.is_block_covariance(y) ||
      !identical(x[["S"]], y[["S"]]) || !identical(x[["K"]], y[["K"]])) {
    stop("Block covariances can only be added to matching block covariances.",
         call. = FALSE)
  }
  blocks <- if (identical(x[["blocks"]], y[["blocks"]])) {
    x[["blocks"]]
  } else {
    .block_covariance_join_blocks(x[["blocks"]], y[["blocks"]], x[["K"]])
  }
  left  <- .block_covariance_coarsen(x, blocks)
  right <- .block_covariance_coarsen(y, blocks)

  .block_covariance(
    blocks = blocks,
    values = lapply(seq_along(blocks), function(index) {
      a <- left[["values"]][[index]]
      b <- right[["values"]][[index]]
      if (identical(dim(a)[[1L]], dim(b)[[1L]])) {
        return(a + b)
      }
      draws    <- max(dim(a)[[1L]], dim(b)[[1L]])
      expanded <- function(value) {
        if (dim(value)[[1L]] == draws) {
          return(value)
        }
        array(rep(value, each = draws), dim = c(draws, dim(value)[-1L]))
      }
      expanded(a) + expanded(b)
    }),
    S = x[["S"]], K = x[["K"]]
  )
}


# Add an S x K matrix of per-row variances to the diagonal.
.block_covariance_add_diagonal <- function(x, diagonal) {

  diagonal <- as.matrix(diagonal)
  if (!identical(dim(diagonal), c(x[["S"]], x[["K"]]))) {
    stop("A block covariance diagonal update must be draw x row.", call. = FALSE)
  }
  .block_covariance(
    blocks = x[["blocks"]],
    values = lapply(seq_along(x[["blocks"]]), function(index) {
      rows  <- x[["blocks"]][[index]]
      value <- x[["values"]][[index]]
      n     <- length(rows)
      if (dim(value)[[1L]] != x[["S"]]) {
        value <- array(rep(value, each = x[["S"]]), dim = c(x[["S"]], n, n))
      }
      for (position in seq_len(n)) {
        value[, position, position] <- value[, position, position] +
          diagonal[, rows[[position]]]
      }
      value
    }),
    S = x[["S"]], K = x[["K"]]
  )
}


# The S x K matrix of per-row variances, as `cube[, k, k]` gave it.
.block_covariance_diag_matrix <- function(x) {

  out <- matrix(0, x[["S"]], x[["K"]])
  for (index in seq_along(x[["blocks"]])) {
    rows  <- x[["blocks"]][[index]]
    value <- x[["values"]][[index]]
    fixed <- dim(value)[[1L]] != x[["S"]]
    for (position in seq_along(rows)) {
      out[, rows[[position]]] <- if (fixed) {
        rep(value[1L, position, position], x[["S"]])
      } else {
        value[, position, position]
      }
    }
  }
  out
}


# The dense n x n covariance of `rows` in one draw, as
# `cube[draw, rows, rows]` gave it. `rows` may be a block, a union of blocks,
# or a subset; everything between two blocks is zero.
.block_covariance_sub <- function(x, draw, rows) {

  rows <- as.integer(rows)
  n    <- length(rows)
  out  <- matrix(0, n, n)
  owner <- x[["row_block"]][rows]
  for (index in unique(owner)) {
    inside   <- which(owner == index)
    position <- x[["row_position"]][rows[inside]]
    value    <- x[["values"]][[index]]
    source   <- if (dim(value)[[1L]] == x[["S"]]) draw else 1L
    out[inside, inside] <- matrix(
      value[source, position, position], length(inside), length(inside)
    )
  }
  out
}


# The dense K x K covariance of one draw. Callers that factorize the matrix use
# this: a blocked LAPACK factorization of the assembled matrix and of its blocks
# do not agree bit for bit, and their results must.
.block_covariance_dense <- function(x, draw) {

  .block_covariance_sub(x, draw, seq_len(x[["K"]]))
}


# Whether every element of every draw is exactly zero, as `all(cube == 0)`.
.block_covariance_is_zero <- function(x) {

  all(x[["zero"]])
}


# Per draw, whether every element is exactly zero, as
# `all(cube[draw, , ] == 0)`.
.block_covariance_zero_draws <- function(x) {

  out <- rep(TRUE, x[["S"]])
  for (index in seq_along(x[["blocks"]])) {
    if (x[["zero"]][[index]]) {
      next
    }
    value <- x[["values"]][[index]]
    n     <- length(x[["blocks"]][[index]])
    if (dim(value)[[1L]] != x[["S"]]) {
      return(rep(FALSE, x[["S"]]))
    }
    nonzero <- rep(FALSE, x[["S"]])
    for (column in seq_len(n)) {
      for (row in seq_len(n)) {
        nonzero <- nonzero | value[, row, column] != 0
      }
    }
    out <- out & !nonzero
  }
  out
}


# Per draw, whether every off-diagonal element is exactly zero, as
# `all(matrix(cube[draw, , ], K, K)[off_diagonal] == 0)`. A singleton-only
# partition answers structurally.
.block_covariance_diagonal_draws <- function(x) {

  out <- rep(TRUE, x[["S"]])
  for (index in seq_along(x[["blocks"]])) {
    n <- length(x[["blocks"]][[index]])
    if (n == 1L || x[["zero"]][[index]]) {
      next
    }
    value <- x[["values"]][[index]]
    fixed <- dim(value)[[1L]] != x[["S"]]
    nonzero <- rep(FALSE, x[["S"]])
    for (column in seq_len(n)) {
      for (row in seq_len(n)) {
        if (row == column) {
          next
        }
        nonzero <- nonzero | if (fixed) {
          rep(value[1L, row, column] != 0, x[["S"]])
        } else {
          value[, row, column] != 0
        }
      }
    }
    out <- out & !nonzero
  }
  out
}


# Whether every stored element is finite, as `all(is.finite(cube))`.
.block_covariance_all_finite <- function(x) {

  for (value in x[["values"]]) {
    if (any(!is.finite(value))) {
      return(FALSE)
    }
  }
  TRUE
}


# The dense product `cube[draw, , ] %*% v`. Block-wise and dense products agree
# bit for bit: the elements the block form leaves out are exactly zero.
.block_covariance_matvec <- function(x, draw, v) {

  out <- numeric(x[["K"]])
  for (index in seq_along(x[["blocks"]])) {
    if (x[["zero"]][[index]]) {
      next
    }
    rows   <- x[["blocks"]][[index]]
    value  <- x[["values"]][[index]]
    n      <- length(rows)
    source <- if (dim(value)[[1L]] == x[["S"]]) draw else 1L
    out[rows] <- as.vector(
      matrix(value[source, , ], n, n) %*% v[rows]
    )
  }
  out
}


# The block-wise solve `cube[draw, , ]^-1 %*% rhs` through each block's
# Cholesky factor.
.block_covariance_chol_solve <- function(x, draw, rhs) {

  out <- numeric(x[["K"]])
  for (index in seq_along(x[["blocks"]])) {
    rows   <- x[["blocks"]][[index]]
    value  <- x[["values"]][[index]]
    n      <- length(rows)
    source <- if (dim(value)[[1L]] == x[["S"]]) draw else 1L
    block  <- matrix(value[source, , ], n, n)
    factor <- tryCatch(chol(block), error = function(error) NULL)
    if (is.null(factor)) {
      stop("A block covariance is not positive definite.", call. = FALSE)
    }
    out[rows] <- backsolve(factor, forwardsolve(t(factor), rhs[rows]))
  }
  out
}


# The same covariance on a subset of the posterior draws, as
# `cube[index, , , drop = FALSE]`.
.block_covariance_draws <- function(x, index) {

  index <- as.integer(index)
  .block_covariance(
    blocks = x[["blocks"]],
    values = lapply(x[["values"]], function(value) {
      if (dim(value)[[1L]] != x[["S"]]) {
        return(value)
      }
      array(value[index, , , drop = FALSE],
            dim = c(length(index), dim(value)[-1L]))
    }),
    S = length(index), K = x[["K"]]
  )
}


# The per-block calling form of the native selected-response RNG: the value
# arrays and the rows they belong to. The kernel reads its dependency blocks
# from its own argument; this partition only has to refine them.
.block_covariance_native_parts <- function(x) {

  list(x[["values"]], lapply(x[["blocks"]], .native_integer_vector))
}
