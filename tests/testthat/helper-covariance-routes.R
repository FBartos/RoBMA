# Covariance blocks that keep the general (dense) route.
#
# Exact recovery represents a dependency block as `diag(d) + U U'` whenever
# that reproduces its supplied entries and the rank stays inside the evidence
# bound, so a fixture that must reach the dense kernels cannot rely on a sign
# pattern or on a small block: every two-row block is exactly rank one. Markov
# correlations on three or more rows are the simplest structure that declines.
# Their only rank-one representation needs a zero residual variance on an
# interior row, and every higher rank is above the bound.
#
# Each constructor re-checks that its block declines, so a fixture cannot drift
# onto the factor routes without the test that uses it saying so.

.dense_route_covariance <- function(variance = c(.02, .03, .04), rho = .8) {

  size <- length(variance)
  if (size < 3L) {
    stop("A dense-route block needs at least three rows.", call. = FALSE)
  }
  block <- .known_v_exact_symmetrize(
    outer(seq_len(size), seq_len(size), function(i, j) rho^abs(i - j)) *
      tcrossprod(sqrt(variance))
  )
  if (!is.null(.covariance_block_constant_factor(block)) ||
      !is.null(.covariance_minimum_rank_factor(block))) {
    stop("A dense-route block must decline exact recovery.", call. = FALSE)
  }
  block
}


# Constant negative correlations decline for a different reason, and their
# normalizer has no monotone covariance envelope, so an anchor error over such
# a block stays unknown. Rank one would need a sign assignment no three rows
# admit, and every higher rank is above the evidence bound.
.dense_route_negative_covariance <- function(variance = c(.04, .09, .05),
                                             rho = -.3) {

  size <- length(variance)
  if (size < 3L || rho >= 0) {
    stop("A negative dense-route block needs at least three rows and rho < 0.",
         call. = FALSE)
  }
  correlation <- matrix(rho, size, size)
  diag(correlation) <- 1
  block <- .known_v_exact_symmetrize(correlation * tcrossprod(sqrt(variance)))
  if (!is.null(.covariance_block_constant_factor(block)) ||
      !is.null(.covariance_minimum_rank_factor(block))) {
    stop("A dense-route block must decline exact recovery.", call. = FALSE)
  }
  block
}


.dense_route_block_diagonal <- function(variance, cluster, rho = .8,
                                        negative = FALSE) {

  out <- diag(variance, length(variance))
  for (group in unique(cluster)) {
    rows <- which(cluster == group)
    out[rows, rows] <- if (length(rows) == 1L) {
      variance[rows]
    } else if (negative) {
      .dense_route_negative_covariance(variance[rows], rho = -abs(rho))
    } else {
      .dense_route_covariance(variance[rows], rho = rho)
    }
  }
  out
}
