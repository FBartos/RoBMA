# Block covariances for fixtures that carry a dense cube.
#
# The joint selection prediction parts carry `RoBMA_block_covariance` objects,
# whose partition the package declares from the fitted dependency structure
# rather than reading it off the values. A fixture that builds its covariance
# by hand has no such declaration, so this helper derives the partition from
# the cube's nonzero pattern. Package code must not do this: the point of the
# representation is that nothing scans for a cross-block zero.

.as_block_covariance <- function(samples) {

  if (.is_block_covariance(samples)) {
    return(samples)
  }
  dimension <- dim(samples)
  if (!identical(length(dimension), 3L)) {
    stop("A fixture covariance must be a draw x row x row array.", call. = FALSE)
  }
  K         <- dimension[[2L]]
  adjacency <- matrix(FALSE, K, K)
  for (column in seq_len(K)) {
    for (row in seq_len(K)) {
      adjacency[row, column] <- any(samples[, row, column] != 0)
    }
  }
  adjacency <- adjacency | t(adjacency)

  .block_covariance_from_array(samples, .known_v_block_indices(adjacency * 1))
}
