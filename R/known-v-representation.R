# Known-V input helpers -----

.new_known_v <- function(fields) {

  if (!is.list(fields) || is.null(names(fields)) || any(!nzchar(names(fields))) ||
      anyDuplicated(names(fields))) {
    stop("Internal error: known-V fields must be a uniquely named list.",
         call. = FALSE)
  }

  known_V <- structure(fields, class = c("brma_known_v", "list"))
  .validate_known_v(known_V)
}


.validate_known_v <- function(known_V) {

  if (!is.list(known_V)) {
    stop("Internal error: known-V representation must be a list.", call. = FALSE)
  }

  required <- c("version", "storage", "K", "diagonal")
  if (!all(required %in% names(known_V))) {
    stop("Internal error: known-V representation is incomplete.", call. = FALSE)
  }
  if (!identical(as.integer(known_V[["version"]]), 2L)) {
    stop("Internal error: unsupported known-V representation version.",
         call. = FALSE)
  }

  storage <- known_V[["storage"]]
  K       <- known_V[["K"]]
  diagonal <- known_V[["diagonal"]]
  if (length(storage) != 1L || !storage %in% c("diagonal", "blocks", "dense") ||
      length(K) != 1L || is.na(K) || K < 1L || K != as.integer(K) ||
      !is.numeric(diagonal) || length(diagonal) != K || anyNA(diagonal) ||
      any(!is.finite(diagonal)) || any(diagonal < 0)) {
    stop("Internal error: invalid known-V representation metadata.",
         call. = FALSE)
  }

  if (identical(storage, "dense")) {
    V <- known_V[["V"]]
    if (!is.matrix(V) || !is.numeric(V) || !identical(dim(V), c(K, K))) {
      stop("Internal error: dense known-V representation is invalid.",
           call. = FALSE)
    }
  } else if (!is.null(known_V[["V"]])) {
    stop("Internal error: compact known-V representation contains dense state.",
         call. = FALSE)
  }

  invisible(known_V)
}


.validate_prepared_known_v <- function(known_V) {

  .validate_known_v(known_V)
  required <- c(
    "parameterization", "parameterization_requested", "effective_backend",
    "correlated", "singular", "residual_variance", "residual_sei", "rank"
  )
  if (!all(required %in% names(known_V))) {
    stop("Internal error: prepared known-V metadata are incomplete.",
         call. = FALSE)
  }

  parameterization  <- known_V[["parameterization"]]
  requested         <- known_V[["parameterization_requested"]]
  backend           <- known_V[["effective_backend"]]
  K                 <- .known_v_nrow(known_V)
  residual_variance <- known_V[["residual_variance"]]
  residual_sei      <- known_V[["residual_sei"]]
  rank              <- known_V[["rank"]]
  if (length(parameterization) != 1L ||
      !parameterization %in% c("latent", "whitened", "block_mvn") ||
      length(requested) != 1L ||
      !requested %in% c("auto", "latent", "whitened", "block_mvn") ||
      length(backend) != 1L ||
      !backend %in% c("diagonal", "latent", "whitened", "block_mvn") ||
      !is.logical(known_V[["correlated"]]) ||
      length(known_V[["correlated"]]) != 1L ||
      is.na(known_V[["correlated"]]) ||
      !is.logical(known_V[["singular"]]) ||
      length(known_V[["singular"]]) != 1L ||
      is.na(known_V[["singular"]]) ||
      !is.numeric(residual_variance) || length(residual_variance) != K ||
      anyNA(residual_variance) || any(!is.finite(residual_variance)) ||
      any(residual_variance < 0) ||
      !is.numeric(residual_sei) || length(residual_sei) != K ||
      anyNA(residual_sei) || any(!is.finite(residual_sei)) ||
      any(residual_sei < 0) ||
      length(rank) != 1L || is.na(rank) || rank < 0L ||
      rank != as.integer(rank)) {
    stop("Internal error: prepared known-V metadata are invalid.",
         call. = FALSE)
  }
  if (!identical(backend, "diagonal")) {
    .known_v_backend_blocks(known_V, backend)
  }

  invisible(known_V)
}


.known_v_update <- function(known_V, fields) {

  if (!is.list(known_V) || !is.list(fields) || is.null(names(fields)) ||
      any(!nzchar(names(fields))) || anyDuplicated(names(fields))) {
    stop("Internal error: invalid known-V update.", call. = FALSE)
  }

  known_V[names(fields)] <- fields
  class(known_V) <- unique(c("brma_known_v", class(known_V), "list"))
  .validate_known_v(known_V)
}


.known_v_parameterization <- function(known_V) {

  known_V[["parameterization"]]
}


.known_v_effective_backend <- function(known_V) {

  known_V[["effective_backend"]]
}


.known_v_requested_parameterization <- function(known_V) {

  known_V[["parameterization_requested"]]
}


.known_v_is_correlated <- function(known_V) {

  isTRUE(known_V[["correlated"]])
}


.known_v_is_singular_representation <- function(known_V) {

  isTRUE(known_V[["singular"]])
}


.known_v_rank <- function(known_V) {

  as.integer(known_V[["rank"]])
}


.known_v_residual_variance <- function(known_V) {

  known_V[["residual_variance"]]
}


.known_v_backend_blocks <- function(
    known_V, backend = c("latent", "whitened", "block_mvn")) {

  backend <- match.arg(backend)
  field <- switch(
    backend,
    latent    = "latent_blocks",
    whitened  = "whitening_blocks",
    block_mvn = "block_mvn_blocks"
  )
  blocks <- known_V[[field]]
  if (!is.list(blocks)) {
    stop(
      "Internal error: current known-V backend metadata are incomplete.",
      call. = FALSE
    )
  }

  blocks
}


.known_v_requested_residual_fraction <- function(known_V) {

  known_V[["residual_fraction_requested"]]
}


.known_v_as_matrix <- function(V, k = NULL, warn_singular = TRUE) {

  V_matrix <- .known_v_as_matrix_structure(V, k = k)

  if (anyNA(V_matrix) || any(!is.finite(V_matrix))) {
    stop("The 'V' argument must contain only finite non-missing values.", call. = FALSE)
  }

  .known_v_check_symmetric(V_matrix, "'V'")

  diagonal <- diag(V_matrix)
  if (any(diagonal <= 0)) {
    stop("The diagonal of 'V' must contain positive variances.", call. = FALSE)
  }

  factorization <- .known_v_correlation_factorization(V_matrix)
  if (isTRUE(factorization[["singular"]])) {
    if (.covariance_is_positive_semidefinite(factorization)) {
      if (isTRUE(warn_singular)) {
        .known_v_warn_singular()
      }
    } else {
      stop("The 'V' argument must be positive semidefinite.", call. = FALSE)
    }
  }

  return(V_matrix)
}


# Convert known-V input without validating covariance values.
.known_v_as_matrix_structure <- function(V, k = NULL) {

  if (is.matrix(V)) {
    V_matrix <- V
  } else if (is.numeric(V) && is.null(dim(V)) && length(V) > 0L) {
    V_matrix <- diag(as.numeric(V), nrow = length(V), ncol = length(V))
  } else if (is.list(V) && length(V) > 0L && all(vapply(V, is.matrix, logical(1)))) {
    V_matrix <- .known_v_blockdiag(V, arg = "V")
  } else {
    stop("The 'V' argument must be a variance vector, a square matrix, or a non-empty list of square matrices.",
         call. = FALSE)
  }

  if (!is.numeric(V_matrix)) {
    stop("The 'V' argument must be numeric.", call. = FALSE)
  }
  if (length(dim(V_matrix)) != 2L || nrow(V_matrix) != ncol(V_matrix)) {
    stop("The 'V' argument must be square.", call. = FALSE)
  }
  if (!is.null(k) && nrow(V_matrix) != k) {
    stop("The dimensions of 'V' must match the length of 'yi'.", call. = FALSE)
  }

  return(V_matrix)
}


# Describe known-V input without materializing block-diagonal storage.
.known_v_input_storage <- function(V, arg = "V") {

  if (is.matrix(V)) {
    if (!is.numeric(V)) {
      stop("The '", arg, "' argument must be numeric.", call. = FALSE)
    }
    if (nrow(V) == 0L || nrow(V) != ncol(V)) {
      stop("The '", arg, "' argument must be a non-empty square matrix.",
           call. = FALSE)
    }
    return("dense")
  }
  if (is.numeric(V) && is.null(dim(V)) && length(V) > 0L) {
    return("diagonal")
  }
  if (is.list(V) && length(V) > 0L && all(vapply(V, is.matrix, logical(1)))) {
    return("blocks")
  }

  stop(
    "The '", arg, "' argument must be a variance vector, a square matrix, ",
    "or a non-empty list of square matrices.",
    call. = FALSE
  )
}


.known_v_input_nrow <- function(V, arg = "V") {

  storage <- .known_v_input_storage(V, arg = arg)
  if (storage == "dense") {
    return(nrow(V))
  }
  if (storage == "diagonal") {
    return(length(V))
  }

  for (block in V) {
    if (!is.numeric(block)) {
      stop("All matrices in the '", arg, "' list must be numeric.", call. = FALSE)
    }
    if (nrow(block) == 0L || nrow(block) != ncol(block)) {
      stop(
        "All matrices in the '", arg,
        "' list must be non-empty square matrices.",
        call. = FALSE
      )
    }
  }

  sum(vapply(V, nrow, integer(1)))
}


.known_v_input_diagonal <- function(V, arg = "V") {

  storage <- .known_v_input_storage(V, arg = arg)
  if (storage == "dense") {
    return(diag(V))
  }
  if (storage == "diagonal") {
    return(as.numeric(V))
  }

  unlist(lapply(V, diag), use.names = FALSE)
}


.known_v_subset_input <- function(V, keep_rows) {

  K <- .known_v_input_nrow(V)
  if (!is.logical(keep_rows) || length(keep_rows) != K || anyNA(keep_rows)) {
    stop("Internal error: invalid known-V row selector.", call. = FALSE)
  }

  storage <- .known_v_input_storage(V)
  if (storage == "dense") {
    return(V[keep_rows, keep_rows, drop = FALSE])
  }
  if (storage == "diagonal") {
    return(as.numeric(V)[keep_rows])
  }

  out   <- list()
  start <- 1L
  for (block in V) {
    index      <- seq.int(start, length.out = nrow(block))
    local_keep <- keep_rows[index]
    if (any(local_keep)) {
      out[[length(out) + 1L]] <- block[local_keep, local_keep, drop = FALSE]
    }
    start <- start + nrow(block)
  }

  out
}


# Canonical known-V accessors.
.known_v_nrow <- function(known_V) {

  as.integer(known_V[["K"]])
}


.known_v_storage <- function(known_V) {

  known_V[["storage"]]
}


.known_v_diagonal <- function(known_V) {

  as.numeric(known_V[["diagonal"]])
}


.known_v_blocks <- function(known_V) {

  storage <- .known_v_storage(known_V)
  if (storage == "diagonal") {
    diagonal <- .known_v_diagonal(known_V)
    return(lapply(seq_len(.known_v_nrow(known_V)), function(i) {
      list(index = i, covariance = matrix(diagonal[[i]], 1L, 1L))
    }))
  }
  if (!is.null(known_V[["blocks"]])) {
    blocks   <- known_V[["blocks"]]
    diagonal <- .known_v_diagonal(known_V)
    for (index in .known_v_independent_indices(known_V)) {
      blocks[[length(blocks) + 1L]] <- list(
        index      = index,
        covariance = matrix(diagonal[[index]], 1L, 1L)
      )
    }
    if (length(blocks) > 1L) {
      blocks <- blocks[order(vapply(blocks, function(x) x[["index"]][[1L]], integer(1)))]
    }
    return(blocks)
  }

  V             <- known_V[["V"]]
  block_indices <- known_V[["block_indices"]]
  if (is.null(block_indices)) {
    block_indices <- .known_v_block_indices(V)
  }
  lapply(block_indices, function(index) {
    covariance <- if (identical(index, seq_len(nrow(V)))) {
      V
    } else {
      V[index, index, drop = FALSE]
    }
    list(index = index, covariance = covariance)
  })
}


.known_v_block_covariance <- function(known_V, block) {

  blocks <- .known_v_blocks(known_V)
  if (length(block) == 1L && is.numeric(block) &&
      block >= 1L && block <= length(blocks)) {
    return(blocks[[block]][["covariance"]])
  }

  match_block <- vapply(blocks, function(x) identical(x[["index"]], block), logical(1))
  if (sum(match_block) != 1L) {
    stop("Known-V dependency block was not found.", call. = FALSE)
  }
  blocks[[which(match_block)]][["covariance"]]
}


.known_v_latent_apply <- function(known_V, z_samples) {

  z_samples <- as.matrix(z_samples)
  K         <- .known_v_nrow(known_V)
  out       <- matrix(0, nrow = nrow(z_samples), ncol = K)

  latent_blocks <- .known_v_backend_blocks(known_V, "latent")
  for (block in latent_blocks) {
    if (block[["rank"]] == 0L) {
      next
    }
    z_index <- seq.int(block[["z_start"]], block[["z_end"]])
    out[, block[["index"]]] <- z_samples[, z_index, drop = FALSE] %*%
      t(block[["B"]])
  }

  out
}


.known_v_materialize <- function(known_V) {

  if (!is.null(known_V[["V"]])) {
    return(known_V[["V"]])
  }

  K   <- .known_v_nrow(known_V)
  out <- matrix(0, nrow = K, ncol = K)
  if (.known_v_storage(known_V) == "diagonal") {
    diag(out) <- .known_v_diagonal(known_V)
    return(out)
  }
  for (block in .known_v_blocks(known_V)) {
    index <- block[["index"]]
    out[index, index] <- block[["covariance"]]
  }
  out
}


.known_v_as_input <- function(known_V) {

  storage <- .known_v_storage(known_V)
  if (storage == "diagonal") {
    return(.known_v_diagonal(known_V))
  }
  if (storage == "dense") {
    return(known_V[["V"]])
  }

  lapply(.known_v_blocks(known_V), `[[`, "covariance")
}

.known_v_check_symmetric <- function(V_matrix, arg) {

  if (any(V_matrix != t(V_matrix))) {
    stop(arg, " must be symmetric.", call. = FALSE)
  }

  return(invisible(TRUE))
}


.known_v_blockdiag <- function(blocks, arg = "V") {

  sizes <- vapply(blocks, nrow, integer(1))
  for (i in seq_along(blocks)) {
    if (!is.numeric(blocks[[i]])) {
      stop("All matrices in the '", arg, "' list must be numeric.",
           call. = FALSE)
    }
    if (nrow(blocks[[i]]) == 0L || nrow(blocks[[i]]) != ncol(blocks[[i]])) {
      stop(
        "All matrices in the '", arg,
        "' list must be non-empty square matrices.",
        call. = FALSE
      )
    }
  }

  total_size <- sum(sizes)
  V_matrix   <- matrix(0, nrow = total_size, ncol = total_size)
  start      <- 1L

  for (i in seq_along(blocks)) {
    end <- start + sizes[[i]] - 1L
    V_matrix[start:end, start:end] <- blocks[[i]]
    start <- end + 1L
  }

  return(V_matrix)
}


.known_v_newdata_prepare <- function(V_new, k) {

  if (.known_v_input_nrow(V_new, arg = "V_new") != k) {
    stop(
      "The dimensions of 'V_new' must match the number of rows in 'newdata'.",
      call. = FALSE
    )
  }

  known_V <- .known_v_canonicalize_newdata(V_new)
  correlated <- length(.known_v_correlated_blocks(known_V)) > 0L
  .known_v_update(known_V, list(
    parameterization  = "block_mvn",
    effective_backend = if (correlated) "block_mvn" else "diagonal",
    correlated        = correlated,
    residual_variance = .known_v_diagonal(known_V),
    residual_sei      = sqrt(pmax(.known_v_diagonal(known_V), 0)),
    rank              = 0L
  ))
}


.known_v_canonicalize_newdata <- function(V_new) {

  storage <- .known_v_input_storage(V_new, arg = "V_new")
  K       <- .known_v_input_nrow(V_new, arg = "V_new")
  if (storage == "diagonal") {
    diagonal <- as.numeric(V_new)
    if (anyNA(diagonal) || any(!is.finite(diagonal))) {
      stop("'V_new' must contain only finite non-missing values.", call. = FALSE)
    }
    if (any(diagonal < 0)) {
      stop("The diagonal of 'V_new' must contain non-negative variances.",
           call. = FALSE)
    }
    return(.new_known_v(list(
      version  = 2L,
      storage  = "diagonal",
      K        = K,
      diagonal = diagonal,
      blocks   = list(),
      singular = any(diagonal == 0)
    )))
  }

  if (storage == "dense") {
    V_new   <- .known_v_validate_newdata_block(V_new)
    indices <- .known_v_block_indices(V_new)
    blocks  <- lapply(indices[lengths(indices) > 1L], function(index) {
      list(index = index, covariance = V_new[index, index, drop = FALSE])
    })
    retain_dense <- length(indices) == 1L &&
      length(indices[[1L]]) == K && K > 1L
    return(.new_known_v(list(
      version  = 2L,
      storage  = if (retain_dense) {
        "dense"
      } else if (length(blocks) == 0L) {
        "diagonal"
      } else {
        "blocks"
      },
      K        = K,
      diagonal = diag(V_new),
      V        = if (retain_dense) V_new else NULL,
      blocks   = if (retain_dense) NULL else blocks,
      block_indices = if (retain_dense) indices else NULL,
      singular = .known_v_newdata_block_is_singular(V_new)
    )))
  }

  blocks   <- list()
  diagonal <- numeric(K)
  singular <- FALSE
  start    <- 1L
  for (input_block in V_new) {
    input_block <- .known_v_validate_newdata_block(input_block)
    input_index <- seq.int(start, length.out = nrow(input_block))
    diagonal[input_index] <- diag(input_block)
    local_indices <- .known_v_block_indices(input_block)
    for (local_index in local_indices[lengths(local_indices) > 1L]) {
      covariance <- input_block[local_index, local_index, drop = FALSE]
      blocks[[length(blocks) + 1L]] <- list(
        index      = input_index[local_index],
        covariance = covariance
      )
    }
    singular <- singular || .known_v_newdata_block_is_singular(input_block)
    start <- start + nrow(input_block)
  }

  .new_known_v(list(
    version  = 2L,
    storage  = "blocks",
    K        = K,
    diagonal = diagonal,
    blocks   = blocks,
    singular = singular
  ))
}


.known_v_validate_newdata_block <- function(V_new) {

  if (!is.numeric(V_new)) {
    stop("'V_new' must be numeric.", call. = FALSE)
  }
  if (nrow(V_new) == 0L || nrow(V_new) != ncol(V_new)) {
    stop("'V_new' must be a non-empty square matrix.", call. = FALSE)
  }
  if (anyNA(V_new) || any(!is.finite(V_new))) {
    stop("'V_new' must contain only finite non-missing values.", call. = FALSE)
  }
  .known_v_check_symmetric(V_new, "'V_new'")
  if (any(diag(V_new) < 0)) {
    stop("The diagonal of 'V_new' must contain non-negative variances.",
         call. = FALSE)
  }
  zero_variance <- diag(V_new) == 0
  if (any(zero_variance) &&
      any(V_new[zero_variance, , drop = FALSE] != 0)) {
    stop("'V_new' must be positive semidefinite.", call. = FALSE)
  }
  factorization <- .known_v_correlation_factorization(V_new)
  if (!is.null(factorization) &&
      !.covariance_is_positive_semidefinite(factorization)) {
    stop("'V_new' must be positive semidefinite.", call. = FALSE)
  }

  V_new
}


.known_v_newdata_block_is_singular <- function(V_new) {

  any(diag(V_new) == 0) ||
    isTRUE(.known_v_correlation_factorization(V_new)[["singular"]])
}


# Factorize the positive-variance correlation block for scale-stable checks.
.known_v_correlation_factorization <- function(V) {

  positive_variance <- diag(V) > 0
  if (!any(positive_variance)) {
    return(NULL)
  }

  correlation <- stats::cov2cor(
    V[positive_variance, positive_variance, drop = FALSE]
  )
  return(.covariance_factorization(correlation))
}
