# Known-V input helpers -----

.known_v_as_matrix <- function(V, k = NULL, warn_singular = TRUE) {

  V_matrix <- .known_v_as_matrix_structure(V, k = k)

  if (anyNA(V_matrix) || any(!is.finite(V_matrix))) {
    stop("The 'V' argument must contain only finite non-missing values.", call. = FALSE)
  }

  .known_v_check_symmetric(V_matrix, "'V'")

  V_matrix <- (V_matrix + t(V_matrix)) / 2
  diagonal <- diag(V_matrix)
  if (any(diagonal <= 0)) {
    stop("The diagonal of 'V' must contain positive variances.", call. = FALSE)
  }

  correlation <- stats::cov2cor(V_matrix)
  eig         <- eigen(correlation, symmetric = TRUE)
  eigenvalues <- eig[["values"]]
  tolerance   <- .Machine$double.eps * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) <= tolerance) {
    semidefinite_tolerance <- sqrt(.Machine$double.eps) *
      max(1, max(abs(eigenvalues)))
    if (min(eigenvalues) >= -semidefinite_tolerance) {
      if (min(eigenvalues) < 0) {
        eigenvalues[eigenvalues < 0] <- 0
        factor <- eig[["vectors"]] %*%
          diag(sqrt(eigenvalues), nrow = length(eigenvalues))
        correlation <- stats::cov2cor(tcrossprod(factor))
        V_matrix <- correlation * tcrossprod(sqrt(diagonal))
        V_matrix <- (V_matrix + t(V_matrix)) / 2
      }
      if (isTRUE(warn_singular)) {
        .known_v_warn_singular()
      }
    } else {
      stop("The 'V' argument must be positive definite.", call. = FALSE)
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


# Canonical known-V accessors. Legacy fields remain readable for cached objects.
.known_v_nrow <- function(known_V) {

  if (!is.null(known_V[["K"]])) {
    return(as.integer(known_V[["K"]]))
  }
  if (!is.null(known_V[["V"]])) {
    return(nrow(known_V[["V"]]))
  }
  if (!is.null(known_V[["diagonal"]])) {
    return(length(known_V[["diagonal"]]))
  }

  stop("Known-V metadata does not record its dimension.", call. = FALSE)
}


.known_v_storage <- function(known_V) {

  if (!is.null(known_V[["storage"]])) {
    return(known_V[["storage"]])
  }
  if (!is.null(known_V[["V"]])) {
    return("dense")
  }

  stop("Known-V metadata does not record its storage representation.",
       call. = FALSE)
}


.known_v_diagonal <- function(known_V) {

  if (!is.null(known_V[["diagonal"]])) {
    return(as.numeric(known_V[["diagonal"]]))
  }
  if (!is.null(known_V[["V"]])) {
    return(diag(known_V[["V"]]))
  }

  stop("Known-V metadata does not contain its diagonal.", call. = FALSE)
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


.known_v_covariance_for_index <- function(known_V, index) {

  if (length(index) == 1L) {
    return(matrix(.known_v_diagonal(known_V)[index], 1L, 1L))
  }

  blocks <- .known_v_correlated_blocks(known_V)
  match_block <- vapply(
    blocks,
    function(block) identical(block[["index"]], index),
    logical(1)
  )
  if (sum(match_block) != 1L) {
    stop("Known-V dependency block was not found.", call. = FALSE)
  }
  blocks[[which(match_block)]][["covariance"]]
}


.known_v_latent_apply <- function(known_V, z_samples) {

  z_samples <- as.matrix(z_samples)
  K         <- .known_v_nrow(known_V)
  out       <- matrix(0, nrow = nrow(z_samples), ncol = K)

  if (!is.null(known_V[["latent_blocks"]])) {
    for (block in known_V[["latent_blocks"]]) {
      if (block[["rank"]] == 0L) {
        next
      }
      z_index <- seq.int(block[["z_start"]], block[["z_end"]])
      out[, block[["index"]]] <- z_samples[, z_index, drop = FALSE] %*%
        t(block[["B"]])
    }
    return(out)
  }
  if (!is.null(known_V[["B"]])) {
    return(z_samples %*% t(known_V[["B"]]))
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

  symmetry_error <- max(abs(V_matrix - t(V_matrix)))
  tolerance      <- sqrt(.Machine$double.eps) * max(1, max(abs(V_matrix)))

  if (symmetry_error > tolerance) {
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
  known_V[["parameterization"]]  <- "block_mvn"
  known_V[["effective_backend"]] <- if (
    length(.known_v_correlated_blocks(known_V)) == 0L
  ) "diagonal" else "block_mvn"
  known_V[["correlated"]]        <-
    length(.known_v_correlated_blocks(known_V)) > 0L
  known_V[["residual_variance"]] <- .known_v_diagonal(known_V)
  known_V[["residual_sei"]]      <- sqrt(pmax(.known_v_diagonal(known_V), 0))
  known_V[["rank"]]              <- 0L

  known_V
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
    return(list(
      version  = 2L,
      storage  = "diagonal",
      K        = K,
      diagonal = diagonal,
      blocks   = list(),
      singular = any(diagonal == 0)
    ))
  }

  if (storage == "dense") {
    V_new   <- .known_v_validate_newdata_block(V_new)
    indices <- .known_v_block_indices(V_new)
    blocks  <- lapply(indices[lengths(indices) > 1L], function(index) {
      list(index = index, covariance = V_new[index, index, drop = FALSE])
    })
    retain_dense <- length(indices) == 1L &&
      length(indices[[1L]]) == K && K > 1L
    return(list(
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
    ))
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

  list(
    version  = 2L,
    storage  = "blocks",
    K        = K,
    diagonal = diagonal,
    blocks   = blocks,
    singular = singular
  )
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
  V_new <- (V_new + t(V_new)) / 2
  if (any(diag(V_new) < 0)) {
    stop("The diagonal of 'V_new' must contain non-negative variances.",
         call. = FALSE)
  }
  eigenvalues <- eigen(V_new, symmetric = TRUE, only.values = TRUE)[["values"]]
  tolerance   <- sqrt(.Machine$double.eps) * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -tolerance) {
    stop("'V_new' must be positive semidefinite.", call. = FALSE)
  }

  V_new
}


.known_v_newdata_block_is_singular <- function(V_new) {

  values    <- eigen(V_new, symmetric = TRUE, only.values = TRUE)[["values"]]
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(values)))
  min(values) <= tolerance
}


.known_v_newdata_as_matrix <- function(V_new, k) {

  if (is.matrix(V_new)) {
    V_matrix <- V_new
  } else if (is.numeric(V_new) && is.null(dim(V_new)) && length(V_new) > 0L) {
    V_matrix <- diag(as.numeric(V_new), nrow = length(V_new), ncol = length(V_new))
  } else if (is.list(V_new) && length(V_new) > 0L &&
             all(vapply(V_new, is.matrix, logical(1)))) {
    V_matrix <- .known_v_blockdiag(V_new, arg = "V_new")
  } else {
    stop(
      "'V_new' must be a variance vector, a square matrix, or a non-empty list of square matrices.",
      call. = FALSE
    )
  }

  if (!is.numeric(V_matrix)) {
    stop("'V_new' must be numeric.", call. = FALSE)
  }
  if (length(dim(V_matrix)) != 2L || nrow(V_matrix) != ncol(V_matrix)) {
    stop("'V_new' must be square.", call. = FALSE)
  }
  if (nrow(V_matrix) != k) {
    stop(
      "The dimensions of 'V_new' must match the number of rows in 'newdata'.",
      call. = FALSE
    )
  }
  if (anyNA(V_matrix) || any(!is.finite(V_matrix))) {
    stop("'V_new' must contain only finite non-missing values.", call. = FALSE)
  }

  .known_v_check_symmetric(V_matrix, "'V_new'")

  V_matrix <- (V_matrix + t(V_matrix)) / 2
  if (any(diag(V_matrix) < 0)) {
    stop("The diagonal of 'V_new' must contain non-negative variances.",
         call. = FALSE)
  }

  eigenvalues <- eigen(V_matrix, symmetric = TRUE, only.values = TRUE)[["values"]]
  tolerance   <- sqrt(.Machine$double.eps) * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -tolerance) {
    stop("'V_new' must be positive semidefinite.", call. = FALSE)
  }

  return(V_matrix)
}

.known_v_prepare <- function(V, keep_rows, known_v_parameterization,
                             known_v_residual_fraction,
                             known_v_residual_fraction_specified = FALSE,
                             known_v_is_scale = FALSE,
                             warn_singular = TRUE) {

  known_v_parameterization <- match.arg(
    known_v_parameterization,
    c("auto", "latent", "whitened", "block_mvn")
  )
  known_v_requested_parameterization <- known_v_parameterization

  if (!is.logical(keep_rows) ||
      length(keep_rows) != .known_v_input_nrow(V)) {
    stop("Internal error: invalid known-V row selector.", call. = FALSE)
  }
  BayesTools::check_bool(known_v_is_scale, "known_v_is_scale")

  if (is.null(known_v_residual_fraction)) {
    known_v_residual_fraction <- 0.10
  }
  .known_v_check_residual_fraction(known_v_residual_fraction)

  V       <- .known_v_subset_input(V, keep_rows)
  known_V <- .known_v_canonicalize(V, warn_singular = warn_singular)
  covariance_blocks <- .known_v_correlated_blocks(known_V)
  block_indices     <- lapply(covariance_blocks, `[[`, "index")
  if (length(.known_v_independent_indices(known_V)) > 0L) {
    block_indices[[length(block_indices) + 1L]] <- 1L
  }
  correlated <- length(covariance_blocks) > 0L
  singular       <- isTRUE(known_V[["singular"]])

  if (known_v_parameterization == "whitened" && known_v_is_scale && correlated) {
    stop(
      "known_v_parameterization = 'whitened' is currently available only without scale regression.",
      call. = FALSE
    )
  }

  if (known_v_parameterization == "auto") {
    known_v_parameterization <- .known_v_auto_parameterization(
      block_indices         = block_indices,
      known_v_is_scale      = known_v_is_scale,
      known_v_is_singular   = singular
    )
  }
  known_v_residual_fraction_metadata <- if (
    known_v_parameterization == "latent" ||
      isTRUE(known_v_residual_fraction_specified)
  ) {
    known_v_residual_fraction
  } else {
    NULL
  }
  if (isTRUE(singular) && known_v_parameterization == "latent") {
    stop(
      "Singular all-correlated known-V matrices cannot use ",
      "known_v_parameterization = 'latent'. Use 'block_mvn', or 'whitened' ",
      "when no scale regression or row-varying marginalized variance is present.",
      call. = FALSE
    )
  }
  known_V[["parameterization"]]           <- known_v_parameterization
  known_V[["parameterization_requested"]] <- known_v_requested_parameterization
  known_V[["effective_backend"]] <- if (correlated) {
    known_v_parameterization
  } else {
    "diagonal"
  }
  known_V[["correlated"]]                  <- correlated
  known_V[["residual_fraction_requested"]] <- known_v_residual_fraction_metadata

  if (known_v_parameterization == "whitened") {
    .known_v_warn_unused_residual_fraction(
      known_v_parameterization,
      known_v_residual_fraction_specified &&
        !identical(known_v_requested_parameterization, "auto")
    )

    whitening <- .known_v_whiten_blocks(known_V)
    return(c(
      known_V,
      list(
        residual_variance = .known_v_diagonal(known_V),
        residual_sei      = sqrt(.known_v_diagonal(known_V)),
        rank              = 0L
      ),
      whitening
    ))
  }

  if (known_v_parameterization == "block_mvn") {
    .known_v_warn_unused_residual_fraction(
      known_v_parameterization,
      known_v_residual_fraction_specified &&
        !identical(known_v_requested_parameterization, "auto")
    )

    block_mvn <- .known_v_block_mvn_blocks(known_V)
    return(c(
      known_V,
      list(
        residual_variance = .known_v_diagonal(known_V),
        residual_sei      = sqrt(.known_v_diagonal(known_V)),
        rank              = 0L
      ),
      block_mvn
    ))
  }

  decomposition <- .known_v_decompose_blocks(
    known_V           = known_V,
    residual_fraction = known_v_residual_fraction
  )

  return(c(known_V, decomposition))
}


.known_v_canonicalize <- function(V, warn_singular = TRUE) {

  storage <- .known_v_input_storage(V)
  K       <- .known_v_input_nrow(V)
  if (K == 0L) {
    stop("The 'V' argument must be non-empty.", call. = FALSE)
  }

  if (storage == "diagonal") {
    diagonal <- as.numeric(V)
    if (anyNA(diagonal) || any(!is.finite(diagonal))) {
      stop("The 'V' argument must contain only finite non-missing values.", call. = FALSE)
    }
    if (any(diagonal <= 0)) {
      stop("The diagonal of 'V' must contain positive variances.", call. = FALSE)
    }
    return(list(
      version  = 2L,
      storage  = "diagonal",
      K        = K,
      diagonal = diagonal,
      blocks   = list(),
      singular = FALSE
    ))
  }

  if (storage == "dense") {
    V       <- .known_v_as_matrix(V, warn_singular = FALSE)
    indices <- .known_v_block_indices(V)
    blocks  <- lapply(indices[lengths(indices) > 1L], function(index) {
      list(index = index, covariance = V[index, index, drop = FALSE])
    })
    singular <- .known_v_is_singular(V)
    if (singular && isTRUE(warn_singular)) {
      .known_v_warn_singular()
    }
    retain_dense <- length(indices) == 1L &&
      length(indices[[1L]]) == K && K > 1L
    return(list(
      version  = 2L,
      storage  = if (retain_dense) {
        "dense"
      } else if (length(blocks) == 0L) {
        "diagonal"
      } else {
        "blocks"
      },
      K        = K,
      diagonal = diag(V),
      V        = if (retain_dense) V else NULL,
      blocks   = if (retain_dense) NULL else blocks,
      block_indices = if (retain_dense) indices else NULL,
      singular = singular
    ))
  }

  blocks   <- list()
  diagonal <- numeric(K)
  singular <- FALSE
  start    <- 1L
  for (input_block in V) {
    input_block <- .known_v_as_matrix(input_block, warn_singular = FALSE)
    input_index <- seq.int(start, length.out = nrow(input_block))
    diagonal[input_index] <- diag(input_block)
    local_indices <- .known_v_block_indices(input_block)
    for (local_index in local_indices[lengths(local_indices) > 1L]) {
      covariance <- input_block[local_index, local_index, drop = FALSE]
      blocks[[length(blocks) + 1L]] <- list(
        index      = input_index[local_index],
        covariance = covariance
      )
      singular <- singular || .known_v_is_singular(covariance)
    }
    start <- start + nrow(input_block)
  }
  if (singular && isTRUE(warn_singular)) {
    .known_v_warn_singular()
  }

  list(
    version  = 2L,
    storage  = "blocks",
    K        = K,
    diagonal = diagonal,
    blocks   = blocks,
    singular = singular
  )
}


.known_v_correlated_blocks <- function(known_V) {

  blocks <- known_V[["blocks"]]
  if (!is.null(blocks)) {
    return(blocks)
  }

  Filter(function(block) length(block[["index"]]) > 1L,
         .known_v_blocks(known_V))
}


.known_v_independent_indices <- function(known_V) {

  K      <- .known_v_nrow(known_V)
  blocks <- .known_v_correlated_blocks(known_V)
  if (length(blocks) == 0L) {
    return(seq_len(K))
  }

  setdiff(seq_len(K), unlist(lapply(blocks, `[[`, "index"), use.names = FALSE))
}

.known_v_sampling_factor <- function(V) {

  chol_factor <- tryCatch(
    chol(V),
    error = function(e) NULL
  )
  if (!is.null(chol_factor)) {
    return(chol_factor)
  }

  eig       <- eigen(V, symmetric = TRUE)
  values    <- eig[["values"]]
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(values)))
  if (min(values) < -tolerance) {
    stop("Known-V sampling covariance is not positive semidefinite.",
         call. = FALSE)
  }

  values[values < 0] <- 0
  return(
    eig[["vectors"]] %*%
      diag(sqrt(values), nrow = length(values)) %*%
      t(eig[["vectors"]])
  )
}

.known_v_is_singular <- function(V) {

  correlation <- stats::cov2cor(V)
  eigenvalues <- eigen(
    correlation,
    symmetric   = TRUE,
    only.values = TRUE
  )[["values"]]
  tolerance   <- .Machine$double.eps * max(1, max(abs(eigenvalues)))

  return(min(eigenvalues) <= tolerance)
}

.known_v_warn_singular <- function() {

  warning(
    "The 'V' argument is positive semidefinite, not positive definite, ",
    "because at least one dependency block has a rank-deficient correlation ",
    "structure. Acceptance is provisional until the fitted priors and ",
    "compiled random effects confirm that every singular direction is ",
    "regularized.",
    call.      = FALSE,
    immediate. = TRUE
  )
}


# Reject singular known-V blocks without integrated conditional variance.
.brma_mv_check_singular_v_regularization <- function(object) {

  data <- object[["data"]]
  if (!.is_data_known_v(data)) {
    return(invisible(TRUE))
  }

  known_V <- .data_known_v_data(data)
  if (!isTRUE(known_V[["singular"]]) &&
      !(!is.null(known_V[["V"]]) && .known_v_is_singular(known_V[["V"]]))) {
    return(invisible(TRUE))
  }

  K                 <- .known_v_nrow(known_V)
  regularized_rows  <- .brma_mv_regularized_variance_rows(object, K = K)
  covariance_blocks <- .known_v_correlated_blocks(known_V)
  invalid_blocks    <- Filter(function(block) {
    covariance <- block[["covariance"]]
    index      <- block[["index"]]
    .known_v_is_singular(covariance) &&
      !.known_v_nullspace_is_regularized(
        covariance       = covariance,
        regularized_rows = regularized_rows[index]
      )
  }, covariance_blocks)

  if (length(invalid_blocks) > 0L) {
    block_labels <- .known_v_block_labels(invalid_blocks)
    stop(
      "Singular 'V' dependency block(s) at retained rows ",
      paste(block_labels, collapse = ", "),
      " are not structurally regularized by integrated conditional variance. ",
      "Sampled random effects change the conditional mean only and do not ",
      "regularize 'V'. Specify strictly positive heterogeneity or a marginalized ",
      "random-effect variance that covers every singular direction, or supply ",
      "a positive-definite 'V'.",
      call. = FALSE
    )
  }

  fixed_variance <- .brma_mv_fixed_integrated_variance(object, K = K)
  if (identical(known_V[["parameterization"]], "block_mvn") &&
      !is.null(fixed_variance)) {
    invalid_numeric_blocks <- Filter(function(block) {
      index      <- block[["index"]]
      covariance <- block[["covariance"]] +
        diag(fixed_variance[index], nrow = length(index))
      !.known_v_is_numerically_positive_definite(covariance)
    }, covariance_blocks)
    if (length(invalid_numeric_blocks) > 0L) {
      block_labels <- .known_v_block_labels(invalid_numeric_blocks)
      stop(
        "Fixed integrated variance does not make singular 'V' dependency ",
        "block(s) at retained rows ", paste(block_labels, collapse = ", "),
        " numerically positive definite for the block-MVN backend. Increase ",
        "the fixed heterogeneity/random-effect SD or supply a ",
        "positive-definite 'V'.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


# Format retained-row labels for known-V dependency blocks.
.known_v_block_labels <- function(blocks) {

  vapply(blocks, function(block) {
    paste0("{", paste(block[["index"]], collapse = ", "), "}")
  }, character(1))
}


# Return fixed integrated scalar variance when it is known before fitting.
.brma_mv_fixed_integrated_variance <- function(object, K) {

  data               <- object[["data"]]
  fixed_scale_values <- if (.is_data_scale(data)) {
    .brma_mv_fixed_scale_values(object, K = K)
  } else {
    list()
  }

  if (.is_data_random(data)) {
    terms <- .data_marginalized_random_effects(data)
    if (length(terms) == 0L) {
      return(NULL)
    }

    design     <- .fitted_formula_design(object, "mu", required = TRUE)
    prior_list <- design[["prior_list"]]
    variance   <- numeric(K)
    for (term in terms) {
      term_variance <- .marginalized_random_term_fixed_variance(
        term               = term,
        data               = data,
        prior_list         = prior_list,
        fixed_scale_values = fixed_scale_values,
        K                  = K
      )
      if (is.null(term_variance)) {
        return(NULL)
      }
      variance <- variance + term_variance
    }
    return(variance)
  }

  if (.is_data_scale(data)) {
    scale_specs <- .data_scale_component_specs(data)
    if (length(scale_specs) != 1L) {
      return(NULL)
    }
    source <- scale_specs[[1L]][["source"]]
    values <- fixed_scale_values[[source]]
    if (is.null(values)) {
      return(NULL)
    }
    return(values^2)
  }

  tau_prior <- object[["priors"]][["outcome"]][["tau"]]
  if (is.null(tau_prior) || BayesTools::is.prior.mixture(tau_prior) ||
      !BayesTools::is.prior.point(tau_prior)) {
    return(NULL)
  }

  location <- tau_prior[["parameters"]][["location"]]
  if (length(location) == 1L) {
    location <- rep(location, K)
  }
  if (length(location) != K || anyNA(location) ||
      any(!is.finite(location)) || any(location < 0)) {
    return(NULL)
  }

  location^2
}


# Return fitted-row SD values for each fully point-fixed scale formula.
.brma_mv_fixed_scale_values <- function(object, K) {

  scale_specs <- .data_scale_component_specs(object[["data"]])
  values      <- lapply(scale_specs, function(scale_spec) {
    design <- .fitted_formula_design(
      object    = object,
      parameter = scale_spec[["parameter"]],
      required  = TRUE
    )
    log_scale <- .formula_design_fixed_values(design)
    if (is.null(log_scale) || length(log_scale) != K) {
      return(NULL)
    }

    scale <- exp(log_scale)
    if (length(scale) != K || anyNA(scale) || any(!is.finite(scale)) ||
        any(scale < 0)) {
      return(NULL)
    }

    as.numeric(scale)
  })
  names(values) <- vapply(scale_specs, `[[`, character(1), "source")

  values
}


# Evaluate a fully point-fixed formula design on its fitted rows.
.formula_design_fixed_values <- function(design) {

  model_matrix <- design[["model_matrix"]]
  prior_list   <- design[["prior_list"]]
  parameter    <- design[["parameter"]]
  if (!is.matrix(model_matrix) || anyNA(model_matrix) ||
      any(!is.finite(model_matrix)) || is.null(prior_list) ||
      is.null(parameter)) {
    return(NULL)
  }

  output         <- numeric(nrow(model_matrix))
  intercept_name <- paste0(parameter, "_intercept")
  if (intercept_name %in% names(prior_list)) {
    prior      <- prior_list[[intercept_name]]
    value      <- .prior_fixed_values(prior)
    multiplier <- .prior_fixed_multiplier(prior)
    if (is.null(value) || length(value) != 1L || is.null(multiplier)) {
      return(NULL)
    }
    if (isTRUE(design[["log_intercept"]])) {
      if (!is.finite(value) || value <= 0) {
        stop(
          "Point-fixed scale intercepts must be strictly positive.",
          call. = FALSE
        )
      }
      value <- log(value)
    }
    output <- output + multiplier * value
  }

  term_names <- setdiff(names(prior_list), intercept_name)
  for (term_name in term_names) {
    prior      <- prior_list[[term_name]]
    values     <- .prior_fixed_values(prior)
    multiplier <- .prior_fixed_multiplier(prior)
    model_term <- sub(paste0("^", parameter, "_"), "", term_name)
    term_index <- match(model_term, design[["model_terms"]])
    if (is.null(values) || is.null(multiplier) || is.na(term_index)) {
      return(NULL)
    }
    columns <- which(design[["assign"]] == (term_index - 1L))
    if (length(columns) == 0L) {
      return(NULL)
    }
    if (BayesTools::is.prior.factor(prior)) {
      if (length(values) == 1L) {
        values <- rep(values, length(columns))
      }
      if (length(values) != length(columns)) {
        return(NULL)
      }
      contribution <- drop(model_matrix[, columns, drop = FALSE] %*% values)
    } else {
      if (length(values) != 1L || length(columns) != 1L) {
        return(NULL)
      }
      contribution <- model_matrix[, columns] * values
    }
    output <- output + multiplier * contribution
  }

  as.numeric(output)
}


# Return a prior's fixed numeric multiplier, if it has one.
.prior_fixed_multiplier <- function(prior) {

  multiplier <- attr(prior, "multiply_by", exact = TRUE)
  if (is.null(multiplier)) {
    return(1)
  }
  if (!is.numeric(multiplier) || length(multiplier) != 1L ||
      is.na(multiplier) || !is.finite(multiplier)) {
    return(NULL)
  }

  as.numeric(multiplier)
}


# Return a marginalized term's row-wise variance when all sources are fixed.
.marginalized_random_term_fixed_variance <- function(
    term, data, prior_list, fixed_scale_values = list(), K) {

  parameter      <- term[["sd_parameter_names"]]
  parameter      <- parameter[!is.na(parameter) & nzchar(parameter)]
  has_allocation <- .marginalized_random_effect_has_allocation(term)
  if (!has_allocation && length(parameter) == 1L) {
    fixed_sd <- .prior_fixed_values(prior_list[[parameter]])
  } else if (has_allocation) {
    binding    <- term[["sd_binding"]]
    allocation <- binding[["allocations"]][[1L]]
    fixed_sd <- .random_sd_source_fixed_values(
      source             = allocation[["source"]],
      data               = data,
      prior_list         = prior_list,
      fixed_scale_values = fixed_scale_values,
      K                  = K
    )
    if (is.null(fixed_sd)) {
      return(NULL)
    }
    for (factor in .marginalized_random_effect_allocation_factors(term)) {
      weights <- .prior_fixed_values(prior_list[[factor[["weight_name"]]]])
      index   <- factor[["index"]]
      if (is.null(weights) || length(weights) < index) {
        return(NULL)
      }
      multiplier <- weights[[index]]
      if (identical(factor[["scale"]], "mean_variance")) {
        multiplier <- factor[["n_targets"]] * multiplier
      } else if (!identical(factor[["scale"]], "total_variance")) {
        return(NULL)
      }
      if (!is.finite(multiplier) || multiplier < 0) {
        return(NULL)
      }
      fixed_sd <- fixed_sd * sqrt(multiplier)
    }
  } else {
    binding <- term[["sd_binding"]]
    source  <- if (is.null(binding)) NULL else binding[["source"]]
    if (is.null(source) && !is.null(binding) &&
        length(binding[["sources_by_column"]]) == 1L) {
      source <- binding[["sources_by_column"]][[1L]]
    }
    fixed_sd <- .random_sd_source_fixed_values(
      source             = source,
      data               = data,
      prior_list         = prior_list,
      fixed_scale_values = fixed_scale_values,
      K                  = K
    )
  }

  if (is.null(fixed_sd) || !length(fixed_sd) %in% c(1L, K)) {
    return(NULL)
  }
  variance <- .marginalized_random_effect_variance_samples(
    term       = term,
    sd_samples = matrix(fixed_sd, nrow = 1L),
    K          = K
  )

  as.numeric(variance[1L, ])
}


# Return fixed SD-source values, or NULL for a sampled/formula source.
.random_sd_source_fixed_values <- function(
    source, data, prior_list, fixed_scale_values = list(), K) {

  if (is.null(source)) {
    return(NULL)
  }

  values <- source[["values"]]
  if (is.null(values)) {
    name <- source[["name"]]
    if (!is.null(name) && name %in% .data_scale_formula_sources(data)) {
      return(fixed_scale_values[[name]])
    }
    if (!is.null(name) && name %in% names(prior_list)) {
      values <- .prior_fixed_values(prior_list[[name]])
    }
  }
  if (is.null(values)) {
    nested_source <- source[["source"]]
    if (!is.null(nested_source) && !identical(nested_source, source)) {
      return(.random_sd_source_fixed_values(
        source             = nested_source,
        data               = data,
        prior_list         = prior_list,
        fixed_scale_values = fixed_scale_values,
        K                  = K
      ))
    }
    return(NULL)
  }

  shape <- source[["shape"]]
  if (identical(shape, "scalar") && length(values) != 1L) {
    return(NULL)
  }
  if (identical(shape, "row") && length(values) == 1L) {
    values <- rep(values, K)
  }
  if (!identical(shape, "scalar") && !identical(shape, "row") &&
      length(values) != 1L) {
    return(NULL)
  }
  if (!length(values) %in% c(1L, K) || anyNA(values) ||
      any(!is.finite(values)) || any(values < 0)) {
    return(NULL)
  }

  as.numeric(values)
}


# Return a point prior's fixed values, including identical point mixtures.
.prior_fixed_values <- function(prior) {

  if (is.null(prior)) {
    return(NULL)
  }
  if (BayesTools::is.prior.mixture(prior)) {
    values <- lapply(prior, .prior_fixed_values)
    if (any(vapply(values, is.null, logical(1))) ||
        !all(vapply(values[-1L], identical, logical(1), values[[1L]]))) {
      return(NULL)
    }
    return(values[[1L]])
  }
  if (!BayesTools::is.prior.point(prior)) {
    return(NULL)
  }

  values <- prior[["parameters"]][["location"]]
  if (!is.numeric(values) || length(values) == 0L || anyNA(values) ||
      any(!is.finite(values))) {
    return(NULL)
  }

  as.numeric(values)
}


# Check positive definiteness at the covariance's relative numerical scale.
.known_v_is_numerically_positive_definite <- function(covariance) {

  values <- eigen(
    (covariance + t(covariance)) / 2,
    symmetric   = TRUE,
    only.values = TRUE
  )[["values"]]
  scale <- max(abs(values))
  if (!is.finite(scale) || scale == 0) {
    return(FALSE)
  }

  tolerance <- .Machine$double.eps * max(1, nrow(covariance)) * scale
  min(values) > tolerance
}


.brma_mv_regularized_variance_rows <- function(object, K) {

  data <- object[["data"]]
  if (!.is_data_random(data)) {
    if (.is_data_scale(data)) {
      return(rep(TRUE, K))
    }

    tau_prior <- object[["priors"]][["outcome"]][["tau"]]
    return(rep(.prior_coordinate_is_structurally_positive(tau_prior), K))
  }

  terms <- .data_marginalized_random_effects(data)
  if (length(terms) == 0L) {
    return(rep(FALSE, K))
  }

  design     <- .fitted_formula_design(object, "mu", required = TRUE)
  prior_list <- design[["prior_list"]]
  support    <- rep(FALSE, K)
  for (term in terms) {
    support <- support | .marginalized_random_term_positive_rows(
      term       = term,
      data       = data,
      prior_list = prior_list,
      K          = K
    )
  }

  support
}


.marginalized_random_term_positive_rows <- function(term, data, prior_list, K) {

  binding       <- term[["sd_binding"]]
  has_allocation <- .marginalized_random_effect_has_allocation(term)
  parameter <- term[["sd_parameter_names"]]
  parameter <- parameter[!is.na(parameter) & nzchar(parameter)]
  if (!has_allocation && length(parameter) == 1L) {
    support <- rep(
      .prior_coordinate_is_structurally_positive(prior_list[[parameter]]),
      K
    )
  } else {
    source  <- if (is.null(binding)) NULL else binding[["source"]]
    if (is.null(source) && !is.null(binding) &&
        length(binding[["sources_by_column"]]) == 1L) {
      source <- binding[["sources_by_column"]][[1L]]
    }
    support <- .random_sd_source_positive_rows(
      source     = source,
      data       = data,
      prior_list = prior_list,
      K          = K
    )

    if (has_allocation) {
      factors <- .marginalized_random_effect_allocation_factors(term)
      for (factor in factors) {
        factor_prior <- prior_list[[factor[["weight_name"]]]]
        factor_positive <- .prior_coordinate_is_structurally_positive(
          prior = factor_prior,
          index = factor[["index"]]
        )
        support <- support & factor_positive
      }
    }
  }

  multiplier <- .marginalized_random_effect_row_multiplier(term, K = K)
  if (!is.null(multiplier)) {
    support <- support & multiplier > 0
  }

  support
}


.random_sd_source_positive_rows <- function(source, data, prior_list, K) {

  if (is.null(source)) {
    return(rep(FALSE, K))
  }

  values <- source[["values"]]
  if (!is.null(values)) {
    if (length(values) == 1L) {
      values <- rep(values, K)
    }
    if (length(values) != K || anyNA(values) || any(!is.finite(values))) {
      stop("Invalid fixed random-effect SD source metadata.", call. = FALSE)
    }
    return(values > 0)
  }

  name <- source[["name"]]
  if (!is.null(name) && name %in% .data_scale_formula_sources(data)) {
    return(rep(TRUE, K))
  }
  if (!is.null(name) && name %in% names(prior_list)) {
    return(rep(
      .prior_coordinate_is_structurally_positive(prior_list[[name]]),
      K
    ))
  }

  nested_source <- source[["source"]]
  if (!is.null(nested_source) && !identical(nested_source, source)) {
    return(.random_sd_source_positive_rows(
      source     = nested_source,
      data       = data,
      prior_list = prior_list,
      K          = K
    ))
  }

  stop(
    "Cannot verify singular-V regularization from random-effect SD source '",
    if (is.null(name)) "<unnamed>" else name,
    "'.",
    call. = FALSE
  )
}


.prior_coordinate_is_structurally_positive <- function(prior, index = 1L) {

  if (is.null(prior)) {
    return(FALSE)
  }
  if (BayesTools::is.prior.mixture(prior)) {
    return(all(vapply(
      prior,
      .prior_coordinate_is_structurally_positive,
      logical(1),
      index = index
    )))
  }
  if (!BayesTools::is.prior.point(prior)) {
    return(TRUE)
  }

  location <- prior[["parameters"]][["location"]]
  if (length(location) == 1L) {
    location <- rep(location, index)
  }

  length(location) >= index && is.finite(location[[index]]) &&
    location[[index]] > 0
}


.known_v_nullspace_is_regularized <- function(covariance, regularized_rows) {

  if (!is.logical(regularized_rows) ||
      length(regularized_rows) != nrow(covariance) ||
      anyNA(regularized_rows)) {
    stop("Internal error: invalid singular-V regularization rows.",
         call. = FALSE)
  }

  eig       <- eigen(covariance, symmetric = TRUE)
  tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(eig[["values"]])))
  null      <- eig[["vectors"]][, eig[["values"]] <= tolerance, drop = FALSE]
  if (ncol(null) == 0L) {
    return(TRUE)
  }
  if (!any(regularized_rows)) {
    return(FALSE)
  }

  qr(null[regularized_rows, , drop = FALSE], tol = sqrt(.Machine$double.eps))[["rank"]] ==
    ncol(null)
}

.known_v_auto_parameterization <- function(block_indices, known_v_is_scale,
                                           known_v_is_singular = FALSE,
                                           max_block_size = NULL) {

  if (!isTRUE(known_v_is_scale)) {
    return("whitened")
  }
  if (isTRUE(known_v_is_singular)) {
    return("block_mvn")
  }

  if (.known_v_block_mvn_auto_feasible(
    block_indices  = block_indices,
    max_block_size = max_block_size
  )) {
    return("block_mvn")
  }

  return("latent")
}

.known_v_block_mvn_auto_feasible <- function(block_indices,
                                             max_block_size = NULL) {

  max_block_size <- .known_v_block_mvn_max_block_size(max_block_size)
  largest_block  <- .known_v_largest_block_size(block_indices)

  is.infinite(max_block_size) || largest_block <= max_block_size
}

.known_v_block_mvn_max_block_size <- function(max_block_size = NULL) {

  if (is.null(max_block_size)) {
    max_block_size <- getOption(
      "RoBMA.known_v_block_mvn_max_block_size",
      128L
    )
  }

  valid <- is.numeric(max_block_size) &&
    length(max_block_size) == 1L &&
    !is.na(max_block_size) &&
    max_block_size > 0 &&
    (is.infinite(max_block_size) || max_block_size == floor(max_block_size))

  if (!valid) {
    stop(
      "'RoBMA.known_v_block_mvn_max_block_size' must be a single positive ",
      "integer or Inf.",
      call. = FALSE
    )
  }

  if (is.infinite(max_block_size)) {
    return(Inf)
  }

  return(as.integer(max_block_size))
}

.known_v_largest_block_size <- function(block_indices) {

  if (length(block_indices) == 0L) {
    return(0L)
  }

  max(vapply(block_indices, length, integer(1)))
}

.known_v_check_residual_fraction <- function(known_v_residual_fraction) {

  BayesTools::check_real(
    known_v_residual_fraction,
    "known_v_residual_fraction",
    check_length = 1,
    lower        = 0,
    upper        = 1,
    allow_bound  = FALSE,
    allow_NA     = FALSE
  )

  return(invisible(TRUE))
}

.known_v_warn_unused_residual_fraction <- function(known_v_parameterization,
                                                   known_v_residual_fraction_specified) {

  if (isTRUE(known_v_residual_fraction_specified)) {
    warning(
      "'known_v_residual_fraction' is only used with ",
      "'known_v_parameterization = \"latent\"' and was disregarded for ",
      "'known_v_parameterization = \"", known_v_parameterization, "\"'.",
      call.      = FALSE,
      immediate. = TRUE
    )
  }

  return(invisible(TRUE))
}

.known_v_block_mvn_blocks <- function(known_V) {

  covariance_blocks <- .known_v_correlated_blocks(known_V)
  diagnostics       <- vector("list", length(covariance_blocks))
  blocks            <- vector("list", length(covariance_blocks))

  for (b in seq_along(covariance_blocks)) {
    idx     <- covariance_blocks[[b]][["index"]]
    V_block <- covariance_blocks[[b]][["covariance"]]
    eig     <- eigen(V_block, symmetric = TRUE, only.values = TRUE)[["values"]]

    blocks[[b]] <- list(
      index   = idx,
      size    = length(idx),
      v_lower = V_block[lower.tri(V_block, diag = TRUE)]
    )

    diagnostics[[b]] <- data.frame(
      block          = b,
      block_size     = length(idx),
      rank           = length(idx),
      min_eigenvalue = min(eig),
      max_eigenvalue = max(eig),
      stringsAsFactors = FALSE
    )
  }

  diagnostics <- .known_v_bind_diagnostics(diagnostics)

  return(list(
    block_mvn_blocks = blocks,
    diagnostics      = diagnostics
  ))
}

.known_v_block_indices <- function(V) {

  K         <- nrow(V)
  adjacency <- V != 0
  diag(adjacency) <- TRUE

  seen   <- rep(FALSE, K)
  blocks <- list()

  for (i in seq_len(K)) {
    if (seen[[i]]) {
      next
    }

    queue <- i
    seen[[i]] <- TRUE
    block <- integer(0)

    while (length(queue) > 0L) {
      current <- queue[[1]]
      queue   <- queue[-1]
      block   <- c(block, current)

      neighbors <- which(adjacency[current, ] & !seen)
      if (length(neighbors) > 0L) {
        seen[neighbors] <- TRUE
        queue <- c(queue, neighbors)
      }
    }

    blocks[[length(blocks) + 1L]] <- sort(block)
  }

  return(blocks)
}

.known_v_check_reconstruction <- function(V, residual_variance, B) {

  reconstruction <- diag(residual_variance, nrow = nrow(V)) + tcrossprod(B)
  error          <- max(abs(reconstruction - V))
  tolerance      <- 1e-8 * max(1, max(abs(V)))

  if (error > tolerance) {
    stop(
      "Known-V decomposition did not reconstruct 'V' within numerical tolerance.",
      call. = FALSE
    )
  }

  return(error)
}

.known_v_whiten_blocks <- function(known_V) {

  covariance_blocks <- .known_v_correlated_blocks(known_V)
  whitening_blocks  <- vector("list", length(covariance_blocks))
  diagnostics       <- vector("list", length(covariance_blocks))

  for (b in seq_along(covariance_blocks)) {
    idx       <- covariance_blocks[[b]][["index"]]
    V_block   <- covariance_blocks[[b]][["covariance"]]
    eig       <- eigen(V_block, symmetric = TRUE)
    values    <- eig[["values"]]
    tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(values)))
    if (min(values) < -tolerance) {
      stop("Known-V whitening covariance is not positive semidefinite.",
           call. = FALSE)
    }
    values[values < 0] <- 0

    whitening_blocks[[b]] <- list(
      index    = idx,
      size     = length(idx),
      rotation = t(eig[["vectors"]]),
      variance = values
    )
    diagnostics[[b]] <- data.frame(
      block                   = b,
      block_size              = length(idx),
      rank                    = length(idx),
      min_whitening_variance  = min(values),
      max_whitening_variance  = max(values),
      stringsAsFactors        = FALSE
    )
  }

  diagnostics <- .known_v_bind_diagnostics(diagnostics)

  return(list(
    whitening_blocks = whitening_blocks,
    diagnostics      = diagnostics
  ))
}

.known_v_decompose_blocks <- function(known_V, residual_fraction) {

  covariance_blocks <- .known_v_correlated_blocks(known_V)
  residual_variance <- .known_v_diagonal(known_V)
  latent_blocks     <- vector("list", length(covariance_blocks))
  diagnostics       <- vector("list", length(covariance_blocks))
  reduction_warning <- FALSE
  rank_total        <- 0L

  for (b in seq_along(covariance_blocks)) {
    idx      <- covariance_blocks[[b]][["index"]]
    V_block  <- covariance_blocks[[b]][["covariance"]]
    decomp   <- .known_v_decompose_block(V_block, residual_fraction)

    residual_variance[idx] <- decomp[["residual_variance"]]
    rank_block <- ncol(decomp[["B"]])
    latent_blocks[[b]] <- list(
      index   = idx,
      size    = length(idx),
      B       = decomp[["B"]],
      rank    = rank_block,
      z_start = rank_total + 1L,
      z_end   = rank_total + rank_block
    )
    rank_total <- rank_total + rank_block

    diagnostics[[b]] <- data.frame(
      block                          = b,
      block_size                     = length(idx),
      requested_residual_fraction    = residual_fraction,
      effective_residual_fraction    = decomp[["effective_residual_fraction"]],
      rank                           = ncol(decomp[["B"]]),
      max_reconstruction_error       = decomp[["max_reconstruction_error"]],
      min_latent_eigenvalue          = decomp[["min_latent_eigenvalue"]],
      min_residual_variance_fraction = min(decomp[["residual_variance"]] / diag(V_block)),
      stringsAsFactors               = FALSE
    )
    reduction_warning <- reduction_warning || decomp[["reduced"]]
  }

  if (reduction_warning) {
    warning(
      "'known_v_residual_fraction' was reduced for at least one known-V block ",
      "to keep V - D positive semidefinite.",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  diagnostics <- .known_v_bind_diagnostics(diagnostics)

  return(list(
    residual_variance = residual_variance,
    residual_sei      = sqrt(residual_variance),
    latent_blocks     = latent_blocks,
    rank              = rank_total,
    diagnostics       = diagnostics
  ))
}


.known_v_bind_diagnostics <- function(diagnostics) {

  if (length(diagnostics) == 0L) {
    return(data.frame())
  }

  out <- do.call(rbind, diagnostics)
  rownames(out) <- NULL
  out
}

.known_v_decompose_block <- function(V_block, residual_fraction) {

  block_size <- nrow(V_block)
  diagonal   <- diag(V_block)

  if (block_size == 1L || all(V_block[row(V_block) != col(V_block)] == 0)) {
    return(list(
      residual_variance            = diagonal,
      B                            = matrix(numeric(0), nrow = block_size, ncol = 0L),
      effective_residual_fraction  = 1,
      max_reconstruction_error     = 0,
      min_latent_eigenvalue        = NA_real_,
      reduced                      = FALSE
    ))
  }

  correlation <- stats::cov2cor(V_block)
  lambda_min  <- min(eigen(correlation, symmetric = TRUE, only.values = TRUE)[["values"]])
  alpha_max   <- 0.99 * lambda_min

  if (alpha_max <= sqrt(.Machine$double.eps)) {
    stop(
      "A known-V block is too close to singular for a positive residual ",
      "D + BB' decomposition.",
      call. = FALSE
    )
  }

  alpha   <- min(residual_fraction, alpha_max)
  reduced <- alpha < residual_fraction

  residual_variance <- alpha * diagonal
  latent_covariance <- V_block - diag(residual_variance, nrow = block_size)
  latent_covariance <- (latent_covariance + t(latent_covariance)) / 2

  eig <- eigen(latent_covariance, symmetric = TRUE)
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(eig[["values"]]))
  keep      <- eig[["values"]] > 0

  if (any(eig[["values"]] < -tolerance)) {
    stop("Known-V decomposition failed; V - D is not positive semidefinite.",
         call. = FALSE)
  }

  if (any(keep)) {
    B <- eig[["vectors"]][, keep, drop = FALSE] %*%
      diag(sqrt(pmax(eig[["values"]][keep], 0)), nrow = sum(keep))
  } else {
    B <- matrix(numeric(0), nrow = block_size, ncol = 0L)
  }

  reconstruction <- diag(residual_variance, nrow = block_size) + tcrossprod(B)

  return(list(
    residual_variance            = residual_variance,
    B                            = B,
    effective_residual_fraction  = alpha,
    max_reconstruction_error     = max(abs(reconstruction - V_block)),
    min_latent_eigenvalue        = min(eig[["values"]]),
    reduced                      = reduced
  ))
}
