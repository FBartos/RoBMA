# Known-V input helpers -----

.known_v_as_matrix <- function(V, k = NULL) {

  if (is.matrix(V)) {
    V_matrix <- V
  } else if (is.numeric(V) && is.null(dim(V)) && length(V) > 0L) {
    V_matrix <- diag(as.numeric(V), nrow = length(V), ncol = length(V))
  } else if (is.list(V) && length(V) > 0L && all(vapply(V, is.matrix, logical(1)))) {
    V_matrix <- .known_v_blockdiag(V)
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
  if (anyNA(V_matrix) || any(!is.finite(V_matrix))) {
    stop("The 'V' argument must contain only finite non-missing values.", call. = FALSE)
  }

  symmetry_error <- max(abs(V_matrix - t(V_matrix)))
  if (symmetry_error > sqrt(.Machine$double.eps)) {
    stop("The 'V' argument must be symmetric.", call. = FALSE)
  }

  V_matrix <- (V_matrix + t(V_matrix)) / 2
  diagonal <- diag(V_matrix)
  if (any(diagonal <= 0)) {
    stop("The diagonal of 'V' must contain positive variances.", call. = FALSE)
  }

  correlation <- stats::cov2cor(V_matrix)
  eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)[["values"]]
  tolerance   <- .Machine$double.eps * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) <= tolerance) {
    stop("The 'V' argument must be positive definite.", call. = FALSE)
  }

  return(V_matrix)
}

.known_v_blockdiag <- function(blocks) {

  sizes <- vapply(blocks, nrow, integer(1))
  for (i in seq_along(blocks)) {
    if (nrow(blocks[[i]]) != ncol(blocks[[i]])) {
      stop("All matrices in the 'V' list must be square.", call. = FALSE)
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

  V <- .known_v_newdata_as_matrix(V_new, k = k)
  eig <- eigen(V, symmetric = TRUE)

  values <- eig[["values"]]
  values[values < 0] <- 0
  sampling_factor <- eig[["vectors"]] %*%
    diag(sqrt(values), nrow = length(values)) %*%
    t(eig[["vectors"]])

  block_indices <- .known_v_block_indices(V)

  return(list(
    V                 = V,
    block_indices     = block_indices,
    parameterization  = "block_mvn",
    correlated        = any(vapply(block_indices, length, integer(1)) > 1L),
    residual_variance = diag(V),
    residual_sei      = sqrt(pmax(diag(V), 0)),
    B                 = matrix(numeric(0), nrow = nrow(V), ncol = 0L),
    rank              = 0L,
    sampling_factor   = sampling_factor
  ))
}


.known_v_newdata_as_matrix <- function(V_new, k) {

  if (is.matrix(V_new)) {
    V_matrix <- V_new
  } else if (is.numeric(V_new) && is.null(dim(V_new)) && length(V_new) > 0L) {
    V_matrix <- diag(as.numeric(V_new), nrow = length(V_new), ncol = length(V_new))
  } else if (is.list(V_new) && length(V_new) > 0L &&
             all(vapply(V_new, is.matrix, logical(1)))) {
    V_matrix <- .known_v_blockdiag(V_new)
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

  symmetry_error <- max(abs(V_matrix - t(V_matrix)))
  if (symmetry_error > sqrt(.Machine$double.eps)) {
    stop("'V_new' must be symmetric.", call. = FALSE)
  }

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
                             known_v_is_multilevel = FALSE) {

  known_v_parameterization <- match.arg(
    known_v_parameterization,
    c("auto", "latent", "whitened", "block_mvn")
  )
  known_v_requested_parameterization <- known_v_parameterization

  if (!is.logical(keep_rows) || length(keep_rows) != nrow(V)) {
    stop("Internal error: invalid known-V row selector.", call. = FALSE)
  }
  BayesTools::check_bool(known_v_is_scale, "known_v_is_scale")
  BayesTools::check_bool(known_v_is_multilevel, "known_v_is_multilevel")

  if (is.null(known_v_residual_fraction)) {
    known_v_residual_fraction <- 0.10
  }
  .known_v_check_residual_fraction(known_v_residual_fraction)

  V <- V[keep_rows, keep_rows, drop = FALSE]
  V <- .known_v_as_matrix(V)

  block_indices <- .known_v_block_indices(V)
  correlated    <- any(vapply(block_indices, length, integer(1)) > 1L)

  if (known_v_parameterization == "auto") {
    known_v_parameterization <- .known_v_auto_parameterization(
      block_indices         = block_indices,
      known_v_is_scale      = known_v_is_scale,
      known_v_is_multilevel = known_v_is_multilevel
    )
  }

  if (known_v_parameterization == "whitened") {
    .known_v_warn_unused_residual_fraction(
      known_v_parameterization,
      known_v_residual_fraction_specified
    )

    whitening <- .known_v_whiten_blocks(
      V             = V,
      block_indices = block_indices
    )

    return(c(
      list(
        V                           = V,
        block_indices               = block_indices,
        parameterization            = known_v_parameterization,
        parameterization_requested  = known_v_requested_parameterization,
        correlated                  = correlated,
        residual_variance           = diag(V),
        residual_sei                = sqrt(diag(V)),
        B                           = matrix(numeric(0), nrow = nrow(V), ncol = 0L),
        rank                        = 0L,
        sampling_factor             = chol(V)
      ),
      whitening
    ))
  }

  if (known_v_parameterization == "block_mvn") {
    .known_v_warn_unused_residual_fraction(
      known_v_parameterization,
      known_v_residual_fraction_specified
    )

    block_mvn <- .known_v_block_mvn_blocks(
      V             = V,
      block_indices = block_indices
    )

    return(c(
      list(
        V                           = V,
        block_indices               = block_indices,
        parameterization            = known_v_parameterization,
        parameterization_requested  = known_v_requested_parameterization,
        correlated                  = correlated,
        residual_variance           = diag(V),
        residual_sei                = sqrt(diag(V)),
        B                           = matrix(numeric(0), nrow = nrow(V), ncol = 0L),
        rank                        = 0L,
        sampling_factor             = chol(V)
      ),
      block_mvn
    ))
  }

  decomposition <- .known_v_decompose_blocks(
    V                 = V,
    block_indices     = block_indices,
    residual_fraction = known_v_residual_fraction
  )
  reconstruction_error <- .known_v_check_reconstruction(
    V                   = V,
    residual_variance   = decomposition[["residual_variance"]],
    B                   = decomposition[["B"]]
  )

  return(c(
    list(
      V                           = V,
      block_indices               = block_indices,
      parameterization            = known_v_parameterization,
      parameterization_requested  = known_v_requested_parameterization,
      correlated                  = correlated,
      residual_fraction_requested = known_v_residual_fraction,
      max_reconstruction_error    = reconstruction_error,
      sampling_factor             = chol(V)
    ),
    decomposition
  ))
}

.known_v_auto_parameterization <- function(block_indices, known_v_is_scale,
                                           known_v_is_multilevel,
                                           max_block_size = NULL) {

  if (!isTRUE(known_v_is_scale) && !isTRUE(known_v_is_multilevel)) {
    return("whitened")
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

.known_v_block_mvn_blocks <- function(V, block_indices) {

  diagnostics <- vector("list", length(block_indices))
  blocks      <- vector("list", length(block_indices))

  for (b in seq_along(block_indices)) {
    idx     <- block_indices[[b]]
    V_block <- V[idx, idx, drop = FALSE]
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

  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL

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

.known_v_whiten_blocks <- function(V, block_indices) {

  K                  <- nrow(V)
  whitening_matrix   <- matrix(0, nrow = K, ncol = K)
  whitening_variance <- numeric(K)
  diagnostics        <- vector("list", length(block_indices))

  for (b in seq_along(block_indices)) {
    idx     <- block_indices[[b]]
    V_block <- V[idx, idx, drop = FALSE]
    eig     <- eigen(V_block, symmetric = TRUE)

    whitening_matrix[idx, idx]   <- t(eig[["vectors"]])
    whitening_variance[idx]      <- eig[["values"]]
    diagnostics[[b]] <- data.frame(
      block                   = b,
      block_size              = length(idx),
      rank                    = length(idx),
      min_whitening_variance  = min(eig[["values"]]),
      max_whitening_variance  = max(eig[["values"]]),
      stringsAsFactors        = FALSE
    )
  }

  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL

  return(list(
    whitening_matrix   = whitening_matrix,
    whitening_variance = whitening_variance,
    diagnostics        = diagnostics
  ))
}

.known_v_decompose_blocks <- function(V, block_indices, residual_fraction) {

  K                 <- nrow(V)
  residual_variance <- numeric(K)
  B                 <- matrix(numeric(0), nrow = K, ncol = 0L)
  diagnostics       <- vector("list", length(block_indices))
  reduction_warning <- FALSE

  for (b in seq_along(block_indices)) {
    idx      <- block_indices[[b]]
    V_block  <- V[idx, idx, drop = FALSE]
    decomp   <- .known_v_decompose_block(V_block, residual_fraction)

    residual_variance[idx] <- decomp[["residual_variance"]]
    if (ncol(decomp[["B"]]) > 0L) {
      B_block <- matrix(0, nrow = K, ncol = ncol(decomp[["B"]]))
      B_block[idx, ] <- decomp[["B"]]
      B <- cbind(B, B_block)
    }

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

  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL

  return(list(
    residual_variance = residual_variance,
    residual_sei      = sqrt(residual_variance),
    B                 = B,
    rank              = ncol(B),
    diagnostics       = diagnostics
  ))
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
