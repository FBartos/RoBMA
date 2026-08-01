# Known-V backend preparation -----

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
  singular       <- .known_v_is_singular_representation(known_V)

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
  exact_rank_one_latent <- isTRUE(singular) &&
    .known_v_singular_blocks_are_exact_rank_one(known_V)
  if (isTRUE(singular) && known_v_parameterization == "latent" &&
      !exact_rank_one_latent) {
    stop(
      "Singular all-correlated known-V matrices cannot use ",
      "known_v_parameterization = 'latent'. Use 'block_mvn', or 'whitened' ",
      "when no scale regression or row-varying marginalized variance is present.",
      call. = FALSE
    )
  }
  effective_backend <- if (!correlated) {
    "diagonal"
  } else if (exact_rank_one_latent) {
    "latent"
  } else {
    known_v_parameterization
  }
  known_V <- .known_v_update(known_V, list(
    parameterization           = known_v_parameterization,
    parameterization_requested = known_v_requested_parameterization,
    effective_backend          = effective_backend,
    correlated                 = correlated,
    residual_fraction_requested = known_v_residual_fraction_metadata
  ))

  if (effective_backend == "whitened") {
    .known_v_warn_unused_residual_fraction(
      known_v_parameterization,
      known_v_residual_fraction_specified &&
        !identical(known_v_requested_parameterization, "auto")
    )

    whitening <- .known_v_whiten_blocks(known_V)
    return(.known_v_update(
      known_V,
      c(list(
        residual_variance = .known_v_diagonal(known_V),
        residual_sei      = sqrt(.known_v_diagonal(known_V)),
        rank              = 0L
      ), whitening)
    ))
  }

  if (effective_backend == "block_mvn") {
    .known_v_warn_unused_residual_fraction(
      known_v_parameterization,
      known_v_residual_fraction_specified &&
        !identical(known_v_requested_parameterization, "auto")
    )

    block_mvn <- .known_v_block_mvn_blocks(known_V)
    return(.known_v_update(
      known_V,
      c(list(
        residual_variance = .known_v_diagonal(known_V),
        residual_sei      = sqrt(.known_v_diagonal(known_V)),
        rank              = 0L
      ), block_mvn)
    ))
  }

  decomposition <- .known_v_decompose_blocks(
    known_V           = known_V,
    residual_fraction = known_v_residual_fraction
  )

  return(.known_v_update(known_V, decomposition))
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
    return(.new_known_v(list(
      version  = 2L,
      storage  = "diagonal",
      K        = K,
      diagonal = diagonal,
      blocks   = list(),
      singular = FALSE
    )))
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
      diagonal = diag(V),
      V        = if (retain_dense) V else NULL,
      blocks   = if (retain_dense) NULL else blocks,
      block_indices = if (retain_dense) indices else NULL,
      singular = singular
    )))
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

  .new_known_v(list(
    version  = 2L,
    storage  = "blocks",
    K        = K,
    diagonal = diagonal,
    blocks   = blocks,
    singular = singular
  ))
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

  factorization <- .covariance_factorization(V)
  factor        <- .covariance_sampling_factor(factorization)
  if (is.null(factor)) {
    stop("Known-V sampling covariance is not positive semidefinite.",
         call. = FALSE)
  }

  factor
}

.known_v_is_singular <- function(V) {

  .known_v_covariance_classification(V)[["singular"]]
}


# Whether every singular dependency block has an exact rank-one factor.
.known_v_singular_blocks_are_exact_rank_one <- function(known_V) {

  singular_blocks <- Filter(function(block) {
    .known_v_is_singular(block[["covariance"]])
  }, .known_v_correlated_blocks(known_V))

  length(singular_blocks) > 0L && all(vapply(
    singular_blocks,
    function(block) {
      !is.null(.covariance_exact_rank_one_factor(block[["covariance"]]))
    },
    logical(1)
  ))
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
    eig     <- .covariance_factorization(V_block)[["eigenvalues"]]

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

.known_v_whiten_blocks <- function(known_V) {

  covariance_blocks <- .known_v_correlated_blocks(known_V)
  whitening_blocks  <- vector("list", length(covariance_blocks))
  diagnostics       <- vector("list", length(covariance_blocks))

  for (b in seq_along(covariance_blocks)) {
    idx       <- covariance_blocks[[b]][["index"]]
    V_block   <- covariance_blocks[[b]][["covariance"]]
    factorization <- .covariance_factorization(V_block)
    values        <- factorization[["decomposition_values"]]
    if (!.covariance_is_positive_semidefinite(factorization)) {
      stop("Known-V whitening covariance is not positive semidefinite.",
           call. = FALSE)
    }
    if (any(values < 0)) {
      stop("Known-V whitening covariance has negative eigenvalues.",
           call. = FALSE)
    }

    whitening_blocks[[b]] <- list(
      index    = idx,
      size     = length(idx),
      rotation = t(factorization[["eigenvectors"]]),
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

  rank_one_factor <- .covariance_exact_rank_one_factor(V_block)
  if (!is.null(rank_one_factor) && block_size > 1L) {
    return(list(
      residual_variance            = numeric(block_size),
      B                            = matrix(rank_one_factor, ncol = 1L),
      effective_residual_fraction  = 0,
      max_reconstruction_error     = 0,
      min_latent_eigenvalue        = 0,
      reduced                      = FALSE
    ))
  }

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
  lambda_min  <- min(.covariance_factorization(correlation)[["eigenvalues"]])
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
  eig  <- .covariance_factorization(latent_covariance)
  keep <- eig[["decomposition_values"]] > 0

  if (!.covariance_is_positive_semidefinite(eig)) {
    stop("Known-V decomposition failed; V - D is not positive semidefinite.",
         call. = FALSE)
  }
  if (any(eig[["decomposition_values"]] < 0)) {
    stop("Known-V decomposition produced negative eigenvalues.", call. = FALSE)
  }

  if (any(keep)) {
    B <- eig[["eigenvectors"]][, keep, drop = FALSE] %*%
      diag(sqrt(eig[["decomposition_values"]][keep]), nrow = sum(keep))
  } else {
    B <- matrix(numeric(0), nrow = block_size, ncol = 0L)
  }

  reconstruction <- diag(residual_variance, nrow = block_size) + tcrossprod(B)

  return(list(
    residual_variance            = residual_variance,
    B                            = B,
    effective_residual_fraction  = alpha,
    max_reconstruction_error     = max(abs(reconstruction - V_block)),
    min_latent_eigenvalue        = min(eig[["decomposition_values"]]),
    reduced                      = reduced
  ))
}
