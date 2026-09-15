# Known-V backend preparation -----

.selection_bind_groups <- function(model, row_index, input_data,
                                   cluster = NULL,
                                   allow_singletons = FALSE) {

  if (!inherits(model, "selection_model") || !is.list(model)) {
    stop("Internal error: a validated selection model is required for group binding.", call. = FALSE)
  }
  if (!is.numeric(row_index) || !length(row_index) || anyNA(row_index) ||
      any(!is.finite(row_index)) || any(row_index < 1L) ||
      any(row_index != as.integer(row_index)) || anyDuplicated(row_index)) {
    stop("Internal error: the selection row map is invalid.", call. = FALSE)
  }
  row_index <- as.integer(row_index)
  BayesTools::check_bool(allow_singletons, "allow_singletons")
  requested <- model[["group"]]
  if (identical(model[["weight_rule"]], "product")) {
    # Native selection interfaces still need a row partition. These singleton
    # identities carry no publication meaning and never determine dependencies.
    requested <- NULL
    values <- row_index
    provenance <- "inactive"
  } else if (!is.null(requested)) {
    if (!is.character(requested) || length(requested) != 1L ||
        is.na(requested) || !nzchar(requested)) {
      stop("The selection model 'group' must be a column name or NULL.", call. = FALSE)
    }
    if (!is.data.frame(input_data) || !requested %in% names(input_data)) {
      stop("The 'group' column '", requested, "' was not found in 'data'.", call. = FALSE)
    }
    if (max(row_index) > nrow(input_data)) {
      stop("The 'group' column does not match the original data row map.", call. = FALSE)
    }
    values <- input_data[[requested]][row_index]
    provenance <- "explicit"
  } else if (!is.null(cluster)) {
    if (max(row_index) > length(cluster)) {
      stop("The 'cluster' identifiers do not match the original data row map.", call. = FALSE)
    }
    values <- cluster[row_index]
    provenance <- "cluster"
  } else if (allow_singletons) {
    values <- row_index
    provenance <- "singleton"
  } else {
    stop(
      "Publication groups are unavailable for this input. Specify 'group' in 'selection_model()'.",
      call. = FALSE
    )
  }
  if (!is.atomic(values) || !is.null(dim(values)) ||
      length(values) != length(row_index)) {
    stop("Publication group identifiers must have one value per retained data row.", call. = FALSE)
  }
  if (anyNA(values)) {
    stop("Publication group identifiers must not be missing among retained data rows.", call. = FALSE)
  }
  labels <- unique(values)
  group_index <- match(values, labels)
  list(
    version     = 1L,
    requested   = requested,
    row_index   = row_index,
    row_labels  = values,
    group_index = group_index,
    row_blocks  = unname(split(seq_along(group_index), group_index)),
    labels      = labels,
    provenance  = provenance
  )
}

.known_v_prepare <- function(V, keep_rows, known_v_parameterization,
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

  metadata <- .known_v_input_metadata(V)
  metadata <- .known_v_subset_selection_metadata(metadata, which(keep_rows))
  V       <- .known_v_subset_input(V, keep_rows)
  known_V <- .known_v_canonicalize(V, warn_singular = warn_singular)
  known_V <- .known_v_update(known_V, list(selection_metadata = metadata))
  declared_factor <- identical(.known_v_storage(known_V), "factor")
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
  residual_fraction_requested <- if (declared_factor) {
    NULL
  } else if (known_v_parameterization == "latent") {
    0.10
  } else {
    NULL
  }
  exact_declared_latent <- declared_factor ||
    (isTRUE(singular) &&
       .known_v_singular_blocks_are_exact_rank_one(known_V))
  if (isTRUE(singular) && known_v_parameterization == "latent" &&
      !exact_declared_latent) {
    stop(
      "Singular all-correlated known-V matrices cannot use ",
      "known_v_parameterization = 'latent'. Use 'block_mvn', or 'whitened' ",
      "when no scale regression or row-varying marginalized variance is present.",
      call. = FALSE
    )
  }
  effective_backend <- if (!correlated) {
    "diagonal"
  } else if (isTRUE(singular) && exact_declared_latent) {
    "latent"
  } else {
    known_v_parameterization
  }
  known_V <- .known_v_update(known_V, list(
    parameterization           = known_v_parameterization,
    parameterization_requested = known_v_requested_parameterization,
    effective_backend          = effective_backend,
    correlated                 = correlated,
    residual_fraction_requested = residual_fraction_requested
  ))

  if (effective_backend == "whitened") {

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

  if (effective_backend == "diagonal") {
    return(.known_v_update(known_V, list(
      residual_variance = .known_v_diagonal(known_V),
      residual_sei      = sqrt(.known_v_diagonal(known_V)),
      latent_blocks     = list(),
      rank              = 0L,
      diagnostics       = data.frame()
    )))
  }

  decomposition <- if (declared_factor) {
    .known_v_decompose_declared_factor(known_V)
  } else {
    .known_v_decompose_blocks(known_V)
  }

  return(.known_v_update(known_V, decomposition))
}


.known_v_canonicalize <- function(V, warn_singular = TRUE) {

  storage <- .known_v_input_storage(V)
  K       <- .known_v_input_nrow(V)
  metadata <- .known_v_input_metadata(V)
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
      selection_metadata = metadata,
      storage  = "diagonal",
      K        = K,
      diagonal = diagonal,
      blocks   = list(),
      singular = FALSE
    )))
  }

  if (storage == "factor") {
    components <- .known_v_factor_components(V)
    diagonal   <- components[["diagonal"]] +
      rowSums(components[["loading"]]^2)
    if (any(diagonal <= 0)) {
      stop("The diagonal of 'V' must contain positive variances.",
           call. = FALSE)
    }
    covariance <- .known_v_factor_covariance(
      components[["diagonal"]],
      components[["loading"]]
    )
    if (any(!is.finite(covariance))) {
      stop("The 'V' argument must contain only finite non-missing values.",
           call. = FALSE)
    }
    # Statistical dependence belongs to V, not to overlapping auxiliary
    # loading columns whose cross-products can cancel exactly.
    block_indices <- .known_v_block_indices(covariance)
    singular <- .known_v_is_singular(covariance)
    if (singular && isTRUE(warn_singular)) {
      .known_v_warn_singular()
    }
    return(.new_known_v(list(
      version         = 2L,
      selection_metadata = metadata,
      storage         = "factor",
      K               = K,
      diagonal        = diagonal,
      factor_diagonal = components[["diagonal"]],
      factor_loading  = components[["loading"]],
      blocks          = NULL,
      block_indices   = block_indices,
      singular        = singular
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
    return(.known_v_attach_certified_factor(.new_known_v(list(
      version  = 2L,
      selection_metadata = metadata,
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
    ))))
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

  .known_v_attach_certified_factor(.new_known_v(list(
    version  = 2L,
    selection_metadata = metadata,
    storage  = "blocks",
    K        = K,
    diagonal = diagonal,
    blocks   = blocks,
    singular = singular
  )))
}


# Attach a certified block-constant factor when every correlated dependency
# block recovers exactly. The supplied entries stay authoritative; the factor
# is a derived computational representation of the same matrix, used only
# where a declared factor is already consumed.
#
# Recovery is per block: a block that declines keeps the supplied entries and
# its dense route, and the rest of the matrix still reaches the factor rules.
.known_v_attach_certified_factor <- function(known_V) {

  blocks <- .known_v_blocks(known_V)
  correlated <- which(vapply(blocks, function(block) {
    length(block[["index"]]) > 1L
  }, logical(1)))
  if (length(correlated) == 0L) {
    return(known_V)
  }

  diagonal   <- .known_v_diagonal(known_V)
  certified  <- list()
  dense_rows <- integer(0)
  rank       <- 0L
  for (position in correlated) {
    index  <- blocks[[position]][["index"]]
    factor <- .covariance_block_constant_factor(
      blocks[[position]][["covariance"]]
    )
    if (is.null(factor)) {
      dense_rows <- c(dense_rows, index)
      next
    }
    diagonal[index] <- factor[["diagonal"]]
    rank <- rank + factor[["rank"]]
    certified[[length(certified) + 1L]] <- list(
      index    = index,
      loading  = factor[["loading"]],
      levels   = factor[["levels"]],
      supports = lapply(factor[["supports"]], function(rows) index[rows]),
      support  = factor[["support"]],
      depth    = factor[["depth"]],
      residual = factor[["residual"]]
    )
  }
  if (length(certified) == 0L) {
    return(known_V)
  }

  .known_v_update(known_V, list(certified_factor = list(
    status     = "recovered_block_constant",
    diagonal   = diagonal,
    blocks     = certified,
    dense_rows = sort(dense_rows),
    rank       = rank
  )))
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
    "structure.",
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
    values        <- factorization[["spectral_values"]]
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

.known_v_decompose_blocks <- function(known_V) {

  covariance_blocks <- .known_v_correlated_blocks(known_V)
  residual_variance <- .known_v_diagonal(known_V)
  latent_blocks     <- vector("list", length(covariance_blocks))
  diagnostics       <- vector("list", length(covariance_blocks))
  rank_total        <- 0L

  for (b in seq_along(covariance_blocks)) {
    idx      <- covariance_blocks[[b]][["index"]]
    V_block  <- covariance_blocks[[b]][["covariance"]]
    decomp   <- .known_v_decompose_block(V_block)

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
      requested_residual_fraction    = 0.10,
      effective_residual_fraction    = decomp[["effective_residual_fraction"]],
      rank                           = ncol(decomp[["B"]]),
      max_reconstruction_error       = decomp[["max_reconstruction_error"]],
      min_latent_eigenvalue          = decomp[["min_latent_eigenvalue"]],
      min_residual_variance_fraction = min(decomp[["residual_variance"]] / diag(V_block)),
      stringsAsFactors               = FALSE
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


.known_v_decompose_declared_factor <- function(known_V) {

  diagonal         <- known_V[["factor_diagonal"]]
  loading          <- known_V[["factor_loading"]]
  source_ids       <- .known_v_selection_metadata(known_V)[["source_ids"]]
  total_diagonal   <- .known_v_diagonal(known_V)
  covariance_blocks <- Filter(function(block) {
    any(loading[block[["index"]], , drop = FALSE] != 0)
  }, .known_v_blocks(known_V))
  latent_blocks    <- vector("list", length(covariance_blocks))
  diagnostics      <- vector("list", length(covariance_blocks))
  rank_total       <- 0L

  for (b in seq_along(covariance_blocks)) {
    index <- covariance_blocks[[b]][["index"]]
    support_size <- if (ncol(loading) == 0L) {
      integer()
    } else {
      colSums(loading[index, , drop = FALSE] != 0)
    }
    factor_columns <- which(support_size > 0L)
    B          <- loading[index, factor_columns, drop = FALSE]
    rank_block <- ncol(B)
    latent_blocks[[b]] <- list(
      index   = index,
      size    = length(index),
      B       = B,
      source_ids = source_ids[factor_columns],
      rank    = rank_block,
      z_start = rank_total + 1L,
      z_end   = rank_total + rank_block
    )
    rank_total <- rank_total + rank_block
    diagnostics[[b]] <- data.frame(
      block                          = b,
      block_size                     = length(index),
      requested_residual_fraction    = NA_real_,
      effective_residual_fraction    = NA_real_,
      rank                           = rank_block,
      max_reconstruction_error       = 0,
      min_latent_eigenvalue          = NA_real_,
      min_residual_variance_fraction = min(
        diagonal[index] / total_diagonal[index]
      ),
      stringsAsFactors = FALSE
    )
  }

  list(
    residual_variance = diagonal,
    residual_sei      = sqrt(diagonal),
    latent_blocks     = latent_blocks,
    rank              = rank_total,
    diagnostics       = .known_v_bind_diagnostics(diagnostics)
  )
}


.known_v_resolve_selection_structure <- function(known_V) {

  if (!is.null(known_V[["selection_structure"]])) {
    .known_v_check_selection_structure(known_V)
    return(known_V)
  }
  metadata   <- .known_v_selection_metadata(known_V)
  K          <- .known_v_nrow(known_V)
  blocks     <- .known_v_blocks(known_V)
  rank_total <- 0L
  for (b in seq_along(blocks)) {
    rows <- blocks[[b]][["index"]]
    if (identical(.known_v_storage(known_V), "factor")) {
      diagonal         <- known_V[["factor_diagonal"]][rows]
      loading          <- known_V[["factor_loading"]][rows, , drop = FALSE]
      loading          <- loading[, colSums(abs(loading)) > 0, drop = FALSE]
      independent      <- which(diagonal > 0)
      diagonal_loading <- matrix(0, length(rows), length(independent))
      if (length(independent)) {
        diagonal_loading[cbind(independent, seq_along(independent))] <-
          sqrt(diagonal[independent])
      }
      B <- cbind(diagonal_loading, loading)
    } else {
      covariance <- blocks[[b]][["covariance"]]
      if (all(covariance == 0)) {
        B <- matrix(0, length(rows), 0L)
      } else {
        root <- tryCatch(chol(covariance), error = function(e) NULL)
        if (is.null(root)) {
          # Pivoted Cholesky handles an already validated PSD covariance.
          # tol = 0 preserves positive pivots rather than dropping small SDs.
          root <- suppressWarnings(chol(covariance, pivot = TRUE, tol = 0))
          rank <- attr(root, "rank")
          B    <- t(root[seq_len(rank), order(attr(root, "pivot")), drop = FALSE])
        } else {
          B <- t(root)
        }
      }
    }
    rank <- ncol(B)
    blocks[[b]] <- list(
      index = rows, size = length(rows), B = unname(B), rank = rank,
      z_start = rank_total + 1L, z_end = rank_total + rank
    )
    rank_total <- rank_total + rank
  }
  blocks    <- Filter(function(block) block[["rank"]] > 0L, blocks)
  structure <- list(
    version           = 2L,
    metadata_hash     = metadata[["hash"]],
    row_index         = metadata[["row_index"]],
    residual_variance = numeric(K),
    residual_sei      = numeric(K),
    latent_blocks     = blocks,
    rank              = rank_total,
    source_ids        = "sampling_error",
    provenance        = list(
      input_origin       = metadata[["origin"]],
      factor_status      = metadata[["factor_status"]],
      policy             = "whole_sampling_error",
      policy_version     = 1L
    )
  )
  structure[["hash"]] <- rlang::hash(structure)
  .known_v_update(known_V, list(selection_structure = structure))
}


.selection_sampling_structure <- function(data) {

  known_V <- .data_known_v_data(data)
  if (is.null(known_V)) {
    K         <- nrow(data[["outcome"]])
    row_index <- attr(data, "selection_binding", exact = TRUE)[["row_index"]]
    if (is.null(row_index)) {
      row_index <- attr(data, "selection_model", exact = TRUE)[["groups"]][["row_index"]]
    }
    if (is.null(row_index)) {
      row_index <- seq_len(K)
    }
    sei      <- data[["outcome"]][["sei"]]
    positive <- which(sei > 0)
    blocks   <- lapply(seq_along(positive), function(index) {

      row <- positive[[index]]
      list(index = row, size = 1L, B = matrix(sei[[row]], 1L, 1L),
           rank = 1L, z_start = index, z_end = index)
    })
    return(list(
      version = 2L, row_index = row_index,
      residual_variance = numeric(K), residual_sei = numeric(K),
      latent_blocks = blocks, rank = length(positive), source_ids = "sampling_error",
      provenance = list(input_origin = "vi/sei", policy = "whole_sampling_error",
        policy_version = 1L)
    ))
  }
  .known_v_check_selection_structure(known_V)
}


.known_v_check_selection_structure <- function(known_V) {

  structure <- known_V[["selection_structure"]]
  metadata  <- .known_v_selection_metadata(known_V)
  if (!is.list(structure) || !identical(structure[["version"]], 2L)) {
    stop("The conditional sampling structure is missing or invalid.", call. = FALSE)
  }
  payload <- structure
  payload[["hash"]] <- NULL
  if (!identical(structure[["hash"]], rlang::hash(payload)) ||
      !identical(structure[["metadata_hash"]], metadata[["hash"]]) ||
      !identical(structure[["row_index"]], metadata[["row_index"]])) {
    stop("The conditional sampling structure is missing or no longer matches the retained rows.", call. = FALSE)
  }
  structure
}


.known_v_bind_diagnostics <- function(diagnostics) {

  if (length(diagnostics) == 0L) {
    return(data.frame())
  }

  out <- do.call(rbind, diagnostics)
  rownames(out) <- NULL
  out
}

.known_v_decompose_block <- function(V_block) {

  block_size <- nrow(V_block)
  diagonal   <- diag(V_block)

  rank_one_factor <- .covariance_exact_rank_one_factor(V_block)
  if (!is.null(rank_one_factor) && block_size > 1L) {
    return(list(
      residual_variance            = numeric(block_size),
      B                            = matrix(rank_one_factor, ncol = 1L),
      effective_residual_fraction  = 0,
      max_reconstruction_error     = 0,
      min_latent_eigenvalue        = 0
    ))
  }

  if (block_size == 1L || all(V_block[row(V_block) != col(V_block)] == 0)) {
    return(list(
      residual_variance            = diagonal,
      B                            = matrix(numeric(0), nrow = block_size, ncol = 0L),
      effective_residual_fraction  = 1,
      max_reconstruction_error     = 0,
      min_latent_eigenvalue        = NA_real_
    ))
  }

  # Standardization is a symmetric covariance congruence. Avoid separately
  # rounded row/column scaling while retaining the supplied V unchanged.
  correlation <- V_block / tcrossprod(sqrt(diagonal))
  diag(correlation) <- 1
  lambda_min  <- min(.covariance_factorization(correlation)[["eigenvalues"]])
  alpha_max   <- 0.99 * lambda_min

  if (alpha_max <= sqrt(.Machine$double.eps)) {
    stop(
      "A known-V block is too close to singular for a positive residual ",
      "D + BB' decomposition.",
      call. = FALSE
    )
  }

  alpha <- min(0.10, alpha_max)

  residual_variance <- alpha * diagonal
  latent_covariance <- V_block - diag(residual_variance, nrow = block_size)
  eig  <- .covariance_factorization(latent_covariance)
  values <- eig[["spectral_values"]]
  keep   <- values > 0

  if (!.covariance_is_positive_semidefinite(eig)) {
    stop("Known-V decomposition failed; V - D is not positive semidefinite.",
         call. = FALSE)
  }
  if (any(values < 0)) {
    stop("Known-V decomposition produced negative eigenvalues.", call. = FALSE)
  }
  if (any(keep)) {
    B <- eig[["eigenvectors"]][, keep, drop = FALSE] %*%
      diag(sqrt(values[keep]), nrow = sum(keep))
  } else {
    B <- matrix(numeric(0), nrow = block_size, ncol = 0L)
  }

  reconstruction <- diag(residual_variance, nrow = block_size) + tcrossprod(B)

  return(list(
    residual_variance            = residual_variance,
    B                            = B,
    effective_residual_fraction  = alpha,
    max_reconstruction_error     = max(abs(reconstruction - V_block)),
    min_latent_eigenvalue        = min(values)
  ))
}
