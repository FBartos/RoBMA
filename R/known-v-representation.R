# Known-V input helpers -----

#' Declare a diagonal-plus-factor sampling covariance
#'
#' @description
#' Constructs an exact known sampling covariance representation
#' \deqn{V = \mathrm{diag}(d) + U U^\mathsf{T}.}
#' The supplied diagonal and loading matrix define `V`; no numerical rank is
#' inferred from a materialized covariance matrix. This representation permits
#' likelihood implementations to use the declared factor dimension when that
#' route is supported. Gaussian and selection models retain the same target
#' as a conventional dense `V` input with the same source settings and, for
#' best-p-value selection, publication groups.
#'
#' With `selection_model(known_sampling_variance = "condition")`, the declared
#' covariance defines the complete retained sampling-error vector, including
#' both diagonal and factor variation. With `known_sampling_variance = "integrate"`
#' (the default), the full sampling error is integrated before selection
#' normalization. The diagonal-plus-factor split is a numerical representation;
#' it does not define which sampling sources are conditioned upon.
#'
#' Ordinary covariance matrices, including results from [metafor::vcalc()],
#' can also be supplied directly. Covariance inputs do not specify publication
#' groups. Only `weight_rule = "best"` uses publication groups; declare them
#' separately with `selection_model(weight_rule = "best", group = ...)`.
#' The default product rule requires no publication grouping.
#'
#' Matrix inputs use base R's numerical symmetry tolerance. Accepted matrices
#' are stored symmetrically by averaging corresponding off-diagonal entries,
#' leaving the supplied R object and its diagonal unchanged. Positive
#' semidefiniteness is judged from the standardized eigenvalues, with values
#' inside the eigensolver's roundoff envelope treated as zero. Exactly singular
#' designs are therefore accepted and factorized to their true rank, including
#' [metafor::vcalc()] output with `rho = 1`, whose diagonal and off-diagonal
#' reach working precision along different rounding paths. A singular
#' covariance still requires model structure covering its null space; see
#' [brma.mv()]. The stored covariance uses this symmetric representation and
#' retains input type and row mapping; it does not include a separate copy of
#' unequal input triangles.
#'
#' @param diagonal finite non-negative numeric vector `d`.
#' @param loading finite numeric matrix `U` with one row per element of
#'   `diagonal`. A zero-column matrix is allowed and represents a diagonal
#'   covariance.
#'
#' @return A `RoBMA_known_v_factor` object accepted as the `V` or `V_new`
#'   argument of multivariate-model interfaces.
#'
#' @export
known_v_factor <- function(diagonal, loading) {

  if (!is.numeric(diagonal) || !is.null(dim(diagonal)) ||
      length(diagonal) == 0L || anyNA(diagonal) ||
      any(!is.finite(diagonal)) || any(diagonal < 0)) {
    stop("'diagonal' must be a finite non-negative numeric vector.",
         call. = FALSE)
  }
  if (!is.numeric(loading) || !is.matrix(loading) ||
      nrow(loading) != length(diagonal) || anyNA(loading) ||
      any(!is.finite(loading))) {
    stop(
      "'loading' must be a finite numeric matrix with one row per ",
      "element of 'diagonal'.",
      call. = FALSE
    )
  }
  total_diagonal <- as.numeric(diagonal) + rowSums(loading^2)
  if (any(!is.finite(total_diagonal))) {
    stop("The declared sampling covariance must contain finite variances.",
         call. = FALSE)
  }

  factor <- structure(
    list(
      diagonal = as.numeric(diagonal),
      loading  = unname(loading)
    ),
    class = c("RoBMA_known_v_factor", "list")
  )
  .known_v_stamp_factor(factor, list(
    version     = 1L,
    row_index   = seq_along(diagonal),
    input_nrow  = length(diagonal),
    source_ids  = if (ncol(loading)) paste0("sampling_factor_", seq_len(ncol(loading))) else character()
  ))
}


.known_v_stamp_factor <- function(factor, metadata) {

  metadata[["hash"]] <- NULL
  metadata[["hash"]] <- rlang::hash(list(
    diagonal = factor[["diagonal"]], loading = factor[["loading"]], metadata = metadata
  ))
  attr(factor, "RoBMA_factor_metadata") <- metadata
  factor
}


.known_v_factor_metadata <- function(factor, arg = "V") {

  metadata <- attr(factor, "RoBMA_factor_metadata", exact = TRUE)
  if (!is.list(metadata) || !identical(metadata[["version"]], 1L)) {
    stop("The '", arg, "' factor metadata are invalid.", call. = FALSE)
  }
  payload <- metadata
  payload[["hash"]] <- NULL
  if (!identical(metadata[["hash"]], rlang::hash(list(
        diagonal = factor[["diagonal"]], loading = factor[["loading"]], metadata = payload))) ||
      !.known_v_valid_row_map(metadata[["row_index"]], metadata[["input_nrow"]],
                             length(factor[["diagonal"]])) ||
      !is.character(metadata[["source_ids"]]) ||
      length(metadata[["source_ids"]]) != ncol(factor[["loading"]]) ||
      anyNA(metadata[["source_ids"]]) || any(!nzchar(metadata[["source_ids"]])) ||
      anyDuplicated(metadata[["source_ids"]])) {
    stop("The '", arg, "' factor metadata no longer match its declaration.", call. = FALSE)
  }
  metadata
}


.known_v_factor_components <- function(V, arg = "V") {

  factor <- V
  if (!inherits(factor, "RoBMA_known_v_factor") || !is.list(factor) ||
      !identical(names(factor), c("diagonal", "loading"))) {
    stop("The '", arg, "' factor representation is invalid.", call. = FALSE)
  }
  components <- tryCatch(
    known_v_factor(factor[["diagonal"]], factor[["loading"]]),
    error = function(e) {
      stop("The '", arg, "' factor representation is invalid: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  metadata <- .known_v_factor_metadata(factor, arg = arg)
  .known_v_stamp_factor(components, metadata)
}


.known_v_has_declared_factor <- function(V, arg = "V") {

  inherits(V, "RoBMA_known_v_factor")
}


# A certified computational factor for the stored covariance, declared by the
# user through known_v_factor() or recovered exactly from the supplied entries.
# V stays authoritative: this is a representation of the same matrix, never a
# source-role or conditioning change.
#
# Recovered factors are stored per dependency block. Their columns never cross
# blocks, so block-sparse storage keeps the representation compact for inputs
# with many small blocks; consumers that need one global loading materialize
# it through .known_v_certified_factor_loading().
.known_v_certified_factor <- function(known_V) {

  if (identical(.known_v_storage(known_V), "factor")) {
    return(list(status = "declared"))
  }

  known_V[["certified_factor"]]
}


.known_v_has_certified_factor <- function(known_V) {

  !is.null(.known_v_certified_factor(known_V))
}


.known_v_certified_factor_status <- function(known_V) {

  factor <- .known_v_certified_factor(known_V)
  if (is.null(factor)) {
    return("undeclared")
  }

  factor[["status"]]
}


# The residual diagonal of the certified factor, one value per row.
.known_v_certified_factor_diagonal <- function(known_V) {

  if (identical(.known_v_storage(known_V), "factor")) {
    return(as.numeric(known_V[["factor_diagonal"]]))
  }

  as.numeric(known_V[["certified_factor"]][["diagonal"]])
}


# The certified loading as one K x p matrix, ordered by dependency block.
.known_v_certified_factor_loading <- function(known_V) {

  if (identical(.known_v_storage(known_V), "factor")) {
    return(unname(known_V[["factor_loading"]]))
  }

  factor  <- known_V[["certified_factor"]]
  K       <- .known_v_nrow(known_V)
  loading <- matrix(0, nrow = K, ncol = factor[["rank"]])
  offset  <- 0L
  for (block in factor[["blocks"]]) {
    columns <- offset + seq_len(ncol(block[["loading"]]))
    loading[block[["index"]], columns] <- block[["loading"]]
    offset <- offset + ncol(block[["loading"]])
  }

  loading
}


# Correlated rows a certified factor does not represent. Their blocks keep the
# supplied entries and the dense route.
.known_v_certified_factor_dense_rows <- function(known_V) {

  if (identical(.known_v_storage(known_V), "factor")) {
    return(integer(0))
  }

  as.integer(known_V[["certified_factor"]][["dense_rows"]])
}


# The certified factor of each correlated dependency block, block-local.
.known_v_certified_factor_blocks <- function(known_V) {

  if (identical(.known_v_storage(known_V), "factor")) {
    loading  <- known_V[["factor_loading"]]
    diagonal <- as.numeric(known_V[["factor_diagonal"]])
    return(lapply(known_V[["block_indices"]], function(index) {
      block_loading <- loading[index, , drop = FALSE]
      list(
        index    = index,
        diagonal = diagonal[index],
        loading  = unname(block_loading[, colSums(abs(block_loading)) > 0, drop = FALSE])
      )
    }))
  }

  factor   <- known_V[["certified_factor"]]
  diagonal <- as.numeric(factor[["diagonal"]])
  lapply(factor[["blocks"]], function(block) {
    list(
      index    = block[["index"]],
      diagonal = diagonal[block[["index"]]],
      loading  = block[["loading"]]
    )
  })
}


# Validate a recovered certified factor against the stored covariance. The
# reconstruction identity is re-checked here so that a representation carried
# through subsetting, updating, or caching cannot drift from its entries.
.known_v_validate_certified_factor <- function(known_V) {

  factor <- known_V[["certified_factor"]]
  if (is.null(factor)) {
    return(invisible(known_V))
  }

  K <- .known_v_nrow(known_V)
  if (identical(.known_v_storage(known_V), "factor")) {
    stop("Internal error: declared factor storage cannot carry a recovered factor.",
         call. = FALSE)
  }
  diagonal   <- factor[["diagonal"]]
  blocks     <- factor[["blocks"]]
  dense_rows <- factor[["dense_rows"]]
  if (!is.list(factor) || !identical(factor[["status"]], "recovered_block_constant") ||
      !is.numeric(diagonal) || !is.null(dim(diagonal)) || length(diagonal) != K ||
      anyNA(diagonal) || any(!is.finite(diagonal)) || any(diagonal < 0) ||
      !is.integer(dense_rows) || anyNA(dense_rows) || anyDuplicated(dense_rows) ||
      any(dense_rows < 1L) || any(dense_rows > K) ||
      !is.list(blocks) || length(blocks) == 0L ||
      !identical(factor[["rank"]], sum(vapply(blocks, function(block) {
        ncol(block[["loading"]])
      }, integer(1))))) {
    stop("Internal error: recovered known-V factor metadata are invalid.",
         call. = FALSE)
  }

  # The supplied blocks own the partition; a corrupted one must fail there.
  partition <- lapply(.known_v_blocks(known_V), `[[`, "index")
  .known_v_validate_dependency_blocks(partition, K)
  covered <- integer(0)

  for (block in blocks) {
    index   <- block[["index"]]
    loading <- block[["loading"]]
    if (!is.numeric(loading) || !is.matrix(loading) ||
        nrow(loading) != length(index) || ncol(loading) == 0L ||
        anyNA(loading) || any(!is.finite(loading)) ||
        !any(vapply(partition, function(rows) identical(rows, index), logical(1)))) {
      stop("Internal error: recovered known-V factor metadata are invalid.",
           call. = FALSE)
    }
    if (any(colSums(loading != 0) < 2L)) {
      stop("Internal error: a recovered known-V factor column has no dependence.",
           call. = FALSE)
    }
    reconstruction <- tcrossprod(loading)
    diag(reconstruction) <- diag(reconstruction) + diagonal[index]
    supplied <- .known_v_block_covariance(known_V, index)
    scale    <- max(abs(supplied))
    residual <- max(abs(reconstruction - supplied)) / scale
    if (!is.finite(residual) ||
        residual > 8 * length(index) * .Machine$double.eps) {
      stop("Internal error: a recovered known-V factor no longer reproduces 'V'.",
           call. = FALSE)
    }
    covered <- c(covered, index)
  }

  if (anyDuplicated(covered) || length(intersect(covered, dense_rows)) > 0L) {
    stop("Internal error: recovered known-V factor blocks overlap.", call. = FALSE)
  }
  independent <- setdiff(seq_len(K), covered)
  if (!identical(as.numeric(diagonal[independent]),
                 as.numeric(.known_v_diagonal(known_V)[independent]))) {
    stop("Internal error: a recovered known-V factor changed an independent variance.",
         call. = FALSE)
  }

  invisible(known_V)
}


# The supplied entries of one dependency block, without materializing V.
.known_v_block_covariance <- function(known_V, index) {

  if (length(index) == 1L) {
    return(matrix(.known_v_diagonal(known_V)[[index]], 1L, 1L))
  }
  V <- known_V[["V"]]
  if (!is.null(V)) {
    return(V[index, index, drop = FALSE])
  }
  for (block in known_V[["blocks"]]) {
    if (identical(as.integer(block[["index"]]), as.integer(index))) {
      return(block[["covariance"]])
    }
  }

  stop("Internal error: unknown known-V dependency block.", call. = FALSE)
}


.known_v_factor_covariance <- function(diagonal, loading, index = NULL) {

  if (is.null(index)) {
    index <- seq_along(diagonal)
  }
  out <- tcrossprod(loading[index, , drop = FALSE])
  diag(out) <- diag(out) + diagonal[index]
  out
}


.known_v_input_metadata <- function(V, arg = "V") {

  factor <- if (inherits(V, "RoBMA_known_v_factor")) V else NULL
  factor_metadata <- if (!is.null(factor)) .known_v_factor_metadata(factor, arg) else NULL
  K <- .known_v_input_nrow(V, arg = arg)
  metadata <- list(
    version             = 2L,
    origin              = if (!is.null(factor)) "known_v_factor" else "matrix",
    row_index           = if (!is.null(factor_metadata)) factor_metadata[["row_index"]] else seq_len(K),
    input_nrow          = if (!is.null(factor_metadata)) factor_metadata[["input_nrow"]] else K,
    factor_status       = if (!is.null(factor)) "declared" else "undeclared",
    source_ids          = factor_metadata[["source_ids"]]
  )
  .known_v_stamp_selection_metadata(metadata)
}


.known_v_stamp_selection_metadata <- function(metadata) {

  metadata[["hash"]] <- NULL
  metadata[["hash"]] <- rlang::hash(metadata)
  metadata
}


.known_v_subset_selection_metadata <- function(metadata, row_index) {

  metadata[["row_index"]] <- metadata[["row_index"]][row_index]
  .known_v_stamp_selection_metadata(metadata)
}


.known_v_selection_metadata <- function(known_V) {

  metadata <- known_V[["selection_metadata"]]
  .known_v_validate_selection_metadata(metadata, .known_v_nrow(known_V))
  metadata
}


.known_v_validate_selection_metadata <- function(metadata, K) {

  if (!is.list(metadata) || !identical(metadata[["version"]], 2L)) {
    stop("The known-V selection metadata are missing or invalid.", call. = FALSE)
  }
  payload <- metadata
  payload[["hash"]] <- NULL
  if (!identical(metadata[["hash"]], rlang::hash(payload)) ||
      !.known_v_valid_row_map(metadata[["row_index"]], metadata[["input_nrow"]], K)) {
    stop("The known-V selection metadata no longer match their declared rows.", call. = FALSE)
  }
  invisible(metadata)
}


.known_v_valid_row_map <- function(row_index, input_nrow, K) {

  is.numeric(input_nrow) && length(input_nrow) == 1L &&
    !is.na(input_nrow) && is.finite(input_nrow) && input_nrow >= K &&
    input_nrow == floor(input_nrow) &&
    is.numeric(row_index) && is.null(dim(row_index)) && length(row_index) == K &&
    !anyNA(row_index) && all(is.finite(row_index)) &&
    !anyDuplicated(row_index) && all(row_index >= 1L) &&
    all(row_index <= input_nrow) && all(row_index == floor(row_index))
}

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
  if (length(storage) != 1L ||
      !storage %in% c("diagonal", "blocks", "dense", "factor") ||
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
  if (identical(storage, "factor")) {
    factor_diagonal <- known_V[["factor_diagonal"]]
    factor_loading  <- known_V[["factor_loading"]]
    valid_factor <- is.numeric(factor_diagonal) &&
      is.null(dim(factor_diagonal)) && length(factor_diagonal) == K &&
      !anyNA(factor_diagonal) && all(is.finite(factor_diagonal)) &&
      all(factor_diagonal >= 0) && is.numeric(factor_loading) &&
      is.matrix(factor_loading) && nrow(factor_loading) == K &&
      !anyNA(factor_loading) && all(is.finite(factor_loading))
    if (!valid_factor ||
        !identical(
          as.numeric(diagonal),
          as.numeric(factor_diagonal) + rowSums(factor_loading^2)
        )) {
      stop("Internal error: factor known-V representation is invalid.",
           call. = FALSE)
    }
  } else if (!is.null(known_V[["factor_diagonal"]]) ||
             !is.null(known_V[["factor_loading"]])) {
    stop("Internal error: non-factor known-V representation contains factors.",
         call. = FALSE)
  }
  .known_v_validate_certified_factor(known_V)

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


.known_v_as_matrix <- function(V, k = NULL, warn_singular = TRUE) {

  V_matrix <- .known_v_as_matrix_structure(V, k = k)

  if (anyNA(V_matrix) || any(!is.finite(V_matrix))) {
    stop("The 'V' argument must contain only finite non-missing values.", call. = FALSE)
  }

  .known_v_check_symmetric(V_matrix, "'V'")
  V_matrix <- .known_v_exact_symmetrize(V_matrix)

  diagonal <- diag(V_matrix)
  if (any(diagonal <= 0)) {
    stop("The diagonal of 'V' must contain positive variances.", call. = FALSE)
  }

  classification <- .known_v_covariance_classification(V_matrix)
  if (!isTRUE(classification[["positive_semidefinite"]])) {
    stop("The 'V' argument must be positive semidefinite.", call. = FALSE)
  }
  if (isTRUE(classification[["singular"]]) && isTRUE(warn_singular)) {
    .known_v_warn_singular()
  }

  return(V_matrix)
}


# Convert known-V input without validating covariance values.
.known_v_as_matrix_structure <- function(V, k = NULL) {

  if (.known_v_has_declared_factor(V)) {
    components <- .known_v_factor_components(V)
    V_matrix <- .known_v_factor_covariance(
      components[["diagonal"]],
      components[["loading"]]
    )
  } else if (is.matrix(V)) {
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

  if (.known_v_has_declared_factor(V, arg = arg)) {
    .known_v_factor_components(V, arg = arg)
    return("factor")
  }
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
    "a non-empty list of square matrices, or a known_v_factor() object.",
    call. = FALSE
  )
}


.known_v_input_nrow <- function(V, arg = "V") {

  storage <- .known_v_input_storage(V, arg = arg)
  if (storage == "factor") {
    return(length(.known_v_factor_components(V, arg = arg)[["diagonal"]]))
  }
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
  if (storage == "factor") {
    components <- .known_v_factor_components(V, arg = arg)
    return(components[["diagonal"]] + rowSums(components[["loading"]]^2))
  }
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
  if (is.logical(keep_rows) && length(keep_rows) == K && !anyNA(keep_rows)) {
    row_index <- which(keep_rows)
  } else if (is.numeric(keep_rows) && !anyNA(keep_rows) &&
             all(is.finite(keep_rows)) && all(keep_rows == as.integer(keep_rows)) &&
             all(keep_rows >= 1L & keep_rows <= K) && !anyDuplicated(keep_rows)) {
    row_index <- as.integer(keep_rows)
  } else {
    stop("Internal error: invalid known-V row selector.", call. = FALSE)
  }

  storage <- .known_v_input_storage(V)
  if (storage == "factor") {
    components <- .known_v_factor_components(V)
    metadata <- .known_v_factor_metadata(components)
    metadata[["row_index"]] <- metadata[["row_index"]][row_index]
    return(.known_v_stamp_factor(known_v_factor(
      diagonal = components[["diagonal"]][row_index],
      loading  = components[["loading"]][row_index, , drop = FALSE]
    ), metadata))
  }
  if (storage == "dense") {
    return(V[row_index, row_index, drop = FALSE])
  }
  if (storage == "diagonal") {
    return(as.numeric(V)[row_index])
  }

  if (is.unsorted(row_index)) {
    return(.known_v_blockdiag(V)[row_index, row_index, drop = FALSE])
  }
  keep_rows <- seq_len(K) %in% row_index

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
  if (storage == "factor") {
    diagonal      <- known_V[["factor_diagonal"]]
    loading       <- known_V[["factor_loading"]]
    block_indices <- known_V[["block_indices"]]
    return(lapply(block_indices, function(index) {
      list(
        index      = index,
        covariance = .known_v_factor_covariance(diagonal, loading, index)
      )
    }))
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


.known_v_covariance_matrix <- function(known_V) {

  K          <- .known_v_nrow(known_V)
  covariance <- matrix(0, nrow = K, ncol = K)
  for (block in .known_v_blocks(known_V)) {
    index <- block[["index"]]
    covariance[index, index] <- block[["covariance"]]
  }

  covariance
}


.known_v_dependency_covariance <- function(
    data, sampling_latent_marginalized = FALSE) {

  known_V <- .data_known_v_data(data)
  backend <- .data_known_v_effective_backend(data)
  if (backend == "latent" && !sampling_latent_marginalized) {
    return(diag(
      .known_v_residual_variance(known_V),
      nrow = .known_v_nrow(known_V)
    ))
  }

  .known_v_covariance_matrix(known_V)
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
  if (.known_v_storage(known_V) == "factor") {
    return(.known_v_factor_covariance(
      known_V[["factor_diagonal"]],
      known_V[["factor_loading"]]
    ))
  }
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
  if (storage == "factor") {
    factor <- known_v_factor(
      known_V[["factor_diagonal"]],
      known_V[["factor_loading"]]
    )
    metadata <- .known_v_selection_metadata(known_V)
    return(.known_v_stamp_factor(factor, list(
      version    = 1L,
      row_index  = metadata[["row_index"]],
      input_nrow = metadata[["input_nrow"]],
      source_ids = metadata[["source_ids"]]
    )))
  }
  if (storage == "diagonal") {
    return(.known_v_diagonal(known_V))
  }
  if (storage == "dense") {
    return(known_V[["V"]])
  }

  lapply(.known_v_blocks(known_V), `[[`, "covariance")
}

.known_v_check_symmetric <- function(V_matrix, arg) {

  if (!isSymmetric(V_matrix, check.attributes = FALSE)) {
    stop(arg, " must be symmetric.", call. = FALSE)
  }

  return(invisible(TRUE))
}


.known_v_exact_symmetrize <- function(V_matrix) {

  if (identical(V_matrix, t(V_matrix))) {
    return(V_matrix)
  }

  upper <- upper.tri(V_matrix)
  midpoint <- rowMeans(cbind(
    V_matrix[upper],
    t(V_matrix)[upper]
  ))
  V_matrix[upper] <- midpoint
  V_matrix[lower.tri(V_matrix)] <- t(V_matrix)[lower.tri(V_matrix)]

  V_matrix
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

  known_V           <- .known_v_canonicalize_newdata(V_new)
  correlated        <- length(.known_v_correlated_blocks(known_V)) > 0L
  residual_variance <- .known_v_diagonal(known_V)
  .known_v_update(known_V, list(
    parameterization  = "block_mvn",
    effective_backend = if (correlated) "block_mvn" else "diagonal",
    correlated        = correlated,
    residual_variance = residual_variance,
    residual_sei      = sqrt(residual_variance),
    rank              = 0L
  ))
}


.known_v_canonicalize_newdata <- function(V_new) {

  storage <- .known_v_input_storage(V_new, arg = "V_new")
  K       <- .known_v_input_nrow(V_new, arg = "V_new")
  metadata <- .known_v_input_metadata(V_new, arg = "V_new")
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
      selection_metadata = metadata,
      storage  = "diagonal",
      K        = K,
      diagonal = diagonal,
      blocks   = list(),
      singular = any(diagonal == 0)
    )))
  }

  if (storage == "factor") {
    components <- .known_v_factor_components(V_new, arg = "V_new")
    diagonal   <- components[["diagonal"]] +
      rowSums(components[["loading"]]^2)
    covariance <- .known_v_factor_covariance(
      components[["diagonal"]],
      components[["loading"]]
    )
    if (any(!is.finite(covariance))) {
      stop("'V_new' must contain only finite non-missing values.",
           call. = FALSE)
    }
    block_indices <- .known_v_block_indices(covariance)
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
      singular        = .known_v_newdata_block_is_singular(covariance)
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
      selection_metadata = metadata,
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
    selection_metadata = metadata,
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
  V_new <- .known_v_exact_symmetrize(V_new)
  if (any(diag(V_new) < 0)) {
    stop("The diagonal of 'V_new' must contain non-negative variances.",
         call. = FALSE)
  }
  zero_variance <- diag(V_new) == 0
  if (any(zero_variance) &&
      any(V_new[zero_variance, , drop = FALSE] != 0)) {
    stop("'V_new' must be positive semidefinite.", call. = FALSE)
  }
  classification <- .known_v_covariance_classification(V_new)
  if (!isTRUE(classification[["positive_semidefinite"]])) {
    stop("'V_new' must be positive semidefinite.", call. = FALSE)
  }

  V_new
}


.known_v_newdata_block_is_singular <- function(V_new) {

  any(diag(V_new) == 0) ||
    isTRUE(.known_v_covariance_classification(V_new)[["singular"]])
}


# Classify dependency blocks without modifying the supplied covariance.
.known_v_covariance_classification <- function(V) {

  positive_variance <- diag(V) > 0
  if (!any(positive_variance)) {
    return(list(positive_semidefinite = TRUE, singular = TRUE))
  }

  covariance <- V[positive_variance, positive_variance, drop = FALSE]
  indices    <- .known_v_block_indices(covariance)
  singular   <- any(!positive_variance)

  # The factorization is the single authority on block validity. A pairwise
  # bound |cov_ij| <= sd_i sd_j adds nothing: an excess e makes the 2x2
  # principal submatrix carry eigenvalue -e, so by Cauchy interlacing the block
  # has lambda_min <= -e and the factorization already rejects every excess
  # above its own roundoff tolerance. Screening the pairs separately only
  # applied a second, stricter convention inside that tolerance, which refused
  # exactly-singular designs such as metafor::vcalc(rho = 1) whose off-diagonal
  # and diagonal reached working precision along different rounding paths.
  for (index in indices) {
    block <- covariance[index, index, drop = FALSE]

    factorization <- .covariance_factorization(block)
    if (!.covariance_is_positive_semidefinite(factorization)) {
      return(list(positive_semidefinite = FALSE, singular = TRUE))
    }
    singular <- singular || isTRUE(factorization[["singular"]])
  }

  return(list(positive_semidefinite = TRUE, singular = singular))
}
