# ============================================================================ #
# brma.unit_level.R
# ============================================================================ #
#
# Shared helpers for diagnostics that need to distinguish the output unit from
# the conditioning depth.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# .normalize_unit
# ---------------------------------------------------------------------------- #
#
# @param unit character; output/deletion unit.
#
# @return normalized unit value.
#
# ---------------------------------------------------------------------------- #
.normalize_unit <- function(unit) {

  return(match.arg(unit, c("estimate", "cluster")))
}


# ---------------------------------------------------------------------------- #
# .normalize_conditioning_depth
# ---------------------------------------------------------------------------- #
#
# @param conditioning_depth character; conditioning depth.
#
# @return normalized conditioning depth value.
#
# ---------------------------------------------------------------------------- #
.normalize_conditioning_depth <- function(conditioning_depth) {

  return(match.arg(conditioning_depth, c("marginal", "cluster", "estimate")))
}


# ---------------------------------------------------------------------------- #
# .check_unit_conditioning_depth
# ---------------------------------------------------------------------------- #
#
# Validate combinations of unit and conditioning depth for diagnostics.
#
# @param object             brma object.
# @param unit               character; output/deletion unit.
# @param conditioning_depth character; conditioning depth.
# @param caller             character; caller name for error messages.
#
# @return invisible NULL.
#
# ---------------------------------------------------------------------------- #
.check_unit_conditioning_depth <- function(object, unit, conditioning_depth, caller) {

  is_multilevel <- .is_multilevel(object)

  if (unit == "cluster" && !is_multilevel) {
    stop(caller, " with unit = 'cluster' is only available for multilevel models.",
         call. = FALSE)
  }

  if (conditioning_depth == "cluster" && !is_multilevel) {
    stop(caller, " with conditioning_depth = 'cluster' is only available for multilevel models.",
         call. = FALSE)
  }

  if (unit == "cluster" && conditioning_depth == "estimate") {
    stop(caller, " does not support unit = 'cluster' with conditioning_depth = 'estimate'.",
         call. = FALSE)
  }

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .loo_conditioning_depth_from_unit
# ---------------------------------------------------------------------------- #
#
# LOO exposes only the deletion unit. The implied conditioning depth is stored
# as metadata.
#
# @param unit character; output/deletion unit.
#
# @return character; implied conditioning depth.
#
# ---------------------------------------------------------------------------- #
.loo_conditioning_depth_from_unit <- function(unit) {

  unit <- .normalize_unit(unit)

  if (unit == "estimate") {
    return("estimate")
  } else {
    return("cluster")
  }
}


# ---------------------------------------------------------------------------- #
# .get_cluster_labels
# ---------------------------------------------------------------------------- #
#
# @param object brma object.
#
# @return named character vector mapping numeric cluster indices to labels.
#
# ---------------------------------------------------------------------------- #
.get_cluster_labels <- function(object) {

  if (!.is_multilevel(object)) {
    return(NULL)
  }

  outcome_data <- object[["data"]][["outcome"]]
  cluster      <- outcome_data[["cluster"]]

  if ("cluster_label" %in% names(outcome_data)) {
    labels <- tapply(outcome_data[["cluster_label"]], cluster, function(x) x[1])
  } else {
    labels <- tapply(cluster, cluster, function(x) as.character(x[1]))
  }
  label_names <- names(labels)
  labels      <- as.character(labels)
  names(labels) <- label_names

  return(labels)
}


# ---------------------------------------------------------------------------- #
# .get_estimate_labels
# ---------------------------------------------------------------------------- #
#
# @param object brma object.
#
# @return character vector of estimate labels.
#
# ---------------------------------------------------------------------------- #
.get_estimate_labels <- function(object) {

  outcome_data <- object[["data"]][["outcome"]]

  if ("slab" %in% names(outcome_data)) {
    labels <- as.character(outcome_data[["slab"]])
  } else {
    labels <- as.character(seq_len(nrow(outcome_data)))
  }

  return(labels)
}


# ---------------------------------------------------------------------------- #
# .get_outcome_hash
# ---------------------------------------------------------------------------- #
#
# @param object brma object.
#
# @return deterministic hash of the outcome target.
#
# ---------------------------------------------------------------------------- #
.get_outcome_hash <- function(object) {

  outcome_type <- .outcome_type(object)
  outcome_data <- object[["data"]][["outcome"]]

  payload <- list(outcome_type = outcome_type)

  if (outcome_type == "norm") {
    payload[["outcome"]] <- list(
      yi  = unname(as.numeric(outcome_data[["yi"]])),
      sei = unname(as.numeric(outcome_data[["sei"]]))
    )
  } else if (outcome_type == "bin") {
    payload[["outcome"]] <- list(
      ai  = unname(as.integer(outcome_data[["ai"]])),
      ci  = unname(as.integer(outcome_data[["ci"]])),
      n1i = unname(as.integer(outcome_data[["n1i"]])),
      n2i = unname(as.integer(outcome_data[["n2i"]]))
    )
  } else if (outcome_type == "pois") {
    payload[["outcome"]] <- list(
      x1i = unname(as.integer(outcome_data[["x1i"]])),
      x2i = unname(as.integer(outcome_data[["x2i"]])),
      t1i = unname(as.numeric(outcome_data[["t1i"]])),
      t2i = unname(as.numeric(outcome_data[["t2i"]]))
    )
  }

  if ("weights" %in% names(outcome_data)) {
    payload[["weights"]] <- unname(as.numeric(outcome_data[["weights"]]))
  }

  if (.is_data_known_v(object[["data"]])) {
    known_V <- .data_known_v_data(object[["data"]])
    blocks  <- .known_v_correlated_blocks(known_V)
    diagonal <- .normalize_hash_zero(.known_v_diagonal(known_V))
    payload[["known_V"]] <- list(
      version  = 2L,
      K        = .known_v_nrow(known_V),
      diagonal = unname(diagonal),
      blocks   = lapply(blocks, function(block) {
        list(
          index      = unname(as.integer(block[["index"]])),
          covariance = unname(.normalize_hash_zero(
            as.numeric(block[["covariance"]])
          ))
        )
      })
    )
  }

  bytes <- as.integer(.outcome_hash_bytes(payload))
  hash1 <- 5381
  hash2 <- 0

  for (byte in bytes) {
    hash1 <- (hash1 * 33 + byte) %% 2147483647
    hash2 <- (hash2 * 65599 + byte) %% 2147483629
  }

  return(paste0(
    "v1:",
    sprintf("%08x", as.integer(hash1)),
    sprintf("%08x", as.integer(hash2))
  ))
}


# Encode outcome-hash payloads independently of the R serialization version.
.outcome_hash_bytes <- function(x) {

  if (length(x) > .Machine$integer.max) {
    stop("Internal error: outcome-hash payload is too long.",
         call. = FALSE)
  }

  encode_length <- function(x) {
    writeBin(as.integer(x), raw(), size = 4L, endian = "big")
  }
  encode_character <- function(x) {
    value_bytes <- lapply(x, function(value) {
      if (is.na(value)) {
        return(encode_length(-1L))
      } else {
        bytes <- charToRaw(enc2utf8(value))
        return(c(encode_length(length(bytes)), bytes))
      }
    })
    do.call(c, c(
      list(charToRaw("c"), encode_length(length(x))),
      value_bytes
    ))
  }

  if (is.null(x)) {
    return(charToRaw("n"))
  }
  if (is.list(x)) {
    element_bytes <- lapply(x, .outcome_hash_bytes)
    return(do.call(c, c(
      list(
        charToRaw("l"),
        encode_length(length(x)),
        .outcome_hash_bytes(names(x))
      ),
      element_bytes
    )))
  }
  if (is.double(x)) {
    x[!is.na(x) & x == 0] <- 0
    return(c(
      charToRaw("d"),
      encode_length(length(x)),
      writeBin(x, raw(), size = 8L, endian = "big")
    ))
  }
  if (is.integer(x)) {
    return(c(
      charToRaw("i"),
      encode_length(length(x)),
      writeBin(x, raw(), size = 4L, endian = "big")
    ))
  }
  if (is.character(x)) {
    return(encode_character(x))
  }

  stop("Internal error: unsupported outcome-hash payload type.",
       call. = FALSE)
}


# ---------------------------------------------------------------------------- #
# .check_brma_compare_targets
# ---------------------------------------------------------------------------- #
#
# Reject Bayes factor/model-probability comparisons across different data.
#
# @param objects list of brma objects.
# @param caller  character; caller name for error messages.
#
# @return invisible NULL.
#
# ---------------------------------------------------------------------------- #
.check_brma_compare_targets <- function(objects, caller) {

  data_hashes <- vapply(objects, .get_outcome_hash, character(1))

  if (length(unique(data_hashes)) > 1) {
    stop(caller, " requires models fitted to the same outcome data.",
         call. = FALSE)
  }

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .check_cluster_unit_deferred
# ---------------------------------------------------------------------------- #
#
# Error used by diagnostics whose cluster-unit implementation is intentionally
# deferred.
#
# @param caller character; caller name for error messages.
#
# @return stops.
#
# ---------------------------------------------------------------------------- #
.check_cluster_unit_deferred <- function(caller, argument = "unit") {
  # Cluster residual diagnostics need a separate Mahalanobis/chi-square design
  stop(
    caller, " with ", argument, " = 'cluster' is not implemented currently.",
    call. = FALSE
  )
}


.normalize_hash_zero <- function(x) {

  x[x == 0] <- 0
  x
}


# ---------------------------------------------------------------------------- #
# .check_random_formula_postfit_deferred
# ---------------------------------------------------------------------------- #
#
# Error used by post-fit methods whose random-formula semantics still need a
# dedicated RoBMA design.
#
# @param object brma object.
# @param caller character; caller name for error messages.
#
# @return invisible NULL or stops.
#
# ---------------------------------------------------------------------------- #
.check_random_formula_postfit_deferred <- function(object, caller) {

  if (.is_random(object)) {
    stop(
      caller,
      " is not implemented for brma.mv() random-formula models yet.",
      call. = FALSE
    )
  }

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .check_log_lik_target_available
# ---------------------------------------------------------------------------- #
#
# Shared availability gate for pointwise log-likelihood, LOO, and WAIC.
#
# @param object brma object.
# @param unit   character; normalized or raw output/deletion unit.
# @param caller character; caller name for error messages.
#
# @return invisible NULL or stops.
#
# ---------------------------------------------------------------------------- #
.check_log_lik_target_available <- function(object, unit, caller) {

  unit       <- .normalize_unit(unit)
  data       <- object[["data"]]
  is_brma_mv <- inherits(object, "brma.mv")

  if (unit == "cluster" && is_brma_mv) {
    mv_scope <- if (.is_data_known_v(data)) {
      "brma.mv() known-V models"
    } else {
      "brma.mv() models"
    }
    stop(
      caller,
      " with unit = 'cluster' is not implemented for ", mv_scope,
      " yet. Use unit = 'estimate'.",
      call. = FALSE
    )
  }

  if (unit == "cluster" && !.is_multilevel(object)) {
    stop(caller, " with unit = 'cluster' is only available for multilevel models.",
         call. = FALSE)
  }

  if (.is_random(object) && !(is_brma_mv && .is_data_known_v(data))) {
    .check_random_formula_postfit_deferred(object, caller)
  }

  invisible(TRUE)
}


# ---------------------------------------------------------------------------- #
# .check_legacy_level_arg
# ---------------------------------------------------------------------------- #
#
# Reject the old residual-conditioning argument name when it arrives through
# `...`.
#
# @param dots   list; captured dots.
# @param caller character; caller name for error messages.
#
# @return invisible NULL.
#
# ---------------------------------------------------------------------------- #
.check_legacy_level_arg <- function(dots, caller) {

  if ("level" %in% names(dots)) {
    stop(caller, " uses 'conditioning_depth' for residual conditioning; ",
         "do not set 'level'.", call. = FALSE)
  }

  return(invisible(NULL))
}


# ---------------------------------------------------------------------------- #
# .add_loo_target_metadata
# ---------------------------------------------------------------------------- #
#
# Store the LOO/WAIC target on the returned loo object.
#
# @param object             loo or waic object.
# @param unit               character; output/deletion unit.
# @param conditioning_depth character; implied conditioning depth.
# @param targets            character; target labels.
# @param data_hash          character; hash of the outcome target.
#
# @return object with RoBMA target metadata.
#
# ---------------------------------------------------------------------------- #
.add_loo_target_metadata <- function(object, unit, conditioning_depth, targets,
                                     data_hash, metadata = NULL) {

  target <- list(
    unit               = unit,
    conditioning_depth = conditioning_depth,
    n                  = length(targets),
    targets            = targets,
    data_hash          = data_hash
  )
  if (!is.null(metadata)) {
    extra  <- metadata[setdiff(names(metadata), names(target))]
    target <- c(target, extra)
  }

  attr(object, "RoBMA_target") <- target

  return(object)
}


# ---------------------------------------------------------------------------- #
# .get_loo_target_metadata
# ---------------------------------------------------------------------------- #
#
# @param object loo or waic object.
#
# @return target metadata list or NULL.
#
# ---------------------------------------------------------------------------- #
.get_loo_target_metadata <- function(object) {

  return(attr(object, "RoBMA_target", exact = TRUE))
}


# ---------------------------------------------------------------------------- #
# .current_predictive_target_key
# ---------------------------------------------------------------------------- #
#
# Build the current LOO/WAIC compatibility key without evaluating draws.
#
# ---------------------------------------------------------------------------- #
.current_predictive_target_key <- function(object, unit) {

  unit      <- .normalize_unit(unit)
  data_hash <- .get_outcome_hash(object)
  if (unit == "estimate") {
    target <- if (.known_v_estimate_target_uses_schur_conditioning(
        object[["data"]])) {
      "known_v_estimate"
    } else {
      "factorized_estimate"
    }
  } else {
    target <- "cluster_joint"
  }

  list(
    unit               = unit,
    conditioning_depth = .loo_conditioning_depth_from_unit(unit),
    target             = target,
    data_hash          = data_hash
  )
}


# ---------------------------------------------------------------------------- #
# .predictive_target_fingerprint
# ---------------------------------------------------------------------------- #
#
# Collapse the public LOO/WAIC compatibility key to a stable scalar.
#
# ---------------------------------------------------------------------------- #
.predictive_target_fingerprint <- function(metadata) {

  fields <- c("unit", "conditioning_depth", "target", "data_hash")
  if (is.null(metadata) || !all(fields %in% names(metadata))) {
    return(NULL)
  }
  key <- metadata[fields]
  if (any(vapply(key, function(x) length(x) != 1L || is.na(x), logical(1)))) {
    return(NULL)
  }

  paste(vapply(key, as.character, character(1)), collapse = "|")
}


# ---------------------------------------------------------------------------- #
# .check_cached_predictive_target
# ---------------------------------------------------------------------------- #
#
# Reject cached predictive diagnostics that no longer match the fitted target.
#
# ---------------------------------------------------------------------------- #
.check_cached_predictive_target <- function(object, metadata, unit, method) {

  recompute <- if (identical(method, "LOO")) "add_loo()" else "add_waic()"
  stored_fingerprint <- .predictive_target_fingerprint(metadata)
  if (is.null(stored_fingerprint)) {
    stop(
      "Stored ", method, " has incomplete RoBMA target metadata. Recompute ",
      "with ", recompute, ".",
      call. = FALSE
    )
  }
  if (!identical(metadata[["unit"]], unit)) {
    stop(
      "Stored ", method, " was computed with unit = '", metadata[["unit"]],
      "'. Recompute with ", recompute, " using unit = '", unit, "'.",
      call. = FALSE
    )
  }

  current <- .current_predictive_target_key(object, unit)
  if (!identical(metadata[["data_hash"]], current[["data_hash"]])) {
    stop(
      "Stored ", method, " does not match the current outcome data. ",
      "Recompute with ", recompute, ".",
      call. = FALSE
    )
  }
  if (!identical(
      stored_fingerprint,
      .predictive_target_fingerprint(current))) {
    stop(
      "Stored ", method, " does not match the current likelihood target. ",
      "Recompute with ", recompute, ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# ---------------------------------------------------------------------------- #
# .check_loo_target
# ---------------------------------------------------------------------------- #
#
# Check that a stored LOO object matches the requested unit.
#
# @param object brma object.
# @param unit   character; requested unit.
#
# @return stored loo object.
#
# ---------------------------------------------------------------------------- #
.check_loo_target <- function(object, unit) {

  unit       <- .normalize_unit(unit)
  .check_log_lik_target_available(object, unit, "loo()")
  loo_store  <- object[["loo"]]

  if (is.null(loo_store)) {
    stop("LOO has not been computed. Call 'object <- add_loo(object)' first.",
         call. = FALSE)
  }

  loo_result <- loo_store[[unit]]

  if (is.null(loo_result)) {
    stop(
      "LOO with unit = '", unit, "' has not been computed. ",
      "Call 'object <- add_loo(object, unit = \"", unit, "\")' first.",
      call. = FALSE
    )
  }

  metadata <- .get_loo_target_metadata(loo_result)
  .check_cached_predictive_target(object, metadata, unit, "LOO")

  return(loo_result)
}


# ---------------------------------------------------------------------------- #
# .check_waic_target
# ---------------------------------------------------------------------------- #
#
# Check that a stored WAIC object matches the requested unit.
#
# @param object brma object.
# @param unit   character; requested unit.
#
# @return stored waic object.
#
# ---------------------------------------------------------------------------- #
.check_waic_target <- function(object, unit) {

  unit       <- .normalize_unit(unit)
  .check_log_lik_target_available(object, unit, "waic()")
  waic_store <- object[["waic"]]

  if (is.null(waic_store)) {
    stop("WAIC has not been computed. Call 'object <- add_waic(object)' first.",
         call. = FALSE)
  }

  waic_result <- waic_store[[unit]]

  if (is.null(waic_result)) {
    stop(
      "WAIC with unit = '", unit, "' has not been computed. ",
      "Call 'object <- add_waic(object, unit = \"", unit, "\")' first.",
      call. = FALSE
    )
  }

  metadata <- .get_loo_target_metadata(waic_result)
  .check_cached_predictive_target(object, metadata, unit, "WAIC")

  return(waic_result)
}


# ---------------------------------------------------------------------------- #
# .get_target_conditioning_depth
# ---------------------------------------------------------------------------- #
#
# Extract the conditioning-depth metadata, accepting older cached objects that
# used `level`.
#
# @param metadata list; RoBMA target metadata.
#
# @return character scalar.
#
# ---------------------------------------------------------------------------- #
.get_target_conditioning_depth <- function(metadata) {

  if (!is.null(metadata[["conditioning_depth"]])) {
    return(metadata[["conditioning_depth"]])
  }
  if (!is.null(metadata[["level"]])) {
    return(metadata[["level"]])
  }

  return(NA_character_)
}


# ---------------------------------------------------------------------------- #
# .check_loo_compare_targets
# ---------------------------------------------------------------------------- #
#
# Reject comparisons across different LOO targets.
#
# @param loo_objects list of loo objects.
#
# @return invisible NULL.
#
# ---------------------------------------------------------------------------- #
.check_loo_compare_targets <- function(loo_objects) {

  metadata <- lapply(loo_objects, .get_loo_target_metadata)
  missing  <- vapply(metadata, is.null, logical(1))

  if (length(metadata) <= 1) {
    return(invisible(NULL))
  }

  if (any(missing)) {
    stop("LOO/WAIC objects without RoBMA target metadata cannot be compared.",
         call. = FALSE)
  }

  criteria <- vapply(loo_objects, function(object) {
    if (inherits(object, "waic")) "WAIC" else "LOO"
  }, character(1))
  if (length(unique(criteria)) > 1L) {
    stop("LOO and WAIC objects cannot be compared in the same table.",
         call. = FALSE)
  }

  missing_hash <- vapply(metadata, function(x) is.null(x[["data_hash"]]), logical(1))
  if (any(missing_hash)) {
    stop("LOO/WAIC objects without RoBMA data hashes cannot be compared.",
         call. = FALSE)
  }

  units                <- vapply(metadata, `[[`, character(1), "unit")
  conditioning_depths  <- vapply(metadata, .get_target_conditioning_depth, character(1))
  data_hashes          <- vapply(metadata, `[[`, character(1), "data_hash")

  if (length(unique(units)) > 1 ||
      length(unique(conditioning_depths)) > 1 ||
      length(unique(data_hashes)) > 1) {
    stop("LOO/WAIC objects with different data, unit, or conditioning-depth targets cannot be compared.",
         call. = FALSE)
  }

  target_kinds <- vapply(metadata, function(x) {
    if (is.null(x[["target"]])) "" else as.character(x[["target"]])
  }, character(1))
  if (any(!nzchar(target_kinds))) {
    stop("LOO/WAIC objects without likelihood target labels cannot be compared.",
         call. = FALSE)
  }
  if (length(unique(target_kinds)) > 1) {
    stop("LOO/WAIC objects with different likelihood targets cannot be compared.",
         call. = FALSE)
  }

  return(invisible(NULL))
}
