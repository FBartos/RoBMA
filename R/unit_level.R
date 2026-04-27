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
# .get_cluster_indices
# ---------------------------------------------------------------------------- #
#
# @param object brma object.
#
# @return named list of integer vectors, one per cluster.
#
# ---------------------------------------------------------------------------- #
.get_cluster_indices <- function(object) {

  if (!.is_multilevel(object)) {
    return(NULL)
  }

  outcome_data <- object[["data"]][["outcome"]]
  cluster      <- outcome_data[["cluster"]]
  indices      <- split(seq_along(cluster), cluster)
  labels       <- .get_cluster_labels(object)

  names(indices) <- labels[names(indices)]

  return(indices)
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
# @return deterministic hash of the yi/sei target.
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

  if ("cluster" %in% names(outcome_data)) {
    payload[["cluster"]] <- unname(as.integer(outcome_data[["cluster"]]))
  }

  if ("cluster_label" %in% names(outcome_data)) {
    payload[["cluster_label"]] <- unname(as.character(outcome_data[["cluster_label"]]))
  }

  bytes <- as.integer(serialize(payload, NULL, version = 3))
  hash1 <- 5381
  hash2 <- 0

  for (byte in bytes) {
    hash1 <- (hash1 * 33 + byte) %% 2147483647
    hash2 <- (hash2 * 65599 + byte) %% 2147483629
  }

  return(paste0(
    sprintf("%08x", as.integer(hash1)),
    sprintf("%08x", as.integer(hash2))
  ))
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
    stop(caller, " requires models fitted to the same yi/sei data.",
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
.check_cluster_unit_deferred <- function(caller) {
  # Cluster residual diagnostics need a separate Mahalanobis/chi-square design
  stop(caller, " with unit = 'cluster' is not implemented currently.", call. = FALSE)
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
# @param data_hash          character; hash of the yi/sei target.
#
# @return object with RoBMA target metadata.
#
# ---------------------------------------------------------------------------- #
.add_loo_target_metadata <- function(object, unit, conditioning_depth, targets,
                                     data_hash) {

  attr(object, "RoBMA_target") <- list(
    unit               = unit,
    conditioning_depth = conditioning_depth,
    n                  = length(targets),
    targets            = targets,
    data_hash          = data_hash
  )

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

  if (!is.null(metadata) && metadata[["unit"]] != unit) {
    stop(
      "Stored LOO was computed with unit = '", metadata[["unit"]],
      "'. Recompute with add_loo(object, unit = '", unit, "').",
      call. = FALSE
    )
  }

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

  if (!is.null(metadata) && metadata[["unit"]] != unit) {
    stop(
      "Stored WAIC was computed with unit = '", metadata[["unit"]],
      "'. Recompute with add_waic(object, unit = '", unit, "').",
      call. = FALSE
    )
  }

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

  return(invisible(NULL))
}
