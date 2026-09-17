# ============================================================================ #
# IWMDE Estimate Cache Helpers
# ============================================================================ #

.iwmde_estimate_cache <- function() {

  return(list(
    estimates   = new.env(parent = emptyenv()),
    diagnostics = .iwmde_diagnostic_cache()
  ))
}


.iwmde_estimate_cache_env <- function(cache) {

  if (is.list(cache) && is.environment(cache[["estimates"]])) {
    return(cache[["estimates"]])
  }

  return(NULL)
}


.iwmde_estimate_diagnostic_cache <- function(cache) {

  if (is.list(cache) && is.environment(cache[["diagnostics"]])) {
    return(cache[["diagnostics"]])
  }

  return(NULL)
}


.iwmde_estimate_cache_has <- function(cache, key) {

  cache_env <- .iwmde_estimate_cache_env(cache)
  return(!is.null(cache_env) && exists(key, envir = cache_env, inherits = FALSE))
}


.iwmde_estimate_cache_get <- function(cache, key) {

  cache_env <- .iwmde_estimate_cache_env(cache)
  return(get(key, envir = cache_env, inherits = FALSE))
}


.iwmde_estimate_cache_set <- function(cache, key, value) {

  cache_env <- .iwmde_estimate_cache_env(cache)
  if (!is.null(cache_env)) {
    assign(key, value, envir = cache_env)
  }

  return(invisible(value))
}
.iwmde_diagnostic_cache <- function() {

  cache <- new.env(parent = emptyenv())

  return(cache)
}


.iwmde_target_key <- function(parameter, parameter_spec) {

  condition_key <- .iwmde_parameter_condition_key(parameter_spec)
  structure_key <- unlist(lapply(c(
    target    = "target_columns",
    auxiliary = "auxiliary_columns",
    exclude   = "conditioning_exclude"
  ), function(field) {
    values <- parameter_spec[[field]]
    if (is.null(values) || length(values) == 0L) {
      return(character())
    }

    paste(sort(as.character(values)), collapse = ",")
  }), use.names = TRUE)
  if (length(structure_key) > 0L) {
    structure_key <- paste0(names(structure_key), "=", structure_key)
  }
  direction <- parameter_spec[["direction"]]
  if (!is.null(direction)) {
    direction <- direction[order(names(direction))]
    structure_key <- c(structure_key,
      paste0("direction:", names(direction), "=", .iwmde_key_number(direction)),
      paste0("chart:", parameter_spec[["conditioning_chart"]]))
  }
  if (!is.null(parameter_spec[["gate_metadata"]])) {
    structure_key <- c(structure_key, .iwmde_hash(
      "random_inclusion", parameter_spec[["gate_metadata"]]
    ))
  }

  if (identical(parameter_spec[["type"]], "linear")) {
    weights <- parameter_spec[["weights"]]
    weights <- weights[order(names(weights))]
    parts   <- paste0(
      names(weights),
      "=",
      .iwmde_key_number(weights)
    )

    return(paste(c("linear", parts, structure_key, condition_key), collapse = "|"))
  }

  if (identical(parameter_spec[["type"]], "simplex_pair")) {
    return(paste(c(
      "simplex_pair",
      parameter_spec[["parameter"]],
      parameter_spec[["index"]],
      parameter_spec[["n_targets"]],
      structure_key,
      condition_key
    ), collapse = "|"))
  }

  if (identical(parameter_spec[["type"]], "random_component_sd")) {
    factors <- vapply(parameter_spec[["factors"]], function(factor) {
      paste(
        factor[["weight_name"]],
        factor[["index"]],
        factor[["scale"]],
        factor[["n_targets"]],
        sep = ":"
      )
    }, character(1))

    return(paste(
      c(
        "random_component_sd",
        parameter_spec[["source_parameter"]],
        factors,
        structure_key,
        condition_key
      ),
      collapse = "|"
    ))
  }

  if (identical(parameter_spec[["type"]], "primitive")) {
    return(paste(
      c("primitive", parameter, structure_key, condition_key),
      collapse = "|"
    ))
  }

  if (identical(parameter_spec[["status"]], "unsupported")) {
    return(paste0("unsupported|", parameter))
  }

  return(paste0("primitive|", parameter))
}


.iwmde_parameter_condition_key <- function(parameter_spec) {

  condition_key <- parameter_spec[["condition_key"]]
  if (!is.null(condition_key) && length(condition_key) > 0L) {
    condition_key <- as.character(condition_key[[1L]])
    if (!is.na(condition_key) && nzchar(condition_key)) {
      return(condition_key)
    }
  }

  conditional <- parameter_spec[["conditional"]]
  if (is.null(conditional) || length(conditional) == 0L) {
    return(NULL)
  }

  conditional <- sort(unique(as.character(conditional)))
  rule        <- parameter_spec[["conditional_rule"]]
  if (is.null(rule) || length(rule) == 0L) {
    rule <- "OR"
  }

  return(paste0("conditional=", rule, ":", paste(conditional, collapse = ",")))
}


.iwmde_key_number <- function(x) {

  x <- as.numeric(x)
  out <- character(length(x))
  if (length(x) == 0L) {
    return(out)
  }

  # Chen proposal keys are built for every posterior row of a plan, so the
  # represented-coordinate encoding runs over the whole vector at once instead
  # of once per value. The encoded strings are the same.
  nan_values <- is.nan(x)
  na_values  <- is.na(x) & !nan_values
  infinite   <- is.infinite(x)
  out[nan_values] <- "NaN"
  out[na_values]  <- "NA"
  out[infinite & x > 0] <- "Inf"
  out[infinite & x < 0] <- "-Inf"

  finite <- !nan_values & !na_values & !infinite
  if (any(finite)) {
    values <- x[finite]
    # Negative zero is the same coordinate as zero.
    values[values == 0] <- 0
    bytes <- as.integer(writeBin(values, raw(), size = 8L, endian = "big"))
    hex   <- sprintf("%02x", bytes)
    dim(hex) <- c(8L, length(values))
    out[finite] <- do.call(
      paste0,
      lapply(seq_len(8L), function(byte) hex[byte, ])
    )
  }

  return(out)
}


.iwmde_cache_has <- function(cache, key) {

  return(!is.null(cache) && exists(key, envir = cache, inherits = FALSE))
}


.iwmde_cache_get <- function(cache, key) {

  return(get(key, envir = cache, inherits = FALSE))
}


.iwmde_cache_set <- function(cache, key, value) {

  if (!is.null(cache)) {
    assign(key, value, envir = cache)
  }

  return(invisible(value))
}
