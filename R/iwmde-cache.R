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

  if (identical(parameter_spec[["type"]], "linear")) {
    weights <- parameter_spec[["weights"]]
    weights <- weights[order(names(weights))]
    parts   <- paste0(
      names(weights),
      "=",
      .iwmde_key_number(weights)
    )

    return(paste(c("linear", parts, condition_key), collapse = "|"))
  }

  if (identical(parameter_spec[["type"]], "primitive")) {
    return(paste(c("primitive", parameter, condition_key), collapse = "|"))
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

  return(vapply(x, function(value) {

    if (is.nan(value)) {
      return("NaN")
    }
    if (is.na(value)) {
      return("NA")
    }
    if (is.infinite(value)) {
      return(if (value > 0) "Inf" else "-Inf")
    }
    if (value == 0) {
      value <- 0
    }

    bytes <- writeBin(value, raw(), size = 8L, endian = "big")
    return(paste(sprintf("%02x", as.integer(bytes)), collapse = ""))
  }, character(1), USE.NAMES = FALSE))
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


.iwmde_select_active_rows <- function(rows, max_samples) {

  if (length(rows) <= max_samples) {
    return(rows)
  }

  selected <- .thin_sample_rows(length(rows), max_samples)

  return(rows[selected])
}
