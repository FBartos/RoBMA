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

  if (identical(parameter_spec[["type"]], "scale_log_intercept")) {
    return(paste(c(
      "scale_log_intercept",
      parameter,
      condition_key
    ), collapse = "|"))
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

  x <- zapsmall(as.numeric(x), digits = 14)
  x[abs(x) < sqrt(.Machine$double.eps)] <- 0

  return(formatC(x, digits = 15, format = "fg", flag = "#"))
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


.iwmde_select_active_rows <- function(rows, max_samples, context = NULL) {

  if (length(rows) <= max_samples) {
    return(rows)
  }

  selected <- NULL
  if (!is.null(context) && length(context[["indicator_names"]]) > 0L) {
    samples <- context[["posterior_samples"]]
    group <- vapply(rows, function(row) {
      .iwmde_active_key(context, samples[row, ])
    }, character(1))
    selected <- .thin_sample_rows_by_group(group, max_samples)
  }
  if (is.null(selected)) {
    selected <- .thin_sample_rows(length(rows), max_samples)
  }

  return(rows[selected])
}
