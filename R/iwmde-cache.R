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
