#' @title Update a brma Fit
#'
#' @description
#' Extends an existing fitted \code{brma} object by additional MCMC samples,
#' updates study labels, and optionally recomputes cached fit-dependent
#' quantities.
#'
#' @param object a fitted \code{brma} object.
#' @param formula. unsupported; included for compatibility with
#'   \code{\link[stats]{update}}.
#' @param ... unsupported additional arguments.
#' @param sample_extend integer. Number of additional samples per chain.
#' @param slab optional character vector of study labels. Updating labels does
#'   not refit or extend the model.
#' @param autofit_control list of autofit control settings. Values are merged
#'   with the existing settings before extending.
#' @param convergence_checks list of convergence check settings. Values are
#'   merged with the existing settings and used to re-check the fit.
#' @param recompute whether cached \code{loo}, \code{waic}, and
#'   \code{marglik} values already stored in \code{object} should be recomputed
#'   after extension (\code{"all"}) or dropped with a warning (\code{"drop"}).
#' @param parallel logical. Whether to extend chains in parallel.
#' @param cores integer. Number of cores to use when \code{parallel = TRUE}.
#' @param silent logical. Whether to suppress JAGS output during extension.
#' @param seed optional seed used before extending.
#' @param evaluate unsupported; included for compatibility with
#'   \code{\link[stats]{update}}.
#'
#' @return The updated \code{brma} object.
#'
#' @details
#' Extending a fit adds posterior samples only. It does not rerun adaptation or
#' burn-in. Prior, data, and model-structure updates are intentionally not
#' supported by this method. For \code{brma.mv()} objects, the stored \code{V},
#' \code{R}, known-\code{V} backend, and marginalized random-effect metadata are
#' preserved because the fitted model structure is reused.
#'
#' @examples \dontrun{
#' fit <- update(fit, sample_extend = 1000)
#' fit <- update(fit, slab = paste("Study", seq_len(nobs(fit))))
#' }
#'
#' @export
update.brma <- function(
    object, formula. = NULL, ...,
    sample_extend = NULL, slab = NULL,
    autofit_control = NULL, convergence_checks = NULL,
    recompute = c("all", "drop"),
    parallel = NULL, cores = NULL, silent = NULL, seed = NULL,
    evaluate = TRUE) {

  if (!inherits(object, "brma")) {
    stop("'object' must be a 'brma' object.", call. = FALSE)
  }
  if (!is.null(formula.)) {
    stop("Updating formulas is not supported by update.brma().", call. = FALSE)
  }
  if (!isTRUE(evaluate)) {
    stop("update.brma() does not support 'evaluate = FALSE'.", call. = FALSE)
  }

  dots <- list(...)
  .check_unused_dots(
    dots    = dots,
    allowed = character(0),
    caller  = "update.brma()"
  )

  recompute <- match.arg(recompute)

  if (!is.null(sample_extend)) {
    BayesTools::check_int(
      sample_extend,
      "sample_extend",
      lower     = 1,
      allow_NA  = FALSE
    )
  }

  if (!is.null(slab)) {
    object <- .update_brma_slab(object = object, slab = slab)
  }

  object[["fit_control"]] <- .update_fit_control(
    old_fit_control = object[["fit_control"]],
    chains          = NULL,
    adapt           = NULL,
    burnin          = NULL,
    sample          = NULL,
    thin            = NULL,
    autofit         = NULL,
    parallel        = parallel,
    cores           = cores,
    silent          = silent,
    seed            = seed
  )
  object[["autofit_control"]] <- .update_autofit_control(
    old_autofit_control = object[["autofit_control"]],
    autofit_control     = .update_brma_autofit_input(
      autofit_control = autofit_control,
      sample_extend   = sample_extend
    )
  )
  object[["convergence_checks"]] <- .update_convergence_checks(
    old_convergence_checks = object[["convergence_checks"]],
    convergence_checks     = convergence_checks
  )

  fit_extended <- !is.null(sample_extend)
  cached       <- .collect_brma_fit_cache(object)

  if (fit_extended) {
    object <- .extend_brma_fit_once(object)
    object[["summary"]]      <- .object_summary(object)
    object[["coefficients"]] <- .object_coefficients(object)
    object <- .refresh_brma_fit_cache(
      object    = object,
      cached    = cached,
      recompute = recompute
    )
  } else if (!is.null(convergence_checks)) {
    object[["fit"]] <- .recheck_brma_fit(object)
  }

  return(object)
}


# Update study labels without touching the posterior fit.
.update_brma_slab <- function(object, slab) {

  outcome <- object[["data"]][["outcome"]]

  if (length(slab) != nrow(outcome)) {
    stop(
      "The 'slab' argument must have length ",
      nrow(outcome),
      " (same as the fitted data).",
      call. = FALSE
    )
  }

  object[["data"]][["outcome"]][["slab"]] <- as.character(slab)
  attr(object[["data"]], "slab")          <- TRUE

  return(object)
}


# Merge update-time autofit settings before validation.
.update_brma_autofit_input <- function(autofit_control, sample_extend) {

  if (is.null(autofit_control)) {
    autofit_control <- list()
  }
  if (!is.null(sample_extend)) {
    autofit_control[["sample_extend"]] <- sample_extend
  }

  return(autofit_control)
}


# Extend exactly one sample_extend chunk, even when stored max_extend is larger.
.extend_brma_fit_once <- function(object) {

  if (is.null(object[["fit"]]) || length(object[["fit"]]) == 0L) {
    stop("'object' does not contain a fitted model to extend.", call. = FALSE)
  }

  stored_autofit_control <- object[["autofit_control"]]
  extend_autofit_control <- stored_autofit_control
  extend_autofit_control[["max_extend"]] <- 1
  extend_autofit_control <- BayesTools::JAGS_check_and_list_autofit_settings(
    autofit_control = extend_autofit_control
  )

  object[["autofit_control"]] <- extend_autofit_control
  object[["fit"]]             <- .fit(object, extend = TRUE)
  object[["autofit_control"]] <- stored_autofit_control
  .stop_fit_errors(object[["fit"]])

  return(object)
}


# Re-run convergence checks on the current posterior samples.
.recheck_brma_fit <- function(object) {

  fit <- object[["fit"]]

  if (is.null(fit) || length(fit) == 0L) {
    stop("'object' does not contain a fitted model to check.", call. = FALSE)
  }

  prior_list <- attr(fit, "prior_list")
  if (is.null(prior_list)) {
    prior_list <- .create_fit_priors(
      data   = object[["data"]],
      priors = object[["priors"]]
    )
  }

  check_fit <- BayesTools::JAGS_check_convergence(
    fit            = fit,
    prior_list     = prior_list,
    add_parameters = .convergence_structural_parameters(object[["priors"]]),
    max_Rhat             = object[["convergence_checks"]][["max_Rhat"]],
    min_ESS              = object[["convergence_checks"]][["min_ESS"]],
    max_error            = object[["convergence_checks"]][["max_error"]],
    max_SD_error         = object[["convergence_checks"]][["max_SD_error"]],
    check_indicators     = isTRUE(object[["convergence_checks"]][["check_indicators"]]),
    monitor              = object[["convergence_checks"]][["monitor"]],
    allow_not_assessable = isTRUE(object[["convergence_checks"]][["allow_not_assessable"]])
  )

  fit[["converged"]]     <- check_fit
  fit[["has_posterior"]] <- TRUE
  fit[["warnings"]]      <- c(attr(fit, "warnings"), attr(check_fit, "errors"))

  return(fit)
}


# Record which fit-dependent quantities need cache handling after extension.
.collect_brma_fit_cache <- function(object) {

  cached <- list(
    loo     = object[["loo"]],
    waic    = object[["waic"]],
    marglik = object[["marglik"]]
  )

  return(cached)
}


# Recompute or drop fit-dependent cached quantities after extension.
.refresh_brma_fit_cache <- function(object, cached, recompute) {

  cached_names <- .brma_cached_names(cached)

  if (length(cached_names) == 0L) {
    return(object)
  }

  object <- .drop_brma_fit_cache(object)

  if (recompute == "drop") {
    warning(
      "Dropping cached ",
      paste(cached_names, collapse = ", "),
      " because the fit was extended.",
      call. = FALSE
    )
    return(object)
  }

  object <- .recompute_brma_fit_cache(
    object = object,
    cached = cached
  )

  return(object)
}


# Names of cached fit-dependent quantities present in the original object.
.brma_cached_names <- function(cached) {

  cached_names <- character(0)

  if (!is.null(cached[["loo"]])) {
    cached_names <- c(cached_names, "loo")
  }
  if (!is.null(cached[["waic"]])) {
    cached_names <- c(cached_names, "waic")
  }
  if (!is.null(cached[["marglik"]])) {
    cached_names <- c(cached_names, "marglik")
  }

  return(cached_names)
}


# Remove cached fit-dependent quantities.
.drop_brma_fit_cache <- function(object) {

  object[["loo"]]     <- NULL
  object[["waic"]]    <- NULL
  object[["marglik"]] <- NULL

  return(object)
}


# Recompute the same cached quantities and units that existed before extension.
.recompute_brma_fit_cache <- function(object, cached) {

  loo_units <- .brma_cache_units(cached[["loo"]])
  for (unit in loo_units) {
    object <- add_loo(object, unit = unit)
  }

  waic_units <- .brma_cache_units(cached[["waic"]])
  for (unit in waic_units) {
    object <- add_waic(object, unit = unit)
  }

  if (!is.null(cached[["marglik"]]) && !inherits(object, "RoBMA")) {
    object <- add_marglik(object)
  }

  return(object)
}


# Extract cached LOO/WAIC unit names, accepting a single unnamed legacy cache.
.brma_cache_units <- function(cache) {

  if (is.null(cache)) {
    return(character(0))
  }

  units <- names(cache)
  if (is.null(units)) {
    if (length(cache) == 1L) {
      return("estimate")
    }
    stop("Cached LOO/WAIC objects must be named by unit.", call. = FALSE)
  }

  units <- units[nzchar(units)]

  return(units)
}
