#' @title Fitting specification
#' @name fitting_specification
#'
#' @description
#' The `brma` family of functions uses the following arguments to specify the
#' MCMC sampling and fitting settings.
#'
#' @param sample numeric. Number of MCMC samples to save. Defaults to `5000`.
#' @param burnin numeric. Number of burn-in iterations. Defaults to `2000`.
#' @param adapt numeric. Number of adaptation iterations. Defaults to `500`.
#' @param chains numeric. Number of MCMC chains. Defaults to `3`.
#' @param thin numeric. Thinning interval. Defaults to `1`.
#' @param parallel logical. Whether to run MCMC chains in parallel. Defaults to `FALSE`.
#' @param autofit logical. Whether to automatically extend the MCMC chains if convergence is not met.
#' Defaults to `FALSE`.
#' @param autofit_control list of autofit control settings. See [set_autofit_control()] for details.
#' @param convergence_checks list of convergence check settings. See [set_convergence_checks()] for details.
#' @param seed numeric. Random seed for reproducibility. Defaults to `NULL`.
#' @param silent logical. Whether to suppress output. When omitted,
#' constructors use `RoBMA.get_option("silent")`.
#' @param ... additional advanced arguments. Fitting functions reject unused
#' arguments; currently recognized internal arguments include `only_data`,
#' `only_priors`, `is_JASP`, and `is_JASP_prefix`.
#'
#' @seealso \code{\link{brma}}, \code{\link{set_autofit_control}}, \code{\link{set_convergence_checks}}
#' @aliases fitting_specification
NULL

### creates basic RoBMA object with
# - fit_control
# - autofit_control
# - convergence_checks
#
# @param dots additional arguments originally passed as ...
# @param chains Number of MCMC chains
# @param adapt Number of adaptation iterations
# @param burnin Number of burnin iterations
# @param sample Number of sampling iterations
# @param thin Thinning interval
# @param autofit Whether to use autofit
# @param parallel Whether to run chains in parallel
# @param silent Whether to suppress output
# @param seed Random seed
# @param autofit_control Autofit control settings
# @param convergence_checks Convergence check settings
#
# @return A list containing fit_control, autofit_control, convergence_checks
.createObject <- function(
    dots, class,
    chains, adapt, burnin, sample, thin,
    autofit, parallel, silent, seed,
    autofit_control, convergence_checks) {

  object       <- NULL

  ### input global settings if unspecified
  if (missing(silent)) {
    silent <- RoBMA.get_option("silent")
  }

  ### check and store MCMC settings
  object$fit_control <- BayesTools::JAGS_check_and_list_fit_settings(
    chains = chains, adapt = adapt, burnin = burnin, sample = sample,
    thin = thin, autofit = autofit, parallel = parallel, cores = chains,
    silent = silent, seed = seed
  )
  object$autofit_control    <- BayesTools::JAGS_check_and_list_autofit_settings(autofit_control = autofit_control)
  object$convergence_checks <- .check_and_list_convergence_checks(convergence_checks = convergence_checks)


  ### include JASP indicators for progress bars
  if (!is.null(dots[["is_JASP"]])) {
    object[["is_JASP"]]        <- dots[["is_JASP"]]
    object[["is_JASP_prefix"]] <- dots[["is_JASP_prefix"]]
  }


  ### add class
  class(object) <- class

  return(object)
}

.set_only_priors_class <- function(object) {

  class(object) <- unique(c("only_priors.brma", class(object)))

  return(object)
}

.validate_constructor_dots <- function(dots, caller, allowed_extra = character()) {

  allowed <- c(
    "only_data", "only_priors", "is_JASP", "is_JASP_prefix",
    allowed_extra
  )
  .check_unused_dots(
    dots    = dots,
    allowed = allowed,
    caller  = caller
  )

  bool_arguments <- intersect(
    c("only_data", "only_priors", "is_JASP"),
    names(dots)
  )
  for (argument in bool_arguments) {
    BayesTools::check_bool(
      dots[[argument]],
      argument,
      allow_NA = FALSE
    )
  }

  if ("is_JASP_prefix" %in% names(dots)) {
    BayesTools::check_char(
      dots[["is_JASP_prefix"]],
      "is_JASP_prefix",
      check_length = 1,
      allow_NA     = FALSE
    )
  }

  return(dots)
}

.check_unused_dots <- function(dots, allowed, caller) {

  if (length(dots) == 0L) {
    return(invisible(TRUE))
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }

  unused <- dot_names[!nzchar(dot_names) | !dot_names %in% allowed]
  if (length(unused) == 0L) {
    return(invisible(TRUE))
  }

  unused[!nzchar(unused)] <- "<unnamed>"
  stop(
    "Unused argument", if (length(unused) > 1L) "s" else "",
    " in ", caller, ": ",
    paste0("'", unused, "'", collapse = ", "),
    call. = FALSE
  )
}

.warn_unused_dots <- function(dots, allowed, caller) {

  if (length(dots) == 0L) {
    return(invisible(TRUE))
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }

  unused <- dot_names[!nzchar(dot_names) | !dot_names %in% allowed]
  if (length(unused) == 0L) {
    return(invisible(TRUE))
  }

  unused[!nzchar(unused)] <- "<unnamed>"
  warning(
    "Unused argument", if (length(unused) > 1L) "s" else "",
    " in ", caller, ": ",
    paste0("'", unused, "'", collapse = ", "),
    call. = FALSE
  )

  return(invisible(TRUE))
}

.keep_allowed_dots <- function(dots, allowed) {

  if (length(dots) == 0L) {
    return(list())
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    return(list())
  }

  keep <- nzchar(dot_names) & dot_names %in% allowed

  return(dots[keep])
}

.autocompute_brma <- function(object, loo = TRUE, waic = TRUE,
                              marglik = !inherits(object, "RoBMA")) {

  if (loo && RoBMA.get_option("autocompute.loo")) {
    object <- add_loo(object)
  }
  if (waic && RoBMA.get_option("autocompute.waic")) {
    object <- add_waic(object)
  }
  if (marglik && RoBMA.get_option("autocompute.marglik")) {
    object <- add_marglik(object)
  }

  return(object)
}


# Validate model indicators without changing their represented values.
.as_exact_model_indicator <- function(value, name, scalar = FALSE) {

  if ((!is.numeric(value) && !is.integer(value)) ||
      (scalar && length(value) != 1L) ||
      any(!is.finite(value))) {
    stop("'", name, "' must contain finite model indicators.", call. = FALSE)
  }
  if (any(value != trunc(value))) {
    stop("'", name, "' must be integer-valued.", call. = FALSE)
  }
  if (any(abs(value) > .Machine$integer.max)) {
    stop("'", name, "' exceeds R's supported integer range.", call. = FALSE)
  }

  return(as.integer(value))
}


.extract_posterior_indicator <- function(posterior_samples, parameter,
                                         prior = NULL, column = NULL) {

  if (is.null(column)) {
    column <- paste0(parameter, "_indicator")
  }
  if (!column %in% colnames(posterior_samples)) {
    stop("Missing posterior model indicator: '", column, "'.",
         call. = FALSE)
  }

  indicator <- tryCatch(
    .as_exact_model_indicator(posterior_samples[, column], column),
    error = function(error) {
      stop(
        "Invalid posterior model indicator: '", column, "'. ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  if (!is.null(prior)) {
    valid_values <- seq_len(length(prior))
    if (any(!indicator %in% valid_values)) {
      stop("Invalid posterior model indicator range: '", column, "'.",
           call. = FALSE)
    }
  }

  return(indicator)
}


### object tools options
# add simple summary and model coefficients to the object
# (this differ from more customizable user facing summary function)
.object_summary      <- function(object) {

  # provide a simple summary
  estimates <- BayesTools::JAGS_estimates_table(
    fit               = object[["fit"]],
    transform_factors = TRUE,
    transform_scaled  = TRUE,
    remove_spike_0    = FALSE,
    remove_parameters = c(
      "theta", # remove random-effects (estimate-level)
      "gamma", # remove random-effects (cluster-level)
      "sampling_z", # remove known-V sampling dependency factors
      "pi",    # remove baserate for OR models
      "phi"    # remove lograte for IRR models
    ),
    random_effects_summary = "standard"
  )

  return(estimates)
}
.object_coefficients <- function(object) {

  estimates        <- object[["summary"]][,"Mean"]
  names(estimates) <- rownames(object[["summary"]])

  return(estimates)
}
