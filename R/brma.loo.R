# ============================================================================ #
# brma.loo.R
# ============================================================================ #
#
# This file implements LOO-PSIS (Leave-One-Out Cross-Validation via Pareto
# Smoothed Importance Sampling) diagnostics for brma class objects.
#
# The LOO is computed at the estimate level: each effect size estimate (or
# pair of counts for binomial/Poisson models) is treated as one observation.
# For multilevel models, the LOO conditions on the fitted study-level random
# effects.
#
# Important: LOO-PSIS evaluates how well the model predicts new observations.
# This is different from evaluating how well the model fits the observed data.
# Model comparison via loo_compare() selects models based on their expected
# out-of-sample predictive performance, not based on in-sample fit.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# add_loo generic and method
# ---------------------------------------------------------------------------- #

#' @export
add_loo <- function(object, ...) UseMethod("add_loo")


#' @title Add LOO-PSIS to brma Objects
#'
#' @description Compute approximate leave-one-out cross-validation (LOO-CV)
#' using Pareto smoothed importance sampling (PSIS) for brma model objects
#' and store the result in the object.
#'
#' @param object a brma model object.
#' @param r_eff optional vector of relative effective sample sizes. If not
#' provided, it is computed from the log-likelihood values.
#' @param parallel Logical. If \code{TRUE} computations are parallelized and
#' the number of cores is taken from \code{RoBMA.get_option("max_cores")}. If
#' \code{FALSE} computations are run on a single core.
#'
#' @details
#' LOO-CV is computed at the estimate level: for binomial and Poisson models,
#' each pair of counts (ai/ci or x1i/x2i) that defines a single effect size
#' estimate is treated as one observation.
#'
#' For multilevel models, the LOO conditions on the fitted study-level random
#' effects (gamma). This means the LOO evaluates prediction of new estimates
#' within studies, not prediction of new studies.
#'
#' For selection models, the LOO evaluates the weighted likelihood, conditioning
#' on the posterior omega samples.
#' 
#' The PSIS object is essential for model comparison via
#' \code{\link[loo]{loo_compare}} and is automatically saved in the loo result.
#'
#' \strong{Important for model comparison:} When comparing models via
#' \code{\link[loo]{loo_compare}}, the selection is based on expected
#' out-of-sample predictive performance. This evaluates how well models predict
#' \emph{new} observations, not how well they fit the observed data.
#'
#' @return The brma object with the LOO result stored in \code{object[["loo"]]}.
#'
#' @seealso \code{\link{loo.brma}}, \code{\link[loo]{loo}},
#' \code{\link[loo]{loo_compare}}, \code{\link[loo]{pareto_k_ids}}
#'
#' @examples \dontrun{
#' # Fit a brma model
#' fit <- brma(yi = yi, sei = sei, data = dat)
#'
#' # Compute and store LOO-PSIS
#' fit <- add_loo(fit)
#'
#' # Extract LOO object
#' loo_fit <- loo(fit)
#' print(loo_fit)
#'
#' # Check Pareto k diagnostics
#' plot(loo_fit)
#' }
#'
#' @references
#' \insertCite{vehtari2017practical}{RoBMA}
#' \insertCite{vehtari2024pareto}{RoBMA}
#'
#' @export
add_loo.brma <- function(object, r_eff = NULL, parallel = FALSE) {
  BayesTools::check_bool(parallel, "parallel")

  # compute the log-likelihood matrix (S x K)
  log_lik <- .pdf.brma(object)

  # determine number of cores based on `parallel` and package options
  cores <- if (parallel) max(1, RoBMA.get_option("max_cores")) else 1

  # compute relative effective sample sizes if not provided
  if (is.null(r_eff)) {
    # loo::relative_eff expects exp(log_lik) with chain_id for matrix input
    # Get chain / iteration information directly from the runjags fit
    n_chains <- length(object[["fit"]][["mcmc"]])
    n_iter <- object[["fit"]][["sample"]]

    chain_id <- rep(seq_len(n_chains), each = n_iter)

    r_eff <- loo::relative_eff(exp(log_lik), chain_id = chain_id, cores = cores)
  }

  # call loo on the log-likelihood matrix
  loo_result <- loo::loo(log_lik, r_eff = r_eff, save_psis = TRUE, cores = cores)

  # store in object and return

  object[["loo"]] <- loo_result
  return(object)
}


# ---------------------------------------------------------------------------- #
# add_waic generic and method
# ---------------------------------------------------------------------------- #

#' @export
add_waic <- function(object, ...) UseMethod("add_waic")


#' @title Add WAIC to brma Objects
#'
#' @description Compute the Widely Applicable Information Criterion (WAIC)
#' for brma model objects and store the result in the object.
#'
#' @param object a brma model object.
#' @param ... additional arguments passed to \code{\link[loo]{waic}}.
#'
#' @details
#' WAIC is an alternative to LOO-CV for estimating out-of-sample predictive
#' accuracy. Like LOO, it evaluates expected predictive performance for new
#' observations.
#'
#' In most cases, LOO-PSIS (via \code{\link{add_loo}}) is preferred over WAIC
#' because it provides better estimates and includes diagnostics (Pareto k
#' values) that indicate when the approximation may be unreliable.
#'
#' @return The brma object with the WAIC result stored in \code{object[["waic"]]}.
#'
#' @seealso \code{\link{waic.brma}}, \code{\link{add_loo}}, \code{\link[loo]{waic}}
#'
#' @export
add_waic.brma <- function(object, ...) {
  # compute the log-likelihood matrix (S x K)
  log_lik <- .pdf.brma(object)

  # call waic on the log-likelihood matrix
  waic_result <- loo::waic(log_lik, ...)

  # store in object and return
  object[["waic"]] <- waic_result
  return(object)
}


# ---------------------------------------------------------------------------- #
# loo generic and extraction method
# ---------------------------------------------------------------------------- #

#' @export
loo <- function(x, ...) UseMethod("loo")


#' @title LOO-PSIS for brma Objects
#'
#' @description Extract the LOO-PSIS object from a brma model object.
#' The LOO must first be computed using \code{\link{add_loo}}.
#'
#' @param x a brma model object.
#' @param ... additional arguments (currently unused).
#'
#' @details
#' This function extracts the LOO object that was previously computed and
#' stored using \code{\link{add_loo}}. If LOO has not been computed,
#' an error is thrown.
#'
#' @return An object of class \code{c("psis_loo", "loo")} as returned by
#' \code{\link[loo]{loo}}.
#'
#' @seealso \code{\link{add_loo}}, \code{\link[loo]{loo}},
#' \code{\link[loo]{loo_compare}}, \code{\link[loo]{pareto_k_ids}}
#'
#' @examples \dontrun{
#' # Fit a brma model
#' fit <- brma(yi = yi, sei = sei, data = dat)
#'
#' # Compute LOO-PSIS (must be called first)
#' fit <- add_loo(fit)
#'
#' # Extract LOO object
#' loo_fit <- loo(fit)
#' print(loo_fit)
#' }
#'
#' @export
loo.brma <- function(x, ...) {
  if (is.null(x[["loo"]])) {
    stop("LOO has not been computed. Call 'object <- add_loo(object)' first.",
      call. = FALSE
    )
  }

  return(x[["loo"]])
}


#' @title Extract Log-Likelihood Matrix from brma Object
#'
#' @description Extract the pointwise log-likelihood matrix from a brma model
#' object. This is an S x K matrix where S is the number of posterior samples
#' and K is the number of observations. This method implements the S3
#' \'logLik\' generic for \code{brma} objects and returns the matrix of
#' pointwise log-likelihoods (one column per observation, one row per sample).
#'
#' @param object a brma model object.
#' @param ... currently unused.
#'
#' @details
#' The log-likelihood is computed for each observation at each posterior sample.
#' For binomial and Poisson models, each observation consists of a pair of
#' counts (ai/ci or x1i/x2i) that together define a single effect size estimate.
#'
#' @return An S x K matrix of log-likelihood values.
#'
#' @seealso \code{\link{loo.brma}}
#'
#' @export
logLik.brma <- function(object, ...) {
  out <- .pdf.brma(object)
  class(out) <- c("logLik.brma", class(out))
  return(out)
}


#' @export
print.logLik.brma <- function(x, ...) {
  S <- if (is.matrix(x)) nrow(x) else length(x)
  K <- if (is.matrix(x)) ncol(x) else 1
  cat(sprintf("%d*%d pointwise log-likelihood matrix\n", S, K))
  invisible(x)
}


# ---------------------------------------------------------------------------- #
# loo_compare generic and methods
# ---------------------------------------------------------------------------- #

#' @export
loo_compare <- function(x, ...) UseMethod("loo_compare")


#' @title Compare brma Models Using LOO
#'
#' @description Compare multiple brma models using LOO-PSIS cross-validation.
#' This is a convenience wrapper around \code{\link[loo]{loo_compare}}.
#'
#' @param x a brma model object (the first model to compare).
#' @param ... additional brma model objects or \code{loo} objects to compare.
#'
#' @details
#' This function compares models based on their expected out-of-sample
#' predictive performance (ELPD).
#'
#' \strong{Important for model comparison:} When comparing models via
#' \code{\link[loo]{loo_compare}}, the selection is based on expected
#' out-of-sample predictive performance. This evaluates how well models predict
#' \emph{new} observations, not how well they fit the observed data.
#'
#' @return A matrix of class \code{"compare.loo"} as returned by
#' \code{\link[loo]{loo_compare}}.
#'
#' @seealso \code{\link{loo.brma}}, \code{\link[loo]{loo_compare}}
#'
#' @examples \dontrun{
#' # Fit models with and without publication bias adjustment
#' fit_bias <- brma(yi = yi, sei = sei, data = dat)
#' fit_nobias <- brma(yi = yi, sei = sei, data = dat, priors_bias = NULL)
#'
#' # Add LOO to both models
#' fit_bias <- add_loo(fit_bias)
#' fit_nobias <- add_loo(fit_nobias)
#'
#' # Compare models based on expected out-of-sample predictive performance
#' loo_compare(fit_bias, fit_nobias)
#'
#' # Alternatively, compare loo objects directly
#' loo_compare(loo(fit_bias), loo(fit_nobias))
#' }
#'
#' @export
loo_compare.brma <- function(x, ...) {
  # collect all models: x plus any in ...
  models <- c(list(x), list(...))

  if (length(models) < 2) {
    stop("At least two models are required for comparison.", call. = FALSE)
  }

  # convert brma objects to loo objects if necessary
  loo_objects <- lapply(models, function(m) {
    if (inherits(m, "brma")) {
      loo.brma(m)
    } else if (inherits(m, "loo")) {
      m
    } else {
      stop("All arguments must be brma or loo objects.", call. = FALSE)
    }
  })

  # call the 'loo' package's default implementation to avoid dispatch recursion
  loo_compare_fun <- get("loo_compare.default", envir = asNamespace("loo"), inherits = FALSE)
  result <- do.call(loo_compare_fun, loo_objects)

  return(result)
}


#' @title Compare loo Objects Using LOO
#'
#' @description Method for comparing loo objects directly.
#'
#' @param x a loo object (the first model to compare).
#' @param ... additional loo or brma objects to compare.
#'
#' @return A matrix of class \code{"compare.loo"} as returned by
#' \code{\link[loo]{loo_compare}}.
#'
#' @seealso \code{\link{loo.brma}}, \code{\link[loo]{loo_compare}}
#'
#' @export
loo_compare.loo <- function(x, ...) {
  # collect all models: x plus any in ...
  models <- c(list(x), list(...))

  if (length(models) < 2) {
    stop("At least two models are required for comparison.", call. = FALSE)
  }

  # convert brma objects to loo objects if necessary
  loo_objects <- lapply(models, function(m) {
    if (inherits(m, "brma")) {
      loo.brma(m)
    } else if (inherits(m, "loo")) {
      m
    } else {
      stop("All arguments must be brma or loo objects.", call. = FALSE)
    }
  })

  # call the 'loo' package's default implementation to avoid dispatch recursion
  loo_compare_fun <- get("loo_compare.default", envir = asNamespace("loo"), inherits = FALSE)
  result <- do.call(loo_compare_fun, loo_objects)

  return(result)
}


# ---------------------------------------------------------------------------- #
# waic generic and extraction method
# ---------------------------------------------------------------------------- #

#' @export
waic <- function(x, ...) UseMethod("waic")


#' @title WAIC for brma Objects
#'
#' @description Extract the WAIC object from a brma model object.
#' The WAIC must first be computed using \code{\link{add_waic}}.
#'
#' @param x a brma model object.
#' @param ... additional arguments (currently unused).
#'
#' @details
#' This function extracts the WAIC object that was previously computed and
#' stored using \code{\link{add_waic}}. If WAIC has not been computed,
#' an error is thrown.
#'
#' In most cases, LOO-PSIS (via \code{\link{loo.brma}}) is preferred over WAIC
#' because it provides better estimates and includes diagnostics (Pareto k
#' values) that indicate when the approximation may be unreliable.
#'
#' @return An object of class \code{"waic"} as returned by
#' \code{\link[loo]{waic}}.
#'
#' @seealso \code{\link{add_waic}}, \code{\link{loo.brma}}, \code{\link[loo]{waic}}
#'
#' @export
waic.brma <- function(x, ...) {
  if (is.null(x[["waic"]])) {
    stop("WAIC has not been computed. Call 'object <- add_waic(object)' first.",
      call. = FALSE
    )
  }

  return(x[["waic"]])
}
