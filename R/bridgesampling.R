# ============================================================================ #
# brma.bridgesampling.R
# ============================================================================ #
#
# This file implements S3 methods for bridgesampling package generics
# for brma class objects. The marginal likelihood computed via bridge
# sampling is used for Bayesian model comparison via Bayes factors.
#
# The implementation follows the same pattern as brma.loo.R and brma.as_draws.R:
# generics are defined locally so they work regardless of whether the
# bridgesampling package is loaded first, last, or not at all.
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# Generics
# ---------------------------------------------------------------------------- #
# These generics mirror the bridgesampling package functions. They enable
# the S3 dispatch to work regardless of package load order.
# ---------------------------------------------------------------------------- #

#' @export
bridge_sampler <- function(samples, ...) UseMethod("bridge_sampler")

#' @export
logml <- function(x, ...) UseMethod("logml")

#' @export
post_prob <- function(x, ...) UseMethod("post_prob")

#' @export
bf <- function(x1, x2, log = FALSE, ...) UseMethod("bf")

#' @export
bayes_factor <- function(x1, x2, log = FALSE, ...) UseMethod("bayes_factor")


# ---------------------------------------------------------------------------- #
# S3 Methods for brma Objects
# ---------------------------------------------------------------------------- #

#' @title Bridge Sampling for brma Objects
#'
#' @description Extract the marginal likelihood bridge sampling object from
#' a brma model. The marginal likelihood must first be computed using
#' \code{\link{add_marglik}}.
#'
#' @param samples a brma model object.
#' @param ... additional arguments (currently not used).
#'
#' @details
#' This function extracts the bridge sampling object that was previously
#' computed and stored using \code{\link{add_marglik}}. If the marginal
#' likelihood has not been computed, an error is thrown.
#'
#' The returned object can be used for Bayesian model comparison via
#' \code{\link{bf}} and \code{\link{post_prob}}.
#'
#' @return An object of class \code{"bridge"} as returned by
#' \code{\link[bridgesampling]{bridge_sampler}}.
#'
#' @seealso \code{\link{add_marglik}}, \code{\link[bridgesampling]{bridge_sampler}},
#' \code{\link{logml.brma}}, \code{\link{bf.brma}}, \code{\link{post_prob.brma}}
#'
#' @examples
#' \dontrun{
#' fit <- brma(y = d, se = se, data = Bem2011)
#'
#' # Compute marginal likelihood first
#' fit <- add_marglik(fit)
#'
#' # Extract bridge sampling object
#' bridge <- bridge_sampler(fit)
#' print(bridge)
#' }
#'
#' @export
#' @exportS3Method bridgesampling::bridge_sampler
bridge_sampler.brma <- function(samples, ...) {
  if (is.null(samples[["marglik"]])) {
    stop("Marginal likelihood has not been computed. Call 'object <- add_marglik(object)' first.",
      call. = FALSE
    )
  }

  return(samples[["marglik"]])
}


#' @title Log Marginal Likelihood for brma Objects
#'
#' @description Extract the log marginal likelihood from a brma model.
#' The marginal likelihood must first be computed using \code{\link{add_marglik}}.
#'
#' @param x a brma model object.
#' @param ... additional arguments (currently not used).
#'
#' @details
#' This function extracts the log marginal likelihood from the bridge sampling
#' object that was previously computed and stored using \code{\link{add_marglik}}.
#'
#' @return A scalar numeric value representing the log marginal likelihood.
#'
#' @seealso \code{\link{add_marglik}}, \code{\link{bridge_sampler.brma}},
#' \code{\link{bf.brma}}, \code{\link{post_prob.brma}}
#'
#' @examples
#' \dontrun{
#' fit <- brma(y = d, se = se, data = Bem2011)
#'
#' # Compute marginal likelihood first
#' fit <- add_marglik(fit)
#'
#' # Get log marginal likelihood
#' logml(fit)
#' }
#'
#' @export
#' @exportS3Method bridgesampling::logml
logml.brma <- function(x, ...) {
  bridge <- bridge_sampler.brma(x)
  return(bridge$logml)
}


#' @title Posterior Model Probabilities for brma Objects
#'
#' @description Compute posterior model probabilities from marginal
#' likelihoods of brma models.
#'
#' @param x a brma model object.
#' @param ... additional brma model objects.
#' @param prior_prob numeric vector with prior model probabilities.
#' If omitted, a uniform prior is used (i.e., all models are equally
#' likely a priori). The default \code{NULL} corresponds to equal
#' prior model weights.
#' @param model_names character vector with model names. If \code{NULL}
#' (the default), names will be derived from deparsing the call.
#'
#' @details
#' The marginal likelihoods must first be computed using \code{\link{add_marglik}}.
#' All models must be fitted to the same \code{yi}/\code{sei} target.
#'
#' @return A named numeric vector with posterior model probabilities
#' (i.e., which sum to one).
#'
#' @seealso \code{\link{add_marglik}}, \code{\link{bridge_sampler.brma}},
#' \code{\link{bf.brma}}, \code{\link{logml.brma}}
#'
#' @examples
#' \dontrun{
#' fit1 <- brma(y = d, se = se, data = Bem2011)
#' fit2 <- brma(
#'   y = d, se = se, data = Bem2011,
#'   priors = prior(family = "point", value = 0, parameter = "mu")
#' )
#'
#' # Compute marginal likelihoods
#' fit1 <- add_marglik(fit1)
#' fit2 <- add_marglik(fit2)
#'
#' post_prob(fit1, fit2)
#' }
#'
#' @export
#' @exportS3Method bridgesampling::post_prob
post_prob.brma <- function(x, ..., prior_prob = NULL, model_names = NULL) {
  dots <- list(...)
  mc <- match.call()

  # check that all objects in ... are brma objects
  modb <- vapply(dots, inherits, NA, what = "brma")
  if (sum(modb) == 0) {
    stop("Only one object of class 'brma' passed.", call. = FALSE)
  }
  if (sum(modb) != length(dots)) {
    warning("Objects not of class 'brma' are ignored.", call. = FALSE)
  }
  models <- c(list(x), dots[modb])
  .check_brma_compare_targets(models, "post_prob()")

  # get model names
  if (is.null(model_names)) {
    model_names <- c(
      deparse(mc[["x"]]),
      vapply(which(modb), function(i) deparse(mc[[i + 2]]), "")
    )
  }

  # compute log marginal likelihoods (this will error if add_marglik not called)
  logml_values <- c(
    bridge_sampler.brma(x)$logml,
    vapply(models[-1], function(obj) bridge_sampler.brma(obj)$logml, 0)
  )

  # delegate to bridgesampling's implementation
  bridgesampling:::post_prob.default(
    logml_values,
    prior_prob  = prior_prob,
    model_names = model_names
  )
}


#' @title Bayes Factor for brma Objects
#'
#' @description Compute the Bayes factor comparing two brma models.
#'
#' @param x1 a brma model object (numerator).
#' @param x2 a brma model object (denominator).
#' @param log logical; if \code{TRUE}, the log Bayes factor is returned.
#' Default is \code{FALSE}.
#' @param ... additional arguments (currently not used).
#'
#' @details
#' Computes the Bayes factor in favor of the model \code{x1} over the
#' model \code{x2}. The marginal likelihoods must first be computed using
#' \code{\link{add_marglik}}. Both models must be fitted to the same
#' \code{yi}/\code{sei} target.
#'
#' @return A list of class \code{"bf_default"} with components:
#' \itemize{
#'   \item \code{bf}: (scalar) value of the Bayes factor in favor of
#'   \code{x1} over \code{x2}.
#'   \item \code{log}: Boolean indicating whether \code{bf} corresponds
#'   to the log Bayes factor.
#' }
#'
#' @seealso \code{\link{add_marglik}}, \code{\link{bridge_sampler.brma}},
#' \code{\link{post_prob.brma}}, \code{\link{logml.brma}}
#'
#' @examples
#' \dontrun{
#' fit1 <- brma(y = d, se = se, data = Bem2011)
#' fit2 <- brma(
#'   y = d, se = se, data = Bem2011,
#'   priors = prior(family = "point", value = 0, parameter = "mu")
#' )
#'
#' # Compute marginal likelihoods
#' fit1 <- add_marglik(fit1)
#' fit2 <- add_marglik(fit2)
#'
#' bf(fit1, fit2)
#' }
#'
#' @export
#' @exportS3Method bridgesampling::bf
bf.brma <- function(x1, x2, log = FALSE, ...) {
  if (!inherits(x2, "brma")) {
    stop("x2 needs to be of class 'brma'.", call. = FALSE)
  }
  .check_brma_compare_targets(list(x1, x2), "bf()")

  # compute log marginal likelihoods and delegate to bridgesampling
  logml1 <- bridge_sampler.brma(x1)$logml
  logml2 <- bridge_sampler.brma(x2)$logml
  bridgesampling:::bf.default(logml1, logml2, log = log)
}


#' @rdname bf.brma
#' @export
#' @exportS3Method bridgesampling::bayes_factor
bayes_factor.brma <- function(x1, x2, log = FALSE, ...) {
  bf.brma(x1, x2, log = log, ...)
}
