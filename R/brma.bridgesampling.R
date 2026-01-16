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
# S3 Methods for brma_bridge Objects
# ---------------------------------------------------------------------------- #
# These methods handle bridge objects returned by bridge_sampler.brma().
# They delegate to bridgesampling's internal methods.
# ---------------------------------------------------------------------------- #

#' @export
logml.brma_bridge <- function(x, ...) {

  bridgesampling:::logml.bridge(x, ...)
}

#' @export
post_prob.brma_bridge <- function(x, ...) {
  bridgesampling:::post_prob.bridge(x, ...)
}

#' @export
bf.brma_bridge <- function(x1, x2, log = FALSE, ...) {
  bridgesampling:::bf.bridge(x1, x2, log = log, ...)
}

#' @export
bayes_factor.brma_bridge <- function(x1, x2, log = FALSE, ...) {
  bf.brma_bridge(x1, x2, log = log, ...)
}


# ---------------------------------------------------------------------------- #
# S3 Methods for brma Objects
# ---------------------------------------------------------------------------- #

#' @title Bridge Sampling for brma Objects
#'
#' @description Compute the marginal likelihood of a brma model using
#' bridge sampling. The returned object can be used for Bayesian model
#' comparison via \code{\link{bf}} and \code{\link{post_prob}}.
#'
#' @param samples a brma model object.
#' @param ... additional arguments (currently not used).
#'
#' @details
#' The marginal likelihood is computed using the \code{bridgesampling} package
#' via the \code{BayesTools::JAGS_bridgesampling} wrapper.
#'
#' This is the primary method for computing marginal likelihoods from brma
#' objects. The convenience functions \code{\link{logml.brma}},
#' \code{\link{bf.brma}}, and \code{\link{post_prob.brma}} internally call
#' this method.
#'
#' @return An object of class \code{"bridge"} as returned by
#' \code{\link[bridgesampling]{bridge_sampler}}.
#'
#' @seealso \code{\link[bridgesampling]{bridge_sampler}}, \code{\link{logml.brma}},
#' \code{\link{bf.brma}}, \code{\link{post_prob.brma}}
#'
#' @examples
#' \dontrun{
#' fit <- brma(y = d, se = se, data = Bem2011)
#' bridge <- bridge_sampler(fit)
#' print(bridge)
#' }
#'
#' @export
#' @exportS3Method bridgesampling::bridge_sampler
bridge_sampler.brma <- function(samples, ...) {
  out <- .marglik(samples)
  class(out) <- c("brma_bridge", class(out))
  return(out)
}


#' @title Log Marginal Likelihood for brma Objects
#'
#' @description Extract the log marginal likelihood from a brma model
#' computed via bridge sampling.
#'
#' @param x a brma model object.
#' @param ... additional arguments (currently not used).
#'
#' @details
#' The marginal likelihood is computed using \code{\link{bridge_sampler.brma}}.
#'
#' @return A scalar numeric value representing the log marginal likelihood.
#'
#' @seealso \code{\link{bridge_sampler.brma}}, \code{\link{bf.brma}},
#' \code{\link{post_prob.brma}}
#'
#' @examples
#' \dontrun{
#' fit <- brma(y = d, se = se, data = Bem2011)
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
#' The marginal likelihoods are computed using \code{\link{bridge_sampler.brma}}.
#'
#' @return A named numeric vector with posterior model probabilities
#' (i.e., which sum to one).
#'
#' @seealso \code{\link{bridge_sampler.brma}}, \code{\link{bf.brma}},
#' \code{\link{logml.brma}}
#'
#' @examples
#' \dontrun{
#' fit1 <- brma(y = d, se = se, data = Bem2011)
#' fit2 <- brma(y = d, se = se, data = Bem2011,
#'              priors = prior(family = "point", value = 0, parameter = "mu"))
#' post_prob(fit1, fit2)
#' }
#'
#' @export
#' @exportS3Method bridgesampling::post_prob
post_prob.brma <- function(x, ..., prior_prob = NULL, model_names = NULL) {

  dots <- list(...)
  mc   <- match.call()

  # check that all objects in ... are brma objects
  modb <- vapply(dots, inherits, NA, what = "brma")
  if (sum(modb) == 0) {
    stop("Only one object of class 'brma' passed.", call. = FALSE)
  }
  if (sum(modb) != length(dots)) {
    warning("Objects not of class 'brma' are ignored.", call. = FALSE)
  }

  # get model names
  if (is.null(model_names)) {
    model_names <- c(
      deparse(mc[["x"]]),
      vapply(which(modb), function(i) deparse(mc[[i + 2]]), "")
    )
  }

  # compute log marginal likelihoods
  logml_values <- c(
    bridge_sampler.brma(x)$logml,
    vapply(dots[modb], function(obj) bridge_sampler.brma(obj)$logml, 0)
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
#' model \code{x2}. The marginal likelihoods are computed using
#' \code{\link{bridge_sampler.brma}}.
#'
#' @return A list of class \code{"bf_default"} with components:
#' \itemize{
#'   \item \code{bf}: (scalar) value of the Bayes factor in favor of
#'   \code{x1} over \code{x2}.
#'   \item \code{log}: Boolean indicating whether \code{bf} corresponds
#'   to the log Bayes factor.
#' }
#'
#' @seealso \code{\link{bridge_sampler.brma}}, \code{\link{post_prob.brma}},
#' \code{\link{logml.brma}}
#'
#' @examples
#' \dontrun{
#' fit1 <- brma(y = d, se = se, data = Bem2011)
#' fit2 <- brma(y = d, se = se, data = Bem2011,
#'              priors = prior(family = "point", value = 0, parameter = "mu"))
#' bf(fit1, fit2)
#' }
#'
#' @export
#' @exportS3Method bridgesampling::bf
bf.brma <- function(x1, x2, log = FALSE, ...) {

  if (!inherits(x2, "brma")) {
    stop("x2 needs to be of class 'brma'.", call. = FALSE)
  }

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
