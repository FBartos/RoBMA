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
# Re-export bridgesampling generics
# ---------------------------------------------------------------------------- #

#' @importFrom bridgesampling bridge_sampler
#' @export
bridgesampling::bridge_sampler

#' @importFrom bridgesampling logml
#' @export
bridgesampling::logml

#' @importFrom bridgesampling post_prob
#' @export
bridgesampling::post_prob

#' @importFrom bridgesampling bf
#' @export
bridgesampling::bf

#' @importFrom bridgesampling bayes_factor
#' @export
bridgesampling::bayes_factor


# ---------------------------------------------------------------------------- #
# S3 Methods for brma Objects
# ---------------------------------------------------------------------------- #

.log_sum_exp <- function(x) {

  max_x <- max(x)
  return(max_x + log(sum(exp(x - max_x))))
}

.post_prob_from_logml <- function(logml, prior_prob = NULL, model_names = NULL) {

  BayesTools::check_real(logml, "logml", check_length = 0, allow_NA = FALSE)
  if (any(!is.finite(logml))) {
    stop("'logml' must contain only finite values.", call. = FALSE)
  }

  if (length(logml) < 2L) {
    stop("At least two log marginal likelihoods are required.", call. = FALSE)
  }

  if (is.null(prior_prob)) {
    prior_prob <- rep(1 / length(logml), length(logml))
  } else {
    BayesTools::check_real(prior_prob, "prior_prob", lower = 0, check_length = 0, allow_NA = FALSE)
    if (any(!is.finite(prior_prob))) {
      stop("'prior_prob' must contain only finite values.", call. = FALSE)
    }
    if (length(prior_prob) != length(logml)) {
      stop("'prior_prob' must have the same length as the number of models.", call. = FALSE)
    }
    if (sum(prior_prob) <= 0) {
      stop("'prior_prob' must contain at least one positive value.", call. = FALSE)
    }
    prior_prob <- prior_prob / sum(prior_prob)
  }

  if (!is.null(model_names)) {
    BayesTools::check_char(model_names, "model_names", check_length = 0, allow_NA = FALSE)
    if (length(model_names) != length(logml)) {
      stop("'model_names' must have the same length as the number of models.", call. = FALSE)
    }
  }

  log_posterior <- logml + log(prior_prob)
  posterior     <- exp(log_posterior - .log_sum_exp(log_posterior))

  if (!is.null(model_names)) {
    names(posterior) <- model_names
  }

  return(posterior)
}

.bf_from_logml <- function(logml1, logml2, log = FALSE) {

  BayesTools::check_real(logml1, "logml1", check_length = 1, allow_NA = FALSE)
  BayesTools::check_real(logml2, "logml2", check_length = 1, allow_NA = FALSE)
  if (!is.finite(logml1) || !is.finite(logml2)) {
    stop("'logml1' and 'logml2' must be finite.", call. = FALSE)
  }
  BayesTools::check_bool(log, "log")

  log_bf <- logml1 - logml2
  out    <- list(
    bf  = if (log) log_bf else exp(log_bf),
    log = log
  )

  class(out) <- "bf_default"
  return(out)
}

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
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   fit <- add_marglik(fit)
#'
#'   bridge <- bridge_sampler(fit)
#'   print(bridge)
#' }
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
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   fit <- add_marglik(fit)
#'
#'   logml(fit)
#' }
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
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit1 <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit2 <- brma(
#'     yi           = yi,
#'     vi           = vi,
#'     data         = dat.lehmann2018,
#'     measure      = "SMD",
#'     prior_effect = FALSE
#'   )
#'
#'   fit1 <- add_marglik(fit1)
#'   fit2 <- add_marglik(fit2)
#'
#'   post_prob(fit1, fit2)
#' }
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

  return(.post_prob_from_logml(
    logml       = logml_values,
    prior_prob = prior_prob,
    model_names = model_names
  ))
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
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit1 <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit2 <- brma(
#'     yi           = yi,
#'     vi           = vi,
#'     data         = dat.lehmann2018,
#'     measure      = "SMD",
#'     prior_effect = FALSE
#'   )
#'
#'   fit1 <- add_marglik(fit1)
#'   fit2 <- add_marglik(fit2)
#'
#'   bf(fit1, fit2)
#' }
#' }
#'
#' @export
#' @exportS3Method bridgesampling::bf
bf.brma <- function(x1, x2, log = FALSE, ...) {
  if (!inherits(x2, "brma")) {
    stop("x2 needs to be of class 'brma'.", call. = FALSE)
  }
  .check_brma_compare_targets(list(x1, x2), "bf()")

  # compute log marginal likelihoods
  logml1 <- bridge_sampler.brma(x1)$logml
  logml2 <- bridge_sampler.brma(x2)$logml

  return(.bf_from_logml(
    logml1 = logml1,
    logml2 = logml2,
    log    = log
  ))
}


#' @rdname bf.brma
#' @export
#' @exportS3Method bridgesampling::bayes_factor
bayes_factor.brma <- function(x1, x2, log = FALSE, ...) {
  bf.brma(x1, x2, log = log, ...)
}
