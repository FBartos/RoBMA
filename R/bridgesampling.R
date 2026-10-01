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

  prior_prob <- .post_prob_normalize_prior(prior_prob, length(logml))
  .post_prob_check_names(model_names, length(logml))

  log_posterior <- logml + log(prior_prob)
  posterior     <- exp(log_posterior - .log_sum_exp(log_posterior))

  if (!is.null(model_names)) {
    names(posterior) <- model_names
  }

  return(posterior)
}

.post_prob_normalize_prior <- function(prior_prob, n_models) {

  if (is.null(prior_prob)) {
    prior_prob <- rep(1 / n_models, n_models)
  } else {
    BayesTools::check_real(prior_prob, "prior_prob", lower = 0, check_length = 0, allow_NA = FALSE)
    if (any(!is.finite(prior_prob))) {
      stop("'prior_prob' must contain only finite values.", call. = FALSE)
    }
    if (length(prior_prob) != n_models) {
      stop("'prior_prob' must have the same length as the number of models.", call. = FALSE)
    }
    if (sum(prior_prob) <= 0) {
      stop("'prior_prob' must contain at least one positive value.", call. = FALSE)
    }
    prior_prob <- prior_prob / sum(prior_prob)
  }
  return(prior_prob)
}

.post_prob_check_names <- function(model_names, n_models) {

  if (!is.null(model_names)) {
    BayesTools::check_char(model_names, "model_names", check_length = 0, allow_NA = FALSE)
    if (length(model_names) != n_models) {
      stop("'model_names' must have the same length as the number of models.", call. = FALSE)
    }
  }
  return(invisible(NULL))
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
#' This function extracts the raw \pkg{bridgesampling} object underlying the
#' marginal-likelihood result stored by \code{\link{add_marglik}}. If the
#' marginal likelihood has not been computed, an error is thrown.
#' Product-space model-averaging objects (\code{BMA.norm}, \code{BMA.glmm},
#' \code{BMA.mv}, \code{RoBMA}, and \code{RoBMA.mv}) do not expose
#' bridge-sampling marginal likelihoods.
#'
#' Fully fixed models can have a zero-dimensional marginal likelihood evaluated
#' exactly. Such fits have no bridge-sampling object, so this method raises a
#' \code{RoBMA_exact_marglik_no_bridge} error; use \code{\link{logml}} to
#' retrieve their exact log marginal likelihood. For Bayesian model comparison,
#' pass fitted \code{brma} objects directly to \code{\link{bf}} or
#' \code{\link{post_prob}}.
#'
#' @return An object of class \code{"bridge"}, or \code{"bridge_list"} for
#' repeated bridge estimates, as returned by
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
  marglik <- .brma_stored_marglik(samples)
  if (inherits(marglik, c("bridge", "bridge_list"))) {
    return(marglik)
  }

  upstream <- marglik[["diagnostics"]][["upstream"]]
  if (!inherits(upstream, c("bridge", "bridge_list"))) {
    .brma_stop_bridge_unavailable(marglik)
  }

  target <- attr(marglik, "RoBMA_target", exact = TRUE)
  if (!is.null(target)) {
    attr(upstream, "RoBMA_target") <- target
  }

  return(upstream)
}


.brma_stored_marglik <- function(object) {

  if (inherits(object, "RoBMA")) {
    .stop_product_space_marglik()
  }

  marglik <- object[["marglik"]]
  if (is.null(marglik)) {
    stop(
      "Marginal likelihood has not been computed. Call ",
      "'object <- add_marglik(object)' first.",
      call. = FALSE
    )
  }

  return(marglik)
}


.brma_stop_bridge_unavailable <- function(marglik) {

  is_exact <- identical(
    marglik[["aggregation"]][["rule"]],
    "exact_zero_dimensional"
  )
  message <- if (is_exact) {
    paste0(
      "No bridge-sampling object exists because the zero-dimensional ",
      "marginal likelihood was evaluated exactly. Use 'logml()' to retrieve ",
      "the exact log marginal likelihood."
    )
  } else {
    paste0(
      "The stored marginal-likelihood result does not contain a raw ",
      "bridge-sampling object. Recompute it with 'add_marglik()'."
    )
  }
  condition_class <- if (is_exact) {
    c(
      "RoBMA_exact_marglik_no_bridge",
      "RoBMA_bridge_unavailable",
      "error",
      "condition"
    )
  } else {
    c("RoBMA_bridge_unavailable", "error", "condition")
  }
  condition <- structure(
    list(
      message = message,
      call    = NULL,
      reason  = if (is_exact) "exact_zero_dimensional" else "missing_upstream",
      logml   = marglik[["logml"]]
    ),
    class = condition_class
  )

  stop(condition)
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
#' Product-space model-averaging objects (\code{BMA.norm}, \code{BMA.glmm},
#' \code{BMA.mv}, \code{RoBMA}, and \code{RoBMA.mv}) do not expose
#' bridge-sampling marginal likelihoods.
#'
#' For repeated standardized results, this method returns the stored scalar
#' median of the included log marginal likelihoods. Legacy raw bridge lists
#' use the upstream scalar median method. Use \code{\link{bridge_sampler}} to
#' extract the raw bridge repetitions, or \code{\link{bf}} and
#' \code{\link{post_prob}} to compare models across repetitions.
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
  marglik <- .brma_stored_marglik(x)
  if (inherits(marglik, "bridge_list")) {
    return(bridgesampling::logml(marglik))
  }
  return(marglik$logml)
}


.brma_comparison_logml <- function(object) {

  marglik <- .brma_stored_marglik(object)
  repetitions <- marglik[["repetitions"]]
  upstream <- marglik[["diagnostics"]][["upstream"]]
  legacy <- inherits(marglik, "bridge_list")
  if (is.data.frame(repetitions) && nrow(repetitions) > 0L) {
    values <- repetitions[["logml"]]
  } else if (inherits(upstream, "bridge_list")) {
    values <- upstream[["logml"]]
  } else {
    values <- marglik[["logml"]]
  }
  BayesTools::check_real(values, "logml", check_length = 0, allow_NA = TRUE)
  repeated <- legacy || inherits(upstream, "bridge_list") || length(values) > 1L

  if (!legacy && any(!is.finite(values))) {
    if (identical(marglik[["aggregation"]][["nonfinite_policy"]], "drop")) {
      values <- values[is.finite(values)]
    } else {
      stop("Stored log marginal likelihoods must contain only finite values.",
           call. = FALSE)
    }
  }
  if (length(values) == 0L) {
    stop("No included log marginal likelihood repetitions are available.",
         call. = FALSE)
  }
  return(list(logml = values, repeated = repeated))
}


.brma_recycle_logml <- function(logml) {

  repetition_counts <- lengths(logml)
  n_repetitions <- max(repetition_counts)
  if (any(repetition_counts != n_repetitions)) {
    warning("Not all objects provide ", n_repetitions,
            " logmls. Some values are recycled.", call. = FALSE)
    logml <- lapply(logml, rep, length.out = n_repetitions)
  }
  return(logml)
}


#' @title Posterior Model Probabilities for brma Objects
#'
#' @description Compute posterior model probabilities from marginal
#' likelihoods of brma models.
#'
#' @param x a brma model object.
#' @param ... additional brma model objects.
#' @param prior_prob numeric vector with prior model probabilities or weights.
#' Values must be finite, nonnegative, have the same length as the retained
#' models, and have a positive total. If omitted, a uniform prior is used.
#' Supplied values are normalized internally.
#' @param model_names character vector with model names. If \code{NULL}
#' (the default), names will be derived from deparsing the call.
#'
#' @details
#' The marginal likelihoods must first be computed using \code{\link{add_marglik}}.
#' \code{x} and at least one additional \code{brma} model must be supplied.
#' Non-\code{brma} objects in \code{...} are ignored with a warning. All retained
#' models must be fitted to the same outcome target/data, including outcome type
#' and, when present, likelihood weights. Cluster identifiers define model
#' structure and are not part of the outcome-data identity.
#'
#' If any model contains repeated bridge estimates, comparisons retain those
#' repetitions regardless of model order. Models with fewer repetitions are
#' recycled to the largest repetition count with a warning, following
#' \pkg{bridgesampling}. Prior weights are normalized before each comparison.
#' A stored explicit non-finite-repetition drop policy is honored; otherwise
#' standardized marginal-likelihood results must be finite. Legacy raw
#' \code{"bridge_list"} results retain non-finite repetitions and use upstream
#' missing-value propagation and warnings.
#'
#' @return A named numeric vector with posterior model probabilities for single
#' estimates, or a repetition-by-model numeric matrix when any model has repeated
#' estimates. Each successfully evaluated comparison sums to one.
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
  comparisons <- lapply(models, .brma_comparison_logml)
  logml_values <- lapply(comparisons, "[[", "logml")
  if (any(vapply(comparisons, "[[", logical(1), "repeated"))) {
    prior_prob <- .post_prob_normalize_prior(prior_prob, length(models))
    .post_prob_check_names(model_names, length(models))
    logml_values <- .brma_recycle_logml(logml_values)
    posterior <- vapply(seq_along(logml_values[[1L]]), function(repetition) {

      bridgesampling::post_prob(
        vapply(logml_values, "[[", numeric(1), repetition),
        prior_prob = prior_prob, model_names = model_names
      )
    }, numeric(length(models)))
    return(t(posterior))
  }

  return(.post_prob_from_logml(
    logml       = unlist(logml_values, use.names = FALSE),
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
#' outcome target/data, including outcome type and, when present, likelihood
#' weights. Cluster identifiers define model structure and are not part of the
#' outcome-data identity.
#'
#' When either model contains repeated bridge estimates, the shorter sequence
#' is recycled with a warning, following \pkg{bridgesampling}. Repetitions are
#' compared in order. The median-based estimate is computed from the separate
#' medians of the log marginal likelihoods before recycling, rather than the
#' median of the repetition-wise Bayes factors. A recorded explicit drop policy
#' excludes non-finite repetitions from standardized results. Legacy raw
#' \code{"bridge_list"} results retain failed repetitions in the returned
#' estimates, with upstream missing-value behavior when printed.
#'
#' @return For single estimates, a list of class \code{"bf_default"} with components:
#' \itemize{
#'   \item \code{bf}: (scalar) value of the Bayes factor in favor of
#'   \code{x1} over \code{x2}.
#'   \item \code{log}: Boolean indicating whether \code{bf} corresponds
#'   to the log Bayes factor.
#' }
#' For repeated estimates, a list of class \code{c("bf.brma", "bf_bridge_list")}
#' with a vector \code{bf}, scalar \code{bf_median_based}, and \code{log}. It
#' uses the upstream repeated-Bayes-factor print method. Both
#' \code{as.data.frame()} and \code{data.frame()} return a long data frame with
#' \code{component}, \code{parameter}, and full-precision numeric \code{value}
#' columns, containing every repetition and the displayed median-based, range,
#' and interquartile-range summaries.
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
  comparison1 <- .brma_comparison_logml(x1)
  comparison2 <- .brma_comparison_logml(x2)
  logml1 <- comparison1[["logml"]]
  logml2 <- comparison2[["logml"]]
  if (comparison1[["repeated"]] || comparison2[["repeated"]]) {
    BayesTools::check_bool(log, "log")
    logml_values <- .brma_recycle_logml(list(logml1, logml2))
    values <- logml_values[[1L]] - logml_values[[2L]]
    median_based <- stats::median(logml1, na.rm = TRUE) -
      stats::median(logml2, na.rm = TRUE)
    out <- list(
      bf = if (log) values else exp(values),
      bf_median_based = if (log) median_based else exp(median_based),
      log = log
    )
    class(out) <- c("bf.brma", "bf_bridge_list")
    mc <- match.call()
    attr(out, "model_names") <- c(
      paste(deparse(mc[["x1"]]), collapse = ""),
      paste(deparse(mc[["x2"]]), collapse = "")
    )
    return(out)
  }

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
  out <- bf.brma(x1, x2, log = log, ...)
  if (inherits(out, "bf.brma")) {
    mc <- match.call()
    attr(out, "model_names") <- c(
      paste(deparse(mc[["x1"]]), collapse = ""),
      paste(deparse(mc[["x2"]]), collapse = "")
    )
  }
  return(out)
}


#' @export
as.data.frame.bf.brma <- function(x, row.names = NULL, optional = FALSE, ...) {

  values <- x[["bf"]]
  limits <- range(values, na.rm = TRUE)
  scale <- if (x[["log"]]) "log_BF" else "BF"
  out <- data.frame(
    component = c(rep("repetitions", length(values)), rep("summary", 4L)),
    parameter = c(
      paste0(scale, "[", seq_along(values), "]"),
      paste0(scale, c("_median_based", "_minimum", "_maximum", "_IQR"))
    ),
    value = c(values, x[["bf_median_based"]], limits, stats::IQR(values, na.rm = TRUE)),
    row.names = row.names,
    stringsAsFactors = FALSE
  )
  return(out)
}
