# ============================================================================ #
# brma.loo.R
# ============================================================================ #
#
# This file implements LOO-PSIS (Leave-One-Out Cross-Validation via Pareto
# Smoothed Importance Sampling) diagnostics for brma class objects.
#
# LOO can be computed at the estimate unit or, for multilevel models, at the
# cluster unit.
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
#' @param unit output/deletion unit. \code{"estimate"} computes one contribution
#' per effect-size estimate. \code{"cluster"} computes one contribution per
#' cluster and is available only for multilevel models.
#' @param r_eff optional vector of relative effective sample sizes. If not
#' provided, it is computed from the log-likelihood values.
#' @param parallel Logical. If \code{TRUE}, \code{loo::relative_eff()} and
#' \code{loo::loo()} use \code{RoBMA.get_option("max_cores")}. Log-likelihood
#' construction is unchanged. If \code{FALSE}, those computations use one core.
#' @param ... additional arguments (currently ignored).
#'
#' @details
#' With \code{unit = "estimate"}, LOO-CV is computed with one contribution per
#' effect-size estimate. For binomial and Poisson models, each pair of counts
#' (ai/ci or x1i/x2i) that defines a single effect size estimate is treated as
#' one contribution.
#'
#' With \code{unit = "cluster"}, LOO-CV is computed with one joint contribution
#' per cluster. For unweighted normal models without selection this uses the
#' analytic cluster block covariance. Selection, data-weighted normal, and
#' GLMM models integrate the held-out cluster effect with Gauss-Hermite
#' quadrature.
#'
#' For selection models, the LOO evaluates the weighted likelihood, conditioning
#' on the posterior omega samples.
#'
#' For correlated known-\code{V} \code{brma.mv()} models, estimate-unit LOO uses
#' conditional predictive contributions \eqn{p(y_i \mid y_{-i}, \theta)} within
#' known-\code{V} dependency blocks. This is an existing-estimate diagnostic, not
#' an independent new-estimate prediction target.
#'
#' For \code{brma.mv()} models with known random-effect group covariance
#' \code{R}, estimate-unit LOO keeps sampled random effects at the estimate
#' conditioning depth. The known \code{R} matrix shapes the posterior and prior
#' for those sampled random effects, but it is not added again as a marginal
#' \eqn{ZGZ'} covariance term in this target.
#'
#' The PSIS object is essential for model comparison via
#' \code{\link[loo]{loo_compare}} and is automatically saved in the loo result.
#' RoBMA stores target metadata so comparisons can reject mismatched data,
#' unit, or conditioning-depth targets.
#'
#' \strong{Important for model comparison:} When comparing models via
#' \code{\link[loo]{loo_compare}}, the selection is based on expected
#' out-of-sample predictive performance. This evaluates how well models predict
#' \emph{new} observations, not how well they fit the observed data.
#'
#' @return The brma object with the LOO result stored in
#' \code{object[["loo"]][[unit]]}.
#'
#' @seealso \code{\link{loo.brma}}, \code{\link[loo]{loo}},
#' \code{\link[loo]{loo_compare}}, \code{\link[loo]{pareto_k_ids}}
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   fit <- add_loo(fit)
#'   loo_fit <- loo(fit)
#'   print(loo_fit)
#'   plot(loo_fit)
#' }
#' }
#'
#' @references
#' \insertCite{vehtari2017practical}{RoBMA}
#' \insertCite{vehtari2024pareto}{RoBMA}
#'
#' @aliases add_loo
#' @export
add_loo.brma <- function(object, unit = "estimate", r_eff = NULL, parallel = FALSE, ...) {
  unit <- .normalize_unit(unit)
  BayesTools::check_bool(parallel, "parallel")

  conditioning_depth <- .loo_conditioning_depth_from_unit(unit)
  .check_unit_conditioning_depth(
    object             = object,
    unit               = unit,
    conditioning_depth = conditioning_depth,
    caller             = "add_loo()"
  )

  # compute the log-likelihood matrix (S x K)
  log_lik <- .log_lik.brma(object, unit = unit, caller = "add_loo()")
  target  <- attr(log_lik, "RoBMA_target", exact = TRUE)
  .warn_known_v_schur_log_score(target, "LOO")

  # determine number of cores based on `parallel` and package options
  cores <- if (parallel) max(1, RoBMA.get_option("max_cores")) else 1

  # compute relative effective sample sizes if not provided
  if (is.null(r_eff)) {
    # loo::relative_eff expects exp(log_lik) with chain_id for matrix input
    chain_id <- .loo_chain_id(object[["fit"]], n_samples = nrow(log_lik))

    r_eff <- loo::relative_eff(exp(log_lik), chain_id = chain_id, cores = cores)
  }

  # call loo on the log-likelihood matrix
  loo_result <- loo::loo(log_lik, r_eff = r_eff, save_psis = TRUE, cores = cores)
  loo_result <- .add_loo_target_metadata(
    object             = loo_result,
    unit               = target[["unit"]],
    conditioning_depth = .get_target_conditioning_depth(target),
    targets            = target[["targets"]],
    data_hash          = target[["data_hash"]],
    metadata           = target
  )

  # store in object and return
  if (is.null(object[["loo"]])) {
    object[["loo"]] <- list()
  }
  object[["loo"]][[unit]] <- loo_result
  return(object)
}


# ---------------------------------------------------------------------------- #
# .loo_chain_id
# ---------------------------------------------------------------------------- #
#
# Derive chain IDs from the retained MCMC samples.
#
# @param fit       runjags fit object.
# @param n_samples integer; expected number of posterior rows.
#
# @return integer vector of length n_samples.
#
# ---------------------------------------------------------------------------- #
.loo_chain_id <- function(fit, n_samples) {

  if (is.null(fit[["mcmc"]])) {
    stop(
      "Cannot infer LOO chain IDs because fitted MCMC draws are missing. ",
      "Supply 'r_eff' explicitly.",
      call. = FALSE
    )
  }

  chain_lengths <- vapply(fit[["mcmc"]], NROW, integer(1))
  chain_id      <- rep(seq_along(chain_lengths), times = chain_lengths)

  if (length(chain_id) != n_samples) {
    stop(
      "Could not derive valid chain IDs for relative_eff(): expected ",
      n_samples, " posterior rows but got ", length(chain_id), ". ",
      "Supply 'r_eff' explicitly.",
      call. = FALSE
    )
  }

  return(chain_id)
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
#' @param unit output/deletion unit. See \code{\link{add_loo}}; the same
#' accepted values and multilevel constraint apply.
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
#' For correlated known-\code{V} \code{brma.mv()} models, estimate-unit WAIC uses
#' the same conditional log-likelihood matrix as estimate-unit LOO,
#' \eqn{p(y_i \mid y_{-i}, \theta)}. It therefore has the same interpretation as
#' an existing-estimate diagnostic rather than an independent new-estimate
#' prediction target.
#'
#' For \code{brma.mv()} models with known random-effect group covariance
#' \code{R}, estimate-unit WAIC has the same conditioning convention as
#' estimate-unit LOO: sampled random effects are conditioned on, and known
#' \code{R} is not added again as a marginal \eqn{ZGZ'} covariance term.
#'
#' @return The brma object with the WAIC result stored in
#' \code{object[["waic"]][[unit]]}.
#'
#' @seealso \code{\link{waic.brma}}, \code{\link{add_loo}}, \code{\link[loo]{waic}}
#'
#' @aliases add_waic
#' @export
add_waic.brma <- function(object, unit = "estimate", ...) {
  unit <- .normalize_unit(unit)
  conditioning_depth <- .loo_conditioning_depth_from_unit(unit)
  .check_unit_conditioning_depth(
    object             = object,
    unit               = unit,
    conditioning_depth = conditioning_depth,
    caller             = "add_waic()"
  )

  # compute the log-likelihood matrix (S x K)
  log_lik <- .log_lik.brma(object, unit = unit, caller = "add_waic()")
  target  <- attr(log_lik, "RoBMA_target", exact = TRUE)
  .warn_known_v_schur_log_score(target, "WAIC")

  # call waic on the log-likelihood matrix
  waic_result <- loo::waic(log_lik, ...)
  waic_result <- .add_loo_target_metadata(
    object             = waic_result,
    unit               = target[["unit"]],
    conditioning_depth = .get_target_conditioning_depth(target),
    targets            = target[["targets"]],
    data_hash          = target[["data_hash"]],
    metadata           = target
  )

  # store in object and return
  if (is.null(object[["waic"]])) {
    object[["waic"]] <- list()
  }
  object[["waic"]][[unit]] <- waic_result
  return(object)
}


.warn_known_v_schur_log_score <- function(target, method) {

  if (!isTRUE(target[["known_v_schur"]])) {
    return(invisible(FALSE))
  }

  target_row <- .brma_mv_target_row(if (identical(method, "WAIC")) {
    "add_waic()/waic()"
  } else {
    "add_loo()/loo()"
  })

  warning(
    "Estimate-unit ", method, " for brma.mv() known-V models uses conditional ",
    "predictive contributions p(y_i | y_-i, theta). Interpret it as an ",
    "existing-estimate diagnostic, not independent new-estimate prediction. ",
    "Target: ", target_row[["target"]], ".",
    call. = FALSE
  )

  return(invisible(TRUE))
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
#' @param unit output/deletion unit. See \code{\link{add_loo}}.
#' @param ... additional arguments (currently unused).
#'
#' @details
#' This function extracts the LOO object that was previously computed and
#' stored using \code{object <- add_loo(object, unit = unit)}. If LOO has not
#' been computed for the requested unit, an error is thrown.
#'
#' This is the RoBMA S3 generic and \code{brma} method. Use
#' \code{\link[loo]{loo}} directly for raw log-likelihood arrays or matrices.
#'
#' @return An object of class \code{c("psis_loo", "loo")} as returned by
#' \code{\link[loo]{loo}}.
#'
#' @seealso \code{\link{add_loo}}, \code{\link[loo]{loo}},
#' \code{\link[loo]{loo_compare}}, \code{\link[loo]{pareto_k_ids}}
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit <- add_loo(fit)
#'
#'   loo_fit <- loo(fit)
#'   print(loo_fit)
#' }
#' }
#'
#' @aliases loo
#' @export
loo.brma <- function(x, unit = "estimate", ...) {
  return(.check_loo_target(x, unit = unit))
}


#' @title Extract Log-Likelihood Matrix from brma Object
#'
#' @description Extract the pointwise log-likelihood matrix from a brma model
#' object. This is an S x K or S x G matrix where S is the number of posterior
#' samples, K is the number of estimates, and G is the number of clusters.
#' This method implements the S3
#' \'logLik\' generic for \code{brma} objects and returns the matrix of
#' pointwise log-likelihoods (one column per observation, one row per sample).
#'
#' @param object a brma model object.
#' @param unit output unit. See \code{\link{add_loo}}.
#' @param ... currently unused.
#'
#' @details
#' The log-likelihood is computed for each observation at each posterior sample.
#' For binomial and Poisson models, each observation consists of a pair of
#' counts (ai/ci or x1i/x2i) that together define a single effect size estimate.
#'
#' @return An S x K or S x G matrix of log-likelihood values.
#'
#' @seealso \code{\link{loo.brma}}
#'
#' @export
logLik.brma <- function(object, unit = "estimate", ...) {
  out <- .log_lik.brma(object, unit = unit, caller = "logLik()")
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
#' @param unit output/deletion unit used when extracting LOO from brma objects.
#'
#' @details
#' This function compares models based on their expected out-of-sample
#' predictive performance (ELPD).
#'
#' \strong{Important for model comparison:} When comparing models via
#' \code{\link[loo]{loo_compare}}, the selection is based on expected
#' out-of-sample predictive performance. This evaluates how well models predict
#' \emph{new} observations, not how well they fit the observed data.
#' RoBMA rejects comparisons with different outcome targets/data, \code{unit},
#' or implied \code{conditioning_depth}.
#'
#' @return A matrix of class \code{"compare.loo"} as returned by
#' \code{\link[loo]{loo_compare}}.
#'
#' @seealso \code{\link{loo.brma}}, \code{\link[loo]{loo_compare}}
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit_bias <- RoBMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit_nobias <- BMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   fit_bias <- add_loo(fit_bias)
#'   fit_nobias <- add_loo(fit_nobias)
#'
#'   loo_compare(fit_bias, fit_nobias)
#'   loo_compare(loo(fit_bias), loo(fit_nobias))
#' }
#' }
#'
#' @aliases loo_compare
#' @export
loo_compare.brma <- function(x, ..., unit = "estimate") {
  unit <- .normalize_unit(unit)

  # collect all models: x plus any in ...
  models <- c(list(x), list(...))

  if (length(models) < 2) {
    stop("At least two models are required for comparison.", call. = FALSE)
  }

  # convert brma objects to loo objects if necessary
  loo_objects <- lapply(models, function(m) {
    if (inherits(m, "brma")) {
      loo.brma(m, unit = unit)
    } else if (inherits(m, "loo")) {
      m
    } else {
      stop("All arguments must be brma or loo objects.", call. = FALSE)
    }
  })

  # call the 'loo' package's default implementation to avoid dispatch recursion
  .check_loo_compare_targets(loo_objects)
  loo_compare_fun <- get("loo_compare.default", envir = asNamespace("loo"), inherits = FALSE)
  result <- do.call(loo_compare_fun, loo_objects)

  return(result)
}


#' @title Compare loo Objects Using LOO
#'
#' @description Method for comparing RoBMA-targeted \code{loo} or \code{waic}
#' objects directly.
#'
#' @param x a RoBMA-targeted \code{loo} or \code{waic} object (the first model
#' to compare).
#' @param ... additional RoBMA-targeted \code{loo}/\code{waic} or \code{brma}
#' objects to compare.
#' @param unit output/deletion unit used when extracting LOO from brma objects.
#'
#' @return A matrix of class \code{"compare.loo"} as returned by
#' \code{\link[loo]{loo_compare}}.
#'
#' @seealso \code{\link{loo.brma}}, \code{\link[loo]{loo_compare}}
#'
#' @export
loo_compare.loo <- function(x, ..., unit = "estimate") {
  unit <- .normalize_unit(unit)

  # collect all models: x plus any in ...
  models <- c(list(x), list(...))

  if (length(models) < 2) {
    stop("At least two models are required for comparison.", call. = FALSE)
  }

  # convert brma objects to loo objects if necessary
  loo_objects <- lapply(models, function(m) {
    if (inherits(m, "brma")) {
      loo.brma(m, unit = unit)
    } else if (inherits(m, "loo")) {
      m
    } else {
      stop("All arguments must be brma or loo objects.", call. = FALSE)
    }
  })

  # call the 'loo' package's default implementation to avoid dispatch recursion
  .check_loo_compare_targets(loo_objects)
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
#' @param unit output/deletion unit. See \code{\link{add_loo}}.
#' @param ... additional arguments (currently unused).
#'
#' @details
#' This function extracts the WAIC object that was previously computed and
#' stored using \code{object <- add_waic(object, unit = unit)}. If WAIC has not
#' been computed for the requested unit, an error is thrown.
#'
#' This is the RoBMA S3 generic and \code{brma} method. The method is also
#' registered for \code{\link[loo]{waic}}, so \code{loo::waic(fit)} extracts
#' the cached WAIC object for \code{brma} fits. Use \code{\link[loo]{waic}}
#' directly for raw log-likelihood arrays or matrices.
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
#' @aliases waic
#' @export
#' @exportS3Method loo::waic
waic.brma <- function(x, unit = "estimate", ...) {
  return(.check_waic_target(x, unit = unit))
}


# ---------------------------------------------------------------------------- #
# loo_weights generic and methods
# ---------------------------------------------------------------------------- #

#' @export
loo_weights <- function(object, ...) UseMethod("loo_weights")


#' @title Extract Normalized PSIS Weights from brma Object
#'
#' @description Extract the normalized Pareto smoothed importance sampling
#' (PSIS) weights from a brma model object.
#'
#' @param object a brma model object.
#' @param unit output/deletion unit. See \code{\link{add_loo}}.
#' @param ... currently unused.
#'
#' @details LOO must first be computed with
#' \code{object <- add_loo(object, unit = unit)}. This method extracts the
#' stored PSIS object and does not compute LOO.
#'
#' @return An \code{S x K} matrix for estimate-unit LOO, or \code{S x G}
#' matrix for cluster-unit LOO, with posterior samples in rows and LOO targets
#' in columns. Columns are normalized to sum to one.
#'
#' @aliases loo_weights
#' @exportS3Method
loo_weights.brma <- function(object, unit = "estimate", ...) {
  # extract loo
  loo_result <- loo.brma(object, unit = unit)

  # extract PSIS object and get normalized weights
  psis_object <- loo_result[["psis_object"]]
  psis_weights <- loo::weights.importance_sampling(psis_object, log = FALSE, normalize = TRUE)

  return(psis_weights)
}


# ---------------------------------------------------------------------------- #
# check_loo generic and methods
# ---------------------------------------------------------------------------- #

#' @export
check_loo <- function(object, ...) UseMethod("check_loo")


#' @title Check LOO Diagnostics for brma Objects
#'
#' @description Check Pareto k diagnostics for a brma model object and
#' warn if any values are high.
#'
#' @param object a brma model object.
#' @param unit output/deletion unit. See \code{\link{add_loo}}.
#' @param ... currently unused.
#'
#' @details LOO must first be computed with
#' \code{object <- add_loo(object, unit = unit)}. The method warns when any
#' Pareto \eqn{k} diagnostic is greater than 0.7.
#'
#' @return NULL (throws warning if diagnostics are unreliable).
#'
#' @aliases check_loo
#' @exportS3Method
check_loo.brma <- function(object, unit = "estimate", ...) {
  # extract loo
  unit       <- .normalize_unit(unit)
  loo_result <- loo.brma(object, unit = unit)

  # check Pareto k diagnostics and warn if unreliable
  pareto_k <- loo_result[["diagnostics"]][["pareto_k"]]
  bad_k <- which(pareto_k > 0.7)
  if (length(bad_k) > 0) {
    warning(
      "Some Pareto k values are high (> 0.7), indicating potentially unreliable ",
      "LOO diagnostics for ", unit, "s: ", paste(bad_k, collapse = ", "), ". ",
      "Inspect the loo fit by using loo(object).",
      call. = FALSE
    )
  }

  return()
}
