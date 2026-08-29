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
# Re-export loo generics
# ---------------------------------------------------------------------------- #

#' @importFrom loo loo_model_weights
#' @export
loo::loo_model_weights


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
#' Estimate-unit GLMM log-likelihoods support point or beta baserate priors and
#' point or normal log-rate priors, including their supported truncations.
#' Other nuisance-prior families are rejected when the model input is validated.
#'
#' With \code{unit = "cluster"}, LOO-CV is computed with one joint contribution
#' per cluster. For unweighted normal models without selection this uses the
#' analytic cluster block covariance. Selection and data-weighted normal models
#' integrate the held-out cluster effect with Gauss-Hermite quadrature.
#' Cluster-unit binomial and Poisson GLMM log-likelihoods are unavailable until
#' certified nested adaptive quadrature is implemented; use
#' \code{unit = "estimate"} for GLMMs.
#'
#' For approximate selection models, LOO evaluates the row-wise selected-normal
#' likelihood conditional on sampled shared effects and posterior omega. For
#' exact selection models, estimate deletion uses the selected Gaussian Schur
#' conditional \eqn{p_E(y_i \mid y_{-i}, \theta)}. The selection weights of the
#' retained estimates cancel; the deleted estimate retains its own weight and
#' conditional selection normalizer. This is the estimate-level finite-vector
#' deletion target, not prediction for an entirely new dependency block.
#'
#' For Gaussian multilevel and known-\code{V} \code{brma.mv()} models,
#' estimate-unit LOO integrates Gaussian local effects and uses
#' \eqn{p(y_i \mid y_{-i}, \theta)} within each fitted dependency block. Thus
#' deleting one estimate retains the other estimates in its cluster or random-
#' effect block; it does not represent prediction for a wholly new group.
#'
#' For \code{brma.mv()} models with known random-effect group covariance
#' \code{R}, the known \code{R} matrix and fitted random-effect covariance enter
#' the marginal \eqn{ZGZ'} covariance through BayesTools' compiled random-effect
#' metadata. The same calculation applies to every supported random structure.
#'
#' Estimate-unit deletion retains the fitted grouping and known-\code{V}
#' dependency structure; it does not represent deletion of an entire new random
#' group. Models with no local effects, sampled local effects, marginalized
#' local effects, or a mixture of sampled and marginalized effects remain
#' comparable when their data hash, deletion unit, retained context, and
#' likelihood target agree. The sampled/marginalized labels are retained as
#' provenance rather than comparison keys.
#'
#' The PSIS object is essential for model comparison via
#' \code{\link[loo]{loo_compare}} and is automatically saved in the loo result.
#' RoBMA stores target metadata so comparisons can reject mismatched data,
#' unit, or retained-context targets.
#' Each finite log-likelihood column observed to be constant across posterior
#' draws uses exact uniform importance ratios. Such columns bypass relative-ESS
#' and generalized-Pareto fitting, including when other columns vary. Their
#' Pareto-k diagnostic is recorded as the sentinel zero, meaning constant
#' ratios with no tail instability rather than an estimated Pareto shape.
#' Finite varying columns are centered by a column-specific constant before
#' relative-ESS and PSIS evaluation. The offset is restored exactly afterward,
#' preserving the LOO target while avoiding probability-scale underflow.
#'
#' \strong{Important for model comparison:} When comparing models via
#' \code{\link[loo]{loo_compare}}, the selection is based on expected
#' predictive performance for the documented deletion target, not in-sample
#' fit. For dependent estimates, deleting one estimate may retain other
#' observations from its fitted dependency block; this is not new-group
#' prediction.
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

  # compute the log-likelihood matrix (S x K)
  log_lik <- .log_lik.brma(object, unit = unit, caller = "add_loo()")
  target  <- attr(log_lik, "RoBMA_target", exact = TRUE)

  # determine number of cores based on `parallel` and package options
  cores <- if (parallel) max(1, RoBMA.get_option("max_cores")) else 1

  deterministic     <- .deterministic_log_lik_columns(log_lik)
  variable          <- !deterministic
  centered_variable <- NULL
  if (any(variable)) {
    variable_log_lik <- if (all(variable)) {
      log_lik
    } else {
      log_lik[, variable, drop = FALSE]
    }
    centered_variable <- .loo_center_finite_columns(variable_log_lik)
  }

  # compute relative effective sample sizes if not provided
  if (is.null(r_eff)) {
    r_eff <- rep(1, ncol(log_lik))
    if (any(variable)) {
      # loo::relative_eff expects exp(log_lik) with chain_id for matrix input
      chain_id <- .loo_chain_id(object[["fit"]], n_samples = nrow(log_lik))
      r_eff[variable] <- loo::relative_eff(
        exp(centered_variable[["log_lik"]]),
        chain_id = chain_id,
        cores    = cores
      )
    }
  }
  r_eff <- .loo_prepare_r_eff(
    r_eff         = r_eff,
    n_units       = ncol(log_lik),
    deterministic = deterministic
  )

  loo_result <- .loo_with_deterministic_columns(
    log_lik           = log_lik,
    r_eff             = r_eff,
    deterministic     = deterministic,
    cores             = cores,
    centered_variable = centered_variable
  )
  loo_result <- .add_loo_target_metadata(
    object           = loo_result,
    unit             = target[["unit"]],
    retained_context = target[["retained_context"]],
    targets          = target[["targets"]],
    data_hash        = target[["data_hash"]],
    metadata         = target
  )

  # store in object and return
  if (is.null(object[["loo"]])) {
    object[["loo"]] <- list()
  }
  object[["loo"]][[unit]] <- loo_result
  return(object)
}


.deterministic_log_lik_columns <- function(log_lik) {

  if (!is.matrix(log_lik) || nrow(log_lik) == 0L) {
    return(rep(FALSE, if (is.matrix(log_lik)) ncol(log_lik) else 0L))
  }

  return(vapply(seq_len(ncol(log_lik)), function(i) {
    column <- log_lik[, i]
    all(is.finite(column)) && isTRUE(all(column == column[[1L]]))
  }, logical(1)))
}


.loo_prepare_r_eff <- function(r_eff, n_units, deterministic) {

  if (length(r_eff) == 1L) {
    r_eff <- rep(r_eff, n_units)
  } else if (length(r_eff) != n_units) {
    stop("'r_eff' must have one value or one value per observation.",
         call. = FALSE)
  }
  r_eff[deterministic] <- 1
  invalid <- !deterministic & (!is.finite(r_eff) | r_eff <= 0)
  if (any(invalid)) {
    stop(
      "'r_eff' must be finite and positive for every varying observation.",
      call. = FALSE
    )
  }

  return(r_eff)
}


.loo_center_finite_columns <- function(log_lik) {

  offsets <- rep(0, ncol(log_lik))
  names(offsets) <- colnames(log_lik)
  finite <- vapply(seq_len(ncol(log_lik)), function(i) {
    all(is.finite(log_lik[, i]))
  }, logical(1))
  if (any(finite)) {
    offsets[finite] <- apply(log_lik[, finite, drop = FALSE], 2L, max)
    log_lik[, finite] <- sweep(
      log_lik[, finite, drop = FALSE],
      MARGIN = 2L,
      STATS  = offsets[finite],
      FUN    = "-"
    )
  }

  return(list(log_lik = log_lik, offsets = offsets))
}


.loo_with_deterministic_columns <- function(log_lik, r_eff, deterministic,
                                            cores, centered_variable = NULL) {

  if (!any(deterministic)) {
    centered <- centered_variable
    if (is.null(centered)) {
      centered <- .loo_center_finite_columns(log_lik)
    }
    result <- loo::loo(
      centered[["log_lik"]],
      r_eff     = r_eff,
      save_psis = TRUE,
      cores     = cores
    )
    return(.loo_restore_log_lik_offsets(result, centered[["offsets"]]))
  }

  variable <- which(!deterministic)
  variable_result <- NULL
  if (length(variable) > 0L) {
    centered <- centered_variable
    if (is.null(centered)) {
      centered <- .loo_center_finite_columns(
        log_lik[, variable, drop = FALSE]
      )
    }
    variable_result <- loo::loo(
      centered[["log_lik"]],
      r_eff     = r_eff[variable],
      save_psis = TRUE,
      cores     = cores
    )
    variable_result <- .loo_restore_log_lik_offsets(
      variable_result,
      centered[["offsets"]]
    )
  }

  .loo_combine_deterministic_columns(
    variable_result = variable_result,
    log_lik         = log_lik,
    r_eff           = r_eff,
    deterministic   = deterministic
  )
}


.loo_restore_log_lik_offsets <- function(result, offsets) {

  if (length(offsets) == 0L || all(offsets == 0)) {
    return(result)
  }

  pointwise <- result[["pointwise"]]
  if (length(offsets) != nrow(pointwise)) {
    stop("LOO offset dimensions are inconsistent.", call. = FALSE)
  }
  pointwise[, "elpd_loo"] <- pointwise[, "elpd_loo"] + offsets
  pointwise[, "looic"]    <- pointwise[, "looic"] - 2 * offsets
  result[["pointwise"]]   <- pointwise

  psis <- result[["psis_object"]]
  psis[["log_weights"]] <- sweep(
    psis[["log_weights"]],
    MARGIN = 2L,
    STATS  = offsets,
    FUN    = "-"
  )
  attr(psis, "norm_const_log") <-
    attr(psis, "norm_const_log", exact = TRUE) - offsets
  result[["psis_object"]] <- psis

  estimates <- .loo_estimates_from_pointwise(pointwise)
  result[["estimates"]]   <- estimates
  result[["elpd_loo"]]    <- estimates["elpd_loo", "Estimate"]
  result[["p_loo"]]       <- estimates["p_loo", "Estimate"]
  result[["looic"]]       <- estimates["looic", "Estimate"]
  result[["se_elpd_loo"]] <- estimates["elpd_loo", "SE"]
  result[["se_p_loo"]]    <- estimates["p_loo", "SE"]
  result[["se_looic"]]    <- estimates["looic", "SE"]

  return(result)
}


.loo_combine_deterministic_columns <- function(variable_result, log_lik,
                                               r_eff, deterministic) {

  n_samples  <- nrow(log_lik)
  n_units    <- ncol(log_lik)
  exact      <- which(deterministic)
  variable   <- which(!deterministic)
  unit_names <- colnames(log_lik)
  pointwise_names <- c(
    "elpd_loo", "mcse_elpd_loo", "p_loo", "looic",
    "influence_pareto_k"
  )

  pointwise <- matrix(
    NA_real_,
    nrow     = n_units,
    ncol     = length(pointwise_names),
    dimnames = list(unit_names, pointwise_names)
  )
  if (length(variable) > 0L) {
    pointwise[variable, ] <- variable_result[["pointwise"]]
  }
  if (length(exact) > 0L) {
    constant <- log_lik[1L, exact]
    pointwise[exact, ] <- cbind(
      elpd_loo           = constant,
      mcse_elpd_loo      = 0,
      p_loo              = 0,
      looic              = -2 * constant,
      influence_pareto_k = 0
    )
  }

  diagnostics <- list(
    pareto_k = rep(0, n_units),
    n_eff    = rep(n_samples, n_units),
    r_eff    = r_eff
  )
  if (length(variable) > 0L) {
    variable_diagnostics <- variable_result[["diagnostics"]]
    diagnostics[["pareto_k"]][variable] <-
      variable_diagnostics[["pareto_k"]]
    diagnostics[["n_eff"]][variable] <- variable_diagnostics[["n_eff"]]
    diagnostics[["r_eff"]][variable] <- variable_diagnostics[["r_eff"]]
  }

  psis <- .loo_combine_psis(
    variable_result = variable_result,
    log_lik         = log_lik,
    diagnostics     = diagnostics,
    deterministic   = deterministic
  )

  estimates <- .loo_estimates_from_pointwise(pointwise)
  result <- list(
    estimates   = estimates,
    pointwise   = pointwise,
    diagnostics = diagnostics,
    psis_object = psis,
    elpd_loo    = estimates["elpd_loo", "Estimate"],
    p_loo       = estimates["p_loo", "Estimate"],
    looic       = estimates["looic", "Estimate"],
    se_elpd_loo = estimates["elpd_loo", "SE"],
    se_p_loo    = estimates["p_loo", "SE"],
    se_looic    = estimates["looic", "SE"]
  )
  attr(result, "dims") <- c(n_samples, n_units)
  class(result) <- c("psis_loo", "importance_sampling_loo", "loo")

  return(result)
}


.loo_combine_psis <- function(variable_result, log_lik, diagnostics,
                              deterministic) {

  n_samples  <- nrow(log_lik)
  n_units    <- ncol(log_lik)
  exact      <- which(deterministic)
  variable   <- which(!deterministic)
  unit_names <- colnames(log_lik)

  log_weights <- matrix(
    NA_real_,
    nrow     = n_samples,
    ncol     = n_units,
    dimnames = list(NULL, unit_names)
  )
  norm_const_log <- rep(NA_real_, n_units)
  tail_len       <- integer(n_units)

  if (length(variable) > 0L) {
    variable_psis <- variable_result[["psis_object"]]
    log_weights[, variable] <- variable_psis[["log_weights"]]
    norm_const_log[variable] <- attr(
      variable_psis,
      "norm_const_log",
      exact = TRUE
    )
    tail_len[variable] <- attr(variable_psis, "tail_len", exact = TRUE)
  }
  if (length(exact) > 0L) {
    log_weights[, exact]  <- 0
    norm_const_log[exact] <- log(n_samples)
  }
  if (!is.null(unit_names)) {
    names(norm_const_log) <- unit_names
    names(tail_len)       <- unit_names
  }

  psis <- list(
    log_weights = log_weights,
    diagnostics = diagnostics
  )
  attr(psis, "norm_const_log") <- norm_const_log
  attr(psis, "tail_len")       <- tail_len
  attr(psis, "r_eff")          <- diagnostics[["r_eff"]]
  attr(psis, "dims")           <- c(n_samples, n_units)
  attr(psis, "method")         <- "psis"
  class(psis) <- c("psis", "importance_sampling", "list")

  return(psis)
}


.loo_estimates_from_pointwise <- function(pointwise) {

  n_units <- nrow(pointwise)
  estimate_names <- c("elpd_loo", "p_loo", "looic")
  estimate <- colSums(pointwise[, estimate_names, drop = FALSE])
  standard_error <- vapply(estimate_names, function(name) {
    sqrt(n_units * stats::var(pointwise[, name]))
  }, numeric(1))

  estimates <- cbind(Estimate = estimate, SE = standard_error)
  rownames(estimates) <- estimate_names

  return(estimates)
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
#' For Gaussian multilevel and known-\code{V} \code{brma.mv()} models,
#' estimate-unit WAIC uses the same integrated conditional log-likelihood matrix
#' as estimate-unit LOO, \eqn{p(y_i \mid y_{-i}, \theta)}. Estimate deletion
#' retains the fitted grouping and dependency structure.
#'
#' For \code{brma.mv()} models with known random-effect group covariance
#' \code{R}, the known \code{R} matrix and fitted random-effect covariance enter
#' through BayesTools' metadata-defined marginal \eqn{ZGZ'} covariance. The same
#' four-field comparison policy applies across sampled, marginalized, and mixed
#' local-effect representations.
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

  # compute the log-likelihood matrix (S x K)
  log_lik <- .log_lik.brma(object, unit = unit, caller = "add_waic()")
  target  <- attr(log_lik, "RoBMA_target", exact = TRUE)

  # call waic on the log-likelihood matrix
  waic_result <- loo::waic(log_lik, ...)
  waic_result <- .add_loo_target_metadata(
    object           = waic_result,
    unit             = target[["unit"]],
    retained_context = target[["retained_context"]],
    targets          = target[["targets"]],
    data_hash        = target[["data_hash"]],
    metadata         = target
  )

  # store in object and return
  if (is.null(object[["waic"]])) {
    object[["waic"]] <- list()
  }
  object[["waic"]][[unit]] <- waic_result
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


#' @title Extract Pointwise Log-Likelihood Draws
#'
#' @description Extract posterior pointwise log-likelihood draws from a fitted
#' \code{brma} object. The result is an \eqn{S \times K} or \eqn{S \times G}
#' matrix, where \eqn{S} is the number of posterior draws, \eqn{K} is the
#' number of estimates, and \eqn{G} is the number of clusters.
#'
#' @param object a brma model object.
#' @param unit output unit. See \code{\link{add_loo}}.
#' @param ... currently unused.
#'
#' @details
#' This function returns posterior log-likelihood draws, not a scalar maximized
#' log-likelihood. RoBMA therefore does not implement \code{stats::logLik()}.
#'
#' The log-likelihood is computed for each target at each posterior sample.
#' For binomial and Poisson models, each observation consists of a pair of
#' counts (ai/ci or x1i/x2i) that together define a single effect size estimate.
#' For correlated known-\code{V} \code{brma.mv()} models, estimate-unit columns
#' are Schur conditional scores \eqn{p(y_i \mid y_{-i}, \theta)}. Their row sum
#' is not presented as a full joint likelihood.
#' Estimate-unit GLMM draws support point or beta baserate priors and point or
#' normal log-rate priors, including their supported truncations. Other
#' nuisance-prior families are rejected when the model input is validated.
#' Cluster-unit GLMM draws are unavailable until certified nested
#' adaptive quadrature is implemented.
#'
#' @return An \eqn{S \times K} or \eqn{S \times G} matrix of pointwise
#' log-likelihood draws with \code{RoBMA_target} metadata.
#'
#' @seealso \code{\link{loo.brma}}, \code{\link{add_waic}}
#'
#' @export
log_lik <- function(object, ...) UseMethod("log_lik")


#' @rdname log_lik
#' @export
log_lik.brma <- function(object, unit = "estimate", ...) {

  return(.log_lik.brma(object, unit = unit, caller = "log_lik()"))
}


#' @title Unsupported Information Criteria for brma Objects
#'
#' @description `brma` objects expose posterior pointwise log-likelihood draws,
#' not a maximized scalar log-likelihood. AIC and BIC are therefore undefined
#' and these methods always stop with an explicit error.
#'
#' @param object a brma model object.
#' @param ... ignored.
#' @param k ignored; retained for the \code{stats::AIC()} method signature.
#'
#' @return These methods do not return a value.
#'
#' @seealso \code{\link{log_lik}}, \code{\link{loo.brma}},
#'   \code{\link{waic.brma}}
#'
#' @rdname information_criteria_brma
#' @export
AIC.brma <- function(object, ..., k = 2) {

  stop(
    "AIC() is not defined for brma objects because RoBMA does not expose a ",
    "scalar maximized likelihood. Use loo() or waic() for predictive model ",
    "assessment.",
    call. = FALSE
  )
}


#' @rdname information_criteria_brma
#' @export
BIC.brma <- function(object, ...) {

  stop(
    "BIC() is not defined for brma objects because RoBMA does not expose a ",
    "scalar maximized likelihood. Use loo() or waic() for predictive model ",
    "assessment.",
    call. = FALSE
  )
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
#' predictive performance for the documented deletion target, not in-sample
#' fit. For dependent estimates, deleting one estimate may retain other
#' observations from its fitted dependency block; this is not new-group
#' prediction.
#' RoBMA rejects comparisons with different outcome targets/data, \code{unit},
#' or retained context. The same compatibility key is applied
#' separately to LOO and WAIC; the two criteria cannot be mixed in one
#' comparison table. Sampled versus marginalized local-effect provenance does
#' not by itself block comparison.
#'
#' @return A numeric matrix of class \code{"compare.loo"}. The matrix retains
#' the comparison columns and printing contract returned by loo versions before
#' 2.10.0.
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
  return(.loo_compare_objects(c(list(x), list(...)), unit))
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
#' @return A numeric matrix of class \code{"compare.loo"}. The matrix retains
#' the comparison columns and printing contract returned by loo versions before
#' 2.10.0.
#'
#' @seealso \code{\link{loo.brma}}, \code{\link[loo]{loo_compare}}
#'
#' @export
loo_compare.loo <- function(x, ..., unit = "estimate") {
  return(.loo_compare_objects(c(list(x), list(...)), unit))
}


.loo_compare_objects <- function(models, unit) {

  unit <- .normalize_unit(unit)

  if (length(models) < 2) {
    stop("At least two models are required for comparison.", call. = FALSE)
  }

  loo_objects <- lapply(models, function(m) {
    if (inherits(m, "brma")) {
      loo.brma(m, unit = unit)
    } else if (inherits(m, "loo")) {
      m
    } else {
      stop("All arguments must be brma or loo objects.", call. = FALSE)
    }
  })

  .check_loo_compare_targets(loo_objects)
  loo_compare_fun <- get("loo_compare.default", envir = asNamespace("loo"), inherits = FALSE)
  result <- .as_legacy_loo_compare(do.call(loo_compare_fun, loo_objects))

  return(result)
}


# Preserve RoBMA's released numeric loo-comparison contract across loo 2.10.
.as_legacy_loo_compare <- function(x) {

  if (!is.data.frame(x)) {
    stop("Internal error: unsupported loo comparison table.", call. = FALSE)
  }
  metric_columns <- list(
    loo = c(
      "elpd_loo", "se_elpd_loo", "p_loo", "se_p_loo",
      "looic", "se_looic"
    ),
    waic = c(
      "elpd_waic", "se_elpd_waic", "p_waic", "se_p_waic",
      "waic", "se_waic"
    )
  )
  matched_metric <- names(metric_columns)[vapply(
    metric_columns,
    function(columns) all(columns %in% colnames(x)),
    logical(1)
  )]
  if (length(matched_metric) != 1L || !"model" %in% colnames(x)) {
    stop("Internal error: unsupported loo comparison table.", call. = FALSE)
  }
  legacy_columns <- c("elpd_diff", "se_diff", metric_columns[[matched_metric]])
  if (!all(legacy_columns %in% colnames(x))) {
    stop("Internal error: incomplete loo comparison table.", call. = FALSE)
  }

  out <- as.matrix(x[, legacy_columns, drop = FALSE])
  if (!is.numeric(out) || anyNA(x[["model"]]) ||
      any(!nzchar(as.character(x[["model"]])))) {
    stop("Internal error: invalid loo comparison table.", call. = FALSE)
  }
  rownames(out) <- as.character(x[["model"]])
  class(out)    <- c("compare.loo", "matrix", "array")

  return(out)
}


.print_compare_loo <- function(x, ..., digits = 1, simplify = TRUE) {

  if (is.data.frame(x)) {
    print_method <- get(
      "print.compare.loo",
      envir    = asNamespace("loo"),
      inherits = FALSE
    )
    return(print_method(x, ..., digits = digits))
  }

  x_copy <- x
  if (ncol(x_copy) >= 2L && isTRUE(simplify)) {
    x_copy <- x_copy[, c("elpd_diff", "se_diff"), drop = FALSE]
  }
  print(
    format(round(x_copy, digits), nsmall = digits),
    quote = FALSE
  )

  invisible(x)
}


# ---------------------------------------------------------------------------- #
# loo_model_weights method
# ---------------------------------------------------------------------------- #

#' @title Predictive Model Weights for brma Objects
#'
#' @description Compute stacking or pseudo-BMA weights from LOO results already
#' stored in multiple \code{brma} objects.
#'
#' @param x a \code{brma} model object (the first model to weight).
#' @param ... additional \code{brma} model objects or RoBMA-targeted
#'   \code{psis_loo} objects.
#' @param unit output/deletion unit used when extracting stored LOO results from
#'   \code{brma} objects.
#' @param method weighting method; \code{"stacking"} or \code{"pseudobma"}.
#' @param optim_method,optim_control optimization method and controls for
#'   stacking. See \code{\link[loo]{loo_model_weights}}.
#' @param BB,BB_n,alpha Bayesian-bootstrap controls for pseudo-BMA weighting.
#'   See \code{\link[loo]{loo_model_weights}}.
#' @param r_eff_list optional relative effective sample sizes. Retained for the
#'   upstream interface and ignored when stored \code{psis_loo} objects are
#'   supplied.
#' @param cores number of cores passed to \code{loo::loo_model_weights()}.
#'
#' @details Each \code{brma} object must first be updated with
#' \code{object <- add_loo(object, unit = unit)}. All models must use the same
#' outcome data, deletion \code{unit}, and likelihood target. RoBMA validates
#' these targets before delegating to \code{loo::loo_model_weights()}.
#'
#' Model names are derived from the supplied model expressions. See
#' \code{\link[loo]{loo_model_weights}} for details about stacking, pseudo-BMA,
#' and pseudo-BMA+.
#'
#' @return A named numeric vector containing one predictive weight per model.
#'
#' @seealso \code{\link{add_loo}}, \code{\link{loo.brma}},
#'   \code{\link{loo_compare.brma}}, \code{\link[loo]{loo_model_weights}}
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'
#'   fit_bias <- RoBMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'   fit_nobias <- BMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   fit_bias   <- add_loo(fit_bias)
#'   fit_nobias <- add_loo(fit_nobias)
#'
#'   loo_model_weights(fit_bias, fit_nobias)
#' }
#' }
#'
#' @export
#' @exportS3Method loo::loo_model_weights
loo_model_weights.brma <- function(
    x, ..., unit = "estimate", method = c("stacking", "pseudobma"),
    optim_method = "BFGS", optim_control = list(), BB = TRUE, BB_n = 1000,
    alpha = 1, r_eff_list = NULL, cores = getOption("mc.cores", 1)) {

  unit   <- .normalize_unit(unit)
  dots   <- list(...)
  models <- c(list(x), dots)

  if (length(models) < 2L) {
    stop("At least two models are required for predictive model weights.",
         call. = FALSE)
  }

  loo_objects <- lapply(models, function(model) {
    if (inherits(model, "brma")) {
      return(loo.brma(model, unit = unit))
    }
    if (inherits(model, "psis_loo")) {
      return(model)
    }
    stop("All arguments must be brma or psis_loo objects.", call. = FALSE)
  })

  model_call        <- match.call(expand.dots = FALSE)
  model_expressions <- c(
    list(model_call[["x"]]),
    as.list(model_call[["..."]])
  )
  model_names <- vapply(model_expressions, deparse1, character(1))
  dot_names   <- names(dots)
  named_dots  <- which(!is.na(dot_names) & nzchar(dot_names))
  if (length(named_dots) > 0L) {
    model_names[named_dots + 1L] <- dot_names[named_dots]
  }
  names(loo_objects) <- model_names

  .check_loo_compare_targets(loo_objects)

  draw_counts   <- vapply(
    loo_objects, function(object) dim(object)[[1L]], integer(1)
  )
  target_counts <- vapply(
    loo_objects, function(object) dim(object)[[2L]], integer(1)
  )

  # TODO: Once https://github.com/stan-dev/loo/issues/389 is fixed and the
  # corresponding loo release is required, remove this unequal-draw workaround
  # and delegate all cases to loo_model_weights.default().
  if (length(unique(draw_counts)) > 1L &&
      length(unique(target_counts)) == 1L) {
    method <- match.arg(method)
    lpd_point <- vapply(
      loo_objects,
      function(object) object[["pointwise"]][, "elpd_loo"],
      numeric(target_counts[[1L]])
    )

    if (method == "stacking") {
      weights <- loo::stacking_weights(
        lpd_point     = lpd_point,
        optim_method  = optim_method,
        optim_control = optim_control
      )
    } else {
      weights <- loo::pseudobma_weights(
        lpd_point = lpd_point,
        BB        = BB,
        BB_n      = BB_n,
        alpha     = alpha
      )
    }
    names(weights) <- names(loo_objects)

    return(weights)
  }

  loo_model_weights_fun <- get(
    "loo_model_weights.default",
    envir    = asNamespace("loo"),
    inherits = FALSE
  )

  return(loo_model_weights_fun(
    x             = loo_objects,
    method        = method,
    optim_method  = optim_method,
    optim_control = optim_control,
    BB            = BB,
    BB_n          = BB_n,
    alpha         = alpha,
    r_eff_list    = r_eff_list,
    cores         = cores
  ))
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
#' stored PSIS object and does not compute LOO. It does not compute weights
#' across models. For predictive model weights, use
#' \code{loo_model_weights(model_1, model_2)}; for a model comparison table,
#' use \code{loo_compare(model_1, model_2)}.
#'
#' @return An \code{S x K} matrix for estimate-unit LOO, or \code{S x G}
#' matrix for cluster-unit LOO, with posterior samples in rows and LOO targets
#' in columns. Columns are normalized to sum to one.
#'
#' @seealso \code{\link{loo.brma}}, \code{\link{loo_compare.brma}},
#'   \code{\link{loo_model_weights.brma}}
#'
#' @aliases loo_weights
#' @exportS3Method
loo_weights.brma <- function(object, unit = "estimate", ...) {

  if (inherits(unit, c("brma", "loo"))) {
    stop(
      "loo_weights() extracts PSIS importance weights from one brma model; ",
      "it does not compute weights across models. Use ",
      "loo_model_weights(model_1, model_2) for predictive ",
      "model weights, or loo_compare(model_1, model_2) to compare models.",
      call. = FALSE
    )
  }

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
