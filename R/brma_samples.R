# ============================================================================ #
# brma_samples.R
# ============================================================================ #
#
# This file defines the brma_samples class for posterior samples returned by
# predict.brma() and related wrapper functions. The class:
#
# - Stores posterior samples as a matrix (S x K)
# - Retains MCMC metadata (nchains, niter) for posterior package integration
# - Provides print/summary methods using BayesTools::ensemble_estimates_table
# - Supports as_draws conversion for posterior package compatibility
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# Constructor
# ---------------------------------------------------------------------------- #

#' @title Create a brma_samples Object
#'
#' @description Internal constructor for creating \code{brma_samples} objects
#' that store posterior samples with metadata for printing and conversion
#' to \pkg{posterior} draws formats.
#'
#' @param samples a matrix of posterior samples with dimensions S x K
#' (samples x parameters/observations)
#' @param n_chains number of MCMC chains used in sampling
#' @param n_iter number of iterations per chain (after thinning)
#' @param title title for the summary table output
#' @param probs default quantiles for credible intervals. Defaults to
#' \code{c(.025, .975)}
#' @param data optional data associated with the predictions (e.g., for
#' non-aggregated predictions)
#' @param effect_transform optional effect-size transformation metadata
#'
#' @return An object of class \code{brma_samples} which inherits from
#' \code{matrix}.
#'
#' @keywords internal
.new_brma_samples <- function(samples, n_chains, n_iter, title,
                              probs = c(.025, .975), data = NULL,
                              effect_transform = NULL) {

  # ensure samples is a matrix with proper column names
  if (!is.matrix(samples)) {
    samples <- as.matrix(samples)
  }

  # add attributes for MCMC structure
  attr(samples, "nchains") <- n_chains
  attr(samples, "niter")   <- n_iter

  # add attributes for display
  attr(samples, "title")    <- title
  attr(samples, "probs")    <- probs
  attr(samples, "data")     <- data

  if (!is.null(effect_transform)) {
    attr(samples, "effect_transform") <- effect_transform
    attr(samples, "footnotes")        <- effect_transform[["note"]]
  }

  # set class (inherits from matrix for backward compatibility)
  class(samples) <- c("brma_samples", "matrix", "array")

  return(samples)
}


# ---------------------------------------------------------------------------- #
# Print method
# ---------------------------------------------------------------------------- #

#' @title Print brma_samples Object
#'
#' @description Prints a summary table of posterior samples using
#' \code{BayesTools::ensemble_estimates_table}. Returns the samples
#' invisibly for method chaining.
#'
#' @param x a \code{brma_samples} object
#' @param probs quantiles for credible intervals. If \code{NULL}, uses
#' the default stored in the object (typically \code{c(.025, .975)})
#' @param ... additional arguments passed to
#' \code{BayesTools::ensemble_estimates_table}
#'
#' @return Returns \code{x} invisibly.
#'
#' @export
print.brma_samples <- function(x, probs = NULL, ...) {

  summary <- summary.brma_samples(x, probs = probs, ...)
  print.summary.brma_samples(summary)

  return(invisible(x))
}


# ---------------------------------------------------------------------------- #
# Summary method
# ---------------------------------------------------------------------------- #

#' @title Summarize brma_samples Object
#'
#' @description Creates, prints, and returns a summary table of posterior samples
#' using \code{BayesTools::ensemble_estimates_table}.
#'
#' @param object a \code{brma_samples} object
#' @param probs quantiles for credible intervals. If \code{NULL}, uses
#' the default stored in the object (typically \code{c(.025, .975)})
#' @param ... additional arguments passed to
#' \code{BayesTools::ensemble_estimates_table}
#'
#' @return A \code{BayesTools_table} object containing the summary statistics.
#'
#' @export
summary.brma_samples <- function(object, probs = NULL, ...) {

  # use stored probs if not provided
  if (is.null(probs)) {
    probs <- attr(object, "probs")
  }

  dots <- list(...)
  if (is.null(dots[["footnotes"]]) && !is.null(attr(object, "footnotes"))) {
    dots[["footnotes"]] <- attr(object, "footnotes")
  }

  # create summary table
  summary_table <- do.call(
    BayesTools::ensemble_estimates_table,
    c(
      list(
        samples    = asplit(object, 2),
        parameters = colnames(object),
        probs      = probs,
        title      = attr(object, "title")
      ),
      dots
    )
  )

  class(summary_table) <- c("summary.brma_samples", class(summary_table))
 
  return(summary_table)
}


# ---------------------------------------------------------------------------- #
# Print summary method
# ---------------------------------------------------------------------------- #

#' @title Print summary.brma_samples Object
#'
#' @description Prints a summary table of posterior samples using
#' \code{BayesTools::ensemble_estimates_table}. Returns the summary table
#' invisibly.
#'
#' @param x a \code{summary.brma_samples} object
#' @param probs quantiles for credible intervals. If \code{NULL}, uses
#' the default stored in the object (typically \code{c(.025, .975)})
#' @param ... additional arguments passed to
#' \code{BayesTools::ensemble_estimates_table}
#'
#' @return Returns the summary table invisibly.
#'
#' @export
print.summary.brma_samples <- function(x, probs = NULL, ...) {

  class(x) <- setdiff(class(x), "summary.brma_samples")

  cat("\n")
  print(x)
  cat("\n")

  return(invisible(x))
}


# ---------------------------------------------------------------------------- #
# as.matrix method (for backward compatibility)
# ---------------------------------------------------------------------------- #

#' @title Convert brma_samples to Matrix
#'
#' @description Converts a \code{brma_samples} object to a plain matrix,
#' removing all brma_samples-specific attributes.
#'
#' @param x a \code{brma_samples} object
#' @param ... additional arguments (ignored)
#'
#' @return A plain matrix of posterior samples.
#'
#' @export
as.matrix.brma_samples <- function(x, ...) {
  # remove class and custom attributes
  attributes(x) <- list(
    dim      = dim(x),
    dimnames = dimnames(x)
  )
  return(x)
}
