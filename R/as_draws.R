# ============================================================================ #
# brma.as_draws.R
# ============================================================================ #
#
# This file provides an interface to the posterior package for brma objects.
# It defines:
#
# - Generic S3 methods for as_draws family (for use when posterior is not loaded)
# - Methods for brma class that convert MCMC samples to posterior draws formats
# - Helper function to extract mcmc.list from brma objects
#
# ============================================================================ #


# ---------------------------------------------------------------------------- #
# Documentation
# ---------------------------------------------------------------------------- #

#' @title Convert brma Objects to posterior Draws Formats
#'
#' @description Provides an interface to the \pkg{posterior} package
#' for \code{brma} objects. These functions convert the MCMC samples
#' from a fitted \code{brma} model to various draws formats supported
#' by the \pkg{posterior} package, enabling the use of posterior's
#' rich set of diagnostics and summary functions.
#'
#' @param x a fitted \code{brma} object.
#' @param ... additional arguments passed to the corresponding
#' \pkg{posterior} function.
#'
#' @details
#' These functions are S3 generics. Their default methods forward to the
#' corresponding \pkg{posterior} generics so attaching \pkg{RoBMA} preserves
#' the usual \pkg{posterior} behavior for non-\code{brma} objects.
#'
#' The following conversion functions are available:
#' \itemize{
#'   \item \code{as_draws}: converts to the default draws format
#'   \item \code{as_draws_array}: converts to a 3-D array (iteration x chain x variable)
#'   \item \code{as_draws_df}: converts to a data frame with columns for iterations, chains, and variables
#'   \item \code{as_draws_list}: converts to a list of lists
#'   \item \code{as_draws_matrix}: converts to a 2-D matrix (draw x variable)
#'   \item \code{as_draws_rvars}: converts to random variable objects
#' }
#'
#' These methods require the \pkg{posterior} package to be installed.
#' The conversion is performed by first extracting the MCMC samples
#' as a \code{mcmc.list} object and then using the corresponding
#' \pkg{posterior} conversion function.
#'
#' @return An object of the corresponding \pkg{posterior} draws class.
#'
#' @seealso \code{\link[posterior]{draws}}, \code{\link{brma}},
#' \code{\link{as_draws.brma_samples}}
#'
#' @name as_draws.brma
NULL


# ---------------------------------------------------------------------------- #
# Helper: Extract mcmc.list from brma objects
# ---------------------------------------------------------------------------- #

#' @title Extract mcmc.list from brma Objects
#'
#' @description Internal helper to extract the MCMC samples from a fitted
#' \code{brma} object as a \code{mcmc.list} object suitable for conversion
#' to \pkg{posterior} draws formats.
#'
#' @param x a fitted \code{brma} object
#'
#' @return A \code{mcmc.list} object containing the MCMC samples.
#'
#' @noRd
.brma_to_mcmc.list <- function(x) {

  if (!inherits(x, "brma")) {
    stop("'x' must be a 'brma' object", call. = FALSE)
  }

  if (is.null(x[["fit"]]) || length(x[["fit"]]) == 0) {
    stop("The 'brma' object does not contain a valid fit", call. = FALSE)
  }

  return(coda::as.mcmc.list(x[["fit"]]))
}


# ---------------------------------------------------------------------------- #
# Helper: Check posterior package availability
# ---------------------------------------------------------------------------- #

.check_posterior_package <- function() {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "Package 'posterior' is required for as_draws conversion.\n",
      "Install it with: install.packages('posterior')",
      call. = FALSE
    )
  }
}

.posterior_default_method <- function(generic) {

  return(getS3method(
    f     = generic,
    class = "default",
    envir = asNamespace("posterior")
  ))
}

#' @rdname as_draws.brma
#' @export
as_draws <- function(x, ...) {
  UseMethod("as_draws")
}

#' @rdname as_draws.brma
#' @export
as_draws_array <- function(x, ...) {
  UseMethod("as_draws_array")
}

#' @rdname as_draws.brma
#' @export
as_draws_df <- function(x, ...) {
  UseMethod("as_draws_df")
}

#' @rdname as_draws.brma
#' @export
as_draws_list <- function(x, ...) {
  UseMethod("as_draws_list")
}

#' @rdname as_draws.brma
#' @export
as_draws_matrix <- function(x, ...) {
  UseMethod("as_draws_matrix")
}

#' @rdname as_draws.brma
#' @export
as_draws_rvars <- function(x, ...) {
  UseMethod("as_draws_rvars")
}

#' @rdname as_draws.brma
#' @export
as_draws.default <- function(x, ...) {

  .check_posterior_package()

  return(.posterior_default_method("as_draws")(x, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_array.default <- function(x, ...) {

  .check_posterior_package()

  return(.posterior_default_method("as_draws_array")(x, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_df.default <- function(x, ...) {

  .check_posterior_package()

  return(.posterior_default_method("as_draws_df")(x, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_list.default <- function(x, ...) {

  .check_posterior_package()

  return(.posterior_default_method("as_draws_list")(x, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_matrix.default <- function(x, ...) {

  .check_posterior_package()

  return(.posterior_default_method("as_draws_matrix")(x, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_rvars.default <- function(x, ...) {

  .check_posterior_package()

  return(.posterior_default_method("as_draws_rvars")(x, ...))
}


# ---------------------------------------------------------------------------- #
# as_draws methods for brma class
# ---------------------------------------------------------------------------- #

#' @rdname as_draws.brma
#' @export
as_draws.brma <- function(x, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(x)

  return(posterior::as_draws(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_array.brma <- function(x, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(x)

  return(posterior::as_draws_array(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_df.brma <- function(x, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(x)

  return(posterior::as_draws_df(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_list.brma <- function(x, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(x)

  return(posterior::as_draws_list(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_matrix.brma <- function(x, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(x)

  return(posterior::as_draws_matrix(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_rvars.brma <- function(x, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(x)

  return(posterior::as_draws_rvars(mcmc.list, ...))
}
