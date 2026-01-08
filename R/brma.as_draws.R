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
#' @seealso \code{\link[posterior]{draws}}, \code{\link{brma}}
#'
#' @name as_draws.brma
NULL

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
as_draws.brma <- function(x, ...) {

  require("posterior")

  # extract the mcmc samples
  mcmc.list <- coda::as.mcmc.list(x[["fit"]])

  return(posterior::as_draws(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_array.brma <- function(x, ...) {

  require("posterior")

  # extract the mcmc samples
  mcmc.list <- coda::as.mcmc.list(x[["fit"]])

  return(posterior::as_draws_array(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_df.brma <- function(x, ...) {

  require("posterior")

  # extract the mcmc samples
  mcmc.list <- coda::as.mcmc.list(x[["fit"]])

  return(posterior::as_draws_df(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_list.brma <- function(x, ...) {

  require("posterior")

  # extract the mcmc samples
  mcmc.list <- coda::as.mcmc.list(x[["fit"]])

  return(posterior::as_draws_list(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_matrix.brma <- function(x, ...) {

  require("posterior")

  # extract the mcmc samples
  mcmc.list <- coda::as.mcmc.list(x[["fit"]])

  return(posterior::as_draws_matrix(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_rvars.brma <- function(x, ...) {

  require("posterior")

  # extract the mcmc samples
  mcmc.list <- coda::as.mcmc.list(x[["fit"]])

  return(posterior::as_draws_rvars(mcmc.list, ...))
}
