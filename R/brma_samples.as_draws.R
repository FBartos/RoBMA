# ---------------------------------------------------------------------------- #
# as_draws methods for posterior package integration
# ---------------------------------------------------------------------------- #

#' @title Convert brma_samples to posterior Draws Formats
#'
#' @description Provides an interface to the \pkg{posterior} package
#' for \code{brma_samples} objects. These functions convert the posterior
#' samples to various draws formats supported by the \pkg{posterior} package.
#'
#' @param x a \code{brma_samples} object
#' @param ... additional arguments passed to the corresponding
#' \pkg{posterior} function
#'
#' @details
#' The conversion reconstructs the MCMC chain structure from the stored
#' \code{n_chains} and \code{n_iter} attributes. The samples are assumed
#' to be ordered with chains concatenated (i.e., all iterations from chain 1,
#' then all from chain 2, etc.).
#'
#' @return An object of the corresponding \pkg{posterior} draws class.
#'
#' @seealso \code{\link[posterior]{draws}}, \code{\link{as_draws.brma}}
#'
#' @name as_draws.brma_samples
NULL

#' @rdname as_draws.brma_samples
#' @export
as_draws.brma_samples <- function(x, ...) {

  .check_posterior_package()

  # reconstruct mcmc.list from samples, convert to draws_array first
  # (to properly set chain dimensions), then to draws_matrix
  # this matches posterior::as_draws.matrix behavior but preserves chain info
  mcmc.list <- .brma_samples_to_mcmc.list(x)

  return(posterior::as_draws_matrix(mcmc.list, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_array.brma_samples <- function(x, ...) {

  .check_posterior_package()

  # reconstruct mcmc.list from samples
  mcmc.list <- .brma_samples_to_mcmc.list(x)

  return(posterior::as_draws_array(mcmc.list, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_df.brma_samples <- function(x, ...) {

  .check_posterior_package()

  # reconstruct mcmc.list from samples
  mcmc.list <- .brma_samples_to_mcmc.list(x)

  return(posterior::as_draws_df(mcmc.list, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_list.brma_samples <- function(x, ...) {

  .check_posterior_package()

  # reconstruct mcmc.list from samples
  mcmc.list <- .brma_samples_to_mcmc.list(x)

  return(posterior::as_draws_list(mcmc.list, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_matrix.brma_samples <- function(x, ...) {

  .check_posterior_package()

  # reconstruct mcmc.list from samples
  mcmc.list <- .brma_samples_to_mcmc.list(x)

  return(posterior::as_draws_matrix(mcmc.list, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_rvars.brma_samples <- function(x, ...) {

  .check_posterior_package()

  # reconstruct mcmc.list from samples
  mcmc.list <- .brma_samples_to_mcmc.list(x)

  return(posterior::as_draws_rvars(mcmc.list, ...))
}


# ---------------------------------------------------------------------------- #
# Helper: Convert brma_samples to mcmc.list
# ---------------------------------------------------------------------------- #

#' @title Convert brma_samples to mcmc.list
#'
#' @description Internal helper to reconstruct an \code{mcmc.list} object
#' from a \code{brma_samples} matrix by splitting samples back into chains.
#'
#' @param x a \code{brma_samples} object
#'
#' @return An \code{mcmc.list} object suitable for conversion to posterior
#' draws formats.
#'
#' @keywords internal
.brma_samples_to_mcmc.list <- function(x) {

  n_chains <- attr(x, "nchains")
  n_iter   <- attr(x, "niter")

  # convert to plain matrix
  samples_mat <- as.matrix(x)

  # split into chains
  # samples are ordered: chain1_iter1, ..., chain1_iterN, chain2_iter1, ..., etc.
  chains <- vector("list", n_chains)
  for (i in seq_len(n_chains)) {
    start_idx <- (i - 1) * n_iter + 1
    end_idx   <- i * n_iter
    chains[[i]] <- coda::mcmc(samples_mat[start_idx:end_idx, , drop = FALSE])
  }

  return(coda::mcmc.list(chains))
}
