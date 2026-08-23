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
#' \code{nchains} and \code{niter} attributes. The samples are assumed
#' to be ordered with chains concatenated (i.e., all iterations from chain 1,
#' then all from chain 2, etc.). Conditional RoBMA samples are intentionally
#' stored as one flattened chain because conditioning subsets posterior rows
#' across chains.
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
  return(.brma_samples_as_draws(x, posterior::as_draws_matrix, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_array.brma_samples <- function(x, ...) {
  return(.brma_samples_as_draws(x, posterior::as_draws_array, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_df.brma_samples <- function(x, ...) {
  return(.brma_samples_as_draws(x, posterior::as_draws_df, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_list.brma_samples <- function(x, ...) {
  return(.brma_samples_as_draws(x, posterior::as_draws_list, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_matrix.brma_samples <- function(x, ...) {
  return(.brma_samples_as_draws(x, posterior::as_draws_matrix, ...))
}

#' @rdname as_draws.brma_samples
#' @export
as_draws_rvars.brma_samples <- function(x, ...) {
  return(.brma_samples_as_draws(x, posterior::as_draws_rvars, ...))
}


.brma_samples_as_draws <- function(x, converter, ...) {

  .check_posterior_package()
  mcmc.list <- .brma_samples_to_mcmc.list(x)

  return(converter(mcmc.list, ...))
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
#' @noRd
.brma_samples_to_mcmc.list <- function(x) {

  n_chains <- attr(x, "nchains")
  n_iter   <- attr(x, "niter")

  # convert to plain matrix
  samples_mat <- as.matrix(x)

  if (length(n_chains) != 1L || length(n_iter) != 1L ||
      is.na(n_chains) || is.na(n_iter) ||
      n_chains < 1L || n_iter < 1L ||
      n_chains * n_iter != nrow(samples_mat)) {
    stop("Invalid brma_samples chain metadata: 'nchains * niter' must equal the number of sample rows.",
         call. = FALSE)
  }

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
