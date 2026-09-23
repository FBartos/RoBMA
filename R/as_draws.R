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
#' @param x an object to convert. The \code{brma} methods expect a fitted
#' \code{brma} object; default methods forward non-\code{brma} objects to the
#' corresponding \pkg{posterior} conversion function.
#' @param include_auxiliary logical; whether to include raw backend auxiliary
#' variables in addition to the stable model-parameter schema. Formula-random
#' parameters retain their RoBMA semantic names in either case.
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
#' For \code{brma} methods, conversion is performed by first extracting the MCMC
#' samples as a \code{mcmc.list} object and then using the corresponding
#' \pkg{posterior} conversion function. By default, backend-private latent
#' effects, known-\code{V} dependency factors, random-effect simulation and
#' Cholesky factors, and prior-parameterization variables are omitted. Stable
#' fitted coordinates such as model indicators and likelihood parameters
#' remain available. Formula-random coordinates are replaced by their semantic
#' RoBMA quantities using the same names as summaries, plots, and hypotheses,
#' such as `tau_total`, `tau_common`, `rho(...)`, and `tau2_prop(...)`.
#' Random-effect SD switches are exposed as exact columns named after the
#' standard deviation they act on, with an `_indicator` suffix, pairing with
#' the semantic SD column the same way `mu_x_indicator` pairs with `mu_x`. A
#' heterogeneity-component inclusion gate is 0/1, such as `study: tau_indicator`
#' and `observation: tau_indicator`, and carries the inclusion state the
#' summary reports under `Heterogeneity:`. An allocation-SD prior indicator,
#' such as `tau_total_indicator`, records which mixture component allocated the
#' total standard deviation and therefore takes one value per component,
#' matching the order of that prior's components.
#' Structural point-prior coordinates are returned as exact constant columns,
#' with chain timing taken from the fitted draw geometry. Private backend
#' anchors are never exposed. Use `include_auxiliary = TRUE` to skip RoBMA
#' auxiliary filtering; added backend variables retain their coordinate names,
#' while public random-effect parameters retain their RoBMA names.
#' \code{brma_samples} objects have
#' separate methods documented at \code{\link{as_draws.brma_samples}}.
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
#' @param include_auxiliary logical; whether to retain raw backend auxiliary
#' variables.
#'
#' @return A \code{mcmc.list} object containing the MCMC samples.
#'
#' @noRd
.brma_to_mcmc.list <- function(x, include_auxiliary = FALSE) {

  if (!inherits(x, "brma")) {
    stop("'x' must be a 'brma' object", call. = FALSE)
  }
  BayesTools::check_bool(include_auxiliary, "include_auxiliary")

  if (is.null(x[["fit"]]) || length(x[["fit"]]) == 0) {
    stop("The 'brma' object does not contain a valid fit", call. = FALSE)
  }

  .brma_validate_fit_contract(
    x,
    requires = c("parameter_map", "draw_geometry")
  )
  mcmc_list <- BayesTools::JAGS_materialize_draws(x[["fit"]])
  if (include_auxiliary) {
    if (.is_random(x)) {
      return(.brma_replace_random_coordinates(x, mcmc_list))
    }
    return(mcmc_list)
  }

  mcmc_list <- .brma_filter_auxiliary_variables(x, mcmc_list)
  if (.is_random(x)) {
    return(.brma_replace_random_coordinates(x, mcmc_list))
  }
  return(mcmc_list)
}


# Replace fitted formula-random coordinates by RoBMA semantic quantities.
.brma_replace_random_coordinates <- function(x, mcmc_list) {

  coordinates <- BayesTools::parameter_coordinates(x[["fit"]])
  random_coordinates <- coordinates[
    coordinates[["role"]] %in% c("random_sd", "random_correlation"),
    "coordinate_name",
    drop = TRUE
  ]
  bundle          <- .brma_random_parameter_bundle(x, chains = TRUE)
  semantic_chains <- bundle[["samples"]]
  semantic_names  <- bundle[["specs"]][["label"]]
  # Heterogeneity-component inclusion gates live on private coordinates that
  # JAGS_materialize_draws() hides. Re-expose them beside the SD they switch on.
  gate_chains     <- .brma_random_inclusion_indicator_chains(x)
  if (length(semantic_names) != ncol(semantic_chains[[1L]])) {
    stop("Semantic random-effect draw names are inconsistent.", call. = FALSE)
  }

  chains <- lapply(seq_along(mcmc_list), function(chain_i) {
    coordinate_chain <- mcmc_list[[chain_i]]
    semantic_chain   <- semantic_chains[[chain_i]]
    if (nrow(coordinate_chain) != nrow(semantic_chain) ||
        !identical(coda::mcpar(coordinate_chain), coda::mcpar(semantic_chain))) {
      stop(
        "Semantic random-effect draws disagree with the fitted draw geometry.",
        call. = FALSE
      )
    }

    coordinate_values <- as.matrix(coordinate_chain)
    coordinate_values <- coordinate_values[
      , !colnames(coordinate_values) %in% random_coordinates,
      drop = FALSE
    ]
    semantic_values <- as.matrix(semantic_chain)
    colnames(semantic_values) <- semantic_names
    if (length(gate_chains) > 0L) {
      gate_values <- gate_chains[[chain_i]]
      if (nrow(gate_values) != nrow(semantic_values)) {
        stop(
          "Random-effect inclusion draws disagree with the fitted draw geometry.",
          call. = FALSE
        )
      }
      semantic_values <- cbind(semantic_values, gate_values)
    }
    duplicate_names <- intersect(
      colnames(coordinate_values),
      colnames(semantic_values)
    )
    if (length(duplicate_names) > 0L) {
      stop(
        "RoBMA semantic random-effect names conflict with fitted draw names: ",
        paste0("'", duplicate_names, "'", collapse = ", "), ".",
        call. = FALSE
      )
    }

    values <- cbind(coordinate_values, semantic_values)
    mcpar  <- coda::mcpar(coordinate_chain)
    coda::mcmc(
      values,
      start = mcpar[[1L]],
      end   = mcpar[[2L]],
      thin  = mcpar[[3L]]
    )
  })

  coda::mcmc.list(chains)
}


# Exact indicator draws of every random-effect SD switch, one matrix per chain,
# named after the SD it acts on: 0/1 for a component inclusion gate, and the
# mixture component index for an allocation-SD prior. They sit on private
# coordinates, so they must be materialized with the internal coordinates
# included. The values are passed through unchanged; re-encoding them would
# assert a meaning the prior does not carry.
.brma_random_inclusion_indicator_chains <- function(x) {

  gate_names <- c(.random_inclusion_sd_names(x), .random_slab_sd_names(x))
  if (length(gate_names) == 0L) {
    return(list())
  }

  chains <- BayesTools::JAGS_materialize_draws(
    x[["fit"]],
    parameters       = unname(names(gate_names)),
    include_internal = TRUE
  )
  available <- intersect(names(gate_names), colnames(as.matrix(chains[[1L]])))
  if (length(available) == 0L) {
    return(list())
  }

  lapply(chains, function(chain) {

    values <- as.matrix(chain)[, available, drop = FALSE]
    if (any(!is.finite(values)) || !all(values == trunc(values))) {
      stop("Random-effect indicator draws are not whole numbers.",
           call. = FALSE)
    }
    colnames(values) <- paste0(unname(gate_names[available]), "_indicator")
    values
  })
}


# Catalog backend-private variables present in a fitted brma object.
.brma_auxiliary_variable_catalog <- function(x, variables) {

  if (!inherits(x, "brma")) {
    stop("'x' must be a 'brma' object", call. = FALSE)
  }
  if (is.null(variables)) {
    variables <- character(0)
  }
  if (!is.character(variables) || anyNA(variables)) {
    stop("'variables' must be a character vector without missing values.",
         call. = FALSE)
  }

  variable_base <- sub("\\[.*$", "", variables)
  category      <- rep(NA_character_, length(variables))

  category[variable_base %in% c("theta", "gamma")] <- "latent_effect"
  category[variable_base == "sampling_z"] <- "sampling_dependency"
  category[endsWith(variable_base, "_xRE_Zx")] <- "random_effect_latent"
  category[endsWith(variable_base, "_xRE_CORx_L")] <- "random_correlation_factor"
  category[endsWith(variable_base, "_xRE_CORx_lkj_u")] <- "random_correlation_coordinate"
  category[endsWith(variable_base, "_xRE_CORx_lkj_cpc")] <- "random_correlation_coordinate"
  category[startsWith(variable_base, "prior_par_eta_")] <- "prior_parameterization"
  category[startsWith(variable_base, "inv_")] <- "transformed_prior"
  category[variable_base %in% .brma_spike_and_slab_auxiliary_bases(x)] <-
    "spike_and_slab_latent"
  category[variable_base %in% c(
    "eta",
    "log_omega",
    "alpha",
    "pi_null"
  )] <- "selection_latent"

  keep <- !is.na(category)
  return(data.frame(
    variable         = variables[keep],
    category         = category[keep],
    stringsAsFactors = FALSE
  ))
}


# Identify exact private coordinates generated by spike-and-slab priors.
.brma_spike_and_slab_auxiliary_bases <- function(x) {

  prior_list <- attr(x[["fit"]], "prior_list", exact = TRUE)
  if (is.null(prior_list) || length(prior_list) == 0L ||
      is.null(names(prior_list))) {
    return(character(0))
  }

  is_spike_and_slab <- vapply(prior_list, function(prior) {
    isTRUE(BayesTools::is.prior.spike_and_slab(prior))
  }, logical(1))
  parameters <- names(prior_list)[is_spike_and_slab & nzchar(names(prior_list))]

  return(c(
    paste0(parameters, "_variable"),
    paste0(parameters, "_inclusion")
  ))
}


# Remove auxiliary variables while retaining chain timing and raw values.
.brma_filter_auxiliary_variables <- function(x, mcmc_list) {

  chain_variables <- lapply(mcmc_list, function(chain) {
    colnames(chain)
  })
  variables <- chain_variables[[1L]]
  if (is.null(variables)) {
    variables <- character(0)
    chain_variables <- lapply(chain_variables, function(x) {
      if (is.null(x)) character(0) else x
    })
  }
  if (!all(vapply(
    chain_variables,
    identical,
    logical(1),
    y = variables
  ))) {
    stop("MCMC chains contain inconsistent variable schemas.", call. = FALSE)
  }
  if (length(variables) == 0L) {
    return(mcmc_list)
  }

  catalog <- .brma_auxiliary_variable_catalog(x, variables)
  keep    <- !variables %in% catalog[["variable"]]
  chains  <- lapply(mcmc_list, function(chain) {
    chain[, keep, drop = FALSE]
  })

  return(coda::mcmc.list(chains))
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
as_draws.brma <- function(x, include_auxiliary = FALSE, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(
    x,
    include_auxiliary = include_auxiliary
  )

  return(posterior::as_draws(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_array.brma <- function(x, include_auxiliary = FALSE, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(
    x,
    include_auxiliary = include_auxiliary
  )

  return(posterior::as_draws_array(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_df.brma <- function(x, include_auxiliary = FALSE, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(
    x,
    include_auxiliary = include_auxiliary
  )

  return(posterior::as_draws_df(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_list.brma <- function(x, include_auxiliary = FALSE, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(
    x,
    include_auxiliary = include_auxiliary
  )

  return(posterior::as_draws_list(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_matrix.brma <- function(x, include_auxiliary = FALSE, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(
    x,
    include_auxiliary = include_auxiliary
  )

  return(posterior::as_draws_matrix(mcmc.list, ...))
}

#' @rdname as_draws.brma
#' @export
as_draws_rvars.brma <- function(x, include_auxiliary = FALSE, ...) {

  .check_posterior_package()

  mcmc.list <- .brma_to_mcmc.list(
    x,
    include_auxiliary = include_auxiliary
  )

  return(posterior::as_draws_rvars(mcmc.list, ...))
}
