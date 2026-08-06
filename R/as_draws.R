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
#' variables. Defaults to \code{FALSE}, returning the stable model-parameter
#' schema. Set to \code{TRUE} to expose the non-private coordinates recorded
#' by the fitted BayesTools parameter registry without adding RoBMA-derived
#' quantities.
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
#' Cholesky factors, and prior-parameterization variables are omitted. Public
#' covariance parameters, derived correlation matrices, allocation weights,
#' model indicators, and likelihood parameters remain available.
#' Scalar-structure random-effect correlation matrices are reconstructed from
#' compact rho draws on demand; those compact internal coordinates are never
#' returned. Conversion fails before allocating the dense matrices when
#' their combined draw-by-matrix size exceeds
#' `getOption("RoBMA.max_derived_correlation_cells", 2e7)`. Increase this option
#' explicitly when the dense public matrices are required. Use
#' Structural point-prior coordinates are returned as exact constant columns,
#' with chain timing taken from the fitted draw geometry. Private backend
#' anchors are never exposed. Use `include_auxiliary = TRUE` to skip RoBMA
#' filtering and derivation while retaining BayesTools' public/private
#' boundary. \code{brma_samples} objects have
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
    requires = c("parameter_registry", "draw_geometry")
  )
  mcmc_list <- BayesTools::JAGS_materialize_draws(x[["fit"]])
  if (include_auxiliary) {
    return(mcmc_list)
  }

  mcmc_list <- .brma_filter_auxiliary_variables(x, mcmc_list)
  return(.brma_append_derived_random_correlations(x, mcmc_list))
}


# Reconstruct public scalar-structure correlation matrices from compact rho draws.
.brma_append_derived_random_correlations <- function(x, mcmc_list) {

  formula_design <- attr(x[["fit"]], "formula_design", exact = TRUE)
  if (is.null(formula_design) || !is.list(formula_design)) {
    return(mcmc_list)
  }

  random_terms <- unlist(lapply(formula_design, function(design) {
    if (is.null(design[["random_effects"]])) {
      return(list())
    }
    design[["random_effects"]]
  }), recursive = FALSE)
  compact_names <- .brma_random_correlation_coordinate_names(random_terms)
  mcmc_list <- .brma_materialize_random_correlation_coordinates(
    x            = x,
    mcmc_list    = mcmc_list,
    random_terms = random_terms
  )
  .brma_check_derived_random_correlation_budget(
    random_terms = random_terms,
    mcmc_list    = mcmc_list
  )
  return(.brma_append_derived_random_correlation_terms(
    mcmc_list      = mcmc_list,
    random_terms   = random_terms,
    omit_variables = compact_names
  ))
}


# Return the compact coordinates used to derive public scalar correlations.
.brma_random_correlation_coordinate_names <- function(random_terms) {

  names <- unlist(lapply(random_terms, function(random_term) {
    spec <- .brma_derived_random_correlation_spec(random_term)
    if (is.null(spec)) {
      return(character())
    }
    correlation <- spec[["random_term"]][["correlation"]]
    if (!is.list(correlation)) {
      return(character())
    }
    c(correlation[["rho_name"]], correlation[["sample_name"]])
  }), use.names = FALSE)
  names <- names[is.character(names) & !is.na(names) & nzchar(names)]

  return(unique(names))
}


# Materialize only missing sampled rho coordinates needed for reconstruction.
.brma_materialize_random_correlation_coordinates <- function(x, mcmc_list,
                                                              random_terms) {

  variables <- colnames(as.matrix(mcmc_list[[1L]]))
  parameters <- unlist(lapply(random_terms, function(random_term) {
    spec <- .brma_derived_random_correlation_spec(random_term)
    if (is.null(spec) || spec[["n_columns"]] == 1L) {
      return(character())
    }
    correlation <- spec[["random_term"]][["correlation"]]
    if (!is.list(correlation) || !is.null(correlation[["sample_fixed"]])) {
      return(character())
    }
    candidates <- unique(c(
      correlation[["rho_name"]],
      correlation[["sample_name"]]
    ))
    candidates <- candidates[
      is.character(candidates) & !is.na(candidates) & nzchar(candidates)
    ]
    if (any(candidates %in% variables)) {
      return(character())
    }
    correlation[["rho_name"]]
  }), use.names = FALSE)
  parameters <- unique(parameters[
    is.character(parameters) & !is.na(parameters) & nzchar(parameters)
  ])
  if (length(parameters) == 0L) {
    return(mcmc_list)
  }

  internal <- BayesTools::JAGS_materialize_draws(
    x[["fit"]],
    parameters       = parameters,
    include_internal = TRUE
  )
  if (length(internal) != length(mcmc_list)) {
    stop(
      "Internal random-correlation draws disagree with the public draw geometry.",
      call. = FALSE
    )
  }
  chains <- lapply(seq_along(mcmc_list), function(chain_i) {
    public_chain   <- mcmc_list[[chain_i]]
    internal_chain <- internal[[chain_i]]
    if (nrow(public_chain) != nrow(internal_chain) ||
        !identical(coda::mcpar(public_chain), coda::mcpar(internal_chain))) {
      stop(
        "Internal random-correlation draws disagree with the public draw geometry.",
        call. = FALSE
      )
    }
    values <- cbind(as.matrix(public_chain), as.matrix(internal_chain))
    mcpar  <- coda::mcpar(public_chain)
    coda::mcmc(
      values,
      start = mcpar[[1L]],
      end   = mcpar[[2L]],
      thin  = mcpar[[3L]]
    )
  })

  return(coda::mcmc.list(chains))
}


# Guard the total dense output allocation before deriving any matrices.
.brma_check_derived_random_correlation_budget <- function(random_terms,
                                                          mcmc_list) {

  max_cells <- getOption("RoBMA.max_derived_correlation_cells", 2e7)
  if (!is.numeric(max_cells) || length(max_cells) != 1L ||
      is.na(max_cells) || max_cells < 1) {
    stop(
      "Option 'RoBMA.max_derived_correlation_cells' must be a positive numeric scalar.",
      call. = FALSE
    )
  }
  n_draws <- sum(vapply(mcmc_list, nrow, integer(1)))
  required_cells <- sum(vapply(random_terms, function(random_term) {
    spec <- .brma_derived_random_correlation_spec(random_term)
    if (is.null(spec)) {
      return(0)
    }
    as.numeric(spec[["n_columns"]])^2 * n_draws
  }, numeric(1)))
  if (required_cells > max_cells) {
    stop(
      "Default draw conversion would derive ",
      format(required_cells, scientific = FALSE, trim = TRUE),
      " dense random-correlation cells, exceeding option ",
      "'RoBMA.max_derived_correlation_cells' (",
      format(max_cells, scientific = FALSE, trim = TRUE),
      "). Increase the option or use 'include_auxiliary = TRUE' to return ",
      "raw backend draws without derivation.",
      call. = FALSE
    )
  }

  return(invisible(required_cells))
}


# Validate metadata needed to reconstruct one scalar correlation structure.
.brma_derived_random_correlation_spec <- function(random_term) {

  structure <- random_term[["structure"]]
  if (!is.character(structure) || length(structure) != 1L ||
      is.na(structure) ||
      !tolower(structure) %in% c("cs", "hcs", "ar1", "car", "har")) {
    return(NULL)
  }
  monitor <- random_term[["monitor"]]
  if (!is.null(monitor) && !isTRUE(monitor[["correlation"]])) {
    return(NULL)
  }

  parameter_stem <- random_term[["parameter_stem"]]
  n_columns      <- random_term[["n_columns"]]
  if (!is.character(parameter_stem) || length(parameter_stem) != 1L ||
      is.na(parameter_stem) || !nzchar(parameter_stem) ||
      !is.numeric(n_columns) || length(n_columns) != 1L ||
      is.na(n_columns) || !is.finite(n_columns) ||
      n_columns != floor(n_columns) || n_columns < 1 ||
      n_columns > .Machine$integer.max) {
    stop("Invalid scalar random-correlation metadata in fitted object.",
         call. = FALSE)
  }
  n_columns        <- as.integer(n_columns)
  correlation_base <- paste0(parameter_stem, "_xRE_CORx_R")
  correlation      <- random_term[["correlation"]]
  anchor_names    <- unique(c(
    random_term[["sd_parameter_names"]],
    if (is.list(correlation)) correlation[["rho_name"]] else NULL,
    if (is.list(correlation)) correlation[["sample_name"]] else NULL
  ))
  anchor_names <- anchor_names[
    is.character(anchor_names) & !is.na(anchor_names) & nzchar(anchor_names)
  ]

  return(list(
    random_term      = random_term,
    block_name       = random_term[["block_name"]],
    parameter_stem   = parameter_stem,
    n_columns        = n_columns,
    correlation_base = correlation_base,
    anchor_names     = anchor_names
  ))
}


# Generate canonical column names for one derived correlation matrix.
.brma_derived_random_correlation_names <- function(spec) {

  return(unlist(lapply(seq_len(spec[["n_columns"]]), function(column) {
    paste0(
      spec[["correlation_base"]], "[",
      seq_len(spec[["n_columns"]]), ",", column, "]"
    )
  }), use.names = FALSE))
}


# Locate the legacy block-local insertion point for derived matrix columns.
.brma_derived_random_correlation_anchors <- function(variables, specs) {

  anchors <- vapply(specs, function(spec) {
    indexes <- unlist(lapply(spec[["anchor_names"]], function(anchor) {
      which(variables == anchor | startsWith(variables, paste0(anchor, "[")))
    }), use.names = FALSE)
    if (length(indexes) == 0L) {
      return(NA_integer_)
    }
    max(indexes)
  }, integer(1))

  for (i in which(is.na(anchors))) {
    next_anchors     <- anchors[seq_along(anchors) > i & !is.na(anchors)]
    previous_anchors <- anchors[seq_along(anchors) < i & !is.na(anchors)]
    if (length(next_anchors) > 0L) {
      anchors[[i]] <- min(next_anchors) - 1L
    } else if (length(previous_anchors) > 0L) {
      anchors[[i]] <- max(previous_anchors)
    } else {
      anchors[[i]] <- length(variables)
    }
  }

  return(anchors)
}


# Interleave appended matrices at their legacy block-local monitor positions.
.brma_derived_random_correlation_order <- function(variables, specs,
                                                   derived) {

  anchors <- .brma_derived_random_correlation_anchors(variables, specs)
  after   <- vector("list", length(variables) + 1L)
  offset  <- length(variables)
  for (i in seq_along(specs)) {
    indexes <- offset + seq_len(ncol(derived[[i]]))
    after[[anchors[[i]] + 1L]] <- c(after[[anchors[[i]] + 1L]], indexes)
    offset <- offset + ncol(derived[[i]])
  }
  order <- vector("list", length(variables) + 1L)
  order[[1L]] <- after[[1L]]
  for (i in seq_along(variables)) {
    order[[i + 1L]] <- c(i, after[[i + 1L]])
  }

  return(unlist(order, use.names = FALSE))
}


# Reconstruct scalar correlation matrices with one copy per chain.
.brma_append_derived_random_correlation_terms <- function(
    mcmc_list, random_terms, omit_variables = character()) {

  specs <- lapply(random_terms, .brma_derived_random_correlation_spec)
  specs <- specs[!vapply(specs, is.null, logical(1))]
  if (length(specs) == 0L) {
    return(mcmc_list)
  }

  chains <- lapply(mcmc_list, function(chain) {
    values  <- as.matrix(chain)
    derived <- lapply(specs, function(spec) {
      matrix(
        BayesTools::random_effects_correlation_draws(
          random_term       = spec[["random_term"]],
          posterior_samples = values
        ),
        nrow     = nrow(values),
        dimnames = list(NULL, .brma_derived_random_correlation_names(spec))
      )
    })
    combined <- do.call(cbind, c(list(values), derived))
    combined <- combined[, .brma_derived_random_correlation_order(
      variables = colnames(values),
      specs     = specs,
      derived   = derived
    ), drop = FALSE]
    combined <- combined[
      , !colnames(combined) %in% omit_variables,
      drop = FALSE
    ]
    mcpar <- coda::mcpar(chain)
    coda::mcmc(
      combined,
      start = mcpar[[1L]],
      end   = mcpar[[2L]],
      thin  = mcpar[[3L]]
    )
  })

  return(coda::mcmc.list(chains))
}


# Reconstruct one CS/HCS/AR1/HAR/CAR correlation matrix per posterior draw.
.brma_append_derived_random_correlation <- function(mcmc_list, random_term) {

  return(.brma_append_derived_random_correlation_terms(
    mcmc_list    = mcmc_list,
    random_terms = list(random_term)
  ))
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
