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
#' @param component public component name used by data-frame coercion
#' @param probs default quantiles for credible intervals. Defaults to
#' \code{c(.025, .975)}
#' @param data optional data associated with the predictions (e.g., for
#' non-aggregated predictions)
#' @param effect_transform optional effect-size transformation metadata
#' @param prediction_mixture optional posterior predictive distribution of one
#' new value per column: a list of \code{location} and \code{scale} matrices
#' with the dimensions of \code{samples}, on the fitted (untransformed) scale.
#' Each posterior draw contributes a normal component (a point mass where the
#' scale is zero). When supplied, \code{summary()} appends the exact quantiles
#' of this mixture, transformed by \code{effect_transform}, as
#' prediction-interval columns.
#'
#' @return An object of class \code{brma_samples} which inherits from
#' \code{matrix}.
#'
#' @noRd
.new_brma_samples <- function(samples, n_chains, n_iter, title,
                              component = "result",
                              probs = c(.025, .975), data = NULL,
                              effect_transform = NULL,
                              prediction_mixture = NULL) {

  # ensure samples is a matrix with proper column names
  if (!is.matrix(samples)) {
    samples <- as.matrix(samples)
  }

  n_chains <- as.integer(n_chains)
  n_iter   <- as.integer(n_iter)

  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !nzchar(component)) {
    stop("'component' must be one non-empty character value.", call. = FALSE)
  }
  if (length(n_chains) != 1L || length(n_iter) != 1L ||
      is.na(n_chains) || is.na(n_iter) ||
      n_chains < 1L || n_iter < 1L ||
      n_chains * n_iter != nrow(samples)) {
    stop("Invalid brma_samples chain metadata: 'n_chains * n_iter' must equal the number of sample rows.",
         call. = FALSE)
  }

  if (!is.null(prediction_mixture)) {
    prediction_mixture <- .new_prediction_mixture(
      location = prediction_mixture[["location"]],
      scale    = prediction_mixture[["scale"]],
      samples  = samples
    )
  }

  # add attributes for MCMC structure
  attr(samples, "nchains") <- n_chains
  attr(samples, "niter")   <- n_iter

  # add attributes for display
  attr(samples, "title")     <- title
  attr(samples, "component") <- component
  attr(samples, "probs")     <- probs
  attr(samples, "data")      <- data

  if (!is.null(prediction_mixture)) {
    attr(samples, "prediction_mixture") <- prediction_mixture
  }

  if (!is.null(effect_transform)) {
    attr(samples, "effect_transform") <- effect_transform
    attr(samples, "footnotes")        <- effect_transform[["note"]]
  }

  # set class (inherits from matrix for backward compatibility)
  class(samples) <- c("brma_samples", "matrix", "array")

  return(samples)
}


# Add a result class to a list of brma_samples objects so list-level methods can
# preserve its component structure.
.new_brma_samples_list <- function(x) {

  if (!is.list(x) || is.matrix(x)) {
    stop("Internal error: expected a list of brma_samples results.",
         call. = FALSE)
  }

  class(x) <- unique(c("brma_samples_list", class(x), "list"))

  return(x)
}


# Derive valid chain metadata for brma_samples objects.
.brma_samples_chain_info <- function(fit = NULL, n_samples) {

  n_samples <- as.integer(n_samples)

  if (!is.null(fit) && !is.null(fit[["mcmc"]]) && length(fit[["mcmc"]]) > 0L) {
    chain_lengths <- vapply(fit[["mcmc"]], function(x) NROW(x), integer(1))
    if (sum(chain_lengths) == n_samples && length(unique(chain_lengths)) == 1L) {
      return(list(
        n_chains = length(chain_lengths),
        n_iter   = chain_lengths[[1]]
      ))
    }
  }

  return(list(
    n_chains = 1L,
    n_iter   = n_samples
  ))
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
#' @description Creates and returns a summary table of posterior samples
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
        samples    = .brma_samples_columns(object),
        parameters = colnames(object),
        probs      = probs,
        title      = attr(object, "title")
      ),
      dots
    )
  )

  prediction_mixture <- attr(object, "prediction_mixture", exact = TRUE)
  if (!is.null(prediction_mixture)) {
    # transforming the quantiles is exact: every effect transformation is
    # strictly increasing
    prediction_quantiles <- .transform_effect_matrix(
      .prediction_mixture_quantiles(prediction_mixture, probs),
      attr(object, "effect_transform", exact = TRUE)
    )
    for (i in seq_along(probs)) {
      summary_table[[paste("PI", probs[i])]] <- unname(prediction_quantiles[i, ])
    }
    attr(summary_table, "type") <- c(
      attr(summary_table, "type"),
      rep("estimate", length(probs))
    )
  }

  class(summary_table) <- c("summary.brma_samples", class(summary_table))
  attr(summary_table, "component") <- .brma_samples_component(object)

  return(summary_table)
}


# The columns of posterior samples as a named list of plain numeric vectors,
# the form BayesTools::ensemble_estimates_table() summarizes.
.brma_samples_columns <- function(samples) {

  stats::setNames(
    lapply(seq_len(ncol(samples)), function(i) as.numeric(samples[, i])),
    colnames(samples)
  )
}


# ---------------------------------------------------------------------------- #
# Posterior predictive mixtures
# ---------------------------------------------------------------------------- #

# The posterior predictive distribution of one new value is the equally
# weighted normal mixture over the posterior draws,
#   F(x) = mean_s pnorm(x, location_s, scale_s),
# with a point mass where scale_s = 0 (pnorm(sd = 0) is its right-continuous
# step). Its quantiles are computed from the draws directly rather than from
# one simulated predictive draw per posterior draw, which removes that
# simulation error and needs no random numbers.
.new_prediction_mixture <- function(location, scale, samples) {

  location <- as.matrix(location)
  scale    <- as.matrix(scale)
  if (!is.numeric(location) || !is.numeric(scale) ||
      !identical(dim(location), dim(samples)) ||
      !identical(dim(scale), dim(samples)) ||
      any(!is.finite(location)) || any(!is.finite(scale)) || any(scale < 0)) {
    stop(
      "The predictive mixture needs finite location and non-negative finite ",
      "scale matrices with the dimensions of the samples.",
      call. = FALSE
    )
  }
  dimnames(location) <- dimnames(scale) <- list(NULL, colnames(samples))

  return(list(location = location, scale = scale))
}

# Quantiles (rows) of each column's predictive mixture, on the fitted scale.
.prediction_mixture_quantiles <- function(prediction_mixture, probs) {

  location <- prediction_mixture[["location"]]
  scale    <- prediction_mixture[["scale"]]
  out <- vapply(seq_len(ncol(location)), function(j) {
    vapply(probs, .prediction_mixture_quantile, numeric(1),
           location = location[, j], scale = scale[, j])
  }, numeric(length(probs)))

  return(matrix(out, nrow = length(probs), ncol = ncol(location),
                dimnames = list(probs, colnames(location))))
}

# The generalized inverse inf{x : F(x) >= p}. Point masses alone are the
# empirical distribution of their locations. Otherwise the solve runs in
# standardized coordinates (the center and extent of the mixture), so that the
# bracket and tolerance follow the distribution's own scale. F(x) - p is
# assembled from an exact count and the small tails pnorm(-|z|) of the normal
# components, so its sign stays accurate where F rounds to p (low-density
# gaps); where every tail underflows, the log tails decide the sign. Brent's
# method narrows the bracket; where F is flat at p (between point masses) a
# root is not unique, and bisection on F(x) >= p finds the left end. A quantile
# in the jump of a point mass is that location exactly.
.prediction_mixture_quantile <- function(p, location, scale) {

  continuous <- scale > 0
  if (p <= 0) {
    return(if (any(continuous)) -Inf else min(location))
  }
  if (p >= 1) {
    return(if (any(continuous)) Inf else max(location))
  }

  n <- length(location)
  if (!any(continuous)) {
    # smallest k with k / n >= p (one step corrects the rounding of p * n)
    k <- min(max(ceiling(p * n), 1), n)
    if (k > 1 && (k - 1) / n >= p) {
      k <- k - 1
    } else if (k < n && k / n < p) {
      k <- k + 1
    }
    return(sort(location)[[k]])
  }

  center <- min(location) / 2 + max(location) / 2
  extent <- max(scale, abs(location - center))
  normal_location <- (location[continuous] - center) / extent
  normal_scale    <- scale[continuous] / extent
  atom_location   <- location[!continuous]
  atom_position   <- (atom_location - center) / extent
  distinct_atoms <- unique(atom_location)
  if (any(normal_scale == 0) ||
      anyDuplicated((distinct_atoms - center) / extent) > 0L) {
    stop("The predictive mixture spans too wide a range relative to its ",
         "smallest scale or the spacing of its point masses to compute its ",
         "quantiles in double precision.",
         call. = FALSE)
  }

  log_sum_exp <- function(x) {
    if (length(x) == 0L) {
      return(-Inf)
    }
    top <- max(x)
    top + log(sum(exp(x - top)))
  }

  # F(z) - p in standardized coordinates
  residual <- function(z) {
    standardized <- (z - normal_location) / normal_scale
    tail         <- stats::pnorm(-abs(standardized))
    above        <- standardized >= 0
    value <- (sum(above) + sum(atom_position <= z)) / n - p +
      (sum(tail[!above]) - sum(tail[above])) / n
    if (value == 0) {
      # only the sign is needed: compare the tails on the log scale
      log_tail <- stats::pnorm(-abs(standardized), log.p = TRUE)
      value    <- sign(log_sum_exp(log_tail[!above]) - log_sum_exp(log_tail[above])) *
        .Machine$double.xmin
      if (is.nan(value)) {
        value <- 0
      }
    }
    value
  }

  # F(lower) < p <= F(upper): every component lies below p at lower (above
  # 1 - p at upper), and lower is strictly below every point mass
  width <- 2 - stats::qnorm(min(p, 1 - p))
  lower <- min(c(normal_location - width * normal_scale, atom_position)) - 1e-8
  upper <- max(c(normal_location + width * normal_scale, atom_position))
  # near machine precision of the standardized coordinate, so that nearly
  # equal probabilities (0.025 and (1 - 0.95) / 2) give nearly equal quantiles
  tol   <- 1e-14
  root  <- stats::uniroot(residual, lower = lower, upper = upper, tol = tol)[["root"]]
  if (residual(root) >= 0) {
    upper <- root
    if (residual(root - tol) < 0) {
      lower <- root - tol
    }
  } else {
    lower <- root
    if (residual(root + tol) >= 0) {
      upper <- root + tol
    }
  }
  while (upper - lower > tol) {
    middle <- lower + (upper - lower) / 2
    if (middle <= lower || middle >= upper) {
      break
    }
    if (residual(middle) >= 0) {
      upper <- middle
    } else {
      lower <- middle
    }
  }

  in_bracket <- atom_position > lower & atom_position <= upper
  atoms      <- unique(atom_location[in_bracket])
  atoms      <- atoms[vapply((atoms - center) / extent, residual, numeric(1)) >= 0]
  if (length(atoms) > 0L) {
    return(min(atoms))
  }

  quantile <- center + extent * upper
  if (!is.finite(quantile)) {
    # extent * upper can overflow while the sum is representable
    quantile <- extent * (center / extent + upper)
  }
  # (a finite mixture can have a quantile beyond the double range, e.g. scales
  # near .Machine$double.xmax; report it instead of returning Inf)
  if (!is.finite(quantile)) {
    stop("The ", p, " prediction quantile is outside the double-precision range.",
         call. = FALSE)
  }

  return(quantile)
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
# Data-frame coercion
# ---------------------------------------------------------------------------- #

#' @title Convert brma_samples Results to Data Frames
#'
#' @description \code{as.data.frame()} and \code{data.frame()} of a
#' \code{brma_samples} result return its posterior draws (one row per draw, one
#' column per parameter), as for the underlying matrix. The posterior summary
#' table displayed by \code{print()} is \code{as.data.frame(summary(x))}: a
#' long data frame with component and parameter identifiers (multi-component
#' paths use \code{/} as a separator), which is also the result for a
#' multi-component \code{brma_samples_list}. Printed quantile labels
#' such as \code{0.025} and \code{PI 0.025} are returned as the syntactic column
#' names \code{CI_0.025} and \code{PI_0.025}, respectively.
#'
#' @param x a \code{brma_samples}, \code{summary.brma_samples}, or
#' multi-component \code{brma_samples_list} result
#' @param row.names \code{NULL} or a character vector giving the row names for
#' the resulting data frame. Custom row names are unsupported when
#' \code{format = "list"}.
#' @param optional logical; passed to the final data-frame coercion method.
#' @param format for multi-component results, whether to return a single
#' \code{"long"} data frame or a named, possibly nested \code{"list"} with a
#' plain data frame at each leaf. Defaults to \code{"long"}.
#' @param stringsAsFactors logical; passed to the output data frame.
#' @param ... for a \code{brma_samples} result, additional arguments passed to
#' the matrix coercion; for a \code{brma_samples_list}, additional arguments
#' passed to \code{summary.brma_samples()}
#'
#' @return For an individual \code{brma_samples} object, a \code{data.frame}
#' of the posterior draws. For a \code{summary.brma_samples} object, a plain
#' long \code{data.frame} containing the displayed posterior summary statistics
#' and leading \code{component} and \code{parameter} columns.
#' For a \code{brma_samples_list}, a long data frame by default or, with
#' \code{format = "list"}, a named, possibly nested list of such data frames.
#'
#' @export
as.data.frame.brma_samples <- function(x, row.names = NULL, optional = FALSE, ...) {

  return(as.data.frame(
    as.matrix(x),
    row.names = row.names,
    optional  = optional,
    ...
  ))
}


# The posterior summary table of a brma_samples result as a plain long data
# frame (what as.data.frame() of its summary returns).
.brma_samples_summary_data_frame <- function(x, row.names = NULL,
                                             optional = FALSE,
                                             stringsAsFactors = FALSE, ...) {

  summary_table <- summary.brma_samples(x, ...)
  output <- .output_table_as_long_data_frame(
    table            = summary_table,
    component        = .brma_samples_component(x),
    row.names        = row.names,
    optional         = optional,
    stringsAsFactors = stringsAsFactors
  )

  return(output)
}


#' @rdname as.data.frame.brma_samples
#' @export
as.data.frame.summary.brma_samples <- function(
    x, row.names = NULL, optional = FALSE, stringsAsFactors = FALSE, ...) {

  output <- .output_table_as_long_data_frame(
    table            = x,
    component        = .brma_samples_component(x),
    row.names        = row.names,
    optional         = optional,
    stringsAsFactors = stringsAsFactors
  )

  return(output)
}


#' @rdname as.data.frame.brma_samples
#' @export
as.data.frame.brma_samples_list <- function(
    x, row.names = NULL, optional = FALSE, format = c("long", "list"),
    stringsAsFactors = FALSE, ...) {

  format <- match.arg(format)
  if (identical(format, "list")) {
    if (!is.null(row.names)) {
      stop(
        "'row.names' is unsupported when format = 'list'.",
        call. = FALSE
      )
    }
    return(.brma_samples_list_as_data_frames(
      x,
      stringsAsFactors = stringsAsFactors,
      ...
    ))
  }

  output <- .brma_samples_list_as_long_data_frame(
    x,
    stringsAsFactors = stringsAsFactors,
    ...
  )
  output <- as.data.frame(
    output,
    row.names = row.names,
    optional  = optional
  )

  return(output)
}


.output_data_frame_names <- function(x) {

  column_names <- names(x)
  is_prediction_interval <- startsWith(column_names, "PI ")
  quantile_labels <- ifelse(
    is_prediction_interval,
    substring(column_names, 4L),
    column_names
  )
  is_credible_interval <- !is_prediction_interval &
    !is.na(suppressWarnings(as.numeric(quantile_labels)))

  interval_columns <- is_credible_interval | is_prediction_interval
  interval_names <- paste0(
    ifelse(is_prediction_interval, "PI_", "CI_"),
    quantile_labels
  )
  column_names[interval_columns] <- make.names(
    interval_names[interval_columns],
    unique = TRUE
  )
  names(x) <- column_names

  return(x)
}


.output_plain_data_frame <- function(x) {

  if (is.data.frame(x)) {
    class(x) <- "data.frame"
    for(column in names(x)){
      if(inherits(x[[column]], "BayesTools_BF")){
        x[[column]] <- BayesTools::BF_carrier_drop(x[[column]])
      }
    }
    return(x)
  }

  return(as.data.frame(x, optional = TRUE))
}


.output_table_as_long_data_frame <- function(
    table, component, parameter = NULL, row.names = NULL, optional = FALSE,
    stringsAsFactors = FALSE) {

  table <- .output_plain_data_frame(table)
  table <- .output_data_frame_names(table)
  if (is.null(parameter)) {
    parameter <- rownames(table)
  }
  if (is.null(parameter)) {
    parameter <- as.character(seq_len(nrow(table)))
  }
  if (length(parameter) != nrow(table)) {
    stop("Internal error: output parameter labels do not match table rows.",
         call. = FALSE)
  }
  if (length(component) == 1L) {
    component <- rep(component, nrow(table))
  }
  if (!is.character(component) || length(component) != nrow(table) ||
      anyNA(component) || any(!nzchar(component))) {
    stop("Internal error: invalid output component labels.", call. = FALSE)
  }

  output <- data.frame(
    component = component,
    parameter = as.character(parameter),
    table,
    row.names        = NULL,
    check.names      = FALSE,
    stringsAsFactors = stringsAsFactors
  )
  output <- as.data.frame(
    output,
    row.names = row.names,
    optional  = optional
  )

  return(output)
}


.output_bind_long_data_frames <- function(
    tables, row.names = NULL, optional = FALSE) {

  tables <- Filter(Negate(is.null), tables)
  if (length(tables) == 0L) {
    output <- data.frame(
      component = character(),
      parameter = character()
    )
  } else {
    column_names <- unique(unlist(lapply(tables, names), use.names = FALSE))
    aligned_tables <- lapply(tables, function(table) {
      missing_names <- setdiff(column_names, names(table))
      for (name in missing_names) {
        table[[name]] <- rep(NA, nrow(table))
      }
      table[column_names]
    })
    output <- do.call(rbind, aligned_tables)
    row_counts <- vapply(tables, nrow, integer(1L))
    for (column in column_names) {
      output[[column]] <- .output_bind_bf_attributes(
        output[[column]], lapply(tables, `[[`, column), row_counts
      )
    }
    rownames(output) <- NULL
  }

  output <- as.data.frame(
    output,
    row.names = row.names,
    optional  = optional
  )

  return(output)
}


.output_bind_bf_attributes <- function(output, sources, row_counts) {

  # Plain exports retain their visible BF values and existing bound metadata.
  # A first-row carrier left behind by rbind cannot describe the bound column.
  for(attribute in c("canonical_log_BF", "producer_values", "carrier_logBF",
                     "carrier_BF01", "numerical_diagnostics")) attr(output, attribute) <- NULL
  has_bounds <- vapply(sources, function(source) {
    !is.null(attr(source, "bound_operator", exact = TRUE))
  }, logical(1L))
  is_bf <- vapply(sources, inherits, logical(1L), what = "BayesTools_BF")
  if (!any(has_bounds | is_bf)) {
    return(output)
  }

  # rbind() retains only the first column's attributes, even when a later
  # table supplies its first Bayes factor. Keep its scale and presentation
  # metadata, and bind row-level bounds with the same missing-cell padding.
  source <- sources[[which(has_bounds | is_bf)[[1L]]]]
  for (attribute in c("name", "logBF", "BF01")) {
    attr(output, attribute) <- attr(source, attribute, exact = TRUE)
  }
  bounds <- lapply(seq_along(sources), function(i) {
    .output_bf_bound_operators(
      attr(sources[[i]], "bound_operator", exact = TRUE), row_counts[[i]]
    )
  })
  attr(output, "bound_operator") <- unlist(bounds, use.names = FALSE)
  if (any(!is.na(attr(output, "bound_operator", exact = TRUE)))) {
    class(output) <- unique(c("BayesTools_BF", class(output)))
  }
  return(output)
}


.output_bf_bound_operators <- function(bounds, n_rows) {

  if (is.null(bounds)) {
    return(rep(NA_character_, n_rows))
  }
  # BayesTools accepts a scalar as one operator shared by every BF. Other
  # mismatched lengths can be stale row metadata and must not be recycled.
  if (length(bounds) == 1L) {
    return(rep(bounds, n_rows))
  }
  if (length(bounds) != n_rows) {
    stop("Internal error: Bayes factor bounds do not match table rows.",
         call. = FALSE)
  }
  return(bounds)
}


.brma_samples_component <- function(x) {

  component <- attr(x, "component", exact = TRUE)
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !nzchar(component)) {
    return("result")
  }

  return(component)
}


.brma_samples_list_leaf_component <- function(x, path) {

  component <- .brma_samples_component(x)
  if (!identical(component, "result")) {
    return(component)
  }

  return(paste(path, collapse = "/"))
}


.brma_samples_list_as_long_data_frame <- function(x, ...) {

  tables <- .brma_samples_list_flatten_data_frames(x, ...)
  if (length(tables) == 0L) {
    return(data.frame(
      component = character(),
      parameter = character()
    ))
  }

  return(.output_bind_long_data_frames(tables))
}


.brma_samples_list_flatten_data_frames <- function(x, ..., path = character()) {

  output          <- list()
  component_names <- names(x)
  if (is.null(component_names)) {
    component_names <- rep("", length(x))
  }

  for (i in seq_along(x)) {
    component_name <- component_names[[i]]
    if (!nzchar(component_name)) {
      component_name <- as.character(i)
    }
    component_path <- c(path, component_name)
    component      <- x[[i]]

    if (inherits(component, "brma_samples")) {
      table <- .brma_samples_summary_data_frame(component, ...)
      table[["component"]] <- rep(
        .brma_samples_list_leaf_component(component, component_path),
        nrow(table)
      )
      table <- table[c(
        "component", "parameter",
        setdiff(names(table), c("component", "parameter"))
      )]
      rownames(table) <- NULL
      output[[length(output) + 1L]] <- table
      next
    }

    if (is.list(component) && !is.matrix(component)) {
      nested <- .brma_samples_list_flatten_data_frames(
        component,
        ...,
        path = component_path
      )
      output <- c(output, nested)
      next
    }

    stop("Internal error: invalid component in a brma_samples result list.",
         call. = FALSE)
  }

  return(output)
}


.brma_samples_list_as_data_frames <- function(
    x, ..., path = character()) {

  component_names <- names(x)
  if (is.null(component_names)) {
    component_names <- rep("", length(x))
  }
  output <- lapply(seq_along(x), function(i) {
    component_name <- component_names[[i]]
    if (!nzchar(component_name)) {
      component_name <- as.character(i)
    }
    component_path <- c(path, component_name)
    component      <- x[[i]]
    if (inherits(component, "brma_samples")) {
      table <- .brma_samples_summary_data_frame(component, ...)
      table[["component"]] <- rep(
        .brma_samples_list_leaf_component(component, component_path),
        nrow(table)
      )
      table <- table[c(
        "component", "parameter",
        setdiff(names(table), c("component", "parameter"))
      )]
      return(table)
    }
    if (is.list(component) && !is.matrix(component)) {
      return(.brma_samples_list_as_data_frames(
        component,
        ...,
        path = component_path
      ))
    }
    stop("Internal error: invalid component in a brma_samples result list.",
         call. = FALSE)
  })
  names(output) <- component_names

  return(output)
}


# Keep list printing unchanged while suppressing the implementation class.
#' @export
#' @noRd
print.brma_samples_list <- function(x, ...) {

  print(unclass(x), ...)

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
