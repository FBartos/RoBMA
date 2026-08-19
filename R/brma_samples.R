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
#' @param prediction_samples optional posterior predictive samples with the
#' same dimensions as \code{samples}. When supplied, their posterior quantiles
#' are appended as prediction-interval columns by \code{summary()}.
#'
#' @return An object of class \code{brma_samples} which inherits from
#' \code{matrix}.
#'
#' @noRd
.new_brma_samples <- function(samples, n_chains, n_iter, title,
                              component = "result",
                              probs = c(.025, .975), data = NULL,
                              effect_transform = NULL,
                              prediction_samples = NULL) {

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

  if (!is.null(prediction_samples)) {
    if (!is.matrix(prediction_samples)) {
      prediction_samples <- as.matrix(prediction_samples)
    }
    if (!is.numeric(prediction_samples) ||
        !identical(dim(prediction_samples), dim(samples)) ||
        any(!is.finite(prediction_samples))) {
      stop(
        "'prediction_samples' must be a finite numeric matrix with the same dimensions as 'samples'.",
        call. = FALSE
      )
    }
    colnames(prediction_samples) <- colnames(samples)
  }

  # add attributes for MCMC structure
  attr(samples, "nchains") <- n_chains
  attr(samples, "niter")   <- n_iter

  # add attributes for display
  attr(samples, "title")     <- title
  attr(samples, "component") <- component
  attr(samples, "probs")     <- probs
  attr(samples, "data")      <- data

  if (!is.null(prediction_samples)) {
    attr(samples, "prediction_samples") <- prediction_samples
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
        samples    = asplit(object, 2),
        parameters = colnames(object),
        probs      = probs,
        title      = attr(object, "title")
      ),
      dots
    )
  )

  prediction_samples <- attr(object, "prediction_samples", exact = TRUE)
  if (!is.null(prediction_samples)) {
    prediction_table <- BayesTools::ensemble_estimates_table(
      samples    = asplit(prediction_samples, 2),
      parameters = colnames(object),
      probs      = probs,
      title      = attr(object, "title")
    )
    quantile_columns <- utils::tail(
      colnames(prediction_table),
      length(probs)
    )
    for (column in quantile_columns) {
      summary_table[[paste("PI", column)]] <- prediction_table[[column]]
    }
    attr(summary_table, "type") <- c(
      attr(summary_table, "type"),
      rep("estimate", length(quantile_columns))
    )
  }

  class(summary_table) <- c("summary.brma_samples", class(summary_table))
  attr(summary_table, "component") <- .brma_samples_component(object)
 
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
# Data-frame coercion
# ---------------------------------------------------------------------------- #

#' @title Convert brma_samples Results to Data Frames
#'
#' @description Converts a \code{brma_samples} result to the posterior summary
#' table displayed by \code{print()}. Results are returned in long form with
#' component and parameter identifiers; multi-component paths use \code{/} as
#' a separator. Printed quantile labels
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
#' @param ... additional arguments passed to \code{summary.brma_samples()}
#'
#' @return For an individual \code{brma_samples} object, a plain long
#' \code{data.frame} containing the displayed posterior summary statistics and
#' leading \code{component} and \code{parameter} columns.
#' For a \code{brma_samples_list}, a long data frame by default or, with
#' \code{format = "list"}, a named, possibly nested list of such data frames.
#'
#' @export
as.data.frame.brma_samples <- function(x, row.names = NULL, optional = FALSE,
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
    tables <- lapply(tables, function(table) {
      missing_names <- setdiff(column_names, names(table))
      for (name in missing_names) {
        table[[name]] <- rep(NA, nrow(table))
      }
      table[column_names]
    })
    output <- do.call(rbind, tables)
    rownames(output) <- NULL
  }

  output <- as.data.frame(
    output,
    row.names = row.names,
    optional  = optional
  )

  return(output)
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

  return(do.call(rbind, tables))
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
      table <- as.data.frame(component, ...)
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
      table <- as.data.frame(component, ...)
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
