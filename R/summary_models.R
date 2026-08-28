#' @title Summarize Model-Averaged Component Weights
#'
#' @description Creates marginal or individual model-weight summaries for
#' RoBMA-class product-space objects.
#'
#' @param object a fitted RoBMA-class product-space object, including
#' \code{RoBMA}, \code{BMA}/\code{BMA.norm}, \code{BMA.glmm}, and
#' \code{BMA.mv}.
#' @param type whether to summarize marginal component prior distributions
#'   (\code{"marginal"}) or individual model combinations (\code{"individual"}).
#' @param include_mcmc_diagnostics whether to include Bayes factor MCMC
#'   diagnostics in the output. Defaults to \code{TRUE}.
#' @param ... additional arguments
#'
#' @details Only mixture-prior components are summarized; non-mixture
#' components are omitted.
#'
#' @return A list of class \code{summary_models.RoBMA} with elements
#' \code{name} and \code{type}. For \code{type = "marginal"}, element
#' \code{marginal} contains component tables with columns such as
#' \code{prior_prob}, \code{post_prob}, and \code{inclusion_BF}. For
#' \code{type = "individual"}, element \code{individual} contains individual
#' model combinations and posterior probabilities.
#'
#' @examples \dontrun{
#' if (requireNamespace("metadat", quietly = TRUE)) {
#'   data(dat.lehmann2018, package = "metadat")
#'   fit <- RoBMA(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")
#'
#'   summary_models(fit)
#'   summary_models(fit, type = "individual")
#' }
#' }
#'
#' @export
summary_models <- function(object, ...) {

  UseMethod("summary_models")
}

#' @export
summary_models.default <- function(object, ...) {

  stop("'summary_models' is available only for RoBMA objects.", call. = FALSE)
}

#' @rdname summary_models
#' @export
summary_models.RoBMA <- function(object, type = "marginal",
                                 include_mcmc_diagnostics = TRUE, ...) {

  BayesTools::check_char(type, "type")
  BayesTools::check_bool(include_mcmc_diagnostics, "include_mcmc_diagnostics")
  type <- match.arg(type, c("marginal", "individual"))

  if (is.null(object[["fit"]]) || length(object[["fit"]]) == 0) {
    stop("'summary_models' requires a fitted RoBMA object.", call. = FALSE)
  }

  components             <- .summary_models_components(object)
  posterior_samples      <- .get_posterior_samples(object[["fit"]])
  posterior_samples_list <- if (include_mcmc_diagnostics) {
    .summary_models_posterior_samples_list(object[["fit"]])
  } else {
    NULL
  }

  if (length(components) == 0L) {
    stop("No mixture-prior model components found.", call. = FALSE)
  }

  output <- list(
    name = .summary.brma_model_names(object),
    type = type
  )

  if (type == "marginal") {
    output[["marginal"]] <- .summary_models_marginal_tables(
      components             = components,
      posterior_samples      = posterior_samples,
      posterior_samples_list = posterior_samples_list
    )
  } else {
    output[["individual"]] <- .summary_models_individual_table(
      components             = components,
      posterior_samples      = posterior_samples,
      posterior_samples_list = posterior_samples_list
    )
  }

  class(output) <- "summary_models.RoBMA"

  return(output)
}

#' @rdname summary_models
#' @param x a `summary_models.RoBMA` object.
#' @export
print.summary_models.RoBMA <- function(x, ...) {

  cat("\n")

  if (x[["type"]] == "marginal") {
    for (i in seq_along(x[["marginal"]])) {
      if (i > 1L) {
        cat("\n")
      }
      table <- x[["marginal"]][[i]]
      print(table)
    }
  } else {
    print(x[["individual"]])
  }

  cat("\n")

  return(invisible(x))
}


#' @title Convert Model-Weight Summaries to a Data Frame
#'
#' @description Converts every table displayed by a
#' \code{summary_models.RoBMA} object to one component-aware long data frame.
#'
#' @param x a \code{summary_models.RoBMA} object.
#' @param row.names \code{NULL} or a character vector giving the row names.
#' @param optional logical; passed to the final data-frame coercion.
#' @param stringsAsFactors accepted for compatibility with \code{data.frame()}.
#' @param ... unused additional arguments.
#'
#' @return A plain \code{data.frame} with leading \code{component} and
#' \code{parameter} columns.
#'
#' @export
as.data.frame.summary_models.RoBMA <- function(
    x, row.names = NULL, optional = FALSE, stringsAsFactors = FALSE, ...) {

  if (identical(x[["type"]], "marginal")) {
    tables <- lapply(names(x[["marginal"]]), function(component) {
      .output_table_as_long_data_frame(
        table            = x[["marginal"]][[component]],
        component        = component,
        stringsAsFactors = stringsAsFactors
      )
    })
  } else {
    table <- x[["individual"]]
    tables <- list(.output_table_as_long_data_frame(
      table            = table,
      component        = "individual models",
      parameter        = paste("model", seq_len(nrow(table))),
      stringsAsFactors = stringsAsFactors
    ))
  }

  output <- .output_bind_long_data_frames(
    tables    = tables,
    row.names = row.names,
    optional  = optional
  )

  return(output)
}


# Extract top-level model components used by summary_models().
.summary_models_components <- function(object) {

  prior_list <- attr(object[["fit"]], "prior_list")
  components <- list()

  effect_parameter <- if ("mu" %in% names(prior_list)) "mu" else "mu_intercept"
  components <- .summary_models_add_component(
    components = components,
    prior_list = prior_list,
    component  = "Effect",
    parameter  = effect_parameter
  )

  location_parameters <- grep("^mu_", names(prior_list), value = TRUE)
  location_parameters <- location_parameters[location_parameters != "mu_intercept"]
  location_parameters <- setdiff(
    location_parameters,
    .random_slab_prior_parameters(prior_list)
  )
  for (parameter in location_parameters) {
    component <- paste0(
      "Location: ",
      .summary_parameter_label(sub("^mu_", "", parameter))
    )
    components <- .summary_models_add_component(
      components = components,
      prior_list = prior_list,
      component  = component,
      parameter  = parameter
    )
  }

  heterogeneity_parameter <- if ("tau" %in% names(prior_list)) "tau" else "log_tau_intercept"
  components <- .summary_models_add_component(
    components = components,
    prior_list = prior_list,
    component  = "Heterogeneity",
    parameter  = heterogeneity_parameter
  )
  components <- .summary_models_add_random_slab_components(
    components = components,
    prior_list = prior_list
  )

  scale_parameters <- grep("^log_tau_", names(prior_list), value = TRUE)
  scale_parameters <- scale_parameters[scale_parameters != "log_tau_intercept"]
  scale_parameters <- setdiff(
    scale_parameters,
    .random_slab_prior_parameters(prior_list)
  )
  for (parameter in scale_parameters) {
    component <- paste0(
      "Scale: ",
      .summary_parameter_label(sub("^log_tau_", "", parameter))
    )
    components <- .summary_models_add_component(
      components = components,
      prior_list = prior_list,
      component  = component,
      parameter  = parameter
    )
  }

  components <- .summary_models_add_component(
    components = components,
    prior_list = prior_list,
    component  = "Publication Bias",
    parameter  = "bias"
  )
  components <- .summary_models_add_random_components(
    components = components,
    object     = object,
    prior_list = prior_list
  )

  return(components)
}


.summary_models_add_random_slab_components <- function(components, prior_list) {

  parameters <- .random_slab_prior_parameters(prior_list)
  for (parameter in parameters) {
    allocation <- attr(
      prior_list[[parameter]], "random_allocation", exact = TRUE
    )
    component <- if (length(parameters) == 1L) {
      "Heterogeneity Slab"
    } else {
      paste0("Random Slab: ", allocation)
    }
    components <- .summary_models_add_component(
      components = components,
      prior_list = prior_list,
      component  = component,
      parameter  = parameter
    )
  }

  components
}

.summary_models_add_component <- function(components, prior_list, component, parameter) {

  if (!parameter %in% names(prior_list)) {
    return(components)
  }
  if (!BayesTools::is.prior.mixture(prior_list[[parameter]])) {
    return(components)
  }

  prior <- prior_list[[parameter]]
  components[[component]] <- list(
    component        = component,
    parameter        = parameter,
    prior            = prior,
    indicator        = NULL,
    indicator_offset = 0L,
    names            = .summary_models_prior_names(prior),
    hypothesis       = .summary_models_hypothesis(prior),
    prior_probs      = .summary_models_prior_probs(prior)
  )

  return(components)
}

.summary_models_add_random_components <- function(components, object,
                                                  prior_list) {

  quantities <- BayesTools::parameter_catalog(object[["fit"]])[["quantities"]]
  random_inclusion <- which(
    quantities[["role"]] == "random_inclusion" &
      !quantities[["internal"]] &
      quantities[["status"]] != "unavailable"
  )
  for (i in random_inclusion) {
    key       <- quantities[["extraction_key"]][[i]]
    indicator <- key[["source_parameter"]]
    gate_prior <- .random_allocation_inclusion_prior(
      prior_list     = prior_list,
      indicator_name = indicator,
      required       = TRUE
    )
    prior_probability <- mean(gate_prior)
    if (prior_probability %in% c(0, 1)) {
      next
    }
    component_name <- paste(quantities[["arguments"]][[i]], collapse = ", ")
    component      <- paste0("Random: ", component_name)
    components[[component]] <- list(
      component        = component,
      parameter        = sub("_indicator$", "", indicator),
      prior            = vector("list", 2L),
      indicator        = indicator,
      indicator_offset = 1L,
      names            = c("Excluded", "Included"),
      hypothesis       = c("Null", "Alternative"),
      prior_probs      = c(1 - prior_probability, prior_probability)
    )
  }

  components
}

.summary_models_marginal_tables <- function(components, posterior_samples,
                                            posterior_samples_list = NULL) {

  tables <- list()

  for (component in names(components)) {

    info       <- components[[component]]
    indicators <- .summary_models_indicators(
      posterior_samples = posterior_samples,
      parameter         = info[["parameter"]],
      prior             = info[["prior"]],
      column            = info[["indicator"]],
      offset            = info[["indicator_offset"]]
    )
    post_probs <- vapply(seq_along(info[["names"]]), function(i) {

      mean(indicators == i)
    }, numeric(1))
    inclusion_BF <- .summary_models_inclusion_BF(
      prior_probs = info[["prior_probs"]],
      post_probs  = post_probs
    )
    BF_error_percent <- .summary_models_BF_error_percent(
      indicator_matrices = .summary_models_marginal_indicator_matrices(
        posterior_samples_list = posterior_samples_list,
        parameter              = info[["parameter"]],
        prior                  = info[["prior"]],
        column                 = info[["indicator"]],
        offset                 = info[["indicator_offset"]]
      ),
      prior_probs        = info[["prior_probs"]],
      post_probs         = post_probs,
      inclusion_BF       = inclusion_BF
    )

    table <- data.frame(
      Hypothesis   = info[["hypothesis"]],
      prior_prob   = info[["prior_probs"]],
      post_prob    = post_probs,
      inclusion_BF = inclusion_BF,
      row.names    = info[["names"]],
      check.names  = FALSE
    )

    table[["inclusion_BF"]] <- BayesTools::format_BF(
      table[["inclusion_BF"]],
      inclusion = TRUE
    )
    table <- .summary_models_add_BF_diagnostics(
      table            = table,
      BF_error_percent = BF_error_percent
    )

    class(table)             <- c("BayesTools_table", "RoBMA_summary_models_marginal", class(table))
    attr(table, "type")      <- c(
      "string_left", "prior_prob", "post_prob", "inclusion_BF",
      .summary_models_BF_diagnostic_types(BF_error_percent)
    )
    attr(table, "rownames")  <- TRUE
    attr(table, "title")     <- component
    attr(table, "footnotes") <- NULL
    attr(table, "warnings")  <- NULL

    tables[[component]] <- table
  }

  return(tables)
}

.summary_models_individual_table <- function(components, posterior_samples,
                                             posterior_samples_list = NULL) {

  component_names <- names(components)
  component_grid  <- expand.grid(
    lapply(components, function(info) seq_along(info[["names"]])),
    KEEP.OUT.ATTRS   = FALSE,
    stringsAsFactors = FALSE
  )

  prior_columns <- data.frame(
    lapply(seq_along(components), function(i) {

      info <- components[[i]]
      info[["names"]][component_grid[[i]]]
    }),
    check.names = FALSE
  )
  names(prior_columns) <- component_names

  prior_probs <- vapply(seq_len(nrow(component_grid)), function(i) {

    prod(vapply(seq_along(components), function(j) {

      components[[j]][["prior_probs"]][component_grid[i, j]]
    }, numeric(1)))
  }, numeric(1))

  indicators <- lapply(components, function(info) {

    .summary_models_indicators(
      posterior_samples = posterior_samples,
      parameter         = info[["parameter"]],
      prior             = info[["prior"]],
      column            = info[["indicator"]],
      offset            = info[["indicator_offset"]]
    )
  })

  post_probs <- vapply(seq_len(nrow(component_grid)), function(i) {

    selected <- rep(TRUE, nrow(posterior_samples))
    for (j in seq_along(components)) {
      selected <- selected & indicators[[j]] == component_grid[i, j]
    }
    mean(selected)
  }, numeric(1))

  inclusion_BF <- .summary_models_inclusion_BF(
    prior_probs = prior_probs,
    post_probs  = post_probs
  )
  BF_error_percent <- .summary_models_BF_error_percent(
    indicator_matrices = .summary_models_individual_indicator_matrices(
      components             = components,
      component_grid         = component_grid,
      posterior_samples_list = posterior_samples_list
    ),
    prior_probs        = prior_probs,
    post_probs         = post_probs,
    inclusion_BF       = inclusion_BF
  )

  table <- cbind(
    prior_columns,
    data.frame(
      prior_prob   = prior_probs,
      post_prob    = post_probs,
      inclusion_BF = BayesTools::format_BF(inclusion_BF, inclusion = TRUE),
      check.names  = FALSE
    )
  )
  table <- .summary_models_add_BF_diagnostics(
    table            = table,
    BF_error_percent = BF_error_percent
  )

  class(table)             <- c("BayesTools_table", "RoBMA_summary_models_individual", class(table))
  attr(table, "type")      <- c(
    rep("string_left", length(components)),
    "prior_prob", "post_prob", "inclusion_BF",
    .summary_models_BF_diagnostic_types(BF_error_percent)
  )
  attr(table, "rownames")  <- FALSE
  attr(table, "title")     <- "Individual Models"
  attr(table, "footnotes") <- NULL
  attr(table, "warnings")  <- NULL

  return(table)
}

.summary_models_prior_names <- function(prior) {

  prior_names <- names(prior)
  if (is.null(prior_names)) {
    prior_names <- rep("", length(prior))
  }

  missing_names <- !nzchar(prior_names)
  if (any(missing_names)) {
    prior_names[missing_names] <- vapply(prior[missing_names], function(p) {

      prior_name <- paste(capture.output(print(p)), collapse = " ")
      prior_name <- sub(" ~.*$", "", prior_name)
      trimws(prior_name)
    }, character(1))
  }

  return(make.unique(prior_names))
}

.summary_models_hypothesis <- function(prior) {

  components <- attr(prior, "components")
  hypothesis <- ifelse(components == "null", "Null", "Alternative")

  return(hypothesis)
}

.summary_models_prior_probs <- function(prior) {

  prior_probs <- attr(prior, "prior_weights")
  prior_probs <- prior_probs / sum(prior_probs)

  return(prior_probs)
}

.summary_models_indicators <- function(posterior_samples, parameter, prior,
                                       column = NULL, offset = 0L) {

  indicators <- .extract_posterior_indicator(
    posterior_samples = posterior_samples,
    parameter         = parameter,
    prior             = NULL,
    column            = column
  )
  indicators <- indicators + offset
  if (any(!indicators %in% seq_len(length(prior)))) {
    indicator_name <- if (is.null(column)) {
      paste0(parameter, "_indicator")
    } else {
      column
    }
    stop(
      "Invalid posterior model indicator range: '", indicator_name, "'.",
      call. = FALSE
    )
  }

  return(indicators)
}

.summary_models_posterior_samples_list <- function(fit) {

  posterior_samples_list <- suppressWarnings(coda::as.mcmc.list(fit))
  posterior_samples_list <- lapply(posterior_samples_list, as.matrix)

  return(posterior_samples_list)
}

.summary_models_indicators_list <- function(posterior_samples_list, parameter,
                                            prior, column = NULL, offset = 0L) {

  if (is.null(posterior_samples_list)) {
    return(NULL)
  }

  indicators_list <- lapply(posterior_samples_list, function(samples) {

    .summary_models_indicators(
      posterior_samples = samples,
      parameter         = parameter,
      prior             = prior,
      column            = column,
      offset            = offset
    )
  })

  return(indicators_list)
}

.summary_models_marginal_indicator_matrices <- function(posterior_samples_list,
                                                        parameter, prior,
                                                        column = NULL,
                                                        offset = 0L) {

  indicators_list <- .summary_models_indicators_list(
    posterior_samples_list = posterior_samples_list,
    parameter              = parameter,
    prior                  = prior,
    column                 = column,
    offset                 = offset
  )
  if (is.null(indicators_list)) {
    return(NULL)
  }

  indicator_matrices <- lapply(indicators_list, function(indicators) {

    out <- matrix(
      as.numeric(outer(indicators, seq_len(length(prior)), "==")),
      nrow = length(indicators),
      ncol = length(prior)
    )
    colnames(out) <- paste0("component_", seq_len(length(prior)))

    return(out)
  })

  return(indicator_matrices)
}

.summary_models_individual_indicator_matrices <- function(components,
                                                          component_grid,
                                                          posterior_samples_list) {

  if (is.null(posterior_samples_list)) {
    return(NULL)
  }

  indicators_list <- lapply(components, function(info) {

    .summary_models_indicators_list(
      posterior_samples_list = posterior_samples_list,
      parameter              = info[["parameter"]],
      prior                  = info[["prior"]],
      column                 = info[["indicator"]],
      offset                 = info[["indicator_offset"]]
    )
  })

  indicator_matrices <- lapply(seq_along(posterior_samples_list), function(i) {

    out <- matrix(
      0,
      nrow = nrow(posterior_samples_list[[i]]),
      ncol = nrow(component_grid)
    )
    for (j in seq_len(nrow(component_grid))) {
      selected <- rep(TRUE, nrow(posterior_samples_list[[i]]))
      for (k in seq_along(components)) {
        selected <- selected & indicators_list[[k]][[i]] == component_grid[j, k]
      }
      out[, j] <- as.numeric(selected)
    }
    colnames(out) <- paste0("model_", seq_len(nrow(component_grid)))

    return(out)
  })

  return(indicator_matrices)
}

.summary_models_BF_error_percent <- function(indicator_matrices, prior_probs,
                                             post_probs, inclusion_BF) {

  if (is.null(indicator_matrices)) {
    return(NULL)
  }

  indicator_mcmc <- coda::mcmc.list(lapply(indicator_matrices, coda::mcmc))
  mcmc_summary   <- summary(indicator_mcmc, quantiles = NULL)[["statistics"]]
  if (is.null(dim(mcmc_summary))) {
    mcmc_summary <- t(mcmc_summary)
  }

  MCMC_error        <- mcmc_summary[, "Time-series SE"]
  BF_error_percent <- 100 * MCMC_error / (post_probs * (1 - post_probs))
  invalid          <- post_probs <= 0 | post_probs >= 1 |
    prior_probs <= 0 | prior_probs >= 1 |
    !is.finite(inclusion_BF) | !is.finite(MCMC_error)
  BF_error_percent[invalid] <- NA_real_

  return(BF_error_percent)
}

.summary_models_add_BF_diagnostics <- function(table, BF_error_percent) {

  if (is.null(BF_error_percent)) {
    return(table)
  }

  table[["BF_error_percent"]] <- BF_error_percent
  attr(table[["BF_error_percent"]], "name") <- "error%(Inclusion BF)"

  return(table)
}

.summary_models_BF_diagnostic_types <- function(BF_error_percent) {

  if (is.null(BF_error_percent)) {
    return(character())
  }

  return("BF_error")
}

.summary_models_inclusion_BF <- function(prior_probs, post_probs) {

  inclusion_BF <- vapply(seq_along(prior_probs), function(i) {

    BayesTools::inclusion_BF(
      prior_probs = prior_probs,
      post_probs  = post_probs,
      is_null     = seq_along(prior_probs) != i
    )
  }, numeric(1))

  return(inclusion_BF)
}
