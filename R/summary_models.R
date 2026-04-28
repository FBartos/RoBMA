#' @title Summarize RoBMA Model Weights
#'
#' @description Creates marginal or individual model-weight summaries for
#' RoBMA product-space objects.
#'
#' @param object a fitted RoBMA object
#' @param type whether to summarize marginal component prior distributions
#'   (\code{"marginal"}) or individual model combinations (\code{"individual"}).
#' @param ... additional arguments
#'
#' @return A list of class \code{summary_models.RoBMA}.
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
summary_models.RoBMA <- function(object, type = "marginal", ...) {

  BayesTools::check_char(type, "type")
  type <- match.arg(type, c("marginal", "individual"))

  if (is.null(object[["fit"]]) || length(object[["fit"]]) == 0) {
    stop("'summary_models' requires a fitted RoBMA object.", call. = FALSE)
  }

  components        <- .summary_models_components(object)
  posterior_samples <- .get_posterior_samples(object[["fit"]])

  if (length(components) == 0L) {
    stop("No mixture-prior model components found.", call. = FALSE)
  }

  output <- list(
    name = .summary.brma_model_names(object),
    type = type
  )

  if (type == "marginal") {
    output[["marginal"]] <- .summary_models_marginal_tables(
      components        = components,
      posterior_samples = posterior_samples
    )
  } else {
    output[["individual"]] <- .summary_models_individual_table(
      components        = components,
      posterior_samples = posterior_samples
    )
  }

  class(output) <- "summary_models.RoBMA"

  return(output)
}

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

  scale_parameters <- grep("^log_tau_", names(prior_list), value = TRUE)
  scale_parameters <- scale_parameters[scale_parameters != "log_tau_intercept"]
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

  return(components)
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
    component   = component,
    parameter   = parameter,
    prior       = prior,
    names       = .summary_models_prior_names(prior),
    hypothesis  = .summary_models_hypothesis(prior),
    prior_probs = .summary_models_prior_probs(prior)
  )

  return(components)
}

.summary_models_marginal_tables <- function(components, posterior_samples) {

  tables <- list()

  for (component in names(components)) {

    info       <- components[[component]]
    indicators <- .summary_models_indicators(
      posterior_samples = posterior_samples,
      parameter         = info[["parameter"]]
    )
    post_probs <- vapply(seq_along(info[["names"]]), function(i) {

      mean(indicators == i)
    }, numeric(1))
    inclusion_BF <- .summary_models_inclusion_BF(
      prior_probs = info[["prior_probs"]],
      post_probs  = post_probs
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

    class(table)             <- c("BayesTools_table", "RoBMA_summary_models_marginal", class(table))
    attr(table, "type")      <- c("string_left", "prior_prob", "post_prob", "inclusion_BF")
    attr(table, "rownames")  <- TRUE
    attr(table, "title")     <- component
    attr(table, "footnotes") <- NULL
    attr(table, "warnings")  <- NULL

    tables[[component]] <- table
  }

  return(tables)
}

.summary_models_individual_table <- function(components, posterior_samples) {

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
      parameter         = info[["parameter"]]
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

  table <- cbind(
    prior_columns,
    data.frame(
      prior_prob   = prior_probs,
      post_prob    = post_probs,
      inclusion_BF = BayesTools::format_BF(inclusion_BF, inclusion = TRUE),
      check.names  = FALSE
    )
  )

  class(table)             <- c("BayesTools_table", "RoBMA_summary_models_individual", class(table))
  attr(table, "type")      <- c(
    rep("string_left", length(components)),
    "prior_prob", "post_prob", "inclusion_BF"
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

.summary_models_indicators <- function(posterior_samples, parameter) {

  column <- paste0(parameter, "_indicator")
  if (!column %in% colnames(posterior_samples)) {
    stop("Missing posterior model indicator: '", column, "'.", call. = FALSE)
  }

  return(as.integer(posterior_samples[, column]))
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
