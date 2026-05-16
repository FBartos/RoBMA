# ============================================================================ #
# IWMDE Parameter Specs, Components, and Priors
# ============================================================================ #

.iwmde_parameter_spec <- function(context, parameter, parameter_spec = NULL) {

  samples <- context[["posterior_samples"]]

  if (is.null(parameter_spec)) {
    if (!parameter %in% colnames(samples)) {
      return(list(status = "unsupported", reason = "posterior column is missing"))
    }
    if (.iwmde_parameter_is_indicator(parameter)) {
      return(list(status = "unsupported", reason = "indicator columns are discrete"))
    }
    if (.iwmde_parameter_is_weightfunction_coordinate(parameter)) {
      return(list(
        status = "unsupported",
        reason = "weightfunction coordinates need a joint omega/eta replacement map"
      ))
    }

    return(list(type = "primitive", parameter = parameter, status = "ok"))
  }

  if (!identical(parameter_spec[["type"]], "linear")) {
    return(list(status = "unsupported", reason = "unknown IWMDE parameter spec"))
  }

  weights <- .iwmde_linear_weights(parameter_spec[["weights"]])
  if (is.null(weights)) {
    return(list(status = "unsupported", reason = "linear weights are unavailable"))
  }
  if (length(weights) == 0L) {
    return(list(status = "unsupported", reason = "linear weights are all zero"))
  }

  missing <- setdiff(names(weights), colnames(samples))
  if (length(missing) > 0L) {
    return(list(
      status = "unsupported",
      reason = paste0(
        "linear weight columns are missing: ",
        paste(missing, collapse = ", ")
      )
    ))
  }

  bad <- names(weights)[
    vapply(names(weights), function(name) {
      .iwmde_parameter_is_indicator(name) ||
        .iwmde_parameter_is_weightfunction_coordinate(name)
    }, logical(1))
  ]
  if (length(bad) > 0L) {
    return(list(
      status = "unsupported",
      reason = paste0(
        "linear weights include unsupported columns: ",
        paste(bad, collapse = ", ")
      )
    ))
  }

  parameter_spec[["weights"]] <- weights
  parameter_spec[["status"]]  <- "ok"

  return(parameter_spec)
}


.iwmde_linear_weights <- function(weights) {

  if (is.null(weights)) {
    return(NULL)
  }
  if (is.data.frame(weights)) {
    weights <- as.matrix(weights)
  }
  if (is.matrix(weights)) {
    if (nrow(weights) != 1L) {
      return(NULL)
    }
    weights <- weights[1, ]
  }
  weight_names <- names(weights)
  weights <- as.numeric(weights)
  names(weights) <- weight_names
  if (is.null(names(weights)) || any(names(weights) == "")) {
    return(NULL)
  }

  keep <- is.finite(weights) & abs(weights) > sqrt(.Machine$double.eps)
  weights <- weights[keep]

  return(weights)
}


.iwmde_parameter_values <- function(context, parameter, parameter_spec) {

  samples <- context[["posterior_samples"]]
  if (identical(parameter_spec[["type"]], "linear")) {
    return(.iwmde_linear_values(samples, parameter_spec[["weights"]]))
  }

  return(as.numeric(samples[, parameter]))
}


.iwmde_linear_values <- function(samples, weights) {

  columns <- names(weights)

  return(as.numeric(as.matrix(samples[, columns, drop = FALSE]) %*% weights))
}


.iwmde_parameter_is_indicator <- function(parameter) {

  return(parameter == "bias_indicator" || grepl("(^|_)indicator$", parameter))
}


.iwmde_parameter_is_weightfunction_coordinate <- function(parameter) {

  return(grepl("^(omega|log_omega|eta)(\\[|$)", parameter))
}


.iwmde_parameter_is_eta <- function(parameter) {

  return(grepl("^eta(\\[|$)", parameter))
}


.iwmde_unsupported <- function(parameter, reason) {

  out <- list(
    parameter = parameter,
    status    = "unsupported",
    reason    = reason
  )
  class(out) <- c("iwmde_parameter_diagnostic", "list")
  return(out)
}


.iwmde_point_only_diagnostic <- function(parameter, samples, component) {

  out <- list(
    parameter    = parameter,
    status       = "point_only",
    samples      = samples,
    point_masses = component[["point_masses"]],
    reason       = "parameter has no continuous active posterior component"
  )
  class(out) <- c("iwmde_parameter_diagnostic", "list")
  return(out)
}


.iwmde_max_or_na <- function(x) {

  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }

  return(max(x))
}


.iwmde_parameter_components <- function(context, parameter,
                                        parameter_spec = NULL) {

  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "linear")) {
    return(.iwmde_linear_components(context, parameter_spec))
  }

  samples <- context[["posterior_samples"]]
  n       <- nrow(samples)

  static <- .iwmde_static_parameter_components(context, parameter, n)
  if (!is.null(static)) {
    return(static)
  }

  active  <- rep(TRUE, n)
  points  <- data.frame(x = numeric(), mass = numeric())

  point_location <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    state <- .iwmde_focal_prior_state(context, parameter, samples[i, ])
    if (identical(state[["status"]], "unsupported")) {
      active[i] <- FALSE
    } else if (identical(state[["status"]], "point")) {
      active[i]         <- FALSE
      point_location[i] <- state[["location"]]
    } else {
      active[i] <- TRUE
    }
  }

  point_location <- point_location[is.finite(point_location)]
  if (length(point_location) > 0L) {
    point_table <- table(signif(point_location, 12))
    points <- data.frame(
      x    = as.numeric(names(point_table)),
      mass = as.numeric(point_table) / n,
      row.names = NULL
    )
  }

  return(list(active = active, point_masses = points))
}


.iwmde_linear_components <- function(context, parameter_spec) {

  samples        <- context[["posterior_samples"]]
  weights        <- parameter_spec[["weights"]]
  n              <- nrow(samples)
  active         <- rep(FALSE, n)
  point_location <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    row         <- samples[i, ]
    point_sum   <- 0
    unsupported <- FALSE

    for (parameter in names(weights)) {
      state <- .iwmde_focal_prior_state(context, parameter, row)

      if (identical(state[["status"]], "unsupported")) {
        unsupported <- TRUE
        break
      } else if (identical(state[["status"]], "point")) {
        point_sum <- point_sum + weights[[parameter]] * state[["location"]]
      } else if (identical(state[["status"]], "continuous")) {
        active[i] <- TRUE
      }
    }

    if (!active[i] && !unsupported) {
      point_location[i] <- point_sum
    }
  }

  points         <- data.frame(x = numeric(), mass = numeric())
  point_location <- point_location[is.finite(point_location)]
  if (length(point_location) > 0L) {
    point_table <- table(signif(point_location, 12))
    points <- data.frame(
      x    = as.numeric(names(point_table)),
      mass = as.numeric(point_table) / n,
      row.names = NULL
    )
  }

  return(list(active = active, point_masses = points))
}


.iwmde_static_parameter_components <- function(context, parameter, n) {

  prior_name <- .iwmde_parameter_prior_name(context, parameter)
  flat_prior <- context[["flat_prior_list"]]
  points     <- data.frame(x = numeric(), mass = numeric())

  if (is.null(prior_name) || !prior_name %in% names(flat_prior)) {
    return(list(active = rep(TRUE, n), point_masses = points))
  }

  prior <- flat_prior[[prior_name]]
  if (BayesTools::is.prior.mixture(prior)) {
    return(NULL)
  }
  if (identical(prior_name, "bias")) {
    if (parameter == "PET" && !BayesTools::is.prior.PET(prior)) {
      return(.iwmde_static_point_components(n, 0))
    }
    if (parameter == "PEESE" && !BayesTools::is.prior.PEESE(prior)) {
      return(.iwmde_static_point_components(n, 0))
    }
  }
  if (BayesTools::is.prior.none(prior)) {
    return(.iwmde_static_point_components(n, 0))
  }
  if (BayesTools::is.prior.point(prior)) {
    return(.iwmde_static_point_components(
      n        = n,
      location = prior[["parameters"]][["location"]]
    ))
  }
  if (BayesTools::is.prior.weightfunction(prior)) {
    return(NULL)
  }

  return(list(active = rep(TRUE, n), point_masses = points))
}


.iwmde_static_point_components <- function(n, location) {

  return(list(
    active = rep(FALSE, n),
    point_masses = data.frame(
      x         = location,
      mass      = 1,
      row.names = NULL
    )
  ))
}


.iwmde_focal_prior_state <- function(context, parameter, row) {

  key <- paste(parameter, .iwmde_active_key(context, row), sep = "|")
  if (exists(key, envir = context[["focal_prior_cache"]], inherits = FALSE)) {
    return(get(key, envir = context[["focal_prior_cache"]]))
  }

  prior <- .iwmde_focal_prior(context, parameter, row)

  if (is.null(prior)) {
    out <- list(status = "continuous", prior = NULL)
    assign(key, out, envir = context[["focal_prior_cache"]])
    return(out)
  }
  if (BayesTools::is.prior.none(prior)) {
    out <- list(status = "point", location = 0)
    assign(key, out, envir = context[["focal_prior_cache"]])
    return(out)
  }
  if (BayesTools::is.prior.point(prior)) {
    out <- list(status = "point", location = prior[["parameters"]][["location"]])
    assign(key, out, envir = context[["focal_prior_cache"]])
    return(out)
  }
  if (BayesTools::is.prior.weightfunction(prior) ||
      grepl("^(omega|log_omega)(\\[|$)", parameter)) {
    out <- list(status = "unsupported")
    assign(key, out, envir = context[["focal_prior_cache"]])
    return(out)
  }

  out <- list(status = "continuous", prior = prior)
  assign(key, out, envir = context[["focal_prior_cache"]])
  return(out)
}


.iwmde_focal_prior <- function(context, parameter, row) {

  flat_prior <- context[["flat_prior_list"]]
  prior_name <- .iwmde_parameter_prior_name(context, parameter)
  if (is.null(prior_name) || !prior_name %in% names(flat_prior)) {
    return(NULL)
  }

  prior <- flat_prior[[prior_name]]
  if (BayesTools::is.prior.mixture(prior)) {
    indicator <- .iwmde_indicator_name(prior_name)
    if (!indicator %in% names(row)) {
      return(NULL)
    }
    index <- as.integer(round(row[[indicator]]))
    if (!is.finite(index) || index < 1L || index > length(prior)) {
      return(NULL)
    }
    prior <- prior[[index]]
  }

  if (identical(prior_name, "bias")) {
    if (parameter == "PET" && !BayesTools::is.prior.PET(prior)) {
      return(BayesTools::prior("point", list(location = 0)))
    }
    if (parameter == "PEESE" && !BayesTools::is.prior.PEESE(prior)) {
      return(BayesTools::prior("point", list(location = 0)))
    }
  }

  return(prior)
}


.iwmde_parameter_prior_name <- function(context, parameter) {

  flat_prior <- context[["flat_prior_list"]]
  if (parameter %in% names(flat_prior)) {
    return(parameter)
  }

  indexed <- regexec("^(.+)\\[[0-9]+\\]$", parameter)
  match   <- regmatches(parameter, indexed)[[1]]
  if (length(match) == 2L && match[2] %in% names(flat_prior)) {
    return(match[2])
  }

  if (parameter %in% c("PET", "PEESE") && "bias" %in% names(flat_prior)) {
    return("bias")
  }

  return(NULL)
}


.iwmde_indicator_name <- function(parameter) {

  if (identical(parameter, "bias")) {
    return("bias_indicator")
  }

  return(paste0(parameter, "_indicator"))
}


.iwmde_parameter_support <- function(context, parameter, rows,
                                     parameter_spec = NULL) {

  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "linear")) {
    return(.iwmde_linear_support(context, rows, parameter_spec[["weights"]]))
  }

  if (.iwmde_parameter_is_eta(parameter)) {
    return(c(0, Inf))
  }

  samples <- context[["posterior_samples"]]
  lower   <- Inf
  upper   <- -Inf
  keys    <- vapply(rows, function(row) {
    .iwmde_active_key(context, samples[row, ])
  }, character(1))
  rows    <- rows[!duplicated(keys)]

  for (row in rows) {
    prior <- .iwmde_focal_prior_cached(
      context    = context,
      parameter  = parameter,
      row        = samples[row, ],
      active_key = .iwmde_active_key(context, samples[row, ])
    )
    prior_support <- .iwmde_prior_support(prior)
    lower <- min(lower, prior_support[1])
    upper <- max(upper, prior_support[2])
  }

  if (!is.finite(lower)) lower <- -Inf
  if (!is.finite(upper)) upper <- Inf

  return(c(lower, upper))
}


.iwmde_linear_support <- function(context, rows, weights) {

  samples <- context[["posterior_samples"]]
  lower   <- Inf
  upper   <- -Inf
  keys    <- vapply(rows, function(row_index) {
    .iwmde_active_key(context, samples[row_index, ])
  }, character(1))
  rows    <- rows[!duplicated(keys)]

  for (row_index in rows) {
    row       <- samples[row_index, ]
    row_lower <- 0
    row_upper <- 0

    for (parameter in names(weights)) {
      weight <- weights[[parameter]]
      state  <- .iwmde_focal_prior_state(context, parameter, row)
      if (identical(state[["status"]], "point")) {
        support <- rep(state[["location"]], 2L)
      } else {
        support <- .iwmde_prior_support(state[["prior"]])
      }

      if (weight >= 0) {
        row_lower <- row_lower + weight * support[1]
        row_upper <- row_upper + weight * support[2]
      } else {
        row_lower <- row_lower + weight * support[2]
        row_upper <- row_upper + weight * support[1]
      }
    }

    lower <- min(lower, row_lower)
    upper <- max(upper, row_upper)
  }

  if (!is.finite(lower)) lower <- -Inf
  if (!is.finite(upper)) upper <- Inf

  return(c(lower, upper))
}


.iwmde_focal_prior_cached <- function(context, parameter, row,
                                      active_key = NULL) {

  if (is.null(active_key)) {
    active_key <- .iwmde_active_key(context, row)
  }
  key <- paste("support_prior", parameter, active_key, sep = "|")
  if (exists(key, envir = context[["support_cache"]], inherits = FALSE)) {
    return(get(key, envir = context[["support_cache"]], inherits = FALSE))
  }

  prior <- .iwmde_focal_prior(context, parameter, row)
  assign(key, prior, envir = context[["support_cache"]])

  return(prior)
}


.iwmde_prior_support <- function(prior) {

  if (is.null(prior) || BayesTools::is.prior.none(prior)) {
    return(c(-Inf, Inf))
  }
  if (!is.null(prior[["truncation"]])) {
    return(c(prior[["truncation"]][["lower"]], prior[["truncation"]][["upper"]]))
  }

  return(c(-Inf, Inf))
}


