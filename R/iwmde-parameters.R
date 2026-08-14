# ============================================================================ #
# IWMDE Parameter Specs, Components, and Priors
# ============================================================================ #

.iwmde_parameter_spec <- function(context, parameter, parameter_spec = NULL) {

  samples <- context[["posterior_samples"]]

  if (is.null(parameter_spec)) {
    if (!.iwmde_parameter_has_resolved_values(context, parameter)) {
      return(list(status = "unsupported", reason = "posterior column is missing"))
    }
    if (.iwmde_parameter_is_indicator(parameter)) {
      return(list(status = "unsupported", reason = "indicator columns are discrete"))
    }
    if (.iwmde_parameter_is_weightfunction_coordinate(
      parameter = parameter,
      context   = context
    )) {
      return(list(
        status = "unsupported",
        reason = "weightfunction coordinates need a joint omega/eta replacement map"
      ))
    }

    return(list(type = "primitive", parameter = parameter, status = "ok"))
  }

  if (identical(parameter_spec[["type"]], "primitive")) {
    if (!.iwmde_parameter_has_resolved_values(context, parameter)) {
      return(list(status = "unsupported", reason = "posterior column is missing"))
    }
    if (.iwmde_parameter_is_indicator(parameter)) {
      return(list(status = "unsupported", reason = "indicator columns are discrete"))
    }
    if (.iwmde_parameter_is_weightfunction_coordinate(
      parameter = parameter,
      context   = context
    )) {
      return(list(
        status = "unsupported",
        reason = "weightfunction coordinates need a joint omega/eta replacement map"
      ))
    }

    parameter_spec[["parameter"]] <- parameter
    parameter_spec[["status"]]    <- "ok"

    return(parameter_spec)
  }

  if (identical(parameter_spec[["type"]], "simplex_pair")) {
    source  <- parameter_spec[["parameter"]]
    index   <- parameter_spec[["index"]]
    columns <- paste0(source, "[", 1:2, "]")
    if (!is.character(source) || length(source) != 1L || is.na(source) ||
        !nzchar(source) || !is.numeric(index) || length(index) != 1L ||
        is.na(index) || !index %in% 1:2 ||
        !all(columns %in% colnames(samples))) {
      return(list(
        status = "unsupported",
        reason = "two-component simplex coordinates are unavailable"
      ))
    }
    parameter_spec[["index"]]     <- as.integer(index)
    parameter_spec[["n_targets"]] <- 2L
    parameter_spec[["status"]]    <- "ok"

    return(parameter_spec)
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
  missing <- missing[!vapply(missing, function(parameter) {
    .iwmde_parameter_has_resolved_values(context, parameter)
  }, logical(1))]
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
        .iwmde_parameter_is_weightfunction_coordinate(
          parameter = name,
          context   = context
        )
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
    if (!all(vapply(weights, is.numeric, logical(1)))) {
      stop("Linear target weights must be numeric.", call. = FALSE)
    }
    weights <- as.matrix(weights)
  }
  if (is.matrix(weights)) {
    if (nrow(weights) != 1L) {
      stop("Linear target weights must contain exactly one row.", call. = FALSE)
    }
    weights <- weights[1, ]
  }
  if (!is.numeric(weights)) {
    stop("Linear target weights must be numeric.", call. = FALSE)
  }
  weight_names <- names(weights)
  weights <- as.numeric(weights)
  names(weights) <- weight_names
  if (is.null(names(weights)) || anyNA(names(weights)) ||
      any(names(weights) == "")) {
    stop("Linear target weights must be fully named.", call. = FALSE)
  }
  duplicate_names <- unique(names(weights)[duplicated(names(weights))])
  if (length(duplicate_names) > 0L) {
    stop(
      "Linear target weights contain duplicate parameter name(s): ",
      paste(duplicate_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(weights))) {
    stop("Linear target weights must all be finite.", call. = FALSE)
  }

  keep <- weights != 0
  weights <- weights[keep]

  return(weights)
}


.iwmde_parameter_has_resolved_values <- function(context, parameter) {

  samples <- context[["posterior_samples"]]
  if (parameter %in% colnames(samples)) {
    return(TRUE)
  }

  return(.iwmde_parameter_is_point_only(context, parameter))
}


.iwmde_parameter_is_point_only <- function(context, parameter) {

  samples <- context[["posterior_samples"]]
  if (is.null(.iwmde_parameter_prior_name(context, parameter)) ||
      is.null(samples) ||
      nrow(samples) == 0L) {
    return(FALSE)
  }

  return(all(vapply(seq_len(nrow(samples)), function(i) {
    state <- .iwmde_focal_prior_state(context, parameter, samples[i, ])

    identical(state[["status"]], "point")
  }, logical(1))))
}


.iwmde_parameter_values <- function(context, parameter, parameter_spec) {

  samples <- context[["posterior_samples"]]
  if (identical(parameter_spec[["type"]], "linear")) {
    return(.iwmde_linear_values(context, samples, parameter_spec[["weights"]]))
  }

  return(.iwmde_parameter_column_values(context, samples, parameter))
}


.iwmde_parameter_column_values <- function(context, samples, parameter) {

  return(vapply(seq_len(nrow(samples)), function(i) {
    .iwmde_parameter_value_row(context, samples[i, ], parameter)
  }, numeric(1)))
}


.iwmde_parameter_value_row <- function(context, row, parameter) {

  state <- .iwmde_focal_prior_state(context, parameter, row)
  if (identical(state[["status"]], "point")) {
    return(as.numeric(state[["location"]]))
  }

  if (parameter %in% names(row)) {
    return(as.numeric(row[[parameter]]))
  }

  return(NA_real_)
}


.iwmde_linear_values <- function(context, samples, weights) {

  columns <- names(weights)
  values <- vapply(columns, function(parameter) {
    .iwmde_parameter_column_values(context, samples, parameter)
  }, numeric(nrow(samples)))

  return(as.numeric(values %*% weights))
}


.iwmde_parameter_is_indicator <- function(parameter) {

  return(parameter == "bias_indicator" || grepl("(^|_)indicator$", parameter))
}


.iwmde_selection_coordinate_names <- function(context = NULL,
                                              include_private = TRUE) {

  selection_spec <- NULL
  if (is.list(context) && !is.null(context[["selection_spec"]])) {
    selection_spec <- context[["selection_spec"]]
  } else if (is.list(context) && !is.null(context[["jags_omega"]])) {
    selection_spec <- context
  }

  if (!is.null(selection_spec)) {
    coordinates <- selection_spec[["jags_omega"]]
  } else {
    coordinates <- character()
  }
  if (isTRUE(include_private)) {
    coordinates <- c(coordinates, "eta")
  }

  coordinates <- unique(coordinates[!is.na(coordinates) & nzchar(coordinates)])
  return(coordinates)
}


.iwmde_parameter_matches_coordinate <- function(parameter, coordinates) {

  coordinates <- unique(coordinates[!is.na(coordinates) & nzchar(coordinates)])
  if (length(coordinates) == 0L) {
    return(FALSE)
  }

  matches <- parameter %in% coordinates
  for (coordinate in coordinates) {
    matches <- matches | BayesTools::JAGS_indexed_parameter_columns(
      columns   = parameter,
      parameter = coordinate
    )
  }

  return(any(matches))
}


.iwmde_parameter_is_weightfunction_coordinate <- function(parameter, context = NULL,
                                                         include_private = TRUE) {

  return(.iwmde_parameter_matches_coordinate(
    parameter   = parameter,
    coordinates = .iwmde_selection_coordinate_names(
      context         = context,
      include_private = include_private
    )
  ))
}


.iwmde_parameter_is_eta <- function(parameter) {

  return(grepl("^eta(\\[|$)", parameter))
}


.iwmde_unsupported <- function(parameter, reason) {

  return(.iwmde_new_diagnostic(list(
    parameter = parameter,
    status    = "unsupported",
    reason    = reason
  )))
}


.iwmde_point_only_diagnostic <- function(parameter, samples, component) {

  return(.iwmde_new_diagnostic(list(
    parameter    = parameter,
    status       = "point_only",
    samples      = samples,
    point_masses = component[["point_masses"]],
    reason       = "parameter has no continuous active posterior component"
  )))
}


.iwmde_max_or_na <- function(x) {

  x <- x[!is.na(x)]
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

  points <- .iwmde_point_mass_table(point_location, denominator = n)

  return(list(
    active         = active,
    point_masses   = points,
    point_location = point_location
  ))
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

  points <- .iwmde_point_mass_table(point_location, denominator = n)

  return(list(
    active         = active,
    point_masses   = points,
    point_location = point_location
  ))
}


.iwmde_static_parameter_components <- function(context, parameter, n) {

  prior_name <- .iwmde_parameter_prior_name(context, parameter)
  flat_prior <- context[["flat_prior_list"]]
  points     <- data.frame(x = numeric(), mass = numeric())

  if (is.null(prior_name) || !prior_name %in% names(flat_prior)) {
    return(list(
      active         = rep(TRUE, n),
      point_masses   = points,
      point_location = rep(NA_real_, n)
    ))
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

  return(list(
    active         = rep(TRUE, n),
    point_masses   = points,
    point_location = rep(NA_real_, n)
  ))
}


.iwmde_static_point_components <- function(n, location) {

  return(list(
    active = rep(FALSE, n),
    point_masses = data.frame(
      x         = location,
      mass      = 1,
      row.names = NULL
    ),
    point_location = rep(location, n)
  ))
}


.iwmde_parameter_condition_rows <- function(context, parameter_spec) {

  n           <- nrow(context[["posterior_samples"]])
  conditional <- parameter_spec[["conditional"]]
  if (is.null(conditional) || length(conditional) == 0L) {
    return(rep(TRUE, n))
  }

  conditional <- unique(as.character(conditional))
  conditional <- conditional[!is.na(conditional) & nzchar(conditional)]
  if (length(conditional) == 0L) {
    return(rep(TRUE, n))
  }

  rule <- parameter_spec[["conditional_rule"]]
  if (is.null(rule) || length(rule) == 0L) {
    rule <- "OR"
  }
  rule <- match.arg(rule, c("AND", "OR"))

  rows <- matrix(FALSE, nrow = n, ncol = length(conditional))
  for (i in seq_along(conditional)) {
    rows[, i] <- .iwmde_parameter_active_rows(context, conditional[[i]])
  }

  if (identical(rule, "AND")) {
    return(apply(rows, 1L, all))
  }

  return(apply(rows, 1L, any))
}


.iwmde_parameter_active_rows <- function(context, parameter) {

  samples <- context[["posterior_samples"]]
  if (!parameter %in% colnames(samples) &&
      is.null(.iwmde_parameter_prior_name(context, parameter))) {
    return(rep(FALSE, nrow(samples)))
  }

  return(vapply(seq_len(nrow(samples)), function(i) {
    state <- .iwmde_focal_prior_state(context, parameter, samples[i, ])
    identical(state[["status"]], "continuous")
  }, logical(1)))
}


.iwmde_restrict_parameter_component <- function(component, rows) {

  component[["active"]] <- component[["active"]] & rows

  point_location <- component[["point_location"]]
  if (is.null(point_location)) {
    point_location <- rep(NA_real_, length(rows))
  }
  point_location[!rows] <- NA_real_
  component[["point_location"]] <- point_location

  component[["point_masses"]] <- .iwmde_point_mass_table(
    locations   = point_location,
    denominator = sum(rows)
  )

  return(component)
}


.iwmde_point_mass_table <- function(locations, denominator) {

  locations <- locations[is.finite(locations)]
  if (length(locations) == 0L) {
    return(data.frame(x = numeric(), mass = numeric()))
  }
  if (!is.numeric(denominator) || length(denominator) != 1L ||
      !is.finite(denominator) || denominator <= 0) {
    stop("Point-mass denominator must be finite and positive.", call. = FALSE)
  }

  values <- sort(unique(locations))
  counts <- vapply(values, function(value) {
    sum(locations == value)
  }, integer(1L))

  return(data.frame(
    x         = values,
    mass      = counts / denominator,
    row.names = NULL
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
      .iwmde_parameter_is_weightfunction_coordinate(
        parameter = parameter,
        context   = context
      )) {
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
    index <- .iwmde_indicator_index(
      row[[indicator]],
      indicator,
      max_index = length(prior)
    )
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


.iwmde_parameter_controls_sampled_random_sd <- function(context, parameter) {

  if (is.null(parameter) || length(parameter) != 1L ||
      is.na(parameter) || !nzchar(parameter)) {
    return(FALSE)
  }

  data <- context[["data"]]
  if (is.null(data) || !.is_data_random(data)) {
    return(FALSE)
  }

  terms <- .iwmde_sampled_random_effect_terms(context)
  if (length(terms) == 0L) {
    return(FALSE)
  }

  prior_name <- .iwmde_parameter_prior_name(context, parameter)
  candidates <- unique(c(parameter, prior_name))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  sd_columns <- unique(unlist(
    lapply(terms, .iwmde_random_term_sd_columns),
    use.names = FALSE
  ))
  sd_columns <- sd_columns[!is.na(sd_columns) & nzchar(sd_columns)]
  if (length(intersect(candidates, sd_columns)) > 0L) {
    return(TRUE)
  }

  if (.is_data_scale(data) &&
      length(intersect(.data_scale_formula_sources(data), sd_columns)) > 0L) {
    scale_parameters <- .data_scale_formula_parameters(data)
    for (scale_parameter in scale_parameters) {
      if (any(candidates == scale_parameter) ||
          any(startsWith(candidates, paste0(scale_parameter, "_")))) {
        return(TRUE)
      }
    }
  }

  return(FALSE)
}


.iwmde_sampled_random_effect_terms <- function(context) {

  object <- context[["object"]]
  design <- if (!is.null(object)) {
    .fitted_formula_design(object, "mu", required = FALSE)
  } else {
    NULL
  }

  if (!is.null(design[["random_effects"]])) {
    terms <- design[["random_effects"]]
  } else {
    terms <- .data_effective_sampled_random_effect_terms(context[["data"]])
  }

  if (length(terms) == 0L) {
    return(list())
  }

  terms[vapply(terms, function(term) {
    identical(.random_effect_term_compile_mode(term), "sampled")
  }, logical(1))]
}


.iwmde_random_term_sd_columns <- function(term) {

  columns <- term[["sd_parameter_names"]]

  sources <- .predict_known_v_marginalized_sd_sources(term)
  for (source in sources) {
    columns <- c(columns, source[["name"]])
  }

  if (.marginalized_random_effect_has_allocation(term)) {
    factors <- .marginalized_random_effect_allocation_factors(term)
    columns <- c(columns, vapply(factors, function(factor) {
      paste0(factor[["weight_name"]], "[", factor[["index"]], "]")
    }, character(1)))
  }

  unique(columns[!is.na(columns) & nzchar(columns)])
}


.iwmde_indicator_name <- function(parameter) {

  if (identical(parameter, "bias")) {
    return("bias_indicator")
  }

  return(paste0(parameter, "_indicator"))
}


.iwmde_parameter_support <- function(context, parameter, rows,
                                     parameter_spec = NULL) {

  supports <- .iwmde_parameter_row_supports(
    context        = context,
    parameter      = parameter,
    rows           = rows,
    parameter_spec = parameter_spec
  )
  if (nrow(supports) == 0L) {
    return(c(-Inf, Inf))
  }

  lower <- min(supports[, 1L], na.rm = TRUE)
  upper <- max(supports[, 2L], na.rm = TRUE)
  if (!is.finite(lower)) lower <- -Inf
  if (!is.finite(upper)) upper <- Inf

  return(c(lower, upper))
}


.iwmde_parameter_row_supports <- function(context, parameter, rows,
                                          parameter_spec = NULL) {

  if (length(rows) == 0L) {
    return(matrix(numeric(), nrow = 0L, ncol = 2L))
  }

  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "linear")) {
    return(.iwmde_linear_row_supports(
      context = context,
      rows    = rows,
      weights = parameter_spec[["weights"]]
    ))
  }
  if (!is.null(parameter_spec) &&
      identical(parameter_spec[["type"]], "simplex_pair")) {
    return(matrix(
      c(0, 1),
      nrow  = length(rows),
      ncol  = 2L,
      byrow = TRUE
    ))
  }

  if (.iwmde_parameter_is_eta(parameter)) {
    return(matrix(
      c(0, Inf),
      nrow = length(rows),
      ncol = 2L,
      byrow = TRUE
    ))
  }

  samples <- context[["posterior_samples"]]
  out     <- matrix(NA_real_, nrow = length(rows), ncol = 2L)

  for (i in seq_along(rows)) {
    row <- rows[[i]]
    prior <- .iwmde_focal_prior_cached(
      context    = context,
      parameter  = parameter,
      row        = samples[row, ],
      active_key = .iwmde_active_key(context, samples[row, ])
    )
    out[i, ] <- .iwmde_prior_support(prior)
  }

  return(out)
}


.iwmde_linear_row_supports <- function(context, rows, weights) {

  samples <- context[["posterior_samples"]]
  out     <- matrix(NA_real_, nrow = length(rows), ncol = 2L)

  for (i in seq_along(rows)) {
    row          <- samples[rows[[i]], ]
    current      <- .iwmde_linear_value_row(context, row, weights)
    active_names <- .iwmde_linear_active_columns(context, row, weights)

    if (!is.finite(current)) {
      out[i, ] <- c(NA_real_, NA_real_)
      next
    }
    if (length(active_names) == 0L) {
      out[i, ] <- rep(current, 2L)
      next
    }

    active_weights <- weights[active_names]
    denominator    <- sum(active_weights^2)
    if (!is.finite(denominator) || denominator <= 0) {
      out[i, ] <- c(NA_real_, NA_real_)
      next
    }

    coefficients <- active_weights / denominator
    row_lower    <- -Inf
    row_upper    <- Inf

    for (parameter in active_names) {
      coefficient <- coefficients[[parameter]]
      state  <- .iwmde_focal_prior_state(context, parameter, row)
      support <- .iwmde_prior_support(state[["prior"]])
      value   <- .iwmde_parameter_value_row(context, row, parameter)

      if (!is.finite(coefficient) || coefficient == 0) {
        next
      }

      if (coefficient > 0) {
        lower <- current + (support[1] - value) / coefficient
        upper <- current + (support[2] - value) / coefficient
      } else {
        lower <- current + (support[2] - value) / coefficient
        upper <- current + (support[1] - value) / coefficient
      }

      row_lower <- max(row_lower, lower)
      row_upper <- min(row_upper, upper)
    }

    out[i, ] <- c(row_lower, row_upper)
  }

  return(out)
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


# Classify exact requested ordinates from deterministic prior provenance.
.iwmde_prepare_prior_ordinates <- function(context, parameter, parameter_spec,
                                           values) {

  if (is.null(parameter_spec)) {
    parameter_spec <- list(type = "primitive")
  }
  if (!is.null(parameter_spec[["prior_ordinates"]])) {
    requested <- .iwmde_unique_ordinate_values(values)
    prepared  <- parameter_spec[["prior_ordinates"]]
    valid <- is.list(prepared) && length(prepared) == length(requested) &&
      (length(prepared) == 0L || identical(names(prepared), names(requested)))
    if (valid && length(prepared) > 0L) {
      for (i in seq_along(prepared)) {
        valid <- isTRUE(tryCatch({
          .iwmde_validate_prior_ordinate(prepared[[i]], unname(requested[[i]]))
          TRUE
        }, error = function(error) FALSE))
        if (!valid) {
          break
        }
      }
    }
    if (!isTRUE(valid)) {
      stop("Prepared prior ordinates do not match the current request.",
           call. = FALSE)
    }
    parameter_spec[["prior_density"]] <- NULL
    return(parameter_spec)
  }

  prior_density <- parameter_spec[["prior_density"]]
  conditional   <- parameter_spec[["conditional"]]
  if (is.null(prior_density) &&
      identical(parameter_spec[["type"]], "primitive") &&
      (is.null(conditional) || length(conditional) == 0L)) {
    prior_name <- .iwmde_parameter_prior_name(context, parameter)
    if (!is.null(prior_name)) {
      prior_density <- context[["flat_prior_list"]][[prior_name]]
    }
  }

  classifications <- .iwmde_prior_ordinate_classifications(
    prior_density = prior_density,
    values        = values
  )
  parameter_spec[["prior_density"]] <- NULL
  parameter_spec["prior_ordinates"] <- list(classifications)

  return(parameter_spec)
}


# Classify each exact finite value once through the BayesTools contract.
.iwmde_prior_ordinate_classifications <- function(prior_density, values) {

  values <- .iwmde_unique_ordinate_values(values)
  if (length(values) == 0L) {
    return(list())
  }

  classifications <- lapply(unname(values), function(value) {
    result <- if (is.null(prior_density)) {
      .iwmde_unknown_prior_ordinate(value)
    } else {
      BayesTools::prior_density_ordinate(prior_density, value)
    }

    .iwmde_validate_prior_ordinate(result, value)
    unclass(result)
  })
  names(classifications) <- names(values)

  return(classifications)
}


# Retain unique finite values under exact binary64 keys.
.iwmde_unique_ordinate_values <- function(values) {

  values <- as.numeric(values)
  values <- values[is.finite(values)]
  keys   <- vapply(values, .iwmde_key_number, character(1))
  keep   <- !duplicated(keys)
  values <- values[keep]
  names(values) <- keys[keep]

  return(values)
}


# Represent a target whose deterministic prior density is unavailable.
.iwmde_unknown_prior_ordinate <- function(value) {

  return(list(
    schema_version = "1",
    value          = as.numeric(value),
    behavior       = "unknown",
    log_density    = NA_real_,
    point_mass     = 0,
    exact          = FALSE,
    method         = "unsupported_provenance",
    reason         = "The deterministic target prior density is unavailable.",
    provenance     = list(kind = "missing_target_density")
  ))
}


# Fail closed if the installed BayesTools ordinate schema is incompatible.
.iwmde_validate_prior_ordinate <- function(result, value) {

  fields <- c(
    "schema_version", "value", "behavior", "log_density", "point_mass",
    "exact", "method", "reason", "provenance"
  )
  behaviors <- c(
    "regular", "zero", "infinite", "point_mass", "undefined", "unknown"
  )
  methods <- c(
    "primitive", "point", "finite_mixture", "scalar_affine",
    "linear_normal", "named_transform", "unsupported_provenance"
  )
  valid <- is.list(result) && identical(names(result), fields) &&
    identical(result[["schema_version"]], "1") &&
    is.numeric(result[["value"]]) && length(result[["value"]]) == 1L &&
    identical(as.numeric(result[["value"]]), as.numeric(value)) &&
    is.character(result[["behavior"]]) &&
    length(result[["behavior"]]) == 1L &&
    !is.na(result[["behavior"]]) && result[["behavior"]] %in% behaviors &&
    is.numeric(result[["log_density"]]) &&
    length(result[["log_density"]]) == 1L &&
    is.numeric(result[["point_mass"]]) &&
    length(result[["point_mass"]]) == 1L &&
    is.logical(result[["exact"]]) && length(result[["exact"]]) == 1L &&
    !is.na(result[["exact"]]) &&
    is.character(result[["method"]]) && length(result[["method"]]) == 1L &&
    !is.na(result[["method"]]) && result[["method"]] %in% methods &&
    (is.null(result[["reason"]]) ||
      (is.character(result[["reason"]]) && length(result[["reason"]]) == 1L)) &&
    is.list(result[["provenance"]])
  if (!isTRUE(valid)) {
    stop(
      "BayesTools returned an incompatible prior-density ordinate result.",
      call. = FALSE
    )
  }

  invisible(result)
}


# Warn before estimation when the target prior ordinate is nonregular.
.iwmde_ordinate_prior_warnings <- function(parameter, prior_ordinates) {

  if (length(prior_ordinates) == 0L) {
    return(character())
  }

  descriptions <- c(
    zero       = "tends to zero",
    infinite   = "is singular (tends to infinity)",
    point_mass = "assigns positive point mass",
    undefined  = "is undefined"
  )
  warnings <- unlist(lapply(prior_ordinates, function(ordinate) {
    behavior <- ordinate[["behavior"]]
    if (identical(behavior, "regular")) {
      return(character())
    }
    if (identical(behavior, "unknown")) {
      return(paste0(
        "The qCMDE/IWMDE target prior density for '", parameter, "' at ",
        format(ordinate[["value"]], digits = 17L, trim = TRUE),
        " could not be classified from deterministic provenance. ",
        "Nonregularity cannot be ruled out; the exact requested value is ",
        "retained and the Bayes factor may be unavailable."
      ))
    }

    paste0(
      "The qCMDE/IWMDE ordinate for '", parameter, "' at ",
      format(ordinate[["value"]], digits = 17L, trim = TRUE),
      " is nonregular because the target prior density ",
      descriptions[[behavior]],
      ". The exact requested value is retained; the Bayes factor may be unavailable."
    )
  }), use.names = FALSE)

  return(unique(warnings))
}
