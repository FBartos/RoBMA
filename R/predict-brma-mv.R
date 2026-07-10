# predict-brma-mv.R
# ============================================================================ #
#
# brma.mv-specific prediction helpers. These stay outside predict.R so the
# public prediction dispatch is not carrying known-V covariance construction.
#
# ============================================================================ #


.predict_known_v_newdata_add_vi <- function(newdata, V_new) {

  diag_v <- diag(V_new)
  if ("vi" %in% names(newdata)) {
    .predict_known_v_newdata_check_variance(
      supplied = newdata[["vi"]],
      expected = diag_v,
      label    = "vi"
    )
  }
  if ("sei" %in% names(newdata)) {
    .predict_known_v_newdata_check_variance(
      supplied = newdata[["sei"]]^2,
      expected = diag_v,
      label    = "sei"
    )
  }

  newdata[["vi"]] <- diag_v
  if ("sei" %in% names(newdata)) {
    newdata[["sei"]] <- NULL
  }

  return(newdata)
}


.predict_known_v_newdata_check_variance <- function(supplied, expected, label) {

  if (!is.numeric(supplied) || length(supplied) != length(expected) ||
      anyNA(supplied) || any(!is.finite(supplied))) {
    stop(
      "The '", label, "' column in 'newdata' must contain finite numeric ",
      "values matching diag(V_new).",
      call. = FALSE
    )
  }

  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(expected)))
  if (max(abs(supplied - expected)) > tolerance) {
    stop(
      "The '", label, "' column in 'newdata' must match diag(V_new).",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.predict_known_v_response_covariance <- function(object, data, known_V,
                                                 posterior_samples,
                                                 caller = "known-V response prediction") {

  S <- nrow(posterior_samples)
  K <- nrow(data[["outcome"]])

  .known_v_check_full_covariance_allocation(
    S      = S,
    K      = K,
    caller = caller
  )

  covariance_samples <- array(0, dim = c(S, K, K))
  formula_design     <- .predict_known_v_formula_design_with_row_source_values(
    object = object,
    data   = data
  )
  sampled_blocks     <- .formula_design_sampled_random_effect_blocks(formula_design)

  if (length(sampled_blocks) > 0L) {
    location_priors <- attr(object[["fit"]], "prior_list")
    if (is.null(location_priors)) {
      location_priors <- object[["priors"]][["location"]]
    }

    formula_fit <- .posterior_formula_fit(
      fit               = object[["fit"]],
      posterior_samples = posterior_samples,
      formula_design    = TRUE
    )
    attr(formula_fit, "formula_design") <- list(mu = formula_design)

    random_vcov <- tryCatch(
      BayesTools::random_effects_marginal_vcov(
        fit               = formula_fit,
        parameter         = "mu",
        data              = data[["location"]],
        posterior_samples = posterior_samples,
        prior_list        = location_priors,
        blocks            = sampled_blocks,
        new_levels        = "sample"
      ),
      error = function(e) {
        stop(
          "Cannot compute random-effect covariance for ", caller, ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    covariance_samples <- random_vcov[["samples"]]
    if (length(dim(covariance_samples)) != 3L ||
        dim(covariance_samples)[1L] != S ||
        dim(covariance_samples)[2L] != K ||
        dim(covariance_samples)[3L] != K) {
      stop("Random-effect covariance samples have inconsistent dimensions.",
           call. = FALSE)
    }
  }

  marginalized_variance <- .predict_known_v_newdata_marginalized_variance(
    object            = object,
    data              = data,
    posterior_samples = posterior_samples
  )

  for (s in seq_len(S)) {
    covariance_samples[s, , ] <- covariance_samples[s, , ] +
      diag(marginalized_variance[s, ], nrow = K, ncol = K)
  }

  return(.known_v_add_base_covariance(
    base_covariance    = known_V[["V"]],
    covariance_samples = covariance_samples
  ))
}


.predict_known_v_newdata_response_covariance <- function(object, data,
                                                         known_V_new,
                                                         posterior_samples) {

  .predict_known_v_response_covariance(
    object            = object,
    data              = data,
    known_V           = known_V_new,
    posterior_samples = posterior_samples,
    caller            = "known-V newdata response prediction"
  )
}


.predict_known_v_newdata_marginalized_variance <- function(object, data,
                                                           posterior_samples) {

  S <- nrow(posterior_samples)
  K <- nrow(data[["outcome"]])

  if (!.data_has_marginalized_random_effects(object[["data"]])) {
    return(matrix(0, nrow = S, ncol = K))
  }

  return(.evaluate_marginalized_random_variance(
    data              = object[["data"]],
    posterior_samples = posterior_samples,
    K                 = K,
    source_samples    = .predict_known_v_newdata_marginalized_source_samples(
      object            = object,
      data              = data,
      posterior_samples = posterior_samples
    )
  ))
}


.predict_known_v_marginalized_random_draws <- function(object, data,
                                                       posterior_samples) {

  terms <- .data_marginalized_random_effects(object[["data"]])
  S     <- nrow(posterior_samples)
  K     <- nrow(data[["outcome"]])
  if (length(terms) == 0L) {
    return(matrix(0, nrow = S, ncol = K))
  }

  source_samples <- .predict_known_v_newdata_marginalized_source_samples(
    object            = object,
    data              = data,
    posterior_samples = posterior_samples
  )
  draws    <- matrix(0, nrow = S, ncol = K)
  fitted_K <- nrow(object[["data"]][["outcome"]])

  for (term in terms) {
    sd_samples <- .marginalized_random_effect_sd_samples(
      term              = term,
      posterior_samples = posterior_samples,
      K                 = K,
      source_samples    = source_samples,
      fitted_K          = fitted_K
    )
    sd_samples <- .expand_brma_mv_heterogeneity_samples(
      samples = sd_samples,
      S       = S,
      K       = K
    )
    draws <- draws + .predict_known_v_marginalized_random_term_draws(
      term       = term,
      data       = data,
      sd_samples = sd_samples
    )
  }

  return(draws)
}


.predict_known_v_marginalized_random_term_draws <- function(term, data,
                                                            sd_samples) {

  S          <- nrow(sd_samples)
  K          <- ncol(sd_samples)
  group_keys <- .predict_known_v_marginalized_random_group_keys(term, data, K)
  if (is.null(group_keys)) {
    return(matrix(
      stats::rnorm(length(sd_samples), mean = 0, sd = as.vector(sd_samples)),
      nrow = S,
      ncol = K
    ))
  }

  group_factor <- factor(group_keys, levels = unique(group_keys))
  z <- matrix(
    stats::rnorm(S * nlevels(group_factor)),
    nrow = S,
    ncol = nlevels(group_factor)
  )
  out <- matrix(NA_real_, nrow = S, ncol = K)
  for (level in seq_len(nlevels(group_factor))) {
    rows <- which(group_factor == levels(group_factor)[[level]])
    out[, rows] <- sd_samples[, rows, drop = FALSE] *
      matrix(z[, level], nrow = S, ncol = length(rows))
  }

  return(out)
}


.predict_known_v_marginalized_random_group_keys <- function(term, data, K) {

  location <- data[["location"]]
  if (is.null(location) || nrow(location) != K) {
    return(NULL)
  }

  group_label <- term[["group_label"]]
  if (is.character(group_label) && length(group_label) == 1L &&
      !is.na(group_label) && nzchar(group_label)) {
    if (group_label %in% names(location)) {
      return(as.character(location[[group_label]]))
    }
    variables <- strsplit(group_label, ":", fixed = TRUE)[[1L]]
    if (length(variables) > 1L && all(variables %in% names(location))) {
      parts <- lapply(variables, function(variable) {
        as.character(location[[variable]])
      })
      return(do.call(paste, c(parts, sep = ":")))
    }
  }

  grouping_factor <- attr(term, "grouping_factor", exact = TRUE)
  if (is.character(grouping_factor) && length(grouping_factor) == 1L &&
      !is.na(grouping_factor) && nzchar(grouping_factor) &&
      grouping_factor %in% names(location)) {
    return(as.character(location[[grouping_factor]]))
  }

  return(NULL)
}


.predict_known_v_formula_design_with_row_source_values <- function(object, data) {

  formula_design <- .fitted_formula_design(object, "mu", required = TRUE)
  if (!.is_scale(object)) {
    return(formula_design)
  }

  values <- .predict_known_v_tau_source_values_function(
    object = object,
    data   = data
  )
  formula_design[["random_effects"]] <- lapply(
    formula_design[["random_effects"]],
    .predict_known_v_random_term_with_tau_source_values,
    values = values
  )

  return(formula_design)
}


.predict_known_v_tau_source_values_function <- function(object, data) {

  fit          <- object[["fit"]]
  model_data   <- data
  priors       <- object[["priors"]]
  source_names <- unname(.data_scale_formula_sources(model_data))

  force(fit)
  force(model_data)
  force(priors)

  values <- lapply(source_names, function(source_name) {

    force(source_name)
    function(parameters, data, n_rows) {

      parameter_values <- unlist(parameters, use.names = TRUE)
      posterior_row    <- matrix(as.numeric(parameter_values), nrow = 1L)
      colnames(posterior_row) <- names(parameter_values)

      source_samples <- .predict_known_v_scale_source_samples(
        fit               = fit,
        data              = model_data,
        priors            = priors,
        posterior_samples = posterior_row,
        source_names      = source_name
      )
      values <- source_samples[[source_name]]
      if (ncol(values) != n_rows) {
        stop(
          "Known-V row SD source '", source_name, "' evaluated to ",
          ncol(values), " row value(s), expected ", n_rows, ".",
          call. = FALSE
        )
      }

      as.numeric(values[1L, ])
    }
  })
  names(values) <- source_names

  return(values)
}


.predict_known_v_scale_source_samples <- function(fit, data, priors,
                                                  posterior_samples,
                                                  source_names = NULL) {

  if (!.is_data_scale(data)) {
    return(list())
  }

  available_sources <- unname(.data_scale_formula_sources(data))
  if (is.null(source_names)) {
    source_names <- available_sources
  }
  source_names <- unique(source_names)

  missing_sources <- setdiff(source_names, available_sources)
  if (length(missing_sources) > 0L) {
    stop(
      "Known-V prediction cannot evaluate random-effect row SD source(s): ",
      paste(missing_sources, collapse = ", "),
      call. = FALSE
    )
  }

  K             <- nrow(data[["outcome"]])
  scale_samples <- .evaluate.brma.scale_terms(
    fit               = fit,
    data              = data,
    priors            = priors,
    posterior_samples = posterior_samples,
    as_list           = FALSE
  )

  out <- lapply(source_names, function(source_name) {
    columns <- paste0(source_name, "[", seq_len(K), "]")
    if (!all(columns %in% colnames(scale_samples))) {
      stop(
        "Known-V prediction cannot evaluate random-effect row SD source '",
        source_name, "'.",
        call. = FALSE
      )
    }

    unname(scale_samples[, columns, drop = FALSE])
  })
  names(out) <- source_names

  return(out)
}


.predict_known_v_random_term_with_tau_source_values <- function(term, values) {

  binding <- term[["sd_binding"]]
  if (is.null(binding)) {
    return(term)
  }

  binding <- .predict_known_v_binding_with_tau_source_values(
    binding = binding,
    values  = values
  )
  term[["sd_binding"]] <- binding

  return(term)
}


.predict_known_v_binding_with_tau_source_values <- function(binding, values) {

  if (!is.null(binding[["source"]])) {
    binding[["source"]] <- .predict_known_v_source_with_tau_values(
      source = binding[["source"]],
      values = values
    )
  }

  if (is.list(binding[["sources_by_column"]])) {
    binding[["sources_by_column"]] <- lapply(
      binding[["sources_by_column"]],
      .predict_known_v_source_with_tau_values,
      values = values
    )
  }

  if (is.list(binding[["allocations"]])) {
    binding[["allocations"]] <- lapply(binding[["allocations"]], function(allocation) {
      if (!is.null(allocation[["source"]])) {
        allocation[["source"]] <- .predict_known_v_source_with_tau_values(
          source = allocation[["source"]],
          values = values
        )
      }
      allocation
    })
  }

  return(binding)
}


.predict_known_v_source_with_tau_values <- function(source, values) {

  source_name <- source[["name"]]
  if (is.null(source) ||
      !is.character(source_name) || length(source_name) != 1L ||
      is.na(source_name) || !nzchar(source_name) ||
      !identical(source[["shape"]], "row")) {
    return(source)
  }
  if (!source_name %in% names(values)) {
    return(source)
  }

  BayesTools::random_sd_source(
    BayesTools::parameter_source(
      name   = source_name,
      shape  = "row",
      values = values[[source_name]]
    )
  )
}


.predict_known_v_newdata_marginalized_source_samples <- function(object, data,
                                                                 posterior_samples) {

  if (!.is_scale(object)) {
    return(NULL)
  }

  source_names <- .predict_known_v_marginalized_row_source_names(object)
  if (length(source_names) == 0L) {
    return(NULL)
  }
  .predict_known_v_scale_source_samples(
    fit               = object[["fit"]],
    data              = data,
    priors            = object[["priors"]],
    posterior_samples = posterior_samples,
    source_names      = source_names
  )
}


.predict_known_v_marginalized_row_source_names <- function(object) {

  terms <- .data_marginalized_random_effects(object[["data"]])
  if (length(terms) == 0L) {
    return(character())
  }

  out <- character()
  for (term in terms) {
    sources <- .predict_known_v_marginalized_sd_sources(term)
    for (source in sources) {
      if (identical(source[["shape"]], "row")) {
        out <- c(out, source[["name"]])
      }
    }
  }

  unique(out[!is.na(out) & nzchar(out)])
}


.predict_known_v_marginalized_sd_sources <- function(term) {

  binding <- term[["sd_binding"]]
  if (is.null(binding)) {
    return(list())
  }
  if (.marginalized_random_effect_has_allocation(term)) {
    allocation <- binding[["allocations"]][[1L]]
    return(list(allocation[["source"]]))
  }
  if (!is.null(binding[["source"]])) {
    return(list(binding[["source"]]))
  }
  sources <- binding[["sources_by_column"]]
  if (is.list(sources)) {
    return(Filter(Negate(is.null), sources))
  }

  list()
}


.predict_brma_attach_mv_metadata <- function(samples, object, type, same_data,
                                             aggregate, random_mv,
                                             known_V_new = NULL) {

  if (!inherits(object, "brma.mv")) {
    return(samples)
  }

  metadata <- .predict_brma_mv_metadata(
    object      = object,
    type        = type,
    same_data   = same_data,
    aggregate   = aggregate,
    random_mv   = random_mv,
    known_V_new = known_V_new
  )

  if (is.list(samples) && !is.matrix(samples)) {
    samples <- lapply(samples, function(component_samples) {
      attr(component_samples, "brma_mv_prediction_target") <- metadata
      component_samples
    })
  }
  attr(samples, "brma_mv_prediction_target") <- metadata

  return(samples)
}


.predict_brma_mv_metadata <- function(object, type, same_data, aggregate,
                                      random_mv, known_V_new = NULL) {

  formula_target <- switch(type,
    "terms"       = "fixed",
    "terms.scale" = "scale",
    "cluster"     = "legacy_cluster",
    "estimate"    = if (random_mv) {
      if (aggregate) "marginal" else "conditional"
    } else {
      "fixed_plus_heterogeneity"
    },
    "response"    = if (random_mv) {
      "fixed"
    } else {
      "fixed_plus_heterogeneity"
    },
    NA_character_
  )

  random_effect_target    <- "none"
  random_effects_in_mean <- FALSE
  mean_target             <- if (type %in% c("terms", "response")) {
    "fixed_location"
  } else {
    formula_target
  }
  if (random_mv) {
    if (type == "estimate" && is.null(known_V_new)) {
      random_effect_target    <- if (aggregate) "marginal_sample" else "conditional_mean"
      random_effects_in_mean <- TRUE
      mean_target             <- if (aggregate) {
        "fixed_plus_marginal_random_effects"
      } else {
        "fixed_plus_conditional_random_effects"
      }
    } else if (type == "response") {
      random_effect_target <- "marginal_covariance"
    } else if (type == "terms") {
      random_effect_target <- "excluded"
    }
  }

  response_covariance_target <- NA_character_
  if (type == "response") {
    response_covariance_target <- if (!is.null(known_V_new) && random_mv) {
      "V_new_plus_marginal_random_effect_covariance"
    } else if (!is.null(known_V_new)) {
      "V_new_plus_heterogeneity"
    } else if (random_mv && .is_data_known_v(object[["data"]]) && same_data) {
      "known_V_plus_marginal_random_effect_covariance"
    } else if (random_mv && .data_has_marginalized_random_effects(object[["data"]])) {
      "known_V_plus_marginalized_estimate_level_variance"
    } else if (.is_data_known_v(object[["data"]])) {
      "known_V_plus_heterogeneity"
    } else {
      "sampling_variance_plus_heterogeneity"
    }
  }

  new_levels <- NA_character_
  if (random_mv && type == "response" && !is.null(known_V_new)) {
    new_levels <- "marginal_random_effect_covariance"
  } else if (random_mv && !same_data) {
    new_levels <- if (aggregate) "sample_marginally" else "sample_unseen_levels"
  } else if (random_mv && same_data) {
    new_levels <- "existing_levels_only"
  }

  return(list(
    method                     = "predict.brma",
    class                      = "brma.mv",
    type                       = type,
    unit                       = if (aggregate) "aggregate" else "estimate",
    same_data                  = same_data,
    known_v                    = .is_data_known_v(object[["data"]]) || !is.null(known_V_new),
    random_formula             = random_mv,
    v_new                      = !is.null(known_V_new),
    mean_target                = mean_target,
    covariance_target          = response_covariance_target,
    formula_target             = formula_target,
    random_effect_target       = random_effect_target,
    random_effects_in_mean     = random_effects_in_mean,
    response_covariance_target = response_covariance_target,
    V_new                      = !is.null(known_V_new),
    observed_cross_covariance  = if (!is.null(known_V_new)) FALSE else NA,
    new_levels                 = new_levels
  ))
}
