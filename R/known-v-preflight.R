# Singular known-V preflight -----

.brma_mv_check_singular_v_regularization <- function(object) {

  data <- object[["data"]]
  if (!.is_data_known_v(data)) {
    return(invisible(TRUE))
  }

  known_V <- .data_known_v_data(data)
  if (!.known_v_is_singular_representation(known_V) &&
      !(!is.null(known_V[["V"]]) && .known_v_is_singular(known_V[["V"]]))) {
    return(invisible(TRUE))
  }

  K                 <- .known_v_nrow(known_V)
  regularized_rows  <- .brma_mv_regularized_variance_rows(object, K = K)
  covariance_blocks <- .known_v_correlated_blocks(known_V)
  invalid_blocks    <- Filter(function(block) {
    covariance <- block[["covariance"]]
    index      <- block[["index"]]
    .known_v_is_singular(covariance) &&
      !.known_v_nullspace_is_regularized(
        covariance       = covariance,
        regularized_rows = regularized_rows[index]
      )
  }, covariance_blocks)

  if (length(invalid_blocks) > 0L) {
    block_labels <- .known_v_block_labels(invalid_blocks)
    stop(
      "Singular 'V' dependency block(s) at retained rows ",
      paste(block_labels, collapse = ", "),
      " are not structurally regularized by integrated conditional variance. ",
      "Sampled random effects change the conditional mean only and do not ",
      "regularize 'V'. Specify strictly positive heterogeneity or a marginalized ",
      "random-effect variance that covers every singular direction, or supply ",
      "a positive-definite 'V'.",
      call. = FALSE
    )
  }

  fixed_variance <- .brma_mv_fixed_integrated_variance(object, K = K)
  if (identical(.known_v_effective_backend(known_V), "block_mvn") &&
      !is.null(fixed_variance)) {
    invalid_numeric_blocks <- Filter(function(block) {
      index      <- block[["index"]]
      covariance <- block[["covariance"]] +
        diag(fixed_variance[index], nrow = length(index))
      !.known_v_is_numerically_positive_definite(covariance)
    }, covariance_blocks)
    if (length(invalid_numeric_blocks) > 0L) {
      block_labels <- .known_v_block_labels(invalid_numeric_blocks)
      stop(
        "Fixed integrated variance does not make singular 'V' dependency ",
        "block(s) at retained rows ", paste(block_labels, collapse = ", "),
        " numerically positive definite for the block-MVN backend. Increase ",
        "the fixed heterogeneity/random-effect SD or supply a ",
        "positive-definite 'V'.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


# Explain when non-focal variance leaves a tau-zero known-V ordinate singular.
.iwmde_known_v_tau_zero_boundary_reason <- function(
    context, parameter, parameter_spec, values) {

  data <- context[["data"]]
  if (!identical(parameter, "tau") ||
      !identical(parameter_spec[["type"]], "primitive") ||
      !any(values == 0, na.rm = TRUE) ||
      !.is_data_known_v(data)) {
    return(NULL)
  }

  known_V <- .data_known_v_data(data)
  singular_blocks <- Filter(function(block) {
    .known_v_is_singular(block[["covariance"]])
  }, .known_v_correlated_blocks(known_V))
  if (length(singular_blocks) == 0L) {
    return(NULL)
  }

  K <- .known_v_nrow(known_V)
  regularized_rows <- if (.is_data_random(data)) {
    .brma_mv_regularized_variance_rows(context[["object"]], K = K)
  } else {
    rep(FALSE, K)
  }
  invalid_blocks <- Filter(function(block) {
    index <- block[["index"]]
    !.known_v_nullspace_is_regularized(
      covariance       = block[["covariance"]],
      regularized_rows = regularized_rows[index]
    )
  }, singular_blocks)
  if (length(invalid_blocks) == 0L) {
    return(NULL)
  }

  return(paste0(
    "the tau = 0 likelihood is degenerate for singular known-V dependency ",
    "block(s) at retained rows ",
    paste(.known_v_block_labels(invalid_blocks), collapse = ", "),
    " because no non-focal integrated variance regularizes every null ",
    "direction"
  ))
}


# Format retained-row labels for known-V dependency blocks.
.known_v_block_labels <- function(blocks) {

  vapply(blocks, function(block) {
    paste0("{", paste(block[["index"]], collapse = ", "), "}")
  }, character(1))
}


# Return fixed integrated scalar variance when it is known before fitting.
.brma_mv_fixed_integrated_variance <- function(object, K) {

  data               <- object[["data"]]
  fixed_scale_values <- if (.is_data_scale(data)) {
    .brma_mv_fixed_scale_values(object, K = K)
  } else {
    list()
  }

  if (.is_data_random(data)) {
    terms <- .data_marginalized_random_effects(data)
    if (length(terms) == 0L) {
      return(NULL)
    }

    design     <- .fitted_formula_design(object, "mu", required = TRUE)
    prior_list <- design[["prior_list"]]
    variance   <- numeric(K)
    for (term in terms) {
      term_variance <- .marginalized_random_term_fixed_variance(
        term               = term,
        data               = data,
        prior_list         = prior_list,
        fixed_scale_values = fixed_scale_values,
        K                  = K
      )
      if (is.null(term_variance)) {
        return(NULL)
      }
      variance <- variance + term_variance
    }
    return(variance)
  }

  if (.is_data_scale(data)) {
    scale_specs <- .data_scale_component_specs(data)
    if (length(scale_specs) != 1L) {
      return(NULL)
    }
    source <- scale_specs[[1L]][["source"]]
    values <- fixed_scale_values[[source]]
    if (is.null(values)) {
      return(NULL)
    }
    return(values^2)
  }

  tau_prior <- object[["priors"]][["outcome"]][["tau"]]
  if (is.null(tau_prior) || BayesTools::is.prior.mixture(tau_prior) ||
      !BayesTools::is.prior.point(tau_prior)) {
    return(NULL)
  }

  location <- tau_prior[["parameters"]][["location"]]
  if (length(location) == 1L) {
    location <- rep(location, K)
  }
  if (length(location) != K || anyNA(location) ||
      any(!is.finite(location)) || any(location < 0)) {
    return(NULL)
  }

  location^2
}


# Return fitted-row SD values for each fully point-fixed scale formula.
.brma_mv_fixed_scale_values <- function(object, K) {

  scale_specs <- .data_scale_component_specs(object[["data"]])
  values      <- lapply(scale_specs, function(scale_spec) {
    design <- .fitted_formula_design(
      object    = object,
      parameter = scale_spec[["parameter"]],
      required  = TRUE
    )
    log_scale <- .formula_design_fixed_values(design)
    if (is.null(log_scale) || length(log_scale) != K) {
      return(NULL)
    }

    scale <- exp(log_scale)
    if (length(scale) != K || anyNA(scale) || any(!is.finite(scale)) ||
        any(scale < 0)) {
      return(NULL)
    }

    as.numeric(scale)
  })
  names(values) <- vapply(scale_specs, `[[`, character(1), "source")

  values
}


# Evaluate a fully point-fixed formula design on its fitted rows.
.formula_design_fixed_values <- function(design) {

  model_matrix <- design[["model_matrix"]]
  prior_list   <- design[["prior_list"]]
  parameter    <- design[["parameter"]]
  if (!is.matrix(model_matrix) || anyNA(model_matrix) ||
      any(!is.finite(model_matrix)) || is.null(prior_list) ||
      is.null(parameter)) {
    return(NULL)
  }

  output         <- numeric(nrow(model_matrix))
  intercept_name <- paste0(parameter, "_intercept")
  if (intercept_name %in% names(prior_list)) {
    prior      <- prior_list[[intercept_name]]
    value      <- .prior_fixed_values(prior)
    multiplier <- .prior_fixed_multiplier(prior)
    if (is.null(value) || length(value) != 1L || is.null(multiplier)) {
      return(NULL)
    }
    if (isTRUE(design[["log_intercept"]])) {
      if (!is.finite(value) || value <= 0) {
        stop(
          "Point-fixed scale intercepts must be strictly positive.",
          call. = FALSE
        )
      }
      value <- log(value)
    }
    output <- output + multiplier * value
  }

  term_names <- setdiff(names(prior_list), intercept_name)
  for (term_name in term_names) {
    prior      <- prior_list[[term_name]]
    values     <- .prior_fixed_values(prior)
    multiplier <- .prior_fixed_multiplier(prior)
    model_term <- sub(paste0("^", parameter, "_"), "", term_name)
    term_index <- match(model_term, design[["model_terms"]])
    if (is.null(values) || is.null(multiplier) || is.na(term_index)) {
      return(NULL)
    }
    columns <- which(design[["assign"]] == (term_index - 1L))
    if (length(columns) == 0L) {
      return(NULL)
    }
    if (BayesTools::is.prior.factor(prior)) {
      if (length(values) == 1L) {
        values <- rep(values, length(columns))
      }
      if (length(values) != length(columns)) {
        return(NULL)
      }
      contribution <- drop(model_matrix[, columns, drop = FALSE] %*% values)
    } else {
      if (length(values) != 1L || length(columns) != 1L) {
        return(NULL)
      }
      contribution <- model_matrix[, columns] * values
    }
    output <- output + multiplier * contribution
  }

  as.numeric(output)
}


# Return a prior's fixed numeric multiplier, if it has one.
.prior_fixed_multiplier <- function(prior) {

  multiplier <- attr(prior, "multiply_by", exact = TRUE)
  if (is.null(multiplier)) {
    return(1)
  }
  if (!is.numeric(multiplier) || length(multiplier) != 1L ||
      is.na(multiplier) || !is.finite(multiplier)) {
    return(NULL)
  }

  as.numeric(multiplier)
}


# Return a marginalized term's row-wise variance when all sources are fixed.
.marginalized_random_term_fixed_variance <- function(
    term, data, prior_list, fixed_scale_values = list(), K) {

  parameter      <- term[["sd_parameter_names"]]
  parameter      <- parameter[!is.na(parameter) & nzchar(parameter)]
  has_allocation <- .marginalized_random_effect_has_allocation(term)
  if (!has_allocation && length(parameter) == 1L) {
    fixed_sd <- .prior_fixed_values(prior_list[[parameter]])
  } else if (has_allocation) {
    binding    <- term[["sd_binding"]]
    allocation <- binding[["allocations"]][[1L]]
    fixed_sd <- .random_sd_source_fixed_values(
      source             = allocation[["source"]],
      data               = data,
      prior_list         = prior_list,
      fixed_scale_values = fixed_scale_values,
      K                  = K
    )
    if (is.null(fixed_sd)) {
      return(NULL)
    }
    for (factor in .marginalized_random_effect_allocation_factors(term)) {
      weights <- .prior_fixed_values(prior_list[[factor[["weight_name"]]]])
      index   <- factor[["index"]]
      if (is.null(weights) || length(weights) < index) {
        return(NULL)
      }
      multiplier <- weights[[index]]
      if (identical(factor[["scale"]], "mean_variance")) {
        multiplier <- factor[["n_targets"]] * multiplier
      } else if (!identical(factor[["scale"]], "total_variance")) {
        return(NULL)
      }
      if (!is.finite(multiplier) || multiplier < 0) {
        return(NULL)
      }
      fixed_sd <- fixed_sd * sqrt(multiplier)
    }
  } else {
    binding <- term[["sd_binding"]]
    source  <- if (is.null(binding)) NULL else binding[["source"]]
    if (is.null(source) && !is.null(binding) &&
        length(binding[["sources_by_column"]]) == 1L) {
      source <- binding[["sources_by_column"]][[1L]]
    }
    fixed_sd <- .random_sd_source_fixed_values(
      source             = source,
      data               = data,
      prior_list         = prior_list,
      fixed_scale_values = fixed_scale_values,
      K                  = K
    )
  }

  if (is.null(fixed_sd) || !length(fixed_sd) %in% c(1L, K)) {
    return(NULL)
  }
  variance <- .marginalized_random_effect_variance_samples(
    term       = term,
    sd_samples = matrix(fixed_sd, nrow = 1L),
    K          = K
  )

  as.numeric(variance[1L, ])
}


# Return fixed SD-source values, or NULL for a sampled/formula source.
.random_sd_source_fixed_values <- function(
    source, data, prior_list, fixed_scale_values = list(), K) {

  if (is.null(source)) {
    return(NULL)
  }

  values <- source[["values"]]
  if (is.null(values)) {
    name <- source[["name"]]
    if (!is.null(name) && name %in% .data_scale_formula_sources(data)) {
      return(fixed_scale_values[[name]])
    }
    if (!is.null(name) && name %in% names(prior_list)) {
      values <- .prior_fixed_values(prior_list[[name]])
    }
  }
  if (is.null(values)) {
    nested_source <- source[["source"]]
    if (!is.null(nested_source) && !identical(nested_source, source)) {
      return(.random_sd_source_fixed_values(
        source             = nested_source,
        data               = data,
        prior_list         = prior_list,
        fixed_scale_values = fixed_scale_values,
        K                  = K
      ))
    }
    return(NULL)
  }

  shape <- source[["shape"]]
  if (identical(shape, "scalar") && length(values) != 1L) {
    return(NULL)
  }
  if (identical(shape, "row") && length(values) == 1L) {
    values <- rep(values, K)
  }
  if (!identical(shape, "scalar") && !identical(shape, "row") &&
      length(values) != 1L) {
    return(NULL)
  }
  if (!length(values) %in% c(1L, K) || anyNA(values) ||
      any(!is.finite(values)) || any(values < 0)) {
    return(NULL)
  }

  as.numeric(values)
}


# Return a point prior's fixed values, including identical point mixtures.
.prior_fixed_values <- function(prior) {

  if (is.null(prior)) {
    return(NULL)
  }
  if (BayesTools::is.prior.mixture(prior)) {
    values <- lapply(prior, .prior_fixed_values)
    if (any(vapply(values, is.null, logical(1))) ||
        !all(vapply(values[-1L], identical, logical(1), values[[1L]]))) {
      return(NULL)
    }
    return(values[[1L]])
  }
  if (!BayesTools::is.prior.point(prior)) {
    return(NULL)
  }

  values <- prior[["parameters"]][["location"]]
  if (!is.numeric(values) || length(values) == 0L || anyNA(values) ||
      any(!is.finite(values))) {
    return(NULL)
  }

  as.numeric(values)
}


# Check positive definiteness at the covariance's relative numerical scale.
.known_v_is_numerically_positive_definite <- function(covariance) {

  .covariance_is_numerically_positive_definite(
    .covariance_factorization(covariance)
  )
}


.brma_mv_regularized_variance_rows <- function(object, K) {

  data <- object[["data"]]
  if (!.is_data_random(data)) {
    if (.is_data_scale(data)) {
      return(rep(TRUE, K))
    }

    tau_prior <- object[["priors"]][["outcome"]][["tau"]]
    return(rep(.prior_coordinate_is_structurally_positive(tau_prior), K))
  }

  terms <- .data_marginalized_random_effects(data)
  if (length(terms) == 0L) {
    return(rep(FALSE, K))
  }

  design     <- .fitted_formula_design(object, "mu", required = TRUE)
  prior_list <- design[["prior_list"]]
  support    <- rep(FALSE, K)
  for (term in terms) {
    support <- support | .marginalized_random_term_positive_rows(
      term       = term,
      data       = data,
      prior_list = prior_list,
      K          = K
    )
  }

  support
}


.marginalized_random_term_positive_rows <- function(term, data, prior_list, K) {

  binding       <- term[["sd_binding"]]
  has_allocation <- .marginalized_random_effect_has_allocation(term)
  parameter <- term[["sd_parameter_names"]]
  parameter <- parameter[!is.na(parameter) & nzchar(parameter)]
  if (!has_allocation && length(parameter) == 1L) {
    support <- rep(
      .prior_coordinate_is_structurally_positive(prior_list[[parameter]]),
      K
    )
  } else {
    source  <- if (is.null(binding)) NULL else binding[["source"]]
    if (is.null(source) && !is.null(binding) &&
        length(binding[["sources_by_column"]]) == 1L) {
      source <- binding[["sources_by_column"]][[1L]]
    }
    support <- .random_sd_source_positive_rows(
      source     = source,
      data       = data,
      prior_list = prior_list,
      K          = K
    )

    if (has_allocation) {
      factors <- .marginalized_random_effect_allocation_factors(term)
      for (factor in factors) {
        factor_prior <- prior_list[[factor[["weight_name"]]]]
        factor_positive <- .prior_coordinate_is_structurally_positive(
          prior = factor_prior,
          index = factor[["index"]]
        )
        support <- support & factor_positive
      }
    }
  }

  multiplier <- .marginalized_random_effect_row_multiplier(term, K = K)
  if (!is.null(multiplier)) {
    support <- support & multiplier > 0
  }

  support
}


.random_sd_source_positive_rows <- function(source, data, prior_list, K) {

  if (is.null(source)) {
    return(rep(FALSE, K))
  }

  values <- source[["values"]]
  if (!is.null(values)) {
    if (length(values) == 1L) {
      values <- rep(values, K)
    }
    if (length(values) != K || anyNA(values) || any(!is.finite(values))) {
      stop("Invalid fixed random-effect SD source metadata.", call. = FALSE)
    }
    return(values > 0)
  }

  name <- source[["name"]]
  if (!is.null(name) && name %in% .data_scale_formula_sources(data)) {
    return(rep(TRUE, K))
  }
  if (!is.null(name) && name %in% names(prior_list)) {
    return(rep(
      .prior_coordinate_is_structurally_positive(prior_list[[name]]),
      K
    ))
  }

  nested_source <- source[["source"]]
  if (!is.null(nested_source) && !identical(nested_source, source)) {
    return(.random_sd_source_positive_rows(
      source     = nested_source,
      data       = data,
      prior_list = prior_list,
      K          = K
    ))
  }

  stop(
    "Cannot verify singular-V regularization from random-effect SD source '",
    if (is.null(name)) "<unnamed>" else name,
    "'.",
    call. = FALSE
  )
}


.prior_coordinate_is_structurally_positive <- function(prior, index = 1L) {

  if (is.null(prior)) {
    return(FALSE)
  }
  if (BayesTools::is.prior.mixture(prior)) {
    return(all(vapply(
      prior,
      .prior_coordinate_is_structurally_positive,
      logical(1),
      index = index
    )))
  }
  if (!BayesTools::is.prior.point(prior)) {
    return(TRUE)
  }

  location <- prior[["parameters"]][["location"]]
  if (length(location) == 1L) {
    location <- rep(location, index)
  }

  length(location) >= index && is.finite(location[[index]]) &&
    location[[index]] > 0
}


.known_v_nullspace_is_regularized <- function(covariance, regularized_rows) {

  if (!is.logical(regularized_rows) ||
      length(regularized_rows) != nrow(covariance) ||
      anyNA(regularized_rows)) {
    stop("Internal error: invalid singular-V regularization rows.",
         call. = FALSE)
  }

  factorization <- .covariance_factorization(covariance)
  null <- factorization[["eigenvectors"]][
    , factorization[["decomposition_values"]] <= factorization[["psd_tolerance"]],
    drop = FALSE
  ]
  if (ncol(null) == 0L) {
    return(TRUE)
  }
  if (!any(regularized_rows)) {
    return(FALSE)
  }

  qr(null[regularized_rows, , drop = FALSE], tol = sqrt(.Machine$double.eps))[["rank"]] ==
    ncol(null)
}
