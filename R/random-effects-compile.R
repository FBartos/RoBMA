# ============================================================================ #
# random-effects-compile.R
# ============================================================================ #
#
# RoBMA-side selection and likelihood variance construction for BayesTools
# marginalized random-effect compilation.
#
# ============================================================================ #


.prepare_brma_mv_random_effects_compile <- function(object,
                                                    marginalize_estimate_level) {

  BayesTools::check_bool(
    marginalize_estimate_level,
    "marginalize_estimate_level",
    allow_NA = FALSE
  )

  if (!.is_data_random(object[["data"]]) || !isTRUE(marginalize_estimate_level)) {
    return(object)
  }

  sampled_design <- .object_bayestools_formula_design(
    object                 = object,
    parameter              = "mu",
    source                 = "location",
    random_effects_compile = NULL
  )
  candidates <- .marginalized_random_effect_candidates(
    formula_design = sampled_design,
    data           = object[["data"]]
  )

  if (length(candidates) == 0L) {
    return(object)
  }
  if (length(candidates) > 1L) {
    stop(
      "Multiple random-effect blocks map one-to-one to estimates and could ",
      "be marginalized: ",
      paste(vapply(candidates, `[[`, character(1), "block_name"), collapse = ", "),
      ". Specify a single estimate-level random term or set ",
      "'marginalize_estimate_level = FALSE'.",
      call. = FALSE
    )
  }

  marginalized_block <- candidates[[1L]][["block_name"]]
  compile            <- BayesTools::random_effects_compile(
    marginalized = marginalized_block
  )
  compiled_design <- .object_bayestools_formula_design(
    object                 = object,
    parameter              = "mu",
    source                 = "location",
    random_effects_compile = compile
  )
  marginalized_terms <- .formula_design_random_effects_by_mode(
    compiled_design,
    "marginalized"
  )
  known_r_metadata <- .attach_known_r_marginal_variance_factors(
    formula_design = compiled_design,
    terms          = marginalized_terms
  )
  compiled_design     <- known_r_metadata[["formula_design"]]
  marginalized_terms  <- known_r_metadata[["terms"]]
  object[["data"]] <- .known_v_auto_update_for_marginalized_random(
    data  = object[["data"]],
    terms = marginalized_terms
  )
  .validate_marginalized_random_effects(
    terms = marginalized_terms,
    data  = object[["data"]]
  )

  object[["data"]] <- .set_data_random_effects_compile(
    data                 = object[["data"]],
    compile              = compile,
    marginalized_effects = marginalized_terms
  )

  formula_design <- object[["formula_design"]]
  if (is.null(formula_design)) {
    formula_design <- list()
  }
  formula_design[["mu"]]     <- compiled_design
  object[["formula_design"]] <- formula_design

  object[["random_effects_compile"]] <- list(mu = compile)

  return(object)
}


.known_v_auto_update_for_marginalized_random <- function(data, terms) {

  if (!.is_data_known_v(data) || length(terms) == 0L ||
      !.marginalized_random_effects_row_varying(terms)) {
    return(data)
  }

  known_V <- .data_known_v_data(data)
  if (!identical(known_V[["parameterization_requested"]], "auto") ||
      !identical(known_V[["parameterization"]], "whitened")) {
    return(data)
  }

  known_v_residual_fraction_specified <- !is.null(
    known_V[["residual_fraction_requested"]]
  )
  new_known_V <- .known_v_prepare(
    V                                   = known_V[["V"]],
    keep_rows                           = rep(TRUE, nrow(known_V[["V"]])),
    known_v_parameterization            = "auto",
    known_v_residual_fraction           = known_V[["residual_fraction_requested"]],
    known_v_residual_fraction_specified = known_v_residual_fraction_specified,
    known_v_is_scale                    = TRUE,
    known_v_is_multilevel               = FALSE
  )
  attr(data, "known_V_data") <- new_known_V
  data[["outcome"]][["sei"]] <- sqrt(diag(new_known_V[["V"]]))

  return(data)
}


.marginalized_random_effect_candidates <- function(formula_design, data) {

  random_effects <- formula_design[["random_effects"]]
  if (length(random_effects) == 0L) {
    return(list())
  }

  k <- nrow(data[["outcome"]])
  random_effects[vapply(
    random_effects,
    .is_marginalized_estimate_level_candidate,
    logical(1),
    k              = k,
    formula_design = formula_design
  )]
}


.is_marginalized_estimate_level_candidate <- function(term, k,
                                                      formula_design = NULL) {

  if (!identical(term[["compile_mode"]], "sampled")) {
    return(FALSE)
  }
  if (.random_effect_term_has_known_group_covariance(term)) {
    return(.known_r_marginal_variance_factors_available(
      formula_design = formula_design,
      block          = term[["block_name"]]
    ))
  }

  .is_estimate_level_random_intercept(term, k)
}

.random_effect_term_has_known_group_covariance <- function(term) {

  group_covariance <- term[["group_covariance"]]
  if (is.null(group_covariance)) {
    return(FALSE)
  }

  inherits(group_covariance, "random_group_covariance") ||
    inherits(group_covariance, "random_group_covariance_kernel") ||
    (is.list(group_covariance) && identical(group_covariance[["type"]], "known"))
}


.known_r_marginal_variance_factors_available <- function(formula_design,
                                                         block) {

  if (is.null(formula_design)) {
    return(FALSE)
  }

  .brma_mv_require_marginal_variance_factors_api()
  factors <- tryCatch(
    BayesTools::random_effects_marginal_variance_factors(
      formula_design      = formula_design,
      blocks              = block,
      require_diagonal    = TRUE,
      require_one_to_one  = TRUE
    ),
    error = function(e) NULL
  )

  !is.null(factors)
}


.brma_mv_require_marginal_variance_factors_api <- function() {

  if (!exists(
    "random_effects_marginal_variance_factors",
    envir   = asNamespace("BayesTools"),
    mode    = "function",
    inherits = FALSE
  )) {
    stop(
      "Known-R marginalization requires BayesTools with ",
      "random_effects_marginal_variance_factors().",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.attach_known_r_marginal_variance_factors <- function(formula_design, terms) {

  known_r_blocks <- vapply(
    terms,
    function(term) {
      if (.random_effect_term_has_known_group_covariance(term)) {
        term[["block_name"]]
      } else {
        NA_character_
      }
    },
    character(1)
  )
  known_r_blocks <- known_r_blocks[!is.na(known_r_blocks)]

  if (length(known_r_blocks) == 0L) {
    return(list(formula_design = formula_design, terms = terms))
  }

  .brma_mv_require_marginal_variance_factors_api()
  factors <- BayesTools::random_effects_marginal_variance_factors(
    formula_design      = formula_design,
    blocks              = known_r_blocks,
    require_diagonal    = TRUE,
    require_one_to_one  = TRUE
  )

  for (i in seq_along(terms)) {
    block_name <- terms[[i]][["block_name"]]
    if (!block_name %in% known_r_blocks) {
      next
    }
    terms[[i]] <- .attach_known_r_marginal_variance_factor(
      term   = terms[[i]],
      factor = factors[["blocks"]][[block_name]]
    )
  }

  formula_design[["random_effects"]] <- .replace_formula_design_random_terms(
    random_effects = formula_design[["random_effects"]],
    replacements   = terms
  )

  list(formula_design = formula_design, terms = terms)
}


.attach_known_r_marginal_variance_factor <- function(term, factor) {

  if (is.null(factor)) {
    stop(
      "Known-R marginal variance metadata is missing for random-effect block '",
      term[["block_name"]], "'.",
      call. = FALSE
    )
  }

  row_multiplier <- as.numeric(factor[["row_multiplier"]])
  if (length(row_multiplier) != factor[["n_rows"]] ||
      anyNA(row_multiplier) ||
      any(!is.finite(row_multiplier)) ||
      any(row_multiplier < 0)) {
    stop(
      "Known-R marginal variance multipliers are invalid for random-effect ",
      "block '", term[["block_name"]], "'.",
      call. = FALSE
    )
  }
  names(row_multiplier) <- names(factor[["row_multiplier"]])

  term[["marginal_variance_factor"]] <- factor
  term[["row_multiplier"]]           <- row_multiplier
  term[["row_multiplier_name"]]      <- .marginalized_random_effect_row_multiplier_name(
    term[["block_name"]]
  )

  term
}


.replace_formula_design_random_terms <- function(random_effects, replacements) {

  if (length(replacements) == 0L) {
    return(random_effects)
  }

  replacement_names <- vapply(replacements, `[[`, character(1), "block_name")
  for (i in seq_along(random_effects)) {
    block_name <- random_effects[[i]][["block_name"]]
    match      <- match(block_name, replacement_names)
    if (!is.na(match)) {
      random_effects[[i]] <- replacements[[match]]
    }
  }

  random_effects
}

.is_estimate_level_random_intercept <- function(term, k) {

  if (!identical(as.integer(term[["n_columns"]]), 1L)) {
    return(FALSE)
  }

  structure <- term[["structure"]]
  if (!is.character(structure) || length(structure) != 1L || is.na(structure)) {
    return(FALSE)
  }
  structure <- tolower(structure)
  if (!structure %in% c("id", "diag", "us")) {
    return(FALSE)
  }

  model_matrix <- term[["model_matrix"]]
  if (!is.matrix(model_matrix) || nrow(model_matrix) != k ||
      ncol(model_matrix) != 1L || anyNA(model_matrix) ||
      any(abs(model_matrix[, 1L] - 1) > 1e-12)) {
    return(FALSE)
  }

  group_map <- term[["group_map"]]
  if (!is.numeric(group_map) && !is.integer(group_map)) {
    return(FALSE)
  }
  if (length(group_map) != k || anyNA(group_map)) {
    return(FALSE)
  }
  group_map_int <- as.integer(group_map)
  if (any(!is.finite(group_map)) || anyNA(group_map_int) ||
      any(group_map != group_map_int)) {
    return(FALSE)
  }
  if (!identical(sort(group_map_int), seq_len(k))) {
    return(FALSE)
  }
  if (!identical(as.integer(term[["n_groups"]]), k)) {
    return(FALSE)
  }

  TRUE
}


.validate_marginalized_random_effects <- function(terms, data) {

  if (length(terms) == 0L) {
    stop(
      "Internal error: selected random-effect block was not compiled as ",
      "marginalized.",
      call. = FALSE
    )
  }
  if (length(terms) > 1L) {
    stop(
      "Internal error: brma.mv() currently supports only one marginalized ",
      "estimate-level random-effect block.",
      call. = FALSE
    )
  }

  k <- nrow(data[["outcome"]])
  for (term in terms) {
    if (!identical(term[["compile_mode"]], "marginalized") ||
        !.is_estimate_level_random_intercept(term, k = k)) {
      stop(
        "Internal error: marginalized random-effect block is not a supported ",
        "one-to-one random intercept.",
        call. = FALSE
      )
    }
    .marginalized_random_effect_variance_expression(term, row_index = "i")
  }

  if (.is_data_known_v_parameterization(data, "whitened") &&
      .marginalized_random_effects_row_varying(terms)) {
    stop(
      "known_v_parameterization = 'whitened' cannot be combined with ",
      "row-varying marginalized estimate-level random effects. Use ",
      "known_v_parameterization = 'latent' or 'block_mvn', or set ",
      "'marginalize_estimate_level = FALSE'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.set_data_random_effects_compile <- function(data, compile,
                                             marginalized_effects) {

  location <- data[["location"]]
  attr(location, "random_effects_compile")      <- compile
  attr(location, "marginalized_random_effects") <- marginalized_effects
  data[["location"]] <- location

  data
}


.data_random_effects_compile <- function(data) {

  if (is.null(data[["location"]])) {
    return(NULL)
  }

  attr(data[["location"]], "random_effects_compile", exact = TRUE)
}


.data_marginalized_random_effects <- function(data) {

  if (is.null(data[["location"]])) {
    return(list())
  }

  terms <- attr(data[["location"]], "marginalized_random_effects", exact = TRUE)
  if (is.null(terms)) {
    return(list())
  }

  terms
}


.data_has_marginalized_random_effects <- function(data) {

  length(.data_marginalized_random_effects(data)) > 0L
}


.data_has_marginalized_known_group_covariance <- function(data) {

  terms <- .data_marginalized_random_effects(data)
  if (length(terms) == 0L) {
    return(FALSE)
  }

  any(vapply(
    terms,
    .marginalized_random_effect_has_row_multiplier,
    logical(1)
  ))
}


.data_has_sampled_estimate_level_random_effects <- function(data) {

  terms <- .data_random_effect_terms(data)
  if (length(terms) == 0L) {
    return(FALSE)
  }

  k <- nrow(data[["outcome"]])
  any(vapply(
    terms,
    function(term) {
      identical(term[["compile_mode"]], "sampled") &&
        .is_estimate_level_random_intercept(term, k = k)
    },
    logical(1)
  ))
}


.data_random_effect_terms <- function(data) {

  if (!.is_data_random(data) || is.null(data[["location"]])) {
    return(list())
  }

  random_effects <- attr(data[["location"]], "random_effects", exact = TRUE)
  if (!inherits(random_effects, "BayesTools_random_effects")) {
    return(list())
  }

  terms <- random_effects[["terms"]]
  if (is.null(terms)) {
    return(list())
  }

  return(terms)
}


.random_effect_term_compile_mode <- function(term) {

  mode <- term[["compile_mode"]]
  if (is.null(mode)) {
    mode <- attr(term, "compile_mode", exact = TRUE)
  }
  if (is.null(mode) || length(mode) != 1L || is.na(mode) || !nzchar(mode)) {
    return("sampled")
  }

  return(as.character(mode))
}


.data_sampled_random_effect_blocks <- function(data) {

  terms <- .data_random_effect_terms(data)
  if (length(terms) == 0L) {
    return(character())
  }

  sampled <- vapply(
    terms,
    function(term) identical(.random_effect_term_compile_mode(term), "sampled"),
    logical(1)
  )
  if (!any(sampled)) {
    return(character())
  }

  return(vapply(terms[sampled], `[[`, character(1), "block_name"))
}


.formula_design_sampled_random_effect_blocks <- function(formula_design) {

  terms <- formula_design[["random_effects"]]
  if (length(terms) == 0L) {
    return(character())
  }

  sampled <- vapply(
    terms,
    function(term) identical(.random_effect_term_compile_mode(term), "sampled"),
    logical(1)
  )
  if (!any(sampled)) {
    return(character())
  }

  return(vapply(terms[sampled], `[[`, character(1), "block_name"))
}


.formula_design_blocks_have_known_group_covariance <- function(formula_design,
                                                               blocks = NULL) {

  terms <- formula_design[["random_effects"]]
  if (length(terms) == 0L) {
    return(FALSE)
  }

  if (!is.null(blocks)) {
    block_names <- vapply(terms, `[[`, character(1), "block_name")
    terms <- terms[block_names %in% blocks]
  }
  if (length(terms) == 0L) {
    return(FALSE)
  }

  any(vapply(terms, .random_effect_term_has_known_group_covariance, logical(1)))
}


.object_formula_random_effects_compile <- function(object, source) {

  if (source != "location" || !.is_data_random(object[["data"]])) {
    return(NULL)
  }

  compile <- .data_random_effects_compile(object[["data"]])
  if (!is.null(compile)) {
    return(compile)
  }

  compile_list <- object[["random_effects_compile"]]
  if (!is.null(compile_list[["mu"]])) {
    return(compile_list[["mu"]])
  }

  NULL
}


.formula_design_random_effects_by_mode <- function(formula_design, mode) {

  random_effects <- formula_design[["random_effects"]]
  if (length(random_effects) == 0L) {
    return(list())
  }

  random_effects[vapply(random_effects, function(term) {
    identical(term[["compile_mode"]], mode)
  }, logical(1))]
}


.data_marginalized_random_variance_expression <- function(data,
                                                          row_index = "i") {

  terms <- .data_marginalized_random_effects(data)
  if (length(terms) == 0L) {
    return("0")
  }

  expressions <- vapply(
    terms,
    function(term) {
      .marginalized_random_effect_variance_expression(
        term      = term,
        row_index = row_index
      )[["expression"]]
    },
    character(1)
  )

  paste(expressions, collapse = " + ")
}


.evaluate_marginalized_random_variance <- function(data, posterior_samples,
                                                   K = nrow(data[["outcome"]]),
                                                   source_samples = NULL) {

  posterior_samples <- as.matrix(posterior_samples)
  terms             <- .data_marginalized_random_effects(data)
  variance          <- matrix(0, nrow = nrow(posterior_samples), ncol = K)

  if (length(terms) == 0L) {
    return(variance)
  }

  for (term in terms) {
    sd_samples <- .marginalized_random_effect_sd_samples(
      term              = term,
      posterior_samples = posterior_samples,
      K                 = K,
      source_samples    = source_samples
    )

    term_variance <- .marginalized_random_effect_variance_samples(
      term       = term,
      sd_samples = sd_samples,
      K          = K
    )

    if (ncol(term_variance) == K) {
      variance <- variance + term_variance
    } else {
      stop(
        "Marginalized random-effect variance samples must be row-wise.",
        call. = FALSE
      )
    }
  }

  return(variance)
}


.marginalized_random_effect_variance_samples <- function(term, sd_samples, K) {

  if (ncol(sd_samples) == 1L) {
    variance <- matrix(
      sd_samples[, 1L]^2,
      nrow = nrow(sd_samples),
      ncol = K
    )
  } else if (ncol(sd_samples) == K) {
    variance <- sd_samples^2
  } else {
    stop(
      "Marginalized random-effect SD samples must be scalar or row-wise.",
      call. = FALSE
    )
  }

  multiplier <- .marginalized_random_effect_row_multiplier(term, K = K)
  if (!is.null(multiplier)) {
    variance <- variance * matrix(
      multiplier,
      nrow  = nrow(variance),
      ncol  = K,
      byrow = TRUE
    )
  }

  variance
}


.marginalized_random_effect_row_multiplier <- function(term, K) {

  if (!.marginalized_random_effect_has_row_multiplier(term)) {
    return(NULL)
  }

  multiplier <- term[["row_multiplier"]]
  if (!is.numeric(multiplier) ||
      length(multiplier) != K ||
      anyNA(multiplier) ||
      any(!is.finite(multiplier)) ||
      any(multiplier < 0)) {
    stop(
      "Known-R marginalized random-effect multipliers cannot be reused for ",
      "new or reordered data.",
      call. = FALSE
    )
  }

  multiplier
}


.marginalized_random_effects_row_varying <- function(terms) {

  any(vapply(
    terms,
    function(term) {
      .marginalized_random_effect_variance_expression(
        term      = term,
        row_index = "i"
      )[["row_varying"]]
    },
    logical(1)
  ))
}


.marginalized_random_effect_variance_expression <- function(term,
                                                            row_index = "i") {

  sd_expression <- .marginalized_random_effect_sd_expression(
    term      = term,
    row_index = row_index
  )
  multiplier_expression <- .marginalized_random_effect_row_multiplier_expression(
    term      = term,
    row_index = row_index
  )

  expression <- paste0("pow(", sd_expression[["expression"]], ",2)")
  if (!identical(multiplier_expression[["expression"]], "1")) {
    expression <- paste0(
      expression,
      " * ",
      multiplier_expression[["expression"]]
    )
  }

  return(list(
    expression  = expression,
    row_varying = sd_expression[["row_varying"]] ||
      multiplier_expression[["row_varying"]]
  ))
}


.marginalized_random_effect_row_multiplier_expression <- function(term,
                                                                  row_index) {

  if (!.marginalized_random_effect_has_row_multiplier(term)) {
    return(list(expression = "1", row_varying = FALSE))
  }

  return(list(
    expression  = paste0(term[["row_multiplier_name"]], "[", row_index, "]"),
    row_varying = .marginalized_random_effect_row_multiplier_varying(term)
  ))
}


.marginalized_random_effect_has_row_multiplier <- function(term) {

  !is.null(term[["row_multiplier"]]) &&
    !is.null(term[["row_multiplier_name"]])
}


.marginalized_random_effect_row_multiplier_varying <- function(term) {

  if (!.marginalized_random_effect_has_row_multiplier(term)) {
    return(FALSE)
  }

  multiplier <- term[["row_multiplier"]]
  max(abs(multiplier - multiplier[[1L]])) > sqrt(.Machine$double.eps) *
    max(1, max(abs(multiplier)))
}


.marginalized_random_effect_row_multiplier_name <- function(block_name) {

  paste0("mu__xRE_MVARx_", block_name)
}


.add_marginalized_random_effect_row_multiplier_data <- function(fit_data,
                                                                data) {

  terms <- .data_marginalized_random_effects(data)
  if (length(terms) == 0L) {
    return(fit_data)
  }

  K <- nrow(data[["outcome"]])
  for (term in terms) {
    multiplier <- .marginalized_random_effect_row_multiplier(term, K = K)
    if (is.null(multiplier)) {
      next
    }
    fit_data[[term[["row_multiplier_name"]]]] <- multiplier
  }

  fit_data
}


.marginalized_random_effect_sd_expression <- function(term,
                                                      row_index = "i") {

  binding   <- term[["sd_binding"]]
  parameter <- term[["sd_parameter_names"]]
  if (length(parameter) == 1L && !is.na(parameter) && nzchar(parameter)) {
    return(list(
      expression  = parameter,
      row_varying = FALSE
    ))
  }

  if (.marginalized_random_effect_has_allocation(term)) {
    allocation        <- binding[["allocations"]][[1L]]
    source_expression <- .random_sd_source_sd_expression(
      allocation[["source"]],
      row_index
    )
    factor_expression <- .random_allocation_factors_expression(
      .marginalized_random_effect_allocation_factors(term)
    )
    parts <- c(
      source_expression[["expression"]],
      if (!identical(factor_expression, "1")) factor_expression
    )
    return(list(
      expression  = paste(parts, collapse = " * "),
      row_varying = source_expression[["row_varying"]]
    ))
  }

  source <- if (is.null(binding)) NULL else binding[["source"]]
  if (!is.null(source)) {
    return(.random_sd_source_sd_expression(source, row_index))
  }

  sources_by_column <- if (is.null(binding)) NULL else binding[["sources_by_column"]]
  if (length(sources_by_column) == 1L && !is.null(sources_by_column[[1L]])) {
    return(.random_sd_source_sd_expression(
      sources_by_column[[1L]],
      row_index
    ))
  }

  stop(
    "Cannot construct likelihood variance for marginalized random-effect ",
    "block '", term[["block_name"]], "'.",
    call. = FALSE
  )
}


.marginalized_random_effect_sd_samples <- function(term, posterior_samples, K,
                                                  source_samples = NULL) {

  parameter <- term[["sd_parameter_names"]]
  if (length(parameter) == 1L && !is.na(parameter) && nzchar(parameter)) {
    return(.extract_scalar_posterior_samples(
      posterior_samples = posterior_samples,
      parameter         = parameter
    ))
  }

  binding <- term[["sd_binding"]]
  if (.marginalized_random_effect_has_allocation(term)) {
    return(.marginalized_random_effect_allocated_sd_samples(
      term              = term,
      posterior_samples = posterior_samples,
      K                 = K,
      source_samples    = source_samples
    ))
  }

  source  <- if (is.null(binding)) NULL else binding[["source"]]
  if (!is.null(source)) {
    return(.random_sd_source_samples(
      source            = source,
      posterior_samples = posterior_samples,
      K                 = K,
      source_samples    = source_samples
    ))
  }

  sources_by_column <- if (is.null(binding)) NULL else binding[["sources_by_column"]]
  if (length(sources_by_column) == 1L && !is.null(sources_by_column[[1L]])) {
    return(.random_sd_source_samples(
      source            = sources_by_column[[1L]],
      posterior_samples = posterior_samples,
      K                 = K,
      source_samples    = source_samples
    ))
  }

  stop(
    "Cannot evaluate marginalized random-effect variance for block '",
    term[["block_name"]], "'.",
    call. = FALSE
  )
}


.marginalized_random_effect_has_allocation <- function(term) {

  binding <- term[["sd_binding"]]
  !is.null(binding) &&
    isTRUE(binding[["true_allocation"]]) &&
    is.list(binding[["allocations"]]) &&
    length(binding[["allocations"]]) > 0L
}


.marginalized_random_effect_allocated_sd_samples <- function(term,
                                                             posterior_samples,
                                                             K,
                                                             source_samples = NULL) {

  binding    <- term[["sd_binding"]]
  allocation <- binding[["allocations"]][[1L]]
  base <- .random_sd_source_samples(
    source            = allocation[["source"]],
    posterior_samples = posterior_samples,
    K                 = K,
    source_samples    = source_samples
  )

  factors <- .marginalized_random_effect_allocation_factors(term)
  for (factor in factors) {
    multiplier <- .random_allocation_factor_samples(
      factor            = factor,
      posterior_samples = posterior_samples
    )
    base <- .recycle_rows(base, multiplier)
  }

  return(base)
}


.marginalized_random_effect_allocation_factors <- function(term) {

  binding    <- term[["sd_binding"]]
  allocation <- binding[["allocations"]][[1L]]
  target     <- allocation[["target"]]

  if (identical(target, "block")) {
    factors <- allocation[["factors"]]
    if (is.null(factors)) {
      factors <- binding[["factors"]]
    }
    return(factors)
  }

  if (identical(target, "sd_component")) {
    factors    <- allocation[["parent_factors"]]
    leaf_index <- allocation[["leaf_index_by_column"]]
    if (length(leaf_index) != 1L || is.na(leaf_index)) {
      stop(
        "Cannot evaluate marginalized random-effect SD-component allocation ",
        "for block '", term[["block_name"]], "'.",
        call. = FALSE
      )
    }
    leaf_factor <- list(
      weight_name = allocation[["weight_name"]],
      index       = leaf_index[[1L]],
      scale       = allocation[["scale"]],
      n_targets   = allocation[["n_targets"]]
    )
    return(c(factors, list(leaf_factor)))
  }

  stop(
    "Cannot evaluate marginalized random-effect allocation target '",
    target, "'.",
    call. = FALSE
  )
}


.random_allocation_factor_samples <- function(factor, posterior_samples) {

  weight_name <- factor[["weight_name"]]
  index       <- factor[["index"]]
  scale       <- factor[["scale"]]
  n_targets   <- factor[["n_targets"]]
  column      <- paste0(weight_name, "[", index, "]")

  if (!column %in% colnames(posterior_samples)) {
    stop(
      "Cannot evaluate marginalized random-effect allocation; missing ",
      "posterior column: ", column,
      call. = FALSE
    )
  }

  weights <- posterior_samples[, column]
  if (any(!is.finite(weights) | weights < 0)) {
    stop(
      "Cannot evaluate marginalized random-effect allocation with invalid ",
      "weights in posterior column: ", column,
      call. = FALSE
    )
  }

  if (identical(scale, "mean_variance")) {
    weights <- n_targets * weights
  } else if (!identical(scale, "total_variance")) {
    stop("Random-effect allocation factor scale is unsupported.",
         call. = FALSE)
  }

  matrix(sqrt(weights), ncol = 1L)
}


.recycle_rows <- function(values, multiplier) {

  if (nrow(values) != nrow(multiplier) || ncol(multiplier) != 1L) {
    stop("Cannot recycle allocation multiplier over marginalized SD samples.",
         call. = FALSE)
  }

  values * matrix(multiplier[, 1L], nrow = nrow(values), ncol = ncol(values))
}


.random_sd_source_samples <- function(source, posterior_samples, K,
                                      source_samples = NULL) {

  name <- source[["name"]]
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop(
      "Cannot evaluate marginalized random-effect variance from an unnamed ",
      "SD source.",
      call. = FALSE
    )
  }

  if (!is.null(source_samples) && name %in% names(source_samples)) {
    return(.validate_random_sd_source_samples(
      samples = source_samples[[name]],
      name    = name,
      S       = nrow(posterior_samples),
      K       = K
    ))
  }

  shape <- source[["shape"]]
  if (identical(shape, "row")) {
    return(.extract_posterior_matrix(
      posterior_samples = posterior_samples,
      parameter         = name,
      K                 = K
    ))
  }
  if (identical(shape, "scalar")) {
    return(.extract_scalar_posterior_samples(
      posterior_samples = posterior_samples,
      parameter         = name
    ))
  }

  stop(
    "Cannot evaluate marginalized random-effect variance from SD source '",
    name, "' with shape '", shape, "'.",
    call. = FALSE
  )
}


.extract_scalar_posterior_samples <- function(posterior_samples, parameter) {

  if (!parameter %in% colnames(posterior_samples)) {
    stop("Missing posterior column: ", parameter, call. = FALSE)
  }

  matrix(posterior_samples[, parameter], ncol = 1L)
}


.validate_random_sd_source_samples <- function(samples, name, S, K) {

  if (!is.matrix(samples) || nrow(samples) != S ||
      !ncol(samples) %in% c(1L, K) ||
      anyNA(samples) || any(!is.finite(samples))) {
    stop(
      "External SD source samples for '", name, "' must be a finite matrix ",
      "with ", S, " row(s) and either 1 or ", K, " column(s).",
      call. = FALSE
    )
  }

  samples
}


.random_allocation_factors_expression <- function(factors) {

  if (length(factors) == 0L) {
    return("1")
  }

  expressions <- vapply(
    factors,
    .random_allocation_factor_expression,
    character(1)
  )
  paste(expressions, collapse = " * ")
}


.random_allocation_factor_expression <- function(factor) {

  weight_name <- factor[["weight_name"]]
  index       <- factor[["index"]]
  scale       <- factor[["scale"]]
  n_targets   <- factor[["n_targets"]]

  if (!is.character(weight_name) || length(weight_name) != 1L ||
      is.na(weight_name) || !nzchar(weight_name) ||
      !is.numeric(index) || length(index) != 1L || is.na(index)) {
    stop("Random-effect allocation factor metadata are incomplete.",
         call. = FALSE)
  }

  if (identical(scale, "mean_variance")) {
    multiplier <- paste0(n_targets, " * ", weight_name, "[", index, "]")
  } else if (identical(scale, "total_variance")) {
    multiplier <- paste0(weight_name, "[", index, "]")
  } else {
    stop("Random-effect allocation factor scale is unsupported.",
         call. = FALSE)
  }

  paste0("sqrt(", multiplier, ")")
}


.random_sd_source_sd_expression <- function(source, row_index) {

  name <- source[["name"]]
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop(
      "Cannot construct marginalized random-effect variance from an unnamed ",
      "SD source.",
      call. = FALSE
    )
  }

  shape <- source[["shape"]]
  if (identical(shape, "row")) {
    return(list(
      expression  = paste0(name, "[", row_index, "]"),
      row_varying = TRUE
    ))
  }
  if (identical(shape, "scalar")) {
    return(list(
      expression  = name,
      row_varying = FALSE
    ))
  }

  stop(
    "Cannot construct marginalized random-effect variance from SD source '",
    name, "' with shape '", shape, "'.",
    call. = FALSE
  )
}
