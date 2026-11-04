# Ordered scalar targets use fitted primitive recipes and source coordinates.
# Projection values and measures belong to BayesTools; this file only chooses
# conditional charts and evaluates their complete conditional kernels.

.iwmde_ordered_stop <- function(parameter, cause, detail, value = NULL) {

  condition <- .iwmde_unavailable(paste0("Ordered density for '", parameter,
    "' is unavailable: ", detail, ". Inspect the fitted target and its density diagnostics."),
    "ordered_offset")
  condition[["class"]] <- c(paste0("RoBMA_ordered_", cause), condition[["class"]])
  condition[["cause"]] <- cause
  condition[["value"]] <- value
  stop(structure(list(message = condition[["reason"]], call = NULL, cause = cause,
    value = value, value_key = if (is.null(value)) NULL else .iwmde_key_number(value)),
    class = c(condition[["class"]], "error", "condition")))
}


.iwmde_ordered_refit <- function() {

  tryCatch(.stop_refit_required(
    "Ordered source draws are unavailable. Refit the model with the current RoBMA/BayesTools build."),
    RoBMA_refit_required = function(error) {
      class(error) <- c("RoBMA_ordered_draws_unavailable",
        "BayesTools_ordered_coordinates_unavailable", "BayesTools_ordered_unavailable", class(error))
      stop(error)
    })
}


.iwmde_require_target_type <- function(spec, replacement = FALSE) {

  if (is.null(spec)) return(invisible(TRUE))
  type <- spec[["type"]]
  known <- if (replacement) c("primitive", "scalar", "indexed", "fallback", "linear",
    "simplex_pair", "random_component_sd", "ordered", "ordered_source") else
    c("primitive", "linear", "simplex_pair", "random_component_sd", "ordered")
  if (is.null(type) && identical(spec[["status"]], "unsupported")) return(invisible(TRUE))
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% known) {
    .iwmde_stop_unavailable(.iwmde_unavailable(
      "Density estimation is unavailable for an unknown target specification. Inspect the fitted target metadata.",
      "unknown_target"))
  }
  invisible(TRUE)
}


.iwmde_ordered_prior_support <- function(prior) {

  if (BayesTools::is.prior.mixture(prior)) {
    intervals <- vapply(prior, .iwmde_ordered_prior_support, numeric(2))
    return(c(min(intervals[1L, ]), max(intervals[2L, ])))
  }
  if (BayesTools::is.prior.point(prior)) {
    location <- prior[["parameters"]][["location"]]
    if (!is.numeric(location) || length(location) != 1L || !is.finite(location)) {
      .iwmde_ordered_stop("ordered target", "support_unavailable", "the primitive support is not numeric")
    }
    return(rep(location, 2L))
  }
  if (BayesTools::is.prior.discrete(prior)) {
    if (!identical(prior[["distribution"]], "bernoulli")) {
      .iwmde_ordered_stop("ordered target", "support_unavailable", "the discrete support is unavailable")
    }
    support <- c(0, 1)
    bounds <- .iwmde_prior_support(prior)
    support <- support[support >= bounds[[1L]] & support <= bounds[[2L]] &
      is.finite(BayesTools::lpdf(prior, support))]
    return(range(support))
  }
  .iwmde_prior_support(prior)
}


.iwmde_ordered_support <- function(context, parameter_spec) {

  projection <- parameter_spec[["ordered_projection"]]
  support <- c(0, 0)
  multiply <- function(left, right) {
    products <- vapply(left, function(x) vapply(right, function(y) {
      if (x == 0 || y == 0) 0 else x * y
    }, numeric(1)), numeric(2))
    range(products)
  }
  for (term in projection[["contractions"]]) {
    coefficients <- term[["tensor"]][["coefficients"]]
    if (all(coefficients == 0)) next
    prior <- parameter_spec[["ordered_specs"]][[term[["parameter"]]]][["total_prior"]]
    support <- support + multiply(.iwmde_ordered_prior_support(prior), range(coefficients))
  }
  for (name in names(projection[["ordinary_weights"]])) {
    owner <- .iwmde_parameter_prior_name(context, name)
    if (is.null(owner)) .iwmde_ordered_refit()
    prior <- context[["flat_prior_list"]][[owner]]
    support <- support + sort(projection[["ordinary_weights"]][[name]] * .iwmde_ordered_prior_support(prior))
  }
  support
}


.iwmde_ordered_offset_descriptor <- function(context, parameter_spec) {

  projection <- parameter_spec[["ordered_projection"]]
  terms <- projection[["contractions"]]
  terms <- terms[!vapply(terms, function(term) all(term[["tensor"]][["coefficients"]] == 0), logical(1))]
  if (!length(terms) || any(vapply(terms, function(term) length(term[["tensor"]][["records"]]) != 1L, logical(1)))) return(NULL)
  offset <- 0
  for (name in names(projection[["ordinary_weights"]])) {
    owner <- .iwmde_parameter_prior_name(context, name)
    prior <- context[["flat_prior_list"]][[owner]]
    location <- prior[["parameters"]][["location"]]
    if (!BayesTools::is.prior.point(prior) || !is.numeric(location) ||
        length(location) != 1L || !is.finite(location)) return(NULL)
    offset <- offset + projection[["ordinary_weights"]][[name]] * location
  }
  # Several declared totals can multiply the same grouped share. Their affine
  # image is W = a*T + b(U); population support, rather than sampled signs,
  # determines whether the reciprocal supplier is admissible.
  continuous <- which(vapply(terms, function(term) {
    prior <- parameter_spec[["ordered_specs"]][[term[["parameter"]]]][["total_prior"]]
    !BayesTools::is.prior.point(prior) && !BayesTools::is.prior.discrete(prior)
  }, logical(1)))
  pivot <- if (length(continuous)) continuous[[1L]] else 1L
  term <- terms[[pivot]]
  coefficients <- term[["tensor"]][["coefficients"]]
  record <- term[["tensor"]][["records"]][[1L]]
  if (!identical(record[["spec"]][["type"]], "dirichlet")) return(NULL)
  J <- which(coefficients != 0)
  interior <- min(coefficients) < 0 && max(coefficients) > 0
  endpoint <- any(coefficients == 0) && !interior
  if (!interior && !endpoint) return(NULL)
  pure <- length(unique(coefficients)) == 2L && endpoint
  k <- if (pure) coefficients[J][[1L]] else NA_real_
  affine <- list()
  shift_support <- c(0, 0)
  if (length(terms) > 1L) {
    if (!pure) return(NULL)
    for (i in seq_along(terms)) {
      other <- terms[[i]]
      other_record <- other[["tensor"]][["records"]][[1L]]
      values <- other[["tensor"]][["coefficients"]]
      if (!identical(other_record[["key"]], record[["key"]]) ||
          length(values) != length(coefficients) || any(values[-J] != 0) ||
          length(unique(values[J])) != 1L) return(NULL)
      ratio <- values[J][[1L]] / k
      affine[[length(affine) + 1L]] <- list(term = other[["parameter"]], slice = other[["slice"]], a = ratio)
      if (i != pivot) {
        prior <- parameter_spec[["ordered_specs"]][[other[["parameter"]]]][["total_prior"]]
        shift_support <- shift_support + sort(ratio * .iwmde_ordered_prior_support(prior))
      }
    }
  } else affine <- list(list(term = term[["parameter"]], slice = term[["slice"]], a = 1))
  list(offset = offset, term = term[["parameter"]], slice = term[["slice"]],
    record = record, J = J, alpha = sum(record[["spec"]][["alpha"]][J]),
    pure = pure, interior = interior, k = k, a = 1, affine_terms = affine,
    shift_support = shift_support, shifted = any(shift_support != 0))
}


.iwmde_ordered_offset_status <- function(context, parameter_spec) {

  descriptor <- .iwmde_ordered_offset_descriptor(context, parameter_spec)
  if (is.null(descriptor)) return(NULL)
  prior <- parameter_spec[["ordered_specs"]][[descriptor[["term"]]]][["total_prior"]]
  components <- if (BayesTools::is.prior.mixture(prior)) unclass(prior) else list(prior)
  statuses <- lapply(components, function(component) {
    if (.iwmde_retained_location_dynamic_prior(component)) return(list(status = "capability",
      cause = "supplier_provenance", detail = "the active total has a dependent or expression-valued prior"))
    if (BayesTools::is.prior.point(component)) {
      location <- component[["parameters"]][["location"]]
      if (!is.numeric(location) || length(location) != 1L || !is.finite(location)) {
        return(list(status = "capability", cause = "point_supplier", detail = "the point supplier is not a declared numeric constant"))
      }
      if (descriptor[["shifted"]]) {
        shifts <- descriptor[["shift_support"]]
        if (shifts[[1L]] != shifts[[2L]]) return(list(status = "capability",
          cause = "shifted_point_supplier", detail = "the focal point total does not certify an atom or a supplier for its random affine multiplier"))
        location <- location + shifts[[1L]]
      }
      if (location == 0) return(list(status = "atom"))
      if (isTRUE(descriptor[["interior"]])) {
        if (descriptor[["alpha"]] > 1) return(list(status = "raw", supplier = "share", exponent = Inf))
        return(list(status = "capability", cause = "interior_zero_supplier",
          detail = "the interior share limit has no supplied finite or nonfinite certificate"))
      }
      if (descriptor[["alpha"]] < 1) return(list(status = "infinite", cause = "nonfinite_offset",
        detail = "the declared fixed-total share limit is infinite"))
      if (descriptor[["alpha"]] > 1) return(list(status = "zero", cause = "structural_zero_offset"))
      return(list(status = "raw", supplier = "share", exponent = Inf))
    }
    if (BayesTools::is.prior.discrete(component) && identical(component[["distribution"]], "bernoulli")) {
      if (descriptor[["shifted"]]) return(list(status = "capability", cause = "shifted_discrete_supplier",
        detail = "the discrete focal total does not certify the affine projected component measure"))
      if (!is.finite(BayesTools::lpdf(component, 1))) return(list(status = "atom"))
      if (isTRUE(descriptor[["interior"]])) {
        if (descriptor[["alpha"]] > 1) return(list(status = "raw", supplier = "share", exponent = Inf))
        return(list(status = "capability", cause = "interior_zero_supplier",
          detail = "the interior share limit has no supplied finite or nonfinite certificate"))
      }
      if (descriptor[["alpha"]] < 1) return(list(status = "infinite", cause = "nonfinite_offset",
        detail = "the declared nonzero Bernoulli share limit is infinite"))
      if (descriptor[["alpha"]] > 1) return(list(status = "zero", cause = "structural_zero_offset"))
      return(list(status = "raw", supplier = "share", exponent = Inf))
    }
    distribution <- component[["distribution"]]
    support <- .iwmde_prior_support(component)
    shape <- component[["parameters"]][["shape"]]
    if (isTRUE(descriptor[["pure"]]) && distribution == "gamma" && support[[1L]] == 0 && is.finite(support[[2L]]) &&
        descriptor[["alpha"]] == 1 && is.numeric(shape) && length(shape) == 1L &&
        shape > 1 && shape <= 4 && (!descriptor[["shifted"]] || descriptor[["shift_support"]][[1L]] == 0)) {
      return(list(status = "capability", cause = "truncated_gamma_inverse_supplier",
        detail = "the finite-upper-truncated Gamma share supplier has no supplied inverse-prior block route"))
    }
    if (descriptor[["shifted"]]) {
      shifts <- descriptor[["shift_support"]]
      if (distribution == "normal" && all(is.infinite(support)) && descriptor[["alpha"]] > 1) {
        return(list(status = if (descriptor[["alpha"]] <= 4) "block" else "raw", supplier = "total", exponent = descriptor[["alpha"]]))
      }
      if (distribution == "normal" && all(is.infinite(support)) && descriptor[["alpha"]] <= 1) {
        return(list(status = "infinite", cause = "nonfinite_offset", detail = "the declared affine total has positive local density at zero and the inverse share moment is nonintegrable"))
      }
      shape <- component[["parameters"]][["shape"]]
      if (distribution == "gamma" && support[[1L]] == 0 && shifts[[1L]] >= 0 &&
          descriptor[["alpha"]] == 1 && is.numeric(shape) &&
          (shape > 1 || shape >= 1 && shifts[[1L]] > 0)) {
        return(list(status = if (shifts[[1L]] > 0 || shape > 4) "raw" else "block", supplier = "share",
          exponent = if (shifts[[1L]] > 0) Inf else shape))
      }
      return(list(status = "capability", cause = "shifted_supplier_population",
        detail = "the declared affine shift has no supplied uniform supplier and reciprocal-moment certificate"))
    }
    density_zero <- if (0 >= support[[1L]] && 0 <= support[[2L]]) BayesTools::lpdf(component, 0) else -Inf
    if (descriptor[["interior"]]) {
      if (is.finite(density_zero)) return(list(status = "infinite", cause = "nonfinite_offset",
        detail = "the declared interior contraction has a positive total density at zero"))
      return(list(status = "capability", cause = "interior_zero_supplier",
        detail = "the interior-zero supplier is unavailable from the declared prior route"))
    }
    alpha <- descriptor[["alpha"]]
    if (distribution == "normal" && is.finite(density_zero)) {
      if (alpha <= 1) return(list(status = "infinite", cause = "nonfinite_offset",
        detail = "the declared share has a nonintegrable inverse first moment"))
      if (!descriptor[["pure"]] && alpha <= 4) return(list(status = "capability", cause = "nonpure_block",
        detail = "the nonpure endpoint requires an unavailable companion block"))
      return(list(status = if (alpha <= 4) "block" else "raw", supplier = "total", exponent = alpha))
    }
    if (distribution == "gamma" && support[[1L]] == 0) {
      shape <- component[["parameters"]][["shape"]]
      if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape)) {
        return(list(status = "capability", cause = "point_precision", detail = "the Gamma boundary shape is not numeric"))
      }
      if (shape < 1 || alpha < 1 || (shape == 1 && alpha == 1)) {
        return(list(status = "infinite", cause = "nonfinite_offset", detail = "the declared product limit is infinite"))
      }
      if (shape > 1 && alpha > 1) return(list(status = "zero", cause = "structural_zero_offset"))
      if (!descriptor[["pure"]]) return(list(status = "capability", cause = "nonpure_supplier",
        detail = "the nonpure endpoint supplier is unavailable from the declared prior route"))
      exponent <- if (alpha == 1) shape else alpha
      return(list(status = if (exponent <= 4) "block" else "raw",
        supplier = if (alpha == 1) "share" else "total", exponent = exponent))
    }
    gap <- support[[1L]] > 0 || support[[2L]] < 0
    if (descriptor[["pure"]] && gap && alpha == 1 && distribution %in% c("normal", "gamma")) {
      return(list(status = "raw", supplier = "share", exponent = Inf))
    }
    list(status = "capability", cause = "point_precision",
      detail = "the supplier and its population precision certificate are unavailable")
  })
  if (any(vapply(statuses, function(status) identical(status[["status"]], "infinite"), logical(1)))) {
    status <- statuses[[which(vapply(statuses, function(status) identical(status[["status"]], "infinite"), logical(1)))[[1L]]]]
  } else if (length(statuses) == 1L) status <- statuses[[1L]] else if (
      all(vapply(statuses, function(status) status[["status"]] %in% c("block", "raw", "zero", "atom"), logical(1))) &&
      any(vapply(statuses, function(status) identical(status[["status"]], "block"), logical(1)))) {
    status <- list(status = "block", component_statuses = statuses)
  } else if (
      all(vapply(statuses, function(status) status[["status"]] %in% c("raw", "zero", "atom"), logical(1)))) {
    status <- list(status = if (any(vapply(statuses, function(status) identical(status[["status"]], "raw"), logical(1)))) "raw" else "zero")
  } else {
    status <- list(status = "capability", cause = "mixed_supplier",
      detail = "the component suppliers require separate normalized point plans")
  }
  if (length(statuses) > 1L && is.null(status[["component_statuses"]])) status[["component_statuses"]] <- statuses
  c(status, list(descriptor = descriptor, value = descriptor[["offset"]],
    value_key = .iwmde_key_number(descriptor[["offset"]])))
}


.iwmde_ordered_source_normalizer <- function(context, state) {

  gaussian <- state[["source_gaussian"]]
  prior <- state[["ordered_chart"]][["prior"]]
  if (is.null(gaussian) || !prior[["distribution"]] %in% c("normal", "gamma", "beta")) return(NULL)
  unit <- if (.is_data_multilevel(context[["data"]])) "cluster" else "estimate"
  if (!identical(unit, "estimate") || .is_data_known_v(context[["data"]])) return(NULL)
  row <- matrix(state[["row"]], nrow = 1L, dimnames = list(NULL, names(state[["row"]])))
  replayed <- .iwmde_ordered_replay(context, row, state[["ordered_specs"]])
  setup <- .iwmde_log_lik_posterior_setup_active_branch(context, replayed, state[["active_setup"]], unit)
  variance <- as.numeric(setup[["tau_within"]])^2 + setup[["sei"]]^2
  weights <- if (is.null(setup[["weights"]])) rep(1, length(variance)) else setup[["weights"]]
  gaussian_L <- sum(weights * stats::dnorm(setup[["yi"]], as.numeric(setup[["mu"]]), sqrt(variance), log = TRUE))
  prior_constant <- state[["prior_evaluator"]](replayed) + state[["ordered_kernel"]](replayed) -
    BayesTools::lpdf(prior, state[["ordered_chart"]][["current"]])
  current <- state[["ordered_chart"]][["current"]]
  A <- gaussian[["quadratic"]]
  H <- gaussian[["linear"]] + A * current
  C <- gaussian_L + prior_constant - gaussian[["linear"]] * current - .5 * A * current^2
  counts <- new.env(parent = emptyenv())
  for (name in c("candidate_rows", "native_nodes", "inner_calls", "outer_nodes", "maximum_chunk")) counts[[name]] <- 0L
  counts[["max_log_scale"]] <- max(1, abs(C), abs(prior_constant), abs(gaussian_L))
  block <- list(source_state = state, prior = prior, normal = .iwmde_qcmde_normal_prior_parameters(prior),
    normal_probability = .iwmde_ordered_probability_prior(prior), gaussian = list(A = A, H = H, C = C),
    counts = counts, selected = gaussian[["selected"]], log_correction_lower = 0, log_correction_upper = 0)
  if (.iwmde_ordered_source_residual_available(context, state)) {
    block[["gaussian_anchor"]] <- list(A = A, t_ref = current, H_ref = gaussian[["linear"]],
      C_ref = gaussian_L + prior_constant, piece_scale = abs(prior_constant) + abs(gaussian_L),
      evaluate = function(value) .iwmde_ordered_source_gaussian_log(context, state, value, prior_constant))
    counts[["max_log_scale"]] <- max(1, abs(prior_constant), abs(gaussian_L))
  }
  if (block[["selected"]]) {
    selection <- state[["active_setup"]][["priors"]][["outcome"]][["bias"]]
    if (!BayesTools::is.prior.weightfunction(selection) || !identical(selection[["weights"]][["type"]], "fixed") ||
        any(selection[["weights"]][["omega"]] <= 0)) return(NULL)
    selection_context <- .iwmde_selection_context_active_branch(context, state[["active_setup"]], replayed)
    if (is.null(selection_context) || any(selection_context[["kernel_mode"]] != SELKERNEL_STEP)) return(NULL)
    norm <- .selection_step_log_norm_matrix(setup[["mu"]], matrix(sqrt(variance), 1L), setup[["selection_sei"]], selection_context)
    selected_L <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(context, replayed, state[["active_setup"]], unit)
    observed <- .iwmde_ordered_selection_observed_constant(context, state, setup, selection_context)
    if (is.null(observed)) observed <- selected_L - gaussian_L + sum(weights * norm)
    omega <- selection[["weights"]][["omega"]]
    block[["log_correction_lower"]] <- observed - sum(weights) * log(max(omega))
    block[["log_correction_upper"]] <- observed - sum(weights) * log(min(omega))
  }
  untruncated <- !is.null(block[["normal"]]) && identical(.iwmde_prior_support(prior), c(-Inf, Inf))
  if (block[["selected"]] && untruncated) {
    results <- list()
    attempts <- list()
    for (order in c(15L, 31L, 63L, 127L, 255L, 511L, 1023L)) {
      results[[length(results) + 1L]] <- .iwmde_ordered_block_inner(context, block, 0, order)
      if (length(results) < 3L) next
      use <- tail(results, 3L)
      changes <- abs(expm1(diff(vapply(use, `[[`, numeric(1), "log_value"))))
      log_error <- max(.iwmde_ordered_log_difference(use[[1L]][["log_value"]], use[[2L]][["log_value"]]),
        .iwmde_ordered_log_difference(use[[2L]][["log_value"]], use[[3L]][["log_value"]]))
      final <- use[[3L]]
      fraction <- exp(.iwmde_ordered_log_add(log_error, final[["log_exterior"]]) - final[["log_value"]])
      attempts[[length(attempts) + 1L]] <- list(order = order, changes = changes, relative_error = fraction)
      if (all(changes <= .025) && fraction <= .025 / 1.025) {
        final[["log_error"]] <- log_error
        final[["relative_error"]] <- fraction
        final[["attempts"]] <- attempts
        break
      }
      if (order == 1023L) .iwmde_ordered_numeric_stop("source_ladder_exhausted",
        "the selected source normalizer exhausted order 1023", attempts)
    }
  } else final <- .iwmde_ordered_block_inner(context, block, 0)
  list(log_normalizer = final[["log_value"]], discretization_relative_error = exp(final[["log_error"]] - final[["log_value"]]),
    exterior_relative_error = exp(final[["log_exterior"]] - final[["log_value"]]),
    max_log_scale = counts[["max_log_scale"]], method = final[["method"]], proof = final)
}


.iwmde_ordered_affine_shift <- function(parameter_spec, row, source_term = NULL, source_slice = NULL) {

  descriptor <- parameter_spec[["ordered_point_plan"]][["descriptor"]]
  if (is.null(descriptor)) return(0)
  if (is.null(source_term)) source_term <- descriptor[["term"]]
  if (is.null(source_slice)) source_slice <- descriptor[["slice"]]
  shift <- 0
  for (term in descriptor[["affine_terms"]]) {
    if (identical(term[["term"]], source_term) && term[["slice"]] == source_slice) next
    spec <- parameter_spec[["ordered_specs"]][[term[["term"]]]]
    component <- .iwmde_ordered_total_component(spec, row)
    total <- if (BayesTools::is.prior.point(component[["prior"]])) component[["prior"]][["parameters"]][["location"]] else
      row[[component[["coordinates"]][[term[["slice"]]]]]]
    shift <- shift + term[["a"]] * total
  }
  unname(shift)
}


.iwmde_ordered_point_certificate <- function(context, parameter_spec, row_states) {

  unavailable <- function(cause, detail) list(available = FALSE, cause = cause, detail = detail)
  if (all(vapply(row_states, function(state) identical(state[["ordered_chart"]][["kind"]], "ordinary"), logical(1)))) {
    return(list(available = TRUE, kind = "existing_ordinary_pivot"))
  }
  block_point <- identical(parameter_spec[["ordered_point_plan"]][["status"]], "block")
  if (block_point && (.is_data_random(context[["data"]]) || .is_data_multilevel(context[["data"]]))) {
    return(unavailable("covariance_certificate", "the offset block requires the supplied invariant diagonal Gaussian covariance"))
  }
  if (block_point) {
    descriptor <- parameter_spec[["ordered_point_plan"]][["descriptor"]]
    spec <- parameter_spec[["ordered_specs"]][[descriptor[["term"]]]]
    supported <- vapply(row_states, function(state) {
      component <- .iwmde_ordered_total_component(spec, state[["row"]])
      statuses <- parameter_spec[["ordered_point_plan"]][["component_statuses"]]
      if (length(statuses) && !identical(statuses[[component[["component"]]]][["status"]], "block")) return(TRUE)
      prior <- component[["prior"]]
      prior[["distribution"]] %in% c("normal", "gamma") && !.iwmde_retained_location_dynamic_prior(prior)
    }, logical(1))
    if (!all(supported)) return(unavailable("offset_block_capability_unavailable",
      "the source prior is outside the supplied Gaussian block implementation"))
  }
  if (!identical(.data_outcome_type(context[["data"]]), "norm") ||
      .iwmde_context_uses_local_likelihood(context) || .is_data_known_v(context[["data"]]) ||
      any(vapply(row_states, function(state) is.null(state[["source_gaussian"]]), logical(1)))) {
    return(unavailable("covariance_certificate", "the invariant Gaussian source and fixed positive measurement covariance are not certified"))
  }
  object <- context[["object"]]
  object[["data"]] <- context[["data"]]
  sei <- .outcome_data_sei(object)
  if (!length(sei) || any(!is.finite(sei)) || any(sei <= 0)) {
    return(unavailable("covariance_floor", "the measurement covariance has no declared positive diagonal floor"))
  }
  numeric_moments <- function(prior) {
    if (BayesTools::is.prior.ordered(prior)) return(numeric_moments(prior[["total"]]))
    if (BayesTools::is.prior.mixture(prior)) return(all(vapply(prior, numeric_moments, logical(1))))
    if (.iwmde_retained_location_dynamic_prior(prior)) return(FALSE)
    if (BayesTools::is.prior.point(prior) || BayesTools::is.prior.discrete(prior)) return(TRUE)
    prior[["distribution"]] %in% c("normal", "gamma", "lnorm", "lognormal", "beta", "uniform")
  }
  means <- BayesTools::JAGS_formula_design(context[["object"]][["fit"]], parameter = "mu")[["prior_list"]]
  coordinates <- BayesTools::parameter_coordinates(context[["object"]][["fit"]])[["coordinate_name"]]
  dependencies <- BayesTools::JAGS_formula_coordinate_dependencies(context[["object"]][["fit"]], coordinates)
  if (any(dependencies[["formula_parameter"]] == "mu" & dependencies[["dependency_type"]] != "coefficient")) {
    return(unavailable("mean_dependency_certificate", "a mean dependency has no supplied independent numeric primitive certificate"))
  }
  if (!all(vapply(means, numeric_moments, logical(1)))) {
    return(unavailable("point_precision", "a retained mean source lacks an independent numeric all-polynomial-moment certificate"))
  }
  for (state in row_states) {
    prior <- state[["ordered_chart"]][["prior"]]
    if (identical(prior[["distribution"]], "normal")) {
      support <- .iwmde_prior_support(prior)
      full <- identical(support, c(-Inf, Inf))
      half <- identical(support, c(0, Inf)) && identical(prior[["parameters"]][["mean"]], 0)
      if (!full && !half) return(unavailable("point_precision", "the truncated Normal source has no supplied population ordinate-moment certificate"))
    }
    if (identical(prior[["distribution"]], "gamma") && prior[["parameters"]][["shape"]] < 1 ||
        identical(prior[["distribution"]], "beta") &&
        any(unlist(prior[["parameters"]][c("alpha", "beta")], use.names = FALSE) < 1) ||
        !prior[["distribution"]] %in% c("normal", "gamma", "beta")) {
      return(unavailable("point_precision", "the moving source family has no supplied root-weighted or local precision certificate"))
    }
    if (isTRUE(state[["active_setup"]][["is_weightfunction"]])) {
      selection <- state[["active_setup"]][["priors"]][["outcome"]][["bias"]]
      plan <- .data_selection_execution_plan(context[["data"]])
      model <- .data_selection_model(context[["data"]])
      branches <- model[["branches"]][model[["active_branches"]]]
      if (!BayesTools::is.prior.weightfunction(selection) ||
          !identical(selection[["weights"]][["type"]], "fixed") ||
          any(!is.finite(selection[["weights"]][["omega"]])) || any(selection[["weights"]][["omega"]] <= 0) ||
          !length(plan[["block_methods"]]) || any(plan[["block_methods"]] != "singleton") ||
          !all(vapply(branches, function(branch) identical(branch[["weight_rule"]], "product"), logical(1))) ||
          .selection_retains_sampling(context[["data"]]) || .selection_retains_any_random(context[["data"]])) {
        return(unavailable("selection_precision", "the selection route lacks fixed positive partition weights and invariant singleton product integration"))
      }
    }
  }
  if (is.null(parameter_spec[["ordered_point_plan"]])) {
    contractions <- parameter_spec[["ordered_projection"]][["contractions"]]
    bounded <- vapply(contractions, function(term) {
      coefficients <- term[["tensor"]][["coefficients"]]
      all(coefficients == 0) || all(coefficients > 0) || all(coefficients < 0)
    }, logical(1))
    if (!all(bounded)) return(unavailable("point_precision", "the target has no supplied supplier or uniformly separated allocation multiplier certificate"))
  }
  list(available = TRUE, kind = "independent_numeric_gaussian_core",
    measurement_variance_floor = min(sei^2), moment_basis = "declared_primitive_priors")
}


.iwmde_ordered_constant_certificate <- function(context, parameter_spec, row_states) {

  if (.is_data_random(context[["data"]]) || .is_data_scale(context[["data"]]) ||
      .is_data_known_v(context[["data"]]) || .is_data_multilevel(context[["data"]]) ||
      .iwmde_context_uses_local_likelihood(context)) return(FALSE)
  charts <- lapply(row_states, `[[`, "ordered_chart")
  if (!length(charts) || !all(vapply(charts, function(chart) chart[["kind"]] %in% c("total", "ordinary"), logical(1)))) return(FALSE)
  source <- unique(vapply(charts, `[[`, character(1), "source"))
  if (length(source) != 1L) return(FALSE)
  focal <- .iwmde_parameter_prior_name(context, source)
  focal_ordered <- unique(vapply(charts, function(chart) if (is.null(chart[["term"]])) "" else chart[["term"]], character(1)))
  for (name in names(context[["flat_prior_list"]])) {
    prior <- context[["flat_prior_list"]][[name]]
    if (name %in% focal_ordered) {
      if (length(parameter_spec[["ordered_specs"]][[name]][["total_names"]]) != 1L ||
          BayesTools::is.prior.mixture(prior[["total"]]) ||
          any(vapply(parameter_spec[["ordered_specs"]][[name]][["allocations"]], function(record) {
            !identical(record[["spec"]][["type"]], "fixed") && record[["dim"]] > 1L
          }, logical(1)))) return(FALSE)
      next
    }
    if (identical(name, focal)) next
    if (BayesTools::is.prior.ordered(prior)) {
      if (!BayesTools::is.prior.point(prior[["total"]]) ||
          !is.numeric(prior[["total"]][["parameters"]][["location"]]) ||
          prior[["total"]][["parameters"]][["location"]] != 0) return(FALSE)
      next
    }
    if (BayesTools::is.prior.weightfunction(prior)) {
      if (!identical(prior[["weights"]][["type"]], "fixed")) return(FALSE)
      next
    }
    if (!BayesTools::is.prior.point(prior) || .iwmde_retained_location_dynamic_prior(prior)) return(FALSE)
  }
  TRUE
}


.iwmde_ordered_target_spec <- function(context, parameter, parameter_spec) {

  if (identical(parameter_spec[["type"]], "ordered")) return(parameter_spec)
  if (identical(parameter_spec[["status"]], "unsupported")) return(parameter_spec)
  if (!is.null(parameter_spec) &&
      !parameter_spec[["type"]] %in% c("primitive", "linear")) return(parameter_spec)
  fit <- context[["object"]][["fit"]]
  priors <- context[["flat_prior_list"]]
  ordered <- names(priors)[vapply(priors, BayesTools::is.prior.ordered, logical(1))]
  if (!length(ordered) || is.null(fit)) return(parameter_spec)
  weights <- if (identical(parameter_spec[["type"]], "linear")) {
    .iwmde_linear_weights(parameter_spec[["weights"]])
  } else stats::setNames(1, parameter)
  if (!length(weights)) return(parameter_spec)
  prior_owners <- vapply(names(weights), function(name) {
    owner <- .iwmde_parameter_prior_name(context, name)
    if (is.null(owner)) "" else owner
  }, character(1))
  known_ordinary <- nzchar(prior_owners) & !prior_owners %in% ordered
  if (all(known_ordinary)) return(parameter_spec)

  # Select owners from persisted coordinate/catalog rows before asking for any
  # ordered provenance. An ordinary request never validates unrelated terms.
  coordinates <- BayesTools::parameter_coordinates(fit)
  selected <- coordinates[coordinates[["coordinate_name"]] %in% names(weights), , drop = FALSE]
  quantities <- BayesTools::parameter_catalog(fit)[["quantities"]]
  selected_quantities <- quantities[quantities[["canonical_name"]] %in% names(weights), , drop = FALSE]
  owner_rows <- rbind(selected[, c("term", "formula_parameter"), drop = FALSE],
    selected_quantities[, c("term", "formula_parameter"), drop = FALSE])
  owners <- unique(unlist(lapply(seq_len(nrow(owner_rows)), function(i) {
    term <- owner_rows[["term"]][[i]]
    formula <- owner_rows[["formula_parameter"]][[i]]
    if (is.na(term) || !nzchar(term) || is.na(formula) || !nzchar(formula)) return(character())
    BayesTools::JAGS_parameter_names(term, formula_parameter = formula)
  }), use.names = FALSE))
  owners <- intersect(owners, ordered)
  unresolved <- names(weights)[!known_ordinary]
  resolved_coefficients <- coordinates[["coordinate_name"]][
    coordinates[["role"]] == "fixed_coefficient" & coordinates[["coordinate_name"]] %in% unresolved]
  unresolved <- setdiff(unresolved, resolved_coefficients)
  unresolved <- intersect(unresolved, coordinates[["coordinate_name"]])
  if (length(unresolved)) {
    nodes <- .iwmde_deterministic_nodes(context)
    candidates <- nodes[["family"]] == "ordered_coefficient" & nodes[["parameter"]] %in% ordered
    matched <- vapply(seq_len(nrow(nodes)), function(i) candidates[[i]] &&
      any(unresolved %in% c(nodes[["coordinates"]][[i]], nodes[["dependencies"]][[i]])), logical(1))
    owners <- unique(c(owners, nodes[["parameter"]][matched]))
  }
  if (!length(owners)) return(parameter_spec)
  selected_specs <- BayesTools::JAGS_ordered_parameter_spec(fit, parameters = owners)
  affected <- unlist(lapply(selected_specs, function(spec) {
    c(spec[["coefficient_names"]], spec[["total_names"]])
  }), use.names = FALSE)
  if (!any(names(weights) %in% affected)) return(parameter_spec)

  # The complete prior/likelihood requires every ordered consumer, including
  # terms sharing the moved allocation key. Keep the original bound wrappers.
  specs <- BayesTools::JAGS_ordered_parameter_spec(fit, parameters = ordered)
  draws <- .iwmde_ordered_source_draws(context, specs, coordinates)
  projection <- BayesTools::JAGS_ordered_parameter_spec(fit, parameters = ordered,
    weights = weights, draws = draws)
  if (any(projection[["state"]] == "unavailable")) stop(projection[["reason"]])
  if (is.null(parameter_spec)) parameter_spec <- list()
  parameter_spec[["type"]] <- "ordered"
  parameter_spec[["parameter"]] <- parameter
  parameter_spec[["weights"]] <- weights
  parameter_spec[["ordered_specs"]] <- specs
  parameter_spec[["ordered_projection"]] <- projection
  parameter_spec[["ordered_source_key"]] <- .iwmde_hash("ordered_sources", list(
    names(specs), weights, .iwmde_source_fingerprint(context)))
  parameter_spec[["target_columns"]] <- names(weights)
  parameter_spec[["conditioning_exclude"]] <- unique(unlist(lapply(specs, function(spec) {
    c(spec[["coefficient_names"]], spec[["total_names"]], spec[["source_coordinates"]],
      unlist(lapply(spec[["allocations"]], `[[`, "coordinates"), use.names = FALSE))
  }), use.names = FALSE))
  parameter_spec[["status"]] <- "ok"
  parameter_spec[["ordered_point_plan"]] <- .iwmde_ordered_offset_status(context, parameter_spec)
  context[["evaluator_cache"]][[parameter_spec[["ordered_source_key"]]]] <- draws
  parameter_spec
}


.iwmde_ordered_source_draws <- function(context, specs, coordinates) {

  candidates <- unique(unlist(lapply(specs, function(spec) {
    c(spec[["source_coordinates"]], spec[["total_auxiliary_coordinates"]])
  }), use.names = FALSE))
  available <- coordinates[["coordinate_name"]][
    coordinates[["monitor_status"]] %in% c("sampled", "structural")]
  requested <- intersect(candidates, available)
  draws <- context[["posterior_samples"]]
  if (length(requested)) {
    sources <- as.matrix(BayesTools::JAGS_materialize_draws(context[["object"]][["fit"]],
      parameters = requested, include_internal = TRUE))
    if (nrow(sources) != nrow(draws)) .iwmde_ordered_refit()
    draws <- cbind(draws[, !colnames(draws) %in% colnames(sources), drop = FALSE], sources)
  }
  required_gamma <- unique(unlist(lapply(specs, function(spec) {
    unlist(lapply(spec[["allocations"]], `[[`, "gamma_coordinates"), use.names = FALSE)
  }), use.names = FALSE))
  if (!all(required_gamma %in% colnames(draws))) .iwmde_ordered_refit()
  # Optional unmonitored mixture sources can use the monitored fitted total
  # snapshot. A declared real source is never overridden by that snapshot.
  for (spec in specs) {
    node <- spec[["total_node"]]
    if (is.null(node)) next
    indicator <- node[["spec"]][["indicator"]]
    if (!indicator %in% colnames(draws)) .iwmde_ordered_refit()
    components <- node[["spec"]][["components"]]
    for (i in seq_along(components)) {
      component <- components[[i]]
      if (!is.null(component[["location"]])) next
      source <- component[["coordinates"]]
      absent <- setdiff(source, colnames(draws))
      if (!length(absent)) next
      if (any(absent %in% available) || !all(spec[["total_names"]] %in% colnames(draws))) {
        .iwmde_ordered_refit()
      }
      active <- if (node[["spec"]][["kind"]] %in% c("spike_and_slab", "ordered_spike_and_slab")) {
        draws[, indicator] == 1
      } else draws[, indicator] == i
      added <- matrix(0, nrow(draws), length(absent), dimnames = list(NULL, absent))
      slices <- match(absent, source)
      added[active, ] <- draws[active, spec[["total_names"]][slices], drop = FALSE]
      draws <- cbind(draws, added)
    }
  }
  draws
}


.iwmde_ordered_components <- function(parameter_spec) {

  projection <- parameter_spec[["ordered_projection"]]
  list(active = projection[["state"]] == "continuous",
    point_location = projection[["atom"]],
    point_masses = .iwmde_point_mass_table(projection[["atom"]], length(projection[["values"]])))
}


.iwmde_ordered_total_component <- function(spec, row) {

  prior <- spec[["total_prior"]]
  node <- spec[["total_node"]]
  if (is.null(node)) return(list(prior = prior, component = 1L,
    coordinates = spec[["total_names"]]))
  indicator <- row[[node[["spec"]][["indicator"]]]]
  if (node[["spec"]][["kind"]] %in% c("spike_and_slab", "ordered_spike_and_slab")) {
    null <- which(vapply(prior, function(part) BayesTools::is.prior.point(part) &&
      is.numeric(part[["parameters"]][["location"]]) && part[["parameters"]][["location"]] == 0,
      logical(1)))
    if (length(null) != 1L || !indicator %in% c(0, 1)) {
      .iwmde_ordered_stop(spec[["parameter"]], "invalid_state", "the declared total state is invalid")
    }
    component <- if (indicator == 0) null else setdiff(seq_along(prior), null)
    coordinates <- node[["spec"]][["components"]][[1L]][["coordinates"]]
  } else {
    if (!indicator %in% seq_along(prior)) {
      .iwmde_ordered_stop(spec[["parameter"]], "invalid_state", "the declared total state is invalid")
    }
    component <- as.integer(indicator)
    coordinates <- node[["spec"]][["components"]][[component]][["coordinates"]]
  }
  list(prior = prior[[component]], component = component, coordinates = coordinates)
}


.iwmde_ordered_row_chart <- function(context, parameter_spec, row_index) {

  draws <- context[["evaluator_cache"]][[parameter_spec[["ordered_source_key"]]]]
  row <- draws[row_index, ]
  projection <- parameter_spec[["ordered_projection"]]
  parameter <- parameter_spec[["parameter"]]
  ordinary <- projection[["ordinary_weights"]]
  for (name in names(ordinary)) {
    state <- .iwmde_focal_prior_state(context, name, row)
    if (identical(state[["status"]], "continuous")) {
      return(list(kind = "ordinary", source = name, prior = state[["prior"]],
        current = row[[name]], slope = ordinary[[name]],
        offset = projection[["values"]][[row_index]] - ordinary[[name]] * row[[name]],
        component = 1L))
    }
  }
  terms <- projection[["contractions"]]
  terms <- terms[!vapply(terms, function(term) all(term[["tensor"]][["coefficients"]] == 0), logical(1))]
  point_plan <- parameter_spec[["ordered_point_plan"]]
  descriptor <- point_plan[["descriptor"]]
  chart_component <- 1L
  if (length(point_plan[["component_statuses"]])) {
    focal_spec <- parameter_spec[["ordered_specs"]][[descriptor[["term"]]]]
    chart_component <- .iwmde_ordered_total_component(focal_spec, row)[["component"]]
    point_plan <- point_plan[["component_statuses"]][[chart_component]]
  }
  if (!is.null(descriptor) && !isTRUE(descriptor[["interior"]]) && identical(point_plan[["supplier"]], "share")) {
    return(.iwmde_ordered_allocation_chart(context, parameter_spec, projection, row, row_index,
      descriptor[["record"]], J = descriptor[["J"]], component = chart_component))
  }
  for (term in terms) {
    spec <- parameter_spec[["ordered_specs"]][[term[["parameter"]]]]
    total <- .iwmde_ordered_total_component(spec, row)
    if (BayesTools::is.prior.point(total[["prior"]]) ||
        BayesTools::is.prior.discrete(total[["prior"]])) next
    tensor <- term[["tensor"]]
    # Supplier precedes total priority for the pure Gamma/Beta endpoint case.
    records <- tensor[["records"]]
    coefficients <- tensor[["coefficients"]]
    if (length(terms) == 1L && length(records) == 1L &&
        length(unique(coefficients)) == 2L && any(coefficients == 0)) {
      record <- records[[1L]]
      J <- which(coefficients != 0)
      alpha <- sum(record[["spec"]][["alpha"]][J])
      if (identical(total[["prior"]][["distribution"]], "gamma") && alpha == 1 &&
          total[["prior"]][["parameters"]][["shape"]] > 1) {
        return(.iwmde_ordered_allocation_chart(context, parameter_spec, projection, row, row_index,
          record, J = J, component = total[["component"]]))
      }
    }
    source <- total[["coordinates"]][[term[["slice"]]]]
    if (!source %in% names(row)) .iwmde_ordered_refit()
    slope <- projection[["total_contraction_weights"]][[term[["parameter"]]]][row_index, term[["slice"]]]
    if (!is.finite(slope) || slope == 0) {
      .iwmde_stop_construction_failure(context[["ordered_requested_method"]], parameter, row_index,
        "source-chart construction", "the continuous source slope is zero or unrepresentable")
    }
    return(list(kind = "total", source = source, total = spec[["total_names"]][[term[["slice"]]]],
      term = term[["parameter"]], slice = term[["slice"]], prior = total[["prior"]],
      component = total[["component"]], current = row[[source]], slope = slope,
      offset = if (!is.null(descriptor) && isTRUE(descriptor[["pure"]])) {
        a <- descriptor[["affine_terms"]][[which(vapply(descriptor[["affine_terms"]], function(part)
          identical(part[["term"]], term[["parameter"]]) && part[["slice"]] == term[["slice"]], logical(1)))[[1L]]]][["a"]]
        unname(descriptor[["offset"]] + slope / a * .iwmde_ordered_affine_shift(parameter_spec, row, term[["parameter"]], term[["slice"]]))
      } else projection[["values"]][[row_index]] - slope * row[[source]]))
  }
  for (term in terms) {
    for (record in term[["tensor"]][["records"]]) {
      coefficients <- term[["tensor"]][["coefficients"]]
      if (length(term[["tensor"]][["records"]]) == 1L && length(unique(coefficients)) == 2L) {
        J <- if (any(coefficients == 0)) which(coefficients != 0) else
          which(coefficients == coefficients[[1L]])
        return(.iwmde_ordered_allocation_chart(context, parameter_spec, projection, row, row_index, record, J = J))
      }
      contraction <- projection[["allocation_contractions"]][[record[["key"]]]][row_index, ]
      pair <- which(contraction != contraction[[1L]])
      if (length(pair)) return(.iwmde_ordered_allocation_chart(context, parameter_spec, projection, row,
        row_index, record, indices = c(1L, pair[[1L]])))
    }
  }
  .iwmde_ordered_stop(parameter, "chart_unavailable", "no declared affine continuous source chart is available")
}


.iwmde_ordered_allocation_chart <- function(context, parameter_spec, projection, row,
                                            row_index, record, J = NULL,
                                            indices = NULL, component = 1L) {

  eta <- row[record[["gamma_coordinates"]]]
  H <- sum(eta)
  contraction <- projection[["allocation_contractions"]][[record[["key"]]]][row_index, ]
  if (!is.null(J)) {
    K <- setdiff(seq_along(eta), J)
    r_J <- eta[J] / sum(eta[J])
    r_K <- eta[K] / sum(eta[K])
    group_value <- function(values, ratios) {
      if (all(values == values[[1L]])) values[[1L]] else sum(ratios * values)
    }
    slope <- group_value(contraction[J], r_J) - group_value(contraction[K], r_K)
    current <- sum(eta[J]) / H
    shapes <- c(sum(record[["spec"]][["alpha"]][J]), sum(record[["spec"]][["alpha"]][K]))
    chart <- list(kind = "group", J = J, K = K, H = H, r_J = r_J, r_K = r_K)
  } else {
    h <- sum(eta[indices])
    slope <- h / H * (contraction[indices[[1L]]] - contraction[indices[[2L]]])
    current <- eta[indices[[1L]]] / h
    shapes <- record[["spec"]][["alpha"]][indices]
    chart <- list(kind = "pair", indices = indices, h = h)
  }
  if (!is.finite(slope) || slope == 0 || !is.finite(current)) {
    .iwmde_stop_construction_failure(context[["ordered_requested_method"]], parameter_spec[["parameter"]],
      row_index, "source-chart construction", "the continuous source slope is zero or unrepresentable")
  }
  descriptor <- parameter_spec[["ordered_point_plan"]][["descriptor"]]
  if (is.null(descriptor)) descriptor <- .iwmde_ordered_offset_descriptor(context, parameter_spec)
  offset <- if (!is.null(descriptor) && !isTRUE(descriptor[["interior"]]) && identical(record[["key"]], descriptor[["record"]][["key"]]) &&
      !is.null(J) && identical(as.integer(J), as.integer(descriptor[["J"]]))) {
    descriptor[["offset"]]
  } else projection[["values"]][[row_index]] - slope * current
  c(chart, list(key = record[["key"]], gamma_coordinates = record[["gamma_coordinates"]],
    prior = BayesTools::prior("beta", as.list(shapes)), component = component,
    current = current, slope = slope,
    offset = offset))
}


.iwmde_ordered_replay <- function(context, samples, specs) {

  nodes <- .iwmde_deterministic_nodes(context)
  selected <- nodes[["node"]][nodes[["family"]] %in% c("ordered_allocation", "ordered_coefficient")]
  changed <- unique(unlist(lapply(specs, function(spec) {
    c(spec[["source_coordinates"]], spec[["coefficient_names"]])
  }), use.names = FALSE))
  repeat {
    dependent <- vapply(nodes[["dependencies"]], function(x) any(x %in% changed), logical(1))
    expanded <- unique(c(selected, nodes[["node"]][dependent]))
    if (identical(expanded, selected)) break
    selected <- expanded
    changed <- unique(c(changed, unlist(nodes[["coordinates"]][nodes[["node"]] %in% selected], use.names = FALSE)))
  }
  evaluated <- .iwmde_deterministic_evaluator(context, selected)(samples)
  cbind(samples[, !colnames(samples) %in% colnames(evaluated), drop = FALSE], evaluated)
}


.iwmde_ordered_replace_source <- function(context, state, values) {

  chart <- state[["ordered_chart"]]
  samples <- matrix(rep(state[["row"]], each = length(values)), nrow = length(values),
    dimnames = list(NULL, names(state[["row"]])))
  support <- .iwmde_prior_support(chart[["prior"]])
  valid <- is.finite(values) & values >= support[[1L]] & values <= support[[2L]]
  if (chart[["kind"]] %in% c("total", "ordinary")) {
    samples[, chart[["source"]]] <- values
  } else if (identical(chart[["kind"]], "group")) {
    samples[, chart[["gamma_coordinates"]][chart[["J"]]]] <-
      outer(chart[["H"]] * values, chart[["r_J"]])
    samples[, chart[["gamma_coordinates"]][chart[["K"]]]] <-
      outer(chart[["H"]] * (1 - values), chart[["r_K"]])
  } else if (identical(chart[["kind"]], "pair")) {
    samples[, chart[["gamma_coordinates"]][chart[["indices"]][[1L]]]] <- chart[["h"]] * values
    samples[, chart[["gamma_coordinates"]][chart[["indices"]][[2L]]]] <- chart[["h"]] * (1 - values)
  } else .iwmde_ordered_stop(state[["parameter"]], "chart_unavailable", "the source chart kind is unknown")
  if (any(valid)) {
    replayed <- .iwmde_ordered_replay(context, samples[valid, , drop = FALSE], state[["ordered_specs"]])
    missing <- setdiff(colnames(replayed), colnames(samples))
    if (length(missing)) samples <- cbind(samples, matrix(NA_real_, nrow(samples), length(missing),
      dimnames = list(NULL, missing)))
    samples[valid, colnames(replayed)] <- replayed
  }
  list(samples = samples, valid = valid)
}


.iwmde_ordered_row_state <- function(context, parameter_spec, row_index) {

  draws <- context[["evaluator_cache"]][[parameter_spec[["ordered_source_key"]]]]
  chart <- .iwmde_ordered_row_chart(context, parameter_spec, row_index)
  row <- draws[row_index, ]
  active_key <- .iwmde_active_keys(context)[[row_index]]
  active_setup <- .iwmde_active_setup(context, row, active_key)
  state_scope <- .iwmde_state_scope(parameter_spec[["parameter"]], parameter_spec, context)
  prior_list <- .iwmde_active_flat_prior_list(context, row,
    parameter_spec[["parameter"]], state_scope = state_scope)
  ordered <- parameter_spec[["ordered_specs"]]
  generic <- prior_list[!vapply(prior_list, BayesTools::is.prior.ordered, logical(1))]
  allocation_chart <- if (chart[["kind"]] %in% c("group", "pair")) {
    stats::setNames(list(chart[intersect(c("kind", "J", "indices"), names(chart))]), chart[["key"]])
  } else NULL
  state <- list(row = row, row_index = row_index, active_key = active_key,
    active_setup = active_setup, parameter = parameter_spec[["parameter"]],
    ordered_source_key = parameter_spec[["ordered_source_key"]],
    ordered_specs = ordered, ordered_chart = chart,
    ordered_kernel = BayesTools::JAGS_ordered_density_kernel(lapply(ordered, `[[`, "prior"), allocation_chart),
    prior_list = generic,
    prior_evaluator = BayesTools::JAGS_marglik_priors_rows_evaluator(generic),
    focal_prior = chart[["prior"]], likelihood_mode = .iwmde_likelihood_mode(parameter_spec[["parameter"]], parameter_spec, context),
    state_scope = state_scope,
    use_focal_prior_delta = FALSE)
  state[["baseline_log_q"]] <- .iwmde_ordered_log_q_source(context, state, chart[["current"]])[[1L]]
  state[["source_gaussian"]] <- .iwmde_ordered_source_gaussian(context, state)
  state
}


.iwmde_ordered_log_q_source <- function(context, state, values) {

  if (length(values) > 1024L) {
    chunks <- split(seq_along(values), ceiling(seq_along(values) / 1024L))
    evaluated <- lapply(chunks, function(take) .iwmde_ordered_log_q_source(context, state, values[take]))
    out <- unlist(evaluated, use.names = FALSE)
    attr(out, "max_accumulated_logterm") <- unlist(lapply(evaluated, attr, "max_accumulated_logterm", exact = TRUE), use.names = FALSE)
    return(out)
  }

  replaced <- .iwmde_ordered_replace_source(context, state, values)
  out <- rep(-Inf, length(values))
  valid <- replaced[["valid"]]
  if (!any(valid)) {
    attr(out, "max_accumulated_logterm") <- numeric(length(values))
    return(out)
  }
  samples <- replaced[["samples"]][valid, , drop = FALSE]
  unit <- if (.is_data_multilevel(context[["data"]])) "cluster" else "estimate"
  likelihood <- if (isTRUE(.data_outcome_type(context[["data"]]) %in% c("bin", "pois")) &&
      identical(state[["likelihood_mode"]], "conditional")) {
    .iwmde_glmm_conditional_log_likelihood_samples(context, samples, state[["active_setup"]])
  } else .iwmde_log_lik_from_posterior_samples_sum_active_branch(context,
    samples, state[["active_setup"]], unit)
  ordinary <- state[["prior_evaluator"]](samples)
  ordered <- state[["ordered_kernel"]](samples)
  out[valid] <- likelihood + ordinary + ordered
  scales <- numeric(length(values))
  scales[valid] <- pmax(abs(likelihood), abs(ordinary), abs(ordered), abs(out[valid]))
  attr(out, "max_accumulated_logterm") <- scales
  out
}


.iwmde_ordered_source_grid <- function(context, parameter, values, row_states) {

  if (!length(values) || !length(row_states)) {
    return(matrix(numeric(), nrow = length(values), ncol = length(row_states)))
  }
  out <- vapply(row_states, function(state) .iwmde_ordered_log_q_source(context, state, values),
    numeric(length(values)))
  out <- matrix(out, nrow = length(values), ncol = length(row_states))
  gaussian <- lapply(row_states, `[[`, "source_gaussian")
  if (all(vapply(gaussian, function(kernel) !is.null(kernel) && !kernel[["selected"]], logical(1)))) {
    attr(out, "gaussian_kernel") <- list(
      current = vapply(row_states, function(state) state[["ordered_chart"]][["current"]], numeric(1)),
      linear = vapply(gaussian, `[[`, numeric(1), "linear"),
      quadratic = vapply(gaussian, `[[`, numeric(1), "quadratic"), prior_route = rep("focal", length(gaussian)))
  }
  out
}


.iwmde_ordered_source_gaussian <- function(context, state) {

  if (is.null(context[["data"]]) || is.null(state[["active_setup"]][["priors"]]) ||
      !identical(.data_outcome_type(context[["data"]]), "norm") ||
      .iwmde_retained_location_dynamic_prior(state[["prior_list"]])) return(NULL)
  chart <- state[["ordered_chart"]]
  fit <- context[["object"]][["fit"]]
  directions <- numeric()
  draws <- matrix(state[["row"]], nrow = 1L, dimnames = list(NULL, names(state[["row"]])))
  projection_key <- .iwmde_hash("ordered_coefficient_projections", list(
    names(state[["ordered_specs"]]), .iwmde_source_fingerprint(context)))
  projections <- context[["evaluator_cache"]][[projection_key]]
  if (is.null(projections)) {
    source_draws <- context[["evaluator_cache"]][[state[["ordered_source_key"]]]]
    coefficients <- unlist(lapply(state[["ordered_specs"]], `[[`, "coefficient_names"), use.names = FALSE)
    projections <- lapply(coefficients, function(coefficient) {
      BayesTools::JAGS_ordered_parameter_spec(fit, parameters = names(state[["ordered_specs"]]),
        weights = stats::setNames(1, coefficient), draws = source_draws)
    })
    names(projections) <- coefficients
    context[["evaluator_cache"]][[projection_key]] <- projections
  }
  index <- state[["row_index"]]
  for (spec in state[["ordered_specs"]]) {
    for (coefficient in spec[["coefficient_names"]]) {
      projection <- projections[[coefficient]]
      derivative <- 0
      if (identical(chart[["kind"]], "total")) {
        component <- .iwmde_ordered_total_component(spec, state[["row"]])
        slices <- which(component[["coordinates"]] == chart[["source"]])
        if (length(slices)) derivative <- sum(projection[["total_contraction_weights"]][[spec[["parameter"]]]][index, slices])
      } else if (chart[["kind"]] %in% c("group", "pair")) {
        contraction <- projection[["allocation_contractions"]][[chart[["key"]]]]
        if (!is.null(contraction)) {
          if (identical(chart[["kind"]], "group")) {
            derivative <- sum(chart[["r_J"]] * contraction[index, chart[["J"]]]) -
              sum(chart[["r_K"]] * contraction[index, chart[["K"]]])
          } else {
            H <- sum(state[["row"]][chart[["gamma_coordinates"]]])
            derivative <- chart[["h"]] / H *
              (contraction[index, chart[["indices"]][[1L]]] - contraction[index, chart[["indices"]][[2L]]])
          }
        }
      }
      if (derivative != 0) directions[[coefficient]] <- derivative
    }
  }
  if (identical(chart[["kind"]], "ordinary")) directions[[chart[["source"]]]] <- 1
  if (!length(directions)) return(NULL)
  moving <- if (chart[["kind"]] %in% c("total", "ordinary")) chart[["source"]] else chart[["gamma_coordinates"]]
  nodes <- .iwmde_deterministic_nodes(context)
  repeat {
    selected <- vapply(nodes[["dependencies"]], function(x) any(x %in% moving), logical(1))
    expanded <- unique(c(moving, unlist(nodes[["coordinates"]][selected], use.names = FALSE)))
    if (identical(expanded, moving)) break
    moving <- expanded
  }
  if (any(!nodes[["family"]][selected] %in% c("ordered_allocation", "ordered_coefficient",
    "prior_mixture", "linear_predictor"))) return(NULL)
  basis <- BayesTools::JAGS_formula_predictor_basis(fit, directions, posterior_samples = draws)
  if (!identical(basis[["status"]], "affine") || !identical(basis[["coordinate"]], "identity") ||
      !identical(basis[["parameter"]], "mu")) return(NULL)
  replayed <- .iwmde_ordered_replay(context, draws, state[["ordered_specs"]])
  unit <- if (.is_data_multilevel(context[["data"]])) "cluster" else "estimate"
  setup <- .iwmde_log_lik_posterior_setup_active_branch(context, replayed, state[["active_setup"]], unit)
  change <- .iwmde_normal_location_likelihood_change(context, setup, list(mu_basis = basis[["basis"]]))
  if (is.null(change) || length(change[["linear"]]) != 1L || length(change[["quadratic"]]) != 1L ||
      any(!is.finite(c(change[["linear"]], change[["quadratic"]]))) || change[["quadratic"]] < 0) return(NULL)
  c(change, list(selected = isTRUE(state[["active_setup"]][["is_weightfunction"]]),
    directions = directions, basis = basis[["basis"]], mean = as.numeric(setup[["mu"]])))
}


.iwmde_ordered_source_residual_available <- function(context, state) {

  data <- context[["data"]]
  !is.null(state[["source_gaussian"]][["basis"]]) && !is.null(data) &&
    identical(.data_outcome_type(data), "norm") && !.is_data_random(data) &&
    !.is_data_multilevel(data) && !.is_data_known_v(data) && !.iwmde_context_uses_local_likelihood(context)
}


.iwmde_ordered_source_gaussian_log <- function(context, state, values, prior_constant) {

  replaced <- .iwmde_ordered_replace_source(context, state, values)
  samples <- replaced[["samples"]]
  # Gaussian completion can have its unconstrained mean outside a truncated
  # Normal prior. Replay its polynomial predictor there; the two prior masses
  # still integrate only over the declared support.
  if (!all(replaced[["valid"]])) samples <- .iwmde_ordered_replay(context, samples, state[["ordered_specs"]])
  setup <- .iwmde_log_lik_posterior_setup_active_branch(context, samples, state[["active_setup"]], "estimate")
  mean <- setup[["mu"]]
  variance <- setup[["tau_within"]]^2 + matrix(setup[["sei"]]^2, nrow(mean), ncol(mean), byrow = TRUE)
  weights <- if (is.null(setup[["weights"]])) rep(1, ncol(mean)) else setup[["weights"]]
  pieces <- sweep(stats::dnorm(matrix(setup[["yi"]], nrow(mean), ncol(mean), byrow = TRUE),
    mean, sqrt(variance), log = TRUE), 2L, weights, "*")
  kernel <- state[["source_gaussian"]]
  mean_scale <- matrix(abs(kernel[["mean"]]), nrow(mean), ncol(mean), byrow = TRUE) +
    outer(abs(values) + abs(state[["ordered_chart"]][["current"]]), abs(as.numeric(kernel[["basis"]])))
  scales <- vapply(seq_len(nrow(mean)), function(i) {
    residual_block <- list(setup = setup, weights = weights, variance = variance[i, ])
    .iwmde_ordered_gaussian_residual_scale(residual_block, mean[i, ], mean_scale[i, ])
  }, numeric(1))
  list(value = rowSums(pieces) + prior_constant,
    max_log_scale = pmax(1, rowSums(abs(pieces)), abs(prior_constant), scales))
}


.iwmde_ordered_normal_source_mass <- function(state, context = NULL) {

  kernel <- state[["source_gaussian"]]
  normal <- .iwmde_qcmde_normal_prior_parameters(state[["ordered_chart"]][["prior"]])
  if (is.null(kernel) || kernel[["selected"]] || is.null(normal)) return(NULL)
  current <- state[["ordered_chart"]][["current"]]
  A <- kernel[["quadratic"]]
  precision <- A + 1 / normal[["sd"]]^2
  mean <- current + (kernel[["linear"]] + (normal[["mean"]] - current) / normal[["sd"]]^2) / precision
  baseline <- state[["baseline_log_q"]] - BayesTools::lpdf(state[["ordered_chart"]][["prior"]], current)
  log_L <- baseline + kernel[["linear"]] * (mean - current) - .5 * A * (mean - current)^2
  residual_scale <- 0
  if (!is.null(context) && .iwmde_ordered_source_residual_available(context, state)) {
    row <- matrix(state[["row"]], 1L, dimnames = list(NULL, names(state[["row"]])))
    replayed <- .iwmde_ordered_replay(context, row, state[["ordered_specs"]])
    prior_constant <- state[["prior_evaluator"]](replayed) + state[["ordered_kernel"]](replayed) -
      BayesTools::lpdf(state[["ordered_chart"]][["prior"]], current)
    evaluated <- .iwmde_ordered_source_gaussian_log(context, state, mean, prior_constant)
    log_L <- evaluated[["value"]]
    residual_scale <- evaluated[["max_log_scale"]]
  }
  log_mass <- log_L - .5 * ((mean - normal[["mean"]]) / normal[["sd"]])^2 -
    .5 * log(normal[["sd"]]^2 * precision)
  posterior_mass <- .iwmde_log_normal_interval_prob(normal[["lower"]], normal[["upper"]],
    mean, 1 / sqrt(precision))
  original_mass <- .iwmde_log_normal_interval_prob(normal[["lower"]], normal[["upper"]],
    normal[["mean"]], normal[["sd"]])
  mean_scale <- abs(current) + (abs(kernel[["linear"]]) +
    (abs(normal[["mean"]]) + abs(current)) / normal[["sd"]]^2) / precision
  residual_log_scale <- max(1, abs(log_L), abs(log_mass), abs(posterior_mass), abs(original_mass),
      residual_scale, .5 * ((mean - normal[["mean"]]) / normal[["sd"]])^2,
      abs(.5 * log(normal[["sd"]]^2 * precision)),
      if (is.null(context) || !.iwmde_ordered_source_residual_available(context, state))
        c(abs(kernel[["linear"]] * (mean - current)), abs(.5 * A * (mean - current)^2)) else
        (abs(current) + (abs(kernel[["linear"]]) + (abs(normal[["mean"]]) + abs(current)) / normal[["sd"]]^2) / precision) *
          abs(mean - normal[["mean"]]) / normal[["sd"]]^2)
  prior_scale <- .iwmde_normal_interval_conditioning_scale(normal[["lower"]], normal[["upper"]],
    normal[["mean"]], normal[["sd"]])
  posterior_scale <- .iwmde_normal_interval_conditioning_scale(normal[["lower"]], normal[["upper"]],
    mean, 1 / sqrt(precision), mean_scale)
  list(log_normalizer = log_mass + posterior_mass - original_mass,
    residual_log_mass = log_mass, posterior_interval_log_mass = posterior_mass,
    original_prior_interval_log_mass = original_mass, method = "completed_gaussian_source",
    residual_log_scale = residual_log_scale,
    max_log_scale = max(residual_log_scale, prior_scale, posterior_scale))
}


.iwmde_ordered_normalization_error <- function(parameter, rows, proof) {

  message <- paste0("Ordered density for '", parameter,
    "' was rejected by diagnostics: source integration error ",
    format(max(proof[["ratio_relative_error"]]), digits = 8L),
    ". Inspect the condition's 'proof' field or use 'density_method = \"KDE\"'.")
  stop(structure(list(message = message, call = NULL, parameter = parameter,
    rows = rows, stage = "ordered source normalization", proof = proof),
    class = c("RoBMA_ordered_normalization_error", "iwmde_ordered_normalization_error",
      "iwmde_ordinate_numerical_error", "error", "condition")))
}


.iwmde_ordered_log_terms <- function(context, parameter_spec, row_states,
                                     evaluation_values, control) {

  parameter <- parameter_spec[["parameter"]]
  n <- length(row_states)
  log_terms <- matrix(NA_real_, length(evaluation_values), n)
  proofs <- vector("list", n)
  keys <- vapply(row_states, function(state) {
    chart <- state[["ordered_chart"]]
    .iwmde_hash("ordered_source_group", chart[intersect(c("kind", "source", "key", "J", "indices",
      "term", "slice", "component", "prior"), names(chart))])
  }, character(1))
  for (key in unique(keys)) {
    positions <- which(keys == key)
    states <- row_states[positions]
    rows <- vapply(states, `[[`, integer(1), "row_index")
    prior <- states[[1L]][["ordered_chart"]][["prior"]]
    support <- .iwmde_prior_support(prior)
    transform <- .iwmde_parameter_transform(support)
    current <- vapply(states, function(state) state[["ordered_chart"]][["current"]], numeric(1))
    masses <- lapply(states, .iwmde_ordered_normal_source_mass, context = context)
    analytic <- all(!vapply(masses, is.null, logical(1)))
    certified <- if (!analytic) lapply(states, function(state) .iwmde_ordered_source_normalizer(context, state)) else NULL
    source_inner <- !analytic && all(!vapply(certified, is.null, logical(1)))
    resolution <- numeric(length(states))
    if (analytic) {
      log_Z <- vapply(masses, `[[`, numeric(1), "log_normalizer")
      discretization <- exterior <- rep(0, length(states))
      summary <- list(normalization_range = support, normalization_points = 0L)
    } else if (source_inner) {
      masses <- certified
      log_Z <- vapply(masses, `[[`, numeric(1), "log_normalizer")
      discretization <- vapply(masses, `[[`, numeric(1), "discretization_relative_error")
      exterior <- vapply(masses, `[[`, numeric(1), "exterior_relative_error")
      summary <- list(normalization_range = support, normalization_points = 0L)
    } else {
      starts <- c(current, BayesTools::quant(prior, c(.005, .995)))
      grid <- .iwmde_normalization_grid(starts, numeric(), support, transform,
        control[["normalization_points"]], control[["normalization_prob"]])
      if (is.null(grid)) {
        .iwmde_stop_construction_failure(context[["ordered_requested_method"]], parameter, rows,
          "source normalization grid", "the source support has no representable normalization grid")
      }
      normalization <- .iwmde_qcmde_normalization(context, parameter,
        display_grid = numeric(), normalization_grid = grid, transform = transform,
        normalization_prob = control[["normalization_prob"]], row_states = states,
        replacement = list(type = "ordered_source"), estimator_rows = rows,
        active_mass = 1, n_candidate_rows = length(states))
      summary <- normalization[["summary"]]
      log_Z <- summary[["log_normalizer"]]
      discretization <- abs(expm1(summary[["validation_log_normalizer"]] - log_Z))
      exterior <- rep(summary[["truncation"]][["max"]], length(states))
      resolution <- summary[["resolution"]][["relative_error"]]
      if (!is.numeric(resolution) || length(resolution) != length(states)) resolution <- rep(NA_real_, length(states))
    }
    n_observations <- nrow(context[["data"]][["outcome"]])
    if (is.null(n_observations)) n_observations <- 1L
    floor <- 64 * .Machine$double.eps * max(1L, n_observations) *
      pmax(1, abs(log_Z), abs(vapply(states, `[[`, numeric(1), "baseline_log_q")),
        if (analytic || source_inner) vapply(masses, `[[`, numeric(1), "max_log_scale") else 1)
    denominator_error <- discretization + exterior + resolution + expm1(floor)
    ratio_error <- (1 + expm1(floor)) / (1 - denominator_error) - 1
    proof <- list(internal_method = "ordered_source_conditional_cmde",
      working_distribution = "normalized_conditional_source", normalization_basis = "source_chart",
      error_status = "estimate", numerator_relative_error = expm1(floor),
      denominator_relative_error = denominator_error,
      ratio_relative_error = ratio_error,
      source_support = support, source_grid = summary[["normalization_range"]],
      source_grid_points = summary[["normalization_points"]],
      exterior_relative_error = exterior, floating_relative_error = expm1(floor),
      lattice_resolution_relative_error = resolution,
      candidate_rows = rows)
    if (any(!is.finite(ratio_error)) || any(denominator_error >= 1) || any(ratio_error > .025)) {
      .iwmde_ordered_normalization_error(parameter, rows, proof)
    }
    for (i in seq_along(states)) {
      state <- states[[i]]
      chart <- state[["ordered_chart"]]
      source <- (evaluation_values - chart[["offset"]]) / chart[["slope"]]
      log_q <- .iwmde_ordered_log_q_source(context, state, source)
      finite <- is.finite(log_q)
      evaluation_floor <- numeric(length(evaluation_values))
      evaluation_floor[finite] <- 64 * .Machine$double.eps * max(1L, n_observations) *
        pmax(1, abs(log_Z[[i]]), attr(log_q, "max_accumulated_logterm", exact = TRUE)[finite])
      evaluation_error <- (1 + expm1(evaluation_floor)) / (1 - denominator_error[[i]]) - 1
      if (any(!is.finite(evaluation_error[finite])) || any(evaluation_error[finite] > .025)) {
        current_proof <- proof
        current_proof[["ratio_relative_error"]] <- evaluation_error[finite]
        .iwmde_ordered_normalization_error(parameter, rows[[i]], current_proof)
      }
      log_terms[, positions[[i]]] <- log_q - log_Z[[i]] - log(abs(chart[["slope"]]))
      row_proof <- proof
      for (field in c("numerator_relative_error", "denominator_relative_error", "ratio_relative_error",
        "exterior_relative_error", "floating_relative_error", "lattice_resolution_relative_error")) row_proof[[field]] <- proof[[field]][[i]]
      row_proof[["candidate_rows"]] <- rows[[i]]
      row_proof[["chart"]] <- chart[intersect(c("kind", "source", "key", "J", "indices", "component"), names(chart))]
      row_proof[["log_normalizer"]] <- log_Z[[i]]
      row_proof[["evaluation_relative_error"]] <- evaluation_error
      row_proof[["ratio_relative_error"]] <- max(c(row_proof[["ratio_relative_error"]], evaluation_error[finite]))
      row_proof[["normalizer"]] <- if (analytic || source_inner) masses[[i]] else list(method = "checked_source_grid")
      proofs[[positions[[i]]]] <- row_proof
    }
  }
  if (anyNA(log_terms) || any(log_terms == Inf)) {
    .iwmde_stop_construction_failure(context[["ordered_requested_method"]], parameter,
      vapply(row_states, `[[`, integer(1), "row_index"), "normalized contribution validation",
      "the ordered log-term producer returned an undefined or infinite contribution")
  }
  list(log_terms = log_terms, proofs = proofs,
    internal_method = "ordered_source_conditional_cmde",
    working_distribution = "normalized_conditional_source", normalization_basis = "source_chart")
}


.iwmde_ordered_integral_error <- function(log_terms, x, density_terms,
                                          execution, active_mass,
                                          conditioned_rows = NULL,
                                          conditioned_chain_id = NULL) {

  if (length(x) < 2L) return(list(mcse = NA_real_, relative_mcse = NA_real_))
  dx <- diff(x)
  log_weights <- c(log(dx[[1L]]) - log(2),
    if (length(dx) > 1L) vapply(seq_len(length(dx) - 1L), function(i) {
      .iwmde_log_sum_exp(log(dx[c(i, i + 1L)])) - log(2)
    }, numeric(1)) else numeric(), log(tail(dx, 1L)) - log(2))
  integrals <- vapply(seq_len(ncol(log_terms)), function(i) {
    log(active_mass) + .iwmde_log_sum_exp(log_terms[, i] + log_weights)
  }, numeric(1))
  scale <- max(integrals)
  if (!is.finite(scale)) return(list(mcse = NA_real_, relative_mcse = NA_real_))
  contributions <- matrix(exp(integrals - scale), nrow = 1L)
  for (field in setdiff(names(attributes(density_terms[["contributions"]])), c("dim", "dimnames", "target", "finite_terms"))) {
    attr(contributions, field) <- attr(density_terms[["contributions"]], field, exact = TRUE)
  }
  attr(contributions, "target") <- mean(contributions)
  attr(contributions, "finite_terms") <- sum(is.finite(integrals))
  mcmc <- .iwmde_mcmc_contributions(contributions, execution[["sampling_population_rows"]],
    conditioned_rows, conditioned_chain_id, active_mass)
  error <- .iwmde_mixture_mcse(contributions, mcmc, density_terms[["active_mass_error"]], active_mass)
  log_mcse <- log(error[["mcse"]][[1L]]) + scale
  absolute <- exp(log_mcse)
  if (is.finite(log_mcse) && absolute == 0) absolute <- NA_real_
  list(mcse = absolute, relative_mcse = error[["relative_mcse"]][[1L]],
    log_mcse = log_mcse, log_scale = scale,
    evaluation_scale_underflow = any(is.finite(log_terms) &
      sweep(log_terms, 2L, integrals, "-") < log(.Machine$double.xmin * .Machine$double.eps)))
}


.iwmde_log_sum_exp <- function(values) {

  scale <- max(values)
  if (scale == -Inf) return(-Inf)
  scale + log(sum(exp(values - scale)))
}


.iwmde_ordered_density <- function(context, plan, execution, display_grid,
                                    density_output) {

  context[["ordered_requested_method"]] <- plan[["method"]]
  offset <- plan[["execution_spec"]][["ordered_point_plan"]]
  certificate <- .iwmde_ordered_point_certificate(context, plan[["execution_spec"]], execution[["row_states"]])
  context[["ordered_point_certificate"]] <- certificate
  omitted <- list()
  exact <- if (is.null(offset)) logical(length(display_grid)) else display_grid == offset[["value"]]
  unavailable_block <- any(exact) && identical(offset[["status"]], "block") && !certificate[["available"]]
  if (any(exact) && (offset[["status"]] %in% c("infinite", "capability") || unavailable_block)) {
    cause <- if (unavailable_block || identical(offset[["status"]], "infinite") && !certificate[["available"]]) certificate[["cause"]] else offset[["cause"]]
    detail <- if (unavailable_block || identical(offset[["status"]], "infinite") && !certificate[["available"]]) certificate[["detail"]] else offset[["detail"]]
    if (!density_output) .iwmde_ordered_stop(plan[["target"]][["parameter"]], cause, detail, offset[["value"]])
    omitted <- list(list(value = offset[["value"]], value_key = offset[["value_key"]], cause = cause, reason = detail))
    display_grid <- display_grid[!exact]
    if (!length(display_grid)) .iwmde_ordered_stop(plan[["target"]][["parameter"]], cause, detail, offset[["value"]])
  }
  if (!density_output && !certificate[["available"]]) {
    .iwmde_ordered_stop(plan[["target"]][["parameter"]], certificate[["cause"]], certificate[["detail"]])
  }
  exact <- if (is.null(offset)) logical(length(display_grid)) else display_grid == offset[["value"]]
  block_point <- any(exact) && identical(offset[["status"]], "block")
  if (block_point) {
    block <- .iwmde_ordered_block_log_terms(context, plan[["execution_spec"]], execution[["row_states"]], offset, plan[["control"]])
    if (all(exact)) {
      produced <- block
    } else {
      produced <- .iwmde_ordered_log_terms(context, plan[["execution_spec"]],
        execution[["row_states"]], display_grid[!exact], plan[["control"]])
      combined <- matrix(NA_real_, length(display_grid), length(execution[["row_states"]]))
      combined[!exact, ] <- produced[["log_terms"]]
      combined[exact, ] <- block[["log_terms"]]
      produced[["log_terms"]] <- combined
      produced[["source_subplan"]] <- produced[["proofs"]]
      produced[["block_subplan"]] <- block[["proofs"]]
      for (i in seq_along(produced[["proofs"]])) {
        produced[["proofs"]][[i]][["ratio_relative_error"]] <- max(
          produced[["proofs"]][[i]][["ratio_relative_error"]], block[["proofs"]][[i]][["ratio_relative_error"]])
      }
      produced[["internal_method"]] <- "ordered_source_with_offset_block_cmde"
      produced[["normalization_basis"]] <- "source_chart_with_offset_block"
    }
  } else produced <- .iwmde_ordered_log_terms(context, plan[["execution_spec"]],
    execution[["row_states"]], display_grid, plan[["control"]])
  log_terms <- produced[["log_terms"]]
  structural_zero <- if (is.null(offset) || !identical(offset[["status"]], "zero")) {
    rep(FALSE, length(display_grid))
  } else display_grid == offset[["value"]]
  log_terms[structural_zero, ] <- -Inf
  scales <- apply(log_terms, 1L, max)
  scales[scales == -Inf] <- 0
  scaled <- sweep(log_terms, 1L, scales, "-")
  active_mass <- plan[["rows"]][["active_mass"]]
  mixed <- plan[["rows"]][["point_mass_total"]] > 0
  density_terms <- .iwmde_density_aggregate(scaled, active_mass, ncol(log_terms),
    contribution_rows = execution[["active_rows"]],
    sampling_population_rows = execution[["sampling_population_rows"]],
    chain_id = execution[["chain_id"]], expected_chain_ids = execution[["expected_chain_ids"]],
    conditioned_rows = if (mixed) execution[["population_rows"]] else NULL,
    conditioned_chain_id = if (mixed) execution[["conditioned_chain_id"]] else NULL,
    evaluation_values = display_grid)
  uncertainty <- .iwmde_mixture_mcse(density_terms[["contributions"]],
    density_terms[["mcmc_contributions"]], density_terms[["active_mass_error"]], active_mass)
  constant <- .iwmde_ordered_constant_certificate(context, plan[["execution_spec"]], execution[["row_states"]])
  zero_variance <- is.finite(uncertainty[["mcse"]]) & uncertainty[["mcse"]] == 0 & !structural_zero
  if (any(zero_variance) && !constant) {
    for (field in c("mcse", "relative_mcse", "active_branch_mcse", "active_branch_relative_mcse",
      "active_mass_component_mcse", "ess")) uncertainty[[field]][zero_variance] <- NA_real_
    if (all(is.na(uncertainty[["ess"]]))) {
      uncertainty[["uncertainty_scope"]] <- "unavailable_insufficient_chain_batches"
      uncertainty[["uncertainty_status"]] <- "unavailable"
    } else uncertainty[["uncertainty_status"]] <- "partial"
    uncertainty[["uncertainty_reason"]] <- "zero empirical contribution variance lacks a declared dependency certificate"
  }
  if (constant) {
    for (field in c("mcse", "relative_mcse", "active_branch_mcse", "active_branch_relative_mcse",
      "active_mass_component_mcse")) uncertainty[[field]] <- rep(0, length(display_grid))
    uncertainty[["ess"]] <- rep(ncol(log_terms), length(display_grid))
    uncertainty[["uncertainty_scope"]] <- "declared_constant_conditional_contribution"
    uncertainty[["uncertainty_status"]] <- "available"
    uncertainty[["uncertainty_reason"]] <- NULL
  }
  log_y <- density_terms[["log_y"]] + scales
  y <- exp(log_y)
  if (any(is.infinite(y))) .iwmde_ordered_numeric_stop("ordinate_overflow",
    "the finite log ordinate cannot be represented on the output density scale", list(log_ordinate = log_y))
  underflow <- is.finite(log_y) & y == 0
  log_mcse <- log(uncertainty[["mcse"]]) + scales
  rescale <- function(value) {
    out <- exp(log(value) + scales)
    out[is.finite(log(value)) & out == 0] <- NA_real_
    out
  }
  per_row_error <- vapply(produced[["proofs"]], `[[`, numeric(1), "ratio_relative_error")
  positive <- exp(scaled)
  numerator <- rowSums(sweep(positive, 2L, per_row_error / (1 - per_row_error), "*"))
  denominator <- rowSums(sweep(positive, 2L, 1 / (1 + per_row_error), "*"))
  numerical_error <- rep(0, length(display_grid))
  available <- denominator > 0
  numerical_error[available] <- numerator[available] / denominator[available]
  resolution <- vapply(produced[["proofs"]], function(proof) {
    if (identical(proof[["internal_method"]], "ordered_offset_block_cmde")) return(0)
    value <- proof[["lattice_resolution_relative_error"]]
    if (is.numeric(value) && length(value) == 1L) value else NA_real_
  }, numeric(1))
  resolution_error <- .iwmde_qcmde_resolution_ordinate_error(log_terms, rep(0, ncol(log_terms)),
    list(relative_error = resolution))
  integral <- if (density_output) .iwmde_ordered_integral_error(log_terms, display_grid,
    density_terms, execution, active_mass,
    conditioned_rows = if (mixed) execution[["population_rows"]] else NULL,
    conditioned_chain_id = if (mixed) execution[["conditioned_chain_id"]] else NULL
  ) else list(mcse = NA_real_, relative_mcse = NA_real_)
  source_points <- max(vapply(produced[["proofs"]], `[[`, numeric(1), "source_grid_points"))
  source_ranges <- do.call(rbind, lapply(produced[["proofs"]], `[[`, "source_grid"))
  source_range <- if (source_points == 0L) c(NA_real_, NA_real_) else range(source_ranges)
  normalizers <- vapply(produced[["proofs"]], `[[`, numeric(1), "log_normalizer")
  list(x = display_grid, y = y, log_y = log_y,
    finite_terms = density_terms[["finite_terms"]], max_log_ratio = density_terms[["max_log_ratio"]],
    ess = uncertainty[["ess"]], max_weight_share = density_terms[["max_weight_share"]],
    mcse = rescale(uncertainty[["mcse"]]), relative_mcse = uncertainty[["relative_mcse"]],
    active_branch_mcse = rescale(uncertainty[["active_branch_mcse"]]),
    active_branch_relative_mcse = uncertainty[["active_branch_relative_mcse"]],
    active_mass_mcse = uncertainty[["active_mass_mcse"]],
    active_mass_relative_mcse = uncertainty[["active_mass_relative_mcse"]],
    active_mass_component_mcse = rescale(uncertainty[["active_mass_component_mcse"]]),
    mixture_mcse_type = uncertainty[["mixture_mcse_type"]],
    sampling_mcse = rescale(density_terms[["sampling_mcse"]]),
    sampling_relative_mcse = density_terms[["sampling_relative_mcse"]],
    sampling_fraction = density_terms[["sampling_fraction"]],
    sampling_uncertainty_type = density_terms[["sampling_uncertainty_type"]],
    mcmc_uncertainty_scope = uncertainty[["uncertainty_scope"]],
    mcmc_uncertainty_status = uncertainty[["uncertainty_status"]],
    mcmc_uncertainty_reason = uncertainty[["uncertainty_reason"]],
    log_normalizer = normalizers, pilot_log_normalizer = normalizers,
    conditional_normalization = list(rows = seq_along(normalizers),
      methods = vapply(produced[["proofs"]], function(proof) proof[["normalizer"]][["method"]], character(1))),
    pilot_gate_stopped = FALSE, pilot_bulk_ess = numeric(),
    n_candidate_rows = ncol(log_terms), n_evaluated_rows = ncol(log_terms),
    normalization_points = source_points, normalization_range = source_range,
    normalization_initial_points = source_points, normalization_initial_range = source_range,
    pilot_normalization_integral = active_mass, final_normalization_integral = active_mass,
    normalization_relative_error = max(numerical_error), normalization_scale = produced[["normalization_basis"]],
    normalization_mass_ratio = 1, normalization_truncation = max(vapply(produced[["proofs"]], `[[`, numeric(1), "exterior_relative_error")),
    normalization_truncation_status = "estimate", truncation_ordinate_bound = 0,
    lattice_resolution_relative_error = resolution_error,
    pilot_y = y, validation_y = y, ordinate_relative_change = numerical_error,
    ordinate_log_change = log1p(numerical_error), pilot_ordinate_relative_change = numerical_error,
    pilot_ordinate_log_change = log1p(numerical_error),
    max_normalizer_relative_change = max(per_row_error), max_quadrature_relative_change = max(per_row_error),
    p95_normalizer_relative_change = unname(stats::quantile(per_row_error, .95)),
    median_normalizer_relative_change = stats::median(per_row_error),
    normalization_refined_points = source_points, normalization_refined_range = source_range,
    n_rescued_normalizer = 0L, n_initial_dropped_normalizer = 0L, n_refinement_steps = 0L,
    integral_mcse = integral[["mcse"]], integral_relative_mcse = integral[["relative_mcse"]],
    batch_size = uncertainty[["batch_size"]], n_batches = uncertainty[["n_batches"]],
    estimator = plan[["method"]], weight_method = "normalized_conditional_source",
    internal_method = produced[["internal_method"]], working_distribution = produced[["working_distribution"]],
    normalization_basis = produced[["normalization_basis"]], ordered_proofs = produced[["proofs"]],
    ordered_source_subplan = produced[["source_subplan"]], ordered_block_subplan = if (all(exact) && block_point) produced[["proofs"]] else produced[["block_subplan"]],
    contribution_log_scale = scales, scaled_contributions = density_terms[["contributions"]],
    log_mcse = log_mcse, underflow = underflow, joint_integral_proof = integral,
    omitted_ordered_points = omitted, structural_zero = structural_zero,
    ordered_point_plan = offset, ordered_point_certificate = certificate,
    constant_dependency_certificate = constant)
}
