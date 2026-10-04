# Positive conditional block integration for declared regular ordered offsets.

.iwmde_ordered_numeric_stop <- function(cause, detail, proof = list()) {

  stop(structure(list(message = paste0("Ordered density was rejected by diagnostics: ",
    detail, ". Inspect 'density_diagnostics()' and the block integration proof."),
    call = NULL, cause = cause, proof = proof),
    class = c(paste0("RoBMA_ordered_", cause), "iwmde_ordinate_numerical_error", "error", "condition")))
}


.iwmde_ordered_poly_trim <- function(p) {

  while (length(p) > 1L && tail(p, 1L) == 0) p <- utils::head(p, -1L)
  p
}


.iwmde_ordered_poly_derivative <- function(p) {

  if (length(p) == 1L) 0 else p[-1L] * seq_len(length(p) - 1L)
}


.iwmde_ordered_poly_multiply <- function(a, b) {

  out <- numeric(length(a) + length(b) - 1L)
  for (i in seq_along(a)) for (j in seq_along(b)) out[[i + j - 1L]] <- out[[i + j - 1L]] + a[[i]] * b[[j]]
  out
}


.iwmde_ordered_poly_add <- function(...) {

  arguments <- list(...)
  out <- numeric(max(lengths(arguments)))
  for (argument in arguments) out[seq_along(argument)] <- out[seq_along(argument)] + argument
  .iwmde_ordered_poly_trim(out)
}


.iwmde_ordered_poly_value <- function(p, x) {

  out <- rep(p[[length(p)]], length(x))
  if (length(p) > 1L) for (i in (length(p) - 1L):1L) out <- out * x + p[[i]]
  out
}


.iwmde_ordered_poly_sign <- function(p, x) {

  q <- p / max(abs(p))
  n <- length(p) - 1L
  value <- .iwmde_ordered_poly_value(q, x)
  scale <- .iwmde_ordered_poly_value(abs(q), abs(x))
  gamma <- n * .Machine$double.eps / (1 - n * .Machine$double.eps)
  error <- 2 * gamma * scale
  list(value = value, error = error,
    sign = if (value > error) 1L else if (value < -error) -1L else 0L,
    margin = abs(value) - error)
}


.iwmde_ordered_poly_interval <- function(p, lower, upper) {

  q <- p / max(abs(p))
  coefficients <- cbind(q - 2 * .Machine$double.eps * abs(q), q + 2 * .Machine$double.eps * abs(q))
  out <- coefficients[nrow(coefficients), ]
  if (length(q) > 1L) for (i in (length(q) - 1L):1L) {
    products <- c(out[[1L]] * lower, out[[1L]] * upper, out[[2L]] * lower, out[[2L]] * upper)
    bound <- 4 * .Machine$double.eps * max(abs(products))
    out <- c(min(products) - bound, max(products) + bound) + coefficients[i, ]
    out <- out + c(-1, 1) * 4 * .Machine$double.eps * max(abs(out))
  }
  out
}


.iwmde_ordered_poly_roots <- function(p, depth = 0L) {

  p <- .iwmde_ordered_poly_trim(p)
  refuse <- function(reason, trace = list()) .iwmde_ordered_numeric_stop("root_isolation_unresolved", reason, trace)
  if (any(!is.finite(p)) || length(p) > 6L) refuse("the nominal root polynomial is invalid or has degree above five")
  if (all(p == 0)) return(list(status = "structural_constant", roots = matrix(numeric(), 0L, 2L),
    multiplicity = integer(), derivative = list(), trace = list()))
  k <- 0L
  q <- p
  while (length(q) > 1L && q[[1L]] == 0) {
    k <- k + 1L
    q <- q[-1L]
  }
  endpoint <- if (k) matrix(c(0, 0), 1L, 2L) else matrix(numeric(), 0L, 2L)
  if (length(q) == 1L) return(list(status = "accepted", roots = endpoint,
    multiplicity = if (k) k else integer(), derivative = list(), trace = list(), factor = k, q = q))
  derivative <- .iwmde_ordered_poly_roots(.iwmde_ordered_poly_derivative(q), depth + 1L)
  # Endpoint-touching nonstructural derivative bands remain critical. Removing
  # them would incorrectly accept the original close endpoint-root pair.
  critical <- derivative[["roots"]][!(derivative[["roots"]][, 1L] == 0 &
    derivative[["roots"]][, 2L] == 0), , drop = FALSE]
  trace <- list()
  if (nrow(critical)) for (i in seq_len(nrow(critical))) {
    interval <- .iwmde_ordered_poly_interval(q, critical[i, 1L], critical[i, 2L])
    trace[[length(trace) + 1L]] <- list(type = "critical", bracket = critical[i, ], p_interval = interval)
    if (interval[[1L]] <= 0 && interval[[2L]] >= 0) {
      refuse("a repeated or near-tangent root remains possible in a critical bracket", list(depth = depth, trace = trace, derivative = derivative))
    }
  }
  lower_cuts <- c(0, if (nrow(critical)) critical[, 2L] else numeric())
  upper_cuts <- c(if (nrow(critical)) critical[, 1L] else numeric(), 1)
  roots <- endpoint
  multiplicity <- if (k) k else integer()
  for (i in seq_along(lower_cuts)) {
    lower <- lower_cuts[[i]]
    upper <- upper_cuts[[i]]
    if (lower == upper) next
    left <- .iwmde_ordered_poly_sign(q, lower)
    right <- .iwmde_ordered_poly_sign(q, upper)
    trace[[length(trace) + 1L]] <- list(type = "segment", bracket = c(lower, upper), left = left, right = right)
    if (!left[["sign"]] || !right[["sign"]]) {
      refuse("a segment or endpoint sign is uncertified", list(depth = depth, trace = trace, derivative = derivative))
    }
    if (left[["sign"]] == right[["sign"]]) next
    method <- "bisection"
    iterations <- 0L
    repeat {
      tolerance <- max(1e-12, 64 * .Machine$double.eps * max(abs(c(lower, upper))))
      if (upper - lower <= tolerance) break
      midpoint <- (lower + upper) / 2
      middle <- .iwmde_ordered_poly_sign(q, midpoint)
      iterations <- iterations + 1L
      if (!middle[["sign"]]) {
        interval <- .iwmde_ordered_poly_interval(.iwmde_ordered_poly_derivative(q), lower, upper)
        minimum <- if (interval[[1L]] > 0) interval[[1L]] else if (interval[[2L]] < 0) -interval[[2L]] else 0
        minimum <- minimum * max(abs(.iwmde_ordered_poly_derivative(q))) / max(abs(q))
        if (minimum <= 0) refuse("the uncertain midpoint has no positive derivative lower bound",
          list(depth = depth, trace = trace, derivative = derivative, bracket = c(lower, upper),
            midpoint = middle, derivative_interval = interval))
        radius <- (abs(middle[["value"]]) + middle[["error"]]) / minimum
        candidate <- c(max(lower, midpoint - radius * (1 + 8 * .Machine$double.eps)),
          min(upper, midpoint + radius * (1 + 8 * .Machine$double.eps)))
        candidate <- candidate + c(-1, 1) * 8 * .Machine$double.eps * max(abs(candidate))
        candidate <- c(max(lower, candidate[[1L]]), min(upper, candidate[[2L]]))
        if (diff(candidate) > tolerance) refuse("the certified midpoint radius exceeds the width criterion",
          list(depth = depth, trace = trace, derivative = derivative, bracket = c(lower, upper), midpoint = middle, radius = radius))
        lower <- candidate[[1L]]
        upper <- candidate[[2L]]
        method <- "MVT"
        break
      }
      if (middle[["sign"]] == left[["sign"]]) {
        lower <- midpoint
        left <- middle
      } else {
        upper <- midpoint
        right <- middle
      }
      if (iterations > 200L) refuse("the root bisection work bound was exhausted", trace)
    }
    roots <- rbind(roots, c(lower, upper))
    multiplicity <- c(multiplicity, 1L)
    trace[[length(trace) + 1L]] <- list(type = "root", bracket = c(lower, upper), left = left, right = right,
      method = method, iterations = iterations)
  }
  list(status = "accepted", roots = roots, multiplicity = multiplicity,
    derivative = derivative, trace = trace, factor = k, q = q)
}


.iwmde_ordered_log_add <- function(a, b) {

  scale <- pmax(a, b)
  out <- scale + log(exp(a - scale) + exp(b - scale))
  out[scale == -Inf] <- -Inf
  out
}


.iwmde_ordered_log_difference <- function(a, b) {

  large <- pmax(a, b)
  small <- pmin(a, b)
  out <- large + log(-expm1(small - large))
  out[a == b] <- -Inf
  out
}


.iwmde_ordered_probability_prior <- function(prior) {

  distribution <- prior[["distribution"]]
  parameters <- prior[["parameters"]]
  support <- .iwmde_prior_support(prior)
  if (.iwmde_retained_location_dynamic_prior(prior) ||
      !distribution %in% c("normal", "gamma", "beta")) {
    .iwmde_ordered_stop("ordered block", "proper_cdf_prior", "the source has no supplied proper-CDF panel family")
  }
  base_cdf <- switch(distribution,
    normal = function(value, lower) stats::pnorm(value, parameters[["mean"]], parameters[["sd"]], lower.tail = lower, log.p = TRUE),
    gamma = function(value, lower) stats::pgamma(value, parameters[["shape"]], rate = parameters[["rate"]], lower.tail = lower, log.p = TRUE),
    beta = function(value, lower) stats::pbeta(value, parameters[["alpha"]], parameters[["beta"]], lower.tail = lower, log.p = TRUE))
  base_quantile <- switch(distribution,
    normal = function(value, lower) stats::qnorm(value, parameters[["mean"]], parameters[["sd"]], lower.tail = lower, log.p = TRUE),
    gamma = function(value, lower) stats::qgamma(value, parameters[["shape"]], rate = parameters[["rate"]], lower.tail = lower, log.p = TRUE),
    beta = function(value, lower) stats::qbeta(value, parameters[["alpha"]], parameters[["beta"]], lower.tail = lower, log.p = TRUE))
  raw_panel <- function(left, right) {
    if (any(is.na(left)) || any(is.na(right)) || any(right <= left)) {
      .iwmde_ordered_numeric_stop("probability_panel_unrepresentable", "a probability panel has invalid support", list(left = left, right = right))
    }
    lower <- base_cdf(right, TRUE) <= log(.5)
    lower[!lower & base_cdf(left, FALSE) > log(.5)] <- TRUE
    large <- ifelse(lower, base_cdf(right, TRUE), base_cdf(left, FALSE))
    small <- ifelse(lower, base_cdf(left, TRUE), base_cdf(right, FALSE))
    mass <- large + log(-expm1(small - large))
    if (any(!is.finite(mass))) {
      .iwmde_ordered_numeric_stop("probability_panel_unrepresentable", "a chosen support-tail interval cannot be represented",
        list(stage = "tail_probability_difference", left = left, right = right, large = large, small = small))
    }
    list(lower = left, upper = right, lower_tail = lower, log_base_mass = mass,
      log_anchor = ifelse(lower, base_cdf(left, TRUE), base_cdf(right, FALSE)))
  }
  normalization <- raw_panel(support[[1L]], support[[2L]])[["log_base_mass"]]
  panel <- function(left, right) {
    out <- raw_panel(left, right)
    out[["log_mass"]] <- out[["log_base_mass"]] - normalization
    out
  }
  conditional <- function(panel, fraction) {
    probability <- .iwmde_ordered_log_add(panel[["log_anchor"]], panel[["log_base_mass"]] +
      if (panel[["lower_tail"]]) log(fraction) else log1p(-fraction))
    if (any(probability >= 0) || any(!is.finite(probability))) {
      .iwmde_ordered_numeric_stop("probability_panel_unrepresentable", "a conditional node probability rounded outside open support",
        list(panel = panel, log_probability = probability))
    }
    out <- base_quantile(probability, panel[["lower_tail"]])
    if (any(!is.finite(out)) || any(out <= panel[["lower"]] | out >= panel[["upper"]])) {
      .iwmde_ordered_numeric_stop("probability_panel_unrepresentable", "a conditional quantile lies outside its open panel",
        list(panel = panel, node = out))
    }
    out
  }
  conditional_many <- function(panels, fraction) {
    n <- length(fraction)
    lower_tail <- rep(vapply(panels, `[[`, logical(1), "lower_tail"), each = n)
    anchor <- rep(vapply(panels, `[[`, numeric(1), "log_anchor"), each = n)
    mass <- rep(vapply(panels, `[[`, numeric(1), "log_base_mass"), each = n)
    fractions <- rep(fraction, length(panels))
    probability <- .iwmde_ordered_log_add(anchor, mass + ifelse(lower_tail, log(fractions), log1p(-fractions)))
    if (any(probability >= 0) || any(!is.finite(probability))) {
      .iwmde_ordered_numeric_stop("probability_panel_unrepresentable", "a conditional node probability rounded outside open support",
        list(panels = panels, log_probability = probability))
    }
    out <- numeric(length(probability))
    if (any(lower_tail)) out[lower_tail] <- base_quantile(probability[lower_tail], TRUE)
    if (any(!lower_tail)) out[!lower_tail] <- base_quantile(probability[!lower_tail], FALSE)
    lower <- rep(vapply(panels, `[[`, numeric(1), "lower"), each = n)
    upper <- rep(vapply(panels, `[[`, numeric(1), "upper"), each = n)
    if (any(!is.finite(out)) || any(out <= lower | out >= upper)) {
      .iwmde_ordered_numeric_stop("probability_panel_unrepresentable", "a conditional quantile lies outside its open panel",
        list(panels = panels, node = out))
    }
    out
  }
  quantile <- function(probabilities) {
    out <- numeric(length(probabilities))
    out[probabilities == 0] <- support[[1L]]
    out[probabilities == 1] <- support[[2L]]
    interior <- which(probabilities > 0 & probabilities < 1)
    if (length(interior)) out[interior] <- conditional(panel(support[[1L]], support[[2L]]), probabilities[interior])
    out
  }
  list(prior = prior, support = support, panel = panel, conditional = conditional,
    conditional_many = conditional_many, quantile = quantile,
    log_normalization = normalization)
}


.iwmde_ordered_probability_panels <- function(cuts, prior) {

  stored <- attr(cuts, "ordered_probability_panels", exact = TRUE)
  if (!is.null(stored) && identical(attr(cuts, "ordered_probability_prior", exact = TRUE), prior[["prior"]])) return(stored)

  lapply(seq_len(length(cuts) - 1L), function(i) prior[["panel"]](cuts[[i]], cuts[[i + 1L]]))
}


.iwmde_ordered_probability_refine <- function(cuts, prior) {

  panels <- .iwmde_ordered_probability_panels(cuts, prior)
  .iwmde_ordered_merge_cuts(cuts, vapply(panels, function(panel) prior[["conditional"]](panel, .5), numeric(1)),
    "partition_refinement", prior)
}


.iwmde_ordered_panel_rule <- function(cuts, prior, order, evaluate) {

  rule <- .gauss_legendre_nodes(order)
  fractions <- rule[["nodes"]]
  weights <- rule[["weights"]]
  panels <- .iwmde_ordered_probability_panels(cuts, prior)
  node_cache <- attr(cuts, "ordered_probability_nodes", exact = TRUE)
  node_key <- as.character(order)
  nodes <- if (is.environment(node_cache)) node_cache[[node_key]] else NULL
  if (is.null(nodes)) {
    nodes <- prior[["conditional_many"]](panels, fractions)
    if (is.environment(node_cache)) node_cache[[node_key]] <- nodes
  }
  result <- evaluate(nodes)
  values <- lapply(seq_along(panels), function(i) {
    take <- seq.int((i - 1L) * order + 1L, i * order)
    panel <- panels[[i]]
    list(log_value = .iwmde_log_sum_exp(result[["log_value"]][take] + log(weights)) + panel[["log_mass"]],
      log_inner_error = if (is.null(result[["log_error"]])) -Inf else
        .iwmde_log_sum_exp(result[["log_error"]][take] + log(weights)) + panel[["log_mass"]],
      log_inner_exterior = if (is.null(result[["log_exterior"]])) -Inf else
        .iwmde_log_sum_exp(result[["log_exterior"]][take] + log(weights)) + panel[["log_mass"]])
  })
  panel_values <- vapply(values, `[[`, numeric(1), "log_value")
  list(log_value = .iwmde_log_sum_exp(panel_values), panels = panel_values,
    log_inner_error = .iwmde_log_sum_exp(vapply(values, `[[`, numeric(1), "log_inner_error")),
    log_inner_exterior = .iwmde_log_sum_exp(vapply(values, `[[`, numeric(1), "log_inner_exterior")))
}


.iwmde_ordered_role_cuts <- function(records, prior, history = records,
                                      coalesced = list(), recompute = NULL) {

  if (is.null(records[["index"]])) records[["index"]] <- seq_len(nrow(records))
  history <- history
  noncritical <- c("fixed_prior_probability", "roundtrip_quantile", "source_geometry")
  repeat {
    cuts <- sort(unique(records[["value"]]))
    failure <- NULL
    panels <- vector("list", length(cuts) - 1L)
    for (i in seq_len(length(cuts) - 1L)) {
      checked <- tryCatch(prior[["panel"]](cuts[[i]], cuts[[i + 1L]]), error = identity)
      if (inherits(checked, "error")) {
        failure <- list(index = i, error = checked)
        break
      }
      panels[[i]] <- checked
    }
    if (is.null(failure)) break
    i <- failure[["index"]]
    error <- failure[["error"]]
    left <- records[["role"]][records[["value"]] == cuts[[i]]]
    right <- records[["role"]][records[["value"]] == cuts[[i + 1L]]]
    left_ok <- length(left) && all(left %in% noncritical)
    right_ok <- length(right) && all(right %in% noncritical)
    if (!inherits(error, "RoBMA_ordered_probability_panel_unrepresentable") ||
        !identical(error[["proof"]][["stage"]], "tail_probability_difference") || (!left_ok && !right_ok)) {
      error[["proof"]][["anchor_roles"]] <- list(left = left, right = right)
      stop(error)
    }
    candidates <- c(if (left_ok) i, if (right_ok) i + 1L)
    ranks <- vapply(candidates, function(j) min(match(records[["role"]][records[["value"]] == cuts[[j]]], noncritical)), integer(1))
    indices <- vapply(candidates, function(j) min(records[["index"]][records[["value"]] == cuts[[j]]]), integer(1))
    remove <- candidates[order(ranks, indices)[[1L]]]
    omitted <- cuts[[remove]]
    retained <- cuts[[if (remove == i) i + 1L else i]]
    before <- cuts
    records <- records[records[["value"]] != omitted, , drop = FALSE]
    if (!identical(range(records[["value"]]), prior[["support"]])) {
      .iwmde_ordered_numeric_stop("probability_panel_unrepresentable", "bookkeeping coalescence changed the full prior support")
    }
    if (!is.null(recompute)) recompute()
    coalesced[[length(coalesced) + 1L]] <- list(omitted = omitted, retained = retained,
      before = before, after = sort(unique(records[["value"]])), left_roles = left, right_roles = right,
      original_anchors = history, geometry_recomputed = !is.null(recompute), reason = error[["cause"]])
  }
  attr(cuts, "ordered_anchor_records") <- records
  attr(cuts, "ordered_anchor_history") <- history
  attr(cuts, "ordered_coalesced_cuts") <- coalesced
  attr(cuts, "ordered_probability_panels") <- panels
  attr(cuts, "ordered_probability_prior") <- prior[["prior"]]
  attr(cuts, "ordered_probability_nodes") <- new.env(parent = emptyenv())
  cuts
}


.iwmde_ordered_base_cuts <- function(prior, anchors = numeric(), roles = "required_geometry",
                                      recompute = NULL) {

  probabilities <- c(0, .001, .01, .1, .25, .5, .75, .9, .99, .999, 1)
  support <- prior[["support"]]
  anchors <- vapply(anchors, function(value) max(support[[1L]], min(support[[2L]], value)), numeric(1))
  records <- rbind(data.frame(value = prior[["quantile"]](probabilities), role = "fixed_prior_probability"),
    data.frame(value = support, role = "support_endpoint"),
    data.frame(value = anchors, role = rep(roles, length.out = length(anchors))))
  .iwmde_ordered_role_cuts(records, prior, recompute = recompute)
}


.iwmde_ordered_merge_cuts <- function(cuts, anchors, roles, prior, recompute = NULL) {

  records <- attr(cuts, "ordered_anchor_records", exact = TRUE)
  history <- attr(cuts, "ordered_anchor_history", exact = TRUE)
  added <- data.frame(value = anchors, role = rep(roles, length.out = length(anchors)),
    index = max(records[["index"]]) + seq_along(anchors))
  .iwmde_ordered_role_cuts(rbind(records, added), prior, history = rbind(history, added),
    coalesced = attr(cuts, "ordered_coalesced_cuts", exact = TRUE), recompute = recompute)
}


.iwmde_ordered_block_samples <- function(context, block, total, share) {

  n <- max(length(total), length(share))
  total <- rep_len(total, n)
  share <- rep_len(share, n)
  row <- block[["state"]][["row"]]
  samples <- matrix(rep(row, each = n), nrow = n, dimnames = list(NULL, names(row)))
  chart <- block[["chart"]]
  samples[, block[["source"]]] <- total
  samples[, chart[["gamma_coordinates"]][chart[["J"]]]] <- outer(chart[["H"]] * share, chart[["r_J"]])
  samples[, chart[["gamma_coordinates"]][chart[["K"]]]] <- outer(chart[["H"]] * (1 - share), chart[["r_K"]])
  .iwmde_ordered_replay(context, samples, block[["state"]][["ordered_specs"]])
}


.iwmde_ordered_selection_partition <- function(context, parameter_spec, state) {

  data <- context[["data"]]
  if (!isTRUE(state[["active_setup"]][["is_weightfunction"]]) ||
      .is_data_scale(data) || .is_data_random(data) || .is_data_multilevel(data) ||
      .is_data_known_v(data) || .iwmde_context_uses_local_likelihood(context)) return(NULL)
  certificate <- context[["ordered_point_certificate"]]
  if (is.null(certificate)) certificate <- .iwmde_ordered_point_certificate(context, parameter_spec, list(state))
  if (!isTRUE(certificate[["available"]]) ||
      !identical(certificate[["kind"]], "independent_numeric_gaussian_core")) return(NULL)
  selection <- state[["active_setup"]][["priors"]][["outcome"]][["bias"]]
  if (!BayesTools::is.prior.weightfunction(selection) ||
      !identical(selection[["weights"]][["type"]], "fixed") ||
      !is.numeric(selection[["weights"]][["omega"]]) ||
      any(!is.finite(selection[["weights"]][["omega"]])) ||
      any(selection[["weights"]][["omega"]] <= 0)) return(NULL)
  plan <- .data_selection_execution_plan(data)
  model <- .data_selection_model(data)
  branches <- model[["branches"]][model[["active_branches"]]]
  K <- nrow(data[["outcome"]])
  if (!length(plan[["block_methods"]]) || any(plan[["block_methods"]] != "singleton") ||
      length(plan[["dependent_blocks"]]) || any(plan[["block_sizes"]] != 1L) ||
      length(plan[["sampling"]][["diagonal"]]) != K ||
      !is.matrix(plan[["sampling"]][["loading"]]) || ncol(plan[["sampling"]][["loading"]]) != 0L ||
      !length(branches) || !all(vapply(branches, function(branch) identical(branch[["weight_rule"]], "product"), logical(1))) ||
      .selection_retains_sampling(data) || .selection_retains_any_random(data)) return(NULL)
  design <- BayesTools::JAGS_formula_design(context[["object"]][["fit"]], parameter = "mu")
  if (!is.matrix(design[["model_matrix"]]) || nrow(design[["model_matrix"]]) != K ||
      !is.null(design[["formula_scale"]]) || length(design[["transformed_terms"]]) ||
      length(design[["expression_specs"]]) || length(design[["expression_data"]]) ||
      length(design[["random_effects"]]) || length(design[["random_allocations"]])) return(NULL)
  selection_spec <- .selection_spec(state[["active_setup"]][["priors"]], data[["outcome"]][["yi"]],
    data[["outcome"]][["sei"]], .data_effect_direction(data), signed_data = FALSE)
  if (is.null(selection_spec) || !identical(selection_spec[["mode"]], "step") ||
      isTRUE(selection_spec[["has_phack"]])) return(NULL)
  # Row labels and powers remain original identities. Equality uses every
  # persisted observation and mean-operator input, with exact serialized bits.
  row_inputs <- list(outcome = data[["outcome"]][setdiff(names(data[["outcome"]]), c("slab", "weights"))],
    mods = data[["mods"]], location = data[["location"]], model_frame = design[["model_frame"]],
    source_data = design[["source_data"]], model_matrix = design[["model_matrix"]],
    sampling_variance = plan[["sampling"]][["diagonal"]], obs_bin = selection_spec[["obs_bin"]])
  if (!all(vapply(row_inputs, function(x) {
    is.null(x) || if (is.data.frame(x) || is.matrix(x)) nrow(x) == K else length(x) == K
  }, logical(1)))) return(NULL)
  row_value <- function(x, i) {
    if (is.null(x)) return(NULL)
    if (is.data.frame(x)) return(lapply(x, function(column) unname(column[i])))
    if (is.matrix(x)) return(unname(x[i, , drop = TRUE]))
    unname(x[i])
  }
  rows <- lapply(seq_len(K), function(i) serialize(lapply(row_inputs, row_value, i = i), NULL, version = 2L))
  keys <- vapply(rows, rlang::hash, character(1))
  representatives <- which(!duplicated(keys))
  group_index <- match(keys, keys[representatives])
  # A hash selects candidates only; the original bytes prove equality.
  if (length(representatives) == K || !all(vapply(seq_len(K), function(i) {
    identical(rows[[i]], rows[[representatives[group_index[[i]]]]])
  }, logical(1)))) return(NULL)
  list(original_indices = seq_len(K), group_index = group_index,
    representatives = representatives, multiplicity = tabulate(group_index, length(representatives)),
    selection_spec = selection_spec, certificate = certificate,
    structure_identity = .iwmde_hash("ordered_selection_structure", list(data = data, design = design,
      model = model, plan = plan, selection = selection, recipes = parameter_spec[["ordered_specs"]])))
}


.iwmde_ordered_selection_grouping <- function(block, partition) {

  if (is.null(partition)) return(NULL)
  representatives <- partition[["representatives"]]
  group_index <- partition[["group_index"]]
  exact <- function(x) is.numeric(x) && length(x) == length(group_index) &&
    identical(serialize(unname(x), NULL, version = 2L),
      serialize(unname(x[representatives][group_index]), NULL, version = 2L))
  fields <- c("mu0", "mu_a", "mu_b", "mu_c", "variance", "selection_sei")
  if (!all(vapply(block[fields], exact, logical(1))) || !exact(block[["setup"]][["yi"]])) return(NULL)
  selection_context <- block[["selection_context"]]
  expected <- partition[["selection_spec"]]
  fixed_fields <- setdiff(names(expected), "native_cache")
  if (!identical(serialize(selection_context[fixed_fields], NULL, version = 2L),
      serialize(expected[fixed_fields], NULL, version = 2L)) ||
      !exact(selection_context[["obs_bin"]]) || !exact(selection_context[["sei"]]) ||
      !exact(selection_context[["yi"]]) || !identical(as.numeric(selection_context[["omega"]]),
        .selection_fixed_omega_branch(block[["state"]][["active_setup"]][["priors"]][["outcome"]][["bias"]], expected[["p_cuts"]])) ||
      !identical(selection_context[["alpha"]], 0) || !identical(selection_context[["phack_kind"]], 0L) ||
      any(selection_context[["kernel_mode"]] != SELKERNEL_STEP)) return(NULL)
  powers <- vapply(seq_along(representatives), function(group) sum(block[["weights"]][group_index == group]), numeric(1))
  map <- partition[c("original_indices", "group_index", "representatives", "multiplicity", "structure_identity", "certificate")]
  map[["original_powers"]] <- block[["weights"]]
  map[["powers"]] <- powers
  map[["group_count"]] <- length(representatives)
  map[["map_hash"]] <- .iwmde_hash("ordered_selection_partition", map)
  map[["input_identity"]] <- .iwmde_hash("ordered_selection_inputs", list(
    means = block[fields], yi = block[["setup"]][["yi"]],
    selection = selection_context[setdiff(names(selection_context), "native_cache")]))
  map
}


.iwmde_ordered_selection_observed_constant <- function(context, state, setup, selection_context) {

  data <- context[["data"]]
  selection <- state[["active_setup"]][["priors"]][["outcome"]][["bias"]]
  if (!BayesTools::is.prior.weightfunction(selection) ||
      !identical(selection[["weights"]][["type"]], "fixed") ||
      !is.numeric(selection[["weights"]][["omega"]]) ||
      any(!is.finite(selection[["weights"]][["omega"]])) || any(selection[["weights"]][["omega"]] <= 0) ||
      .is_data_random(data) || .is_data_multilevel(data) || .is_data_known_v(data) ||
      .iwmde_context_uses_local_likelihood(context) ||
      .selection_retains_sampling(data) || .selection_retains_any_random(data)) return(NULL)
  plan <- .data_selection_execution_plan(data)
  model <- .data_selection_model(data)
  branches <- model[["branches"]][model[["active_branches"]]]
  K <- length(setup[["yi"]])
  if (!length(plan[["block_methods"]]) || any(plan[["block_methods"]] != "singleton") ||
      length(plan[["dependent_blocks"]]) || any(plan[["block_sizes"]] != 1L) ||
      length(plan[["sampling"]][["diagonal"]]) != K ||
      !is.matrix(plan[["sampling"]][["loading"]]) || ncol(plan[["sampling"]][["loading"]]) != 0L ||
      !length(branches) || !all(vapply(branches, function(branch) identical(branch[["weight_rule"]], "product"), logical(1)))) return(NULL)
  expected <- .selection_spec(state[["active_setup"]][["priors"]], setup[["yi"]],
    setup[["selection_sei"]], .data_effect_direction(data), signed_data = FALSE)
  if (is.null(expected) || !identical(expected[["mode"]], "step") || isTRUE(expected[["has_phack"]])) return(NULL)
  fields <- setdiff(names(expected), "native_cache")
  omega <- as.numeric(selection_context[["omega"]])
  if (!identical(serialize(selection_context[fields], NULL, version = 2L),
      serialize(expected[fields], NULL, version = 2L)) ||
      !identical(omega, .selection_fixed_omega_branch(selection, expected[["p_cuts"]])) ||
      any(!is.finite(omega)) || any(omega <= 0) ||
      !identical(selection_context[["alpha"]], 0) || !identical(selection_context[["phack_kind"]], 0L) ||
      any(selection_context[["kernel_mode"]] != SELKERNEL_STEP) ||
      !identical(as.numeric(selection_context[["yi"]]), as.numeric(setup[["yi"]])) ||
      !identical(as.numeric(selection_context[["sei"]]), as.numeric(setup[["selection_sei"]])) ||
      !identical(as.numeric(setup[["sei"]]), as.numeric(setup[["selection_sei"]])) ||
      !is.matrix(setup[["mu"]]) || nrow(setup[["mu"]]) != 1L || ncol(setup[["mu"]]) != K ||
      any(!is.finite(setup[["mu"]]))) {
    .iwmde_ordered_numeric_stop("selection_route", "the fixed step-selection inputs do not match their declared context")
  }
  # On this declared route .outcome_pdf_sum.selnorm() uses the same mean and
  # total SD for numerator and normalization; retained-source routes are excluded.
  sd <- .root_sum_squares(setup[["tau_within"]], matrix(setup[["sei"]], 1L))
  if (length(sd) != K || any(!is.finite(sd)) || any(sd <= 0)) {
    .iwmde_ordered_numeric_stop("selection_route", "the fixed step-selection covariance is unrepresentable")
  }
  weights <- if (is.null(setup[["weights"]])) rep(1, K) else setup[["weights"]]
  sum(weights * log(omega[selection_context[["obs_bin"]]]))
}


.iwmde_ordered_block_geometry <- function(context, parameter_spec, state,
                                         observation_partition = .iwmde_ordered_selection_partition(context, parameter_spec, state)) {

  point <- parameter_spec[["ordered_point_plan"]]
  descriptor <- point[["descriptor"]]
  if (is.null(descriptor) || !isTRUE(descriptor[["pure"]])) {
    .iwmde_ordered_stop(parameter_spec[["parameter"]], "offset_block_capability_unavailable", "the block has no declared pure grouped target")
  }
  spec <- parameter_spec[["ordered_specs"]][[descriptor[["term"]]]]
  component <- .iwmde_ordered_total_component(spec, state[["row"]])
  prior <- component[["prior"]]
  if (!prior[["distribution"]] %in% c("normal", "gamma") || .iwmde_retained_location_dynamic_prior(prior) ||
      .is_data_known_v(context[["data"]]) ||
      .is_data_random(context[["data"]]) || .is_data_multilevel(context[["data"]])) {
    .iwmde_ordered_stop(parameter_spec[["parameter"]], "offset_block_capability_unavailable",
      "the source prior or covariance is outside the supplied Gaussian block implementation")
  }
  chart <- .iwmde_ordered_allocation_chart(context, parameter_spec, parameter_spec[["ordered_projection"]],
    state[["row"]], state[["row_index"]], descriptor[["record"]], J = descriptor[["J"]], component = component[["component"]])
  block <- list(state = state, chart = chart, source = component[["coordinates"]][[descriptor[["slice"]]]],
    prior = prior, beta = chart[["prior"]], descriptor = descriptor,
    shift = .iwmde_ordered_affine_shift(parameter_spec, state[["row"]]),
    counts = new.env(parent = emptyenv()))
  block[["counts"]][["candidate_rows"]] <- 0L
  block[["counts"]][["native_nodes"]] <- 0L
  block[["counts"]][["inner_calls"]] <- 0L
  block[["counts"]][["outer_nodes"]] <- 0L
  block[["counts"]][["maximum_chunk"]] <- 0L
  anchor_total <- BayesTools::quant(prior, .5)
  anchor <- .iwmde_ordered_block_samples(context, block, anchor_total, .5)
  setup <- .iwmde_log_lik_posterior_setup_active_branch(context, anchor, state[["active_setup"]], "estimate")
  if (!is.matrix(setup[["mu"]]) || nrow(setup[["mu"]]) != 1L) {
    .iwmde_ordered_stop(parameter_spec[["parameter"]], "offset_block_capability_unavailable", "the Gaussian predictor setup is unavailable")
  }
  K <- ncol(setup[["mu"]])
  corners <- .iwmde_ordered_block_samples(context, block, c(0, 0), c(0, 1))
  directions_a <- directions_b <- directions_c <- numeric()
  fit <- context[["object"]][["fit"]]
  for (owner in state[["ordered_specs"]]) {
    for (coefficient in owner[["coefficient_names"]]) {
      projection <- BayesTools::JAGS_ordered_parameter_spec(fit, parameters = names(state[["ordered_specs"]]),
        weights = stats::setNames(1, coefficient), draws = corners)
      total_weights <- projection[["total_contraction_weights"]][[descriptor[["term"]]]][, descriptor[["slice"]]]
      if (total_weights[[1L]] != 0) directions_a[[coefficient]] <- total_weights[[1L]]
      difference <- total_weights[[2L]] - total_weights[[1L]]
      if (difference != 0) directions_b[[coefficient]] <- difference
      contractions <- projection[["allocation_contractions"]][[chart[["key"]]]]
      if (!is.null(contractions)) {
        derivative <- sum(chart[["r_J"]] * contractions[1L, chart[["J"]]]) -
          sum(chart[["r_K"]] * contractions[1L, chart[["K"]]])
        if (derivative != 0) directions_c[[coefficient]] <- derivative
      }
    }
  }
  basis <- function(directions) {
    if (!length(directions)) return(matrix(0, 1L, K))
    result <- BayesTools::JAGS_formula_predictor_basis(fit, directions, posterior_samples = corners[1L, , drop = FALSE])
    if (!identical(result[["status"]], "affine") || !identical(result[["coordinate"]], "identity") ||
        !identical(result[["parameter"]], "mu")) {
      .iwmde_ordered_stop(parameter_spec[["parameter"]], "offset_block_capability_unavailable", "the mean is not a declared bilinear Gaussian predictor")
    }
    result[["basis"]]
  }
  a <- basis(directions_a)
  b <- basis(directions_b)
  c <- basis(directions_c)
  mu0 <- setup[["mu"]] - anchor_total * (a + .5 * b) - .5 * c
  baseline_setup <- setup
  baseline_setup[["mu"]] <- mu0
  change <- function(direction) .iwmde_normal_location_likelihood_change(context, baseline_setup, list(mu_basis = direction))
  A_a <- change(a)
  A_b <- change(b)
  A_c <- change(c)
  cross_ab <- (change(a + b)[["quadratic"]] - A_a[["quadratic"]] - A_b[["quadratic"]]) / 2
  cross_ac <- (change(a + c)[["quadratic"]] - A_a[["quadratic"]] - A_c[["quadratic"]]) / 2
  cross_bc <- (change(b + c)[["quadratic"]] - A_b[["quadratic"]] - A_c[["quadratic"]]) / 2
  variance <- as.numeric(setup[["tau_within"]])^2 + setup[["sei"]]^2
  if (length(variance) != K || any(!is.finite(variance)) || any(variance <= 0)) {
    .iwmde_ordered_stop(parameter_spec[["parameter"]], "offset_block_capability_unavailable", "the invariant diagonal Gaussian covariance is unavailable")
  }
  weights <- if (is.null(setup[["weights"]])) rep(1, K) else setup[["weights"]]
  log_L0 <- sum(weights * stats::dnorm(setup[["yi"]], as.numeric(mu0), sqrt(variance), log = TRUE))
  block[["gaussian"]] <- list(A = c(A_a[["quadratic"]], 2 * cross_ab, A_b[["quadratic"]]),
    H = c(A_a[["linear"]], A_b[["linear"]] - cross_ac, -cross_bc),
    C = c(log_L0, A_c[["linear"]], -.5 * A_c[["quadratic"]]))
  block[["mu0"]] <- as.numeric(mu0)
  block[["mu_a"]] <- as.numeric(a)
  block[["mu_b"]] <- as.numeric(b)
  block[["mu_c"]] <- as.numeric(c)
  block[["setup"]] <- setup
  block[["variance"]] <- variance
  block[["weights"]] <- weights
  block[["selected"]] <- isTRUE(state[["active_setup"]][["is_weightfunction"]])
  block[["log_correction_lower"]] <- block[["log_correction_upper"]] <- 0
  if (block[["selected"]]) {
    selection_context <- .iwmde_selection_context_active_branch(context, state[["active_setup"]], setup[["posterior_samples"]])
    if (is.null(selection_context) || any(selection_context[["kernel_mode"]] != SELKERNEL_STEP)) {
      .iwmde_ordered_stop(parameter_spec[["parameter"]], "selection_route", "the selected block is outside the positive step-partition allowlist")
    }
    log_norm <- .selection_step_log_norm_matrix(setup[["mu"]], matrix(sqrt(variance), 1L),
      setup[["selection_sei"]], selection_context)
    original_L <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(context, anchor, state[["active_setup"]], "estimate")
    gaussian_L <- sum(weights * stats::dnorm(setup[["yi"]], as.numeric(setup[["mu"]]), sqrt(variance), log = TRUE))
    observed_constant <- .iwmde_ordered_selection_observed_constant(context, state, setup, selection_context)
    if (is.null(observed_constant)) observed_constant <- original_L - gaussian_L + sum(weights * log_norm)
    block[["selection_context"]] <- selection_context
    block[["selection_sei"]] <- setup[["selection_sei"]]
    block[["selection_observed_constant"]] <- observed_constant
    omega <- state[["active_setup"]][["priors"]][["outcome"]][["bias"]][["weights"]][["omega"]]
    block[["log_correction_lower"]] <- observed_constant - sum(weights) * log(max(omega))
    block[["log_correction_upper"]] <- observed_constant - sum(weights) * log(min(omega))
    block[["counts"]][["native_nodes"]] <- 2L
    grouping <- .iwmde_ordered_selection_grouping(block, observation_partition)
    block[["selection_grouping"]] <- grouping
    if (!is.null(grouping)) {
      block[["normalizer_context"]] <- BayesTools::selection_context_subset_observations(selection_context, grouping[["representatives"]])
    }
  }
  block[["normal"]] <- .iwmde_qcmde_normal_prior_parameters(prior)
  block[["normal_probability"]] <- .iwmde_ordered_probability_prior(prior)
  block[["beta_probability"]] <- .iwmde_ordered_probability_prior(block[["beta"]])
  block[["closed_total_numerator"]] <- !isTRUE(descriptor[["shifted"]]) && all(vapply(state[["ordered_specs"]], function(consumer) {
    consumes <- chart[["key"]] %in% names(consumer[["allocations"]])
    if (!consumes) return(TRUE)
    if (identical(consumer[["parameter"]], descriptor[["term"]]) && length(consumer[["total_names"]]) == 1L) return(TRUE)
    total <- consumer[["total_prior"]]
    BayesTools::is.prior.point(total) && is.numeric(total[["parameters"]][["location"]]) && total[["parameters"]][["location"]] == 0
  }, logical(1)))
  block
}


.iwmde_ordered_gaussian_descriptor <- function(block, share) {

  anchor <- block[["gaussian_anchor"]]
  if (!is.null(anchor)) return(anchor)
  # Polynomial-only descriptions retain their established geometry. Production
  # bilinear blocks carry the complete original observations below.
  if (is.null(block[["setup"]])) {
    polynomial <- block[["gaussian"]]
    return(list(A = .iwmde_ordered_poly_value(polynomial[["A"]], share), t_ref = 0,
      H_ref = .iwmde_ordered_poly_value(polynomial[["H"]], share),
      C_ref = .iwmde_ordered_poly_value(polynomial[["C"]], share)))
  }
  a <- block[["mu0"]] + share * block[["mu_c"]]
  b <- block[["mu_a"]] + share * block[["mu_b"]]
  yi <- block[["setup"]][["yi"]]
  p <- block[["weights"]]
  v <- block[["variance"]]
  A <- sum(p * b^2 / v)
  H <- sum(p * b * (yi - a) / v)
  if (!is.finite(A) || A < 0 || !is.finite(H) ||
      (A == 0 && any(p > 0 & b != 0))) {
    .iwmde_ordered_numeric_stop("gaussian_completion_unrepresentable", "the Gaussian reference precision is unrepresentable",
      list(share = share, A = A, H = H))
  }
  t_ref <- if (A > 0) H / A else 0
  mean_ref <- a + b * t_ref
  residual <- yi - mean_ref
  H_ref <- sum(p * b * residual / v)
  pieces <- p * stats::dnorm(yi, mean_ref, sqrt(v), log = TRUE)
  C_ref <- sum(pieces)
  if (any(!is.finite(c(t_ref, mean_ref, H_ref, C_ref)))) {
    .iwmde_ordered_numeric_stop("gaussian_completion_unrepresentable", "the Gaussian reference is unrepresentable",
      list(share = share, A = A, H = H, t_ref = t_ref, H_ref = H_ref, C_ref = C_ref))
  }
  list(A = A, t_ref = t_ref, H_ref = H_ref, C_ref = C_ref, a = a, b = b,
    piece_scale = sum(abs(pieces)), residual_scale = .iwmde_ordered_gaussian_residual_scale(block, mean_ref,
      abs(block[["mu0"]]) + abs(share * block[["mu_c"]]) + abs(b * t_ref)))
}


.iwmde_ordered_gaussian_residual_scale <- function(block, mean, mean_scale) {

  yi <- block[["setup"]][["yi"]]
  residual <- yi - mean
  epsilon <- 64 * .Machine$double.eps * length(block[["weights"]])
  sum(block[["weights"]] * (abs(residual) * (abs(yi) + mean_scale) / block[["variance"]] +
    epsilon * (abs(yi) + mean_scale)^2 / block[["variance"]]))
}


.iwmde_ordered_gaussian_descriptors <- function(block, share) {

  keys <- .iwmde_key_number(share)
  # Descriptors reuse only identical IEEE source values within this call.
  keys[share == 0 & 1 / share < 0] <- "-0"
  representatives <- which(!duplicated(keys))
  list(descriptors = lapply(share[representatives], function(s) .iwmde_ordered_gaussian_descriptor(block, s)),
    index = match(keys, keys[representatives]))
}


.iwmde_ordered_block_gaussian_log <- function(block, total, share) {

  n <- max(length(total), length(share))
  total <- rep_len(total, n)
  share <- rep_len(share, n)
  parts <- .iwmde_ordered_gaussian_descriptors(block, share)
  out <- numeric(n)
  scale <- 1
  for (i in seq_along(parts[["descriptors"]])) {
    descriptor <- parts[["descriptors"]][[i]]
    take <- which(parts[["index"]] == i)
    if (!is.null(descriptor[["evaluate"]])) {
      evaluated <- descriptor[["evaluate"]](total[take])
      out[take] <- evaluated[["value"]]
      scale <- max(scale, evaluated[["max_log_scale"]])
      next
    }
    delta <- total[take] - descriptor[["t_ref"]]
    linear <- descriptor[["H_ref"]] * delta
    quadratic <- .5 * descriptor[["A"]] * delta^2
    out[take] <- descriptor[["C_ref"]] + linear - quadratic
    scale <- max(scale, abs(out[take]), abs(linear), abs(quadratic),
      descriptor[["piece_scale"]], descriptor[["residual_scale"]])
  }
  if (!is.null(block[["counts"]]) && (!is.null(block[["setup"]]) || !is.null(block[["gaussian_anchor"]]))) {
    block[["counts"]][["max_log_scale"]] <- max(block[["counts"]][["max_log_scale"]], scale)
  }
  out
}


.iwmde_ordered_block_log_likelihood <- function(context, block, total, share) {

  n <- max(length(total), length(share))
  total <- rep_len(total, n)
  share <- rep_len(share, n)
  block[["counts"]][["candidate_rows"]] <- block[["counts"]][["candidate_rows"]] + n
  if (!is.null(block[["source_state"]])) {
    return(as.numeric(.iwmde_ordered_log_q_source(context, block[["source_state"]], total)) -
      BayesTools::lpdf(block[["prior"]], total))
  }
  if (!block[["selected"]]) return(.iwmde_ordered_block_gaussian_log(block, total, share))
  grouping <- block[["selection_grouping"]]
  representatives <- if (is.null(grouping)) seq_along(block[["mu0"]]) else grouping[["representatives"]]
  powers <- if (is.null(grouping)) block[["weights"]] else grouping[["powers"]]
  selection_context_base <- if (is.null(grouping)) block[["selection_context"]] else block[["normalizer_context"]]
  out <- numeric(n)
  for (start in seq.int(1L, n, by = 1024L)) {
    take <- seq.int(start, min(n, start + 1023L))
    count <- length(take)
    mu <- matrix(block[["mu0"]][representatives], count, length(representatives), byrow = TRUE) +
      outer(share[take], block[["mu_c"]][representatives]) + outer(total[take], block[["mu_a"]][representatives]) +
      outer(total[take] * share[take], block[["mu_b"]][representatives])
    # The fixed positive singleton selection route has an invariant observation
    # factor and covariance. The declared bilinear predictor and original native
    # normalizer evaluate every source node without recompiling its setup.
    selection_context <- selection_context_base
    if (nrow(selection_context[["omega"]]) == 1L && count > 1L) {
      selection_context[["omega"]] <- matrix(rep(selection_context[["omega"]][1L, ], each = count), nrow = count)
    }
    log_norm <- .selection_step_log_norm_matrix(mu,
      matrix(rep(sqrt(block[["variance"]][representatives]), each = count), nrow = count),
      block[["selection_sei"]][representatives], selection_context)
    out[take] <- .iwmde_ordered_block_gaussian_log(block, total[take], share[take]) +
      block[["selection_observed_constant"]] - as.numeric(log_norm %*% powers)
    block[["counts"]][["native_nodes"]] <- block[["counts"]][["native_nodes"]] + count
    block[["counts"]][["maximum_chunk"]] <- max(block[["counts"]][["maximum_chunk"]], count)
  }
  out
}


.iwmde_ordered_normal_completion <- function(block, share) {

  normal <- block[["normal"]]
  parts <- .iwmde_ordered_gaussian_descriptors(block, share)
  descriptors <- parts[["descriptors"]][parts[["index"]]]
  A <- vapply(descriptors, `[[`, numeric(1), "A")
  t_ref <- vapply(descriptors, `[[`, numeric(1), "t_ref")
  H_ref <- vapply(descriptors, `[[`, numeric(1), "H_ref")
  precision <- A + 1 / normal[["sd"]]^2
  mean <- t_ref + (H_ref + (normal[["mean"]] - t_ref) / normal[["sd"]]^2) / precision
  sd <- 1 / sqrt(precision)
  scales <- numeric(length(share))
  log_L <- vapply(seq_along(share), function(i) {
    descriptor <- descriptors[[i]]
    if (!is.null(descriptor[["evaluate"]])) {
      evaluated <- descriptor[["evaluate"]](mean[[i]])
      scales[[i]] <<- evaluated[["max_log_scale"]]
      return(evaluated[["value"]])
    }
    if (!is.null(descriptor[["a"]])) {
      mu <- descriptor[["a"]] + descriptor[["b"]] * mean[[i]]
      pieces <- block[["weights"]] * stats::dnorm(block[["setup"]][["yi"]], mu, sqrt(block[["variance"]]), log = TRUE)
      scales[[i]] <<- max(sum(abs(pieces)), .iwmde_ordered_gaussian_residual_scale(block, mu,
        abs(block[["mu0"]]) + abs(share[[i]] * block[["mu_c"]]) + abs(descriptor[["b"]] * mean[[i]])))
      return(sum(pieces))
    }
    delta <- mean[[i]] - descriptor[["t_ref"]]
    scales[[i]] <<- max(abs(descriptor[["H_ref"]] * delta), abs(.5 * descriptor[["A"]] * delta^2))
    descriptor[["C_ref"]] + descriptor[["H_ref"]] * delta - .5 * descriptor[["A"]] * delta^2
  }, numeric(1))
  prior_square <- .5 * ((mean - normal[["mean"]]) / normal[["sd"]])^2
  log_precision <- .5 * log(normal[["sd"]]^2 * precision)
  residual <- log_L - prior_square - log_precision
  posterior_mass <- .iwmde_qcmde_log_normal_mass(normal[["lower"]], normal[["upper"]], mean, sd)
  original_mass <- .iwmde_qcmde_log_normal_mass(normal[["lower"]], normal[["upper"]], normal[["mean"]], normal[["sd"]])
  log_mass <- residual + posterior_mass - original_mass
  if (any(!is.finite(log_mass)) || any(!is.finite(mean)) || any(!is.finite(sd)) || any(sd <= 0)) {
    .iwmde_ordered_numeric_stop("gaussian_completion_unrepresentable", "the completed Gaussian mass is unrepresentable",
      list(log_mass = log_mass, mean = mean, sd = sd, residual = residual,
        posterior_interval_log_mass = posterior_mass, original_prior_interval_log_mass = original_mass))
  }
  list(log_mass = log_mass, mean = mean, sd = sd,
    residual_log_mass = residual, posterior_interval_log_mass = posterior_mass,
    original_prior_interval_log_mass = original_mass,
    max_log_scale = pmax(1, abs(log_L), abs(residual), abs(posterior_mass), abs(original_mass),
      scales, prior_square, abs(log_precision),
      (abs(t_ref) + (abs(H_ref) + (abs(normal[["mean"]]) + abs(t_ref)) / normal[["sd"]]^2) / precision) *
        abs(mean - normal[["mean"]]) / normal[["sd"]]^2))
}


.iwmde_ordered_profile_geometry <- function(block) {

  gaussian <- block[["gaussian"]]
  normal_mass <- !is.null(block[["normal"]]) && all(is.infinite(c(block[["normal"]][["lower"]], block[["normal"]][["upper"]])))
  r <- gaussian[["A"]]
  h <- gaussian[["H"]]
  c <- gaussian[["C"]]
  if (normal_mass) {
    normal <- block[["normal"]]
    r[[1L]] <- r[[1L]] + 1 / normal[["sd"]]^2
    h[[1L]] <- h[[1L]] + normal[["mean"]] / normal[["sd"]]^2
    c[[1L]] <- c[[1L]] - .5 * (normal[["mean"]] / normal[["sd"]])^2 - log(normal[["sd"]])
  }
  moving_constant <- all(block[["mu_a"]] == 0) && all(block[["mu_b"]] == 0)
  if (!normal_mass && moving_constant) {
    roots <- .iwmde_ordered_poly_roots(.iwmde_ordered_poly_derivative(c))
    points <- if (nrow(roots[["roots"]])) rowMeans(roots[["roots"]]) else numeric()
    profile <- function(s) .iwmde_ordered_block_gaussian_log(block, 0, s)
    first <- function(s) .iwmde_ordered_poly_value(.iwmde_ordered_poly_derivative(c), s)
    second <- function(s) .iwmde_ordered_poly_value(.iwmde_ordered_poly_derivative(.iwmde_ordered_poly_derivative(c)), s)
    polynomial <- .iwmde_ordered_poly_derivative(c)
    records <- list(stationary = roots)
    switches <- numeric()
  } else {
    candidates <- c(0, 1)
    if (length(r) == 3L && r[[3L]] > 0) {
      vertex <- -r[[2L]] / (2 * r[[3L]])
      if (vertex > 0 && vertex < 1) candidates <- c(candidates, vertex)
    }
    if (any(.iwmde_ordered_poly_value(r, candidates) <= 0)) {
      .iwmde_ordered_numeric_stop("proper_profile_geometry_unresolved", "the Gaussian profile precision is not strictly positive on source support",
        list(precision = r, candidates = candidates))
    }
    multiply <- .iwmde_ordered_poly_multiply
    derivative <- .iwmde_ordered_poly_derivative
    polynomial <- .iwmde_ordered_poly_add(2 * multiply(derivative(c), multiply(r, r)),
      if (normal_mass) -multiply(r, derivative(r)) else 0,
      2 * multiply(multiply(h, derivative(h)), r), -multiply(multiply(h, h), derivative(r)))
    stationary <- .iwmde_ordered_poly_roots(polynomial)
    points <- if (nrow(stationary[["roots"]])) rowMeans(stationary[["roots"]]) else numeric()
    records <- list(stationary = stationary)
    switches <- numeric()
    support <- block[["normal_probability"]][["support"]]
    if (!normal_mass) {
      for (endpoint in support[is.finite(support)]) {
        switch_record <- .iwmde_ordered_poly_roots(.iwmde_ordered_poly_add(h, -endpoint * r))
        records[[paste0("switch_", endpoint)]] <- switch_record
        if (nrow(switch_record[["roots"]])) switches <- c(switches, rowMeans(switch_record[["roots"]]))
        constrained <- .iwmde_ordered_poly_add(c, endpoint * h, -.5 * endpoint^2 * r)
        constrained_record <- .iwmde_ordered_poly_roots(derivative(constrained))
        records[[paste0("constrained_", endpoint)]] <- constrained_record
        if (nrow(constrained_record[["roots"]])) {
          positions <- rowMeans(constrained_record[["roots"]])
          mode <- .iwmde_ordered_poly_value(h, positions) / .iwmde_ordered_poly_value(r, positions)
          positions <- positions[if (endpoint == support[[1L]]) mode <= endpoint else mode >= endpoint]
          points <- c(points, positions)
        }
      }
    }
    profile <- function(s) {
      if (normal_mass) return(.iwmde_ordered_normal_completion(block, s)[["log_mass"]])
      parts <- .iwmde_ordered_gaussian_descriptors(block, s)
      descriptors <- parts[["descriptors"]][parts[["index"]]]
      mode <- vapply(descriptors, function(d) d[["t_ref"]] + d[["H_ref"]] / d[["A"]], numeric(1))
      source <- pmax(support[[1L]], pmin(support[[2L]], mode))
      .iwmde_ordered_block_gaussian_log(block, source, s)
    }
    first <- function(s) {
      rr <- .iwmde_ordered_poly_value(r, s)
      hh <- .iwmde_ordered_poly_value(h, s)
      if (normal_mass) return(.iwmde_ordered_poly_value(polynomial, s) / (2 * rr^2))
      source <- max(support[[1L]], min(support[[2L]], hh / rr))
      .iwmde_ordered_poly_value(derivative(c), s) + .iwmde_ordered_poly_value(derivative(h), s) * source -
        .5 * .iwmde_ordered_poly_value(derivative(r), s) * source^2
    }
    second <- function(s) {
      rr <- .iwmde_ordered_poly_value(r, s)
      r1 <- .iwmde_ordered_poly_value(derivative(r), s)
      if (normal_mass) return(.iwmde_ordered_poly_value(derivative(polynomial), s) / (2 * rr^2) -
        .iwmde_ordered_poly_value(polynomial, s) * r1 / rr^3)
      hh <- .iwmde_ordered_poly_value(h, s)
      source <- max(support[[1L]], min(support[[2L]], hh / rr))
      if (source != hh / rr) return(.iwmde_ordered_poly_value(derivative(derivative(c)), s) +
        source * .iwmde_ordered_poly_value(derivative(derivative(h)), s) -
        .5 * source^2 * .iwmde_ordered_poly_value(derivative(derivative(r)), s))
      h1 <- .iwmde_ordered_poly_value(derivative(h), s)
      .iwmde_ordered_poly_value(derivative(derivative(c)), s) +
        (h1^2 + hh * .iwmde_ordered_poly_value(derivative(derivative(h)), s)) / rr -
        2 * hh * h1 * r1 / rr^2 - .5 * hh^2 * .iwmde_ordered_poly_value(derivative(derivative(r)), s) / rr^2 +
        hh^2 * r1^2 / rr^3
    }
  }
  root_points <- sort(unique(c(points, switches)))
  centers <- widths <- numeric()
  types <- character()
  for (s in points[points > 0 & points < 1]) {
    curvature <- second(s)
    if (curvature < 0) {
      centers <- c(centers, s)
      widths <- c(widths, 1 / sqrt(-curvature))
      types <- c(types, "interior")
    }
  }
  for (s in c(0, 1)) {
    slope <- first(s)
    curvature <- second(s)
    if ((s == 0 && slope <= 0) || (s == 1 && slope >= 0)) {
      choices <- c(if (slope != 0) 1 / abs(slope), if (curvature < 0) 1 / sqrt(-curvature))
      if (!length(choices) && normal_mass && s == 0 && any(polynomial != 0)) {
        power <- which(polynomial != 0)[[1L]]
        choices <- (2 * r[[1L]]^2 * power / abs(polynomial[[power]]))^(1 / power)
      }
      if (!length(choices) && any(polynomial != 0)) {
        .iwmde_ordered_numeric_stop("profile_width_unresolved", "the endpoint profile width is uncertified", list(endpoint = s, polynomial = polynomial))
      }
      centers <- c(centers, s)
      widths <- c(widths, if (length(choices)) min(choices) else 1)
      types <- c(types, "endpoint")
    }
  }
  if (!length(centers) || any(!is.finite(widths)) || any(widths <= 0)) {
    .iwmde_ordered_numeric_stop("profile_width_unresolved", "no valid profile peak band was obtained", list(centers = centers, widths = widths))
  }
  extrema <- function(left, right) {
    relevant <- unlist(lapply(records, function(record) {
      brackets <- record[["roots"]]
      brackets <- brackets[brackets[, 2L] >= left & brackets[, 1L] <= right, , drop = FALSE]
      if (!nrow(brackets)) numeric() else c(brackets[, 1L], rowMeans(brackets), brackets[, 2L])
    }), use.names = FALSE)
    range(profile(c(left, right, relevant[relevant > left & relevant < right])))
  }
  list(normal_mass = normal_mass, polynomial = polynomial, root_records = records,
    roots = points, switches = switches, centers = centers, widths = widths, types = types,
    profile = profile, first = first, second = second, extrema = extrema)
}


.iwmde_ordered_bands <- function(centers, widths, radius) {

  bands <- cbind(pmax(0, centers - radius * widths), pmin(1, centers + radius * widths))
  bands <- bands[order(bands[, 1L]), , drop = FALSE]
  out <- bands[1L, , drop = FALSE]
  if (nrow(bands) > 1L) for (i in 2:nrow(bands)) {
    last <- nrow(out)
    if (bands[i, 1L] <= out[last, 2L]) out[last, 2L] <- max(out[last, 2L], bands[i, 2L]) else out <- rbind(out, bands[i, ])
  }
  out
}


.iwmde_ordered_band_complement <- function(bands) {

  intervals <- cbind(c(0, bands[, 2L]), c(bands[, 1L], 1))
  intervals[intervals[, 2L] > intervals[, 1L], , drop = FALSE]
}


.iwmde_ordered_inner_guard <- function(block, share, prior = block[["normal_probability"]]) {

  gaussian <- block[["gaussian"]]
  A <- .iwmde_ordered_poly_value(gaussian[["A"]], share)
  H <- .iwmde_ordered_poly_value(gaussian[["H"]], share)
  if (A == 0 && H == 0 && !block[["selected"]]) return(list(constant = TRUE))
  if (!is.finite(A) || A <= 0 || !is.finite(H)) {
    .iwmde_ordered_numeric_stop("inner_geometry_unresolved", "the proper-prior source precision is not positive", list(A = A, H = H))
  }
  support <- prior[["support"]]
  intersection <- function(value) max(support[[1L]], min(support[[2L]], value))
  unconstrained <- H / A
  center <- intersection(unconstrained)
  width <- 1 / sqrt(A)
  if (center != unconstrained) width <- min(width, 1 / abs(H - A * center))
  gaussian_log <- function(total) .iwmde_ordered_block_gaussian_log(block, total, share)
  fixed <- .iwmde_ordered_base_cuts(prior, c(center, center - width, center + width),
    recompute = function() gaussian_log(center))
  trace <- list()
  for (j in 0:30) {
    radius <- 2^j
    band <- vapply(center + c(-1, 1) * radius * width, intersection, numeric(1))
    cuts <- .iwmde_ordered_merge_cuts(fixed, band, "peak_band", prior)
    panels <- .iwmde_ordered_probability_panels(cuts, prior)
    lower <- .iwmde_log_sum_exp(vapply(panels, function(panel) {
      ends <- c(panel[["lower"]], panel[["upper"]])
      if (any(!is.finite(ends))) return(-Inf)
      panel[["log_mass"]] + min(gaussian_log(ends)) + block[["log_correction_lower"]]
    }, numeric(1)))
    if (!is.finite(lower)) .iwmde_ordered_numeric_stop("inner_lower_nonpositive", "finite probability panels have no positive lower mass")
    outside <- rbind(c(support[[1L]], band[[1L]]), c(band[[2L]], support[[2L]]))
    outside <- outside[outside[, 2L] > outside[, 1L], , drop = FALSE]
    upper <- if (!nrow(outside)) -Inf else .iwmde_log_sum_exp(vapply(seq_len(nrow(outside)), function(i) {
      panel <- prior[["panel"]](outside[i, 1L], outside[i, 2L])
      maximizer <- max(outside[i, 1L], min(outside[i, 2L], unconstrained))
      panel[["log_mass"]] + gaussian_log(maximizer) + block[["log_correction_upper"]]
    }, numeric(1)))
    fraction <- exp(upper - lower)
    trace[[length(trace) + 1L]] <- list(radius = radius, log_lower = lower, log_upper = upper, fraction = fraction)
    if (fraction <= .025 / 16 || !nrow(outside)) break
    if (j == 30L) .iwmde_ordered_numeric_stop("inner_guard_unresolved", "the source missing-mass guard exhausted thirty doublings", trace)
  }
  geometric <- unlist(lapply(0:j, function(i) center + c(-1, 1) * 2^i * width), use.names = FALSE)
  cuts <- .iwmde_ordered_base_cuts(prior, c(center, geometric), recompute = function() gaussian_log(center))
  list(constant = FALSE, cuts = cuts, log_lower = lower, log_upper = upper,
    fraction = fraction, center = center, width = width, radius = radius, trace = trace)
}


.iwmde_ordered_positive_integral <- function(prior, cuts, evaluate, log_exterior = -Inf,
                                             target = .025 / 8, label = "source", evaluate_error = evaluate) {

  refinements <- 0L
  repeat {
    if (2L * (length(cuts) - 1L) > 100L) {
      .iwmde_ordered_numeric_stop("panel_budget_exceeded", paste0("the ", label, " partition exceeds one hundred fine panels"),
        list(cuts = cuts, refinements = refinements))
    }
    fine <- .iwmde_ordered_probability_refine(cuts, prior)
    coarse <- .iwmde_ordered_panel_rule(cuts, prior, 8L, evaluate)
    primary <- .iwmde_ordered_panel_rule(cuts, prior, 16L, evaluate)
    final <- .iwmde_ordered_panel_rule(fine, prior, 16L, evaluate_error)
    gl_error <- .iwmde_ordered_log_difference(coarse[["log_value"]], primary[["log_value"]])
    partition_error <- .iwmde_ordered_log_difference(primary[["log_value"]], final[["log_value"]])
    inner_error <- final[["log_inner_error"]]
    inner_exterior <- final[["log_inner_exterior"]]
    errors <- c(outer_gl = gl_error, partition = partition_error, weighted_inner = inner_error,
      inner_exterior = inner_exterior, outer_exterior = log_exterior)
    fraction <- exp(.iwmde_log_sum_exp(errors) - final[["log_value"]])
    if (is.finite(fraction) && fraction <= target) break
    if (length(fine) - 1L >= 100L) {
      .iwmde_ordered_numeric_stop("panel_budget_exceeded", paste0("the ", label, " refinement exhausted one hundred fine panels"),
        list(log_value = final[["log_value"]], log_errors = errors, relative_error = fraction,
          fine_cuts = fine, refinements = refinements))
    }
    split <- vapply(seq_along(primary[["panels"]]), function(i) {
      .iwmde_log_sum_exp(final[["panels"]][c(2L * i - 1L, 2L * i)])
    }, numeric(1))
    panel_error <- pmax(.iwmde_ordered_log_difference(coarse[["panels"]], primary[["panels"]]),
      .iwmde_ordered_log_difference(primary[["panels"]], split))
    index <- which.max(panel_error)
    panel <- prior[["panel"]](cuts[[index]], cuts[[index + 1L]])
    cuts <- .iwmde_ordered_merge_cuts(cuts, prior[["conditional"]](panel, .5), "adaptive_refinement", prior)
    refinements <- refinements + 1L
  }
  list(log_value = final[["log_value"]], log_error = .iwmde_log_sum_exp(errors[1:3]),
    log_exterior = .iwmde_log_sum_exp(errors[4:5]), relative_error = fraction,
    components = exp(errors - final[["log_value"]]), fine_panels = length(fine) - 1L,
    cuts = cuts, fine_cuts = fine, refinements = refinements,
    coalesced_cuts = attr(cuts, "ordered_coalesced_cuts", exact = TRUE))
}


.iwmde_ordered_gamma_completion <- function(context, block, share, prior) {

  result <- list(eligible = FALSE, reason = "unresolved", method = "proper_cdf_source",
    effective_prior = prior[["prior"]], z = NA_real_, condition = NA_real_,
    max_log_scale = NA_real_, prospective_floor = NA_real_)
  fallback <- function(reason) {
    result[["reason"]] <- reason
    result
  }
  if (!is.null(block[["source_state"]])) return(fallback("source_state"))
  if (!identical(block[["selected"]], FALSE)) return(fallback("selected"))
  supplied <- prior[["prior"]]
  if (!identical(supplied[["distribution"]], "gamma")) return(fallback("source_family"))
  if (.iwmde_retained_location_dynamic_prior(supplied)) return(fallback("dynamic_prior"))
  if (!identical(prior[["support"]], c(0, Inf)) ||
      !identical(.iwmde_prior_support(supplied), c(0, Inf))) return(fallback("truncated_prior"))
  shape <- supplied[["parameters"]][["shape"]]
  rate <- supplied[["parameters"]][["rate"]]
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) ||
      !shape %in% c(1, 2)) return(fallback("shape_outside_scope"))
  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate) || rate <= 0) return(fallback("invalid_rate"))
  polynomial <- block[["gaussian"]]
  if (!all(c("A", "H", "C") %in% names(polynomial)) ||
      !all(vapply(polynomial[c("A", "H", "C")], function(value)
        is.numeric(value) && length(value) > 0L && all(is.finite(value)), logical(1)))) return(fallback("invalid_geometry"))
  A <- .iwmde_ordered_poly_value(polynomial[["A"]], share)
  H <- .iwmde_ordered_poly_value(polynomial[["H"]], share)
  C <- .iwmde_ordered_poly_value(polynomial[["C"]], share)
  if (length(A) != 1L || length(H) != 1L || length(C) != 1L ||
      any(!is.finite(c(A, H, C))) || A < 0) return(fallback("invalid_geometry"))
  if (A == 0) return(fallback("existing_zero_precision_guard"))
  beta <- H - rate
  z <- beta / sqrt(A)
  result[["z"]] <- z
  if (any(!is.finite(c(beta, z, log(A))))) return(fallback("required_log_unrepresentable"))
  # Reject before the cancellation-prone log Mills-ratio subtraction.
  if (z < 0 && !(z^2 < 63.5)) return(fallback("stable_negative_bound"))
  log_phi <- stats::dnorm(z, log = TRUE)
  log_Phi <- stats::pnorm(z, log.p = TRUE)
  if (any(!is.finite(c(log_phi, log_Phi)))) return(fallback("required_log_unrepresentable"))
  condition <- 1
  subtraction <- NA_real_
  if (z < 0) {
    u <- log(-z) + log_Phi - log_phi
    subtraction <- -expm1(u)
    if (!is.finite(u) || !is.finite(subtraction) || subtraction <= 0) return(fallback("nonpositive_subtraction"))
    condition <- (1 + exp(u)) / subtraction
  }
  result[["condition"]] <- condition
  if (!is.finite(condition) || condition > 128) return(fallback("condition_over_128"))
  log_I <- if (shape == 1) {
    .5 * log(2 * pi) - .5 * log(A) + .5 * z^2 + log_Phi
  } else if (z < 0) {
    -log(A) + log(subtraction)
  } else if (z == 0) {
    -log(A)
  } else {
    term <- log(z) + log_Phi - log_phi
    big <- max(0, term)
    -log(A) + big + log(exp(-big) + exp(term - big))
  }
  log_mass <- C + shape * log(rate) - lgamma(shape) + log_I
  result[["log_mass"]] <- log_mass
  if (any(!is.finite(c(log_I, log_mass)))) return(fallback("result_unrepresentable"))
  t_bound <- abs(beta) / A + sqrt(shape) / sqrt(A)
  scales <- c(constant = 1, C = abs(C), mean = (abs(H) + rate) * t_bound,
    condition = if (z < 0) condition * (1 + abs(log(-z)) + abs(log_Phi) + abs(log_phi)) else
      if (z == 0) 1 + abs(log_Phi) + abs(log_phi) else 1 + z^2 / 2 + abs(log_Phi),
    log_A = abs(log(A)), log_rate = shape * abs(log(rate)), log_gamma = abs(lgamma(shape)),
    log_I = abs(log_I), log_mass = abs(log_mass), quadratic = .5 * (sqrt(A) * t_bound)^2)
  n_observations <- nrow(context[["data"]][["outcome"]])
  if (is.null(n_observations)) n_observations <- 1L
  result[["t_bound"]] <- t_bound
  result[["scales"]] <- scales
  result[["max_log_scale"]] <- max(scales)
  result[["prospective_floor"]] <- 64 * .Machine$double.eps * max(1L, n_observations) * max(scales)
  if (any(!is.finite(c(t_bound, scales, result[["prospective_floor"]])))) return(fallback("scale_unrepresentable"))
  if (result[["prospective_floor"]] > log1p(.025 / 16)) return(fallback("prospective_floor"))
  result[["eligible"]] <- TRUE
  result[["reason"]] <- "completion"
  result[["method"]] <- "completed_gamma_source"
  result
}


.iwmde_ordered_gamma_record <- function(block, completion) {

  ledger <- block[["counts"]][["gamma_completion"]]
  if (is.null(ledger)) ledger <- list(method = "completed_gamma_source", accepted = 0L, fallback = 0L,
    fallback_reasons = integer(), effective_priors = list(), z_range = c(Inf, -Inf),
    condition_range = c(Inf, -Inf), max_log_scale = 0, max_prospective_floor = 0,
    scale_updates = 0L, max_scales = numeric())
  prior <- completion[["effective_prior"]]
  if (!any(vapply(ledger[["effective_priors"]], identical, logical(1), prior))) {
    ledger[["effective_priors"]][[length(ledger[["effective_priors"]]) + 1L]] <- prior
  }
  if (completion[["eligible"]]) {
    ledger[["accepted"]] <- ledger[["accepted"]] + 1L
    ledger[["z_range"]] <- range(ledger[["z_range"]], completion[["z"]], finite = TRUE)
    ledger[["condition_range"]] <- range(ledger[["condition_range"]], completion[["condition"]], finite = TRUE)
    ledger[["max_log_scale"]] <- max(ledger[["max_log_scale"]], completion[["max_log_scale"]])
    ledger[["max_prospective_floor"]] <- max(ledger[["max_prospective_floor"]], completion[["prospective_floor"]])
    ledger[["max_scales"]] <- if (!length(ledger[["max_scales"]])) completion[["scales"]] else
      pmax(ledger[["max_scales"]], completion[["scales"]])
    ledger[["scale_updates"]] <- ledger[["scale_updates"]] + 1L
  } else {
    ledger[["fallback"]] <- ledger[["fallback"]] + 1L
    reason <- completion[["reason"]]
    previous <- ledger[["fallback_reasons"]][reason]
    ledger[["fallback_reasons"]][reason] <- if (is.na(previous)) 1L else previous + 1L
  }
  block[["counts"]][["gamma_completion"]] <- ledger
}


.iwmde_ordered_block_counts <- function(block) {

  counts <- as.list(block[["counts"]])
  completion <- counts[["gamma_completion"]]
  if (!is.null(completion) && completion[["accepted"]] == 0L) {
    completion[["z_range"]] <- completion[["condition_range"]] <- c(NA_real_, NA_real_)
    counts[["gamma_completion"]] <- completion
  }
  counts
}


.iwmde_ordered_block_inner <- function(context, block, share, order = NULL,
                                      prior = block[["normal_probability"]]) {

  block[["counts"]][["inner_calls"]] <- block[["counts"]][["inner_calls"]] + 1L
  if (!is.null(block[["normal"]]) && identical(prior, block[["normal_probability"]]) && !block[["selected"]]) {
    completion <- .iwmde_ordered_normal_completion(block, share)
    block[["counts"]][["max_log_scale"]] <- max(block[["counts"]][["max_log_scale"]], completion[["max_log_scale"]])
    return(list(log_value = completion[["log_mass"]], log_error = -Inf, log_exterior = -Inf,
      fine_panels = 0L, method = "completed_gaussian_source", completion = completion))
  }
  if (!is.null(order)) {
    completion <- .iwmde_ordered_normal_completion(block, share)
    rule <- .gauss_hermite_nodes(order)
    total <- completion[["mean"]] + completion[["sd"]] * rule[["nodes"]]
    actual <- .iwmde_ordered_block_log_likelihood(context, block, total, share)
    correction <- actual - .iwmde_ordered_block_gaussian_log(block, total, share)
    value <- completion[["log_mass"]] + .iwmde_log_sum_exp(log(rule[["weights"]]) + correction)
    tail <- completion[["log_mass"]] + block[["log_correction_upper"]] + log(2) +
      stats::pnorm(max(abs(rule[["nodes"]])), lower.tail = FALSE, log.p = TRUE)
    block[["counts"]][["max_log_scale"]] <- max(block[["counts"]][["max_log_scale"]], completion[["max_log_scale"]], abs(actual), abs(correction))
    return(list(log_value = value, log_error = -Inf, log_exterior = tail,
      fine_panels = 0L, method = "completed_gaussian_gh", order = order))
  }
  gaussian_A <- block[["gaussian"]][["A"]]
  zero_precision <- is.numeric(gaussian_A) && length(gaussian_A) > 0L &&
    length(share) == 1L && is.finite(.iwmde_ordered_poly_value(gaussian_A, share)) &&
    .iwmde_ordered_poly_value(gaussian_A, share) == 0
  if (identical(prior[["prior"]][["distribution"]], "gamma") && !zero_precision) {
    completion <- .iwmde_ordered_gamma_completion(context, block, share, prior)
    if (!identical(completion[["reason"]], "existing_zero_precision_guard")) .iwmde_ordered_gamma_record(block, completion)
    if (completion[["eligible"]]) {
      block[["counts"]][["max_log_scale"]] <- max(block[["counts"]][["max_log_scale"]], completion[["max_log_scale"]])
      return(list(log_value = completion[["log_mass"]], log_error = -Inf, log_exterior = -Inf,
        fine_panels = 0L, method = "completed_gamma_source", completion = completion))
    }
  }
  guard <- .iwmde_ordered_inner_guard(block, share, prior)
  if (guard[["constant"]]) {
    return(list(log_value = .iwmde_ordered_block_gaussian_log(block, 0, share),
      log_error = -Inf, log_exterior = -Inf, fine_panels = 0L, method = "declared_constant_source"))
  }
  evaluate <- function(total) {
    value <- .iwmde_ordered_block_log_likelihood(context, block, total, share)
    block[["counts"]][["max_log_scale"]] <- max(block[["counts"]][["max_log_scale"]], abs(value))
    list(log_value = value)
  }
  out <- .iwmde_ordered_positive_integral(prior, guard[["cuts"]], evaluate, guard[["log_upper"]], label = "inner source")
  out[["guard"]] <- guard
  out[["method"]] <- "proper_cdf_source"
  out
}


.iwmde_ordered_block_gh_many <- function(context, block, shares, order) {

  completion <- .iwmde_ordered_normal_completion(block, shares)
  rule <- .gauss_hermite_nodes(order)
  totals <- rep(completion[["mean"]], each = order) +
    rep(completion[["sd"]], each = order) * rep(rule[["nodes"]], length(shares))
  source_shares <- rep(shares, each = order)
  actual <- .iwmde_ordered_block_log_likelihood(context, block, totals, source_shares)
  correction <- actual - .iwmde_ordered_block_gaussian_log(block, totals, source_shares)
  block[["counts"]][["inner_calls"]] <- block[["counts"]][["inner_calls"]] + length(shares)
  block[["counts"]][["max_log_scale"]] <- max(block[["counts"]][["max_log_scale"]],
    completion[["max_log_scale"]], abs(actual), abs(correction))
  log_tail <- block[["log_correction_upper"]] + log(2) +
    stats::pnorm(max(abs(rule[["nodes"]])), lower.tail = FALSE, log.p = TRUE)
  lapply(seq_along(shares), function(i) {
    take <- seq.int((i - 1L) * order + 1L, i * order)
    list(log_value = completion[["log_mass"]][[i]] +
      .iwmde_log_sum_exp(log(rule[["weights"]]) + correction[take]),
      log_error = -Inf, log_exterior = completion[["log_mass"]][[i]] + log_tail,
      fine_panels = 0L, method = "completed_gaussian_gh", order = order)
  })
}


.iwmde_ordered_block_outer_guard <- function(block, prior = block[["beta_probability"]]) {

  geometry <- .iwmde_ordered_profile_geometry(block)
  normal_mass <- geometry[["normal_mass"]]
  roots <- unique(unlist(lapply(geometry[["root_records"]], function(record) {
    if (nrow(record[["roots"]])) rowMeans(record[["roots"]]) else numeric()
  }), use.names = FALSE))
  anchors <- c(roots, geometry[["switches"]], geometry[["centers"]],
    geometry[["centers"]] - geometry[["widths"]], geometry[["centers"]] + geometry[["widths"]])
  fixed <- .iwmde_ordered_base_cuts(prior, anchors, recompute = function() geometry[["extrema"]](0, 1))
  rectangle <- block[["normal_probability"]][["quantile"]](c(.25, .75))
  trace <- list()
  for (j in 0:30) {
    radius <- 2^j
    bands <- .iwmde_ordered_bands(geometry[["centers"]], geometry[["widths"]], radius)
    outside <- .iwmde_ordered_band_complement(bands)
    cuts <- .iwmde_ordered_merge_cuts(fixed, as.vector(bands), "peak_band", prior)
    panels <- .iwmde_ordered_probability_panels(cuts, prior)
    lower <- .iwmde_log_sum_exp(vapply(panels, function(panel) {
      ends <- c(panel[["lower"]], panel[["upper"]])
      if (normal_mass) {
        minimum <- geometry[["extrema"]](ends[[1L]], ends[[2L]])[[1L]]
      } else {
        corners <- expand.grid(total = rectangle, share = ends)
        minimum <- log(.5) + min(.iwmde_ordered_block_gaussian_log(block, corners[["total"]], corners[["share"]]))
      }
      panel[["log_mass"]] + minimum + block[["log_correction_lower"]]
    }, numeric(1)))
    if (!is.finite(lower)) .iwmde_ordered_numeric_stop("outer_lower_nonpositive", "the positive profile or four-corner lower mass is unavailable")
    upper <- if (!nrow(outside)) -Inf else .iwmde_log_sum_exp(vapply(seq_len(nrow(outside)), function(i) {
      prior[["panel"]](outside[i, 1L], outside[i, 2L])[["log_mass"]] +
        geometry[["extrema"]](outside[i, 1L], outside[i, 2L])[[2L]] + block[["log_correction_upper"]]
    }, numeric(1)))
    fraction <- exp(upper - lower)
    trace[[length(trace) + 1L]] <- list(radius = radius, log_lower = lower, log_upper = upper, fraction = fraction)
    if (fraction <= .025 / 16 || !nrow(outside)) break
    if (j == 30L) .iwmde_ordered_numeric_stop("outer_guard_unresolved", "the share missing-mass guard exhausted thirty doublings", trace)
  }
  # Normal E uses only its final-band cuts. The proper-prior protocol keeps the
  # specified geometric sequence through the final radius.
  if (!normal_mass) {
    geometric <- unlist(lapply(0:j, function(i) c(geometry[["centers"]] - 2^i * geometry[["widths"]],
      geometry[["centers"]] + 2^i * geometry[["widths"]])), use.names = FALSE)
    cuts <- .iwmde_ordered_base_cuts(prior, c(roots, geometry[["switches"]], geometry[["centers"]], geometric),
      recompute = function() geometry[["extrema"]](0, 1))
  }
  list(cuts = cuts, log_lower = lower, log_upper = upper, fraction = fraction,
    bands = bands, outside = outside, radius = radius, geometry = geometry, trace = trace)
}


.iwmde_ordered_block_denominator <- function(context, block) {

  guard <- .iwmde_ordered_block_outer_guard(block)
  selected_gh <- block[["selected"]] && guard[["geometry"]][["normal_mass"]]
  orders <- if (selected_gh) c(15L, 31L, 63L, 127L, 255L, 511L, 1023L) else NA_integer_
  cache <- new.env(parent = emptyenv())
  evaluate <- function(share, order) {
    block[["counts"]][["outer_nodes"]] <- block[["counts"]][["outer_nodes"]] + length(share)
    if (!block[["selected"]] && !is.null(block[["normal"]])) {
      completion <- .iwmde_ordered_normal_completion(block, share)
      block[["counts"]][["inner_calls"]] <- block[["counts"]][["inner_calls"]] + length(share)
      block[["counts"]][["max_log_scale"]] <- max(block[["counts"]][["max_log_scale"]], completion[["max_log_scale"]])
      return(list(log_value = completion[["log_mass"]], log_error = rep(-Inf, length(share)),
        log_exterior = rep(-Inf, length(share))))
    }
    keys <- paste(vapply(share, .iwmde_key_number, character(1)), order, sep = ":")
    missing <- which(!vapply(keys, exists, logical(1), envir = cache, inherits = FALSE))
    if (length(missing)) {
      added <- if (!is.na(order)) .iwmde_ordered_block_gh_many(context, block, share[missing], order) else
        lapply(share[missing], function(s) .iwmde_ordered_block_inner(context, block, s))
      for (i in seq_along(missing)) cache[[keys[[missing[[i]]]]]] <- added[[i]]
    }
    values <- lapply(keys, get, envir = cache, inherits = FALSE)
    list(log_value = vapply(values, `[[`, numeric(1), "log_value"),
      log_error = vapply(values, `[[`, numeric(1), "log_error"),
      log_exterior = vapply(values, `[[`, numeric(1), "log_exterior"))
  }
  if (!selected_gh) {
    out <- .iwmde_ordered_positive_integral(block[["beta_probability"]], guard[["cuts"]],
      function(share) evaluate(share, NA_integer_), guard[["log_upper"]], target = .025 / 1.025,
      label = "outer share")
    out[["orders"]] <- integer()
    out[["guard"]] <- guard
    return(out)
  }
  # First eligible ladder member has two successive completed-mass changes.
  # Weighted absolute inner differences enter the same denominator envelope.
  cuts <- guard[["cuts"]]
  fine <- .iwmde_ordered_probability_refine(cuts, block[["beta_probability"]])
  attempts <- list()
  previous <- list()
  for (order in orders) {
    current <- .iwmde_ordered_panel_rule(fine, block[["beta_probability"]], 16L,
      function(share) evaluate(share, order))
    previous[[length(previous) + 1L]] <- current
    if (length(previous) < 3L) next
    use <- tail(orders[orders <= order], 3L)
    changes <- abs(expm1(diff(vapply(tail(previous, 3L), `[[`, numeric(1), "log_value"))))
    weighted <- vapply(seq_len(2L), function(i) {
      pair <- use[c(i, i + 1L)]
      .iwmde_ordered_panel_rule(fine, block[["beta_probability"]], 16L, function(share) {
        a <- evaluate(share, pair[[1L]])[["log_value"]]
        b <- evaluate(share, pair[[2L]])[["log_value"]]
        list(log_value = .iwmde_ordered_log_difference(a, b))
      })[["log_value"]]
    }, numeric(1))
    inner_error <- max(weighted)
    outside <- .iwmde_ordered_log_add(current[["log_inner_exterior"]], guard[["log_upper"]])
    preliminary <- exp(.iwmde_ordered_log_add(inner_error, outside) - current[["log_value"]])
    attempts[[length(attempts) + 1L]] <- list(order = order, changes = changes, preliminary_error = preliminary)
    if (any(changes > .025) || preliminary > .025 / 1.025) next
    try_integral <- tryCatch(.iwmde_ordered_positive_integral(block[["beta_probability"]], cuts,
      function(share) evaluate(share, order),
      guard[["log_upper"]], target = .025 / 1.025, label = "selected outer share", evaluate_error =
      function(share) {
        value <- evaluate(share, order)
        errors <- lapply(seq_len(2L), function(i) .iwmde_ordered_log_difference(
          evaluate(share, use[[i]])[["log_value"]], evaluate(share, use[[i + 1L]])[["log_value"]]))
        value[["log_error"]] <- pmax(errors[[1L]], errors[[2L]])
        value
      }),
      RoBMA_ordered_panel_budget_exceeded = identity)
    if (!inherits(try_integral, "error")) {
      try_integral[["orders"]] <- use
      try_integral[["inner_changes"]] <- changes
      try_integral[["attempts"]] <- attempts
      try_integral[["guard"]] <- guard
      return(try_integral)
    }
    attempts[[length(attempts)]][["outer_failure"]] <- try_integral[["proof"]]
  }
  .iwmde_ordered_numeric_stop("inner_ladder_exhausted", "the selected source ladder exhausted order 1023 without two passing changes and a complete .025 ratio envelope",
    list(orders = orders, attempts = attempts, guard = guard, counts = .iwmde_ordered_block_counts(block)))
}


.iwmde_ordered_block_numerator <- function(context, block, supplier) {

  descriptor <- block[["descriptor"]]
  k <- abs(descriptor[["k"]])
  shift <- if (is.null(block[["shift"]])) 0 else block[["shift"]]
  zero_total <- -shift / descriptor[["a"]]
  if (is.null(descriptor[["a"]])) zero_total <- -shift
  if (identical(supplier, "total")) {
    alpha <- block[["beta"]][["parameters"]][["alpha"]]
    beta <- block[["beta"]][["parameters"]][["beta"]]
    factor <- BayesTools::lpdf(block[["prior"]], zero_total) + log((alpha + beta - 1) / (alpha - 1)) - log(k)
    if (block[["closed_total_numerator"]]) {
      value <- .iwmde_ordered_block_log_likelihood(context, block, zero_total, .5)
      return(list(log_value = factor + value, log_error = -Inf, log_exterior = -Inf,
        relative_error = 0, fine_panels = 0L, method = "closed_total_supplier"))
    }
    shifted <- BayesTools::prior("beta", list(alpha - 1, beta))
    prior <- .iwmde_ordered_probability_prior(shifted)
    numerator_block <- block
    numerator_block[["normal"]] <- NULL
    numerator_block[["mu0"]] <- block[["mu0"]] + zero_total * block[["mu_a"]]
    numerator_block[["mu_c"]] <- block[["mu_c"]] + zero_total * block[["mu_b"]]
    numerator_block[["mu_a"]] <- numerator_block[["mu_b"]] <- 0
    numerator_block[["gaussian"]] <- list(A = 0, H = 0,
      C = .iwmde_ordered_poly_add(block[["gaussian"]][["C"]], zero_total * block[["gaussian"]][["H"]],
        -.5 * zero_total^2 * block[["gaussian"]][["A"]]))
    guard <- .iwmde_ordered_block_outer_guard(numerator_block, prior)
    out <- .iwmde_ordered_positive_integral(prior, guard[["cuts"]], function(share) {
      list(log_value = .iwmde_ordered_block_log_likelihood(context, block, zero_total, share))
    }, guard[["log_upper"]], label = "shifted-Beta numerator")
    out[["guard"]] <- guard
    out[["log_value"]] <- out[["log_value"]] + factor
    out[["log_error"]] <- out[["log_error"]] + factor
    out[["log_exterior"]] <- out[["log_exterior"]] + factor
    out[["method"]] <- "shifted_beta_total_supplier"
    return(out)
  }
  if (!identical(supplier, "share")) .iwmde_ordered_numeric_stop("supplier_unresolved", "the finite offset supplier is unknown")
  prior <- block[["prior"]]
  shape <- prior[["parameters"]][["shape"]]
  rate <- prior[["parameters"]][["rate"]]
  support <- .iwmde_prior_support(prior)
  if (identical(prior[["distribution"]], "gamma") && shape > 1 && shift == 0 && identical(support, c(0, Inf))) {
    inverse <- BayesTools::prior("gamma", list(shape - 1, rate))
    factor <- BayesTools::lpdf(block[["beta"]], 0) + log(rate / (shape - 1)) - log(k)
    out <- .iwmde_ordered_block_inner(context, block, 0, prior = .iwmde_ordered_probability_prior(inverse))
    out[["log_value"]] <- out[["log_value"]] + factor
    out[["log_error"]] <- out[["log_error"]] + factor
    out[["log_exterior"]] <- out[["log_exterior"]] + factor
    out[["relative_error"]] <- exp(.iwmde_ordered_log_add(out[["log_error"]], out[["log_exterior"]]) - out[["log_value"]])
    out[["method"]] <- "unshifted_gamma_inverse_supplier"
    return(out)
  }
  if (support[[1L]] + shift > 0) {
    probability <- block[["normal_probability"]]
    guard <- .iwmde_ordered_inner_guard(block, 0, probability)
    if (guard[["constant"]]) {
      cuts <- .iwmde_ordered_base_cuts(probability)
      exterior <- -Inf
    } else {
      cuts <- guard[["cuts"]]
      exterior <- guard[["log_upper"]] - log(support[[1L]] + shift)
    }
    out <- .iwmde_ordered_positive_integral(probability, cuts, function(total) {
      list(log_value = .iwmde_ordered_block_log_likelihood(context, block, total, 0) - log(abs(total + shift)))
    }, exterior, label = "true affine-prior inverse numerator")
    factor <- BayesTools::lpdf(block[["beta"]], 0) - log(k)
    for (field in c("log_value", "log_error", "log_exterior")) out[[field]] <- out[[field]] + factor
    out[["method"]] <- "true_affine_prior_inverse_supplier"
    return(out)
  }
  .iwmde_ordered_stop(block[["state"]][["parameter"]], "offset_block_capability_unavailable",
    "the affine inverse supplier has no positive declared gap or finite unshifted reciprocal plan")
}


.iwmde_ordered_block_log_terms <- function(context, parameter_spec, row_states, point_plan, control = NULL) {

  n <- length(row_states)
  log_terms <- matrix(NA_real_, 1L, n)
  proofs <- vector("list", n)
  n_observations <- nrow(context[["data"]][["outcome"]])
  if (is.null(n_observations)) n_observations <- 1L
  active_keys <- vapply(row_states, function(state) {
    key <- state[["active_key"]]
    if (is.null(key)) "" else key
  }, character(1))
  partition <- if (all(nzchar(active_keys)) && length(unique(active_keys)) == 1L) {
    .iwmde_ordered_selection_partition(context, parameter_spec, row_states[[1L]])
  } else NULL
  for (i in seq_len(n)) {
    row_plan <- point_plan
    if (length(point_plan[["component_statuses"]])) {
      spec <- parameter_spec[["ordered_specs"]][[point_plan[["descriptor"]][["term"]]]]
      component <- .iwmde_ordered_total_component(spec, row_states[[i]][["row"]])[["component"]]
      row_plan[names(point_plan[["component_statuses"]][[component]])] <- point_plan[["component_statuses"]][[component]]
      if (identical(row_plan[["status"]], "raw")) {
        produced <- .iwmde_ordered_log_terms(context, parameter_spec, row_states[i], point_plan[["value"]],
          if (is.null(control)) .iwmde_density_control_resolve("qCMDE", NULL, "ordinate") else control)
        log_terms[1L, i] <- produced[["log_terms"]][1L, 1L]
        proofs[[i]] <- produced[["proofs"]][[1L]]
        next
      }
      if (identical(row_plan[["status"]], "zero")) {
        log_terms[1L, i] <- -Inf
        proofs[[i]] <- list(internal_method = "ordered_offset_block_cmde", normalization_basis = "ordered_offset_block",
          error_status = "exact", numerator_relative_error = 0, denominator_relative_error = 0,
          ratio_relative_error = 0, exterior_relative_error = 0, floating_relative_error = 0,
          source_grid_points = 0L, source_grid = c(NA_real_, NA_real_), log_normalizer = 0,
          normalizer = list(method = "declared_zero_component_limit"), counts = list(native_nodes = 0),
          candidate_rows = row_states[[i]][["row_index"]])
        next
      }
    }
    block <- .iwmde_ordered_block_geometry(context, parameter_spec, row_states[[i]], observation_partition = partition)
    cache <- list(cache = .iwmde_qcmde_cache_env(context), key = .iwmde_hash("ordered_normalized_block_row", list(
      algorithm = .iwmde_algorithm_version(), source = .iwmde_source_fingerprint(context),
      weights = parameter_spec[["weights"]], recipes = parameter_spec[["ordered_specs"]],
      point_plan = point_plan, row = row_states[[i]][["row"]], active = row_states[[i]][["active_key"]],
      chart = row_states[[i]][["ordered_chart"]], internal_method = "ordered_offset_block_cmde",
      working_distribution = "normalized_conditional_source", normalization_basis = "ordered_offset_block",
      integration_control = context[["integration_control"]], selection_grouping = block[["selection_grouping"]])))
    if (is.null(cache[["cache"]])) cache <- NULL
    cached <- .iwmde_qcmde_cache_get(cache)
    if (!is.null(cached)) {
      log_terms[1L, i] <- cached[["log_term"]]
      proofs[[i]] <- cached[["proof"]]
      proofs[[i]][["normalizer_cache_reused"]] <- TRUE
      next
    }
    block[["counts"]][["max_log_scale"]] <- 1
    denominator <- tryCatch(.iwmde_ordered_block_denominator(context, block), error = function(error) {
      error[["proof"]][["candidate_row"]] <- row_states[[i]][["row_index"]]
      stop(error)
    })
    numerator <- .iwmde_ordered_block_numerator(context, block, row_plan[["supplier"]])
    floor <- 64 * .Machine$double.eps * max(1L, n_observations) * max(1,
      block[["counts"]][["max_log_scale"]], abs(numerator[["log_value"]]), abs(denominator[["log_value"]]))
    e_N <- numerator[["relative_error"]] + expm1(floor)
    e_Z <- denominator[["relative_error"]] + expm1(floor)
    ratio <- if (e_Z < 1) (1 + e_N) / (1 - e_Z) - 1 else Inf
    proof <- list(internal_method = "ordered_offset_block_cmde", working_distribution = "normalized_conditional_source",
      normalization_basis = "ordered_offset_block", error_status = "estimate",
      numerator_relative_error = e_N, denominator_relative_error = e_Z, ratio_relative_error = ratio,
      outer_relative_error = denominator[["components"]][["outer_gl"]],
      partition_relative_error = denominator[["components"]][["partition"]],
      inner_absolute_relative_error = denominator[["components"]][["weighted_inner"]],
      inner_exterior_relative_bound = denominator[["components"]][["inner_exterior"]],
      outer_exterior_relative_bound = denominator[["components"]][["outer_exterior"]],
      exterior_relative_error = sum(denominator[["components"]][c("inner_exterior", "outer_exterior")]),
      floating_relative_error = expm1(floor), root_status = "accepted_nominal_polynomial",
      orders = denominator[["orders"]], panels = denominator[["fine_panels"]],
      candidate_rows = row_states[[i]][["row_index"]], log_normalizer = denominator[["log_value"]],
      source_support = .iwmde_prior_support(block[["prior"]]), source_grid = c(NA_real_, NA_real_),
      source_grid_points = 0L, normalizer = list(method = "positive_ordered_offset_block"),
      supplier = row_plan[["supplier"]], value = point_plan[["value"]], value_key = point_plan[["value_key"]],
      counts = .iwmde_ordered_block_counts(block), selection_grouping = block[["selection_grouping"]],
      source_identity = .iwmde_source_fingerprint(context), prior_identity = .iwmde_hash("ordered_block_priors", list(block[["prior"]], block[["beta"]])),
      control_identity = .iwmde_hash("ordered_block_controls", list(control, context[["integration_control"]], n_observations)),
      numerator = numerator, denominator = denominator)
    if (!is.finite(ratio) || ratio > .025) {
      .iwmde_ordered_numeric_stop("block_integration_error", "the complete block integration error exceeds .025", proof)
    }
    log_terms[1L, i] <- numerator[["log_value"]] - denominator[["log_value"]]
    proofs[[i]] <- proof
    .iwmde_qcmde_cache_set(cache, list(log_term = log_terms[1L, i], proof = proof))
  }
  list(log_terms = log_terms, proofs = proofs, internal_method = "ordered_offset_block_cmde",
    working_distribution = "normalized_conditional_source", normalization_basis = "ordered_offset_block")
}
