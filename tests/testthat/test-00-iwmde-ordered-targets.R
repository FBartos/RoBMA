context("Ordered primitive conditional density targets")

.ordered_target_fixture <- function(total = BayesTools::prior("normal", list(0, 1)),
                                     allocation = NULL, shared = FALSE, total_shared = NULL) {

  data <- data.frame(f = base::ordered(rep(c("lo", "mid", "hi", "top"), 10L),
    levels = c("lo", "mid", "hi", "top")), x = seq(-1, 1, length.out = 40L))
  normal <- BayesTools::prior("normal", list(0, 1))
  ordered <- BayesTools::prior_ordered(total, allocation = allocation,
    id = if (shared) "shared" else NULL)
  priors <- list(intercept = BayesTools::prior("point", list(0)), f = ordered)
  if (shared) priors <- c(priors, list(x = normal, `f:x` = if (is.null(total_shared)) ordered else
    BayesTools::prior_ordered(total_shared, allocation = allocation, id = "shared")))
  formula <- if (shared) ~ f*x else ~ f
  compiled <- BayesTools::JAGS_formula(formula, "mu", data, priors)
  draws <- matrix(0, 40L, 1L, dimnames = list(NULL, "mu_intercept"))
  if (shared) draws <- cbind(draws, mu_x = stats::qnorm(seq(.05, .95, length.out = 40L)))
  bound <- compiled[["prior_list"]]
  ordered_priors <- bound[vapply(bound, BayesTools::is.prior.ordered, logical(1))]
  source_specs <- Map(BayesTools:::.bt_ordered_spec, names(ordered_priors), ordered_priors)
  for (spec in source_specs) {
    source_prior <- spec[["total_prior"]]
    if (!is.null(spec[["total_node"]])) {
      node <- spec[["total_node"]]
      spike <- node[["spec"]][["kind"]] %in% c("spike_and_slab", "ordered_spike_and_slab")
      indicator <- rep(if (spike) c(0, 1) else seq_along(source_prior), length.out = 40L)
      draws <- cbind(draws, stats::setNames(data.frame(indicator), node[["spec"]][["indicator"]]))
      variable <- if (spike) source_prior[!vapply(source_prior, BayesTools::is.prior.point, logical(1))] else source_prior
      for (i in seq_along(node[["spec"]][["components"]])) {
        component <- node[["spec"]][["components"]][[i]]
        if (!is.null(component[["location"]])) next
        values <- BayesTools::quant(variable[[i]], seq(.05, .95, length.out = 40L))
        draws <- cbind(draws, stats::setNames(data.frame(values), component[["coordinates"]]))
      }
      source <- as.numeric(BayesTools:::.bt_deterministic_node_evaluate(node, BayesTools:::.bt_deterministic_lookup(as.matrix(draws))))
    } else {
    source <- if (BayesTools::is.prior.point(source_prior)) rep(source_prior[["parameters"]][["location"]], 40L) else
      BayesTools::quant(source_prior, seq(.05, .95, length.out = 40L))
    }
    draws <- cbind(draws, stats::setNames(data.frame(source), spec[["total_names"]]))
    for (record in spec[["allocations"]]) {
      gamma_names <- record[["gamma_coordinates"]]
      if (!length(gamma_names) || all(gamma_names %in% colnames(draws))) next
      eta <- cbind((seq_len(40L) + 10) / 10, (rev(seq_len(40L)) + 5) / 10,
        (seq_len(40L) %% 7L + 2) / 10)
      colnames(eta) <- gamma_names
      draws <- cbind(draws, eta)
    }
  }
  draws <- as.matrix(draws)
  for (spec in source_specs) {
    coefficients <- BayesTools:::.bt_deterministic_node_evaluate(
      BayesTools:::.bt_dnode_ordered_coefficients(spec), BayesTools:::.bt_deterministic_lookup(draws))
    draws <- cbind(draws, coefficients)
  }
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(draws)), sample = 40L),
    class = c("runjags", "BayesTools_fit", "list"))
  attr(fit, "prior_list") <- compiled[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = compiled[["formula_design"]])
  attr(fit, "parameter_map") <- BayesTools:::.bt_build_parameter_map(colnames(draws),
    prior_list = compiled[["prior_list"]], formula_design = list(mu = compiled[["formula_design"]]))
  fit <- BayesTools:::.bt_attach_fit_contract(BayesTools:::.bt_attach_draw_geometry(fit))
  context <- .iwmde_context_ensure_caches(list(object = list(fit = fit),
    posterior_samples = draws, flat_prior_list = compiled[["prior_list"]], chain_id = rep(1L, 40L)))
  list(context = context, compiled = compiled, fit = fit, data = data)
}


test_that("mixture and spike sources mutate real mapped primitives and replay stale totals", {

  priors <- list(BayesTools::prior_mixture(list(BayesTools::prior("normal", list(.3, 1)),
    BayesTools::prior("normal", list(-.5, 2)))),
    BayesTools::prior_spike_and_slab(BayesTools::prior("normal", list(0, 1))))
  for (total in priors) {
    fixture <- .ordered_target_fixture(total = total)
    context <- fixture[["context"]]
    context[["data"]] <- list(outcome = data.frame(yi = .3, sei = .4))
    spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
    target <- .iwmde_ordered_target_spec(context, "cell", list(type = "linear",
      weights = stats::setNames(1, spec[["coefficient_names"]][[1L]])))
    local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
      .iwmde_active_setup = function(...) list(),
      .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
      .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(context, posterior_samples, ...) {
        stats::dnorm(.3, posterior_samples[, spec[["coefficient_names"]][[1L]]], .4, log = TRUE)
      }, .is_data_multilevel = function(...) FALSE, .is_data_joint_selection = function(...) FALSE,
      .data_outcome_type = function(...) "norm", .package = "RoBMA")
    values <- .iwmde_parameter_values(context, "cell", target)
    if (BayesTools::is.prior.spike_and_slab(total)) expect_identical(values[seq.int(1L, 40L, by = 2L)], rep(0, 20L))
    for (row in if (BayesTools::is.prior.spike_and_slab(total)) c(2L, 4L) else c(1L, 2L)) {
      state <- .iwmde_ordered_row_state(context, target, row)
      primitive <- state[["ordered_chart"]][["source"]]
      expect_true(primitive %in% spec[["total_node"]][["dependencies"]])
      expect_false(primitive %in% spec[["total_names"]])
      state[["row"]][spec[["total_names"]]] <- state[["row"]][spec[["total_names"]]] + 1
      replacement <- .iwmde_ordered_replace_source(context, state, c(-.4, .8))[["samples"]]
      expect_identical(as.numeric(replacement[, primitive]), c(-.4, .8))
      expect_identical(as.numeric(replacement[, spec[["total_names"]]]), c(-.4, .8))
      inactive <- setdiff(spec[["total_node"]][["dependencies"]], c(primitive, spec[["total_node"]][["spec"]][["indicator"]]))
      if (length(inactive)) expect_identical(replacement[, inactive, drop = FALSE],
        matrix(rep(state[["row"]][inactive], each = 2L), nrow = 2L, dimnames = list(NULL, inactive)))
    }
  }
})


test_that("shared affine supplier classification uses population support and exact declared offset", {

  fixture <- .ordered_target_fixture(shared = TRUE, allocation = BayesTools::prior("dirichlet", list(c(1.5, 1.5, 1))))
  context <- fixture[["context"]]
  specs <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])
  weights <- stats::setNames(c(-2.5, -1.7), c(specs[["mu_f"]][["coefficient_names"]][[1L]],
    specs[["mu_f__xXx__x"]][["coefficient_names"]][[1L]]))
  target <- .iwmde_ordered_target_spec(context, "affine", list(type = "linear", weights = weights))
  descriptor <- target[["ordered_point_plan"]][["descriptor"]]
  expect_true(descriptor[["shifted"]])
  expect_identical(descriptor[["shift_support"]], c(-Inf, Inf))
  expect_identical(descriptor[["offset"]], 0)
  expect_identical(target[["ordered_point_plan"]][["status"]], "block")
  row <- context[["evaluator_cache"]][[target[["ordered_source_key"]]]][1L, ]
  shift <- .iwmde_ordered_affine_shift(target, row)
  expect_equal(shift, (-1.7 / -2.5) * row[[specs[["mu_f__xXx__x"]][["total_names"]]]], tolerance = 0)
  chart <- .iwmde_ordered_row_chart(context, target, 1L)
  expect_equal(unname(chart[["offset"]]), unname(chart[["slope"]] * shift), tolerance = 0)
  expect_identical(.iwmde_key_number(descriptor[["offset"]]), "0000000000000000")
  gamma <- .ordered_target_fixture(total = BayesTools::prior("gamma", list(2, 1)), shared = TRUE)
  gamma_target <- .iwmde_ordered_target_spec(gamma[["context"]], "affine", list(type = "linear", weights = -weights))
  expect_identical(gamma_target[["ordered_point_plan"]][["descriptor"]][["shift_support"]], c(0, Inf))
  expect_identical(gamma_target[["ordered_point_plan"]][["supplier"]], "share")
  expect_identical(.iwmde_ordered_row_chart(gamma[["context"]], gamma_target, 1L)[["kind"]], "group")
})


test_that("a uniform Gamma multiplier gap admits shape one and preserves source side", {

  fixture <- .ordered_target_fixture(total = BayesTools::prior("gamma", list(1, 1)), shared = TRUE,
    total_shared = BayesTools::prior("point", list(.7)))
  context <- fixture[["context"]]
  specs <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])
  names <- c(specs[["mu_f"]][["coefficient_names"]][[1L]], specs[["mu_f__xXx__x"]][["coefficient_names"]][[1L]])
  target <- .iwmde_ordered_target_spec(context, "gap", list(type = "linear", weights = stats::setNames(rep(-2.5, 2L), names)))
  expect_identical(target[["ordered_point_plan"]][["status"]], "raw")
  expect_identical(target[["ordered_point_plan"]][["supplier"]], "share")
  expect_identical(target[["ordered_point_plan"]][["exponent"]], Inf)
  expect_identical(target[["ordered_point_plan"]][["descriptor"]][["shift_support"]], c(.7, .7))
  chart <- .iwmde_ordered_row_chart(context, target, 1L)
  row <- context[["evaluator_cache"]][[target[["ordered_source_key"]]]][1L, ]
  expect_identical(chart[["kind"]], "group")
  expect_identical(chart[["slope"]], -2.5 * (row[[specs[["mu_f"]][["total_names"]]]] + .7))
  expect_lt(chart[["slope"]], 0)
  expect_equal(.iwmde_ordered_support(context, target), c(-Inf, 0), tolerance = 0)
  plain <- .ordered_target_fixture(total = BayesTools::prior("gamma", list(1, 1)))
  plain_spec <- BayesTools::JAGS_ordered_parameter_spec(plain[["fit"]])[[1L]]
  nongap <- .iwmde_ordered_target_spec(plain[["context"]], "nongap", list(type = "linear",
    weights = stats::setNames(1, plain_spec[["coefficient_names"]][[1L]])))
  expect_identical(nongap[["ordered_point_plan"]][["status"]], "infinite")
})


test_that("a zero focal point with a random shared affine total is not an atom", {

  spike <- BayesTools::prior_spike_and_slab(BayesTools::prior("normal", list(0, 1)))
  fixture <- .ordered_target_fixture(total = spike, shared = TRUE, total_shared = BayesTools::prior("normal", list(0, 1)),
    allocation = BayesTools::prior("dirichlet", list(c(1.5, 1.5, 1))))
  context <- fixture[["context"]]
  specs <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])
  weights <- stats::setNames(c(1, 2), c(specs[["mu_f"]][["coefficient_names"]][[1L]], specs[["mu_f__xXx__x"]][["coefficient_names"]][[1L]]))
  target <- .iwmde_ordered_target_spec(context, "affine", list(type = "linear", weights = weights))
  expect_identical(target[["ordered_point_plan"]][["status"]], "capability")
  expect_identical(target[["ordered_point_plan"]][["cause"]], "mixed_supplier")
  expect_true(all(.iwmde_ordered_components(target)[["active"]]))
  chart <- .iwmde_ordered_row_chart(context, target, 1L)
  expect_identical(chart[["term"]], "mu_f__xXx__x")
  expect_identical(chart[["offset"]], 0)
  expect_true(is.finite(chart[["slope"]]) && chart[["slope"]] != 0)
})


test_that("grouped charts keep literal contractions and tiny nonzero static weights", {

  fixture <- .ordered_target_fixture(total = BayesTools::prior("point", list(2)))
  context <- fixture[["context"]]
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  tiny <- 1e-250
  target <- .iwmde_ordered_target_spec(context, "tiny", list(type = "linear", weights = stats::setNames(c(tiny, tiny), spec[["coefficient_names"]][1:2])))
  chart <- .iwmde_ordered_row_chart(context, target, 1L)
  expect_identical(chart[["slope"]], 2 * tiny)
  expect_identical(chart[["offset"]], 0)
  expect_identical(unname(target[["weights"]]), c(tiny, tiny))
  expect_identical(target[["ordered_point_plan"]][["descriptor"]][["J"]], 1:2)
})


test_that("density completeness accepts only metadata-declared zero reference cells", {

  samples <- coda::mcmc(matrix(0, 20L, 1L, dimnames = list(NULL, "cell")))
  expect_false(.plot_brma_declared_zero_cell(samples, c(a = 0)))
  BayesTools::posterior_metadata(samples, "atoms") <- BayesTools::posterior_atom_attribute(
    data.frame(x = 0, mass = 1), source = "zero_factor_design")
  expect_true(.plot_brma_declared_zero_cell(samples, c(a = 0, b = 0)))
  expect_false(.plot_brma_declared_zero_cell(samples, NULL))
  expect_false(.plot_brma_declared_zero_cell(samples, c(a = 0, b = 1e-250)))
  expect_false(.plot_brma_declared_zero_cell(samples, c(a = NA_real_)))
})


test_that("exact normalized block rows replace the point before one common aggregation", {

  fixture <- .ordered_target_fixture()
  context <- fixture[["context"]]
  context[["data"]] <- list(outcome = data.frame(yi = .3, sei = .4))
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  parameter_spec <- list(type = "linear", weights = stats::setNames(c(1, 1), spec[["coefficient_names"]][1:2]))
  aggregate <- .iwmde_density_aggregate
  aggregation_calls <- 0L
  source_values <- NULL
  make_proofs <- function(states, method) lapply(states, function(state) list(
    ratio_relative_error = 0, numerator_relative_error = 0, denominator_relative_error = 0,
    exterior_relative_error = 0, floating_relative_error = 0, error_status = "exact",
    source_grid_points = 0L, source_grid = c(-Inf, Inf), log_normalizer = 0,
    normalizer = list(method = method), candidate_rows = state[["row_index"]]))
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(...) 0,
    .is_data_multilevel = function(...) FALSE, .is_data_joint_selection = function(...) FALSE,
    .data_outcome_type = function(...) "norm",
    .iwmde_ordered_point_certificate = function(...) list(available = TRUE),
    .iwmde_ordered_log_terms = function(context, parameter_spec, row_states, evaluation_values, control) {
      source_values <<- evaluation_values
      list(log_terms = matrix(-1, length(evaluation_values), length(row_states)),
        proofs = make_proofs(row_states, "source"), internal_method = "ordered_source_conditional_cmde",
        working_distribution = "normalized_conditional_source", normalization_basis = "source_chart")
    }, .iwmde_ordered_block_log_terms = function(context, parameter_spec, row_states, point_plan, control) {
      list(log_terms = matrix(log(seq_along(row_states) / 10), 1L),
        proofs = make_proofs(row_states, "block"), internal_method = "ordered_offset_block_cmde",
        working_distribution = "normalized_conditional_source", normalization_basis = "ordered_offset_block")
    }, .iwmde_density_aggregate = function(...) {
      aggregation_calls <<- aggregation_calls + 1L
      aggregate(...)
    }, .package = "RoBMA")
  plan <- .iwmde_plan(context, "cell", "qCMDE", list(samples = 20L),
    outputs = "ordinate", values = 0, parameter_spec = parameter_spec)
  execution <- .iwmde_plan_row_execution(context, plan)
  expect_null(plan[["grids"]][["normalization_grid"]])
  density <- .iwmde_ordered_density(context, plan, execution, c(-.1, 0, .1), TRUE)
  expect_identical(source_values, c(-.1, .1))
  expect_identical(aggregation_calls, 1L)
  expect_equal(density[["y"]][[2L]], mean(seq_len(20L) / 10), tolerance = 1e-14)
  expect_identical(density[["y"]][c(1L, 3L)], rep(exp(-1), 2L))
  expect_identical(attr(density[["scaled_contributions"]], "contribution_rows", exact = TRUE), execution[["active_rows"]])
  expect_identical(density[["internal_method"]], "ordered_source_with_offset_block_cmde")
  expect_length(density[["ordered_source_subplan"]], 20L)
  expect_length(density[["ordered_block_subplan"]], 20L)
})


test_that("ordered factory uses exact projection measures and avoids ordinary provenance reads", {

  fixture <- .ordered_target_fixture(allocation = c(.1, .2, .7))
  context <- fixture[["context"]]
  specs <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])
  weights <- stats::setNames(rep(1, 3L), specs[[1L]][["coefficient_names"]])
  target <- .iwmde_ordered_target_spec(context, "last", list(type = "linear", weights = weights))
  expect_identical(target[["type"]], "ordered")
  expect_identical(target[["ordered_projection"]][["values"]],
    as.numeric(context[["posterior_samples"]][, specs[[1L]][["total_names"]]]))
  direct <- .iwmde_ordered_target_spec(context, specs[[1L]][["total_names"]], NULL)
  expect_identical(direct[["ordered_projection"]][["values"]], target[["ordered_projection"]][["values"]])
  testthat::local_mocked_bindings(JAGS_ordered_parameter_spec = function(...) stop("unrelated provenance read"),
    .package = "BayesTools")
  local_mocked_bindings(.iwmde_deterministic_nodes = function(...) stop("unrelated node-table read"),
    .package = "RoBMA")
  expect_null(.iwmde_ordered_target_spec(context, "mu_intercept", NULL))
})


test_that("unknown ordered dispatch cannot reuse primitive or retained-location fallbacks", {

  unknown <- list(type = "unknown_ordered_stub")
  calls <- list(function() .iwmde_parameter_spec(list(), "x", unknown),
    function() .iwmde_parameter_values(list(), "x", unknown),
    function() .iwmde_parameter_components(list(), "x", unknown),
    function() .iwmde_parameter_row_supports(list(), "x", 1L, unknown),
    function() .iwmde_replacement_spec(list(), "x", unknown),
    function() .iwmde_target_key("x", unknown),
    function() .iwmde_likelihood_mode("x", unknown),
    function() .iwmde_predictor_update_basis(list(), "x", list(), unknown, list()),
    function() .iwmde_log_q_grid(list(), "x", 0, list(), unknown),
    function() .iwmde_qcmde_row_laws(list(), list(), unknown, NULL),
    function() .iwmde_chen_log_weight(list(), "x", unknown, integer(), numeric(), integer(), numeric(), c(-Inf, Inf)),
    function() .iwmde_estimate(list(), "x", "qCMDE", NULL, parameter_spec = unknown))
  for (call in calls) expect_error(call(), class = "RoBMA_density_method_unknown_target")
})


test_that("total source maps and replay replace stale coefficients in every shared consumer", {

  fixture <- .ordered_target_fixture(shared = TRUE)
  context <- fixture[["context"]]
  specs <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])
  spec <- specs[["mu_f"]]
  target <- .iwmde_ordered_target_spec(context, "first",
    list(type = "linear", weights = stats::setNames(-2.5, spec[["coefficient_names"]][[1L]])))
  chart <- .iwmde_ordered_row_chart(context, target, 1L)
  expect_identical(chart[["kind"]], "total")
  state <- list(row = context[["evaluator_cache"]][[target[["ordered_source_key"]]]][1L, ],
    ordered_chart = chart, ordered_specs = target[["ordered_specs"]])
  state[["row"]][unlist(lapply(specs, `[[`, "coefficient_names"), use.names = FALSE)] <- 999
  candidates <- c(-.4, .8)
  replaced <- .iwmde_ordered_replace_source(context, state, candidates)[["samples"]]
  gamma <- spec[["allocations"]][[1L]][["gamma_coordinates"]]
  share <- as.numeric(state[["row"]][gamma]) / sum(state[["row"]][gamma])
  expect_equal(unname(replaced[, spec[["coefficient_names"]], drop = FALSE]),
    outer(candidates, share), tolerance = 1e-15)
  for (consumer in specs) {
    total <- if (identical(consumer[["parameter"]], "mu_f")) candidates else
      rep(state[["row"]][[consumer[["total_names"]]]], length(candidates))
    expect_equal(unname(replaced[, consumer[["coefficient_names"]], drop = FALSE]),
      outer(total, share), tolerance = 1e-15)
  }
  projected <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]], weights = target[["weights"]],
    draws = replaced)
  expect_equal(projected[["values"]], chart[["offset"]] + chart[["slope"]] * candidates,
    tolerance = 1e-15)
})


test_that("group and pair sources preserve Gamma invariants and exact nonunit maps", {

  fixture <- .ordered_target_fixture(total = BayesTools::prior("point", list(2)))
  context <- fixture[["context"]]
  specs <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])
  spec <- specs[[1L]]
  for (weights in list(c(-2.5, -2.5, 0), c(-2.5, 1, 3))) {
    names(weights) <- spec[["coefficient_names"]]
    target <- .iwmde_ordered_target_spec(context, "cell", list(type = "linear", weights = weights))
    chart <- .iwmde_ordered_row_chart(context, target, 1L)
    expect_identical(chart[["kind"]], if (length(unique(weights)) == 2L) "group" else "pair")
    row <- context[["evaluator_cache"]][[target[["ordered_source_key"]]]][1L, ]
    state <- list(row = row, ordered_chart = chart, ordered_specs = specs)
    values <- c(0, .2, .7, 1)
    samples <- .iwmde_ordered_replace_source(context, state, values)[["samples"]]
    eta <- samples[, chart[["gamma_coordinates"]], drop = FALSE]
    expect_equal(rowSums(eta), rep(sum(row[chart[["gamma_coordinates"]]]), length(values)), tolerance = 1e-15)
    projected <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]], weights = weights, draws = samples)
    expect_equal(projected[["values"]], chart[["offset"]] + chart[["slope"]] * values, tolerance = 2e-15)
    if (identical(chart[["kind"]], "pair")) {
      other <- setdiff(seq_len(ncol(eta)), chart[["indices"]])
      expect_identical(as.numeric(eta[, other]), unname(rep(row[chart[["gamma_coordinates"]][other]], length(values))))
    }
  }
})


test_that("normalized source producer agrees with completed Normal conditionals before aggregation", {

  fixture <- .ordered_target_fixture()
  context <- fixture[["context"]]
  context[["data"]] <- list(outcome = data.frame(yi = .3, sei = .4))
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  target <- .iwmde_ordered_target_spec(context, "cell", list(type = "linear",
    weights = stats::setNames(-2.5, spec[["coefficient_names"]][[1L]])))
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(context, posterior_samples, ...) {
      stats::dnorm(.3, posterior_samples[, spec[["coefficient_names"]][[1L]]], .4, log = TRUE)
    }, .is_data_multilevel = function(...) FALSE,
    .is_data_joint_selection = function(...) FALSE,
    .data_outcome_type = function(...) "norm",
    .iwmde_density_aggregate = function(...) stop("producer aggregated independently"), .package = "RoBMA")
  states <- lapply(c(2L, 7L, 21L, 39L), function(row) .iwmde_ordered_row_state(context, target, row))
  values <- c(-.1, .03, .3)
  control <- .iwmde_density_control_resolve("qCMDE", list(normalization_points = 100L), "ordinate")
  result <- .iwmde_ordered_log_terms(context, target, states, values, control)
  expected <- vapply(states, function(state) {
    slope <- state[["ordered_chart"]][["slope"]]
    a <- slope / -2.5
    precision <- 1 + a^2 / .4^2
    mean <- (.3 * a / .4^2) / precision
    stats::dnorm(values / slope, mean, 1 / sqrt(precision), log = TRUE) - log(abs(slope))
  }, numeric(length(values)))
  error <- abs(expm1(result[["log_terms"]] - expected))
  reported <- vapply(result[["proofs"]], `[[`, numeric(1), "ratio_relative_error")
  expect_true(all(error <= matrix(rep(reported, each = length(values)), nrow = length(values))))
  expect_identical(result[["normalization_basis"]], "source_chart")
  expect_identical(result[["working_distribution"]], "normalized_conditional_source")
})


test_that("source kernels count original totals and shared Gamma priors once", {

  fixture <- .ordered_target_fixture(shared = TRUE)
  context <- fixture[["context"]]
  specs <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])
  spec <- specs[["mu_f"]]
  target <- .iwmde_ordered_target_spec(context, "cell", list(type = "linear",
    weights = stats::setNames(1, spec[["coefficient_names"]][[1L]])))
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(context, posterior_samples, ...) {
      # The original likelihood has both the main effect and its shared-shape
      # interaction. Its candidate coefficients must both have been replayed.
      means <- posterior_samples[, spec[["coefficient_names"]][[1L]]] +
        1.7 * posterior_samples[, specs[["mu_f__xXx__x"]][["coefficient_names"]][[2L]]]
      stats::dnorm(.3, means, .4, log = TRUE)
    }, .is_data_multilevel = function(...) FALSE, .package = "RoBMA")
  state <- .iwmde_ordered_row_state(context, target, 1L)
  expect_false(any(vapply(state[["prior_list"]], BayesTools::is.prior.ordered, logical(1))))
  values <- c(-.4, .8)
  row <- state[["row"]]
  gamma <- spec[["allocations"]][[1L]][["gamma_coordinates"]]
  eta <- row[gamma]
  share <- eta / sum(eta)
  interaction <- specs[["mu_f__xXx__x"]]
  original_mean <- values * share[[1L]] + 1.7 * row[[interaction[["total_names"]]]] * share[[2L]]
  original <- stats::dnorm(.3, original_mean, .4, log = TRUE) + stats::dnorm(values, log = TRUE) +
    stats::dnorm(row[[interaction[["total_names"]]]], log = TRUE) +
    sum(stats::dgamma(eta, 1, 1, log = TRUE)) + stats::dnorm(row[["mu_x"]], log = TRUE)
  expect_equal(as.numeric(.iwmde_ordered_log_q_source(context, state, values)), original, tolerance = 2e-14)
  expect_equal(state[["baseline_log_q"]], .iwmde_ordered_log_q_source(context, state,
    state[["ordered_chart"]][["current"]])[[1L]], tolerance = 0)
})


test_that("completed Gaussian source mass retains residual mass and original truncation", {

  prior <- BayesTools::prior("normal", list(.3, 1.2), truncation = list(-.5, 1.4))
  current <- .2
  A <- 2.1
  H <- .5
  log_L <- function(value) .4 + H * value - .5 * A * value^2
  state <- list(ordered_chart = list(prior = prior, current = current),
    baseline_log_q = log_L(current) + BayesTools::lpdf(prior, current),
    source_gaussian = list(quadratic = A, linear = H - A * current, selected = FALSE))
  result <- .iwmde_ordered_normal_source_mass(state)
  independent <- stats::integrate(function(value) exp(log_L(value)) *
    stats::dnorm(value, .3, 1.2) / (stats::pnorm(1.4, .3, 1.2) - stats::pnorm(-.5, .3, 1.2)),
    -.5, 1.4, rel.tol = 1e-12, abs.tol = 0)
  expect_lte(abs(exp(result[["log_normalizer"]]) - independent[["value"]]),
    independent[["abs.error"]] + 64 * .Machine$double.eps * independent[["value"]])
  expect_lt(result[["original_prior_interval_log_mass"]], 0)
  expect_false(isTRUE(all.equal(result[["log_normalizer"]],
    result[["residual_log_mass"]] + result[["posterior_interval_log_mass"]])))
})


test_that("source numerical ratio gate rejects .03 independently of generic thresholds", {

  fixture <- .ordered_target_fixture()
  context <- fixture[["context"]]
  context[["data"]] <- list(outcome = data.frame(yi = .3, sei = .4))
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  target <- .iwmde_ordered_target_spec(context, "cell", list(type = "linear",
    weights = stats::setNames(1, spec[["coefficient_names"]][[1L]])))
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(context, posterior_samples, ...) {
      stats::dnorm(.3, posterior_samples[, spec[["coefficient_names"]][[1L]]], .4, log = TRUE)
    }, .is_data_multilevel = function(...) FALSE,
    .data_outcome_type = function(...) "norm", .is_data_joint_selection = function(...) FALSE,
    .iwmde_qcmde_normalization = function(...) list(summary = list(log_normalizer = 0,
      validation_log_normalizer = log(1.03), truncation = list(max = 0),
      normalization_range = c(-3, 3), normalization_points = 100L)), .package = "RoBMA")
  state <- .iwmde_ordered_row_state(context, target, 1L)
  control <- .iwmde_density_control_resolve("qCMDE", list(normalization_points = 100L), "ordinate")
  error <- tryCatch(.iwmde_ordered_log_terms(context, target, list(state), .1, control), error = identity)
  expect_s3_class(error, "iwmde_ordered_normalization_error")
  expect_s3_class(error, "iwmde_ordinate_numerical_error")
  expect_gt(error[["proof"]][["ratio_relative_error"]], .03)
})


test_that("ordered plans execute requested methods with common rows and retain underflow uncertainty", {

  fixture <- .ordered_target_fixture()
  context <- fixture[["context"]]
  context[["data"]] <- list(outcome = data.frame(yi = .3, sei = .4))
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  parameter_spec <- list(type = "linear", weights = stats::setNames(-2.5,
    spec[["coefficient_names"]][[1L]]))
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(context, posterior_samples, ...) {
      stats::dnorm(.3, posterior_samples[, spec[["coefficient_names"]][[1L]]], .4, log = TRUE)
    }, .is_data_multilevel = function(...) FALSE, .is_data_joint_selection = function(...) FALSE,
    .data_outcome_type = function(...) "norm",
    .iwmde_ordered_point_certificate = function(...) list(available = TRUE, kind = "unit_source_geometry"),
    .iwmde_plan_proposal_weight = function(...) stop("ordered engine constructed a fitted proposal"),
    .package = "RoBMA")
  plans <- lapply(c("qCMDE", "IWMDE"), function(method) .iwmde_plan(context, "cell", method,
    list(samples = 20L, normalization_points = 100L), outputs = "ordinate", values = c(.03, 50),
    parameter_spec = parameter_spec))
  expect_false(identical(plans[[1L]][["plan_key"]], plans[[2L]][["plan_key"]]))
  for (plan in plans) {
    execution <- .iwmde_plan_row_execution(context, plan)
    expect_null(execution[["proposal_weight"]])
    density <- .iwmde_execute_plan_estimator(context, plan, "ordinate", execution)
    expect_identical(density[["estimator"]], plan[["method"]])
    expect_identical(density[["internal_method"]], "ordered_source_conditional_cmde")
    expect_identical(density[["n_evaluated_rows"]], 20L)
    expect_identical(attr(density[["scaled_contributions"]], "contribution_rows", exact = TRUE),
      execution[["active_rows"]])
    expect_true(all(is.finite(density[["log_y"]])))
    expect_true(density[["underflow"]][[2L]])
    expect_identical(density[["y"]][[2L]], 0)
    expect_true(is.finite(density[["log_mcse"]][[2L]]))
    expect_true(is.finite(density[["relative_mcse"]][[2L]]))
    expect_true(is.na(density[["mcse"]][[2L]]))
  }
})


test_that("offset certificates match the supplied block covariance and prior scope", {

  fixture <- .ordered_target_fixture()
  context <- fixture[["context"]]
  context[["data"]] <- structure(list(outcome = data.frame(yi = .3, sei = .4)), outcome_type = "norm")
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  target <- .iwmde_ordered_target_spec(context, "cell", list(type = "linear",
    weights = stats::setNames(c(1, 1), spec[["coefficient_names"]][1:2])))
  state <- list(row = context[["evaluator_cache"]][[target[["ordered_source_key"]]]][1L, ],
    ordered_chart = .iwmde_ordered_row_chart(context, target, 1L), source_gaussian = list(quadratic = 1),
    active_setup = list())
  expect_identical(target[["ordered_point_plan"]][["status"]], "block")
  expect_true(.iwmde_ordered_point_certificate(context, target, list(state))[["available"]])
  for (flag in c("random", "cluster", "known_V")) {
    scoped <- context
    attr(scoped[["data"]], flag) <- TRUE
    certificate <- .iwmde_ordered_point_certificate(scoped, target, list(state))
    expect_false(certificate[["available"]])
    expect_identical(certificate[["cause"]], "covariance_certificate")
    ordinary <- state
    ordinary[["ordered_chart"]][["kind"]] <- "ordinary"
    expect_true(.iwmde_ordered_point_certificate(scoped, target, list(ordinary))[["available"]])
  }
  target[["ordered_specs"]][[1L]][["total_prior"]] <- BayesTools::prior("beta", list(2, 1))
  expect_identical(.iwmde_ordered_point_certificate(context, target, list(state))[["cause"]],
    "offset_block_capability_unavailable")
})


test_that("finite-upper Gamma inverse suppliers fail early and curves retain finite neighbors", {

  fixture <- .ordered_target_fixture(total = BayesTools::prior("gamma", list(2, 1), truncation = list(0, 10)))
  context <- fixture[["context"]]
  context[["data"]] <- structure(list(outcome = data.frame(yi = .3, sei = .4)), outcome_type = "norm")
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  parameter_spec <- list(type = "linear", weights = stats::setNames(1, spec[["coefficient_names"]][[1L]]))
  target <- .iwmde_ordered_target_spec(context, "cell", parameter_spec)
  expect_identical(target[["ordered_point_plan"]][["status"]], "capability")
  expect_identical(target[["ordered_point_plan"]][["cause"]], "truncated_gamma_inverse_supplier")
  evaluated <- NULL
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(...) 0,
    .iwmde_ordered_point_certificate = function(...) list(available = TRUE),
    .iwmde_ordered_block_geometry = function(...) stop("unsupported offset ran a denominator"),
    .iwmde_ordered_log_terms = function(context, parameter_spec, row_states, evaluation_values, control) {
      evaluated <<- evaluation_values
      list(log_terms = matrix(-1, length(evaluation_values), length(row_states)),
        proofs = lapply(row_states, function(state) list(ratio_relative_error = 0,
          numerator_relative_error = 0, denominator_relative_error = 0, exterior_relative_error = 0,
          floating_relative_error = 0, error_status = "exact", source_grid_points = 0L,
          source_grid = c(0, 10), log_normalizer = 0, normalizer = list(method = "declared_fixture"),
          candidate_rows = state[["row_index"]])),
        internal_method = "ordered_source_conditional_cmde", working_distribution = "normalized_conditional_source",
        normalization_basis = "source_chart")
    }, .package = "RoBMA")
  plan <- .iwmde_plan(context, "cell", "qCMDE", list(samples = 20L), outputs = "ordinate",
    values = 0, parameter_spec = parameter_spec)
  execution <- .iwmde_plan_row_execution(context, plan)
  expect_error(.iwmde_ordered_density(context, plan, execution, 0, FALSE),
    class = "RoBMA_ordered_truncated_gamma_inverse_supplier")
  curve <- .iwmde_ordered_density(context, plan, execution, c(0, .1, .2), TRUE)
  expect_identical(evaluated, c(.1, .2))
  expect_identical(curve[["x"]], c(.1, .2))
  expect_identical(curve[["omitted_ordered_points"]][[1L]][["cause"]], "truncated_gamma_inverse_supplier")
})


test_that("all-invalid source chunks retain a full arithmetic-scale vector", {

  state <- list(prior_evaluator = function(samples) rep(0, nrow(samples)),
    ordered_kernel = function(samples) rep(0, nrow(samples)))
  local_mocked_bindings(.iwmde_ordered_replace_source = function(context, state, values) {
    valid <- values > 0
    list(samples = matrix(values, ncol = 1L), valid = valid)
  }, .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(context, samples, ...) samples[, 1L],
    .package = "RoBMA")
  values <- c(rep(-1, 1024L), seq_len(25L), rep(-1, 1024L))
  result <- .iwmde_ordered_log_q_source(list(), state, values)
  scale <- attr(result, "max_accumulated_logterm", exact = TRUE)
  expect_length(scale, length(values))
  expect_identical(scale[values < 0], rep(0, sum(values < 0)))
  expect_identical(scale[values > 0], values[values > 0])
  state[["prior_evaluator"]] <- function(samples) rep(-Inf, nrow(samples))
  nonrepresentable <- .iwmde_ordered_log_q_source(list(), state, 1)
  expect_identical(attr(nonrepresentable, "max_accumulated_logterm", exact = TRUE), Inf)
})


test_that("construction failures use requested or neutral ordered method labels", {

  for (method in c("qCMDE", "IWMDE", "ordered_source_conditional_cmde")) {
    failure <- tryCatch(.iwmde_stop_construction_failure(method, "cell", 1L, "source chart", "unrepresentable"), error = identity)
    label <- if (method == "ordered_source_conditional_cmde") "Ordered density" else method
    expect_identical(conditionMessage(failure), paste0(label,
      " construction failed for target 'cell' at posterior row 1 during source chart: unrepresentable."))
  }
})


test_that("uncertified zero variance affects only its evaluation", {

  fixture <- .ordered_target_fixture()
  context <- fixture[["context"]]
  context[["data"]] <- structure(list(outcome = data.frame(yi = .3, sei = .4)), outcome_type = "norm")
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  parameter_spec <- list(type = "linear", weights = stats::setNames(c(1, 1), spec[["coefficient_names"]][1:2]))
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(...) 0,
    .iwmde_ordered_point_certificate = function(...) list(available = TRUE),
    .iwmde_ordered_constant_certificate = function(...) FALSE,
    .iwmde_ordered_log_terms = function(context, parameter_spec, row_states, evaluation_values, control) {
      n <- length(row_states)
      list(log_terms = rbind(rep(-1, n), log(seq_len(n) / n), log(rev(seq_len(n)) / (n + 1))),
        proofs = lapply(row_states, function(state) list(ratio_relative_error = 0,
          numerator_relative_error = 0, denominator_relative_error = 0, exterior_relative_error = 0,
          floating_relative_error = 0, error_status = "exact", source_grid_points = 0L,
          source_grid = c(-Inf, Inf), log_normalizer = 0, normalizer = list(method = "declared_fixture"),
          candidate_rows = state[["row_index"]])),
        internal_method = "ordered_source_conditional_cmde", working_distribution = "normalized_conditional_source",
        normalization_basis = "source_chart")
    }, .package = "RoBMA")
  plan <- .iwmde_plan(context, "cell", "qCMDE", list(samples = 20L), outputs = "ordinate",
    values = .1, parameter_spec = parameter_spec)
  execution <- .iwmde_plan_row_execution(context, plan)
  result <- .iwmde_ordered_density(context, plan, execution, c(-.2, -.1, .1), TRUE)
  expect_true(is.na(result[["mcse"]][[1L]]))
  expect_true(all(is.finite(result[["mcse"]][2:3])))
  expect_true(all(is.finite(result[["ess"]][2:3])))
  expect_identical(result[["mcmc_uncertainty_scope"]], "unavailable_insufficient_chain_batches")
})


test_that("analytic source normalizers bypass an unused normalization grid", {

  prior <- BayesTools::prior("normal", list(0, 1))
  state <- list(row_index = 1L, ordered_chart = list(kind = "total", source = "total",
    prior = prior, current = .2, slope = 1, offset = 0),
    source_gaussian = list(quadratic = 2, linear = -.4, selected = FALSE),
    baseline_log_q = stats::dnorm(.2, log = TRUE) - .2^2)
  local_mocked_bindings(.iwmde_normalization_grid = function(...) stop("unused grid constructed"),
    .iwmde_ordered_log_q_source = function(context, state, values) {
      result <- stats::dnorm(values, log = TRUE) - values^2
      attr(result, "max_accumulated_logterm") <- abs(result)
      result
    }, .package = "RoBMA")
  result <- .iwmde_ordered_log_terms(list(data = list(outcome = data.frame(yi = 0))),
    list(parameter = "cell"), list(state), c(-.1, .1), .iwmde_density_control_resolve("qCMDE", NULL, "ordinate"))
  expect_equal(as.numeric(result[["log_terms"]]), stats::dnorm(c(-.1, .1), 0, 1 / sqrt(3), log = TRUE), tolerance = 1e-14)
  expect_identical(result[["proofs"]][[1L]][["source_grid_points"]], 0L)
})
