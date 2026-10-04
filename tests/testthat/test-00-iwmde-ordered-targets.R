context("Ordered primitive conditional density targets")

.ordered_target_fixture <- function(total = BayesTools::prior("normal", list(0, 1)),
                                     allocation = NULL, shared = FALSE, total_shared = NULL, data = NULL) {

  if (is.null(data)) data <- data.frame(f = base::ordered(rep(c("lo", "mid", "hi", "top"), 10L),
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


test_that("short ordered contribution sequences retain ESS without inventing MCMC errors", {

  fixture <- .ordered_target_fixture()
  context <- fixture[["context"]]
  draws <- context[["posterior_samples"]]
  fit <- fixture[["fit"]]
  fit[["mcmc"]] <- coda::mcmc.list(coda::mcmc(draws[1:20, , drop = FALSE]),
    coda::mcmc(draws[21:40, , drop = FALSE]))
  fit[["sample"]] <- 20L
  fit <- BayesTools:::.bt_attach_fit_contract(BayesTools:::.bt_attach_draw_geometry(fit))
  context[["object"]][["fit"]] <- fit
  context[["chain_id"]] <- rep(1:2, each = 20L)
  context[["data"]] <- structure(list(outcome = data.frame(yi = .3, sei = .4)), outcome_type = "norm")
  context <- .iwmde_context_ensure_caches(context, reset = TRUE)
  spec <- BayesTools::JAGS_ordered_parameter_spec(fit)[[1L]]
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(context, posterior_samples, ...) {
      stats::dnorm(.3, posterior_samples[, spec[["coefficient_names"]][[1L]]], .4, log = TRUE)
    }, .is_data_multilevel = function(...) FALSE, .is_data_joint_selection = function(...) FALSE,
    .iwmde_ordered_point_certificate = function(...) list(available = TRUE), .package = "RoBMA")
  parameter_spec <- list(type = "linear", weights = stats::setNames(-2.5, spec[["coefficient_names"]][[1L]]))
  selected <- c(2L, 18L, 38L)
  for (method in c("qCMDE", "IWMDE")) {
    for (output in c("ordinate", "density")) {
      plan <- .iwmde_plan(context, "cell", method, list(samples = 20L, normalization_points = 100L),
        outputs = output, values = if (output == "ordinate") c(.03, .1, .2) else NULL,
        parameter_spec = parameter_spec)
      expect_identical(plan[["status"]], "ok")
      if (output == "density") plan[["grids"]][["display_grid"]] <- c(.03, .1, .2)
      values <- .iwmde_parameter_values(context, "cell", plan[["execution_spec"]])
      restrict <- function(plan, rows) {
        finite_rows <- seq_len(40L) %in% rows
        baseline <- .iwmde_plan_baseline_contract(context, plan, rows, values[rows])
        component <- .iwmde_restrict_parameter_component(
          .iwmde_parameter_components(context, "cell", plan[["execution_spec"]]), finite_rows)
        plan[["rows"]] <- .iwmde_plan_rows(values, finite_rows, component,
          continuous_rows = rows, candidate_rows = rows, candidate_values = values[rows],
          continuous_values = values[rows], baseline_contract = baseline,
          chain_coverage = .iwmde_plan_chain_coverage(context, rows))
        plan
      }
      # The public planner requires >=20 continuous rows. This explicitly
      # restricted internal population proves producer/schema robustness only.
      plan <- restrict(plan, selected)
      execution <- .iwmde_plan_row_execution(context, plan)
      expect_identical(execution[["chain_id"]], c(1L, 1L, 2L))
      expect_identical(execution[["expected_chain_ids"]], 1:2)
      result <- .iwmde_execute_plan_estimator(context, plan, output, execution)
      expect_identical(result[["n_evaluated_rows"]], 3L)
      expect_identical(result[["n_batches"]], 0L)
      expect_true(all(apply(result[["scaled_contributions"]], 1L, function(row) length(unique(row)) > 1L)))
      for (field in c("mcse", "relative_mcse", "log_mcse")) expect_true(all(is.na(result[[field]])))
      expect_true(all(is.finite(result[["ess"]]) & result[["ess"]] >= 0 &
        result[["ess"]] <= result[["finite_terms"]]))
      expect_identical(result[["mcmc_uncertainty_status"]], "partial")
      expect_identical(result[["mcmc_uncertainty_scope"]], "selected_continuous_rows_only")
      expect_identical(result[["mcmc_uncertainty_reason"]], "fewer than four MCMC contribution rows")
      expect_identical(.iwmde_new_density_result(result)[["ess"]], result[["ess"]])
      diagnostic <- .iwmde_plan_diagnostic_result(plan, output, execution, result)
      if (output == "ordinate") {
        point <- .iwmde_select_ordinate_diagnostic(diagnostic, .1)
        ordinate <- .iwmde_posterior_ordinate_attribute(point, method,
          metadata = list(parameter = "cell"), density_control = plan[["control"]])
        expect_true(.iwmde_posterior_ordinate_supports_bf(ordinate))
        design <- .iwmde_ordinate_sampling_design(ordinate, plan[["control"]], 40L, 3L)
        expect_false(design[["precision_target_met"]])
        expect_false(design[["target_met"]])
        public <- .iwmde_public_density_diagnostic_rows(ordinate)
        expect_true(is.finite(public[["ess"]]))
        expect_match(paste(capture.output(print(public)), collapse = " "), "fewer than four MCMC contribution rows")
      } else {
        expect_identical(.iwmde_diagnostics_density_failure_reason(diagnostic[["diagnostics"]]),
          "density plot-scale MCSE is NA. Inspect the returned condition's 'density_diagnostics' field; for a plot set 'density_method = \"KDE\"'")
        expect_null(.iwmde_posterior_density_attribute(diagnostic, method,
          metadata = list(parameter = "cell"), density_control = plan[["control"]]))
        samples <- structure(as.numeric(plan[["rows"]][["samples"]]), iwmde_diagnostics = list(diagnostic))
        plot_error <- .plot_brma_iwmde_unavailable_error(samples, method, "cell")
        expect_s3_class(plot_error, "RoBMA_density_plot_error")
        expect_identical(plot_error[["density_diagnostics"]][[1L]][["diagnostics"]][["min_ess"]], min(result[["ess"]]))
      }
      for (field in c("mcse", "relative_mcse", "log_mcse")) {
        invalid <- result
        invalid[[field]][[1L]] <- 0
        expect_error(.iwmde_new_density_result(invalid), "inconsistent partial uncertainty")
        for (value in c(NaN, Inf)) {
          invalid[[field]][[1L]] <- value
          expect_error(.iwmde_new_density_result(invalid), "invalid|inconsistent partial uncertainty")
        }
      }
      for (field in c("mcmc_uncertainty_status", "mcmc_uncertainty_reason")) {
        invalid <- result
        invalid[field] <- list(if (field == "mcmc_uncertainty_status") "available" else NULL)
        expect_error(.iwmde_new_density_result(invalid), "inconsistent partial uncertainty")
      }
      invalid <- result
      invalid[["ess"]][[1L]] <- result[["finite_terms"]][[1L]] + 1
      expect_error(.iwmde_new_density_result(invalid), "invalid effective sample sizes")
      for (rows in list(c(2L, 18L, 22L, 38L), c(2L, 6L, 12L, 18L, 22L, 26L, 32L, 38L))) {
        control <- restrict(plan, rows)
        control_execution <- .iwmde_plan_row_execution(context, control)
        control_result <- .iwmde_execute_plan_estimator(context, control, output, control_execution)
        expect_identical(control_result[["mcmc_uncertainty_status"]], if (length(rows) == 4L) "unavailable" else "available")
        for (field in c("mcse", "relative_mcse", "log_mcse", "ess")) {
          if (length(rows) == 4L) expect_true(all(is.na(control_result[[field]]))) else
            expect_true(all(is.finite(control_result[[field]])))
        }
      }
    }
  }
  contributions <- rbind(c(1, 2, 4), c(3, 1, 5))
  attr(contributions, "chain_id") <- c(1L, 1L, 2L)
  attr(contributions, "expected_chain_ids") <- 1:2
  expect_identical(.iwmde_batch_mcse(contributions)[["uncertainty_status"]], "partial")
  attr(contributions, "chain_id") <- rep(1L, 3L)
  missing <- .iwmde_batch_mcse(contributions)
  expect_identical(missing[["uncertainty_scope"]], "unavailable_missing_selected_chain")
  expect_true(all(is.na(missing[["ess"]])))
})


test_that("consistent conditioned chains retain actual active-mass uncertainty", {

  for (chain_length in c(4L, 8L, 16L)) {
    conditioned_rows <- seq_len(2L * chain_length)
    conditioned_chain_id <- rep(1:2, each = chain_length)
    for (all_active in c(FALSE, TRUE)) {
      active_per_chain <- if (all_active || chain_length == 4L) chain_length else
        as.integer(3L * chain_length / 4L)
      active_rows <- c(seq_len(active_per_chain), chain_length + seq_len(active_per_chain))
      selected_per_chain <- max(4L, as.integer(chain_length / 2L))
      selected_rows <- c(seq_len(selected_per_chain), chain_length + seq_len(selected_per_chain))
      active_mass <- length(active_rows) / length(conditioned_rows)
      values <- rbind(seq(.2, 1.3, length.out = length(selected_rows)),
        seq(1.7, .4, length.out = length(selected_rows)))
      aggregate <- .iwmde_density_aggregate(log(values), active_mass, length(selected_rows),
        contribution_rows = selected_rows, sampling_population_rows = active_rows,
        chain_id = conditioned_chain_id[selected_rows], expected_chain_ids = 1:2,
        conditioned_rows = conditioned_rows, conditioned_chain_id = conditioned_chain_id)
      expect_true(all(selected_rows %in% active_rows))
      expect_true(all(active_rows %in% conditioned_rows))
      expect_true(all(table(factor(conditioned_chain_id[selected_rows], levels = 1:2)) >= 4L))
      expect_true(all(apply(aggregate[["contributions"]], 1L, function(row) length(unique(row)) > 1L)))
      active <- .iwmde_batch_mcse(aggregate[["contributions"]])
      expect_identical(active[["uncertainty_status"]], "available")
      expect_true(all(is.finite(c(active[["mcse"]], active[["relative_mcse"]], active[["ess"]]))))
      mass <- .iwmde_active_mass_mcse(active_rows, conditioned_rows, conditioned_chain_id, 1:2, active_mass)
      expect_identical(aggregate[["active_mass_error"]], mass)
      expect_identical(mass[["uncertainty_status"]], "available")
      expect_true(all(is.finite(c(mass[["mcse"]], mass[["relative_mcse"]]))))
      if (all_active || chain_length == 4L) expect_identical(mass[["mcse"]], 0)
      primary <- .iwmde_mixture_mcse(aggregate[["contributions"]],
        aggregate[["mcmc_contributions"]], mass, active_mass)
      expect_identical(primary[["uncertainty_status"]], "available")
      expect_true(all(is.finite(c(primary[["mcse"]], primary[["relative_mcse"]], primary[["ess"]]))))
      expect_identical(primary[["active_mass_mcse"]], mass[["mcse"]])
      expect_null(primary[["uncertainty_reason"]])
    }
  }
  conditioned_rows <- 1:16
  conditioned_chain_id <- rep(1:2, each = 8L)
  active_rows <- seq(1L, 16L, by = 2L)
  selected_rows <- active_rows[1:4]
  active_mass <- length(active_rows) / length(conditioned_rows)
  aggregate <- .iwmde_density_aggregate(log(rbind(c(.2, .4, .7, 1.1), c(1.3, .9, .5, .3))),
    active_mass, length(selected_rows), contribution_rows = selected_rows,
    sampling_population_rows = active_rows, chain_id = conditioned_chain_id[selected_rows],
    expected_chain_ids = 1:2, conditioned_rows = conditioned_rows,
    conditioned_chain_id = conditioned_chain_id)
  expect_identical(aggregate[["active_mass_error"]][["uncertainty_status"]], "available")
  missing <- .iwmde_batch_mcse(aggregate[["contributions"]])
  expect_identical(missing[["uncertainty_scope"]], "unavailable_missing_selected_chain")
  primary <- .iwmde_mixture_mcse(aggregate[["contributions"]], aggregate[["mcmc_contributions"]],
    aggregate[["active_mass_error"]], active_mass)
  expect_identical(primary[["uncertainty_status"]], "unavailable")
  for (field in c("mcse", "relative_mcse", "ess")) expect_true(all(is.na(primary[[field]])))
  expect_match(primary[["uncertainty_reason"]], "no rows from fitted chain.*2")
})


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
  expect_error(.iwmde_ordered_block_geometry(context, target, state),
    class = "RoBMA_ordered_offset_block_capability_unavailable")
})


test_that("exact ordered selection grouping preserves powers and the complete likelihood", {

  dat <- data.frame(f = base::ordered(rep(c("lo", "mid", "hi", "top"), 10L), levels = c("lo", "mid", "hi", "top")))
  fixture <- .ordered_target_fixture(data = dat)
  dat[["yi"]] <- rep(c(-.1, .1, .2, .3), 10L)
  dat[["sei"]] <- rep(c(.2, .3, .25, .35), 10L)
  prior_bias <- BayesTools::prior_weightfunction("one-sided", c(.025, .05), BayesTools::wf_fixed(c(1, 3, .2)))
  object <- bselmodel(yi = yi, sei = sei, mods = ~f, data = dat, measure = "SMD", only_priors = TRUE,
    prior_effect = BayesTools::prior("point", list(0)), prior_bias = prior_bias,
    prior_mods = list(f = BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)))))
  context <- fixture[["context"]]
  context[["data"]] <- object[["data"]]
  context[["ordered_point_certificate"]] <- list(available = TRUE, kind = "independent_numeric_gaussian_core")
  state <- list(active_setup = list(is_weightfunction = TRUE, priors = object[["priors"]]))
  target <- list(ordered_specs = BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]]))
  partition <- .iwmde_ordered_selection_partition(context, target, state)
  expect_identical(partition[["representatives"]], 1:4)
  expect_identical(partition[["group_index"]], rep(1:4, 10L))
  expect_identical(partition[["multiplicity"]], rep(10L, 4L))
  selection_context <- partition[["selection_spec"]]
  selection_context[["omega"]] <- matrix(.selection_fixed_omega_branch(prior_bias, selection_context[["p_cuts"]]), 1L)
  selection_context[["alpha"]] <- 0
  selection_context[["phack_kind"]] <- 0L
  selection_context[["yi"]] <- dat[["yi"]]
  selection_context[["sei"]] <- dat[["sei"]]
  selection_context <- BayesTools::selection_context_validate(selection_context, n_samples = 1L)
  block <- list(state = state, setup = list(yi = dat[["yi"]]), selected = TRUE,
    mu0 = rep(c(.03, .04, .05, .06), 10L), mu_a = rep(c(0, 1, 1, 1), 10L),
    mu_b = rep(c(0, 0, 1, 1), 10L), mu_c = rep(c(0, .1, .2, .3), 10L),
    variance = .1^2 + dat[["sei"]]^2, selection_sei = dat[["sei"]],
    weights = rep(c(1, 3, .2), length.out = 40L), selection_context = selection_context)
  grouping <- .iwmde_ordered_selection_grouping(block, partition)
  expect_identical(grouping[["original_powers"]], block[["weights"]])
  expect_identical(grouping[["powers"]], vapply(1:4, function(g) sum(block[["weights"]][seq(g, 40L, 4L)]), numeric(1)))
  expect_identical(grouping[["original_indices"]], 1:40)
  expect_identical(grouping[["group_count"]], 4L)
  block[["selection_grouping"]] <- grouping
  block[["normalizer_context"]] <- BayesTools::selection_context_subset_observations(selection_context, 1:4)
  block[["counts"]] <- new.env(parent = emptyenv())
  for (name in c("candidate_rows", "native_nodes", "maximum_chunk")) block[["counts"]][[name]] <- 0L
  powers <- block[["weights"]] / block[["variance"]]
  residual <- dat[["yi"]] - block[["mu0"]]
  block[["gaussian"]] <- list(
    A = c(sum(powers * block[["mu_a"]]^2), 2 * sum(powers * block[["mu_a"]] * block[["mu_b"]]), sum(powers * block[["mu_b"]]^2)),
    H = c(sum(powers * residual * block[["mu_a"]]), sum(powers * (residual * block[["mu_b"]] - block[["mu_c"]] * block[["mu_a"]])),
      -sum(powers * block[["mu_c"]] * block[["mu_b"]])),
    C = c(sum(block[["weights"]] * stats::dnorm(dat[["yi"]], block[["mu0"]], sqrt(block[["variance"]]), log = TRUE)),
      sum(powers * residual * block[["mu_c"]]), -.5 * sum(powers * block[["mu_c"]]^2)))
  block[["selection_observed_constant"]] <- sum(block[["weights"]] * log(selection_context[["omega"]][1L, selection_context[["obs_bin"]]]))
  candidates <- expand.grid(total = c(-.4, 0, .7), share = c(0, .5, 1))
  mu <- matrix(block[["mu0"]], 9L, 40L, byrow = TRUE) + outer(candidates[["total"]], block[["mu_a"]]) +
    outer(candidates[["total"]] * candidates[["share"]], block[["mu_b"]]) + outer(candidates[["share"]], block[["mu_c"]])
  reference <- rowSums(.test_step_reference(dat[["yi"]], mu, matrix(sqrt(block[["variance"]]), 9L, 40L, byrow = TRUE),
    dat[["sei"]], selection_context[["omega"]][rep(1L, 9L), , drop = FALSE], selection_context, block[["weights"]]))
  actual <- .iwmde_ordered_block_log_likelihood(context, block, candidates[["total"]], candidates[["share"]])
  floor <- 64 * .Machine$double.eps * 40 * pmax(1, abs(reference), abs(actual))
  expect_true(all(abs(actual - reference) <= floor))
  expect_identical(block[["counts"]][["native_nodes"]], 9L)
  for (field in c("mu0", "mu_a", "mu_b", "mu_c", "variance", "selection_sei")) {
    changed <- block
    changed[[field]][[5L]] <- changed[[field]][[5L]] + .Machine$double.eps
    expect_null(.iwmde_ordered_selection_grouping(changed, partition))
  }
  changed <- block
  changed[["mu_a"]][c(1L, 5L)] <- c(0, -0)
  expect_null(.iwmde_ordered_selection_grouping(changed, partition))
  for (field in c("z_lower", "omega", "alpha", "obs_bin")) {
    changed <- block
    changed[["selection_context"]][[field]][[1L]] <- changed[["selection_context"]][[field]][[1L]] + .1
    expect_null(.iwmde_ordered_selection_grouping(changed, partition))
  }
  changed <- context
  changed[["data"]][["outcome"]][["yi"]][[5L]] <- changed[["data"]][["outcome"]][["yi"]][[5L]] + .Machine$double.eps
  expect_identical(length(.iwmde_ordered_selection_partition(changed, target, state)[["representatives"]]), 5L)
  for (flag in c("scale", "random", "cluster", "known_V")) {
    changed <- context
    attr(changed[["data"]], flag) <- TRUE
    expect_null(.iwmde_ordered_selection_partition(changed, target, state))
  }
  changed <- context
  changed[["ordered_point_certificate"]][["kind"]] <- "existing_ordinary_pivot"
  expect_null(.iwmde_ordered_selection_partition(changed, target, state))
  changed <- state
  changed[["active_setup"]][["priors"]][["outcome"]][["bias"]][["weights"]][["omega"]][[2L]] <- 0
  expect_null(.iwmde_ordered_selection_partition(context, target, changed))
  full_target <- .iwmde_ordered_target_spec(context, "cell", list(type = "linear",
    weights = stats::setNames(c(1, 1), target[["ordered_specs"]][[1L]][["coefficient_names"]][1:2])))
  changed[["row"]] <- context[["posterior_samples"]][1L, ]
  changed[["ordered_chart"]] <- .iwmde_ordered_row_chart(context, full_target, 1L)
  changed[["source_gaussian"]] <- list(quadratic = 1)
  refusal <- .iwmde_ordered_point_certificate(context, full_target, list(changed))
  expect_false(refusal[["available"]])
  expect_identical(refusal[["cause"]], "selection_precision")
  changed <- context
  attr(changed[["data"]], "selection_model")[["known_sampling_variance"]] <- "condition"
  expect_null(.iwmde_ordered_selection_partition(changed, target, state))
  design <- BayesTools::JAGS_formula_design(fixture[["fit"]], "mu")
  local_mocked_bindings(JAGS_formula_design = function(...) design, .package = "BayesTools")
  for (field in c("transformed_terms", "expression_specs", "expression_data", "random_effects")) {
    design[[field]] <- list(unsupported = TRUE)
    expect_null(.iwmde_ordered_selection_partition(context, target, state))
    design[[field]] <- list()
  }
})


test_that("point-total mixture rows share one normalized block aggregate", {

  # Condition on the declared within-group ratio r_J=(1/2,1/2). Integrating
  # T and S physically gives proper component evidence and posterior weights.
  # The component weights below make the two posterior branches equally likely.
  Q <- function(share) 10 / 2^2 * (1 + 1.25 * share^2)
  point_Z <- stats::integrate(function(share) stats::dbeta(share, 3, 1) * exp(-Q(share) / 2),
    0, 1, rel.tol = 1e-11)
  normal_Z <- stats::integrate(function(share) stats::dbeta(share, 3, 1) / sqrt(1 + Q(share)),
    0, 1, rel.tol = 1e-11)
  prior_weights <- c(1 / point_Z[["value"]], 1 / normal_Z[["value"]])
  prior_weights <- prior_weights / sum(prior_weights)
  total <- BayesTools::prior_mixture(list(BayesTools::prior("point", list(1), prior_weights = prior_weights[[1L]]),
    BayesTools::prior("normal", list(0, 1), prior_weights = prior_weights[[2L]])))
  fixture <- .ordered_target_fixture(total = total,
    allocation = BayesTools::prior("dirichlet", list(c(1.5, 1.5, 1))))
  context <- fixture[["context"]]
  context[["data"]] <- structure(list(outcome = data.frame(yi = rep(0, 40L), sei = rep(2, 40L))), outcome_type = "norm")
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  gamma <- spec[["allocations"]][[1L]][["gamma_coordinates"]]
  draws <- context[["posterior_samples"]]
  draws[, gamma[1:2]] <- rowSums(draws[, gamma[1:2]]) / 2
  draws <- .iwmde_ordered_replay(context, draws, list(mu_f = spec))
  context[["posterior_samples"]] <- draws
  context[["object"]][["fit"]][["mcmc"]][[1L]] <- coda::mcmc(draws)
  design <- outer(as.integer(fixture[["data"]][["f"]]), 1:3, ">")
  setup <- function(context, samples, ...) list(mu = samples[, spec[["coefficient_names"]], drop = FALSE] %*% t(design),
    yi = rep(0, 40L), sei = rep(2, 40L), tau_within = matrix(0, nrow(samples), 40L),
    effect_direction = "positive", priors = list(outcome = list()))
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(priors = list(outcome = list())),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_posterior_setup_active_branch = setup,
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(context, samples, ...) {
      mu <- setup(context, samples)[["mu"]]
      rowSums(stats::dnorm(mu, 0, 2, log = TRUE))
    }, .package = "RoBMA")
  parameter_spec <- list(type = "linear", weights = stats::setNames(c(1, 1), spec[["coefficient_names"]][1:2]))
  posterior_weights <- prior_weights * c(point_Z[["value"]], normal_Z[["value"]])
  posterior_weights <- posterior_weights / sum(posterior_weights)
  expect_equal(posterior_weights, c(.5, .5), tolerance = 1e-14)
  numerator <- stats::dnorm(0) * 3 / 2
  prior_ordinate <- prior_weights[[2L]] * numerator
  posterior_ordinate <- posterior_weights[[2L]] * numerator / normal_Z[["value"]]
  withr::defer(.iwmde_qcmde_cache_clear())
  for (method in c("qCMDE", "IWMDE")) {
    .iwmde_qcmde_cache_clear()
    plan <- .iwmde_plan(context, "cell", method, list(samples = 20L), outputs = "ordinate",
      values = 0, parameter_spec = parameter_spec)
    expect_identical(plan[["execution_spec"]][["ordered_point_plan"]][["status"]], "block")
    expect_identical(vapply(plan[["execution_spec"]][["ordered_point_plan"]][["component_statuses"]], `[[`, character(1), "status"),
      c("zero", "block"))
    expect_identical(plan[["rows"]][["point_mass_total"]], 0)
    expect_identical(plan[["rows"]][["active_mass"]], 1)
    execution <- .iwmde_plan_row_execution(context, plan)
    result <- .iwmde_execute_plan_estimator(context, plan, "ordinate", execution)
    expect_true(result[["ordered_point_certificate"]][["available"]])
    components <- vapply(execution[["row_states"]], function(state)
      .iwmde_ordered_total_component(spec, state[["row"]])[["component"]], integer(1))
    expect_identical(result[["finite_terms"]], as.numeric(sum(components == 2L)))
    expect_identical(result[["n_evaluated_rows"]], 20L)
    expect_identical(result[["sampling_fraction"]], .5)
    expect_identical(attr(result[["scaled_contributions"]], "contribution_rows", exact = TRUE), execution[["active_rows"]])
    expect_identical(attr(result[["scaled_contributions"]], "chain_id", exact = TRUE), execution[["chain_id"]])
    expect_identical(vapply(result[["ordered_proofs"]], `[[`, integer(1), "candidate_rows"), execution[["active_rows"]])
    zero <- components == 1L
    expect_identical(as.numeric(result[["scaled_contributions"]][, zero]), rep(0, sum(zero)))
    expect_true(all(vapply(result[["ordered_proofs"]][zero], function(proof)
      proof[["normalizer"]][["method"]] == "declared_zero_component_limit", logical(1))))
    block <- result[["ordered_proofs"]][[which(!zero)[[1L]]]]
    expect_false(isTRUE(block[["normalizer_cache_reused"]]))
    L0 <- 40 * stats::dnorm(0, 0, 2, log = TRUE)
    expect_lte(abs(exp(block[["log_normalizer"]] - L0) - normal_Z[["value"]]),
      normal_Z[["abs.error"]] + block[["denominator_relative_error"]] * normal_Z[["value"]])
    selected_weight <- mean(components == 2L)
    expected <- selected_weight * numerator / normal_Z[["value"]]
    expect_lte(abs(result[["y"]] - expected), expected * result[["normalization_relative_error"]] +
      normal_Z[["abs.error"]] / normal_Z[["value"]] * expected)
    expect_lte(abs(result[["y"]] - posterior_ordinate), 4 * result[["sampling_mcse"]] +
      posterior_ordinate * result[["normalization_relative_error"]])
    expect_true(is.finite(prior_ordinate) && prior_ordinate > 0)
  }
})


test_that("raw components in mixed block plans retain their conditional source route", {

  total <- BayesTools::prior_mixture(list(BayesTools::prior("point", list(1)), BayesTools::prior("gamma", list(2, 1))))
  fixture <- .ordered_target_fixture(total = total)
  context <- fixture[["context"]]
  context[["data"]] <- structure(list(outcome = data.frame(yi = .3, sei = .4)), outcome_type = "norm")
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  target <- .iwmde_ordered_target_spec(context, "cell", list(type = "linear",
    weights = stats::setNames(1, spec[["coefficient_names"]][[1L]])))
  expect_identical(vapply(target[["ordered_point_plan"]][["component_statuses"]], `[[`, character(1), "status"), c("raw", "block"))
  rows <- context[["evaluator_cache"]][[target[["ordered_source_key"]]]]
  states <- lapply(1:2, function(row) list(row = rows[row, ], row_index = row,
    ordered_chart = .iwmde_ordered_row_chart(context, target, row),
    source_gaussian = list(quadratic = 1), active_setup = list()))
  expect_true(.iwmde_ordered_point_certificate(context, target, states)[["available"]])
  called <- integer()
  local_mocked_bindings(.iwmde_ordered_log_terms = function(context, parameter_spec, row_states, ...) {
    called <<- vapply(row_states, `[[`, integer(1), "row_index")
    list(log_terms = matrix(log(2), 1L, 1L), proofs = list(list(normalizer = list(method = "raw_fixture"))))
  }, .iwmde_ordered_block_geometry = function(...) stop("raw component reached block geometry"), .package = "RoBMA")
  result <- .iwmde_ordered_block_log_terms(context, target, states[1L], target[["ordered_point_plan"]])
  expect_identical(called, 1L)
  expect_identical(result[["log_terms"]], matrix(log(2), 1L, 1L))
  expect_identical(result[["proofs"]][[1L]][["normalizer"]][["method"]], "raw_fixture")
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


test_that("requested method reaches source-chart failures through actual plans", {

  fixture <- .ordered_target_fixture()
  context <- fixture[["context"]]
  context[["data"]] <- structure(list(outcome = data.frame(yi = .3, sei = .4)), outcome_type = "norm")
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  parameter_spec <- list(type = "linear", weights = stats::setNames(c(1, 1), spec[["coefficient_names"]][1:2]))
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(...) 0,
    .package = "RoBMA")
  parameter_spec <- .iwmde_ordered_target_spec(context, "cell", parameter_spec)
  parameter_spec[["ordered_projection"]][["total_contraction_weights"]][["mu_f"]][, ] <- 0
  for (method in c("qCMDE", "IWMDE")) {
    failure <- tryCatch(.iwmde_plan(context, "cell", method, list(samples = 20L), outputs = "ordinate",
      values = .1, parameter_spec = parameter_spec), error = identity)
    expect_s3_class(failure, "iwmde_construction_error")
    expect_match(conditionMessage(failure), paste0("^", method, " construction failed for target 'cell'.*source-chart construction"))
  }
})


test_that("ordered partial uncertainty survives result assembly and public consumers", {

  fixture <- .ordered_target_fixture()
  context <- fixture[["context"]]
  context[["data"]] <- structure(list(outcome = data.frame(yi = .3, sei = .4)), outcome_type = "norm")
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  parameter_spec <- list(type = "linear", weights = stats::setNames(c(1, 1), spec[["coefficient_names"]][1:2]))
  all_constant <- FALSE
  first_constant <- TRUE
  first_log_scale <- -1
  second_log_scale <- NULL
  local_mocked_bindings(.iwmde_active_keys = function(...) rep("one", 40L),
    .iwmde_active_setup = function(...) list(),
    .iwmde_active_flat_prior_list = function(...) context[["flat_prior_list"]],
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = function(...) 0,
    .iwmde_ordered_point_certificate = function(...) list(available = TRUE),
    .iwmde_ordered_constant_certificate = function(...) FALSE,
    .iwmde_ordered_log_terms = function(context, parameter_spec, row_states, evaluation_values, control) {
      n <- length(row_states)
      terms <- if (all_constant) matrix(-1, length(evaluation_values), n) else
        rbind(if (first_constant) rep(first_log_scale, n) else log((seq_len(n) + 1L) / (n + 1L)),
          if (is.null(second_log_scale)) log(seq_len(n) / n) else
          second_log_scale + log1p(seq_len(n) / 1e5), log(rev(seq_len(n)) / (n + 1)))
      list(log_terms = terms,
        proofs = lapply(row_states, function(state) list(ratio_relative_error = 0,
          numerator_relative_error = 0, denominator_relative_error = 0, exterior_relative_error = 0,
          floating_relative_error = 0, error_status = "exact", source_grid_points = 0L,
          source_grid = c(-Inf, Inf), log_normalizer = 0, normalizer = list(method = "declared_fixture"),
          candidate_rows = state[["row_index"]])),
        internal_method = "ordered_source_conditional_cmde", working_distribution = "normalized_conditional_source",
        normalization_basis = "source_chart")
    }, .package = "RoBMA")
  for (method in c("qCMDE", "IWMDE")) {
    all_constant <- FALSE
    plan <- .iwmde_plan(context, "cell", method, list(samples = 20L), outputs = "ordinate",
      values = c(-.2, -.1, .1), parameter_spec = parameter_spec)
    execution <- .iwmde_plan_row_execution(context, plan)
    result <- .iwmde_execute_plan_estimator(context, plan, "ordinate", execution)
    expect_true(is.na(result[["mcse"]][[1L]]))
    expect_true(all(is.finite(result[["mcse"]][2:3])))
    expect_true(all(is.finite(result[["ess"]][2:3])))
    expect_identical(result[["mcmc_uncertainty_scope"]], "selected_continuous_rows_only")
    expect_identical(result[["mcmc_uncertainty_status"]], "partial")
    expect_identical(result[["mcmc_uncertainty_reason"]],
      "zero empirical contribution variance lacks a declared dependency certificate")
    expect_identical(is.na(result[["mcse"]]), is.na(result[["relative_mcse"]]))
    expect_identical(is.na(result[["mcse"]]), is.na(result[["ess"]]))
    diagnostic <- .iwmde_plan_diagnostic_result(plan, "ordinate", execution, result)
    expect_match(.iwmde_diagnostics_density_sample_failure_reason(diagnostic[["diagnostics"]]),
      "MCMC uncertainty is unavailable")
    selected <- .iwmde_select_ordinate_diagnostic(diagnostic, -.2)
    ordinate <- .iwmde_posterior_ordinate_attribute(selected, method, density_control = plan[["control"]],
      metadata = list(parameter = "cell"), allow_rejected = TRUE)
    expect_true(is.na(ordinate[["diagnostics"]][["relative_mcse"]]))
    design <- .iwmde_ordinate_sampling_design(ordinate, plan[["control"]], 40L, 20L)
    expect_false(design[["precision_target_met"]])
    expect_false(design[["target_met"]])
    expect_identical(result[["sampling_mcse"]][[1L]], 0)
    expect_identical(result[["sampling_relative_mcse"]][[1L]], 0)
    estimate <- .iwmde_estimate_fixed_ordinate(context, "cell", method, list(samples = 20L),
      values = c(-.2, -.1, .1), parameter_spec = parameter_spec, metadata = list(parameter = "cell"))
    entries <- c(.iwmde_posterior_ordinate_entries(estimate[["posterior_ordinate"]]),
      .iwmde_posterior_ordinate_entries(estimate[["rejected_posterior_ordinate"]]))
    ordinate <- entries[[which(vapply(entries, function(entry) entry[["value"]] == -.2, logical(1)))]]
    expect_false(ordinate[["diagnostics"]][["precision_target_met"]])
    expect_match(paste(.iwmde_posterior_ordinate_warnings(ordinate), collapse = " "),
      "MCSE.*unavailable")
    rows <- .iwmde_public_density_diagnostic_rows(ordinate)
    condition <- structure(list(density_diagnostics = rows), class = "RoBMA_density_ordinate_error")
    public <- density_diagnostics(condition)
    expect_true(is.na(public[["relative_mcse"]]))
    expect_true(is.na(public[["ess"]]))
    expect_identical(public[["mcmc_uncertainty_scope"]], "selected_continuous_rows_only")
    expect_match(paste(capture.output(print(public)), collapse = " "),
      "Ordered MCMC uncertainty is partial: zero empirical contribution variance lacks a declared dependency certificate")
    for (frame in list(as.data.frame(public), data.frame(public))) {
      expect_identical(frame[["relative_mcse"]], public[["relative_mcse"]])
      expect_identical(frame[["ess"]], public[["ess"]])
      expect_identical(frame[["ordinate"]], public[["ordinate"]])
    }
    # The finite partial point remains numerically accepted under the ordinary
    # BF policy. Use actual ordinary-scale underflow to exercise rejection.
    expect_length(estimate[["ordinate_failures"]][["ordinate"]], 0L)
    first_log_scale <- -900
    rejected <- .iwmde_estimate_fixed_ordinate(context, "cell", method, list(samples = 20L),
      values = c(-.2, -.1, .1), parameter_spec = parameter_spec, metadata = list(parameter = "cell"))
    rejection <- tryCatch(.iwmde_stop_ordinate_unavailable("Fixture ordinate is unavailable", rejected), error = identity)
    rejected_rows <- density_diagnostics(rejection)
    expect_s3_class(rejected_rows, "BayesTools_table")
    expect_identical(rejected_rows[["numerical_status"]], "underflow")
    expect_match(paste(capture.output(print(rejected_rows)), collapse = " "),
      "Ordered MCMC uncertainty is partial: zero empirical contribution variance lacks a declared dependency certificate")
    first_constant <- FALSE
    ordinary <- .iwmde_estimate_fixed_ordinate(context, "cell", method, list(samples = 20L),
      values = c(-.2, -.1, .1), parameter_spec = parameter_spec, metadata = list(parameter = "cell"))
    ordinary_rows <- .iwmde_public_density_diagnostic_rows(ordinary[["posterior_ordinate"]])
    expect_identical(class(ordinary_rows), c("RoBMA_density_diagnostics", "data.frame"))
    expect_null(attr(ordinary_rows, "footnotes", exact = TRUE))
    finite_row <- ordinary_rows[2L, , drop = FALSE]
    expect_true(all(is.finite(c(finite_row[["ordinate"]], finite_row[["relative_mcse"]], finite_row[["ess"]]))))
    make_bf <- function(label, diagnostics) {
      out <- data.frame(BF = 1.125, row.names = label)
      class(out) <- c("BayesTools_hypothesis_BF", "BayesTools_table", "data.frame")
      attr(out, "raw_BF") <- 1.125
      attr(out, "hypothesis_ast") <- structure(list(statements = list(list())), class = "BayesTools_hypothesis_ast")
      attr(out, "density_diagnostics") <- diagnostics
      out
    }
    for (frames in list(list(rejected_rows, finite_row), list(finite_row, rejected_rows))) {
      mixed <- .iwmde_bind_density_diagnostics(frames)
      expect_identical(attr(mixed, "footnotes"), attr(rejected_rows, "footnotes"))
      expect_null(attr(mixed, "type", exact = TRUE))
      expect_identical(mixed[["requested_value"]], unlist(lapply(frames, `[[`, "requested_value"), use.names = FALSE))
      for (field in names(finite_row)) {
        expect_identical(mixed[[field]], c(frames[[1L]][[field]], frames[[2L]][[field]]))
      }
      expect_match(paste(capture.output(print(mixed)), collapse = " "), "Ordered MCMC uncertainty is partial")
      for (converted in list(as.data.frame(mixed), data.frame(mixed))) {
        expect_identical(names(converted), names(mixed))
        for (field in names(mixed)) expect_identical(converted[[field]], mixed[[field]])
      }
      subset <- mixed[2:1, c("ordinate", "relative_mcse", "ess"), drop = FALSE]
      expect_identical(attr(subset, "footnotes"), attr(mixed, "footnotes"))
      expect_identical(subset[["ordinate"]], mixed[["ordinate"]][2:1])
      combined <- .hypothesis_brma_bind_parameter_results(
        list(make_bf("one", frames[[1L]]), make_bf("two", frames[[2L]])), list(1L, 2L),
        structure(list(statements = list(list(), list())), class = "BayesTools_hypothesis_ast"))
      expect_identical(density_diagnostics(combined)[["ordinate"]], mixed[["ordinate"]])
      expect_identical(attr(density_diagnostics(combined), "footnotes"), attr(mixed, "footnotes"))
      expect_match(paste(capture.output(print(density_diagnostics(combined))), collapse = " "),
        "Ordered MCMC uncertainty is partial")
    }
    expect_identical(attr(.iwmde_bind_density_diagnostics(list(rejected_rows, rejected_rows)), "footnotes"),
      attr(rejected_rows, "footnotes"))
    empty <- .iwmde_empty_public_density_diagnostics()
    attr(empty, "footnotes") <- c("Existing empty-frame note.", attr(rejected_rows, "footnotes"))
    expect_identical(attr(.iwmde_bind_density_diagnostics(list(finite_row, rejected_rows, empty)), "footnotes"),
      c(attr(rejected_rows, "footnotes"), "Existing empty-frame note."))
    numerical_error <- tryCatch(.iwmde_stop_ordinate_numerical_failure(c(0, .525), c(-900, log(.525)),
      values = c(-.2, -.1), finite_terms = c(20, 20), max_weight_share = c(.05, .1)), error = identity)
    enriched <- tryCatch(.iwmde_enrich_ordinate_error(numerical_error, rejected[["plan"]], context, "cell",
      method, rejected[["plan"]][["control"]], c(-.2, -.1), metadata = list(parameter = "cell")), error = identity)
    expect_identical(density_diagnostics(enriched)[["numerical_status"]], c("underflow", "finite"))
    expect_identical(class(density_diagnostics(enriched)), c("RoBMA_density_diagnostics", "data.frame"))
    first_constant <- TRUE
    first_log_scale <- -1
    available <- .iwmde_select_ordinate_diagnostic(diagnostic, -.1)
    available <- .iwmde_posterior_ordinate_attribute(available, method, density_control = plan[["control"]],
      metadata = list(parameter = "cell"), allow_rejected = TRUE)
    expect_identical(available[["diagnostics"]][["relative_mcse"]], result[["relative_mcse"]][[2L]])
    expect_identical(available[["diagnostics"]][["ess"]], result[["ess"]][[2L]])
    underflow <- result
    underflow[["y"]][[2L]] <- 0
    underflow[["log_y"]][[2L]] <- -900
    underflow[["mcse"]][[2L]] <- NA_real_
    underflow[["log_mcse"]][[2L]] <- -950
    underflow[["underflow"]][[2L]] <- TRUE
    expect_identical(.iwmde_new_density_result(underflow)[["relative_mcse"]], result[["relative_mcse"]])
    subnormal <- underflow
    subnormal[["y"]][[2L]] <- exp(-740)
    subnormal[["log_y"]][[2L]] <- -740
    subnormal[["underflow"]][[2L]] <- FALSE
    expect_gt(subnormal[["y"]][[2L]], 0)
    expect_identical(.iwmde_new_density_result(subnormal)[["ess"]], result[["ess"]])
    malformed <- subnormal
    malformed[["log_mcse"]][[2L]] <- -700
    expect_gt(exp(malformed[["log_mcse"]][[2L]]), 0)
    expect_error(.iwmde_new_density_result(malformed), "inconsistent partial uncertainty")
    for (field in c("mcmc_uncertainty_scope", "mcmc_uncertainty_status", "mcmc_uncertainty_reason", "internal_method")) {
      for (value in list(character(), NA_character_, c("one", "two"))) {
        invalid <- result
        invalid[field] <- list(value)
        expect_error(.iwmde_new_density_result(invalid), "invalid uncertainty metadata")
      }
    }
    invalid <- result
    invalid["mcmc_uncertainty_scope"] <- list(NULL)
    expect_error(.iwmde_new_density_result(invalid), "invalid uncertainty metadata")
    invalid <- result
    invalid["internal_method"] <- list(NULL)
    expect_error(.iwmde_new_density_result(invalid), "invalid uncertainty metadata")
    invalid <- result
    invalid[["mcmc_uncertainty_scope"]] <- NULL
    expect_error(.iwmde_new_density_result(invalid), "missing required field")
    # Execute the actual density-output assembly and the curve consumers.
    curve_plan <- .iwmde_plan(context, "cell", method,
      list(samples = 20L), outputs = "density", parameter_spec = parameter_spec)
    curve_plan[["grids"]][["display_grid"]] <- c(-.2, -.1, .1)
    curve_execution <- .iwmde_plan_row_execution(context, curve_plan)
    curve_result <- .iwmde_execute_plan_estimator(context, curve_plan, "density", curve_execution)
    curve_diagnostic <- .iwmde_plan_diagnostic_result(curve_plan, "density", curve_execution, curve_result)
    expect_identical(curve_result[["mcmc_uncertainty_status"]], "partial")
    expect_identical(curve_result[["mcse"]][2:3], result[["mcse"]][2:3])
    expect_identical(curve_result[["ess"]][2:3], result[["ess"]][2:3])
    expect_match(.iwmde_diagnostics_density_failure_reason(curve_diagnostic[["diagnostics"]]),
      "MCMC uncertainty is unavailable")
    expect_null(.iwmde_posterior_density_attribute(curve_diagnostic, method,
      metadata = list(parameter = "cell"), density_control = curve_plan[["control"]]))
    samples <- structure(as.numeric(curve_plan[["rows"]][["samples"]]), iwmde_diagnostics = list(curve_diagnostic))
    unavailable_curve <- .plot_brma_iwmde_unavailable_reason(samples)
    expect_true(unavailable_curve[["rejected"]])
    expect_match(unavailable_curve[["reason"]], "MCMC uncertainty is unavailable")
    expect_identical(.iwmde_density_sample_rejection_action(curve_diagnostic[["diagnostics"]]),
      "Inspect the returned condition's 'density_diagnostics' field; for a plot set 'density_method = \"KDE\"'")
    expected_plain_reason <- paste0(
      "density MCMC uncertainty is unavailable for some evaluation values. ",
      "For a plot set 'density_method = \"KDE\"'")
    expect_identical(.iwmde_diagnostics_density_failure_reason(curve_diagnostic[["diagnostics"]],
      has_diagnostics_field = FALSE), expected_plain_reason)
    # Controlled consumer coverage with the actually produced partial record;
    # this does not establish a public ordered-partial marginal-means route.
    marginal_samples <- with_draw_metadata(as.numeric(curve_plan[["rows"]][["samples"]]),
      linear_weights = parameter_spec[["weights"]])
    marginal_object <- list(
      inference = list(averaged = list(mu_f = list(mid = marginal_samples))),
      parameters = "mu_f",
      term_map = data.frame(term = "f", parameter = "mu_f", label = "f", stringsAsFactors = FALSE),
      density_diagnostics = list(averaged = list(`f: mid` = curve_diagnostic)))
    expect_identical(.marginal_means_iwmde_density_rejection_reasons(marginal_object,
      selected = list(term = "f"), type = "averaged", levels = "mid"), expected_plain_reason)
    ordinary_diagnostics <- curve_diagnostic[["diagnostics"]]
    ordinary_diagnostics[["internal_method"]] <- "ordinary"
    for (all_rows in c(TRUE, FALSE)) {
      ordinary_diagnostics[["all_rows_used"]] <- all_rows
      ordinary_diagnostics[["sampling_fraction"]] <- if (all_rows) 1 else .5
      expected_action <- if (all_rows) "Try fitting the model with more posterior draws" else
        "Try a larger 'samples' value in the 'density_control' argument"
      for (has_field in c(TRUE, FALSE)) {
        expect_identical(.iwmde_density_sample_rejection_action(ordinary_diagnostics,
          has_diagnostics_field = has_field), expected_action)
      }
      expect_identical(.iwmde_diagnostics_density_failure_reason(ordinary_diagnostics,
        has_diagnostics_field = FALSE), .iwmde_diagnostics_density_failure_reason(ordinary_diagnostics))
    }
    plot_error <- .plot_brma_iwmde_unavailable_error(samples, method, "cell")
    expect_s3_class(plot_error, "RoBMA_density_plot_error")
    expect_identical(conditionMessage(plot_error), paste0(method,
      " density for 'cell' was rejected by diagnostics: density MCMC uncertainty is unavailable for some evaluation values. ",
      "Inspect the returned condition's 'density_diagnostics' field; for a plot set 'density_method = \"KDE\"'."))
    expect_identical(plot_error[["density_diagnostics"]][[1L]][["diagnostics"]], curve_diagnostic[["diagnostics"]])
    for (scale in c(-740, -900)) {
      second_log_scale <- scale
      represented <- .iwmde_execute_plan_estimator(context, curve_plan, "density", curve_execution)
      expect_identical(represented[["mcmc_uncertainty_status"]], "partial")
      expect_true(is.na(represented[["mcse"]][[2L]]))
      expect_true(is.finite(represented[["log_mcse"]][[2L]]))
      expect_identical(exp(represented[["log_mcse"]][[2L]]), 0)
      expect_true(is.finite(represented[["relative_mcse"]][[2L]]))
      expect_true(is.finite(represented[["ess"]][[2L]]))
      if (scale == -740) expect_gt(represented[["y"]][[2L]], 0) else expect_identical(represented[["y"]][[2L]], 0)
      expect_identical(represented[["underflow"]][[2L]], scale == -900)
    }
    second_log_scale <- NULL
    zero_plan <- plan
    zero_plan[["execution_spec"]][["ordered_point_plan"]][["status"]] <- "zero"
    zero_plan[["execution_spec"]][["ordered_point_plan"]][["value"]] <- -.1
    zero_result <- .iwmde_execute_plan_estimator(context, zero_plan, "ordinate", execution)
    expect_identical(zero_result[["log_y"]][[2L]], -Inf)
    expect_identical(zero_result[["mcse"]][[2L]], 0)
    expect_identical(zero_result[["ess"]][[2L]], 0)
    expect_true(is.na(zero_result[["relative_mcse"]][[2L]]))
    expect_identical(zero_result[["mcmc_uncertainty_status"]], "partial")
    for (field in c("mcse", "relative_mcse", "ess")) {
      invalid <- result
      invalid[[field]][[2L]] <- Inf
      expect_error(.iwmde_new_density_result(invalid), "invalid")
      invalid[[field]][[2L]] <- NaN
      expect_error(.iwmde_new_density_result(invalid), "invalid")
      invalid <- result
      invalid[[field]][[1L]] <- 1
      expect_error(.iwmde_new_density_result(invalid), "inconsistent partial uncertainty")
    }
    for (field in c("mcmc_uncertainty_status", "mcmc_uncertainty_scope", "mcmc_uncertainty_reason", "internal_method")) {
      invalid <- result
      invalid[field] <- list(switch(field, mcmc_uncertainty_status = "available",
        mcmc_uncertainty_scope = "unavailable_insufficient_chain_batches",
        mcmc_uncertainty_reason = NULL, internal_method = "ordinary"))
      expect_error(.iwmde_new_density_result(invalid), "invalid (effective sample sizes|uncertainty metadata)")
    }
    all_constant <- TRUE
    unavailable <- .iwmde_execute_plan_estimator(context, plan, "ordinate", execution)
    expect_true(all(is.na(unavailable[["mcse"]])))
    expect_true(all(is.na(unavailable[["relative_mcse"]])))
    expect_true(all(is.na(unavailable[["ess"]])))
    expect_identical(unavailable[["mcmc_uncertainty_scope"]], "unavailable_insufficient_chain_batches")
    expect_identical(unavailable[["mcmc_uncertainty_status"]], "unavailable")
  }
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


.ordered_narrow_gaussian_fixture <- function(repetitions = 1L, selected = FALSE, direction = "positive",
                                             total = BayesTools::prior("normal", list(0, 1)), current = .2) {

  sign <- if (direction == "negative") -1 else 1
  dat <- data.frame(f = base::ordered(rep(c("lo", "mid", "top"), repetitions),
    levels = c("lo", "mid", "hi", "top")),
    yi = sign * rep(c(0, 1e-4, 1), repetitions), sei = rep(c(1e-4, 1e-6, 1e-4), repetitions))
  allocation <- BayesTools::prior("dirichlet", list(c(2, 2, 1)))
  fixture <- .ordered_target_fixture(total = total, data = dat, allocation = allocation)
  priors <- list(intercept = BayesTools::prior("point", list(0)),
    f = BayesTools::prior_ordered(total, allocation = allocation))
  arguments <- list(yi = dat[["yi"]], sei = dat[["sei"]], mods = ~f, data = dat, measure = "SMD",
    only_priors = TRUE, prior_effect = BayesTools::prior("point", list(0)),
    prior_heterogeneity = BayesTools::prior("point", list(0)), prior_mods = priors)
  if (selected) {
    arguments[["effect_direction"]] <- direction
    arguments[["prior_bias"]] <- BayesTools::prior_weightfunction("one-sided", c(.025, .05), BayesTools::wf_fixed(c(1, .7, .3)))
  }
  object <- do.call(if (selected) bselmodel else brma, arguments)
  context <- fixture[["context"]]
  object[["fit"]] <- fixture[["fit"]]
  context[["object"]] <- object
  context[["data"]] <- object[["data"]]
  context[["priors"]] <- object[["priors"]]
  context[["posterior_samples"]] <- cbind(context[["posterior_samples"]], tau = 0)
  spec <- BayesTools::JAGS_ordered_parameter_spec(fixture[["fit"]])[[1L]]
  context[["posterior_samples"]][1L, spec[["total_names"]]] <- sign * current
  context[["posterior_samples"]][1L, spec[["allocations"]][[1L]][["gamma_coordinates"]]] <- c(1e-4, .5 * (1 - 1e-4), .5 * (1 - 1e-4))
  context[["posterior_samples"]] <- .iwmde_ordered_replay(context, context[["posterior_samples"]], list(mu_f = spec))
  context[["object"]][["fit"]][["mcmc"]][[1L]] <- coda::mcmc(context[["posterior_samples"]])
  context <- .iwmde_context_ensure_caches(context, reset = TRUE)
  target <- .iwmde_ordered_target_spec(context, "cell", list(type = "linear",
    weights = stats::setNames(1, spec[["coefficient_names"]][[1L]])))
  state <- .iwmde_ordered_row_state(context, target, 1L)
  list(context = context, target = target, state = state, data = dat,
    block = .iwmde_ordered_block_geometry(context, target, state))
}


test_that("conditioned Gaussian references retain every original node and witness", {

  references <- readRDS(testthat::test_path("fixtures", "ordered-gaussian-reference.rds"))
  cases <- c("normal-1-positive", "normal-1-negative", "normal-30-positive", "normal-30-negative",
    "selected-1-positive", "selected-1-negative", "selected-30-positive", "selected-30-negative",
    "boundary-1", "boundary-2", "boundary-3")
  expect_identical(sort(names(references)), sort(cases))
  expect_identical(unname(vapply(references[cases], function(x) length(x[["total"]]), integer(1))),
    c(rep(31L, 8L), rep(3L, 3L)))
  witnesses <- vapply(references, function(x) sum(abs(x[["rounded_replay"]] - x[["declared"]]) > x[["allowance"]]), integer(1))
  expect_identical(sum(witnesses), 38L)
  expect_identical(sum(witnesses[startsWith(names(witnesses), "boundary-")]), 2L)
})


test_that("one-ulp alternate replay inputs retain the portable conditioning relation", {

  for (case in c("normal-1-positive", "selected-1-positive", "boundary-1")) {
    reference <- .ordered_gaussian_reference(case)
    expect_length(reference[["alternate_inputs"]], 4L)
    for (alternate in reference[["alternate_inputs"]]) {
      node <- alternate[["node"]]
      point <- reference
      point[["total"]] <- reference[["total"]][node]
      for (field in c("declared_mean", "declared_mean_radius", "replay_mean", "replay_sd")) {
        point[[field]] <- point[[field]][node, , drop = FALSE]
      }
      mean <- alternate[["mean"]]
      sd <- alternate[["sd"]]
      expect_false(.ordered_gaussian_replay_matches(point, mean, sd))
      input <- reference[["inputs"]]
      weights <- input[["powers"]]
      selection <- reference[["selection"]]
      if (is.null(selection)) {
        evaluated <- sum(weights * stats::dnorm(input[["yi"]], mean, sd, log = TRUE))
      } else {
        selection[["n_bins"]] <- length(selection[["z_lower"]])
        evaluated <- sum(.test_step_reference(input[["yi"]], mean, sd, input[["selection_sei"]],
          selection[["omega"]], selection, weights))
      }
      if (!startsWith(case, "boundary-")) evaluated <- evaluated + stats::dnorm(reference[["total"]][[node]], log = TRUE)
      allowance <- reference[["allowance"]][[node]]
      expect_lte(abs(evaluated - alternate[["oracle"]]), allowance)
      bound <- .ordered_gaussian_rounding_bound(point, weights, input[["variance"]], mean, sd, selection)
      # The MP120 shift and analytic perturbation proof are independent of this
      # double-precision bound and the replay's arithmetic.
      expect_lte(alternate[["actual_shift"]], alternate[["analytic_bound"]])
      expect_lte(alternate[["analytic_bound"]], bound)
      expect_lte(abs(evaluated - reference[["declared"]][[node]]), bound + allowance)
      changed_mean <- mean
      changed_mean[[length(changed_mean)]] <- changed_mean[[length(changed_mean)]] + .001
      changed_sd <- sd * 1.01
      expect_true(any(.ordered_gaussian_replay_caps(point, changed_mean, sd)[["mean_deviation"]] >
        .ordered_gaussian_replay_caps(point, changed_mean, sd)[["mean_cap"]]))
      expect_true(any(.ordered_gaussian_replay_caps(point, mean, changed_sd)[["variance_deviation"]] >
        .ordered_gaussian_replay_caps(point, mean, changed_sd)[["variance_cap"]]))
    }
  }
})


test_that("near-fit ordered Gaussian and observed step factors match full original replay", {

  grid <- expand.grid(total = c(-1, 0, 1), share = c(1e-6, 1e-4, .001, .5, 1))
  for (selected in c(FALSE, TRUE)) for (repetitions in c(1L, 30L)) for (direction in c("positive", "negative")) {
    fixture <- .ordered_narrow_gaussian_fixture(repetitions, selected, direction)
    context <- fixture[["context"]]
    block <- fixture[["block"]]
    state <- fixture[["state"]]
    samples <- .iwmde_ordered_block_samples(context, block, grid[["total"]], grid[["share"]])
    original <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(context, samples, state[["active_setup"]], "estimate")
    actual <- .iwmde_ordered_block_log_likelihood(context, block, grid[["total"]], grid[["share"]])
    allowance <- 64 * .Machine$double.eps * nrow(fixture[["data"]]) * pmax(1, abs(original), abs(actual))
    expect_true(all(abs(actual - original) <= allowance))
    if (selected) {
      selection <- block[["selection_context"]]
      observed <- sum(block[["weights"]] * log(selection[["omega"]][1L, selection[["obs_bin"]]]))
      expect_identical(block[["selection_observed_constant"]], observed)
      expect_identical(.iwmde_ordered_selection_observed_constant(context, state, block[["setup"]], selection), observed)
      for (field in c("omega", "obs_bin", "alpha", "phack_kind", "kernel_mode")) {
        changed <- selection
        changed[[field]][[1L]] <- changed[[field]][[1L]] + 1
        error <- tryCatch(.iwmde_ordered_selection_observed_constant(context, state, block[["setup"]], changed), error = identity)
        expect_identical(class(error)[[1L]], "RoBMA_ordered_selection_route")
        expect_s3_class(error, "RoBMA_density_method_unavailable")
        expect_false(inherits(error, "iwmde_ordinate_numerical_error"))
      }
      changed <- block[["setup"]]
      changed[["selection_sei"]][[1L]] <- 2 * changed[["selection_sei"]][[1L]]
      error <- tryCatch(.iwmde_ordered_selection_observed_constant(context, state, changed, selection), error = identity)
      expect_identical(class(error)[[1L]], "RoBMA_ordered_selection_route")
      expect_s3_class(error, "RoBMA_density_method_unavailable")
      expect_false(inherits(error, "iwmde_ordinate_numerical_error"))
      changed <- block[["setup"]]
      changed[["tau_within"]][[1L]] <- Inf
      error <- tryCatch(.iwmde_ordered_selection_observed_constant(context, state, changed, selection), error = identity)
      expect_identical(class(error)[[1L]], "RoBMA_ordered_selection_route")
      expect_s3_class(error, "iwmde_ordinate_numerical_error")
      expect_false(inherits(error, "RoBMA_density_method_unavailable"))
    }
    case <- paste(if (selected) "selected" else "normal", repetitions, direction, sep = "-")
    node_reference <- .ordered_gaussian_reference(case)
    .ordered_gaussian_reference_inputs(node_reference, block)
    expect_identical(node_reference[["share"]], 1e-4)
    expect_identical(block[["prior"]], node_reference[["prior"]])
    expect_identical(context[["priors"]], node_reference[["priors"]])
    expect_identical(fixture[["data"]], node_reference[["data"]])
    expect_identical(block[["chart"]][names(node_reference[["chart"]])], node_reference[["chart"]])
    nodes <- node_reference[["total"]]
    expect_length(nodes, 31L)
    node_samples <- .iwmde_ordered_block_samples(context, block, nodes, 1e-4)
    expect_identical(node_samples[, block[["chart"]][["gamma_coordinates"]], drop = FALSE], node_reference[["gamma"]])
    node_original <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(context, node_samples, state[["active_setup"]], "estimate")
    node_actual <- .iwmde_ordered_block_log_likelihood(context, block, nodes, 1e-4)
    h_original <- node_original + stats::dnorm(nodes, log = TRUE)
    h_actual <- node_actual + BayesTools::lpdf(block[["prior"]], nodes)
    node_allowance <- node_reference[["allowance"]]
    expect_true(all(abs(h_actual - node_reference[["declared"]]) <= node_allowance))
    replay_setup <- .iwmde_log_lik_posterior_setup_active_branch(context, node_samples, state[["active_setup"]], "estimate")
    replay_sd <- .root_sum_squares(replay_setup[["tau_within"]],
      matrix(replay_setup[["sei"]], length(nodes), nrow(fixture[["data"]]), byrow = TRUE))
    if (.ordered_gaussian_replay_matches(node_reference, replay_setup[["mu"]], replay_sd)) {
      expect_true(all(abs(h_original - node_reference[["rounded_replay"]]) <= node_allowance))
    }
    selection_fields <- names(node_reference[["selection"]])
    selection <- if (selected) block[["selection_context"]] else NULL
    if (selected) expect_identical(selection[selection_fields], node_reference[["selection"]])
    rounding_bound <- .ordered_gaussian_rounding_bound(node_reference, block[["weights"]], block[["variance"]],
      replay_setup[["mu"]], replay_sd, selection)
    if (.ordered_gaussian_replay_matches(node_reference, replay_setup[["mu"]], replay_sd)) {
      expect_true(all(node_reference[["forward_bound"]] <= rounding_bound))
    }
    expect_true(all(abs(h_original - node_reference[["declared"]]) <= rounding_bound + node_allowance))
  }
})


test_that("actual unselected source completion retains nuisance and truncation masses", {

  cases <- list(list(mean = 0, sd = 1, lower = -Inf, upper = Inf, current = .2),
    list(mean = .3, sd = 1.2, lower = .999, upper = .9995, current = .99925))
  names(cases) <- c("untruncated", "tight-truncated")
  cases[["tail-upper-0.9992"]] <- list(mean = 0, sd = 1, lower = -Inf, upper = .9992, current = .999)
  references <- readRDS(testthat::test_path("fixtures", "ordered-normal-source-reference.rds"))
  for (name in names(cases)) {
    case <- cases[[name]]
    oracle <- references[[name]]
    expect_identical(case, oracle[["case"]])
    prior <- BayesTools::prior("normal", list(case[["mean"]], case[["sd"]]),
      truncation = list(case[["lower"]], case[["upper"]]))
    fixture <- .ordered_narrow_gaussian_fixture(total = prior, current = case[["current"]])
    context <- fixture[["context"]]
    state <- fixture[["state"]]
    dat <- fixture[["data"]]
    expect_identical(dat, oracle[["data"]])
    expect_identical(dat[["yi"]], oracle[["yi"]])
    expect_identical(dat[["sei"]], oracle[["sei"]])
    expect_false(state[["source_gaussian"]][["selected"]])
    expect_true(.iwmde_ordered_source_residual_available(context, state))
    expect_identical(state[["ordered_chart"]][["prior"]], prior)
    expect_identical(prior, oracle[["prior"]])
    expect_identical(state[["ordered_chart"]][["kind"]], "total")
    expect_identical(state[["ordered_chart"]][["current"]], case[["current"]])
    eta <- state[["row"]][state[["ordered_specs"]][[1L]][["allocations"]][[1L]][["gamma_coordinates"]]]
    expect_identical(unname(eta), oracle[["eta"]])
    share <- eta[[1L]] / sum(eta)
    prior_constant <- sum(stats::dgamma(eta, shape = c(2, 2, 1), rate = 1, log = TRUE))
    expect_lte(abs(prior_constant - oracle[["nuisance"]]),
      64 * .Machine$double.eps * max(1, abs(oracle[["nuisance"]])))
    original_reference <- stats::integrate(function(total) stats::dnorm(total, case[["mean"]], case[["sd"]]),
      case[["lower"]], case[["upper"]], rel.tol = 1e-12, abs.tol = 0)
    original_mass <- original_reference[["value"]]
    original_log_error <- -log1p(-original_reference[["abs.error"]] / original_mass)
    source <- function(total) {
      mu <- outer(total, c(0, share, 1))
      rowSums(stats::dnorm(matrix(dat[["yi"]], length(total), nrow(dat), byrow = TRUE),
        mu, matrix(dat[["sei"]], length(total), nrow(dat), byrow = TRUE), log = TRUE)) +
        stats::dnorm(total, case[["mean"]], case[["sd"]], log = TRUE) - log(original_mass) + prior_constant
    }
    precision <- share^2 / 1e-12 + 1 / 1e-8 + 1 / case[["sd"]]^2
    center <- (share * 1e-4 / 1e-12 + 1 / 1e-8 + case[["mean"]] / case[["sd"]]^2) / precision
    sd <- 1 / sqrt(precision)
    anchor_point <- if (is.finite(case[["upper"]])) case[["upper"]] else center
    anchor <- source(anchor_point)
    reference <- stats::integrate(function(z) exp(source(center + sd * z) - anchor),
      (case[["lower"]] - center) / sd, (case[["upper"]] - center) / sd,
      rel.tol = 1e-10, abs.tol = 0, subdivisions = 1000L)
    expected <- anchor + log(sd) + log(reference[["value"]])
    actual <- .iwmde_ordered_normal_source_mass(state, context)
    expect_identical(actual[["method"]], "completed_gaussian_source")
    expect_identical(reference[["message"]], "OK")
    prior_interval <- .ordered_normal_interval_arithmetic_allowance(case[["lower"]], case[["upper"]],
      case[["mean"]], case[["sd"]])
    # Account for the actual completed interval's inputs, without requiring
    # derived means or BLAS predictor values to match a captured platform.
    kernel <- state[["source_gaussian"]]
    completed_precision <- kernel[["quadratic"]] + 1 / case[["sd"]]^2
    completed_mean <- case[["current"]] +
      (kernel[["linear"]] + (case[["mean"]] - case[["current"]]) / case[["sd"]]^2) / completed_precision
    completed_interval <- .ordered_normal_interval_arithmetic_allowance(case[["lower"]], case[["upper"]],
      completed_mean, 1 / sqrt(completed_precision))
    coordinate_error <- 0
    if (is.finite(case[["lower"]]) && is.finite(case[["upper"]])) {
      gamma2 <- 2 * .Machine$double.eps / (1 - 2 * .Machine$double.eps)
      delta_total <- gamma2 * (abs(center) + max(abs(case[["lower"]] - center), abs(case[["upper"]] - center)))
      basis <- c(0, share, 1)
      derivative <- sum(abs(basis) * pmax(abs(dat[["yi"]] - basis * case[["lower"]]),
        abs(dat[["yi"]] - basis * case[["upper"]])) / dat[["sei"]]^2) +
        max(abs(case[["lower"]] - case[["mean"]]), abs(case[["upper"]] - case[["mean"]])) / case[["sd"]]^2
      coordinate_error <- delta_total * derivative + .5 * precision * delta_total^2
    }
    residual_error <- 64 * .Machine$double.eps * nrow(dat) *
      max(1, actual[["residual_log_scale"]], abs(expected), abs(actual[["log_normalizer"]]))
    producer_error <- prior_interval[["log_error"]] + completed_interval[["log_error"]] + residual_error
    reference_error <- coordinate_error + reference[["abs.error"]] / reference[["value"]] +
      original_log_error + residual_error
    terms <- c(prior_interval = prior_interval[["log_error"]], completed_interval = completed_interval[["log_error"]],
      reference_coordinate = coordinate_error, quadrature = reference[["abs.error"]] / reference[["value"]],
      prior_integral = original_log_error, residual = residual_error)
    expect_true(all(is.finite(terms) & terms >= 0))
    allowance <- sum(terms)
    expect_lte(abs(actual[["log_normalizer"]] - expected), allowance)
    expect_lte(abs(actual[["log_normalizer"]] - oracle[["source"]]), producer_error)
    floating_floor <- 64 * .Machine$double.eps * nrow(dat) * actual[["max_log_scale"]]
    expect_lte(abs(actual[["log_normalizer"]] - oracle[["source"]]), floating_floor)
    expect_true(is.finite(floating_floor))
    expect_lt(floating_floor, .025)
    expect_lte(abs(expected - oracle[["source"]]), reference_error)
    expect_lte(abs(actual[["posterior_interval_log_mass"]] - oracle[["posterior_mass"]]),
      completed_interval[["log_error"]] + residual_error)
    prior_residual_error <- 64 * .Machine$double.eps * max(1, abs(log(original_mass)))
    expect_lte(abs(actual[["original_prior_interval_log_mass"]] - log(original_mass)),
      prior_interval[["log_error"]] + original_log_error + prior_residual_error)
    expect_lte(abs(actual[["original_prior_interval_log_mass"]] - oracle[["prior_mass"]]),
      prior_interval[["log_error"]] + prior_residual_error)
    expect_lte(abs(log(original_mass) - oracle[["prior_mass"]]), original_log_error + prior_residual_error)
    expect_identical(original_reference[["message"]], "OK")
    row <- matrix(state[["row"]], 1L, dimnames = list(NULL, names(state[["row"]])))
    replayed <- .iwmde_ordered_replay(context, row, state[["ordered_specs"]])
    nuisance <- state[["prior_evaluator"]](replayed) + state[["ordered_kernel"]](replayed) -
      BayesTools::lpdf(prior, state[["ordered_chart"]][["current"]])
    expect_lte(abs(nuisance - prior_constant), 64 * .Machine$double.eps * max(1, abs(prior_constant)))
    if (is.finite(case[["upper"]])) {
      if (identical(name, "tight-truncated")) expect_lt(allowance, 1e-9)
      expect_gt(center, case[["upper"]])
      expect_false(any(.iwmde_ordered_replace_source(context, state, center)[["valid"]]))
      expect_lt(actual[["posterior_interval_log_mass"]], 0)
      expect_lt(actual[["original_prior_interval_log_mass"]], 0)
    }
  }
})


test_that("actual selected source normalizers retain the anchored prior constant", {

  for (repetitions in c(1L, 30L)) {
    fixture <- .ordered_narrow_gaussian_fixture(repetitions, selected = TRUE)
    context <- fixture[["context"]]
    state <- fixture[["state"]]
    dat <- fixture[["data"]]
    selection <- fixture[["block"]][["selection_context"]]
    eta <- state[["row"]][state[["ordered_specs"]][[1L]][["allocations"]][[1L]][["gamma_coordinates"]]]
    # Retained primitive Gamma priors are counted once, independently of the
    # repeated outcome rows and of the focal total's Normal density.
    prior_constant <- sum(stats::dgamma(eta, shape = c(2, 2, 1), rate = 1, log = TRUE))
    share <- eta[[1L]] / sum(eta)
    source <- function(total) {
      mu <- outer(total, rep(c(0, share, 1), repetitions))
      sd <- matrix(dat[["sei"]], length(total), nrow(dat), byrow = TRUE)
      rowSums(.test_step_reference(dat[["yi"]], mu, sd, dat[["sei"]],
        selection[["omega"]][rep(1L, length(total)), , drop = FALSE], selection)) +
        stats::dnorm(total, log = TRUE) + prior_constant
    }
    # Independently complete the physical residual quadratic solely to choose
    # a coordinate for the full original selected-density integral.
    precision <- repetitions * (share^2 / 1e-12 + 1 / 1e-8) + 1
    center <- repetitions * (share * 1e-4 / 1e-12 + 1 / 1e-8) / precision
    sd <- 1 / sqrt(precision)
    anchor <- source(center)
    reference <- stats::integrate(function(z) exp(source(center + sd * z) - anchor),
      -Inf, Inf, rel.tol = 1e-9, abs.tol = 1e-11, subdivisions = 1000L)
    expected <- anchor + log(sd) + log(reference[["value"]])
    actual <- .iwmde_ordered_source_normalizer(context, state)
    expect_identical(actual[["method"]], "completed_gaussian_gh")
    allowance <- reference[["abs.error"]] / reference[["value"]] +
      64 * .Machine$double.eps * nrow(dat) * max(1, abs(expected), abs(actual[["log_normalizer"]]))
    expect_lte(abs(actual[["log_normalizer"]] - expected), allowance)
    expect_identical(reference[["message"]], "OK")
  }
})
