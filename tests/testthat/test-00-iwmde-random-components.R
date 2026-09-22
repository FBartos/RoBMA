test_that("shared random inclusion preserves SD atoms and defined proportions", {

  samples <- cbind(tau = 1:4, gate = c(0, 1, 0, 1),
                   `weight[1]` = .25, `weight[2]` = .75,
                   `prior_par_eta_weight[1]` = 1,
                   `prior_par_eta_weight[2]` = 3)
  gate <- list(weight_name = NULL, index = NA_integer_,
               scale = "total_variance", n_targets = 1L,
               inclusion_name = "gate")
  weight <- list(weight_name = "weight", index = 1L,
                 scale = "total_variance", n_targets = 2L)
  allocation <- list(
    source = list(name = "tau", shape = "scalar"),
    scale = "total_variance", target = "block", n_targets = 2L,
    weight_name = "weight", inclusion = list(),
    parent_factors = list(gate), factors = list(gate, weight)
  )
  term <- list(
    block_name = "study", sd_component_terms = "intercept",
    sd_binding = list(true_allocation = TRUE, allocations = list(allocation))
  )
  fit <- structure(list(), prior_list = list(
    weight = BayesTools::prior("dirichlet", list(alpha = c(1, 1)))
  ), formula_design = list(mu = list(parameter = "mu", random_effects = list(term))))
  selected <- list(
    spec = list(source_type = "composite", source_parameter = "",
                source_transform = "identity", quantity = "sd_total",
                evaluator = "allocation_sd", formula_parameter = "mu",
                block = "study", random_component = "intercept",
                allocation_derived = TRUE, allocation_index = 1L),
    allocation_definition = allocation
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .get_posterior_samples = function(...) samples,
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    random_effects_marginal_update_plan = function(...) list(family = "unsupported"),
    .package = "BayesTools"
  )
  context <- list(posterior_samples = samples, flat_prior_list = list())
  expected <- list(sd_total = c(0, 2, 0, 4), sd = c(0, 1, 0, 2),
                   var_prop = c(NA, .25, NA, .25))

  for (quantity in names(expected)) {
    selected[["spec"]][["quantity"]] <- quantity
    selected[["spec"]][["evaluator"]] <- switch(
      quantity, sd_total = "allocation_sd", sd = "sd", var_prop = "allocation"
    )
    selected[["spec"]][["source_transform"]] <- if (quantity == "var_prop") {
      "var_prop"
    } else "identity"
    target <- .brma_random_parameter_density_target(list(fit = fit), quantity)
    expect_null(target[["reason"]], info = quantity)
    spec <- .iwmde_parameter_spec(context, target[["parameter"]], target[["parameter_spec"]])
    expect_identical(spec[["status"]], "ok")
    expect_equal(.iwmde_parameter_values(context, target[["parameter"]], spec),
                 expected[[quantity]], info = quantity)
    component <- .iwmde_parameter_components(context, target[["parameter"]], spec)
    expect_identical(component[["active"]], c(FALSE, TRUE, FALSE, TRUE))
    expect_identical(.iwmde_parameter_condition_rows(context, spec),
                     if (quantity == "var_prop") c(FALSE, TRUE, FALSE, TRUE) else rep(TRUE, 4))
    if (quantity == "var_prop") {
      expect_equal(nrow(component[["point_masses"]]), 0L)
    } else {
      expect_equal(component[["point_masses"]], data.frame(x = 0, mass = .5))
    }
    expect_identical(.iwmde_plan_parameter_spec(spec)[["gate_metadata"]],
                     spec[["gate_metadata"]])
    ungated <- spec
    ungated[["gate_metadata"]] <- NULL
    expect_false(identical(.iwmde_target_key(target[["parameter"]], spec),
                           .iwmde_target_key(target[["parameter"]], ungated)))
  }
})

test_that("shared-gate allocation grids use declared continuous covariance plans", {

  mods <- data.frame(study = factor(c("a", "a", "b", "b")),
                      esid = factor(1:4))
  result <- BayesTools::JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag") +
      random(1 | esid, name = "esid", covariance = "diag"),
    parameter = "mu",
    data = mods,
    prior_list = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(allocation = list(
      BayesTools::random_variance_allocation(
        name = "root", terms = c(component = "component_gated"),
        sd = BayesTools::prior("gamma", list(2, 2)),
        inclusion = list(component = BayesTools::prior("spike", list(location = .5)))
      ),
      BayesTools::random_variance_allocation(
        name = "split", terms = c(study = "study", esid = "esid"),
        parent = BayesTools::allocation_ref("root", "component"),
        weights = BayesTools::prior("dirichlet", list(alpha = c(1, 1)))
      )
    ))
  )
  design <- result$formula_design
  root   <- design$random_allocations$root
  split  <- design$random_allocations$split
  source <- root$source_node
  gate   <- root$inclusion$component$indicator_name
  weight <- split$weight_name
  samples <- matrix(c(.1, .4, 1, .2, .8, -.1, .6, 0, .7, .3,
                       .2, .8, 1, .6, .4), nrow = 3L, byrow = TRUE,
                    dimnames = list(NULL, c("mu_intercept", source, gate,
                                            paste0(weight, "[", 1:2, "]"))))
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = design)
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)
  catalog <- BayesTools::parameter_catalog(fit)
  prior_selection <- BayesTools::parameter_catalog_resolve(
    catalog, "(mu) split: var_prop(study)"
  )
  selected <- list(
    spec = list(quantity = "var_prop", allocation_index = 1L,
                label = "split: tau2_prop(study)"),
    entry = list(parameter = "(mu) split: tau2_prop(study)",
                 selection = prior_selection),
    allocation_definition = split,
    samples = as.matrix(BayesTools::parameter_draws(fit, prior_selection))
  )
  withr::local_seed(1)
  prior_seed <- .Random.seed
  prior_samples <- .brma_random_parameter_mixed_posterior(
    list(fit = fit), "split: tau2_prop(study)", prior = TRUE,
    selected = selected
  )[[1L]]
  prior_density <- attr(prior_samples, "prior_density", exact = TRUE)
  expect_s3_class(prior_density, "prior_linear_density")
  expect_equal(vapply(c(0, .5, 1), function(value) {
    BayesTools::prior_density_ordinate(prior_density, value)$log_density
  }, numeric(1)), rep(0, 3L))
  expect_equal(as.numeric(prior_samples), c(.2, .6))
  expect_identical(.Random.seed, prior_seed)
  targets <- list(
    list(selector = "(mu) split: sd_total", index = NA_integer_, type = "scalar"),
    list(selector = "(mu) esid: sd(intercept)", index = 2L, type = "random_component_sd"),
    list(selector = "(mu) split: var_prop(study)", index = 1L, type = "simplex_pair"),
    list(selector = "(mu) split: var_prop(esid)", index = 2L, type = "simplex_pair")
  )
  values <- c(0, .4, 1)
  rows   <- c(3L, 1L)
  states <- lapply(rows, function(row) list(row_index = row))
  yi     <- c(.1, -.2, .3, .05)
  context <- list(posterior_samples = samples,
                   data = list(outcome = data.frame(yi = yi)))
  testthat::local_mocked_bindings(
    .iwmde_uses_known_v_random_marginal_likelihood = function(...) TRUE,
    .is_data_weights = function(...) FALSE,
    .fitted_formula_design = function(...) design,
    .iwmde_known_v_random_marginal_setup = function(...) list(
      sampling_covariance = sampling,
      dependency_blocks = list(1:2, 3:4),
      group_iid_plan_cache = new.env(parent = emptyenv())
    ),
    .iwmde_predictor_evaluate_fixed_mu = function(context, active_setup, samples) {

      matrix(samples[, "mu_intercept"], nrow = nrow(samples), ncol = 4L)
    },
    .iwmde_predictor_log_prior = function(context, parameter, values, row_states, ...) {

      numeric(length(values) * length(row_states))
    },
    .data_effect_direction = function(...) "positive",
    .package = "RoBMA"
  )

  for (target in targets) {
    selection <- BayesTools::parameter_catalog_resolve(
      catalog, target$selector, "mu"
    )
    update <- BayesTools::random_effects_marginal_update_plan(fit, selection)
    parameter <- if (identical(target$type, "simplex_pair")) {
      paste0(weight, "[", target$index, "]")
    } else source
    replacement <- list(type = target$type, covariance_update = update)
    if (identical(target$type, "simplex_pair")) {
      replacement <- c(replacement, list(parameter = weight,
        index = target$index, n_targets = 2L))
    } else if (identical(target$type, "random_component_sd")) {
      replacement <- c(replacement, list(source_parameter = source,
        factors = list(list(weight_name = weight, index = 2L,
          n_targets = 2L, scale = "total_variance"))))
    }
    plan <- .iwmde_known_v_random_group_iid_plan(context, parameter, replacement)
    expect_identical(plan$required_active, gate)

    for (correlated in c(FALSE, TRUE)) {
      sampling <- diag(c(.04, .05, .06, .07))
      if (correlated) {
        sampling[1, 2] <- sampling[2, 1] <- .01
        sampling[3, 4] <- sampling[4, 3] <- .02
      }
      actual <- .iwmde_log_q_grid_known_v_random_group_iid(
        context, parameter, values, states, replacement, list()
      )
      expected <- matrix(NA_real_, nrow = length(values), ncol = length(rows))
      for (j in seq_along(rows)) {
        row <- rows[[j]]
        for (i in seq_along(values)) {
          sd <- if (target$type == "scalar") values[[i]] else samples[row, source]
          p  <- samples[row, paste0(weight, "[1]")]
          if (target$type == "simplex_pair") {
            p <- if (target$index == 1L) values[[i]] else 1 - values[[i]]
          } else if (target$type == "random_component_sd") {
            sd <- values[[i]] / sqrt(1 - p)
          }
          covariance <- sampling + sd^2 * ((1 - p) * diag(4L) +
            p * outer(mods$study, mods$study, "=="))
          residual <- yi - samples[row, "mu_intercept"]
          factor   <- chol(covariance)
          solved   <- forwardsolve(t(factor), residual)
          expected[i, j] <- -.5 * (4 * log(2 * pi) +
            2 * sum(log(diag(factor))) + sum(solved^2))
        }
      }
      expect_equal(actual, expected, tolerance = 1e-12)
    }
    expect_equal(.iwmde_log_q_grid_known_v_random_group_iid(
      context, parameter, values, states[1L], replacement, list()
    ), actual[, 1L, drop = FALSE], tolerance = 1e-12)
    expect_null(.iwmde_log_q_grid_known_v_random_group_iid(
      context, parameter, values, list(list(row_index = 2L)), replacement, list()
    ))
    undeclared <- replacement
    undeclared$covariance_update$conditional <- NULL
    expect_null(.iwmde_known_v_random_group_iid_plan(context, parameter, undeclared))
  }
  missing_gate <- context
  missing_gate$posterior_samples <- samples[, colnames(samples) != gate, drop = FALSE]
  expect_null(.iwmde_log_q_grid_known_v_random_group_iid(
    missing_gate, parameter, values, states, replacement, list()
  ))
  replacement$covariance_update$conditional$required_active <- c(gate, "another_gate")
  expect_null(.iwmde_known_v_random_group_iid_plan(context, parameter, replacement))
})


test_that("simplex endpoint replacement uses its conditional beta prior", {

  parameter <- "rho"
  eta <- .iwmde_simplex_auxiliary_columns(parameter, 2L)
  samples <- matrix(
    c(0, 0, 4, 1, 2, 0),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("theta", eta))
  )
  prior_list <- list(
    theta = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    rho   = BayesTools::prior("dirichlet", list(alpha = c(1, 1)))
  )
  replacement <- list(
    type              = "simplex_pair",
    parameter         = parameter,
    index             = 1L,
    auxiliary_columns = eta
  )

  actual <- .iwmde_replacement_log_prior_rows(
    samples     = samples,
    prior_list  = prior_list,
    replacement = replacement
  )
  expected <- stats::dnorm(samples[, "theta"], log = TRUE) +
    stats::dbeta(c(0, 1), 1, 1, log = TRUE)

  expect_equal(actual, expected, tolerance = 0)
  expect_true(all(is.finite(actual)))
})

test_that("simplex component replacement supports more than two targets", {

  parameter <- "weight"
  columns   <- paste0(parameter, "[", 1:3, "]")
  eta       <- .iwmde_simplex_auxiliary_columns(parameter, 3L)
  row <- c(theta = 0.2, stats::setNames(c(.2, .3, .5), columns),
           stats::setNames(c(2, 3, 5), eta))
  replacement <- list(
    type              = "simplex_pair",
    parameter         = parameter,
    index             = 2L,
    n_targets         = 3L,
    auxiliary_columns = eta
  )

  replaced <- .iwmde_replace_row_for_value(
    context     = list(object = NULL),
    state       = list(row = row),
    parameter   = columns[[2L]],
    value       = 0.6,
    replacement = replacement
  )

  expect_true(replaced[["valid"]])
  expect_equal(unname(replaced[["row"]][columns]), c(.114285714, .6, .285714286))
  expect_equal(sum(replaced[["row"]][eta]), sum(row[eta]))
  expect_equal(
    replaced[["row"]][eta[[1L]]] / replaced[["row"]][eta[[3L]]],
    row[eta[[1L]]] / row[eta[[3L]]]
  )

  prior_list <- list(
    theta  = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    weight = BayesTools::prior("dirichlet", list(alpha = c(1, 2, 3)))
  )
  samples <- matrix(
    replaced[["row"]],
    nrow = 1L,
    dimnames = list(NULL, names(replaced[["row"]]))
  )
  expect_true(is.finite(.iwmde_replacement_log_prior_rows(
    samples,
    prior_list,
    replacement
  )))
})

test_that("simplex conditional kernels include the allocation Jacobian", {

  weights <- c(.05, .2, .6, .95)
  cases <- list(
    list(alpha = c(2, 3), density = function(w) 12 * w * (1 - w)^2),
    list(alpha = rep(1, 4L), density = function(w) 3 * (1 - w)^2),
    list(alpha = rep(.5, 3L), density = function(w) .5 / sqrt(w))
  )
  for (case in cases) {
    alpha <- case[["alpha"]]
    n     <- length(alpha)
    eta   <- .iwmde_simplex_auxiliary_columns("weight", n)
    other <- seq_len(n - 1L) / sum(seq_len(n - 1L))
    samples <- 4 * cbind(weights, outer(1 - weights, other))
    colnames(samples) <- eta
    prior_list <- list(weight = BayesTools::prior("dirichlet", list(alpha)))
    replacement <- list(type = "simplex_pair", parameter = "weight",
                         index = 1L, auxiliary_columns = eta)
    actual <- exp(.iwmde_replacement_log_prior_rows(
      samples, prior_list, replacement
    ))
    expect_equal(actual, case[["density"]](weights), tolerance = 1e-12)

    # Independently transform the joint gamma density. For two components the
    # Jacobian is constant, so its normalized conditional density is unchanged.
    gamma_kernel <- function(w) {

      coordinates <- 4 * cbind(w, outer(1 - w, other))
      log_density <- rowSums(vapply(seq_len(n), function(j) {
        stats::dgamma(coordinates[, j], alpha[[j]], log = TRUE)
      }, numeric(length(w))))
      exp(log_density) * (1 - w)^(n - 2L)
    }
    normalizer <- stats::integrate(gamma_kernel, 0, 1, rel.tol = 1e-10)$value
    expect_equal(actual, gamma_kernel(weights) / normalizer, tolerance = 1e-10)
  }

  samples <- rbind(c(0, 1, 3), c(4, 0, 0))
  colnames(samples) <- eta
  expect_equal(.iwmde_replacement_log_prior_rows(
    samples, prior_list, replacement
  ), c(Inf, log(.5)), tolerance = 1e-14)
})

test_that("simplex baselines use the same conditional prior as replacements", {

  eta <- .iwmde_simplex_auxiliary_columns("weight", 3L)
  samples <- cbind(center = c(.2, -.4), rbind(c(1, 2, 1), c(3, 2, 1)))
  colnames(samples)[-1L] <- eta
  priors <- list(
    center = BayesTools::prior("normal", list(0, 1)),
    weight = BayesTools::prior("dirichlet", list(alpha = c(2, 1, 1)))
  )
  context <- .iwmde_context_ensure_caches(structure(
    list(posterior_samples = samples, flat_prior_list = priors),
    class = "iwmde_context"
  ))
  spec <- list(type = "simplex_pair", parameter = "weight", index = 1L,
               auxiliary_columns = eta)
  testthat::local_mocked_bindings(
    .iwmde_active_setup = function(...) list(),
    .iwmde_row_parameters = function(...) list(),
    .iwmde_baseline_log_likelihood = function(...) 0,
    .iwmde_log_lik_from_posterior_samples_sum_active_branch =
      function(context, posterior_samples, ...) numeric(nrow(posterior_samples)),
    .package = "RoBMA"
  )
  grouped <- .iwmde_row_states_grouped_marginal(context, 1:2, "weight[1]", spec)
  scalar  <- .iwmde_row_states(context, 1:2, "weight[1]", spec)
  w <- samples[, eta[[1L]]] / rowSums(samples[, eta, drop = FALSE])
  expected <- stats::dnorm(samples[, "center"], log = TRUE) + log(6 * w * (1 - w))
  expect_equal(vapply(grouped, `[[`, numeric(1), "baseline_log_q"), expected)
  expect_equal(vapply(scalar, `[[`, numeric(1), "baseline_log_q"), expected)
  expect_equal(.iwmde_replacement_log_prior_rows(samples, priors, spec), expected)
})

test_that("retained selection baselines preserve local states and scalar diagnostics", {

  withr::local_seed(17)
  original_rng_kind <- RNGkind()
  withr::defer(do.call(RNGkind, as.list(original_rng_kind)))

  samples <- cbind(mu = c(.2, -.1, .4, .5, -.3),
                   gamma = c(.1, .3, -.2, .7, -.5),
                   mu_indicator = c(1, 2, 1, 2, 1))
  priors <- list(
    mu = BayesTools::prior_mixture(list(
      BayesTools::prior("normal", list(0, 1)),
      BayesTools::prior("normal", list(0, 2))
    )),
    gamma = BayesTools::prior("normal", list(0, 1))
  )
  data <- structure(list(), selection_model = structure(list(
    schema_version = 3L, estimate_random_effects = "integrate",
    other_random_effects = "condition", known_sampling_variance = "integrate",
    applicability = list(estimate_random_effects = FALSE,
                         other_random_effects = TRUE,
                         known_sampling_variance = FALSE)
  ), class = "RoBMA_selection_model"))
  new_context <- function() {

    .iwmde_context_ensure_caches(structure(list(
      data = data, posterior_samples = samples, flat_prior_list = priors,
      indicator_names = "mu_indicator"
    ), class = "iwmde_context"))
  }
  mode <- "ok"
  evaluations <- list()
  likelihood <- function(context, posterior_samples, active_setup, ...) {

    evaluations[[length(evaluations) + 1L]] <<-
      match(posterior_samples[, "mu"], samples[, "mu"])
    if (nrow(posterior_samples) > 1L && mode != "ok") {
      RNGkind("Wichmann-Hill")
      stats::runif(3L)
      if (mode == "error") stop("Synthetic batch failure.", call. = FALSE)
      if (mode == "warning") warning("Synthetic batch warning.", call. = FALSE)
      if (mode == "nonfinite") return(rep(NA_real_, nrow(posterior_samples)))
      if (mode == "shape") return(matrix(0, nrow(posterior_samples), 1L))
    }
    if (nrow(posterior_samples) == 1L && posterior_samples[1L, "mu"] == -.1) {
      if (mode == "warning") warning("Synthetic row warning.", call. = FALSE)
      if (mode == "error") stop("Synthetic row failure.", call. = FALSE)
    }
    stats::dnorm(.15, posterior_samples[, "mu"] + posterior_samples[, "gamma"],
                 .3 + .1 * active_setup[["branch"]], log = TRUE)
  }
  testthat::local_mocked_bindings(
    .iwmde_active_setup = function(context, row, ...) {
      list(branch = row[["mu_indicator"]])
    },
    .iwmde_row_parameters = function(context, row, ...) {
      list(mu = row[["mu"]], gamma = row[["gamma"]])
    },
    .iwmde_log_lik_from_posterior_samples_sum_active_branch = likelihood,
    .package = "RoBMA"
  )
  rows <- c(5L, 1L, 4L, 2L, 1L)
  evaluate <- function(context, batch) {

    warnings <- character()
    value <- withCallingHandlers(tryCatch({
      if (batch) {
        .iwmde_row_states(context, rows, "mu", estimator = "q_grid_cmde")
      } else {
        lapply(rows, function(row) {
          .iwmde_row_states(context, row, "mu", estimator = "q_grid_cmde")[[1L]]
        })
      }
    }, error = identity), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings, kind = RNGkind(), seed = .Random.seed)
  }
  for (mode in c("ok", "warning", "error", "nonfinite", "shape")) {
    scalar_context <- new_context()
    batch_context  <- new_context()
    do.call(RNGkind, as.list(original_rng_kind))
    set.seed(17)
    scalar <- evaluate(scalar_context, FALSE)
    evaluations <- list()
    do.call(RNGkind, as.list(original_rng_kind))
    set.seed(17)
    batched <- evaluate(batch_context, TRUE)
    expect_identical(batched[["warnings"]], scalar[["warnings"]], info = mode)
    expect_identical(batched[["kind"]], scalar[["kind"]], info = mode)
    expect_identical(batched[["seed"]], scalar[["seed"]], info = mode)
    if (mode == "error") {
      expected <- paste0("qCMDE construction failed for target 'mu' at posterior row 2 ",
        "during baseline joint-density evaluation: Synthetic row failure.")
      expect_identical(conditionMessage(batched[["value"]]), expected)
      expect_identical(conditionMessage(scalar[["value"]]), expected)
      expect_null(conditionCall(batched[["value"]]))
    } else {
      fields <- c("row_index", "row", "parameters", "prior_list",
        "baseline_log_lik", "baseline_log_prior", "baseline_focal_log_prior",
        "use_focal_prior_delta", "baseline_log_q", "likelihood_mode", "state_scope")
      expect_equal(lapply(batched[["value"]], `[`, fields),
                   lapply(scalar[["value"]], `[`, fields), tolerance = 0, info = mode)
      expected_lik <- stats::dnorm(.15, samples[rows, "mu"] + samples[rows, "gamma"],
        .3 + .1 * samples[rows, "mu_indicator"], log = TRUE)
      expected_prior <- stats::dnorm(samples[rows, "mu"], 0,
        samples[rows, "mu_indicator"], log = TRUE) +
        stats::dnorm(samples[rows, "gamma"], log = TRUE)
      expect_equal(vapply(batched[["value"]], `[[`, numeric(1), "baseline_log_q"),
                   expected_lik + expected_prior, tolerance = 1e-14)
      expect_true(all(vapply(batched[["value"]], function(state) {
        identical(state[["likelihood_mode"]], "conditional") &&
          identical(state[["state_scope"]], "local")
      }, logical(1))))
      if (mode == "ok") {
        expected_batches <- list(c(5L, 1L), c(4L, 2L))
        expect_identical(evaluations, expected_batches)
        repeated <- evaluate(batch_context, TRUE)
        expect_identical(evaluations, expected_batches)
        expect_equal(lapply(repeated[["value"]], `[`, fields),
                     lapply(batched[["value"]], `[`, fields), tolerance = 0)
      }
    }
  }
})

test_that("component allocations apply their own leaf weights", {

  posterior <- cbind(
    tau    = c(.8, 1.2),
    `weight[1]` = c(.25, .4),
    `weight[2]` = c(.75, .6)
  )
  term <- list(
    block_name = "study",
    sd_binding = list(
      true_allocation = TRUE,
      allocations = list(list(
        target               = "sd_component",
        source               = list(name = "tau", shape = "scalar"),
        parent_factors       = list(),
        weight_name          = "weight",
        scale                = "mean_variance",
        n_targets            = 2L,
        leaf_index_by_column = 1:2
      ))
    )
  )

  actual <- .marginalized_random_effect_allocated_sd_samples(
    term,
    posterior,
    K = 2L
  )

  expect_equal(
    unname(actual),
    unname(posterior[, "tau"] *
      sqrt(2 * posterior[, paste0("weight[", 1:2, "]")]))
  )
  expect_equal(
    vapply(
      .marginalized_random_effect_allocation_factors(term, all = TRUE),
      `[[`,
      integer(1),
      "index"
    ),
    1:2
  )
})

test_that("semantic density transformations apply their Jacobians", {

  source_x <- c(-.5, 0, .5)
  density <- list(
    x = source_x,
    y = c(.2, .4, .2),
    point_masses = data.frame(x = .25, mass = .1)
  )
  transform <- list(type = "tanh")
  actual <- .plot_brma_transform_iwmde_density(density, transform)

  expect_equal(actual[["x"]], tanh(source_x))
  expect_equal(actual[["y"]], density[["y"]] / (1 - tanh(source_x)^2))
  expect_equal(actual[["point_masses"]][["x"]], tanh(.25))

  ordinate <- list(
    value            = .5,
    evaluation_value = .5,
    ordinate         = .3,
    diagnostics      = list(mcse = .03, relative_mcse = .1),
    iwmde_provenance = list(value = .5, evaluation_value = .5)
  )
  transformed <- .hypothesis_brma_transform_iwmde_ordinate(
    ordinate,
    transform
  )
  jacobian <- 1 - tanh(.5)^2

  expect_equal(transformed[["value"]], tanh(.5))
  expect_equal(transformed[["ordinate"]], .3 / jacobian)
  expect_equal(transformed[["diagnostics"]][["mcse"]], .03 / jacobian)
  expect_equal(transformed[["diagnostics"]][["relative_mcse"]], .1)
})

test_that("diagonal allocation grid equals the full marginal covariance", {

  yi <- c(0.1, -0.2, 0.3, 0.05)
  vi <- c(0.04, 0.05, 0.06, 0.07)
  cluster_map <- c(1L, 1L, 2L, 2L)
  posterior_samples <- matrix(
    c(0.2, 0.45),
    ncol = 1L,
    dimnames = list(NULL, "tau")
  )
  mu_samples <- rbind(
    c(0.02, 0.02, -0.01, -0.01),
    c(-0.03, 0.01, 0.04, 0.02)
  )
  values <- c(0, 0.4, 1)
  row_states <- list(list(row_index = 1L), list(row_index = 2L))
  context <- list(
    posterior_samples = posterior_samples,
    data = list(outcome = data.frame(yi = yi))
  )
  testthat::local_mocked_bindings(
    .iwmde_known_v_random_group_iid_plan = function(...) {
      list(
        source_parameter   = "tau",
        target_mode        = "proportion",
        target_index       = 2L,
        unique_factor      = 1L,
        group_maps         = list(seq_along(cluster_map), cluster_map),
        cluster_weight_index = 2L,
        source_to_total_sd = 1,
        factor_indices     = 1:2
      )
    },
    .iwmde_known_v_random_marginal_setup = function(...) {
      list(
        sampling_covariance = diag(vi),
        dependency_blocks   = unname(split(seq_along(cluster_map), cluster_map)),
        group_iid_plan_cache = new.env(parent = emptyenv())
      )
    },
    .iwmde_predictor_evaluate_fixed_mu = function(...) mu_samples,
    .iwmde_predictor_log_prior = function(...) {
      numeric(length(values) * length(row_states))
    },
    .data_effect_direction = function(...) "positive",
    .package = "RoBMA"
  )

  actual <- .iwmde_log_q_grid_known_v_random_group_iid(
    context      = context,
    parameter    = "rho[2]",
    values       = values,
    row_states   = row_states,
    replacement  = list(),
    active_setup = list()
  )
  expected <- matrix(NA_real_, nrow = length(values), ncol = nrow(mu_samples))
  same_cluster <- outer(cluster_map, cluster_map, "==")
  for (value_i in seq_along(values)) {
    for (draw in seq_len(nrow(mu_samples))) {
      tau <- posterior_samples[draw, "tau"]
      covariance <- diag(vi + tau^2 * (1 - values[[value_i]])) +
        tau^2 * values[[value_i]] * same_cluster
      expected[value_i, draw] <- .marglik_mvn_log_density(
        y          = yi,
        mean       = mu_samples[draw, ],
        covariance = covariance
      )
    }
  }

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("semantic random qCMDE hypotheses use the plotting density target", {

  selected <- list(
    entry = list(term = "study variance proportion"),
    spec = list(
      quantity     = "var_prop",
      source_parameter = "rho",
      label            = "tau2_prop(study)"
    ),
    samples      = matrix(seq(0.01, 0.99, length.out = 100L), ncol = 1L),
    prior        = BayesTools::prior("beta", list(alpha = 1, beta = 1)),
    source_prior = BayesTools::prior(
      "dirichlet",
      list(alpha = c(1, 1))
    )
  )
  attached_values <- numeric()
  target_spec <- list(
    type                 = "simplex_pair",
    parameter            = "rho",
    index                = 2L,
    n_targets            = 2L,
    target_columns       = paste0("rho[", 1:2, "]"),
    auxiliary_columns    = paste0("rho_eta[", 1:2, "]"),
    conditioning_exclude = paste0("rho[", 1:2, "]")
  )
  used_density_method <- NULL
  attachment_calls    <- 0L
  reused_selected     <- NULL
  reused_prior        <- NULL
  semantic_prior_density <- BayesTools::prior(
    "uniform",
    list(a = 0, b = 1)
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .brma_random_parameter_density_target = function(...) {
      list(parameter = "rho[2]", parameter_spec = target_spec)
    },
    .brma_random_parameter_mixed_posterior = function(
        ..., selected = NULL, prior_selected = NULL) {
      reused_selected <<- selected
      reused_prior    <<- prior_selected
      values <- 1:3
      attr(values, "prior_density") <- semantic_prior_density
      list(theta = values)
    },
    .iwmde_context = function(...) list(),
    .iwmde_estimate_cache = function(...) new.env(parent = emptyenv()),
    .hypothesis_brma_attach_iwmde_scalar = function(
        posterior, raw_posterior, value, parameter, parameter_spec, ...) {
      attachment_calls <<- attachment_calls + 1L
      attached_values <<- c(attached_values, value)
      expect_identical(parameter, "rho[2]")
      expect_identical(parameter_spec, target_spec)
      expect_identical(
        attr(raw_posterior, "prior_density", exact = TRUE),
        semantic_prior_density
      )
      posterior
    },
    .hypothesis_brma_append_iwmde_warnings = function(table, ...) table,
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    marginal_posterior = function(...) seq(0.01, 0.99, length.out = 100L),
    hypothesis_BF = function(..., density_method) {
      used_density_method <<- density_method
      structure(data.frame(BF = 1), class = c(
        "BayesTools_hypothesis_BF",
        "data.frame"
      ))
    },
    .package = "BayesTools"
  )

  out <- .hypothesis_brma_random(
    object                    = list(),
    parameter                 = "theta",
    hypothesis                = BayesTools::hypothesis_parse(c(
      "theta != 0.701406683025 vs theta = 0.701406683025",
      "theta != 1 vs theta = 1"
    )),
    standardized_coefficients = FALSE,
    conditional               = FALSE,
    logBF                     = FALSE,
    BF01                      = FALSE,
    seed                      = 1,
    density_method            = "qCMDE",
    density_control           = list(
      n_points             = 20L,
      samples              = 50L,
      target_relative_mcse = 0.05,
      normalization_points = 50L,
      normalization_prob   = 0.999
    ),
    n_samples = 100L,
    columns   = "default"
  )

  expect_s3_class(out, "BayesTools_hypothesis_BF")
  expect_identical(attachment_calls, 1L)
  expect_identical(attached_values, c(0.701406683025, 1))
  expect_identical(used_density_method, "precomputed")
  expect_identical(reused_selected, selected)
  expect_identical(reused_prior, selected)
})


test_that("semantic random point hypotheses reject singular display boundaries", {

  selected <- list(
    entry = list(term = "tau2_common"),
    spec = list(
      quantity         = "var_common",
      source_parameter = "tau",
      label            = "tau2_common"
    ),
    samples      = matrix(seq(0.01, 0.99, length.out = 100L), ncol = 1L),
    prior        = NULL,
    source_prior = BayesTools::prior("uniform", list(a = 0, b = 1))
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .brma_random_parameter_density_target = function(...) list(
      parameter         = "tau",
      parameter_spec    = list(type = "primitive"),
      display_transform = list(type = "square")
    ),
    .package = "RoBMA"
  )

  expect_error(
    .hypothesis_brma_random(
      object                    = list(),
      parameter                 = "tau2_common",
      hypothesis                = BayesTools::hypothesis_parse(
        "tau2_common = 0"
      ),
      standardized_coefficients = FALSE,
      conditional               = FALSE,
      logBF                     = FALSE,
      BF01                      = FALSE,
      seed                      = 1,
      density_method            = "qCMDE",
      density_control           = list(),
      n_samples                 = 100L,
      columns                   = "default"
    ),
    "public transformation is singular at that support boundary",
    fixed = TRUE
  )
})


test_that("allocation-derived SD zero tests name the omission coordinate", {

  allocation <- list(
    target          = "block",
    scale           = "total_variance",
    index           = 1L,
    component_names = c("study", "observation")
  )
  term <- list(
    block_name  = "study",
    group_label = "study_id",
    sd_binding  = list(
      true_allocation = TRUE,
      allocations     = list(allocation)
    )
  )
  formula_design <- list(mu = list(
    parameter      = "mu",
    random_effects = list(term)
  ))
  object <- list(fit = structure(list(), formula_design = formula_design))
  selected <- list(
    entry = list(term = "study: tau"),
    spec = list(
      label              = "study: tau",
      quantity           = "sd",
      allocation_derived = TRUE,
      formula_parameter  = "mu",
      block              = "study",
      grouping           = "study_id",
      random_component   = "intercept"
    ),
    samples      = matrix(seq(.1, 1, length.out = 20L), ncol = 1L),
    prior        = NULL,
    source_prior = NULL
  )
  testthat::local_mocked_bindings(
    .brma_random_parameter_select = function(...) selected,
    .package = "RoBMA"
  )

  expect_error(
    .hypothesis_brma_random(
      object                    = object,
      parameter                 = "theta",
      hypothesis                = BayesTools::hypothesis_parse("theta = 0"),
      standardized_coefficients = FALSE,
      conditional               = FALSE,
      logBF                     = FALSE,
      BF01                      = FALSE,
      seed                      = 1,
      density_method            = "qCMDE",
      density_control           = list(),
      n_samples                 = 100L,
      columns                   = "default"
    ),
    paste0(
      "Point-null Bayes factors are unavailable for allocation-derived ",
      "random-effect quantity 'study: tau' at 0 because zero is a ",
      "nonregular product boundary of the common scale and allocation ",
      "weight. Test 'tau2_prop(study) = 0' to compare omission of this ",
      "component."
    ),
    fixed = TRUE
  )

  mean_design <- formula_design
  mean_design[["mu"]][["random_effects"]][[1L]][["sd_binding"]][[
    "allocations"
  ]][[1L]][["scale"]] <- "mean_variance"
  mean_object <- list(fit = structure(list(), formula_design = mean_design))
  expect_error(
    .hypothesis_brma_random(
      object                    = mean_object,
      parameter                 = "theta",
      hypothesis                = BayesTools::hypothesis_parse("theta = 0"),
      standardized_coefficients = FALSE,
      conditional               = FALSE,
      logBF                     = FALSE,
      BF01                      = FALSE,
      seed                      = 1,
      density_method            = "qCMDE",
      density_control           = list(),
      n_samples                 = 100L,
      columns                   = "default"
    ),
    paste0(
      "Point-null Bayes factors are unavailable for allocation-derived ",
      "random-effect quantity 'study: tau' at 0 because zero is a ",
      "nonregular product boundary of the common scale and allocation ",
      "weight. Test 'tau2_mult(study) = 0' to compare omission of this ",
      "component."
    ),
    fixed = TRUE
  )
})


test_that("batched ordinates retain value-specific diagnostics", {

  values <- c(0, 1)
  prior_ordinates <- .iwmde_prior_ordinate_classifications(
    prior_density = BayesTools::prior("beta", list(alpha = 1, beta = 1)),
    values        = values
  )
  density <- list(
    x                                = values,
    evaluation_x                     = values,
    y                                = c(0.25, 0.5),
    pilot_y                          = c(0.24, 0.48),
    validation_y                     = c(0.26, 0.52),
    ordinate_relative_change         = c(0.01, 0.02),
    ordinate_log_change              = c(0.01, 0.02),
    pilot_ordinate_relative_change   = c(0.01, 0.02),
    pilot_ordinate_log_change        = c(0.01, 0.02),
    mcse                             = c(0.01, 0.02),
    relative_mcse                    = c(0.04, 0.08),
    active_branch_mcse               = c(0.01, 0.02),
    active_branch_relative_mcse      = c(0.04, 0.08),
    active_mass_component_mcse       = c(0, 0),
    sampling_mcse                    = c(0.03, 0.04),
    sampling_relative_mcse           = c(0.12, 0.16),
    finite_terms                     = c(500L, 450L),
    ess                              = c(300, 40),
    max_weight_share                 = c(0.02, 0.25),
    max_log_ratio                    = c(1, 2),
    estimator                        = "q_grid_cmde"
  )
  diagnostics <- .mock_iwmde_good_diagnostics(
    estimator  = "q_grid_cmde",
    rows       = 500L,
    value      = 0,
    include_bf = TRUE
  )
  diagnostics[["prior_ordinates"]]   <- prior_ordinates
  diagnostics[["ordinate_warnings"]] <- character()
  diagnostic <- list(
    status       = "ok",
    parameter    = "rho[2]",
    target_key   = "simplex_pair|rho|2",
    point_masses = .mock_iwmde_empty_point_masses(),
    iwmde        = density,
    diagnostics  = diagnostics,
    plan = list(
      plan_key           = "shared-plan",
      source_fingerprint = list(source = "test"),
      prior_ordinates    = prior_ordinates,
      target = list(
        target_key = "simplex_pair|rho|2",
        metadata   = list(parameter = "rho[2]")
      )
    )
  )

  attributes <- .iwmde_posterior_ordinate_attributes(
    diagnostic      = diagnostic,
    density_method  = "qCMDE",
    density_control = list(samples = 500L)
  )
  entries <- .iwmde_posterior_ordinate_entries(attributes[["accepted"]])

  expect_equal(vapply(entries, `[[`, numeric(1), "value"), values)
  expect_equal(vapply(entries, `[[`, numeric(1), "ordinate"), density[["y"]])
  expect_equal(vapply(entries, function(entry) {
    entry[["diagnostics"]][["relative_mcse"]]
  }, numeric(1)), density[["relative_mcse"]])
  expect_equal(vapply(entries, function(entry) {
    entry[["diagnostics"]][["sampling_relative_mcse"]]
  }, numeric(1)), density[["sampling_relative_mcse"]])
  expect_equal(vapply(entries, function(entry) {
    entry[["diagnostics"]][["ess"]]
  }, numeric(1)), density[["ess"]])
  expect_true(all(vapply(entries, function(entry) {
    length(entry[["iwmde_provenance"]][["prior_ordinates"]]) == 1L
  }, logical(1))))
  expect_identical(attributes[["rejected"]], NULL)
})

test_that("separate random targets share one hypothesis result", {

  hypothesis <- BayesTools::hypothesis_parse(c(
    "`tau2_prop(study)` != 0 vs `tau2_prop(study)` = 0",
    "tau_total = 0",
    "`tau2_prop(study)` != 1 vs `tau2_prop(study)` = 1"
  ))
  selections <- list(
    list(component = "random", parameter = "rho[2]"),
    list(component = "random", parameter = "tau"),
    list(component = "random", parameter = "rho[2]")
  )
  calls <- list()
  testthat::local_mocked_bindings(
    hypothesis.brma = function(object, hypothesis, ...) {

      n <- length(hypothesis)
      is_allocation <- grepl("tau2_prop", hypothesis[[1L]], fixed = TRUE)
      out <- BayesTools::hypothesis_BF(
        posterior  = if (is_allocation) {
          seq(-1, 2, length.out = 101L)
        } else {
          seq(-2, 1, length.out = 101L)
        },
        prior      = seq(-2, 2, length.out = 101L),
        hypothesis = paste("theta >", seq_len(n) / 10),
        parameter  = "theta"
      )
      attr(out, "warnings") <- stats::setNames(
        paste("warning", seq_len(n)),
        rownames(out)
      )
      calls[[length(calls) + 1L]] <<- list(
        hypothesis = hypothesis,
        result     = out
      )
      out
    },
    .package = "RoBMA"
  )

  out <- .hypothesis_brma_multiple_parameters(
    object                    = list(),
    hypothesis                = hypothesis,
    selections                = selections,
    standardized_coefficients = FALSE,
    conditional               = FALSE,
    conditional_omitted       = TRUE,
    logBF                     = FALSE,
    BF01                      = FALSE,
    seed                      = 1,
    density_method            = "qCMDE",
    density_control           = list(),
    n_samples                 = 100L,
    columns                   = "default"
  )

  expect_length(calls, 2L)
  expect_identical(lengths(lapply(calls, `[[`, "hypothesis")), c(2L, 1L))
  allocation_BF <- attr(calls[[1L]][["result"]], "raw_BF", exact = TRUE)
  total_BF      <- attr(calls[[2L]][["result"]], "raw_BF", exact = TRUE)
  expect_identical(
    attr(out, "raw_BF"),
    c(allocation_BF[[1L]], total_BF, allocation_BF[[2L]])
  )
  expect_identical(
    names(attr(out, "warnings")),
    c("theta", "theta2", "theta1")
  )
  expect_identical(attr(out, "hypothesis_ast"), hypothesis)
  expect_s3_class(out, "BayesTools_hypothesis_BF")
})
