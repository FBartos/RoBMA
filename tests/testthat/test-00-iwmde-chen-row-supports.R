context("(00) IWMDE Chen proposals for row-varying bounded supports")

# The linear route's support is the interval along which the other coordinates
# keep to their own supports, so it differs from draw to draw. The proposal of
# a row has to integrate to one on that row's support: the rows are standardized
# by their own bounds before the proposal family is fitted.


# Two coordinates with uniform(0, 1) priors, and the target a + b. Its support
# on the row (a, b) is [|a - b|, 2 - |a - b|]: moving the target moves both
# coordinates together. The rows are a low-discrepancy sequence of the prior.
.chen_linear_fixture <- function(n) {

  i <- seq_len(n)
  samples <- cbind(
    a = (.5 + i * 0.7548776662466927) %% 1,
    b = (.5 + i * 0.5698402909980532) %% 1
  )
  context <- .iwmde_context_ensure_caches(list(
    posterior_samples = samples,
    flat_prior_list   = list(a = BayesTools::prior("uniform", list(0, 1)),
                             b = BayesTools::prior("uniform", list(0, 1))),
    indicator_names   = character()
  ))

  list(context = context, samples = samples, rows = i,
       spec = list(type = "linear", weights = c(a = 1, b = 1)),
       values = rowSums(samples),
       supports = cbind(abs(samples[, 1L] - samples[, 2L]),
                        2 - abs(samples[, 1L] - samples[, 2L])))
}


.chen_linear_weight <- function(fixture, active_rows = fixture$rows,
                                active_values = fixture$values) {

  .iwmde_chen_log_weight(
    context        = fixture$context,
    parameter      = "a+b",
    parameter_spec = fixture$spec,
    active_rows    = active_rows,
    active_values  = active_values,
    weight_rows    = fixture$rows,
    weight_values  = fixture$values,
    support        = c(0, 2)
  )
}


test_that("a bounded linear target has a finite proposal at every row", {

  fixture <- .chen_linear_fixture(40L)
  supports <- .iwmde_parameter_row_supports(fixture$context, "a+b", fixture$rows,
                                            fixture$spec)
  # The supports are the closed forms, so every row has a support of its own.
  expect_lt(max(abs(supports - fixture$supports)), 1e-12)
  expect_gt(length(unique(.iwmde_chen_support_keys(supports))), 30L)

  weight <- .chen_linear_weight(fixture)
  expect_true(all(is.finite(weight[["log_weight"]])))
  expect_length(weight[["partitions"]], 1L)
  expect_identical(weight[["partitions"]][[1L]][["key"]], "both")
  expect_equal(weight[["partitions"]][[1L]][["support"]], c(0, 1))
  expect_identical(weight[["partitions"]][[1L]][["n_eval_rows"]], 40L)
  expect_identical(weight[["partitions"]][[1L]][["n_finite_terms"]], 40L)
})


test_that("the proposal of every row integrates to one on that row's support", {

  fixture <- .chen_linear_fixture(40L)
  for (row in c(1L, 7L, 19L, 33L, 40L)) {
    proposal <- function(x) {
      exp(.chen_linear_weight(fixture, active_rows = rep(row, length(x)),
                              active_values = x)[["log_weight"]])
    }
    integral <- stats::integrate(proposal, fixture$supports[row, 1L],
                                 fixture$supports[row, 2L], rel.tol = 1e-8,
                                 subdivisions = 1000L)[["value"]]
    expect_equal(integral, 1, tolerance = 1e-6, info = paste("row", row))
  }
})


test_that("the importance-weighted estimate of a + b is the triangular density", {

  # With the prior alone, q is uniform on each row's support, so the estimate
  # at psi is the mean of the proposal over the rows whose support contains it.
  fixture <- .chen_linear_fixture(2000L)
  weight  <- exp(.chen_linear_weight(fixture)[["log_weight"]])
  estimate <- function(psi) {
    mean(weight * (fixture$supports[, 1L] <= psi & psi <= fixture$supports[, 2L]))
  }
  psi <- c(.25, .5, 1, 1.5, 1.75)
  triangular <- ifelse(psi <= 1, psi, 2 - psi)
  expect_lt(max(abs(vapply(psi, estimate, numeric(1)) - triangular)), 1e-2)
})


test_that("one-sided and unbounded rows are standardized by the bounds they have", {

  values   <- c(.4, 1.5, 3, -2, .7, 1)
  supports <- rbind(c(0, 1), c(1, Inf), c(-Inf, 5), c(-Inf, Inf), c(2, 1), c(NA, 2))
  chart    <- .iwmde_chen_standardize(values, supports)

  expect_identical(chart[["key"]], c("both", "lower", "upper", "none", NA, NA))
  expect_equal(chart[["value"]], c(.4, .5, 2, -2, NA, NA))
  expect_equal(chart[["log_jacobian"]], c(0, 0, 0, 0, 0, 0))
  wide <- .iwmde_chen_standardize(.5, rbind(c(-1, 3)))
  expect_equal(wide[["value"]], .375)
  expect_equal(wide[["log_jacobian"]], -log(4))
  expect_identical(.iwmde_chen_standard_support("both"), c(0, 1))
  expect_identical(.iwmde_chen_standard_support("lower"), c(0, Inf))
  expect_identical(.iwmde_chen_standard_support("upper"), c(0, Inf))
  expect_identical(.iwmde_chen_standard_support("none"), c(-Inf, Inf))
})


test_that("a row-varying target is partitioned by the pattern of finite bounds, a constant one by its bounds", {

  # Six rows: two with both bounds, two with only a lower bound, two with only
  # an upper bound. The families are those of the pattern, each fitted to
  # standardized values.
  supports <- rbind(c(0, 2), c(.5, 3), c(1, Inf), c(2, Inf), c(-Inf, 4), c(-Inf, 5))
  values   <- c(1, 2, 2.5, 4, 1, 2)
  rownames(supports) <- as.character(seq_len(6L))
  calls <- list()
  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {
      context[["supports"]][as.character(rows), , drop = FALSE]
    },
    .iwmde_chen_logit_conditional_normal_log_weight = function(
        context, parameter, parameter_spec, active_rows, active_values,
        weight_rows, weight_values, support) {
      calls[[length(calls) + 1L]] <<- list(values = active_values, support = support)
      list(log_weight = rep(0, length(active_values)), method = "bounded")
    },
    .iwmde_chen_compact_conditional_normal_log_weight = function(
        context, parameter, parameter_spec, active_rows, active_values,
        weight_rows, weight_values, support) {
      calls[[length(calls) + 1L]] <<- list(values = active_values, support = support)
      list(log_weight = rep(0, length(active_values)), method = "half_bounded")
    },
    .package = "RoBMA"
  )
  run <- function(type) {
    calls <<- list()
    weight <- .iwmde_chen_log_weight(
      context        = list(supports = supports),
      parameter      = "theta",
      parameter_spec = list(type = type, weights = c(theta = 1)),
      active_rows    = 1:6,
      active_values  = values,
      weight_rows    = 1:6,
      weight_values  = values,
      support        = c(-Inf, Inf)
    )
    list(weight = weight, calls = calls)
  }

  linear <- run("linear")
  expect_equal(vapply(linear$weight$partitions, `[[`, character(1), "key"),
               c("both", "lower", "upper"))
  # Both bounds: the unit interval, with the Jacobian -log(upper - lower);
  # one bound: the distance from it.
  expect_equal(linear$calls[[1L]]$values, c(1 / 2, (2 - .5) / 2.5))
  expect_equal(linear$calls[[1L]]$support, c(0, 1))
  expect_equal(linear$calls[[2L]]$values, c(1.5, 2))
  expect_equal(linear$calls[[2L]]$support, c(0, Inf))
  expect_equal(linear$calls[[3L]]$values, c(3, 3))
  expect_equal(linear$calls[[3L]]$support, c(0, Inf))
  expect_equal(linear$weight$log_weight[1:2], -log(c(2, 2.5)))
  expect_equal(linear$weight$log_weight[3:6], rep(0, 4L))

  # A primitive target keeps one partition for each exact pair of bounds.
  primitive <- run("primitive")
  expect_length(primitive$weight$partitions, 6L)
  expect_equal(primitive$calls[[1L]]$support, c(0, 2))
  expect_equal(primitive$weight$log_weight, rep(0, 6L))
})


test_that("a linear target with one common support keeps the weights of the exact-support fit", {

  # Every row has the support (0, 2): the standardized fit on the unit interval
  # with the Jacobian -log(2) is the beta fit of the raw values on (0, 2), which
  # is what the exact-bounds partition gave (the target has no conditioning
  # columns, so the logit-conditional fit falls back to the beta).
  fixture  <- .chen_linear_fixture(60L)
  supports <- matrix(c(0, 2), nrow = 60L, ncol = 2L, byrow = TRUE)
  rownames(supports) <- as.character(seq_len(60L))
  testthat::local_mocked_bindings(
    .iwmde_parameter_row_supports = function(context, parameter, rows,
                                             parameter_spec) {
      supports[as.character(rows), , drop = FALSE]
    },
    .package = "RoBMA"
  )
  linear <- .chen_linear_weight(fixture)
  exact  <- .iwmde_chen_beta_log_weight(fixture$values, fixture$values, c(0, 2))

  expect_identical(linear[["method"]], "chen_beta")
  expect_equal(linear[["log_weight"]], exact[["log_weight"]], tolerance = 1e-12)
})


# ---------------------------------------------------------------------------- #
# Public route
# ---------------------------------------------------------------------------- #

# A prior-only model with a three-level factor (two treatment-contrast
# coefficients with uniform(-1, 1) priors) and the replay of its fit:
# constructor, compiled formula declarations, and draws.
.chen_factor_object <- function(draws) {

  data <- data.frame(yi = seq(-.1, .2, length.out = 6L), sei = .1,
                     g = factor(rep(c("A", "B", "C"), each = 2L)))
  object <- brma(
    yi = yi, sei = sei, mods = ~ 0 + g, data = data, measure = "GEN",
    prior_mods = list(g = BayesTools::prior_factor("uniform", list(-1, 1),
                                                   contrast = "treatment")),
    prior_effect = BayesTools::prior("point", list(location = 0)),
    prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
    only_priors = TRUE, silent = TRUE
  )
  args <- .create_jags_formula_args(object[["data"]], object[["priors"]])
  compiled <- lapply(names(args[["formula_list"]]), function(parameter) {
    BayesTools::JAGS_formula(
      formula = args[["formula_list"]][[parameter]], parameter = parameter,
      data = args[["formula_data_list"]][[parameter]],
      prior_list = args[["formula_prior_list"]][[parameter]],
      formula_scale = args[["formula_scale_list"]][[parameter]],
      prior_random = args[["formula_random_prior_list"]][[parameter]],
      random_effects_compile =
        args[["formula_random_effects_compile_list"]][[parameter]]
    )
  })
  names(compiled) <- names(args[["formula_list"]])
  object[["formula_design"]] <- lapply(compiled, `[[`, "formula_design")
  map <- compiled[["mu"]][["formula_design"]][["name_map"]]
  g_name <- map[["jags_name"]][map[["kind"]] == "fixed" & map[["term"]] == "g"]
  samples <- draws
  colnames(samples) <- paste0(g_name, "[", 1:2, "]")
  samples <- cbind(samples, mu_intercept = 0, tau = .1)
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "list")
  )
  attr(fit, "formula_design") <- object[["formula_design"]]
  scales <- lapply(compiled, `[[`, "formula_scale")
  if (any(lengths(scales) > 0L)) {
    attr(fit, "formula_scale") <- scales
  }
  priors <- list()
  for (declared in compiled) {
    priors <- c(priors, declared[["prior_list"]])
  }
  attr(fit, "prior_list") <- c(priors,
    .create_fit_priors(object[["data"]], object[["priors"]]))
  object[["fit"]] <- as_bayestools_fit(fit)

  object
}


test_that("IWMDE of a difference of factor coefficients succeeds and agrees with qCMDE and the exact posterior", {

  .iwmde_qcmde_cache_clear()
  withr::defer(.iwmde_qcmde_cache_clear())
  # The posterior of the two coefficients is independent normal (the group
  # means .05 and .17, SD .1) up to a truncation at +-1 of mass 1e-20; the
  # replay draws are a low-discrepancy sample of it, so the importance weights
  # are those of a posterior.
  n  <- 1000L
  i  <- seq_len(n)
  draws <- cbind(
    .05 + .1 * stats::qnorm((.5 + i * 0.7548776662466927) %% 1),
    .17 + .1 * stats::qnorm((.5 + i * 0.5698402909980532) %% 1)
  )
  object <- .chen_factor_object(draws)
  statement <- "g[B] - g[C] = 0.1"

  # The posterior of the difference is normal(-.12, sqrt(.02)); its prior is
  # triangular on (-2, 2), with ordinate (2 - .1) / 4 at .1.
  posterior <- stats::dnorm(.1, -.12, sqrt(.02))
  exact_bf  <- ((2 - .1) / 4) / posterior

  iwmde <- hypothesis(object, statement, density_method = "IWMDE", columns = "all")
  qcmde <- hypothesis(object, statement, density_method = "qCMDE", columns = "all")
  expect_equal(unname(attr(iwmde, "raw_BF")), exact_bf, tolerance = .01)
  expect_equal(unname(attr(qcmde, "raw_BF")), exact_bf, tolerance = .01)
  expect_equal(unname(attr(iwmde, "raw_BF")), unname(attr(qcmde, "raw_BF")),
               tolerance = .01)
  diagnostics <- density_diagnostics(iwmde)
  expect_identical(diagnostics[["estimator"]], "iwmde")
  expect_lt(diagnostics[["relative_mcse"]], .01)
})


# Real public constructor/compiler/fit contract, with explicitly synthetic
# source rows; these tests verify the proposal seam, not fitted BF calibration.
.zero_i <- seq_len(80L)
.zero_object <- .chen_factor_object(cbind(
  .05 + .1 * stats::qnorm((.5 + .zero_i * 0.7548776662466927) %% 1),
  .17 + .1 * stats::qnorm((.5 + .zero_i * 0.5698402909980532) %% 1)))
.zero_context <- .iwmde_context(.zero_object)
.zero_parameter <- colnames(.zero_context[["posterior_samples"]])[[1L]]
.zero_plan <- .iwmde_plan(.zero_context, .zero_parameter, "IWMDE", list(samples = 40L),
  outputs = "ordinate", values = .1, parameter_spec = list(type = "primitive"))
stopifnot(identical(.zero_plan[["status"]], "ok"))
.zero_execution <- .iwmde_plan_row_execution(.zero_context, .zero_plan)
.chen_execute_zero_proposal <- function(proposal) {
  execution <- .zero_execution
  execution$proposal_weight <- proposal
  .iwmde_execute_plan_estimator(.zero_context, .zero_plan, "ordinate", execution)
}

testthat::test_that("certified zero weights retain the exact row population and denominator", {
  proposal <- .zero_execution$proposal_weight
  n <- length(proposal$log_weight)
  mask <- setNames(seq_len(n) %% 5L == 0L, as.character(.zero_execution$active_rows))
  proposal$log_weight[mask] <- -Inf
  proposal$zero_support <- mask
  result <- tryCatch(.chen_execute_zero_proposal(proposal), error = identity)
  testthat::expect_false(inherits(result, "error"))
  if (!inherits(result, "error")) {
    q <- .iwmde_log_q_grid(.zero_context, .zero_parameter, .1,
      .zero_execution$row_states, .zero_plan$replacement)
    baseline <- vapply(.zero_execution$row_states, `[[`, numeric(1), "baseline_log_q")
    expected <- .zero_plan$rows$active_mass * sum(exp(as.numeric(q) - baseline + proposal$log_weight)) / n
    testthat::expect_lte(abs(result$y - expected), 64 * .Machine$double.eps * max(1, expected))
    testthat::expect_identical(result$n_candidate_rows, as.integer(n))
    testthat::expect_identical(result$n_evaluated_rows, as.integer(n))
    testthat::expect_identical(as.integer(result$finite_terms), as.integer(sum(!mask)))
    testthat::expect_equal(result$sampling_fraction, n / length(.zero_execution$sampling_population_rows))
    testthat::expect_null(result$zero_support)
    testthat::expect_false(any(vapply(result$weight_partitions, function(partition) "zero_support" %in% names(partition), logical(1))))
  }
})

testthat::test_that("malformed or contradictory zero certificates fail before density evaluation", {
  proposal <- .zero_execution$proposal_weight
  n <- length(proposal$log_weight)
  valid <- setNames(rep(FALSE, n), as.character(.zero_execution$active_rows))
  malformed <- list(as.integer(valid), valid[-1L], replace(valid, 1L, NA),
    setNames(valid, rev(names(valid))), unname(valid), replace(valid, 1L, TRUE))
  for (mask in malformed) {
    weight <- proposal
    weight$zero_support <- mask
    error <- tryCatch(.chen_execute_zero_proposal(weight), error = identity)
    testthat::expect_s3_class(error, "iwmde_construction_error")
  }
  for (value in c(-Inf, Inf, NA_real_, NaN)) {
    weight <- proposal
    weight$log_weight[1L] <- value
    weight$zero_support <- valid
    testthat::expect_error(.chen_execute_zero_proposal(weight), class = "iwmde_construction_error")
  }
})

testthat::test_that("an all-zero certified proposal remains a diagnostic failure", {
  proposal <- .zero_execution$proposal_weight
  proposal$log_weight[] <- -Inf
  proposal$zero_support <- setNames(rep(TRUE, length(proposal$log_weight)), as.character(.zero_execution$active_rows))
  result <- tryCatch(.chen_execute_zero_proposal(proposal), error = identity)
  testthat::expect_false(inherits(result, "iwmde_construction_error"))
  if (!inherits(result, "error")) {
    testthat::expect_identical(result$n_evaluated_rows, .zero_execution$n_candidate_rows)
    testthat::expect_identical(as.numeric(result$y), 0)
    testthat::expect_false(any(result$structural_zero))
    reason <- .iwmde_posterior_ordinate_bf_failure_reason(
      list(ordinate = result$y[1L], diagnostics = list(estimator = "iwmde")))
    testthat::expect_identical(reason, "posterior ordinate is zero or non-finite")
  }
})
test_that("compiled allocated SD proposals normalize on their own declared intervals", {

  compact <- NULL
  compact_evaluator <- .iwmde_chen_compact_normal_log_weight
  # Observe the genuine callback parameters without changing its result. Its
  # compact endpoints are quadrature breakpoints, not target/prior support.
  testthat::local_mocked_bindings(.iwmde_chen_compact_normal_log_weight = function(
      distance, mean, sd, focal_scale, method) {
    compact <<- list(mean = mean, sd = sd, focal_scale = focal_scale)
    compact_evaluator(distance, mean, sd, focal_scale, method)
  }, .package = "RoBMA")
  u <- seq_len(32L) / 33
  cases <- list(
    uniform = list(prior = BayesTools::prior("uniform", list(0, 1)), values = u, lower = rep(0, 32), upper = rep(1, 32)),
    truncated = list(prior = BayesTools::prior("normal", list(1, .7), truncation = list(.1, 2)),
      values = .1 + 1.9 * u, lower = rep(.1, 32), upper = rep(2, 32)),
    positive = list(prior = BayesTools::prior("gamma", list(2, 2), truncation = list(.4, Inf)),
      values = .4 + u, lower = rep(.4, 32), upper = rep(Inf, 32)),
    mixture = list(prior = BayesTools::prior_mixture(list(
      BayesTools::prior("uniform", list(0, 1)), BayesTools::prior("uniform", list(.2, 2))), is_null = c(FALSE, FALSE)),
      values = c(1:16 / 17, .2 + 1.8 * (1:16 / 17)), indicators = rep(1:2, each = 16),
      lower = rep(c(0, .2), each = 16), upper = rep(c(1, 2), each = 16)),
    full_positive = list(prior = BayesTools::prior("gamma", list(2, 2)), values = u,
      lower = rep(0, 32), upper = rep(Inf, 32))
  )
  for (name in names(cases)) {
    case <- cases[[name]]
    scale <- if (is.null(case$scale)) "total_variance" else case$scale
    fixture <- review_allocated_source_replay("uniform", source_prior = case$prior,
      source_values = case$values, source_indicators = case$indicators, scale = scale)
    object <- fixture$object
    target <- .brma_random_parameter_density_target(object, "(mu) study: tau")
    context <- .iwmde_context(object)
    spec <- .iwmde_parameter_spec(context, target$parameter, target$parameter_spec)
    values <- .iwmde_parameter_values(context, target$parameter, spec)
    multiplier <- sqrt(fixture$samples[, fixture$weights[1L]] * if (scale == "mean_variance") 2 else 1)
    expected <- cbind(case$lower * multiplier, case$upper * multiplier)
    supports <- .iwmde_parameter_row_supports(context, target$parameter, 1:32, spec)
    expect_equal(unname(supports), unname(expected), tolerance = 1e-14, info = name)
    set.seed(1401)
    plan <- .iwmde_plan(context, target$parameter, density_method = "IWMDE",
      density_control = list(samples = 24L), outputs = "ordinate", values = .6, parameter_spec = target$parameter_spec)
    expect_identical(plan$status, "ok", info = name)
    execution <- .iwmde_plan_row_execution(context, plan)
    fit_rows <- plan$rows$continuous_rows
    eval_rows <- execution$active_rows
    expect_length(fit_rows, 32L)
    expect_length(eval_rows, 24L)
    expect_identical(as.numeric(execution$active_values), as.numeric(values[eval_rows]))
    proposal <- function(rows, x) .iwmde_chen_log_weight(context, target$parameter, spec,
      active_rows = rows, active_values = x, weight_rows = fit_rows,
      weight_values = values[fit_rows], support = plan$support$support)
    actual <- proposal(eval_rows, values[eval_rows])
    expect_identical(actual, execution$proposal_weight)
    expect_length(actual$partitions, 1L)
    # All 32 fitting rows in the original bounded fixture, including rows that
    # are not estimator rows; representative rows of the additional laws.
    rows <- if (name == "uniform") 1:32 else c(3L, 11L, 23L)
    for (row in rows) {
      callback <- function(x) exp(proposal(rep(row, length(x)), x)$log_weight)
      compact <- NULL
      proposal(row, values[row])
      captured <- compact
      pieces <- expected[row, ]
      if (!is.null(captured)) {
        log_tail <- stats::pnorm(0, captured$mean, captured$sd, lower.tail = FALSE, log.p = TRUE)
        endpoints <- expected[row, 1L] + captured$focal_scale * stats::qnorm(
          log_tail + log(c(.975, .025)), captured$mean, captured$sd, lower.tail = FALSE, log.p = TRUE)
        expect_true(all(is.finite(endpoints)) && endpoints[1L] > expected[row, 1L] &&
          endpoints[2L] < expected[row, 2L])
        expect_equal(callback(c(mean(c(expected[row, 1L], endpoints[1L])),
          endpoints[2L] + max(1, diff(endpoints)))), c(0, 0))
        pieces <- c(expected[row, 1L], endpoints, expected[row, 2L])
        if (name == "positive" && row == 23L) {
          plain <- stats::integrate(callback, expected[row, 1L], expected[row, 2L],
            rel.tol = 1e-9, abs.tol = 1e-10, subdivisions = 1000L)
          cat("PRESERVED_PLAIN_HALF_LINE_QUADRATURE positive23 value=", format(plain$value, digits = 17),
            " abs.error=", format(plain$abs.error, digits = 17), " message=", plain$message, "\n", sep = "")
        }
      }
      integrals <- lapply(seq_len(length(pieces) - 1L), function(i) stats::integrate(callback,
        pieces[i], pieces[i + 1L], rel.tol = 1e-9, abs.tol = 1e-10, subdivisions = 1000L))
      value <- sum(vapply(integrals, `[[`, numeric(1), "value"))
      error <- sum(vapply(integrals, `[[`, numeric(1), "abs.error"))
      expect_true(abs(value - 1) < 1e-7 + error, info = paste(name, row, value, error))
    }
  }
})

test_that("missing allocated bounds and all-invalid row populations fail without widening", {

  fixture <- review_allocated_source_replay("uniform")
  target <- .brma_random_parameter_density_target(fixture$object, "(mu) study: tau")
  context <- .iwmde_context(fixture$object)
  context$flat_prior_list[[target$parameter]] <- NULL
  expect_error(.iwmde_parameter_row_supports(context, target$parameter, 1:2, target$parameter_spec),
    class = "iwmde_construction_error")
  linear <- .chen_linear_fixture(4L)
  linear$context$posterior_samples[,] <- NA_real_
  warnings <- list()
  outcome <- withCallingHandlers(tryCatch(.iwmde_parameter_support(linear$context, "a+b", 1:4,
    linear$spec), error = identity), warning = function(w) warnings[[length(warnings) + 1L]] <<- w)
  expect_s3_class(outcome, "iwmde_construction_error")
  expect_length(warnings, 0L)
  expect_equal(.iwmde_parameter_support(linear$context, "a+b", integer(), linear$spec), c(-Inf, Inf))
  mixed <- .chen_linear_fixture(4L)
  expected <- .iwmde_parameter_support(mixed$context, "a+b", c(1L, 3L, 4L), mixed$spec)
  mixed$context$posterior_samples[2L, ] <- NA_real_
  expect_equal(.iwmde_parameter_support(mixed$context, "a+b", 1:4, mixed$spec), expected)
  # The existing genuine compiler fixture has a gated root and a child split.
  nested <- shared_gate_random_object(root_gate = rep(1, 4))
  selected <- .brma_random_parameter_density_target(nested, "(mu) study: tau(intercept)")
  nested_context <- .iwmde_context_ensure_caches(list(object = nested,
    posterior_samples = as.matrix(nested$fit[[1L]]),
    flat_prior_list = attr(nested$fit, "prior_list"), indicator_names = character()))
  expect_equal(.iwmde_parameter_row_supports(nested_context, selected$parameter,
    c(4L, 1L), selected$parameter_spec), matrix(c(0, Inf), 2L, 2L, byrow = TRUE))
})

test_that("admitted SD-component allocations retain their nonunit source multiplier", {

  weight <- ((seq_len(32) * 11) %% 37 + .5) / 37
  shares <- cbind(weight, 1 - weight)
  source <- .2 + 1.8 * seq_len(32) / 33
  object <- sd_component_allocation_object(list(source = source, weight = shares, eta = 2 * shares),
    source_prior = BayesTools::prior("uniform", list(.2, 2)))
  target <- .brma_random_parameter_density_target(object, "tau(intercept)")
  samples <- as.matrix(object$fit[[1L]])
  context <- .iwmde_context_ensure_caches(list(object = object, posterior_samples = samples,
    flat_prior_list = attr(object$fit, "prior_list"), indicator_names = character()))
  spec <- .iwmde_parameter_spec(context, target$parameter, target$parameter_spec)
  multiplier <- sqrt(2 * weight)
  expected <- cbind(.2 * multiplier, 2 * multiplier)
  expect_equal(.iwmde_parameter_row_supports(context, target$parameter, 1:32, spec), expected, tolerance = 1e-14)
  values <- .iwmde_parameter_values(context, target$parameter, spec)
  expect_lt(max(abs(values - source * multiplier)), 1e-14)
  for (row in c(3L, 11L, 23L)) {
    proposal <- function(x) exp(.iwmde_chen_log_weight(context, target$parameter, spec,
      rep(row, length(x)), x, 1:32, values, range(expected))$log_weight)
    integral <- stats::integrate(proposal, expected[row, 1L], expected[row, 2L], rel.tol = 1e-9, abs.tol = 1e-10)
    expect_lt(abs(integral$value - 1), 1e-7 + integral$abs.error)
  }
  absent_node <- spec
  absent_node$node <- NULL
  expect_error(.iwmde_parameter_row_supports(context, target$parameter, 1:2, absent_node),
    "multipliers need the deterministic nodes")
})
