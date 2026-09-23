.retained_location_test_states <- function(prior, current, coefficients, sd) {

  change <- .iwmde_retained_location_statistics(coefficients, sd)
  baseline <- .iwmde_focal_log_prior_values(prior, current, "target")
  lapply(seq_along(current), function(row) {
    .iwmde_new_row_state(list(row_index = row, active_key = "all", current = current[row],
      baseline_log_q = baseline[row], baseline_log_lik = 0,
      baseline_log_prior = baseline[row], baseline_focal_log_prior = baseline[row],
      focal_prior = prior, use_focal_prior_delta = TRUE,
      gaussian_change = list(linear = change$linear[row], quadratic = change$quadratic[row]),
      conditioning_transform = list(type = "retained_gaussian_intercept", parameter = "target")))
  })
}

test_that("retained-location q kernels give the independent Gaussian conditional law", {

  current <- c(-.3, .7)
  coefficients <- rbind(c(.2, -.4, .8), c(-.1, .5, -.7))
  sd <- c(.6, .9)
  prior <- BayesTools::prior("normal", list(mean = .15, sd = 1.2))
  states <- .retained_location_test_states(prior, current, coefficients, sd)
  values <- c(-1, -.2, .1, .8, 1.7)
  observed <- .iwmde_log_q_grid_retained_location(list(), "target", values, states,
    list(type = "scalar"))
  m <- coefficients + current
  J <- ncol(m)
  precision <- 1 / 1.2^2 + J / sd^2
  conditional_mean <- (.15 / 1.2^2 + rowSums(m) / sd^2) / precision
  conditional_sd <- sqrt(1 / precision)
  for (row in seq_along(current)) {
    expected <- vapply(values, function(value) {
      stats::dnorm(value, .15, 1.2, log = TRUE) +
        sum(stats::dnorm(m[row, ], value, sd[row], log = TRUE)) -
        sum(stats::dnorm(coefficients[row, ], 0, sd[row], log = TRUE))
    }, numeric(1L))
    expect_equal(observed[, row], expected, tolerance = 1e-13)
    normalizer <- stats::integrate(function(value) {
      exp(.iwmde_log_q_grid_retained_location(list(), "target", value, states[row],
        list(type = "scalar"))[, 1L])
    }, -Inf, Inf, rel.tol = 1e-10)$value
    expect_equal(exp(observed[, row]) / normalizer,
      stats::dnorm(values, conditional_mean[row], conditional_sd[row]), tolerance = 1e-9)
  }
})

test_that("retained-location conditioning preserves nonnormal priors and bounded support", {

  current <- .2
  coefficients <- matrix(c(.4, -.3, .1), 1L)
  sd <- .7
  prior <- BayesTools::prior("t", list(.3, .6, 3),
    truncation = list(lower = -1, upper = 2))
  states <- .retained_location_test_states(prior, current, coefficients, sd)
  values <- c(-2, -.8, -.1, .4, 1.5, 3)
  observed <- as.numeric(.iwmde_log_q_grid_retained_location(list(), "target", values,
    states, list(type = "scalar")))
  m <- as.numeric(coefficients) + current
  prior_mass <- stats::pt((2 - .3) / .6, df = 3) - stats::pt((-1 - .3) / .6, df = 3)
  reference <- function(value) {
    vapply(value, function(x) {
      if (x < -1 || x > 2) return(0)
      stats::dt((x - .3) / .6, df = 3) / (.6 * prior_mass) *
        exp(sum(stats::dnorm(m, x, sd, log = TRUE)))
    }, numeric(1L))
  }
  normalizer <- stats::integrate(reference, -1, 2, rel.tol = 1e-11)$value
  kernel_normalizer <- stats::integrate(function(value) {
    exp(.iwmde_log_q_grid_retained_location(list(), "target", value, states,
      list(type = "scalar"))[, 1L])
  }, -1, 2, rel.tol = 1e-11)$value
  expect_equal(exp(observed) / kernel_normalizer, reference(values) / normalizer,
    tolerance = 1e-10)
  expect_identical(observed[c(1L, 6L)], c(-Inf, -Inf))
})

test_that("row normalizers resolve narrow Gaussian kernels missed by a shared grid", {

  # The full-V regression had conditional SD .0009108 and a spuriously
  # converged shared-grid normalizer. Include that width between grid knots.
  prior <- BayesTools::prior("normal", list(0, sqrt(.5)))
  current <- c(.031, .1)
  conditional_mean <- c(.03284134, .17)
  conditional_sd <- c(.0009108, .075)
  precision <- 1 / conditional_sd^2 - 2
  linear <- conditional_mean / conditional_sd^2 - precision * current
  states <- .retained_location_test_states(prior, current,
    matrix(linear / precision, 2L, 1L), 1 / sqrt(precision))
  expected_log_mass <- log(conditional_sd / sqrt(.5)) - current^2 +
    .5 * ((linear - 2 * current) * conditional_sd)^2
  results <- lapply(states, .iwmde_retained_location_normalizer)
  expect_equal(vapply(results, `[[`, numeric(1L), "log_normalizer"),
    expected_log_mass, tolerance = 1e-11)
  expect_equal(vapply(results, `[[`, numeric(1L), "posterior_mean"),
    conditional_mean, tolerance = 1e-14)
  expect_equal(vapply(results, `[[`, numeric(1L), "posterior_sd"),
    conditional_sd, tolerance = 1e-14)

  transform <- .iwmde_parameter_transform(c(-Inf, Inf))
  grid <- .iwmde_qcmde_grid_from_z(seq(-.5, .5, length.out = 30L), transform)
  plan <- .iwmde_qcmde_normalizer_plan(grid, transform)
  values <- c(conditional_mean, .2)
  evaluation <- .iwmde_qcmde_evaluate_grid_sequence(list(), "target", values,
    plan, states, list(type = "scalar"), 1:2, active_mass = 1, denominator = 2L)
  for (normalizer in evaluation$log_normalizer_sequence) {
    expect_equal(normalizer, expected_log_mass, tolerance = 1e-11)
  }
  expect_identical(evaluation$conditional_normalization$methods,
    rep("normal_product", 2L))
  density <- .iwmde_qcmde_density_from_normalizer(evaluation$log_q_display,
    evaluation$log_normalizer_sequence[[1L]], active_mass = 1, denominator = 2L)
  expected <- rowMeans(vapply(seq_along(current), function(row) {
    stats::dnorm(values, conditional_mean[row], conditional_sd[row])
  }, numeric(length(values))))
  expect_equal(density, expected, tolerance = 1e-10)
  shared_log_mass <- .iwmde_log_trapz_columns(grid$z,
    .iwmde_log_q_grid_retained_location(list(), "target", grid$x, states,
      list(type = "scalar")))
  expect_gt(abs(shared_log_mass[1L] - expected_log_mass[1L]), 1)
})

test_that("conditional normalization preserves truncation and scaled t-prior mass", {

  current <- .032
  kernel_center <- .03284
  kernel_sd <- .0009108
  coefficient <- matrix(kernel_center - current, 1L)
  bounds <- c(.0325, .034)
  normal <- BayesTools::prior("normal", list(.1, .7),
    truncation = list(lower = bounds[1L], upper = bounds[2L]))
  states <- .retained_location_test_states(normal, current = .033,
    coefficients = matrix(kernel_center - .033, 1L), sd = kernel_sd)
  result <- .iwmde_retained_location_normalizer(states[[1L]])
  posterior_precision <- 1 / .7^2 + 1 / kernel_sd^2
  posterior_sd <- 1 / sqrt(posterior_precision)
  posterior_mean <- (.1 / .7^2 + kernel_center / kernel_sd^2) / posterior_precision
  prior_mass <- diff(stats::pnorm(bounds, .1, .7))
  posterior_mass <- diff(stats::pnorm(bounds, posterior_mean, posterior_sd))
  values <- c(bounds[1L] - .001, bounds, mean(bounds), bounds[2L] + .001)
  observed <- exp(.iwmde_log_q_grid_retained_location(list(), "target", values,
    states, list(type = "scalar"))[, 1L] - result$log_normalizer)
  expected <- stats::dnorm(values, posterior_mean, posterior_sd) / posterior_mass
  expected[values < bounds[1L] | values > bounds[2L]] <- 0
  expect_equal(observed, expected, tolerance = 1e-9)
  expect_true(is.finite(prior_mass) && prior_mass > 0)

  for (support in list(c(-Inf, Inf), c(-1, 2))) {
    prior <- BayesTools::prior("t", list(.3, .6, 3),
      truncation = list(lower = support[1L], upper = support[2L]))
    states <- .retained_location_test_states(prior, current, coefficient, kernel_sd)
    result <- .iwmde_retained_location_normalizer(states[[1L]])
    prior_mass <- diff(stats::pt((support - .3) / .6, 3))
    # Independent integral of the known t density against a standard normal;
    # no package q-grid, prior evaluator, or normalizer is used in this oracle.
    reference <- stats::integrate(function(z) {
      value <- kernel_center + kernel_sd * z
      stats::dt((value - .3) / .6, 3) / (.6 * prior_mass) * stats::dnorm(z) *
        (value >= support[1L] & value <= support[2L])
    }, -12, 12, rel.tol = 1e-11)$value
    expected <- .5 * ((kernel_center - current) / kernel_sd)^2 +
      log(kernel_sd * sqrt(2 * pi) * reference)
    expect_equal(result$log_normalizer, expected, tolerance = 1e-7)
    expect_lte(abs(expm1(expected - result$log_normalizer)), result$relative_error + 1e-10)
    expect_identical(result$method, "scaled_gaussian_integral")
  }
  invalid <- states[[1L]]
  invalid$gaussian_change$quadratic <- 0
  expect_error(.iwmde_retained_location_normalizer(invalid),
    "Retained-location normalizer inputs are invalid.", fixed = TRUE)
})

test_that("scaled normalization retains finite tail mass beyond ordinary probability range", {

  prior <- BayesTools::prior("uniform", list(10, 11))
  state <- .retained_location_test_states(prior, current = 10.1,
    coefficients = matrix(-10.1, 1L), sd = .01)[[1L]]
  result <- .iwmde_retained_location_normalizer(state)
  log_tails <- stats::pnorm(c(10, 11), sd = .01, lower.tail = FALSE, log.p = TRUE)
  log_mass <- log_tails[1L] + log1p(-exp(log_tails[2L] - log_tails[1L]))
  expected <- .5 * (10.1 / .01)^2 + log(.01 * sqrt(2 * pi)) + log_mass
  expect_true(all(is.finite(c(result$log_normalizer, result$relative_error))))
  expect_lt(abs(result$log_normalizer - expected), 1e-6)
  expect_lte(result$relative_error, .iwmde_qcmde_refinement_target())
})

test_that("zero retained scale falls back and invalid numeric inputs remain visible", {

  coefficients <- rbind(c(.2, -.1), c(.4, .3))
  expect_null(.iwmde_retained_location_statistics(coefficients, c(.5, 0)))
  expect_error(.iwmde_retained_location_statistics(coefficients, c(.5, NaN)),
    "Retained-location Gaussian conditional inputs are invalid.", fixed = TRUE)
  coefficients[2L, 1L] <- Inf
  expect_error(.iwmde_retained_location_statistics(coefficients, c(.5, .8)),
    "Retained-location Gaussian conditional inputs are invalid.", fixed = TRUE)
})

test_that("the fixed information gate retains ordinary rows and records its conditioning policy", {

  V <- diag(c(.04, .09, .16))
  sampling_plan <- structure(list(schema_version = 5L, sampling = list(
    representation = "diagonal_factor", diagonal = diag(V), loading = matrix(0, 3L, 0L))),
    class = c("RoBMA_selection_execution_plan", "list"))
  data <- structure(list(outcome = data.frame(yi = c(.1, -.2, .3), sei = sqrt(diag(V)))),
    selection_execution_plan = sampling_plan)
  prior <- BayesTools::prior("normal", list(0, 1))
  samples <- cbind(target = c(-.3, .2, .5, .7), sd = c(.5, .08, 0, .12),
    b1 = c(.2, -.1, .3, .4), b2 = c(-.1, .6, .2, -.2))
  context <- .iwmde_context_ensure_caches(list(data = data,
    object = list(fit = NULL, data = data, priors = list()), priors = list(),
    posterior_samples = samples, flat_prior_list = list(target = prior)))
  plan <- list(method = "q_grid_cmde", target = list(parameter = "target"),
    execution_spec = list(type = "primitive"))
  metadata <- list(parameter = "target", block = "retained", term = list(),
    group_map = c(1L, 1L, 2L), group_rows = c(1L, 3L), n_groups = 2L, prior = prior)
  ordinary_calls <- list()
  testthat::local_mocked_bindings(
    .iwmde_retained_location_plan = function(...) metadata,
    .random_effect_term_sd_samples = function(term, posterior_samples, K) {
      matrix(posterior_samples[, "sd"], nrow(posterior_samples), K)
    },
    .evaluate.brma.random_effects = function(..., posterior_samples) {
      posterior_samples[, c("b1", "b1", "b2"), drop = FALSE]
    },
    .iwmde_row_states_grouped_marginal = function(...) NULL,
    .iwmde_row_states = function(context, rows, ...) {
      ordinary_calls[[length(ordinary_calls) + 1L]] <<- rows
      lapply(rows, function(row) .iwmde_new_row_state(list(row_index = row,
        baseline_log_q = 0, active_key = "all", likelihood_mode = "conditional")))
    },
    .package = "RoBMA"
  )
  states <- .iwmde_retained_location_row_states(context, plan, 1:4)
  transformed <- function(states) vapply(states, function(state) {
    !is.null(state$conditioning_transform)
  }, logical(1L))
  expect_identical(vapply(states, `[[`, integer(1L), "row_index"), 1:4)
  expect_identical(transformed(states), c(TRUE, FALSE, FALSE, TRUE))
  expect_identical(ordinary_calls, list(2:3))
  policy <- attr(states, "conditioning_policy")
  expect_identical(policy$information_multiplier, 4L)
  expect_equal(policy$information_cap, 4 * sum(1 / diag(V)), tolerance = 1e-14)
  expect_identical(policy$transformed_rows, 2L)
  expect_identical(policy$ordinary_rows, 2L)

  shifted <- context
  shifted$posterior_samples[, c("target", "b1", "b2")] <-
    shifted$posterior_samples[, c("target", "b1", "b2")] + 10
  shifted_states <- .iwmde_retained_location_row_states(shifted, plan, 1:4)
  expect_identical(attr(shifted_states, "conditioning_policy"), policy)
  smaller <- .iwmde_retained_location_row_states(context, plan, c(1L, 4L))
  expect_identical(transformed(smaller), c(TRUE, FALSE))
  expect_identical(attr(smaller, "conditioning_policy")$information_multiplier, 2L)

  plan$method <- "iwmde"
  expect_null(.iwmde_retained_location_row_states(context, plan, 1:4))
  plan$method <- "q_grid_cmde"
  invalid <- context
  invalid$posterior_samples[1L, "b1"] <- NaN
  expect_error(.iwmde_retained_location_row_states(invalid, plan, 1:4),
    "Retained-location sampled coefficients are invalid.", fixed = TRUE)
})

test_that("constructor metadata routes iid allocations and excludes known-R conditioning", {

  dat <- data.frame(yi = c(.1, -.2, .3, .05), vi = c(.04, .09, .16, .25),
    study = factor(c("a", "a", "b", "b")), esid = factor(1:4))
  S <- 24L
  bias <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .5)), model = selection_model(group = "study"))
  for (sampling_kind in c("full", "diagonal", "known_R")) {
    V <- diag(dat$vi)
    if (sampling_kind != "diagonal") V[cbind(c(1L, 2L, 3L, 4L), c(2L, 1L, 4L, 3L))] <- .01
    arguments <- list(yi = dat$yi, random = ~1 | study/esid, data = dat,
      measure = "GEN", prior_unit_information_sd = 1, prior_bias = bias,
      prior_effect = BayesTools::prior("t", list(0, .7, 3)),
      only_priors = TRUE, silent = TRUE, effect_direction = "positive")
    if (sampling_kind != "diagonal") arguments$V <- V else arguments$vi <- dat$vi
    if (sampling_kind == "known_R") {
      known_R <- matrix(c(1, .4, .4, 1), 2L,
        dimnames = list(levels(dat$study), levels(dat$study)))
      arguments$R <- list(study = known_R)
      arguments$Rscale <- "none"
    }
    object <- do.call(bselmodel.mv, arguments)
    design <- object$formula_design$mu
    target <- design$name_map$jags_name[design$name_map$kind == "fixed" &
      design$name_map$term == "intercept"]
    stopifnot(length(target) == 1L, length(design$random_allocations) == 1L)
    information <- sum(forwardsolve(t(chol(V)), rep(1, 4L))^2)
    # Borderline rows switch from transformed at S=24 to ordinary at S=20.
    # Every injected default Dirichlet/Gamma coordinate must be positive.
    # Exact zero-scale fallback is tested separately above, without pretending
    # that a zero auxiliary is a draw from this continuous allocation graph.
    retained_sd <- rep(c(.5, sqrt(2 / (22 * information)), .04, .005), S / 4L)
    samples <- matrix(seq(-.3, .4, length.out = S), S, 1L,
      dimnames = list(NULL, target))
    source_names <- vapply(.data_selection_model(object$data)$sources$random,
      `[[`, character(1L), "name")
    source_retained <- vapply(.data_selection_model(object$data)$sources$random,
      function(source) isTRUE(source$retained), logical(1L))
    block_sd <- list()
    for (term in design$random_effects) {
      is_retained <- source_retained[match(term$block_name, source_names)]
      scale <- if (is_retained) retained_sd else rep(.2, S)
      block_sd[[term$block_name]] <- scale
      samples <- cbind(samples, scale)
      colnames(samples)[ncol(samples)] <- term$sd_parameter_names
      latent_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
        term, n_groups = term$n_groups, n_columns = 1L))
      for (group in seq_along(latent_names)) {
        samples <- cbind(samples, seq(-.2, .3, length.out = S) + group / 10)
        colnames(samples)[ncol(samples)] <- latent_names[[group]]
      }
    }
    allocation <- design$random_allocations[[1L]]
    variances <- do.call(cbind, lapply(unname(allocation$terms), function(block) block_sd[[block]]^2))
    total <- rowSums(variances)
    samples <- cbind(samples, sqrt(total))
    colnames(samples)[ncol(samples)] <- allocation$scale_name
    eta_name <- BayesTools:::.JAGS_prior_dirichlet_eta_name(allocation$weight_name)
    for (index in seq_len(ncol(variances))) {
      weight <- variances[, index] / total
      samples <- cbind(samples, weight, 2 * weight)
      colnames(samples)[ncol(samples) + c(-1L, 0L)] <- c(
        paste0(allocation$weight_name, "[", index, "]"), paste0(eta_name, "[", index, "]"))
    }
    # Materialize the real constructor's map/design around injected draws;
    # no JAGS model is fitted, compiled, adapted, or sampled.
    fit <- coda::mcmc.list(coda::mcmc(samples))
    class(fit) <- c("BayesTools_fit", class(fit))
    attr(fit, "formula_design") <- object$formula_design
    attr(fit, "prior_list") <- c(design$prior_list, .create_fit_priors(object$data, object$priors))
    fit <- BayesTools:::.bt_attach_parameter_map(fit)
    fit <- BayesTools:::.bt_attach_draw_geometry(fit)
    fit <- BayesTools:::.bt_attach_fit_contract(fit)
    object$fit <- fit
    context <- .iwmde_context(object)
    expect_false(any(nzchar(.random_allocation_inclusion_indicators(context$flat_prior_list))))
    if (sampling_kind == "known_R") {
      known_terms <- Filter(.random_effect_term_has_known_group_covariance, design$random_effects)
      expect_length(known_terms, 1L)
      term <- known_terms[[1L]]
      expect_equal(term$group_covariance$kernel, known_R, tolerance = 0)
      expect_true(source_retained[match(term$block_name, source_names)])
      # For this correlated two-group prior, 1'R^-1 1 = 2/(1+.4), not J=2.
      expect_equal(sum(term$group_covariance$precision), 2 / 1.4, tolerance = 1e-14)
      expect_false(.brma_mv_random_term_is_pure_intercept(term))
      expect_null(.iwmde_retained_location_plan(context, target, list(type = "primitive")))
      next
    }
    set.seed(9)
    plan <- .iwmde_plan(context, target, "qCMDE", list(samples = S, n_points = 20L),
      outputs = "ordinate", values = 0)
    expect_identical(plan$status, "ok")
    states <- plan$rows$row_states
    selected_rows <- plan$rows$estimator_rows
    transformed <- vapply(states, function(state) !is.null(state$conditioning_transform), logical(1L))
    expect_identical(selected_rows, seq_len(S))
    expect_identical(transformed, retained_sd > 0 & 2 / retained_sd^2 <= S * information)
    expect_true(any(transformed))
    expect_true(any(!transformed))
    expect_equal(plan$rows$conditioning_policy$information_cap, S * information, tolerance = 1e-12)
    expect_identical(.iwmde_plan_rows_provenance(plan)$conditioning_policy,
      plan$rows$conditioning_policy)
    values <- c(-.15, 0, .25)
    combined <- .iwmde_log_q_grid(context, target, values, states, plan$replacement)
    for (row in seq_along(states)) {
      separate <- .iwmde_log_q_grid(context, target, values, states[row], plan$replacement)
      expect_equal(combined[, row], as.numeric(separate), tolerance = 1e-10)
    }
    expect_equal(vapply(states[transformed], `[[`, numeric(1L), "baseline_log_q"),
      .iwmde_focal_log_prior_values(states[[which(transformed)[1L]]]$focal_prior,
        samples[transformed, target], target), tolerance = 0)
    set.seed(9)
    smaller <- .iwmde_plan(context, target, "qCMDE", list(samples = 20L, n_points = 20L),
      outputs = "ordinate", values = 0)
    expect_identical(smaller$rows$conditioning_policy$information_multiplier, 20L)
    expect_false(identical(smaller$plan_key, plan$plan_key))
    smaller_rows <- smaller$rows$estimator_rows
    expect_identical(vapply(smaller$rows$row_states,
      function(state) !is.null(state$conditioning_transform), logical(1L)),
      retained_sd[smaller_rows] > 0 & 2 / retained_sd[smaller_rows]^2 <= 20 * information)
    plan$method <- "iwmde"
    expect_null(.iwmde_retained_location_row_states(context, plan, selected_rows))
  }
})

.retained_location_mixture_fixture <- function(sign) {

  dat <- data.frame(yi = c(.1, -.2, .3, .05), study = factor(c("a", "a", "b", "b")))
  V <- diag(c(.04, .09, .16, .25))
  V[cbind(1:4, c(2L, 1L, 4L, 3L))] <- .01
  bias <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .5)), model = selection_model())
  object <- RoBMA.mv(yi = yi, V = V, random = ~1 | study, data = dat,
    measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = list(BayesTools::prior("normal", list(0, .7)),
      BayesTools::prior("normal", list(.3, .4))),
    prior_bias = list(bias, BayesTools::prior_PET("normal", list(0, 1)),
      BayesTools::prior_PEESE("normal", list(0, 1))),
    prior_bias_null = BayesTools::prior_none(),
    only_priors = TRUE, silent = TRUE, effect_direction = sign)
  design <- object$formula_design$mu
  target <- design$name_map$jags_name[design$name_map$kind == "fixed" &
    design$name_map$term == "intercept"]
  prior_list <- c(design$prior_list, .create_fit_priors(object$data, object$priors))
  focal <- prior_list[[target]]
  point <- which(vapply(focal, BayesTools::is.prior.point, logical(1L)))
  continuous <- which(!vapply(focal, BayesTools::is.prior.point, logical(1L)))
  # Include a spike as the first posterior row: candidate rows must select
  # their own focal priors instead of inheriting that first row's prior.
  branches <- expand.grid(effect = c(point, continuous), bias = seq_along(prior_list$bias),
    gate = c(0, 1), replicate = 1:2)
  S <- nrow(branches)
  samples <- matrix(seq(-.2, .4, length.out = S), S, 1L,
    dimnames = list(NULL, target))
  samples[branches$effect == point, target] <- 0
  append_column <- function(name, values) {
    samples <<- cbind(samples, rep(values, length.out = S))
    colnames(samples)[ncol(samples)] <<- name
  }
  append_column(paste0(target, "_indicator"), branches$effect)
  append_column("bias_indicator", branches$bias)
  allocation <- design$random_allocations[[1L]]
  append_column(allocation$source_node, .5)
  gate <- allocation$inclusion[[1L]]$indicator_name
  append_column(gate, branches$gate)
  term <- design$random_effects[[1L]]
  append_column(term$sd_parameter_names, .5 * branches$gate)
  latent <- as.vector(BayesTools:::.bt_random_effect_latent_names(term,
    n_groups = term$n_groups, n_columns = 1L))
  for (index in seq_along(latent)) {
    append_column(latent[[index]], seq(-.3, .2, length.out = S) + index / 10)
  }
  append_column("omega[1]", 1)
  append_column("omega[2]", .5)
  append_column("PET", .1 * vapply(prior_list$bias, BayesTools::is.prior.PET,
    logical(1L))[branches$bias])
  append_column("PEESE", .1 * vapply(prior_list$bias, BayesTools::is.prior.PEESE,
    logical(1L))[branches$bias])
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "formula_design") <- object$formula_design
  attr(fit, "prior_list") <- prior_list
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)
  object$fit <- fit
  list(context = .iwmde_context(object), target = target, latent = latent,
    gate = gate, V = V)
}

test_that("retained-location mixtures preserve branches, allocation gates and scalar priors", {

  for (sign in c("positive", "negative")) {
    fixture <- .retained_location_mixture_fixture(sign)
    context <- fixture$context
    target <- fixture$target
    samples <- context$posterior_samples
    expect_null(.iwmde_retained_location_plan(context, target, list(type = "primitive")))
    plan <- .iwmde_plan(context, target, "qCMDE", list(samples = Inf, n_points = 20L),
      outputs = "density")
    expect_identical(plan$status, "ok")
    rows <- plan$rows$estimator_rows
    states <- plan$rows$row_states
    transformed <- vapply(states, function(state) !is.null(state$conditioning_transform), logical(1L))
    expect_identical(transformed, samples[rows, fixture$gate] == 1)
    expect_true(any(transformed))
    expect_true(any(!transformed))
    expect_equal(plan$rows$point_masses$mass, 1 / 3)
    expect_identical(vapply(states, `[[`, character(1L), "active_key"),
      .iwmde_active_keys(context)[rows])
    expect_identical(plan$rows$conditioning_policy$information_multiplier, length(rows))
    expect_equal(plan$rows$conditioning_policy$information_cap,
      length(rows) * sum(solve(fixture$V)), tolerance = 1e-12)
    values <- c(-.15, 0, .25)
    combined <- .iwmde_log_q_grid(context, target, values, states, plan$replacement)
    for (position in which(transformed)) {
      state <- states[[position]]
      prior <- state$focal_prior
      row <- samples[rows[[position]], ]
      latent_means <- row[[target]] + .5 * row[fixture$latent]
      precision <- 1 / prior$parameters$sd^2 + length(fixture$latent) / .5^2
      conditional_mean <- (prior$parameters$mean / prior$parameters$sd^2 +
        sum(latent_means) / .5^2) / precision
      normalizer <- .iwmde_retained_location_normalizer(state)
      expect_equal(exp(combined[, position] - normalizer$log_normalizer),
        stats::dnorm(values, conditional_mean, 1 / sqrt(precision)), tolerance = 1e-12)
    }
    # The full likelihood itself stays unchanged under the translation, for
    # active selection, ordinary normal, PET and PEESE branches in either sign.
    for (key in unique(.iwmde_active_keys(context)[rows[transformed]])) {
      position <- which(transformed & vapply(states, `[[`, character(1L), "active_key") == key)[[1L]]
      row <- samples[rows[[position]], ]
      candidates <- matrix(rep(row, each = 2L), 2L, dimnames = list(NULL, names(row)))
      candidates[2L, target] <- row[[target]] + .2
      candidates[2L, fixture$latent] <- row[fixture$latent] - .2 / .5
      likelihood <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(context,
        candidates, .iwmde_active_setup(context, row, key), unit = "estimate")
      expect_equal(likelihood[1L], likelihood[2L], tolerance = 1e-11)
    }
    dynamic <- context
    dynamic$flat_prior_list[[target]][[2L]]$parameters$mean <- as.name(target)
    conditioning_row <- samples[rows[which(transformed)[[1L]]], ]
    expect_null(.iwmde_retained_location_plan(dynamic, target, list(type = "primitive"),
      row = conditioning_row))
    # The translation moves the fitted coordinate itself, so the accessor's
    # affine verdict is usable only in the fitted coordinate. A logged
    # intercept is affine in its own logarithm and keeps the generic route.
    expect_false(is.null(.iwmde_retained_location_plan(context, target,
      list(type = "primitive"), row = conditioning_row)))
    local({
      original_basis <- BayesTools::JAGS_formula_predictor_basis
      testthat::local_mocked_bindings(
        JAGS_formula_predictor_basis = function(...) {
          result <- original_basis(...)
          if (identical(result$status, "affine")) result$coordinate <- "log"
          result
        }, .package = "BayesTools")
      expect_null(.iwmde_retained_location_plan(context, target,
        list(type = "primitive"), row = conditioning_row))
    })
  }
})
