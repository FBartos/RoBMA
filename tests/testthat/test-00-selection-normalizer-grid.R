# Deterministic posterior rows exercise the real compiled likelihood and priors;
# constructing these objects does not fit or sample a model.
.normalizer_grid_test_fixture <- function(sign = "positive", weights = c(1, .5),
                                          covariance = .01) {

  dat <- data.frame(yi = c(-.1, .2, .05, .3, -.2, .1),
    study = factor(rep(letters[1:3], each = 2)),
    group = factor(rep(letters[1:3], each = 2)))
  V <- diag(c(.04, .09, .05, .08, .06, .1))
  V[cbind(1:6, c(2, 1, 4, 3, 6, 5))] <- covariance
  bias <- BayesTools::prior_weightfunction("one-sided",
    if (length(weights) == 2L) .025 else c(.025, .05),
    BayesTools::wf_fixed(weights), model = selection_model(group = "study",
      other_random_effects = "condition", known_sampling_variance = "integrate"))
  object <- bselmodel.mv(yi = yi, V = V, mods = ~ group, random = ~1 | study,
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("t", list(0, .7, 4)),
    prior_mods = BayesTools::prior_factor("t", list(0, .6, 4), contrast = "treatment"),
    prior_heterogeneity = BayesTools::prior("point", list(location = .2)),
    prior_bias = bias, effect_direction = sign, only_priors = TRUE, silent = TRUE)
  design <- object$formula_design$mu
  map <- design$name_map
  intercept <- map$jags_name[map$kind == "fixed" & map$term == "intercept"]
  factor_name <- map$jags_name[map$kind == "fixed" & map$term == "group"]
  stopifnot(length(intercept) == 1L, length(factor_name) == 1L)
  samples <- cbind(c(.1, -.2), c(.15, .22), c(-.25, .18))
  colnames(samples) <- c(intercept, paste0(factor_name, "[", 1:2, "]"))
  term <- design$random_effects[[1L]]
  samples <- cbind(samples, rep(.2, 2L))
  colnames(samples)[ncol(samples)] <- term$sd_parameter_names
  latent <- as.vector(BayesTools:::.bt_random_effect_latent_names(term,
    n_groups = term$n_groups, n_columns = 1L))
  for (index in seq_along(latent)) {
    samples <- cbind(samples, c(-.1, .2) + index / 20)
    colnames(samples)[ncol(samples)] <- latent[[index]]
  }
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "formula_design") <- object$formula_design
  attr(fit, "prior_list") <- c(design$prior_list, .create_fit_priors(object$data, object$priors))
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)
  object$fit <- fit
  context <- .iwmde_context(object)
  spec <- .iwmde_linear_conditioning_spec(context,
    list(type = "primitive", parameter = intercept), "qCMDE")
  execution <- .iwmde_plan_execution_spec(spec)
  list(context = context, parameter = intercept, execution = execution,
    replacement = .iwmde_replacement_spec(context, intercept, execution),
    states = .iwmde_row_states(context, 1L, intercept, execution, "q_grid_cmde"))
}

.normalizer_grid_test_joint <- function(fixture, values, geometry_limit = NULL) {

  context <- fixture$context
  state <- fixture$states[[1L]]
  weights <- fixture$execution$weights
  direction <- fixture$execution$direction
  if (is.null(direction)) direction <- weights / sum(weights^2)
  current <- sum(state$row[names(weights)] * weights)
  # Independent full-coordinate line construction, followed by the complete
  # scalar prior evaluator and the original, uninterpolated likelihood.
  candidates <- matrix(rep(state$row, each = length(values)), length(values),
    dimnames = list(NULL, names(state$row)))
  for (i in seq_along(values)) {
    candidates[i, names(direction)] <- state$row[names(direction)] +
      (values[[i]] - current) * direction
  }
  prior_rows <- .resolve_fixed_prior_sample_columns(candidates, state$prior_list)
  prior <- vapply(seq_len(nrow(prior_rows)), function(i) {
    BayesTools::JAGS_marglik_priors(prior_rows[i, ], state$prior_list)
  }, numeric(1L))
  reference <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(context,
    candidates, state$active_setup, unit = "estimate") + prior
  context$normalizer_grid <- .selection_normalizer_grid(context, values, 1L)
  stopifnot(is.environment(context$normalizer_grid))
  if (!is.null(geometry_limit)) context$normalizer_grid$geometry_limit <- geometry_limit
  actual <- as.numeric(.iwmde_log_q_grid(context, fixture$parameter, values,
    fixture$states, fixture$replacement))
  list(actual = actual, reference = reference, shared = context$normalizer_grid,
    diagnostic = .selection_normalizer_grid_diagnostics(context$normalizer_grid))
}

test_that("requested-node interpolation preserves full joint densities for both signs", {

  # Uneven spacing and a close pair exercise actual-node geometry, not the
  # ideal Chebyshev-root remainder. More than 80 points activate interpolation.
  values <- sort(unique(c(seq(-.4, .4, length.out = 161L), .031, .031 + 1e-14)))
  for (sign in c("positive", "negative")) {
    fixture <- .normalizer_grid_test_fixture(sign)
    result <- .normalizer_grid_test_joint(fixture, values)
    diagnostic <- result$diagnostic
    expect_true(diagnostic$used)
    expect_gt(diagnostic$interpolated_points, 0)
    expect_gt(diagnostic$constant_evaluations, 0)
    expect_true(all(is.finite(result$actual)))
    expect_lte(max(abs(result$actual - result$reference)),
      diagnostic$max_log_likelihood_error + 1e-9)
    expect_lte(diagnostic$anchor_evaluations + diagnostic$fallback_evaluations +
      diagnostic$constant_evaluations, diagnostic$unique_requested_normalizers)
    expect_false(diagnostic$untracked_path)
    expect_identical(diagnostic$unknown_error_points, 0L)
    expect_match(diagnostic$error_scope, "finite requested ordinates", fixed = TRUE)
  }
})

test_that("scalar fixed coefficients retain checked full-event normalizer grids", {

  values <- seq(-.4, .4, length.out = 161L)
  for (sign in c("positive", "negative")) {
    fixture <- .normalizer_grid_test_fixture(sign)
    context <- fixture$context
    parameter <- fixture$parameter
    spec <- list(type = "primitive", parameter = parameter)
    states <- .iwmde_row_states(context, 1L, parameter, spec, "q_grid_cmde")
    state <- states[[1L]]
    candidates <- matrix(rep(state$row, each = length(values)), length(values),
      dimnames = list(NULL, names(state$row)))
    candidates[, parameter] <- values
    prior_rows <- .resolve_fixed_prior_sample_columns(candidates, state$prior_list)
    prior <- vapply(seq_len(nrow(prior_rows)), function(i) {
      BayesTools::JAGS_marglik_priors(prior_rows[i, ], state$prior_list)
    }, numeric(1L))
    reference <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(context,
      candidates, state$active_setup, unit = "estimate") + prior
    context$normalizer_grid <- .selection_normalizer_grid(context, values, 1L)
    actual <- .iwmde_log_q_grid(context, parameter, values, states,
      list(type = "scalar"))
    diagnostic <- .selection_normalizer_grid_diagnostics(context$normalizer_grid)
    expect_true(diagnostic$used)
    expect_gt(diagnostic$interpolated_points, 0)
    expect_false(diagnostic$untracked_path)
    expect_lte(max(abs(as.numeric(actual) - reference)),
      diagnostic$max_log_likelihood_error + 1e-9)
  }
})

test_that("normalizer geometry has bounded storage and unsafe-node fallback", {

  fixture <- .normalizer_grid_test_fixture()
  values <- seq(-.4, .4, length.out = 81L)
  result <- .normalizer_grid_test_joint(fixture, values, geometry_limit = 0)
  expect_false(result$diagnostic$used)
  expect_gt(result$diagnostic$fallback_evaluations, 0)
  expect_equal(result$actual, result$reference, tolerance = 1e-10)
  expect_identical(result$shared$geometry_bytes, 0)
  expect_lte(result$diagnostic$fallback_evaluations + result$diagnostic$constant_evaluations,
    result$diagnostic$unique_requested_normalizers)
  # Distinct finite nodes whose difference is subnormal must not be treated
  # as a zero interpolation remainder. Only this unsafe query is declined.
  tiny <- .Machine$double.xmin
  values <- sort(unique(c(values, tiny, tiny * (1 + .Machine$double.eps))))
  shared <- .selection_normalizer_grid(fixture$context, values, 1L)
  leaves <- .selection_normalizer_grid_geometry(shared, seq_along(values))
  covered <- unique(unlist(lapply(leaves, `[[`, "ids"), use.names = FALSE))
  expect_false(match(tiny * (1 + .Machine$double.eps), values) %in% covered)
  expect_lte(shared$geometry_bytes, shared$geometry_limit)
  previous <- options(RoBMA.known_v_covariance_max_bytes = 1)
  on.exit(options(previous), add = TRUE)
  expect_null(.selection_normalizer_grid(fixture$context, values, 1L))
})

test_that("unavailable anchor errors remain unknown through direct evaluation", {

  # Negative sampling covariances use the ordinary native route.
  # Its QMC MCSE is not an absolute normalizer-error bound for interpolation.
  fixture <- .normalizer_grid_test_fixture(weights = c(1, .98, .99), covariance = -.01)
  result <- .normalizer_grid_test_joint(fixture, seq(-.4, .4, length.out = 81L))
  expect_false(result$diagnostic$used)
  expect_gt(result$diagnostic$unknown_error_points, 0)
  expect_true(is.na(result$diagnostic$max_log_likelihood_error))
  expect_equal(result$actual, result$reference, tolerance = 1e-10)
  empty <- .selection_normalizer_grid(fixture$context, seq(-1, 1, length.out = 81L), 1L)
  expect_true(is.na(.selection_normalizer_grid_diagnostics(empty)$max_log_likelihood_error))
  result$shared$untracked <- TRUE
  diagnostic <- .selection_normalizer_grid_diagnostics(result$shared)
  expect_true(diagnostic$untracked_path)
  expect_true(is.na(diagnostic$conditional_relative_error))
})

test_that("normalizer interpolation leaves Gaussian and constant-weight laws exact", {

  values <- seq(-.3, .3, length.out = 81L)
  size <- length(values)
  yi <- c(.1, -.2)
  sei <- c(.2, .3)
  covariance <- matrix(c(.04, .01, .01, .09), 2L)
  means <- cbind(.1 + values, -.05 + .5 * values)
  packed <- matrix(covariance[lower.tri(covariance, diag = TRUE)], size, 3L, byrow = TRUE)
  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .5)), model = selection_model(group = "paper",
      other_random_effects = "integrate", known_sampling_variance = "integrate"))
  selection <- .selection_spec(list(outcome = list(bias = prior)), yi, sei, "positive")
  selection$omega <- matrix(c(1, .5), size, 2L, byrow = TRUE)
  selection$kernel_mode <- rep(SELKERNEL_STEP, size)
  selection$vector_rule <- integer(size)
  selection$alpha <- numeric(size)
  selection$phack_kind <- integer(size)
  selection$use_normal <- rep(FALSE, size)
  plan <- list(points_per_scramble = 512L, scrambles = 8L, relative_tolerance = .005,
    designs = list("2" = BayesTools::selection_qmc_design(4L, 512L, 8L, seed = 173L)),
    quadrature = .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS))
  # The shared state is not consumed by either analytic bypass. No fake
  # normalizer or likelihood callback is involved in this test.
  metadata <- list(state = list(shared = new.env(parent = emptyenv())))
  reference <- vapply(seq_len(size), function(i) {
    mvtnorm::dmvnorm(yi, means[i, ], covariance, log = TRUE)
  }, numeric(1L))
  for (mode in c("gaussian", "constant")) {
    branch <- selection
    if (mode == "gaussian") {
      branch$kernel_mode[] <- SELKERNEL_NORMAL
      branch$use_normal[] <- TRUE
    } else {
      branch$omega[,] <- 2
    }
    expect_null(.selection_normalizer_grid_loglik(yi, means, packed, sei,
      branch, plan, 2L, metadata))
    actual <- .selection_joint_dense_loglik_block(yi, means, packed, sei, branch, plan, 2L)
    expect_equal(actual, reference, tolerance = 1e-11)
  }
})
