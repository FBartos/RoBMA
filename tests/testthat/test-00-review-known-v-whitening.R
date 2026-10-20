.review_whitened_known_v <- function(mixed = FALSE) {

  L <- rbind(c(1, 0), c(0, 1), c(1, 1))
  V <- tcrossprod(L)
  if (mixed) {
    full <- diag(c(1, 1, 2, 2))
    full[1:3, 1:3] <- V
    V <- full
  }
  .known_v_prepare(V, rep(TRUE, nrow(V)), "whitened", warn_singular = FALSE)
}

test_that("accepted whitening zeros reach the sampling factor adapter", {

  for (mixed in c(FALSE, TRUE)) {
    known <- .review_whitened_known_v(mixed)
    bytes <- serialize(known, NULL)
    factor <- .known_v_sampling_factor_plan(known)
    expect_false(is.null(factor))
    if (is.null(factor)) next
    expect_identical(serialize(known, NULL), bytes)
    expect_equal(ncol(factor$factor_plan$model_matrix), 2L)
    reconstruction <- factor$sampling_covariance + tcrossprod(factor$factor_plan$model_matrix)
    expect_lt(max(abs(reconstruction - .known_v_covariance_matrix(known))), 5e-15)
    if (mixed) expect_identical(factor$sampling_covariance[4, 4], 2)
    K <- .known_v_nrow(known)
    for (y in list(rep(0, K), c(1e-8 * c(1, 1, -1), if (mixed) 0))) {
      actual <- .marglik_covariance_plan_loglik(NULL, y, rep(0, K),
        factor$sampling_covariance, list(factor$factor_plan), list(factor$factor_state),
        .known_v_block_indices(.known_v_covariance_matrix(known)), rep(1e-18, K))
      expected <- if (all(y == 0)) 17.417144092998338 else -132.58285590700166
      if (mixed) expected <- expected + dnorm(0, 0, sqrt(2), log = TRUE)
      expect_lt(abs(actual - expected), 1e-11)
    }
  }
  pd <- .known_v_prepare(matrix(c(1, .25, .25, 1), 2), rep(TRUE, 2), "whitened")
  expect_null(.known_v_sampling_factor_plan(pd))
  damaged <- .review_whitened_known_v()
  damaged$whitening_blocks[[1]]$index <- c(1L, 1L, 3L)
  expect_error(.known_v_sampling_factor_plan(damaged), "whitening")
})

# A replay of a known-V model with one integrated random intercept per estimate
# (twelve estimates). The sampling covariance has four blocks of the rank-two
# 3 x 3 matrix tcrossprod(L), so its accepted whitening has exact zeros;
# 'definite' adds a diagonal that makes it positive definite (the control).
.review_whitened_replay <- function(definite = FALSE) {

  loadings <- rbind(c(1, 0), c(0, 1), c(1, 1)) * .25
  K <- 12L
  V <- matrix(0, K, K)
  for (block in 0:3) {
    V[block * 3L + 1:3, block * 3L + 1:3] <- tcrossprod(loadings)
  }
  if (definite) V <- V + diag(.02, K)
  dat <- data.frame(yi = c(.31, -.12, .44, .05, .52, .27, -.2, .38, .14, .61, .09, .33),
    id = factor(seq_len(K)))
  setup <- review_compile_object(suppressWarnings(brma.mv(yi = yi, V = V, data = dat,
    random = ~ 1 | id, measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE)))
  n <- 24L
  quantiles <- stats::qnorm((seq_len(n) - .5) / n)
  samples <- cbind(mu_intercept = .3 + .15 * quantiles,
    mu__xREx__id_intercept = .25 + .08 * rev(quantiles))
  list(object = review_attach_replay(setup, samples), V = V, yi = dat$yi,
       samples = samples)
}

test_that("IWMDE joint densities of a known-V model with exact whitening zeros match the dense normal", {

  singular <- .review_whitened_replay()
  definite <- .review_whitened_replay(definite = TRUE)
  # Dense multivariate-normal log density of y with mean mu and covariance
  # V + sigma^2 I (the random intercepts integrated out), by Cholesky.
  dense <- function(fixture, mu, sigma) {
    factor <- chol(fixture$V + sigma^2 * diag(length(fixture$yi)))
    z <- backsolve(factor, fixture$yi - mu, transpose = TRUE)
    -.5 * (length(z) * log(2 * pi) + 2 * sum(log(diag(factor))) + sum(z^2))
  }
  target <- function(fixture, parameter) {
    context <- .iwmde_context(fixture$object)
    spec    <- .iwmde_parameter_spec(context, parameter)
    list(context = context, parameter = parameter, spec = spec,
         states = .iwmde_row_states(context, 1:6, parameter = parameter, parameter_spec = spec),
         replacement = .iwmde_replacement_spec(context, parameter, spec))
  }
  sd_target <- function(fixture) {
    quantities <- BayesTools::parameter_catalog(fixture$object$fit)$quantities
    selected   <- quantities[quantities$role == "random_sd", , drop = FALSE]
    stopifnot(nrow(selected) == 1L)
    density_target <- .brma_random_parameter_density_target(fixture$object,
      .brma_random_parameter_io_labels(selected, "selector"))
    context <- .iwmde_context(fixture$object)
    spec    <- .iwmde_parameter_spec(context, density_target$parameter,
                                     density_target$parameter_spec)
    list(context = context, parameter = density_target$parameter, spec = spec,
         states = .iwmde_row_states(context, 1:6, parameter = density_target$parameter,
                                    parameter_spec = spec),
         replacement = .iwmde_replacement_spec(context, density_target$parameter, spec))
  }

  # The factor-plan branches are the ones under test: the counters show that
  # the whitened-zero likelihood (sampling-factor sum for a scale target, the
  # location quadratic for a location target) is the one that produced q.
  calls <- new.env(parent = emptyenv())
  counted <- function(name) {
    original <- get(name, envir = asNamespace("RoBMA"), inherits = FALSE)
    calls[[name]] <- 0L
    function(...) {
      calls[[name]] <- calls[[name]] + 1L
      original(...)
    }
  }
  testthat::local_mocked_bindings(
    .known_v_marginal_factor_plan = counted(".known_v_marginal_factor_plan"),
    .iwmde_log_lik_known_v_random_marginal_sum_from_samples =
      counted(".iwmde_log_lik_known_v_random_marginal_sum_from_samples"),
    .iwmde_normal_location_likelihood_change_known_v_random =
      counted(".iwmde_normal_location_likelihood_change_known_v_random"),
    .package = "RoBMA"
  )

  # Agreement with the dense density is to rounding of a 12 x 12 Cholesky solve
  # whose condition number is about 20 here (measured 6e-15); 1e-10 on the log
  # scale leaves room for the platform's BLAS. Everything that does not depend
  # on the value (the other priors of the row) is the baseline offset.
  rows <- 1:6
  offsets <- function(states) {
    vapply(states, `[[`, numeric(1), "baseline_log_prior") -
      vapply(states, `[[`, numeric(1), "baseline_focal_log_prior")
  }
  check_dense <- function(fixture, case, values, mean_at, sd_at) {
    log_q <- .iwmde_log_q_grid(case$context, case$parameter, values, case$states,
                               case$replacement)
    focal <- .iwmde_focal_log_prior_values(case$states[[1L]][["focal_prior"]], values,
                                           case$parameter)
    reference <- outer(seq_along(values), seq_along(rows), Vectorize(function(g, s) {
      dense(fixture, mean_at(values[g], rows[s]), sd_at(values[g], rows[s]))
    }))
    expect_identical(dim(log_q), dim(reference))
    expect_lt(max(abs(sweep(log_q - unname(focal), 2L, offsets(case$states)) - reference)), 1e-10)
  }

  expect_false(is.null(.iwmde_known_v_random_marginal_setup(
    .iwmde_context(singular$object))[["sampling_factors"]]))
  expect_null(.iwmde_known_v_random_marginal_setup(
    .iwmde_context(definite$object))[["sampling_factors"]])

  location <- target(singular, "mu_intercept")
  check_dense(singular, location, c(.1, .3, .45),
              mean_at = function(value, row) value,
              sd_at   = function(value, row) singular$samples[row, "mu__xREx__id_intercept"])
  expect_gt(calls[[".iwmde_normal_location_likelihood_change_known_v_random"]], 0L)

  scale <- sd_target(singular)
  scale_values <- c(.15, .25, .4)
  calls[[".iwmde_log_lik_known_v_random_marginal_sum_from_samples"]] <- 0L
  check_dense(singular, scale, scale_values,
              mean_at = function(value, row) singular$samples[row, "mu_intercept"],
              sd_at   = function(value, row) value)
  expect_gt(calls[[".iwmde_log_lik_known_v_random_marginal_sum_from_samples"]], 0L)
  expect_gt(calls[[".known_v_marginal_factor_plan"]], 0L)

  # The shortcuts that assume a positive-definite sampling covariance decline
  # on the whitened-zero model and take over on the positive-definite control.
  active <- scale$states[[1L]]$active_setup
  decline <- lapply(list(
    group_iid = .iwmde_log_q_grid_known_v_random_group_iid,
    affine    = .iwmde_log_q_grid_known_v_random_affine,
    factor    = .iwmde_log_q_grid_known_v_random_factor),
    function(shortcut) shortcut(scale$context, scale$parameter, scale_values,
                                scale$states, scale$replacement, active))
  expect_true(all(vapply(decline, is.null, logical(1L))))
  control <- sd_target(definite)
  control_active <- control$states[[1L]]$active_setup
  for (shortcut in list(.iwmde_log_q_grid_known_v_random_group_iid,
                        .iwmde_log_q_grid_known_v_random_affine)) {
    expect_true(is.matrix(shortcut(control$context, control$parameter, scale_values,
                                   control$states, control$replacement, control_active)))
  }
})
