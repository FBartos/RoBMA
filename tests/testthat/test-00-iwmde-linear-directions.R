test_that("fixed linear directions retain Gaussian conditional and marginal laws", {

  row <- c(a = 0.2, b = -0.3, c = 0.4)
  weights <- c(a = 1, b = 1)
  context <- list(posterior_samples = matrix(row, 1L, dimnames = list(NULL, names(row))),
    flat_prior_list = lapply(row, function(x) BayesTools::prior("normal", list(mean = 0, sd = 1))),
    indicator_names = character(), row_cache = new.env(parent = emptyenv()),
    focal_prior_cache = new.env(parent = emptyenv()))
  covariance <- matrix(c(1.2, 0.2, -0.1, 0.2, 0.8, 0.15, -0.1, 0.15, 0.7), 3L)
  precision <- solve(covariance)
  center <- c(-0.1, 0.4, 0.2)
  root <- chol(covariance)
  full_weights <- c(1, 1, 0)
  current <- sum(row[names(weights)] * weights)
  directions <- list(c(a = 0.5, b = 0.5), c(b = 1), c(a = 1, c = -0.75))
  for (direction in directions) {
    resolved <- RoBMA:::.iwmde_linear_replacement_state(context,
      list(row = row, row_index = 1L), list(type = "linear", weights = weights, direction = direction))
    expect_true(resolved$valid)
    expect_identical(resolved$coefficients, direction)
    r <- c(a = 0, b = 0, c = 0)
    r[names(direction)] <- direction
    q <- as.numeric(crossprod(r, precision %*% r))
    conditional_sd <- 1 / sqrt(q)
    conditional_mean <- current - as.numeric(crossprod(r, precision %*% (row - center))) / q
    log_joint <- function(value) {
      vapply(value, function(t) {
        candidate <- row
        candidate[resolved$active_columns] <- candidate[resolved$active_columns] +
          (t - resolved$current) * resolved$coefficients
        whitened <- forwardsolve(t(root), candidate - center)
        -0.5 * sum(whitened^2) - sum(log(diag(root))) - 1.5 * log(2 * pi)
      }, numeric(1L))
    }
    peak <- log_joint(conditional_mean)
    integral <- stats::integrate(function(t) exp(log_joint(t) - peak), -Inf, Inf,
      rel.tol = 1e-11, abs.tol = 1e-11)$value
    grid <- conditional_mean + c(-2, -0.5, 0, 1.5) * conditional_sd
    expect_equal(exp(log_joint(grid) - peak) / integral,
      stats::dnorm(grid, conditional_mean, conditional_sd), tolerance = 1e-9)
    # Total variance of the conditional mean plus its constant conditional
    # variance equals the independent analytic marginal variance w' Sigma w.
    nuisance <- full_weights - as.numeric(precision %*% r) / q
    expect_equal(as.numeric(crossprod(nuisance, covariance %*% nuisance)) + 1 / q,
      as.numeric(crossprod(full_weights, covariance %*% full_weights)), tolerance = 1e-12)
  }
  # A second chart for the same row/target must not collide with the first.
  again <- RoBMA:::.iwmde_linear_replacement_state(context,
    list(row = row, row_index = 1L), list(type = "linear", weights = weights, direction = directions[[1L]]))
  expect_identical(again$coefficients, directions[[1L]])
  spec <- list(type = "primitive", parameter = "a", weights = c(a = 1),
    direction = c(a = 1, b = -1), conditioning_chart = "indicator_reference")
  expect_identical(RoBMA:::.iwmde_plan_execution_spec(spec)$type, "linear")
  expect_false(identical(RoBMA:::.iwmde_target_key("a", spec),
    RoBMA:::.iwmde_target_key("a", list(type = "primitive", parameter = "a"))))
})

test_that("linear support intersections follow every moved coordinate and its sign", {

  samples <- rbind(c(a = 0.2, b = 0.3), c(a = -0.1, b = 1.5))
  context <- list(posterior_samples = samples,
    flat_prior_list = list(a = BayesTools::prior("uniform", list(a = -1, b = 1)),
      b = BayesTools::prior("uniform", list(a = 0, b = 2))),
    indicator_names = character(), row_cache = new.env(parent = emptyenv()),
    focal_prior_cache = new.env(parent = emptyenv()))
  negative <- RoBMA:::.iwmde_linear_row_supports(context, 1:2, c(a = 1), c(a = 1, b = -2))
  positive <- RoBMA:::.iwmde_linear_row_supports(context, 1:2, c(a = 1), c(a = 1, b = 2))
  ordinary <- RoBMA:::.iwmde_linear_row_supports(context, 1:2, c(a = 1))
  expect_equal(unname(negative), rbind(c(-0.65, 0.35), c(-0.35, 0.65)), tolerance = 1e-14)
  expect_equal(unname(positive), rbind(c(0.05, 1), c(-0.85, 0.15)), tolerance = 1e-14)
  expect_equal(unname(ordinary), rbind(c(-1, 1), c(-1, 1)), tolerance = 1e-14)
})

test_that("affine joint grids reuse unchanged blocks and preserve the original log joint", {

  original_sum <- .log_lik_estimate_sum_from_setup
  captured_mu <- NULL
  testthat::local_mocked_bindings(.log_lik_estimate_sum_from_setup = function(setup) {
    captured_mu <<- setup$mu
    original_sum(setup)
  }, .package = "RoBMA")
  dat <- data.frame(yi = c(-.1, .2, .05, .3, -.2, .1),
    study = factor(rep(letters[1:3], each = 2)), group = factor(rep(letters[1:3], each = 2)))
  V <- diag(c(.04, .09, .05, .08, .06, .1))
  V[cbind(c(1, 2, 3, 4, 5, 6), c(2, 1, 4, 3, 6, 5))] <- .01
  bias <- BayesTools::prior_weightfunction("one-sided", .025, BayesTools::wf_fixed(c(1, .5)),
    model = selection_model(group = "study", other_random_effects = "condition",
      known_sampling_variance = "integrate"))
  for (sign in c("positive", "negative")) {
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
    stopifnot(length(intercept) == 1L, length(factor_name) == 1L,
      length(design$random_effects) == 1L)
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
    expect_false(is.null(spec$direction))
    execution <- .iwmde_plan_execution_spec(spec)
    replacement <- .iwmde_replacement_spec(context, intercept, execution)
    states <- .iwmde_row_states(context, 1L, intercept, execution, "q_grid_cmde")
    values <- samples[1L, intercept] + c(0, .1, -.2)
    candidates <- .iwmde_build_replacement_samples(context, intercept, values, states, replacement)
    batch <- list(candidates = candidates, valid_positions = which(candidates$valid), row_states = states)
    observed <- .iwmde_joint_affine_log_likelihood(context, intercept, values,
      candidates$samples, states[[1L]]$active_setup, batch, replacement)
    expect_false(is.null(observed))
    baseline <- .iwmde_predictor_setup(context, states, states[[1L]]$active_setup, "estimate")
    unchanged <- dat$group != levels(dat$group)[[1L]]
    expect_identical(unname(captured_mu[, unchanged, drop = FALSE]),
      matrix(rep(baseline$mu[1L, unchanged], each = length(values)), nrow = length(values)))
    original <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(context,
      candidates$samples, states[[1L]]$active_setup, unit = "estimate")
    expect_equal(observed, original, tolerance = 1e-10)
    prior_rows <- .resolve_fixed_prior_sample_columns(candidates$samples, states[[1L]]$prior_list)
    prior <- vapply(seq_len(nrow(prior_rows)), function(i) {
      BayesTools::JAGS_marglik_priors(prior_rows[i, ], states[[1L]]$prior_list)
    }, numeric(1L))
    observed_joint <- as.numeric(.iwmde_log_q_grid(context, intercept, values, states, replacement))
    expect_equal(observed_joint, original + prior, tolerance = 1e-10)
    ordinary <- list(type = "linear", weights = spec$weights)
    ordinary_candidates <- .iwmde_build_replacement_samples(
      context, intercept, values, states, ordinary)
    ordinary_batch <- list(candidates = ordinary_candidates,
      valid_positions = which(ordinary_candidates$valid), row_states = states)
    ordinary_result <- .iwmde_joint_affine_log_likelihood(context, intercept, values,
      ordinary_candidates$samples, states[[1L]]$active_setup, ordinary_batch, ordinary)
    ordinary_reference <- .iwmde_log_lik_from_posterior_samples_sum_active_branch(
      context, ordinary_candidates$samples, states[[1L]]$active_setup, unit = "estimate")
    expect_equal(ordinary_result, ordinary_reference, tolerance = 1e-10)
  }
})
