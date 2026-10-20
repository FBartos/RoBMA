# Genuine compiler/fit-contract supplied-row replays; these are transport and
# declared-target regressions, not MCMC or posterior calibration fixtures.
.review_h_contribution <- function(multiplier = 2) {

  dat <- data.frame(yi = c(-.2, .3, 0, .4, .1, .2), sei = .5,
    g = factor(rep(c("A", "B"), 3)))
  slope <- BayesTools::prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(slope, "multiply_by") <- multiplier
  setup <- review_compile_object(brma(yi = yi, sei = sei, mods = ~g, data = dat,
    measure = "GEN", prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_mods = list(g = slope), prior_heterogeneity = if (is.character(multiplier)) {
      BayesTools::prior("normal", list(0, .5), list(0, Inf))
    } else BayesTools::prior("point", list(.1)),
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE))
  n <- 240L
  quantiles <- stats::qnorm((seq_len(n) - .5) / n)
  samples <- cbind(mu_intercept = .1 + .12 * quantiles,
    mu_g = .2 + .15 * quantiles[rev(seq_len(n))], tau = .1)
  review_attach_replay(setup, samples)
}

.review_h_empty_condition <- function() {

  dat <- data.frame(yi = c(.1, .2, .15, .25), sei = .2, x = c(-1, 1, -1, 1))
  setup <- review_compile_object(BMA(yi = yi, sei = sei, mods = ~ x, data = dat,
    measure = "GEN", prior_mods = list(intercept = BayesTools::prior("normal", list(.4, .2))),
    prior_heterogeneity = BayesTools::prior("point", list(0)), prior_heterogeneity_null = NULL,
    prior_unit_information_sd = 1, standardize_continuous_predictors = FALSE,
    only_priors = TRUE, silent = TRUE))
  slope <- setup$compiled$mu$prior_list$mu_x
  components <- vapply(slope, function(prior) attr(prior, "component"), character(1L))
  alternative <- which(components == "alternative")
  null <- which(components == "null")
  fit_priors <- .create_fit_priors(setup$object$data, setup$object$priors)
  stopifnot(BayesTools::is.prior.mixture(fit_priors$tau), length(fit_priors$tau) == 1L)
  tau_component <- seq_along(fit_priors$tau)
  samples <- cbind(mu_intercept = .2 + .025 * stats::qnorm((seq_len(120) - .5) / 120),
    mu_x = c(rep(0, 60), seq(.01, .15, length.out = 60)),
    mu_x_indicator = rep(c(null, alternative), each = 60), tau = 0,
    tau_indicator = rep(tau_component, 120L))
  monitors <- unique(c(BayesTools::JAGS_to_monitor(c(setup$compiled$mu$prior_list, fit_priors)),
    .create_jags_formula_args(setup$object$data, setup$object$priors)$add_parameters))
  stopifnot(length(setdiff(monitors, colnames(samples))) == 0L)
  review_attach_replay(setup, samples)
}

.review_h_random <- function(kind = c("fisher", "scale", "common"), structure = "cs") {

  kind <- match.arg(kind)
  dat <- data.frame(yi = c(.1, .2, .3, .4), study = factor(c("s1", "s1", "s2", "s2")),
    outcome = factor(c("a", "b", "a", "b")), x = c(-1, 1, -1, 1))
  half <- BayesTools::prior("normal", list(0, .5), list(0, Inf))
  random <- if (kind == "common") ~diag(1 + x | study) else if (structure == "ar1") {
    ~ar1(outcome | study)
  } else ~cs(outcome | study)
  args <- list(yi = dat$yi, V = diag(.04, 4), random = random,
    data = dat, measure = "GEN",
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = if (kind == "fisher") {
      BayesTools::prior_random(sd = half, cor = BayesTools::prior("normal", list(0, .5)))
    } else half, prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  if (kind == "scale") args$scale <- ~x
  object <- do.call(brma.mv, args)
  setup <- review_compile_object(object)
  term <- setup$object$formula_design$mu$random_effects[[1L]]
  priors <- unlist(lapply(setup$compiled, `[[`, "prior_list"), recursive = FALSE)
  names(priors) <- unlist(lapply(setup$compiled, function(compiled) names(compiled$prior_list)))
  priors <- c(priors, .create_fit_priors(object$data, object$priors))
  columns <- review_prior_columns(priors)
  n <- 120L
  samples <- matrix(0, n, length(columns), dimnames = list(NULL, columns))
  samples[, "mu_intercept"] <- seq(.1, .3, length.out = n)
  if (kind == "common") {
    allocation <- setup$object$formula_design$mu$random_allocations[[1L]]
    samples[, allocation$source_node] <- seq(.2, .4, length.out = n)
    weights <- c(.4, .6)
    weight_names <- paste0(allocation$weight_name, "[", 1:2, "]")
    samples[, weight_names] <- matrix(rep(weights, each = n), n)
    eta <- samples[, weight_names, drop = FALSE]
    colnames(eta) <- paste0(BayesTools:::.JAGS_prior_dirichlet_eta_name(allocation$weight_name), "[", 1:2, "]")
    sds <- vapply(seq_len(2), function(i) samples[, allocation$source_node] * sqrt(2 * weights[[i]]), numeric(n))
    colnames(sds) <- term$sd_binding$sd_component_names
    samples <- cbind(samples, eta, sds)
  } else {
    samples[, term$correlation$sample_name] <- seq(-.2, .3, length.out = n)
    if (kind == "fisher") {
      samples[, term$sd_binding$source$name] <- seq(.2, .4, length.out = n)
      rho <- matrix(tanh(samples[, term$correlation$sample_name]), ncol = 1,
        dimnames = list(NULL, term$correlation$rho_name))
      samples <- cbind(samples, rho)
    } else {
      samples[, "log_tau_intercept"] <- seq(.15, .35, length.out = n)
      samples[, "log_tau_x"] <- seq(-.1, .1, length.out = n)
    }
  }
  latent <- grep("_xRE_Zx$", setup$compiled$mu$add_parameters, value = TRUE)
  indices <- expand.grid(group = seq_len(term$n_groups), column = seq_len(term$n_columns))
  values <- vapply(seq_len(nrow(indices)), function(i) .1 * sin(seq_len(n) / (i + 1)), numeric(n))
  colnames(values) <- paste0(latent, "[", indices$group, ",", indices$column, "]")
  review_attach_replay(setup, cbind(samples, values))
}

.review_h_ordinate_refusal <- function() {

  dat <- data.frame(yi = rep(c(.1, 10), each = 3), sei = .2,
    g = factor(rep(c("A", "B"), each = 3)))
  setup <- review_compile_object(brma(yi = yi, sei = sei, mods = ~0+g, data = dat,
    measure = "GEN", prior_effect = BayesTools::prior("point", list(0)),
    prior_mods = list(g = BayesTools::prior_factor("normal", list(0, 1), contrast = "independent")),
    prior_heterogeneity = BayesTools::prior("point", list(.1)), prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE))
  name <- review_coefficient_name(setup, "g")
  columns <- BayesTools:::.JAGS_prior_factor_names(name, setup$compiled$mu$prior_list[[name]])
  n <- 240L
  quantiles <- stats::qnorm((seq_len(n) - .5) / n)
  samples <- cbind(.1 + .12 * quantiles, 9.84 + .12 * rev(quantiles), mu_intercept = 0, tau = .1)
  colnames(samples)[1:2] <- columns
  review_attach_replay(setup, samples)
}

.review_h_atomic <- function() {

  dat <- data.frame(yi = c(-.1, 0, .1, .2), sei = .2, x = c(-1, 0, 1, 2))
  point <- function(value) BayesTools::prior("point", list(value))
  setup <- review_compile_object(BMA.norm(yi = yi, sei = sei, mods = ~x, data = dat,
    measure = "GEN", prior_effect = list(point(-.2), point(.2)), prior_effect_null = NULL,
    prior_mods = list(x = point(0)), prior_mods_null = NULL, prior_heterogeneity = point(.1),
    prior_heterogeneity_null = NULL, standardize_continuous_predictors = FALSE,
    only_priors = TRUE, silent = TRUE))
  samples <- cbind(mu_intercept = rep(c(-.2, .2), c(40, 60)),
    mu_intercept_indicator = rep(1:2, c(40, 60)), mu_x = 0, tau = .1)
  map <- setup$object$formula_design$mu$name_map
  missing <- setdiff(map$jags_name[map$role == "_indicator"], colnames(samples))
  if (length(missing)) samples <- cbind(samples,
    matrix(1, nrow(samples), length(missing), dimnames = list(NULL, missing)))
  review_attach_replay(setup, samples)
}

test_that("E2-1 declared empty marginal cells retain their actual event", {

  object <- .review_h_empty_condition()
  means <- marginal_means(object, bf = FALSE, density_method = "KDE", n_samples = 128)
  cells <- means$inference$conditional$mu_x
  for (parameter in c("intercept", "x")) {
    specs <- .iwmde_marginal_means_specs(means, parameter, "conditional", NULL)
    expect_gt(length(specs), 0L)
    for (spec in specs) {
      cell <- spec$marginal_samples
      declared <- BayesTools::posterior_metadata(cell, "condition")
      expect_identical(spec$condition_key, declared$condition_key)
      expect_identical(spec$condition_event, declared$condition_event)
      expect_identical(spec$prior_density_context, BayesTools::posterior_metadata(cell, "prior_context"))
      if (isTRUE(declared$averaged)) expect_length(spec$conditional, 0L)
    }
  }
  empty <- means$inference$conditional$mu_intercept$intercept
  for (defect in c("missing", "conflict")) {
    changed <- empty
    metadata <- BayesTools::posterior_metadata(changed, "condition")
    if (defect == "missing") metadata$condition_event <- NULL else metadata$condition_key <- "stale"
    BayesTools::posterior_metadata(changed, "condition") <- metadata
    expect_error(.iwmde_sample_condition_metadata(changed, require_child_condition = TRUE), "metadata")
  }
  nonempty <- cells[[which(vapply(cells, function(cell)
    length(BayesTools::posterior_metadata(cell, "condition")$effective_conditional) > 0L, logical(1L)))[[1L]]]]
  metadata <- BayesTools::posterior_metadata(nonempty, "condition")
  metadata$averaged <- TRUE
  BayesTools::posterior_metadata(nonempty, "condition") <- metadata
  expect_error(.iwmde_sample_condition_metadata(nonempty, require_child_condition = TRUE), "metadata")
  # A public point reaches estimation; a subsequent diagnostic refusal is kept.
  result <- tryCatch(hypothesis(means, "intercept = 0", density_method = "qCMDE", seed = 19), error = identity)
  expect_false(inherits(result, "error") && grepl("metadata is incomplete", conditionMessage(result), fixed = TRUE))
})

test_that("E2-1 complete declared monitors admit public empty-cell estimation", {

  object <- .review_h_empty_condition()
  samples <- as.matrix(object$fit$mcmc)
  priors <- attr(object$fit, "prior_list")
  monitors <- unique(c(BayesTools::JAGS_to_monitor(priors),
    .create_jags_formula_args(object$data, object$priors)$add_parameters))
  expect_length(setdiff(monitors, colnames(samples)), 0L)
  expect_identical(unique(samples[, "tau_indicator"]), as.numeric(seq_along(priors$tau)))
  means <- marginal_means(object, bf = FALSE, density_method = "qCMDE", n_samples = 128)
  expect_s3_class(means, "marginal_means.brma")
  densities <- unlist(lapply(means$inference$averaged, function(cells)
    lapply(cells, function(cell) BayesTools::posterior_metadata(cell, "posterior_density"))),
    recursive = FALSE)
  expect_length(densities, 4L)
  expect_true(all(vapply(densities, is.null, logical(1L))))
  failures <- vapply(means$density_diagnostics$averaged, function(record)
    .iwmde_diagnostics_density_failure_reason(record$diagnostics), character(1L))
  expect_true(all(grepl("120 estimator rows across.*2 active states", failures)))
  captured <- NULL
  compute <- BayesTools::hypothesis_BF
  testthat::local_mocked_bindings(
    hypothesis_BF = function(posterior, ...) {
      captured <<- posterior
      compute(posterior = posterior, ...)
    },
    .package = "BayesTools"
  )
  point <- hypothesis(means, "intercept = .1", density_method = "qCMDE", seed = 19)
  child <- captured[["conditional"]][["mu_intercept"]][["intercept"]]
  .expect_iwmde_child_identity(child, "mu_intercept", "intercept", .1)
  expect_identical(BayesTools::posterior_metadata(child, "condition"),
    BayesTools::posterior_metadata(means$inference$conditional$mu_intercept$intercept, "condition"))
  expect_true(is.finite(as.numeric(point$BF)))
  diagnostics <- density_diagnostics(point)
  expect_identical(diagnostics$requested_samples, 500)
  expect_identical(diagnostics$eligible_rows, 120L)
  expect_identical(diagnostics$evaluated_rows, 120L)
  expect_true(diagnostics$all_rows_used)
  # Savage-Dickey value of this conjugate fixture. The design column x sums to
  # zero, so the posterior of the intercept given any retained row (tau = 0) is
  # the same normal: precision 1 / .2^2 + 4 / .2^2 and mean
  # (.4 / .2^2 + sum(yi) / .2^2) / precision. The qCMDE ordinate is the mean of
  # these exact conditional ordinates, and the prior ordinate is that of
  # normal(.4, .2). The allowance is the numerical error reported by the
  # diagnostics (an error e of the ordinate moves the log Bayes factor by
  # log1p(e) at most), plus rounding.
  precision <- 125
  expected <- stats::dnorm(.1, .4, .2, log = TRUE) -
    stats::dnorm(.1, (.4 / .2^2 + sum(c(.1, .2, .15, .25)) / .2^2) / precision,
                 sqrt(1 / precision), log = TRUE)
  reported <- diagnostics$normalization_relative_error + diagnostics$stability_relative_error +
    diagnostics$lattice_resolution_relative_error + diagnostics$truncation_ordinate_bound
  expect_lte(abs(attr(point, "raw_log_BF") - expected),
             log1p(reported) + 64 * .Machine$double.eps)
})

# Draws of the exact posterior of a two-group treatment-contrast model (normal
# priors, tau = .1, six observations), so that the qCMDE row weights are
# well behaved and the density diagnostics of the cells pass.
.review_h_posterior_replay <- function(n = 240L) {

  dat <- data.frame(yi = seq(-.1, .2, length.out = 6L), sei = .1,
    g = factor(rep(c("A", "B"), each = 3L)))
  setup <- review_compile_object(brma(yi = yi, sei = sei, mods = ~ g, data = dat,
    measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    prior_mods = list(g = BayesTools::prior_factor("normal", list(mean = 0, sd = 1),
      contrast = "treatment")),
    prior_heterogeneity = BayesTools::prior("point", list(location = .1)),
    only_priors = TRUE, silent = TRUE))
  design     <- cbind(1, as.numeric(dat$g == "B"))
  variance   <- .1^2 + .1^2
  covariance <- solve(crossprod(design) / variance + diag(2L))
  location   <- covariance %*% crossprod(design, dat$yi) / variance
  withr::local_seed(6123)
  draws <- matrix(stats::rnorm(2L * n), n, 2L) %*% chol(covariance) +
    matrix(location, n, 2L, byrow = TRUE)
  g_name <- review_coefficient_name(setup, "g")
  slope  <- draws[, 2L, drop = FALSE]
  colnames(slope) <- BayesTools:::.JAGS_prior_factor_names(
    g_name, setup[["compiled"]][["mu"]][["prior_list"]][[g_name]]
  )
  samples <- cbind(slope, draws[, 1L], tau = rep(.1, n))
  colnames(samples)[2L] <- review_coefficient_name(setup, "intercept")

  review_attach_replay(setup, samples)
}

test_that("E2-2 actual attachment retains requested curve precision", {

  object <- .review_h_posterior_replay()
  for (target in c(.1, .05)) {
    result <- marginal_means(object, bf = FALSE, parameter = "g", density_method = "qCMDE",
      density_control = list(target_relative_mcse = target, samples = 120L))
    expect_identical(result$density_settings$target_relative_mcse, target)
    cells <- result$inference$averaged$mu_g
    expect_named(cells, c("A", "B"))
    for (cell in cells) {
      density <- BayesTools::posterior_metadata(cell, "posterior_density")
      expect_false(is.null(density))
      expect_identical(density$iwmde_provenance$density_control$target_relative_mcse, target)
    }
  }
})

test_that("E2-3 marginal specs use the actual fitted contribution target", {

  for (multiplier in c(2, 1, -2)) {
    object <- .review_h_contribution(multiplier)
    means <- marginal_means(object, bf = FALSE, density_method = "KDE")
    spec <- .iwmde_marginal_means_specs(means, "g", "conditional", "B")[[1L]]
    weights <- .iwmde_linear_weights(spec$weights)
    draws <- as.matrix(coda::as.mcmc(object$fit))
    projection <- as.numeric(draws[, names(weights), drop = FALSE] %*% weights)
    expect_equal(projection, as.numeric(means$inference$conditional$mu_g$B), tolerance = 1e-14)
    expect_identical(weights, c(mu_intercept = 1, mu_g = multiplier))
    combination <- .hypothesis_plan_marginal_means_build(means, "1 * g[B] = 0.1")$targets[[1L]]$spec
    spec$weights <- weights
    expect_identical(weights[order(names(weights))], combination$weights[order(names(combination$weights))])
  }
  object <- .review_h_contribution("tau")
  means <- marginal_means(object, bf = FALSE, density_method = "KDE")
  expect_error(.iwmde_marginal_means_specs(means, "g", "conditional", "B"),
    class = "RoBMA_density_method_original_scale")
})

test_that("E2-4 method switches resolve fresh controls and same methods inherit", {

  object <- .review_h_contribution(1)
  for (method in c("qCMDE", "IWMDE")) {
    means <- marginal_means(object, bf = FALSE, parameter = "g", density_method = method,
      density_control = list(n_points = 80L, samples = 120L))
    other <- setdiff(c("qCMDE", "IWMDE"), method)
    fresh <- .hypothesis_marginal_means_density_control(means, other, NULL)
    expect_identical(fresh$samples, .density_control_normalize(other, purpose = "ordinate")$samples)
    expect_identical(.hypothesis_marginal_means_density_control(means, method, NULL)$samples, 120L)
    expect_identical(.hypothesis_marginal_means_density_control(means, other, list(samples = 321L))$samples, 321L)
    means$ordinate_settings$samples <- 119L
    expect_identical(.hypothesis_marginal_means_density_control(means, method, NULL)$samples, 119L)
    means$ordinate_settings <- NULL
    expect_identical(.hypothesis_marginal_means_density_control(means, method, NULL)$samples, 120L)
  }
})

test_that("E3-3 plots forward stored curve controls separately from point settings", {

  means <- marginal_means(.review_h_contribution(1), bf = FALSE, parameter = "g",
    density_method = "qCMDE", density_control = list(n_points = 80L, samples = 120L))
  means$ordinate_settings$n_points <- 77L
  seen <- NULL
  testthat::local_mocked_bindings(.marginal_means_plot_prepare_iwmde = function(
      x, selected, type, density_method, density_control, force) {
    seen <<- list(control = density_control, force = force)
    stop("Captured plot controls.", call. = FALSE)
  }, .package = "RoBMA")
  expect_error(plot(means, parameter = "g"), "Captured plot controls.", fixed = TRUE)
  expect_identical(seen$control$n_points, 80L)
  expect_false(seen$force)
  expect_error(plot(means, parameter = "g", density_control = list(n_points = 80L)),
    "Captured plot controls.", fixed = TRUE)
  expect_true(seen$force)
  expect_error(plot(means, parameter = "g", density_method = "IWMDE"), "Captured plot controls.", fixed = TRUE)
  expect_identical(seen$control$samples, 1000L)
})

test_that("E2-5 transformed source ordinates keep their declared source prior", {

  object <- .review_h_random("fisher")
  warnings <- character()
  result <- withCallingHandlers(hypothesis(object, "`(mu) rho` = 0.8", density_method = "qCMDE", seed = 19),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
  expect_s3_class(result, "BayesTools_hypothesis_BF")
  expect_false(any(grepl("target prior density tends to zero", warnings, fixed = TRUE)))
  diagnostics <- density_diagnostics(result)
  expect_false(any(grepl("target prior density tends to zero", diagnostics$warnings, fixed = TRUE)))
  target <- .brma_random_parameter_density_target(object, "(mu) rho", operation = "point hypotheses")
  context <- .iwmde_context(object)
  spec <- .iwmde_prepare_prior_ordinates(context, target$parameter, target$parameter_spec,
    BayesTools::parameter_transform_inverse(.8, target$display_transform))
  expect_identical(spec$prior_ordinates[[1L]]$behavior, "regular")
})

test_that("E2-6 interpretation names the actual pooled common quantity", {

  object <- .review_h_random("common")
  pooled <- pooled_heterogeneity(object)
  expect_identical(colnames(pooled), "tau_common")
  result <- interpret(object, scope = "estimates", central = "mean", probs = c(.025, .975), digits = 10)
  records <- attr(result, "records")
  records <- records[grepl("heterogeneity", records$record_id, fixed = TRUE), , drop = FALSE]
  expect_identical(records$parameter, "tau_common")
  expect_equal(records$central_value, mean(pooled), tolerance = 1e-14)
  expect_equal(c(records$lower_value, records$upper_value), unname(stats::quantile(pooled, c(.025, .975))), tolerance = 1e-14)
  expect_true(any(grepl("mean tau_common", result, fixed = TRUE)))
})

test_that("E2-7 genuine precomputed numerical refusals retain exact row diagnostics", {

  object <- .review_h_ordinate_refusal()
  means <- marginal_means(object, bf = TRUE, parameter = "g", density_method = "qCMDE",
    n_samples = 256, density_control = list(samples = 100L))
  cells <- means$inference$conditional$mu_g
  bf <- means$inference$inference$mu_g
  records <- attr(cells$B, "ordinate_failures", exact = TRUE)
  expect_identical(records$numerical_status, "underflow")
  expect_identical(records$requested_value, 0)
  expect_true(is.finite(as.numeric(bf$A)))
  expect_true(is.na(bf$B))
  condition <- attr(bf$B, "numerical_diagnostics", exact = TRUE)
  expect_s3_class(condition, "RoBMA_density_ordinate_error")
  expect_identical(attr(bf$B, "posterior_density_source", exact = TRUE), "numerical_unavailable")
  if (inherits(condition, "condition")) expect_identical(condition$density_diagnostics$log_ordinate, records$log_ordinate)
  stale_request <- .marginal_means_ordinate_request_provenance(means, "mu_g", 0, "IWMDE")
  stale <- .marginal_means_iwmde_bf(cells, 0, stale_request, "IWMDE")
  expect_false(grepl("underflow", attr(stale$B, "warnings"), fixed = TRUE))
  expect_null(attr(stale$B, "numerical_diagnostics", exact = TRUE))
  other_request <- .marginal_means_ordinate_request_provenance(means, "mu_g", .1, "qCMDE")
  other <- .marginal_means_iwmde_bf(cells, .1, other_request, "qCMDE")
  expect_false(grepl("underflow", attr(other$B, "warnings"), fixed = TRUE))
  expect_null(attr(other$B, "numerical_diagnostics", exact = TRUE))
  changed <- means
  changed$ordinate_settings$samples <- 101L
  changed_request <- .marginal_means_ordinate_request_provenance(changed, "mu_g", 0, "qCMDE")
  changed_bf <- .marginal_means_iwmde_bf(cells, 0, changed_request, "qCMDE")
  expect_false(grepl("underflow", attr(changed_bf$B, "warnings"), fixed = TRUE))
  expect_null(attr(changed_bf$B, "numerical_diagnostics", exact = TRUE))
  legacy <- cells
  attr(attr(legacy$B, "ordinate_failures"), "iwmde_provenance") <- NULL
  expected <- .marginal_means_ordinate_request_provenance(means, "mu_g", 0, "qCMDE")
  legacy_bf <- .marginal_means_iwmde_bf(legacy, 0, expected, "qCMDE")
  expect_null(attr(legacy_bf$B, "numerical_diagnostics", exact = TRUE))
  expect_false(grepl("underflow", attr(legacy_bf$B, "warnings"), fixed = TRUE))
  all_failed <- marginal_means(object, null_hypothesis = 100, bf = TRUE, parameter = "g",
    density_method = "qCMDE", n_samples = 256, density_control = list(samples = 100L))
  for (item in all_failed$inference$inference$mu_g) {
    expect_true(is.na(item))
    expect_s3_class(attr(item, "numerical_diagnostics"), "RoBMA_density_ordinate_error")
  }
  table <- summary(means, bf = TRUE)
  diagnostics <- attr(table, "numerical_diagnostics", exact = TRUE)
  expect_true(any(vapply(diagnostics, inherits, logical(1), "RoBMA_density_ordinate_error")))
  for (frame in list(as.data.frame(table), data.frame(table))) {
    expect_identical(frame$parameter, c("g[A]", "g[B]"))
    expect_equal(as.numeric(frame$inclusion_BF), c(as.numeric(bf$A), NA_real_))
    expect_null(attr(frame$inclusion_BF, "numerical_diagnostics", exact = TRUE))
  }
})

test_that("E3-1 scale-bound random correlation summaries remain available", {

  for (structure in c("cs", "ar1")) {
    object <- .review_h_random("scale", structure)
    catalog <- BayesTools::parameter_catalog(object$fit)$quantities
    expect_true(any(catalog$role == "random_correlation"))
    direct <- BayesTools::JAGS_estimates_table(object$fit, random_effects_summary = "standard",
      keep_parameters = "random", remove_diagnostics = TRUE, diagnostic_columns = "none",
      simplify_names = TRUE, transform_scaled = TRUE, formula_prefix = FALSE,
      random_effects_metadata = TRUE, remove_spike_0 = FALSE, remove_inclusion = FALSE)
    result <- summary(object, include_mcmc_diagnostics = FALSE)
    expect_true(.summary_random_components_enabled(object))
    expect_true(is.data.frame(result$estimates_random) && nrow(result$estimates_random) > 0L)
    if (is.data.frame(result$estimates_random) && nrow(result$estimates_random) > 0L)
      expect_equal(result$estimates_random$Mean, direct$Mean, tolerance = 1e-14)
    for (frame in list(as.data.frame(result), data.frame(result))) {
      expect_true(any(grepl("rho", frame$parameter, fixed = TRUE)))
    }
    expect_true(any(grepl("rho", capture.output(print(result)), fixed = TRUE)))
  }
})

test_that("E3-2 inclusion warnings follow actual relabel and row selection", {

  original <- c("(mu) intercept", "(mu) inclusion(study)")
  table <- data.frame(prior_prob = .5, post_prob = c(.75, .8),
    inclusion_BF = c(3, 4), row.names = original)
  table$inclusion_BF <- BayesTools::format_BF(table$inclusion_BF, inclusion = TRUE)
  class(table) <- c("BayesTools_table", "data.frame")
  attr(table, "type") <- c("prior_prob", "post_prob", "inclusion_BF")
  attr(table, "parameters") <- original
  attr(table, "rownames") <- TRUE
  attr(table, "logBF") <- FALSE
  attr(table, "BF01") <- FALSE
  attr(table, "warnings") <- stats::setNames(c("Synthetic effect marker E.", "Synthetic random marker R.",
    "Synthetic global marker G."), c(original, ""))
  core <- .summary.inclusion_subtable(table, 1L, "Effect", "Component Inclusion")
  random <- .summary.inclusion_subtable(table, 2L, "study: tau", "Random-Effect Inclusion")
  result <- structure(list(inclusion_components = core, inclusion_random = random,
    name = "Synthetic warning transport"), class = "summary.brma")
  merged <- .summary_brma_component_inclusion(result)
  warnings <- attr(merged, "warnings", exact = TRUE)
  expect_true(all(c("Effect", "Heterogeneity: study: tau") %in% names(warnings)))
  expect_true(any(grepl("marker G", warnings, fixed = TRUE)))
  selected <- merged["Effect", , drop = FALSE]
  expect_true(any(grepl("marker E", attr(selected, "warnings"), fixed = TRUE)))
  expect_false(any(grepl("marker R", attr(selected, "warnings"), fixed = TRUE)))
  expect_equal(as.numeric(merged$inclusion_BF), c(3, 4))
  for (field in c("mods", "scale")) attr(result, field) <- FALSE
  attr(result, "random") <- TRUE
  printed <- capture.output(print(result))
  expect_equal(sum(grepl("marker E", printed, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("marker R", printed, fixed = TRUE)), 1L)
})

test_that("E1-1 multiple declared atomic locations are not one fixed scalar", {

  means <- marginal_means(.review_h_atomic(), bf = FALSE, density_method = "KDE", n_samples = 128)
  cell <- means$inference$averaged$mu_intercept$intercept
  expect_false(.hypothesis_plan_draws_fixed(cell))
  result <- hypothesis(means, "intercept > 0", density_method = "KDE", seed = 19)
  direct <- BayesTools::hypothesis_BF(means$inference$averaged$mu_intercept,
    hypothesis = "intercept > 0", parameter = "intercept", density_method = "KDE")
  expect_equal(as.numeric(result$BF), as.numeric(direct$BF), tolerance = 1e-14)
  expect_equal(as.numeric(result$BF), 1.5, tolerance = 1e-14)
  duplicate <- BayesTools::posterior_atom_attribute(data.frame(x = c(.2, .2), mass = c(.4, .6)))
  fixed <- with_draw_metadata(rep(.2, 3), atoms = duplicate)
  expect_true(.hypothesis_plan_draws_fixed(fixed))
  expect_error(BayesTools::posterior_atom_attribute(data.frame(x = c(-.2, .2), mass = c(0, 1))))
  expect_false(.hypothesis_plan_draws_fixed(numeric()))
})

test_that("E1-2 indexed scalar labels do not acquire factor suffixes", {

  dat <- data.frame(yi = 0, study = factor(rep(c("a", "b"), each = 2)),
    out = factor(rep(c("x", "y"), 2)))
  setup <- review_compile_object(brma.mv(yi = yi, V = diag(4), data = dat, measure = "GEN",
    random = ~us(0 + out | study), prior_effect = BayesTools::prior("point", list(0)),
    prior_heterogeneity = BayesTools::prior_random(sd = BayesTools::prior("gamma", list(2, 2)),
      study = BayesTools::random_block(contrasts = c(out = "independent"))),
    only_priors = TRUE, silent = TRUE))
  term <- setup$object$formula_design$mu$random_effects[[1L]]
  columns <- unique(c(review_prior_columns(setup$compiled$mu$prior_list), term$correlation$primitive_names))
  samples <- matrix(.5, 100, length(columns), dimnames = list(NULL, columns))
  samples[, "mu_intercept"] <- 0
  for (column in setdiff(columns, "mu_intercept")) samples[, column] <- seq(.2, .8, length.out = 100)
  object <- review_attach_replay(setup, samples)
  for (parameter in c("tau(out[x])", "rho(out[x],out[y])")) {
    statement <- paste0("`", parameter, "` > 0.5")
    result <- hypothesis(object, rep(statement, 2), component = "random", density_method = "KDE", seed = 19)
    expect_identical(rownames(result), paste0(parameter, " (", 1:2, ")"))
    expect_identical(as.data.frame(result)$parameter, rownames(result))
  }
})

test_that("E1-3 whole-term and explicit-level modes form compatible groups", {

  object <- review_factor_replay(c("A", "B", "C"), "meandif", independent = TRUE, n = 240L)
  means <- marginal_means(object, bf = FALSE, density_method = "KDE")
  for (input in list(object, means)) {
    for (operator in c(">", "=")) {
      statements <- paste(c("g[B]", "g"), operator, "0")
      run <- function(statements) hypothesis(input, statements, density_method = "KDE", seed = 19)
      for (order in list(1:2, 2:1)) {
        combined <- run(statements[order])
        separate <- lapply(statements[order], run)
        expect_identical(nrow(combined), 4L)
        expect_equal(attr(combined, "raw_log_BF"), unlist(lapply(separate, attr, which = "raw_log_BF"), use.names = FALSE), tolerance = 1e-12)
      }
    }
  }
})

test_that("E1-4 model-average request limits count actual expanded targets", {

  means <- marginal_means(.review_h_empty_condition(), bf = FALSE, density_method = "KDE", n_samples = 128)
  plan <- .hypothesis_plan_marginal_means(means, "x = 0")
  expect_gt(length(plan$targets), 1L)
  refusal <- .hypothesis_plan_marginal_means_request_refusal(list(plan), means)
  expect_true(!is.null(refusal))
  expect_error(hypothesis(means, "x = 0", density_method = "KDE", seed = 19),
    "spanning multiple marginal-means levels", class = "RoBMA_hypothesis_statement")
  one <- hypothesis(means, "intercept = 0", density_method = "KDE", seed = 19)
  alias <- hypothesis(means, "mu_intercept = 0", density_method = "KDE", seed = 19)
  expect_equal(attr(one, "raw_log_BF"), attr(alias, "raw_log_BF"), tolerance = 1e-14)
  repeated <- lapply(c("intercept = 0", "mu_intercept = 0.1"), function(statement)
    .hypothesis_plan_marginal_means(means, statement))
  expect_null(.hypothesis_plan_marginal_means_request_refusal(repeated, means))
})

test_that("E2-7 actual structural point-mass conditions keep their source identity", {

  means <- marginal_means(.review_h_atomic(), bf = FALSE, density_method = "KDE", n_samples = 128)
  posterior <- means$inference$conditional$mu_intercept
  result <- .marginal_means_inclusion_bf(posterior, -.2, TRUE)$intercept
  original <- attr(result, "numerical_diagnostics")
  expect_s3_class(original, "BayesTools_posterior_point_mass_at_null")
  expect_true(is.na(result))
  expect_true(is.character(conditionMessage(original)))
  expect_identical(attr(result, "posterior_density_source"), "null_point_mass")
  expect_error(.marginal_means_row_bf(stop("Programming failure.", call. = FALSE)), "Programming failure.", fixed = TRUE)
})

test_that("synthetic accepted diagnostics still require declared lattice accuracy", {

  diagnostics <- list(estimator = "q_grid_cmde", ordinate_relative_change = 0,
    lattice_resolution_relative_error = 0)
  expect_null(.iwmde_diagnostics_bf_failure_reason(diagnostics))
  diagnostics$lattice_resolution_relative_error <- NULL
  expect_match(.iwmde_diagnostics_bf_failure_reason(diagnostics), "normalization diagnostics are unavailable", fixed = TRUE)
})


# Append to test-00-review-hypothesis-summary.R, whose genuine replay helpers
# this block reuses.
test_that("E1-c1 default named and positional hypotheses share the actual carrier", {

  object <- review_factor_replay(c("A", "B", "C"), "meandif",
    independent = TRUE, distribution = "mnormal", n = 128L)
  posterior <- BayesTools::marginal_posterior(
    .brma_as_mixed_posteriors(object, parameters = "mu_intercept"),
    parameter = "mu_intercept", prior_samples = TRUE, use_formula = FALSE,
    n_samples = 128L
  )
  original <- posterior
  named <- hypothesis(object = posterior, hypothesis = "mu_intercept > 0",
    parameter = "mu_intercept", seed = 17)
  positional <- hypothesis(posterior, "mu_intercept > 0", parameter = "mu_intercept", seed = 17)
  expect_identical(positional, named)
  expect_identical(posterior, original)
})

test_that("E2-c2 warning ownership distinguishes exact references and expressions", {

  statements <- c("tau = 0", "tau_total = 0", "2 * tau = 0",
    "g[A] = 0", "g[AA] = 0")
  ast <- BayesTools::hypothesis_parse(statements)
  table <- data.frame(
    Alternative = statements,
    Null = c("tau != 0", "tau_total != 0", "2 * tau != 0",
      "g[A] != 0", "g[AA] != 0"),
    row.names = c("tau", "tau_total", "scaled_tau", "g[A]", "g[AA]"),
    check.names = FALSE
  )
  attr(table, "hypothesis_ast") <- ast
  record <- function(parameter, label = parameter, level = NA_character_, value = 0) {
    data.frame(parameter = parameter, label = label, level = level, value = value,
      warning = "Declared ordinate warning.", stringsAsFactors = FALSE)
  }
  expect_identical(.hypothesis_brma_warning_rows(table, record("tau")), "tau")
  expect_identical(.hypothesis_brma_warning_rows(table, record("tau_total")), "tau_total")
  expect_identical(.hypothesis_brma_warning_rows(table, record("2 * tau")), "scaled_tau")
  expect_identical(.hypothesis_brma_warning_rows(table, record("g", "g[A]", "A")), "g[A]")
  expect_identical(.hypothesis_brma_warning_rows(table, record("tau", value = .1)), character())
  expanded <- table[rep(4L, 4L), , drop = FALSE]
  rownames(expanded) <- c("g[A] (1)", "g[B] (1)", "g[A] (2)", "g[B] (2)")
  expanded$Alternative <- "g = 0"
  expanded$Null <- "g != 0"
  attr(expanded, "hypothesis_ast") <- BayesTools::hypothesis_parse(c("g = 0", "g = 0"))
  expect_identical(.hypothesis_brma_warning_rows(expanded, record("g", "g[A]", "A")),
    c("g[A] (1)", "g[A] (2)"))
})

test_that("E1-c1 unsupported random comparisons refuse before SD-pair work", {

  testthat::local_mocked_bindings(
    .brma_random_parameter_sd_pair = function(...) stop("SD-pair lookup was evaluated."),
    .package = "RoBMA"
  )
  for (statement in c("tau2 != 0 vs tau2 > 0", "tau2 = -1 vs tau2 = 0",
    "tau2 = 0 vs tau2 = .1")) {
    ast <- BayesTools::hypothesis_parse(statement)
    plan <- .hypothesis_plan_new(ast, ast, "tau2", "tau2", "random")
    result <- .hypothesis_plan_random_evaluation(plan, object = NULL, cache = NULL,
      selected = list(spec = list(label = "tau2")))
    expect_identical(result$refusal$class, "RoBMA_hypothesis_statement")
    expect_match(result$refusal$reason, "Unsupported hypothesis comparison.", fixed = TRUE)
  }
})

test_that("E1-c1 a variance refusal preserves SD evaluation evidence", {

  diagnostics <- data.frame(requested_value = .3, evaluation_value = .3,
    estimator = "q_grid_cmde", density_method = "qCMDE", status = "unavailable",
    numerical_status = "underflow", failure_reason = "Ordinary-scale ordinate underflow.")
  original <- .iwmde_ordinate_unavailable_condition("Original SD message.",
    list(ordinate_failures = diagnostics))
  testthat::local_mocked_bindings(
    .hypothesis_plan_random_evaluate = function(...) stop(original),
    .package = "RoBMA"
  )
  plan <- list(label = "50% study: tau2",
    evaluation = list(parameter = "tau", samples = list(tau = c(.2, .3))))
  result <- tryCatch(.hypothesis_plan_random_sd_points(plan,
    BayesTools::hypothesis_parse("tau2 = .09"), values = .09,
    density_method = "qCMDE"), error = identity)
  expect_s3_class(result, "RoBMA_density_ordinate_error")
  expect_match(conditionMessage(result), "50% study: tau2 = 0.09", fixed = TRUE)
  expect_identical(result$density_diagnostics$requested_value, .09)
  expect_identical(result$density_diagnostics$evaluation_value, diagnostics$evaluation_value)
  expect_identical(result$density_diagnostics$estimator, diagnostics$estimator)
  expect_identical(result$density_diagnostics$failure_reason, diagnostics$failure_reason)
  expect_identical(result$density_diagnostics$evaluation_scale, "sd")
  expect_identical(result$density_diagnostics$requested_scale, "variance")
})

test_that("E3-c1 combined keyed notes follow exact colliding row labels", {

  make_table <- function(value) {
    table <- data.frame(Mean = value, row.names = "tau")
    class(table) <- c("BayesTools_table", "data.frame")
    attr(table, "type") <- "estimate"
    attr(table, "parameters") <- "tau-source"
    attr(table, "footnotes") <- stats::setNames(
      c("tau: Defined-draw explanation.", "Ordinary global explanation."),
      c("tau-source", ""))
    table
  }
  output <- .summary_brma_combine_tables(list(make_table(.1), make_table(.2)), "Estimates")
  notes <- attr(output, "footnotes", exact = TRUE)
  expect_identical(rownames(output), c("tau", "tau.1"))
  expect_identical(unname(notes[names(notes) == "tau"]), "tau: Defined-draw explanation.")
  expect_identical(unname(notes[names(notes) == "tau.1"]), "tau.1: Defined-draw explanation.")
  expect_equal(sum(notes == "Ordinary global explanation."), 1L)
  expect_identical(output$Mean, c(.1, .2))
})

test_that("E3-c2 unused summary arguments warn and lines rejects explicit add", {

  object <- .review_h_contribution()
  expect_warning(summary(object, include_mcmc_diagnostics = FALSE,
    unused_cleanup_argument = TRUE),
    "Unused argument in summary.brma(): 'unused_cleanup_argument'", fixed = TRUE)
  models <- .review_h_empty_condition()
  expect_warning(summary_models(models, include_mcmc_diagnostics = FALSE,
    unused_cleanup_argument = TRUE),
    "Unused argument in summary_models(): 'unused_cleanup_argument'", fixed = TRUE)
  expect_error(lines(object, add = FALSE),
    "'lines()' always adds to an existing plot. Omit 'add'; use 'plot()' to create a new plot.",
    fixed = TRUE)
})

test_that("E3-c2 model tables support both plain data-frame coercions", {

  object <- .review_h_empty_condition()
  for (type in c("marginal", "individual")) {
    result <- summary_models(object, type = type, include_mcmc_diagnostics = FALSE)
    frame <- as.data.frame(result)
    expect_identical(data.frame(result, check.names = FALSE), frame)
    expect_identical(names(data.frame(result)), make.names(names(frame), unique = TRUE))
    expect_identical(class(frame), "data.frame")
    expect_identical(names(frame)[1:2], c("component", "parameter"))
    expect_true(all(vapply(frame[c("prior_prob", "post_prob", "inclusion_BF")], is.numeric, logical(1L))))
    if (type == "marginal") {
      expected <- unlist(lapply(result$marginal, `[[`, "post_prob"), use.names = FALSE)
      expect_identical(frame$post_prob, expected)
      expect_identical(unique(frame$component), names(result$marginal))
    } else {
      expect_identical(frame$post_prob, result$individual$post_prob)
      expect_identical(unique(frame$component), "individual models")
    }
  }
})

test_that("E3-c3 a hidden BF summary skips refresh and retains its source evidence", {

  means <- marginal_means(.review_h_contribution(), bf = FALSE,
    density_method = "KDE", n_samples = 128L)
  original <- means
  testthat::local_mocked_bindings(
    .marginal_means_refresh_iwmde_bf = function(...) stop("Hidden BF refresh was evaluated."),
    .package = "RoBMA"
  )
  result <- summary(means, bf = FALSE)
  expect_s3_class(result, "summary.marginal_means.brma")
  expect_false(any(grepl("BF", names(result), fixed = TRUE)))
  expect_identical(means, original)
})

test_that("E3-c3 precomputed marginal means skip the preliminary KDE BF", {

  count <- 0L
  compute <- .marginal_means_inclusion_bf_one
  testthat::local_mocked_bindings(
    .marginal_means_inclusion_bf_one = function(posterior, null_hypothesis) {
      count <<- count + 1L
      compute(posterior, null_hypothesis)
    },
    .package = "RoBMA"
  )
  result <- marginal_means(.review_h_contribution(), bf = TRUE,
    density_method = "qCMDE", n_samples = 128L)
  expect_s3_class(result, "marginal_means.brma")
  expect_identical(count, 0L)
  expect_identical(result$density_method, "qCMDE")
  expect_true(is.list(result$ordinate_diagnostics))
})

test_that("A5-5 genuine BMA gate eligibility controls undefined-draw notes", {

  fixture <- review_bma_gate_replay()
  object <- fixture$object
  catalog <- BayesTools::parameter_catalog(object$fit)$quantities
  keep <- startsWith(catalog$role, "random_") & !catalog$internal &
    catalog$quantity %in% .brma_random_parameter_supported_quantities()
  specs <- .brma_random_parameter_specs(catalog[keep, , drop = FALSE])
  target <- specs$parameter[specs$quantity == "var_prop" & specs$random_component == "study"]
  expect_length(target, 1L)
  selected <- .brma_random_parameter_select(object, target)
  states <- BayesTools::parameter_gate_states(object$fit, selected$entry$selection)
  expect_identical(states$event, c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(states$defined, c(FALSE, TRUE, TRUE, TRUE))
  conditional <- hypothesis(object, paste0("`", target, "` > 0.1"), conditional = TRUE,
    seed = 17, n_samples = 128L)
  averaged <- hypothesis(object, paste0("`", target, "` > 0.1"), conditional = FALSE,
    seed = 17, n_samples = 128L)
  expect_false(any(grepl("posterior draws where the quantity is defined",
    attr(conditional, "footnotes", exact = TRUE), fixed = TRUE)))
  expect_match(attr(averaged, "footnotes", exact = TRUE),
    "computed from 3 of 4 posterior draws where the quantity is defined.", fixed = TRUE)
  expect_identical(as.numeric(conditional$BF), Inf)
  expect_equal(as.numeric(averaged$BF), 1.1578947368421058, tolerance = 1e-14)
  parent <- review_heterogeneity_parent_replay("total_variance")$object
  parent_catalog <- BayesTools::parameter_catalog(parent$fit)$quantities
  parent_keep <- startsWith(parent_catalog$role, "random_") & !parent_catalog$internal &
    parent_catalog$quantity %in% .brma_random_parameter_supported_quantities()
  parent_specs <- .brma_random_parameter_specs(parent_catalog[parent_keep, , drop = FALSE])
  child <- parent_specs$parameter[parent_specs$quantity == "var_prop" &
    parent_specs$random_component == "time[1]"]
  expect_length(child, 1L)
  undefined <- hypothesis(parent, paste0("`", child, "` > 0.2"),
    seed = 17, n_samples = 128L)
  expect_match(attr(undefined, "footnotes", exact = TRUE),
    "computed from 2 of 4 posterior draws", fixed = TRUE)
  expect_equal(attr(undefined, "raw_log_BF"), -Inf)
})

test_that("E1-c1 genuine percent labels retain the point-mass refusal", {

  for (label in c("study", "50% study")) {
    fixture <- review_bma_gate_replay(label)
    catalog <- BayesTools::parameter_catalog(fixture$object$fit)$quantities
    keep <- startsWith(catalog$role, "random_") & !catalog$internal &
      catalog$quantity %in% .brma_random_parameter_supported_quantities()
    specs <- .brma_random_parameter_specs(catalog[keep, , drop = FALSE])
    target <- specs$parameter[specs$quantity == "sd" & specs$owner_name == label]
    expect_length(target, 1L)
    selected <- .brma_random_parameter_select(fixture$object, target)
    expect_identical(as.numeric(selected$samples), c(0, .5, 0, .5))
    result <- tryCatch(hypothesis(fixture$object, paste0("`", target, "` = 0")), error = identity)
    expect_s3_class(result, "BayesTools_point_mass_at_null")
    expected <- paste0("Point-null Bayes factors at 0 are unavailable for random-effect quantity '",
      label, ": tau' because its prior has a point mass at 0. Use a region or directional hypothesis; ",
      "the Component Inclusion table from 'summary(object)' or 'summary_models(object)' compares the exclusion and inclusion of gated components.")
    expect_identical(conditionMessage(result), expected)
  }
})

test_that("E3-c4 fixed-zero formulas retain actual mods and location owners", {

  for (owner in c("mods", "location")) {
    object <- if (owner == "mods") .review_h_ordinate_refusal() else review_fixed_zero_location_replay()
    original <- object
    source <- .fitted_formula_source("mu", object$data)
    design <- .fitted_formula_design(object, "mu", required = TRUE)
    formula <- .fitted_formula_evaluable(object, design, source)
    terms <- .marginal_means_terms(design)
    parameters <- BayesTools::JAGS_parameter_names(terms, formula_parameter = "mu")
    dropped <- .marginal_means_drop_fixed_zero_intercept(object, formula, terms, parameters)
    expect_identical(source, owner)
    expect_length(dropped$formula, 2L)
    expect_identical(environment(dropped$formula), environment(formula))
    expect_identical(dropped$terms, "g")
    expect_identical(dropped$parameters, "mu_g")
    means <- marginal_means(object, bf = FALSE, density_method = "KDE")
    expect_identical(names(means$inference$averaged$mu_g), c("A", "B"))
    draws <- as.matrix(coda::as.mcmc(object$fit))
    for (level in c("A", "B")) {
      cell <- means$inference$averaged$mu_g[[level]]
      weights <- .iwmde_linear_weights(BayesTools::posterior_metadata(cell, "linear_weights"))
      expect_equal(as.numeric(cell), as.numeric(draws[, names(weights), drop = FALSE] %*% weights),
        tolerance = 1e-14)
    }
    expect_identical(object, original)
    expect_identical(.fitted_formula_design(object, "mu", required = TRUE), design)
  }
})

test_that("E2-c2 actual scaled-expression warnings retain their rewritten identity", {

  records <- list()
  append <- .hypothesis_brma_append_iwmde_warnings
  testthat::local_mocked_bindings(.hypothesis_brma_append_iwmde_warnings = function(table, posterior, ...) {
    records[[length(records) + 1L]] <<- list(
      record = .iwmde_collect_posterior_ordinate_warning_records(posterior),
      table = table, ast = attr(table, "hypothesis_ast"))
    append(table, posterior, ...)
  }, .package = "RoBMA")
  fixture <- .review_h_contribution(2)
  result <- hypothesis(fixture, c("2 * g[B] = 0.1", "2 * g[B] = 0.2"),
    density_method = "qCMDE", n_samples = 128L, seed = 19)
  expect_length(records, 1L)
  record <- records[[1L]]$record
  expect_identical(record$parameter, rep(".BayesTools_linear_target", 2L))
  expect_identical(record$label, record$parameter)
  expect_identical(record$value, c(.1, .2))
  expect_true(all(is.na(record$level)))
  # Savage-Dickey values of this conjugate fixture. With tau = .1 the data
  # variance is .5^2 + .1^2, and the posterior of g[B] given the intercept of a
  # retained row is normal with precision 1 + sum((2 * z)^2) / variance and mean
  # sum(2 * z * (yi - intercept)) / (variance * precision), where z indicates
  # level B and 2 is the coefficient multiplier. The qCMDE ordinate of
  # 2 * g[B] is the mean of these exact conditional ordinates over the retained
  # rows; its prior ordinate is that of normal(0, 2). The allowance is the
  # numerical error reported by the diagnostics plus rounding; the sampling
  # warnings above concern the replayed rows, not this comparison.
  intercept <- as.matrix(fixture$fit$mcmc)[, "mu_intercept"]
  yi        <- c(-.2, .3, 0, .4, .1, .2)
  z         <- rep(c(0, 1), 3L)
  variance  <- .5^2 + .1^2
  precision <- 1 + sum((2 * z)^2) / variance
  conditional_mean <- vapply(intercept, function(b) {
    sum(2 * z * (yi - b)) / (variance * precision)
  }, numeric(1))
  expected <- vapply(c(.1, .2), function(target) {
    stats::dnorm(target, 0, 2) /
      mean(stats::dnorm(target, 2 * conditional_mean, 2 / sqrt(precision)))
  }, numeric(1))
  diagnostics <- density_diagnostics(result)
  reported <- diagnostics$normalization_relative_error + diagnostics$stability_relative_error +
    diagnostics$lattice_resolution_relative_error + diagnostics$truncation_ordinate_bound
  expect_true(all(abs(as.numeric(result$BF) / expected - 1) <=
                    reported + 64 * .Machine$double.eps))
  warnings <- attr(result, "warnings", exact = TRUE)
  expect_identical(names(warnings), rownames(result))
  expect_identical(unname(warnings), record$warning)
  for (i in 1:2) {
    expect_identical(.hypothesis_brma_warning_rows(records[[1L]]$table, record[i, , drop = FALSE]),
      rownames(records[[1L]]$table)[i])
  }
})

test_that("marginal-means requests count the genuine singleton intercept target", {

  means <- marginal_means(.review_h_empty_condition(), bf = FALSE,
    density_method = "KDE", n_samples = 128L)
  expect_identical(names(means$inference$averaged$mu_intercept), "intercept")
  expect_identical(names(means$inference$conditional$mu_intercept), "intercept")
  original <- means$inference
  level <- names(means$inference$averaged$mu_x)[[1L]]
  slope <- paste0("`x[", level, "]` = 0.1")
  refusal <- function(statements) {
    plans <- lapply(statements, function(statement)
      .hypothesis_plan_marginal_means(means, statement))
    .hypothesis_plan_marginal_means_request_refusal(plans, means)
  }
  reason <- paste0("Point hypotheses spanning multiple marginal-means levels are not ",
    "supported because the levels need not share one conditioning event.")
  other <- paste0("`x[", tail(names(means$inference$averaged$mu_x), 1L), "]` = 0.1")
  for (statements in list(c("intercept = 0.2", slope), c(slope, "intercept = 0.2"),
                           c(slope, other))) {
    expect_identical(refusal(statements)$reason, reason)
    expect_identical(refusal(statements)$class, "RoBMA_hypothesis_statement")
    expect_error(hypothesis(means, statements, seed = 913),
      "spanning multiple marginal-means levels", class = "RoBMA_hypothesis_statement")
  }
  repeated <- c("intercept = 0.2", "mu_intercept = 0.3")
  expect_null(refusal("intercept = 0.2"))
  expect_null(refusal(repeated))
  result <- hypothesis(means, repeated, seed = 913)
  expect_identical(rownames(result),
    c("mu_intercept[intercept] (1)", "mu_intercept[intercept] (2)"))
  expect_identical(refusal(c("intercept = 0.2", "x > 0"))$class,
    "RoBMA_hypothesis_statement")
  expect_null(refusal(c("intercept > 0", "x > 0")))
  expect_identical(means$inference, original)

  single <- marginal_means(review_factor_replay(c("A", "B"), "treatment", distribution = "normal"),
    bf = FALSE, density_method = "KDE", n_samples = 128L)
  statements <- c("intercept = 0.2", "g[A] = 0.1")
  plans <- lapply(statements, function(statement)
    .hypothesis_plan_marginal_means(single, statement))
  expect_null(.hypothesis_plan_marginal_means_request_refusal(plans, single))
})
