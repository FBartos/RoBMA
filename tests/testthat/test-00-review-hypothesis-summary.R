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
  samples <- cbind(mu_intercept = .2 + .025 * stats::qnorm((seq_len(120) - .5) / 120),
    mu_x = c(rep(0, 60), seq(.01, .15, length.out = 60)),
    mu_x_indicator = rep(c(null, alternative), each = 60), tau = 0)
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
  # A public point reaches estimation; a subsequent diagnostic refusal is kept.
  result <- tryCatch(hypothesis(means, "intercept = 0", density_method = "qCMDE", seed = 19), error = identity)
  expect_false(inherits(result, "error") && grepl("metadata is incomplete", conditionMessage(result), fixed = TRUE))
})

test_that("E2-2 actual attachment retains requested curve precision", {

  object <- review_factor_replay(c("A", "B"), "treatment", independent = TRUE,
    distribution = "normal", n = 240L)
  result <- marginal_means(object, bf = FALSE, parameter = "g", density_method = "qCMDE",
    density_control = list(target_relative_mcse = .1, samples = 120L))
  expect_identical(result$density_settings$target_relative_mcse, .1)
  cells <- result$inference$averaged$mu_g
  for (cell in cells) {
    density <- BayesTools::posterior_metadata(cell, "posterior_density")
    if (!is.null(density)) expect_identical(density$iwmde_provenance$density_control$target_relative_mcse, .1)
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
