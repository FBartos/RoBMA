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
