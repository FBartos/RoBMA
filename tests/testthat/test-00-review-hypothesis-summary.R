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
