source(testthat::test_path("common-functions.R"))

skip_if_missing_fits("RoBMA.mv_exact_product_space")
fit_robma_mv <- load_fit("RoBMA.mv_exact_product_space", validate = FALSE)


.robma_mv_random_gate_names <- function(fit) {

  allocation <- fit[["formula_design"]][["mu"]][["random_allocations"]][[1L]]
  vapply(
    allocation[["inclusion"]],
    `[[`,
    character(1),
    "indicator_name"
  )
}


test_that("RoBMA.mv summaries average every product-space component", {

  samples <- .get_posterior_samples(fit_robma_mv[["fit"]])
  out <- summary(
    fit_robma_mv,
    conditional              = TRUE,
    include_mcmc_diagnostics = FALSE
  )
  models <- summary_models(
    fit_robma_mv,
    type                     = "marginal",
    include_mcmc_diagnostics = FALSE
  )[["marginal"]]

  expect_named(
    models,
    c(
      "Effect", "Location: x", "Publication Bias", "Random: study",
      "Random: observation"
    )
  )
  bias_frequency <- tabulate(
    samples[, "bias_indicator"],
    nbins = 4L
  ) / nrow(samples)
  expect_equal(
    unname(models[["Publication Bias"]][["post_prob"]]),
    bias_frequency
  )
  expect_equal(
    out[["inclusion_components"]]["Publication Bias", "post_prob"],
    sum(bias_frequency[-1L])
  )

  gate_names <- .robma_mv_random_gate_names(fit_robma_mv)
  x_scale <- attr(
    fit_robma_mv[["fit"]],
    "formula_scale"
  )[["mu"]][["mu_x"]]
  expected_intercept <-
    samples[, "mu_intercept"] *
      (samples[, "mu_intercept_indicator"] == 2L) -
    x_scale[["mean"]] / x_scale[["sd"]] * samples[, "mu_x"] *
      (samples[, "mu_x_indicator"] == 2L)
  expect_equal(
    unname(out[["inclusion_random"]][["post_prob"]]),
    unname(colMeans(samples[, gate_names, drop = FALSE]))
  )
  expect_equal(
    out[["estimates_mods"]]["intercept", "Mean"],
    mean(expected_intercept),
    tolerance = 5e-4
  )
  expect_equal(
    out[["estimates_bias"]]["PET", "Mean"],
    mean(samples[, "PET"]),
    tolerance = 5e-4
  )
  expect_equal(
    out[["estimates_bias_conditional"]]["PET", "Mean"],
    mean(samples[samples[, "bias_indicator"] == 3L, "PET"]),
    tolerance = 5e-4
  )
  expect_equal(
    out[["estimates_bias_conditional"]]["PEESE", "Mean"],
    mean(samples[samples[, "bias_indicator"] == 4L, "PEESE"]),
    tolerance = 5e-4
  )

  summary_frame <- as.data.frame(out)
  expect_true(all(c(
    "inclusion", "inclusion random", "location", "bias", "random"
  ) %in% summary_frame[["component"]]))
})


test_that("RoBMA.mv model tables enumerate the combined product space", {

  individual <- summary_models(
    fit_robma_mv,
    type                     = "individual",
    include_mcmc_diagnostics = FALSE
  )[["individual"]]

  expect_equal(nrow(individual), 64L)
  expect_equal(sum(individual[["prior_prob"]]), 1, tolerance = 1e-12)
  expect_equal(sum(individual[["post_prob"]]), 1, tolerance = 1e-12)
  expect_setequal(
    individual[["Publication Bias"]],
    c("None", "omega[one-sided: .025]", "PET", "PEESE")
  )
  expect_setequal(individual[["Random: study"]], c("Excluded", "Included"))
  expect_setequal(
    individual[["Random: observation"]],
    c("Excluded", "Included")
  )
})


test_that("RoBMA.mv predictive and scoring targets remain multivariate", {

  n <- nobs(fit_robma_mv)
  terms <- predict(
    fit_robma_mv,
    type               = "terms",
    conditioning_depth = "marginal"
  )
  estimate <- predict(
    fit_robma_mv,
    type               = "estimate",
    conditioning_depth = "marginal"
  )
  response <- predict(
    fit_robma_mv,
    type               = "response",
    conditioning_depth = "marginal"
  )

  expect_brma_samples_matrix(terms, n, "RoBMA.mv terms")
  expect_brma_samples_matrix(estimate, n, "RoBMA.mv estimate")
  expect_brma_samples_matrix(response, n, "RoBMA.mv response")
  expect_length(fitted(
    fit_robma_mv,
    type               = "estimate",
    conditioning_depth = "estimate"
  ), n)
  expect_brma_samples_matrix(
    pooled_effect(fit_robma_mv),
    1L,
    "RoBMA.mv pooled effect"
  )
  expect_brma_samples_matrix(
    pooled_heterogeneity(fit_robma_mv, component = "total"),
    1L,
    "RoBMA.mv pooled heterogeneity"
  )

  log_likelihood <- log_lik(fit_robma_mv)
  expect_equal(dim(log_likelihood), c(nrow(terms), n))
  expect_s3_class(loo(fit_robma_mv), "loo")
  weights <- loo_weights(fit_robma_mv)
  expect_equal(dim(weights), dim(log_likelihood))
  expect_equal(colSums(weights), rep(1, n), tolerance = 1e-10)

  fit_waic <- fit_robma_mv
  fit_waic[["waic"]] <- NULL
  fit_waic <- suppressWarnings(add_waic(fit_waic))
  expect_s3_class(waic(fit_waic), "waic")
})


test_that("RoBMA.mv formula and ensemble post-processing remains available", {

  n <- nobs(fit_robma_mv)
  expect_named(ranef(fit_robma_mv, simplify = FALSE), c("study", "observation"))
  expect_brma_samples_matrix(
    ranef(fit_robma_mv, component = "total", expand = TRUE),
    n,
    "RoBMA.mv total random effect"
  )
  expect_brma_samples_matrix(true_effects(fit_robma_mv), n, "RoBMA.mv effects")
  expect_brma_samples_matrix(blup(fit_robma_mv), n, "RoBMA.mv BLUP")

  expect_s3_class(RoBMA::as_draws(fit_robma_mv), "draws")
  expect_s3_class(RoBMA::as_draws_array(fit_robma_mv), "draws_array")
  expect_s3_class(RoBMA::as_draws_df(fit_robma_mv), "draws_df")
  expect_s3_class(RoBMA::as_draws_matrix(fit_robma_mv), "draws_matrix")
  expect_s3_class(RoBMA::as_draws_rvars(fit_robma_mv), "draws_rvars")
  expect_s3_class(
    marginal_means(fit_robma_mv, n_samples = 100L),
    "marginal_means.brma"
  )
  expect_s3_class(vif(fit_robma_mv), "vif.brma")
  expect_s3_class(interpret(fit_robma_mv), "interpret.brma")
  expect_s3_class(
    summary_heterogeneity(fit_robma_mv),
    "summary_heterogeneity.brma_list"
  )
})


test_that("RoBMA.mv retains established selection diagnostic semantics", {

  n <- nobs(fit_robma_mv)
  expect_length(residuals(fit_robma_mv), n)
  expect_equal(nrow(rstudent(fit_robma_mv)), n)
  expect_equal(nrow(dfbetas(fit_robma_mv)), n)
  expect_length(covratio(fit_robma_mv), n)

  unavailable <- list(
    rstandard = function() rstandard(fit_robma_mv),
    hatvalues = function() hatvalues(fit_robma_mv),
    dffits = function() dffits(fit_robma_mv),
    cooks.distance = function() cooks.distance(fit_robma_mv),
    influence = function() influence(fit_robma_mv)
  )
  for (method in names(unavailable)) {
    expect_error(
      unavailable[[method]](),
      paste0(method, " is not available for selection models (weightfunction)."),
      fixed = TRUE
    )
  }
})


test_that("RoBMA.mv plots and hypotheses use the combined ensemble", {

  expect_s3_class(as_metafor_forest(fit_robma_mv, addpred = TRUE),
                  "metafor_forest.brma")
  expect_s3_class(
    suppressWarnings(funnel(
      fit_robma_mv,
      type      = "rstudent",
      plot_type = "ggplot"
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(qqnorm(
      fit_robma_mv,
      type      = "rstudent",
      reps      = 50L,
      plot_type = "ggplot"
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(regplot(
      fit_robma_mv,
      mod       = "x",
      plot_type = "ggplot"
    )),
    "ggplot"
  )
  expect_s3_class(as_zplot(fit_robma_mv), "zplot_brma")
  expect_s3_class(
    plot(fit_robma_mv, "study: sd", plot_type = "ggplot"),
    "ggplot"
  )
  expect_s3_class(
    hypothesis(
      fit_robma_mv,
      "study: sd < 0.1",
      density_method = "KDE"
    ),
    "BayesTools_hypothesis_BF"
  )
})


test_that("RoBMA.mv preserves product-space limitations and update", {

  expect_error(add_marglik(fit_robma_mv), "model-averaging objects", fixed = TRUE)
  expect_error(
    bridgesampling::bridge_sampler(fit_robma_mv),
    "model-averaging objects",
    fixed = TRUE
  )
  expect_error(AIC(fit_robma_mv), "not defined", fixed = TRUE)
  expect_error(BIC(fit_robma_mv), "not defined", fixed = TRUE)

  n_before <- nrow(.get_posterior_samples(fit_robma_mv[["fit"]]))
  extended <- suppressWarnings(update(
    fit_robma_mv,
    sample_extend      = 20L,
    recompute          = "drop",
    silent             = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  ))

  expect_s3_class(extended, "RoBMA.mv")
  expect_equal(
    nrow(.get_posterior_samples(extended[["fit"]])),
    n_before + 40L
  )
  expect_null(extended[["loo"]])
  expect_identical(extended[["selection_likelihood"]][["type"]], "exact")
})
