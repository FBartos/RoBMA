source(testthat::test_path("common-functions.R"))

skip_if_missing_fits("BMA.mv_random_components")
fit_bma_mv <- load_fit("BMA.mv_random_components", validate = FALSE)


.bma_mv_random_gate_names <- function(fit) {

  allocation <- fit[["formula_design"]][["mu"]][["random_allocations"]][[1L]]
  vapply(
    allocation[["inclusion"]],
    `[[`,
    character(1),
    "indicator_name"
  )
}


.bma_mv_parameter_draws <- function(fit, parameter) {

  catalog <- BayesTools::parameter_catalog(fit[["fit"]])
  selection <- BayesTools::parameter_catalog_resolve(catalog, parameter)
  as.numeric(as.matrix(BayesTools::parameter_draws(
    fit[["fit"]],
    selection
  ))[, 1L])
}


test_that("BMA.mv summary reports exact random-component inclusion states", {

  samples    <- .get_posterior_samples(fit_bma_mv[["fit"]])
  gate_names <- .bma_mv_random_gate_names(fit_bma_mv)
  out <- summary(
    fit_bma_mv,
    conditional              = TRUE,
    include_mcmc_diagnostics = FALSE
  )
  inclusion <- out[["inclusion_random"]]

  expect_equal(rownames(inclusion), names(gate_names))
  expect_equal(unname(inclusion[["prior_prob"]]), c(0.5, 0.5))
  expect_equal(
    unname(inclusion[["post_prob"]]),
    unname(colMeans(samples[, gate_names, drop = FALSE]))
  )

  summary_frame <- as.data.frame(out)
  expect_equal(
    summary_frame[["parameter"]][summary_frame[["component"]] == "inclusion random"],
    names(gate_names)
  )
  expect_false(any(grepl(
    "__xRE_",
    rownames(out[["inclusion_mods"]]),
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "__xRE_",
    c(
      rownames(out[["estimates_random"]]),
      rownames(out[["estimates_random_conditional"]])
    ),
    fixed = TRUE
  )))
  catalog <- BayesTools::parameter_catalog(fit_bma_mv[["fit"]])
  expect_false(any(grepl(
    "__xRE_",
    catalog[["quantities"]][["canonical_name"]],
    fixed = TRUE
  )))
  expect_match(
    attr(out[["estimates_random"]], "footnotes"),
    "before independent component gates",
    fixed = TRUE
  )
  expect_match(
    attr(out[["estimates_random_conditional"]], "footnotes"),
    "condition on their own inclusion gate",
    fixed = TRUE
  )

  for (component in names(gate_names)) {
    sd_draws <- .bma_mv_parameter_draws(
      fit_bma_mv,
      paste0("(mu) ", component, ": sd(intercept)")
    )
    gate <- samples[, gate_names[[component]]]
    expect_true(all(sd_draws[gate == 0] == 0), info = component)
    expect_equal(
      out[["estimates_random_conditional"]][paste0(component, ": sd"), "Mean"],
      mean(sd_draws[gate == 1]),
      tolerance = 5e-4,
      info = component
    )
  }
})


test_that("BMA.mv component SDs use gated non-renormalized slab allocations", {

  samples    <- .get_posterior_samples(fit_bma_mv[["fit"]])
  gate_names <- .bma_mv_random_gate_names(fit_bma_mv)
  sd_total   <- .bma_mv_parameter_draws(fit_bma_mv, "(mu) sd_total")

  for (component in names(gate_names)) {
    actual <- .bma_mv_parameter_draws(
      fit_bma_mv,
      paste0("(mu) ", component, ": sd(intercept)")
    )
    allocation <- .bma_mv_parameter_draws(
      fit_bma_mv,
      paste0("(mu) var_prop(", component, ")")
    )
    expected <- sd_total * samples[, gate_names[[component]]] * sqrt(allocation)
    expect_equal(actual, as.numeric(expected), tolerance = 1e-12, info = component)
  }

  both_off <- rowSums(samples[, gate_names, drop = FALSE]) == 0
  study_allocation <- .bma_mv_parameter_draws(
    fit_bma_mv,
    "(mu) var_prop(study)"
  )
  expect_true(any(both_off))
  expect_true(all(study_allocation[both_off] > 0))
  expect_true(all(study_allocation[both_off] < 1))
})


test_that("BMA.mv model tables include independent random gates", {

  marginal <- summary_models(
    fit_bma_mv,
    type                     = "marginal",
    include_mcmc_diagnostics = FALSE
  )
  individual <- summary_models(
    fit_bma_mv,
    type                     = "individual",
    include_mcmc_diagnostics = FALSE
  )[["individual"]]

  expect_named(
    marginal[["marginal"]],
    c(
      "Effect", "Location: x", "Heterogeneity Slab",
      "Random: study", "Random: observation"
    )
  )
  expect_equal(
    marginal[["marginal"]][["Heterogeneity Slab"]][["prior_prob"]],
    c(0.6, 0.4)
  )
  expect_equal(nrow(individual), 32L)
  expect_equal(
    sort(unique(individual[["prior_prob"]])),
    c(0.025, 0.0375)
  )
  expect_equal(sum(individual[["post_prob"]]), 1, tolerance = 1e-12)
  expect_equal(
    sort(unique(individual[["Random: study"]])),
    c("Excluded", "Included")
  )
  expect_equal(
    sort(unique(individual[["Random: observation"]])),
    c("Excluded", "Included")
  )
})


test_that("BMA.mv model tables omit fixed random gates", {

  prior_list <- attr(fit_bma_mv[["fit"]], "prior_list", exact = TRUE)
  gate_names <- .bma_mv_random_gate_names(fit_bma_mv)
  gate_prior_index <- which(vapply(prior_list, function(prior) {
    identical(
      attr(prior, "random_allocation_indicator", exact = TRUE),
      gate_names[["study"]]
    )
  }, logical(1)))
  expect_length(gate_prior_index, 1L)

  for (fixed_value in c(0, 1)) {
    fixed_prior <- BayesTools::prior(
      "spike",
      parameters = list(location = fixed_value)
    )
    attr(fixed_prior, "random_allocation_indicator") <- gate_names[["study"]]
    fixed_prior_list <- prior_list
    fixed_prior_list[[gate_prior_index]] <- fixed_prior

    components <- .summary_models_add_random_components(
      components = list(),
      object      = fit_bma_mv,
      prior_list  = fixed_prior_list
    )
    expect_named(components, "Random: observation")
  }
})


test_that("BMA.mv prediction and diagnostic targets remain available", {

  n <- nobs(fit_bma_mv)
  terms <- predict(
    fit_bma_mv,
    type               = "terms",
    conditioning_depth = "marginal"
  )
  estimate <- predict(
    fit_bma_mv,
    type               = "estimate",
    conditioning_depth = "marginal"
  )
  response <- predict(
    fit_bma_mv,
    type               = "response",
    conditioning_depth = "marginal"
  )
  scale_conditional <- predict(
    fit_bma_mv,
    type        = "terms.scale",
    conditional = TRUE,
    quiet       = TRUE
  )
  fitted_estimate <- fitted(
    fit_bma_mv,
    type               = "estimate",
    conditioning_depth = "estimate",
    summary             = FALSE
  )

  expect_brma_samples_matrix(terms, n, "BMA.mv terms")
  expect_brma_samples_matrix(estimate, n, "BMA.mv estimate")
  expect_brma_samples_matrix(response, n, "BMA.mv response")
  expect_named(scale_conditional, names(.bma_mv_random_gate_names(fit_bma_mv)))
  posterior_samples <- .get_posterior_samples(fit_bma_mv[["fit"]])
  for (component in names(scale_conditional)) {
    gate <- .bma_mv_random_gate_names(fit_bma_mv)[[component]]
    expect_equal(
      nrow(scale_conditional[[component]]),
      sum(posterior_samples[, gate] == 1),
      info = component
    )
  }
  expect_length(fitted_estimate, n)
  expect_true(all(is.finite(fitted_estimate)))
  expect_brma_samples_matrix(pooled_effect(fit_bma_mv), 1L, "BMA.mv pooled effect")
  expect_brma_samples_matrix(
    pooled_heterogeneity(fit_bma_mv, component = "total"),
    1L,
    "BMA.mv pooled heterogeneity"
  )
  pooled_total_conditional <- pooled_heterogeneity(
    fit_bma_mv,
    conditional = TRUE,
    component   = "total"
  )
  expect_equal(
    nrow(pooled_total_conditional),
    sum(rowSums(posterior_samples[, .bma_mv_random_gate_names(fit_bma_mv),
                                  drop = FALSE]) > 0)
  )

  expect_equal(dim(log_lik(fit_bma_mv)), c(nrow(terms), n))
  expect_length(residuals(fit_bma_mv), n)
  expect_equal(nrow(rstandard(fit_bma_mv)), n)
  expect_equal(nrow(rstudent(fit_bma_mv)), n)
  expect_length(hatvalues(fit_bma_mv), n)
  expect_equal(nrow(dfbetas(fit_bma_mv)), n)
  expect_length(dffits(fit_bma_mv), n)
  expect_length(cooks.distance(fit_bma_mv), n)
  expect_length(covratio(fit_bma_mv), n)
  expect_s3_class(suppressWarnings(influence(fit_bma_mv)), "infl.brma")
})


test_that("BMA.mv formula and ensemble post-processing remains available", {

  random_effects <- ranef(fit_bma_mv, simplify = FALSE)
  expect_named(random_effects, c("study", "observation"))
  expect_brma_samples_matrix(
    ranef(fit_bma_mv, component = "total", expand = TRUE),
    nobs(fit_bma_mv),
    "BMA.mv total random effect"
  )
  expect_brma_samples_matrix(
    true_effects(fit_bma_mv),
    nobs(fit_bma_mv),
    "BMA.mv true effects"
  )
  expect_brma_samples_matrix(
    blup(fit_bma_mv),
    nobs(fit_bma_mv),
    "BMA.mv BLUP"
  )

  draws <- RoBMA::as_draws(fit_bma_mv)
  expect_s3_class(draws, "draws")
  expect_s3_class(RoBMA::as_draws_array(fit_bma_mv), "draws_array")
  expect_s3_class(RoBMA::as_draws_df(fit_bma_mv), "draws_df")
  expect_s3_class(RoBMA::as_draws_matrix(fit_bma_mv), "draws_matrix")
  expect_s3_class(RoBMA::as_draws_rvars(fit_bma_mv), "draws_rvars")
  expect_true(all(
    .bma_mv_random_gate_names(fit_bma_mv) %in% posterior::variables(draws)
  ))
  slab_indicators <- paste0(
    .random_slab_prior_parameters(
      attr(fit_bma_mv[["fit"]], "prior_list", exact = TRUE)
    ),
    "_indicator"
  )
  expect_true(all(slab_indicators %in% posterior::variables(draws)))
  expect_s3_class(marginal_means(fit_bma_mv, n_samples = 100L), "marginal_means.brma")
  expect_s3_class(vif(fit_bma_mv), "vif.brma")
  expect_s3_class(interpret(fit_bma_mv), "interpret.brma")
  expect_s3_class(summary_heterogeneity(fit_bma_mv), "summary_heterogeneity.brma_list")
})


test_that("BMA.mv LOO and WAIC use the model-averaged predictive density", {

  log_likelihood <- log_lik(fit_bma_mv)
  loo_result      <- loo(fit_bma_mv)
  weights         <- loo_weights(fit_bma_mv)

  expect_s3_class(loo_result, "loo")
  expect_equal(dim(weights), dim(log_likelihood))
  expect_equal(colSums(weights), rep(1, nobs(fit_bma_mv)), tolerance = 1e-10)
  expect_no_error(suppressWarnings(check_loo(fit_bma_mv)))

  fit_waic <- fit_bma_mv
  fit_waic[["waic"]] <- NULL
  fit_waic <- suppressWarnings(add_waic(fit_waic))
  expect_s3_class(waic(fit_waic), "waic")
})


test_that("BMA.mv forest and diagnostic plot data use multivariate targets", {

  forest_data <- as_metafor_forest(fit_bma_mv, addpred = TRUE)
  forest_conditional <- as_metafor_forest(
    fit_bma_mv,
    addpred     = TRUE,
    conditional = TRUE
  )
  expect_s3_class(forest_data, "metafor_forest.brma")
  expect_s3_class(forest_conditional, "metafor_forest.brma")
  expect_s3_class(
    suppressWarnings(funnel(
      fit_bma_mv,
      type      = "rstandard",
      plot_type = "ggplot"
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(qqnorm(
      fit_bma_mv,
      type      = "rstandard",
      reps      = 50L,
      plot_type = "ggplot"
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(regplot(fit_bma_mv, mod = "x", plot_type = "ggplot")),
    "ggplot"
  )
  expect_s3_class(as_zplot(fit_bma_mv), "zplot_brma")
})


test_that("BMA.mv random plots and hypotheses respect component gates", {

  expect_s3_class(
    plot(fit_bma_mv, "study: sd", plot_type = "ggplot"),
    "ggplot"
  )
  expect_s3_class(
    plot(
      fit_bma_mv,
      "study: sd",
      conditional = TRUE,
      plot_type    = "ggplot"
    ),
    "ggplot"
  )
  expect_s3_class(
    hypothesis(
      fit_bma_mv,
      "study: sd < 0.1",
      density_method = "KDE"
    ),
    "BayesTools_hypothesis_BF"
  )
  expect_s3_class(
    hypothesis(
      fit_bma_mv,
      "study: sd < 0.1",
      conditional    = TRUE,
      density_method = "KDE"
    ),
    "BayesTools_hypothesis_BF"
  )
  expect_error(
    hypothesis(
      fit_bma_mv,
      "study: sd != 0 vs study: sd = 0",
      density_method = "KDE"
    ),
    "Random-Effect Inclusion table",
    fixed = TRUE
  )
})


test_that("BMA.mv preserves established product-space limitations", {

  expect_error(add_marglik(fit_bma_mv), "model-averaging objects", fixed = TRUE)
  expect_error(
    bridgesampling::bridge_sampler(fit_bma_mv),
    "model-averaging objects",
    fixed = TRUE
  )
  expect_error(radial(fit_bma_mv), "models that contain moderators", fixed = TRUE)
  expect_error(AIC(fit_bma_mv), "not defined", fixed = TRUE)
  expect_error(BIC(fit_bma_mv), "not defined", fixed = TRUE)
})


test_that("BMA.mv fits can be extended without rebuilding the product space", {

  n_before <- nrow(.get_posterior_samples(fit_bma_mv[["fit"]]))
  extended <- suppressWarnings(update(
    fit_bma_mv,
    sample_extend     = 20L,
    recompute         = "drop",
    silent            = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  ))

  expect_s3_class(extended, "BMA.mv")
  expect_equal(
    nrow(.get_posterior_samples(extended[["fit"]])),
    n_before + 40L
  )
  expect_null(extended[["loo"]])
  expect_length(
    summary(extended, include_mcmc_diagnostics = FALSE)[["inclusion_random"]][["post_prob"]],
    2L
  )
})
