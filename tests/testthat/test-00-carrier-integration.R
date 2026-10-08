.carrier_integration_fit_cache <- new.env(parent = emptyenv())

.carrier_integration_fits <- function(){

  if(!is.null(.carrier_integration_fit_cache$fits)) return(.carrier_integration_fit_cache$fits)
  set.seed(1)
  k <- 48L
  data <- data.frame(g1 = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
    g2 = factor(rep(1:4, length.out = k)[sample.int(k)], levels = 1:4), sei = stats::runif(k, .1, .3))
  data$yi <- stats::rnorm(k, c(0, .2, .4)[as.integer(data$g1)] + c(0, .1, .3, .5)[as.integer(data$g2)], data$sei)
  fit <- function(mods, contrast = "treatment", prior_effect = NULL, prior_mods = NULL, averaged = FALSE){
    args <- list(yi = quote(yi), sei = quote(sei), mods = mods, data = data, measure = "SMD",
      set_contrast_factor_predictors = contrast, chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1, silent = TRUE)
    if(!is.null(prior_effect)) args$prior_effect <- prior_effect
    if(!is.null(prior_mods)) args$prior_mods <- prior_mods
    expect_warning(do.call(if(averaged) BMA.norm else brma, args), "Only one chain was run", fixed = TRUE)
  }
  fits <- list(treatment = fit(~g1-1), meandif = fit(~g1-1, "meandif"),
    nonzero = fit(~g1, prior_effect = BayesTools::prior("point", list(1))), averaged = fit(~g1, averaged = TRUE),
    ordinary = fit(~g1, "meandif"), gated = fit(~g1, prior_effect = BayesTools::prior("point", list(1)),
      prior_mods = list(g1 = BayesTools::prior_ordered(BayesTools::prior_spike_and_slab(BayesTools::prior("gamma", list(2, 2))), allocation = c(.25, .75)))))
  .carrier_integration_fit_cache$fits <- fits
  fits
}

test_that("whole marginal-means regions retain structural fixed-cell refusal", {
  skip_on_cran()
  means <- marginal_means(.carrier_integration_fits()$nonzero, bf = FALSE)
  expect_error(hypothesis(means, "g1 > 0", density_method = "KDE"), class = "RoBMA_hypothesis_fixed")
  expect_error(hypothesis(means, "g1[5] > 0", density_method = "KDE"), class = "RoBMA_hypothesis_fixed")
  expect_error(hypothesis(means, "g1 = 0", density_method = "KDE"), class = "RoBMA_hypothesis_fixed")
  expect_s3_class(hypothesis(means, "g1[10] > 0", density_method = "KDE"), "BayesTools_hypothesis_BF")
})

test_that("expanded public factor groups retain every statement-major result row", {
  skip_on_cran()
  fit <- .carrier_integration_fits()$ordinary
  statements <- c("g1 > 0", "mu > 0", "g1 > 0")
  separate <- lapply(statements, function(statement) hypothesis(fit, statement, density_method = "KDE", columns = "all", seed = 1, standardized_coefficients = TRUE))
  combined <- hypothesis(fit, statements, density_method = "KDE", columns = "all", seed = 1, standardized_coefficients = TRUE)
  expect_identical(nrow(combined), 7L)
  expected <- unlist(lapply(separate, function(result) attr(result, "raw_log_BF")), use.names = FALSE)
  expect_identical(attr(combined, "raw_log_BF"), expected)
  expect_length(attr(combined, "numerical_diagnostics"), 7L)
  expect_length(attr(combined, "prior_numerical_diagnostics"), 7L)
  expect_length(attr(combined, "hypothesis_ast")$statements, 3L)
  for(logBF in c(FALSE, TRUE)) for(BF01 in c(FALSE, TRUE)){
    viewed <- update(combined, logBF = logBF, BF01 = BF01)
    expect_identical(attr(viewed, "raw_log_BF"), expected)
    expected_view <- unlist(lapply(separate, function(result) as.numeric(update(result, logBF = logBF, BF01 = BF01)$BF)), use.names = FALSE)
    expect_identical(as.numeric(viewed$BF), expected_view)
  }
  plain <- .output_table_as_long_data_frame(combined, "location")
  expect_identical(nrow(plain), 7L)
  expect_null(attr(plain$BF, "canonical_log_BF"))
})

test_that("point and region marginal groups select their own execution method", {
  skip_on_cran()
  means <- marginal_means(.carrier_integration_fits()$ordinary, bf = FALSE)
  statements <- c("g1[10] - g1[20] = 0", "g1[10] > 0")
  for(method in c("qCMDE", "IWMDE")){
    separate <- lapply(statements, function(statement) hypothesis(means, statement, density_method = method, columns = "all", seed = 1))
    combined <- hypothesis(means, statements, density_method = method, columns = "all", seed = 1)
    reversed <- hypothesis(means, rev(statements), density_method = method, columns = "all", seed = 1)
    expected <- unlist(lapply(separate, function(result) attr(result, "raw_log_BF")), use.names = FALSE)
    expect_equal(attr(combined, "raw_log_BF"), expected, tolerance = 1e-12)
    expect_equal(attr(reversed, "raw_log_BF"), rev(expected), tolerance = 1e-12)
    expect_length(attr(combined, "numerical_diagnostics"), 2L)
  }
})

test_that("derived public OR targets carry their effective declaration to scalar attachment", {
  skip_on_cran()
  means <- marginal_means(.carrier_integration_fits()$averaged, bf = FALSE)
  plan <- .hypothesis_plan_marginal_means(means, "2 * g1[10] = 0")
  declared <- .iwmde_sample_condition_metadata(plan$linear_target$posterior)
  expect_identical(declared$conditional_rule, "OR")
  expect_length(declared$conditional, 2L)
  captured <- NULL
  local_mocked_bindings(.iwmde_estimate = function(context, parameter_spec, ...){
    captured <<- list(context = context, spec = parameter_spec)
    stop("Captured actual public scalar attach", call. = FALSE)
  }, .package = "RoBMA")
  expect_error(hypothesis(means, "2 * g1[10] = 0", density_method = "qCMDE", seed = 1), "Captured actual public scalar attach", fixed = TRUE)
  for(field in names(declared)) expect_identical(captured$spec[[field]], declared[[field]])
  expected_rows <- .iwmde_parameter_condition_rows(captured$context, declared)
  expect_identical(.iwmde_parameter_condition_rows(captured$context, captured$spec), expected_rows)
  and <- declared; and$conditional_rule <- "AND"
  expect_gt(sum(expected_rows), sum(.iwmde_parameter_condition_rows(captured$context, and)))
})

test_that("ordered global totals participate in marginal inclusion conditioning", {
  skip_on_cran()
  fit <- .carrier_integration_fits()$gated
  prior <- attr(fit$fit, "prior_list")$mu_g1
  expect_true(.marginal_means_prior_is_conditional(prior))
  expect_false(.marginal_means_prior_is_conditional(BayesTools::prior_ordered(BayesTools::prior("gamma", list(2, 2)))))
  means <- marginal_means(fit, bf = FALSE)
  conditional <- means$inference$conditional$mu_g1[[3L]]
  condition <- BayesTools::posterior_metadata(conditional, "condition")
  expect_identical(condition$effective_conditional, "mu_g1")
  expect_false(condition$averaged)
  expect_lt(length(conditional), length(means$inference$averaged$mu_g1[[3L]]))
})

test_that("declared stochastic contribution refusals are planned for levels and combinations", {
  skip_on_cran()
  data <- data.frame(g = factor(rep(c("A", "B"), 3L)), yi = c(-.2, .3, 0, .4, .1, .2), sei = rep(.5, 6L))
  for(multiplier in list(2, "tau")){
    slope <- BayesTools::prior_factor("normal", list(0, 1), contrast = "treatment")
    attr(slope, "multiply_by") <- multiplier
    fit <- expect_warning(brma(yi = yi, sei = sei, mods = ~g, data = data, measure = "GEN",
      prior_effect = BayesTools::prior("normal", list(0, 1)), prior_mods = list(g = slope),
      prior_heterogeneity = BayesTools::prior("normal", list(0, 1), truncation = list(lower = 0)),
      chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 281, silent = TRUE), "Only one chain was run", fixed = TRUE)
    means <- marginal_means(fit, bf = FALSE)
    for(statement in c("g[B] = 0.1", "g[B] - g[A] = 0.1")){
      plan <- .hypothesis_plan_marginal_means(means, statement)
      expect_null(plan$methods$KDE)
      if(is.numeric(multiplier)){
        expect_null(plan$methods$qCMDE)
        expect_null(plan$methods$IWMDE)
      }else{
        for(method in c("qCMDE", "IWMDE")){
          expect_match(plan$methods[[method]]$reason, "stochastic formula contribution products", fixed = TRUE)
          expect_true("RoBMA_density_method_original_scale" %in% plan$methods[[method]]$class)
          expect_error(hypothesis(means, statement, density_method = method), class = "RoBMA_hypothesis_method")
        }
      }
    }
  }
})

test_that("random groups select point-owned attachment metadata independently of order", {
  # Existing compiled stand-in: routing only. The saved actual known-V fit
  # separately verifies this public path; the estimator is not executed here.
  object <- gated_random_object()
  attr(object$data, "known_V") <- TRUE
  attr(object$data, "measure") <- "GEN"
  statements <- paste0("`(mu) study: tau(intercept)`", c(" > 0.1", " = 0.2"))
  target <- .hypothesis_plans(object, statements)[[2L]]$density_target
  captured <- list()
  local_mocked_bindings(.hypothesis_plan_attach_scalar = function(parameter, spec, display_transform = NULL, ...){
    captured[[length(captured) + 1L]] <<- list(parameter = parameter, spec = spec, display_transform = display_transform)
    stop("Captured random point attachment", call. = FALSE)
  }, .package = "RoBMA")
  for(hypothesis in list(statements, rev(statements), statements[2L])){
    expect_error(hypothesis(object, hypothesis, density_method = "qCMDE"), "Captured random point attachment", fixed = TRUE)
  }
  for(record in captured){
    expect_identical(record$parameter, target$parameter)
    expect_identical(record$spec, target$parameter_spec)
    expect_identical(record$display_transform, target$display_transform)
  }
})
