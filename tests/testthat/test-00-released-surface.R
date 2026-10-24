# Argument order, defaults, and entry points of RoBMA 4.0.0 that later
# development must keep. New arguments follow the released ones.

.released_arguments_norm <- c(
  "yi", "vi", "sei", "weights", "ni", "mods", "scale", "cluster", "data",
  "slab", "subset", "measure"
)

test_that("selection arguments follow the released RoBMA() and bselmodel() arguments", {

  released_RoBMA <- c(
    .released_arguments_norm, "effect_direction",
    "prior_effect", "prior_heterogeneity", "prior_mods", "prior_scale",
    "prior_heterogeneity_allocation", "prior_bias",
    "prior_effect_null", "prior_heterogeneity_null", "prior_mods_null",
    "prior_scale_null", "prior_heterogeneity_allocation_null", "prior_bias_null",
    "standardize_continuous_predictors", "set_contrast_factor_predictors",
    "prior_unit_information_sd", "rescale_priors",
    "prior_informed_field", "prior_informed_subfield", "model_type",
    "sample", "burnin", "adapt", "chains", "thin", "parallel", "autofit",
    "autofit_control", "convergence_checks", "seed", "silent"
  )
  released_bselmodel <- c(
    .released_arguments_norm,
    "prior_effect", "prior_heterogeneity", "prior_mods", "prior_scale",
    "prior_heterogeneity_allocation", "prior_bias",
    "standardize_continuous_predictors", "set_contrast_factor_predictors",
    "prior_unit_information_sd", "rescale_priors",
    "prior_informed_field", "prior_informed_subfield",
    "effect_direction", "steps",
    "sample", "burnin", "adapt", "chains", "thin", "parallel", "autofit",
    "autofit_control", "convergence_checks", "seed", "silent"
  )

  for (case in list(
    list(name = "RoBMA",     fun = RoBMA,     released = released_RoBMA),
    list(name = "bselmodel", fun = bselmodel, released = released_bselmodel)
  )) {
    arguments <- names(formals(case[["fun"]]))
    expect_identical(arguments, c(case[["released"]], "selection", "selection_control", "..."),
                     info = case[["name"]])
  }
})

test_that("parameter selectors of plot(), plot_prior(), and print_prior() have no default", {

  no_default <- function(fun, name) identical(formals(fun)[[name]], quote(expr = ))

  for (name in c("parameter", "parameter_mods", "parameter_scale")) {
    expect_true(no_default(plot.brma, name), info = paste("plot.brma", name))
    expect_true(no_default(print_prior.brma, name), info = paste("print_prior.brma", name))
  }
  expect_identical(formals(plot_prior.brma)[["parameter"]], "mu")
  expect_true(no_default(plot_prior.brma, "parameter_mods"))
  expect_true(no_default(plot_prior.brma, "parameter_scale"))

  prior_object <- RoBMA(
    yi = c(.1, .2, .3), sei = c(.1, .1, .1), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  expect_no_error(utils::capture.output(print_prior(prior_object)))
  expect_no_error(utils::capture.output(print_prior(prior_object, parameter = "tau")))
  expect_no_error(plot_prior(prior_object))
})


# Small fits shared by the tests below (two chains of 100 draws, built once).
.surface_fit_env <- new.env()

.surface_data <- function() {
  data.frame(
    yi      = c(.12, .30, .05, .41, .22, .18, .35, .10),
    sei     = c(.10, .12, .09, .15, .11, .10, .13, .12),
    x       = c(1, 2, 3, 4, 5, 6, 7, 8) / 8,
    cluster = rep(1:4, each = 2)
  )
}

.surface_fit <- function(name) {

  if (!exists(name, envir = .surface_fit_env, inherits = FALSE)) {
    arguments <- list(
      yi = quote(yi), sei = quote(sei), data = .surface_data(), measure = "GEN",
      prior_unit_information_sd = 1,
      chains = 2, adapt = 50, burnin = 100, sample = 100,
      seed = 1, silent = TRUE
    )
    arguments <- switch(
      name,
      ordinary   = arguments,
      mods       = c(arguments, list(mods = ~ x)),
      multilevel = c(arguments, list(cluster = quote(cluster)))
    )
    fit <- suppressWarnings(do.call(brma, arguments))
    assign(name, fit, envir = .surface_fit_env)
  }

  return(get(name, envir = .surface_fit_env, inherits = FALSE))
}

test_that("predict(newdata = TRUE) returns the aggregated prediction", {

  mods <- .surface_fit("mods")

  terms <- predict(mods, newdata = TRUE, type = "terms", quiet = TRUE)
  fixed <- predict(mods, type = "terms")
  expect_s3_class(terms, "brma_samples")
  expect_identical(colnames(terms), "mu")
  expect_identical(attr(terms, "title"), "Aggregated Location Term Posterior Prediction:")
  expect_equal(as.numeric(terms), unname(rowMeans(unclass(fixed))), tolerance = 1e-12)
  expect_message(
    predict(mods, newdata = TRUE, type = "terms"),
    "Aggregated prediction averages mu across the moderator model matrix (K = 8 observations)",
    fixed = TRUE
  )
  expect_equal(predict(mods, newdata = TRUE, type = "marginal", quiet = TRUE), terms)

  scale <- predict(mods, newdata = TRUE, type = "terms.scale")
  expect_identical(colnames(scale), "tau")
  expect_identical(attr(scale, "title"), "Aggregated Scale Term Posterior Prediction:")
  expect_equal(unclass(scale)[, "tau"], unclass(pooled_heterogeneity(mods))[, "tau"],
               ignore_attr = TRUE)

  set.seed(1)
  estimate <- predict(mods, newdata = TRUE, type = "estimate", quiet = TRUE)
  expect_identical(colnames(estimate), "theta")
  expect_identical(attr(estimate, "title"), "Aggregated True Effect Posterior Prediction:")
  expect_equal(unname(colMeans(unclass(estimate))), unname(colMeans(unclass(terms))),
               tolerance = 1)
  expect_gt(sd(estimate[, "theta"]), sd(terms[, "mu"]))
  set.seed(1)
  expect_equal(predict(mods, newdata = TRUE, type = "blup", quiet = TRUE), estimate)

  expect_error(
    predict(mods, newdata = TRUE, type = "response"),
    "Aggregated predictions (newdata = TRUE) are not available for type = 'response'",
    fixed = TRUE
  )
  expect_error(
    predict(mods, newdata = TRUE, type = "cluster"),
    "type = 'cluster' is only available for multilevel (3-level) models",
    fixed = TRUE
  )
  expect_error(
    predict(mods, newdata = TRUE, conditioning_depth = "estimate"),
    "available only with newdata = NULL"
  )
  expect_error(predict(mods, newdata = FALSE), "'newdata' must be NULL, TRUE, a data.frame")
})

test_that("predict(newdata = TRUE, type = 'cluster') adds a new cluster effect", {

  multilevel <- .surface_fit("multilevel")

  set.seed(1)
  cluster <- predict(multilevel, newdata = TRUE, type = "cluster")
  expect_identical(colnames(cluster), "mu_cluster")
  expect_identical(attr(cluster, "title"), "Aggregated Cluster-Level Posterior Prediction:")
  terms <- predict(multilevel, newdata = TRUE, type = "terms")
  expect_gt(sd(cluster[, "mu_cluster"]), sd(terms[, "mu"]))
})

test_that("marginal_means() accepts normal_approximation at its released position", {

  expect_identical(
    head(names(formals(marginal_means.brma)), 6L),
    c("object", "null_hypothesis", "normal_approximation", "n_samples",
      "output_measure", "transform")
  )

  mods <- .surface_fit("mods")
  kde    <- marginal_means(mods, bf = TRUE, n_samples = 1000)
  normal <- marginal_means(mods, normal_approximation = TRUE, bf = TRUE, n_samples = 1000)

  expect_false(kde[["normal_approximation"]])
  expect_true(normal[["normal_approximation"]])
  expect_no_warning(summary(normal))
  expect_s3_class(summary(normal), "summary.marginal_means.brma")
  bf <- function(x) vapply(x[["inference"]][["inference"]], function(row) as.numeric(row[[1]]), numeric(1))
  expect_true(all(is.finite(bf(normal))))
  expect_false(isTRUE(all.equal(bf(kde), bf(normal))))
  expect_error(
    marginal_means(mods, normal_approximation = NA),
    "normal_approximation"
  )
})

test_that("update(seed = ) extends reproducibly", {

  fit <- .surface_fit("ordinary")

  extend <- function(seed) {
    as.numeric(as_draws_matrix(update(fit, sample_extend = 20, seed = seed))[, "mu"])
  }
  first  <- extend(1)
  second <- extend(1)
  other  <- extend(2)

  expect_length(first, 2 * (100 + 20))
  expect_identical(first, second)
  expect_false(identical(first, other))
  expect_identical(update(fit, sample_extend = 20, seed = 3)[["fit_control"]][["seed"]], 3)
})
