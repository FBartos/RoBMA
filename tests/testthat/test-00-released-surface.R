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
    z       = c(.2, .1, .5, .3, .9, .4, .8, .6),
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
      mods2      = c(arguments, list(mods = ~ x + z)),
      scale      = c(arguments, list(mods = ~ x, scale = ~ x)),
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

test_that("update() has no seed argument and the extension continues the stored chain states", {

  expect_false("seed" %in% names(formals(update.brma)))
  expect_identical(
    names(formals(update.brma)),
    c("object", "formula.", "...", "sample_extend", "slab", "autofit_control",
      "convergence_checks", "recompute", "parallel", "cores", "silent", "evaluate")
  )

  fit <- .surface_fit("ordinary")
  expect_error(
    update(fit, sample_extend = 20, seed = 1),
    "Unused argument in update.brma(): 'seed'", fixed = TRUE
  )

  # the chains continue their stored generator states, so a repeated extension is identical
  extend <- function() as.numeric(as_draws_matrix(update(fit, sample_extend = 20))[, "mu"])
  first  <- extend()
  expect_length(first, 2 * (100 + 20))
  expect_identical(first, extend())
  expect_identical(
    update(fit, sample_extend = 20)[["fit_control"]][["seed"]],
    fit[["fit_control"]][["seed"]]
  )
})

test_that("unused arguments in the dots of released methods warn and never stop", {

  ordinary   <- .surface_fit("ordinary")
  multilevel <- .surface_fit("multilevel")

  calls <- list(
    predict              = function(x) predict(x, bogus_zz = 1),
    fitted               = function(x) fitted(x, bogus_zz = 1),
    blup                 = function(x) blup(x, bogus_zz = 1),
    true_effects         = function(x) true_effects(x, bogus_zz = 1),
    ranef                = function(x) ranef(x, bogus_zz = 1),
    pooled_effect        = function(x) pooled_effect(x, bogus_zz = 1),
    pooled_heterogeneity = function(x) pooled_heterogeneity(x, bogus_zz = 1),
    summary              = function(x) summary(x, bogus_zz = 1),
    marginal_means       = function(x) marginal_means(x, bogus_zz = 1),
    plot                 = function(x) plot(x, bogus_zz = 1)
  )
  fits <- list(
    ordinary = ordinary, multilevel = multilevel,
    mods = .surface_fit("mods")
  )
  for (name in names(calls)) {
    fit <- if (name == "marginal_means") fits[["mods"]] else if (name == "ranef") multilevel else ordinary
    if (name == "plot") {
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
    }
    expect_warning(
      suppressMessages(calls[[name]](fit)),
      "Unused argument in .*'bogus_zz'",
      info = name
    )
  }
})

test_that("pooled_effect() and pooled_heterogeneity() forward the dots to predict()", {

  mods <- .surface_fit("mods")

  # as_measure is an argument of predict() without a role for the pooled effect
  expect_no_warning(pooled_effect(mods, as_measure = FALSE))
  expect_error(
    pooled_effect(mods, conditioning_depth = "estimate"),
    "available only with newdata = NULL"
  )
  expect_error(
    pooled_heterogeneity(mods, conditioning_depth = "cluster"),
    "available only with newdata = NULL"
  )
  expect_warning(
    pooled_effect(mods, newdata = TRUE),
    "Unused argument in pooled_effect.brma(): 'newdata'",
    fixed = TRUE
  )

  pooled    <- pooled_effect(mods)
  predicted <- predict(mods, newdata = TRUE, type = "terms", bias_adjusted = TRUE, quiet = TRUE)
  expect_identical(attr(pooled, "title"), "Pooled Effect Size")
  expect_equal(unclass(pooled)[, "mu"], unclass(predicted)[, "mu"], ignore_attr = TRUE)
  expect_false(is.null(attr(pooled, "prediction_samples")))
  expect_null(attr(predicted, "prediction_samples"))
  expect_identical(attr(pooled_heterogeneity(mods), "title"), "Pooled Heterogeneity")
})

test_that("logLik() returns the pointwise log-likelihood matrix and AIC() and BIC() refuse", {

  fit <- .surface_fit("multilevel")

  log_likelihood <- stats::logLik(fit)
  expect_s3_class(log_likelihood, "logLik.brma")
  expect_identical(dim(log_likelihood), c(200L, 8L))
  expect_equal(unclass(log_likelihood), log_lik(fit), ignore_attr = TRUE)
  expect_output(print(log_likelihood), "200*8 pointwise log-likelihood matrix", fixed = TRUE)
  expect_identical(ncol(stats::logLik(fit, unit = "cluster")), 4L)

  expect_error(stats::AIC(fit), "AIC() is not defined for brma", fixed = TRUE)
  expect_error(stats::BIC(fit), "BIC() is not defined for brma", fixed = TRUE)
})

test_that("brma_samples coerce to the draws and as_draws() returns a draws_matrix", {

  fit  <- .surface_fit("ordinary")
  pred <- predict(fit, type = "estimate")

  draws <- as.data.frame(pred)
  expect_s3_class(draws, "data.frame")
  expect_identical(dim(draws), c(200L, 8L))
  expect_identical(names(draws), colnames(pred))
  expect_equal(unname(as.matrix(draws)), unname(as.matrix(pred)))
  expect_equal(unname(as.matrix(data.frame(pred))), unname(as.matrix(draws)))

  table <- as.data.frame(summary(pred))
  expect_identical(names(table)[1:2], c("component", "parameter"))
  expect_identical(nrow(table), 8L)

  generic <- as_draws(pred)
  expect_s3_class(generic, "draws_matrix")
  expect_identical(posterior::nchains(generic), 2L)
  expect_s3_class(as_draws_array(pred), "draws_array")
})

test_that("print_prior() lists the location priors before tau and the scale priors", {

  expect_identical(
    names(print_prior(.surface_fit("mods"), silent = TRUE)),
    c("mu_intercept", "mu_x", "tau")
  )
  expect_identical(
    names(print_prior(.surface_fit("scale"), silent = TRUE)),
    c("mu_intercept", "mu_x", "log_tau_intercept", "log_tau_x")
  )
  expect_identical(
    names(print_prior(.surface_fit("ordinary"), silent = TRUE)),
    c("mu", "tau")
  )
})

test_that("the location-scale summary carries the released footnote", {

  fit <- .surface_fit("scale")
  footnote <- paste0(
    "exp(Intercept) corresponds to the between-study heterogeneity tau; ",
    "the meta-regression coefficients correspond to the multiplicative ",
    "effects on log-scale."
  )

  expect_identical(.summary_scale_footnotes(fit), footnote)
  expect_identical(
    attr(summary(fit, include_mcmc_diagnostics = FALSE)[["estimates_scale"]], "footnotes"),
    footnote
  )
})

test_that("ranef() expands by default and blup() labels the true effects theta[i]", {

  multilevel <- .surface_fit("multilevel")

  random <- ranef(multilevel)
  expect_identical(names(random), c("cluster", "estimate"))
  expect_identical(ncol(random[["cluster"]]), 8L)
  expect_true(all(grepl("|", colnames(random[["cluster"]]), fixed = TRUE)))
  expect_true(all(startsWith(colnames(random[["cluster"]]), "u_cluster[")))
  expect_identical(ncol(ranef(multilevel, expand = FALSE)[["cluster"]]), 4L)
  expect_identical(
    ranef(multilevel, expand = TRUE)[["cluster"]],
    random[["cluster"]]
  )

  for (effects in list(blup(multilevel), true_effects(multilevel))) {
    expect_identical(colnames(effects), paste0("theta[", 1:8, "]"))
    expect_identical(attr(effects, "title"), "True Effects (BLUP Means)")
  }
  # the new depth-aware prediction keeps its own title
  estimate <- predict(multilevel, type = "estimate", conditioning_depth = "estimate")
  expect_false(identical(attr(estimate, "title"), "True Effects (BLUP Means)"))
})

test_that("vif() plugs in the posterior mean heterogeneity of a simple model", {

  fit <- .surface_fit("mods2")

  X   <- .get_model_matrix(fit)
  vi  <- .surface_data()[["sei"]]^2
  tau <- mean(as_draws_matrix(fit)[, "tau"])
  covariance  <- solve(crossprod(X, X / (vi + tau^2)))
  correlation <- stats::cov2cor(covariance)[-1L, -1L]
  expected    <- diag(solve(correlation))

  actual <- vif(fit, posterior_correlation = FALSE)[["vif"]]
  expect_equal(actual[["GVIF"]], unname(expected), tolerance = 1e-8)
})

test_that("cooks.distance() divides by the rank of the model matrix", {

  fit <- suppressWarnings(add_loo(.surface_fit("mods2")))
  rank <- qr(.get_model_matrix(fit))[["rank"]]
  expect_identical(rank, 3L)

  distance <- suppressWarnings(cooks.distance(fit))
  expect_identical(.cooks_distance_rank(fit), rank)
  unscaled <- suppressWarnings(
    .cooks.distance_internal(
      .influence_fit_samples(fit),
      .diagnostic_psis_context(fit)[["psis_weights"]]
    )
  )
  expect_equal(unname(distance), unname(unscaled) / rank, tolerance = 1e-10)
  expect_equal(
    unname(suppressWarnings(influence(fit))[["inf"]][["cook.d"]]),
    unname(distance),
    tolerance = 1e-10
  )
})

.surface_glmm_fit <- function() {

  if (!exists("glmm", envir = .surface_fit_env, inherits = FALSE)) {
    data("dat.bcg", package = "metadat", envir = environment())
    fit <- suppressWarnings(brma.glmm(
      ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg[1:8, ],
      measure = "OR", chains = 2, adapt = 100, burnin = 150, sample = 150,
      seed = 1, silent = TRUE
    ))
    assign("glmm", suppressWarnings(add_loo(fit)), envir = .surface_fit_env)
  }

  return(get("glmm", envir = .surface_fit_env, inherits = FALSE))
}

test_that("GLMM LOO-PIT residuals are the labelled normal approximation", {

  skip_if_not_installed("metadat")
  fit <- .surface_glmm_fit()

  count_messages <- function(expr) {
    messages <- character()
    value <- withCallingHandlers(expr, message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
    list(value = value, messages = messages)
  }

  table <- count_messages(suppressWarnings(rstudent(fit)))
  expect_length(table$messages, 1L)
  expect_match(table$messages, "descriptive normal approximation", fixed = TRUE)
  expect_identical(dim(table$value), c(8L, 3L))
  expect_true(all(is.finite(unlist(table$value))))

  vector <- count_messages(suppressWarnings(residuals(fit, type = "LOO-PIT")))
  expect_length(vector$messages, 1L)
  expect_equal(vector$value, table$value[["z"]], tolerance = 1e-12)
  expect_length(count_messages(residuals(fit, type = "outcome"))$messages, 0L)

  # independent oracle: the PSIS-weighted normal CDF of the approximate
  # effect sizes and standard errors
  setup   <- .estimate_likelihood_setup.brma(fit)
  weights <- .diagnostic_psis_context(fit)[["psis_weights"]]
  S       <- nrow(setup[["mu"]])
  yi      <- matrix(setup[["yi"]], S, 8L, byrow = TRUE)
  sei     <- matrix(setup[["sei"]], S, 8L, byrow = TRUE)
  cdf     <- stats::pnorm(yi, setup[["mu"]], sqrt(setup[["tau_within"]]^2 + sei^2))
  expect_equal(table$value[["z"]], stats::qnorm(colSums(weights * cdf)), tolerance = 1e-10)

  # plots label the approximation instead of messaging
  qq <- count_messages(suppressWarnings(qqnorm(fit, envelope = FALSE, as_data = TRUE)))
  expect_length(qq$messages, 0L)
  expect_identical(qq$value[["ylab"]], "Sample Quantiles (Descriptive Normal Approximation)")
  funnel_data <- count_messages(suppressWarnings(funnel(fit, residual = TRUE, as_data = TRUE)))
  expect_length(funnel_data$messages, 0L)
  expect_identical(funnel_data$value[["xlab"]], "Residual (Descriptive Normal Approximation)")

  # influence() reports the approximate rstudent and the PIT-free columns
  influence_result <- count_messages(suppressWarnings(influence(fit)))
  expect_length(influence_result$messages, 1L)
  expect_named(influence_result$value[["inf"]], c("rstudent", "cov.r", "tau.del"))
})

test_that("as_draws() returns the parameters shown by summary() and hides the intermediates", {

  catalog <- .brma_auxiliary_variable_catalog(
    structure(list(fit = NULL), class = "brma"),
    c("mu", "tau", "pi[1]", "phi[2]", "theta[3]", "gamma[4]", "eta[1]", "omega[1]")
  )
  expect_identical(catalog[["variable"]], c("pi[1]", "phi[2]", "theta[3]", "gamma[4]", "eta[1]"))

  skip_if_not_installed("metadat")
  skip_if_not_installed("posterior")

  glmm <- .surface_glmm_fit()
  shown <- rownames(summary(glmm, include_mcmc_diagnostics = FALSE)[["estimates"]])
  expect_identical(shown, c("mu", "tau"))
  expect_identical(colnames(as_draws_matrix(glmm)), shown)
  auxiliary <- colnames(as_draws_matrix(glmm, include_auxiliary = TRUE))
  expect_identical(head(auxiliary, 2L), shown)
  expect_true(all(c("pi[1]", "theta[1]") %in% auxiliary))
  expect_setequal(
    sub("[[].*$", "", setdiff(auxiliary, shown)), c("pi", "theta")
  )

  multilevel <- .surface_fit("multilevel")
  shown <- rownames(summary(multilevel, include_mcmc_diagnostics = FALSE)[["estimates"]])
  expect_identical(shown, c("mu", "tau", "rho"))
  expect_identical(colnames(as_draws_matrix(multilevel)), shown)
  expect_true(all(
    paste0("gamma[", 1:4, "]") %in%
      colnames(as_draws_matrix(multilevel, include_auxiliary = TRUE))
  ))

  ordinary <- .surface_fit("ordinary")
  expect_identical(colnames(as_draws_matrix(ordinary)), c("mu", "tau"))
  expect_identical(colnames(as_draws_matrix(ordinary, include_auxiliary = TRUE)), c("mu", "tau"))
})

test_that("parameter lookups of plot(), plot_prior(), and print_prior() use the released wording", {

  ordinary   <- .surface_fit("ordinary")
  mods       <- .surface_fit("mods")
  scale      <- .surface_fit("scale")
  multilevel <- .surface_fit("multilevel")

  lookup <- function(expr, message, class = NULL) {
    error <- tryCatch(expr, error = identity)
    expect_identical(conditionMessage(error), message)
    if (!is.null(class)) {
      expect_s3_class(error, class)
    }
    invisible(error)
  }
  not_found <- "BayesTools_parameter_not_found"

  # plot() describes a fitted object by its fitted structure
  lookup(plot(ordinary, parameter = "bogus"),
    "Unknown parameter 'bogus'. Valid base parameters are: 'mu', 'tau'.", not_found)
  lookup(plot(multilevel, parameter = "bogus"),
    "Unknown parameter 'bogus'. Valid base parameters are: 'mu', 'tau', 'rho'.")
  lookup(plot(ordinary, parameter = "rho"),
    "The 'rho' parameter is only available for multilevel models.", not_found)
  lookup(plot(ordinary, parameter = "PET"),
    "The 'PET' parameter is only available for PET models.")
  lookup(plot(ordinary, parameter = "PEESE"),
    "The 'PEESE' parameter is only available for PEESE models.")
  lookup(plot(ordinary, parameter = "weightfunction"),
    "The 'omega'/'weightfunction' parameter is only available for selection models.")
  lookup(plot(ordinary, parameter_mods = "x"),
    "The 'parameter_mods' argument can only be used for models with moderators.", not_found)
  lookup(plot(mods, parameter_scale = "x"),
    "The 'parameter_scale' argument can only be used for location-scale models.")
  lookup(plot(mods, parameter_mods = "bogus"),
    "The specified 'parameter_mods' term 'bogus' is not available. Available terms are: 'intercept', 'x'.")
  lookup(plot(scale, parameter_scale = "bogus"),
    "The specified 'parameter_scale' term 'bogus' is not available. Available terms are: 'intercept', 'x'.")

  # the prior functions describe an object by its stored priors
  unfitted <- brma(yi = yi, sei = sei, data = .surface_data(), measure = "GEN",
    mods = ~ x, scale = ~ x, prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  pet <- bPET(yi = yi, sei = sei, data = .surface_data(), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  ordinary_priors <- brma(yi = yi, sei = sei, data = .surface_data(), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  for (object in list(ordinary, ordinary_priors)) {
    lookup(print_prior(object, parameter = "bogus", silent = TRUE),
      "Unknown prior parameter 'bogus'. Available parameters are: 'mu', 'tau'.")
    lookup(plot_prior(object, parameter = "weightfunction"),
      "Unknown prior parameter 'omega'. Available parameters are: 'mu', 'tau'.")
  }
  lookup(print_prior(mods, parameter = "bogus", silent = TRUE),
    paste0("Unknown prior parameter 'bogus'. Available parameters are: 'tau', ",
           "'parameter_mods = 'intercept'', 'parameter_mods = 'x''."))
  lookup(plot_prior(mods, parameter_scale = "x"),
    "The 'parameter_scale' argument can only be used when the object contains the corresponding priors.")
  lookup(print_prior(mods, parameter_mods = "bogus", silent = TRUE),
    "The specified 'parameter_mods' term 'bogus' is not available. Available terms are: 'intercept', 'x'.")
  for (object in list(scale, unfitted)) {
    lookup(plot_prior(object, parameter = "x"),
      paste0("The term 'x' is available in both moderator and scale priors. ",
             "Use 'parameter_mods' or 'parameter_scale'."))
    lookup(print_prior(object, parameter = "intercept", silent = TRUE),
      paste0("The term 'intercept' is available in both moderator and scale priors. ",
             "Use 'parameter_mods' or 'parameter_scale'."))
  }
  lookup(print_prior(scale, parameter = "x", silent = TRUE),
    paste0("The term 'x' is available in both moderator and scale priors. ",
           "Use 'parameter_mods' or 'parameter_scale'."), "BayesTools_parameter_ambiguous")
  lookup(print_prior(pet, parameter = "omega", silent = TRUE),
    "The publication-bias prior does not contain a 'omega' component.")
  lookup(plot_prior(pet, parameter = "PEESE"),
    "The publication-bias prior does not contain a 'PEESE' component.")
  lookup(plot(unfitted, parameter = "bogus"),
    paste0("Unknown prior parameter 'bogus'. Available parameters are: ",
           "'parameter_mods = 'intercept'', 'parameter_mods = 'x'', ",
           "'parameter_scale = 'intercept'', 'parameter_scale = 'x''."))

  # the hidden auxiliary priors keep their own message
  glmm <- brma.glmm(ai = 1:3, ci = 2:4, n1i = rep(10, 3), n2i = rep(12, 3),
    measure = "OR", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  expect_error(print_prior(glmm, parameter = "theta", silent = TRUE),
    "The specified parameter 'theta' is not available.", fixed = TRUE)
})

test_that("a stored LOO or WAIC of another unit names the released remedy", {

  metadata <- list(unit = "cluster", retained_context = "none", target = "loo", data_hash = "hash")
  expect_error(
    RoBMA:::.check_cached_predictive_target(NULL, metadata, "estimate", "LOO"),
    "Stored LOO was computed with unit = 'cluster'. Recompute with add_loo(object, unit = 'estimate').",
    fixed = TRUE
  )
  expect_error(
    RoBMA:::.check_cached_predictive_target(NULL, metadata, "estimate", "WAIC"),
    "Stored WAIC was computed with unit = 'cluster'. Recompute with add_waic(object, unit = 'estimate').",
    fixed = TRUE
  )
})
