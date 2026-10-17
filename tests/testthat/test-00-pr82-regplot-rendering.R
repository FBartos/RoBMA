test_that("regplot clips bands without discarding vertices", {

  testthat::skip_if_not_installed("ggplot2")
  band <- data.frame(
    x = c(-2, 0, 2, 2, 0, -2), y = c(-1, 0, 1, 3, 2, 1),
    lower = c(-1, 0, 1, 1, 0, -1), upper = c(1, 2, 3, 3, 2, 1),
    group = "all", group_id = 1
  )
  data <- list(
    points = data.frame(x = 0, y = 1, size = 1, group = "all", group_id = 1),
    pred = NULL, ci = band, pi = NULL, si = NULL, refline = NULL,
    mod_type = "continuous", groups = "all", xlim = c(-1, 1),
    ylim = c(0, 2), xlab = "x", ylab = "y"
  )
  plot <- .regplot_plot_ggplot(data, .set_dots_regplot())
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[1L]]$x, band$x)
  expect_equal(built$data[[1L]]$y, band$y)
  expect_identical(plot$coordinates$limits$x, data$xlim)
  expect_identical(plot$coordinates$limits$y, data$ylim)
})

test_that("regplot supports zero standard errors in cosmetic point sizes", {

  expect_equal(.regplot_normalized_precision(c(0, 0.1, 0.2)), c(1, 0, 0))
  expect_equal(.regplot_normalized_precision(c(0, 0)), c(0, 0))
  expect_error(.regplot_normalized_precision(c(-0.1, 0.2)), "non-negative", fixed = TRUE)
})

test_that("explicit grouped line colors override the automatic palette", {

  dots <- .set_dots_regplot(lcol = "red")
  expect_identical(
    .regplot_palette(c("a", "b"), dots$lcol, attr(dots, "lcol_supplied")),
    c(a = "red", b = "red")
  )
  defaults <- .set_dots_regplot()
  expect_length(unique(.regplot_palette(c("a", "b"), defaults$lcol,
                                         attr(defaults, "lcol_supplied"))), 2L)
})

test_that("new regplot bands refuse fitted known-R multipliers at every grid size", {

  dat <- data.frame(yi = c(0.1, 0.2, 0.3), x = c(0, 1, 2), estimate = factor(c("e1", "e2", "e3")))
  kernel <- diag(c(4, 9, 16))
  dimnames(kernel) <- list(levels(dat$estimate), levels(dat$estimate))
  normal <- BayesTools::prior("normal", list(0, 1))
  object <- brma.mv(
    yi = yi, V = diag(0.01, 3), mods = ~ x, random = ~ 1 | estimate,
    R = kernel, Rscale = "none", data = dat, measure = "GEN",
    prior_unit_information_sd = 1, prior_effect = normal, prior_mods = normal,
    prior_heterogeneity = BayesTools::prior("point", list(0.2)),
    standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE
  )
  args <- .create_jags_formula_args(object[["data"]], object[["priors"]])
  compiled <- BayesTools::JAGS_formula(
    formula = args$formula_list$mu, parameter = "mu", data = args$formula_data_list$mu,
    prior_list = args$formula_prior_list$mu, formula_scale = args$formula_scale_list$mu,
    prior_random = args$formula_random_prior_list$mu,
    random_effects_compile = args$formula_random_effects_compile_list$mu
  )
  columns <- unique(unlist(lapply(names(compiled$prior_list), function(parameter) {
    BayesTools:::.prior_linear_prior_columns(parameter, compiled$prior_list[[parameter]])
  })))
  samples <- matrix(0.2, 24, length(columns), dimnames = list(NULL, columns))
  samples[, "mu_intercept"] <- seq(-0.1, 0.1, length.out = 24)
  samples[, "mu_x"] <- seq(0.05, 0.15, length.out = 24)
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = 24L), class = c("runjags", "list"))
  attr(fit, "formula_design") <- list(mu = compiled$formula_design)
  attr(fit, "formula_scale") <- list(mu = compiled$formula_scale)
  attr(fit, "prior_list") <- c(compiled$prior_list, .create_fit_priors(object$data, object$priors))
  object[["fit"]] <- as_bayestools_fit(fit)
  expect_identical(class(object), c("only_priors.brma", "brma.mv", "brma.norm", "brma"))
  terms <- .data_marginalized_random_effects(object$data)
  variance <- .evaluate_marginalized_random_variance(object$data, samples,
    sd_evaluator = .marginalized_random_sd_evaluator(object$fit, terms))
  expect_equal(unname(variance[1L, ]), 0.2 ^ 2 * c(4, 9, 16), tolerance = 1e-14)
  message <- paste0(
    "Prediction and sampling bands ('pi'/'si') are unavailable for marginalized known-R blocks ",
    "because their row multipliers are tied to the fitted layout. Support for a supplied 'R_new' ",
    "is deferred. Set 'pi = FALSE' and 'si = FALSE'."
  )
  for (n in c(3L, 4L)) {
    at <- seq(0, by = 0.5, length.out = n)
    expect_type(regplot(object, mod = "x", at = at, pi = FALSE, si = FALSE, as_data = TRUE), "list")
    expect_error(regplot(object, mod = "x", at = at, pi = TRUE, as_data = TRUE), message, fixed = TRUE)
    expect_error(regplot(object, mod = "x", at = at, si = TRUE, as_data = TRUE), message, fixed = TRUE)
  }
  expect_error(predict(object, newdata = dat, type = "estimate"), "Support for a supplied 'R_new' is deferred.", fixed = TRUE)
})


test_that("regplot jitter preserves the caller's random stream", {

  withr::local_seed(200)
  before <- .Random.seed
  kind <- RNGkind()
  first <- .regplot_jitter_values(10, 0.2)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  expect_identical(.regplot_jitter_values(10, 0.2), first)
  rm(".Random.seed", envir = .GlobalEnv)
  .regplot_jitter_values(10, 0.2)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
