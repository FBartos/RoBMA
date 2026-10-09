context("Cached ordinary-lane regression fits")

source(testthat::test_path("common-functions.R"))

.lane_random_correlation_builder <- function(heterogeneity = NULL) {

  set.seed(2)
  k   <- 30L
  dat <- data.frame(
    study = rep(sprintf("s%02d", seq_len(10L)), each = 3L),
    x     = stats::rnorm(k)
  )
  dat[["yi"]] <- 0.2 + 0.1 * dat[["x"]] +
    stats::rnorm(10L, 0, 0.1)[as.integer(factor(dat[["study"]]))] +
    stats::rnorm(k, 0, 0.15)
  args <- list(
    yi = quote(yi), V = diag(rep(0.0225, k)), random = ~ us(1 + x | study),
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    chains = 1, sample = 300, burnin = 100, adapt = 100, seed = 1,
    silent = TRUE
  )
  args[["prior_heterogeneity"]] <- heterogeneity
  fit <- suppressWarnings(do.call(brma.mv, args))
  list(fit = fit, info = list(data = dat, arguments = args))
}

test_that("carrier integration fits are cached", {

  skip_on_cran()
  keys <- c("lane_carrier_treatment", "lane_carrier_meandif", "lane_carrier_nonzero", "lane_carrier_averaged", "lane_carrier_ordinary", "lane_carrier_gated")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

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
  save_fit("lane_carrier_treatment", fits$treatment, info = list(data = data))
  save_fit("lane_carrier_meandif", fits$meandif, info = list(data = data))
  save_fit("lane_carrier_nonzero", fits$nonzero, info = list(data = data))
  save_fit("lane_carrier_averaged", fits$averaged, info = list(data = data))
  save_fit("lane_carrier_ordinary", fits$ordinary, info = list(data = data))
  save_fit("lane_carrier_gated", fits$gated, info = list(data = data))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("public contribution multiplier fits are cached", {

  skip_on_cran()
  keys <- c("lane_contribution_constant2", "lane_contribution_constant1", "lane_contribution_tau")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  fits <- list()
  data <- data.frame(g = factor(rep(c("A", "B"), 3)),
    yi = c(-.2, .3, 0, .4, .1, .2), sei = rep(.5, 6))
  for(multiplier in list(2, 1, "tau")){
    slope <- BayesTools::prior_factor("normal", list(0, 1), contrast = "treatment")
    attr(slope, "multiply_by") <- multiplier
    fit <- expect_warning(brma(yi = yi, sei = sei, mods = ~ g, data = data, measure = "GEN",
      prior_effect = BayesTools::prior("normal", list(0, 1)), prior_mods = list(g = slope),
      prior_heterogeneity = BayesTools::prior("normal", list(0, 1), truncation = list(lower = 0)),
      chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 281, silent = TRUE),
      "Only one chain was run", fixed = TRUE)
    key <- if (is.numeric(multiplier)) paste0("lane_contribution_constant", multiplier) else "lane_contribution_tau"
    fits[[key]] <- list(fit = fit, info = list(data = data, slope = slope, multiplier = multiplier))
  }

  save_fit("lane_contribution_constant2", fits[["lane_contribution_constant2"]]$fit,
    info = fits[["lane_contribution_constant2"]]$info)
  save_fit("lane_contribution_constant1", fits[["lane_contribution_constant1"]]$fit,
    info = fits[["lane_contribution_constant1"]]$info)
  save_fit("lane_contribution_tau", fits[["lane_contribution_tau"]]$fit,
    info = fits[["lane_contribution_tau"]]$info)
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("actual brma coefficients and static-method routes use the compiled owner fit is cached", {

  skip_on_cran()
  keys <- c("lane_contribution_point")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  slope <- BayesTools::prior("point", list(2))
  attr(slope, "multiply_by") <- "tau"
  data <- data.frame(x = c(10, 20, 30), yi = c(5, 5, 5), sei = rep(1, 3))
  fit <- expect_warning(brma(yi = yi, sei = sei, mods = ~ x, data = data, measure = "GEN",
    prior_effect = BayesTools::prior("point", list(5)), prior_mods = list(x = slope),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1), truncation = list(lower = 0)),
    chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 281, silent = TRUE),
    "Only one chain was run", fixed = TRUE)
  save_fit("lane_contribution_point", fit, info = list(data = data, slope = slope))
  expect_true(is_cached_fit_valid("lane_contribution_point"))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("scalar discovery retains declared refusal rows and usable sibling quantities fit is cached", {

  skip_on_cran()
  keys <- c("lane_contribution_discovery")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  slope <- BayesTools::prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "tau"
  data <- data.frame(x = c(10, 20, 30), y = c(0, 1, 2), yi = c(.1, .2, .3), sei = .5)
  object <- expect_warning(brma(yi = yi, sei = sei, mods = ~ x * y, data = data, measure = "GEN",
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_mods = list(x = slope, y = BayesTools::prior("normal", list(0, 1)), `x:y` = BayesTools::prior("normal", list(0, 1))),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1), list(lower = 0)),
    chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 281, silent = TRUE),
    "Only one chain was run", fixed = TRUE)
  save_fit("lane_contribution_discovery", object, info = list(data = data, slope = slope))
  expect_true(is_cached_fit_valid("lane_contribution_discovery"))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("exp-affine coefficients refuse qCMDE/IWMDE with the original-scale classes fit is cached", {

  skip_on_cran()
  keys <- c("lane_exp_affine")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  # The original-scale scale intercept of a scale formula with a standardized
  # continuous predictor is an exp(affine) map of the fitted coefficients.
  # plot() and hypothesis() refuse qCMDE/IWMDE for it with the same
  # density-method classes; hypothesis() keeps its own classes first.
  set.seed(1)
  k   <- 30L
  dat <- data.frame(x = stats::rnorm(k), sei = stats::runif(k, 0.1, 0.3))
  dat[["yi"]] <- stats::rnorm(k, 0.2 + 0.1 * dat[["x"]], dat[["sei"]])
  fit <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ x, scale = ~ x, data = dat,
    measure = "SMD", chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))
  save_fit("lane_exp_affine", fit, info = list(data = dat))
  expect_true(is_cached_fit_valid("lane_exp_affine"))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("fixed reference-level statements keep one reason under qCMDE for model-averaged fits fit is cached", {

  skip_on_cran()
  keys <- c("lane_factor_averaged_reference")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  set.seed(1)
  k    <- 48L
  data <- data.frame(
    g1  = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
    sei = stats::runif(k, 0.1, 0.3)
  )
  data[["yi"]] <- stats::rnorm(
    k, c(0, 0.2, 0.4)[as.integer(data[["g1"]])], data[["sei"]]
  )
  fit <- suppressWarnings(BMA.norm(
    yi = yi, sei = sei, mods = ~ g1, data = data, measure = "SMD",
    set_contrast_factor_predictors = "treatment",
    chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))
  save_fit("lane_factor_averaged_reference", fit, info = list(data = data))
  expect_true(is_cached_fit_valid("lane_factor_averaged_reference"))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("ordered factor labels have cached original fits", {

  skip_on_cran()
  keys <- c("lane_factor_ordered_labels", "lane_factor_numeric_labels")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  set.seed(3)
  k    <- 48L
  data <- data.frame(
    g   = factor(rep(c("lo", "mid", "hi"), length.out = k), levels = c("lo", "mid", "hi")),
    sei = stats::runif(k, 0.1, 0.3)
  )
  data[["yi"]] <- stats::rnorm(k, c(0, 0.2, 0.3)[as.integer(data[["g"]])], data[["sei"]])
  ordered <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ g, data = data, measure = "SMD",
    prior_mods = list(g = BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)))),
    chains = 1, sample = 500, burnin = 100, adapt = 100, seed = 1, silent = TRUE
  ))
  save_fit("lane_factor_ordered_labels", ordered, info = list(data = data))
  data[["g"]] <- factor(
    c("1", "2", "3")[as.integer(data[["g"]])],
    levels = c("1", "2", "3")
  )
  numeric_labels <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ g, data = data, measure = "SMD",
    prior_mods = list(g = BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)))),
    chains = 1, sample = 500, burnin = 100, adapt = 100, seed = 1, silent = TRUE
  ))
  save_fit("lane_factor_numeric_labels", numeric_labels, info = list(data = data))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("scalar eligibility fits are cached", {

  skip_on_cran()
  keys <- c("lane_scalar_standardized", "lane_scalar_raw", "lane_scalar_convolution", "lane_scalar_averaged", "lane_scalar_averaged_default")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  set.seed(1)
  k    <- 40L
  data <- data.frame(
    x   = stats::rnorm(k, 1, 2),
    sei = stats::runif(k, 0.1, 0.3)
  )
  data[["yi"]] <- stats::rnorm(k, 0.1 * data[["x"]], data[["sei"]])
  cauchy <- BayesTools::prior("cauchy", list(0, 0.5))
  fits <- lapply(c(standardized = TRUE, raw = FALSE), function(standardize) {
    suppressWarnings(brma(
      yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
      prior_effect = cauchy, prior_mods = list(x = cauchy),
      standardize_continuous_predictors = standardize,
      chains = 1, sample = 1000, burnin = 200, adapt = 100,
      seed = 1, silent = TRUE
    ))
  })
  fits[["convolution"]] <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
    prior_mods = list(x = cauchy),
    chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))
  fits[["averaged"]] <- suppressWarnings(BMA.norm(
    yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
    prior_effect = cauchy, prior_mods = list(x = cauchy),
    chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))
  fits[["averaged_default"]] <- suppressWarnings(BMA.norm(
    yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
    chains = 1, sample = 1000, burnin = 200, adapt = 100,
    seed = 1, silent = TRUE
  ))
  save_fit("lane_scalar_standardized", fits[["standardized"]], info = list(data = data, cauchy = cauchy))
  save_fit("lane_scalar_raw", fits[["raw"]], info = list(data = data, cauchy = cauchy))
  save_fit("lane_scalar_convolution", fits[["convolution"]], info = list(data = data, cauchy = cauchy))
  save_fit("lane_scalar_averaged", fits[["averaged"]], info = list(data = data, cauchy = cauchy))
  save_fit("lane_scalar_averaged_default", fits[["averaged_default"]], info = list(data = data, cauchy = cauchy))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("fitted known-V row SD sources read the coordinates of factor scale terms fit is cached", {

  skip_on_cran()
  keys <- c("lane_known_v_factor_scale")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  # A scale factor with two design columns is sampled as the coordinates
  # 'log_tau_g[1]' and 'log_tau_g[2]' of the node 'log_tau_g'; the values
  # functions read the coordinates, as the draws name them.
  set.seed(3)
  k   <- 12L
  dat <- data.frame(
    yi  = stats::rnorm(k, 0.2, 0.3),
    obs = factor(sprintf("e%02d", seq_len(k))),
    x   = stats::rnorm(k),
    g   = factor(rep(c("a", "b", "c"), length.out = k))
  )
  fit <- suppressWarnings(brma.mv(
    yi = yi, V = diag(rep(0.04, k)), random = ~ 1 | obs, scale = ~ x + g,
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    chains = 1, sample = 200, burnin = 100, adapt = 100, seed = 1,
    silent = TRUE
  ))
  save_fit("lane_known_v_factor_scale", fit, info = list(data = dat))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("random correlation fits are cached", {

  skip_on_cran()
  keys <- c("lane_random_correlation", "lane_random_correlation_direct", "lane_random_correlation_gated")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  default <- .lane_random_correlation_builder()
  save_fit("lane_random_correlation", default$fit, info = default$info)
  direct <- .lane_random_correlation_builder(BayesTools::prior_random(
    sd = BayesTools::prior("normal", list(0, .5), list(0, Inf))
  ))
  save_fit("lane_random_correlation_direct", direct$fit, info = direct$info)
  set.seed(2)
  k   <- 30L
  dat <- data.frame(
    study = rep(sprintf("s%02d", seq_len(10L)), each = 3L),
    x     = stats::rnorm(k)
  )
  dat[["yi"]] <- 0.2 + 0.1 * dat[["x"]] + stats::rnorm(k, 0, 0.15)
  fit <- suppressWarnings(BMA.mv(
    yi = yi, V = diag(rep(0.0225, k)), random = ~ us(1 + x | study),
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    chains = 1, sample = 300, burnin = 100, adapt = 100, seed = 1,
    silent = TRUE
  ))
  save_fit("lane_random_correlation_gated", fit, info = list(data = dat))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("nested gated allocation fit is cached", {

  skip_on_cran()
  keys <- c("lane_nested_total")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  set.seed(3)
  k   <- 12L
  dat <- data.frame(
    study = factor(rep(sprintf("s%02d", seq_len(6L)), each = 2L)),
    esid  = factor(seq_len(k))
  )
  dat[["yi"]] <- 0.2 + stats::rnorm(6L, 0, 0.12)[as.integer(dat[["study"]])] +
    stats::rnorm(k, 0, 0.1)
  fit <- suppressWarnings(BMA.mv(
    yi = yi, V = diag(rep(0.0225, k)), random = ~ 1 | study / esid,
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    chains = 1, sample = 400, burnin = 150, adapt = 100, seed = 1,
    silent = TRUE
  ))
  save_fit("lane_nested_total", fit, info = list(data = dat))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("selection reference weights are structural convergence parameters fit is cached", {

  skip_on_cran()
  keys <- c("lane_selection_reference")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  prior_bias <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_cumulative(c(1, 1, 1))
  )

  fit <- suppressWarnings(bselmodel(
    yi                 = c(.1, .2, .3, .05, .25),
    sei                = c(.1, .1, .1, .12, .08),
    measure            = "SMD",
    prior_bias         = prior_bias,
    chains             = 2,
    sample             = 200,
    burnin             = 100,
    adapt              = 100,
    seed               = 1,
    silent             = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = NULL
    )
  ))
  save_fit("lane_selection_reference", fit, info = list(prior_bias = prior_bias))
  expect_true(is_cached_fit_valid("lane_selection_reference"))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("a fitted point-baserate GLMM supports log-likelihood criteria fit is cached", {

  skip_on_cran()
  keys <- c("lane_glmm_point_baserate")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  prior_baserate <- BayesTools::prior_factor(
    "point",
    list(location = 0.05),
    contrast = "independent"
  )
  expect_warning(
    fit <- brma.glmm(
      ai = c(4L, 6L, 3L), bi = c(119L, 300L, 228L),
      ci = c(11L, 29L, 11L), di = c(128L, 274L, 209L),
      measure = "OR", prior_baserate = prior_baserate,
      chains = 1L, sample = 100L, burnin = 100L, adapt = 100L,
      seed = 42, silent = TRUE
    ),
    paste0("Only one chain was run. R-hat cannot be computed; ",
           "checking the remaining enabled convergence criteria."),
    fixed = TRUE
  )
  save_fit("lane_glmm_point_baserate", fit, info = list(prior_baserate = prior_baserate))
  expect_true(is_cached_fit_valid("lane_glmm_point_baserate"))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("boundary-concentrated beta priors support likelihood criteria fit is cached", {

  skip_on_cran()
  keys <- c("lane_glmm_boundary_beta")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  expect_warning(
    fit <- brma.glmm(
      ai             = c(0L, 0L, 1L),
      bi             = c(1L, 10L, 9L),
      ci             = c(0L, 1L, 0L),
      di             = c(1L, 9L, 10L),
      measure        = "OR",
      prior_baserate = BayesTools::prior("beta", list(0.1, 0.1)),
      chains         = 1L,
      sample         = 100L,
      burnin         = 100L,
      adapt          = 100L,
      seed           = 42L,
      silent         = TRUE
    ),
    paste0("Only one chain was run. R-hat cannot be computed; ",
           "checking the remaining enabled convergence criteria."),
    fixed = TRUE
  )
  save_fit("lane_glmm_boundary_beta", fit, info = NULL)
  expect_true(is_cached_fit_valid("lane_glmm_boundary_beta"))
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})

test_that("fitted truncated-nuisance GLMMs support likelihood criteria fits are cached", {

  skip_on_cran()
  keys <- c("lane_glmm_truncated_beta", "lane_glmm_truncated_lograte")
  if (!any(keys %in% active_fit_catalog()$name)) skip("No ordinary-lane fits in this group are active.")
  refit <- is_true_env("ROBMA_TEST_FORCE_REFIT") || is_false_env("ROBMA_TEST_SKIP_REFIT")
  if (!refit && all(vapply(keys, is_cached_fit_valid, logical(1)))) skip("Valid cached fits already exist.")

  expect_warning(
    fit_bin <- brma.glmm(
      ai = c(4L, 6L, 3L), bi = c(119L, 300L, 228L),
      ci = c(11L, 29L, 11L), di = c(128L, 274L, 209L),
      measure = "OR",
      prior_baserate = BayesTools::prior(
        "beta", list(1.5, 2.5), truncation = list(0.1, 0.9)
      ),
      chains = 1L, sample = 100L, burnin = 100L, adapt = 100L,
      seed = 123L, silent = TRUE
    ),
    paste0("Only one chain was run. R-hat cannot be computed; ",
           "checking the remaining enabled convergence criteria."),
    fixed = TRUE
  )
  expect_warning(
    fit_pois <- brma.glmm(
      x1i = c(0L, 3L, 7L), x2i = c(0L, 4L, 11L),
      t1i = c(20, 15, 18), t2i = c(22, 19, 21),
      measure = "IRR",
      prior_lograte = BayesTools::prior(
        "normal", list(-1, 1.3), truncation = list(-3, 2)
      ),
      chains = 1L, sample = 100L, burnin = 100L, adapt = 100L,
      seed = 123L, silent = TRUE
    ),
    paste0("Only one chain was run. R-hat cannot be computed; ",
           "checking the remaining enabled convergence criteria."),
    fixed = TRUE
  )
  save_fit("lane_glmm_truncated_beta", fit_bin)
  save_fit("lane_glmm_truncated_lograte", fit_pois)
  expect_true(all(vapply(keys, is_cached_fit_valid, logical(1))))
})
