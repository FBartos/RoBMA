context("Exact marginal likelihood for fixed-zero random effects")

skip_on_cran()
skip_if_not_installed("bridgesampling")

.fixed_zero_marglik_settings <- function() {

  list(
    chains             = 2,
    sample             = 100,
    burnin             = 50,
    adapt              = 100,
    seed               = 1,
    silent             = TRUE,
    convergence_checks = set_convergence_checks(
      max_Rhat = NULL,
      min_ESS  = NULL
    )
  )
}

.expect_exact_fixed_zero_marglik <- function(fit, expected) {

  fit    <- add_marglik(fit)
  bridge <- bridge_sampler(fit)

  expect_identical(bridge[["aggregation"]][["rule"]], "exact_zero_dimensional")
  expect_identical(bridge[["aggregation"]][["n_repetitions"]], 0L)
  expect_equal(logml(fit), expected, tolerance = 1e-12)
}

test_that("fixed-zero cluster and random formulas have exact marginal likelihoods", {

  settings <- .fixed_zero_marglik_settings()
  point_zero <- BayesTools::prior(
    "spike",
    parameters = list(location = 0)
  )
  dat <- data.frame(
    yi        = c(0.10, 0.20, 0.30, 0.40),
    study     = c("s1", "s1", "s2", "s2"),
    criterion = c("c1", "c2", "c1", "c2"),
    estimate  = paste0("e", 1:4)
  )
  sei      <- rep(0.20, nrow(dat))
  expected <- sum(stats::dnorm(dat[["yi"]], mean = 0, sd = sei, log = TRUE))

  cluster_fit <- brma(
    yi                        = dat[["yi"]],
    sei                       = sei,
    cluster                   = dat[["study"]],
    measure                   = "GEN",
    prior_effect              = point_zero,
    prior_heterogeneity       = point_zero,
    prior_unit_information_sd = 1,
    chains                    = settings[["chains"]],
    sample                    = settings[["sample"]],
    burnin                    = settings[["burnin"]],
    adapt                     = settings[["adapt"]],
    seed                      = settings[["seed"]],
    silent                    = settings[["silent"]],
    convergence_checks        = settings[["convergence_checks"]]
  )
  .expect_exact_fixed_zero_marglik(cluster_fit, expected)

  for (random in list(
    ~ 1 | study,
    list(~ 1 | study, ~ 1 | criterion)
  )) {
    random_fit <- brma.mv(
      yi                        = yi,
      V                         = diag(sei^2),
      random                    = random,
      data                      = dat,
      known_v_parameterization  = "latent",
      measure                   = "GEN",
      prior_effect              = point_zero,
      prior_heterogeneity       = point_zero,
      prior_unit_information_sd = 1,
      chains                    = settings[["chains"]],
      sample                    = settings[["sample"]],
      burnin                    = settings[["burnin"]],
      adapt                     = settings[["adapt"]],
      seed                      = settings[["seed"]],
      silent                    = settings[["silent"]],
      convergence_checks        = settings[["convergence_checks"]]
    )
    .expect_exact_fixed_zero_marglik(random_fit, expected)
  }

  for (random in list(
    ~ 1 | estimate,
    ~ 1 | study / estimate
  )) {
    marginalized_fit <- brma.mv(
      yi                        = yi,
      V                         = diag(sei^2),
      random                    = random,
      data                      = dat,
      known_v_parameterization  = "block_mvn",
      measure                   = "GEN",
      prior_effect              = point_zero,
      prior_heterogeneity       = point_zero,
      prior_unit_information_sd = 1,
      chains                    = settings[["chains"]],
      sample                    = settings[["sample"]],
      burnin                    = settings[["burnin"]],
      adapt                     = settings[["adapt"]],
      seed                      = settings[["seed"]],
      silent                    = settings[["silent"]],
      convergence_checks        = settings[["convergence_checks"]]
    )
    .expect_exact_fixed_zero_marglik(marginalized_fit, expected)
  }
})


test_that("fixed-zero marginalized random effects preserve sampled fixed effects", {

  settings <- .fixed_zero_marglik_settings()
  point_zero <- BayesTools::prior(
    "spike",
    parameters = list(location = 0)
  )
  normal_effect <- BayesTools::prior(
    "normal",
    parameters = list(mean = 0, sd = 1)
  )
  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30, 0.40),
    estimate = paste0("e", 1:4)
  )

  fit <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.20^2, nrow(dat))),
    random                    = ~ 1 | estimate,
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_effect              = normal_effect,
    prior_heterogeneity       = point_zero,
    prior_unit_information_sd = 1,
    chains                    = settings[["chains"]],
    sample                    = settings[["sample"]],
    burnin                    = settings[["burnin"]],
    adapt                     = settings[["adapt"]],
    seed                      = settings[["seed"]],
    silent                    = settings[["silent"]],
    convergence_checks        = settings[["convergence_checks"]]
  )

  expect_silent(fit <- add_marglik(fit))
  expect_true(is.finite(logml(fit)))
  expect_identical(
    bridge_sampler(fit)[["aggregation"]][["rule"]],
    "median_finite_logml"
  )
})


test_that("only exact structural zero bypasses marginalized random variance", {

  near_zero <- 1e-8
  point_near_zero <- BayesTools::prior(
    "spike",
    parameters = list(location = near_zero)
  )
  dat <- data.frame(
    yi       = c(0.10, 0.20),
    estimate = c("e1", "e2")
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.20^2, nrow(dat))),
    random                    = ~ 1 | estimate,
    data                      = dat,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_heterogeneity       = point_near_zero,
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  design <- .fitted_formula_design(object, "mu", required = TRUE)
  parameter <- design[["random_effects"]][[1L]][["sd_parameter_names"]][[1L]]
  parameters <- stats::setNames(list(near_zero), parameter)

  expect_false(.marglik_formula_random_effects_fixed_zero(
    design = design,
    data   = object[["data"]]
  ))
  expect_equal(
    .marglik_known_v_extra_variance(
      parameters         = parameters,
      model_data         = object[["data"]],
      bridge_context     = NULL,
      tau_within_samples = matrix(0, nrow = 1L, ncol = nrow(dat)),
      is_random          = TRUE,
      K                  = nrow(dat),
      fixed_zero_random  = FALSE
    ),
    matrix(near_zero^2, nrow = 1L, ncol = nrow(dat))
  )
})
