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
    criterion = c("c1", "c2", "c1", "c2")
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
})
