context("Prior input handling for brma.glmm")

skip_on_cran()
source(testthat::test_path("helper-contracts.R"))

test_data_bin <- data.frame(
  ai = c(4L, 6L, 8L),
  bi = c(16L, 14L, 12L),
  ci = c(3L, 5L, 7L),
  di = c(17L, 15L, 13L)
)

test_data_pois <- data.frame(
  x1i = c(4L, 6L, 8L),
  x2i = c(3L, 5L, 7L),
  t1i = c(20, 25, 30),
  t2i = c(21, 24, 31)
)


test_that("Binomial GLMM baserate priors are assigned", {

  result_default <- brma.glmm(
    ai = ai, bi = bi, ci = ci, di = di,
    data = test_data_bin, measure = "OR",
    prior_baserate = NULL, only_priors = TRUE
  )[["priors"]]

  expect_equal(result_default$outcome$pi$distribution, "beta")
  expect_equal(result_default$outcome$pi$parameters$alpha, 1)
  expect_equal(result_default$outcome$pi$parameters$beta,  1)

  custom_prior <- BayesTools::prior("beta", parameters = list(alpha = 2, beta = 3))
  result_custom <- brma.glmm(
    ai = ai, bi = bi, ci = ci, di = di,
    data = test_data_bin, measure = "OR",
    prior_baserate = custom_prior, only_priors = TRUE
  )[["priors"]]

  expect_equal(result_custom$outcome$pi$parameters$alpha, 2)
  expect_equal(result_custom$outcome$pi$parameters$beta,  3)
})


test_that("Binomial GLMM baserate point priors exclude endpoints", {

  for (location in c(0, 1)) {
    prior_baserate <- BayesTools::prior(
      "spike",
      parameters = list(location = location)
    )

    expect_error(
      brma.glmm(
        ai = ai, bi = bi, ci = ci, di = di,
        data = test_data_bin, measure = "OR",
        prior_baserate = prior_baserate, only_priors = TRUE
      ),
      "baserate.*strictly within \\(0, 1\\)"
    )
  }

  expect_no_error(.assign_prior.baserate(BayesTools::prior(
    "spike",
    parameters = list(location = .Machine$double.eps)
  )))
  expect_no_error(.assign_prior.baserate(BayesTools::prior(
    "spike",
    parameters = list(location = 1 - .Machine$double.eps)
  )))
})


test_that("GLMM nuisance priors reject non-point discrete distributions only", {

  discrete <- BayesTools::prior("bernoulli", parameters = list(probability = 0.5))

  expect_error(
    .assign_prior.baserate(discrete),
    "Non-point discrete priors are not supported for 'prior_baserate'; use a point or continuous prior"
  )
  expect_error(
    .assign_prior.lograte(discrete, test_data_pois),
    "Non-point discrete priors are not supported for 'prior_lograte'; use a point or continuous prior"
  )
})


test_that("GLMM nuisance priors accept the released continuous families", {

  # 4.0.0 fitted any continuous base-rate and log-rate prior; a base-rate
  # prior is restricted to [0, 1]
  baserate <- list(
    normal    = BayesTools::prior("normal", list(0.5, 0.1), truncation = list(0, 1)),
    uniform   = BayesTools::prior("uniform", list(0, 1)),
    lognormal = BayesTools::prior("lognormal", list(0, 1)),
    gamma     = BayesTools::prior("gamma", list(shape = 2, rate = 1)),
    beta      = BayesTools::prior("beta", list(2, 2))
  )
  for (name in names(baserate)) {
    fitted_prior <- suppressWarnings(brma.glmm(
      ai = ai, bi = bi, ci = ci, di = di,
      data = test_data_bin, measure = "OR",
      prior_baserate = baserate[[name]], only_priors = TRUE
    ))[["priors"]][["outcome"]][["pi"]]
    expect_identical(fitted_prior[["distribution"]], baserate[[name]][["distribution"]],
                     info = name)
    expect_lte(fitted_prior[["truncation"]][["upper"]], 1)
    expect_gte(fitted_prior[["truncation"]][["lower"]], 0)
  }

  lograte <- list(
    normal    = BayesTools::prior("normal", list(0, 1)),
    cauchy    = BayesTools::prior("cauchy", list(0, 1)),
    t         = BayesTools::prior("t", list(0, 1, 3)),
    lognormal = BayesTools::prior("lognormal", list(0, 1))
  )
  for (name in names(lograte)) {
    fitted_prior <- BMA.glmm(
      x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
      data = test_data_pois, measure = "IRR",
      prior_lograte = lograte[[name]], only_priors = TRUE
    )[["priors"]][["outcome"]][["phi"]]
    expect_identical(fitted_prior[["distribution"]], lograte[[name]][["distribution"]],
                     info = name)
  }
})


# Independent oracles of the estimate-unit GLMM likelihood: the nuisance
# parameter and the standard-normal latent effect are integrated by nested
# stats::integrate() calls (the base rate on the logit scale, dp = p (1 - p) du).
.glmm_oracle_binom <- function(a, n1, c, n2, mu, tau, density, lower = 0, upper = 1) {

  inner <- function(theta) {
    vapply(theta, function(th) {
      effect <- mu + tau * th
      stats::integrate(function(u) {
        p <- stats::plogis(u)
        stats::dbinom(a, n1, stats::plogis(u + 0.5 * effect)) *
          stats::dbinom(c, n2, stats::plogis(u - 0.5 * effect)) *
          density(p) * p * (1 - p)
      }, stats::qlogis(lower), stats::qlogis(upper), subdivisions = 5000L,
      rel.tol = 1e-10, abs.tol = 0)[["value"]] * stats::dnorm(th)
    }, numeric(1))
  }

  # the standard-normal mass outside (-10, 10) is 1.5e-23; an infinite range would
  # evaluate the inner integral at effects that underflow every likelihood factor
  return(log(stats::integrate(inner, -10, 10, subdivisions = 5000L,
                              rel.tol = 1e-9, abs.tol = 0)[["value"]]))
}

.glmm_oracle_pois <- function(x1, t1, x2, t2, mu, tau, density) {

  inner <- function(theta) {
    vapply(theta, function(th) {
      effect <- mu + tau * th
      stats::integrate(function(phi) {
        stats::dpois(x1, exp(phi + log(t1) + 0.5 * effect)) *
          stats::dpois(x2, exp(phi + log(t2) - 0.5 * effect)) * density(phi)
      }, -Inf, Inf, subdivisions = 5000L, rel.tol = 1e-9, abs.tol = 0)[["value"]] *
        stats::dnorm(th)
    }, numeric(1))
  }

  # the standard-normal mass outside (-10, 10) is 1.5e-23; an infinite range would
  # evaluate the inner integral at effects that underflow every likelihood factor
  return(log(stats::integrate(inner, -10, 10, subdivisions = 5000L,
                              rel.tol = 1e-9, abs.tol = 0)[["value"]]))
}

# The refined prior-CDF grid stops at a relative change of 1e-6 between
# successive rules (.glmm_grid_control()), and the oracle integrates to 1e-9;
# a relative tolerance of 1e-6 therefore separates agreement from a wrong grid.
.glmm_oracle_tolerance <- 1e-6

.glmm_released_families <- function() {

  list(
    binomial = list(
      uniform = list(
        prior = BayesTools::prior("uniform", list(0, 1)),
        density = function(p) stats::dunif(p), lower = 0, upper = 1),
      `truncated normal` = list(
        prior = BayesTools::prior("normal", list(0.3, 0.2), truncation = list(0.05, 0.95)),
        density = function(p) stats::dnorm(p, 0.3, 0.2) /
          (stats::pnorm(0.95, 0.3, 0.2) - stats::pnorm(0.05, 0.3, 0.2)),
        lower = 0.05, upper = 0.95)
    ),
    poisson = list(
      cauchy = list(
        prior = BayesTools::prior("cauchy", list(0, 1)),
        density = function(x) stats::dcauchy(x, 0, 1)),
      t = list(
        prior = BayesTools::prior("t", list(0, 1, 3)),
        density = function(x) stats::dt(x, df = 3))
    )
  )
}

test_that("the estimate-unit GLMM likelihood of released nuisance families equals its integral", {

  skip_on_cran()
  families <- .glmm_released_families()

  for (name in names(families[["binomial"]])) {
    family <- families[["binomial"]][[name]]
    stored <- brma.glmm(
      ai = ai, bi = bi, ci = ci, di = di,
      data = test_data_bin, measure = "OR",
      prior_baserate = family[["prior"]], only_priors = TRUE
    )[["priors"]][["outcome"]][["pi"]]
    value <- .outcome_pdf.binom(4L, 3L, 20L, 20L, matrix(0.2), matrix(0.4), stored)
    oracle <- .glmm_oracle_binom(4L, 20L, 3L, 20L, 0.2, 0.4, family[["density"]],
                                 family[["lower"]], family[["upper"]])
    expect_equal(as.numeric(value), oracle, tolerance = .glmm_oracle_tolerance, info = name)
    expect_equal(
      as.numeric(.outcome_pdf_sum.binom(4L, 3L, 20L, 20L, matrix(0.2), matrix(0.4), stored)),
      oracle, tolerance = .glmm_oracle_tolerance, info = name
    )
  }

  for (name in names(families[["poisson"]])) {
    family <- families[["poisson"]][[name]]
    stored <- BMA.glmm(
      x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
      data = test_data_pois, measure = "IRR",
      prior_lograte = family[["prior"]], only_priors = TRUE
    )[["priors"]][["outcome"]][["phi"]]
    value <- .outcome_pdf.pois(4L, 3L, 20, 21, matrix(0.2), matrix(0.4), stored)
    oracle <- .glmm_oracle_pois(4L, 20, 3L, 21, 0.2, 0.4, family[["density"]])
    expect_equal(as.numeric(value), oracle, tolerance = .glmm_oracle_tolerance, info = name)
    expect_equal(
      as.numeric(.outcome_pdf_sum.pois(4L, 3L, 20, 21, matrix(0.2), matrix(0.4), stored)),
      oracle, tolerance = .glmm_oracle_tolerance, info = name
    )
  }
})

test_that("a gamma log-rate prior on a low-event cell stops with the grid nonconvergence class", {

  skip_on_cran()
  # At 1200 nodes (the end of the refinement ladder) the prior-CDF grid is still
  # about 4e-6 from the nested-QUADPACK value of this cell, so it cannot certify
  # its 1e-6 tolerance. The 30-node grid of 4.0.0 returned a value 0.2 too low.
  # The documented outcome is the typed refusal, which names a remedy.
  stored <- BMA.glmm(
    x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
    data = test_data_pois, measure = "IRR",
    prior_lograte = BayesTools::prior("gamma", list(2, 1)), only_priors = TRUE
  )[["priors"]][["outcome"]][["phi"]]
  condition <- tryCatch(
    .outcome_pdf.pois(4L, 3L, 20, 21, matrix(0.2), matrix(0.4), stored),
    error = identity)
  expect_s3_class(condition, "RoBMA_glmm_grid_nonconvergence")
  expect_s3_class(condition, "RoBMA_glmm_numerical_error")
  expect_match(conditionMessage(condition), "a normal 'prior_lograte'", fixed = TRUE)
})

test_that("fitted GLMMs with released nuisance families give log_lik(), LOO, WAIC, and autocompute", {

  skip_on_cran()
  families <- .glmm_released_families()
  fit_bin <- function(prior) suppressWarnings(brma.glmm(
    ai = c(4L, 6L, 3L, 5L), bi = c(119L, 300L, 228L, 99L),
    ci = c(11L, 29L, 11L, 12L), di = c(128L, 274L, 209L, 90L),
    measure = "OR", prior_baserate = prior, prior_unit_information_sd = 1,
    chains = 2, adapt = 50, burnin = 100, sample = 100, seed = 1,
    silent = TRUE, parallel = FALSE
  ))
  fit_pois <- function(prior) suppressWarnings(brma.glmm(
    x1i = c(0L, 3L, 7L, 4L), x2i = c(0L, 4L, 11L, 6L),
    t1i = c(20, 15, 18, 22), t2i = c(22, 19, 21, 20),
    measure = "IRR", prior_lograte = prior, prior_unit_information_sd = 1,
    chains = 2, adapt = 50, burnin = 100, sample = 100, seed = 1,
    silent = TRUE, parallel = FALSE
  ))
  cells <- list(c(1L, 1L), c(60L, 2L), c(150L, 3L), c(200L, 4L))

  for (type in names(families)) {
    for (name in names(families[[type]])) {
      family <- families[[type]][[name]]
      info   <- paste(type, name)
      fit    <- if (type == "binomial") fit_bin(family[["prior"]]) else fit_pois(family[["prior"]])
      data   <- fit[["data"]][["outcome"]]
      draws  <- posterior::as_draws_matrix(fit)
      mu_draws  <- as.numeric(draws[, "mu"])
      tau_draws <- as.numeric(draws[, "tau"])

      log_lik_matrix <- log_lik(fit)
      expect_identical(dim(log_lik_matrix), c(nrow(draws), 4L), info = info)
      expect_true(all(is.finite(log_lik_matrix)), info = info)
      for (cell in cells) {
        s <- cell[[1L]]
        k <- cell[[2L]]
        oracle <- if (type == "binomial") {
          .glmm_oracle_binom(data[["ai"]][k], data[["n1i"]][k], data[["ci"]][k],
            data[["n2i"]][k], mu_draws[s], tau_draws[s], family[["density"]],
            family[["lower"]], family[["upper"]])
        } else {
          .glmm_oracle_pois(data[["x1i"]][k], data[["t1i"]][k], data[["x2i"]][k],
            data[["t2i"]][k], mu_draws[s], tau_draws[s], family[["density"]])
        }
        expect_equal(unname(log_lik_matrix[s, k]), unname(oracle),
                     tolerance = .glmm_oracle_tolerance, info = paste(info, s, k))
      }

      with_loo <- suppressWarnings(add_loo(fit))
      expect_s3_class(with_loo, "brma")
      expect_true(is.finite(loo(with_loo)[["estimates"]]["elpd_loo", "Estimate"]), info = info)
      with_waic <- suppressWarnings(add_waic(fit))
      expect_true(is.finite(waic(with_waic)[["estimates"]]["elpd_waic", "Estimate"]), info = info)
    }
  }

  # an automatic computation returns the fit, as at the 4.0.0 release
  previous <- RoBMA.get_option("autocompute.loo")
  RoBMA.options(autocompute.loo = TRUE)
  on.exit(RoBMA.options(autocompute.loo = previous), add = TRUE)
  automatic <- fit_bin(families[["binomial"]][["uniform"]][["prior"]])
  expect_s3_class(automatic, "brma")
  expect_false(is.null(automatic[["loo"]]))
  expect_true(is.finite(loo(automatic)[["estimates"]]["elpd_loo", "Estimate"]))
})


test_that("accepted GLMM nuisance priors feed the post-fit likelihood route", {

  bin_priors <- list(
    BayesTools::prior("point", list(location = 0.4)),
    BayesTools::prior("beta", list(1.5, 2.5)),
    BayesTools::prior(
      "beta", list(1.5, 2.5), truncation = list(0.1, 0.9)
    )
  )
  pois_priors <- list(
    BayesTools::prior("point", list(location = -1)),
    BayesTools::prior("normal", list(-1, 1.3)),
    BayesTools::prior(
      "normal", list(-1, 1.3), truncation = list(-3, 2)
    )
  )

  for (prior_baserate in bin_priors) {
    prior_pi <- brma.glmm(
      ai = ai, bi = bi, ci = ci, di = di,
      data = test_data_bin, measure = "OR",
      prior_baserate = prior_baserate, only_priors = TRUE
    )[["priors"]][["outcome"]][["pi"]]
    log_lik <- .outcome_pdf.binom(
      4L, 3L, 20L, 20L, matrix(0.2), matrix(0.4), prior_pi
    )
    expect_true(is.finite(log_lik[1L, 1L]))
  }

  for (prior_lograte in pois_priors) {
    prior_phi <- BMA.glmm(
      x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
      data = test_data_pois, measure = "IRR",
      prior_lograte = prior_lograte, only_priors = TRUE
    )[["priors"]][["outcome"]][["phi"]]
    log_lik <- .outcome_pdf.pois(
      4L, 3L, 20, 21, matrix(0.2), matrix(0.4), prior_phi
    )
    expect_true(is.finite(log_lik[1L, 1L]))
  }
})


test_that("Binomial nuisance quadrature matches beta support boundaries", {

  priors <- list(
    full  = BayesTools::prior("beta", list(0.01, 0.01)),
    lower = BayesTools::prior("beta", list(0.01, 2.5),
                              truncation = list(0, 0.9)),
    upper = BayesTools::prior("beta", list(2.5, 0.01),
                              truncation = list(0.1, 1))
  )

  for (prior_pi in priors) {
    grid <- .glmm_binom_logit_pi_grid(
      ai       = 0L,
      ci       = 0L,
      n1i      = 1L,
      n2i      = 1L,
      prior_pi = prior_pi,
      n_pi     = 75L
    )
    nodes <- stats::plogis(grid[["grid"]][, 1L])

    expect_true(all(nodes > 0 & nodes < 1))
    expect_equal(sum(exp(grid[["log_weights"]][, 1L])), 1,
                 tolerance = 1e-10)
  }

  interior <- BayesTools::prior(
    "beta", list(0.01, 0.01), truncation = list(0.1, 0.9)
  )
  probabilities <- .gauss_legendre_nodes(17L)[["nodes"]]
  expected      <- stats::qlogis(as.numeric(
    BayesTools::quant(interior, probabilities)
  ))
  grid <- .glmm_binom_logit_pi_grid(
    ai       = 0L,
    ci       = 0L,
    n1i      = 1L,
    n2i      = 1L,
    prior_pi = interior,
    n_pi     = 17L
  )

  expect_identical(grid[["grid"]][, 1L], expected)
})


test_that("Poisson GLMM lograte priors are assigned", {

  old_lograte_sd <- RoBMA.get_option("default_lograte.sd")
  on.exit(RoBMA.options(default_lograte.sd = old_lograte_sd), add = TRUE)

  result_default <- brma.glmm(
    x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
    data = test_data_pois, measure = "IRR",
    prior_lograte = NULL, only_priors = TRUE
  )[["priors"]]

  expect_equal(result_default$outcome$phi$distribution, "normal")
  expected_lograte_mean <- log(sum(test_data_pois$x1i + test_data_pois$x2i) / sum(test_data_pois$t1i + test_data_pois$t2i))
  expect_equal(result_default$outcome$phi$parameters$mean, expected_lograte_mean)
  expect_equal(result_default$outcome$phi$parameters$sd, 1)

  exposure_multiplier <- 12
  test_data_pois_scaled <- transform(test_data_pois, t1i = t1i * exposure_multiplier, t2i = t2i * exposure_multiplier)
  result_scaled <- brma.glmm(
    x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
    data = test_data_pois_scaled, measure = "IRR",
    prior_lograte = NULL, only_priors = TRUE
  )[["priors"]]

  expect_equal(result_scaled$outcome$phi$parameters$mean, result_default$outcome$phi$parameters$mean - log(exposure_multiplier))
  expect_equal(result_scaled$outcome$phi$parameters$sd, result_default$outcome$phi$parameters$sd)

  RoBMA.options(default_lograte.sd = 0.75)
  result_option <- brma.glmm(
    x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
    data = test_data_pois, measure = "IRR",
    prior_lograte = NULL, only_priors = TRUE
  )[["priors"]]
  expect_equal(result_option$outcome$phi$parameters$sd, 0.75)

  custom_prior <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 2))
  result_custom <- brma.glmm(
    x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
    data = test_data_pois, measure = "IRR",
    prior_lograte = custom_prior, only_priors = TRUE
  )[["priors"]]

  expect_equal(result_custom$outcome$phi$parameters$sd, 2)
})


test_that("Poisson GLMM default lograte prior requires observed events", {

  expect_error(
    brma.glmm(
      x1i           = c(0L, 0L),
      x2i           = c(0L, 0L),
      t1i           = c(10, 12),
      t2i           = c(11, 13),
      measure       = "IRR",
      prior_lograte = NULL,
      only_priors   = TRUE
    ),
    regexp = "prior_lograte.*observed Poisson event|all-zero"
  )
})


test_that("Poisson GLMM validates custom lograte prior before transformation", {

  expect_error(
    brma.glmm(
      x1i           = x1i,
      x2i           = x2i,
      t1i           = t1i,
      t2i           = t2i,
      data           = test_data_pois,
      measure        = "IRR",
      prior_lograte  = list(distribution = "normal"),
      only_priors    = TRUE
    ),
    regexp = "prior_lograte.*prior distribution"
  )
})


test_that("GLMM nuisance priors must match the outcome measure", {

  baserate_prior <- BayesTools::prior("beta", parameters = list(alpha = 2, beta = 3))
  lograte_prior  <- BayesTools::prior("normal", parameters = list(mean = 0, sd = 2))

  expect_error_cases(list(
    list(
      label  = "single-model OR rejects lograte prior",
      expr   = quote(brma.glmm(
        ai            = ai,
        bi            = bi,
        ci            = ci,
        di            = di,
        data          = test_data_bin,
        measure       = "OR",
        prior_lograte = lograte_prior,
        only_priors   = TRUE
      )),
      regexp = "prior_lograte.*measure = 'IRR'|prior_lograte.*prior_baserate"
    ),
    list(
      label  = "single-model IRR rejects baserate prior",
      expr   = quote(brma.glmm(
        x1i            = x1i,
        x2i            = x2i,
        t1i            = t1i,
        t2i            = t2i,
        data           = test_data_pois,
        measure        = "IRR",
        prior_baserate = baserate_prior,
        only_priors    = TRUE
      )),
      regexp = "prior_baserate.*measure = 'OR'|prior_baserate.*prior_lograte"
    ),
    list(
      label  = "BMA.glmm OR rejects lograte prior",
      expr   = quote(BMA.glmm(
        ai            = ai,
        bi            = bi,
        ci            = ci,
        di            = di,
        data          = test_data_bin,
        measure       = "OR",
        prior_lograte = lograte_prior,
        only_priors   = TRUE
      )),
      regexp = "prior_lograte.*measure = 'IRR'|prior_lograte.*prior_baserate"
    ),
    list(
      label  = "BMA.glmm IRR rejects baserate prior",
      expr   = quote(BMA.glmm(
        x1i            = x1i,
        x2i            = x2i,
        t1i            = t1i,
        t2i            = t2i,
        data           = test_data_pois,
        measure        = "IRR",
        prior_baserate = baserate_prior,
        only_priors    = TRUE
      )),
      regexp = "prior_baserate.*measure = 'OR'|prior_baserate.*prior_lograte"
    )
  ))
})


test_that("GLMM rejects malformed publication-bias priors centrally", {

  object <- brma.glmm(
    ai = ai, bi = bi, ci = ci, di = di,
    data = test_data_bin, measure = "OR",
    only_priors = TRUE
  )
  bad_priors <- object[["priors"]]
  bad_priors[["outcome"]][["bias"]] <- BayesTools::prior_PET(
    "normal",
    parameters = list(mean = 0, sd = 1)
  )

  expect_error(
    .check_glmm_no_bias_priors(object[["data"]], bad_priors),
    "Publication-bias priors are not supported for GLMM outcomes"
  )
  expect_error(
    .create_fit_priors(object[["data"]], bad_priors),
    "Publication-bias priors are not supported for GLMM outcomes"
  )
  expect_error(
    .create_model_syntax(object[["data"]], bad_priors),
    "Publication-bias priors are not supported for GLMM outcomes"
  )
  expect_error(
    .log_lik_posterior_setup(
      fit                  = NULL,
      posterior_samples    = matrix(numeric(), nrow = 2L, ncol = 0L),
      data                 = object[["data"]],
      priors               = bad_priors,
      unit                 = "estimate",
      data_hash            = NULL
    ),
    "Publication-bias priors are not supported for GLMM outcomes"
  )
})


test_that("GLMM estimate-level factor priors declare complete factor levels", {

  for (outcome in c("bin", "pois")) {
    object <- if (outcome == "bin") {
      brma.glmm(
        ai = ai, bi = bi, ci = ci, di = di,
        data = test_data_bin, measure = "OR",
        only_priors = TRUE
      )
    } else {
      brma.glmm(
        x1i = x1i, t1i = t1i, x2i = x2i, t2i = t2i,
        data = test_data_pois, measure = "IRR",
        only_priors = TRUE
      )
    }
    fit_priors <- .create_fit_priors(object[["data"]], object[["priors"]])
    K <- nrow(object[["data"]][["outcome"]])
    for (name in c(if (outcome == "bin") "pi" else "phi", "theta")) {
      expect_identical(
        attr(fit_priors[[name]], "level_names"),
        as.character(seq_len(K)),
        info = paste(outcome, name)
      )
      expect_identical(
        dim(attr(fit_priors[[name]], "factor_design")),
        c(K, K),
        info = paste(outcome, name)
      )
    }
  }
})


test_that("BMA.glmm declares mean-difference factor defaults", {

  dat <- transform(test_data_bin, group = factor(c("a", "b", "a")))
  object <- BMA.glmm(ai = ai, bi = bi, ci = ci, di = di, mods = ~ group,
                     data = dat, measure = "OR", only_priors = TRUE)
  expect_identical(attr(object[["data"]], "set_contrast_factor_predictors"), "meandif")
  expect_s3_class(object[["priors"]][["mods"]][["group"]][[1L]], "prior.meandif")
})

test_that("a uniform(0, 1) base-rate prior is Beta(1, 1) and takes the beta route", {

  skip_on_cran()
  stored_prior <- function(prior) brma.glmm(
    ai = ai, bi = bi, ci = ci, di = di, data = test_data_bin, measure = "OR",
    prior_baserate = prior, only_priors = TRUE
  )[["priors"]][["outcome"]][["pi"]]
  uniform   <- stored_prior(BayesTools::prior("uniform", list(0, 1)))
  explicit  <- stored_prior(BayesTools::prior("uniform", list(0, 1), truncation = list(0, 1)))
  beta      <- stored_prior(BayesTools::prior("beta", list(1, 1)))
  half      <- stored_prior(BayesTools::prior("uniform", list(0, 0.5)))

  expect_identical(.glmm_aghq_prior_spec(uniform, "bin"), c(alpha = 1, beta = 1))
  expect_identical(.glmm_aghq_prior_spec(explicit, "bin"), c(alpha = 1, beta = 1))
  expect_identical(.glmm_aghq_prior_spec(beta, "bin"), c(alpha = 1, beta = 1))
  # a narrower uniform prior is not Beta(1, 1) and keeps the prior-CDF grid
  expect_null(.glmm_aghq_prior_spec(half, "bin"))

  # Large arms: the prior-CDF grid cannot certify the uniform density at this
  # cell, whereas the beta route gives the explicit Beta(1, 1) value exactly.
  value <- function(prior) .outcome_pdf.binom(
    80L, 220L, 2460L, 2780L, matrix(0.2), matrix(0.4), prior)
  expect_identical(value(uniform), value(beta))
  expect_identical(value(explicit), value(beta))

  # The same fitted draws give the same estimate-unit log-likelihood.
  fit <- suppressWarnings(brma.glmm(
    ai = c(4L, 6L, 3L, 5L), bi = c(119L, 300L, 228L, 99L),
    ci = c(11L, 29L, 11L, 12L), di = c(128L, 274L, 209L, 90L),
    measure = "OR", prior_baserate = BayesTools::prior("uniform", list(0, 1)),
    prior_unit_information_sd = 1, chains = 2, adapt = 50, burnin = 100,
    sample = 100, seed = 1, silent = TRUE, parallel = FALSE
  ))
  with_beta <- fit
  with_beta[["priors"]][["outcome"]][["pi"]] <- beta
  expect_identical(log_lik(fit), log_lik(with_beta))
})

test_that("the automatic LOO/WAIC of a GLMM with an uncertified likelihood warns and returns the fit", {

  condition <- .glmm_grid_nonconvergence_condition("pois", 1L, c(2e-5, 3e-6), c(4e-7, 5e-8))
  object <- structure(list(fit = NULL), class = c("brma.glmm", "brma"))
  previous <- c(loo = RoBMA.get_option("autocompute.loo"),
                waic = RoBMA.get_option("autocompute.waic"))
  RoBMA.options(autocompute.loo = TRUE, autocompute.waic = TRUE)
  on.exit(RoBMA.options(autocompute.loo = previous[["loo"]],
                        autocompute.waic = previous[["waic"]]), add = TRUE)

  calls <- character()
  testthat::local_mocked_bindings(
    add_loo  = function(object, ...) { calls <<- c(calls, "loo"); stop(condition) },
    add_waic = function(object, ...) { calls <<- c(calls, "waic"); stop(condition) },
    .package = "RoBMA"
  )
  messages <- character()
  returned <- withCallingHandlers(
    .autocompute_brma(object),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(returned, object)
  # one warning with the condition message; WAIC shares the failed log-likelihood
  expect_identical(messages, paste(conditionMessage(condition),
                                   "The fit is returned without LOO/WAIC."))
  expect_identical(calls, "loo")

  # other errors are not swallowed
  testthat::local_mocked_bindings(
    add_loo = function(object, ...) stop("Another error."),
    .package = "RoBMA"
  )
  expect_error(.autocompute_brma(object), "Another error.", fixed = TRUE)
})
