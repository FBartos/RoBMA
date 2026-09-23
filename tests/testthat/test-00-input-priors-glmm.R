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


test_that("GLMM nuisance priors reject unsupported families before fitting", {

  discrete <- BayesTools::prior("bernoulli", parameters = list(probability = 0.5))

  expect_error(
    .assign_prior.baserate(discrete),
    "Non-point discrete priors are not supported for 'prior_baserate'"
  )
  expect_error(
    .assign_prior.lograte(discrete, test_data_pois),
    "Non-point discrete priors are not supported for 'prior_lograte'"
  )

  unsupported <- list(
    BayesTools::prior("gamma", list(shape = 2, rate = 1)),
    BayesTools::prior("lognormal", list(meanlog = 0, sdlog = 1)),
    BayesTools::prior("invgamma", list(shape = 2, scale = 1))
  )
  for (prior in unsupported) {
    expect_error(
      .assign_prior.baserate(prior),
      paste0("Only point and beta priors.*'", prior[["distribution"]], "'")
    )
    expect_error(
      .assign_prior.lograte(prior, test_data_pois),
      paste0("Only point and normal priors.*'", prior[["distribution"]], "'")
    )
  }

  expect_error(
    brma.glmm(
      ai = ai, bi = bi, ci = ci, di = di,
      data = test_data_bin, measure = "OR",
      prior_baserate = BayesTools::prior(
        "normal", list(0.5, 0.2), truncation = list(0, 1)
      ),
      only_priors = TRUE
    ),
    "Only point and beta priors.*post-fit likelihood diagnostics.*'normal'"
  )
  expect_error(
    BMA.glmm(
      x1i = x1i, x2i = x2i, t1i = t1i, t2i = t2i,
      data = test_data_pois, measure = "IRR",
      prior_lograte = BayesTools::prior("lognormal", list(0, 1)),
      only_priors = TRUE
    ),
    "Only point and normal priors.*post-fit likelihood diagnostics.*'lognormal'"
  )
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
