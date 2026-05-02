context("JAGS module functionality")
skip_on_cran()

test_that("RoBMA JAGS module exposes scalar selected-normal step kernel", {

  RoBMA:::.load_RoBMA_module()

  yi  <- c(.05, .18, .32)
  sei <- c(.10, .12, .15)
  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive",
    signed_data      = TRUE
  )

  model_syntax <- paste0(
    "model{\n",
    "  mu ~ dnorm(0, 1)\n",
    "  for(i in 1:K){\n",
    "    yi[i] ~ dselnorm_step(mu, sigma, sei[i], 1, omega,",
    "sel_z_lower, sel_z_upper, sel_obs_bin[i], sel_sign)\n",
    "  }\n",
    "}\n"
  )
  data <- c(
    list(
      yi    = yi,
      sei   = sei,
      K     = length(yi),
      sigma = .25,
      omega = c(1, .5, .25)
    ),
    spec[["jags_data"]]
  )

  model <- rjags::jags.model(
    file       = textConnection(model_syntax),
    data       = data,
    quiet      = TRUE,
    n.adapt    = 1000
  )
  fit <- rjags::coda.samples(
    model          = model,
    variable.names = "mu",
    n.iter         = 20,
    quiet          = TRUE,
    progress.bar   = "none"
  )

  expect_equal(nrow(fit[[1]]), 20)
  expect_true(all(is.finite(fit[[1]][, "mu"])))
})

test_that("RoBMA JAGS module exposes scalar selected-normal step switch kernel", {

  RoBMA:::.load_RoBMA_module()

  yi  <- c(.05, .18, .32)
  sei <- c(.10, .12, .15)
  prior <- BayesTools::prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = BayesTools::wf_fixed(c(1, .5, .25))
  )
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive",
    signed_data      = TRUE
  )

  model_syntax <- paste0(
    "model{\n",
    "  bias_indicator ~ dcat(bias_prob)\n",
    "  sel_kernel_mode_active = 0 * equals(bias_indicator, 1) + ",
    "1 * equals(bias_indicator, 2)\n",
    "  mu ~ dnorm(0, 1)\n",
    "  for(i in 1:K){\n",
    "    yi[i] ~ dselnorm_step_switch(mu, sigma, sei[i], 1, omega,",
    "sel_z_lower, sel_z_upper, sel_obs_bin[i], sel_sign, ",
    "sel_kernel_mode_active)\n",
    "  }\n",
    "}\n"
  )
  data <- c(
    list(
      yi        = yi,
      sei       = sei,
      K         = length(yi),
      sigma     = .25,
      omega     = c(1, .5, .25),
      bias_prob = c(1, 1)
    ),
    spec[["jags_data"]]
  )

  model <- rjags::jags.model(
    file       = textConnection(model_syntax),
    data       = data,
    quiet      = TRUE,
    n.adapt    = 1000
  )
  fit <- rjags::coda.samples(
    model          = model,
    variable.names = c("mu", "bias_indicator"),
    n.iter         = 20,
    quiet          = TRUE,
    progress.bar   = "none"
  )

  expect_equal(nrow(fit[[1]]), 20)
  expect_true(all(is.finite(fit[[1]][, "mu"])))
  expect_true(all(fit[[1]][, "bias_indicator"] %in% c(1, 2)))
})

test_that("dwbinom has discrete bounded support and can sample", {

  RoBMA:::.load_RoBMA_module()

  con <- textConnection("model{\n  y ~ dwbinom(p, N, weight)\n}\n")
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file    = con,
    data    = list(p = .30, N = 5L, weight = 1),
    inits   = list(y = 2L),
    quiet   = TRUE,
    n.adapt = 0
  )
  fit <- rjags::coda.samples(
    model,
    variable.names = "y",
    n.iter         = 50,
    quiet          = TRUE,
    progress.bar   = "none"
  )
  draws <- as.numeric(as.matrix(fit)[, "y"])

  expect_true(all(draws >= 0 & draws <= 5))
  expect_true(all(draws == floor(draws)))

  con_bad <- textConnection("model{\n  y ~ dwbinom(p, N, weight)\n}\n")
  on.exit(close(con_bad), add = TRUE)
  expect_error(
    rjags::jags.model(
      file    = con_bad,
      data    = list(y = 3L, p = .30, N = 5.5, weight = 1),
      quiet   = TRUE,
      n.adapt = 0
    ),
    "discrete-valued parameters|Invalid parent values"
  )
})

test_that("dwpois has discrete nonnegative support and can sample", {

  RoBMA:::.load_RoBMA_module()

  con <- textConnection("model{\n  y ~ dwpois(lambda, weight)\n}\n")
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file    = con,
    data    = list(lambda = 2.5, weight = 1),
    inits   = list(y = 2L),
    quiet   = TRUE,
    n.adapt = 0
  )
  fit <- rjags::coda.samples(
    model,
    variable.names = "y",
    n.iter         = 50,
    quiet          = TRUE,
    progress.bar   = "none"
  )
  draws <- as.numeric(as.matrix(fit)[, "y"])

  expect_true(all(draws >= 0))
  expect_true(all(draws == floor(draws)))

  con_bad <- textConnection("model{\n  y ~ dwpois(lambda, weight)\n}\n")
  on.exit(close(con_bad), add = TRUE)
  expect_error(
    rjags::jags.model(
      file    = con_bad,
      data    = list(y = 3L, lambda = -1, weight = 1),
      quiet   = TRUE,
      n.adapt = 0
    ),
    "Invalid parent values"
  )
})
