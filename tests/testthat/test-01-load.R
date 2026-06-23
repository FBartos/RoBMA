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
    "sel_z_lower, sel_z_upper, sel_obs_bin[i], sel_sign, ",
    "sel_telescope_probabilities)\n",
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
    "sel_kernel_mode_active, sel_telescope_probabilities)\n",
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

test_that("RoBMA JAGS module exposes known-V multivariate normal kernel", {

  RoBMA:::.load_RoBMA_module()

  y       <- c(0.10, -0.20)
  mu      <- c(0.03, -0.01)
  tau2    <- c(0.01, 0.04)
  V       <- matrix(c(0.04, 0.015, 0.015, 0.09), nrow = 2)
  V_lower <- V[lower.tri(V, diag = TRUE)]

  con <- textConnection(paste0(
    "model{\n",
    "  y[1:2] ~ dknown_v_mnorm(mu[1:2], tau2[1:2], V_lower[1:3])\n",
    "}\n"
  ))
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file     = con,
    data     = list(y = y, mu = mu, tau2 = tau2, V_lower = V_lower),
    quiet    = TRUE,
    n.chains = 2,
    n.adapt  = 0
  )
  dic <- rjags::dic.samples(
    model,
    n.iter       = 5,
    type         = "pD",
    progress.bar = "none"
  )
  Sigma              <- V + diag(tau2)
  Sigma_chol         <- chol(Sigma)
  residual           <- y - mu
  precision_residual <- backsolve(
    Sigma_chol,
    forwardsolve(t(Sigma_chol), residual)
  )
  reference_deviance <- length(y) * log(2 * pi) +
    2 * sum(log(diag(Sigma_chol))) +
    sum(residual * precision_residual)

  expect_equal(
    as.numeric(sum(dic[["deviance"]])),
    reference_deviance,
    tolerance = 1e-8
  )

  se_rank_one      <- c(0.20, 0.30, 0.40)
  y_rank_one       <- c(0.10, 0.20, -0.05)
  mu_rank_one      <- c(0.02, 0.00, 0.03)
  tau2_rank_one    <- rep(0.01, 3)
  V_rank_one       <- tcrossprod(se_rank_one)
  V_rank_one_lower <- V_rank_one[lower.tri(V_rank_one, diag = TRUE)]

  con_rank_one <- textConnection(paste0(
    "model{\n",
    "  y[1:3] ~ dknown_v_mnorm(mu[1:3], tau2[1:3], V_lower[1:6])\n",
    "}\n"
  ))
  on.exit(close(con_rank_one), add = TRUE)

  model_rank_one <- rjags::jags.model(
    file     = con_rank_one,
    data     = list(
      y       = y_rank_one,
      mu      = mu_rank_one,
      tau2    = tau2_rank_one,
      V_lower = V_rank_one_lower
    ),
    quiet    = TRUE,
    n.chains = 2,
    n.adapt  = 0
  )
  dic_rank_one <- rjags::dic.samples(
    model_rank_one,
    n.iter       = 5,
    type         = "pD",
    progress.bar = "none"
  )
  Sigma_rank_one              <- V_rank_one + diag(tau2_rank_one)
  Sigma_rank_one_chol         <- chol(Sigma_rank_one)
  residual_rank_one           <- y_rank_one - mu_rank_one
  precision_residual_rank_one <- backsolve(
    Sigma_rank_one_chol,
    forwardsolve(t(Sigma_rank_one_chol), residual_rank_one)
  )
  reference_rank_one_deviance <- length(y_rank_one) * log(2 * pi) +
    2 * sum(log(diag(Sigma_rank_one_chol))) +
    sum(residual_rank_one * precision_residual_rank_one)

  expect_equal(
    as.numeric(sum(dic_rank_one[["deviance"]])),
    reference_rank_one_deviance,
    tolerance = 1e-8
  )
})

test_that("RoBMA JAGS module can sample known-V multivariate normal nodes", {

  RoBMA:::.load_RoBMA_module()

  mu      <- c(0.03, -0.01)
  tau2    <- c(0.01, 0.04)
  V       <- matrix(c(0.04, 0.015, 0.015, 0.09), nrow = 2)
  V_lower <- V[lower.tri(V, diag = TRUE)]

  con <- textConnection(paste0(
    "model{\n",
    "  y[1:2] ~ dknown_v_mnorm(mu[1:2], tau2[1:2], V_lower[1:3])\n",
    "}\n"
  ))
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file     = con,
    data     = list(mu = mu, tau2 = tau2, V_lower = V_lower),
    quiet    = TRUE,
    n.chains = 2,
    n.adapt  = 0
  )
  fit <- rjags::coda.samples(
    model,
    variable.names = "y",
    n.iter         = 5,
    quiet          = TRUE,
    progress.bar   = "none"
  )

  samples <- as.matrix(fit)
  expect_true(all(is.finite(samples[, c("y[1]", "y[2]")])))
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
