probe_state("test-00-joint-bias-native.R")
skip_on_cran()

.joint_bias_native_fit <- function(extra = "", slopes = "c(pet, 0, 0)", power = 1,
                                   chains = 1L, zero_offset = FALSE, declare = TRUE) {

  setup <- .selection_runtime_setup()
  .JAGS_run_runtime_setup <- getFromNamespace(".JAGS_run_runtime_setup", "BayesTools")
  .JAGS_run_runtime_setup(setup, 1L)
  tag <- if (zero_offset) "robma_joint_bias_moderated" else "robma_joint_bias"
  extra_args <- if (zero_offset) ", y, 0, c(0, 0, 0)" else ""
  declaration <- if (declare) paste0("tag <- ", tag, "(bias, effect, raw, 1, -1, c(2, 1, 0), ",
    slopes, ", c(1, 1, 1, weights, 1, 1, 1)", extra_args, ")") else ""
  syntax <- paste0("model {
    bias ~ dcat(c(1, 2, 1))
    effect ~ dcat(c(3, 1))
    raw ~ dnorm(0.13, 2.5)
    pet <- 1.2
    mu <- -raw * equals(effect, 1)
    for (j in 1:3) {
      active[j] <- weights[j] * equals(bias, 2) + 1 - equals(bias, 2)
    }
    for (i in 1:3) {
      location[i] <- mu + pet * tilt[i] * equals(bias, 1)
      y[i] ~ dselnorm_step_switch(location[i], total_sd[i], se[i], power,
                                  active, z_lower, z_upper, bin[i], 1, equals(bias, 2), 0)
    }
    ", declaration, "
    ", extra, "
  }")
  data <- list(y = c(.4, .1, .7), total_sd = c(.2, .4, .65), se = c(.2, .4, .6),
               tilt = c(.2, 0, .6), power = power, weights = c(1, .8, .35),
               z_lower = c(1.7, .6, -1e300), z_upper = c(1e300, 1.7, .6),
               bin = c(1, 0, 2))
  if (nzchar(extra)) data$extra <- 0
  inits <- lapply(seq_len(chains), function(chain) {
    list(bias = 3L, effect = 2L, raw = .1, .RNG.name = "base::Mersenne-Twister",
         .RNG.seed = 1031L + chain)
  })
  con <- textConnection(syntax)
  on.exit(close(con))
  list(model = rjags::jags.model(con, data = data, inits = inits,
                                n.chains = chains, n.adapt = 0L, quiet = TRUE),
       data = data)
}

test_that("the native joint-bias factory exclusively owns the declared coordinates", {

  skip_if_not_installed("rjags")
  fit <- .joint_bias_native_fit()
  samplers <- rjags::list.samplers(fit$model)
  expect_named(samplers, "RoBMA::IndependentBiasJoint")
  expect_setequal(samplers[[1L]], c("bias", "effect", "raw"))
  draws <- as.matrix(rjags::coda.samples(fit$model, c("effect", "raw", "mu"),
                                        n.iter = 20L, progress.bar = "none"))
  expect_equal(draws[, "mu"], -draws[, "raw"] * (draws[, "effect"] == 1L))

  declined <- list(
    .joint_bias_native_fit(extra = "extra ~ dnorm(raw, 1)"),
    .joint_bias_native_fit(extra = "extra ~ dbt_gate(raw, 1)"),
    .joint_bias_native_fit(extra = "extra ~ dnorm(tag, 1)"),
    .joint_bias_native_fit(slopes = "c(raw, 0, 0)"),
    .joint_bias_native_fit(power = 2)
  )
  for (fit in declined) {
    samplers <- rjags::list.samplers(fit$model)
    expect_false("RoBMA::IndependentBiasJoint" %in% names(samplers))
    expect_setequal(unlist(samplers, use.names = FALSE), c("bias", "effect", "raw"))
  }
})

.moderated_joint_bias_native_fit <- function(coefficients = "c(beta1, beta2)",
                                             outcomes = "y", design = "design", declare = TRUE) {

  setup <- .selection_runtime_setup()
  .JAGS_run_runtime_setup <- getFromNamespace(".JAGS_run_runtime_setup", "BayesTools")
  .JAGS_run_runtime_setup(setup, 1L)
  declaration <- if (declare) paste0("tag <- robma_joint_bias_moderated(bias, effect, raw, 1, -1, c(2, 1, 0),
      c(pet, 0, 0), c(1, 1, 1, weights, 1, 1, 1), ", outcomes, ", ", coefficients, ", ", design, ")") else ""
  syntax <- paste0("model {
    bias ~ dcat(c(1, 2, 1))
    effect ~ dcat(c(3, 1))
    raw ~ dnorm(.13, 2.5)
    beta1 ~ dnorm(0, 1)
    beta2 ~ dnorm(0, 1)
    covariates[1] <- beta1
    covariates[2] <- beta2
    covariates[3] <- 0
    covariates[4] <- 0
    covariates[5] <- beta2
    covariates[6] <- beta1
    mu <- -raw * equals(effect, 1)
    for (j in 1:3) {
      active[j] <- weights[j] * equals(bias, 2) + 1 - equals(bias, 2)
    }
    for (i in 1:3) {
      location[i] <- mu - beta1 * x1[i] - beta2 * x2[i] + pet * tilt[i] * equals(bias, 1)
      y[i] ~ dselnorm_step_switch(location[i], total_sd[i], se[i], 1,
                                  active, z_lower, z_upper, bin[i], 1, equals(bias, 2), 0)
    }
    ", declaration, "
  }")
  x <- cbind(c(-1, .25, .75), c(.5, -.5, .5))
  data <- list(y = c(.4, .4, .7), total_sd = c(.3, .3, .65), se = c(.2, .2, .6),
               tilt = c(.2, 0, .6), pet = 1.2, weights = c(1, .8, .35),
               z_lower = c(1.7, .6, -1e300), z_upper = c(1e300, 1.7, .6),
               bin = c(1, 0, 2), x1 = x[, 1L], x2 = x[, 2L])
  if (declare && design == "design") data$design <- c(t(x))
  con <- textConnection(syntax)
  on.exit(close(con))
  inits <- list(bias = 3L, effect = 2L, raw = .1, beta1 = .2, beta2 = -.25,
                .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2554L)
  rjags::jags.model(con, data = data, inits = inits, n.chains = 1L, n.adapt = 0L, quiet = TRUE)
}

test_that("moderated joint bias owns the original block and reads only unowned coefficients", {

  skip_if_not_installed("rjags")
  fit <- .moderated_joint_bias_native_fit()
  owners <- rjags::list.samplers(fit)
  expect_setequal(owners[["RoBMA::IndependentBiasJoint"]], c("bias", "effect", "raw"))
  expect_setequal(unlist(owners, use.names = FALSE), c("bias", "effect", "raw", "beta1", "beta2"))
  draws <- as.matrix(rjags::coda.samples(fit, c("effect", "raw", "mu"),
                                        n.iter = 20L, progress.bar = "none"))
  expect_equal(draws[, "mu"], -draws[, "raw"] * (draws[, "effect"] == 1L))

  for (declined in list(
    .moderated_joint_bias_native_fit(coefficients = "c(beta1, raw)"),
    .moderated_joint_bias_native_fit(outcomes = "y[1:2]"),
    .moderated_joint_bias_native_fit(outcomes = "c(y[1], y[1], y[3])"),
    .moderated_joint_bias_native_fit(design = "covariates"))) {
    owners <- rjags::list.samplers(declined)
    expect_false("RoBMA::IndependentBiasJoint" %in% names(owners))
    expect_setequal(unlist(owners, use.names = FALSE), c("bias", "effect", "raw", "beta1", "beta2"))
  }
})

test_that("moderated joint bias refreshes private coefficient inputs on each chain", {

  skip_if_not_installed("rjags")
  getFromNamespace(".JAGS_run_runtime_setup", "BayesTools")(.selection_runtime_setup(), 1L)
  run <- function(mode) {
    coefficient <- switch(mode, direct = "beta", private = "negative_beta", zero = "0")
    # negative_beta is used only by the terminal tag. Monitoring it could
    # refresh the node and mask a stale private input to the proposal.
    syntax <- paste0("model {
      bias ~ dcat(c(1, 2))
      effect ~ dcat(c(3, 1))
      raw ~ dnorm(.13, 2.5)
      beta ~ dnorm(0, 1)
      negative_beta <- -beta
      mu <- -raw * equals(effect, 1)
      active[1] <- 1
      active[2] <- 1 - .6 * equals(bias, 2)
      for (i in 1:3) {
        location[i] <- mu - beta * x[i]
        y[i] ~ dselnorm_step_switch(location[i], se[i], se[i], 1,
                                    active, zlo, zhi, bin[i], 1, equals(bias, 2), 0)
      }
      tag <- robma_joint_bias_moderated(bias, effect, raw, 1, -1, c(0, 1),
        c(0, 0), c(1, 1, 1, .4), y, ", coefficient, ", design)
    }")
    x <- c(-1, .25, .75)
    data <- list(x = x, design = if (mode == "private") -x else x, y = c(1.4, -.1, -.4),
                 se = c(.21, .28, .35), zlo = c(1.3, -1e300), zhi = c(1e300, 1.3), bin = c(1L, 0L, 2L))
    inits <- lapply(seq_len(2L), function(chain) {
      list(beta = c(-1, 2)[chain], raw = .1, bias = 1L, effect = 2L,
           .RNG.name = "base::Mersenne-Twister", .RNG.seed = 701L + chain)
    })
    con <- textConnection(syntax)
    on.exit(close(con))
    fit <- rjags::jags.model(con, data, inits, n.chains = 2L, n.threads = 1L,
                            n.adapt = 0L, quiet = TRUE)
    owners <- rjags::list.samplers(fit)
    expect_setequal(owners[["RoBMA::IndependentBiasJoint"]], c("bias", "effect", "raw"))
    draws <- rjags::coda.samples(fit, c("beta", "bias", "effect", "raw", "mu", "location[1]"),
                                n.iter = 64L, progress.bar = "none")
    list(owners = owners, draws = draws, state = fit$state(internal = TRUE))
  }

  direct <- run("direct")
  expect_identical(run("private"), direct)
  zero <- run("zero")
  raw_draws <- function(result) lapply(result$draws, function(chain) chain[, "raw"])
  expect_false(identical(raw_draws(zero), raw_draws(direct)))
  for (chain in direct$draws) {
    beta <- chain[, "beta"]
    expect_gt(length(unique(beta[beta > 0])), 1L) # Distinct negative values of private negative_beta.
    expect_equal(chain[, "mu"], -chain[, "raw"] * (chain[, "effect"] == 1L))
    expect_equal(chain[, "location[1]"], chain[, "mu"] + beta)
  }
})

test_that("a zero moderator offset replays the unmoderated joint update", {

  skip_if_not_installed("rjags")
  fits <- list(.joint_bias_native_fit()$model, .joint_bias_native_fit(zero_offset = TRUE)$model)
  draws <- lapply(fits, function(fit) {
    expect_named(rjags::list.samplers(fit), "RoBMA::IndependentBiasJoint")
    as.matrix(rjags::coda.samples(fit, c("bias", "effect", "raw"), n.iter = 30L, progress.bar = "none"))
  })
  expect_identical(draws[[1L]], draws[[2L]])
})

test_that("native-declined bias tags preserve ordinary owners, draws and RNG state", {

  skip_if_not_installed("rjags")
  variants <- list(
    function(tag) .joint_bias_native_fit(power = 2, declare = tag)$model,
    function(tag) .moderated_joint_bias_native_fit(coefficients = "c(beta1, raw)", declare = tag)
  )
  run <- function(fit) {
    owners <- rjags::list.samplers(fit)
    expect_false("RoBMA::IndependentBiasJoint" %in% names(owners))
    draws <- rjags::coda.samples(fit, c("bias", "effect", "raw", "mu"), n.iter = 30L,
                                progress.bar = "none")
    list(owners = owners, draws = draws, state = fit$state(internal = TRUE))
  }
  for (variant in variants) expect_identical(run(variant(TRUE)), run(variant(FALSE)))
})

test_that("the public independent joint bias update works with or without the nested weight layout", {

  previous <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous), add = TRUE)
  RoBMA.options(sampler.nested.weights = TRUE,
                sampler.probabilities = character(), selection.sampler = "slice")
  data("Bem2011", package = "RoBMA", envir = environment())
  fit <- function(fine_alpha, policy = "auto") {
    RoBMA.options(sampler.joint = policy)
    RoBMA(yi = d, sei = se, data = Bem2011, measure = "SMD", prior_bias = list(
      BayesTools::prior_weightfunction("one-sided", c(.05, .5), BayesTools::wf_cumulative(rep(1, 3))),
      BayesTools::prior_weightfunction("one-sided", c(.025, .05, .5), BayesTools::wf_cumulative(fine_alpha))
    ), chains = 1L, adapt = 500L, burnin = 50L, sample = 100L, seed = 704L,
    convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))[["fit"]]
  }
  declarations <- function(fit) attr(fit, "sampler_declarations", exact = TRUE)
  owned <- function(fit) {
    assignment <- fit[["samplers"]]
    assignment[["Node"]][assignment[["Sampler"]] == "RoBMA::IndependentBiasJoint"]
  }
  declined <- fit(c(2, 1, 1, 1))
  nested <- declarations(declined)[["nested_weights"]]
  expect_false(nested[["realized"]])
  expect_match(nested[["reason"]], "no eligible nested")

  combined <- fit(rep(1, 4))
  expect_true(declarations(combined)[["nested_weights"]][["realized"]])
  for (backend in list(declined, combined)) {
    expect_setequal(owned(backend), c("bias_indicator", "mu_indicator", "mu_component_2"))
    joint <- declarations(backend)[["robma"]][["joint"]]
    expect_identical(joint[["policy"]], "auto")
    expect_identical(joint[["sampler"]], "RoBMA::IndependentBiasJoint")
    expect_identical(joint[["status"]], "active")
  }

  nested_only <- fit(rep(1, 4), policy = "off")
  expect_true(declarations(nested_only)[["nested_weights"]][["realized"]])
  expect_length(owned(nested_only), 0L)
})

test_that("a single bias branch preserves the ordinary public fit under auto", {

  previous <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous), add = TRUE)
  RoBMA.options(sampler.nested.weights = FALSE, sampler.probabilities = character(),
                selection.sampler = "slice.coarse")
  data("Bem2011", package = "RoBMA", envir = environment())
  fit <- function(policy) {
    RoBMA.options(sampler.joint = policy)
    RoBMA(yi = d, sei = se, data = Bem2011, measure = "SMD",
      prior_bias = BayesTools::prior_weightfunction("one-sided", c(.05, .5)),
      prior_bias_null = FALSE, chains = 1L, adapt = 500L, burnin = 50L, sample = 100L,
      seed = 704L, convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))[["fit"]]
  }
  auto <- fit("auto")
  off <- fit("off")
  expect_false("RoBMA::IndependentBiasJoint" %in% auto[["samplers"]][["Sampler"]])
  expect_identical(auto[["samplers"]], off[["samplers"]])
  expect_identical(auto[["mcmc"]], off[["mcmc"]])
})
