skip_on_cran()

.joint_bias_native_fit <- function(extra = "", slopes = "c(pet, 0, 0)", power = 1,
                                   chains = 1L) {

  setup <- .selection_runtime_setup()
  .JAGS_run_runtime_setup <- getFromNamespace(".JAGS_run_runtime_setup", "BayesTools")
  .JAGS_run_runtime_setup(setup, 1L)
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
    marker ~ drobma_joint_bias(bias, effect, raw, 1, -1, c(2, 1, 0), ", slopes, ",
                             c(1, 1, 1, weights, 1, 1, 1))
    ", extra, "
  }")
  data <- list(y = c(.4, .1, .7), total_sd = c(.2, .4, .65), se = c(.2, .4, .6),
               tilt = c(.2, 0, .6), power = power, weights = c(1, .8, .35),
               z_lower = c(1.7, .6, -1e300), z_upper = c(1e300, 1.7, .6),
               bin = c(1, 0, 2), marker = 0)
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
    .joint_bias_native_fit(slopes = "c(raw, 0, 0)"),
    .joint_bias_native_fit(power = 2)
  )
  for (fit in declined) {
    samplers <- rjags::list.samplers(fit$model)
    expect_false("RoBMA::IndependentBiasJoint" %in% names(samplers))
    expect_setequal(unlist(samplers, use.names = FALSE), c("bias", "effect", "raw"))
  }
})
