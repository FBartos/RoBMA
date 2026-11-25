probe_state("test-00-joint-bias-nested-native.R")
skip_on_cran()

.joint_bias_nested_graph <- function(direction = -1, beta_shape = "1 + equals(bias, 3)",
                                     truncation = "", observed = FALSE, extra = "", data = list()) {

  setup <- .selection_runtime_setup()
  .JAGS_run_runtime_setup <- getFromNamespace(".JAGS_run_runtime_setup", "BayesTools")
  .JAGS_run_runtime_setup(setup, 1L)
  syntax <- paste0("model {
    bias ~ dcat(c(1, 1.7, 0.56))
    effect ~ dcat(c(0.7, 1.4))
    raw ~ dnorm(0.13, 2.5)
    tau ~ dnorm(0, 4)T(0,)
    a ~ dbeta(1 + equals(bias, 2) + equals(bias, 3), ", beta_shape, ")", truncation, "
    v ~ dbeta(1, 1)
    u ~ dbeta(1, 1)
    coarse[1:4] <- c(1, 1, a, a * v)
    fine[1:4] <- c(1, a + (1 - a) * u, a, a * v)
    for (j in 1:4) {
      active[j] <- equals(bias, 1) + coarse[j] * equals(bias, 2) + fine[j] * equals(bias, 3)
    }
    mu <- direction * raw * equals(effect, 2)
    for (i in 1:5) {
      sd[i] <- sqrt(pow(se[i], 2) + pow(tau, 2))
      y[i] ~ dselnorm_step_switch(mu, sd[i], se[i], 1, active, zlo, zhi,
                                  bin[i], 1, 1 - equals(bias, 1), 0)
    }
    tag <- robma_joint_bias(bias, effect, raw, 2, direction, c(0, 1, 1),
                               c(0, 0, 0), c(1, 1, 1, 1, coarse, fine))
    ", extra, "
  }")
  se <- c(.16, .25, .32, .50, .81)
  data <- c(list(se = se, y = se * c(2.3, 1.3, .7, .1, -1.3), bin = c(1L, 2L, 3L, 4L, 0L),
                 zlo = c(1.8, 1.15, .55, -1e300), zhi = c(1e300, 1.8, 1.15, .55),
                 direction = direction), data)
  initial <- list(bias = 1L, effect = 2L, raw = -.19, a = .09, v = .28, u = .33, tau = .07,
                  .RNG.name = "base::Mersenne-Twister", .RNG.seed = 317L)
  if (observed) {
    data[["a"]] <- initial[["a"]]
    initial[["a"]] <- NULL
  }
  con <- textConnection(syntax)
  on.exit(close(con))
  rjags::jags.model(con, data = data, inits = list(initial), n.chains = 1L, n.threads = 1L,
                   n.adapt = 0L, quiet = TRUE)
}

test_that("the joint bias sampler can leave the shared Beta weight as a nuisance root", {

  skip_if_not_installed("rjags")
  for (direction in c(-1, 1)) {
    fit <- .joint_bias_nested_graph(direction = direction)
    owners <- rjags::list.samplers(fit)
    expect_setequal(owners[["RoBMA::IndependentBiasJoint"]], c("bias", "effect", "raw"))
    expect_setequal(unlist(owners[names(owners) != "RoBMA::IndependentBiasJoint"], use.names = FALSE),
                    c("a", "v", "u", "tau"))
    draws <- as.matrix(rjags::coda.samples(fit, c("effect", "raw", "mu"), n.iter = 20L,
                                          progress.bar = "none"))
    expect_equal(draws[, "mu"], direction * draws[, "raw"] * (draws[, "effect"] == 2L))
  }
})

test_that("joint bias declines shared-weight graphs with additional dependencies", {

  skip_if_not_installed("rjags")
  controls <- list(
    effect_shape = list(beta_shape = "1 + equals(bias, 3) + equals(effect, 2)"),
    observed_shape = list(beta_shape = "1 + equals(bias, 3) + obs_shape",
                          extra = "obs_shape ~ dnorm(0, 1)", data = list(obs_shape = .5)),
    nuisance_shape = list(beta_shape = "1 + equals(bias, 3) + nuisance",
                          extra = "nuisance ~ dbeta(1, 1)"),
    observed_weight = list(observed = TRUE),
    truncated_weight = list(truncation = "T(.05, .97)"),
    unrelated_weight = list(extra = "other ~ dbeta(1 + equals(bias, 2), 1)\n witness ~ dnorm(other, 1)",
                            data = list(witness = 0)),
    extra_weight_child = list(extra = "child ~ dbeta(1 + a, 1)\n witness ~ dnorm(child, 1)",
                              data = list(witness = 0))
  )
  for (reason in names(controls)) {
    fit <- do.call(.joint_bias_nested_graph, controls[[reason]])
    owners <- rjags::list.samplers(fit)
    expect_false("RoBMA::IndependentBiasJoint" %in% names(owners), info = reason)
    expect_true(all(c("bias", "effect", "raw") %in% unlist(owners, use.names = FALSE)), info = reason)
  }
})
