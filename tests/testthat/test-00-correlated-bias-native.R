skip_on_cran()

.correlated_bias_native_fit <- function(gate, nuisance = "within_sd", declaration = "", chains = 1L) {

  setup <- .selection_runtime_setup()
  getFromNamespace(".JAGS_run_runtime_setup", "BayesTools")(setup, 1L)
  y <- c(.2, -.12)
  loading <- c(.4, .55)
  variance <- c(.35, .47)^2
  V <- diag(variance) + tcrossprod(loading)
  basis <- eigen(V, symmetric = TRUE)
  quadrature <- .selection_joint_cluster_quadrature_rules(c(1L, 3L, 7L, 15L, 31L, 63L))
  qmc <- BayesTools::selection_qmc_design(2L, 512L, 8L, 9049L)
  data <- c(list(y = y, variance = variance, loading = loading, within_sd = .2,
    se = sqrt(diag(V)), z_lower = c(0, -1e300), z_upper = c(1e300, 0), bin = c(1, 2),
    rows = c(y, sqrt(diag(V)), diag(V)), eigenvalues = basis$values,
    eigenvectors = as.double(basis$vectors), gate = gate, declaration = 0,
    qmc = as.numeric(qmc)), quadrature[c("nodes", "log_weights", "orders")])
  if (nzchar(declaration)) data$foreign <- 0

  syntax <- paste0("model {
    bias ~ dcat(c(1, 1))
    gate ~ dcat(c(1, 1))
    raw ~ dnorm(0, 2)
    selected <- equals(bias, 2)
    weights[1] <- 1
    weights[2] <- 1 - .75 * selected
    for (i in 1:2) {
      location[i] <- equals(gate, 2) * raw + .1 * i
      sd[i] <- sqrt(variance[i] + pow(within_sd, 2))
    }
    y[1:2] ~ dselnorm_cluster_step(location, sd, loading, se, weights, z_lower, z_upper,
      bin, 1, 1, selected, nodes, log_weights, orders, qmc, 64, 512, 8, 0.005, 0)
    declaration ~ drobma_joint_correlated(bias, raw, gate, 2, ", nuisance, ",
      c(.1, .2), 0, c(0, 0), c(0, 1), c(0, 0), c(1, 1, 1, .25), rows, eigenvalues,
      eigenvectors, z_lower, z_upper)
    ", declaration, "
  }")
  inits <- lapply(seq_len(chains), function(i) {
    list(bias = 1L, raw = 0, .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2707L + i)
  })
  con <- textConnection(syntax)
  on.exit(close(con))
  rjags::jags.model(con, data = data, inits = inits, n.chains = chains, n.adapt = 0L, quiet = TRUE)
}

test_that("the correlated sampler owns only the label and raw intercept on both gate strata", {

  skip_if_not_installed("rjags")
  for (gate in 1:2) {
    model <- .correlated_bias_native_fit(gate, chains = 2L)
    samplers <- rjags::list.samplers(model)
    expect_named(samplers, "RoBMA::CorrelatedBiasJoint")
    expect_setequal(samplers[[1L]], c("bias", "raw"))
    draws <- rjags::coda.samples(model, c("bias", "raw", "location"), n.iter = 20L,
                               progress.bar = "none")
    for (chain in draws) {
      chain <- as.matrix(chain)
      expect_equal(chain[, "location[1]"], .1 + (gate == 2) * chain[, "raw"])
      expect_true(all(chain[, "bias"] %in% 1:2))
      expect_gt(length(unique(chain[, "raw"])), 10L)
    }
  }
})

test_that("the correlated sampler leaves dependent roles and overlapping blocks to ordinary owners", {

  skip_if_not_installed("rjags")
  declined <- list(
    .correlated_bias_native_fit(2L, nuisance = "within_sd + .1 * raw"),
    .correlated_bias_native_fit(2L, declaration = "foreign ~ dbt_block(raw, 0)")
  )
  for (model in declined) {
    samplers <- rjags::list.samplers(model)
    expect_false("RoBMA::CorrelatedBiasJoint" %in% names(samplers))
    expect_setequal(unlist(samplers, use.names = FALSE), c("bias", "raw"))
  }
})
