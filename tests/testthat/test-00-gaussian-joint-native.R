test_that("Gaussian joint ownership requires the entire independent likelihood graph", {

  skip_on_cran()
  skip_if_not_installed("rjags")

  make_model <- function(location, precision = "1/(pow(se[i], 2) + pow(tau, 2))") {

    model_syntax <- paste0(
      "model {\n",
      "  g0 ~ dcat(c(2, 1))\n",
      "  raw0 ~ dnorm(0.3, 2)\n",
      "  intercept <- raw0 * (g0 == 1) + 0 * (g0 == 2)\n",
      "  g1 ~ dcat(c(1, 2))\n",
      "  raw1 ~ dnorm(-0.2, 1)\n",
      "  slope <- 0 * (g1 == 1) + raw1 * (g1 == 2)\n",
      "  gt ~ dcat(c(1, 2))\n",
      "  rawt ~ dnorm(0, 3) T(0,)\n",
      "  tau <- rawt * (gt == 1) + 0 * (gt == 2)\n",
      "  for (i in 1:K) {\n",
      "    linear[i] <- ", location, "\n",
      "    y[i] ~ dnorm(-linear[i], ", precision, ")\n",
      "  }\n",
      "  marker ~ drobma_gaussian_joint(g0, raw0, g1, raw1, gt, rawt,\n",
      "    active[1:3], x[1:K], y[1:K], se[1:K], -1)\n",
      "}"
    )
    connection <- textConnection(model_syntax)
    on.exit(close(connection), add = TRUE)
    rjags::jags.model(connection,
      data = list(y = c(-0.5, 0.2, 0.1, -0.4), x = c(-1, -.4, .2, .9),
                  se = c(.4, .6, .5, .7), K = 4L, active = c(1L, 2L, 1L), marker = 0),
      inits = list(g0 = 1L, g1 = 2L, gt = 1L, raw0 = .2, raw1 = 0, rawt = .25,
                   .RNG.name = "base::Mersenne-Twister", .RNG.seed = 823L),
      n.chains = 1L, n.adapt = 0L, quiet = TRUE)
  }

  ordinary_location <- "intercept + slope * x[i]"
  owned <- make_model(ordinary_location)
  assigned <- rjags::list.samplers(owned)
  expect_identical(names(assigned), "RoBMA::GaussianJoint")
  expect_setequal(assigned[[1L]], c("g0", "raw0", "g1", "raw1", "gt", "rawt"))

  monitored <- c("g0", "raw0", "g1", "raw1", "gt", "rawt", "intercept", "slope", "tau")
  draws <- as.matrix(rjags::coda.samples(owned, monitored, 200L, progress.bar = "none"))
  expect_equal(draws[, "intercept"], draws[, "raw0"] * (draws[, "g0"] == 1))
  expect_equal(draws[, "slope"], draws[, "raw1"] * (draws[, "g1"] == 2))
  expect_equal(draws[, "tau"], draws[, "rawt"] * (draws[, "gt"] == 1))
  expect_gt(length(unique(draws[draws[, "gt"] == 2, "rawt"])), 1L)

  # The extra term is zero at several ordinary sentinel values but nonlinear.
  nonlinear <- paste0(ordinary_location, " + intercept * (intercept - 1) * (intercept + 0.5)")
  wrong_location <- make_model(nonlinear)
  wrong_variance <- make_model(ordinary_location,
    precision = "1/(pow(se[i], 2) + pow(tau, 2) + pow(tau * (tau - 0.5) * (tau - 1.25), 2))")
  for (fit in list(wrong_location, wrong_variance)) {
    assigned <- rjags::list.samplers(fit)
    expect_false("RoBMA::GaussianJoint" %in% names(assigned))
    expect_setequal(unlist(assigned, use.names = FALSE), c("g0", "raw0", "g1", "raw1", "gt", "rawt"))
  }
})
