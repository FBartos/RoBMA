context("(2) Distribution functions")
skip_on_cran()

### weighted normal distributions ----
test_that("Input checks work", {

  x <- seq(-5, 5, .5)
  expect_error(dwnorm(x, 0, -1, steps = c(.5), omega = c(1, 1), type = "one.sided"), "The 'sd' must be equal or higher than 0.")
  expect_error(dwnorm(x, 0,  1, steps = c(1), omega = c(1, 1), type = "one.sided"), "The 'steps' must be lower than 1.")
  expect_error(dwnorm(x, 0,  1, steps = c(0), omega = c(1, 1), type = "one.sided"), "The 'steps' must be higher than 0.")
  expect_error(dwnorm(x, 0,  1, steps = c(.5), omega = c(2, 1), type = "one.sided"), "The 'omega' must be equal or lower than 1.")
  expect_error(dwnorm(x, 0,  1, steps = c(.5), omega = c(-1, 1), type = "one.sided"), "The 'omega' must be equal or higher than 0.")
  expect_error(dwnorm(x, 0,  1, steps = c(.05, .10), omega = c(1, 1), type = "one.sided"), "'omega' argument must have one more weight than the number of defined steps with 'steps'/'crit_x' argument.")
  expect_error(dwnorm(x, 0,  1, steps = c(.10, .05), omega = c(1, 1, 1), type = "one.sided"), "'steps'/'crit_x' argument must be inreasing.")
  expect_error(dwnorm(x, 0,  1, crit_x = c(0, 1), omega = c(1, 1), type = "one.sided"), "'omega' argument must have one more weight than the number of defined steps with 'steps'/'crit_x' argument.")
  expect_error(dwnorm(x, 0,  1, crit_x = c(1, 0), omega = c(1, 1, 1), type = "one.sided"), "'steps'/'crit_x' argument must be inreasing.")
  expect_error(dwnorm(x, 0,  1, crit_x = c(-1, 1), omega = c(1, 1, 1), type = "two.sided"), "The 'crit_x' must be equal or higher than 0.")

})

test_that("Density function works", {

  x <- seq(-5, 5, .5)

  # verify no-weights against standard normal
  expect_equal(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided"), dnorm(x, 0, 1))
  expect_equal(dwnorm(x, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided"), dnorm(x, 1, 2))
  expect_equal(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided", log = TRUE), dnorm(x, 0, 1, log = TRUE))
  expect_equal(dwnorm(x, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided", log = TRUE), dnorm(x, 1, 2, log = TRUE))

  # verify that steps and crit_x are exchangeable
  expect_equal(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided"), dwnorm(x, 0, 1, crit_x = qnorm(.5 / 2, lower.tail = FALSE), omega = c(1, 1), type = "one.sided"))
  expect_equal(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, 1), type = "two.sided"), dwnorm(x, 0, 1, crit_x = qnorm(.5, lower.tail = FALSE), omega = c(1, 1), type = "two.sided"))

  # verify that the distributions integrate to 1 (up to some tolerance)
  expect_equal(integrate(function(x)dwnorm(x, 0, 1, steps = c(.50), omega = c(1, .5), type = "one.sided"), lower = -Inf, upper = Inf)$value, 1, tolerance = 1e-5)
  expect_equal(integrate(function(x)dwnorm(x, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided"), lower = -Inf, upper = Inf)$value, 1, tolerance = 1e-5)
  expect_equal(integrate(function(x)dwnorm(x, 0, 1, steps = c(.5), omega = c(1, .5), type = "two.sided"), lower = -Inf, upper = Inf)$value, 1, tolerance = 1e-5)
  expect_equal(integrate(function(x)dwnorm(x, 0.5, 1, steps = c(.05, .10), omega = c(1, .5, .1), type = "two.sided"), lower = -Inf, upper = Inf)$value, 1, tolerance = 1e-4)

  # visual inspection
  vdiffr::expect_doppelganger("dwnorm-1", function()curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), from = -3, to = 3))
  vdiffr::expect_doppelganger("dwnorm-2", function()curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided"), from = -3, to = 3))
  vdiffr::expect_doppelganger("dwnorm-3", function()curve(dwnorm(x, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided"), from = -3, to = 3))
  vdiffr::expect_doppelganger("dwnorm-4", function()curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, .5), type = "two.sided"), from = -3, to = 3))
  vdiffr::expect_doppelganger("dwnorm-5", function()curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided"), from = -3, to = 3))
  vdiffr::expect_doppelganger("dwnorm-6", function()curve(dwnorm(x, 0.5, 1, steps = c(.05, .10), omega = c(1, .5, .1), type = "two.sided"), from = -3, to = 3))

  # verify vectorization
  set.seed(1)
  x1     <- rnorm(10)
  mean1  <- rnorm(10)
  sd1    <- runif(10, 0.2, 1)
  steps1 <- cbind(runif(10, 0, .5), runif(10, .5, 1))
  omega1 <- cbind(1, runif(10, 0, .5), runif(10, 0, .3))

  expect_equal(
    dwnorm(x1, mean1, sd1, steps = steps1, omega = omega1, type = "one.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[i], sd1[i], steps = steps1[i,], omega = omega1[i,], type = "one.sided"))))
  expect_equal(
    dwnorm(x1, mean1[1], sd1, steps = steps1, omega = omega1, type = "one.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[1], sd1[i], steps = steps1[i,], omega = omega1[i,], type = "one.sided"))))
  expect_equal(
    dwnorm(x1, mean1, sd1[1], steps = steps1, omega = omega1, type = "one.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[i], sd1[1], steps = steps1[i,], omega = omega1[i,], type = "one.sided"))))
  expect_equal(
    dwnorm(x1, mean1, sd1, steps = steps1[1,], omega = omega1, type = "one.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[i], sd1[i], steps = steps1[1,], omega = omega1[i,], type = "one.sided"))))
  expect_equal(
    dwnorm(x1, mean1, sd1, steps = steps1, omega = omega1[1,], type = "one.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[i], sd1[i], steps = steps1[i,], omega = omega1[1,], type = "one.sided"))))

  expect_equal(
    dwnorm(x1, mean1, sd1, steps = steps1, omega = omega1, type = "two.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[i], sd1[i], steps = steps1[i,], omega = omega1[i,], type = "two.sided"))))
  expect_equal(
    dwnorm(x1, mean1[1], sd1, steps = steps1, omega = omega1, type = "two.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[1], sd1[i], steps = steps1[i,], omega = omega1[i,], type = "two.sided"))))
  expect_equal(
    dwnorm(x1, mean1, sd1[1], steps = steps1, omega = omega1, type = "two.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[i], sd1[1], steps = steps1[i,], omega = omega1[i,], type = "two.sided"))))
  expect_equal(
    dwnorm(x1, mean1, sd1, steps = steps1[1,], omega = omega1, type = "two.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[i], sd1[i], steps = steps1[1,], omega = omega1[i,], type = "two.sided"))))
  expect_equal(
    dwnorm(x1, mean1, sd1, steps = steps1, omega = omega1[1,], type = "two.sided"),
    do.call(c, lapply(1:length(x1), function(i) dwnorm(x1[i], mean1[i], sd1[i], steps = steps1[i,], omega = omega1[1,], type = "two.sided"))))

  # verify fast computation
  expect_equal(
    RoBMA:::.dwnorm_fast.bridge(x1, mean1, sd1, omega1[1,], matrix(c(0.50, 1.50), nrow = 10, ncol = 2, byrow = TRUE), type = "one.sided", log = TRUE),
    dwnorm(x1, mean1, sd1, crit_x = matrix(c(0.50, 1.50), nrow = 10, ncol = 2, byrow = TRUE), omega = omega1[1,], type = "one.sided", log = TRUE)
  )
  expect_equal(
    RoBMA:::.dwnorm_fast.bridge(x1, mean1, sd1, omega1[1,], matrix(c(0.50, 1.50), nrow = 10, ncol = 2, byrow = TRUE), type = "two.sided", log = TRUE),
    dwnorm(x1, mean1, sd1, crit_x = matrix(c(0.50, 1.50), nrow = 10, ncol = 2, byrow = TRUE), omega = omega1[1,], type = "two.sided", log = TRUE)
  )
})

test_that("Random number generator function works", {
  skip_on_os("linux")
  set.seed(1)

  # verify (visually) no-weights against standard normal
  expect_equal(rwnorm(0, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided"), rnorm(0, 0, 1))
  vdiffr::expect_doppelganger("rwnorm-1", function(){
    plot(density(rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided")))
    curve(dnorm(x, 0, 1), from = -10, to = 10, col = "blue", add = TRUE)
  })
  vdiffr::expect_doppelganger("rwnorm-2", function(){
    plot(density(rwnorm(10000, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided")))
    curve(dnorm(x, 1, 2), from = -10, to = 10, col = "blue", add = TRUE)
  })

  # visually verify against pdf
  vdiffr::expect_doppelganger("rwnorm-3", function(){
    hist(rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), from = -10, to = 10, col = "blue", add = TRUE)
  })
  vdiffr::expect_doppelganger("rwnorm-4", function(){
    hist(rwnorm(10000, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided"), from = -10, to = 10, col = "blue", add = TRUE)
    })
  vdiffr::expect_doppelganger("rwnorm-5", function(){
    hist(rwnorm(10000, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided"), from = -10, to = 10, col = "blue", add = TRUE)
    })
  vdiffr::expect_doppelganger("rwnorm-6", function(){
    hist(rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, .5), type = "two.sided"), from = -10, to = 10, col = "blue", add = TRUE)
    })
  vdiffr::expect_doppelganger("rwnorm-7", function(){
    hist(rwnorm(10000, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided"), from = -10, to = 10, col = "blue", add = TRUE)
    })
  vdiffr::expect_doppelganger("rwnorm-8", function(){
    hist(rwnorm(10000, 0.5, 1, steps = c(.05, .10), omega = c(1, .5, .1), type = "two.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 0.5, 1, steps = c(.05, .10), omega = c(1, .5, .1), type = "two.sided"), from = -10, to = 10, col = "blue", add = TRUE)
    })
})

test_that("Distribution function works", {

  q <- seq(-5, 5, .5)

  # verify no-weights against standard normal
  expect_equal(pwnorm(q, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided"), pnorm(q, 0, 1))
  expect_equal(pwnorm(q, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided"), pnorm(q, 1, 2))
  expect_equal(pwnorm(q, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided", log = TRUE), pnorm(q, 0, 1, log = TRUE))
  expect_equal(pwnorm(q, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided", log = TRUE), pnorm(q, 1, 2, log = TRUE))
  expect_equal(pwnorm(q, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided", lower.tail = TRUE), pnorm(q, 0, 1, lower.tail = TRUE))
  expect_equal(pwnorm(q, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided", lower.tail = TRUE), pnorm(q, 1, 2, lower.tail = TRUE))
  expect_equal(pwnorm(q, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided", log = TRUE, lower.tail = TRUE), pnorm(q, 0, 1, log = TRUE, lower.tail = TRUE))
  expect_equal(pwnorm(q, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided", log = TRUE, lower.tail = TRUE), pnorm(q, 1, 2, log = TRUE, lower.tail = TRUE))

  # verify against rng
  set.seed(1)
  expect_equal({
    samples <- rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided")
    sapply(1:length(q), function(i) mean(samples < q[i]))
  }, pwnorm(q, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), tolerance = 1e-2)
  expect_equal({
    samples <- rwnorm(10000, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided")
    sapply(1:length(q), function(i) mean(samples < q[i]))
  }, pwnorm(q, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided"), tolerance = 1e-2)
  expect_equal({
    samples <- rwnorm(10000, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided")
    sapply(1:length(q), function(i) mean(samples < q[i]))
  }, pwnorm(q, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided"), tolerance = 1e-2)
  expect_equal({
    samples <- rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, .5), type = "two.sided")
    sapply(1:length(q), function(i) mean(samples < q[i]))
  }, pwnorm(q, 0, 1, steps = c(.5), omega = c(1, .5), type = "two.sided"), tolerance = 1e-2)
  expect_equal({
    samples <- rwnorm(10000, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided")
    sapply(1:length(q), function(i) mean(samples < q[i]))
  }, pwnorm(q, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided"), tolerance = 1e-2)
  expect_equal({
    samples <- rwnorm(10000, 0.5, 1, steps = c(.05, .10), omega = c(1, .5, .1), type = "two.sided")
    sapply(1:length(q), function(i) mean(samples < q[i]))
  }, pwnorm(q, 0.5, 1, steps = c(.05, .10), omega = c(1, .5, .1), type = "two.sided"), tolerance = 1e-2)

})

test_that("Quantile function work", {

  p <- seq(0, 1, .1)

  # verify no-weights against standard normal
  expect_equal(qwnorm(p, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided"), qnorm(p, 0, 1), tolerance = 1e-5)
  expect_equal(qwnorm(p, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided"), qnorm(p, 1, 2), tolerance = 1e-5)
  expect_equal(qwnorm(log(p), 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided", log.p = TRUE), qnorm(log(p), 0, 1, log.p = TRUE), tolerance = 1e-5)
  expect_equal(qwnorm(log(p), 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided", log.p = TRUE), qnorm(log(p), 1, 2, log.p = TRUE), tolerance = 1e-5)
  expect_equal(qwnorm(p, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided", lower.tail = TRUE), qnorm(p, 0, 1), tolerance = 1e-5, lower.tail = TRUE)
  expect_equal(qwnorm(p, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided", lower.tail = TRUE), qnorm(p, 1, 2), tolerance = 1e-5, lower.tail = TRUE)
  expect_equal(qwnorm(log(p), 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided", log.p = TRUE, lower.tail = TRUE), qnorm(log(p), 0, 1, log.p = TRUE), tolerance = 1e-5, lower.tail = TRUE)
  expect_equal(qwnorm(log(p), 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided", log.p = TRUE, lower.tail = TRUE), qnorm(log(p), 1, 2, log.p = TRUE), tolerance = 1e-5, lower.tail = TRUE)

  # verify against rng
  p <- seq(0.1, 0.9, 0.1)
  set.seed(1)
  expect_equal({
    samples <- rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided")
    unname(quantile(samples, p))
  }, qwnorm(p, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), tolerance = 0.05)
  expect_equal({
    samples <- rwnorm(10000, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided")
    unname(quantile(samples, p))
  }, qwnorm(p, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided"), tolerance = 0.05)
  expect_equal({
    samples <- rwnorm(10000, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided")
    unname(quantile(samples, p))
  }, qwnorm(p, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided"), tolerance = 0.05)
  expect_equal({
    samples <- rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, .5), type = "two.sided")
    unname(quantile(samples, p))
  }, qwnorm(p, 0, 1, steps = c(.5), omega = c(1, .5), type = "two.sided"), tolerance = 0.05)
  expect_equal({
    samples <- rwnorm(10000, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided")
    unname(quantile(samples, p))
  }, qwnorm(p, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided"), tolerance = 0.05)
  expect_equal({
    samples <- rwnorm(10000, 0.5, 1, steps = c(.05, .10), omega = c(1, .5, .1), type = "two.sided")
    unname(quantile(samples, p))
  }, qwnorm(p, 0.5, 1, steps = c(.05, .10), omega = c(1, .5, .1), type = "two.sided"), tolerance = 0.05)

})

### multivariate weighted normal distributions ----

test_that("Density function works", {

  set.seed(1)
  mu    <- c(0.2, 0.5, 0.15)
  sigma <- matrix(c(
    1.2, 0.3, 0.1,
    0.3, 0.8, 0.2,
    0.1, 0.2, 1.1), nrow = 3, ncol = 3)
  omega <- c(1, 1)
  x     <- mvtnorm::rmvnorm(10, mu, sigma)


  # verify no-weights against standard normal
  expect_equal(
    sapply(1:nrow(x), function(i){
      RoBMA:::.dwmnorm_fast(x = x[i,], mean = mu, sigma = sigma, omega = omega, crit_x = matrix(1, nrow = 1, ncol = 3), type = "one.sided", log = TRUE)
    }),
    mvtnorm::dmvnorm(x = x, mean = mu, sigma = sigma, log = TRUE),
    tolerance = 1e-4
  )

  expect_equal(
    sapply(1:nrow(x), function(i){
      RoBMA:::.dwmnorm_fast(x = x[i,], mean = mu, sigma = sigma, omega = omega, crit_x = matrix(1, nrow = 1, ncol = 3), type = "two.sided", log = TRUE)
    }),
    mvtnorm::dmvnorm(x = x, mean = mu, sigma = sigma, log = TRUE),
    tolerance = 1e-4
  )

  # verify against independent univariate
  sigma  <- diag(diag(sigma), ncol(sigma))
  omega  <- c(0.25, 1)
  crit_x <- matrix(c(
    1.25,
    1.30,
    0.80), nrow = 1, ncol = 3)

  expect_equal(
    sapply(1:nrow(x), function(i){
      RoBMA:::.dwmnorm_fast(x = x[i,], mean = mu, sigma = sigma, omega = omega, crit_x = crit_x, type = "one.sided", log = TRUE)
    }),
    sapply(1:nrow(x), function(i){
      sum(RoBMA:::.dwnorm_fast.bridge(x = x[i,], mean = mu, sd = sqrt(diag(sigma)), omega = omega, crit_x = t(crit_x), type = "one.sided", log = TRUE))
    }),
    tolerance = 1e-4
  )

  omega  <- c(0.25, .50,  1)
  crit_x <- matrix(c(
    1.25, 1.50,
    1.30, 1.45,
    0.80, 1.00), nrow = 2, ncol = 3)

  expect_equal(
    sapply(1:nrow(x), function(i){
      RoBMA:::.dwmnorm_fast(x = x[i,], mean = mu, sigma = sigma, omega = omega, crit_x = crit_x, type = "one.sided", log = TRUE)
    }),
    sapply(1:nrow(x), function(i){
      sum(RoBMA:::.dwnorm_fast.bridge(x = x[i,], mean = mu, sd = sqrt(diag(sigma)), omega = omega, crit_x = t(crit_x), type = "one.sided", log = TRUE))
    }),
    tolerance = 1e-4
  )

  expect_equal(
    sapply(1:nrow(x), function(i){
      RoBMA:::.dwmnorm_fast(x = x[i,], mean = mu, sigma = sigma, omega = omega, crit_x = crit_x, type = "two.sided", log = TRUE)
    }),
    sapply(1:nrow(x), function(i){
      sum(RoBMA:::.dwnorm_fast.bridge(x = x[i,], mean = mu, sd = sqrt(diag(sigma)), omega = omega, crit_x = t(crit_x), type = "two.sided", log = TRUE))
    }),
    tolerance = 1e-4
  )

})

test_that("R and JAGS density is consistent", {

  # re-load the module
  RoBMA:::.load_RoBMA_module()

  ### one sided
  # univariate (simple)
  set.seed(1)
  model_syntax <-
    "model
    {
      x ~ dnorm(0, 1)

      mu ~ dnorm(0, pow(0.30, -2))
      sigma ~ dunif(0.1, 1.5)

      omega[1] ~ dunif(0, 1)
      omega[2] ~ dunif(0, 1)

      log_lik = wnorm_1s_lpdf(x, mu, sigma, crit_x, omega)
    }"

  data <- list(
    crit_x = matrix(c(0.25), nrow = 1, ncol = 1)
  )

  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
  fit   <- rjags::coda.samples(model = model, variable.names = c("x", "omega", "mu", "sigma", "log_lik"), n.iter = 10, quiet = TRUE, progress.bar = "none")


  expect_equal(as.vector(fit[[1]][,"log_lik"]), unname(sapply(1:10, function(i){
    RoBMA:::.dwnorm_fast.bridge(x = fit[[1]][i,"x"], mean = fit[[1]][i,"mu"], sd = fit[[1]][i,"sigma"], omega = fit[[1]][i,c("omega[1]", "omega[2]")], crit_x = data$crit_x, type = "one.sided", log = TRUE)
  })), tolerance = 1e-3)

  # univariate (more complex)
  set.seed(1)
  model_syntax <-
    "model
    {
      x ~ dnorm(0, 1)

      mu ~ dnorm(0, pow(0.30, -2))
      sigma ~ dunif(0.1, 1.5)

      omega[1] ~ dunif(0, 1)
      omega[2] ~ dunif(0, 1)
      omega[3] ~ dunif(0, 1)

      log_lik = wnorm_1s_lpdf(x, mu, sigma, crit_x, omega)
    }"

  data <- list(
    crit_x = matrix(c(
      0.25,
      0.70), nrow = 1, ncol = 2)
  )

  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
  fit   <- rjags::coda.samples(model = model, variable.names = c("x", "omega", "mu", "sigma", "log_lik"), n.iter = 10, quiet = TRUE, progress.bar = "none")


  expect_equal(as.vector(fit[[1]][,"log_lik"]), unname(sapply(1:10, function(i){
    RoBMA:::.dwnorm_fast.bridge(x = fit[[1]][i,"x"], mean = fit[[1]][i,"mu"], sd = fit[[1]][i,"sigma"], omega = fit[[1]][i,c("omega[1]", "omega[2]", "omega[3]")], crit_x = data$crit_x, type = "one.sided", log = TRUE)
  })), tolerance = 1e-3)

  # multivariate
  set.seed(1)
  model_syntax <-
    "model
    {
      x[1] ~ dnorm(0, 1)
      x[2] ~ dnorm(0, 1)
      x[3] ~ dnorm(0, 1)

      mu[1] ~ dnorm(0, pow(0.30, -2))
      mu[2] ~ dnorm(0, pow(0.30, -2))
      mu[3] ~ dnorm(0, pow(0.30, -2))

      omega[1] ~ dunif(0, 1)
      omega[2] ~ dunif(0, 1)
      omega[3] ~ dunif(0, 1)

      log_lik = wmnorm_1s_lpdf(x, mu, sigma, crit_x, omega)
    }"

  data <- list(
    sigma  = matrix(c(
      1.5, 1.0, 0.5,
      1.0, 1.8, 0.7,
      0.5, 0.7, 1.2), nrow = 3, ncol = 3),
    crit_x = matrix(c(
      1.25, 1.96,
      1.30, 2.05,
      1.10, 1.50), nrow = 2, ncol = 3)
  )

  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
  fit   <- rjags::coda.samples(model = model, variable.names = c("x", "omega", "mu", "log_lik"), n.iter = 10, quiet = TRUE, progress.bar = "none")

  expect_equal(as.vector(fit[[1]][,"log_lik"]), sapply(1:10, function(i){
    RoBMA:::.dwmnorm_fast(x = fit[[1]][i,c("x[1]", "x[2]", "x[3]")], mean = fit[[1]][i,c("mu[1]", "mu[2]", "mu[3]")], sigma = data$sigma, omega = fit[[1]][i,c("omega[1]", "omega[2]", "omega[3]")], crit_x = data$crit_x, type = "one.sided", log = TRUE)
  }), tolerance = 1e-3)


  ### two sided
  # univariate
  set.seed(1)
  model_syntax <-
    "model
    {
      x ~ dnorm(0, 1)

      mu ~ dnorm(0, pow(0.30, -2))
      sigma ~ dunif(0.1, 1.5)

      omega[1] ~ dunif(0, 1)
      omega[2] ~ dunif(0, 1)
      omega[3] ~ dunif(0, 1)

      log_lik = wnorm_2s_lpdf(x, mu, sigma, crit_x, omega)
    }"

  data <- list(
    crit_x = matrix(c(
      0.25,
      0.70), nrow = 1, ncol = 2)
  )

  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
  fit   <- rjags::coda.samples(model = model, variable.names = c("x", "omega", "mu", "sigma", "log_lik"), n.iter = 10, quiet = TRUE, progress.bar = "none")


  expect_equal(as.vector(fit[[1]][,"log_lik"]), unname(sapply(1:10, function(i){
    RoBMA:::.dwnorm_fast.bridge(x = fit[[1]][i,"x"], mean = fit[[1]][i,"mu"], sd = fit[[1]][i,"sigma"], omega = fit[[1]][i,c("omega[1]", "omega[2]", "omega[3]")], crit_x = data$crit_x, type = "two.sided", log = TRUE)
  })), tolerance = 1e-3)

  # multivariate
  set.seed(1)
  model_syntax <-
    "model
    {
      x[1] ~ dnorm(0, 1)
      x[2] ~ dnorm(0, 1)
      x[3] ~ dnorm(0, 1)
      x[4] ~ dnorm(0, 1)

      mu[1] ~ dnorm(0, pow(0.30, -2))
      mu[2] ~ dnorm(0, pow(0.30, -2))
      mu[3] ~ dnorm(0, pow(0.30, -2))
      mu[4] ~ dnorm(0, pow(0.30, -2))

      omega[1] ~ dunif(0, 1)
      omega[2] ~ dunif(0, 1)

      log_lik = wmnorm_2s_lpdf(x, mu, sigma, crit_x, omega)
    }"

  data <- list(
    sigma  = matrix(c(
      1.5, 1.0, 0.5, 0.2,
      1.0, 1.8, 0.7, 0.5,
      0.5, 0.7, 1.2, 0.4,
      0.2, 0.5, 0.4, 0.8), nrow = 4, ncol = 4),
    crit_x = matrix(c(
      1.25,
      1.30,
      1.10,
      0.80), nrow = 1, ncol = 4)
  )

  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
  fit   <- rjags::coda.samples(model = model, variable.names = c("x", "omega", "mu", "log_lik"), n.iter = 10, quiet = TRUE, progress.bar = "none")

  expect_equal(as.vector(fit[[1]][,"log_lik"]), sapply(1:10, function(i){
    RoBMA:::.dwmnorm_fast(x = fit[[1]][i,c("x[1]", "x[2]", "x[3]", "x[4]")], mean = fit[[1]][i,c("mu[1]", "mu[2]", "mu[3]", "mu[4]")], sigma = data$sigma, omega = fit[[1]][i,c("omega[1]", "omega[2]")], crit_x = data$crit_x, type = "two.sided", log = TRUE)
  }), tolerance = 1e-3)


  ### two sided: normal vs. as vector
  set.seed(1)
  model_syntax <-
    "model
    {
      omega[1] ~ dunif(0, 1)
      omega[2] = 1

      log_lik = wmnorm_2s_v_lpdf(x, mu, se2, tau2, rho2, crit_x, omega, indx)
    }"

  data <- list(
    x    = rnorm(5),
    mu   = c(1.1, 1.2, 2.1, 2.2, 2.3),
    se2  = c(.11, .12, .21, .22, .23),
    tau2 = .05,
    rho2 = .30^2,
    crit_x = matrix(c(
      1.01,
      1.02,
      2.01,
      2.02,
      2.03), nrow = 1, ncol = 5),
    indx = c(2,5)
  )

  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE, n.adapt = 10, inits = list(.RNG.seed = 1, .RNG.name = "base::Super-Duper"))
  fit   <- rjags::coda.samples(model = model, variable.names = "log_lik", n.iter = 10, quiet = TRUE, progress.bar = "none")


  set.seed(1)
  model_syntax2 <-
    "model
    {
      omega[1] ~ dunif(0, 1)
      omega[2] = 1

      log_lik[1] = wmnorm_2s_lpdf(x1, mu1, sigma1, crit_x1, omega)
      log_lik[2] = wmnorm_2s_lpdf(x2, mu2, sigma2, crit_x2, omega)

    }"

  data2 <- list(
    x1      = data$x[1:2],
    x2      = data$x[3:5],
    mu1     = data$mu[1:2],
    mu2     = data$mu[3:5],
    sigma1  = diag(data$se2[1:2], 2) + diag(data$tau2 * (1-data$rho2), 2) + matrix(data$tau2 * data$rho2, 2, 2),
    sigma2  = diag(data$se2[3:5], 3) + diag(data$tau2 * (1-data$rho2), 3) + matrix(data$tau2 * data$rho2, 3, 3),
    crit_x1 = data$crit_x[,1:2, drop = FALSE],
    crit_x2 = data$crit_x[,3:5, drop = FALSE]
  )

  model2 <- rjags::jags.model(file = textConnection(model_syntax2), data = data2, quiet = TRUE, n.adapt = 10, inits = list(.RNG.seed = 1, .RNG.name = "base::Super-Duper"))
  fit2   <- rjags::coda.samples(model = model2, variable.names = "log_lik", n.iter = 10, quiet = TRUE, progress.bar = "none")

  expect_equal(apply(fit2[[1]], 1, sum), as.vector(fit[[1]]), tolerance = 1e-4)


  ### one sided: normal vs. as vector
  set.seed(1)
  model_syntax <-
    "model
    {
      omega[1] ~ dunif(0.0, 0.5)
      omega[2] ~ dunif(0.5, 1.0)
      omega[3] = 1

      log_lik = wmnorm_1s_v_lpdf(x, mu, se2, tau2, rho2, crit_x, omega, indx)
    }"

  data <- list(
    x    = rnorm(5),
    mu   = c(1.1, 1.2, 2.1, 2.2, 2.3),
    se2  = c(.11, .12, .21, .22, .23),
    tau2 = .05,
    rho2 = .30^2,
    crit_x = matrix(c(
      1.25, 1.96,
      1.30, 2.05,
      1.10, 1.50,
      0.10, 0.50,
      0.50, 1.00), nrow = 2, ncol = 5),
    indx = c(2,5)
  )


  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE, n.adapt = 10, inits = list(.RNG.seed = 1, .RNG.name = "base::Super-Duper"))
  fit   <- rjags::coda.samples(model = model, variable.names = "log_lik", n.iter = 10, quiet = TRUE, progress.bar = "none")


  set.seed(1)
  model_syntax2 <-
    "model
    {
      omega[1] ~ dunif(0.0, 0.5)
      omega[2] ~ dunif(0.5, 1.0)
      omega[3] = 1

      log_lik[1] = wmnorm_1s_lpdf(x1, mu1, sigma1, crit_x1, omega)
      log_lik[2] = wmnorm_1s_lpdf(x2, mu2, sigma2, crit_x2, omega)

    }"

  data2 <- list(
    x1      = data$x[1:2],
    x2      = data$x[3:5],
    mu1     = data$mu[1:2],
    mu2     = data$mu[3:5],
    sigma1  = diag(data$se2[1:2], 2) + diag(data$tau2 * (1-data$rho2), 2) + matrix(data$tau2 * data$rho2, 2, 2),
    sigma2  = diag(data$se2[3:5], 3) + diag(data$tau2 * (1-data$rho2), 3) + matrix(data$tau2 * data$rho2, 3, 3),
    crit_x1 = data$crit_x[,1:2, drop = FALSE],
    crit_x2 = data$crit_x[,3:5, drop = FALSE]
  )

  model2 <- rjags::jags.model(file = textConnection(model_syntax2), data = data2, quiet = TRUE, n.adapt = 10, inits = list(.RNG.seed = 1, .RNG.name = "base::Super-Duper"))
  fit2   <- rjags::coda.samples(model = model2, variable.names = "log_lik", n.iter = 10, quiet = TRUE, progress.bar = "none")

  expect_equal(apply(fit2[[1]], 1, sum), as.vector(fit[[1]]), tolerance = 1e-4)


})

### other ----

test_that("Test re-sampling and inverse transformation fast version of weighted-normal distribution", {

  # Test parameters
  n        <- 5
  mean_val <- c(0, 0.3, 0.6, 0.9, 1.2)
  sd_val   <- rep(0.2, 5)
  omega    <- matrix(c(0.3, 0.7, 1.0), byrow = TRUE, nrow = n, ncol = 3)
  crit_x   <- c(0.2, 0.6)

  # Generate samples using .rwnorm_fast.ss
  samples_fast <- replicate(5000, .rwnorm_fast.ss(
    mean   = mean_val,
    sd     = sd_val,
    omega  = omega,
    crit_x = crit_x
  ))

  samples_fast2 <- replicate(5000, .rwnorm_fast.ss2(
    mean   = mean_val,
    sd     = sd_val,
    omega  = omega,
    crit_x = crit_x
  ))

  vdiffr::expect_doppelganger(paste0("plot_fast_rng"),        function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(2,3))
    for(i in 1:5){
      common_breaks <- seq(min(c(samples_fast[i,], samples_fast2[i,])), max(c(samples_fast[i,], samples_fast2[i,])), length.out = 20)
      common_breaks <- sort(unique(c(common_breaks, crit_x)))
      hist(samples_fast[i,],  breaks = common_breaks)
      hist(samples_fast2[i,], breaks = common_breaks, col = scales::alpha("orange", 10), add = TRUE)
    }
  })
})

### fast spike-and-slab weighted normal distributions ----

test_that("Fast spike-and-slab output preformated distribution functions", {

  # Test parameters (single row for visual checks)
  mean_val <- 0.2
  sd_val   <- 0.3
  omega    <- matrix(c(0.1, 0.5, 1.0), nrow = 1, ncol = 3)
  crit_x   <- c(0.4, 0.8)

  # Visual check 1: .dwnorm_fast.ss density curve
  vdiffr::expect_doppelganger("dwnorm_fast_ss-1", function(){
    x_vals <- seq(-2, 3, 0.01)
    y_vals <- sapply(x_vals, function(x) {
      RoBMA:::.dwnorm_fast.ss(x, mean_val, sd_val, omega, crit_x, log = FALSE)
    })
    # Compare with exported function
    y_exported <- dwnorm(x_vals, mean_val, sd_val, crit_x = crit_x, omega = as.vector(omega), type = "one.sided")

    plot(x_vals, y_vals, type = "l", main = "dwnorm_fast.ss vs exported dwnorm",
         xlab = "x", ylab = "density", col = "blue", lwd = 2)
    lines(x_vals, y_exported, col = "red", lty = 2, lwd = 2)
    abline(v = crit_x, col = "gray", lty = 3)
    legend("topright", c("fast.ss", "exported"), col = c("blue", "red"), lty = c(1, 2))
  })  # Visual check 2: .dwnorm_fast.ss vs normal (when omega = c(1,1,1))
  omega_uniform <- matrix(c(1.0, 1.0, 1.0), nrow = 1, ncol = 3)
  vdiffr::expect_doppelganger("dwnorm_fast_ss-2", function(){
    x_vals <- seq(-2, 3, 0.01)
    y_fast <- sapply(x_vals, function(x) {
      RoBMA:::.dwnorm_fast.ss(x, mean_val, sd_val, omega_uniform, crit_x, log = FALSE)
    })
    y_norm <- dnorm(x_vals, mean_val, sd_val)
    plot(x_vals, y_fast, type = "l", col = "blue", main = "dwnorm_fast.ss vs normal (uniform weights)",
         xlab = "x", ylab = "density")
    lines(x_vals, y_norm, col = "red", lty = 2)
    legend("topright", c("fast.ss", "normal"), col = c("blue", "red"), lty = c(1, 2))
  })

  # Visual check 3: .rwnorm_fast.ss histogram vs .dwnorm_fast.ss
  set.seed(123)
  vdiffr::expect_doppelganger("rwnorm_fast_ss-1", function(){
    samples_fast <- replicate(5000, RoBMA:::.rwnorm_fast.ss(mean_val, sd_val, omega, crit_x))
    samples_exported <- rwnorm(5000, mean_val, sd_val, crit_x = crit_x, omega = as.vector(omega), type = "one.sided")

    common_range <- range(c(samples_fast, samples_exported))
    common_breaks <- seq(common_range[1], common_range[2], length.out = 40)

    hist(samples_fast, breaks = common_breaks, freq = FALSE, col = scales::alpha("blue", 0.6),
         main = "rwnorm_fast.ss vs exported rwnorm", xlab = "x")
    hist(samples_exported, breaks = common_breaks, freq = FALSE, col = scales::alpha("red", 0.6), add = TRUE)

    # Add theoretical density
    x_vals <- seq(common_range[1], common_range[2], 0.01)
    y_vals <- sapply(x_vals, function(x) {
      RoBMA:::.dwnorm_fast.ss(x, mean_val, sd_val, omega, crit_x, log = FALSE)
    })
    lines(x_vals, y_vals, col = "black", lwd = 2)
    abline(v = crit_x, col = "gray", lty = 3)
    legend("topright", c("fast.ss", "exported", "theory"),
           fill = c("blue", "red", NA), border = c("black", "black", NA),
           lty = c(NA, NA, 1), col = c(NA, NA, "black"))
  })

  # Visual check 4: .rwnorm_fast.ss vs .rwnorm_fast.ss2 comparison
  set.seed(123)
  vdiffr::expect_doppelganger("rwnorm_fast_ss-2", function(){
    samples1 <- replicate(2000, RoBMA:::.rwnorm_fast.ss(mean_val, sd_val, omega, crit_x))
    samples2 <- replicate(2000, RoBMA:::.rwnorm_fast.ss2(mean_val, sd_val, omega, crit_x))

    common_range <- range(c(samples1, samples2))
    common_breaks <- seq(common_range[1], common_range[2], length.out = 30)

    hist(samples1, breaks = common_breaks, freq = FALSE, col = scales::alpha("blue", 0.5),
         main = "rwnorm_fast.ss vs rwnorm_fast.ss2", xlab = "x")
    hist(samples2, breaks = common_breaks, freq = FALSE, col = scales::alpha("red", 0.5), add = TRUE)
    legend("topright", c("fast.ss", "fast.ss2"), fill = c("blue", "red"))
    abline(v = crit_x, col = "black", lty = 2)
  })

  # Visual check 5: .pwnorm_fast.ss cumulative distribution
  vdiffr::expect_doppelganger("pwnorm_fast_ss-1", function(){
    q_vals <- seq(-2, 3, 0.05)
    p_fast <- sapply(q_vals, function(q) {
      exp(RoBMA:::.pwnorm_fast.ss(q, mean_val, sd_val, omega, crit_x))
    })
    p_exported <- pwnorm(q_vals, mean_val, sd_val, crit_x = crit_x, omega = as.vector(omega), type = "one.sided")

    plot(q_vals, p_fast, type = "l", main = "pwnorm_fast.ss vs exported pwnorm",
         xlab = "q", ylab = "P(X <= q)", col = "blue", lwd = 2)
    lines(q_vals, p_exported, col = "red", lty = 2, lwd = 2)
    abline(v = crit_x, col = "gray", lty = 3)
    abline(h = c(0, 1), col = "gray", lty = 3)
    legend("bottomright", c("fast.ss", "exported"), col = c("blue", "red"), lty = c(1, 2))
  })

  # Visual check 6: .qwnorm_fast.ss quantile function
  vdiffr::expect_doppelganger("qwnorm_fast_ss-1", function(){
    p_vals <- seq(0.01, 0.99, 0.01)
    q_fast <- sapply(p_vals, function(p) {
      RoBMA:::.qwnorm_fast.ss(p, mean_val, sd_val, omega, crit_x)
    })
    q_exported <- qwnorm(p_vals, mean_val, sd_val, crit_x = crit_x, omega = as.vector(omega), type = "one.sided")

    plot(p_vals, q_fast, type = "l", main = "qwnorm_fast.ss vs exported qwnorm",
         xlab = "p", ylab = "quantile", col = "blue", lwd = 2)
    lines(p_vals, q_exported, col = "red", lty = 2, lwd = 2)
    abline(h = crit_x, col = "gray", lty = 3)
    legend("topleft", c("fast.ss", "exported"), col = c("blue", "red"), lty = c(1, 2))
  })

  # Visual check 7: Different weight patterns effect on density
  omega_low_high <- matrix(c(0.2, 0.5, 1.0), nrow = 1, ncol = 3)  # increasing weights
  omega_high_low <- matrix(c(1.0, 0.5, 0.2), nrow = 1, ncol = 3)  # decreasing weights

  vdiffr::expect_doppelganger("dwnorm_fast_ss-3", function(){
    x_vals <- seq(-2, 3, 0.01)

    # Fast.ss versions
    y_fast_low_high <- sapply(x_vals, function(x) {
      RoBMA:::.dwnorm_fast.ss(x, mean_val, sd_val, omega_low_high, crit_x, log = FALSE)
    })
    y_fast_high_low <- sapply(x_vals, function(x) {
      RoBMA:::.dwnorm_fast.ss(x, mean_val, sd_val, omega_high_low, crit_x, log = FALSE)
    })

    # Exported versions for comparison
    y_exp_low_high <- dwnorm(x_vals, mean_val, sd_val, crit_x = crit_x, omega = as.vector(omega_low_high), type = "one.sided")
    y_exp_high_low <- dwnorm(x_vals, mean_val, sd_val, crit_x = crit_x, omega = as.vector(omega_high_low), type = "one.sided")

    plot(x_vals, y_fast_low_high, type = "l", col = "blue", lwd = 2, main = "Weight patterns: fast.ss vs exported",
         xlab = "x", ylab = "density")
    lines(x_vals, y_fast_high_low, col = "blue", lwd = 2)
    lines(x_vals, y_exp_low_high, col = "red", lty = 2, lwd = 1)
    lines(x_vals, y_exp_high_low, col = "red", lty = 2, lwd = 1)
    abline(v = crit_x, col = "gray", lty = 3)
    legend("topright", c("low→high (fast.ss)", "high→low (fast.ss)", "low→high (exported)", "high→low (exported)"),
           col = c("blue", "blue", "red", "red"), lty = c(1, 1, 2, 2), lwd = c(2, 2, 1, 1))
  })

  # Visual check 8: .rwnorm_true_fast.ss for hierarchical sampling
  set.seed(123)
  vdiffr::expect_doppelganger("rwnorm_true_fast_ss-1", function(){
    tau_val <- 0.15
    se_val  <- 0.25

    true_effects <- replicate(3000, RoBMA:::.rwnorm_true_fast.ss(mean_val, tau_val, se_val, omega, crit_x))

    hist(true_effects, freq = FALSE, breaks = 30,
         main = "rwnorm_true_fast.ss (hierarchical sampling)",
         xlab = "true effects")

    # Overlay the prior distribution of true effects
    x_vals <- seq(min(true_effects) - 0.5, max(true_effects) + 0.5, 0.01)
    prior_dens <- dnorm(x_vals, mean_val, tau_val)
    lines(x_vals, prior_dens, col = "blue", lwd = 2, lty = 2)

    abline(v = crit_x, col = "red", lty = 2)
    legend("topright", c("sampled true effects", "prior distribution"),
           col = c("black", "blue"), lty = c(1, 2))
  })

  # Visual check 9: Consistency check - p and q functions are inverses
  vdiffr::expect_doppelganger("pqwnorm_fast_ss-consistency", function(){
    # Test if p and q are inverses for both fast.ss and exported
    p_vals <- seq(0.05, 0.95, 0.05)

    # Fast.ss versions
    q_fast <- sapply(p_vals, function(p) RoBMA:::.qwnorm_fast.ss(p, mean_val, sd_val, omega, crit_x))
    p_recovered_fast <- sapply(q_fast, function(q) exp(RoBMA:::.pwnorm_fast.ss(q, mean_val, sd_val, omega, crit_x)))

    # Exported versions
    q_exported <- qwnorm(p_vals, mean_val, sd_val, crit_x = crit_x, omega = as.vector(omega), type = "one.sided")
    p_recovered_exported <- pwnorm(q_exported, mean_val, sd_val, crit_x = crit_x, omega = as.vector(omega), type = "one.sided")

    plot(p_vals, p_recovered_fast, main = "P-Q consistency: fast.ss vs exported",
         xlab = "original p", ylab = "recovered p",
         xlim = c(0, 1), ylim = c(0, 1), pch = 16, col = "blue")
    points(p_vals, p_recovered_exported, pch = 17, col = "red")
    abline(a = 0, b = 1, col = "gray", lty = 2)
    legend("bottomright", c("fast.ss", "exported", "perfect"),
           pch = c(16, 17, NA), col = c("blue", "red", "gray"),
           lty = c(NA, NA, 2))
  })

  # Visual check 10: Multiple critical values
  crit_x_multi <- c(0.2, 0.6, 1.2)
  omega_multi <- matrix(c(0.3, 0.6, 0.1, 1.0), nrow = 1, ncol = 4)

  vdiffr::expect_doppelganger("dwnorm_fast_ss-multi", function(){
    x_vals <- seq(-2, 3, 0.01)
    y_fast <- sapply(x_vals, function(x) {
      RoBMA:::.dwnorm_fast.ss(x, mean_val, sd_val, omega_multi, crit_x_multi, log = FALSE)
    })
    y_exported <- dwnorm(x_vals, mean_val, sd_val, crit_x = crit_x_multi, omega = as.vector(omega_multi), type = "one.sided")

    plot(x_vals, y_fast, type = "l", main = "Multiple critical values: fast.ss vs exported",
         xlab = "x", ylab = "density", col = "blue", lwd = 2)
    lines(x_vals, y_exported, col = "red", lty = 2, lwd = 2)
    abline(v = crit_x_multi, col = "gray", lty = 3)

    # Add weight labels
    text(x = c(-1, 0.4, 0.9, 1.8), y = max(y_fast) * 0.8,
         labels = paste("ω =", omega_multi[1,]), col = "darkgreen")
    legend("topright", c("fast.ss", "exported"), col = c("blue", "red"), lty = c(1, 2), lwd = 2)
  })
})
