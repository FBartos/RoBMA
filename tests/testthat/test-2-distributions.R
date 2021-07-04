context("(2) Distribution functions")
skip_on_cran()

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
  expect_doppelganger("dwnorm-1", function()curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), from = -3, to = 3))
  expect_doppelganger("dwnorm-2", function()curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided"), from = -3, to = 3))
  expect_doppelganger("dwnorm-3", function()curve(dwnorm(x, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided"), from = -3, to = 3))
  expect_doppelganger("dwnorm-4", function()curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, .5), type = "two.sided"), from = -3, to = 3))
  expect_doppelganger("dwnorm-5", function()curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided"), from = -3, to = 3))
  expect_doppelganger("dwnorm-6", function()curve(dwnorm(x, 0.5, 1, steps = c(.05, .10), omega = c(1, .5, .1), type = "two.sided"), from = -3, to = 3))

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
    RoBMA:::.dwnorm_fast(x1, mean1, sd1, omega1[1,], matrix(c(0.50, 1.50), nrow = 10, ncol = 2, byrow = TRUE), type = "one.sided", log = TRUE),
    dwnorm(x1, mean1, sd1, crit_x = matrix(c(0.50, 1.50), nrow = 10, ncol = 2, byrow = TRUE), omega = omega1[1,], type = "one.sided", log = TRUE)
  )
  expect_equal(
    RoBMA:::.dwnorm_fast(x1, mean1, sd1, omega1[1,], matrix(c(0.50, 1.50), nrow = 10, ncol = 2, byrow = TRUE), type = "two.sided", log = TRUE),
    dwnorm(x1, mean1, sd1, crit_x = matrix(c(0.50, 1.50), nrow = 10, ncol = 2, byrow = TRUE), omega = omega1[1,], type = "two.sided", log = TRUE)
  )
})

test_that("Random number generator function works", {

  set.seed(1)

  # verify (visually) no-weights against standard normal
  expect_equal(rwnorm(0, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided"), rnorm(0, 0, 1))
  expect_doppelganger("rwnorm-1", function(){
    plot(density(rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, 1), type = "one.sided")))
    curve(dnorm(x, 0, 1), from = -10, to = 10, col = "blue", add = TRUE)
  })
  expect_doppelganger("rwnorm-2", function(){
    plot(density(rwnorm(10000, 1, 2, steps = c(.5), omega = c(1, 1), type = "two.sided")))
    curve(dnorm(x, 1, 2), from = -10, to = 10, col = "blue", add = TRUE)
  })

  # visually verify against pdf
  expect_doppelganger("rwnorm-3", function(){
    hist(rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), from = -10, to = 10, col = "blue", add = TRUE)
  })
  expect_doppelganger("rwnorm-4", function(){
    hist(rwnorm(10000, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(.5, 1), type = "one.sided"), from = -10, to = 10, col = "blue", add = TRUE)
    })
  expect_doppelganger("rwnorm-5", function(){
    hist(rwnorm(10000, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 1, 2, steps = c(.10, .50), omega = c(1, .2, .5), type = "one.sided"), from = -10, to = 10, col = "blue", add = TRUE)
    })
  expect_doppelganger("rwnorm-6", function(){
    hist(rwnorm(10000, 0, 1, steps = c(.5), omega = c(1, .5), type = "one.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(1, .5), type = "two.sided"), from = -10, to = 10, col = "blue", add = TRUE)
    })
  expect_doppelganger("rwnorm-7", function(){
    hist(rwnorm(10000, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided"), freq = FALSE, breaks = 50)
    curve(dwnorm(x, 0, 1, steps = c(.5), omega = c(.5, 1), type = "two.sided"), from = -10, to = 10, col = "blue", add = TRUE)
    })
  expect_doppelganger("rwnorm-8", function(){
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
