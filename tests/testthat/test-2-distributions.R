context("(2) distribution functions")
skip_on_cran()


x <- seq(-5, 5, .5)
x_d0 <- dwt(x, 15, -1, steps = c(.5), omega = c(1, 1))
x_d1 <- dwt(x, 15, 0, steps = c(.5, .30), omega = c(.15, .50, 1), type = "one.sided")
x_d2 <- dwt(x, 15, 0, steps = c(.5, .30), omega = c(.15, .50, 1), type = "two.sided")

x_d1a   <- dwt(x, 15, 0, steps = c(.5, .30), omega = c(.15, .50, 1), type = "one.sided")
x_d1b   <- dwt(x, 5,  1, steps = c(.5, .30), omega = c(.15, .50, 1), type = "one.sided")
x_d1c   <- dwt(x, 15, 2, steps = c(.5, .30), omega = c(.15, .50, 1), type = "one.sided")
x_d1abc <- dwt(rep(x, 3), c(rep(15, length(x)), rep(5, length(x)), rep(15, length(x))),
               c(rep(0, length(x)), rep(1, length(x)), rep(2, length(x))), steps = c(.5, .30), omega = c(.15, .50, 1), type = "one.sided")

test_that("Density functions works", {

  # plot(x, x_d1, "l")
  # plot(x, x_d2, "l")

  expect_equal(x_d0, dt(x, 15, -1))
  expect_equal(round(x_d1,10), c(0.0000484535, 0.0001332085, 0.0003723157, 0.0010444143, 0.0028847949, 0.0076371866, 0.0186971784, 0.0405034563, 0.0739341388, 0.1085540748, 0.4130033293, 0.3618469161, 0.4928942587, 0.2700230419, 0.1246478561, 0.0509145777, 0.0192319662, 0.0069627623, 0.0024821045, 0.0008880566, 0.0003230232))
 expect_equal(round(x_d2, 10), c(0.0003230232, 0.0008880566, 0.0024821045, 0.0069627623, 0.0192319662, 0.0509145777, 0.1246478561, 0.2700230419, 0.2464471294, 0.1085540748, 0.1239009988, 0.1085540748, 0.2464471294, 0.2700230419, 0.1246478561, 0.0509145777, 0.0192319662, 0.0069627623, 0.0024821045, 0.0008880566, 0.0003230232))
  expect_equal(x_d1abc, c(x_d1a, x_d1b, x_d1c))

})


x_p0 <- pwt(x, 15, -1, steps = c(.5), omega = c(1, 1))
x_p1 <- pwt(x, 15, 0, steps = c(.5, .30), omega = c(.15, .50, 1), type = "one.sided")
x_p2 <- pwt(x, 15, 0, steps = c(.5, .30), omega = c(.15, .50, 1), type = "two.sided")

test_that("Probability density functions works", {

  expect_equal(x_p0, pt(x, 15, -1))
  expect_equal(round(x_p1, 10), c(0.0000250057, 0.0000668365, 0.0001830500, 0.0005089786, 0.0014167480, 0.0038693374, 0.0100965801, 0.0243736832, 0.0526058109, 0.0985784390, 0.1578947368, 0.3556157297, 0.6492945938, 0.8375087785, 0.9326894660, 0.9742044176, 0.9905550132, 0.9966068095, 0.9987796665, 0.9995544234, 0.9998332952))
  expect_equal(round(x_p2, 10), c(0.0001667048, 0.0004455766, 0.0012203335, 0.0033931905, 0.0094449868, 0.0257955824, 0.0673105340, 0.1624912215, 0.3332474400, 0.4406837021, 0.5000000000, 0.5593162979, 0.6667525600, 0.8375087785, 0.9326894660, 0.9742044176, 0.9905550132, 0.9966068095, 0.9987796665, 0.9995544234, 0.9998332952))

})


# pretty imprecise
x_q0 <- qwt(x_p0, 15, -1, steps = c(.5), omega = c(1, 1))
x_q1 <- qwt(x_p1, 15, 0, steps = c(.5, .30), omega = c(.15, .50, 1), type = "one.sided")
x_q2 <- qwt(x_p2, 15, 0, steps = c(.5, .30), omega = c(.15, .50, 1), type = "two.sided")

test_that("Quantile functions works", {

  expect_equal(round(x_q0, 5), round(qt(x_p0, 15, -1), 5))
  expect_equal(round(x_q1, 5), round(x, 5))
  expect_equal(round(x_q2, 5), round(x, 5))

})


set.seed(1)
r0 <- rwt(100000, 15, -1, steps = c(.5), omega = c(1, 1))
r1 <- rwt(100000, 15, 0,  steps = c(.5, .30), omega = c(.15, .50, 1), type = "one.sided")
r2 <- rwt(100000, 15, 0,  steps = c(.5, .30), omega = c(.15, .50, 1), type = "two.sided")

test_that("Random number generator functions works", {

  expect_doppelganger("rwt0", ggplot2::ggplot(data = data.frame(x = r0)) + ggplot2::geom_density(ggplot2::aes(x)))
  expect_doppelganger("rwt_one.sided", ggplot2::ggplot(data = data.frame(x = r1)) + ggplot2::geom_density(ggplot2::aes(x)))
  expect_doppelganger("rwt_two.sided", ggplot2::ggplot(data = data.frame(x = r2)) + ggplot2::geom_density(ggplot2::aes(x)))

})
