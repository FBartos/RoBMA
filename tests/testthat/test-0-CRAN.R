context("(0) Basic tests for CRAN")
# These are just a very rudimentary tests that don't require time or saved files.
# The full range of tests is run locally.

# test the prior specification function
p  <- prior("normal", list(0, 1))
p1 <- prior("normal", list(0, 1), list(-Inf, 0))

test_that("print prior function (1)", {
  expect_equal(
    capture.output(p),
    "Normal(0, 1)[-Inf, Inf]"
  )
})

test_that("print prior function ¨(2)", {
  expect_equal(
    capture.output(p1),
    "Normal(0, 1)[-Inf, 0]"
  )
})


# fit a default model
t <- c(1,2,3)
n <- c(20, 15, 10)
fit_default <- RoBMA(t = t, n = n, chains = 2, burnin = 1000, iter = 4000, control = list(silent = TRUE), seed = 666)
test_that("Default model can be fit",
  expect_equal(TRUE, is.RoBMA(fit_default)))

# check the S3 methods
test_that("print function", {
  expect_equal(
    capture.output(fit_default),
    c("Call:"                                                                                      ,
      "RoBMA(t = t, n = n, chains = 2, iter = 4000, burnin = 1000, control = list(silent = TRUE), ",
      "    seed = 666)"                                                                            ,
      ""                                                                                           ,
      "Estimates:"                                                                                 ,
      "             mu             tau   omega[0,0.05] omega[0.05,0.1]    omega[0.1,1] "           ,
      "      0.5999943       0.1475925       1.0000000       0.8019722       0.6671523 " )
  )
})


test_that("summary function", {
  expect_equal(
    capture.output(summary(fit_default)),
    c("Call:"                                                                                      ,
      "RoBMA(t = t, n = n, chains = 2, iter = 4000, burnin = 1000, control = list(silent = TRUE), ",
      "    seed = 666)"                                                                            ,
      ""                                                                                           ,
      "Robust Bayesian Meta-Analysis"                                                              ,
      "              Models Prior prob. Post. prob. Incl. BF"                                      ,
      "Effect          6/12       0.500       0.807    4.184"                                      ,
      "Heterogeneity   6/12       0.500       0.477    0.912"                                      ,
      "Pub. bias       8/12       0.500       0.545    1.197"                                      ,
      ""                                                                                           ,
      "Model-averaged estimates"                                                                   ,
      "                 Mean Median 0.025 0.975"                                                   ,
      "mu              0.600  0.640 0.000 1.383"                                                   ,
      "tau             0.148  0.000 0.000 0.946"                                                   ,
      "omega[0,0.05]   1.000  1.000 1.000 1.000"                                                   ,
      "omega[0.05,0.1] 0.802  0.952 0.171 1.000"                                                   ,
      "omega[0.1,1]    0.667  0.832 0.043 1.000"                                                   ,
      "(Estimated omegas correspond to two-sided p-values)"  )
  )
})
