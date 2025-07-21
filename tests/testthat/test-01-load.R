context("(1) JAGS module functionality")
skip_on_cran()

### one-sided weight function
test_that("Module can be loaded and the one-sided normal distribution works", {

  module_location <- gsub('/$','', file.path(.libPaths(), "RoBMA", 'libs', if(.Platform$r_arch!="") .Platform$r_arch else ""))
  sapply(module_location, function(path) rjags::load.module("RoBMA", path = path))

  model_syntax <- "
  model{
    eta[1] ~ dgamma(1, 1)
    eta[2] ~ dgamma(1, 1)
    for(j in 1:2){
      std_eta[j]  = eta[j] / sum(eta)
      omega[j]    = sum(std_eta[1:j])
    }
    for(i in 1:K){
      y[i] ~ dwnorm_1s(0, 1, crit_y[i,], omega)
    }
  }"

  data <- list(
    y      = rnorm(10),
    K      = 10,
    crit_y = matrix(1, ncol = 1, nrow = 10)
  )

  RoBMA:::.load_RoBMA_module()
  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
  fit   <- rjags::jags.samples(model = model, variable.names = "omega", n.iter = 100, quiet = TRUE, progress.bar = "none")



  expect_equal(summary(fit)[1], "200")
  expect_equal(summary(fit)[2], "mcarray")
  expect_equal(summary(fit)[3], "numeric")
  expect_equal(unname(dim(fit$omega)), c(2, 100, 1))
  expect_equal(fit$omega[2], 1)
  expect_equal(fit$omega[1] > 0 & fit$omega[1] < 1, TRUE)
})

### two-sided weight function
test_that("Module can be loaded and the two-sided weighted-t distribution works", {

  module_location <- gsub('/$','', file.path(.libPaths(), "RoBMA", 'libs', if(.Platform$r_arch!="") .Platform$r_arch else ""))
  sapply(module_location, function(path) rjags::load.module("RoBMA", path = path))

  model_syntax <- "
  model{
    eta[1] ~ dgamma(1, 1)
    eta[2] ~ dgamma(1, 1)
    eta[3] ~ dgamma(1, 1)
    for(j in 1:3){
      std_eta[j]  = eta[j] / sum(eta)
      omega[j]    = sum(std_eta[1:j])
    }
    for(i in 1:K){
      y[i] ~ dwnorm_2s(0, 1, crit_y[i,], omega)
    }
  }"

  data <- list(
    y      = rnorm(10),
    K      = 10,
    crit_y = matrix(c(1, 2), ncol = 2, nrow = 10, byrow = TRUE)
  )

  model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
  fit   <- rjags::jags.samples(model = model, variable.names = "omega", n.iter = 100, quiet = TRUE, progress.bar = "none")

  expect_equal(summary(fit)[1], "300")
  expect_equal(summary(fit)[2], "mcarray")
  expect_equal(summary(fit)[3], "numeric")
  expect_equal(unname(dim(fit$omega)), c(3, 100, 1))
  expect_equal(fit$omega[3], 1)
  expect_equal(fit$omega[1] < fit$omega[2],         TRUE)
  expect_equal(fit$omega[1] > 0 & fit$omega[1] < 1, TRUE)
  expect_equal(fit$omega[2] > 0 & fit$omega[2] < 1, TRUE)
})
