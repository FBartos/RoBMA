context("(1) JAGS module functionality")

hereIsTheModule <- sub("/DESCRIPTION", '', sub("/Meta.*", '', attr(packageDescription("RoBMA"), "file")))
rjags::load.module("RoBMA", path = paste0(hereIsTheModule, "/libs", Sys.getenv("R_ARCH")) )

### one-sided weight function
model_syntax <- '
model{
for(j in 1:2){
  eta[j] ~ dgamma(1, 1)
}
for(j in 1:2){
  std_eta[j]  = eta[j] / sum(eta)
  omega[j]    = sum(std_eta[1:j])
}
for(i in 1:10){
  t[i] ~ dwt_1s(10, 0, crit_t[i,], omega)
}
}'

data <- list(
  t      = rt(10, 10),
  crit_t = matrix(1.96, ncol = 1, nrow = 10)
)

model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
fit   <- rjags::jags.samples(model = model, variable.names = "omega", n.iter = 100, quiet = TRUE, progress.bar = "none")


test_that("Module can be loaded and the one-sided weighted-t distribution works", {
  expect_equal(summary(fit)[1], "200")
  expect_equal(summary(fit)[2], "mcarray")
  expect_equal(summary(fit)[3], "numeric")
  expect_equal(unname(dim(fit$omega)), c(2, 100, 1))
  expect_equal(fit$omega[2], 1)
  expect_equal(fit$omega[1] > 0 & fit$omega[1] < 1, TRUE)
})



### two-sided weight function
model_syntax <- '
model{
for(j in 1:3){
  eta[j] ~ dgamma(1, 1)
}
for(j in 1:3){
  std_eta[j]  = eta[j] / sum(eta)
  omega[j]    = sum(std_eta[1:j])
}
for(i in 1:10){
  t[i] ~ dwt_2s(10, 0, crit_t[i,], omega)
}
}'

data <- list(
  t      = rt(10, 10),
  crit_t = matrix(c(1.96, 3), ncol = 2, nrow = 10, byrow = TRUE)
)

model <- rjags::jags.model(file = textConnection(model_syntax), data = data, quiet = TRUE)
fit   <- rjags::jags.samples(model = model, variable.names = "omega", n.iter = 100, quiet = TRUE, progress.bar = "none")


test_that("Module can be loaded and the two-sided weighted-t distribution works", {
  expect_equal(summary(fit)[1], "300")
  expect_equal(summary(fit)[2], "mcarray")
  expect_equal(summary(fit)[3], "numeric")
  expect_equal(unname(dim(fit$omega)), c(3, 100, 1))
  expect_equal(fit$omega[3], 1)
  expect_equal(fit$omega[1] < fit$omega[2],         TRUE)
  expect_equal(fit$omega[1] > 0 & fit$omega[1] < 1, TRUE)
  expect_equal(fit$omega[2] > 0 & fit$omega[2] < 1, TRUE)
})
