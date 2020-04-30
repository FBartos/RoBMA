#Ok so this test is pretty stupid but should at least tell me it can (a) load the jags-module and (b) get the same output as before
#Where is the module?
hereIsTheModule <- sub("/DESCRIPTION", '', sub("/Meta.*", '', attr(packageDescription("RoBMA"), "file")))

rjags::load.module("RoBMA", path = paste0(hereIsTheModule, "/libs", Sys.getenv("R_ARCH")) )

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

model <- rjags::jags.model(file = textConnection(model_syntax), data = data)
fit   <- rjags::jags.samples(model = model, variable.names = "omega", n.iter = 100)

print("summary(fit):")
print(summary(fit))

#Check the summary because I do not know anything about how this works
test_that("Module loaded and summary of fit is same as before", {
  expect_equal(summary(fit)[1], "200")
  expect_equal(summary(fit)[2], "mcarray")
  expect_equal(summary(fit)[3], "numeric")
})