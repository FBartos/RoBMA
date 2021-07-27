model_syntax1 <- "
  model{
    omega[1] = .1
    omega[2] ~ dunif(0.20,0.201);
    omega[3] = 1

    y[] ~ dwmnorm_1s(mu[], sigma[,], omega[], crit_y)
  }"

data1 <- list(
  y      = c(1.5, 0, 3),
  mu     = c(0.2, 0.5, -0.1),
  sigma  = matrix(c(
    1.5, 1.0, 0.5,
    1.0, 1.8, 0.7,
    0.5, 0.7, 1.2), nrow = 3, ncol = 3),
  crit_y = matrix(c(
    1.25, 1.96,
    1.30, 2.05,
    1.10, 1.50), nrow = 3, ncol = 2)
)

model_syntax2 <- "
  model{
    omega[1] ~ dunif(0.10,0.101);
    omega[2] = 1

    y[] ~ dwmnorm_1s(mu[], sigma[,], omega[], crit_y)
  }"

data2 <- list(
  y      = c(-1, -1.5, -1.1,- 1.2),
  mu     = c(0.2, 0.5, -0.1, 1.1),
  sigma  = matrix(c(
    1.5, 1.0, 0.5, 0.2,
    1.0, 1.8, 0.7, 0.6,
    0.5, 0.7, 1.2, 0.8,
    0.2, 0.6, 0.8, 3.6), nrow = 4, ncol = 4),
  crit_y = matrix(c(
    1.25,
    1.30,
    1.05,
    2.05), nrow = 4, ncol = 1)
)

model1 <- rjags::jags.model(file = textConnection(model_syntax1), data = data1, quiet = TRUE, n.adapt=1)
fit1   <- rjags::jags.samples(model = model1, variable.names = "omega", n.iter = 2, quiet = TRUE, progress.bar = "none")


model2 <- rjags::jags.model(file = textConnection(model_syntax2), data = data2, quiet = TRUE, n.adapt=1)
fit2   <- rjags::jags.samples(model = model2, variable.names = "omega", n.iter = 1, quiet = TRUE, progress.bar = "none")

corr[upper.tri(corr)]

K <- 5
m <- matrix(0:(K*K-1), ncol = K, nrow = K, byrow = T)
m
m[lower.tri(m)]

for(i in 0:(K-1))
  print(K * i + i)


TT <- function(n) n * (n + 1) / 2

K = 4
for(i in 0:(K-1)){
  for(j in 0:(K-1)){

    print(sigma[K * i + j + 1])
    print(sigma[K * j + i + 1])
    print("xxx")
  }
}

