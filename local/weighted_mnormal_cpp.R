model_syntax1 <- "
  model{
    omega[1] ~ dunif(0.20, 1);
    omega[2] ~ dunif(0.20, 1);
    omega[3] = 1

    y[] ~ dwmnorm_2s(mu[], sigma[,], crit_y, omega[])
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
    1.10, 1.50), nrow = 2, ncol = 3)
)

model_syntax2 <- "
  model{
    omega[1] ~ dunif(0.10,0.101);
    omega[2] = 1

    y[] ~ dwmnorm_2s(mu[], sigma[,], crit_y, omega[])
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
    2.05), nrow = 1, ncol = 4)
)

model1 <- rjags::jags.model(file = textConnection(model_syntax1), data = data1, quiet = TRUE, n.adapt=1)
data1$crit_y

fit1   <- rjags::jags.samples(model = model1, variable.names = "omega", n.iter = 50)
, quiet = TRUE, progress.bar = "none"
mvtnorm::dmvnorm(
  x      = c(1.5, 0, 3),
  mean   = c(0.2, 0.5, -0.1),
  sigma  = matrix(c(
    1.5, 1.0, 0.5,
    1.0, 1.8, 0.7,
    0.5, 0.7, 1.2), nrow = 3, ncol = 3),
  log = T)


model2 <- rjags::jags.model(file = textConnection(model_syntax2), data = data2, quiet = TRUE, n.adapt=10)
fit2   <- rjags::jags.samples(model = model2, variable.names = "omega", n.iter = 2, quiet = TRUE, progress.bar = "none")

mvtnorm::dmvnorm(
  x     = c(-1, -1.5, -1.1,- 1.2),
  mean  = c(0.2, 0.5, -0.1, 1.1),
  sigma = matrix(c(
    1.5, 1.0, 0.5, 0.2,
    1.0, 1.8, 0.7, 0.6,
    0.5, 0.7, 1.2, 0.8,
    0.2, 0.6, 0.8, 3.6), nrow = 4, ncol = 4),
  log = T)

y      = c(1.5, 0, 3)
mu     = c(0.2, 0.5, -0.1)
sigma  = matrix(c(
  1.5, 1.0, 0.5,
  1.0, 1.8, 0.7,
  0.5, 0.7, 1.2), nrow = 3, ncol = 3)


chol_decomp <- chol(sigma)
chol_decomp
rooti       <- solve(chol_decomp)
rootisum    <- sum(log(diag(rooti)))
constants   <- -nrow(sigma)/2.0 * log(2 * pi)
other_terms <- rootisum + constants

z          <- (y - mu)
z          <- inplace_tri_mat_mult(z, rooti)
 - 0.5 * sum(z * z) + other_terms

mvtnorm::dmvnorm(
  x     = y,
  mean  = mu,
  sigma = sigma, log = T)

inplace_tri_mat_mult <-  function(x, mat){
  n = ncol(mat)

  for(j in n:1){
    tmp <- 0
    for(i in 1:j){
      tmp = tmp + mat[i, j] * x[i]
      #print(paste0("[", i, ",", j, "] = ", mat[i, j]))
      #print(paste0("z = ", x[i]))
    }
    #print(paste0("z[", j, "] = ", tmp))

    x[j] = tmp
  }

  return(x)
}

#           [,1]       [,2]       [,3]
# [1,] 0.8164966 -0.6262243 -0.1230100
# [2,] 0.0000000  0.9393364 -0.3382774
# [3,] 0.0000000  0.0000000  1.0455848


# chol_inv[i,]:  0.81649
# chol_inv[i,]: -0.62622 	0.939336
# chol_inv[i,]: -0.12301 -0.338277	1.04558
