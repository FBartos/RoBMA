# the individual variance components
sigma_study   <- 0.5
sigma_effect  <- 0.3
sigma_error   <- c(0.25, 0.15, 0.08)

# study heterogeneity matrix
E_study <- diag(3)
E_study[1:2, 1:2] <- sigma_study^2
E_study[3,   3]   <- sigma_study^2

e_study <- mvtnorm::rmvnorm(10000, rep(0, 3), E_study)
apply(e_study, 2, sd)
cor(e_study)

# effect heterogeneity matrix
E_effect <- diag(sigma_effect^2, 3)

e_effect <- mvtnorm::rmvnorm(10000, rep(0, 3), E_effect)
apply(e_effect, 2, sd)
cor(e_effect)

# residual error
E_error <- diag(sigma_error^2, 3)

e_error <- mvtnorm::rmvnorm(10000, rep(0, 3), E_error)
apply(e_error, 2, sd)
cor(e_error)

# or sampled together:
E_all <- E_study + E_effect + E_error

e_all <- mvtnorm::rmvnorm(10000, rep(0, 3), E_all)

apply(e_all, 2, sd)
cor(e_all)

apply(e_study + e_effect + e_error, 2, sd)
cor(e_study + e_effect + e_error)


lower[isNInf(lower)] <- 0
upper[isInf(upper)] <- 0
error <- 0
value <- 0
inform <- 0
.C(C_mvtdst, N = as.integer(n), NU = as.integer(df), LOWER = as.double(lower),
   UPPER = as.double(upper), INFIN = as.integer(infin),
   CORREL = as.double(corrF), DELTA = as.double(delta),
   MAXPTS = as.integer(x$maxpts), ABSEPS = as.double(x$abseps),
   RELEPS = as.double(x$releps), error = as.double(error),
   value = as.double(value), inform = as.integer(inform),
   RND = as.integer(1))

library(mvtnorm)
mvtnorm::pmvnorm(
  lower = rep(-Inf, 2),
  upper = rep(1.96, 2),
  mean  = c(0, 0),,
  sigma = E_all[1:2, 1:2])
