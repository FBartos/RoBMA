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
  lower = c(-Inf, .5),
  upper = c(1, Inf),
  mean  = c(1.2, 1.5),
  sigma = E_all[1:2, 1:2])

### pre-processing
lower <- (lower - mean)/sqrt(diag(sigma))
upper <- (upper - mean)/sqrt(diag(sigma))
delta <- rep(0, length(lower))
corr  <- cov2cor(sigma)
RET  <- mvt(lower = lower, upper = upper, df = 0,
            corr  = corr, delta = mean, algorithm = algorithm,
            ...)

### return 0 on the same upper and lower bounds
any(abs(lower - upper) < sqrt(.Machine$double.eps) * (abs(lower) + abs(upper)) | lower == upper)
return(0)

### settings of inif argument:
# =  1  if the upper bound is Inf
# =  0  if the lower bound is -Inf
# = -1  if the lower bound is -Inf and the upper bound is InF
infin[isInf(upper)]  <- 1
infin[isNInf(lower)] <- 0
infin[isNInf(lower) & isInf(upper)] <- -1
debug(pmvnorm)

### return 1 when all range is unrestricted
all(infin < 0)

### prepare the correlation matrix input
corr <- matrix(as.vector(corr), ncol = n, byrow = TRUE)
corr <- corr[upper.tri(corr)]

ret <- probval(algorithm, n, df, lower, upper, infin, corr,
               delta)


### inside of probval
if (isInf(df))
   df <- 0
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
