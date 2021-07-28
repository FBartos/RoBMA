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
fit1   <- rjags::jags.samples(model = model1, variable.names = "omega", n.iter = 1, quiet = TRUE, progress.bar = "none")


model2 <- rjags::jags.model(file = textConnection(model_syntax2), data = data2, quiet = TRUE, n.adapt=2)
fit2   <- rjags::jags.samples(model = model2, variable.names = "omega", n.iter = 2, quiet = TRUE, progress.bar = "none")


y      = c(1.5, 0, 3)
mu     = c(0.2, 0.5, -0.1)
sigma  = matrix(c(
  1.5, 1.0, 0.5,
  1.0, 1.8, 0.7,
  0.5, 0.7, 1.2), nrow = 3, ncol = 3)

solve(chol_decomp)
solve(t(chol_decomp))

chol_decomp <- chol(sigma)
chol_decomp
rooti       <- solve(chol_decomp)
rootisum    <- sum(log(diag(rooti)))
constants   <- -nrow(sigma)/2.0 * log(2 * pi)
other_terms <- rootisum + constants

z          <- (y - mu)
z          <- inplace_tri_mat_mult(z, rooti)
other_terms - 0.5 * sum(z * z)

mvtnorm::dmvnorm(
  x     = y,
  mean  = mu,
  sigma = sigma, log = T)

inplace_tri_mat_mult <-  function(x, mat){
  n = ncol(mat)

  for(j in n:1){
    tmp <- 0
    for(i in 1:j)
      tmp = tmp + mat[i, j] * x[i]
      x[j] = tmp
  }

  return(x)
}
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
      x[j] = tmp;
  }
}

arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
double const rootisum = arma::sum(log(rooti.diag())),
constants = -(double)xdim/2.0 * log2pi,
other_terms = rootisum + constants;

arma::rowvec z;
for (uword i = 0; i < n; i++) {
  z = (x.row(i) - mean);
  inplace_tri_mat_mult(z, rooti);
  out(i) = other_terms - 0.5 * arma::dot(z, z);
}
