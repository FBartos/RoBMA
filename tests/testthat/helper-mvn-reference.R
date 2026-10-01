# Independent dense multivariate-normal oracle for covariance-kernel tests.
.marglik_mvn_log_density <- function(y, mean, covariance) {

  as.numeric(mvtnorm::dmvnorm(
    x     = y,
    mean  = mean,
    sigma = covariance,
    log   = TRUE
  ))
}
