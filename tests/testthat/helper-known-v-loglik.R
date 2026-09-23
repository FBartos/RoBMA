.test_log_lik_known_v_joint_sum_from_evaluated_predictors <- function(
    fit, data, priors, mu_samples, tau_within_samples,
    tau_between_samples = NULL, posterior_samples = NULL, unit = "estimate",
    data_hash = NULL, random_effects_conditioning = "none") {

  unit  <- .normalize_unit(unit)
  setup <- .log_lik_evaluated_setup(
    fit                         = fit,
    data                        = data,
    priors                      = priors,
    unit                        = unit,
    data_hash                   = data_hash,
    mu_samples                  = mu_samples,
    tau_within_samples          = tau_within_samples,
    tau_between_samples         = tau_between_samples,
    posterior_samples           = posterior_samples,
    random_effects_conditioning = random_effects_conditioning
  )

  return(.log_lik_known_v_joint_sum_from_setup(setup))
}

.known_v_conditional_loglik_reference <- function(yi, mu, covariance) {

  precision          <- solve(covariance)
  precision_diagonal <- diag(precision)
  variance           <- 1 / precision_diagonal
  residual           <- as.vector(precision %*% (yi - mu)) /
    precision_diagonal

  -0.5 * (log(2 * pi * variance) + residual^2 / variance)
}

.known_v_gls_projection_reference <- function(X, y, covariance) {

  covariance_factor <- t(chol(covariance))
  W                  <- solve(covariance)
  WX                 <- W %*% X
  XtWX_inv           <- solve(crossprod(X, WX))
  H                  <- X %*% XtWX_inv %*% t(WX)
  beta_hat           <- as.vector(XtWX_inv %*% crossprod(X, W %*% y))

  list(
    covariance        = covariance,
    covariance_factor = covariance_factor,
    W                 = W,
    WX                = WX,
    XtWX_inv          = XtWX_inv,
    H                 = H,
    beta_hat          = beta_hat,
    residual          = y - as.vector(X %*% beta_hat)
  )
}
