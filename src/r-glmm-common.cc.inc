// Native GLMM evaluators reproduce the R likelihood contracts in log space.
// theta is the standard-normal heterogeneity node; logit_pi/log_phi are nuisance
// integration grids. Observation weights form L(y | theta, nuisance)^weight,
// while quadrature weights remain unscaled additive log weights. Boundary
// probabilities and rates retain their mathematical limits without clamping.

static int matrix_nrow_glmm(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[0];
}

static int matrix_ncol_glmm(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[1];
}

static void check_real_glmm(SEXP x, const char *name)
{
  if (TYPEOF(x) != REALSXP) {
    Rf_error("'%s' must be numeric.", name);
  }
}

static void check_integer_glmm(SEXP x, const char *name)
{
  if (TYPEOF(x) != INTSXP) {
    Rf_error("'%s' must be integer.", name);
  }
}

static double buffered_logsumexp_glmm(const std::vector<double> &terms)
{
  double term_max = -INFINITY;
  for (std::vector<double>::const_iterator it = terms.begin(); it != terms.end(); ++it) {
    if (std::isnan(*it)) {
      return NA_REAL;
    }
    if (*it > term_max) {
      term_max = *it;
    }
  }

  if (std::isinf(term_max)) {
    return term_max;
  }

  double term_sum = 0;
  for (std::vector<double>::const_iterator it = terms.begin(); it != terms.end(); ++it) {
    term_sum += std::exp(*it - term_max);
  }

  return term_max + std::log(term_sum);
}

static void log_inv_logit_pair_glmm(double eta, double *log_p, double *log_q)
{
  if (eta >= 0) {
    const double log_denom = std::log1p(std::exp(-eta));
    *log_p = -log_denom;
    *log_q = -eta - log_denom;
  } else {
    const double log_denom = std::log1p(std::exp(eta));
    *log_p = eta - log_denom;
    *log_q = -log_denom;
  }
}

static double poisson_log_kernel_glmm(int x, double log_lambda)
{
  if (std::isnan(log_lambda)) {
    return NA_REAL;
  }
  if (log_lambda == INFINITY) {
    return -INFINITY;
  }
  if (log_lambda == -INFINITY) {
    return x == 0 ? 0.0 : -INFINITY;
  }

  const double lambda = std::exp(log_lambda);
  if (std::isinf(lambda)) {
    return -INFINITY;
  }

  return (x == 0 ? 0.0 : x * log_lambda) - lambda;
}

static void validate_common_glmm(SEXP mu_samples, SEXP tau_within,
                                 SEXP theta_grid, SEXP log_theta_weights,
                                 int *S, int *K, int *n_theta)
{
  check_real_glmm(mu_samples, "mu_samples");
  check_real_glmm(tau_within, "tau_within");
  check_real_glmm(theta_grid, "theta_grid");
  check_real_glmm(log_theta_weights, "log_theta_weights");

  *S       = matrix_nrow_glmm(mu_samples, "mu_samples");
  *K       = matrix_ncol_glmm(mu_samples, "mu_samples");
  *n_theta = Rf_length(theta_grid);

  if (matrix_nrow_glmm(tau_within, "tau_within") != *S ||
      matrix_ncol_glmm(tau_within, "tau_within") != *K) {
    Rf_error("'tau_within' dimensions must match 'mu_samples'.");
  }
  if (Rf_length(log_theta_weights) != *n_theta) {
    Rf_error("'theta_grid' and 'log_theta_weights' lengths must match.");
  }
}

static void validate_weights_glmm(SEXP weights, int K)
{
  if (weights == R_NilValue) {
    return;
  }

  check_real_glmm(weights, "weights");
  if (Rf_length(weights) != K) {
    Rf_error("'weights' must have one value per observation.");
  }

  const double *weights_p = REAL(weights);
  for (int k = 0; k < K; ++k) {
    if (!std::isfinite(weights_p[k]) || weights_p[k] <= 0.0) {
      Rf_error("'weights' must contain finite positive values.");
    }
  }
}

static void validate_conditional_common_glmm(SEXP mu_samples, SEXP nuisance,
                                             SEXP weights, int *S, int *K)
{
  check_real_glmm(mu_samples, "mu_samples");
  check_real_glmm(nuisance, "nuisance");

  *S = matrix_nrow_glmm(mu_samples, "mu_samples");
  *K = matrix_ncol_glmm(mu_samples, "mu_samples");

  if (matrix_nrow_glmm(nuisance, "nuisance") != *S ||
      matrix_ncol_glmm(nuisance, "nuisance") != *K) {
    Rf_error("'nuisance' dimensions must match 'mu_samples'.");
  }

  validate_weights_glmm(weights, *K);
}

static double binom_marginal_loglik_scalar(
  int ai, int ci, int n1i, int n2i, double log_coef,
  double mu, double tau, double weight, const double *theta, const double *logit_pi,
  const double *grid_weights, int n_theta, int n_pi,
  std::vector<double> &terms, std::vector<double> &half_effect)
{
  int grid_index = 0;

  for (int t = 0; t < n_theta; ++t) {
    half_effect[t] = 0.5 * (mu + tau * theta[t]);
  }

  for (int p = 0; p < n_pi; ++p) {
    const double lpi = logit_pi[p];

    for (int t = 0; t < n_theta; ++t) {
      const double eta_1 = lpi + half_effect[t];
      const double eta_2 = lpi - half_effect[t];
      double log_p_1, log_q_1, log_p_2, log_q_2;
      log_inv_logit_pair_glmm(eta_1, &log_p_1, &log_q_1);
      log_inv_logit_pair_glmm(eta_2, &log_p_2, &log_q_2);

      const double log_lik = glmm_binomial_log_likelihood(
        ai, ci, n1i, n2i, log_coef,
        log_p_1, log_q_1, log_p_2, log_q_2
      );

      terms[grid_index] = weight * log_lik + grid_weights[grid_index];
      ++grid_index;
    }
  }

  return buffered_logsumexp_glmm(terms);
}

static double pois_marginal_loglik_scalar(
  int x1i, int x2i, double log_t1i, double log_t2i,
  double log_factorial, double mu, double tau, double weight, const double *theta,
  const double *log_phi, const double *grid_weights, int n_theta, int n_phi,
  std::vector<double> &terms, std::vector<double> &half_effect)
{
  int grid_index = 0;

  for (int t = 0; t < n_theta; ++t) {
    half_effect[t] = 0.5 * (mu + tau * theta[t]);
  }

  for (int p = 0; p < n_phi; ++p) {
    const double base_log_lambda_1 = log_phi[p] + log_t1i;
    const double base_log_lambda_2 = log_phi[p] + log_t2i;
    for (int t = 0; t < n_theta; ++t) {
      const double log_lambda_1 = base_log_lambda_1 + half_effect[t];
      const double log_lambda_2 = base_log_lambda_2 - half_effect[t];

      const double log_lik =
        poisson_log_kernel_glmm(x1i, log_lambda_1) +
        poisson_log_kernel_glmm(x2i, log_lambda_2) -
        log_factorial;

      terms[grid_index] = weight * log_lik + grid_weights[grid_index];
      ++grid_index;
    }
  }

  return buffered_logsumexp_glmm(terms);
}
