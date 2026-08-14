extern "C" SEXP RoBMA_glmm_binom_marginal_loglik(SEXP ai, SEXP ci,
                                                  SEXP n1i, SEXP n2i,
                                                  SEXP mu_samples,
                                                  SEXP tau_within,
                                                  SEXP weights,
                                                  SEXP theta_grid,
                                                  SEXP log_theta_weights,
                                                  SEXP logit_pi_grid,
                                                  SEXP log_pi_weights)
{
  check_integer_glmm(ai, "ai");
  check_integer_glmm(ci, "ci");
  check_integer_glmm(n1i, "n1i");
  check_integer_glmm(n2i, "n2i");
  check_real_glmm(logit_pi_grid, "logit_pi_grid");
  check_real_glmm(log_pi_weights, "log_pi_weights");

  int S, K, n_theta;
  validate_common_glmm(
    mu_samples, tau_within, theta_grid, log_theta_weights,
    &S, &K, &n_theta
  );
  validate_weights_glmm(weights, K);

  const int n_pi = matrix_nrow_glmm(logit_pi_grid, "logit_pi_grid");
  if (matrix_ncol_glmm(logit_pi_grid, "logit_pi_grid") != K ||
      matrix_nrow_glmm(log_pi_weights, "log_pi_weights") != n_pi ||
      matrix_ncol_glmm(log_pi_weights, "log_pi_weights") != K) {
    Rf_error("'logit_pi_grid' and 'log_pi_weights' must be n_pi x K matrices.");
  }
  if (Rf_length(ai) != K || Rf_length(ci) != K ||
      Rf_length(n1i) != K || Rf_length(n2i) != K) {
    Rf_error("binomial outcome vectors must have length K.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const int *ai_p             = INTEGER(ai);
  const int *ci_p             = INTEGER(ci);
  const int *n1i_p            = INTEGER(n1i);
  const int *n2i_p            = INTEGER(n2i);
  const double *mu_p          = REAL(mu_samples);
  const double *tau_p         = REAL(tau_within);
  const double *theta_p       = REAL(theta_grid);
  const double *log_theta_p   = REAL(log_theta_weights);
  const double *logit_pi_p    = REAL(logit_pi_grid);
  const double *log_pi_w_p    = REAL(log_pi_weights);
  const double *weights_p      = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p               = REAL(out);
  std::vector<double> terms(static_cast<size_t>(n_pi) * static_cast<size_t>(n_theta));
  std::vector<double> grid_weights(static_cast<size_t>(n_pi) * static_cast<size_t>(n_theta));
  std::vector<double> half_effect(static_cast<size_t>(n_theta));

  for (int k = 0; k < K; ++k) {
    const int ai_k  = ai_p[k];
    const int ci_k  = ci_p[k];
    const int n1_k  = n1i_p[k];
    const int n2_k  = n2i_p[k];
    const int bi_k  = n1_k - ai_k;
    const int di_k  = n2_k - ci_k;

    const double log_coef =
      std::lgamma(n1_k + 1.0) - std::lgamma(ai_k + 1.0) - std::lgamma(bi_k + 1.0) +
      std::lgamma(n2_k + 1.0) - std::lgamma(ci_k + 1.0) - std::lgamma(di_k + 1.0);

    int grid_weight_index = 0;
    for (int p = 0; p < n_pi; ++p) {
      const double log_pi_w = log_pi_w_p[p + n_pi * k];
      for (int t = 0; t < n_theta; ++t) {
        grid_weights[grid_weight_index] = log_theta_p[t] + log_pi_w;
        ++grid_weight_index;
      }
    }

    for (int s = 0; s < S; ++s) {
      const double mu_s  = mu_p[s + S * k];
      const double tau_s = tau_p[s + S * k];
      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];
      int grid_index     = 0;

      for (int t = 0; t < n_theta; ++t) {
        half_effect[t] = 0.5 * (mu_s + tau_s * theta_p[t]);
      }

      for (int p = 0; p < n_pi; ++p) {
        const double lpi = logit_pi_p[p + n_pi * k];

        for (int t = 0; t < n_theta; ++t) {
          const double eta_1 = lpi + half_effect[t];
          const double eta_2 = lpi - half_effect[t];
          double log_p_1, log_q_1, log_p_2, log_q_2;
          log_inv_logit_pair_glmm(eta_1, &log_p_1, &log_q_1);
          log_inv_logit_pair_glmm(eta_2, &log_p_2, &log_q_2);

          const double log_lik = glmm_binomial_log_likelihood(
            ai_k, ci_k, n1_k, n2_k, log_coef,
            log_p_1, log_q_1, log_p_2, log_q_2
          );

          terms[grid_index] = weight_k * log_lik + grid_weights[grid_index];
          ++grid_index;
        }
      }

      out_p[s + S * k] = buffered_logsumexp_glmm(terms);
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_glmm_pois_marginal_loglik(SEXP x1i, SEXP x2i,
                                                 SEXP t1i, SEXP t2i,
                                                 SEXP mu_samples,
                                                 SEXP tau_within,
                                                 SEXP weights,
                                                 SEXP theta_grid,
                                                 SEXP log_theta_weights,
                                                 SEXP log_phi_grid,
                                                 SEXP log_phi_weights)
{
  check_integer_glmm(x1i, "x1i");
  check_integer_glmm(x2i, "x2i");
  check_real_glmm(t1i, "t1i");
  check_real_glmm(t2i, "t2i");
  check_real_glmm(log_phi_grid, "log_phi_grid");
  check_real_glmm(log_phi_weights, "log_phi_weights");

  int S, K, n_theta;
  validate_common_glmm(
    mu_samples, tau_within, theta_grid, log_theta_weights,
    &S, &K, &n_theta
  );
  validate_weights_glmm(weights, K);

  const int n_phi = matrix_nrow_glmm(log_phi_grid, "log_phi_grid");
  if (matrix_ncol_glmm(log_phi_grid, "log_phi_grid") != K ||
      matrix_nrow_glmm(log_phi_weights, "log_phi_weights") != n_phi ||
      matrix_ncol_glmm(log_phi_weights, "log_phi_weights") != K) {
    Rf_error("'log_phi_grid' and 'log_phi_weights' must be n_phi x K matrices.");
  }
  if (Rf_length(x1i) != K || Rf_length(x2i) != K ||
      Rf_length(t1i) != K || Rf_length(t2i) != K) {
    Rf_error("Poisson outcome vectors must have length K.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  const int *x1i_p            = INTEGER(x1i);
  const int *x2i_p            = INTEGER(x2i);
  const double *t1i_p         = REAL(t1i);
  const double *t2i_p         = REAL(t2i);
  const double *mu_p          = REAL(mu_samples);
  const double *tau_p         = REAL(tau_within);
  const double *theta_p       = REAL(theta_grid);
  const double *log_theta_p   = REAL(log_theta_weights);
  const double *log_phi_p     = REAL(log_phi_grid);
  const double *log_phi_w_p   = REAL(log_phi_weights);
  const double *weights_p      = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p               = REAL(out);
  std::vector<double> terms(static_cast<size_t>(n_phi) * static_cast<size_t>(n_theta));
  std::vector<double> grid_weights(static_cast<size_t>(n_phi) * static_cast<size_t>(n_theta));
  std::vector<double> half_effect(static_cast<size_t>(n_theta));

  for (int k = 0; k < K; ++k) {
    const int x1_k = x1i_p[k];
    const int x2_k = x2i_p[k];
    const double log_t1 = std::log(t1i_p[k]);
    const double log_t2 = std::log(t2i_p[k]);
    const double log_factorial = std::lgamma(x1_k + 1.0) + std::lgamma(x2_k + 1.0);

    int grid_weight_index = 0;
    for (int p = 0; p < n_phi; ++p) {
      const double log_phi_w = log_phi_w_p[p + n_phi * k];
      for (int t = 0; t < n_theta; ++t) {
        grid_weights[grid_weight_index] = log_theta_p[t] + log_phi_w;
        ++grid_weight_index;
      }
    }

    for (int s = 0; s < S; ++s) {
      const double mu_s  = mu_p[s + S * k];
      const double tau_s = tau_p[s + S * k];
      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];
      int grid_index     = 0;

      for (int t = 0; t < n_theta; ++t) {
        half_effect[t] = 0.5 * (mu_s + tau_s * theta_p[t]);
      }

      for (int p = 0; p < n_phi; ++p) {
        const double base_log_lambda_1 = log_phi_p[p + n_phi * k] + log_t1;
        const double base_log_lambda_2 = log_phi_p[p + n_phi * k] + log_t2;
        for (int t = 0; t < n_theta; ++t) {
          const double log_lambda_1 = base_log_lambda_1 + half_effect[t];
          const double log_lambda_2 = base_log_lambda_2 - half_effect[t];

          const double log_lik =
            poisson_log_kernel_glmm(x1_k, log_lambda_1) +
            poisson_log_kernel_glmm(x2_k, log_lambda_2) -
            log_factorial;

          terms[grid_index] = weight_k * log_lik + grid_weights[grid_index];
          ++grid_index;
        }
      }

      out_p[s + S * k] = buffered_logsumexp_glmm(terms);
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_glmm_binom_marginal_loglik_row_sum(
  SEXP ai, SEXP ci, SEXP n1i, SEXP n2i, SEXP mu_samples,
  SEXP tau_within, SEXP weights, SEXP theta_grid,
  SEXP log_theta_weights, SEXP logit_pi_grid, SEXP log_pi_weights)
{
  check_integer_glmm(ai, "ai");
  check_integer_glmm(ci, "ci");
  check_integer_glmm(n1i, "n1i");
  check_integer_glmm(n2i, "n2i");
  check_real_glmm(logit_pi_grid, "logit_pi_grid");
  check_real_glmm(log_pi_weights, "log_pi_weights");

  int S, K, n_theta;
  validate_common_glmm(
    mu_samples, tau_within, theta_grid, log_theta_weights,
    &S, &K, &n_theta
  );
  validate_weights_glmm(weights, K);

  const int n_pi = matrix_nrow_glmm(logit_pi_grid, "logit_pi_grid");
  if (matrix_ncol_glmm(logit_pi_grid, "logit_pi_grid") != K ||
      matrix_nrow_glmm(log_pi_weights, "log_pi_weights") != n_pi ||
      matrix_ncol_glmm(log_pi_weights, "log_pi_weights") != K) {
    Rf_error("'logit_pi_grid' and 'log_pi_weights' must be n_pi x K matrices.");
  }
  if (Rf_length(ai) != K || Rf_length(ci) != K ||
      Rf_length(n1i) != K || Rf_length(n2i) != K) {
    Rf_error("binomial outcome vectors must have length K.");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, S));

  const int *ai_p           = INTEGER(ai);
  const int *ci_p           = INTEGER(ci);
  const int *n1i_p          = INTEGER(n1i);
  const int *n2i_p          = INTEGER(n2i);
  const double *mu_p        = REAL(mu_samples);
  const double *tau_p       = REAL(tau_within);
  const double *theta_p     = REAL(theta_grid);
  const double *log_theta_p = REAL(log_theta_weights);
  const double *logit_pi_p  = REAL(logit_pi_grid);
  const double *log_pi_w_p  = REAL(log_pi_weights);
  const double *weights_p   = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p             = REAL(out);

  const int n_grid = n_pi * n_theta;
  std::vector<double> terms(static_cast<size_t>(n_grid));
  std::vector<double> half_effect(static_cast<size_t>(n_theta));
  std::vector<double> log_coef(static_cast<size_t>(K));
  std::vector<double> grid_weights(static_cast<size_t>(K) * static_cast<size_t>(n_grid));

  for (int k = 0; k < K; ++k) {
    const int bi_k = n1i_p[k] - ai_p[k];
    const int di_k = n2i_p[k] - ci_p[k];
    log_coef[k] =
      std::lgamma(n1i_p[k] + 1.0) - std::lgamma(ai_p[k] + 1.0) - std::lgamma(bi_k + 1.0) +
      std::lgamma(n2i_p[k] + 1.0) - std::lgamma(ci_p[k] + 1.0) - std::lgamma(di_k + 1.0);

    double *grid_weights_k = grid_weights.data() + static_cast<size_t>(k) * static_cast<size_t>(n_grid);
    int grid_index = 0;
    for (int p = 0; p < n_pi; ++p) {
      const double log_pi_w = log_pi_w_p[p + n_pi * k];
      for (int t = 0; t < n_theta; ++t) {
        grid_weights_k[grid_index] = log_theta_p[t] + log_pi_w;
        ++grid_index;
      }
    }
  }

  for (int s = 0; s < S; ++s) {
    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];
      row_sum += binom_marginal_loglik_scalar(
        ai_p[k], ci_p[k], n1i_p[k], n2i_p[k], log_coef[k],
        mu_p[s + S * k], tau_p[s + S * k], weight_k, theta_p,
        logit_pi_p + n_pi * k,
        grid_weights.data() + static_cast<size_t>(k) * static_cast<size_t>(n_grid),
        n_theta, n_pi, terms, half_effect
      );
    }
    out_p[s] = row_sum;
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_glmm_pois_marginal_loglik_row_sum(
  SEXP x1i, SEXP x2i, SEXP t1i, SEXP t2i, SEXP mu_samples,
  SEXP tau_within, SEXP weights, SEXP theta_grid,
  SEXP log_theta_weights, SEXP log_phi_grid, SEXP log_phi_weights)
{
  check_integer_glmm(x1i, "x1i");
  check_integer_glmm(x2i, "x2i");
  check_real_glmm(t1i, "t1i");
  check_real_glmm(t2i, "t2i");
  check_real_glmm(log_phi_grid, "log_phi_grid");
  check_real_glmm(log_phi_weights, "log_phi_weights");

  int S, K, n_theta;
  validate_common_glmm(
    mu_samples, tau_within, theta_grid, log_theta_weights,
    &S, &K, &n_theta
  );
  validate_weights_glmm(weights, K);

  const int n_phi = matrix_nrow_glmm(log_phi_grid, "log_phi_grid");
  if (matrix_ncol_glmm(log_phi_grid, "log_phi_grid") != K ||
      matrix_nrow_glmm(log_phi_weights, "log_phi_weights") != n_phi ||
      matrix_ncol_glmm(log_phi_weights, "log_phi_weights") != K) {
    Rf_error("'log_phi_grid' and 'log_phi_weights' must be n_phi x K matrices.");
  }
  if (Rf_length(x1i) != K || Rf_length(x2i) != K ||
      Rf_length(t1i) != K || Rf_length(t2i) != K) {
    Rf_error("Poisson outcome vectors must have length K.");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, S));

  const int *x1i_p          = INTEGER(x1i);
  const int *x2i_p          = INTEGER(x2i);
  const double *t1i_p       = REAL(t1i);
  const double *t2i_p       = REAL(t2i);
  const double *mu_p        = REAL(mu_samples);
  const double *tau_p       = REAL(tau_within);
  const double *theta_p     = REAL(theta_grid);
  const double *log_theta_p = REAL(log_theta_weights);
  const double *log_phi_p   = REAL(log_phi_grid);
  const double *log_phi_w_p = REAL(log_phi_weights);
  const double *weights_p   = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p             = REAL(out);

  const int n_grid = n_phi * n_theta;
  std::vector<double> terms(static_cast<size_t>(n_grid));
  std::vector<double> half_effect(static_cast<size_t>(n_theta));
  std::vector<double> log_t1(static_cast<size_t>(K));
  std::vector<double> log_t2(static_cast<size_t>(K));
  std::vector<double> log_factorial(static_cast<size_t>(K));
  std::vector<double> grid_weights(static_cast<size_t>(K) * static_cast<size_t>(n_grid));

  for (int k = 0; k < K; ++k) {
    log_t1[k]        = std::log(t1i_p[k]);
    log_t2[k]        = std::log(t2i_p[k]);
    log_factorial[k] = std::lgamma(x1i_p[k] + 1.0) + std::lgamma(x2i_p[k] + 1.0);

    double *grid_weights_k = grid_weights.data() + static_cast<size_t>(k) * static_cast<size_t>(n_grid);
    int grid_index = 0;
    for (int p = 0; p < n_phi; ++p) {
      const double log_phi_w = log_phi_w_p[p + n_phi * k];
      for (int t = 0; t < n_theta; ++t) {
        grid_weights_k[grid_index] = log_theta_p[t] + log_phi_w;
        ++grid_index;
      }
    }
  }

  for (int s = 0; s < S; ++s) {
    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];
      row_sum += pois_marginal_loglik_scalar(
        x1i_p[k], x2i_p[k], log_t1[k], log_t2[k], log_factorial[k],
        mu_p[s + S * k], tau_p[s + S * k], weight_k, theta_p,
        log_phi_p + n_phi * k,
        grid_weights.data() + static_cast<size_t>(k) * static_cast<size_t>(n_grid),
        n_theta, n_phi, terms, half_effect
      );
    }
    out_p[s] = row_sum;
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_glmm_binom_conditional_loglik_sum(
  SEXP ai, SEXP ci, SEXP n1i, SEXP n2i, SEXP mu_samples,
  SEXP logit_baserate, SEXP weights)
{
  check_integer_glmm(ai, "ai");
  check_integer_glmm(ci, "ci");
  check_integer_glmm(n1i, "n1i");
  check_integer_glmm(n2i, "n2i");

  int S, K;
  validate_conditional_common_glmm(mu_samples, logit_baserate, weights, &S, &K);

  if (Rf_length(ai) != K || Rf_length(ci) != K ||
      Rf_length(n1i) != K || Rf_length(n2i) != K) {
    Rf_error("binomial outcome vectors must have length K.");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, S));

  const int *ai_p          = INTEGER(ai);
  const int *ci_p          = INTEGER(ci);
  const int *n1i_p         = INTEGER(n1i);
  const int *n2i_p         = INTEGER(n2i);
  const double *mu_p       = REAL(mu_samples);
  const double *base_p     = REAL(logit_baserate);
  const double *weights_p  = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p            = REAL(out);
  std::vector<double> log_coef(static_cast<size_t>(K));

  for (int k = 0; k < K; ++k) {
    const int bi_k = n1i_p[k] - ai_p[k];
    const int di_k = n2i_p[k] - ci_p[k];
    log_coef[k] =
      std::lgamma(n1i_p[k] + 1.0) - std::lgamma(ai_p[k] + 1.0) - std::lgamma(bi_k + 1.0) +
      std::lgamma(n2i_p[k] + 1.0) - std::lgamma(ci_p[k] + 1.0) - std::lgamma(di_k + 1.0);
  }

  for (int s = 0; s < S; ++s) {
    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      const int cell = s + S * k;
      const double half_effect = 0.5 * mu_p[cell];
      const double eta_1 = base_p[cell] + half_effect;
      const double eta_2 = base_p[cell] - half_effect;
      double log_p_1, log_q_1, log_p_2, log_q_2;
      log_inv_logit_pair_glmm(eta_1, &log_p_1, &log_q_1);
      log_inv_logit_pair_glmm(eta_2, &log_p_2, &log_q_2);

      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];
      const double log_lik = glmm_binomial_log_likelihood(
        ai_p[k], ci_p[k], n1i_p[k], n2i_p[k], log_coef[k],
        log_p_1, log_q_1, log_p_2, log_q_2
      );
      row_sum += weight_k * log_lik;
    }
    out_p[s] = row_sum;
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_glmm_pois_conditional_loglik_sum(
  SEXP x1i, SEXP x2i, SEXP t1i, SEXP t2i, SEXP mu_samples,
  SEXP log_phi, SEXP weights)
{
  check_integer_glmm(x1i, "x1i");
  check_integer_glmm(x2i, "x2i");
  check_real_glmm(t1i, "t1i");
  check_real_glmm(t2i, "t2i");

  int S, K;
  validate_conditional_common_glmm(mu_samples, log_phi, weights, &S, &K);

  if (Rf_length(x1i) != K || Rf_length(x2i) != K ||
      Rf_length(t1i) != K || Rf_length(t2i) != K) {
    Rf_error("Poisson outcome vectors must have length K.");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, S));

  const int *x1i_p         = INTEGER(x1i);
  const int *x2i_p         = INTEGER(x2i);
  const double *t1i_p      = REAL(t1i);
  const double *t2i_p      = REAL(t2i);
  const double *mu_p       = REAL(mu_samples);
  const double *log_phi_p  = REAL(log_phi);
  const double *weights_p  = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p            = REAL(out);
  std::vector<double> log_t1(static_cast<size_t>(K));
  std::vector<double> log_t2(static_cast<size_t>(K));
  std::vector<double> log_factorial(static_cast<size_t>(K));

  for (int k = 0; k < K; ++k) {
    log_t1[k]        = std::log(t1i_p[k]);
    log_t2[k]        = std::log(t2i_p[k]);
    log_factorial[k] = std::lgamma(x1i_p[k] + 1.0) + std::lgamma(x2i_p[k] + 1.0);
  }

  for (int s = 0; s < S; ++s) {
    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      const int cell = s + S * k;
      const double half_effect = 0.5 * mu_p[cell];
      const double log_lambda_1 = log_phi_p[cell] + log_t1[k] + half_effect;
      const double log_lambda_2 = log_phi_p[cell] + log_t2[k] - half_effect;
      const double lambda_1     = std::exp(log_lambda_1);
      const double lambda_2     = std::exp(log_lambda_2);

      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];
      const double log_lik =
        x1i_p[k] * log_lambda_1 - lambda_1 +
        x2i_p[k] * log_lambda_2 - lambda_2 -
        log_factorial[k];
      row_sum += weight_k * log_lik;
    }
    out_p[s] = row_sum;
  }

  UNPROTECT(1);
  return out;
}
