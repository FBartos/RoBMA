extern "C" SEXP RoBMA_glmm_binom_cluster_loglik(SEXP ai, SEXP ci,
                                                 SEXP n1i, SEXP n2i,
                                                 SEXP mu_samples,
                                                 SEXP tau_within,
                                                 SEXP tau_between,
                                                 SEXP cluster_index,
                                                 SEXP cluster_size,
                                                 SEXP weights,
                                                 SEXP theta_grid,
                                                 SEXP log_theta_weights,
                                                 SEXP logit_pi_grid,
                                                 SEXP log_pi_weights,
                                                 SEXP gamma_grid,
                                                 SEXP log_gamma_weights)
{
  check_integer_glmm(ai, "ai");
  check_integer_glmm(ci, "ci");
  check_integer_glmm(n1i, "n1i");
  check_integer_glmm(n2i, "n2i");
  check_real_glmm(logit_pi_grid, "logit_pi_grid");
  check_real_glmm(log_pi_weights, "log_pi_weights");

  int S, K, n_theta, n_gamma, G, n_cluster_index;
  validate_cluster_common_glmm(
    mu_samples, tau_within, tau_between, theta_grid, log_theta_weights,
    gamma_grid, log_gamma_weights, cluster_index, cluster_size, weights,
    &S, &K, &n_theta, &n_gamma, &G, &n_cluster_index
  );
  (void)n_cluster_index;

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

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, G));

  const int *ai_p              = INTEGER(ai);
  const int *ci_p              = INTEGER(ci);
  const int *n1i_p             = INTEGER(n1i);
  const int *n2i_p             = INTEGER(n2i);
  const int *cluster_index_p   = INTEGER(cluster_index);
  const int *cluster_size_p    = INTEGER(cluster_size);
  const double *mu_p           = REAL(mu_samples);
  const double *tau_within_p   = REAL(tau_within);
  const double *tau_between_p  = REAL(tau_between);
  const double *theta_p        = REAL(theta_grid);
  const double *log_theta_p    = REAL(log_theta_weights);
  const double *logit_pi_p     = REAL(logit_pi_grid);
  const double *log_pi_w_p     = REAL(log_pi_weights);
  const double *gamma_p        = REAL(gamma_grid);
  const double *log_gamma_p    = REAL(log_gamma_weights);
  const double *weights_p      = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p                = REAL(out);

  const int n_grid = n_pi * n_theta;
  std::vector<double> terms(static_cast<size_t>(n_grid));
  std::vector<double> gamma_terms(static_cast<size_t>(n_gamma));
  std::vector<double> cluster_terms(static_cast<size_t>(S) * static_cast<size_t>(n_gamma));
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

  int cluster_offset = 0;
  for (int g = 0; g < G; ++g) {
    const int cluster_n = cluster_size_p[g];
    std::fill(cluster_terms.begin(), cluster_terms.end(), 0.0);

    for (int m = 0; m < cluster_n; ++m) {
      const int k = cluster_index_p[cluster_offset + m] - 1;
      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];

      for (int j = 0; j < n_gamma; ++j) {
        const double gamma_j = gamma_p[j];

        for (int s = 0; s < S; ++s) {
          const double mu_node = mu_p[s + S * k] +
            gamma_j * tau_between_p[s + S * k];
          const double marginal_loglik = binom_marginal_loglik_scalar(
            ai_p[k], ci_p[k], n1i_p[k], n2i_p[k], log_coef[k],
            mu_node, tau_within_p[s + S * k], weight_k, theta_p,
            logit_pi_p + n_pi * k,
            grid_weights.data() + static_cast<size_t>(k) * static_cast<size_t>(n_grid),
            n_theta, n_pi, terms, half_effect
          );
          cluster_terms[s + S * j] += marginal_loglik;
        }
      }
    }

    for (int s = 0; s < S; ++s) {
      for (int j = 0; j < n_gamma; ++j) {
        gamma_terms[j] = cluster_terms[s + S * j] + log_gamma_p[j];
      }
      out_p[s + S * g] = buffered_logsumexp_glmm(gamma_terms);
    }

    cluster_offset += cluster_n;
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_glmm_pois_cluster_loglik(SEXP x1i, SEXP x2i,
                                                SEXP t1i, SEXP t2i,
                                                SEXP mu_samples,
                                                SEXP tau_within,
                                                SEXP tau_between,
                                                SEXP cluster_index,
                                                SEXP cluster_size,
                                                SEXP weights,
                                                SEXP theta_grid,
                                                SEXP log_theta_weights,
                                                SEXP log_phi_grid,
                                                SEXP log_phi_weights,
                                                SEXP gamma_grid,
                                                SEXP log_gamma_weights)
{
  check_integer_glmm(x1i, "x1i");
  check_integer_glmm(x2i, "x2i");
  check_real_glmm(t1i, "t1i");
  check_real_glmm(t2i, "t2i");
  check_real_glmm(log_phi_grid, "log_phi_grid");
  check_real_glmm(log_phi_weights, "log_phi_weights");

  int S, K, n_theta, n_gamma, G, n_cluster_index;
  validate_cluster_common_glmm(
    mu_samples, tau_within, tau_between, theta_grid, log_theta_weights,
    gamma_grid, log_gamma_weights, cluster_index, cluster_size, weights,
    &S, &K, &n_theta, &n_gamma, &G, &n_cluster_index
  );
  (void)n_cluster_index;

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

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, G));

  const int *x1i_p             = INTEGER(x1i);
  const int *x2i_p             = INTEGER(x2i);
  const int *cluster_index_p   = INTEGER(cluster_index);
  const int *cluster_size_p    = INTEGER(cluster_size);
  const double *t1i_p          = REAL(t1i);
  const double *t2i_p          = REAL(t2i);
  const double *mu_p           = REAL(mu_samples);
  const double *tau_within_p   = REAL(tau_within);
  const double *tau_between_p  = REAL(tau_between);
  const double *theta_p        = REAL(theta_grid);
  const double *log_theta_p    = REAL(log_theta_weights);
  const double *log_phi_p      = REAL(log_phi_grid);
  const double *log_phi_w_p    = REAL(log_phi_weights);
  const double *gamma_p        = REAL(gamma_grid);
  const double *log_gamma_p    = REAL(log_gamma_weights);
  const double *weights_p      = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p                = REAL(out);

  const int n_grid = n_phi * n_theta;
  std::vector<double> terms(static_cast<size_t>(n_grid));
  std::vector<double> gamma_terms(static_cast<size_t>(n_gamma));
  std::vector<double> cluster_terms(static_cast<size_t>(S) * static_cast<size_t>(n_gamma));
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

  int cluster_offset = 0;
  for (int g = 0; g < G; ++g) {
    const int cluster_n = cluster_size_p[g];
    std::fill(cluster_terms.begin(), cluster_terms.end(), 0.0);

    for (int m = 0; m < cluster_n; ++m) {
      const int k = cluster_index_p[cluster_offset + m] - 1;
      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];

      for (int j = 0; j < n_gamma; ++j) {
        const double gamma_j = gamma_p[j];

        for (int s = 0; s < S; ++s) {
          const double mu_node = mu_p[s + S * k] +
            gamma_j * tau_between_p[s + S * k];
          const double marginal_loglik = pois_marginal_loglik_scalar(
            x1i_p[k], x2i_p[k], log_t1[k], log_t2[k], log_factorial[k],
            mu_node, tau_within_p[s + S * k], weight_k, theta_p,
            log_phi_p + n_phi * k,
            grid_weights.data() + static_cast<size_t>(k) * static_cast<size_t>(n_grid),
            n_theta, n_phi, terms, half_effect
          );
          cluster_terms[s + S * j] += marginal_loglik;
        }
      }
    }

    for (int s = 0; s < S; ++s) {
      for (int j = 0; j < n_gamma; ++j) {
        gamma_terms[j] = cluster_terms[s + S * j] + log_gamma_p[j];
      }
      out_p[s + S * g] = buffered_logsumexp_glmm(gamma_terms);
    }

    cluster_offset += cluster_n;
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_glmm_binom_cluster_loglik_row_sum(
  SEXP ai, SEXP ci, SEXP n1i, SEXP n2i, SEXP mu_samples,
  SEXP tau_within, SEXP tau_between, SEXP cluster_index,
  SEXP cluster_size, SEXP weights, SEXP theta_grid,
  SEXP log_theta_weights, SEXP logit_pi_grid, SEXP log_pi_weights,
  SEXP gamma_grid, SEXP log_gamma_weights)
{
  check_integer_glmm(ai, "ai");
  check_integer_glmm(ci, "ci");
  check_integer_glmm(n1i, "n1i");
  check_integer_glmm(n2i, "n2i");
  check_real_glmm(logit_pi_grid, "logit_pi_grid");
  check_real_glmm(log_pi_weights, "log_pi_weights");

  int S, K, n_theta, n_gamma, G, n_cluster_index;
  validate_cluster_common_glmm(
    mu_samples, tau_within, tau_between, theta_grid, log_theta_weights,
    gamma_grid, log_gamma_weights, cluster_index, cluster_size, weights,
    &S, &K, &n_theta, &n_gamma, &G, &n_cluster_index
  );
  (void)n_cluster_index;

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

  const int *ai_p              = INTEGER(ai);
  const int *ci_p              = INTEGER(ci);
  const int *n1i_p             = INTEGER(n1i);
  const int *n2i_p             = INTEGER(n2i);
  const int *cluster_index_p   = INTEGER(cluster_index);
  const int *cluster_size_p    = INTEGER(cluster_size);
  const double *mu_p           = REAL(mu_samples);
  const double *tau_within_p   = REAL(tau_within);
  const double *tau_between_p  = REAL(tau_between);
  const double *theta_p        = REAL(theta_grid);
  const double *log_theta_p    = REAL(log_theta_weights);
  const double *logit_pi_p     = REAL(logit_pi_grid);
  const double *log_pi_w_p     = REAL(log_pi_weights);
  const double *gamma_p        = REAL(gamma_grid);
  const double *log_gamma_p    = REAL(log_gamma_weights);
  const double *weights_p      = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p                = REAL(out);

  std::fill(out_p, out_p + S, 0.0);

  const int n_grid = n_pi * n_theta;
  std::vector<double> terms(static_cast<size_t>(n_grid));
  std::vector<double> gamma_terms(static_cast<size_t>(n_gamma));
  std::vector<double> cluster_terms(static_cast<size_t>(S) * static_cast<size_t>(n_gamma));
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

  int cluster_offset = 0;
  for (int g = 0; g < G; ++g) {
    const int cluster_n = cluster_size_p[g];
    std::fill(cluster_terms.begin(), cluster_terms.end(), 0.0);

    for (int m = 0; m < cluster_n; ++m) {
      const int k = cluster_index_p[cluster_offset + m] - 1;
      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];

      for (int j = 0; j < n_gamma; ++j) {
        const double gamma_j = gamma_p[j];

        for (int s = 0; s < S; ++s) {
          const double mu_node = mu_p[s + S * k] +
            gamma_j * tau_between_p[s + S * k];
          const double marginal_loglik = binom_marginal_loglik_scalar(
            ai_p[k], ci_p[k], n1i_p[k], n2i_p[k], log_coef[k],
            mu_node, tau_within_p[s + S * k], weight_k, theta_p,
            logit_pi_p + n_pi * k,
            grid_weights.data() + static_cast<size_t>(k) * static_cast<size_t>(n_grid),
            n_theta, n_pi, terms, half_effect
          );
          cluster_terms[s + S * j] += marginal_loglik;
        }
      }
    }

    for (int s = 0; s < S; ++s) {
      for (int j = 0; j < n_gamma; ++j) {
        gamma_terms[j] = cluster_terms[s + S * j] + log_gamma_p[j];
      }
      out_p[s] += buffered_logsumexp_glmm(gamma_terms);
    }

    cluster_offset += cluster_n;
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_glmm_pois_cluster_loglik_row_sum(
  SEXP x1i, SEXP x2i, SEXP t1i, SEXP t2i, SEXP mu_samples,
  SEXP tau_within, SEXP tau_between, SEXP cluster_index,
  SEXP cluster_size, SEXP weights, SEXP theta_grid,
  SEXP log_theta_weights, SEXP log_phi_grid, SEXP log_phi_weights,
  SEXP gamma_grid, SEXP log_gamma_weights)
{
  check_integer_glmm(x1i, "x1i");
  check_integer_glmm(x2i, "x2i");
  check_real_glmm(t1i, "t1i");
  check_real_glmm(t2i, "t2i");
  check_real_glmm(log_phi_grid, "log_phi_grid");
  check_real_glmm(log_phi_weights, "log_phi_weights");

  int S, K, n_theta, n_gamma, G, n_cluster_index;
  validate_cluster_common_glmm(
    mu_samples, tau_within, tau_between, theta_grid, log_theta_weights,
    gamma_grid, log_gamma_weights, cluster_index, cluster_size, weights,
    &S, &K, &n_theta, &n_gamma, &G, &n_cluster_index
  );
  (void)n_cluster_index;

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

  const int *x1i_p             = INTEGER(x1i);
  const int *x2i_p             = INTEGER(x2i);
  const int *cluster_index_p   = INTEGER(cluster_index);
  const int *cluster_size_p    = INTEGER(cluster_size);
  const double *t1i_p          = REAL(t1i);
  const double *t2i_p          = REAL(t2i);
  const double *mu_p           = REAL(mu_samples);
  const double *tau_within_p   = REAL(tau_within);
  const double *tau_between_p  = REAL(tau_between);
  const double *theta_p        = REAL(theta_grid);
  const double *log_theta_p    = REAL(log_theta_weights);
  const double *log_phi_p      = REAL(log_phi_grid);
  const double *log_phi_w_p    = REAL(log_phi_weights);
  const double *gamma_p        = REAL(gamma_grid);
  const double *log_gamma_p    = REAL(log_gamma_weights);
  const double *weights_p      = weights == R_NilValue ? 0 : REAL(weights);
  double *out_p                = REAL(out);

  std::fill(out_p, out_p + S, 0.0);

  const int n_grid = n_phi * n_theta;
  std::vector<double> terms(static_cast<size_t>(n_grid));
  std::vector<double> gamma_terms(static_cast<size_t>(n_gamma));
  std::vector<double> cluster_terms(static_cast<size_t>(S) * static_cast<size_t>(n_gamma));
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

  int cluster_offset = 0;
  for (int g = 0; g < G; ++g) {
    const int cluster_n = cluster_size_p[g];
    std::fill(cluster_terms.begin(), cluster_terms.end(), 0.0);

    for (int m = 0; m < cluster_n; ++m) {
      const int k = cluster_index_p[cluster_offset + m] - 1;
      const double weight_k = weights_p == 0 ? 1.0 : weights_p[k];

      for (int j = 0; j < n_gamma; ++j) {
        const double gamma_j = gamma_p[j];

        for (int s = 0; s < S; ++s) {
          const double mu_node = mu_p[s + S * k] +
            gamma_j * tau_between_p[s + S * k];
          const double marginal_loglik = pois_marginal_loglik_scalar(
            x1i_p[k], x2i_p[k], log_t1[k], log_t2[k], log_factorial[k],
            mu_node, tau_within_p[s + S * k], weight_k, theta_p,
            log_phi_p + n_phi * k,
            grid_weights.data() + static_cast<size_t>(k) * static_cast<size_t>(n_grid),
            n_theta, n_phi, terms, half_effect
          );
          cluster_terms[s + S * j] += marginal_loglik;
        }
      }
    }

    for (int s = 0; s < S; ++s) {
      for (int j = 0; j < n_gamma; ++j) {
        gamma_terms[j] = cluster_terms[s + S * j] + log_gamma_p[j];
      }
      out_p[s] += buffered_logsumexp_glmm(gamma_terms);
    }

    cluster_offset += cluster_n;
  }

  UNPROTECT(1);
  return out;
}
