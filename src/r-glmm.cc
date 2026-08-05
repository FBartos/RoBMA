#include <Rinternals.h>
#include <R_ext/Error.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "glmm-binomial-loglik.h"

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

static void validate_cluster_common_glmm(SEXP mu_samples, SEXP tau_within,
                                         SEXP tau_between, SEXP theta_grid,
                                         SEXP log_theta_weights,
                                         SEXP gamma_grid,
                                         SEXP log_gamma_weights,
                                         SEXP cluster_index,
                                         SEXP cluster_size,
                                         SEXP weights,
                                         int *S, int *K, int *n_theta,
                                         int *n_gamma, int *G,
                                         int *n_cluster_index)
{
  validate_common_glmm(
    mu_samples, tau_within, theta_grid, log_theta_weights,
    S, K, n_theta
  );
  check_real_glmm(tau_between, "tau_between");
  check_real_glmm(gamma_grid, "gamma_grid");
  check_real_glmm(log_gamma_weights, "log_gamma_weights");
  check_integer_glmm(cluster_index, "cluster_index");
  check_integer_glmm(cluster_size, "cluster_size");

  if (matrix_nrow_glmm(tau_between, "tau_between") != *S ||
      matrix_ncol_glmm(tau_between, "tau_between") != *K) {
    Rf_error("'tau_between' dimensions must match 'mu_samples'.");
  }

  *n_gamma = Rf_length(gamma_grid);
  if (Rf_length(log_gamma_weights) != *n_gamma) {
    Rf_error("'gamma_grid' and 'log_gamma_weights' lengths must match.");
  }

  *G = Rf_length(cluster_size);
  *n_cluster_index = Rf_length(cluster_index);
  const int *cluster_index_p = INTEGER(cluster_index);
  const int *cluster_size_p  = INTEGER(cluster_size);
  int expected_index = 0;

  for (int g = 0; g < *G; ++g) {
    if (cluster_size_p[g] == NA_INTEGER || cluster_size_p[g] <= 0) {
      Rf_error("'cluster_size' must contain positive integers.");
    }
    expected_index += cluster_size_p[g];
  }
  if (expected_index != *n_cluster_index) {
    Rf_error("'cluster_index' length must equal sum('cluster_size').");
  }
  for (int i = 0; i < *n_cluster_index; ++i) {
    if (cluster_index_p[i] == NA_INTEGER ||
        cluster_index_p[i] < 1 || cluster_index_p[i] > *K) {
      Rf_error("'cluster_index' contains an invalid estimate index.");
    }
  }

  validate_weights_glmm(weights, *K);
}

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
