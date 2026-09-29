// Normal cluster-unit likelihoods integrate one Gaussian gamma shared by all
// observations in a cluster: exactly through the determinant lemma when the
// rows are unweighted, and through the supplied fixed Gaussian rule otherwise.
// Selection models use the joint selection kernels instead.

struct NormClusterArgs {
  int S;
  int K;
  int G;
  int n_gamma;
  const double *yi;
  const double *sei;
  const double *mu;
  const double *tau_within;
  const double *tau_between;
  const int *cluster_index;
  const int *cluster_size;
  const double *weights;
  const double *gamma;
  const double *log_gamma;
};

struct NormClusterAnalyticArgs {
  int S;
  int K;
  int G;
  const double *yi;
  const double *vi;
  const double *mu;
  const double *tau_within;
  const double *tau_between;
  const double *tau_total;
  const int *cluster_index;
  const int *cluster_size;
};

double norm_cluster_logsumexp(const std::vector<double> &terms, int n = -1)
{
  const int n_terms = n < 0 ? static_cast<int>(terms.size()) : n;
  double term_max = -INFINITY;

  for (int j = 0; j < n_terms; ++j) {
    if (std::isnan(terms[j])) {
      return NA_REAL;
    }
    if (terms[j] > term_max) {
      term_max = terms[j];
    }
  }

  if (std::isinf(term_max)) {
    return term_max;
  }

  double term_sum = 0;
  for (int j = 0; j < n_terms; ++j) {
    term_sum += std::exp(terms[j] - term_max);
  }

  return term_max + std::log(term_sum);
}

void validate_cluster_common_data(SEXP yi, SEXP sampling, const char *sampling_name,
                                  SEXP mu_samples, SEXP tau_within,
                                  SEXP tau_between, SEXP cluster_index,
                                  SEXP cluster_size, SEXP weights,
                                  int *S, int *K, int *G)
{
  check_real(yi, "yi");
  check_real(sampling, sampling_name);
  check_real(mu_samples, "mu_samples");
  check_real(tau_within, "tau_within");
  check_real(tau_between, "tau_between");
  check_integer(cluster_index, "cluster_index");
  check_integer(cluster_size, "cluster_size");

  if (weights != R_NilValue) {
    check_real(weights, "weights");
  }

  *S = matrix_nrow(mu_samples, "mu_samples");
  *K = matrix_ncol(mu_samples, "mu_samples");
  *G = Rf_length(cluster_size);

  validate_matrix_dims_named(tau_within, "tau_within", *S, *K, "mu_samples");
  validate_matrix_dims_named(tau_between, "tau_between", *S, *K, "mu_samples");

  if (Rf_length(yi) != *K || Rf_length(sampling) != *K) {
    Rf_error("'yi' and '%s' must have one value per observation.", sampling_name);
  }
  if (weights != R_NilValue && Rf_length(weights) != *K) {
    Rf_error("'weights' must have one value per observation.");
  }

  const int *cluster_index_p = INTEGER(cluster_index);
  const int *cluster_size_p  = INTEGER(cluster_size);
  const int n_cluster_index  = Rf_length(cluster_index);
  int cluster_total = 0;

  for (int g = 0; g < *G; ++g) {
    if (cluster_size_p[g] <= 0) {
      Rf_error("'cluster_size' must contain positive cluster sizes.");
    }
    cluster_total += cluster_size_p[g];
  }
  if (cluster_total != n_cluster_index) {
    Rf_error("'cluster_size' must sum to the length of 'cluster_index'.");
  }
  for (int i = 0; i < n_cluster_index; ++i) {
    if (cluster_index_p[i] < 1 || cluster_index_p[i] > *K) {
      Rf_error("'cluster_index' values must be between 1 and the number of observations.");
    }
  }
}

void validate_cluster_common_norm(SEXP yi, SEXP sei, SEXP mu_samples,
                                  SEXP tau_within, SEXP tau_between,
                                  SEXP cluster_index, SEXP cluster_size,
                                  SEXP weights, SEXP gamma_grid,
                                  SEXP log_gamma_weights,
                                  int *S, int *K, int *G, int *n_gamma)
{
  validate_cluster_common_data(
    yi, sei, "sei", mu_samples, tau_within, tau_between, cluster_index,
    cluster_size, weights, S, K, G
  );
  check_real(gamma_grid, "gamma_grid");
  check_real(log_gamma_weights, "log_gamma_weights");

  *n_gamma = Rf_length(gamma_grid);
  if (*n_gamma < 1 || Rf_length(log_gamma_weights) != *n_gamma) {
    Rf_error("'gamma_grid' and 'log_gamma_weights' must have matching positive lengths.");
  }
}

const double NORM_CLUSTER_LOG_2PI = std::log(2.0 * std::acos(-1.0));

// Integrate a normal cluster effect exactly through the determinant lemma and
// the stable residual-plus-projection form of the Sherman-Morrison quadratic.
double norm_cluster_analytic_loglik(const NormClusterAnalyticArgs &args, int s,
                                    int cluster_offset, int cluster_n,
                                    double rho = NA_REAL)
{
  const bool allocation_grid = args.tau_total != 0;
  const double one_minus_rho = allocation_grid ? 1.0 - rho : NA_REAL;
  const double sqrt_rho      = allocation_grid ? std::sqrt(rho) : NA_REAL;
  double log_det = 0.0;
  double denom   = 1.0;
  double score   = 0.0;

  for (int m = 0; m < cluster_n; ++m) {
    const int k = args.cluster_index[cluster_offset + m] - 1;
    const int cell = s + args.S * k;
    const double tau_total = allocation_grid ?
      args.tau_total[cell] : NA_REAL;
    const double diagonal = allocation_grid ?
      args.vi[k] + tau_total * tau_total * one_minus_rho :
      args.vi[k] + args.tau_within[cell] * args.tau_within[cell];
    const double rank_one = allocation_grid ?
      tau_total * sqrt_rho : args.tau_between[cell];
    const double residual = args.yi[k] - args.mu[cell];
    if (!std::isfinite(diagonal) || diagonal <= 0.0 ||
        !std::isfinite(rank_one) || !std::isfinite(residual) ||
        (allocation_grid &&
         (!std::isfinite(tau_total) || tau_total < 0.0))) {
      return NA_REAL;
    }

    const double inv_diagonal = 1.0 / diagonal;
    log_det += std::log(diagonal);
    denom   += rank_one * rank_one * inv_diagonal;
    score   += rank_one * residual * inv_diagonal;
  }
  if (!std::isfinite(denom) || denom <= 0.0) {
    return NA_REAL;
  }

  const double projection = score / denom;
  double quadratic = projection * projection;
  for (int m = 0; m < cluster_n; ++m) {
    const int k = args.cluster_index[cluster_offset + m] - 1;
    const int cell = s + args.S * k;
    const double tau_total = allocation_grid ?
      args.tau_total[cell] : NA_REAL;
    const double diagonal = allocation_grid ?
      args.vi[k] + tau_total * tau_total * one_minus_rho :
      args.vi[k] + args.tau_within[cell] * args.tau_within[cell];
    const double adjusted_residual = args.yi[k] - args.mu[cell] -
      projection * (allocation_grid ?
        tau_total * sqrt_rho : args.tau_between[cell]);
    quadratic += adjusted_residual * adjusted_residual / diagonal;
  }

  const double out = -0.5 * (
    cluster_n * NORM_CLUSTER_LOG_2PI + log_det + std::log(denom) + quadratic
  );
  return std::isfinite(out) ? out : NA_REAL;
}

SEXP norm_cluster_analytic_accumulate(const NormClusterAnalyticArgs &args,
                                      bool row_sum)
{
  SEXP out = PROTECT(row_sum ?
    Rf_allocVector(REALSXP, args.S) :
    allocate_selnorm_matrix(REALSXP, args.S, args.G));
  double *out_p = REAL(out);
  if (row_sum) {
    std::fill(out_p, out_p + args.S, 0.0);
  }

  int cluster_offset = 0;
  for (int g = 0; g < args.G; ++g) {
    const int cluster_n = args.cluster_size[g];
    for (int s = 0; s < args.S; ++s) {
      const double value = norm_cluster_analytic_loglik(
        args, s, cluster_offset, cluster_n
      );
      if (row_sum) {
        out_p[s] += value;
      } else {
        out_p[s + args.S * g] = value;
      }
    }
    cluster_offset += cluster_n;
  }

  UNPROTECT(1);
  return out;
}

SEXP norm_cluster_analytic_rho_grid_accumulate(
    const NormClusterAnalyticArgs &args, const double *rho, int n_rho)
{
  SEXP out = PROTECT(allocate_selnorm_matrix(REALSXP, n_rho, args.S));
  double *out_p = REAL(out);

  for (int s = 0; s < args.S; ++s) {
    for (int r = 0; r < n_rho; ++r) {
      double row_sum = 0.0;
      int cluster_offset = 0;
      for (int g = 0; g < args.G; ++g) {
        const int cluster_n = args.cluster_size[g];
        row_sum += norm_cluster_analytic_loglik(
          args, s, cluster_offset, cluster_n, rho[r]
        );
        cluster_offset += cluster_n;
      }
      out_p[r + n_rho * s] = row_sum;
    }
  }

  UNPROTECT(1);
  return out;
}

double norm_cluster_cell_loglik(const NormClusterArgs &args, int k,
                                double mu_node, double sd)
{
  const double weight = args.weights == 0 ? 1.0 : args.weights[k];

  return weight * Rf_dnorm4(args.yi[k], mu_node, sd, 1);
}

SEXP norm_cluster_accumulate(const NormClusterArgs &args)
{
  SEXP out = PROTECT(allocate_selnorm_matrix(REALSXP, args.S, args.G));
  double *out_p = REAL(out);

  std::vector<double> gamma_terms(static_cast<size_t>(args.n_gamma));
  std::vector<double> cluster_terms(
    static_cast<size_t>(args.S) * static_cast<size_t>(args.n_gamma)
  );

  int cluster_offset = 0;
  for (int g = 0; g < args.G; ++g) {
    const int cluster_n = args.cluster_size[g];
    std::fill(cluster_terms.begin(), cluster_terms.end(), 0.0);

    for (int m = 0; m < cluster_n; ++m) {
      const int k = args.cluster_index[cluster_offset + m] - 1;

      for (int j = 0; j < args.n_gamma; ++j) {
        const double gamma_j = args.gamma[j];

        for (int s = 0; s < args.S; ++s) {
          const int cell = s + args.S * k;
          const double mu_node = args.mu[cell] +
            gamma_j * args.tau_between[cell];
          const double sd = std::sqrt(
            args.tau_within[cell] * args.tau_within[cell] +
              args.sei[k] * args.sei[k]
          );

          cluster_terms[s + args.S * j] += norm_cluster_cell_loglik(
            args, k, mu_node, sd
          );
        }
      }
    }

    for (int s = 0; s < args.S; ++s) {
      for (int j = 0; j < args.n_gamma; ++j) {
        gamma_terms[j] = cluster_terms[s + args.S * j] + args.log_gamma[j];
      }
      out_p[s + args.S * g] = norm_cluster_logsumexp(gamma_terms);
    }

    cluster_offset += cluster_n;
  }

  UNPROTECT(1);
  return out;
}

SEXP norm_cluster_accumulate_row_sum(const NormClusterArgs &args)
{
  SEXP out = PROTECT(Rf_allocVector(REALSXP, args.S));
  double *out_p = REAL(out);
  std::fill(out_p, out_p + args.S, 0.0);

  std::vector<double> gamma_terms(static_cast<size_t>(args.n_gamma));
  std::vector<double> cluster_terms(
    static_cast<size_t>(args.S) * static_cast<size_t>(args.n_gamma)
  );

  int cluster_offset = 0;
  for (int g = 0; g < args.G; ++g) {
    const int cluster_n = args.cluster_size[g];
    std::fill(cluster_terms.begin(), cluster_terms.end(), 0.0);

    for (int m = 0; m < cluster_n; ++m) {
      const int k = args.cluster_index[cluster_offset + m] - 1;

      for (int j = 0; j < args.n_gamma; ++j) {
        const double gamma_j = args.gamma[j];

        for (int s = 0; s < args.S; ++s) {
          const int cell = s + args.S * k;
          const double mu_node = args.mu[cell] +
            gamma_j * args.tau_between[cell];
          const double sd = std::sqrt(
            args.tau_within[cell] * args.tau_within[cell] +
              args.sei[k] * args.sei[k]
          );

          cluster_terms[s + args.S * j] += norm_cluster_cell_loglik(
            args, k, mu_node, sd
          );
        }
      }
    }

    for (int s = 0; s < args.S; ++s) {
      for (int j = 0; j < args.n_gamma; ++j) {
        gamma_terms[j] = cluster_terms[s + args.S * j] + args.log_gamma[j];
      }
      out_p[s] += norm_cluster_logsumexp(gamma_terms);
    }

    cluster_offset += cluster_n;
  }

  UNPROTECT(1);
  return out;
}
