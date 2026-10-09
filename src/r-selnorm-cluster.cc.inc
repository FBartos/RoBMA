// Normal cluster-unit likelihoods integrate one Gaussian gamma shared by all
// observations in a cluster in closed form, for unit and for general powers:
// the powered Gaussian integral is evaluated by the stable one-factor routine
// of the selected-normal kernels (cpp_gaussian_latent_loglik), whose row
// standard deviations combine the sampling and the within-cluster scale without
// squaring them. Selection models use the joint selection kernels instead.

struct NormClusterAnalyticArgs {
  int S;
  int K;
  int G;
  const double *yi;
  const double *sei;
  const double *mu;
  const double *tau_within;
  const double *tau_between;
  const double *tau_total;
  const double *weights;
  const int *cluster_index;
  const int *cluster_size;
};

// Per-cluster rows gathered for the stable routine: residual, row sd, loading
// and power.
struct NormClusterRows {
  std::vector<double> residual;
  std::vector<double> sd;
  std::vector<double> loading;
  std::vector<double> power;
  std::vector<double> work;

  explicit NormClusterRows(int capacity)
    : residual(static_cast<size_t>(capacity)),
      sd(static_cast<size_t>(capacity)),
      loading(static_cast<size_t>(capacity)),
      power(static_cast<size_t>(capacity)),
      work(static_cast<size_t>(capacity)) {}
};

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

// Row sd sqrt(sei^2 + tau^2): the sum of squares while both operands keep
// their squares representable, std::hypot otherwise.
inline double norm_cluster_row_sd(double sei, double tau)
{
  const double a = std::fabs(sei);
  const double b = std::fabs(tau);
  if (b == 0.0) {
    return a;
  }
  if (a >= 1e-150 && a <= 1e150 && b >= 1e-150 && b <= 1e150) {
    return std::sqrt(a * a + b * b);
  }
  return std::hypot(a, b);
}

// Integrate a normal cluster effect in closed form. Invalid rows (a
// non-finite value, a non-positive row sd, a negative power) give NA, which
// the R wrappers resolve with their reference implementation.
double norm_cluster_analytic_loglik(const NormClusterAnalyticArgs &args, int s,
                                    int cluster_offset, int cluster_n,
                                    NormClusterRows &rows,
                                    double rho = NA_REAL)
{
  const bool allocation_grid = args.tau_total != 0;
  const double within_scale  = allocation_grid ? std::sqrt(1.0 - rho) : NA_REAL;
  const double between_scale = allocation_grid ? std::sqrt(rho) : NA_REAL;

  for (int m = 0; m < cluster_n; ++m) {
    const int k = args.cluster_index[cluster_offset + m] - 1;
    const int cell = s + args.S * k;
    double within;
    double between;
    if (allocation_grid) {
      const double tau_total = args.tau_total[cell];
      if (!std::isfinite(tau_total) || tau_total < 0.0) {
        return NA_REAL;
      }
      within  = tau_total * within_scale;
      between = tau_total * between_scale;
    } else {
      within  = args.tau_within[cell];
      between = args.tau_between[cell];
    }
    rows.residual[m] = args.yi[k] - args.mu[cell];
    rows.sd[m]       = norm_cluster_row_sd(args.sei[k], within);
    rows.loading[m]  = between;
    rows.power[m]    = args.weights == 0 ? 1.0 : args.weights[k];
  }

  const double out = cpp_gaussian_latent_loglik(
    cluster_n, rows.residual.data(), rows.sd.data(), rows.loading.data(),
    args.weights == 0 ? nullptr : rows.power.data(), rows.work.data()
  );
  // An unrepresentable result is -Inf; only an invalid row gives NA.
  return (std::isfinite(out) || out == -std::numeric_limits<double>::infinity()) ?
    out : NA_REAL;
}

int norm_cluster_max_size(const NormClusterAnalyticArgs &args)
{
  int largest = 1;
  for (int g = 0; g < args.G; ++g) {
    largest = std::max(largest, args.cluster_size[g]);
  }
  return largest;
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

  NormClusterRows rows(norm_cluster_max_size(args));
  int cluster_offset = 0;
  for (int g = 0; g < args.G; ++g) {
    const int cluster_n = args.cluster_size[g];
    for (int s = 0; s < args.S; ++s) {
      const double value = norm_cluster_analytic_loglik(
        args, s, cluster_offset, cluster_n, rows
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

  NormClusterRows rows(norm_cluster_max_size(args));
  for (int s = 0; s < args.S; ++s) {
    for (int r = 0; r < n_rho; ++r) {
      double row_sum = 0.0;
      int cluster_offset = 0;
      for (int g = 0; g < args.G; ++g) {
        const int cluster_n = args.cluster_size[g];
        row_sum += norm_cluster_analytic_loglik(
          args, s, cluster_offset, cluster_n, rows, rho[r]
        );
        cluster_offset += cluster_n;
      }
      out_p[r + n_rho * s] = row_sum;
    }
  }

  UNPROTECT(1);
  return out;
}
