#include <Rinternals.h>
#include <R_ext/Error.h>

#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <limits>
#include <vector>

#include "plot-root.h"
#include "selnorm/selnorm.h"

static int regplot_matrix_nrow(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[0];
}

static int regplot_matrix_ncol(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[1];
}

static void regplot_check_real(SEXP x, const char *name)
{
  if (TYPEOF(x) != REALSXP) {
    Rf_error("'%s' must be numeric.", name);
  }
}

static void regplot_check_integer(SEXP x, const char *name)
{
  if (TYPEOF(x) != INTSXP) {
    Rf_error("'%s' must be integer.", name);
  }
}

static void regplot_check_logical(SEXP x, const char *name)
{
  if (TYPEOF(x) != LGLSXP || Rf_length(x) != 1) {
    Rf_error("'%s' must be a logical scalar.", name);
  }
}

static int regplot_scalar_or_row_int(const int *x, int n, int row,
                                     const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (row >= n) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
  return x[row];
}

static double regplot_scalar_or_row_real(const double *x, int n, int row,
                                         const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (row >= n) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
  return x[row];
}

static void regplot_validate_omega_matrix(const double *omega, int S, int B)
{
  for (int b = 0; b < B; ++b) {
    for (int s = 0; s < S; ++s) {
      const double w = omega[s + S * b];
      if (!std::isfinite(w) || w < 0) {
        Rf_error("'omega' must contain finite non-negative values.");
      }
    }
  }
}

static void regplot_set_kernel_data(SelNormKernelData *data, int B,
                                    int n_segments, SEXP sign, SEXP q,
                                    SEXP z_lower, SEXP z_upper,
                                    SEXP phack_z_source,
                                    SEXP phack_z_dest,
                                    SEXP segment_bounds,
                                    SEXP segment_step_bin,
                                    SEXP segment_phack_region,
                                    SEXP telescope_probabilities)
{
  data->n_bins                 = B;
  data->n_segments             = n_segments;
  data->effect_sign            = INTEGER(sign)[0];
  data->q                      = INTEGER(q)[0];
  data->z_lower                = REAL(z_lower);
  data->z_upper                = REAL(z_upper);
  data->phack_z_source         = REAL(phack_z_source);
  data->phack_z_dest           = REAL(phack_z_dest);
  data->segment_bounds         = REAL(segment_bounds);
  data->segment_step_bin       = INTEGER(segment_step_bin);
  data->segment_phack_region   = INTEGER(segment_phack_region);
  data->segment_step_bin_real  = 0;
  data->segment_phack_region_real = 0;
  data->trusted_step_partition = selnorm_is_descending_step_partition(
    data->z_lower, data->z_upper, data->n_bins
  );
  data->telescope_probabilities = data->trusted_step_partition &&
    LOGICAL(telescope_probabilities)[0] == TRUE;
}

static bool regplot_row_has_active_phack(int mode, double alpha, int phack_kind)
{
  return (mode == SELKERNEL_PHACK_POWER ||
    mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
}

static double regplot_probability(double x)
{
  return std::isfinite(x) && x >= 0 && x <= 1 ? x : NA_REAL;
}

static double regplot_empirical_quantile(const double *x, int n, double p)
{
  std::vector<double> sorted(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) {
      return NA_REAL;
    }
    sorted[i] = x[i];
  }
  std::sort(sorted.begin(), sorted.end());

  const int index = std::max(
    0,
    std::min(n - 1, static_cast<int>(std::ceil(static_cast<double>(n) * p)) - 1)
  );
  return sorted[static_cast<size_t>(index)];
}

static double regplot_normal_mixture_cdf(double q, const double *mean,
                                         const double *sd, int S)
{
  double cdf_sum = 0;

  for (int s = 0; s < S; ++s) {
    double cdf_s;
    if (sd[s] == 0) {
      cdf_s = q >= mean[s] ? 1.0 : 0.0;
    } else {
      cdf_s = pnorm(q, mean[s], sd[s], true, false);
    }
    cdf_sum += regplot_probability(cdf_s);
  }

  return cdf_sum / static_cast<double>(S);
}

static double regplot_selnorm_mixture_cdf_cached(
  double q,
  const double *mean,
  const double *sd,
  double se,
  const double *omega,
  const std::vector<double> &row_alpha,
  const std::vector<int> &row_phack,
  const std::vector<int> &row_mode,
  const std::vector<bool> &row_active_phack,
  const std::vector<char> &step_plan_valid,
  const std::vector<int> &step_plan_groups,
  const std::vector<double> &step_lower_score,
  const std::vector<double> &step_upper_score,
  const std::vector<double> &step_log_weight,
  const SelNormKernelData &data,
  int S)
{
  double cdf_sum = 0;

  for (int s = 0; s < S; ++s) {
    double cdf_s;
    if (sd[s] == 0) {
      cdf_s = q >= mean[s] ? 1.0 : 0.0;
    } else if (row_active_phack[static_cast<size_t>(s)]) {
      cdf_s = NA_REAL;
    } else if (step_plan_valid[static_cast<size_t>(s)]) {
      const size_t offset = static_cast<size_t>(s) *
        static_cast<size_t>(data.n_bins);
      cdf_s = cpp_selnorm_step_cdf_from_plan(
        q,
        mean[s],
        sd[s],
        se,
        data,
        step_lower_score.data() + offset,
        step_upper_score.data() + offset,
        step_log_weight.data() + offset,
        step_plan_groups[static_cast<size_t>(s)],
        true
      );
    } else {
      cdf_s = cpp_selnorm_kernel_cdf(
        q,
        mean[s],
        sd[s],
        se,
        omega + s,
        row_alpha[static_cast<size_t>(s)],
        row_phack[static_cast<size_t>(s)],
        row_mode[static_cast<size_t>(s)],
        data,
        S,
        true,
        false
      );
    }
    cdf_sum += regplot_probability(cdf_s);
  }

  return cdf_sum / static_cast<double>(S);
}

template <typename CdfFun>
static double regplot_mixture_quantile(double p, const double *mean,
                                       const double *sd, int S,
                                       bool full_support,
                                       CdfFun cdf_fun,
                                       double previous = NA_REAL,
                                       bool use_previous = false)
{
  bool all_zero_sd = true;
  bool all_positive_sd = true;
  double lower = INFINITY;
  double upper = -INFINITY;
  double step = 0;

  for (int s = 0; s < S; ++s) {
    all_zero_sd     = all_zero_sd && sd[s] == 0;
    all_positive_sd = all_positive_sd && sd[s] > 0;

    const double spread = sd[s];
    lower = std::min(lower, mean[s] - 10 * spread);
    upper = std::max(upper, mean[s] + 10 * spread);
    step  = std::max(step, spread);

  }

  if (all_zero_sd) {
    return regplot_empirical_quantile(mean, S, p);
  }

  if (!std::isfinite(lower) || !std::isfinite(upper) || lower >= upper) {
    return NA_REAL;
  }

  const double global_lower = lower;
  const double global_upper = upper;
  const double global_step  = step;
  bool using_previous = use_previous && std::isfinite(previous);

  if (using_previous) {
    lower = previous - step;
    upper = previous + step;
  }

  const double tolerance = RoBMA::plot_root_tolerance(
    global_lower, global_upper, global_step
  );

  double lower_value = cdf_fun(lower) - p;
  double upper_value = cdf_fun(upper) - p;
  if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    if (!using_previous) {
      return NA_REAL;
    }
    using_previous = false;
    lower = global_lower;
    upper = global_upper;
    step  = global_step;
    lower_value = cdf_fun(lower) - p;
    upper_value = cdf_fun(upper) - p;
    if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
      return NA_REAL;
    }
  }

  for (int i = 0; i < 25; ++i) {
    if (lower_value < 0 && upper_value >= 0) {
      break;
    }
    if (lower_value >= 0) {
      lower -= step;
      lower_value = cdf_fun(lower) - p;
      if (!std::isfinite(lower_value)) {
        if (using_previous) {
          break;
        }
        return NA_REAL;
      }
    }
    if (upper_value < 0) {
      upper += step;
      upper_value = cdf_fun(upper) - p;
      if (!std::isfinite(upper_value)) {
        if (using_previous) {
          break;
        }
        return NA_REAL;
      }
    }
    step *= 2;
  }

  if ((lower_value >= 0 || upper_value < 0 ||
       !std::isfinite(lower_value) || !std::isfinite(upper_value)) &&
      using_previous) {
    lower = global_lower;
    upper = global_upper;
    step  = global_step;
    lower_value = cdf_fun(lower) - p;
    upper_value = cdf_fun(upper) - p;

    for (int i = 0; i < 25; ++i) {
      if (lower_value < 0 && upper_value >= 0) {
        break;
      }
      if (lower_value >= 0) {
        lower -= step;
        lower_value = cdf_fun(lower) - p;
      }
      if (upper_value < 0) {
        upper += step;
        upper_value = cdf_fun(upper) - p;
      }
      step *= 2;
    }
  }

  if (lower_value >= 0 || upper_value < 0 ||
      !std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    return NA_REAL;
  }

  // Full-support continuous mixtures have a unique root and do not need exact
  // generalized-inverse bisection. Atoms and selection gaps retain the exact
  // path.
  if (all_positive_sd && full_support) {
    const double root = RoBMA::plot_brent_root(
      [cdf_fun, p](double q) {
        return cdf_fun(q) - p;
      },
      lower,
      upper,
      lower_value,
      upper_value,
      tolerance
    );
    if (std::isfinite(root)) {
      return root;
    }
  }

  while (true) {
    const double mid = 0.5 * lower + 0.5 * upper;
    if (mid <= lower || mid >= upper) {
      break;
    }
    const double mid_value = cdf_fun(mid) - p;
    if (!std::isfinite(mid_value)) {
      return NA_REAL;
    }

    if (mid_value >= 0) {
      upper = mid;
    } else {
      lower = mid;
    }

  }

  return upper;
}

static double regplot_normal_mixture_quantile(double p, const double *mean,
                                              const double *sd, int S,
                                              double previous = NA_REAL,
                                              bool use_previous = false)
{
  return regplot_mixture_quantile(
    p, mean, sd, S, true,
    [mean, sd, S](double q) {
      return regplot_normal_mixture_cdf(q, mean, sd, S);
    },
    previous, use_previous
  );
}

static bool regplot_selection_mixture_has_full_support(
  const double *omega,
  const std::vector<int> &row_mode,
  const std::vector<bool> &row_active_phack,
  int S,
  int B)
{
  for (int s = 0; s < S; ++s) {
    const int mode = row_mode[static_cast<size_t>(s)];
    if (mode == SELKERNEL_NORMAL ||
        (mode == SELKERNEL_PHACK_POWER &&
         !row_active_phack[static_cast<size_t>(s)])) {
      return true;
    }
  }

  for (int b = 0; b < B; ++b) {
    bool bin_has_support = false;
    for (int s = 0; s < S; ++s) {
      const int mode = row_mode[static_cast<size_t>(s)];
      if ((mode == SELKERNEL_STEP || mode == SELKERNEL_STEP_PHACK_POWER) &&
          omega[s + S * b] > 0) {
        bin_has_support = true;
        break;
      }
    }
    if (!bin_has_support) {
      return false;
    }
  }

  return true;
}

static SEXP regplot_interval_result(int K, std::vector<double> &lower,
                                    std::vector<double> &upper)
{
  SEXP lower_out = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP upper_out = PROTECT(Rf_allocVector(REALSXP, K));

  for (int k = 0; k < K; ++k) {
    if (!std::isfinite(lower[static_cast<size_t>(k)]) ||
        !std::isfinite(upper[static_cast<size_t>(k)]) ||
        lower[static_cast<size_t>(k)] > upper[static_cast<size_t>(k)]) {
      Rf_error(
        "Regression-plot quantiles could not be computed from a valid bracketed CDF."
      );
    }
    REAL(lower_out)[k] = lower[static_cast<size_t>(k)];
    REAL(upper_out)[k] = upper[static_cast<size_t>(k)];
  }

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, lower_out);
  SET_VECTOR_ELT(out, 1, upper_out);

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("lower"));
  SET_STRING_ELT(names, 1, Rf_mkChar("upper"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(4);
  return out;
}

static void regplot_validate_mean_sd(SEXP mean, SEXP sd, int *S, int *K)
{
  regplot_check_real(mean, "mean");
  regplot_check_real(sd, "sd");

  *S = regplot_matrix_nrow(mean, "mean");
  *K = regplot_matrix_ncol(mean, "mean");

  if (regplot_matrix_nrow(sd, "sd") != *S ||
      regplot_matrix_ncol(sd, "sd") != *K) {
    Rf_error("'sd' dimensions must match 'mean'.");
  }
}

static void regplot_validate_probs(SEXP probs)
{
  regplot_check_real(probs, "probs");
  if (Rf_length(probs) != 2) {
    Rf_error("'probs' must have length 2.");
  }
  if (!std::isfinite(REAL(probs)[0]) || !std::isfinite(REAL(probs)[1]) ||
      REAL(probs)[0] <= 0 || REAL(probs)[0] >= 1 ||
      REAL(probs)[1] <= 0 || REAL(probs)[1] >= 1) {
    Rf_error("'probs' must contain probabilities in (0, 1).");
  }
}

extern "C" SEXP RoBMA_regplot_normal_mixture_interval(SEXP mean, SEXP sd,
                                                       SEXP probs)
{
  int S, K;
  regplot_validate_mean_sd(mean, sd, &S, &K);
  regplot_validate_probs(probs);

  const double *mean_p = REAL(mean);
  const double *sd_p   = REAL(sd);
  const double p_lower = REAL(probs)[0];
  const double p_upper = REAL(probs)[1];

  std::vector<double> lower(static_cast<size_t>(K));
  std::vector<double> upper(static_cast<size_t>(K));

  for (int k = 0; k < K; ++k) {
    const double *mean_k = mean_p + static_cast<size_t>(S) * static_cast<size_t>(k);
    const double *sd_k   = sd_p   + static_cast<size_t>(S) * static_cast<size_t>(k);

    lower[static_cast<size_t>(k)] = regplot_normal_mixture_quantile(
      p_lower, mean_k, sd_k, S,
      k > 0 ? lower[static_cast<size_t>(k - 1)] : NA_REAL,
      k > 0
    );
    upper[static_cast<size_t>(k)] = regplot_normal_mixture_quantile(
      p_upper, mean_k, sd_k, S,
      k > 0 ? upper[static_cast<size_t>(k - 1)] : NA_REAL,
      k > 0
    );
  }

  return regplot_interval_result(K, lower, upper);
}

extern "C" SEXP RoBMA_regplot_selnorm_mixture_interval(
  SEXP mean, SEXP sd, SEXP se, SEXP probs, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  int S, K;
  regplot_validate_mean_sd(mean, sd, &S, &K);
  regplot_validate_probs(probs);
  regplot_check_real(se, "se");
  regplot_check_real(omega, "omega");
  regplot_check_real(alpha, "alpha");
  regplot_check_integer(phack_kind, "phack_kind");
  regplot_check_integer(kernel_mode, "kernel_mode");
  regplot_check_real(z_lower, "z_lower");
  regplot_check_real(z_upper, "z_upper");
  regplot_check_integer(sign, "sign");
  regplot_check_integer(q, "q");
  regplot_check_real(phack_z_source, "phack_z_source");
  regplot_check_real(phack_z_dest, "phack_z_dest");
  regplot_check_real(segment_bounds, "segment_bounds");
  regplot_check_integer(segment_step_bin, "segment_step_bin");
  regplot_check_integer(segment_phack_region, "segment_phack_region");
  regplot_check_logical(telescope_probabilities, "telescope_probabilities");

  if ((Rf_length(se) != 1 && Rf_length(se) != K)) {
    Rf_error("'se' must have length 1 or one value per prediction.");
  }
  for (int k = 0; k < Rf_length(se); ++k) {
    if (!(REAL(se)[k] > 0) || !std::isfinite(REAL(se)[k])) {
      Rf_error("'se' must contain finite positive values.");
    }
  }

  const int B = regplot_matrix_ncol(omega, "omega");
  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }
  if (regplot_matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(alpha) != 1 && Rf_length(alpha) != S) {
    Rf_error("'alpha' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(phack_kind) != 1 && Rf_length(phack_kind) != S) {
    Rf_error("'phack_kind' must have length 1 or one value per posterior sample.");
  }
  if (Rf_length(kernel_mode) != 1 && Rf_length(kernel_mode) != S) {
    Rf_error("'kernel_mode' must have length 1 or one value per posterior sample.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  const double *mean_p       = REAL(mean);
  const double *sd_p         = REAL(sd);
  const double *se_p         = REAL(se);
  const double *omega_p      = REAL(omega);
  const double *alpha_p      = REAL(alpha);
  const int *phack_kind_p    = INTEGER(phack_kind);
  const int *kernel_mode_p   = INTEGER(kernel_mode);
  const double p_lower       = REAL(probs)[0];
  const double p_upper       = REAL(probs)[1];

  regplot_validate_omega_matrix(omega_p, S, B);

  SelNormKernelData data;
  regplot_set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region, telescope_probabilities
  );

  std::vector<double> lower(static_cast<size_t>(K));
  std::vector<double> upper(static_cast<size_t>(K));
  std::vector<double> row_alpha(static_cast<size_t>(S));
  std::vector<int> row_phack(static_cast<size_t>(S));
  std::vector<int> row_mode(static_cast<size_t>(S));
  std::vector<bool> row_active_phack(static_cast<size_t>(S));
  std::vector<char> step_plan_valid(static_cast<size_t>(S));
  std::vector<int> step_plan_groups(static_cast<size_t>(S));
  const size_t step_plan_size = static_cast<size_t>(S) *
    static_cast<size_t>(B);
  std::vector<double> step_lower_score(step_plan_size);
  std::vector<double> step_upper_score(step_plan_size);
  std::vector<double> step_log_weight(step_plan_size);

  for (int s = 0; s < S; ++s) {
    row_alpha[static_cast<size_t>(s)] = regplot_scalar_or_row_real(
      alpha_p, Rf_length(alpha), s, "alpha"
    );
    row_phack[static_cast<size_t>(s)] = regplot_scalar_or_row_int(
      phack_kind_p, Rf_length(phack_kind), s, "phack_kind"
    );
    row_mode[static_cast<size_t>(s)] = regplot_scalar_or_row_int(
      kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"
    );
    row_active_phack[static_cast<size_t>(s)] = regplot_row_has_active_phack(
      row_mode[static_cast<size_t>(s)],
      row_alpha[static_cast<size_t>(s)],
      row_phack[static_cast<size_t>(s)]
    );
  }

  const bool full_support = regplot_selection_mixture_has_full_support(
    omega_p, row_mode, row_active_phack, S, B
  );

  for (int k = 0; k < K; ++k) {
    const double *mean_k = mean_p + static_cast<size_t>(S) * static_cast<size_t>(k);
    const double *sd_k   = sd_p   + static_cast<size_t>(S) * static_cast<size_t>(k);
    const double se_k    = se_p[Rf_length(se) == 1 ? 0 : k];

    for (int s = 0; s < S; ++s) {
      const size_t si = static_cast<size_t>(s);
      step_plan_valid[si]  = 0;
      step_plan_groups[si] = 0;
      const int mode = row_mode[si];
      if (sd_k[s] > 0 && !row_active_phack[si] &&
          (mode == SELKERNEL_STEP || mode == SELKERNEL_STEP_PHACK_POWER)) {
        const size_t offset = si * static_cast<size_t>(B);
        int n_groups = 0;
        step_plan_valid[si] = cpp_selnorm_step_cdf_plan(
          mean_k[s],
          sd_k[s],
          se_k,
          omega_p + s,
          data,
          step_lower_score.data() + offset,
          step_upper_score.data() + offset,
          step_log_weight.data() + offset,
          &n_groups,
          S,
          false
        );
        step_plan_groups[si] = n_groups;
      }
    }

    auto cdf_fun = [mean_k, sd_k, se_k, omega_p, &row_alpha, &row_phack,
                    &row_mode, &row_active_phack, &step_plan_valid,
                    &step_plan_groups, &step_lower_score, &step_upper_score,
                    &step_log_weight,
                    &data, S](double q_val) {
      return regplot_selnorm_mixture_cdf_cached(
        q_val,
        mean_k,
        sd_k,
        se_k,
        omega_p,
        row_alpha,
        row_phack,
        row_mode,
        row_active_phack,
        step_plan_valid,
        step_plan_groups,
        step_lower_score,
        step_upper_score,
        step_log_weight,
        data,
        S
      );
    };

    lower[static_cast<size_t>(k)] = regplot_mixture_quantile(
      p_lower, mean_k, sd_k, S, full_support, cdf_fun,
      k > 0 ? lower[static_cast<size_t>(k - 1)] : NA_REAL,
      k > 0
    );
    upper[static_cast<size_t>(k)] = regplot_mixture_quantile(
      p_upper, mean_k, sd_k, S, full_support, cdf_fun,
      k > 0 ? upper[static_cast<size_t>(k - 1)] : NA_REAL,
      k > 0
    );
  }

  return regplot_interval_result(K, lower, upper);
}
