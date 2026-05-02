#include <Rinternals.h>
#include <R_ext/Error.h>

#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <vector>

#include "source/selnorm.h"

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
                                    SEXP segment_phack_region)
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
}

static bool regplot_row_has_active_phack(int mode, double alpha, int phack_kind)
{
  return (mode == SELKERNEL_PHACK_POWER ||
    mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
}

static double regplot_clamp_probability(double x)
{
  if (x < 0) {
    return 0;
  }
  if (x > 1) {
    return 1;
  }
  return x;
}

static double regplot_quantile_type8(const double *x, int n, double p)
{
  std::vector<double> sorted(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) {
      return NA_REAL;
    }
    sorted[i] = x[i];
  }
  std::sort(sorted.begin(), sorted.end());

  const double h = (static_cast<double>(n) + 1.0 / 3.0) * p + 1.0 / 3.0;
  if (h <= 1) {
    return sorted[0];
  }
  if (h >= n) {
    return sorted[static_cast<size_t>(n - 1)];
  }

  const int j = static_cast<int>(std::floor(h));
  const double g = h - static_cast<double>(j);
  return (1 - g) * sorted[static_cast<size_t>(j - 1)] +
    g * sorted[static_cast<size_t>(j)];
}

static double regplot_normal_mixture_cdf(double q, const double *mean,
                                         const double *sd, int S)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);
  double cdf_sum = 0;

  for (int s = 0; s < S; ++s) {
    double cdf_s;
    if (sd[s] < eps_sd) {
      cdf_s = q >= mean[s] ? 1.0 : 0.0;
    } else {
      cdf_s = pnorm(q, mean[s], sd[s], true, false);
    }
    cdf_sum += regplot_clamp_probability(cdf_s);
  }

  return cdf_sum / static_cast<double>(S);
}

static double regplot_selnorm_mixture_cdf(
  double q,
  const double *mean,
  const double *sd,
  double se,
  const double *omega,
  const double *alpha,
  const int *phack_kind,
  const int *kernel_mode,
  const SelNormKernelData &data,
  int S,
  int alpha_len,
  int phack_kind_len,
  int kernel_mode_len)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);
  double cdf_sum = 0;

  for (int s = 0; s < S; ++s) {
    double cdf_s;
    if (sd[s] < eps_sd) {
      cdf_s = q >= mean[s] ? 1.0 : 0.0;
    } else {
      const double alpha_s = regplot_scalar_or_row_real(alpha, alpha_len, s, "alpha");
      const int phack_s = regplot_scalar_or_row_int(
        phack_kind, phack_kind_len, s, "phack_kind"
      );
      const int mode_s = regplot_scalar_or_row_int(
        kernel_mode, kernel_mode_len, s, "kernel_mode"
      );

      if (regplot_row_has_active_phack(mode_s, alpha_s, phack_s)) {
        cdf_s = NA_REAL;
      } else {
        cdf_s = cpp_selnorm_kernel_cdf(
          q,
          mean[s],
          sd[s],
          se,
          omega + s,
          alpha_s,
          phack_s,
          mode_s,
          data,
          S,
          true,
          false
        );
      }
    }
    cdf_sum += regplot_clamp_probability(cdf_s);
  }

  return cdf_sum / static_cast<double>(S);
}

static void regplot_normal_mixture_cdf_pdf(double q, const double *mean,
                                           const double *sd, int S,
                                           double *cdf, double *pdf)
{
  double cdf_sum = 0;
  double pdf_sum = 0;

  for (int s = 0; s < S; ++s) {
    const double cdf_s = pnorm(q, mean[s], sd[s], true, false);
    const double pdf_s = dnorm(q, mean[s], sd[s], false);
    cdf_sum += regplot_clamp_probability(cdf_s);
    pdf_sum += pdf_s;
  }

  *cdf = cdf_sum / static_cast<double>(S);
  *pdf = pdf_sum / static_cast<double>(S);
}

template <typename CdfFun>
static double regplot_grid_quantile(double p, double lower, double upper,
                                    CdfFun cdf_fun)
{
  const int n_grid = 1000;

  if (lower == upper) {
    return upper;
  }

  for (int i = 0; i < n_grid; ++i) {
    const double q = lower +
      (upper - lower) * static_cast<double>(i) / static_cast<double>(n_grid - 1);
    if (cdf_fun(q) >= p) {
      return q;
    }
  }

  return upper;
}

template <typename CdfFun>
static double regplot_mixture_quantile(double p, const double *mean,
                                       const double *sd, int S,
                                       CdfFun cdf_fun)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);
  bool all_zero_sd = true;
  double lower = INFINITY;
  double upper = -INFINITY;
  double step = 0;
  double max_abs_mean = 1;

  for (int s = 0; s < S; ++s) {
    if (!(sd[s] < eps_sd)) {
      all_zero_sd = false;
    }

    const double spread = std::max(sd[s], eps_sd);
    lower = std::min(lower, mean[s] - 10 * spread);
    upper = std::max(upper, mean[s] + 10 * spread);
    step  = std::max(step, spread);

    if (std::isfinite(mean[s])) {
      max_abs_mean = std::max(max_abs_mean, std::fabs(mean[s]));
    }
  }

  if (all_zero_sd) {
    return regplot_quantile_type8(mean, S, p);
  }

  if (!std::isfinite(lower) || !std::isfinite(upper)) {
    return NA_REAL;
  }
  if (lower >= upper) {
    lower -= 1;
    upper += 1;
  }
  if (!std::isfinite(step) || step <= 0) {
    step = max_abs_mean;
  }

  double lower_value = cdf_fun(lower) - p;
  double upper_value = cdf_fun(upper) - p;

  for (int i = 0; i < 25; ++i) {
    if (lower_value <= 0 && upper_value >= 0) {
      break;
    }
    if (lower_value > 0) {
      lower -= step;
      lower_value = cdf_fun(lower) - p;
    }
    if (upper_value < 0) {
      upper += step;
      upper_value = cdf_fun(upper) - p;
    }
    step *= 2;
  }

  if (lower_value > 0 || upper_value < 0) {
    return regplot_grid_quantile(p, lower, upper, cdf_fun);
  }

  for (int i = 0; i < 100; ++i) {
    const double mid = lower + 0.5 * (upper - lower);
    const double mid_value = cdf_fun(mid) - p;

    if (mid_value >= 0) {
      upper = mid;
      upper_value = mid_value;
    } else {
      lower = mid;
      lower_value = mid_value;
    }

    if (std::fabs(upper - lower) <= 1e-6) {
      break;
    }
    if (mid_value == 0) {
      return mid;
    }
  }

  return lower + 0.5 * (upper - lower);
}

static double regplot_normal_mixture_quantile(double p, const double *mean,
                                              const double *sd, int S)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);
  bool all_zero_sd = true;
  bool any_zero_sd = false;
  double lower = INFINITY;
  double upper = -INFINITY;
  double step = 0;
  double mix_mean = 0;
  double mix_second = 0;

  for (int s = 0; s < S; ++s) {
    const bool zero_sd = sd[s] < eps_sd;
    all_zero_sd = all_zero_sd && zero_sd;
    any_zero_sd = any_zero_sd || zero_sd;

    const double spread = std::max(sd[s], eps_sd);
    lower = std::min(lower, mean[s] - 10 * spread);
    upper = std::max(upper, mean[s] + 10 * spread);
    step  = std::max(step, spread);
    mix_mean   += mean[s];
    mix_second += mean[s] * mean[s] + sd[s] * sd[s];
  }

  if (all_zero_sd) {
    return regplot_quantile_type8(mean, S, p);
  }

  if (any_zero_sd) {
    return regplot_mixture_quantile(
      p, mean, sd, S,
      [mean, sd, S](double q) {
        return regplot_normal_mixture_cdf(q, mean, sd, S);
      }
    );
  }

  if (!std::isfinite(lower) || !std::isfinite(upper)) {
    return NA_REAL;
  }
  if (lower >= upper) {
    lower -= 1;
    upper += 1;
  }
  if (!std::isfinite(step) || step <= 0) {
    step = 1;
  }

  double lower_value = regplot_normal_mixture_cdf(lower, mean, sd, S) - p;
  double upper_value = regplot_normal_mixture_cdf(upper, mean, sd, S) - p;

  for (int i = 0; i < 25; ++i) {
    if (lower_value <= 0 && upper_value >= 0) {
      break;
    }
    if (lower_value > 0) {
      lower -= step;
      lower_value = regplot_normal_mixture_cdf(lower, mean, sd, S) - p;
    }
    if (upper_value < 0) {
      upper += step;
      upper_value = regplot_normal_mixture_cdf(upper, mean, sd, S) - p;
    }
    step *= 2;
  }

  if (lower_value > 0 || upper_value < 0) {
    return regplot_grid_quantile(
      p, lower, upper,
      [mean, sd, S](double q) {
        return regplot_normal_mixture_cdf(q, mean, sd, S);
      }
    );
  }

  mix_mean /= static_cast<double>(S);
  mix_second /= static_cast<double>(S);
  const double mix_var = std::max(mix_second - mix_mean * mix_mean, 0.0);
  double q = mix_mean + std::sqrt(mix_var) * qnorm(p, 0, 1, true, false);

  if (!std::isfinite(q) || q <= lower || q >= upper) {
    q = lower + 0.5 * (upper - lower);
  }

  for (int i = 0; i < 40; ++i) {
    double cdf = 0;
    double pdf = 0;
    regplot_normal_mixture_cdf_pdf(q, mean, sd, S, &cdf, &pdf);
    const double value = cdf - p;

    if (value >= 0) {
      upper = q;
      upper_value = value;
    } else {
      lower = q;
      lower_value = value;
    }

    if (std::fabs(value) <= 1e-9 || std::fabs(upper - lower) <= 1e-6) {
      return q;
    }

    double candidate = NA_REAL;
    if (std::isfinite(pdf) && pdf > 0) {
      candidate = q - value / pdf;
    }
    if (!std::isfinite(candidate) || candidate <= lower || candidate >= upper) {
      candidate = lower + 0.5 * (upper - lower);
    }
    q = candidate;
  }

  return lower + 0.5 * (upper - lower);
}

static SEXP regplot_interval_result(int K, std::vector<double> &lower,
                                    std::vector<double> &upper)
{
  SEXP lower_out = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP upper_out = PROTECT(Rf_allocVector(REALSXP, K));

  for (int k = 0; k < K; ++k) {
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
      p_lower, mean_k, sd_k, S
    );
    upper[static_cast<size_t>(k)] = regplot_normal_mixture_quantile(
      p_upper, mean_k, sd_k, S
    );
  }

  return regplot_interval_result(K, lower, upper);
}

extern "C" SEXP RoBMA_regplot_selnorm_mixture_interval(
  SEXP mean, SEXP sd, SEXP se, SEXP probs, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region)
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

  if (Rf_length(se) != 1 || !(REAL(se)[0] > 0)) {
    Rf_error("'se' must be a positive scalar.");
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
  const double se_p          = REAL(se)[0];
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
    segment_step_bin, segment_phack_region
  );

  std::vector<double> lower(static_cast<size_t>(K));
  std::vector<double> upper(static_cast<size_t>(K));

  for (int k = 0; k < K; ++k) {
    const double *mean_k = mean_p + static_cast<size_t>(S) * static_cast<size_t>(k);
    const double *sd_k   = sd_p   + static_cast<size_t>(S) * static_cast<size_t>(k);

    auto cdf_fun = [mean_k, sd_k, se_p, omega_p, alpha_p, phack_kind_p,
                    kernel_mode_p, &data, S, alpha, phack_kind,
                    kernel_mode](double q_val) {
      return regplot_selnorm_mixture_cdf(
        q_val,
        mean_k,
        sd_k,
        se_p,
        omega_p,
        alpha_p,
        phack_kind_p,
        kernel_mode_p,
        data,
        S,
        Rf_length(alpha),
        Rf_length(phack_kind),
        Rf_length(kernel_mode)
      );
    };

    lower[static_cast<size_t>(k)] = regplot_mixture_quantile(
      p_lower, mean_k, sd_k, S, cdf_fun
    );
    upper[static_cast<size_t>(k)] = regplot_mixture_quantile(
      p_upper, mean_k, sd_k, S, cdf_fun
    );
  }

  return regplot_interval_result(K, lower, upper);
}
