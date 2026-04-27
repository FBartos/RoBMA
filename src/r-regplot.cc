#include <Rinternals.h>
#include <R_ext/Error.h>

#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <vector>

#include "source/wmnorm.h"

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
  if (TYPEOF(x) != LGLSXP) {
    Rf_error("'%s' must be logical.", name);
  }
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

static int regplot_branch_active_cuts(int branch, const int *mapping_max,
                                      const double *mapping, int n_branches,
                                      int n_map, int n_crit, int n_omega)
{
  if (branch == NA_INTEGER || branch < 1 || branch > n_branches) {
    Rf_error("'bias_indicator' contains an invalid branch index.");
  }

  const int active_cuts = mapping_max[branch - 1];
  if (active_cuts == NA_INTEGER || active_cuts < 0 ||
      active_cuts > n_map || active_cuts > n_crit) {
    Rf_error("'crit_yi_mapping_max' contains an invalid active cutoff count.");
  }

  for (int j = 0; j < active_cuts; ++j) {
    const double raw_index = mapping[j + n_map * (branch - 1)];
    const int index = static_cast<int>(raw_index);
    if (!std::isfinite(raw_index) || raw_index != static_cast<double>(index) ||
        index < 1 || index > n_crit) {
      Rf_error("'crit_yi_mapping' contains an invalid cutoff index.");
    }
    if (index >= n_omega) {
      Rf_error("'omega' must have one more column than the largest mapped cutoff index.");
    }
  }

  return active_cuts;
}

static void regplot_prepare_weighted_log_norm(
  const double *mean, const double *sd, int S,
  const double *omega, int n_omega,
  const int *bias_indicator, const int *is_weightfunction,
  const double *crit_yi, const double *mapping, const int *mapping_max,
  int n_branches, int n_map, int n_crit, bool negative_direction,
  std::vector<double> &log_norm)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);

  for (int s = 0; s < S; ++s) {
    log_norm[static_cast<size_t>(s)] = 0;
    if (is_weightfunction[s] != TRUE || sd[s] < eps_sd) {
      continue;
    }

    const int branch = bias_indicator[s];
    const int active_cuts = regplot_branch_active_cuts(
      branch, mapping_max, mapping, n_branches, n_map, n_crit, n_omega
    );
    const double mean_eval = negative_direction ? -mean[s] : mean[s];

    log_norm[static_cast<size_t>(s)] = cpp_wnorm_mix_log_norm(
      mean_eval, sd[s], crit_yi, omega + s,
      mapping + n_map * (branch - 1), active_cuts, S
    );
  }
}

static double regplot_weighted_mixture_cdf(
  double q, const double *mean, const double *sd, int S,
  const double *omega, int n_omega,
  const int *bias_indicator, const int *is_weightfunction,
  const double *crit_yi, const double *mapping, const int *mapping_max,
  int n_branches, int n_map, int n_crit, bool negative_direction,
  const std::vector<double> &log_norm)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);
  double cdf_sum = 0;

  for (int s = 0; s < S; ++s) {
    double cdf_s;
    if (sd[s] < eps_sd) {
      cdf_s = q >= mean[s] ? 1.0 : 0.0;
    } else if (is_weightfunction[s] != TRUE) {
      cdf_s = pnorm(q, mean[s], sd[s], true, false);
    } else {
      const int branch = bias_indicator[s];
      const int active_cuts = regplot_branch_active_cuts(
        branch, mapping_max, mapping, n_branches, n_map, n_crit, n_omega
      );
      const double q_eval    = negative_direction ? -q       : q;
      const double mean_eval = negative_direction ? -mean[s] : mean[s];
      const bool lower_tail  = !negative_direction;

      cdf_s = cpp_wnorm_mix_cdf_precomputed(
        q_eval, mean_eval, sd[s], crit_yi, omega + s,
        mapping + n_map * (branch - 1), active_cuts,
        log_norm[static_cast<size_t>(s)], lower_tail, false, S
      );
    }

    cdf_sum += regplot_clamp_probability(cdf_s);
  }

  return cdf_sum / static_cast<double>(S);
}

extern "C" SEXP RoBMA_regplot_weighted_mixture_interval(
  SEXP mean, SEXP sd, SEXP probs, SEXP omega, SEXP bias_indicator,
  SEXP is_weightfunction, SEXP crit_yi, SEXP crit_yi_mapping,
  SEXP crit_yi_mapping_max, SEXP effect_direction)
{
  int S, K;
  regplot_validate_mean_sd(mean, sd, &S, &K);
  regplot_validate_probs(probs);
  regplot_check_real(omega, "omega");
  regplot_check_integer(bias_indicator, "bias_indicator");
  regplot_check_logical(is_weightfunction, "is_weightfunction");
  regplot_check_real(crit_yi, "crit_yi");
  regplot_check_real(crit_yi_mapping, "crit_yi_mapping");
  regplot_check_integer(crit_yi_mapping_max, "crit_yi_mapping_max");

  if (TYPEOF(effect_direction) != STRSXP || Rf_length(effect_direction) != 1) {
    Rf_error("'effect_direction' must be a character scalar.");
  }

  const int n_omega    = regplot_matrix_ncol(omega, "omega");
  const int n_crit     = regplot_matrix_nrow(crit_yi, "crit_yi");
  const int K_crit     = regplot_matrix_ncol(crit_yi, "crit_yi");
  const int n_map      = regplot_matrix_nrow(crit_yi_mapping, "crit_yi_mapping");
  const int n_branches = Rf_length(crit_yi_mapping_max);

  if (regplot_matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(bias_indicator) != S) {
    Rf_error("'bias_indicator' must have one value per posterior sample.");
  }
  if (Rf_length(is_weightfunction) != S) {
    Rf_error("'is_weightfunction' must have one value per posterior sample.");
  }
  if (K_crit != 1 && K_crit != K) {
    Rf_error("'crit_yi' must have one column or one column per prediction.");
  }
  if (regplot_matrix_ncol(crit_yi_mapping, "crit_yi_mapping") != n_branches) {
    Rf_error("'crit_yi_mapping' columns must match 'crit_yi_mapping_max'.");
  }

  const bool negative_direction =
    std::strcmp(CHAR(STRING_ELT(effect_direction, 0)), "negative") == 0;

  const double *mean_p   = REAL(mean);
  const double *sd_p     = REAL(sd);
  const double *omega_p  = REAL(omega);
  const double *crit_p   = REAL(crit_yi);
  const double *mapping  = REAL(crit_yi_mapping);
  const int *mapping_max = INTEGER(crit_yi_mapping_max);
  const int *bias_p      = INTEGER(bias_indicator);
  const int *is_w_p      = LOGICAL(is_weightfunction);
  const double p_lower   = REAL(probs)[0];
  const double p_upper   = REAL(probs)[1];

  std::vector<double> lower(static_cast<size_t>(K));
  std::vector<double> upper(static_cast<size_t>(K));
  std::vector<double> log_norm(static_cast<size_t>(S));

  for (int k = 0; k < K; ++k) {
    const double *mean_k = mean_p + static_cast<size_t>(S) * static_cast<size_t>(k);
    const double *sd_k   = sd_p   + static_cast<size_t>(S) * static_cast<size_t>(k);
    const int crit_col   = K_crit == 1 ? 0 : k;
    const double *crit_k = crit_p + static_cast<size_t>(n_crit) * static_cast<size_t>(crit_col);

    regplot_prepare_weighted_log_norm(
      mean_k, sd_k, S, omega_p, n_omega, bias_p, is_w_p, crit_k,
      mapping, mapping_max, n_branches, n_map, n_crit,
      negative_direction, log_norm
    );

    lower[static_cast<size_t>(k)] = regplot_mixture_quantile(
      p_lower, mean_k, sd_k, S,
      [mean_k, sd_k, S, omega_p, n_omega, bias_p, is_w_p, crit_k,
       mapping, mapping_max, n_branches, n_map, n_crit, negative_direction,
       &log_norm](double q) {
        return regplot_weighted_mixture_cdf(
          q, mean_k, sd_k, S, omega_p, n_omega, bias_p, is_w_p,
          crit_k, mapping, mapping_max, n_branches, n_map, n_crit,
          negative_direction, log_norm
        );
      }
    );

    upper[static_cast<size_t>(k)] = regplot_mixture_quantile(
      p_upper, mean_k, sd_k, S,
      [mean_k, sd_k, S, omega_p, n_omega, bias_p, is_w_p, crit_k,
       mapping, mapping_max, n_branches, n_map, n_crit, negative_direction,
       &log_norm](double q) {
        return regplot_weighted_mixture_cdf(
          q, mean_k, sd_k, S, omega_p, n_omega, bias_p, is_w_p,
          crit_k, mapping, mapping_max, n_branches, n_map, n_crit,
          negative_direction, log_norm
        );
      }
    );
  }

  return regplot_interval_result(K, lower, upper);
}
