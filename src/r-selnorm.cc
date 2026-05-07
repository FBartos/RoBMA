#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Random.h>

#include <Rmath.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <vector>

#include "selnorm/selnorm.h"

extern "C" double Rf_dnorm4(double, double, double, int);

namespace {

int matrix_nrow(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[0];
}

int matrix_ncol(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[1];
}

void check_real(SEXP x, const char *name)
{
  if (TYPEOF(x) != REALSXP) {
    Rf_error("'%s' must be numeric.", name);
  }
}

void check_integer(SEXP x, const char *name)
{
  if (TYPEOF(x) != INTSXP) {
    Rf_error("'%s' must be integer.", name);
  }
}

int scalar_or_row_int(const int *x, int n, int row, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (row >= n) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
  return x[row];
}

double scalar_or_row_real(const double *x, int n, int row, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (row >= n) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
  return x[row];
}

void validate_matrix_dims(SEXP x, const char *name, int S, int K)
{
  if (matrix_nrow(x, name) != S || matrix_ncol(x, name) != K) {
    Rf_error("'%s' dimensions must match 'mu_num'.", name);
  }
}

void validate_matrix_dims_named(SEXP x, const char *name, int S, int K,
                                const char *ref)
{
  if (matrix_nrow(x, name) != S || matrix_ncol(x, name) != K) {
    Rf_error("'%s' dimensions must match '%s'.", name, ref);
  }
}

void check_logical(SEXP x, const char *name)
{
  if (TYPEOF(x) != LGLSXP || Rf_length(x) != 1) {
    Rf_error("'%s' must be a single logical value.", name);
  }
}

double scalar_or_col_real(const double *x, int n, int col, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (col >= n) {
    Rf_error("'%s' must have length 1 or one value per observation.", name);
  }
  return x[col];
}

void validate_omega_matrix(const double *omega, int S, int B)
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

void set_kernel_data(SelNormKernelData *data, int B, int n_segments,
                     SEXP sign, SEXP q, SEXP z_lower, SEXP z_upper,
                     SEXP phack_z_source, SEXP phack_z_dest,
                     SEXP segment_bounds, SEXP segment_step_bin,
                     SEXP segment_phack_region,
                     bool telescope_probabilities)
{
  data->n_bins               = B;
  data->n_segments           = n_segments;
  data->effect_sign          = INTEGER(sign)[0];
  data->q                    = INTEGER(q)[0];
  data->z_lower              = REAL(z_lower);
  data->z_upper              = REAL(z_upper);
  data->phack_z_source       = REAL(phack_z_source);
  data->phack_z_dest         = REAL(phack_z_dest);
  data->segment_bounds       = REAL(segment_bounds);
  data->segment_step_bin     = INTEGER(segment_step_bin);
  data->segment_phack_region = INTEGER(segment_phack_region);
  data->segment_step_bin_real = 0;
  data->segment_phack_region_real = 0;
  data->trusted_step_partition = selnorm_is_descending_step_partition(
    data->z_lower, data->z_upper, B
  );
  data->telescope_probabilities = data->trusted_step_partition &&
    telescope_probabilities;
}

int z_step_bin(double z, const SelNormKernelData &data)
{
  for (int b = 0; b < data.n_bins; ++b) {
    if (z >= data.z_lower[b] - 1e-12 && z <= data.z_upper[b] + 1e-12) {
      return b + 1;
    }
  }

  return z > data.z_upper[0] ? 1 : data.n_bins;
}

bool row_has_active_phack(int mode, double alpha, int phack_kind)
{
  return (mode == SELKERNEL_PHACK_POWER ||
    mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
}

double clamp_probability_value(double x)
{
  if (x < 0) {
    return 0;
  }
  if (x > 1) {
    return 1;
  }
  return x;
}

double funnel_interval_prob(double lower, double upper, double mean, double sd)
{
  if (!(sd > 0) || lower >= upper) {
    return 0;
  }

  if (std::isinf(lower) && lower < 0) {
    if (std::isinf(upper) && upper > 0) {
      return 1;
    }
    return pnorm(upper, mean, sd, true, false);
  }
  if (std::isinf(upper) && upper > 0) {
    return pnorm(lower, mean, sd, false, false);
  }

  double out = pnorm(upper, mean, sd, true, false) -
    pnorm(lower, mean, sd, true, false);
  if (out < 0 && out > -1e-14) {
    out = 0;
  }
  if (out > 1 && out < 1 + 1e-14) {
    out = 1;
  }
  return out;
}

double funnel_quantile_type8(const std::vector<double> &x, double p)
{
  const int n = static_cast<int>(x.size());
  std::vector<double> sorted(x);

  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(sorted[static_cast<size_t>(i)])) {
      return NA_REAL;
    }
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

struct FunnelSeContext {
  int S;
  double se;
  int direction;
  const double *mu;
  const double *tau;
  const double *PET;
  const double *PEESE;
  const int *is_weightfunction;
  const double *omega;
  const std::vector<double> *row_alpha;
  const std::vector<int> *row_phack;
  const std::vector<int> *row_mode;
  const std::vector<bool> *row_active_phack;
  const SelNormKernelData *data;
  std::vector<double> location;
  std::vector<double> total_sd;
  std::vector<double> mean_z;
  std::vector<double> sd_z;
  std::vector<double> log_norm;
  std::vector<double> normalizer;
  std::vector<double> boundary_tail;
  std::vector<char> zero_sd;
  std::vector<char> weighted;
  std::vector<char> fast_step;
};

double funnel_selected_cdf_zero_se(double q, int s, const FunnelSeContext &ctx)
{
  const SelNormKernelData &data = *ctx.data;
  const double *omega_s = ctx.omega + s;
  const double mean = data.effect_sign * ctx.mu[s];
  const double sd = ctx.total_sd[static_cast<size_t>(s)];
  const double q_signed = data.effect_sign * q;
  double denom = 0;
  double mass = 0;

  for (int b = 0; b < data.n_bins; ++b) {
    const double z_lower = data.z_lower[b];
    const double z_upper = data.z_upper[b];
    double lower;
    double upper;

    if (std::isinf(z_lower) && z_lower < 0 &&
        std::isinf(z_upper) && z_upper > 0) {
      lower = -INFINITY;
      upper =  INFINITY;
    } else if (std::isinf(z_upper) && z_upper > 0) {
      lower = 0;
      upper = INFINITY;
    } else if (std::isinf(z_lower) && z_lower < 0) {
      lower = -INFINITY;
      upper = 0;
    } else {
      continue;
    }

    const double bin_mass = funnel_interval_prob(lower, upper, mean, sd);
    const double weight = omega_s[static_cast<size_t>(ctx.S) * b];
    denom += weight * bin_mass;

    const double selected_mass = data.effect_sign == 1 ?
      funnel_interval_prob(lower, std::min(upper, q_signed), mean, sd) :
      funnel_interval_prob(std::max(lower, q_signed), upper, mean, sd);
    mass += weight * selected_mass;
  }

  if (!(denom > 0)) {
    return NA_REAL;
  }

  return clamp_probability_value(mass / denom);
}

void funnel_prepare_se_context(FunnelSeContext *ctx)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);
  const SelNormKernelData &data = *ctx->data;
  const double se2 = ctx->se * ctx->se;
  const int n_boundaries = std::max(data.n_bins - 1, 0);

  ctx->location.assign(static_cast<size_t>(ctx->S), 0);
  ctx->total_sd.assign(static_cast<size_t>(ctx->S), 0);
  ctx->mean_z.assign(static_cast<size_t>(ctx->S), 0);
  ctx->sd_z.assign(static_cast<size_t>(ctx->S), 0);
  ctx->log_norm.assign(static_cast<size_t>(ctx->S), 0);
  ctx->normalizer.assign(static_cast<size_t>(ctx->S), 0);
  ctx->boundary_tail.assign(
    static_cast<size_t>(ctx->S) * static_cast<size_t>(n_boundaries), 0
  );
  ctx->zero_sd.assign(static_cast<size_t>(ctx->S), 0);
  ctx->weighted.assign(static_cast<size_t>(ctx->S), 0);
  ctx->fast_step.assign(static_cast<size_t>(ctx->S), 0);

  for (int s = 0; s < ctx->S; ++s) {
    const size_t si = static_cast<size_t>(s);
    const double location = ctx->mu[s] +
      ctx->direction * ctx->PET[s] * ctx->se +
      ctx->direction * ctx->PEESE[s] * se2;
    const double total_sd = std::sqrt(se2 + ctx->tau[s] * ctx->tau[s]);

    ctx->location[si] = location;
    ctx->total_sd[si] = total_sd;
    ctx->zero_sd[si]  = total_sd < eps_sd;
    ctx->weighted[si] = ctx->is_weightfunction[s] != 0 && !ctx->zero_sd[si];

    if (ctx->weighted[si] && ctx->se > 0) {
      if ((*ctx->row_active_phack)[si]) {
        ctx->log_norm[si] = std::numeric_limits<double>::quiet_NaN();
      } else {
        const int row_mode = (*ctx->row_mode)[si];
        const bool step_mode = row_mode == SELKERNEL_STEP ||
          row_mode == SELKERNEL_STEP_PHACK_POWER;
        const double mean_z = data.effect_sign * location / ctx->se;
        const double sd_z = total_sd / ctx->se;

        ctx->mean_z[si] = mean_z;
        ctx->sd_z[si] = sd_z;

        if (step_mode && data.trusted_step_partition &&
            data.telescope_probabilities) {
          const double *omega_s = ctx->omega + s;
          double normalizer = omega_s[static_cast<size_t>(ctx->S) *
            static_cast<size_t>(data.n_bins - 1)];
          double norm_abs = std::abs(normalizer);
          bool ok = std::isfinite(normalizer);

          for (int b = 0; ok && b < n_boundaries; ++b) {
            const double w_left = omega_s[static_cast<size_t>(ctx->S) *
              static_cast<size_t>(b)];
            const double w_right = omega_s[static_cast<size_t>(ctx->S) *
              static_cast<size_t>(b + 1)];
            if (w_left == w_right) {
              continue;
            }

            const double tail_boundary = pnorm(
              data.z_lower[b], mean_z, sd_z, false, false
            );
            if (!std::isfinite(tail_boundary)) {
              ok = false;
              break;
            }

            ctx->boundary_tail[
              si * static_cast<size_t>(n_boundaries) + static_cast<size_t>(b)
            ] = tail_boundary;

            const double contribution = (w_left - w_right) * tail_boundary;
            normalizer += contribution;
            norm_abs += std::abs(contribution);
          }

          if (ok && std::isfinite(normalizer) && normalizer > 0 &&
              normalizer > std::numeric_limits<double>::min() &&
              normalizer > 1e-8 * norm_abs) {
            ctx->normalizer[si] = normalizer;
            ctx->log_norm[si] = std::log(normalizer);
            ctx->fast_step[si] = 1;
            continue;
          }
        }

        ctx->log_norm[si] = cpp_selnorm_kernel_log_norm(
          location,
          total_sd,
          ctx->se,
          ctx->omega + s,
          (*ctx->row_alpha)[si],
          (*ctx->row_phack)[si],
          (*ctx->row_mode)[si],
          data,
          ctx->S,
          false
        );
      }
    }
  }
}

double funnel_step_selected_cdf_precomputed(double q, int s,
                                            const FunnelSeContext &ctx)
{
  const SelNormKernelData &data = *ctx.data;
  const size_t si = static_cast<size_t>(s);
  const int n_boundaries = std::max(data.n_bins - 1, 0);
  const double q_z = data.effect_sign * q / ctx.se;
  const double tail_q = pnorm(q_z, ctx.mean_z[si], ctx.sd_z[si], false, false);
  const double cdf_q = pnorm(q_z, ctx.mean_z[si], ctx.sd_z[si], true, false);

  if (!std::isfinite(tail_q) || !std::isfinite(cdf_q)) {
    return NA_REAL;
  }

  const double *omega_s = ctx.omega + s;
  const bool signed_lower_tail = data.effect_sign == 1;
  const double normalizer = ctx.normalizer[si];
  double requested = signed_lower_tail ?
    omega_s[static_cast<size_t>(ctx.S) * static_cast<size_t>(data.n_bins - 1)] *
      cdf_q :
    omega_s[static_cast<size_t>(ctx.S) * static_cast<size_t>(data.n_bins - 1)] *
      tail_q;
  double req_abs = std::abs(requested);

  for (int b = 0; b < n_boundaries; ++b) {
    const double w_left = omega_s[static_cast<size_t>(ctx.S) *
      static_cast<size_t>(b)];
    const double w_right = omega_s[static_cast<size_t>(ctx.S) *
      static_cast<size_t>(b + 1)];
    if (w_left == w_right) {
      continue;
    }

    const double diff = w_left - w_right;
    const double tail_boundary = ctx.boundary_tail[
      si * static_cast<size_t>(n_boundaries) + static_cast<size_t>(b)
    ];
    double req_contrib = 0;

    if (signed_lower_tail) {
      if (q_z > data.z_lower[b]) {
        double interval = tail_boundary - tail_q;
        const double scale = std::max(std::abs(tail_boundary), std::abs(tail_q));
        const double tol = 1e-14 * std::max(scale, 1.0);
        if (interval < 0 && interval > -tol) {
          interval = 0;
        } else if (interval < 0) {
          return NA_REAL;
        }
        req_contrib = diff * interval;
      }
    } else {
      req_contrib = diff * (q_z >= data.z_lower[b] ? tail_q : tail_boundary);
    }

    requested += req_contrib;
    req_abs += std::abs(req_contrib);
  }

  const double req_tol = 1e-10 *
    std::max(std::max(req_abs, normalizer), 1.0);
  if (requested < 0 && requested > -req_tol) {
    requested = 0;
  } else if (requested < 0) {
    return NA_REAL;
  }
  if (requested > normalizer && requested < normalizer + req_tol) {
    requested = normalizer;
  } else if (requested > normalizer) {
    return NA_REAL;
  }

  return clamp_probability_value(requested / normalizer);
}

double funnel_model_averaged_cdf_native(double q, const FunnelSeContext &ctx)
{
  double cdf_sum = 0;

  for (int s = 0; s < ctx.S; ++s) {
    const size_t si = static_cast<size_t>(s);
    double cdf_s;

    if (ctx.zero_sd[si]) {
      cdf_s = q >= ctx.location[si] ? 1.0 : 0.0;
    } else if (ctx.weighted[si]) {
      if (ctx.se <= 0) {
        cdf_s = funnel_selected_cdf_zero_se(q, s, ctx);
      } else if ((*ctx.row_active_phack)[si]) {
        cdf_s = NA_REAL;
      } else if (ctx.fast_step[si]) {
        cdf_s = funnel_step_selected_cdf_precomputed(q, s, ctx);
        if (!std::isfinite(cdf_s)) {
          cdf_s = cpp_selnorm_kernel_cdf_with_log_norm(
            q,
            ctx.location[si],
            ctx.total_sd[si],
            ctx.se,
            ctx.omega + s,
            (*ctx.row_alpha)[si],
            (*ctx.row_phack)[si],
            (*ctx.row_mode)[si],
            *ctx.data,
            ctx.log_norm[si],
            ctx.S,
            true,
            false
          );
        }
      } else {
        cdf_s = cpp_selnorm_kernel_cdf_with_log_norm(
          q,
          ctx.location[si],
          ctx.total_sd[si],
          ctx.se,
          ctx.omega + s,
          (*ctx.row_alpha)[si],
          (*ctx.row_phack)[si],
          (*ctx.row_mode)[si],
          *ctx.data,
          ctx.log_norm[si],
          ctx.S,
          true,
          false
        );
      }
    } else {
      cdf_s = pnorm(q, ctx.location[si], ctx.total_sd[si], true, false);
    }

    cdf_sum += clamp_probability_value(cdf_s);
  }

  return cdf_sum / static_cast<double>(ctx.S);
}

double funnel_grid_quantile_native(double p, double lower, double upper,
                                   const FunnelSeContext &ctx)
{
  const int n_grid = 1000;

  if (lower == upper) {
    return upper;
  }

  for (int i = 0; i < n_grid; ++i) {
    const double q = lower +
      (upper - lower) * static_cast<double>(i) / static_cast<double>(n_grid - 1);
    const double cdf = funnel_model_averaged_cdf_native(q, ctx);
    if (std::isfinite(cdf) && cdf >= p) {
      return q;
    }
  }

  return upper;
}

double funnel_model_averaged_quantile_native(double p,
                                             const FunnelSeContext &ctx,
                                             double previous,
                                             bool use_previous)
{
  const double eps_sd = std::sqrt(DBL_EPSILON);
  bool all_zero_sd = true;
  double lower = INFINITY;
  double upper = -INFINITY;
  double step = 0;
  double max_abs_location = 1;

  for (int s = 0; s < ctx.S; ++s) {
    const size_t si = static_cast<size_t>(s);
    all_zero_sd = all_zero_sd && ctx.zero_sd[si];
    const double spread = std::max(ctx.total_sd[si], eps_sd);
    lower = std::min(lower, ctx.location[si] - 10 * spread);
    upper = std::max(upper, ctx.location[si] + 10 * spread);
    step  = std::max(step, spread);

    if (std::isfinite(ctx.location[si])) {
      max_abs_location = std::max(max_abs_location, std::fabs(ctx.location[si]));
    }
  }

  if (all_zero_sd) {
    return funnel_quantile_type8(ctx.location, p);
  }

  if (!std::isfinite(lower) || !std::isfinite(upper)) {
    return NA_REAL;
  }
  if (lower >= upper) {
    lower -= 1;
    upper += 1;
  }
  if (!std::isfinite(step) || step <= 0) {
    step = max_abs_location;
  }

  const double global_lower = lower;
  const double global_upper = upper;
  const double global_step  = step;
  bool using_previous = use_previous && std::isfinite(previous);

  if (using_previous) {
    step = std::max(1e-3, 0.05 * global_step);
    lower = previous - step;
    upper = previous + step;
  }

  double lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
  double upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;

  if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    if (!using_previous) {
      return funnel_grid_quantile_native(p, lower, upper, ctx);
    }
    using_previous = false;
    lower = global_lower;
    upper = global_upper;
    step = global_step;
    lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
    upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;
    if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
      return funnel_grid_quantile_native(p, lower, upper, ctx);
    }
  }

  for (int i = 0; i < 25; ++i) {
    if (lower_value <= 0 && upper_value >= 0) {
      break;
    }
    if (lower_value > 0) {
      lower -= step;
      lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
    }
    if (upper_value < 0) {
      upper += step;
      upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;
    }
    step *= 2;
  }

  if (lower_value > 0 || upper_value < 0 ||
      !std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    if (using_previous) {
      lower = global_lower;
      upper = global_upper;
      step = global_step;
      lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
      upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;

      for (int i = 0; i < 25; ++i) {
        if (lower_value <= 0 && upper_value >= 0) {
          break;
        }
        if (lower_value > 0) {
          lower -= step;
          lower_value = funnel_model_averaged_cdf_native(lower, ctx) - p;
        }
        if (upper_value < 0) {
          upper += step;
          upper_value = funnel_model_averaged_cdf_native(upper, ctx) - p;
        }
        step *= 2;
      }
    }
  }

  if (lower_value > 0 || upper_value < 0 ||
      !std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    return funnel_grid_quantile_native(p, lower, upper, ctx);
  }

  for (int i = 0; i < 100; ++i) {
    const double mid = lower + 0.5 * (upper - lower);
    const double mid_value = funnel_model_averaged_cdf_native(mid, ctx) - p;

    if (!std::isfinite(mid_value)) {
      return funnel_grid_quantile_native(p, lower, upper, ctx);
    }

    if (mid_value >= 0) {
      upper = mid;
      upper_value = mid_value;
    } else {
      lower = mid;
      lower_value = mid_value;
    }

    if (std::fabs(upper - lower) <= 1e-6 || mid_value == 0) {
      break;
    }
  }

  return lower + 0.5 * (upper - lower);
}

}

extern "C" SEXP RoBMA_selnorm_kernel_loglik_matrix(
  SEXP yi, SEXP mu_num, SEXP sigma_num, SEXP mu_norm, SEXP sigma_norm,
  SEXP sei, SEXP weights, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
  SEXP q, SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(yi, "yi");
  check_real(mu_num, "mu_num");
  check_real(sigma_num, "sigma_num");
  check_real(mu_norm, "mu_norm");
  check_real(sigma_norm, "sigma_norm");
  check_real(sei, "sei");
  check_real(weights, "weights");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(obs_bin, "obs_bin");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = matrix_nrow(mu_num, "mu_num");
  const int K = matrix_ncol(mu_num, "mu_num");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims(sigma_num, "sigma_num", S, K);
  validate_matrix_dims(mu_norm, "mu_norm", S, K);
  validate_matrix_dims(sigma_norm, "sigma_norm", S, K);

  if (Rf_length(yi) != K || Rf_length(sei) != K || Rf_length(weights) != K ||
      Rf_length(obs_bin) != K) {
    Rf_error("'yi', 'sei', 'weights', and 'obs_bin' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *yi_p         = REAL(yi);
  const double *mu_num_p     = REAL(mu_num);
  const double *sigma_num_p  = REAL(sigma_num);
  const double *mu_norm_p    = REAL(mu_norm);
  const double *sigma_norm_p = REAL(sigma_norm);
  const double *sei_p        = REAL(sei);
  const double *weights_p    = REAL(weights);
  const double *omega_p      = REAL(omega);
  const double *alpha_p      = REAL(alpha);
  const int *phack_kind_p    = INTEGER(phack_kind);
  const int *kernel_mode_p   = INTEGER(kernel_mode);
  const int *obs_bin_p       = INTEGER(obs_bin);
  double *out_p              = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_lpdf(
        yi_p[k],
        mu_num_p[cell],
        sigma_num_p[cell],
        mu_norm_p[cell],
        sigma_norm_p[cell],
        sei_p[k],
        weights_p[k],
        omega_p + s,
        obs_bin_p[k],
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        S,
        false
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_log_norm_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *sei_p       = REAL(sei);
  const double *omega_p     = REAL(omega);
  const double *alpha_p     = REAL(alpha);
  const int *phack_kind_p   = INTEGER(phack_kind);
  const int *kernel_mode_p  = INTEGER(kernel_mode);
  double *out_p             = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_log_norm(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        S,
        false
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_cdf_matrix(
  SEXP x, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP lower_tail, SEXP telescope_probabilities)
{
  check_real(x, "q");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(lower_tail, "lower.tail");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if ((Rf_length(x) != 1 && Rf_length(x) != K) || Rf_length(sei) != K) {
    Rf_error("'q' must have length 1 or K and 'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *x_p          = REAL(x);
  const double *mean_p       = REAL(mean);
  const double *sd_p         = REAL(sd);
  const double *sei_p        = REAL(sei);
  const double *omega_p      = REAL(omega);
  const double *alpha_p      = REAL(alpha);
  const int *phack_kind_p    = INTEGER(phack_kind);
  const int *kernel_mode_p   = INTEGER(kernel_mode);
  const bool lower_tail_bool = LOGICAL(lower_tail)[0] == TRUE;
  double *out_p              = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_cdf(
        scalar_or_col_real(x_p, Rf_length(x), k, "q"),
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        S,
        lower_tail_bool,
        false
      );
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_moments_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out_mean = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  SEXP out_second = PROTECT(Rf_allocMatrix(REALSXP, S, K));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));

  SET_STRING_ELT(names, 0, Rf_mkChar("mean"));
  SET_STRING_ELT(names, 1, Rf_mkChar("second"));
  SET_VECTOR_ELT(out, 0, out_mean);
  SET_VECTOR_ELT(out, 1, out_second);
  Rf_setAttrib(out, R_NamesSymbol, names);

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  const double *omega_p    = REAL(omega);
  const double *alpha_p    = REAL(alpha);
  const int *phack_kind_p  = INTEGER(phack_kind);
  const int *kernel_mode_p = INTEGER(kernel_mode);
  double *mean_out_p       = REAL(out_mean);
  double *second_out_p     = REAL(out_second);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      cpp_selnorm_kernel_moments(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        mean_out_p + cell,
        second_out_p + cell,
        S,
        false
      );
    }
  }

  UNPROTECT(4);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_rng_matrix(
  SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");

  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
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

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, K));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  const double *omega_p    = REAL(omega);
  const double *alpha_p    = REAL(alpha);
  const int *phack_kind_p  = INTEGER(phack_kind);
  const int *kernel_mode_p = INTEGER(kernel_mode);
  double *out_p            = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  GetRNGstate();
  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      out_p[cell] = cpp_selnorm_kernel_rng(
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        unif_rand(),
        unif_rand(),
        scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha"),
        scalar_or_row_int(phack_kind_p, Rf_length(phack_kind), s, "phack_kind"),
        scalar_or_row_int(kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"),
        data,
        S,
        false
      );
    }
  }
  PutRNGstate();

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_kernel_weighted_summary(
  SEXP yi, SEXP mean, SEXP sd, SEXP sei, SEXP psis_weights, SEXP omega,
  SEXP alpha, SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  check_real(yi, "yi");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(psis_weights, "psis_weights");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  validate_matrix_dims_named(psis_weights, "psis_weights", S, K, "mean");

  if (Rf_length(yi) != K || Rf_length(sei) != K) {
    Rf_error("'yi' and 'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
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

  SEXP out_cdf = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out_mean = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out_second = PROTECT(Rf_allocVector(REALSXP, K));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));

  SET_STRING_ELT(names, 0, Rf_mkChar("cdf"));
  SET_STRING_ELT(names, 1, Rf_mkChar("mean"));
  SET_STRING_ELT(names, 2, Rf_mkChar("second"));
  SET_VECTOR_ELT(out, 0, out_cdf);
  SET_VECTOR_ELT(out, 1, out_mean);
  SET_VECTOR_ELT(out, 2, out_second);
  Rf_setAttrib(out, R_NamesSymbol, names);

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *yi_p        = REAL(yi);
  const double *mean_p      = REAL(mean);
  const double *sd_p        = REAL(sd);
  const double *sei_p       = REAL(sei);
  const double *weights_p   = REAL(psis_weights);
  const double *omega_p     = REAL(omega);
  const double *alpha_p     = REAL(alpha);
  const int *phack_kind_p   = INTEGER(phack_kind);
  const int *kernel_mode_p  = INTEGER(kernel_mode);
  double *cdf_out_p         = REAL(out_cdf);
  double *mean_out_p        = REAL(out_mean);
  double *second_out_p      = REAL(out_second);

  validate_omega_matrix(omega_p, S, B);

  for (int k = 0; k < K; ++k) {
    double cdf_sum = 0;
    double mean_sum = 0;
    double second_sum = 0;

    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;
      const double w = weights_p[cell];
      const double alpha_s = scalar_or_row_real(alpha_p, Rf_length(alpha), s, "alpha");
      const int phack_s = scalar_or_row_int(
        phack_kind_p, Rf_length(phack_kind), s, "phack_kind"
      );
      const int mode_s = scalar_or_row_int(
        kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"
      );
      double cdf = NA_REAL;
      double moment_mean = NA_REAL;
      double moment_second = NA_REAL;

      cpp_selnorm_kernel_summary(
        yi_p[k],
        mean_p[cell],
        sd_p[cell],
        sei_p[k],
        omega_p + s,
        alpha_s,
        phack_s,
        mode_s,
        data,
        &cdf,
        &moment_mean,
        &moment_second,
        S,
        true,
        false
      );

      cdf_sum += w * cdf;
      mean_sum += w * moment_mean;
      second_sum += w * moment_second;
    }

    cdf_out_p[k] = cdf_sum;
    mean_out_p[k] = mean_sum;
    second_out_p[k] = second_sum;
  }

  UNPROTECT(5);
  return out;
}

extern "C" SEXP RoBMA_funnel_model_averaged_quantiles(
  SEXP se_sequence, SEXP mu, SEXP tau, SEXP PET, SEXP PEESE,
  SEXP is_weightfunction, SEXP omega, SEXP alpha, SEXP phack_kind,
  SEXP kernel_mode, SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region, SEXP direction,
  SEXP telescope_probabilities)
{
  check_real(se_sequence, "se_sequence");
  check_real(mu, "mu");
  check_real(tau, "tau");
  check_real(PET, "PET");
  check_real(PEESE, "PEESE");
  check_integer(is_weightfunction, "is_weightfunction");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_integer(direction, "direction");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = Rf_length(mu);
  const int L = Rf_length(se_sequence);
  const int B = matrix_ncol(omega, "omega");

  if (S < 1) {
    Rf_error("'mu' must contain at least one posterior sample.");
  }
  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }
  if (Rf_length(tau) != S || Rf_length(PET) != S ||
      Rf_length(PEESE) != S || Rf_length(is_weightfunction) != S) {
    Rf_error("'mu', 'tau', 'PET', 'PEESE', and 'is_weightfunction' must have equal lengths.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
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
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }
  if (Rf_length(direction) != 1 ||
      (INTEGER(direction)[0] != 1 && INTEGER(direction)[0] != -1)) {
    Rf_error("'direction' must be 1 or -1.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }

  SEXP out_lower = PROTECT(Rf_allocVector(REALSXP, L));
  SEXP out_upper = PROTECT(Rf_allocVector(REALSXP, L));
  SEXP out_mid = PROTECT(Rf_allocVector(REALSXP, L));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));

  SET_STRING_ELT(names, 0, Rf_mkChar("lower"));
  SET_STRING_ELT(names, 1, Rf_mkChar("upper"));
  SET_STRING_ELT(names, 2, Rf_mkChar("mid"));
  SET_VECTOR_ELT(out, 0, out_lower);
  SET_VECTOR_ELT(out, 1, out_upper);
  SET_VECTOR_ELT(out, 2, out_mid);
  Rf_setAttrib(out, R_NamesSymbol, names);

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *se_p      = REAL(se_sequence);
  const double *mu_p      = REAL(mu);
  const double *tau_p     = REAL(tau);
  const double *PET_p     = REAL(PET);
  const double *PEESE_p   = REAL(PEESE);
  const int *weight_p     = INTEGER(is_weightfunction);
  const double *omega_p   = REAL(omega);
  const double *alpha_p   = REAL(alpha);
  const int *phack_p      = INTEGER(phack_kind);
  const int *mode_p       = INTEGER(kernel_mode);
  const int direction_p   = INTEGER(direction)[0];
  double *lower_p         = REAL(out_lower);
  double *upper_p         = REAL(out_upper);
  double *mid_p           = REAL(out_mid);

  validate_omega_matrix(omega_p, S, B);

  std::vector<double> row_alpha(static_cast<size_t>(S));
  std::vector<int> row_phack(static_cast<size_t>(S));
  std::vector<int> row_mode(static_cast<size_t>(S));
  std::vector<bool> row_active_phack(static_cast<size_t>(S));

  for (int s = 0; s < S; ++s) {
    row_alpha[static_cast<size_t>(s)] = scalar_or_row_real(
      alpha_p, Rf_length(alpha), s, "alpha"
    );
    row_phack[static_cast<size_t>(s)] = scalar_or_row_int(
      phack_p, Rf_length(phack_kind), s, "phack_kind"
    );
    row_mode[static_cast<size_t>(s)] = scalar_or_row_int(
      mode_p, Rf_length(kernel_mode), s, "kernel_mode"
    );
    row_active_phack[static_cast<size_t>(s)] = row_has_active_phack(
      row_mode[static_cast<size_t>(s)],
      row_alpha[static_cast<size_t>(s)],
      row_phack[static_cast<size_t>(s)]
    );
  }

  FunnelSeContext ctx;
  ctx.S = S;
  ctx.direction = direction_p;
  ctx.mu = mu_p;
  ctx.tau = tau_p;
  ctx.PET = PET_p;
  ctx.PEESE = PEESE_p;
  ctx.is_weightfunction = weight_p;
  ctx.omega = omega_p;
  ctx.row_alpha = &row_alpha;
  ctx.row_phack = &row_phack;
  ctx.row_mode = &row_mode;
  ctx.row_active_phack = &row_active_phack;
  ctx.data = &data;

  for (int ell = 0; ell < L; ++ell) {
    ctx.se = se_p[ell];
    funnel_prepare_se_context(&ctx);

    const bool has_previous = ell > 0;
    lower_p[ell] = funnel_model_averaged_quantile_native(
      0.025, ctx, has_previous ? lower_p[ell - 1] : NA_REAL, has_previous
    );
    upper_p[ell] = funnel_model_averaged_quantile_native(
      0.975, ctx, has_previous ? upper_p[ell - 1] : NA_REAL, has_previous
    );
    mid_p[ell] = funnel_model_averaged_quantile_native(
      0.5, ctx, has_previous ? mid_p[ell - 1] : NA_REAL, has_previous
    );
  }

  UNPROTECT(5);
  return out;
}

extern "C" SEXP RoBMA_zcurve_normal_density_matrix(
  SEXP z_sequence, SEXP mean, SEXP sd, SEXP sei)
{
  check_real(z_sequence, "z_sequence");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int L = Rf_length(z_sequence);

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, L));

  const double *z_p    = REAL(z_sequence);
  const double *mean_p = REAL(mean);
  const double *sd_p   = REAL(sd);
  const double *sei_p  = REAL(sei);
  double *out_p        = REAL(out);

  for (int ell = 0; ell < L; ++ell) {
    for (int s = 0; s < S; ++s) {
      double dens_sum = 0;
      for (int k = 0; k < K; ++k) {
        const int cell = s + S * k;
        dens_sum += Rf_dnorm4(
          z_p[ell] * sei_p[k], mean_p[cell], sd_p[cell], false
        ) * sei_p[k];
      }
      out_p[s + S * ell] = dens_sum / K;
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP RoBMA_selnorm_zcurve_threshold_summary(
  SEXP z_threshold, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP extrapolate, SEXP telescope_probabilities)
{
  check_real(z_threshold, "z_threshold");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(extrapolate, "extrapolate");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int B = matrix_ncol(omega, "omega");

  if (Rf_length(z_threshold) != 1) {
    Rf_error("'z_threshold' must have length 1.");
  }
  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
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

  SEXP out_edr = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP out_weights = PROTECT(Rf_allocVector(REALSXP, S));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));

  SET_STRING_ELT(names, 0, Rf_mkChar("EDR"));
  SET_STRING_ELT(names, 1, Rf_mkChar("weights"));
  SET_VECTOR_ELT(out, 0, out_edr);
  SET_VECTOR_ELT(out, 1, out_weights);
  Rf_setAttrib(out, R_NamesSymbol, names);

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double threshold = REAL(z_threshold)[0];
  const double *mean_p   = REAL(mean);
  const double *sd_p     = REAL(sd);
  const double *sei_p    = REAL(sei);
  const double *omega_p  = REAL(omega);
  const double *alpha_p  = REAL(alpha);
  const int *phack_p     = INTEGER(phack_kind);
  const int *mode_p      = INTEGER(kernel_mode);
  const bool extrapolate_b = LOGICAL(extrapolate)[0] == TRUE;
  double *edr_p          = REAL(out_edr);
  double *weights_p      = REAL(out_weights);

  validate_omega_matrix(omega_p, S, B);

  std::vector<double> row_alpha(S);
  std::vector<int> row_phack(S);
  std::vector<int> row_mode(S);
  std::vector<bool> row_active_phack(S);
  std::vector<bool> row_use_step(S);

  for (int s = 0; s < S; ++s) {
    row_alpha[static_cast<size_t>(s)] = scalar_or_row_real(
      alpha_p, Rf_length(alpha), s, "alpha"
    );
    row_phack[static_cast<size_t>(s)] = scalar_or_row_int(
      phack_p, Rf_length(phack_kind), s, "phack_kind"
    );
    row_mode[static_cast<size_t>(s)] = scalar_or_row_int(
      mode_p, Rf_length(kernel_mode), s, "kernel_mode"
    );
    row_active_phack[static_cast<size_t>(s)] = row_has_active_phack(
      row_mode[static_cast<size_t>(s)],
      row_alpha[static_cast<size_t>(s)],
      row_phack[static_cast<size_t>(s)]
    );
    row_use_step[static_cast<size_t>(s)] =
      row_mode[static_cast<size_t>(s)] == SELKERNEL_STEP ||
      row_mode[static_cast<size_t>(s)] == SELKERNEL_STEP_PHACK_POWER;
  }

  for (int s = 0; s < S; ++s) {
    double threshold_sum = 0;
    double weight_sum = 0;

    for (int k = 0; k < K; ++k) {
      const int cell = s + S * k;
      double inverse_weight = 1;
      double threshold_prob;

      if (row_active_phack[static_cast<size_t>(s)]) {
        threshold_prob = NA_REAL;
        inverse_weight = NA_REAL;
      } else if (row_use_step[static_cast<size_t>(s)] && extrapolate_b) {
        double unused_inverse;
        threshold_prob = cpp_selnorm_kernel_threshold(
          threshold,
          mean_p[cell],
          sd_p[cell],
          sei_p[k],
          omega_p + s,
          row_alpha[static_cast<size_t>(s)],
          row_phack[static_cast<size_t>(s)],
          SELKERNEL_NORMAL,
          data,
          &unused_inverse,
          S,
          false
        );

        const double log_norm = cpp_selnorm_kernel_log_norm(
          mean_p[cell],
          sd_p[cell],
          sei_p[k],
          omega_p + s,
          row_alpha[static_cast<size_t>(s)],
          row_phack[static_cast<size_t>(s)],
          row_mode[static_cast<size_t>(s)],
          data,
          S,
          false
        );
        inverse_weight = std::isfinite(log_norm) ?
          std::exp(-log_norm) : NA_REAL;
      } else {
        threshold_prob = cpp_selnorm_kernel_threshold(
          threshold,
          mean_p[cell],
          sd_p[cell],
          sei_p[k],
          omega_p + s,
          row_alpha[static_cast<size_t>(s)],
          row_phack[static_cast<size_t>(s)],
          row_mode[static_cast<size_t>(s)],
          data,
          &inverse_weight,
          S,
          false
        );
      }

      if (extrapolate_b) {
        threshold_sum += threshold_prob * inverse_weight;
      } else {
        threshold_sum += threshold_prob;
      }
      weight_sum += inverse_weight;
    }

    edr_p[s] = extrapolate_b ? threshold_sum / weight_sum : threshold_sum / K;
    weights_p[s] = weight_sum / K;
  }

  UNPROTECT(4);
  return out;
}

extern "C" SEXP RoBMA_selnorm_zcurve_density_matrix(
  SEXP z_sequence, SEXP mean, SEXP sd, SEXP sei, SEXP omega, SEXP alpha,
  SEXP phack_kind, SEXP kernel_mode, SEXP z_lower, SEXP z_upper,
  SEXP sign, SEXP q, SEXP phack_z_source, SEXP phack_z_dest,
  SEXP segment_bounds, SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP extrapolate, SEXP telescope_probabilities)
{
  check_real(z_sequence, "z_sequence");
  check_real(mean, "mean");
  check_real(sd, "sd");
  check_real(sei, "sei");
  check_real(omega, "omega");
  check_real(alpha, "alpha");
  check_integer(phack_kind, "phack_kind");
  check_integer(kernel_mode, "kernel_mode");
  check_real(z_lower, "z_lower");
  check_real(z_upper, "z_upper");
  check_integer(sign, "sign");
  check_integer(q, "q");
  check_real(phack_z_source, "phack_z_source");
  check_real(phack_z_dest, "phack_z_dest");
  check_real(segment_bounds, "segment_bounds");
  check_integer(segment_step_bin, "segment_step_bin");
  check_integer(segment_phack_region, "segment_phack_region");
  check_logical(extrapolate, "extrapolate");
  check_logical(telescope_probabilities, "telescope_probabilities");

  const int S = matrix_nrow(mean, "mean");
  const int K = matrix_ncol(mean, "mean");
  const int L = Rf_length(z_sequence);
  const int B = matrix_ncol(omega, "omega");

  if (B < 1) {
    Rf_error("'omega' must have at least one column.");
  }

  validate_matrix_dims_named(sd, "sd", S, K, "mean");
  if (Rf_length(sei) != K) {
    Rf_error("'sei' must have one value per observation.");
  }
  if (matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior sample.");
  }
  if (Rf_length(z_lower) != B || Rf_length(z_upper) != B) {
    Rf_error("'z_lower' and 'z_upper' must match the number of omega columns.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }
  if (Rf_length(sign) != 1 || (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
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

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, S, L));

  SelNormKernelData data;
  set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region,
    LOGICAL(telescope_probabilities)[0] == TRUE
  );

  const double *z_p        = REAL(z_sequence);
  const double *mean_p     = REAL(mean);
  const double *sd_p       = REAL(sd);
  const double *sei_p      = REAL(sei);
  const double *omega_p    = REAL(omega);
  const double *alpha_p    = REAL(alpha);
  const int *phack_kind_p  = INTEGER(phack_kind);
  const int *kernel_mode_p = INTEGER(kernel_mode);
  const bool extrapolate_b = LOGICAL(extrapolate)[0] == TRUE;
  double *out_p            = REAL(out);

  validate_omega_matrix(omega_p, S, B);

  std::vector<double> row_alpha(S);
  std::vector<int> row_phack(S);
  std::vector<int> row_mode(S);
  std::vector<bool> row_active_phack(S);
  std::vector<bool> row_use_step(S);
  bool any_step = false;

  for (int s = 0; s < S; ++s) {
    row_alpha[static_cast<size_t>(s)] = scalar_or_row_real(
      alpha_p, Rf_length(alpha), s, "alpha"
    );
    row_phack[static_cast<size_t>(s)] = scalar_or_row_int(
      phack_kind_p, Rf_length(phack_kind), s, "phack_kind"
    );
    row_mode[static_cast<size_t>(s)] = scalar_or_row_int(
      kernel_mode_p, Rf_length(kernel_mode), s, "kernel_mode"
    );
    row_active_phack[static_cast<size_t>(s)] = row_has_active_phack(
      row_mode[static_cast<size_t>(s)],
      row_alpha[static_cast<size_t>(s)],
      row_phack[static_cast<size_t>(s)]
    );
    row_use_step[static_cast<size_t>(s)] =
      row_mode[static_cast<size_t>(s)] == SELKERNEL_STEP ||
      row_mode[static_cast<size_t>(s)] == SELKERNEL_STEP_PHACK_POWER;
    any_step = any_step || row_use_step[static_cast<size_t>(s)];
  }

  std::vector<double> log_norm(S * K, 0);
  for (int k = 0; k < K; ++k) {
    for (int s = 0; s < S; ++s) {
      const int cell = s + S * k;

      if (row_active_phack[static_cast<size_t>(s)]) {
        log_norm[cell] = std::numeric_limits<double>::quiet_NaN();
      } else if (row_use_step[static_cast<size_t>(s)]) {
        log_norm[cell] = cpp_selnorm_kernel_log_norm(
          mean_p[cell],
          sd_p[cell],
          sei_p[k],
          omega_p + s,
          row_alpha[static_cast<size_t>(s)],
          row_phack[static_cast<size_t>(s)],
          row_mode[static_cast<size_t>(s)],
          data,
          S,
          false
        );
      }
    }
  }

  std::vector<int> z_bin(L, 1);
  if (any_step && !extrapolate_b) {
    for (int ell = 0; ell < L; ++ell) {
      z_bin[static_cast<size_t>(ell)] = z_step_bin(
        data.effect_sign * z_p[ell], data
      );
    }
  }

  for (int ell = 0; ell < L; ++ell) {
    for (int s = 0; s < S; ++s) {
      if (row_active_phack[static_cast<size_t>(s)]) {
        out_p[s + S * ell] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }

      double dens_sum = 0;
      if (row_use_step[static_cast<size_t>(s)]) {
        double log_observed_weight = 0;

        if (!extrapolate_b) {
          const int bin = z_bin[static_cast<size_t>(ell)];
          const double w = omega_p[s + S * (bin - 1)];
          if (!(w > 0)) {
            out_p[s + S * ell] = 0;
            continue;
          }
          log_observed_weight = std::log(w);
        }

        for (int k = 0; k < K; ++k) {
          const int cell = s + S * k;
          double log_density = Rf_dnorm4(
            z_p[ell] * sei_p[k], mean_p[cell], sd_p[cell], true
          );

          log_density -= log_norm[cell];
          if (!extrapolate_b) {
            log_density += log_observed_weight;
          }

          dens_sum += std::exp(log_density) * sei_p[k];
        }
      } else {
        for (int k = 0; k < K; ++k) {
          const int cell = s + S * k;
          dens_sum += Rf_dnorm4(
            z_p[ell] * sei_p[k], mean_p[cell], sd_p[cell], false
          ) * sei_p[k];
        }
      }

      out_p[s + S * ell] = dens_sum / K;
    }
  }

  UNPROTECT(1);
  return out;
}
