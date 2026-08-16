#include <Rinternals.h>
#include <R_ext/Error.h>

#include <Rmath.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "plot-root.h"
#include "selnorm/selnorm.h"

// Funnel and regression plots invert the same posterior-row mixture CDF. This
// file owns their shared normal/selected-normal evaluation and root semantics;
// callers only construct plot-specific mean and SD matrices.
namespace {

int plot_matrix_nrow(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[0];
}

int plot_matrix_ncol(SEXP x, const char *name)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
    Rf_error("'%s' must be a matrix.", name);
  }
  return INTEGER(dim)[1];
}

void plot_check_real(SEXP x, const char *name)
{
  if (TYPEOF(x) != REALSXP) {
    Rf_error("'%s' must be numeric.", name);
  }
}

void plot_check_integer(SEXP x, const char *name)
{
  if (TYPEOF(x) != INTSXP) {
    Rf_error("'%s' must be integer.", name);
  }
}

void plot_check_logical(SEXP x, const char *name)
{
  if (TYPEOF(x) != LGLSXP || Rf_length(x) != 1) {
    Rf_error("'%s' must be a logical scalar.", name);
  }
}

int plot_scalar_or_row_int(const int *x, int n, int row, const char *name)
{
  if (n == 1) {
    return x[0];
  }
  if (row >= n) {
    Rf_error("'%s' must have length 1 or one value per posterior sample.", name);
  }
  return x[row];
}

double plot_scalar_or_row_real(const double *x, int n, int row,
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

void plot_validate_mean_sd(SEXP mean, SEXP sd, int *S, int *K)
{
  plot_check_real(mean, "mean");
  plot_check_real(sd, "sd");

  *S = plot_matrix_nrow(mean, "mean");
  *K = plot_matrix_ncol(mean, "mean");
  if (*S < 1 || *K < 1) {
    Rf_error("'mean' and 'sd' must have positive dimensions.");
  }
  if (plot_matrix_nrow(sd, "sd") != *S ||
      plot_matrix_ncol(sd, "sd") != *K) {
    Rf_error("'sd' dimensions must match 'mean'.");
  }

  const R_xlen_t N = XLENGTH(mean);
  for (R_xlen_t i = 0; i < N; ++i) {
    if (!std::isfinite(REAL(mean)[i])) {
      Rf_error("'mean' must contain finite values.");
    }
    if (!std::isfinite(REAL(sd)[i]) || REAL(sd)[i] < 0) {
      Rf_error("'sd' must contain finite non-negative values.");
    }
  }
}

void plot_validate_probs(SEXP probs)
{
  plot_check_real(probs, "probs");
  if (Rf_length(probs) < 1) {
    Rf_error("'probs' must contain at least one probability.");
  }
  for (int i = 0; i < Rf_length(probs); ++i) {
    if (!std::isfinite(REAL(probs)[i]) ||
        REAL(probs)[i] <= 0 || REAL(probs)[i] >= 1) {
      Rf_error("'probs' must contain probabilities in (0, 1).");
    }
  }
}

double plot_validate_weights(SEXP weights, int S)
{
  plot_check_real(weights, "weights");
  if (Rf_length(weights) != S) {
    Rf_error("'weights' must contain one value per posterior mixture row.");
  }

  double total = 0;
  for (int s = 0; s < S; ++s) {
    const double weight = REAL(weights)[s];
    if (!std::isfinite(weight) || weight < 0) {
      Rf_error("'weights' must contain finite non-negative values.");
    }
    total += weight;
  }
  if (!(total > 0) || !std::isfinite(total)) {
    Rf_error("'weights' must have a finite positive sum.");
  }
  return total;
}

void plot_validate_omega(const double *omega, int S, int B)
{
  for (int b = 0; b < B; ++b) {
    for (int s = 0; s < S; ++s) {
      const double weight = omega[s + S * b];
      if (!std::isfinite(weight) || weight < 0) {
        Rf_error("'omega' must contain finite non-negative values.");
      }
    }
  }
}

void plot_set_kernel_data(SelNormKernelData *data, int B, int n_segments,
                          SEXP sign, SEXP q, SEXP z_lower, SEXP z_upper,
                          SEXP phack_z_source, SEXP phack_z_dest,
                          SEXP segment_bounds, SEXP segment_step_bin,
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

bool plot_row_has_active_phack(int mode, double alpha, int phack_kind)
{
  return (mode == SELKERNEL_PHACK_POWER ||
    mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
}

double plot_probability_or_na(double x)
{
  return std::isfinite(x) && x >= 0 && x <= 1 ? x : NA_REAL;
}

double plot_interval_probability(double lower, double upper,
                                 double mean, double sd)
{
  if (!(sd > 0) || lower >= upper) {
    return 0;
  }
  if (std::isinf(lower) && lower < 0) {
    return std::isinf(upper) && upper > 0 ?
      1 : pnorm(upper, mean, sd, true, false);
  }
  if (std::isinf(upper) && upper > 0) {
    return pnorm(lower, mean, sd, false, false);
  }
  if (lower > mean) {
    return pnorm(lower, mean, sd, false, false) -
      pnorm(upper, mean, sd, false, false);
  }
  return pnorm(upper, mean, sd, true, false) -
    pnorm(lower, mean, sd, true, false);
}

struct PlotMixtureContext {
  int S;
  double se;
  const double *mean;
  const double *sd;
  const double *weights;
  double weight_sum;
  bool full_support;
  const int *selected;
  const double *omega;
  const std::vector<double> *row_alpha;
  const std::vector<int> *row_phack;
  const std::vector<int> *row_mode;
  const std::vector<bool> *row_active_phack;
  const SelNormKernelData *data;
  std::vector<char> telescope_plan_valid;
  std::vector<double> telescope_boundary_tail;
  std::vector<double> telescope_omega_diff;
  std::vector<double> telescope_omega_last;
  std::vector<double> telescope_normalizer;
  std::vector<char> log_plan_valid;
  std::vector<int> log_plan_groups;
  std::vector<double> log_lower_score;
  std::vector<double> log_upper_score;
  std::vector<double> log_weight;
};

bool plot_has_selection(const PlotMixtureContext &ctx)
{
  return ctx.selected != NULL && ctx.data != NULL;
}

void plot_prepare_selection(PlotMixtureContext *ctx)
{
  const int B = ctx->data->n_bins;
  const int boundaries = std::max(0, B - 1);
  const size_t boundary_size = static_cast<size_t>(ctx->S) *
    static_cast<size_t>(boundaries);
  const size_t log_size = static_cast<size_t>(ctx->S) *
    static_cast<size_t>(B);

  ctx->telescope_plan_valid.assign(static_cast<size_t>(ctx->S), 0);
  ctx->telescope_boundary_tail.assign(boundary_size, 0);
  ctx->telescope_omega_diff.assign(boundary_size, 0);
  ctx->telescope_omega_last.assign(static_cast<size_t>(ctx->S), 0);
  ctx->telescope_normalizer.assign(static_cast<size_t>(ctx->S), 0);
  ctx->log_plan_valid.assign(static_cast<size_t>(ctx->S), 0);
  ctx->log_plan_groups.assign(static_cast<size_t>(ctx->S), 0);
  ctx->log_lower_score.assign(log_size, 0);
  ctx->log_upper_score.assign(log_size, 0);
  ctx->log_weight.assign(log_size, 0);

  for (int s = 0; s < ctx->S; ++s) {
    const size_t si = static_cast<size_t>(s);
    const int mode = (*ctx->row_mode)[si];
    if (ctx->selected[s] == 0 || !(ctx->sd[s] > 0) ||
        !(ctx->se > 0) || (*ctx->row_active_phack)[si] ||
        !(mode == SELKERNEL_STEP || mode == SELKERNEL_STEP_PHACK_POWER)) {
      continue;
    }

    const size_t boundary_offset = si * static_cast<size_t>(boundaries);
    double *boundary_tail = boundaries > 0 ?
      ctx->telescope_boundary_tail.data() + boundary_offset : NULL;
    double *omega_diff = boundaries > 0 ?
      ctx->telescope_omega_diff.data() + boundary_offset : NULL;
    ctx->telescope_plan_valid[si] =
      cpp_selnorm_step_cdf_telescope_plan(
        ctx->mean[s], ctx->sd[s], ctx->se, ctx->omega + s, *ctx->data,
        boundary_tail, omega_diff,
        &ctx->telescope_omega_last[si],
        &ctx->telescope_normalizer[si], ctx->S, false
      );
    if (ctx->telescope_plan_valid[si]) {
      continue;
    }

    const size_t log_offset = si * static_cast<size_t>(B);
    int n_groups = 0;
    ctx->log_plan_valid[si] = cpp_selnorm_step_cdf_plan(
      ctx->mean[s], ctx->sd[s], ctx->se, ctx->omega + s, *ctx->data,
      ctx->log_lower_score.data() + log_offset,
      ctx->log_upper_score.data() + log_offset,
      ctx->log_weight.data() + log_offset,
      &n_groups, ctx->S, false
    );
    ctx->log_plan_groups[si] = n_groups;
  }
}

double plot_selected_cdf_zero_se(double q, int s,
                                 const PlotMixtureContext &ctx)
{
  const SelNormKernelData &data = *ctx.data;
  const double *omega_s = ctx.omega + s;
  const double signed_mean = data.effect_sign * ctx.mean[s];
  const double signed_q    = data.effect_sign * q;
  double requested_mass  = 0;
  double complement_mass = 0;

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

    const double weight = omega_s[static_cast<size_t>(ctx.S) * b];
    if (data.effect_sign == 1) {
      requested_mass += weight * plot_interval_probability(
        lower, std::min(upper, signed_q), signed_mean, ctx.sd[s]
      );
      complement_mass += weight * plot_interval_probability(
        std::max(lower, signed_q), upper, signed_mean, ctx.sd[s]
      );
    } else {
      requested_mass += weight * plot_interval_probability(
        std::max(lower, signed_q), upper, signed_mean, ctx.sd[s]
      );
      complement_mass += weight * plot_interval_probability(
        lower, std::min(upper, signed_q), signed_mean, ctx.sd[s]
      );
    }
  }

  const double total_mass = requested_mass + complement_mass;
  if (!(total_mass > 0) || !std::isfinite(total_mass)) {
    return NA_REAL;
  }
  return plot_probability_or_na(requested_mass / total_mass);
}

double plot_selected_cdf(double q, int s, const PlotMixtureContext &ctx)
{
  const size_t si = static_cast<size_t>(s);
  if (ctx.se <= 0) {
    return plot_selected_cdf_zero_se(q, s, ctx);
  }
  if ((*ctx.row_active_phack)[si]) {
    return NA_REAL;
  }

  const int B = ctx.data->n_bins;
  if (ctx.telescope_plan_valid[si]) {
    const int boundaries = std::max(0, B - 1);
    const size_t offset = si * static_cast<size_t>(boundaries);
    const double *boundary_tail = boundaries > 0 ?
      ctx.telescope_boundary_tail.data() + offset : NULL;
    const double *omega_diff = boundaries > 0 ?
      ctx.telescope_omega_diff.data() + offset : NULL;
    const double probability = cpp_selnorm_step_cdf_from_telescope_plan(
      q, ctx.mean[s], ctx.sd[s], ctx.se, *ctx.data,
      boundary_tail, omega_diff,
      ctx.telescope_omega_last[si], ctx.telescope_normalizer[si], true
    );
    if (std::isfinite(probability)) {
      return probability;
    }
  }

  if (ctx.log_plan_valid[si]) {
    const size_t offset = si * static_cast<size_t>(B);
    return cpp_selnorm_step_cdf_from_plan(
      q, ctx.mean[s], ctx.sd[s], ctx.se, *ctx.data,
      ctx.log_lower_score.data() + offset,
      ctx.log_upper_score.data() + offset,
      ctx.log_weight.data() + offset,
      ctx.log_plan_groups[si], true
    );
  }

  return cpp_selnorm_kernel_cdf(
    q, ctx.mean[s], ctx.sd[s], ctx.se, ctx.omega + s,
    (*ctx.row_alpha)[si], (*ctx.row_phack)[si], (*ctx.row_mode)[si],
    *ctx.data, ctx.S, true, false
  );
}

double plot_mixture_cdf(double q, const PlotMixtureContext &ctx)
{
  double cdf_sum = 0;
  for (int s = 0; s < ctx.S; ++s) {
    if (ctx.weights[s] == 0) {
      continue;
    }

    double cdf_s;
    if (ctx.sd[s] == 0) {
      cdf_s = q >= ctx.mean[s] ? 1.0 : 0.0;
    } else if (plot_has_selection(ctx) && ctx.selected[s] != 0) {
      cdf_s = plot_selected_cdf(q, s, ctx);
    } else {
      cdf_s = pnorm(q, ctx.mean[s], ctx.sd[s], true, false);
    }
    cdf_s = plot_probability_or_na(cdf_s);
    if (!std::isfinite(cdf_s)) {
      return NA_REAL;
    }
    cdf_sum += ctx.weights[s] * cdf_s;
  }
  return plot_probability_or_na(cdf_sum / ctx.weight_sum);
}

bool plot_mixture_has_full_support(const PlotMixtureContext &ctx)
{
  if (!plot_has_selection(ctx)) {
    return true;
  }

  for (int s = 0; s < ctx.S; ++s) {
    if (ctx.weights[s] == 0) {
      continue;
    }
    if (ctx.selected[s] == 0) {
      return true;
    }
    const size_t si = static_cast<size_t>(s);
    const int mode = (*ctx.row_mode)[si];
    if (mode == SELKERNEL_NORMAL ||
        (mode == SELKERNEL_PHACK_POWER &&
         !(*ctx.row_active_phack)[si])) {
      return true;
    }
  }

  for (int b = 0; b < ctx.data->n_bins; ++b) {
    bool bin_has_support = false;
    for (int s = 0; s < ctx.S; ++s) {
      if (ctx.weights[s] == 0 || ctx.selected[s] == 0) {
        continue;
      }
      const int mode = (*ctx.row_mode)[static_cast<size_t>(s)];
      if ((mode == SELKERNEL_STEP || mode == SELKERNEL_STEP_PHACK_POWER) &&
          ctx.omega[s + ctx.S * b] > 0) {
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

double plot_weighted_empirical_quantile(const PlotMixtureContext &ctx,
                                        double p)
{
  std::vector<std::pair<double, double> > values;
  values.reserve(static_cast<size_t>(ctx.S));
  for (int s = 0; s < ctx.S; ++s) {
    if (ctx.weights[s] > 0) {
      values.push_back(std::make_pair(ctx.mean[s], ctx.weights[s]));
    }
  }
  std::sort(values.begin(), values.end());

  const double threshold = p * ctx.weight_sum;
  double cumulative = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    cumulative += values[i].second;
    if (cumulative >= threshold) {
      return values[i].first;
    }
  }
  return values.empty() ? NA_REAL : values.back().first;
}

double plot_mixture_quantile(double p, const PlotMixtureContext &ctx,
                             double previous, bool use_previous)
{
  bool all_zero_sd = true;
  bool all_positive_sd = true;
  double lower = INFINITY;
  double upper = -INFINITY;
  double step = 0;

  for (int s = 0; s < ctx.S; ++s) {
    if (ctx.weights[s] == 0) {
      continue;
    }
    all_zero_sd     = all_zero_sd && ctx.sd[s] == 0;
    all_positive_sd = all_positive_sd && ctx.sd[s] > 0;
    lower = std::min(lower, ctx.mean[s] - 10 * ctx.sd[s]);
    upper = std::max(upper, ctx.mean[s] + 10 * ctx.sd[s]);
    step  = std::max(step, ctx.sd[s]);
  }

  if (all_zero_sd) {
    return plot_weighted_empirical_quantile(ctx, p);
  }
  if (!std::isfinite(lower) || !std::isfinite(upper) || lower >= upper ||
      !(step > 0)) {
    return NA_REAL;
  }

  const double global_lower = lower;
  const double global_upper = upper;
  const double global_step  = step;
  bool using_previous = use_previous && std::isfinite(previous);
  if (using_previous) {
    step  = std::max(
      RoBMA::plot_root_tolerance(global_lower, global_upper, global_step),
      0.05 * global_step
    );
    lower = previous - step;
    upper = previous + step;
  }

  double lower_value = plot_mixture_cdf(lower, ctx) - p;
  double upper_value = plot_mixture_cdf(upper, ctx) - p;
  if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    if (!using_previous) {
      return NA_REAL;
    }
    using_previous = false;
    lower = global_lower;
    upper = global_upper;
    step  = global_step;
    lower_value = plot_mixture_cdf(lower, ctx) - p;
    upper_value = plot_mixture_cdf(upper, ctx) - p;
  }

  for (int attempt = 0; attempt < 2; ++attempt) {
    for (int i = 0; i < 25; ++i) {
      if (lower_value < 0 && upper_value >= 0) {
        break;
      }
      if (lower_value >= 0) {
        lower -= step;
        lower_value = plot_mixture_cdf(lower, ctx) - p;
      }
      if (upper_value < 0) {
        upper += step;
        upper_value = plot_mixture_cdf(upper, ctx) - p;
      }
      if (!std::isfinite(lower_value) || !std::isfinite(upper_value)) {
        break;
      }
      step *= 2;
    }

    if (lower_value < 0 && upper_value >= 0 &&
        std::isfinite(lower_value) && std::isfinite(upper_value)) {
      break;
    }
    if (!using_previous) {
      return NA_REAL;
    }
    using_previous = false;
    lower = global_lower;
    upper = global_upper;
    step  = global_step;
    lower_value = plot_mixture_cdf(lower, ctx) - p;
    upper_value = plot_mixture_cdf(upper, ctx) - p;
  }

  if (lower_value >= 0 || upper_value < 0 ||
      !std::isfinite(lower_value) || !std::isfinite(upper_value)) {
    return NA_REAL;
  }

  if (all_positive_sd && ctx.full_support) {
    const double tolerance = RoBMA::plot_root_tolerance(
      global_lower, global_upper, global_step
    );
    const double root = RoBMA::plot_brent_root(
      [&ctx, p](double value) {
        return plot_mixture_cdf(value, ctx) - p;
      },
      lower, upper, lower_value, upper_value, tolerance
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
    const double mid_value = plot_mixture_cdf(mid, ctx) - p;
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

SEXP plot_mixture_quantile_matrix(SEXP mean, SEXP sd, SEXP probs,
                                  SEXP weights, PlotMixtureContext *ctx,
                                  bool prepare_selection,
                                  const double *se = NULL,
                                  int se_length = 0)
{
  const int S = plot_matrix_nrow(mean, "mean");
  const int K = plot_matrix_ncol(mean, "mean");
  const int P = Rf_length(probs);
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, K, P));
  std::vector<double> previous(static_cast<size_t>(P), NA_REAL);

  ctx->S          = S;
  ctx->weights    = REAL(weights);
  ctx->weight_sum = plot_validate_weights(weights, S);
  ctx->full_support = plot_mixture_has_full_support(*ctx);
  for (int k = 0; k < K; ++k) {
    ctx->mean = REAL(mean) + static_cast<size_t>(S) * static_cast<size_t>(k);
    ctx->sd   = REAL(sd)   + static_cast<size_t>(S) * static_cast<size_t>(k);
    if (prepare_selection) {
      ctx->se = se[se_length == 1 ? 0 : k];
      plot_prepare_selection(ctx);
    }

    for (int j = 0; j < P; ++j) {
      const double quantile = plot_mixture_quantile(
        REAL(probs)[j], *ctx, previous[static_cast<size_t>(j)], k > 0
      );
      if (!std::isfinite(quantile)) {
        Rf_error(
          "Plot mixture quantiles could not be computed from a valid bracketed CDF."
        );
      }
      REAL(out)[k + K * j] = quantile;
      previous[static_cast<size_t>(j)] = quantile;
    }
  }

  UNPROTECT(1);
  return out;
}

} // namespace

extern "C" SEXP RoBMA_plot_normal_mixture_quantiles(
  SEXP mean, SEXP sd, SEXP probs, SEXP weights)
{
  int S, K;
  plot_validate_mean_sd(mean, sd, &S, &K);
  plot_validate_probs(probs);
  plot_validate_weights(weights, S);

  PlotMixtureContext ctx;
  ctx.se               = 0;
  ctx.selected         = NULL;
  ctx.omega            = NULL;
  ctx.row_alpha        = NULL;
  ctx.row_phack        = NULL;
  ctx.row_mode         = NULL;
  ctx.row_active_phack = NULL;
  ctx.data             = NULL;
  return plot_mixture_quantile_matrix(
    mean, sd, probs, weights, &ctx, false
  );
}

extern "C" SEXP RoBMA_plot_selnorm_mixture_quantiles(
  SEXP mean, SEXP sd, SEXP se, SEXP probs, SEXP weights, SEXP selected,
  SEXP omega, SEXP alpha, SEXP phack_kind, SEXP kernel_mode,
  SEXP z_lower, SEXP z_upper, SEXP sign, SEXP q,
  SEXP phack_z_source, SEXP phack_z_dest, SEXP segment_bounds,
  SEXP segment_step_bin, SEXP segment_phack_region,
  SEXP telescope_probabilities)
{
  int S, K;
  plot_validate_mean_sd(mean, sd, &S, &K);
  plot_validate_probs(probs);
  plot_validate_weights(weights, S);
  plot_check_real(se, "se");
  plot_check_integer(selected, "selected");
  plot_check_real(omega, "omega");
  plot_check_real(alpha, "alpha");
  plot_check_integer(phack_kind, "phack_kind");
  plot_check_integer(kernel_mode, "kernel_mode");
  plot_check_real(z_lower, "z_lower");
  plot_check_real(z_upper, "z_upper");
  plot_check_integer(sign, "sign");
  plot_check_integer(q, "q");
  plot_check_real(phack_z_source, "phack_z_source");
  plot_check_real(phack_z_dest, "phack_z_dest");
  plot_check_real(segment_bounds, "segment_bounds");
  plot_check_integer(segment_step_bin, "segment_step_bin");
  plot_check_integer(segment_phack_region, "segment_phack_region");
  plot_check_logical(telescope_probabilities, "telescope_probabilities");

  if (Rf_length(se) != 1 && Rf_length(se) != K) {
    Rf_error("'se' must have length 1 or one value per prediction.");
  }
  for (int k = 0; k < Rf_length(se); ++k) {
    if (!std::isfinite(REAL(se)[k]) || REAL(se)[k] < 0) {
      Rf_error("'se' must contain finite non-negative values.");
    }
  }
  if (Rf_length(selected) != S) {
    Rf_error("'selected' must contain one value per posterior mixture row.");
  }
  for (int s = 0; s < S; ++s) {
    if (INTEGER(selected)[s] != 0 && INTEGER(selected)[s] != 1) {
      Rf_error("'selected' must contain only zero or one.");
    }
  }

  const int B = plot_matrix_ncol(omega, "omega");
  if (B < 1 || plot_matrix_nrow(omega, "omega") != S) {
    Rf_error("'omega' must have one row per posterior mixture row.");
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
  if (Rf_length(sign) != 1 ||
      (INTEGER(sign)[0] != 1 && INTEGER(sign)[0] != -1)) {
    Rf_error("'sign' must be 1 or -1.");
  }
  if (Rf_length(q) != 1 || INTEGER(q)[0] < 0 || INTEGER(q)[0] > 2) {
    Rf_error("'q' must be 0, 1, or 2.");
  }
  if (Rf_length(phack_z_source) != 2 || Rf_length(phack_z_dest) != 2) {
    Rf_error("'phack_z_source' and 'phack_z_dest' must have length 2.");
  }

  const int n_segments = Rf_length(segment_step_bin);
  if (Rf_length(segment_phack_region) != n_segments ||
      Rf_length(segment_bounds) != n_segments + 1) {
    Rf_error("Segment arrays have incompatible lengths.");
  }
  plot_validate_omega(REAL(omega), S, B);

  SelNormKernelData data;
  plot_set_kernel_data(
    &data, B, n_segments, sign, q, z_lower, z_upper,
    phack_z_source, phack_z_dest, segment_bounds,
    segment_step_bin, segment_phack_region, telescope_probabilities
  );

  std::vector<double> row_alpha(static_cast<size_t>(S));
  std::vector<int> row_phack(static_cast<size_t>(S));
  std::vector<int> row_mode(static_cast<size_t>(S));
  std::vector<bool> row_active_phack(static_cast<size_t>(S));
  for (int s = 0; s < S; ++s) {
    row_alpha[static_cast<size_t>(s)] = plot_scalar_or_row_real(
      REAL(alpha), Rf_length(alpha), s, "alpha"
    );
    row_phack[static_cast<size_t>(s)] = plot_scalar_or_row_int(
      INTEGER(phack_kind), Rf_length(phack_kind), s, "phack_kind"
    );
    row_mode[static_cast<size_t>(s)] = plot_scalar_or_row_int(
      INTEGER(kernel_mode), Rf_length(kernel_mode), s, "kernel_mode"
    );
    row_active_phack[static_cast<size_t>(s)] = plot_row_has_active_phack(
      row_mode[static_cast<size_t>(s)],
      row_alpha[static_cast<size_t>(s)],
      row_phack[static_cast<size_t>(s)]
    );
  }

  PlotMixtureContext ctx;
  ctx.selected         = INTEGER(selected);
  ctx.omega            = REAL(omega);
  ctx.row_alpha        = &row_alpha;
  ctx.row_phack        = &row_phack;
  ctx.row_mode         = &row_mode;
  ctx.row_active_phack = &row_active_phack;
  ctx.data             = &data;

  return plot_mixture_quantile_matrix(
    mean, sd, probs, weights, &ctx, true, REAL(se), Rf_length(se)
  );
}
