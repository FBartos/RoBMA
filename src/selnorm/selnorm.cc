#include "selnorm.h"

#include <JRmath.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

double interval_prob_cdf(double lower, double upper, double mean, double sd)
{
  if (!(sd > 0) || lower >= upper) {
    return 0;
  }

  if (std::isinf(lower) && lower < 0) {
    if (std::isinf(upper) && upper > 0) {
      return 1;
    }
    if (upper <= -1e290) {
      return 0;
    }
    if (upper >= 1e290) {
      return 1;
    }
    return pnorm(upper, mean, sd, true, false);
  }
  if (std::isinf(upper) && upper > 0) {
    return pnorm(lower, mean, sd, false, false);
  }

  const double p_lower = pnorm(lower, mean, sd, true, false);
  const double p_upper = pnorm(upper, mean, sd, true, false);

  double out = p_upper - p_lower;
  if (out < 0 && out > -1e-14) {
    out = 0;
  }
  if (out > 1 && out < 1 + 1e-14) {
    out = 1;
  }
  return out;
}

double cdf_at(double x, double mean, double sd)
{
  if (std::isinf(x)) {
    return x < 0 ? 0 : 1;
  }
  if (x <= -1e290) {
    return 0;
  }
  if (x >= 1e290) {
    return 1;
  }
  return pnorm(x, mean, sd, true, false);
}

double tail_at(double x, double mean, double sd)
{
  if (std::isinf(x)) {
    return x < 0 ? 1 : 0;
  }
  if (x <= -1e290) {
    return 1;
  }
  if (x >= 1e290) {
    return 0;
  }
  return pnorm(x, mean, sd, false, false);
}

double add_logspace(double log_a, double log_b)
{
  if (log_a == -INFINITY) {
    return log_b;
  }
  if (log_b == -INFINITY) {
    return log_a;
  }
  if (log_b > log_a) {
    std::swap(log_a, log_b);
  }

  return log_a + std::log1p(std::exp(log_b - log_a));
}

double sub_logspace(double log_a, double log_b)
{
  if (log_b == -INFINITY) {
    return log_a;
  }
  if (log_b > log_a) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double ratio = std::exp(log_b - log_a);
  if (ratio >= 1) {
    return -INFINITY;
  }

  return log_a + std::log1p(-ratio);
}

double interval_log_prob(double lower, double upper, double mean, double sd)
{
  if (!(sd > 0) || lower >= upper) {
    return -INFINITY;
  }
  if (std::isinf(lower) && lower < 0 && std::isinf(upper) && upper > 0) {
    return 0;
  }

  if (lower >= mean) {
    const double log_s_lower = std::isinf(lower) && lower > 0 ?
      -INFINITY : pnorm(lower, mean, sd, false, true);
    const double log_s_upper = std::isinf(upper) && upper > 0 ?
      -INFINITY : pnorm(upper, mean, sd, false, true);
    return sub_logspace(log_s_lower, log_s_upper);
  }

  if (upper <= mean) {
    const double log_f_upper = std::isinf(upper) && upper < 0 ?
      -INFINITY : pnorm(upper, mean, sd, true, true);
    const double log_f_lower = std::isinf(lower) && lower < 0 ?
      -INFINITY : pnorm(lower, mean, sd, true, true);
    return sub_logspace(log_f_upper, log_f_lower);
  }

  const double prob = interval_prob_cdf(lower, upper, mean, sd);
  return prob > 0 ? std::log(prob) : -INFINITY;
}

double interval_prob(double lower, double upper, double mean, double sd)
{
  const double log_prob = interval_log_prob(lower, upper, mean, sd);
  if (log_prob == -INFINITY) {
    return 0;
  }
  return std::exp(log_prob);
}

double truncated_normal_quantile(double lower, double upper, double mean,
                                 double sd, double u)
{
  if (!(sd > 0) || lower >= upper || !(u >= 0) || !(u <= 1)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (std::isinf(lower) && lower < 0 && std::isinf(upper) && upper > 0) {
    return qnorm(u, mean, sd, true, false);
  }

  if (lower >= mean) {
    const double log_s_lower = std::isinf(lower) && lower > 0 ?
      -INFINITY : pnorm(lower, mean, sd, false, true);
    const double log_s_upper = std::isinf(upper) && upper > 0 ?
      -INFINITY : pnorm(upper, mean, sd, false, true);

    if (log_s_lower == -INFINITY) {
      return lower;
    }

    double log_s;
    if (log_s_upper == -INFINITY) {
      log_s = log_s_lower + std::log1p(-u);
    } else {
      const double ratio = std::exp(log_s_upper - log_s_lower);
      log_s = log_s_lower + std::log((1 - u) + u * ratio);
    }

    return qnorm(log_s, mean, sd, false, true);
  }

  if (upper <= mean) {
    const double log_f_upper = std::isinf(upper) && upper < 0 ?
      -INFINITY : pnorm(upper, mean, sd, true, true);
    const double log_f_lower = std::isinf(lower) && lower < 0 ?
      -INFINITY : pnorm(lower, mean, sd, true, true);

    if (log_f_upper == -INFINITY) {
      return upper;
    }

    double log_f;
    if (log_f_lower == -INFINITY) {
      log_f = log_f_upper + std::log(u);
    } else {
      const double ratio = std::exp(log_f_lower - log_f_upper);
      log_f = log_f_upper + std::log(ratio + u * (1 - ratio));
    }

    return qnorm(log_f, mean, sd, true, true);
  }

  const double f_lower = cdf_at(lower, mean, sd);
  const double f_upper = cdf_at(upper, mean, sd);
  return qnorm(f_lower + u * (f_upper - f_lower), mean, sd, true, false);
}

double std_phi(double x)
{
  if (!std::isfinite(x)) {
    return 0;
  }
  return dnorm(x, 0, 1, false);
}

void truncated_raw_moments(double lower, double upper, double mean, double sd,
                           double *m0, double *m1, double *m2)
{
  *m0 = interval_prob(lower, upper, mean, sd);

  if (*m0 == 0 || !(sd > 0)) {
    *m1 = 0;
    *m2 = 0;
    return;
  }

  const double a = (lower - mean) / sd;
  const double b = (upper - mean) / sd;
  const double pa = std_phi(a);
  const double pb = std_phi(b);

  double a_pa = a * pa;
  double b_pb = b * pb;
  if (!std::isfinite(a_pa)) {
    a_pa = 0;
  }
  if (!std::isfinite(b_pb)) {
    b_pb = 0;
  }

  const double tail_shift   = pa - pb;
  const double second_shift = a_pa - b_pb;

  *m1 = mean * *m0 + sd * tail_shift;
  *m2 = (mean * mean + sd * sd) * *m0 +
    2 * mean * sd * tail_shift + sd * sd * second_shift;
}

double source_power_moment(double lower, double upper, double mean, double sd,
                           double source_lower, double source_upper, int q)
{
  lower = std::max(lower, source_lower);
  upper = std::min(upper, source_upper);
  if (lower >= upper) {
    return 0;
  }

  double m0, m1, m2;
  truncated_raw_moments(lower, upper, mean, sd, &m0, &m1, &m2);

  const double width = source_upper - source_lower;
  if (!(width > 0)) {
    return 0;
  }

  double out;
  if (q == 1) {
    out = (m1 - source_lower * m0) / width;
  } else if (q == 2) {
    out = (m2 - 2 * source_lower * m1 + source_lower * source_lower * m0) /
      (width * width);
  } else {
    out = 0;
  }

  return out < 0 && out > -1e-14 ? 0 : out;
}

double dest_power_moment(double lower, double upper, double mean, double sd,
                         double dest_lower, double dest_upper, int q)
{
  lower = std::max(lower, dest_lower);
  upper = std::min(upper, dest_upper);
  if (lower >= upper) {
    return 0;
  }

  double m0, m1, m2;
  truncated_raw_moments(lower, upper, mean, sd, &m0, &m1, &m2);

  const double width = dest_upper - dest_lower;
  if (!(width > 0)) {
    return 0;
  }

  double out;
  if (q == 1) {
    out = (dest_upper * m0 - m1) / width;
  } else if (q == 2) {
    out = (dest_upper * dest_upper * m0 - 2 * dest_upper * m1 + m2) /
      (width * width);
  } else {
    out = 0;
  }

  return out < 0 && out > -1e-14 ? 0 : out;
}

double omega_at(const double *omega, int bin, int stride)
{
  return omega[(bin - 1) * stride];
}

bool valid_kernel_mode(int mode)
{
  return mode == SELKERNEL_NORMAL ||
    mode == SELKERNEL_STEP ||
    mode == SELKERNEL_PHACK_POWER ||
    mode == SELKERNEL_STEP_PHACK_POWER;
}

bool valid_phack_kind(int q)
{
  return q == 0 || q == 1 || q == 2;
}

int segment_step_bin_at(const SelNormKernelData &data, int segment)
{
  if (data.segment_step_bin != 0) {
    return data.segment_step_bin[segment];
  }
  return static_cast<int>(data.segment_step_bin_real[segment]);
}

int segment_phack_region_at(const SelNormKernelData &data, int segment)
{
  if (data.segment_phack_region != 0) {
    return data.segment_phack_region[segment];
  }
  return static_cast<int>(data.segment_phack_region_real[segment]);
}

bool validate_omega_values(const double *omega, const SelNormKernelData &data,
                           int omega_stride)
{
  for (int b = 0; b < data.n_bins; ++b) {
    const double w = omega_at(omega, b + 1, omega_stride);
    if (!std::isfinite(w) || w < 0) {
      return false;
    }
  }

  return true;
}

double step_log_normalizer_trusted_telescoping(const double *omega,
                                               const SelNormKernelData &data,
                                               double mean_z, double sd_z,
                                               int omega_stride)
{
  double normalizer = omega_at(omega, data.n_bins, omega_stride);
  double total_abs  = std::abs(normalizer);

  for (int b = 0; b < data.n_bins - 1; ++b) {
    const double w_left  = omega_at(omega, b + 1, omega_stride);
    const double w_right = omega_at(omega, b + 2, omega_stride);
    if (w_left == w_right) {
      continue;
    }

    const double contribution = (w_left - w_right) *
      tail_at(data.z_lower[b], mean_z, sd_z);
    normalizer += contribution;
    total_abs  += std::abs(contribution);
  }

  if (!std::isfinite(normalizer) || !(normalizer > 0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (normalizer <= std::numeric_limits<double>::min()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (normalizer <= 1e-8 * total_abs) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return std::log(normalizer);
}

double step_log_normalizer_trusted_partition(const double *omega,
                                             const SelNormKernelData &data,
                                             double mean_z, double sd_z,
                                             int omega_stride)
{
  double out = -INFINITY;

  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;

    while (last + 1 < data.n_bins &&
           omega_at(omega, last + 2, omega_stride) == w) {
      ++last;
    }

    if (w > 0) {
      const double log_prob = interval_log_prob(
        data.z_lower[last], data.z_upper[b], mean_z, sd_z
      );
      if (log_prob != -INFINITY) {
        out = add_logspace(out, std::log(w) + log_prob);
      }
    }

    b = last + 1;
  }

  return out;
}

double step_log_normalizer(const double *omega, const SelNormKernelData &data,
                           double mean_z, double sd_z, int omega_stride)
{
  if (data.trusted_step_partition) {
    if (data.telescope_probabilities) {
      const double log_normalizer = step_log_normalizer_trusted_telescoping(
        omega, data, mean_z, sd_z, omega_stride
      );
      if (std::isfinite(log_normalizer)) {
        return log_normalizer;
      }
    }

    return step_log_normalizer_trusted_partition(
      omega, data, mean_z, sd_z, omega_stride
    );
  }

  double out = -INFINITY;

  for (int b = 0; b < data.n_bins; ++b) {
    const double w = omega_at(omega, b + 1, omega_stride);
    if (w <= 0) {
      continue;
    }

    const double log_prob = interval_log_prob(
      data.z_lower[b], data.z_upper[b], mean_z, sd_z
    );
    if (log_prob == -INFINITY) {
      continue;
    }

    out = add_logspace(out, std::log(w) + log_prob);
  }

  return out;
}

double clamp_probability(double x)
{
  if (x < 0) {
    return 0;
  }
  if (x > 1) {
    return 1;
  }
  return x;
}

void normal_tail_raw_moments(double lower, double mean, double sd,
                             double *m0, double *m1, double *m2)
{
  *m0 = tail_at(lower, mean, sd);

  if (*m0 == 0 || !(sd > 0)) {
    *m1 = 0;
    *m2 = 0;
    return;
  }

  const double a  = (lower - mean) / sd;
  const double pa = std_phi(a);

  double a_pa = a * pa;
  if (!std::isfinite(a_pa)) {
    a_pa = 0;
  }

  *m1 = mean * *m0 + sd * pa;
  *m2 = (mean * mean + sd * sd) * *m0 +
    2 * mean * sd * pa + sd * sd * a_pa;
}

double trusted_interval_prob(double lower, double upper, double mean,
                             double sd, bool *ok)
{
  if (!(sd > 0) || lower >= upper) {
    return 0;
  }

  const double tail_lower = tail_at(lower, mean, sd);
  const double tail_upper = tail_at(upper, mean, sd);
  double prob = tail_lower - tail_upper;

  if (!std::isfinite(prob)) {
    *ok = false;
    return 0;
  }

  const double scale = std::max(std::abs(tail_lower), std::abs(tail_upper));
  const double tol   = 1e-14 * std::max(scale, 1.0);
  if (prob < 0 && prob > -tol) {
    prob = 0;
  } else if (prob < 0) {
    *ok = false;
    return 0;
  }

  return prob;
}

double step_cdf_trusted_telescoping(const double *omega,
                                    const SelNormKernelData &data,
                                    double q_z, double mean_z, double sd_z,
                                    int omega_stride,
                                    bool signed_lower_tail,
                                    bool *ok)
{
  *ok = false;

  const double tail_q = tail_at(q_z, mean_z, sd_z);
  const double cdf_q  = cdf_at(q_z, mean_z, sd_z);
  if (!std::isfinite(tail_q) || !std::isfinite(cdf_q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double w_last = omega_at(omega, data.n_bins, omega_stride);
  double normalizer   = w_last;
  double norm_abs     = std::abs(normalizer);
  double requested    = signed_lower_tail ? w_last * cdf_q : w_last * tail_q;
  double req_abs      = std::abs(requested);

  for (int b = 0; b < data.n_bins - 1; ++b) {
    const double w_left  = omega_at(omega, b + 1, omega_stride);
    const double w_right = omega_at(omega, b + 2, omega_stride);
    if (w_left == w_right) {
      continue;
    }

    const double diff          = w_left - w_right;
    const double tail_boundary = tail_at(data.z_lower[b], mean_z, sd_z);
    const double norm_contrib  = diff * tail_boundary;
    normalizer += norm_contrib;
    norm_abs   += std::abs(norm_contrib);

    double req_contrib = 0;
    if (signed_lower_tail) {
      if (q_z > data.z_lower[b]) {
        double interval = tail_boundary - tail_q;
        const double scale = std::max(std::abs(tail_boundary), std::abs(tail_q));
        const double tol   = 1e-14 * std::max(scale, 1.0);
        if (interval < 0 && interval > -tol) {
          interval = 0;
        } else if (interval < 0) {
          return std::numeric_limits<double>::quiet_NaN();
        }
        req_contrib = diff * interval;
      }
    } else {
      req_contrib = diff * (q_z >= data.z_lower[b] ? tail_q : tail_boundary);
    }

    requested += req_contrib;
    req_abs   += std::abs(req_contrib);
  }

  if (!std::isfinite(normalizer) || !(normalizer > 0) ||
      normalizer <= std::numeric_limits<double>::min() ||
      normalizer <= 1e-8 * norm_abs ||
      !std::isfinite(requested)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double req_tol = 1e-10 * std::max(std::max(req_abs, normalizer), 1.0);
  if (requested < 0 && requested > -req_tol) {
    requested = 0;
  } else if (requested < 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (requested > normalizer && requested < normalizer + req_tol) {
    requested = normalizer;
  } else if (requested > normalizer) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  *ok = true;
  return clamp_probability(requested / normalizer);
}

double step_cdf_logspace(const double *omega, const SelNormKernelData &data,
                         double q_z, double mean_z, double sd_z,
                         int omega_stride, bool signed_lower_tail,
                         bool trusted_partition)
{
  double log_denom      = -INFINITY;
  double log_lower_mass = -INFINITY;
  double log_upper_mass = -INFINITY;

  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;

    if (trusted_partition) {
      while (last + 1 < data.n_bins &&
             omega_at(omega, last + 2, omega_stride) == w) {
        ++last;
      }
    }

    if (w > 0) {
      const double lower = trusted_partition ? data.z_lower[last] : data.z_lower[b];
      const double upper = data.z_upper[b];
      const double log_w = std::log(w);
      const double log_prob = interval_log_prob(lower, upper, mean_z, sd_z);
      if (log_prob != -INFINITY) {
        log_denom = add_logspace(log_denom, log_w + log_prob);
      }

      const double log_lower_prob = interval_log_prob(
        lower, std::min(upper, q_z), mean_z, sd_z
      );
      if (log_lower_prob != -INFINITY) {
        log_lower_mass = add_logspace(log_lower_mass, log_w + log_lower_prob);
      }

      const double log_upper_prob = interval_log_prob(
        std::max(lower, q_z), upper, mean_z, sd_z
      );
      if (log_upper_prob != -INFINITY) {
        log_upper_mass = add_logspace(log_upper_mass, log_w + log_upper_prob);
      }
    }

    b = trusted_partition ? last + 1 : b + 1;
  }

  if (!std::isfinite(log_denom)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double log_requested_mass = signed_lower_tail ?
    log_lower_mass : log_upper_mass;
  if (log_requested_mass == -INFINITY) {
    return 0;
  }

  return clamp_probability(std::exp(log_requested_mass - log_denom));
}

double step_cdf_logspace_with_log_norm(const double *omega,
                                       const SelNormKernelData &data,
                                       double q_z, double mean_z, double sd_z,
                                       int omega_stride,
                                       bool signed_lower_tail,
                                       bool trusted_partition,
                                       double log_denom)
{
  if (!std::isfinite(log_denom)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double log_lower_mass = -INFINITY;
  double log_upper_mass = -INFINITY;

  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;

    if (trusted_partition) {
      while (last + 1 < data.n_bins &&
             omega_at(omega, last + 2, omega_stride) == w) {
        ++last;
      }
    }

    if (w > 0) {
      const double lower = trusted_partition ? data.z_lower[last] : data.z_lower[b];
      const double upper = data.z_upper[b];
      const double log_w = std::log(w);

      const double log_lower_prob = interval_log_prob(
        lower, std::min(upper, q_z), mean_z, sd_z
      );
      if (log_lower_prob != -INFINITY) {
        log_lower_mass = add_logspace(log_lower_mass, log_w + log_lower_prob);
      }

      const double log_upper_prob = interval_log_prob(
        std::max(lower, q_z), upper, mean_z, sd_z
      );
      if (log_upper_prob != -INFINITY) {
        log_upper_mass = add_logspace(log_upper_mass, log_w + log_upper_prob);
      }
    }

    b = trusted_partition ? last + 1 : b + 1;
  }

  const double log_requested_mass = signed_lower_tail ?
    log_lower_mass : log_upper_mass;
  if (log_requested_mass == -INFINITY) {
    return 0;
  }

  return clamp_probability(std::exp(log_requested_mass - log_denom));
}

bool step_moments_trusted_telescoping(const double *omega,
                                      const SelNormKernelData &data,
                                      double signed_mean, double sd,
                                      double sei, int omega_stride,
                                      double *denom,
                                      double *signed_first,
                                      double *signed_second)
{
  const double w_last = omega_at(omega, data.n_bins, omega_stride);
  *denom        = w_last;
  *signed_first = w_last * signed_mean;
  *signed_second = w_last * (sd * sd + signed_mean * signed_mean);

  double denom_abs = std::abs(*denom);

  for (int b = 0; b < data.n_bins - 1; ++b) {
    const double w_left  = omega_at(omega, b + 1, omega_stride);
    const double w_right = omega_at(omega, b + 2, omega_stride);
    if (w_left == w_right) {
      continue;
    }

    double m0, m1, m2;
    normal_tail_raw_moments(
      data.z_lower[b] * sei, signed_mean, sd, &m0, &m1, &m2
    );

    const double diff = w_left - w_right;
    const double d0   = diff * m0;
    *denom        += d0;
    *signed_first += diff * m1;
    *signed_second += diff * m2;
    denom_abs     += std::abs(d0);
  }

  if (!std::isfinite(*denom) || !std::isfinite(*signed_first) ||
      !std::isfinite(*signed_second) || !(*denom > 0) ||
      *denom <= std::numeric_limits<double>::min() ||
      *denom <= 1e-8 * denom_abs) {
    return false;
  }

  const double second_tol = 1e-10 *
    std::max(std::abs(*signed_second), sd * sd + signed_mean * signed_mean);
  if (*signed_second < 0 && *signed_second > -second_tol) {
    *signed_second = 0;
  } else if (*signed_second < 0) {
    return false;
  }

  return true;
}

bool clamp_requested_mass(double *requested, double normalizer, double req_abs)
{
  const double req_tol = 1e-10 *
    std::max(std::max(req_abs, normalizer), 1.0);

  if (*requested < 0 && *requested > -req_tol) {
    *requested = 0;
  } else if (*requested < 0) {
    return false;
  }
  if (*requested > normalizer && *requested < normalizer + req_tol) {
    *requested = normalizer;
  } else if (*requested > normalizer) {
    return false;
  }

  return true;
}

double step_cdf_trusted_telescoping_with_log_norm(
  const double *omega,
  const SelNormKernelData &data,
  double q_z,
  double mean_z,
  double sd_z,
  int omega_stride,
  bool signed_lower_tail,
  double log_normalizer,
  bool *ok)
{
  *ok = false;

  const double normalizer = std::exp(log_normalizer);
  if (!std::isfinite(normalizer) ||
      normalizer <= std::numeric_limits<double>::min()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double tail_q = tail_at(q_z, mean_z, sd_z);
  const double cdf_q  = cdf_at(q_z, mean_z, sd_z);
  if (!std::isfinite(tail_q) || !std::isfinite(cdf_q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double w_last = omega_at(omega, data.n_bins, omega_stride);
  double requested    = signed_lower_tail ? w_last * cdf_q : w_last * tail_q;
  double req_abs      = std::abs(requested);

  for (int b = 0; b < data.n_bins - 1; ++b) {
    const double w_left  = omega_at(omega, b + 1, omega_stride);
    const double w_right = omega_at(omega, b + 2, omega_stride);
    if (w_left == w_right) {
      continue;
    }

    const double diff = w_left - w_right;
    const double tail_boundary = tail_at(data.z_lower[b], mean_z, sd_z);
    double req_contrib = 0;
    if (signed_lower_tail) {
      if (q_z > data.z_lower[b]) {
        double interval = tail_boundary - tail_q;
        const double scale = std::max(std::abs(tail_boundary), std::abs(tail_q));
        const double tol   = 1e-14 * std::max(scale, 1.0);
        if (interval < 0 && interval > -tol) {
          interval = 0;
        } else if (interval < 0) {
          return std::numeric_limits<double>::quiet_NaN();
        }
        req_contrib = diff * interval;
      }
    } else {
      req_contrib = diff * (q_z >= data.z_lower[b] ? tail_q : tail_boundary);
    }

    requested += req_contrib;
    req_abs   += std::abs(req_contrib);
  }

  if (!std::isfinite(requested) ||
      !clamp_requested_mass(&requested, normalizer, req_abs)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  *ok = true;
  return clamp_probability(requested / normalizer);
}

bool step_summary_trusted_telescoping(const double *omega,
                                      const SelNormKernelData &data,
                                      double q_z, double mean_z, double sd_z,
                                      double signed_mean, double sd,
                                      double sei, int omega_stride,
                                      bool signed_lower_tail,
                                      double *cdf,
                                      double *moment_mean,
                                      double *moment_second)
{
  const double tail_q = tail_at(q_z, mean_z, sd_z);
  const double cdf_q  = cdf_at(q_z, mean_z, sd_z);
  if (!std::isfinite(tail_q) || !std::isfinite(cdf_q)) {
    return false;
  }

  const double w_last = omega_at(omega, data.n_bins, omega_stride);
  double normalizer   = w_last;
  double norm_abs     = std::abs(normalizer);
  double requested    = signed_lower_tail ? w_last * cdf_q : w_last * tail_q;
  double req_abs      = std::abs(requested);
  double signed_first = w_last * signed_mean;
  double signed_second = w_last * (sd * sd + signed_mean * signed_mean);

  for (int b = 0; b < data.n_bins - 1; ++b) {
    const double w_left  = omega_at(omega, b + 1, omega_stride);
    const double w_right = omega_at(omega, b + 2, omega_stride);
    if (w_left == w_right) {
      continue;
    }

    const double diff          = w_left - w_right;
    const double tail_boundary = tail_at(data.z_lower[b], mean_z, sd_z);
    const double norm_contrib  = diff * tail_boundary;
    normalizer += norm_contrib;
    norm_abs   += std::abs(norm_contrib);

    double req_contrib = 0;
    if (signed_lower_tail) {
      if (q_z > data.z_lower[b]) {
        double interval = tail_boundary - tail_q;
        const double scale = std::max(std::abs(tail_boundary), std::abs(tail_q));
        const double tol   = 1e-14 * std::max(scale, 1.0);
        if (interval < 0 && interval > -tol) {
          interval = 0;
        } else if (interval < 0) {
          return false;
        }
        req_contrib = diff * interval;
      }
    } else {
      req_contrib = diff * (q_z >= data.z_lower[b] ? tail_q : tail_boundary);
    }

    requested += req_contrib;
    req_abs   += std::abs(req_contrib);

    double m0, m1, m2;
    normal_tail_raw_moments(
      data.z_lower[b] * sei, signed_mean, sd, &m0, &m1, &m2
    );
    signed_first  += diff * m1;
    signed_second += diff * m2;
  }

  if (!std::isfinite(normalizer) || !(normalizer > 0) ||
      normalizer <= std::numeric_limits<double>::min() ||
      normalizer <= 1e-8 * norm_abs ||
      !std::isfinite(requested) ||
      !std::isfinite(signed_first) ||
      !std::isfinite(signed_second)) {
    return false;
  }

  if (!clamp_requested_mass(&requested, normalizer, req_abs)) {
    return false;
  }

  const double second_tol = 1e-10 *
    std::max(std::abs(signed_second), sd * sd + signed_mean * signed_mean);
  if (signed_second < 0 && signed_second > -second_tol) {
    signed_second = 0;
  } else if (signed_second < 0) {
    return false;
  }

  *cdf           = clamp_probability(requested / normalizer);
  *moment_mean   = data.effect_sign * signed_first / normalizer;
  *moment_second = signed_second / normalizer;
  return true;
}

bool step_cdf_pair_trusted_telescoping(const double *omega,
                                       const SelNormKernelData &data,
                                       double q1_z, bool q1_signed_lower_tail,
                                       double q2_z, bool q2_signed_lower_tail,
                                       double mean_z, double sd_z,
                                       int omega_stride,
                                       double *p1, double *p2,
                                       double *inverse_weight)
{
  const double q_z[2] = {q1_z, q2_z};
  const bool signed_lower_tail[2] = {
    q1_signed_lower_tail, q2_signed_lower_tail
  };
  double tail_q[2];
  double cdf_q[2];

  for (int i = 0; i < 2; ++i) {
    tail_q[i] = tail_at(q_z[i], mean_z, sd_z);
    cdf_q[i]  = cdf_at(q_z[i], mean_z, sd_z);
    if (!std::isfinite(tail_q[i]) || !std::isfinite(cdf_q[i])) {
      return false;
    }
  }

  const double w_last = omega_at(omega, data.n_bins, omega_stride);
  double normalizer   = w_last;
  double norm_abs     = std::abs(normalizer);
  double requested[2] = {
    signed_lower_tail[0] ? w_last * cdf_q[0] : w_last * tail_q[0],
    signed_lower_tail[1] ? w_last * cdf_q[1] : w_last * tail_q[1]
  };
  double req_abs[2] = {std::abs(requested[0]), std::abs(requested[1])};

  for (int b = 0; b < data.n_bins - 1; ++b) {
    const double w_left  = omega_at(omega, b + 1, omega_stride);
    const double w_right = omega_at(omega, b + 2, omega_stride);
    if (w_left == w_right) {
      continue;
    }

    const double diff          = w_left - w_right;
    const double tail_boundary = tail_at(data.z_lower[b], mean_z, sd_z);
    const double norm_contrib  = diff * tail_boundary;
    normalizer += norm_contrib;
    norm_abs   += std::abs(norm_contrib);

    for (int i = 0; i < 2; ++i) {
      double req_contrib = 0;
      if (signed_lower_tail[i]) {
        if (q_z[i] > data.z_lower[b]) {
          double interval = tail_boundary - tail_q[i];
          const double scale = std::max(std::abs(tail_boundary), std::abs(tail_q[i]));
          const double tol   = 1e-14 * std::max(scale, 1.0);
          if (interval < 0 && interval > -tol) {
            interval = 0;
          } else if (interval < 0) {
            return false;
          }
          req_contrib = diff * interval;
        }
      } else {
        req_contrib = diff * (q_z[i] >= data.z_lower[b] ? tail_q[i] : tail_boundary);
      }

      requested[i] += req_contrib;
      req_abs[i]   += std::abs(req_contrib);
    }
  }

  if (!std::isfinite(normalizer) || !(normalizer > 0) ||
      normalizer <= std::numeric_limits<double>::min() ||
      normalizer <= 1e-8 * norm_abs ||
      !std::isfinite(requested[0]) ||
      !std::isfinite(requested[1])) {
    return false;
  }

  if (!clamp_requested_mass(&requested[0], normalizer, req_abs[0]) ||
      !clamp_requested_mass(&requested[1], normalizer, req_abs[1])) {
    return false;
  }

  *p1             = clamp_probability(requested[0] / normalizer);
  *p2             = clamp_probability(requested[1] / normalizer);
  *inverse_weight = 1 / normalizer;
  return true;
}

bool step_cdf_pair_logspace(const double *omega,
                            const SelNormKernelData &data,
                            double q1_z, bool q1_signed_lower_tail,
                            double q2_z, bool q2_signed_lower_tail,
                            double mean_z, double sd_z,
                            int omega_stride,
                            bool trusted_partition,
                            double *p1, double *p2,
                            double *inverse_weight)
{
  const double q_z[2] = {q1_z, q2_z};
  const bool signed_lower_tail[2] = {
    q1_signed_lower_tail, q2_signed_lower_tail
  };
  double log_mass[2] = {-INFINITY, -INFINITY};
  double log_denom   = -INFINITY;

  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;

    if (trusted_partition) {
      while (last + 1 < data.n_bins &&
             omega_at(omega, last + 2, omega_stride) == w) {
        ++last;
      }
    }

    if (w > 0) {
      const double lower = trusted_partition ? data.z_lower[last] : data.z_lower[b];
      const double upper = data.z_upper[b];
      const double log_w = std::log(w);
      const double log_prob = interval_log_prob(lower, upper, mean_z, sd_z);

      if (log_prob != -INFINITY) {
        log_denom = add_logspace(log_denom, log_w + log_prob);
      }

      for (int i = 0; i < 2; ++i) {
        const double mass_lower = signed_lower_tail[i] ?
          lower : std::max(lower, q_z[i]);
        const double mass_upper = signed_lower_tail[i] ?
          std::min(upper, q_z[i]) : upper;
        const double log_requested = interval_log_prob(
          mass_lower, mass_upper, mean_z, sd_z
        );

        if (log_requested != -INFINITY) {
          log_mass[i] = add_logspace(log_mass[i], log_w + log_requested);
        }
      }
    }

    b = trusted_partition ? last + 1 : b + 1;
  }

  if (!std::isfinite(log_denom)) {
    return false;
  }

  *p1 = log_mass[0] == -INFINITY ?
    0 : clamp_probability(std::exp(log_mass[0] - log_denom));
  *p2 = log_mass[1] == -INFINITY ?
    0 : clamp_probability(std::exp(log_mass[1] - log_denom));
  *inverse_weight = std::exp(-log_denom);

  return true;
}

void step_moments_partition(const double *omega,
                            const SelNormKernelData &data,
                            double signed_mean, double sd, double sei,
                            int omega_stride, bool trusted_partition,
                            double *denom, double *signed_first,
                            double *signed_second)
{
  *denom         = 0;
  *signed_first  = 0;
  *signed_second = 0;

  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;

    if (trusted_partition) {
      while (last + 1 < data.n_bins &&
             omega_at(omega, last + 2, omega_stride) == w) {
        ++last;
      }
    }

    if (w > 0) {
      double m0, m1, m2;
      truncated_raw_moments(
        (trusted_partition ? data.z_lower[last] : data.z_lower[b]) * sei,
        data.z_upper[b] * sei,
        signed_mean,
        sd,
        &m0,
        &m1,
        &m2
      );

      *denom         += w * m0;
      *signed_first  += w * m1;
      *signed_second += w * m2;
    }

    b = trusted_partition ? last + 1 : b + 1;
  }
}

double step_rng_trusted_partition(const double *omega,
                                  const SelNormKernelData &data,
                                  double mean_z, double sd_z,
                                  double u_bin, double u_interval,
                                  int omega_stride,
                                  bool *ok)
{
  *ok = false;

  std::vector<double> mass;
  std::vector<double> lower;
  std::vector<double> upper;
  double denom = 0;

  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;

    while (last + 1 < data.n_bins &&
           omega_at(omega, last + 2, omega_stride) == w) {
      ++last;
    }

    if (w > 0) {
      bool interval_ok = true;
      const double group_lower = data.z_lower[last];
      const double group_upper = data.z_upper[b];
      const double prob = trusted_interval_prob(
        group_lower, group_upper, mean_z, sd_z, &interval_ok
      );
      if (!interval_ok) {
        return std::numeric_limits<double>::quiet_NaN();
      }

      const double group_mass = w * prob;
      if (!std::isfinite(group_mass) || group_mass < 0) {
        return std::numeric_limits<double>::quiet_NaN();
      }
      if (group_mass > 0) {
        mass.push_back(group_mass);
        lower.push_back(group_lower);
        upper.push_back(group_upper);
        denom += group_mass;
      }
    }

    b = last + 1;
  }

  if (!std::isfinite(denom) || !(denom > 0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double target = u_bin * denom;
  double cumulative = 0;
  int selected = static_cast<int>(mass.size()) - 1;

  for (int g = 0; g < static_cast<int>(mass.size()); ++g) {
    cumulative += mass[static_cast<size_t>(g)];
    if (target <= cumulative) {
      selected = g;
      break;
    }
  }

  *ok = true;
  return truncated_normal_quantile(
    lower[static_cast<size_t>(selected)],
    upper[static_cast<size_t>(selected)],
    mean_z,
    sd_z,
    u_interval
  );
}

double step_rng_logspace(const double *omega, const SelNormKernelData &data,
                         double mean_z, double sd_z,
                         double u_bin, double u_interval,
                         int omega_stride, bool trusted_partition)
{
  std::vector<double> log_mass;
  std::vector<double> lower;
  std::vector<double> upper;
  double max_log_mass = -INFINITY;

  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;

    if (trusted_partition) {
      while (last + 1 < data.n_bins &&
             omega_at(omega, last + 2, omega_stride) == w) {
        ++last;
      }
    }

    if (w > 0) {
      const double group_lower = trusted_partition ?
        data.z_lower[last] : data.z_lower[b];
      const double group_upper = data.z_upper[b];
      const double lm = std::log(w) + interval_log_prob(
        group_lower, group_upper, mean_z, sd_z
      );
      if (lm != -INFINITY) {
        log_mass.push_back(lm);
        lower.push_back(group_lower);
        upper.push_back(group_upper);
        if (lm > max_log_mass) {
          max_log_mass = lm;
        }
      }
    }

    b = trusted_partition ? last + 1 : b + 1;
  }

  if (max_log_mass == -INFINITY || !std::isfinite(max_log_mass)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double denom = 0;
  for (int g = 0; g < static_cast<int>(log_mass.size()); ++g) {
    denom += std::exp(log_mass[static_cast<size_t>(g)] - max_log_mass);
  }

  if (!(denom > 0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double target = u_bin * denom;
  double cumulative = 0;
  int selected = static_cast<int>(log_mass.size()) - 1;

  for (int g = 0; g < static_cast<int>(log_mass.size()); ++g) {
    cumulative += std::exp(log_mass[static_cast<size_t>(g)] - max_log_mass);
    if (target <= cumulative) {
      selected = g;
      break;
    }
  }

  return truncated_normal_quantile(
    lower[static_cast<size_t>(selected)],
    upper[static_cast<size_t>(selected)],
    mean_z,
    sd_z,
    u_interval
  );
}

double phack_beta(double alpha, int q, const SelNormKernelData &data,
                  double mean_z, double sd_z)
{
  if (!(alpha > 0) || q <= 0) {
    return 0;
  }

  const double source_lower = data.phack_z_source[0];
  const double source_upper = data.phack_z_source[1];
  const double dest_lower   = data.phack_z_dest[0];
  const double dest_upper   = data.phack_z_dest[1];

  const double source_moment = source_power_moment(
    source_lower, source_upper, mean_z, sd_z, source_lower, source_upper, q
  );
  const double dest_moment = dest_power_moment(
    dest_lower, dest_upper, mean_z, sd_z, dest_lower, dest_upper, q
  );

  if (!(dest_moment > 0)) {
    return 0;
  }

  return alpha * source_moment / dest_moment;
}

double phack_kernel_weight(double z, double alpha, double beta, int q,
                           const SelNormKernelData &data)
{
  if (q <= 0 || !(alpha > 0)) {
    return 1;
  }

  const double source_lower = data.phack_z_source[0];
  const double source_upper = data.phack_z_source[1];
  const double dest_lower   = data.phack_z_dest[0];
  const double dest_upper   = data.phack_z_dest[1];

  if (z >= source_lower && z <= source_upper) {
    const double x = (z - source_lower) / (source_upper - source_lower);
    const double depletion = q == 1 ? x : x * x;
    return 1 - alpha * depletion;
  }

  if (z > dest_lower && z <= dest_upper) {
    const double x = (dest_upper - z) / (dest_upper - dest_lower);
    const double boost = q == 1 ? x : x * x;
    return 1 + beta * boost;
  }

  return 1;
}

double combined_normalizer(const double *omega, double alpha, double beta, int q,
                           const SelNormKernelData &data, double mean_z,
                           double sd_z, int omega_stride)
{
  double out = 0;

  for (int s = 0; s < data.n_segments; ++s) {
    const double lower = data.segment_bounds[s];
    const double upper = data.segment_bounds[s + 1];
    const int step_bin = segment_step_bin_at(data, s);

    if (step_bin < 1 || step_bin > data.n_bins) {
      continue;
    }

    const double w = omega_at(omega, step_bin, omega_stride);
    if (w <= 0) {
      continue;
    }

    double contribution = interval_prob(lower, upper, mean_z, sd_z);
    const int region = segment_phack_region_at(data, s);

    if (region == 1 && alpha > 0 && q > 0) {
      contribution -= alpha * source_power_moment(
        lower, upper, mean_z, sd_z, data.phack_z_source[0],
        data.phack_z_source[1], q
      );
    } else if (region == 2 && beta > 0 && q > 0) {
      contribution += beta * dest_power_moment(
        lower, upper, mean_z, sd_z, data.phack_z_dest[0],
        data.phack_z_dest[1], q
      );
    }

    if (contribution < 0 && contribution > -1e-12) {
      contribution = 0;
    }
    out += w * contribution;
  }

  return out;
}

}

bool selnorm_bound_is_negative_infinity(double x)
{
  return x <= -1e290 || (std::isinf(x) && x < 0);
}

bool selnorm_bound_is_positive_infinity(double x)
{
  return x >= 1e290 || (std::isinf(x) && x > 0);
}

bool selnorm_same_finite_bound(double x, double y)
{
  if (!std::isfinite(x) || !std::isfinite(y)) {
    return false;
  }

  const double scale = 1 + std::max(std::abs(x), std::abs(y));
  return std::abs(x - y) <= 1e-12 * scale;
}

bool selnorm_is_descending_step_partition(const double *z_lower,
                                          const double *z_upper,
                                          int n_bins)
{
  if (n_bins < 1 ||
      !selnorm_bound_is_positive_infinity(z_upper[0]) ||
      !selnorm_bound_is_negative_infinity(z_lower[n_bins - 1])) {
    return false;
  }

  for (int b = 0; b < n_bins; ++b) {
    if (!std::isfinite(z_lower[b]) &&
        !selnorm_bound_is_negative_infinity(z_lower[b]) &&
        !selnorm_bound_is_positive_infinity(z_lower[b])) {
      return false;
    }
    if (!std::isfinite(z_upper[b]) &&
        !selnorm_bound_is_negative_infinity(z_upper[b]) &&
        !selnorm_bound_is_positive_infinity(z_upper[b])) {
      return false;
    }
    if (!(z_lower[b] < z_upper[b])) {
      return false;
    }
  }

  for (int b = 1; b < n_bins; ++b) {
    if (!selnorm_same_finite_bound(z_upper[b], z_lower[b - 1])) {
      return false;
    }
  }

  return true;
}

double cpp_selnorm_kernel_lpdf(
  double y,
  double mu_num,
  double sigma_num,
  double mu_norm,
  double sigma_norm,
  double sei,
  double weight,
  const double *omega,
  int obs_bin,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride,
  bool validate_omega)
{
  if (!(sigma_num > 0) || !(sigma_norm > 0) || !(sei > 0) || !(weight > 0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (!valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const int mode = kernel_mode;
  const int q = phack_kind > 0 ? phack_kind : data.q;
  const bool use_step = mode == SELKERNEL_STEP ||
    mode == SELKERNEL_STEP_PHACK_POWER;
  const bool use_phack = (mode == SELKERNEL_PHACK_POWER ||
    mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;

  const double signed_y       = data.effect_sign * y;
  const double signed_mu_num  = data.effect_sign * mu_num;
  double log_lik = dnorm(signed_y, signed_mu_num, sigma_num, true);

  if (mode == SELKERNEL_NORMAL) {
    return weight * log_lik;
  }

  if (use_phack && (!(alpha >= 0) || !(alpha < 1) || !valid_phack_kind(q))) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_mu_norm = data.effect_sign * mu_norm;
  const double mean_z         = signed_mu_norm / sei;
  const double sd_z           = sigma_norm / sei;

  double beta = 0;
  if (use_phack) {
    beta = phack_beta(alpha, q, data, mean_z, sd_z);
    const double z = signed_y / sei;
    const double r = phack_kernel_weight(z, alpha, beta, q, data);
    if (!(r > 0)) {
      return -INFINITY;
    }
    log_lik += std::log(r);
  }

  if (use_step) {
    if (obs_bin < 1 || obs_bin > data.n_bins) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    if (validate_omega) {
      if (!validate_omega_values(omega, data, omega_stride)) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    }

    const double observed_weight = omega_at(omega, obs_bin, omega_stride);
    if (!(observed_weight > 0)) {
      return -INFINITY;
    }
    if (observed_weight != 1) {
      log_lik += std::log(observed_weight);
    }
  }

  double log_normalizer = 0;
  if (mode == SELKERNEL_STEP) {
    log_normalizer = step_log_normalizer(
      omega, data, mean_z, sd_z, omega_stride
    );
  } else if (mode == SELKERNEL_STEP_PHACK_POWER) {
    const double normalizer = combined_normalizer(
      omega, use_phack ? alpha : 0, beta, use_phack ? q : 0,
      data, mean_z, sd_z, omega_stride
    );
    if (!(normalizer > 0)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    log_normalizer = std::log(normalizer);
  }

  if (!std::isfinite(log_normalizer)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  log_lik -= log_normalizer;
  return weight == 1 ? log_lik : weight * log_lik;
}

double cpp_selnorm_kernel_log_norm(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (kernel_mode == SELKERNEL_NORMAL ||
      kernel_mode == SELKERNEL_PHACK_POWER) {
    return 0;
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const int q = phack_kind > 0 ? phack_kind : data.q;
  const bool use_phack = kernel_mode == SELKERNEL_STEP_PHACK_POWER &&
    phack_kind > 0 && alpha > 0;

  if (use_phack && (!(alpha >= 0) || !(alpha < 1) || !valid_phack_kind(q))) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double mean_z = data.effect_sign * mean / sei;
  const double sd_z   = sd / sei;
  double log_normalizer = 0;

  if (kernel_mode == SELKERNEL_STEP) {
    log_normalizer = step_log_normalizer(
      omega, data, mean_z, sd_z, omega_stride
    );
  } else if (kernel_mode == SELKERNEL_STEP_PHACK_POWER) {
    const double beta = use_phack ? phack_beta(alpha, q, data, mean_z, sd_z) : 0;
    const double normalizer = combined_normalizer(
      omega, use_phack ? alpha : 0, beta, use_phack ? q : 0,
      data, mean_z, sd_z, omega_stride
    );
    if (!(normalizer > 0)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    log_normalizer = std::log(normalizer);
  }

  if (!std::isfinite(log_normalizer)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return log_normalizer;
}

double cpp_selnorm_kernel_cdf(
  double q,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride,
  bool lower_tail,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    return pnorm(q, mean, sd, lower_tail, false);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    return pnorm(q, mean, sd, lower_tail, false);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_q = data.effect_sign * q;
  const double q_z      = signed_q / sei;
  const double mean_z   = data.effect_sign * mean / sei;
  const double sd_z     = sd / sei;
  const bool signed_lower_tail = data.effect_sign == 1 ?
    lower_tail : !lower_tail;

  if (data.trusted_step_partition && data.telescope_probabilities) {
    bool ok = false;
    const double out = step_cdf_trusted_telescoping(
      omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail, &ok
    );
    if (ok) {
      return out;
    }
  }

  return step_cdf_logspace(
    omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail,
    data.trusted_step_partition
  );
}

double cpp_selnorm_kernel_cdf_with_log_norm(
  double q,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double log_normalizer,
  int omega_stride,
  bool lower_tail,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    return pnorm(q, mean, sd, lower_tail, false);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    return pnorm(q, mean, sd, lower_tail, false);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double signed_q = data.effect_sign * q;
  const double q_z      = signed_q / sei;
  const double mean_z   = data.effect_sign * mean / sei;
  const double sd_z     = sd / sei;
  const bool signed_lower_tail = data.effect_sign == 1 ?
    lower_tail : !lower_tail;

  if (data.trusted_step_partition && data.telescope_probabilities) {
    bool ok = false;
    const double out = step_cdf_trusted_telescoping_with_log_norm(
      omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail,
      log_normalizer, &ok
    );
    if (ok) {
      return out;
    }
  }

  return step_cdf_logspace_with_log_norm(
    omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail,
    data.trusted_step_partition, log_normalizer
  );
}

void cpp_selnorm_kernel_moments(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *moment_mean,
  double *moment_second,
  int omega_stride,
  bool validate_omega)
{
  *moment_mean   = std::numeric_limits<double>::quiet_NaN();
  *moment_second = std::numeric_limits<double>::quiet_NaN();

  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return;
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    *moment_mean   = mean;
    *moment_second = sd * sd + mean * mean;
    return;
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return;
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    *moment_mean   = mean;
    *moment_second = sd * sd + mean * mean;
    return;
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return;
  }

  const double signed_mean = data.effect_sign * mean;
  double denom;
  double signed_first;
  double signed_second;

  if (!(data.trusted_step_partition && data.telescope_probabilities &&
        step_moments_trusted_telescoping(
          omega, data, signed_mean, sd, sei, omega_stride,
          &denom, &signed_first, &signed_second
        ))) {
    step_moments_partition(
      omega, data, signed_mean, sd, sei, omega_stride,
      data.trusted_step_partition, &denom, &signed_first, &signed_second
    );
  }

  if (!(denom > 0)) {
    return;
  }

  *moment_mean   = data.effect_sign * signed_first / denom;
  *moment_second = signed_second / denom;
}

void cpp_selnorm_kernel_summary(
  double q,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *cdf,
  double *moment_mean,
  double *moment_second,
  int omega_stride,
  bool lower_tail,
  bool validate_omega)
{
  *cdf           = std::numeric_limits<double>::quiet_NaN();
  *moment_mean   = std::numeric_limits<double>::quiet_NaN();
  *moment_second = std::numeric_limits<double>::quiet_NaN();

  if (!(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return;
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    *cdf           = pnorm(q, mean, sd, lower_tail, false);
    *moment_mean   = mean;
    *moment_second = sd * sd + mean * mean;
    return;
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return;
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    *cdf           = pnorm(q, mean, sd, lower_tail, false);
    *moment_mean   = mean;
    *moment_second = sd * sd + mean * mean;
    return;
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return;
  }

  const double signed_q          = data.effect_sign * q;
  const double q_z               = signed_q / sei;
  const double signed_mean       = data.effect_sign * mean;
  const double mean_z            = signed_mean / sei;
  const double sd_z              = sd / sei;
  const bool signed_lower_tail   = data.effect_sign == 1 ?
    lower_tail : !lower_tail;

  if (data.trusted_step_partition && data.telescope_probabilities &&
      step_summary_trusted_telescoping(
        omega, data, q_z, mean_z, sd_z, signed_mean, sd, sei,
        omega_stride, signed_lower_tail, cdf, moment_mean, moment_second
      )) {
    return;
  }

  cpp_selnorm_kernel_moments(
    mean, sd, sei, omega, alpha, phack_kind, kernel_mode, data,
    moment_mean, moment_second, omega_stride, false
  );
  *cdf = cpp_selnorm_kernel_cdf(
    q, mean, sd, sei, omega, alpha, phack_kind, kernel_mode, data,
    omega_stride, lower_tail, false
  );
}

double cpp_selnorm_kernel_threshold(
  double z_threshold,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *inverse_weight,
  int omega_stride,
  bool validate_omega)
{
  *inverse_weight = std::numeric_limits<double>::quiet_NaN();

  if (!(z_threshold >= 0) || !(sd > 0) || !(sei > 0) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double q_upper = z_threshold * sei;
  const double q_lower = -q_upper;

  if (kernel_mode == SELKERNEL_NORMAL) {
    *inverse_weight = 1;
    return pnorm(q_upper, mean, sd, false, false) +
      pnorm(q_lower, mean, sd, true, false);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    *inverse_weight = 1;
    return pnorm(q_upper, mean, sd, false, false) +
      pnorm(q_lower, mean, sd, true, false);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double mean_z = data.effect_sign * mean / sei;
  const double sd_z   = sd / sei;

  if (data.trusted_step_partition && data.telescope_probabilities) {
    double upper_prob;
    double lower_prob;
    const double q_upper_z = data.effect_sign * z_threshold;
    const double q_lower_z = -q_upper_z;
    const bool upper_signed_lower_tail = data.effect_sign == 1 ? false : true;
    const bool lower_signed_lower_tail = data.effect_sign == 1 ? true : false;

    if (step_cdf_pair_trusted_telescoping(
          omega, data,
          q_upper_z, upper_signed_lower_tail,
          q_lower_z, lower_signed_lower_tail,
          mean_z, sd_z, omega_stride,
          &upper_prob, &lower_prob, inverse_weight
        )) {
      return upper_prob + lower_prob;
    }
  }

  double upper_prob;
  double lower_prob;
  if (step_cdf_pair_logspace(
        omega, data,
        data.effect_sign * z_threshold,
        data.effect_sign == 1 ? false : true,
        -data.effect_sign * z_threshold,
        data.effect_sign == 1 ? true : false,
        mean_z, sd_z, omega_stride, data.trusted_step_partition,
        &upper_prob, &lower_prob, inverse_weight
      )) {
    return upper_prob + lower_prob;
  }

  return std::numeric_limits<double>::quiet_NaN();
}

double cpp_selnorm_kernel_rng(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double u_bin,
  double u_interval,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride,
  bool validate_omega)
{
  if (!(sd > 0) || !(sei > 0) ||
      !(u_bin >= 0) || !(u_bin <= 1) ||
      !(u_interval >= 0) || !(u_interval <= 1) ||
      !valid_kernel_mode(kernel_mode) ||
      !(data.effect_sign == 1 || data.effect_sign == -1) ||
      !valid_phack_kind(phack_kind) ||
      !valid_phack_kind(data.q)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (kernel_mode == SELKERNEL_NORMAL) {
    return qnorm(u_interval, mean, sd, true, false);
  }

  const bool use_phack = (kernel_mode == SELKERNEL_PHACK_POWER ||
    kernel_mode == SELKERNEL_STEP_PHACK_POWER) && phack_kind > 0 && alpha > 0;
  if (use_phack) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (kernel_mode == SELKERNEL_PHACK_POWER) {
    return qnorm(u_interval, mean, sd, true, false);
  }

  if (validate_omega && !validate_omega_values(omega, data, omega_stride)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double mean_z = data.effect_sign * mean / sei;
  const double sd_z   = sd / sei;

  if (data.trusted_step_partition && data.telescope_probabilities) {
    bool ok = false;
    const double z = step_rng_trusted_partition(
      omega, data, mean_z, sd_z, u_bin, u_interval, omega_stride, &ok
    );
    if (ok) {
      return data.effect_sign * z * sei;
    }
  }

  const double z = step_rng_logspace(
    omega, data, mean_z, sd_z, u_bin, u_interval, omega_stride,
    data.trusted_step_partition
  );

  return data.effect_sign * z * sei;
}
