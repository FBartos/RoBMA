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

double step_log_normalizer(const double *omega, const SelNormKernelData &data,
                           double mean_z, double sd_z, int omega_stride)
{
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
  const double z              = signed_y / sei;
  const double mean_z         = signed_mu_norm / sei;
  const double sd_z           = sigma_norm / sei;

  double beta = 0;
  if (use_phack) {
    beta = phack_beta(alpha, q, data, mean_z, sd_z);
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
    log_lik += std::log(observed_weight);
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
  return weight * log_lik;
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
  double log_denom      = -INFINITY;
  double log_lower_mass = -INFINITY;
  double log_upper_mass = -INFINITY;

  for (int b = 0; b < data.n_bins; ++b) {
    const double w = omega_at(omega, b + 1, omega_stride);
    if (w <= 0) {
      continue;
    }

    const double lower = data.z_lower[b];
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

  if (!std::isfinite(log_denom)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double log_requested_mass = data.effect_sign == 1 ?
    (lower_tail ? log_lower_mass : log_upper_mass) :
    (lower_tail ? log_upper_mass : log_lower_mass);
  double out = log_requested_mass == -INFINITY ? 0 :
    std::exp(log_requested_mass - log_denom);

  if (out < 0 && out > -1e-12) {
    out = 0;
  }
  if (out > 1 && out < 1 + 1e-12) {
    out = 1;
  }
  return std::min(std::max(out, 0.0), 1.0);
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
  double denom             = 0;
  double signed_first      = 0;
  double signed_second     = 0;

  for (int b = 0; b < data.n_bins; ++b) {
    const double w = omega_at(omega, b + 1, omega_stride);
    if (w <= 0) {
      continue;
    }

    double m0, m1, m2;
    truncated_raw_moments(
      data.z_lower[b] * sei,
      data.z_upper[b] * sei,
      signed_mean,
      sd,
      &m0,
      &m1,
      &m2
    );

    denom         += w * m0;
    signed_first  += w * m1;
    signed_second += w * m2;
  }

  if (!(denom > 0)) {
    return;
  }

  *moment_mean   = data.effect_sign * signed_first / denom;
  *moment_second = signed_second / denom;
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

  std::vector<double> log_mass(data.n_bins);
  double max_log_mass = -INFINITY;

  const double mean_z = data.effect_sign * mean / sei;
  const double sd_z   = sd / sei;

  for (int b = 0; b < data.n_bins; ++b) {
    const double w = omega_at(omega, b + 1, omega_stride);
    if (w <= 0) {
      log_mass[b] = -INFINITY;
      continue;
    }

    log_mass[b] = std::log(w) + interval_log_prob(
      data.z_lower[b], data.z_upper[b], mean_z, sd_z
    );
    if (log_mass[b] > max_log_mass) {
      max_log_mass = log_mass[b];
    }
  }

  if (max_log_mass == -INFINITY || !std::isfinite(max_log_mass)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double denom = 0;
  for (int b = 0; b < data.n_bins; ++b) {
    if (log_mass[b] > -INFINITY) {
      denom += std::exp(log_mass[b] - max_log_mass);
    }
  }

  if (!(denom > 0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double target = u_bin * denom;
  double cumulative   = 0;
  int selected_bin    = data.n_bins - 1;

  for (int b = 0; b < data.n_bins; ++b) {
    if (log_mass[b] == -INFINITY) {
      continue;
    }

    cumulative += std::exp(log_mass[b] - max_log_mass);
    if (target <= cumulative) {
      selected_bin = b;
      break;
    }
  }

  const double z = truncated_normal_quantile(
    data.z_lower[selected_bin],
    data.z_upper[selected_bin],
    mean_z,
    sd_z,
    u_interval
  );

  return data.effect_sign * z * sei;
}
