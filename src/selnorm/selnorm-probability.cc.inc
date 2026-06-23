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

[[maybe_unused]] double source_power_moment(
  double lower, double upper, double mean, double sd,
  double source_lower, double source_upper, int q
)
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

[[maybe_unused]] double dest_power_moment(
  double lower, double upper, double mean, double sd,
  double dest_lower, double dest_upper, int q
)
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
