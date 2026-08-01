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

double step_log_cdf_logspace(const double *omega,
                             const SelNormKernelData &data,
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
    return -INFINITY;
  }

  return log_requested_mass - log_denom;
}

double step_cdf_logspace(const double *omega, const SelNormKernelData &data,
                         double q_z, double mean_z, double sd_z,
                         int omega_stride, bool signed_lower_tail,
                         bool trusted_partition)
{
  const double log_probability = step_log_cdf_logspace(
    omega, data, q_z, mean_z, sd_z, omega_stride, signed_lower_tail,
    trusted_partition
  );
  if (std::isnan(log_probability)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (log_probability == -INFINITY) {
    return 0;
  }

  return clamp_probability(std::exp(log_probability));
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
                            double *signed_second,
                            double *centered_first = NULL,
                            double *centered_second = NULL)
{
  *denom         = 0;
  *signed_first  = 0;
  *signed_second = 0;
  if (centered_first != NULL) {
    *centered_first = 0;
  }
  if (centered_second != NULL) {
    *centered_second = 0;
  }

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
      double m0, m1, m2, centered_m1, centered_m2;
      truncated_raw_moments(
        (trusted_partition ? data.z_lower[last] : data.z_lower[b]) * sei,
        data.z_upper[b] * sei,
        signed_mean,
        sd,
        &m0,
        &m1,
        &m2,
        &centered_m1,
        &centered_m2
      );

      *denom         += w * m0;
      *signed_first  += w * m1;
      *signed_second += w * m2;
      if (centered_first != NULL) {
        *centered_first += w * centered_m1;
      }
      if (centered_second != NULL) {
        *centered_second += w * centered_m2;
      }
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
