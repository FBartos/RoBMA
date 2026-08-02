double probability_from_log_ratio(double log_mass, double log_denom)
{
  if (log_mass == -INFINITY) {
    return 0;
  }
  if (!std::isfinite(log_mass) || !std::isfinite(log_denom) ||
      log_mass > log_denom) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double probability = std::exp(log_mass - log_denom);
  if (!std::isfinite(probability) ||
      probability < 0 || probability > 1) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return probability;
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
    *ok = false;
    return 0;
  }

  const double tail_lower = tail_at(lower, mean, sd);
  const double tail_upper = tail_at(upper, mean, sd);
  const double prob = tail_lower - tail_upper;
  if (!(prob > 0) || prob > 1) {
    *ok = false;
    return 0;
  }

  return prob;
}

bool add_weighted_standard_normal_interval(double lower, double upper,
                                           double log_weight,
                                           double *log_mass)
{
  if (!(lower < upper)) {
    return true;
  }

  const double log_probability = interval_log_prob(lower, upper, 0, 1);
  if (log_probability == -INFINITY) {
    return true;
  }
  if (!std::isfinite(log_probability)) {
    return false;
  }
  const double log_term = log_weight + log_probability;
  if (log_term == -INFINITY) {
    return true;
  }
  if (!std::isfinite(log_term)) {
    return false;
  }
  *log_mass = add_logspace(*log_mass, log_term);
  return true;
}

double step_log_normalizer_score_space(const double *omega,
                                       const SelNormKernelData &data,
                                       double mean, double sd, double sei,
                                       int omega_stride)
{
  double log_mass = -INFINITY;

  for (int b = 0; b < data.n_bins; ) {
    const double w = omega_at(omega, b + 1, omega_stride);
    int last = b;

    if (data.trusted_step_partition) {
      while (last + 1 < data.n_bins &&
             omega_at(omega, last + 2, omega_stride) == w) {
        ++last;
      }
    }

    if (w > 0) {
      const double lower_z = data.trusted_step_partition ?
        data.z_lower[last] : data.z_lower[b];
      const double upper_z = data.z_upper[b];
      const double lower   = stable_affine_score(lower_z, sei, mean, sd);
      const double upper   = stable_affine_score(upper_z, sei, mean, sd);
      if (std::isnan(lower) || std::isnan(upper) || !(lower < upper) ||
          !add_weighted_standard_normal_interval(
            lower, upper, std::log(w), &log_mass
          )) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    } else if (w < 0 || std::isnan(w)) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    b = data.trusted_step_partition ? last + 1 : b + 1;
  }

  return log_mass;
}

bool normalize_log_tail_pair(double log_lower_mass, double log_upper_mass,
                             double *log_lower, double *log_upper)
{
  if (std::isnan(log_lower_mass) || std::isnan(log_upper_mass) ||
      (log_lower_mass == -INFINITY && log_upper_mass == -INFINITY)) {
    return false;
  }
  if (log_lower_mass == -INFINITY) {
    *log_lower = -INFINITY;
    *log_upper = 0;
    return true;
  }
  if (log_upper_mass == -INFINITY) {
    *log_lower = 0;
    *log_upper = -INFINITY;
    return true;
  }

  *log_lower = -log1pexp(log_upper_mass - log_lower_mass);
  *log_upper = -log1pexp(log_lower_mass - log_upper_mass);
  return std::isfinite(*log_lower) && std::isfinite(*log_upper);
}

bool step_log_tail_pair_logspace(const double *omega,
                                 const SelNormKernelData &data,
                                 double q_score, double mean, double sd,
                                 double sei, int omega_stride,
                                 bool trusted_partition,
                                 double *log_lower, double *log_upper)
{
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
      const double lower_z = trusted_partition ?
        data.z_lower[last] : data.z_lower[b];
      const double upper_z = data.z_upper[b];
      const double lower = stable_affine_score(lower_z, sei, mean, sd);
      const double upper = stable_affine_score(upper_z, sei, mean, sd);
      const double log_w = std::log(w);
      if (std::isnan(lower) || std::isnan(upper) || !(lower < upper) ||
          !add_weighted_standard_normal_interval(
            lower, std::min(upper, q_score), log_w, &log_lower_mass
          ) ||
          !add_weighted_standard_normal_interval(
            std::max(lower, q_score), upper, log_w, &log_upper_mass
          )) {
        return false;
      }
    } else if (w < 0 || std::isnan(w)) {
      return false;
    }

    b = trusted_partition ? last + 1 : b + 1;
  }

  return normalize_log_tail_pair(
    log_lower_mass, log_upper_mass, log_lower, log_upper
  );
}

double step_cdf_logspace_with_log_norm(const double *omega,
                                       const SelNormKernelData &data,
                                       double q_score, double mean, double sd,
                                       double sei,
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
      const double lower_z = trusted_partition ?
        data.z_lower[last] : data.z_lower[b];
      const double upper_z = data.z_upper[b];
      const double lower   = stable_affine_score(lower_z, sei, mean, sd);
      const double upper   = stable_affine_score(upper_z, sei, mean, sd);
      const double log_w = std::log(w);
      if (std::isnan(lower) || std::isnan(upper) || !(lower < upper) ||
          !add_weighted_standard_normal_interval(
            lower, std::min(upper, q_score), log_w, &log_lower_mass
          ) ||
          !add_weighted_standard_normal_interval(
            std::max(lower, q_score), upper, log_w, &log_upper_mass
          )) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    } else if (w < 0 || std::isnan(w)) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    b = trusted_partition ? last + 1 : b + 1;
  }

  const double log_requested_mass = signed_lower_tail ?
    log_lower_mass : log_upper_mass;
  return probability_from_log_ratio(log_requested_mass, log_denom);
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

  if (*signed_second < 0) {
    return false;
  }

  return true;
}

bool step_cdf_pair_trusted_telescoping(const double *omega,
                                       const SelNormKernelData &data,
                                       double q1_z, bool q1_signed_lower_tail,
                                       double q2_z, bool q2_signed_lower_tail,
                                       double signed_mean, double sd,
                                       double sei,
                                       int omega_stride,
                                       double *p1, double *p2,
                                       double *inverse_weight)
{
  if (!std::isfinite(step_log_normalizer_score_space(
        omega, data, signed_mean, sd, sei, omega_stride
      ))) {
    return false;
  }

  const double q_z[2] = {q1_z, q2_z};
  const bool signed_lower_tail[2] = {
    q1_signed_lower_tail, q2_signed_lower_tail
  };
  double q_score[2];
  double tail_q[2];
  double cdf_q[2];

  for (int i = 0; i < 2; ++i) {
    q_score[i] = stable_affine_score(q_z[i], sei, signed_mean, sd);
    tail_q[i]  = tail_at(q_score[i], 0, 1);
    cdf_q[i]   = cdf_at(q_score[i], 0, 1);
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
  for (int b = 0; b < data.n_bins - 1; ++b) {
    const double w_left  = omega_at(omega, b + 1, omega_stride);
    const double w_right = omega_at(omega, b + 2, omega_stride);
    if (w_left == w_right) {
      continue;
    }

    const double diff = w_left - w_right;
    const double boundary_score = stable_affine_score(
      data.z_lower[b], sei, signed_mean, sd
    );
    const double tail_boundary = tail_at(boundary_score, 0, 1);
    const double norm_contrib  = diff * tail_boundary;
    normalizer += norm_contrib;
    norm_abs   += std::abs(norm_contrib);

    for (int i = 0; i < 2; ++i) {
      double req_contrib = 0;
      if (signed_lower_tail[i]) {
        if (q_z[i] > data.z_lower[b]) {
          const double interval = tail_boundary - tail_q[i];
          if (!(interval > 0)) {
            return false;
          }
          req_contrib = diff * interval;
        }
      } else {
        req_contrib = diff * (q_z[i] >= data.z_lower[b] ? tail_q[i] : tail_boundary);
      }

      requested[i] += req_contrib;
    }
  }

  if (!std::isfinite(normalizer) || !(normalizer > 0) ||
      normalizer <= std::numeric_limits<double>::min() ||
      normalizer <= 1e-8 * norm_abs ||
      !std::isfinite(requested[0]) ||
      !std::isfinite(requested[1])) {
    return false;
  }

  if (requested[0] < 0 || requested[0] > normalizer ||
      requested[1] < 0 || requested[1] > normalizer) {
    return false;
  }

  *p1             = requested[0] / normalizer;
  *p2             = requested[1] / normalizer;
  *inverse_weight = 1 / normalizer;
  return std::isfinite(*p1) && *p1 >= 0 && *p1 <= 1 &&
    std::isfinite(*p2) && *p2 >= 0 && *p2 <= 1 &&
    std::isfinite(*inverse_weight) && *inverse_weight > 0;
}

bool step_cdf_pair_logspace(const double *omega,
                            const SelNormKernelData &data,
                            double q1_z, bool q1_signed_lower_tail,
                            double q2_z, bool q2_signed_lower_tail,
                            double signed_mean, double sd, double sei,
                            int omega_stride,
                            bool trusted_partition,
                            double *p1, double *p2,
                            double *inverse_weight)
{
  const bool signed_lower_tail[2] = {
    q1_signed_lower_tail, q2_signed_lower_tail
  };
  const double q_score[2] = {
    stable_affine_score(q1_z, sei, signed_mean, sd),
    stable_affine_score(q2_z, sei, signed_mean, sd)
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
      const double lower_z = trusted_partition ?
        data.z_lower[last] : data.z_lower[b];
      const double upper_z = data.z_upper[b];
      const double lower = stable_affine_score(
        lower_z, sei, signed_mean, sd
      );
      const double upper = stable_affine_score(
        upper_z, sei, signed_mean, sd
      );
      if (std::isnan(lower) || std::isnan(upper) || !(lower < upper)) {
        return false;
      }
      const double log_w = std::log(w);
      const double log_prob = interval_log_prob(lower, upper, 0, 1);

      if (log_prob != -INFINITY) {
        log_denom = add_logspace(log_denom, log_w + log_prob);
      }

      for (int i = 0; i < 2; ++i) {
        const double mass_lower = signed_lower_tail[i] ?
          lower : std::max(lower, q_score[i]);
        const double mass_upper = signed_lower_tail[i] ?
          std::min(upper, q_score[i]) : upper;
        const double log_requested = interval_log_prob(
          mass_lower, mass_upper, 0, 1
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

  *p1 = probability_from_log_ratio(log_mass[0], log_denom);
  *p2 = probability_from_log_ratio(log_mass[1], log_denom);
  *inverse_weight = std::exp(-log_denom);

  return std::isfinite(*p1) && *p1 >= 0 && *p1 <= 1 &&
    std::isfinite(*p2) && *p2 >= 0 && *p2 <= 1 &&
    std::isfinite(*inverse_weight) && *inverse_weight > 0;
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
                                  double signed_mean, double sd, double sei,
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
      const double group_lower = stable_affine_score(
        data.z_lower[last], sei, signed_mean, sd
      );
      const double group_upper = stable_affine_score(
        data.z_upper[b], sei, signed_mean, sd
      );
      const double prob = trusted_interval_prob(
        group_lower, group_upper, 0, 1, &interval_ok
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
    0,
    1,
    u_interval
  );
}

double step_rng_logspace(const double *omega, const SelNormKernelData &data,
                         double signed_mean, double sd, double sei,
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
      const double lower_z = trusted_partition ?
        data.z_lower[last] : data.z_lower[b];
      const double upper_z = data.z_upper[b];
      const double group_lower = stable_affine_score(
        lower_z, sei, signed_mean, sd
      );
      const double group_upper = stable_affine_score(
        upper_z, sei, signed_mean, sd
      );
      if (std::isnan(group_lower) || std::isnan(group_upper) ||
          !(group_lower < group_upper)) {
        return std::numeric_limits<double>::quiet_NaN();
      }
      const double lm = std::log(w) + interval_log_prob(
        group_lower, group_upper, 0, 1
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
    0,
    1,
    u_interval
  );
}
