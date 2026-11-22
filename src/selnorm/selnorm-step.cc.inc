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

  const double prob = interval_prob(lower, upper, mean, sd);
  if (!(prob > 0) || !std::isfinite(prob) || prob > 1) {
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
                                       int omega_stride,
                                       double requested_lower = -INFINITY,
                                       double requested_upper = INFINITY,
                                       double *max_log_weight = nullptr)
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
      const double log_weight = std::log(w);
      if (max_log_weight) *max_log_weight = std::max(*max_log_weight, log_weight);
      const double lower_z = data.trusted_step_partition ?
        data.z_lower[last] : data.z_lower[b];
      const double upper_z = data.z_upper[b];
      const double lower   = stable_affine_score(lower_z, sei, mean, sd);
      const double upper   = stable_affine_score(upper_z, sei, mean, sd);
      if (std::isnan(lower) || std::isnan(upper) || !(lower < upper) ||
          !add_weighted_standard_normal_interval(
            std::max(lower, requested_lower), std::min(upper, requested_upper),
            log_weight, &log_mass
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

bool near_step_cut(double query_score, double cut_score)
{
  return std::isfinite(cut_score) && std::abs(query_score - cut_score) <= 1e-5;
}

// Integrate a thin Gaussian interval using its original outcome width.
bool narrow_normal_log_probability(double width, double sd, double lo, double hi,
                                   double &log_p)
{
  if (!(width > 0) || !std::isfinite(width) || !(sd > 0) || !std::isfinite(sd) ||
      !std::isfinite(lo) || !std::isfinite(hi)) return false;
  const double dz = width / sd;
  const double nearest = std::min(std::abs(lo), std::abs(hi));
  if (!(dz <= 1e-5) && !(lo == hi && nearest >= 1e6)) return false;
  const double mid = .5*lo + .5*hi;
  const double slope = mid * dz;
  // Up to this slope the omitted fourth-order term is at most about
  // (1e-3)^4/1920. Above it, dz <= 1e-5 implies |z| >= 100, so the
  // quadratic remainder below fits the rank-one enclosure's log allowance.
  if (std::abs(slope) <= 1e-3) {
    log_p = std::log(width) - std::log(sd) + dnorm(mid, 0, 1, true) +
      std::log1p((slope * slope - dz * dz) / 24);
  } else {
    // On one side of the mode, integrate the exponential slope at the
    // nearer endpoint exactly. The omitted quadratic has relative effect
    // at most dz^2/2, or 1/nearest^2 for a remote one-sided interval.
    const double decay = nearest * dz;
    log_p = dnorm(nearest, 0, 1, true) + std::log(-std::expm1(-decay)) - std::log(nearest);
  }
  return !std::isnan(log_p);
}

bool add_partial_step_mass(double query, double query_score, double cut, double cut_score,
                           double sd, bool above_cut, double log_weight, double *log_mass)
{
  const double lo = above_cut ? cut_score : query_score;
  const double hi = above_cut ? query_score : cut_score;
  double log_p;
  if (narrow_normal_log_probability(std::abs(query - cut), sd, lo, hi, log_p)) {
    *log_mass = add_logspace(*log_mass, log_weight + log_p);
    return !std::isnan(*log_mass);
  }
  return lo < hi && add_weighted_standard_normal_interval(lo, hi, log_weight, log_mass);
}

bool step_log_tail_pair_logspace(const double *omega,
                                 const SelNormKernelData &data,
                                 double query, double q_score, double mean, double sd,
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
      if (!(lower < upper)) return false;
      if (near_step_cut(q_score, lower) || near_step_cut(q_score, upper)) {
        const double raw_lower = lower_z * sei, raw_upper = upper_z * sei;
        if (query <= raw_lower || query >= raw_upper) {
          double *whole = query <= raw_lower ? &log_upper_mass : &log_lower_mass;
          if (!add_weighted_standard_normal_interval(lower, upper, log_w, whole)) return false;
        } else if (!add_partial_step_mass(query, q_score, raw_lower, lower, sd, true,
                                          log_w, &log_lower_mass) ||
                   !add_partial_step_mass(query, q_score, raw_upper, upper, sd, false,
                                          log_w, &log_upper_mass)) return false;
      } else {
        if (!add_weighted_standard_normal_interval(lower, std::min(upper, q_score), log_w,
                                                   &log_lower_mass) ||
            !add_weighted_standard_normal_interval(std::max(lower, q_score), upper, log_w,
                                                   &log_upper_mass)) return false;
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

bool step_log_tail_pair_from_plan(
  double q_score,
  const double *lower_score,
  const double *upper_score,
  const double *log_weight,
  int n_groups,
  double *log_lower,
  double *log_upper)
{
  double log_lower_mass = -INFINITY;
  double log_upper_mass = -INFINITY;

  for (int g = 0; g < n_groups; ++g) {
    const double lower = lower_score[g];
    const double upper = upper_score[g];
    if (near_step_cut(q_score, lower) || near_step_cut(q_score, upper)) return false;
    const double log_w = log_weight[g];
    if (!add_weighted_standard_normal_interval(
          lower, std::min(upper, q_score), log_w, &log_lower_mass
        ) ||
        !add_weighted_standard_normal_interval(
          std::max(lower, q_score), upper, log_w, &log_upper_mass
        )) {
      return false;
    }
  }

  return normalize_log_tail_pair(
    log_lower_mass, log_upper_mass, log_lower, log_upper
  );
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

bool step_positive_score_partition(const double *omega,
                                    const SelNormKernelData &data,
                                    double mean, double sd, double sei,
                                    int omega_stride)
{
  if (!data.trusted_step_partition || data.n_bins < 1 ||
      !std::isfinite(mean) || !(sd > 0) || !std::isfinite(sd) ||
      !(sei > 0) || !std::isfinite(sei)) {
    return false;
  }

  double previous = INFINITY;
  for (int b = 0; b < data.n_bins; ++b) {
    const double weight = omega_at(omega, b + 1, omega_stride);
    if (!(weight > 0) || !std::isfinite(weight)) return false;
    if (b + 1 < data.n_bins &&
        weight != omega_at(omega, b + 2, omega_stride)) {
      const double score = stable_affine_score(data.z_lower[b], sei, mean, sd);
      if (!std::isfinite(score) || !(score < previous)) return false;
      previous = score;
    }
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
                                       double *inverse_weight,
                                       double *log_normalizer = nullptr)
{
  // A full score-space partition with finite positive weights has finite
  // log mass. Check its coalesced boundaries without evaluating that unused
  // mass; retain the original gate for zero weights or collapsed/extreme scores.
  // The skipped log-normalizer evaluation could also fail through rounding-level
  // non-monotonicity of log pnorm; the telescoped value is then used. The
  // review of this gate found that regime unreachable from fitted models.
  if (!step_positive_score_partition(omega, data, signed_mean, sd, sei, omega_stride) &&
      !std::isfinite(step_log_normalizer_score_space(
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
    // The upper-tail query never consumes its lower tail. For a unit normal,
    // finite tail_q also establishes that the unused lower tail is finite.
    cdf_q[i]   = signed_lower_tail[i] ? cdf_at(q_score[i], 0, 1) : 0;
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
    // Significance queries often coincide exactly with a selection boundary.
    // Reuse only the same score, including its zero sign, so every subsequent
    // multiplication and accumulation remains unchanged.
    const bool first_query = boundary_score == q_score[0] &&
      std::signbit(boundary_score) == std::signbit(q_score[0]);
    const bool second_query = boundary_score == q_score[1] &&
      std::signbit(boundary_score) == std::signbit(q_score[1]);
    const double tail_boundary = first_query ? tail_q[0] :
      second_query ? tail_q[1] : tail_at(boundary_score, 0, 1);
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
  if (log_normalizer != nullptr) *log_normalizer = std::log(normalizer);
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
                            double *inverse_weight,
                            double *log_normalizer = nullptr)
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
  if (log_normalizer != nullptr) *log_normalizer = log_denom;

  return std::isfinite(*p1) && *p1 >= 0 && *p1 <= 1 &&
    std::isfinite(*p2) && *p2 >= 0 && *p2 <= 1 &&
    std::isfinite(*inverse_weight) && *inverse_weight > 0;
}

// A half-bin starts at the point nearest the Gaussian mode. In units of
// h = 1/max(1,a) sampling SDs, its density relative to that endpoint is
// exp(-a*h*s - (h*s)^2/2), even when phi(a) itself has underflowed.
struct StepMomentHalf {
  double score;
  double cut;
  double log_mass;
  double anchor;
  double offset;
  double variance;
};

bool step_half_moments(double a, double width_score, double cut, double anchor, double signed_sd,
                       double log_weight, StepMomentHalf &half)
{
  if (!std::isfinite(a) || !(width_score > 0) || !std::isfinite(anchor)) return false;
  const double h = 1 / std::max(1.0, a);
  const double rate = std::min(a, 1.0);
  // Beyond a log-density drop of 50 the omitted exponential second-moment
  // fraction is at most (1 + 50 + 50^2/2)*exp(-50).
  const double end = std::min(width_score / h,
    100 / (rate + std::hypot(rate, 10*h)));
  if (!(end > 0)) return false;

  // Eight-point Gauss-Legendre on panels at most two scaled units wide.
  // Integrate powers of u=s/end, so even a narrow bin has resolved moments.
  static const double nodes[] = {.1834346424956498, .5255324099163290,
                                  .7966664774136267, .9602898564975362};
  static const double weights[] = {.3626837833783620, .3137066458778873,
                                    .2223810344533745, .1012285362903763};
  const int panels = std::max(1, static_cast<int>(std::ceil(end / 2)));
  double mass = 0, first = 0, second = 0;
  for (int p = 0; p < panels; ++p) {
    for (int j = 0; j < 4; ++j) {
      for (int sign = -1; sign <= 1; sign += 2) {
        const double u = (p + .5 + sign * .5 * nodes[j]) / panels;
        const double s = u * end;
        const double z = h * s;
        const double w = weights[j] * std::exp(-rate*s - .5*z*z);
        mass += w;
        first += w*u;
        second += w*u*u;
      }
    }
  }
  const double location = first / mass;
  const double spread = second / mass - location*location;
  const double width = signed_sd * (h * end);
  half.score = a;
  half.cut = cut;
  half.log_mass = log_weight + std::log(h) + std::log(end) + std::log(mass / (2*panels));
  half.anchor = anchor;
  half.offset = width * location;
  const double half_sd = width * std::sqrt(spread);
  half.variance = half_sd * half_sd;
  return std::isfinite(half.log_mass) && std::isfinite(half.offset) && std::isfinite(half.variance);
}

bool step_rare_moments(const double *omega, const SelNormKernelData &data,
                       double mean, double sd, double sei, int stride,
                       double &result_mean, double &result_variance)
{
  std::vector<StepMomentHalf> halves;
  for (int b = 0; b < data.n_bins; ++b) {
    const double w = omega_at(omega, b + 1, stride);
    const double zu = data.z_upper[b];
    while (data.trusted_step_partition && b + 1 < data.n_bins &&
           omega_at(omega, b + 2, stride) == w) ++b;
    if (w == 0) continue;
    const double zl = data.z_lower[b];
    const double lo = stable_affine_score(zl, sei, mean, sd);
    const double hi = stable_affine_score(zu, sei, mean, sd);
    if (!(lo < hi)) return false;
    const double log_w = std::log(w);
    const double width_score = (zu - zl) * sei / sd;
    if (lo < 0) {
      StepMomentHalf half;
      if (!step_half_moments(std::max(-hi, 0.0), hi < 0 ? width_score : -lo,
                             zu, hi < 0 ? zu * sei : mean,
                             -sd, log_w, half)) return false;
      halves.push_back(half);
    }
    if (hi > 0) {
      StepMomentHalf half;
      if (!step_half_moments(std::max(lo, 0.0), lo > 0 ? width_score : hi,
                             zl, lo > 0 ? zl * sei : mean,
                             sd, log_w, half)) return false;
      halves.push_back(half);
    }
  }
  // Factor out the nearest endpoint density before adding log weights. This
  // retains weight ratios even for equally remote tails on opposite sides.
  if (halves.empty()) return false;
  const auto nearest = std::min_element(halves.begin(), halves.end(),
    [](const StepMomentHalf &x, const StepMomentHalf &y) { return x.score < y.score; });
  double peak = -INFINITY, origin = 0;
  for (auto &half : halves) {
    if (nearest->score == 0) {
      half.log_mass -= (.5 * half.score) * half.score;
    } else if (half.cut != nearest->cut) {
      // Factor the difference of squared distances in the supplied cuts.
      // Opposing standardized scores can round to equal magnitudes even
      // when the two peak densities differ appreciably. FastTwoSum keeps
      // the rounding remainder of the midpoint before centering it on mu.
      double x = .5 * half.cut, y = .5 * nearest->cut;
      if (std::abs(x) < std::abs(y)) std::swap(x, y);
      const double mid = x + y;
      const double remainder = (x - mid) + y;
      const double mid_offset = std::fma(mid, sei, -mean) + remainder * sei;
      half.log_mass -= ((half.cut - nearest->cut) * sei / sd) * (mid_offset / sd);
    }
    if (half.log_mass > peak) {
      peak = half.log_mass;
      origin = half.anchor;
    }
  }
  if (!std::isfinite(peak)) return false;
  double total = 0, first = 0;
  for (auto &half : halves) {
    const double mass = std::exp(half.log_mass - peak);
    if (mass == 0) continue;
    half.offset += half.anchor - origin;
    total += mass;
    first += mass * half.offset;
  }
  const double offset = first / total;
  result_mean = origin + offset;
  result_variance = 0;
  for (const auto &half : halves) {
    const double p = std::exp(half.log_mass - peak) / total;
    if (p == 0) continue;
    const double delta = half.offset - offset;
    result_variance += p * half.variance + (p * delta) * delta;
  }
  return std::isfinite(result_mean) && std::isfinite(result_variance);
}

// Multiply before applying the shared power-of-two scale, without losing a
// small weight or overflowing the unscaled product.
double scaled_moment_product(double weight_fraction, int shift, double value)
{
  int exponent = 0;
  const double fraction = std::frexp(value, &exponent);
  return std::scalbn(weight_fraction * fraction, shift + exponent);
}

template<bool scaled_scores = false>
void step_moments_partition(const double *omega,
                            const SelNormKernelData &data,
                            double signed_mean, double sd, double sei,
                            int omega_stride, bool trusted_partition,
                            double *denom, double *signed_first,
                            double *signed_second,
                            double *centered_first = NULL,
                            double *centered_second = NULL, int scale_exponent = 0,
                            double *operand_sum = nullptr)
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
      const double lower = data.z_lower[trusted_partition ? last : b];
      const double upper = data.z_upper[b];
      double m0, m1, m2, centered_m1, centered_m2;
      if constexpr (scaled_scores) {
        const double a = stable_affine_score(lower, sei, signed_mean, sd);
        const double b = stable_affine_score(upper, sei, signed_mean, sd);
        truncated_raw_moments(a, b, 0, 1, &m0, &m1, &m2);
        int exponent = 0;
        const double fraction = std::frexp(w, &exponent);
        const int shift = exponent - scale_exponent;
        *denom         += scaled_moment_product(fraction, shift, m0);
        *signed_first  += scaled_moment_product(fraction, shift, m1);
        *signed_second += scaled_moment_product(fraction, shift, m2);
        // Scale of all the operands in P, phi(a)-phi(b), and
        // P+a*phi(a)-b*phi(b), before cancellation inside this bin.
        double operands = a > 0 ? tail_at(a, 0, 1) : (b < 0 ? cdf_at(b, 0, 1) : m0);
        if (std::isfinite(a)) operands += (1 + std::abs(a)) * std_phi(a);
        if (std::isfinite(b)) operands += (1 + std::abs(b)) * std_phi(b);
        *operand_sum += scaled_moment_product(fraction, shift, operands);
      } else {
        truncated_raw_moments(lower * sei, upper * sei, signed_mean, sd,
                              &m0, &m1, &m2, &centered_m1, &centered_m2);
        *denom         += w * m0;
        *signed_first  += w * m1;
        *signed_second += w * m2;
        if (centered_first != NULL) *centered_first += w * centered_m1;
        if (centered_second != NULL) *centered_second += w * centered_m2;
      }
    }

    b = trusted_partition ? last + 1 : b + 1;
  }
}

bool scaled_step_moments(const double *omega, const SelNormKernelData &data,
                         double signed_mean, double sd, double sei, int stride,
                         double log_normalizer, double *moment_mean,
                         double *moment_second, double *moment_variance)
{
  // Raw weighted sums can overflow even though dividing by their common
  // normalizer gives ordinary values. Retry in standard normal coordinates.
  const double exponent = std::nearbyint(log_normalizer / M_LN2);
  if (exponent < 2 * (std::numeric_limits<double>::min_exponent -
                      std::numeric_limits<double>::digits) ||
      exponent > 2 * std::numeric_limits<double>::max_exponent) return false;
  double denom, signed_first, signed_second, operands = 0;
  step_moments_partition<true>(omega, data, signed_mean, sd, sei, stride,
    data.trusted_step_partition, &denom, &signed_first, &signed_second,
    nullptr, nullptr, static_cast<int>(exponent), &operands);
  if (!(denom > 0) || !std::isfinite(denom)) return false;
  const double log_error = std::abs((std::log(denom) + exponent * M_LN2) - log_normalizer);
  if (log_error > 64 * std::numeric_limits<double>::epsilon() *
                  std::max(1.0, std::abs(log_normalizer))) return false;
  const double shift = signed_first / denom;
  const double second_score = signed_second / denom;
  const double variance_score = std::fma(-shift, shift, second_score);
  // A finite retry must also resolve the variance. A narrow, dominant bin
  // can lose its digits inside the raw moments before the final subtraction.
  // Limit the operand-to-result ratio as in the fused cell summary.
  const double center_scale = 1 + std::abs(shift);
  if (!(variance_score > 0) || !std::isfinite(operands) ||
      (operands / denom) * center_scale * center_scale > 1e4 * variance_score) return false;
  const double location = std::fma(sd, shift, signed_mean);
  const double spread = sd * std::sqrt(variance_score);
  const double variance = spread * spread;
  const double second = std::fma(location, location, variance);
  if (!std::isfinite(location) || !std::isfinite(variance) ||
      (!moment_variance && !std::isfinite(second))) return false;
  *moment_mean = data.effect_sign * location;
  if (moment_variance) *moment_variance = variance;
  if (std::isfinite(second)) *moment_second = second;
  return true;
}

double step_normalizer_trusted_partition(
    const double *omega, const SelNormKernelData &data,
    double signed_mean, double sd, double sei, int omega_stride,
    double *mass, double *lower, double *upper, int *n_groups, bool *ok)
{
  *ok = false;

  double denom = 0;
  *n_groups = 0;

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
        if (mass != NULL) mass[*n_groups] = group_mass;
        if (lower != NULL) lower[*n_groups] = group_lower;
        if (upper != NULL) upper[*n_groups] = group_upper;
        ++(*n_groups);
        denom += group_mass;
      }
    }

    b = last + 1;
  }

  if (!std::isfinite(denom) || !(denom > 0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  *ok = true;
  return denom;
}


bool step_normalizer_trusted_telescoping(
    const double *omega, const SelNormKernelData &data,
    double signed_mean, double sd, double sei, int omega_stride,
    double *boundary_tail, double *omega_diff, double *omega_last,
    double *normalizer)
{
  // For an ordered partition, sum_b omega_b Pr(I_b) equals the last
  // weight plus the boundary-tail probabilities weighted by adjacent
  // differences. Callers retain the interval-based path as a cancellation
  // fallback, so this changes only how the same normalizer is evaluated.
  *omega_last = omega_at(omega, data.n_bins, omega_stride);
  *normalizer = *omega_last;
  double normalizer_abs = std::abs(*normalizer);

  for (int b = 0; b < data.n_bins - 1; ++b) {
    const double difference = omega_at(omega, b + 1, omega_stride) -
      omega_at(omega, b + 2, omega_stride);
    if (difference == 0.0 && boundary_tail == nullptr) {
      if (omega_diff != nullptr) omega_diff[b] = difference;
      continue;
    }
    const double boundary_score = stable_affine_score(
      data.z_lower[b], sei, signed_mean, sd
    );
    // The standard erfc identity avoids the general location/scale CDF call
    // at each quadrature point and retains a direct upper-tail evaluation.
    const double tail = 0.5 * std::erfc(boundary_score * M_SQRT1_2);
    const double contribution = difference * tail;
    if (boundary_tail != NULL) boundary_tail[b] = tail;
    if (omega_diff != NULL) omega_diff[b] = difference;
    *normalizer    += contribution;
    normalizer_abs += std::abs(contribution);
  }

  return std::isfinite(*omega_last) && std::isfinite(*normalizer) &&
    *normalizer > std::numeric_limits<double>::min() &&
    *normalizer > 1e-8 * normalizer_abs;
}


double step_rng_prepared_partition(double u_bin, double u_interval,
    const double *mass, const double *lower, const double *upper,
    int n_groups, double denom)
{
  if (n_groups < 1 || !(denom > 0.0) || !std::isfinite(denom)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double target = u_bin * denom;
  double cumulative = 0;
  int selected = n_groups - 1;
  for (int g = 0; g < n_groups; ++g) {
    cumulative += mass[g];
    if (target <= cumulative) {
      selected = g;
      break;
    }
  }
  return truncated_normal_quantile(
    lower[selected], upper[selected], 0, 1, u_interval
  );
}

double step_rng_trusted_partition(const double *omega,
                                   const SelNormKernelData &data,
                                   double signed_mean, double sd, double sei,
                                  double u_bin, double u_interval,
                                  int omega_stride,
                                  double *mass, double *lower, double *upper,
                                  double *log_normalizer,
                                  bool *ok)
{
  int n_groups = 0;
  const double denom = step_normalizer_trusted_partition(
    omega, data, signed_mean, sd, sei, omega_stride,
    mass, lower, upper, &n_groups, ok
  );
  if (!*ok) return std::numeric_limits<double>::quiet_NaN();

  if (log_normalizer != NULL) {
    *log_normalizer = std::log(denom);
  }

  return step_rng_prepared_partition(
    u_bin, u_interval, mass, lower, upper, n_groups, denom
  );
}

double step_rng_logspace(const double *omega, const SelNormKernelData &data,
                         double signed_mean, double sd, double sei,
                         double u_bin, double u_interval,
                         int omega_stride, bool trusted_partition,
                         double *log_mass, double *lower, double *upper)
{
  double max_log_mass = -INFINITY;
  int n_groups = 0;

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
        log_mass[n_groups] = lm;
        lower[n_groups]    = group_lower;
        upper[n_groups]    = group_upper;
        ++n_groups;
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
  for (int g = 0; g < n_groups; ++g) {
    denom += std::exp(log_mass[g] - max_log_mass);
  }

  if (!(denom > 0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double target = u_bin * denom;
  double cumulative = 0;
  int selected = n_groups - 1;

  for (int g = 0; g < n_groups; ++g) {
    cumulative += std::exp(log_mass[g] - max_log_mass);
    if (target <= cumulative) {
      selected = g;
      break;
    }
  }

  return truncated_normal_quantile(
    lower[selected],
    upper[selected],
    0,
    1,
    u_interval
  );
}
