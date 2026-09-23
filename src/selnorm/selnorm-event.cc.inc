// Shared finite-vector event algebra. These helpers are evaluated under the
// supplied conditional Gaussian law; they never choose a conditioning target.

double vector_observed_log_weight(const double *x, int dimension,
    const double *selection_se, const double *omega, const double *z_lower,
    int n_bins, int effect_sign, int vector_rule, const int *obs_bin)
{
  if (vector_rule == SELVECTOR_PRODUCT) {
    double out = 0.0;
    for (int i = 0; i < dimension; ++i) out += std::log(omega[obs_bin[i] - 1]);
    return out;
  }
  for (int bin = 0; bin < n_bins; ++bin) {
    for (int i = 0; i < dimension; ++i) {
      const double value = vector_rule == SELVECTOR_BEST_TWO_SIDED ?
        std::fabs(x[i]) : effect_sign * x[i];
      if (value >= selection_se[i] * z_lower[bin]) return std::log(omega[bin]);
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

struct EventInterval {
  double lower;
  double upper;
};

// A declared one-dimensional Gaussian support permits exact integration by
// splitting its scalar latent variable at the original selection thresholds.
// The loading is an alternative covariance representation, never inferred from
// a numerical eigendecomposition or a small residual variance.
double rank_one_event_log_mass(const double *mean, const double *loading,
    int dimension, const double *selection_se, const double *omega,
    const SelNormKernelData &selection, int kernel_mode, int vector_rule,
    const double *lower, const double *upper)
{
  const double infinity = std::numeric_limits<double>::infinity();
  double lo = -infinity, hi = infinity;
  for (int i = 0; i < dimension; ++i) {
    const double left = lower == nullptr ? -infinity : lower[i];
    const double right = upper == nullptr ? infinity : upper[i];
    if (left > right) return -infinity;
    if (loading[i] == 0.0) {
      if (mean[i] < left || mean[i] > right) return -infinity;
    } else {
      const double first = (left - mean[i]) / loading[i];
      const double second = (right - mean[i]) / loading[i];
      lo = std::max(lo, std::min(first, second));
      hi = std::min(hi, std::max(first, second));
    }
  }
  if (!(lo < hi)) return -infinity;
  std::vector<double> cuts;
  cuts.push_back(lo);
  cuts.push_back(hi);
  const int sign = vector_rule == SELVECTOR_BEST_TWO_SIDED ? 1 : selection.effect_sign;
  if (kernel_mode != SELKERNEL_NORMAL) {
    for (int i = 0; i < dimension; ++i) {
      if (loading[i] == 0.0) continue;
      for (int bin = 0; bin < selection.n_bins - 1; ++bin) {
        const double cut = (sign * selection_se[i] * selection.z_lower[bin] -
          mean[i]) / loading[i];
        if (cut > lo && cut < hi) cuts.push_back(cut);
        if (vector_rule == SELVECTOR_BEST_TWO_SIDED) {
          const double reflected = (-selection_se[i] * selection.z_lower[bin] -
            mean[i]) / loading[i];
          if (reflected > lo && reflected < hi) cuts.push_back(reflected);
        }
      }
    }
  }
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  std::vector<double> outcome(dimension);
  std::vector<int> bins(dimension);
  double result = -infinity;
  for (std::size_t interval = 1; interval < cuts.size(); ++interval) {
    const double probability = cpp_selnorm_normal_interval_log_prob(
      cuts[interval - 1], cuts[interval], 0.0, 1.0);
    if (!std::isfinite(probability)) continue;
    double log_weight = 0.0;
    if (kernel_mode != SELKERNEL_NORMAL) {
      const double latent = cpp_selnorm_normal_interval_quantile(
        cuts[interval - 1], cuts[interval], 0.0, 1.0, .5);
      for (int i = 0; i < dimension; ++i) {
        outcome[i] = mean[i] + loading[i] * latent;
        bins[i] = 0;
        for (int bin = 0; bin < selection.n_bins; ++bin) {
          if (selection.effect_sign * outcome[i] / selection_se[i] >=
              selection.z_lower[bin]) { bins[i] = bin + 1; break; }
        }
        if (bins[i] == 0) return std::numeric_limits<double>::quiet_NaN();
      }
      log_weight = vector_observed_log_weight(outcome.data(), dimension,
        selection_se, omega, selection.z_lower, selection.n_bins,
        selection.effect_sign, vector_rule, bins.data());
    }
    result = log_add_exp(result, probability + log_weight);
  }
  return result;
}

EventInterval event_interval(double lower, double upper, double se,
    int sign, double box_lower, double box_upper)
{
  EventInterval out;
  out.lower = std::max(box_lower, sign == 1 ? se * lower : -se * upper);
  out.upper = std::min(box_upper, sign == 1 ? se * upper : -se * lower);
  return out;
}

double event_interval_log_probability(const EventInterval &interval,
    double mean, double sd, double *derivative = nullptr)
{
  const double value = cpp_selnorm_normal_interval_log_prob(
    interval.lower, interval.upper, mean, sd);
  if (derivative != nullptr) {
    *derivative = 0.0;
    if (std::isfinite(value)) {
      const double lower_density = std::isfinite(interval.lower) ?
        std::exp(dnorm(interval.lower, mean, sd, true) - value) : 0.0;
      const double upper_density = std::isfinite(interval.upper) ?
        std::exp(dnorm(interval.upper, mean, sd, true) - value) : 0.0;
      *derivative = lower_density - upper_density;
    }
  }
  return value;
}

// For each p bin, partition its event by the first row attaining that bin.
// The recursion sums positive probabilities instead of subtracting nearly
// equal joint CDFs. Thus nonmonotone and zero heights need no telescoping.
double independent_best_log_mass(const double *mean, const double *sd,
    int dimension, const double *selection_se, const double *omega,
    const SelNormKernelData &selection, int vector_rule,
    const double *box_lower = nullptr, const double *box_upper = nullptr,
    const double *loading = nullptr, int rank = 0, double *gradient = nullptr)
{
  const double infinity = std::numeric_limits<double>::infinity();
  double result = -infinity;
  if (gradient != nullptr) std::fill_n(gradient, rank, 0.0);
  std::vector<double> small_gradient(rank), band_gradient(rank), next_gradient(rank);
  for (int bin = 0; bin < selection.n_bins; ++bin) {
    if (!(omega[bin] > 0.0)) continue;
    double lo = selection.z_lower[bin];
    const double hi = selection.z_upper[bin];
    const bool two_sided = vector_rule == SELVECTOR_BEST_TWO_SIDED;
    if (two_sided) {
      if (!(hi > 0.0)) continue;
      lo = std::max(0.0, lo);
    }
    double log_small = 0.0;
    double log_band = -infinity;
    std::fill(small_gradient.begin(), small_gradient.end(), 0.0);
    std::fill(band_gradient.begin(), band_gradient.end(), 0.0);
    for (int i = 0; i < dimension; ++i) {
      const double lower = box_lower == nullptr ? -infinity : box_lower[i];
      const double upper = box_upper == nullptr ? infinity : box_upper[i];
      const int sign = two_sided ? 1 : selection.effect_sign;
      const EventInterval small = event_interval(two_sided ? -lo : -infinity,
        lo, selection_se[i], sign, lower, upper);
      const EventInterval big = event_interval(two_sided ? -hi : -infinity,
        hi, selection_se[i], sign, lower, upper);
      const EventInterval band = event_interval(lo, hi, selection_se[i],
        sign, lower, upper);
      double derivative_small = 0.0, derivative_big = 0.0, derivative_band = 0.0;
      const double p_small = event_interval_log_probability(small, mean[i], sd[i],
        gradient == nullptr ? nullptr : &derivative_small);
      const double p_big = event_interval_log_probability(big, mean[i], sd[i],
        gradient == nullptr ? nullptr : &derivative_big);
      double p_band = event_interval_log_probability(band, mean[i], sd[i],
        gradient == nullptr ? nullptr : &derivative_band);
      if (two_sided) {
        const EventInterval opposite = event_interval(-hi, -lo,
          selection_se[i], 1, lower, upper);
        double derivative_opposite = 0.0;
        const double p_opposite = event_interval_log_probability(opposite,
          mean[i], sd[i], gradient == nullptr ? nullptr : &derivative_opposite);
        const double combined = log_add_exp(p_band, p_opposite);
        if (gradient != nullptr && std::isfinite(combined)) {
          derivative_band = std::exp(p_band - combined) * derivative_band +
            std::exp(p_opposite - combined) * derivative_opposite;
        }
        p_band = combined;
      }
      const double old_band = log_band + p_big;
      const double new_band = log_small + p_band;
      const double combined = log_add_exp(old_band, new_band);
      if (gradient != nullptr && std::isfinite(combined)) {
        const double old_weight = std::exp(old_band - combined);
        const double new_weight = std::exp(new_band - combined);
        for (int factor = 0; factor < rank; ++factor) {
          const double coefficient = loading[i + dimension * factor];
          next_gradient[factor] = old_weight * (band_gradient[factor] +
            coefficient * derivative_big) + new_weight *
            (small_gradient[factor] + coefficient * derivative_band);
        }
        band_gradient = next_gradient;
      }
      if (gradient != nullptr) {
        for (int factor = 0; factor < rank; ++factor) {
          small_gradient[factor] += loading[i + dimension * factor] *
            derivative_small;
        }
      }
      log_band = combined;
      log_small += p_small;
    }
    const double term = std::log(omega[bin]) + log_band;
    const double combined = log_add_exp(result, term);
    if (gradient != nullptr && std::isfinite(combined)) {
      const double old_weight = std::exp(result - combined);
      const double new_weight = std::exp(term - combined);
      for (int factor = 0; factor < rank; ++factor) {
        gradient[factor] = old_weight * gradient[factor] +
          new_weight * band_gradient[factor];
      }
    }
    result = combined;
  }
  return result;
}

double gaussian_rectangle_log_particle(const double *mean,
    const std::vector<double> &cholesky, int dimension,
    const double *lower, const double *upper, const double *qmc,
    int points, int scrambles, int scramble, int point,
    std::vector<double> &latent)
{
  double out = 0.0;
  for (int i = 0; i < dimension; ++i) {
    if (!(lower[i] < upper[i])) return -std::numeric_limits<double>::infinity();
    double conditional_mean = mean[i];
    for (int j = 0; j < i; ++j) {
      conditional_mean += cholesky[i + dimension * j] * latent[j];
    }
    const double sd = cholesky[i + dimension * i];
    const double local = cpp_selnorm_normal_interval_log_prob(
      lower[i], upper[i], conditional_mean, sd);
    if (!std::isfinite(local)) return local;
    out += local;
    bool draw_needed = false;
    for (int j = i + 1; j < dimension; ++j) {
      draw_needed = draw_needed || cholesky[j + dimension * i] != 0.0;
    }
    if (draw_needed) {
      const double u = qmc_value(qmc, scrambles, points, 2 * dimension,
        scramble, point, 2 * i + 1);
      const double draw = cpp_selnorm_normal_interval_quantile(
        lower[i], upper[i], conditional_mean, sd, u);
      if (!std::isfinite(draw)) return std::numeric_limits<double>::quiet_NaN();
      latent[i] = (draw - conditional_mean) / sd;
    } else {
      latent[i] = 0.0;
    }
  }
  return out;
}

struct GaussianEventRectangle {
  std::vector<double> lower;
  std::vector<double> upper;
  double log_weight;
};

std::vector<GaussianEventRectangle> best_event_rectangles(int dimension,
    const double *selection_se, const double *omega,
    const SelNormKernelData &selection, int vector_rule,
    const double *box_lower, const double *box_upper)
{
  const double infinity = std::numeric_limits<double>::infinity();
  const bool two_sided = vector_rule == SELVECTOR_BEST_TWO_SIDED;
  std::vector<GaussianEventRectangle> rectangles;
  for (int bin = 0; bin < selection.n_bins; ++bin) {
    if (!(omega[bin] > 0.0)) continue;
    double lo = selection.z_lower[bin];
    const double hi = selection.z_upper[bin];
    if (two_sided) {
      if (!(hi > 0.0)) continue;
      lo = std::max(lo, 0.0);
    }
    for (int winner = 0; winner < dimension; ++winner) {
      for (int side = 0; side < (two_sided ? 2 : 1); ++side) {
        GaussianEventRectangle rectangle;
        rectangle.lower.resize(dimension);
        rectangle.upper.resize(dimension);
        rectangle.log_weight = std::log(omega[bin]);
        bool nonempty = true;
        for (int i = 0; i < dimension; ++i) {
          double lower = two_sided ? -(i < winner ? lo : hi) : -infinity;
          double upper = i < winner ? lo : hi;
          if (i == winner) {
            lower = side == 0 ? lo : -hi;
            upper = side == 0 ? hi : -lo;
          }
          const EventInterval interval = event_interval(lower, upper,
            selection_se[i], two_sided ? 1 : selection.effect_sign,
            box_lower == nullptr ? -infinity : box_lower[i],
            box_upper == nullptr ? infinity : box_upper[i]);
          rectangle.lower[i] = interval.lower;
          rectangle.upper[i] = interval.upper;
          nonempty = nonempty && interval.lower < interval.upper;
        }
        if (nonempty) rectangles.push_back(rectangle);
      }
    }
  }
  return rectangles;
}

double gaussian_event_log_mass(const double *mean,
    const std::vector<double> &cholesky, int dimension,
    const double *selection_se, const double *omega,
    const SelNormKernelData &selection, int kernel_mode, int vector_rule,
    const double *box_lower, const double *box_upper, const double *qmc,
    int points, int scrambles, double *relative_mcse)
{
  const double infinity = std::numeric_limits<double>::infinity();
  *relative_mcse = 0.0;
  std::vector<double> sd(dimension);
  bool diagonal = true;
  for (int i = 0; i < dimension; ++i) {
    sd[i] = cholesky[i + dimension * i];
    for (int j = 0; j < i; ++j) diagonal = diagonal && cholesky[i + dimension * j] == 0.0;
  }
  std::vector<double> heights(selection.n_bins, 1.0);
  if (kernel_mode != SELKERNEL_NORMAL) {
    std::copy_n(omega, selection.n_bins, heights.begin());
  } else {
    vector_rule = SELVECTOR_PRODUCT;
  }
  bool constant = true;
  for (int bin = 1; bin < selection.n_bins; ++bin) {
    constant = constant && heights[bin] == heights[0];
  }
  if (constant) {
    bool unrestricted = true;
    for (int i = 0; i < dimension; ++i) {
      unrestricted = unrestricted &&
        (box_lower == nullptr || box_lower[i] == -infinity) &&
        (box_upper == nullptr || box_upper[i] == infinity);
    }
    if (unrestricted) {
      return std::log(heights[0]) *
        (vector_rule == SELVECTOR_PRODUCT ? dimension : 1);
    }
  }
  if (diagonal && vector_rule != SELVECTOR_PRODUCT && !constant) {
    return independent_best_log_mass(mean, sd.data(), dimension, selection_se,
      heights.data(), selection, vector_rule, box_lower, box_upper);
  }
  std::vector<GaussianEventRectangle> rectangles;
  if (constant) {
    GaussianEventRectangle rectangle;
    rectangle.lower.resize(dimension, -infinity);
    rectangle.upper.resize(dimension, infinity);
    if (box_lower != nullptr) std::copy_n(box_lower, dimension, rectangle.lower.begin());
    if (box_upper != nullptr) std::copy_n(box_upper, dimension, rectangle.upper.begin());
    rectangle.log_weight = std::log(heights[0]) *
      (vector_rule == SELVECTOR_PRODUCT ? dimension : 1);
    rectangles.push_back(rectangle);
  } else if (vector_rule != SELVECTOR_PRODUCT) {
    rectangles = best_event_rectangles(dimension, selection_se, heights.data(),
      selection, vector_rule, box_lower, box_upper);
  }
  std::vector<double> latent(dimension), log_mass(selection.n_bins);
  std::vector<EventInterval> intervals(selection.n_bins);
  std::vector<double> log_means(scrambles, -infinity);
  const int evaluated_points = diagonal ? 1 : points;
  const int evaluated_scrambles = diagonal ? 1 : scrambles;
  for (int scramble = 0; scramble < evaluated_scrambles; ++scramble) {
    double log_sum = -infinity;
    for (int point = 0; point < evaluated_points; ++point) {
      double particle = -infinity;
      if (!rectangles.empty() || vector_rule != SELVECTOR_PRODUCT) {
        for (std::size_t rectangle = 0; rectangle < rectangles.size(); ++rectangle) {
          const GaussianEventRectangle &r = rectangles[rectangle];
          const double probability = gaussian_rectangle_log_particle(mean,
            cholesky, dimension, r.lower.data(), r.upper.data(), qmc,
            points, scrambles, scramble, point, latent);
          particle = log_add_exp(particle, r.log_weight + probability);
        }
      } else {
        particle = 0.0;
        for (int i = 0; i < dimension; ++i) {
          double conditional_mean = mean[i];
          for (int j = 0; j < i; ++j) {
            conditional_mean += cholesky[i + dimension * j] * latent[j];
          }
          double local = -infinity;
          int last_positive = -1;
          for (int bin = 0; bin < selection.n_bins; ++bin) {
            intervals[bin] = event_interval(selection.z_lower[bin],
              selection.z_upper[bin], selection_se[i], selection.effect_sign,
              box_lower == nullptr ? -infinity : box_lower[i],
              box_upper == nullptr ? infinity : box_upper[i]);
            log_mass[bin] = std::log(heights[bin]) +
              event_interval_log_probability(intervals[bin], conditional_mean, sd[i]);
            local = log_add_exp(local, log_mass[bin]);
            if (std::isfinite(log_mass[bin])) last_positive = bin;
          }
          if (!std::isfinite(local)) { particle = local; break; }
          particle += local;
          if (!diagonal && i + 1 < dimension) {
            const double u = qmc_value(qmc, scrambles, points, 2 * dimension,
              scramble, point, 2 * i);
            double cumulative = 0.0;
            int chosen = last_positive;
            for (int bin = 0; bin < selection.n_bins; ++bin) {
              cumulative += std::exp(log_mass[bin] - local);
              if (u < cumulative) { chosen = bin; break; }
            }
            const double interval_u = qmc_value(qmc, scrambles, points,
              2 * dimension, scramble, point, 2 * i + 1);
            const double draw = cpp_selnorm_normal_interval_quantile(
              intervals[chosen].lower, intervals[chosen].upper,
              conditional_mean, sd[i], interval_u);
            latent[i] = (draw - conditional_mean) / sd[i];
          }
        }
      }
      if (std::isnan(particle)) { *relative_mcse = infinity; return particle; }
      log_sum = log_add_exp(log_sum, particle);
    }
    log_means[scramble] = log_sum - std::log(static_cast<double>(evaluated_points));
  }
  if (diagonal) return log_means[0];
  double result = -infinity;
  for (int s = 0; s < scrambles; ++s) result = log_add_exp(result, log_means[s]);
  result -= std::log(static_cast<double>(scrambles));
  if (result == -infinity) return result;
  double squared = 0.0;
  for (int s = 0; s < scrambles; ++s) {
    const double difference = std::expm1(log_means[s] - result);
    squared += difference * difference;
  }
  *relative_mcse = std::sqrt(squared / (scrambles * (scrambles - 1.0)));
  return result;
}
