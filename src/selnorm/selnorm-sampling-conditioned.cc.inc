// Matheron augmentation for full sampling-error conditioning. The retained
// covariance includes any retained true-effect sources. This fragment shares
// the finite Gaussian event algebra in selnorm-mv.cc.
namespace {

int sampling_conditioned_bin(double x, double se,
    const SelNormKernelData &selection, int vector_rule)
{
  const double z = vector_rule == SELVECTOR_BEST_TWO_SIDED ?
    std::fabs(x) / se : selection.effect_sign * x / se;
  for (int bin = 0; bin < selection.n_bins; ++bin) {
    if (z >= selection.z_lower[bin]) return bin;
  }
  return -1;
}

double sampling_conditioned_diagonal_mass(const double *mean,
    const double *diagonal, int dimension, const double *se,
    const double *omega, const SelNormKernelData &selection, int vector_rule,
    const double *lower, const double *upper)
{
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  if (vector_rule == SELVECTOR_PRODUCT) {
    double value = 0.0;
    for (int row = 0; row < dimension; ++row) {
      const double lo = lower == nullptr ? negative_infinity : lower[row];
      const double hi = upper == nullptr ? -negative_infinity : upper[row];
      if (diagonal[row] == 0.0) {
        if (mean[row] < lo || mean[row] > hi) return negative_infinity;
        const int bin = sampling_conditioned_bin(mean[row], se[row], selection,
          vector_rule);
        if (bin < 0) return std::numeric_limits<double>::quiet_NaN();
        value += std::log(omega[bin]);
      } else if (lo == negative_infinity && hi == -negative_infinity) {
        value += cpp_selnorm_step_log_norm(mean[row], std::sqrt(diagonal[row]),
          se[row], omega, selection, 1, false);
      } else {
        double mass = negative_infinity;
        for (int bin = 0; bin < selection.n_bins; ++bin) {
          const double first = selection.effect_sign * selection.z_lower[bin] * se[row];
          const double second = selection.effect_sign * selection.z_upper[bin] * se[row];
          mass = log_add_exp(mass, std::log(omega[bin]) +
            cpp_selnorm_normal_interval_log_prob(std::max(lo, std::min(first, second)),
              std::min(hi, std::max(first, second)), mean[row], std::sqrt(diagonal[row])));
        }
        value += mass;
      }
    }
    return value;
  }
  // Track the best p-value bin directly. Positive disjoint probabilities avoid
  // subtraction of nearly equal joint CDFs and handle deterministic rows.
  std::vector<double> best(selection.n_bins, negative_infinity);
  std::vector<double> probability(selection.n_bins), next(selection.n_bins);
  for (int row = 0; row < dimension; ++row) {
    const double left = lower == nullptr ? negative_infinity : lower[row];
    const double right = upper == nullptr ? -negative_infinity : upper[row];
    std::fill(probability.begin(), probability.end(), negative_infinity);
    if (diagonal[row] == 0.0) {
      if (mean[row] < left || mean[row] > right) return negative_infinity;
      const int bin = sampling_conditioned_bin(mean[row], se[row], selection,
        vector_rule);
      if (bin < 0) return std::numeric_limits<double>::quiet_NaN();
      probability[bin] = 0.0;
    } else {
      const double sd = std::sqrt(diagonal[row]);
      for (int bin = 0; bin < selection.n_bins; ++bin) {
        double lo = selection.z_lower[bin], hi = selection.z_upper[bin];
        if (vector_rule == SELVECTOR_BEST_TWO_SIDED) {
          if (!(hi > 0.0)) continue;
          lo = std::max(0.0, lo);
          probability[bin] = log_add_exp(
            cpp_selnorm_normal_interval_log_prob(std::max(left, lo * se[row]),
              std::min(right, hi * se[row]),
              mean[row], sd),
            cpp_selnorm_normal_interval_log_prob(std::max(left, -hi * se[row]),
              std::min(right, -lo * se[row]),
              mean[row], sd));
        } else {
          const double first = selection.effect_sign * lo * se[row];
          const double second = selection.effect_sign * hi * se[row];
          probability[bin] = cpp_selnorm_normal_interval_log_prob(
            std::max(left, std::min(first, second)), std::min(right, std::max(first, second)), mean[row], sd);
        }
      }
    }
    if (row == 0) {
      best = probability;
      continue;
    }
    std::fill(next.begin(), next.end(), negative_infinity);
    for (int previous = 0; previous < selection.n_bins; ++previous) {
      for (int current = 0; current < selection.n_bins; ++current) {
        const int bin = std::min(previous, current);
        next[bin] = log_add_exp(next[bin], best[previous] + probability[current]);
      }
    }
    best.swap(next);
  }
  double result = negative_infinity;
  for (int bin = 0; bin < selection.n_bins; ++bin) {
    result = log_add_exp(result, best[bin] + std::log(omega[bin]));
  }
  return result;
}

double sampling_conditioned_factor_mass(const double *mean,
    const double *diagonal, const double *loading, int dimension, int rank,
    const double *se, const double *omega, const SelNormKernelData &selection,
    int vector_rule, const double *qmc, int points, int design_points,
    int scrambles, double *relative_mcse, const double *lower, const double *upper)
{
  const double negative_infinity = -std::numeric_limits<double>::infinity();
  std::vector<double> shifted(dimension), logs(scrambles, negative_infinity);
  // Rao-Blackwellize all independent residual errors. For a factor-only law,
  // integrate the final declared factor exactly by its threshold intervals.
  const bool factor_only = std::all_of(diagonal, diagonal + dimension,
    [](double value) { return value == 0.0; });
  const int sampled_rank = factor_only ? rank - 1 : rank;
  for (int scramble = 0; scramble < scrambles; ++scramble) {
    for (int point = 0; point < points; ++point) {
      std::copy_n(mean, dimension, shifted.begin());
      for (int factor = 0; factor < sampled_rank; ++factor) {
        const double probability = qmc[scramble + scrambles *
          (point + design_points * factor)];
        const double z = qnorm(probability, 0.0, 1.0, true, false);
        for (int row = 0; row < dimension; ++row) {
          shifted[row] += loading[row + dimension * factor] * z;
        }
      }
      const double value = factor_only ? rank_one_event_log_mass(shifted.data(),
        loading + dimension * (rank - 1), dimension, se, omega, selection,
        SELKERNEL_STEP, vector_rule, lower, upper) :
        sampling_conditioned_diagonal_mass(shifted.data(), diagonal, dimension,
          se, omega, selection, vector_rule, lower, upper);
      logs[scramble] = log_add_exp(logs[scramble], value);
    }
    logs[scramble] -= std::log(static_cast<double>(points));
  }
  double result = negative_infinity;
  for (double value : logs) result = log_add_exp(result, value);
  result -= std::log(static_cast<double>(scrambles));
  *relative_mcse = std::numeric_limits<double>::infinity();
  if (std::isfinite(result)) {
    double squared = 0.0;
    for (double value : logs) {
      const double relative = std::expm1(value - result);
      squared += relative * relative;
    }
    *relative_mcse = std::sqrt(squared / (scrambles * (scrambles - 1.0)));
  }
  return result;
}

// With two declared Gaussian factors and no residual error, integrate the
// second factor exactly and apply the existing one-dimensional Gaussian rules
// to the first. Partition at every change in the ordering of affine selection
// boundaries, including jumps from rows unaffected by the second factor.
double sampling_conditioned_rank_two_mass(const double *mean,
    const double *loading, int dimension, const double *se, const double *omega,
    const SelNormKernelData &selection, int vector_rule, const double *lower,
    const double *upper, const SelNormConditionedQuadrature &quadrature,
    double tolerance, double *relative_change)
{
  const double infinity = std::numeric_limits<double>::infinity();
  struct Line { double intercept; double slope; };
  std::vector<Line> lines;
  std::vector<double> cuts = {-infinity, infinity};
  const int sign = vector_rule == SELVECTOR_BEST_TWO_SIDED ? 1 : selection.effect_sign;
  for (int row = 0; row < dimension; ++row) {
    std::vector<double> boundaries;
    for (int bin = 0; bin < selection.n_bins; ++bin) {
      if (std::isfinite(selection.z_lower[bin])) boundaries.push_back(sign * se[row] * selection.z_lower[bin]);
    }
    if (lower != nullptr && std::isfinite(lower[row])) boundaries.push_back(lower[row]);
    if (upper != nullptr && std::isfinite(upper[row])) boundaries.push_back(upper[row]);
    const double first = loading[row], second = loading[row + dimension];
    for (double boundary : boundaries) {
      if (second != 0.0) {
        lines.push_back({(boundary - mean[row]) / second, -first / second});
      } else if (first != 0.0) {
        const double cut = (boundary - mean[row]) / first;
        if (std::isfinite(cut)) cuts.push_back(cut);
      }
    }
  }
  for (std::size_t first = 0; first < lines.size(); ++first) {
    for (std::size_t second = 0; second < first; ++second) {
      if (lines[first].slope == lines[second].slope) continue;
      const double cut = (lines[first].intercept - lines[second].intercept) /
        (lines[second].slope - lines[first].slope);
      if (std::isfinite(cut)) cuts.push_back(cut);
    }
  }
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  std::vector<double> shifted(dimension);
  std::vector<bool> active_interval(cuts.size() - 1, false);
  for (std::size_t interval = 1; interval < cuts.size(); ++interval) {
    const double latent = cpp_selnorm_normal_interval_quantile(
      cuts[interval - 1], cuts[interval], 0.0, 1.0, .5);
    if (!std::isfinite(latent)) return std::numeric_limits<double>::quiet_NaN();
    for (int row = 0; row < dimension; ++row) shifted[row] = mean[row] + loading[row] * latent;
    const double mass = rank_one_event_log_mass(shifted.data(), loading + dimension,
      dimension, se, omega, selection, SELKERNEL_STEP, vector_rule, lower, upper);
    if (std::isnan(mass)) return mass;
    // Between consecutive crossings the boundary ordering and deterministic
    // row bins are fixed. Thus positive conditional support cannot appear or
    // disappear inside an interval; this certifies empty events exactly.
    active_interval[interval - 1] = mass != -infinity;
  }
  if (std::none_of(active_interval.begin(), active_interval.end(), [](bool value) { return value; })) {
    *relative_change = 0.0;
    return -infinity;
  }
  int offset = 0;
  double previous = std::numeric_limits<double>::quiet_NaN(), previous_change = infinity;
  *relative_change = infinity;
  for (int rule = 0; rule < quadrature.cluster_rule_count; ++rule) {
    const int order = static_cast<int>(quadrature.cluster_orders[rule]);
    double current = -infinity;
    for (std::size_t interval = 1; interval < cuts.size(); ++interval) {
      if (!active_interval[interval - 1]) continue;
      const double lo = cuts[interval - 1], hi = cuts[interval];
      const double probability = cpp_selnorm_normal_interval_log_prob(lo, hi, 0.0, 1.0);
      if (probability == -infinity) continue;
      double integral = -infinity;
      for (int point = 0; point < order; ++point) {
        const double node = quadrature.cluster_nodes[offset + point];
        // Use the smaller tail directly, avoiding subtraction from one.
        const double p = pnorm(std::fabs(node), 0.0, 1.0, false, false);
        if (!(p > 0.0)) return std::numeric_limits<double>::quiet_NaN();
        const double latent = node <= 0.0 ?
          cpp_selnorm_normal_interval_quantile(lo, hi, 0.0, 1.0, p) :
          -cpp_selnorm_normal_interval_quantile(-hi, -lo, 0.0, 1.0, p);
        if (!std::isfinite(latent)) return std::numeric_limits<double>::quiet_NaN();
        for (int row = 0; row < dimension; ++row) shifted[row] = mean[row] + loading[row] * latent;
        const double mass = rank_one_event_log_mass(shifted.data(), loading + dimension,
          dimension, se, omega, selection, SELKERNEL_STEP, vector_rule, lower, upper);
        integral = log_add_exp(integral, quadrature.cluster_log_weights[offset + point] + mass);
      }
      current = log_add_exp(current, probability + integral);
    }
    *relative_change = cluster_relative_change(previous, current);
    if (std::isfinite(current) && *relative_change <= tolerance && previous_change <= tolerance) return current;
    previous = current;
    previous_change = *relative_change;
    offset += order;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double sampling_conditioned_normalizer(const double *mean,
    const double *diagonal, const double *loading, int dimension, int rank,
    const double *se, const double *omega, const SelNormKernelData &selection,
    int vector_rule, const double *qmc, int initial_points, int max_points,
    int scrambles, double tolerance, double *relative_mcse,
    double *relative_change, const SelNormConditionedQuadrature &quadrature,
    const double *lower, const double *upper)
{
  *relative_mcse = 0.0;
  *relative_change = 0.0;
  bool constant = true;
  for (int bin = 1; bin < selection.n_bins; ++bin) {
    constant = constant && omega[bin] == omega[0];
  }
  if (constant && lower == nullptr && upper == nullptr) return std::log(omega[0]) *
    (vector_rule == SELVECTOR_PRODUCT ? dimension : 1);

  // Remove only exactly zero declared columns, including inactive mixture
  // branches. Never infer a smaller statistical rank from eigenvalues.
  std::vector<double> active_loading, residual(diagonal, diagonal + dimension);
  for (int factor = 0; factor < rank; ++factor) {
    int support = 0, single_row = -1;
    for (int row = 0; row < dimension; ++row) {
      if (loading[row + dimension * factor] != 0.0) { ++support; single_row = row; }
    }
    if (support == 1) {
      const double value = loading[single_row + dimension * factor];
      residual[single_row] += value * value;
    } else if (support > 1) active_loading.insert(active_loading.end(),
      loading + dimension * factor, loading + dimension * (factor + 1));
  }
  diagonal = residual.data();
  rank = static_cast<int>(active_loading.size()) / dimension;
  if (rank == 0) return sampling_conditioned_diagonal_mass(mean, diagonal,
    dimension, se, omega, selection, vector_rule, lower, upper);
  const bool factor_only = std::all_of(diagonal, diagonal + dimension,
    [](double value) { return value == 0.0; });
  if (factor_only && rank == 1) return rank_one_event_log_mass(mean,
    active_loading.data(), dimension, se, omega, selection, SELKERNEL_STEP,
    vector_rule, lower, upper);
  if (factor_only && rank == 2) {
    const double result = sampling_conditioned_rank_two_mass(mean,
      active_loading.data(), dimension, se, omega, selection, vector_rule,
      lower, upper, quadrature, tolerance, relative_change);
    if (std::isfinite(result) || result == -std::numeric_limits<double>::infinity()) return result;
  }

  const bool positive_residual = std::all_of(diagonal, diagonal + dimension,
    [](double value) { return value > 0.0; });
  if (positive_residual && rank <= std::min(4, dimension) &&
      lower == nullptr && upper == nullptr) {
    std::vector<double> sd(dimension);
    std::vector<int> bins(dimension, 1);
    for (int row = 0; row < dimension; ++row) sd[row] = std::sqrt(diagonal[row]);
    double normalizer = std::numeric_limits<double>::quiet_NaN();
    if (rank == 1) {
      cpp_selnorm_cluster_step_lpdf(mean, mean, sd.data(), active_loading.data(),
        dimension, se, omega, selection.n_bins, selection.z_lower, selection.z_upper,
        bins.data(), selection.effect_sign, selection.telescope_probabilities,
        SELKERNEL_STEP, quadrature.cluster_nodes, quadrature.cluster_log_weights,
        quadrature.cluster_orders, quadrature.cluster_rule_count, qmc,
        initial_points, max_points, scrambles, tolerance, relative_mcse,
        relative_change, nullptr, &normalizer, vector_rule);
    } else {
      cpp_selnorm_factor_step_lpdf(mean, mean, sd.data(), active_loading.data(),
        dimension, rank, se, omega, selection.n_bins, selection.z_lower,
        selection.z_upper, bins.data(), selection.effect_sign,
        selection.telescope_probabilities, SELKERNEL_STEP, quadrature.factor_nodes,
        quadrature.factor_log_weights, quadrature.factor_orders,
        quadrature.factor_rule_count, quadrature.factor_rule_counts,
        qmc, initial_points, max_points, scrambles, tolerance, relative_mcse,
        relative_change, &normalizer, vector_rule);
    }
    return normalizer;
  }

  std::vector<double> covariance(dimension * dimension, 0.0);
  for (int row = 0; row < dimension; ++row) covariance[row + dimension * row] = diagonal[row];
  for (int factor = 0; factor < rank; ++factor) {
    for (int column = 0; column < dimension; ++column) {
      for (int row = 0; row < dimension; ++row) {
        covariance[row + dimension * column] += active_loading[row + dimension * factor] *
          active_loading[column + dimension * factor];
      }
    }
  }
  // Positive residual variances or a complete triangular loading support
  // certify row rank without numerical rank estimation. Otherwise retain the
  // declared factor law even if roundoff makes a singular Gram matrix pass chol.
  std::vector<bool> first_pivot(dimension), last_pivot(dimension);
  for (int row = 0; row < dimension; ++row) {
    first_pivot[row] = last_pivot[row] = diagonal[row] > 0.0;
  }
  for (int factor = 0; factor < rank; ++factor) {
    int first = -1, last = -1;
    for (int row = 0; row < dimension; ++row) {
      if (active_loading[row + dimension * factor] != 0.0) {
        if (first == -1) first = row;
        last = row;
      }
    }
    first_pivot[first] = true;
    last_pivot[last] = true;
  }
  bool use_dense = std::all_of(first_pivot.begin(), first_pivot.end(),
    [](bool value) { return value; }) ||
    std::all_of(last_pivot.begin(), last_pivot.end(), [](bool value) { return value; });
  if (use_dense) {
    int info = 0;
    F77_CALL(dpotrf)("L", &dimension, covariance.data(), &dimension, &info FCONE);
    if (info != 0) use_dense = false;
  }
  auto evaluate = [&](int points, double *error) {
    if (!use_dense) return sampling_conditioned_factor_mass(mean, diagonal,
      active_loading.data(), dimension, rank, se, omega, selection, vector_rule,
      qmc, points, max_points, scrambles, error, lower, upper);
    // The existing event evaluator uses a contiguous design for each budget.
    std::vector<double> prefix(static_cast<std::size_t>(scrambles) * points * 2 * dimension);
    for (int stream = 0; stream < 2 * dimension; ++stream) {
      std::copy_n(qmc + static_cast<std::size_t>(scrambles) * max_points * stream,
        scrambles * points, prefix.data() + static_cast<std::size_t>(scrambles) * points * stream);
    }
    return gaussian_event_log_mass(mean, covariance, dimension, se, omega,
      selection, SELKERNEL_STEP, vector_rule, lower, upper, prefix.data(),
      points, scrambles, error);
  };
  double coarse_error = 0.0;
  double coarse = evaluate(std::max(1, initial_points / 2), &coarse_error);
  int points = initial_points;
  for (;;) {
    const double result = evaluate(points, relative_mcse);
    *relative_change = std::isfinite(result) && std::isfinite(coarse) ?
      std::fabs(std::expm1(coarse - result)) : std::numeric_limits<double>::infinity();
    if ((*relative_mcse <= tolerance && *relative_change <= tolerance) || points == max_points) {
      return result;
    }
    coarse = result;
    points = std::min(max_points, 2 * points);
  }
}

} // namespace

double cpp_selnorm_conditioned_log_normalizer(
    const double *mean, const double *diagonal, const double *loading,
    int dimension, int rank, const double *selection_se, const double *omega,
    int n_bins, const double *z_lower, const double *z_upper, int effect_sign,
    int kernel_mode, bool telescope_probabilities, int vector_rule,
    const int *group_index, const double *qmc, int initial_points, int max_points,
    int scrambles, double tolerance, double *relative_mcse, double *relative_change,
    const SelNormConditionedQuadrature &quadrature,
    const double *lower, const double *upper)
{
  *relative_mcse = 0.0;
  *relative_change = 0.0;
  if (kernel_mode == SELKERNEL_NORMAL && lower == nullptr && upper == nullptr) return 0.0;
  std::vector<double> normal_weights;
  if (kernel_mode == SELKERNEL_NORMAL) {
    normal_weights.assign(n_bins, 1.0);
    omega = normal_weights.data();
    vector_rule = SELVECTOR_PRODUCT;
  }
  SelNormKernelData selection = {};
  selection.n_bins = n_bins;
  selection.effect_sign = effect_sign;
  selection.z_lower = z_lower;
  selection.z_upper = z_upper;
  selection.trusted_step_partition = true;
  selection.telescope_probabilities = telescope_probabilities;
  if (vector_rule == SELVECTOR_PRODUCT) return sampling_conditioned_normalizer(
    mean, diagonal, loading, dimension, rank, selection_se, omega, selection,
    vector_rule, qmc, initial_points, max_points, scrambles, tolerance,
    relative_mcse, relative_change, quadrature, lower, upper);
  std::vector<int> groups;
  for (int row = 0; row < dimension; ++row) {
    if (std::find(groups.begin(), groups.end(), group_index[row]) == groups.end()) {
      groups.push_back(group_index[row]);
    }
    for (int column = 0; column < row; ++column) {
      if (group_index[row] == group_index[column]) continue;
      double covariance = 0.0;
      for (int factor = 0; factor < rank; ++factor) {
        covariance += loading[row + dimension * factor] * loading[column + dimension * factor];
      }
      if (covariance != 0.0) throw std::runtime_error(
        "Best-rule selection is unavailable because integrated sources connect publication groups.");
    }
  }
  double result = 0.0;
  for (int group : groups) {
    std::vector<int> rows;
    for (int row = 0; row < dimension; ++row) if (group_index[row] == group) rows.push_back(row);
    const int size = static_cast<int>(rows.size());
    std::vector<double> m(size), d(size), l(size * rank), se(size), lo(size), hi(size);
    for (int i = 0; i < size; ++i) {
      m[i] = mean[rows[i]];
      d[i] = diagonal[rows[i]];
      se[i] = selection_se[rows[i]];
      if (lower != nullptr) lo[i] = lower[rows[i]];
      if (upper != nullptr) hi[i] = upper[rows[i]];
      for (int factor = 0; factor < rank; ++factor) l[i + size * factor] = loading[rows[i] + dimension * factor];
    }
    double error = 0.0, change = 0.0;
    result += sampling_conditioned_normalizer(m.data(), d.data(), l.data(), size,
      rank, se.data(), omega, selection, vector_rule, qmc, initial_points,
      max_points, scrambles, tolerance, &error, &change, quadrature,
      lower == nullptr ? nullptr : lo.data(), upper == nullptr ? nullptr : hi.data());
    *relative_mcse = std::max(*relative_mcse, error);
    *relative_change = std::max(*relative_change, change);
  }
  return result;
}

double cpp_selnorm_sampling_conditioned_lpdf(
    const double *x, const double *mean, const double *sampling_lower,
    const double *diagonal, const double *loading, int dimension, int rank,
    const double *sampling_auxiliary, const double *random_auxiliary,
    const double *selection_se, const double *omega, int n_bins,
    const double *z_lower, const double *z_upper, const int *obs_bin,
    int effect_sign, int kernel_mode, bool telescope_probabilities,
    int vector_rule, const int *group_index, const double *qmc,
    int initial_points, int max_points, int scrambles, double tolerance,
    double *relative_mcse, double *relative_change, double *sampling_effect,
    double *delta, double *log_normalizer,
    const SelNormConditionedQuadrature &quadrature)
{
  *relative_mcse = 0.0;
  *relative_change = 0.0;
  *log_normalizer = 0.0;
  std::vector<double> sampling(dimension * dimension), total;
  bool deterministic = true;
  int position = 0;
  for (int column = 0; column < dimension; ++column) {
    for (int row = column; row < dimension; ++row) {
      sampling[row + dimension * column] = sampling_lower[position];
      sampling[column + dimension * row] = sampling_lower[position++];
    }
  }
  total = sampling;
  for (int row = 0; row < dimension; ++row) {
    total[row + dimension * row] += diagonal[row];
    deterministic = deterministic && diagonal[row] == 0.0;
    delta[row] = x[row] - mean[row] - sampling_auxiliary[row] - random_auxiliary[row];
  }
  for (int factor = 0; factor < rank; ++factor) {
    for (int row = 0; row < dimension; ++row) {
      deterministic = deterministic && loading[row + dimension * factor] == 0.0;
      for (int column = 0; column < dimension; ++column) {
        total[row + dimension * column] += loading[row + dimension * factor] *
          loading[column + dimension * factor];
      }
    }
  }
  int info = 0, one = 1;
  F77_CALL(dpotrf)("L", &dimension, total.data(), &dimension, &info FCONE);
  if (info != 0) throw std::runtime_error("The total outcome covariance is not positive definite.");
  F77_CALL(dpotrs)("L", &dimension, &one, total.data(), &dimension,
    delta, &dimension, &info FCONE);
  if (info != 0) throw std::runtime_error("The total outcome covariance solve failed.");
  std::vector<double> normalizer_mean(dimension), whitened(dimension);
  double log_density = -.5 * dimension * std::log(2.0 * std::acos(-1.0));
  for (int row = 0; row < dimension; ++row) {
    sampling_effect[row] = deterministic ? x[row] - mean[row] : sampling_auxiliary[row];
    if (!deterministic) {
      for (int column = 0; column < dimension; ++column) {
        sampling_effect[row] += sampling[row + dimension * column] * delta[column];
      }
    }
    normalizer_mean[row] = mean[row] + sampling_effect[row];
    double value = x[row] - mean[row];
    for (int column = 0; column < row; ++column) {
      value -= total[row + dimension * column] * whitened[column];
    }
    whitened[row] = value / total[row + dimension * row];
    log_density -= std::log(total[row + dimension * row]) + .5 * whitened[row] * whitened[row];
  }
  if (kernel_mode == SELKERNEL_NORMAL) return log_density;
  bool hard_zero = false;
  for (int bin = 0; bin < n_bins; ++bin) {
    hard_zero = hard_zero || (omega[bin] == 0.0 &&
      (vector_rule != SELVECTOR_BEST_TWO_SIDED || z_upper[bin] > 0.0));
  }
  if (hard_zero && !deterministic) {
    std::vector<bool> row_varies(dimension);
    for (int row = 0; row < dimension; ++row) {
      row_varies[row] = diagonal[row] > 0.0;
      for (int factor = 0; factor < rank; ++factor) {
        row_varies[row] = row_varies[row] || loading[row + dimension * factor] != 0.0;
      }
      bool group_varies = row_varies[row];
      if (vector_rule != SELVECTOR_PRODUCT) {
        for (int other = 0; other < dimension; ++other) {
          if (group_index[other] != group_index[row]) continue;
          group_varies = group_varies || diagonal[other] > 0.0;
          for (int factor = 0; factor < rank; ++factor) {
            group_varies = group_varies || loading[other + dimension * factor] != 0.0;
          }
        }
      }
      if (!group_varies) throw std::runtime_error(
        "Conditional sampling selection is unavailable because a fully retained outcome has zero acceptance probability. Use strictly positive selection weights or integrate an outcome-generating source.");
    }
  }
  if (deterministic) {
    for (int bin = 0; bin < n_bins; ++bin) {
      if (omega[bin] == 0.0 &&
          (vector_rule != SELVECTOR_BEST_TWO_SIDED || z_upper[bin] > 0.0)) {
        throw std::runtime_error("Conditional sampling selection is unavailable because a fully retained outcome has zero acceptance probability. Use strictly positive selection weights or integrate an outcome-generating source.");
      }
    }
    *log_normalizer = cpp_selnorm_conditioned_log_normalizer(x, diagonal,
      loading, dimension, rank, selection_se, omega, n_bins, z_lower, z_upper,
      effect_sign, kernel_mode, telescope_probabilities, vector_rule, group_index,
      qmc, initial_points, max_points, scrambles, tolerance, relative_mcse,
      relative_change, quadrature);
    return log_density;
  }
  *log_normalizer = cpp_selnorm_conditioned_log_normalizer(normalizer_mean.data(),
    diagonal, loading, dimension, rank, selection_se, omega, n_bins, z_lower,
    z_upper, effect_sign, kernel_mode, telescope_probabilities, vector_rule,
    group_index, qmc, initial_points, max_points, scrambles, tolerance,
    relative_mcse, relative_change, quadrature);
  if (*log_normalizer == -std::numeric_limits<double>::infinity() && *relative_mcse == 0.0) {
    throw std::runtime_error("Conditional sampling selection is unavailable because a retained context has zero acceptance probability. Use strictly positive selection weights or integrate an outcome-generating source.");
  }
  log_density -= *log_normalizer;
  std::vector<int> groups;
  if (vector_rule == SELVECTOR_PRODUCT) groups.push_back(0);
  else {
    for (int row = 0; row < dimension; ++row) {
      if (std::find(groups.begin(), groups.end(), group_index[row]) == groups.end()) {
        groups.push_back(group_index[row]);
      }
    }
  }
  for (int group : groups) {
    std::vector<int> rows;
    for (int row = 0; row < dimension; ++row) {
      if (vector_rule == SELVECTOR_PRODUCT || group_index[row] == group) rows.push_back(row);
    }
    const int size = static_cast<int>(rows.size());
    std::vector<double> local_se(size), local_x(size);
    std::vector<int> local_bin(size);
    for (int i = 0; i < size; ++i) {
      const int row = rows[i];
      local_se[i] = selection_se[row];
      local_x[i] = x[row];
      local_bin[i] = obs_bin[row];
    }
    log_density += vector_observed_log_weight(local_x.data(), size, local_se.data(),
      omega, z_lower, n_bins, effect_sign, vector_rule, local_bin.data());
  }
  return log_density;
}

void cpp_selnorm_factor_projection(
    const double *mean, const double *diagonal, const double *loading,
    int dimension, int rank, const double *selection_se, const double *omega,
    int n_bins, const double *z_lower, const double *z_upper, int effect_sign,
    int kernel_mode, bool telescope_probabilities, int vector_rule,
    const int *group_index, const double *qmc, int initial_points, int max_points,
    int scrambles, double tolerance, const SelNormConditionedQuadrature &quadrature,
    const double *z, int size, int kind, double *density,
    double *log_normalizer, double *relative_error)
{
  const double infinity = std::numeric_limits<double>::infinity();
  std::fill_n(density, size, 0.0);
  *log_normalizer = 0.0;
  *relative_error = 0.0;
  std::vector<double> absolute_error(size, 0.0);
  std::vector<int> groups;
  if (vector_rule == SELVECTOR_PRODUCT || kernel_mode == SELKERNEL_NORMAL) groups.push_back(0);
  else for (int row = 0; row < dimension; ++row) {
    if (std::find(groups.begin(), groups.end(), group_index[row]) == groups.end()) groups.push_back(group_index[row]);
    for (int column = 0; column < row; ++column) {
      if (group_index[row] == group_index[column]) continue;
      double covariance = 0.0;
      for (int factor = 0; factor < rank; ++factor) covariance +=
        loading[row + dimension * factor] * loading[column + dimension * factor];
      if (covariance != 0.0) throw std::runtime_error(
        "Best-rule selection is unavailable because integrated sources connect publication groups.");
    }
  }
  for (int group : groups) {
    std::vector<int> rows;
    for (int row = 0; row < dimension; ++row) {
      if (vector_rule == SELVECTOR_PRODUCT || kernel_mode == SELKERNEL_NORMAL || group_index[row] == group) rows.push_back(row);
    }
    const int k = static_cast<int>(rows.size());
    std::vector<double> m(k), d(k), l(k * rank), se(k);
    std::vector<int> one_group(k, 1);
    for (int row = 0; row < k; ++row) {
      m[row] = mean[rows[row]];
      d[row] = diagonal[rows[row]];
      se[row] = selection_se[rows[row]];
      for (int column = 0; column < rank; ++column) l[row + k * column] = loading[rows[row] + dimension * column];
    }
    double error = 0.0, change = 0.0;
    const double denominator = cpp_selnorm_conditioned_log_normalizer(m.data(),
      d.data(), l.data(), k, rank, se.data(), omega, n_bins, z_lower, z_upper,
      effect_sign, kernel_mode, telescope_probabilities, vector_rule,
      one_group.data(), qmc, initial_points, max_points, scrambles, tolerance,
      &error, &change, quadrature);
    *log_normalizer += denominator;
    *relative_error = std::max(*relative_error, std::max(error, change));
    if (!std::isfinite(denominator)) {
      *relative_error = infinity;
      std::fill_n(density, size, std::numeric_limits<double>::quiet_NaN());
      return;
    }
    // Include only explicitly nonzero Gaussian coordinates. Diagonal errors
    // become independent columns; this changes representation, not covariance.
    std::vector<double> factor;
    for (int row = 0; row < k; ++row) {
      if (d[row] == 0.0) continue;
      const int start = static_cast<int>(factor.size());
      factor.resize(start + k, 0.0);
      factor[start + row] = std::sqrt(d[row]);
    }
    for (int column = 0; column < rank; ++column) {
      bool active = false;
      for (int row = 0; row < k; ++row) active = active || l[row + k * column] != 0.0;
      if (active) factor.insert(factor.end(), l.data() + k * column, l.data() + k * (column + 1));
    }
    const int coordinates = static_cast<int>(factor.size()) / k;
    bool unselected = kernel_mode == SELKERNEL_NORMAL;
    if (!unselected) unselected = std::all_of(omega, omega + n_bins,
      [omega](double value) { return value == omega[0]; });
    for (int focal = 0; focal < k; ++focal) {
      double sd = 0.0;
      int pivot = 0;
      for (int column = 0; column < coordinates; ++column) {
        sd = std::hypot(sd, factor[focal + k * column]);
        if (std::fabs(factor[focal + k * column]) > std::fabs(factor[focal + k * pivot])) pivot = column;
      }
      // Its marginal is a point mass conditional on the retained context. The
      // outer mixture caller adds this row's ordinary retained Gaussian law.
      if (sd == 0.0) continue;
      std::vector<double> conditional_loading, slope(k), conditional_mean(k), zero(k, 0.0);
      if (kind == 0 && !unselected) {
        std::vector<int> columns;
        columns.push_back(pivot);
        for (int column = 0; column < coordinates; ++column) if (column != pivot) columns.push_back(column);
        std::vector<double> reflection(coordinates);
        for (int column = 0; column < coordinates; ++column) reflection[column] = factor[focal + k * columns[column]] / sd;
        reflection[0] += std::copysign(1.0, reflection[0]);
        double norm = 0.0;
        for (double value : reflection) norm = std::hypot(norm, value);
        for (double &value : reflection) value /= norm;
        conditional_loading.resize(k * (coordinates - 1));
        for (int row = 0; row < k; ++row) {
          double dot = 0.0, covariance = 0.0;
          for (int column = 0; column < coordinates; ++column) {
            dot += factor[row + k * columns[column]] * reflection[column];
            covariance += factor[row + k * column] * factor[focal + k * column];
          }
          slope[row] = (covariance / sd) / sd;
          for (int column = 1; column < coordinates; ++column) {
            conditional_loading[row + k * (column - 1)] = row == focal ? 0.0 :
              factor[row + k * columns[column]] - 2.0 * dot * reflection[column];
          }
        }
      }
      std::vector<double> lower(k, -infinity), upper(k, infinity);
      for (int index = 0; index < size; ++index) {
        const double q = z[index] * se[focal];
        double value = 0.0, local_error = 0.0;
        if (kind == 1 && !std::isfinite(z[index])) {
          value = z[index] > 0.0 ? 1.0 : 0.0;
        } else if (kind == 2 && z[index] == 0.0) {
          value = 1.0;
        } else if (unselected) {
          value = kind == 0 ? se[focal] * dnorm(q, m[focal], sd, false) :
            kind == 1 ? pnorm(q, m[focal], sd, true, false) :
            pnorm(-q, m[focal], sd, true, false) + pnorm(q, m[focal], sd, false, false);
        } else if (kind == 0) {
          for (int row = 0; row < k; ++row) conditional_mean[row] = m[row] + slope[row] * (q - m[focal]);
          conditional_mean[focal] = q;
          const double numerator = cpp_selnorm_conditioned_log_normalizer(
            conditional_mean.data(), zero.data(), conditional_loading.data(), k,
            coordinates - 1, se.data(), omega, n_bins, z_lower, z_upper,
            effect_sign, kernel_mode, telescope_probabilities, vector_rule,
            one_group.data(), qmc, initial_points, max_points, scrambles, tolerance,
            &error, &change, quadrature);
          local_error = std::max(error, change);
          value = std::exp(std::log(se[focal]) + dnorm(q, m[focal], sd, true) + numerator - denominator);
        } else {
          upper[focal] = kind == 1 ? q : -q;
          const double first = cpp_selnorm_conditioned_log_normalizer(m.data(),
            d.data(), l.data(), k, rank, se.data(), omega, n_bins, z_lower, z_upper,
            effect_sign, kernel_mode, telescope_probabilities, vector_rule,
            one_group.data(), qmc, initial_points, max_points, scrambles, tolerance,
            &error, &change, quadrature, lower.data(), upper.data());
          local_error = std::max(error, change);
          value = std::exp(first - denominator);
          upper[focal] = infinity;
          if (kind == 2) {
            const double first_error = local_error;
            lower[focal] = q;
            const double second = cpp_selnorm_conditioned_log_normalizer(m.data(),
              d.data(), l.data(), k, rank, se.data(), omega, n_bins, z_lower, z_upper,
              effect_sign, kernel_mode, telescope_probabilities, vector_rule,
              one_group.data(), qmc, initial_points, max_points, scrambles, tolerance,
              &error, &change, quadrature, lower.data(), upper.data());
            local_error = std::max(first_error, std::max(error, change));
            value = std::exp(log_add_exp(first, second) - denominator);
            lower[focal] = -infinity;
          }
        }
        density[index] += value / dimension;
        if (!std::isfinite(local_error) || !std::isfinite(value)) {
          absolute_error[index] = infinity;
        } else {
          absolute_error[index] += value * local_error / dimension;
        }
      }
    }
  }
  double peak = 0.0, max_error = 0.0;
  for (int index = 0; index < size; ++index) {
    peak = std::max(peak, density[index]);
    max_error = std::max(max_error, absolute_error[index]);
  }
  *relative_error = std::max(*relative_error,
    peak > 0.0 ? max_error / peak : max_error == 0.0 ? 0.0 : infinity);
}
