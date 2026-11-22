#ifndef ROBMA_GAUSSIAN_JOINT_KERNEL_H_
#define ROBMA_GAUSSIAN_JOINT_KERNEL_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace robma_gaussian_joint {

struct Observation {
  double y;
  double se;
  double moderator;
};

struct Model {
  std::vector<Observation> rows;
  double effect_sign;
  std::array<double, 2> mean;
  std::array<double, 3> sd;
  // Rows are intercept, moderator and heterogeneity. Within a row the
  // probabilities are in physical order (spike, slab), not JAGS label order.
  std::array<std::array<long double, 2>, 3> log_prior;
};

struct Conditional {
  long double log_likelihood = 0;
  long double intercept_precision_root = 1;
  long double slope_precision_root = 1;
  long double cross = 0;
  long double intercept_rhs = 0;
  long double slope_rhs = 0;
};

struct OrdinateBounds {
  long double normalizer_lower;
  long double normalizer_upper;
  long double half_quadratic_lower;
  long double half_quadratic_upper;
  bool valid;
};

struct IntervalBounds {
  long double lower;
  long double upper;
};

struct OrdinateWork {
  long double normalizer = 0;
  long double sum_abs_normalizer = 0;
  long double twice_quadratic = 0;
  long double design_squared_norm = 2;
  long double rhs_squared_norm = 0;
};

// The augmented weighted design contains either or both identity prior rows.
// Its least singular value is at least one. Guard the normalizing determinant
// and residual norm against a normwise QR perturbation before applying bounds.
inline OrdinateBounds enclose_ordinate(OrdinateWork const &work, std::size_t rows)
{
  // 64 rounded operations per row covers the two Givens rotations and the
  // norm/normalizer accumulation. Allow for double-precision libm rounding as
  // well, including builds where long double aliases double.
  const long double roundoff =
    64 * (static_cast<long double>(rows) + 2) * std::numeric_limits<long double>::epsilon() +
    32 * std::numeric_limits<double>::epsilon();
  const long double rho = roundoff * std::sqrt(work.design_squared_norm);
  if (!(rho < 1) || !std::isfinite(work.normalizer) ||
      !std::isfinite(work.sum_abs_normalizer) || !std::isfinite(work.twice_quadratic) ||
      !std::isfinite(work.rhs_squared_norm)) {
    return {0, 0, 0, 0, false};
  }
  const long double res = std::sqrt(work.twice_quadratic);
  const long double rhs = std::sqrt(work.rhs_squared_norm);
  // Singular-value perturbation bounds the determinant. Perturbing the
  // column-space projector and right-hand side bounds the residual norm.
  const long double normalizer_error =
    roundoff * (1 + work.sum_abs_normalizer) - 2 * std::log1p(-rho);
  const long double residual_error =
    (roundoff * rhs + 2 * rho * rhs) / (1 - rho) + roundoff * (1 + res);
  const long double residual_lower = std::max(0.0L, res - residual_error);
  const long double residual_upper = res + residual_error;
  const long double half_q_lower = 0.5L * residual_lower * residual_lower;
  const long double half_q_upper = 0.5L * residual_upper * residual_upper;
  if (!std::isfinite(normalizer_error) || !std::isfinite(half_q_upper)) {
    return {0, 0, 0, 0, false};
  }
  const long double inf = std::numeric_limits<long double>::infinity();
  const OrdinateBounds out{
    std::nextafter(work.normalizer - normalizer_error, -inf),
    std::nextafter(work.normalizer + normalizer_error, inf),
    std::max(0.0L, std::nextafter(half_q_lower, -inf)),
    std::nextafter(half_q_upper, inf),
    true
  };
  if (!std::isfinite(out.normalizer_lower) || !std::isfinite(out.normalizer_upper) ||
      !std::isfinite(out.half_quadratic_upper)) return {0, 0, 0, 0, false};
  return out;
}

inline bool active(unsigned int cell, unsigned int coefficient)
{
  return (cell >> coefficient) & 1U;
}

inline long double upper_round(long double value)
{
  return std::nextafter(value, std::numeric_limits<long double>::infinity());
}

inline IntervalBounds enclose_interval(OrdinateBounds const &left, OrdinateBounds const &right,
                                       long double left_tau_lower, long double right_tau_upper,
                                       long double min_variance_lower, long double mean_diagonal_upper,
                                       long double n)
{
  // Both the log prefactor and inverse quadratic decrease with tau.
  // At infinity the right prefactor is -Inf and the right quadratic is zero.
  IntervalBounds result{right.normalizer_lower - left.half_quadratic_upper,
                        left.normalizer_upper - right.half_quadratic_lower};
  const long double inf = std::numeric_limits<long double>::infinity();
  const long double left_tau_sq = left_tau_lower * left_tau_lower;
  if (!std::isfinite(left_tau_sq) || !(min_variance_lower > 0)) return result;
  const long double variance_sum = min_variance_lower + std::nextafter(left_tau_sq, 0.0L);
  if (!std::isfinite(variance_sum)) return result;
  const long double dispersion_sq = std::nextafter(variance_sum, 0.0L);
  if (!(dispersion_sq > 0)) return result;
  const long double quadratic = upper_round(2 * left.half_quadratic_upper);
  const long double left_upper = upper_round(left.normalizer_upper - left.half_quadratic_lower);

  if (!std::isfinite(right_tau_upper)) {
    // Jensen gives tr(C^-1) >= n / mean(diag(C)). The ratio of the mean
    // eigenvalue to its lower bound is decreasing in u=tau². If the
    // quadratic times this ratio is below n, log L decreases on the tail.
    const long double ratio = upper_round(upper_round(mean_diagonal_upper +
      upper_round(left_tau_sq)) / dispersion_sq);
    if (std::isfinite(left_upper) && upper_round(quadratic * ratio) < n) {
      result.upper = std::min(result.upper, left_upper);
    }
    return result;
  }

  // In u=tau², (log L)'' lies between -Q(a)/(min_variance+a)^2 and
  // n/[2(min_variance+a)^2]. Those curvature limits enclose the difference
  // from the endpoint secant. Enlarge b-a using the padded physical edges.
  const long double tau_difference = upper_round(right_tau_upper - left_tau_lower);
  const long double tau_sum = upper_round(right_tau_upper + left_tau_lower);
  const long double width = upper_round(upper_round(tau_difference * tau_sum) / dispersion_sq);
  const long double width_sq = upper_round(width * width);
  const long double lower_error = upper_round(upper_round(n / 16) * width_sq);
  const long double upper_error = upper_round(upper_round(quadratic / 8) * width_sq);
  const long double endpoint_lower = std::min(
    std::nextafter(left.normalizer_lower - left.half_quadratic_upper, -inf),
    std::nextafter(right.normalizer_lower - right.half_quadratic_upper, -inf));
  const long double endpoint_upper = std::max(left_upper,
    upper_round(right.normalizer_upper - right.half_quadratic_lower));
  const long double lower = std::nextafter(endpoint_lower - lower_error, -inf);
  const long double upper = upper_round(endpoint_upper + upper_error);
  if (std::isfinite(lower)) result.lower = std::max(result.lower, lower);
  if (std::isfinite(upper)) result.upper = std::min(result.upper, upper);
  return result;
}

// Scale each coefficient by its prior SD. QR starts with the identity prior
// precision and the scaled prior means on its right-hand side, then incorporates
// whitened observation rows by Givens rotations. The residuals give the
// marginal quadratic directly. Keeping the means on the prior rows avoids
// subtracting two prior-scale values to restore a data-scale coefficient.
template <bool enclose = false>
inline Conditional evaluate(Model const &model, unsigned int cell, double tau,
                            OrdinateBounds *bounds = nullptr)
{
  constexpr long double half_log_two_pi = 0.9189385332046727417803297364056176L;
  const long double intercept_scale = active(cell, 0) ? model.effect_sign * model.sd[0] : 0;
  const long double slope_scale = active(cell, 1) ? model.effect_sign * model.sd[1] : 0;
  Conditional out;
  OrdinateWork work;
  out.intercept_rhs = static_cast<long double>(model.mean[0]) / model.sd[0];
  out.slope_rhs = static_cast<long double>(model.mean[1]) / model.sd[1];
  if constexpr (enclose) {
    if (intercept_scale != 0) work.rhs_squared_norm += out.intercept_rhs * out.intercept_rhs;
    if (slope_scale != 0) work.rhs_squared_norm += out.slope_rhs * out.slope_rhs;
  }
  for (Observation const &row : model.rows) {
    const long double sd = std::hypot(static_cast<long double>(row.se),
                                      static_cast<long double>(tau));
    const long double a = intercept_scale / sd;
    long double b = slope_scale * row.moderator / sd;
    long double rhs = row.y / sd;
    if constexpr (enclose) {
      work.design_squared_norm += a * a + b * b;
      work.rhs_squared_norm += rhs * rhs;
    }
    if (a != 0) {
      const long double r = std::hypot(out.intercept_precision_root, a);
      const long double c = out.intercept_precision_root / r;
      const long double s = a / r;
      const long double cross = out.cross;
      const long double intercept_rhs = out.intercept_rhs;
      out.intercept_precision_root = r;
      out.cross = c * cross + s * b;
      b = c * b - s * cross;
      out.intercept_rhs = c * intercept_rhs + s * rhs;
      rhs = c * rhs - s * intercept_rhs;
    }
    if (b != 0) {
      const long double r = std::hypot(out.slope_precision_root, b);
      const long double c = out.slope_precision_root / r;
      const long double s = b / r;
      const long double slope_rhs = out.slope_rhs;
      out.slope_precision_root = r;
      out.slope_rhs = c * slope_rhs + s * rhs;
      rhs = c * rhs - s * slope_rhs;
    }
    const long double normalizer = half_log_two_pi + std::log(sd);
    out.log_likelihood -= normalizer + 0.5L * rhs * rhs;
    if constexpr (enclose) {
      work.normalizer -= normalizer;
      work.sum_abs_normalizer += std::fabs(normalizer);
      work.twice_quadratic += rhs * rhs;
    }
  }
  const long double log_precision = std::log(out.intercept_precision_root) +
                                    std::log(out.slope_precision_root);
  out.log_likelihood -= log_precision;
  if constexpr (enclose) {
    work.normalizer -= log_precision;
    work.sum_abs_normalizer += std::fabs(log_precision);
    *bounds = enclose_ordinate(work, model.rows.size());
  }
  return out;
}

inline std::array<double, 2> restore_coefficients(Model const &model,
                                                   Conditional const &at,
                                                   double intercept_normal,
                                                   double slope_normal)
{
  const long double slope = (at.slope_rhs + slope_normal) / at.slope_precision_root;
  const long double intercept = (at.intercept_rhs + intercept_normal - at.cross * slope)
                                / at.intercept_precision_root;
  return {
    static_cast<double>(model.sd[0] * intercept),
    static_cast<double>(model.sd[1] * slope)
  };
}

class Proposal {
public:
  static constexpr unsigned int slab_bins = 64;
  static constexpr unsigned int scale_bins = slab_bins + 1;
  static constexpr unsigned int coefficient_cells = 4;
  static constexpr unsigned int cell_count = scale_bins * coefficient_cells;
  struct Coverage {
    long double log_bound;
    unsigned int subintervals_per_bin;
    bool certified;
  };

  struct State {
    unsigned int coefficients;
    unsigned int bin;
    double tau;
    Conditional conditional;
    long double log_ratio;
  };

  // 'quantile' takes a half-normal upper-tail probability in (0, 1).
  // This avoids subtraction from one when a draw falls in the tail.
  template <typename Quantile>
  bool prepare(Model const &model, Quantile quantile)
  {
    for (unsigned int i = 0; i < _cuts.size(); ++i) {
      _cuts[i] = quantile(static_cast<double>(slab_bins - i - 1) / slab_bins);
      if (!std::isfinite(_cuts[i]) || !(_cuts[i] > (i == 0 ? 0 : _cuts[i - 1]))) {
        return false;
      }
    }

    constexpr long double prior_share = 0.01L;
    std::array<long double, cell_count> log_mass;
    std::array<long double, cell_count> prior_mass;
    long double largest = -std::numeric_limits<long double>::infinity();
    for (unsigned int bin = 0; bin < scale_bins; ++bin) {
      const double tau = bin == 0 ? 0 : quantile((slab_bins - bin + 0.5) / slab_bins);
      if (!std::isfinite(tau) || (bin != 0 && bucket(tau) != bin)) return false;
      for (unsigned int cell = 0; cell < coefficient_cells; ++cell) {
        const unsigned int i = index(cell, bin);
        const long double log_prior_mass = log_prior(model, cell, bin)
                                            - log_frequency(bin);
        prior_mass[i] = std::exp(log_prior_mass);
        log_mass[i] = evaluate(model, cell, tau).log_likelihood + log_prior_mass;
        if (!std::isfinite(log_mass[i]) || !(prior_mass[i] > 0)) return false;
        largest = std::max(largest, log_mass[i]);
      }
    }

    long double likelihood_total = 0;
    for (long double &mass : log_mass) {
      mass = std::exp(mass - largest);
      likelihood_total += mass;
    }

    std::array<long double, cell_count> proposal_mass;
    long double total = 0;
    for (unsigned int i = 0; i < cell_count; ++i) {
      const long double mass = (1 - prior_share) * log_mass[i] / likelihood_total
                               + prior_share * prior_mass[i];
      proposal_mass[i] = mass;
      total += mass;
    }
    // Building the cumulative probabilities in ascending mass order retains
    // small prior components in floating-point arithmetic.
    for (unsigned int i = 0; i < cell_count; ++i) _order[i] = i;
    std::stable_sort(_order.begin(), _order.end(),
      [&](unsigned int i, unsigned int j) { return proposal_mass[i] < proposal_mass[j]; });
    long double cumulative = 0;
    for (unsigned int i = 0; i < cell_count; ++i) {
      const unsigned int j = _order[i];
      const long double next = i == cell_count - 1 ? 1 : cumulative + proposal_mass[j] / total;
      if (!(next > cumulative) || next > 1) return false;
      _cdf[i] = next;
      // Score and verify the actual CDF interval, including rounding and
      // the last cell's remaining probability.
      _log_density[j] = std::log(next - cumulative) + log_frequency(j % scale_bins);
      cumulative = next;
    }
    return true;
  }

  unsigned int bucket(double tau) const
  {
    return 1 + static_cast<unsigned int>(
      std::upper_bound(_cuts.begin(), _cuts.end(), tau) - _cuts.begin());
  }

  State state(Model const &model, unsigned int coefficients, unsigned int bin, double tau) const
  {
    const Conditional at = evaluate(model, coefficients, tau);
    return {coefficients, bin, tau, at,
      at.log_likelihood + log_prior(model, coefficients, bin)
                        - _log_density[index(coefficients, bin)]};
  }

  template <typename Quantile, typename RNG>
  State draw(Model const &model, Quantile quantile, RNG &rng) const
  {
    const auto selected = std::upper_bound(_cdf.begin(), _cdf.end(), rng.uniform());
    const unsigned int i = _order[selected - _cdf.begin()];
    const unsigned int bin = i % scale_bins;
    const double tau = bin == 0 ? 0 : quantile((slab_bins - bin + rng.uniform()) / slab_bins);
    return state(model, i / scale_bins, bin, tau);
  }

  long double log_density(unsigned int coefficients, unsigned int bin) const
  {
    return _log_density[index(coefficients, bin)];
  }

  static long double coverage_limit()
  {
    // P(no posterior refresh after 100 sweeps) <= .001 from every start.
    return -std::log(-std::expm1(std::log(0.001L) / 100));
  }

  template <typename Quantile>
  Coverage coverage(Model const &model, Quantile quantile) const
  {
    const long double limit = coverage_limit();
    const long double coarse = coverage_bound(model, quantile, 1);
    if (!std::isfinite(coarse)) return {coarse, 1, false};
    if (coarse <= limit) return {coarse, 1, true};
    const long double refined = coverage_bound(model, quantile, 4);
    return {refined, 4, std::isfinite(refined) && refined <= limit};
  }

  // Exposed to the standalone numerical verification for the two stages.
  template <typename Quantile>
  long double coverage_bound(Model const &model, Quantile quantile,
                             unsigned int subintervals_per_bin) const
  {
    const long double inf = std::numeric_limits<long double>::infinity();
    if (subintervals_per_bin != 1 && subintervals_per_bin != 4) return inf;
    const auto n_intervals = slab_bins * subintervals_per_bin;
    const long double log_width = -std::log(static_cast<long double>(n_intervals));
    const auto smallest_se = std::min_element(model.rows.begin(), model.rows.end(),
      [](Observation const &a, Observation const &b) { return a.se < b.se; });
    if (smallest_se == model.rows.end() || !(smallest_se->se > 0)) return inf;
    const long double n = model.rows.size();
    const long double min_se = smallest_se->se;
    const long double min_variance = min_se * min_se;
    const long double min_variance_lower =
      std::isfinite(min_variance) ? std::nextafter(min_variance, 0.0L) : 0;
    long double sampling_trace = 0;
    long double slope_trace = 0;
    for (const Observation &row : model.rows) {
      const long double se = row.se;
      const long double slope_sd = upper_round(std::fabs(row.moderator) *
                                               static_cast<long double>(model.sd[1]));
      sampling_trace = upper_round(sampling_trace + upper_round(se * se));
      slope_trace = upper_round(slope_trace + upper_round(slope_sd * slope_sd));
    }
    const long double intercept_sd = model.sd[0];
    const long double sampling_mean = upper_round(sampling_trace / n);
    const long double intercept_mean = upper_round(intercept_sd * intercept_sd);
    const long double slope_mean = upper_round(slope_trace / n);
    std::vector<RankEdge> edges(n_intervals - 1);
    for (unsigned int j = 1; j < n_intervals; ++j) {
      const double survival = static_cast<double>(n_intervals - j) / n_intervals;
      if (!rank_edge(model.sd[2], smallest_se->se, quantile, survival, edges[j - 1])) return inf;
    }
    long double log_evidence_lower = -inf;
    long double log_weight_upper = -inf;
    for (unsigned int cell = 0; cell < coefficient_cells; ++cell) {
      OrdinateBounds atom;
      evaluate<true>(model, cell, 0, &atom);
      if (!atom.valid) return inf;
      const long double log_atom = log_prior(model, cell, 0);
      log_evidence_lower = add_log_mass(log_evidence_lower,
        log_atom + atom.normalizer_lower - atom.half_quadratic_upper);
      log_weight_upper = std::max(log_weight_upper,
        log_atom + atom.normalizer_upper - atom.half_quadratic_lower - log_density(cell, 0));
      const long double log_slab = log_prior(model, cell, 1);
      OrdinateBounds left = atom;
      long double left_tau_lower = 0;
      const long double mean_diagonal_upper = upper_round(sampling_mean +
        upper_round((active(cell, 0) ? intercept_mean : 0) + (active(cell, 1) ? slope_mean : 0)));
      for (unsigned int j = 1; j <= n_intervals; ++j) {
        OrdinateBounds right = {-inf, -inf, 0, 0, true};
        long double right_tau_upper = inf;
        if (j < n_intervals) {
          right = rank_edge_bounds(model, cell, edges[j - 1]);
          if (!right.valid) return inf;
          right_tau_upper = edges[j - 1].upper_tau;
        }
        const unsigned int bin = 1 + (j - 1) / subintervals_per_bin;
        const auto interval = enclose_interval(left, right, left_tau_lower, right_tau_upper,
                                               min_variance_lower, mean_diagonal_upper, n);
        log_evidence_lower = add_log_mass(log_evidence_lower, log_slab + interval.lower + log_width);
        log_weight_upper = std::max(log_weight_upper, log_slab + interval.upper - log_density(cell, bin));
        left = right;
        if (j < n_intervals) left_tau_lower = edges[j - 1].lower_tau;
      }
    }
    if (!std::isfinite(log_evidence_lower) || !std::isfinite(log_weight_upper)) return inf;
    // Round the assembled log ratio outwards as well. Include cancellation
    // between the ordinate and the evidence when setting the margin.
    const long double error = 64 * n_intervals * std::numeric_limits<long double>::epsilon() *
      (1 + std::fabs(log_evidence_lower) + std::fabs(log_weight_upper));
    const long double bound = std::nextafter(log_weight_upper - log_evidence_lower + error, inf);
    return bound >= 0 ? bound : inf;
  }

private:
  std::array<double, slab_bins - 1> _cuts;
  std::array<unsigned int, cell_count> _order;
  std::array<long double, cell_count> _cdf;
  std::array<long double, cell_count> _log_density;

  struct RankEdge {
    double tau;
    double lower_tau;
    double upper_tau;
    long double left_dispersion_ratio;
    long double right_dispersion_ratio;
  };

  static long double add_log_mass(long double a, long double b)
  {
    const long double lower = std::min(a, b);
    const long double upper = std::max(a, b);
    if (!std::isfinite(lower)) return upper;
    return upper + std::log1p(std::exp(lower - upper));
  }

  template <typename Quantile>
  static bool rank_edge(double tau_sd, double min_se, Quantile quantile,
                        double survival, RankEdge &edge)
  {
    // Evaluate outward of the nominal rank edge, and independently check
    // that the physical interval encloses it with a numerical margin.
    constexpr double pad = 64 * std::numeric_limits<double>::epsilon();
    const double lo = quantile(survival + pad);
    const double tau = quantile(survival);
    const double hi = quantile(survival - pad);
    if (!(lo > 0 && lo < tau && tau < hi) || !std::isfinite(hi)) return false;
    const long double scale = tau_sd * std::sqrt(2.0L);
    if (!(std::erfc(lo / scale) > survival + pad / 2 &&
          std::erfc(hi / scale) < survival - pad / 2)) return false;
    const auto dispersion = [=](double at) {
      return std::hypot(static_cast<long double>(min_se), static_cast<long double>(at));
    };
    constexpr long double margin = 1 + 8 * std::numeric_limits<long double>::epsilon();
    const long double dispersion_at_tau = dispersion(tau);
    edge = {tau, lo, hi, margin * dispersion_at_tau / dispersion(lo),
                 margin * dispersion(hi) / dispersion_at_tau};
    return true;
  }

  static OrdinateBounds rank_edge_bounds(Model const &model, unsigned int cell, RankEdge const &edge)
  {
    OrdinateBounds out;
    evaluate<true>(model, cell, edge.tau, &out);
    if (!out.valid) return out;
    // C(t) = C(0) + t²I. Its eigenvalues have relative changes no larger
    // than the variance changes at the smallest sampling SE. Enclose the
    // log determinant and inverse quadratic across the padded edge using
    // this ratio, with just one QR evaluation for the edge and cell.
    const long double inf = std::numeric_limits<long double>::infinity();
    const long double n = model.rows.size();
    const long double left = edge.left_dispersion_ratio;
    const long double right = edge.right_dispersion_ratio;
    out.normalizer_lower = std::nextafter(out.normalizer_lower - n * std::log(right), -inf);
    out.normalizer_upper = std::nextafter(out.normalizer_upper + n * std::log(left), inf);
    out.half_quadratic_lower = std::max(0.0L,
      std::nextafter(out.half_quadratic_lower / (right * right), -inf));
    out.half_quadratic_upper = std::nextafter(out.half_quadratic_upper * left * left, inf);
    out.valid = std::isfinite(out.normalizer_lower) && std::isfinite(out.normalizer_upper) &&
                std::isfinite(out.half_quadratic_upper);
    return out;
  }

  static unsigned int index(unsigned int coefficients, unsigned int bin)
  {
    return coefficients * scale_bins + bin;
  }

  static long double log_prior(Model const &model, unsigned int cell, unsigned int bin)
  {
    return model.log_prior[0][active(cell, 0)] + model.log_prior[1][active(cell, 1)]
           + model.log_prior[2][bin != 0];
  }

  static long double log_frequency(unsigned int bin)
  {
    return bin == 0 ? 0 : std::log(static_cast<long double>(slab_bins));
  }
};

} // namespace robma_gaussian_joint

#endif
