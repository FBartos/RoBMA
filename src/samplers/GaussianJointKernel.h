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

inline bool active(unsigned int cell, unsigned int coefficient)
{
  return (cell >> coefficient) & 1U;
}

// Scale each coefficient by its prior SD. QR starts with the identity prior
// precision and the scaled prior means on its right-hand side, then incorporates
// whitened observation rows by Givens rotations. The residuals give the
// marginal quadratic directly. Keeping the means on the prior rows avoids
// subtracting two prior-scale values to restore a data-scale coefficient.
inline Conditional evaluate(Model const &model, unsigned int cell, double tau)
{
  constexpr long double half_log_two_pi = 0.9189385332046727417803297364056176L;
  const long double intercept_scale = active(cell, 0) ? model.effect_sign * model.sd[0] : 0;
  const long double slope_scale = active(cell, 1) ? model.effect_sign * model.sd[1] : 0;
  Conditional out;
  out.intercept_rhs = static_cast<long double>(model.mean[0]) / model.sd[0];
  out.slope_rhs = static_cast<long double>(model.mean[1]) / model.sd[1];
  for (Observation const &row : model.rows) {
    const long double sd = std::hypot(static_cast<long double>(row.se),
                                      static_cast<long double>(tau));
    const long double a = intercept_scale / sd;
    long double b = slope_scale * row.moderator / sd;
    long double rhs = row.y / sd;
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
    out.log_likelihood -= half_log_two_pi + std::log(sd) + 0.5L * rhs * rhs;
  }
  out.log_likelihood -= std::log(out.intercept_precision_root)
                       + std::log(out.slope_precision_root);
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
      cumulative += proposal_mass[j] / total;
      _cdf[i] = cumulative;
      _log_density[j] = std::log(proposal_mass[j] / total) + log_frequency(j % scale_bins);
    }
    _cdf.back() = 1;
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

private:
  std::array<double, slab_bins - 1> _cuts;
  std::array<unsigned int, cell_count> _order;
  std::array<long double, cell_count> _cdf;
  std::array<long double, cell_count> _log_density;

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
