#ifndef ROBMA_INDEPENDENT_BIAS_PROPOSAL_H
#define ROBMA_INDEPENDENT_BIAS_PROPOSAL_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "../selnorm/selnorm-tail.h"

namespace jags::RoBMA::joint_bias {

constexpr double log_sqrt_two_pi = 0.91893853320467274178;
constexpr double prior_mixture = 0.03;

enum class Family { unadjusted, selection, pet, peese };

// All locations are in the likelihood's coordinates. An exempt row has
// bin zero and a zero PET/PEESE predictor, but contributes its normal density.
struct Row {
  double outcome;
  double selection_se;
  double tilt_se;
  int bin;
  int selection_sign;
};

struct Geometry {
  std::vector<Family> families;
  std::vector<Row> rows;
  // Descending standard normal cutoffs, including +Inf and -Inf.
  std::vector<double> cuts;
  // Cells are in original bias-label order, spike cells first, then slab.
  std::vector<double> log_prior;
  double prior_mean;
  double prior_precision;
  double prior_sd;
};

inline double log_normal(double x, double mean, double sd)
{
  const double z = (x - mean) / sd;
  return -log_sqrt_two_pi - std::log(sd) - 0.5 * z * z;
}

struct Proposal {
  std::vector<double> mean;
  std::vector<double> sd;
  std::vector<double> log_mass;
  bool local = false;

  explicit Proposal(std::size_t cells) : mean(cells), sd(cells), log_mass(cells) {}

  double log_density(const Geometry &geometry, std::size_t cell, double value) const
  {
    const double prior = geometry.log_prior[cell] +
      log_normal(value, geometry.prior_mean, geometry.prior_sd);
    if (!local) return prior;
    const double a = std::log(prior_mixture) + prior;
    const double b = std::log1p(-prior_mixture) + log_mass[cell] +
      log_normal(value, mean[cell], sd[cell]);
    const double top = std::max(a, b);
    return std::isfinite(top) ? top + std::log1p(std::exp(-std::abs(a - b))) : top;
  }

  std::size_t cell(const Geometry &geometry, bool use_prior, double uniform) const
  {
    const auto &mass = use_prior ? geometry.log_prior : log_mass;
    std::size_t last = 0;
    for (std::size_t i = 0; i < mass.size(); ++i) {
      const double probability = std::exp(mass[i]);
      if (!(probability > 0)) continue;
      last = i;
      uniform -= probability;
      if (uniform < 0) break;
    }
    return last;
  }
};

struct NormalBoundary {
  double score;
  double tail;
  double density;
  double score_density;
};

inline NormalBoundary boundary(double score)
{
  if (std::isinf(score)) return {score, 0, 0, 0};
  const double density = std::exp(-0.5 * score * score - log_sqrt_two_pi);
  return {score, selnorm_tail::scalar_kernel::upper_tail(std::abs(score)),
          density, score * density};
}

inline double interval_area(const NormalBoundary &lower, const NormalBoundary &upper)
{
  if (lower.score >= 0) return lower.tail - upper.tail;
  if (upper.score <= 0) return upper.tail - lower.tail;
  return 1 - (lower.tail + upper.tail);
}

struct Bin {
  double central;
  double spike;
  double first;
  double second;
};

// One pass over the common partition supplies every branch's normalizer and
// its first two local derivatives. The moment identities also underlie the
// selected-normal kernel's shared tail and moment calculation.
inline void normal_bins(const Geometry &geometry, const Row &row, double center,
                        double sd, std::vector<Bin> &bins)
{
  NormalBoundary upper = boundary(geometry.cuts.front());
  NormalBoundary upper_spike = upper;
  for (std::size_t i = 0; i < bins.size(); ++i) {
    const double cut = geometry.cuts[i + 1] * row.selection_se;
    const auto lower = boundary((cut - row.selection_sign * center) / sd);
    const auto lower_spike = boundary(cut / sd);
    const double central = interval_area(lower, upper);
    bins[i] = {central, interval_area(lower_spike, upper_spike),
               lower.density - upper.density,
               central + lower.score_density - upper.score_density};
    upper = lower;
    upper_spike = lower_spike;
  }
}

struct Quadratic {
  double spike = 0;
  double central = 0;
  double gradient = 0;
  double curvature = 0;
  bool possible = true;
};

struct Scratch {
  std::vector<Bin> bins;
  std::vector<Quadratic> terms;

  explicit Scratch(const Geometry &geometry)
    : bins(geometry.cuts.size() - 1), terms(geometry.families.size()) {}
};

inline void add_gaussian(Quadratic &term, double residual, double center,
                         double log_height, double precision)
{
  const double error = residual - center;
  term.spike += log_height - 0.5 * residual * (residual * precision);
  term.central += log_height - 0.5 * error * (error * precision);
  term.gradient += error * precision;
  term.curvature += precision;
}

// The expansion is a proposal. A failed expansion uses precisely the original
// joint prior for the same fixed nuisance state in both directions of MH.
inline void propose(const Geometry &geometry, const std::vector<double> &total_sd,
                    const double *slopes, const double *omega, Proposal &out, Scratch &scratch)
{
  out.local = false;
  double precision = geometry.prior_precision;
  double weighted_center = precision * geometry.prior_mean;
  for (std::size_t i = 0; i < geometry.rows.size(); ++i) {
    if (!std::isfinite(total_sd[i]) || !(total_sd[i] > 0)) return;
    const double inverse_sd = 1 / total_sd[i];
    const double row_precision = inverse_sd * inverse_sd;
    precision += row_precision;
    weighted_center += row_precision * geometry.rows[i].outcome;
  }
  const double center = weighted_center / precision;
  if (!std::isfinite(center) || !std::isfinite(precision) || !(precision > 0)) return;

  auto &bins = scratch.bins;
  auto &terms = scratch.terms;
  std::fill(terms.begin(), terms.end(), Quadratic{});
  Quadratic gaussian;
  for (std::size_t i = 0; i < geometry.rows.size(); ++i) {
    const auto &row = geometry.rows[i];
    const double sd = total_sd[i];
    const double inverse_sd = 1 / sd;
    const double row_precision = inverse_sd * inverse_sd;
    const double log_height = -log_sqrt_two_pi - std::log(sd);
    add_gaussian(gaussian, row.outcome, center, log_height, row_precision);
    if (row.bin > 0) normal_bins(geometry, row, center, sd, bins);

    for (std::size_t b = 0; b < terms.size(); ++b) {
      auto &term = terms[b];
      switch (geometry.families[b]) {
        case Family::unadjusted:
          continue;
        case Family::pet:
          add_gaussian(term, row.outcome - slopes[b] * row.tilt_se, center,
                       log_height, row_precision);
          continue;
        case Family::peese:
          add_gaussian(term, row.outcome - slopes[b] * (row.tilt_se * row.tilt_se),
                       center, log_height, row_precision);
          continue;
        case Family::selection:
          break;
      }
      if (row.bin == 0 || !term.possible) continue;
      const double *weight = omega + b * bins.size();
      const double observed = weight[row.bin - 1];
      if (observed == 0) {
        term.possible = false;
        continue;
      }
      double z = 0, z_spike = 0, first = 0, second = 0;
      for (std::size_t j = 0; j < bins.size(); ++j) {
        z += weight[j] * bins[j].central;
        z_spike += weight[j] * bins[j].spike;
        first += weight[j] * bins[j].first;
        second += weight[j] * bins[j].second;
      }
      if (!(observed > 0 && z > 0 && z_spike > 0)) return;
      const double mean_score = first / z;
      term.spike += std::log(observed) - std::log(z_spike);
      term.central += std::log(observed) - std::log(z);
      term.gradient -= row.selection_sign * mean_score * inverse_sd;
      term.curvature += (second / z - mean_score * mean_score - 1) * row_precision;
    }
  }

  const double prior_center = log_normal(center, geometry.prior_mean, geometry.prior_sd);
  const double prior_gradient = (geometry.prior_mean - center) * geometry.prior_precision;
  for (std::size_t b = 0; b < terms.size(); ++b) {
    auto term = terms[b];
    if (geometry.families[b] == Family::unadjusted) term = gaussian;
    if (geometry.families[b] == Family::selection) {
      term.spike += gaussian.spike;
      term.central += gaussian.central;
      term.gradient += gaussian.gradient;
      term.curvature += gaussian.curvature;
    }
    const std::size_t active = b + terms.size();
    out.mean[b] = geometry.prior_mean;
    out.sd[b] = geometry.prior_sd;
    if (!term.possible) {
      out.mean[active] = geometry.prior_mean;
      out.sd[active] = geometry.prior_sd;
      out.log_mass[b] = -std::numeric_limits<double>::infinity();
      out.log_mass[active] = out.log_mass[b];
      continue;
    }

    const double curvature = term.curvature + geometry.prior_precision;
    const double gradient = term.gradient + prior_gradient;
    if (!std::isfinite(curvature) || !(curvature > 0)) return;
    out.mean[active] = center + gradient / curvature;
    out.sd[active] = 1 / std::sqrt(curvature);
    out.log_mass[b] = term.spike + geometry.log_prior[b];
    out.log_mass[active] = geometry.log_prior[active] + term.central + prior_center +
      0.5 * gradient * (gradient / curvature) + log_sqrt_two_pi + std::log(out.sd[active]);
    if (!std::isfinite(out.mean[active]) || !std::isfinite(out.sd[active]) ||
        !std::isfinite(out.log_mass[b]) || !std::isfinite(out.log_mass[active])) return;
  }

  const double top = *std::max_element(out.log_mass.begin(), out.log_mass.end());
  if (!std::isfinite(top)) return;
  double sum = 0;
  for (double mass : out.log_mass) sum += std::exp(mass - top);
  const double log_sum = std::log(sum);
  for (double &mass : out.log_mass) mass = (mass - top) - log_sum;
  out.local = true;
}

} // namespace jags::RoBMA::joint_bias

#endif
