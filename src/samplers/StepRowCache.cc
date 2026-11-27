#include "StepRowCache.h"
#include "../distributions/DSELNORMSTEPSWITCH.h"
#include "../distributions/selnorm-jags-bounds.h"
#include "../selnorm/selnorm.h"
#include "../selnorm/selnorm-tail.h"

#include <distribution/Distribution.h>

#include <JRmath.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace jags {
namespace RoBMA {
namespace {

const double negative_infinity = -std::numeric_limits<double>::infinity();
const double infinity = std::numeric_limits<double>::infinity();

// Parameter positions of dselnorm_step_switch.
enum Parameter { p_mean = 0, p_sd = 1, p_se = 2, p_power = 3, p_omega = 4, p_lower = 5, p_upper = 6,
                 p_bin = 7, p_sign = 8, p_mode = 9, p_telescope = 10, p_count = 11 };

std::mutex test_mutex;
StepRowCacheTestFlags test_flags;

// Q(x) = P(Z > x) by the package's fast upper tail, infinite arguments included.
double upper_tail(double x)
{
  if (x == infinity) return 0.0;
  if (x == -infinity) return 1.0;
  return selnorm_tail::scalar_kernel::upper_tail(x);
}

}

bool is_step_switch_row(StochasticNode const *child)
{
  return dynamic_cast<DSELNORMSTEPSWITCH const *>(child->distribution()) != nullptr &&
    child->parents().size() >= p_count && child->length() == 1;
}

bool step_row_cache_mismatch(long double log_ratio, double exact_old, double exact_new)
{
  const double scale = std::max(1.0, std::max(std::fabs(exact_old), std::fabs(exact_new)));
  return !(std::fabs(log_ratio) <= step_row_cache_tolerance * scale);
}

StepRowCache::StepRowCache(std::vector<StochasticNode *> const &rows, unsigned int chain)
  : _chain(chain)
{
  const StepRowCacheTestFlags flags = step_row_cache_test_flags();
  _mutant = flags.mutant;
  _perturbation = flags.perturbation;
  std::size_t offset = 0;
  for (StochasticNode const *node : rows) {
    if (!is_step_switch_row(node)) throw std::logic_error("A cached row is not a step-selection row.");
    std::vector<Node const *> const &parents = node->parents();
    Row row;
    row.node = node;
    row.x = node->value(chain);
    row.mean = parents[p_mean]->value(chain);
    row.sd = parents[p_sd]->value(chain);
    row.omega = parents[p_omega]->value(chain);
    row.mode = parents[p_mode]->value(chain);
    row.bins = static_cast<int>(parents[p_omega]->length());
    row.offset = offset;
    // The parameters other than the mean, the SD, the weights and the kernel
    // mode are fixed, so their checks are made once: the mode-independent ones
    // of DSELNORMSTEPSWITCH::checkParameterValue (a finite selection SE too)
    // and, for step mode, the observed bin and a descending step partition.
    // A row that fails them, or has a parameter that is not fixed, is
    // evaluated by its node.
    bool fixed = row.bins >= 1 && parents[p_lower]->length() == parents[p_omega]->length() &&
      parents[p_upper]->length() == parents[p_omega]->length();
    for (unsigned int i = 0; i < p_count; ++i) {
      if (i == p_mean || i == p_sd || i == p_omega || i == p_mode) continue;
      fixed = fixed && parents[i]->isFixed();
    }
    if (fixed) {
      const double se = parents[p_se]->value(chain)[0], power = parents[p_power]->value(chain)[0];
      const double direction = parents[p_sign]->value(chain)[0], telescope = parents[p_telescope]->value(chain)[0];
      row.cacheable = se > 0 && std::isfinite(se) && power > 0 && std::isfinite(power) &&
        (direction == 1 || direction == -1) && (telescope == 0 || telescope == 1);
      row.se = se;
      row.power = power;
      row.effect_sign = direction;
      double const *lower = parents[p_lower]->value(chain), *upper = parents[p_upper]->value(chain);
      const double bin = parents[p_bin]->value(chain)[0];
      row.step_valid = selnorm_jags_integer_in_range(bin, 0, row.bins);
      for (int b = 0; b < row.bins; ++b) {
        _bounds.push_back(selnorm_jags_bound_to_ieee(lower[b]));
      }
      for (int b = 0; b < row.bins; ++b) {
        _bounds.push_back(selnorm_jags_bound_to_ieee(upper[b]));
      }
      double const *decoded = _bounds.data() + 2 * offset;
      for (int b = 0; b < row.bins; ++b) row.step_valid = row.step_valid && decoded[b] < decoded[row.bins + b];
      row.step_valid = row.step_valid && selnorm_is_descending_step_partition(decoded, decoded + row.bins, row.bins);
      row.observed = row.step_valid ? static_cast<int>(bin) : 0;
    } else {
      _bounds.insert(_bounds.end(), 2 * static_cast<std::size_t>(std::max(row.bins, 0)), 0.0);
    }
    offset += static_cast<std::size_t>(std::max(row.bins, 0));
    _rows.push_back(row);
  }
  _p.assign(offset, 0.0);
  _rounding.assign(offset, 0.0);
}

// The Gaussian terms of normal and step mode (DSELNORMSTEPSWITCH::logDensity
// and cpp_selnorm_kernel_lpdf) at the row's current observation, mean and SD,
// which become its cached state; the bin probabilities are left to
// fill_bins().
void StepRowCache::fill_gaussian(Row &row)
{
  const double x = *row.x, mean = *row.mean, sd = *row.sd;
  row.filled = true;
  row.binned = false;
  row.c_x = x;
  row.c_mean = mean;
  row.c_sd = sd;
  row.gauss_direct = true;
  if (!std::isfinite(mean) || !(sd > 0) || !std::isfinite(sd)) return;
  const double log_lik = dnorm(x, mean, sd, true);
  row.normal_term = row.power == 1 ? log_lik : row.power * log_lik;
  row.gauss_step = cpp_selnorm_normal_lpdf(row.effect_sign * x, row.effect_sign * mean, sd);
  row.gauss_direct = std::isnan(row.normal_term) || !std::isfinite(row.gauss_step);
}

// The bin probabilities in linear space at the cached mean and SD, and a bound
// of the rounding of each: the tails that the subtraction cancels, and the
// rounding of each standardized bound s (about eps (3 |s| + |shift / sd|),
// which moves the probability by the normal density there, at most
// Q(|s|) (1 + |s|)), at four units of rounding for the fast tail, and the
// smallest normal number for tails that underflow. Adjacent bins share a
// boundary, whose tail is evaluated once.
void StepRowCache::fill_bins(Row &row)
{
  row.binned = true;
  row.bins_direct = true;
  const double mean = row.c_mean, sd = row.c_sd;
  if (row.gauss_direct) return;
  const int bins = row.bins;
  double const *lower = _bounds.data() + 2 * row.offset, *upper = lower + bins;
  double *p = _p.data() + row.offset;
  double *rounding = _rounding.data() + row.offset;
  const double shift = row.effect_sign * mean;
  const double relative_shift = std::fabs(shift / sd);
  const double unit = 4 * std::numeric_limits<double>::epsilon();
  double last_score = std::numeric_limits<double>::quiet_NaN(), last_tail = 0;
  auto tail = [&](double score) {
    if (!(score == last_score)) {
      last_score = score;
      last_tail = upper_tail(std::fabs(score));
    }
    return last_tail;
  };
  // the rounding of the probability through a finite standardized bound
  auto bound_rounding = [&](double score, double score_tail) {
    if (!std::isfinite(score) || score_tail == 0) return 0.0;
    const double magnitude = std::fabs(score);
    return score_tail * (1 + magnitude) * (3 * magnitude + relative_shift);
  };
  double largest = 0;
  bool finite = true;
  for (int b = 0; b < bins; ++b) {
    const double up = upper[b] == infinity ? infinity : (upper[b] * row.se - shift) / sd;
    const double up_tail = tail(up);
    const double low = lower[b] == negative_infinity ? negative_infinity : (lower[b] * row.se - shift) / sd;
    const double low_tail = tail(low);
    // without cancellation: in the upper tail Q(low) - Q(up), in the lower
    // tail Phi(up) - Phi(low), across zero 1 - Q(up) - Phi(low)
    double probability, cancelled;
    if (low >= 0) {
      probability = low_tail - up_tail;
      cancelled = low_tail + up_tail;
    } else if (up <= 0) {
      probability = up_tail - low_tail;
      cancelled = up_tail + low_tail;
    } else {
      probability = 1.0 - up_tail - low_tail;
      cancelled = 1.0 + up_tail + low_tail;
    }
    if (_perturbation != 0) probability *= 1.0 + _perturbation * (b + 1) / bins;
    p[b] = probability;
    rounding[b] = unit * (cancelled + bound_rounding(up, up_tail) + bound_rounding(low, low_tail)) +
      std::numeric_limits<double>::min();
    finite = finite && std::isfinite(probability) && std::isfinite(rounding[b]);
    largest = std::max(largest, probability);
  }
  row.bins_direct = !finite || !(largest > step_row_cache_min_probability);
}

void StepRowCache::refresh()
{
  step_row_cache_counters().refreshes.fetch_add(1, std::memory_order_relaxed);
  for (Row &row : _rows) {
    if (!row.cacheable) continue;
    // the stale-cache mutant always reads the update's start
    if (_mutant == 3 || !row.filled || !(*row.x == row.c_x && *row.mean == row.c_mean && *row.sd == row.c_sd)) {
      fill_gaussian(row);
    }
  }
}

// One row's term; `direct` is set when the node evaluated it. checked_omega
// and omega_valid carry the validity of the last weights vector checked in
// this evaluation (rows usually share one).
double StepRowCache::row_value(Row &row, double const *&checked_omega, bool &omega_valid, bool &direct)
{
  direct = true;
  if (!row.cacheable) return row.node->logDensity(_chain, PDF_LIKELIHOOD);
  // the cache of this row's current state (the stale-cache mutant keeps the
  // state of the refresh)
  if (!row.filled ||
      (_mutant != 3 && !(*row.x == row.c_x && *row.mean == row.c_mean && *row.sd == row.c_sd))) {
    fill_gaussian(row);
  }
  if (row.gauss_direct) return row.node->logDensity(_chain, PDF_LIKELIHOOD);
  direct = false;
  const double mode = *row.mode;
  if (mode == SELKERNEL_NORMAL) return row.normal_term;
  if (mode != SELKERNEL_STEP) return negative_infinity;
  if (row.omega != checked_omega) {
    checked_omega = row.omega;
    omega_valid = true;
    for (int b = 0; b < row.bins && omega_valid; ++b) {
      omega_valid = std::isfinite(row.omega[b]) && row.omega[b] >= 0;
    }
  }
  if (!row.step_valid || !omega_valid) return negative_infinity;
  if (row.observed == 0) return row.normal_term;
  const double observed = row.omega[row.observed - 1];
  if (!(observed > 0)) return negative_infinity;
  if (!row.binned) fill_bins(row);
  if (row.bins_direct) {
    direct = true;
    return row.node->logDensity(_chain, PDF_LIKELIHOOD);
  }
  double const *p = _p.data() + row.offset;
  double const *rounding = _rounding.data() + row.offset;
  double normalizer = 0, normalizer_rounding = 0;
  for (int b = 0; b < row.bins; ++b) {
    normalizer += row.omega[b] * p[b];
    normalizer_rounding += row.omega[b] * rounding[b];
  }
  // a normalizer that the linear sum cannot carry or resolve goes to the node
  if (!(normalizer > step_row_cache_min_probability) || !std::isfinite(normalizer) ||
      !(normalizer_rounding <= step_row_cache_precision * normalizer)) {
    direct = true;
    return row.node->logDensity(_chain, PDF_LIKELIHOOD);
  }
  // the order of cpp_selnorm_kernel_lpdf
  double out = row.gauss_step;
  if (observed != 1 && _mutant != 2) out += std::log(observed);
  out -= std::log(normalizer);
  return row.power == 1 || _mutant == 1 ? out : row.power * out;
}

double StepRowCache::value(StochasticNode const **failed)
{
  double out = 0;
  double const *checked_omega = nullptr;
  bool omega_valid = false;
  std::uint64_t direct_rows = 0;
  for (Row &row : _rows) {
    bool direct;
    const double term = row_value(row, checked_omega, omega_valid, direct);
    direct_rows += direct;
    if (std::isnan(term)) {
      if (failed != nullptr) *failed = row.node;
      return term;
    }
    out += term;
  }
  StepRowCacheCounters &counters = step_row_cache_counters();
  counters.direct.fetch_add(direct_rows, std::memory_order_relaxed);
  counters.cached.fetch_add(_rows.size() - direct_rows, std::memory_order_relaxed);
  return out;
}

void step_row_cache_test_set(StepRowCacheTestFlags const &flags)
{
  const std::lock_guard<std::mutex> lock(test_mutex);
  test_flags = flags;
}

StepRowCacheTestFlags step_row_cache_test_flags()
{
  const std::lock_guard<std::mutex> lock(test_mutex);
  return test_flags;
}

StepRowCacheCounters &step_row_cache_counters()
{
  // Never destroyed: samplers may count late.
  static StepRowCacheCounters *counters = new StepRowCacheCounters();
  return *counters;
}

}
}
