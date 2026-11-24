#ifndef ROBMA_STEP_ROW_CACHE_H
#define ROBMA_STEP_ROW_CACHE_H

// The row-cache surrogate of singleton step-selection rows
// (dselnorm_step_switch children) for RoBMA's corrected samplers of the
// publication-bias side: the gated weight block (the weights of one label
// move) and the screened indicator (the label moves).
//
// In step mode with observed bin o > 0, row i's log density is
//   pw_i [log phi(sign y_i; sign m_i, sd_i) + log omega_o - log A_i],
//   A_i = sum_b omega_b p_ib,  p_ib = P(lower_b < sign y / se_i < upper_b),
// y ~ N(m_i, sd_i^2), and pw_i (log phi(y_i; m_i, sd_i)) in normal mode or
// with o = 0 (DSELNORMSTEPSWITCH::logDensity). While only the weights or the
// label move, the rows' observations, means, SDs, selection SEs, powers, bin
// bounds, observed bins, effect directions and telescoping are constant, so
// p_ib and both Gaussian terms are cached once per update (refresh) and an
// evaluation is a sum over the rows of a short dot product and a logarithm.
// The kernel mode, the weights and their validity (finite, nonnegative, a
// positive observed weight) are read at every evaluation; telescoping only
// changes the route of the exact normalizer, not its value.
//
// Exactness does not rest on the cache: it is the surrogate of samplers that
// correct every update with the exact density through the graph. A row is
// evaluated by its own node (the exact density) instead of the cache when
//   - at the refresh, its parameters are not ordinary (a parameter check of
//     the distribution would fail at the refresh values, an infinite SD or
//     selection SE, a parameter other than the mean, the SD, the weights and
//     the kernel mode is not fixed), or no bin probability is at least
//     `step_row_cache_min_probability` (the probabilities cannot carry the
//     normalizer);
//   - at an evaluation, its mean or SD differs from the refresh values (a
//     label with PET or PEESE moves the mean), or the cached normalizer
//     sum_b omega_b p_ib is not above `step_row_cache_min_probability`, the
//     threshold below which StepRowNormalizer leaves its linear sum for the
//     log scale (the sum would underflow or lose its leading digits).
//
// Evidence: a corrected proposal whose surrogate includes cached rows counts as
// evidence when |log acceptance ratio| exceeds
// step_row_cache_tolerance * max(1, |exact log densities|): exact up to
// rounding, the cache never sends a sampler to its fallback by itself.

#include <graph/StochasticNode.h>

#include <atomic>
#include <cstdint>
#include <vector>

namespace jags {
namespace RoBMA {

// A cached normalizer (and the largest bin probability of a row) must exceed
// this value; the same threshold as StepRowNormalizer's linear path.
const double step_row_cache_min_probability = 1e-280;
// Relative tolerance of the evidence comparison (see above).
const double step_row_cache_tolerance = 1e-9;

// Whether a stochastic child is a singleton step-selection row.
bool is_step_switch_row(StochasticNode const *child);

// Whether the evidence comparison finds a difference between the surrogate and
// exact changes of a proposal: |log_ratio| above the tolerance relative to the
// largest exact log density.
bool step_row_cache_mismatch(long double log_ratio, double exact_old, double exact_new);

class StepRowCache {
public:
  StepRowCache() = default;
  // The rows of one chain (each must satisfy is_step_switch_row).
  StepRowCache(std::vector<StochasticNode *> const &rows, unsigned int chain);

  bool empty() const { return _rows.empty(); }
  std::size_t size() const { return _rows.size(); }

  // Reads every row at the chain's current values and caches its bin
  // probabilities and Gaussian terms (once per update).
  void refresh();
  // The sum of the rows' log densities (PDF_LIKELIHOOD) at the chain's
  // current values. A NaN sum names the first row with a NaN term in `failed`.
  double value(StochasticNode const **failed) const;

private:
  struct Row {
    StochasticNode const *node = nullptr;
    // the observation and the parameters that may move (stable value pointers)
    double const *x = nullptr, *mean = nullptr, *sd = nullptr, *omega = nullptr, *mode = nullptr;
    int bins = 0;
    // fixed parameters, checked once: cacheable when every parameter other
    // than the mean, SD, weights and mode is fixed and passes its check
    bool cacheable = false;
    bool step_valid = false;  // the step-mode checks that do not involve the weights
    int observed = 0;
    double se = 1, power = 1, effect_sign = 1;
    std::size_t offset = 0;   // the row's probabilities in _p, its bounds at 2 offset in _bounds
    // refresh state
    bool direct = true;       // evaluated by its node at every evaluation of this update
    double c_mean = 0, c_sd = 0;
    double normal_term = 0;   // pw log phi(y; m, sd) as in normal mode
    double gauss_step = 0;    // log phi(sign y; sign m, sd) as in step mode
  };
  std::vector<Row> _rows;
  std::vector<double> _p;
  std::vector<double> _bounds;  // per row its decoded lower, then upper bounds
  unsigned int _chain = 0;
  int _mutant = 0;
  double _perturbation = 0;

  double row_value(Row const &row, double const *&checked_omega, bool &omega_valid, bool &direct) const;
};

// Test hooks (process-wide, read when a model's samplers are built). Mutants
// of the cache (not exact surrogates; the samplers stay exact through their
// correction): 1 the power omitted, 2 the observed-bin term omitted, 3 a stale
// cache (the mean and SD comparison skipped). perturbation: each cached bin
// probability p_b of a row is multiplied by 1 + perturbation (b + 1) / bins (a
// deliberately poor surrogate for the exactness fixtures). disable: the
// factories ignore step-selection rows (the eligibility before the cache).
// Counters: refreshes, rows evaluated from the cache, and rows evaluated by
// their node (all reasons).
struct StepRowCacheTestFlags {
  int mutant = 0;
  double perturbation = 0;
  bool disable = false;
};
void step_row_cache_test_set(StepRowCacheTestFlags const &flags);
StepRowCacheTestFlags step_row_cache_test_flags();
struct StepRowCacheCounters {
  std::atomic<std::uint64_t> refreshes{0};
  std::atomic<std::uint64_t> cached{0};
  std::atomic<std::uint64_t> direct{0};
};
StepRowCacheCounters &step_row_cache_counters();

}
}
#endif
