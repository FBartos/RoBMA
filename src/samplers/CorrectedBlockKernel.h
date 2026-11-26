#ifndef ROBMA_CORRECTED_BLOCK_KERNEL_H
#define ROBMA_CORRECTED_BLOCK_KERNEL_H

// The corrected block kernel shared by RoBMA's block samplers
// (RoBMA::CoarseCorrectedBlock for declared blocks, RoBMA::GatedCorrectedBlock
// for the gated selection weights of one bias branch).
//
// One corrected update of the members of a GraphView: both densities at block
// entry (the members' full densities, the exact terms of the other children
// and, for the surrogate, each coarse-capable selection child's surrogate), a
// fair forward or backward order of the directions, one slice update of the
// surrogate per direction with chained values (BayesTools' header
// factor_slice_kernel.h, LinkingTo), a fresh draw of every gamma-group total
// from its exact law, one full evaluation and the long-double correction;
// rejection and exceptions restore the block-entry values, an exhausted
// interval the direction-entry values. The random order makes the sweep a
// reversible proposal on the surrogate, so one correction over the same child
// set at both states is exact; a fixed order with one correction is not.
//
// Coordinates: each member outside a gamma group is held on BayesTools'
// transformed scale (identity, log distance from a bound, logit of an
// interval); the members of a gamma group (GammaGroup.h) by the m - 1
// log-ratio coordinates of their proportions, with the group total T drawn
// from its Gamma(sum a, rate) law once per update. In the corrected update the
// total is drawn after the surrogate sweep and before the full evaluation, as
// an independence proposal from its exact marginal: the correction includes
// log Gamma(T) - log Gamma(T'), and the full densities include the Jacobian
// (m - 1) log T + sum_i log p_i at their own totals, so the total's terms
// cancel and the correction is the one of the proportions alone (in exact
// arithmetic). Drawing it there leaves the last full evaluation at the state
// the update ends in, which the selection nodes' memo then holds. The exact
// kernel draws the totals first (a Gibbs step) and then sweeps.
//
// Adaptation: while adapting, a rejected correction is followed by one exact
// sweep (the full conditional, the current axes, the entry totals) from the
// restored entry state, and the axes adapt to the state the update leaves.
// Until the axes are first learned (100 adaptive updates) the coordinate axes
// are used with the initial widths of the options (1 by default).
//
// Fallback: at the first adapt_off() the kernel switches for good to the exact
// kernel (BayesTools' fixed-order sweep of the full conditional with the
// learned axes, the totals drawn first) by the evidence rule of
// CorrectedTally, judged on acceptance^(1/d), d the number of directions.

#include "CoarseSurrogate.h"
#include "CorrectedSampler.h"
#include "GammaGroup.h"
#include "StepRowCache.h"

#include <BayesTools/factor_slice_kernel.h>

#include <sampler/GraphView.h>

#include <cstdint>
#include <string>
#include <vector>

namespace jags {
struct RNG;
namespace RoBMA {

struct CorrectedBlockOptions {
  // the start of every error message ("Corrected slice (block)")
  std::string label;
  // the corrected kind the sampler counts as (CorrectedSampler.h)
  int kind = corrected_kind_block;
  // test mutants (not exact kernels): the directions in their fixed order;
  // the gamma groups without the sum_i log p_i term of their Jacobian
  bool fixed_order = false;
  bool drop_sum_jacobian = false;
  std::vector<GammaGroup> groups;
  // widths of the coordinate axes until the axes are learned, one per
  // direction (plain members in block order, then each group's m - 1
  // log-ratio coordinates); empty for widths 1
  std::vector<double> initial_widths;
  // singleton step-selection children enter the surrogate through the row
  // cache (StepRowCache.h), refreshed at the start of each corrected update
  bool row_cache = false;
};

class CorrectedBlockKernel : public CorrectedSamplerWindow {
public:
  CorrectedBlockKernel(GraphView const *graph, unsigned int chain, SelNormCoarseSettings const &settings,
                       double min_acceptance, CorrectedBlockOptions const &options);
  ~CorrectedBlockKernel() override;
  CorrectedBlockKernel(CorrectedBlockKernel const &) = delete;
  CorrectedBlockKernel &operator=(CorrectedBlockKernel const &) = delete;

  // One update from the current member values (corrected, or exact after a
  // fallback); any exception restores the block-entry values.
  void update(RNG *rng);
  // The first call ends adaptation and decides the kernel.
  void adapt_off();
  bool adapting() const { return _tally.adapting(); }
  bool exact() const { return _exact; }
  // the number of directions
  unsigned int dimension() const { return _du; }
  // the number of adaptive updates the axes learned from
  unsigned int axis_updates() const { return _axes.updates(); }

  // Test hook: chain, the number of directions and the current widths.
  void append_widths(std::vector<double> &out) const;
  // Test hook of the row cache (frozen states of weights moves): refreshed at
  // the current state, then at that state and at points moved by -0.7 and
  // +0.7 along each coordinate axis, per point (chain, cache value, sum of the
  // rows' graph log densities); the current state is restored.
  void append_row_probe(std::vector<double> &out);

  std::uint64_t windowProposed() const override { return _tally.window_proposed(); }
  std::uint64_t windowAccepted() const override { return _tally.window_accepted(); }
  void clearWindow() override { _tally.clear_window(); }
  // A block move composes d direction moves: its rate is judged as rate^(1/d).
  double rateExponent() const override { return 1.0 / static_cast<double>(_du); }

private:
  GraphView const *_graph;
  const unsigned int _chain;
  const unsigned int _d;
  const SelNormCoarseSettings _settings;
  const double _min_acceptance;
  const CorrectedBlockOptions _options;
  // members outside the gamma groups (block positions) and their coordinates
  std::vector<unsigned int> _plain;
  std::vector<::BayesTools::factor_slice::Coordinate> _coordinates;
  unsigned int _du;
  // the current totals of the gamma groups
  std::vector<double> _totals;
  ::BayesTools::factor_slice::Axes _axes;
  // whether a refresh of the axes has replaced their initial state (a refresh
  // keeps the axes when an eigenvalue of the covariance is not positive)
  bool _learned = false;
  std::vector<double> _x;  // scratch: member values of a point
  std::vector<SelectedChild> _selected;
  std::vector<int> _child_selected;
  bool _exact;
  // whether a surrogate of the current update was not exact (evidence)
  bool _surrogate_not_exact;
  CorrectedTally _tally;
  // the step-selection rows of the surrogate (options.row_cache)
  std::vector<StochasticNode *> _row_nodes;
  StepRowCache _rows;

  void density_failure(Node const *node, char const *where, char const *cause) const;
  void adapt_axes(std::vector<double> const &u_start, std::vector<double> const &u);
  void check_sum(double out, char const *where) const;
  double member_terms(char const *where) const;
  double log_full(char const *where, LikelihoodMagnitude *magnitude = nullptr) const;
  double log_full_current() const;
  double log_surrogate(char const *where);
  void trace(int kind, double value) const;
  double width(unsigned int k) const;
  bool set_point(std::vector<double> const &u);
  double log_jacobian(std::vector<double> const &u) const;
  double surrogate_at(std::vector<double> const &u);
  double full_at(std::vector<double> const &u);
  std::vector<double> current_u();
  long double draw_totals(RNG *rng);
  std::vector<double> exact_sweep(RNG *rng);
  void exact_update(RNG *rng);
  void corrected_update(RNG *rng, std::vector<double> const &entry);
  void leave_live_blocks();
};

// Test hooks of both block samplers, not package interface. fixed_order:
// block samplers built from now on sweep the directions in their fixed order in
// corrected mode (a kernel that is not exact, for the power of the
// certification fixtures). trace: every density evaluation of a block sampler
// is recorded as (kind 0 surrogate / 1 full, chain, d, member values, log
// density without the Jacobian). The exhausted count is the number of slice
// intervals exhausted in floating point in corrected mode since the last take.
void corrected_block_test_set(bool fixed_order, bool trace);
bool corrected_block_test_fixed_order();
bool corrected_block_test_trace();
std::vector<double> corrected_block_test_take_trace();
std::uint64_t corrected_block_test_take_exhausted();
// The current widths of the live block samplers, in construction order: per
// sampler its chain, the number of directions and their widths.
std::vector<double> corrected_block_test_widths();
// The row probes (append_row_probe) of the live block samplers with a row
// cache, in construction order.
std::vector<double> corrected_block_test_row_probe();

}
}
#endif
