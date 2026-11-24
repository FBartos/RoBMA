#include "CorrectedBlockKernel.h"


#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>

#include <JRmath.h>

namespace jags {
namespace RoBMA {
namespace {

// The univariate slice update, the coordinate transforms and the adaptation of
// the axes are BayesTools' factor slice kernel (installed header, LinkingTo).
namespace kernel = ::BayesTools::factor_slice;

const double negative_infinity = -std::numeric_limits<double>::infinity();

// ---- test hooks

std::atomic<bool> test_fixed_order(false);
std::atomic<bool> test_trace(false);
std::atomic<std::uint64_t> test_exhausted(0);

struct TraceStore {
  std::mutex mutex;
  std::vector<double> records;
};

// Never destroyed, like the sampler registry: samplers may run late.
TraceStore &trace_store()
{
  static TraceStore *store = new TraceStore();
  return *store;
}

// The live block kernels, for the test hook that reads their widths. Never
// destroyed, like the sampler registry: samplers may be destroyed late.
struct LiveBlocks {
  std::mutex mutex;
  std::vector<CorrectedBlockKernel const *> live;
};

LiveBlocks &live_blocks()
{
  static LiveBlocks *blocks = new LiveBlocks();
  return *blocks;
}

std::vector<bool> grouped_members(unsigned int d, std::vector<GammaGroup> const &groups)
{
  std::vector<bool> grouped(d, false);
  for (GammaGroup const &group : groups) {
    if (group.members.size() < 2) throw std::logic_error("A gamma group needs two members.");
    for (unsigned int index : group.members) {
      if (index >= d || grouped[index]) throw std::logic_error("Invalid gamma group of a block.");
      grouped[index] = true;
    }
  }
  return grouped;
}

}

CorrectedBlockKernel::CorrectedBlockKernel(GraphView const *graph, unsigned int chain,
                                           SelNormCoarseSettings const &settings, double min_acceptance,
                                           CorrectedBlockOptions const &options)
  : _graph(graph), _chain(chain), _d(graph->length()), _settings(settings),
    _min_acceptance(min_acceptance), _options(options), _du(0), _totals(options.groups.size(), 0.0),
    _axes(1), _x(graph->length()), _exact(false), _surrogate_not_exact(false), _tally(options.kind)
{
  std::vector<StochasticNode *> const &nodes = graph->nodes();
  const std::vector<bool> grouped = grouped_members(_d, options.groups);
  for (unsigned int i = 0; i < _d; ++i) {
    if (grouped[i]) continue;
    // The support is fixed (checked by the factory), so it is read once.
    double lower, upper;
    nodes[i]->support(&lower, &upper, 1, chain);
    _plain.push_back(i);
    _coordinates.push_back(kernel::Coordinate(lower, upper));
  }
  _du = static_cast<unsigned int>(_plain.size());
  for (GammaGroup const &group : options.groups) {
    _du += static_cast<unsigned int>(group.members.size()) - 1;
  }
  if (_du == 0) throw std::logic_error("A corrected block needs a direction.");
  if (!options.initial_widths.empty() && options.initial_widths.size() != _du) {
    throw std::logic_error("The initial widths do not match the block's directions.");
  }
  _axes = kernel::Axes(_du);
  for (StochasticNode *child : graph->stochasticChildren()) {
    const ChildKind kind = child_kind(child, settings);
    if (kind == ChildKind::other && options.row_cache && is_step_switch_row(child)) {
      // a row of the row cache (-2)
      _child_selected.push_back(-2);
      _row_nodes.push_back(child);
    } else if (kind == ChildKind::other) {
      _child_selected.push_back(-1);
    } else {
      _child_selected.push_back(static_cast<int>(_selected.size()));
      _selected.emplace_back(child, kind);
    }
  }
  if (_selected.empty() && _row_nodes.empty()) throw std::logic_error(options.label + " has no selection child.");
  _rows = StepRowCache(_row_nodes, chain);
  graph->checkFinite(chain);
  // The one limit of the sampler's scale that the bounds alone decide (as in
  // BayesTools' factor slice sampler): a two-sided member whose span
  // upper - lower is not finite cannot be held through its distance from a
  // bound. Other precision limits (bounds very far from where the posterior
  // concentrates) are not detected; see the kernel header.
  for (std::size_t k = 0; k < _plain.size(); ++k) {
    kernel::Coordinate const &coordinate = _coordinates[k];
    if (coordinate.transform() == kernel::Transform::interval &&
        !std::isfinite(coordinate.upper() - coordinate.lower())) {
      throwNodeError(nodes[_plain[k]], _options.label +
        ": the support bounds of this node are too wide for the sampler's scale "
        "(their difference exceeds the largest double)");
    }
  }
  // Joins the live list first and the sampler registry last: a constructor
  // that throws runs no destructor, so a failed registration leaves the live
  // list again and neither list keeps a pointer to this kernel.
  {
    LiveBlocks &blocks = live_blocks();
    const std::lock_guard<std::mutex> lock(blocks.mutex);
    blocks.live.push_back(this);
  }
  try {
    register_corrected_sampler(this, options.kind);
  } catch (...) {
    leave_live_blocks();
    throw;
  }
}

// A sampler destroyed before the report is read still counts for the worst
// rate of its run.
CorrectedBlockKernel::~CorrectedBlockKernel()
{
  unregister_corrected_sampler(this);
  leave_live_blocks();
}

void CorrectedBlockKernel::leave_live_blocks()
{
  LiveBlocks &blocks = live_blocks();
  const std::lock_guard<std::mutex> lock(blocks.mutex);
  std::vector<CorrectedBlockKernel const *> &live = blocks.live;
  for (std::size_t i = 0; i < live.size(); ++i) {
    if (live[i] == this) {
      live.erase(live.begin() + static_cast<std::ptrdiff_t>(i));
      break;
    }
  }
}

void CorrectedBlockKernel::append_widths(std::vector<double> &out) const
{
  out.push_back(_chain);
  out.push_back(_du);
  for (unsigned int k = 0; k < _du; ++k) out.push_back(width(k));
}

void CorrectedBlockKernel::append_row_probe(std::vector<double> &out)
{
  if (_rows.empty()) return;
  std::vector<double> entry(_d);
  _graph->getValue(entry, _chain);
  const std::vector<double> u = current_u();
  _rows.refresh();
  auto record = [&]() {
    double graph = 0;
    for (StochasticNode const *row : _row_nodes) graph += row->logDensity(_chain, PDF_LIKELIHOOD);
    out.push_back(_chain);
    out.push_back(_rows.value(nullptr));
    out.push_back(graph);
  };
  record();
  for (unsigned int k = 0; k < _du; ++k) {
    for (double offset : {-0.7, 0.7}) {
      std::vector<double> moved(u);
      moved[k] += offset;
      if (set_point(moved)) record();
    }
  }
  _graph->setValue(entry, _chain);
}

void CorrectedBlockKernel::density_failure(Node const *node, char const *where, char const *cause) const
{
  throwNodeError(node, _options.label + ": failure to calculate the log density of the block at " +
                 where + ": " + cause);
}

void CorrectedBlockKernel::check_sum(double out, char const *where) const
{
  if (std::isnan(out) || out > std::numeric_limits<double>::max()) {
    density_failure(_graph->nodes()[0], where, std::isnan(out) ?
      "the summed log density of the block whose first member is this node is NaN" :
      "the summed log density of the block whose first member is this node is +Inf");
  }
}

// The log density terms of the members: every member's full log density (a
// member may be a stochastic child of another member). The same terms, in
// the same order, as BayesTools' factor slice sampler.
double CorrectedBlockKernel::member_terms(char const *where) const
{
  double out = 0.0;
  for (StochasticNode const *node : _graph->nodes()) {
    const double term = node->logDensity(_chain, PDF_FULL);
    if (std::isnan(term)) density_failure(node, where, "the log density term of this node is NaN");
    out += term;
  }
  return out;
}

// log pi(x | rest) at the current node values (without the Jacobian): the
// members' terms and the likelihood terms of every stochastic child of the
// block, in the order of BayesTools' factor slice sampler. A NaN term, a NaN
// sum or a +Inf sum is a node error.
double CorrectedBlockKernel::log_full(char const *where) const
{
  double out = member_terms(where);
  for (StochasticNode const *child : _graph->stochasticChildren()) {
    const double term = child->logDensity(_chain, PDF_LIKELIHOOD);
    if (std::isnan(term)) density_failure(child, where, "the log density term of this node is NaN");
    out += term;
  }
  check_sum(out, where);
  trace(1, out);
  return out;
}

// log q(x | rest) at the current node values (without the Jacobian): the full
// conditional with each coarse-capable selection child's term replaced by its
// surrogate (SelectedChild::logDensity), which also flags a surrogate that is
// not exact. -Inf is outside the support; NaN and +Inf are node errors.
double CorrectedBlockKernel::log_surrogate(char const *where)
{
  double out = member_terms(where);
  std::vector<StochasticNode *> const &children = _graph->stochasticChildren();
  for (std::size_t i = 0; i < children.size(); ++i) {
    const int selected = _child_selected[i];
    if (selected == -2) continue;
    const double term = selected < 0 ? children[i]->logDensity(_chain, PDF_LIKELIHOOD) :
      _selected[selected].logDensity(_chain, _settings, &_surrogate_not_exact);
    if (std::isnan(term)) density_failure(children[i], where, "the surrogate log density term of this node is NaN");
    out += term;
  }
  if (!_rows.empty()) {
    StochasticNode const *failed = nullptr;
    const double rows = _rows.value(&failed);
    if (std::isnan(rows)) density_failure(failed, where, "the surrogate log density term of this node is NaN");
    out += rows;
  }
  check_sum(out, where);
  trace(0, out);
  return out;
}

void CorrectedBlockKernel::trace(int kind, double value) const
{
  if (!test_trace.load(std::memory_order_relaxed)) return;
  std::vector<double> x(_d);
  _graph->getValue(x, _chain);
  TraceStore &store = trace_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  store.records.push_back(kind);
  store.records.push_back(_chain);
  store.records.push_back(_d);
  store.records.insert(store.records.end(), x.begin(), x.end());
  store.records.push_back(value);
}

// The width of direction k: the learned one once the axes were refreshed, the
// initial width of the coordinate axis before.
double CorrectedBlockKernel::width(unsigned int k) const
{
  if (_axes.updates() >= kernel::first_refresh || _options.initial_widths.empty()) return _axes.width(k);
  return _options.initial_widths[k];
}

// Sets the node values to the point u (the gamma groups at their current
// totals) when every member is representable strictly inside its support.
bool CorrectedBlockKernel::set_point(std::vector<double> const &u)
{
  for (std::size_t k = 0; k < _plain.size(); ++k) {
    if (!std::isfinite(u[k]) || !_coordinates[k].to_x(u[k], &_x[_plain[k]])) return false;
  }
  std::size_t k = _plain.size();
  for (std::size_t g = 0; g < _options.groups.size(); ++g) {
    std::vector<unsigned int> const &members = _options.groups[g].members;
    const std::size_t m = members.size();
    double largest = 0.0;
    for (std::size_t i = 0; i + 1 < m; ++i) {
      if (!std::isfinite(u[k + i])) return false;
      largest = std::max(largest, u[k + i]);
    }
    double total = std::exp(-largest);
    for (std::size_t i = 0; i + 1 < m; ++i) total += std::exp(u[k + i] - largest);
    const double log_total = std::log(total);
    for (std::size_t i = 0; i < m; ++i) {
      const double log_p = (i + 1 < m ? u[k + i] : 0.0) - largest - log_total;
      const double x = _totals[g] * std::exp(log_p);
      if (!std::isfinite(x) || !(x > 0.0)) return false;
      _x[members[i]] = x;
    }
    k += m - 1;
  }
  _graph->setValue(_x, _chain);
  return true;
}

// log |dx/du| at u: the plain coordinates' terms, and for each gamma group
// (m - 1) log T + sum_i log p_i at its current total (the test mutant drops the
// sum over the proportions).
double CorrectedBlockKernel::log_jacobian(std::vector<double> const &u) const
{
  double out = 0.0;
  for (std::size_t k = 0; k < _plain.size(); ++k) out += _coordinates[k].log_jacobian(u[k]);
  std::size_t k = _plain.size();
  for (std::size_t g = 0; g < _options.groups.size(); ++g) {
    const std::size_t m = _options.groups[g].members.size();
    out += static_cast<double>(m - 1) * std::log(_totals[g]);
    if (!_options.drop_sum_jacobian) {
      double largest = 0.0;
      for (std::size_t i = 0; i + 1 < m; ++i) largest = std::max(largest, u[k + i]);
      double total = std::exp(-largest);
      for (std::size_t i = 0; i + 1 < m; ++i) total += std::exp(u[k + i] - largest);
      const double log_total = std::log(total);
      for (std::size_t i = 0; i + 1 < m; ++i) out += u[k + i] - largest - log_total;
      out += -largest - log_total;
    }
    k += m - 1;
  }
  return out;
}

// log q_u and log pi_u of a point u (the Jacobian included); -Inf when a
// coordinate does not back-transform (nothing is clamped).
double CorrectedBlockKernel::surrogate_at(std::vector<double> const &u)
{
  if (!set_point(u)) return negative_infinity;
  return log_surrogate("a slice point") + log_jacobian(u);
}

double CorrectedBlockKernel::full_at(std::vector<double> const &u)
{
  if (!set_point(u)) return negative_infinity;
  return log_full("a proposal") + log_jacobian(u);
}

// u of the current member values, and the totals of the gamma groups; a value
// on a bound of its support is a node error.
std::vector<double> CorrectedBlockKernel::current_u()
{
  _graph->getValue(_x, _chain);
  std::vector<StochasticNode *> const &nodes = _graph->nodes();
  std::vector<double> u(_du);
  for (std::size_t k = 0; k < _plain.size(); ++k) {
    const double x = _x[_plain[k]];
    u[k] = _coordinates[k].to_u(x);
    if (!std::isfinite(u[k]) || !(x > _coordinates[k].lower() && x < _coordinates[k].upper())) {
      throwNodeError(nodes[_plain[k]], _options.label + ": current value lies on a bound of its support");
    }
  }
  std::size_t k = _plain.size();
  for (std::size_t g = 0; g < _options.groups.size(); ++g) {
    std::vector<unsigned int> const &members = _options.groups[g].members;
    const std::size_t m = members.size();
    double total = 0.0;
    for (unsigned int index : members) {
      if (!std::isfinite(_x[index]) || !(_x[index] > 0.0)) {
        throwNodeError(nodes[index], _options.label + ": current value lies on a bound of its support");
      }
      total += _x[index];
    }
    const double log_last = std::log(_x[members[m - 1]]);
    for (std::size_t i = 0; i + 1 < m; ++i) u[k + i] = std::log(_x[members[i]]) - log_last;
    if (!std::isfinite(total)) {
      throwNodeError(nodes[members[0]], _options.label + ": the sum of the gamma group of this node is not finite");
    }
    _totals[g] = total;
    k += m - 1;
  }
  return u;
}

// New totals of the gamma groups from their Gamma(sum a, rate) laws; returns
// the log density ratio log Gamma(T) - log Gamma(T') of the old to the new.
long double CorrectedBlockKernel::draw_totals(RNG *rng)
{
  long double out = 0.0L;
  for (std::size_t g = 0; g < _options.groups.size(); ++g) {
    GammaGroup const &group = _options.groups[g];
    const double scale = 1.0 / group.rate;
    const double drawn = rgamma(group.shape_sum, scale, rng);
    out += static_cast<long double>(dgamma(_totals[g], group.shape_sum, scale, 1)) -
      dgamma(drawn, group.shape_sum, scale, 1);
    _totals[g] = drawn;
  }
  return out;
}

// One fixed-order sweep of the full conditional with the current axes and
// totals: BayesTools' factor slice update (BTFactorSlice.cc,
// FactorSlice::update), the same calls in the same order. Returns the u of the
// state it leaves.
std::vector<double> CorrectedBlockKernel::exact_sweep(RNG *rng)
{
  std::vector<double> u = current_u();
  std::vector<double> x_start(_d);
  for (unsigned int k = 0; k < _du; ++k) {
    _graph->getValue(x_start, _chain);
    const double lp0 = log_full("the current state") + log_jacobian(u);
    if (!std::isfinite(lp0)) {
      throwNodeError(_graph->nodes()[0], _options.label +
        " (exact kernel): current value of the block has no finite density");
    }
    const kernel::SliceResult result = kernel::slice_along(
      u, _axes.direction(k), width(k), lp0,
      [this](std::vector<double> const &point) { return full_at(point); },
      [rng]() { return rng->uniform(); },
      [rng]() { return rng->exponential(); });
    if (!result.accepted) _graph->setValue(x_start, _chain);
  }
  return u;
}

// The exact kernel after a fallback: the gamma-group totals from their law (a
// Gibbs step; a total whose proportions are not representable keeps the old
// one), then the exact sweep.
void CorrectedBlockKernel::exact_update(RNG *rng)
{
  if (!_options.groups.empty()) {
    const std::vector<double> u = current_u();
    const std::vector<double> totals(_totals);
    draw_totals(rng);
    if (!set_point(u)) {
      _totals = totals;
      if (!set_point(u)) {
        throwNodeError(_graph->nodes()[0], _options.label + ": the block's current value is not representable");
      }
    }
  }
  exact_sweep(rng);
}

void CorrectedBlockKernel::corrected_update(RNG *rng, std::vector<double> const &entry)
{
  std::vector<double> u = current_u();
  const std::vector<double> u_start(u);
  const std::vector<double> totals_start(_totals);
  const double jacobian_old = log_jacobian(u);
  // the row cache of this update, at the entry state
  if (!_rows.empty()) _rows.refresh();
  // 1. both densities recomputed at block entry
  const double full_old = log_full("the current state") + jacobian_old;
  if (!std::isfinite(full_old)) {
    throwNodeError(_graph->nodes()[0], _options.label + " requires a finite current full density");
  }
  const double surrogate_old = log_surrogate("the current state") + jacobian_old;
  if (!std::isfinite(surrogate_old)) {
    throwNodeError(_graph->nodes()[0], _options.label + " requires a finite current surrogate density");
  }
  // 2. a fair draw of the sweep order
  const bool backward = !_options.fixed_order && rng->uniform() < 0.5;
  // 3. one slice update of q per direction, chaining the accepted values
  double surrogate = surrogate_old;
  std::vector<double> x_start(_d);
  for (unsigned int j = 0; j < _du; ++j) {
    const unsigned int k = backward ? _du - 1 - j : j;
    _graph->getValue(x_start, _chain);
    const kernel::SliceResult result = kernel::slice_along(
      u, _axes.direction(k), width(k), surrogate,
      [this](std::vector<double> const &point) { return surrogate_at(point); },
      [rng]() { return rng->uniform(); },
      [rng]() { return rng->exponential(); });
    if (result.accepted) {
      surrogate = result.log_density;
    } else {
      // floating-point exhaustion: a self-transition of this direction
      _graph->setValue(x_start, _chain);
      test_exhausted.fetch_add(1, std::memory_order_relaxed);
    }
  }
  // 4. new gamma-group totals from their law (an independence proposal), then
  // the full density at the end of the sweep
  long double log_total_ratio = 0.0L;
  bool representable = true;
  if (!_options.groups.empty()) {
    log_total_ratio = draw_totals(rng);
    representable = set_point(u);
  }
  const double full_new = representable ?
    log_full("a proposal") + log_jacobian(u) : negative_infinity;
  bool accepted;
  // evidence of the row cache: the surrogate and exact changes differ beyond
  // rounding (StepRowCache.h)
  bool row_evidence = false;
  if (full_new == negative_infinity) {
    accepted = false;
    // a finite surrogate where the exact density is zero
    row_evidence = !_rows.empty() && representable;
  } else {
    // 5. the correction
    const long double log_ratio =
      (static_cast<long double>(full_new) - full_old) +
      (static_cast<long double>(surrogate_old) - surrogate) + log_total_ratio;
    if (std::isnan(log_ratio)) {
      throwNodeError(_graph->nodes()[0], _options.label + " proposed an undefined acceptance ratio");
    }
    accepted = !(log_ratio < 0 && std::log(static_cast<long double>(rng->uniform())) > log_ratio);
    row_evidence = !_rows.empty() && step_row_cache_mismatch(log_ratio, full_old, full_new);
  }
  if (!accepted) {
    // the block-entry snapshot, never through the inverse transform
    _graph->setValue(entry, _chain);
    u = u_start;
    _totals = totals_start;
    // While adapting, a rejected correction is followed by one exact sweep
    // (the full conditional, the current axes, the entry totals) from the
    // restored entry state. A chain that starts where the surrogate is poor
    // would otherwise reject over and over and not move, and axes learned from
    // repeated states collapse, freezing the chain after adaptation too. The
    // adaptive phase need not leave the target invariant (its draws are
    // discarded). The proposal still counts as rejected, but the adaptive chain
    // leaves regions where the surrogate is poor faster than the corrected
    // kernel alone would, so the adaptive rate that the low-acceptance rule
    // reads is optimistic. A block kept on that reading is still exact (the
    // rule guards efficiency only), and the per-fit warning judges the rate
    // after adaptation.
    if (_tally.adapting()) u = exact_sweep(rng);
  }
  const bool adapting = _tally.adapting();
  _tally.count(accepted, _surrogate_not_exact || row_evidence);
  if (adapting) _axes.adapt(u_start, u);
}

void CorrectedBlockKernel::update(RNG *rng)
{
  std::vector<double> entry(_d);
  _graph->getValue(entry, _chain);
  _surrogate_not_exact = false;
  // Any exception restores the block-entry values before it propagates.
  try {
    if (_exact) exact_update(rng);
    else corrected_update(rng, entry);
  }
  catch (...) {
    _graph->setValue(entry, _chain);
    throw;
  }
}

void CorrectedBlockKernel::adapt_off()
{
  if (_tally.adapt_off(_min_acceptance, rateExponent())) _exact = true;
}

void corrected_block_test_set(bool fixed_order, bool trace)
{
  test_fixed_order.store(fixed_order, std::memory_order_relaxed);
  test_trace.store(trace, std::memory_order_relaxed);
}

bool corrected_block_test_fixed_order()
{
  return test_fixed_order.load(std::memory_order_relaxed);
}

bool corrected_block_test_trace()
{
  return test_trace.load(std::memory_order_relaxed);
}

std::vector<double> corrected_block_test_take_trace()
{
  TraceStore &store = trace_store();
  const std::lock_guard<std::mutex> lock(store.mutex);
  std::vector<double> out;
  out.swap(store.records);
  return out;
}

std::vector<double> corrected_block_test_widths()
{
  LiveBlocks &blocks = live_blocks();
  const std::lock_guard<std::mutex> lock(blocks.mutex);
  std::vector<double> out;
  for (CorrectedBlockKernel const *block : blocks.live) block->append_widths(out);
  return out;
}

std::vector<double> corrected_block_test_row_probe()
{
  LiveBlocks &blocks = live_blocks();
  const std::lock_guard<std::mutex> lock(blocks.mutex);
  std::vector<double> out;
  for (CorrectedBlockKernel const *block : blocks.live) {
    // the live list holds const pointers; the probe restores the state it moves
    const_cast<CorrectedBlockKernel *>(block)->append_row_probe(out);
  }
  return out;
}

std::uint64_t corrected_block_test_take_exhausted()
{
  return test_exhausted.exchange(0, std::memory_order_relaxed);
}

}
}
