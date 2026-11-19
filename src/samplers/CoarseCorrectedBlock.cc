#include "CoarseCorrectedBlock.h"
#include "CoarseSurrogate.h"
#include "CorrectedSampler.h"

#include <BayesTools/factor_slice_kernel.h>

#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/MutableSampleMethod.h>
#include <sampler/MutableSampler.h>
#include <sampler/Sampler.h>

#include <atomic>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

// The univariate slice update, the coordinate transforms and the adaptation of
// the axes are BayesTools' factor slice kernel (installed header, LinkingTo).
namespace kernel = ::BayesTools::factor_slice;

const char *const block_sampler_name = "RoBMA::CoarseCorrectedBlock";
const double negative_infinity = -std::numeric_limits<double>::infinity();

// ---- the declarations of BayesTools::JAGS_fit() (src/samplers/
// BTSamplerDeclarations.cc of BayesTools, read the same way): a marker declares
// the node it is a direct stochastic child of only when it is observed with
// value 0 and its first parent is that node itself.

std::vector<StochasticNode const *> markers_of(StochasticNode *node, char const *name)
{
  std::vector<StochasticNode const *> out;
  for (StochasticNode const *child : *node->stochasticChildren()) {
    if (child->distribution()->name() != name) continue;
    std::vector<Node const *> const &parents = child->parents();
    if (parents.size() < 2 || parents[0] != node) continue;
    out.push_back(child);
  }
  return out;
}

bool observed_zero(StochasticNode const *marker)
{
  return isObserved(marker) && marker->length() == 1 && marker->value(0)[0] == 0.0;
}

bool has_gate_marker(StochasticNode *node)
{
  return !markers_of(node, "dbt_gate").empty();
}

// The id of the single valid dbt_block marker of the node (a fixed positive
// integer up to 1e9), or 0.
unsigned int block_id(StochasticNode *node)
{
  std::vector<StochasticNode const *> markers = markers_of(node, "dbt_block");
  if (markers.size() != 1 || !observed_zero(markers[0])) return 0;
  Node const *id = markers[0]->parents()[1];
  if (!id->isFixed() || id->length() != 1) return 0;
  const double value = id->value(0)[0];
  if (!std::isfinite(value) || value < 1.0 || value != std::floor(value) || value > 1.0e9) return 0;
  return static_cast<unsigned int>(value);
}

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

class CoarseCorrectedBlock : public MutableSampleMethod, public CorrectedSamplerWindow {
  GraphView const *_graph;
  const unsigned int _chain;
  const unsigned int _d;
  const SelNormCoarseSettings _settings;
  const double _min_acceptance;
  // Test hook: the direction order is not drawn (not an exact kernel).
  const bool _fixed_order;
  std::vector<kernel::Coordinate> _coordinates;
  kernel::Axes _axes;
  std::vector<double> _x;  // scratch: member values of a point
  // The coarse-capable selection children, and for every stochastic child of
  // the view (in the view's order) its index among them or -1 (an exact term).
  std::vector<SelectedChild> _selected;
  std::vector<int> _child_selected;
  // JAGS gives a sampler no accessor for the model's adaptation: true from
  // construction until the first adaptOff().
  bool _adapting;
  // The kernel chosen at the first adaptOff(): the exact fixed-order sweep of
  // the full conditional (BayesTools' factor slice update) for the rest of this
  // object's life.
  bool _exact;
  // Adaptive proposals and accepted corrections, and whether any surrogate
  // evaluated while adapting used a normalizer that is not the exact one.
  std::uint64_t _adaptive_proposed;
  std::uint64_t _adaptive_accepted;
  bool _surrogate_not_exact;
  // Proposals and accepted corrections after adaptOff() in corrected mode since
  // the last clearing read of the report (the window of one run).
  std::uint64_t _window_proposed;
  std::uint64_t _window_accepted;
  // Updates after adaptation (for the periodic resolvability check), and a
  // failure of the check at adaptOff(), raised at the next update().
  std::uint64_t _updates_after_adaptation;
  std::exception_ptr _deferred_error;

  static std::vector<kernel::Coordinate> coordinates(GraphView const *graph, unsigned int chain)
  {
    std::vector<kernel::Coordinate> out;
    for (StochasticNode *node : graph->nodes()) {
      // The support is fixed (checked by the factory), so it is read once.
      double lower, upper;
      node->support(&lower, &upper, 1, chain);
      out.push_back(kernel::Coordinate(lower, upper));
    }
    return out;
  }

  static void density_failure(Node const *node, char const *where, char const *cause)
  {
    throwNodeError(node, std::string(
      "Corrected slice (block): failure to calculate the log density of the block at ") + where + ": " + cause);
  }

  void check_sum(double out, char const *where) const
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
  double member_terms(char const *where) const
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
  double log_full(char const *where) const
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
  // surrogate (SelectedChild::logDensity). -Inf is outside the support; NaN and
  // +Inf are node errors.
  double log_surrogate(char const *where)
  {
    double out = member_terms(where);
    std::vector<StochasticNode *> const &children = _graph->stochasticChildren();
    for (std::size_t i = 0; i < children.size(); ++i) {
      const int selected = _child_selected[i];
      const double term = selected < 0 ? children[i]->logDensity(_chain, PDF_LIKELIHOOD) :
        _selected[selected].logDensity(_chain, _settings, &_surrogate_not_exact);
      if (std::isnan(term)) density_failure(children[i], where, "the surrogate log density term of this node is NaN");
      out += term;
    }
    check_sum(out, where);
    trace(0, out);
    return out;
  }

  void trace(int kind, double value) const
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

  // Sets the node values to g^-1(u) when every coordinate is representable
  // strictly inside its support.
  bool set_point(std::vector<double> const &u)
  {
    for (unsigned int i = 0; i < _d; ++i) {
      if (!std::isfinite(u[i]) || !_coordinates[i].to_x(u[i], &_x[i])) return false;
    }
    _graph->setValue(_x, _chain);
    return true;
  }

  // log q_u and log pi_u of a point u (the Jacobian included); -Inf when a
  // coordinate does not back-transform (nothing is clamped).
  double surrogate_at(std::vector<double> const &u)
  {
    if (!set_point(u)) return negative_infinity;
    return log_surrogate("a slice point") + kernel::log_jacobian(_coordinates, u);
  }

  double full_at(std::vector<double> const &u)
  {
    if (!set_point(u)) return negative_infinity;
    return log_full("a proposal") + kernel::log_jacobian(_coordinates, u);
  }

  // u = g(x) of the current member values; a value on a bound of its support
  // is a node error.
  std::vector<double> current_u()
  {
    _graph->getValue(_x, _chain);
    std::vector<double> u(_d);
    for (unsigned int i = 0; i < _d; ++i) {
      u[i] = _coordinates[i].to_u(_x[i]);
      if (!std::isfinite(u[i]) || !(_x[i] > _coordinates[i].lower() && _x[i] < _coordinates[i].upper())) {
        throwNodeError(_graph->nodes()[i], "Corrected slice (block): current value lies on a bound of its support");
      }
    }
    return u;
  }

  // Stops the fit when a transformed member's coordinate cannot resolve its
  // posterior at the current state (kernel::resolvable() on the full
  // conditional, as BayesTools' factor slice sampler does): the coordinate
  // reaches a value only to a relative precision of about |u| eps of its
  // distance from the anchoring bound, so with bounds far from where the
  // posterior concentrates the sampler would move on a grid coarser than the
  // posterior. Runs at construction, at every refresh of the axes, at the
  // first adaptOff() and every kernel::resolvability_interval updates after it, in
  // both kernels; it draws no random numbers and restores the node values, so
  // it changes no draw. A member on a bound is left to update().
  void check_resolvable(char const *site)
  {
    std::vector<StochasticNode *> const &nodes = _graph->nodes();
    std::vector<double> x(_d);
    _graph->getValue(x, _chain);
    std::vector<double> trial(x);
    for (unsigned int i = 0; i < _d; ++i) {
      kernel::Coordinate const &coordinate = _coordinates[i];
      if (coordinate.transform() == kernel::Transform::identity ||
          !(x[i] > coordinate.lower() && x[i] < coordinate.upper())) continue;
      bool ok;
      try {
        ok = kernel::resolvable(coordinate, x[i], [&](double value) {
          trial[i] = value;
          _graph->setValue(trial, _chain);
          return log_full("a resolvability check");
        });
      }
      catch (...) {
        _graph->setValue(x, _chain);
        throw;
      }
      trial[i] = x[i];
      _graph->setValue(x, _chain);
      if (!ok) {
        throwNodeError(nodes[i], std::string(
          "Corrected slice (block): the support bounds of this node are too far from where its posterior "
          "concentrates for the sampler's scale (found ") + site + ")");
      }
    }
  }

  // One completed proposal of the corrected update.
  void countOutcome(bool accepted)
  {
    if (_adapting) {
      ++_adaptive_proposed;
      if (accepted) ++_adaptive_accepted;
      return;
    }
    ++_window_proposed;
    corrected_sampler_counters.report_sampling_proposed.fetch_add(1, std::memory_order_relaxed);
    if (accepted) {
      ++_window_accepted;
      corrected_sampler_counters.report_sampling_accepted.fetch_add(1, std::memory_order_relaxed);
    }
  }

  // The exact kernel: BayesTools' factor slice update (BTFactorSlice.cc,
  // FactorSlice::update), the same calls in the same order, on the full
  // conditional with the axes learned while adapting.
  void exact_update(RNG *rng)
  {
    std::vector<double> u = current_u();
    for (unsigned int k = 0; k < _d; ++k) {
      std::vector<double> x_start(_d);
      _graph->getValue(x_start, _chain);
      const double lp0 = log_full("the current state") + kernel::log_jacobian(_coordinates, u);
      if (!std::isfinite(lp0)) {
        throwNodeError(_graph->nodes()[0],
          "Corrected slice (block, exact kernel): current value of the block has no finite density");
      }
      const kernel::SliceResult result = kernel::slice_along(
        u, _axes.direction(k), _axes.width(k), lp0,
        [this](std::vector<double> const &point) { return full_at(point); },
        [rng]() { return rng->uniform(); },
        [rng]() { return rng->exponential(); });
      if (!result.accepted) _graph->setValue(x_start, _chain);
    }
  }

  // The corrected update (plan: Kernel, steps 1-5), then the adaptation from
  // the retained state (a repeated state after a rejection) while adapting.
  void corrected_update(RNG *rng, std::vector<double> const &entry)
  {
    std::vector<double> u = current_u();
    const std::vector<double> u_start(u);
    const double jacobian_old = kernel::log_jacobian(_coordinates, u);
    // 1. both densities recomputed at block entry
    const double full_old = log_full("the current state") + jacobian_old;
    if (!std::isfinite(full_old)) {
      throwNodeError(_graph->nodes()[0], "Corrected slice (block) requires a finite current full density");
    }
    const double surrogate_old = log_surrogate("the current state") + jacobian_old;
    if (!std::isfinite(surrogate_old)) {
      throwNodeError(_graph->nodes()[0], "Corrected slice (block) requires a finite current surrogate density");
    }
    // 2. a fair draw of the sweep order
    const bool backward = !_fixed_order && rng->uniform() < 0.5;
    // 3. one slice update of q per direction, chaining the accepted values
    double surrogate = surrogate_old;
    std::vector<double> x_start(_d);
    for (unsigned int j = 0; j < _d; ++j) {
      const unsigned int k = backward ? _d - 1 - j : j;
      _graph->getValue(x_start, _chain);
      const kernel::SliceResult result = kernel::slice_along(
        u, _axes.direction(k), _axes.width(k), surrogate,
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
    corrected_sampler_counters.coarse_proposed.fetch_add(1, std::memory_order_relaxed);
    // 4. the full density at the end of the sweep (the node values hold u)
    const double full_new = log_full("a proposal") + kernel::log_jacobian(_coordinates, u);
    bool accepted;
    if (full_new == negative_infinity) {
      accepted = false;
    } else {
      // 5. the correction
      const long double log_ratio =
        (static_cast<long double>(full_new) - full_old) +
        (static_cast<long double>(surrogate_old) - surrogate);
      if (std::isnan(log_ratio)) {
        throwNodeError(_graph->nodes()[0], "Corrected slice (block) proposed an undefined acceptance ratio");
      }
      accepted = !(log_ratio < 0 && std::log(static_cast<long double>(rng->uniform())) > log_ratio);
    }
    if (accepted) {
      corrected_sampler_counters.coarse_accepted.fetch_add(1, std::memory_order_relaxed);
    } else {
      // the block-entry snapshot, never through the inverse transform
      _graph->setValue(entry, _chain);
      u = u_start;
      corrected_sampler_counters.coarse_correction_rejected.fetch_add(1, std::memory_order_relaxed);
    }
    countOutcome(accepted);
    if (_adapting && _axes.adapt(u_start, u)) check_resolvable("at a refresh of the axes during adaptation");
  }

public:
  CoarseCorrectedBlock(GraphView const *graph, unsigned int chain, SelNormCoarseSettings const &settings,
                       double min_acceptance, bool fixed_order)
    : _graph(graph), _chain(chain), _d(graph->length()), _settings(settings),
      _min_acceptance(min_acceptance), _fixed_order(fixed_order),
      _coordinates(coordinates(graph, chain)), _axes(_d), _x(_d), _adapting(true), _exact(false),
      _adaptive_proposed(0), _adaptive_accepted(0), _surrogate_not_exact(false),
      _window_proposed(0), _window_accepted(0), _updates_after_adaptation(0)
  {
    for (StochasticNode *child : graph->stochasticChildren()) {
      const ChildKind kind = child_kind(child, settings);
      if (kind == ChildKind::other) {
        _child_selected.push_back(-1);
      } else {
        _child_selected.push_back(static_cast<int>(_selected.size()));
        _selected.emplace_back(child, kind);
      }
    }
    if (_selected.empty()) throw std::logic_error("Corrected slice (block) has no selection child.");
    graph->checkFinite(chain);
    check_resolvable("at the start of the fit");
    register_corrected_sampler(this);
  }

  // A sampler destroyed before the report is read still counts for the worst
  // rate of its run.
  ~CoarseCorrectedBlock() override
  {
    unregister_corrected_sampler(this);
  }

  // Window counters of the report (read with the registry mutex held).
  std::uint64_t windowProposed() const override { return _window_proposed; }
  std::uint64_t windowAccepted() const override { return _window_accepted; }
  void clearWindow() override
  {
    _window_proposed = 0;
    _window_accepted = 0;
  }
  // A block move composes d direction moves: its rate is judged as rate^(1/d).
  double rateExponent() const override { return 1.0 / static_cast<double>(_d); }

  void update(RNG *rng) override
  {
    if (_deferred_error) std::rethrow_exception(_deferred_error);
    std::vector<double> entry(_d);
    _graph->getValue(entry, _chain);
    // Any exception restores the block-entry values before it propagates.
    try {
      if (_exact) exact_update(rng);
      else corrected_update(rng, entry);
    }
    catch (...) {
      _graph->setValue(entry, _chain);
      throw;
    }
    if (!_adapting && ++_updates_after_adaptation % kernel::resolvability_interval == 0) {
      check_resolvable("during sampling");
    }
  }

  bool isAdaptive() const override
  {
    return true;
  }

  // The kernel is chosen once, at the first adaptOff(), from the adaptive
  // updates: with the fallback enabled, the exact kernel when (a) fewer than 50
  // proposals were made, (c) no surrogate differed from the exact density by
  // construction, or (b) the per-direction acceptance rate^(1/d) is below
  // min_acceptance; the first rule that holds is the recorded reason. Later
  // calls change nothing (JAGS calls adaptOff() repeatedly). JAGS reports an
  // exception thrown here only as a message (rjags' adapt() returns normally),
  // so a failed resolvability check is raised at the next update().
  void adaptOff() override
  {
    if (!_adapting) return;
    _adapting = false;
    try {
      check_resolvable("at the end of adaptation");
    }
    catch (...) {
      _deferred_error = std::current_exception();
    }
    if (!(_min_acceptance > 0.0)) return;
    std::atomic<std::uint64_t> *reason = nullptr;
    if (_adaptive_proposed < corrected_min_adaptive_proposals) {
      reason = &corrected_sampler_counters.report_fallback_few_proposals;
    } else if (!_surrogate_not_exact) {
      reason = &corrected_sampler_counters.report_fallback_surrogate_exact;
    } else if (std::pow(static_cast<double>(_adaptive_accepted) / static_cast<double>(_adaptive_proposed),
                        rateExponent()) < _min_acceptance) {
      reason = &corrected_sampler_counters.report_fallback_low_acceptance;
    }
    if (reason != nullptr) {
      _exact = true;
      reason->fetch_add(1, std::memory_order_relaxed);
    }
  }

  // As BayesTools' factor slice sampler: incomplete until one adaptive update.
  bool checkAdaptation() const override
  {
    return _axes.updates() > 0;
  }
};

bool has_selection_child(GraphView const &view, SelNormCoarseSettings const &settings)
{
  for (StochasticNode *child : view.stochasticChildren()) {
    if (child_kind(child, settings) != ChildKind::other) return true;
  }
  return false;
}

}

std::vector<Sampler *> CoarseCorrectedBlockFactory::makeSamplers(
    std::list<StochasticNode *> const &nodes, Graph const &graph) const
{
  std::vector<Sampler *> samplers;
  const CorrectedSamplerConfig config = corrected_sampler_configuration();
  if (!config.enabled) return samplers;
  const bool fixed_order = test_fixed_order.load(std::memory_order_relaxed);

  // The members of BayesTools' factor slice factory (BTFactorSlice.cc): scalar,
  // continuous, unobserved nodes with df > 0 and a fixed support with
  // lower < upper, without a gate marker, with one valid block marker.
  std::map<unsigned int, std::vector<StochasticNode *> > blocks;
  for (StochasticNode *node : nodes) {
    if (node->length() != 1 || node->isDiscreteValued() || node->df() == 0 ||
        isObserved(node) || !isSupportFixed(node)) continue;
    if (has_gate_marker(node)) continue;
    const unsigned int id = block_id(node);
    if (id == 0) continue;
    double lower, upper;
    node->support(&lower, &upper, 1, 0);
    if (!(lower < upper)) continue;
    blocks[id].push_back(node);
  }

  const std::string sampler_name(block_sampler_name);
  try {
    for (auto const &block : blocks) {
      if (block.second.size() < 2) continue;
      std::unique_ptr<GraphView> view(new GraphView(block.second, graph, true));
      // a group without a coarse-capable selection child stays with BayesTools
      if (!has_selection_child(*view, config.settings)) continue;
      std::vector<std::unique_ptr<MutableSampleMethod> > owned_methods;
      std::vector<MutableSampleMethod *> methods;
      for (unsigned int chain = 0; chain < block.second[0]->nchain(); ++chain) {
        owned_methods.emplace_back(new CoarseCorrectedBlock(view.get(), chain, config.settings,
                                                            config.min_acceptance, fixed_order));
        methods.push_back(owned_methods.back().get());
      }
      // As in the scalar factory: the Sampler base adopts the graph, and
      // owned_methods keeps the methods until construction succeeds.
      std::unique_ptr<Sampler> sampler(new MutableSampler(view.release(), methods, sampler_name));
      for (std::unique_ptr<MutableSampleMethod> &method : owned_methods) method.release();
      samplers.push_back(sampler.get());
      sampler.release();
    }
  }
  catch (...) {
    for (Sampler *sampler : samplers) delete sampler;
    throw;
  }
  return samplers;
}

std::string CoarseCorrectedBlockFactory::name() const
{
  return block_sampler_name;
}

CoarseCorrectedBlockFactory *coarse_corrected_block_factory()
{
  static CoarseCorrectedBlockFactory factory;
  return &factory;
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

std::uint64_t corrected_block_test_take_exhausted()
{
  return test_exhausted.exchange(0, std::memory_order_relaxed);
}

}
}
