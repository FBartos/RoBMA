#include "CoarseCorrectedSlice.h"
#include "../distributions/DSELNORMMVSTEP.h"
#include "../distributions/DSELNORMCLUSTERSTEP.h"
#include "../distributions/DSELNORMFACTORSTEP.h"
#include "../selnorm/selnorm.h"

#include <sampler/Slicer.h>
#include <sampler/SingletonGraphView.h>
#include <sampler/MutableSampler.h>
#include <sampler/Sampler.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>

#include <cmath>
#include <atomic>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

std::atomic<std::uint64_t> coarse_proposed(0);
std::atomic<std::uint64_t> coarse_accepted(0);
std::atomic<std::uint64_t> coarse_correction_rejected(0);

// Per-fit counters (CoarseCorrectedSliceReport).
std::atomic<std::uint64_t> report_samplers(0);
std::atomic<std::uint64_t> report_fallback_few_proposals(0);
std::atomic<std::uint64_t> report_fallback_surrogate_exact(0);
std::atomic<std::uint64_t> report_fallback_low_acceptance(0);
std::atomic<std::uint64_t> report_sampling_proposed(0);
std::atomic<std::uint64_t> report_sampling_accepted(0);

// Fallback rules decided at adaptOff(): fewer adaptive proposals than this
// cannot be judged. The worst-sampler rate counts samplers with at least
// min_rate_proposals proposals since the last clearing read of the report.
const std::uint64_t min_adaptive_proposals = 50;
const std::uint64_t min_rate_proposals = 1000;

// The rate of one sampler's window and whether it counts for the worst rate.
bool window_rate(std::uint64_t proposed, std::uint64_t accepted, double *rate)
{
  if (proposed < min_rate_proposals) return false;
  *rate = static_cast<double>(accepted) / static_cast<double>(proposed);
  return true;
}

// How a stochastic child of the sampled node enters the surrogate density: a
// selection block with a surrogate evaluator (the dense multivariate block, or a
// cluster or factor block with a coarse rung), or an exact term.
enum class ChildKind { other, mv, cluster, factor };

// Parameter positions of the cluster and factor distributions (they share a
// layout up to the loading): kernel mode, vector rule, quadrature orders.
const unsigned int cluster_factor_mode = 10;
const unsigned int cluster_factor_orders = 13;
const unsigned int cluster_factor_rule = 19;

// A cluster or factor child is coarse-capable when its fixed controls allow a
// coarse rule: a vector rule not fixed at a best-weight value, a kernel mode not
// fixed at NORMAL, and, for a factor block, rule min(max_rules, rules) of its
// schedule affordable at the uncompacted rank (a compacted block or a forest
// only costs less). In model-averaged fits the vector rule is a node of the bias
// indicator that is product in every state of a product-only model, so a
// non-fixed rule counts as capable; in a state with a best-weight rule the
// child's surrogate is its exact density. The decision reads fixed parents only,
// so it is the same in every chain. Any other cluster or factor child stays an
// exact term.
ChildKind child_kind(StochasticNode const *child, SelNormCoarseSettings const &settings)
{
  Distribution const *distribution = child->distribution();
  if (dynamic_cast<DSELNORMMVSTEP const *>(distribution) != nullptr) return ChildKind::mv;
  const bool cluster = dynamic_cast<DSELNORMCLUSTERSTEP const *>(distribution) != nullptr;
  const bool factor = !cluster && dynamic_cast<DSELNORMFACTORSTEP const *>(distribution) != nullptr;
  if (!cluster && !factor) return ChildKind::other;
  std::vector<Node const *> const &parents = child->parents();
  if (parents.size() < distribution->npar()) return ChildKind::other;
  Node const *rule = parents[cluster_factor_rule];
  if (rule->isFixed() && rule->value(0)[0] != SELVECTOR_PRODUCT) return ChildKind::other;
  Node const *mode = parents[cluster_factor_mode];
  if (mode->isFixed() && mode->value(0)[0] == SELKERNEL_NORMAL) return ChildKind::other;
  Node const *orders = parents[cluster_factor_orders];
  if (!orders->isFixed() || orders->length() < 1) return ChildKind::other;
  if (factor) {
    const unsigned int rank = parents[2]->length() / parents[0]->length();
    const unsigned int rung = std::min(settings.max_rules, orders->length());
    if (rung < 1 || !selnorm_factor_rule_affordable(
          static_cast<int>(orders->value(0)[rung - 1]), static_cast<int>(rank))) {
      return ChildKind::other;
    }
  }
  return cluster ? ChildKind::cluster : ChildKind::factor;
}

struct SelectedChild {
  StochasticNode *node;
  ChildKind kind;
  // The distribution of the child, of the class the kind names (checked once).
  Distribution const *distribution;
  std::vector<Node const *> parents;
  std::vector<unsigned int> lengths;
  mutable std::vector<double const *> parameters;
  bool controls_fixed = true;
  mutable bool controls_checked = false;
  mutable bool controls_valid = false;

  // The parents that may change during sampling: mean, covariance (residual SD
  // and loadings), weights, kernel mode and, for cluster and factor children,
  // the vector rule (a mixture of weight rules composes it from the bias
  // indicator). Everything else is a control.
  static bool dynamic_parent(ChildKind kind, unsigned int i)
  {
    if (kind == ChildKind::mv) return i == 0 || i == 1 || i == 3 || i == 9;
    return i == 0 || i == 1 || i == 2 || i == 4 || i == cluster_factor_mode ||
      i == cluster_factor_rule;
  }

  SelectedChild(StochasticNode *child, ChildKind kind_of_child)
    : node(child), kind(kind_of_child), distribution(child->distribution())
  {
    const bool consistent =
      (kind == ChildKind::mv && dynamic_cast<DSELNORMMVSTEP const *>(distribution) != nullptr) ||
      (kind == ChildKind::cluster && dynamic_cast<DSELNORMCLUSTERSTEP const *>(distribution) != nullptr) ||
      (kind == ChildKind::factor && dynamic_cast<DSELNORMFACTORSTEP const *>(distribution) != nullptr);
    if (!consistent) {
      throw std::logic_error("Selection child has no surrogate evaluator.");
    }
    const unsigned int count = distribution->npar();
    if (child->parents().size() < count) {
      throw std::logic_error("Selection child is missing distribution parameters.");
    }
    // StochasticNode appends truncation bounds after its distribution parents.
    parents.assign(child->parents().begin(), child->parents().begin() + count);
    lengths.resize(count);
    parameters.resize(count);
    for (unsigned int i = 0; i < count; ++i) {
      lengths[i] = parents[i]->length();
      if (!dynamic_parent(kind, i) && !parents[i]->isFixed()) {
        controls_fixed = false;
      }
    }
  }

  bool checkControl() const
  {
    switch (kind) {
    case ChildKind::mv:
      return static_cast<DSELNORMMVSTEP const *>(distribution)->checkControlParameterValue(parameters, lengths);
    case ChildKind::cluster:
      return static_cast<DSELNORMCLUSTERSTEP const *>(distribution)->checkControlParameterValue(parameters, lengths);
    default:
      return static_cast<DSELNORMFACTORSTEP const *>(distribution)->checkControlParameterValue(parameters, lengths);
    }
  }

  bool checkDynamic() const
  {
    switch (kind) {
    case ChildKind::mv:
      return static_cast<DSELNORMMVSTEP const *>(distribution)->checkDynamicParameterValue(parameters, lengths);
    case ChildKind::cluster:
      return static_cast<DSELNORMCLUSTERSTEP const *>(distribution)->checkDynamicParameterValue(parameters, lengths);
    default:
      return static_cast<DSELNORMFACTORSTEP const *>(distribution)->checkDynamicParameterValue(parameters, lengths);
    }
  }

  double surrogate(unsigned int chain, SelNormCoarseSettings const &settings,
                   bool *not_exact) const
  {
    double const *x = node->value(chain);
    double const *lower = node->lowerLimit(chain);
    double const *upper = node->upperLimit(chain);
    switch (kind) {
    case ChildKind::mv:
      return static_cast<DSELNORMMVSTEP const *>(distribution)->surrogateLogDensity(
        x, node->length(), PDF_LIKELIHOOD, parameters, lengths, lower, upper, settings,
        not_exact);
    case ChildKind::cluster:
      return static_cast<DSELNORMCLUSTERSTEP const *>(distribution)->surrogateLogDensity(
        x, node->length(), PDF_LIKELIHOOD, parameters, lengths, lower, upper, settings,
        not_exact);
    default:
      return static_cast<DSELNORMFACTORSTEP const *>(distribution)->surrogateLogDensity(
        x, node->length(), PDF_LIKELIHOOD, parameters, lengths, lower, upper, settings,
        not_exact);
    }
  }

  // not_exact is set when the surrogate used a normalizer that is not exact.
  double logDensity(unsigned int chain, SelNormCoarseSettings const &settings,
                    bool *not_exact) const
  {
    for (unsigned int i = 0; i < parents.size(); ++i) {
      parameters[i] = parents[i]->value(chain);
    }
    bool valid;
    if (controls_fixed) {
      // The first surrogate call follows initialized-graph checkFinite().
      // This certificate belongs to one SelectedChild in one chain's sampler.
      if (!controls_checked) {
        controls_valid = checkControl();
        controls_checked = true;
      }
      valid = controls_valid && checkDynamic();
    } else {
      // Any nonfixed control retains the complete original node gate.
      valid = node->checkParentValues(chain);
    }
    if (!valid) return -std::numeric_limits<double>::infinity();
    return surrogate(chain, settings, not_exact);
  }
};

class CoarseCorrectedSlice;

// The live corrected samplers of this process, and the lowest window rate of
// those destroyed since the last clearing read of the report. Registration,
// destruction and the report's read take the mutex; window counters change only
// inside updates, which JAGS runs on the thread that also reads the report
// between runs.
struct SamplerRegistry {
  std::mutex mutex;
  std::vector<CoarseCorrectedSlice *> live;
  bool has_retired = false;
  double retired_min_rate = 0.0;
};

SamplerRegistry &sampler_registry()
{
  static SamplerRegistry registry;
  return registry;
}

class CoarseCorrectedSlice : public Slicer {
  SingletonGraphView const *_graph;
  const unsigned int _chain;
  const SelNormCoarseSettings _settings;
  const double _min_acceptance;
  std::vector<SelectedChild> _selected;
  std::vector<StochasticNode *> _other;
  // Valid only between setValue calls inside one scalar update. Other sampled
  // coordinates may change between updates, so update() always invalidates it.
  mutable bool _surrogate_valid;
  mutable double _surrogate_value;
  // JAGS gives a sampler no accessor for the model's adaptation, so the sampler
  // keeps its own: true from construction until its first adaptOff().
  bool _adapting;
  // The kernel chosen at the first adaptOff(): the exact slice kernel replaces
  // the corrected update for the rest of this object's life.
  bool _exact;
  // Adaptive proposals and accepted corrections, and whether any surrogate
  // evaluated while adapting used a normalizer that is not the exact one.
  std::uint64_t _adaptive_proposed;
  std::uint64_t _adaptive_accepted;
  mutable bool _surrogate_not_exact;
  // Proposals and accepted corrections after adaptOff() in corrected mode since
  // the last clearing read of the report (the window of one run).
  std::uint64_t _window_proposed;
  std::uint64_t _window_accepted;

  // One completed proposal of the corrected update.
  void countOutcome(bool accepted)
  {
    if (_adapting) {
      ++_adaptive_proposed;
      if (accepted) ++_adaptive_accepted;
      return;
    }
    ++_window_proposed;
    report_sampling_proposed.fetch_add(1, std::memory_order_relaxed);
    if (accepted) {
      ++_window_accepted;
      report_sampling_accepted.fetch_add(1, std::memory_order_relaxed);
    }
  }

  // The exact kernel's target: the full conditional, -Inf allowed (outside the
  // support), NaN and +Inf refused as in the corrected update.
  double exactLogDensity() const
  {
    const double full = _graph->logFullConditional(_chain);
    if (std::isnan(full) || full == std::numeric_limits<double>::infinity()) {
      throwNodeError(_graph->node(), "Corrected slice encountered a non-finite full density");
    }
    return full;
  }
public:
  CoarseCorrectedSlice(SingletonGraphView const *graph, unsigned int chain,
                      SelNormCoarseSettings const &settings, double min_acceptance)
    : Slicer(1.0, 10), _graph(graph), _chain(chain), _settings(settings),
      _min_acceptance(min_acceptance), _surrogate_valid(false), _surrogate_value(0.0),
      _adapting(true), _exact(false), _adaptive_proposed(0), _adaptive_accepted(0),
      _surrogate_not_exact(false), _window_proposed(0), _window_accepted(0)
  {
    for (StochasticNode *child : graph->stochasticChildren()) {
      const ChildKind kind = child_kind(child, settings);
      if (kind != ChildKind::other) _selected.emplace_back(child, kind);
      else _other.push_back(child);
    }
    if (_selected.empty()) throw std::logic_error("Corrected slice has no selection child.");
    graph->checkFinite(chain);
    SamplerRegistry &registry = sampler_registry();
    const std::lock_guard<std::mutex> lock(registry.mutex);
    registry.live.push_back(this);
    report_samplers.fetch_add(1, std::memory_order_relaxed);
  }

  // A sampler destroyed before the report is read still counts for the worst
  // rate of its run.
  ~CoarseCorrectedSlice()
  {
    SamplerRegistry &registry = sampler_registry();
    const std::lock_guard<std::mutex> lock(registry.mutex);
    std::vector<CoarseCorrectedSlice *> &live = registry.live;
    for (std::size_t i = 0; i < live.size(); ++i) {
      if (live[i] == this) {
        live[i] = live.back();
        live.pop_back();
        break;
      }
    }
    double rate;
    if (window_rate(_window_proposed, _window_accepted, &rate) &&
        (!registry.has_retired || rate < registry.retired_min_rate)) {
      registry.retired_min_rate = rate;
      registry.has_retired = true;
    }
  }

  // Window counters of the report (read with the registry mutex held).
  std::uint64_t windowProposed() const { return _window_proposed; }
  std::uint64_t windowAccepted() const { return _window_accepted; }
  void clearWindow()
  {
    _window_proposed = 0;
    _window_accepted = 0;
  }

  double value() const
  {
    return _graph->node()->value(_chain)[0];
  }

  void setValue(double value)
  {
    _surrogate_valid = false;
    _graph->setValue(&value, 1, _chain);
  }

  void getLimits(double *lower, double *upper) const
  {
    _graph->node()->support(lower, upper, 1, _chain);
  }

  // The density the inherited slice transition follows: the surrogate in
  // corrected mode, the full conditional in exact mode.
  double logDensity() const
  {
    if (_exact) return exactLogDensity();
    if (_surrogate_valid) return _surrogate_value;
    // The singleton GraphView decomposition preserves every original prior and
    // every non-selected child. Only selected children's log A is approximated.
    double result = _graph->logPrior(_chain);
    for (StochasticNode *child : _other) result += child->logDensity(_chain, PDF_LIKELIHOOD);
    for (SelectedChild const &child : _selected) {
      result += child.logDensity(_chain, _settings, &_surrogate_not_exact);
    }
    if (std::isnan(result) || result == std::numeric_limits<double>::infinity()) {
      throwNodeError(_graph->node(), "Corrected slice encountered a non-finite surrogate density");
    }
    _surrogate_value = result;
    _surrogate_valid = true;
    return result;
  }

  // The kernel is chosen once, at the first adaptOff(), from the adaptive
  // draws: with the fallback enabled, the exact kernel when (a) fewer than 50
  // proposals were made, (c) no surrogate differed from the exact density by
  // construction, or (b) the acceptance rate is below min_acceptance; the first
  // rule that holds is the recorded reason. Both kernels leave the exact target
  // invariant, so the chain after adaptation is one homogeneous exact kernel.
  void adaptOff() override
  {
    Slicer::adaptOff();
    if (!_adapting) return;
    _adapting = false;
    if (!(_min_acceptance > 0.0)) return;
    std::atomic<std::uint64_t> *reason = nullptr;
    if (_adaptive_proposed < min_adaptive_proposals) {
      reason = &report_fallback_few_proposals;
    } else if (!_surrogate_not_exact) {
      reason = &report_fallback_surrogate_exact;
    } else if (static_cast<double>(_adaptive_accepted) /
               static_cast<double>(_adaptive_proposed) < _min_acceptance) {
      reason = &report_fallback_low_acceptance;
    }
    if (reason != nullptr) {
      _exact = true;
      _surrogate_valid = false;
      reason->fetch_add(1, std::memory_order_relaxed);
    }
  }

  void update(RNG *rng)
  {
    _surrogate_valid = false;
    const double old_value = value();
    if (_exact) {
      // One ordinary slice update of the full conditional: no surrogate, no
      // correction, no counting.
      try {
        if (!updateStep(rng)) {
          throwNodeError(_graph->node(), state() == SLICER_POSINF ?
            "Corrected slice (exact kernel) is at a value with infinite density" :
            "Corrected slice (exact kernel) has a current value inconsistent with the data");
        }
      }
      catch (...) {
        setValue(old_value);
        throw;
      }
      return;
    }
    try {
      // Never reuse a full conditional from an earlier coordinate update.
      const double full_old = _graph->logFullConditional(_chain);
      if (!std::isfinite(full_old)) {
        throwNodeError(_graph->node(), "Corrected slice requires a finite current full density");
      }
      const double surrogate_old = logDensity();
      if (!std::isfinite(surrogate_old)) {
        throwNodeError(_graph->node(), "Corrected slice requires a finite current surrogate density");
      }
      // Exactly one scalar reversible surrogate transition. The inherited width
      // adaptation is used only until JAGS invokes adaptOff().
      if (!updateStep(rng)) {
        throwNodeError(_graph->node(), "Surrogate slice transition failed");
      }
      coarse_proposed.fetch_add(1, std::memory_order_relaxed);
      const double full_new = _graph->logFullConditional(_chain);
      if (std::isnan(full_new) || full_new == std::numeric_limits<double>::infinity()) {
        throwNodeError(_graph->node(), "Corrected slice proposed a non-finite full density");
      }
      if (full_new == -std::numeric_limits<double>::infinity()) {
        setValue(old_value);
        coarse_correction_rejected.fetch_add(1, std::memory_order_relaxed);
        countOutcome(false);
        return;
      }
      const double surrogate_new = logDensity();
      if (!std::isfinite(surrogate_new)) {
        throwNodeError(_graph->node(), "Corrected slice proposed a non-finite surrogate density");
      }
      const long double log_ratio =
        (static_cast<long double>(full_new) - full_old) +
        (static_cast<long double>(surrogate_old) - surrogate_new);
      if (std::isnan(log_ratio)) {
        throwNodeError(_graph->node(), "Corrected slice proposed an undefined acceptance ratio");
      }
      if (log_ratio < 0 &&
          std::log(static_cast<long double>(rng->uniform())) > log_ratio) {
        setValue(old_value);
        coarse_correction_rejected.fetch_add(1, std::memory_order_relaxed);
        countOutcome(false);
      }
      else {
        coarse_accepted.fetch_add(1, std::memory_order_relaxed);
        countOutcome(true);
      }
    }
    catch (...) {
      setValue(old_value);
      throw;
    }
  }
};

// A node is eligible when at least one child has a surrogate evaluator, so a
// node whose blocks would all fall back to the exact density keeps its stock
// sampler instead of paying two exact evaluations per update.
bool hasSelectionChild(SingletonGraphView const &graph, SelNormCoarseSettings const &settings)
{
  for (StochasticNode *child : graph.stochasticChildren()) {
    if (child_kind(child, settings) != ChildKind::other) return true;
  }
  return false;
}

// A node with a BayesTools sampler declaration (a 'dbt_block' or 'dbt_gate'
// marker child, written by BayesTools::JAGS_fit()) is left to the samplers of
// the BayesTools module, whatever the order in which the factories run.
bool hasSamplerDeclaration(StochasticNode *node)
{
  for (StochasticNode const *child : *node->stochasticChildren()) {
    std::string const &name = child->distribution()->name();
    if (name == "dbt_block" || name == "dbt_gate") return true;
  }
  return false;
}

}

CoarseCorrectedSliceFactory::CoarseCorrectedSliceFactory(SelNormCoarseSettings const &settings)
{
  configure(false, settings);
}

void CoarseCorrectedSliceFactory::configure(bool enabled, SelNormCoarseSettings const &settings,
                                            double min_acceptance)
{
  if (!(settings.mean_step >= 0) || !std::isfinite(settings.mean_step) ||
      !(settings.diagonal_step >= 0) || !std::isfinite(settings.diagonal_step) ||
      !(settings.log_weight_step >= 0) || !std::isfinite(settings.log_weight_step) ||
      settings.max_rules < 3 || !(min_acceptance >= 0.0 && min_acceptance <= 1.0)) {
    throw std::invalid_argument("Invalid corrected selection surrogate settings.");
  }
  const std::lock_guard<std::mutex> lock(_configuration_mutex);
  _configuration.enabled = enabled;
  _configuration.settings = settings;
  _configuration.min_acceptance = min_acceptance;
}

CoarseCorrectedSliceConfig CoarseCorrectedSliceFactory::configuration() const
{
  const std::lock_guard<std::mutex> lock(_configuration_mutex);
  return _configuration;
}

std::vector<Sampler *> CoarseCorrectedSliceFactory::makeSamplers(
    std::list<StochasticNode *> const &nodes, Graph const &graph) const
{
  std::vector<Sampler *> samplers;
  const CoarseCorrectedSliceConfig config = configuration();
  if (!config.enabled) return samplers;
  const std::string sampler_name("RoBMA::CoarseCorrectedSlice");
  try {
    for (StochasticNode *node : nodes) {
      // Multivariate, discrete, and structural nodes retain their stock sampler.
      if (node->length() != 1 || node->isDiscreteValued() || node->df() == 0 ||
          isObserved(node) || hasSamplerDeclaration(node)) continue;
      std::unique_ptr<SingletonGraphView> view(new SingletonGraphView(node, graph));
      if (!hasSelectionChild(*view, config.settings)) continue;
      std::vector<std::unique_ptr<MutableSampleMethod> > owned_methods;
      std::vector<MutableSampleMethod *> methods;
      for (unsigned int chain = 0; chain < node->nchain(); ++chain) {
        owned_methods.emplace_back(
          std::make_unique<CoarseCorrectedSlice>(view.get(), chain, config.settings,
                                               config.min_acceptance));
        methods.push_back(owned_methods.back().get());
      }
      // The Sampler base adopts the graph before MutableSampler copies members.
      // Passing view.get() would double-delete if a later constructor member
      // throws. C++17 sequences allocation before evaluating view.release().
      // sampler_name is already constructed, so no argument conversion can
      // allocate after release but before the base constructor adopts the graph.
      // owned_methods retains method ownership until construction succeeds.
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

std::string CoarseCorrectedSliceFactory::name() const
{
  return "RoBMA::CoarseCorrectedSlice";
}

CoarseCorrectedSliceFactory *coarse_corrected_slice_factory()
{
  static CoarseCorrectedSliceFactory factory;
  return &factory;
}

CoarseCorrectedSliceStats coarse_corrected_slice_stats(bool reset)
{
  if (reset) {
    return {
      coarse_proposed.exchange(0, std::memory_order_relaxed),
      coarse_accepted.exchange(0, std::memory_order_relaxed),
      coarse_correction_rejected.exchange(0, std::memory_order_relaxed)
    };
  }
  return {
    coarse_proposed.load(std::memory_order_relaxed),
    coarse_accepted.load(std::memory_order_relaxed),
    coarse_correction_rejected.load(std::memory_order_relaxed)
  };
}

CoarseCorrectedSliceReport coarse_corrected_slice_report(bool clear)
{
  // Each counter is read, and cleared in the same operation when requested.
  const auto take = [clear](std::atomic<std::uint64_t> &counter) {
    return clear ? counter.exchange(0, std::memory_order_relaxed) :
      counter.load(std::memory_order_relaxed);
  };
  CoarseCorrectedSliceReport report;
  report.samplers = take(report_samplers);
  report.fallback_few_proposals = take(report_fallback_few_proposals);
  report.fallback_surrogate_exact = take(report_fallback_surrogate_exact);
  report.fallback_low_acceptance = take(report_fallback_low_acceptance);
  report.sampling_proposed = take(report_sampling_proposed);
  report.sampling_accepted = take(report_sampling_accepted);
  // The lowest rate of one sampler over its window: the destroyed samplers'
  // minimum and the live samplers with enough proposals (fallen-back samplers
  // have none).
  SamplerRegistry &registry = sampler_registry();
  const std::lock_guard<std::mutex> lock(registry.mutex);
  report.has_min_rate = registry.has_retired;
  report.min_sampler_rate = registry.has_retired ? registry.retired_min_rate : 0.0;
  for (CoarseCorrectedSlice *sampler : registry.live) {
    double rate;
    if (window_rate(sampler->windowProposed(), sampler->windowAccepted(), &rate) &&
        (!report.has_min_rate || rate < report.min_sampler_rate)) {
      report.min_sampler_rate = rate;
      report.has_min_rate = true;
    }
    if (clear) sampler->clearWindow();
  }
  if (clear) {
    registry.has_retired = false;
    registry.retired_min_rate = 0.0;
  }
  return report;
}

}
}
