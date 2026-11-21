#include "CoarseCorrectedSlice.h"
#include "CoarseSurrogate.h"
#include "CorrectedSampler.h"

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
#include <stdexcept>
#include <string>
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

class CoarseCorrectedSlice : public Slicer, public CorrectedSamplerWindow {
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
    corrected_sampler_counters.report_sampling_proposed.fetch_add(1, std::memory_order_relaxed);
    if (accepted) {
      ++_window_accepted;
      corrected_sampler_counters.report_sampling_accepted.fetch_add(1, std::memory_order_relaxed);
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
    register_corrected_sampler(this);
  }

  // A sampler destroyed before the report is read still counts for the worst
  // rate of its run.
  ~CoarseCorrectedSlice() override
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
  // A scalar sampler's acceptance rate is judged as it is.
  double rateExponent() const override { return 1.0; }

  double value() const override
  {
    return _graph->node()->value(_chain)[0];
  }

  void setValue(double value) override
  {
    _surrogate_valid = false;
    _graph->setValue(&value, 1, _chain);
  }

  void getLimits(double *lower, double *upper) const override
  {
    _graph->node()->support(lower, upper, 1, _chain);
  }

  // The density the inherited slice transition follows: the surrogate in
  // corrected mode, the full conditional in exact mode.
  double logDensity() const override
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
    if (_adaptive_proposed < corrected_min_adaptive_proposals) {
      reason = &corrected_sampler_counters.report_fallback_few_proposals;
    } else if (!_surrogate_not_exact) {
      reason = &corrected_sampler_counters.report_fallback_surrogate_exact;
    } else if (static_cast<double>(_adaptive_accepted) /
               static_cast<double>(_adaptive_proposed) < _min_acceptance) {
      reason = &corrected_sampler_counters.report_fallback_low_acceptance;
    }
    if (reason != nullptr) {
      _exact = true;
      _surrogate_valid = false;
      reason->fetch_add(1, std::memory_order_relaxed);
    }
  }

  void update(RNG *rng) override
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
      corrected_sampler_counters.coarse_proposed.fetch_add(1, std::memory_order_relaxed);
      const double full_new = _graph->logFullConditional(_chain);
      if (std::isnan(full_new) || full_new == std::numeric_limits<double>::infinity()) {
        throwNodeError(_graph->node(), "Corrected slice proposed a non-finite full density");
      }
      if (full_new == -std::numeric_limits<double>::infinity()) {
        setValue(old_value);
        corrected_sampler_counters.coarse_correction_rejected.fetch_add(1, std::memory_order_relaxed);
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
        corrected_sampler_counters.coarse_correction_rejected.fetch_add(1, std::memory_order_relaxed);
        countOutcome(false);
      }
      else {
        corrected_sampler_counters.coarse_accepted.fetch_add(1, std::memory_order_relaxed);
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

std::vector<Sampler *> CoarseCorrectedSliceFactory::makeSamplers(
    std::list<StochasticNode *> const &nodes, Graph const &graph) const
{
  std::vector<Sampler *> samplers;
  const CorrectedSamplerConfig config = corrected_sampler_configuration();
  if (!config.enabled) return samplers;
  const std::string sampler_name("RoBMA::CoarseCorrectedSlice");
  try {
    for (StochasticNode *node : nodes) {
      // Multivariate, discrete, and structural nodes retain their stock sampler.
      if (node->length() != 1 || node->isDiscreteValued() || !node->fullRank() ||
          hasSamplerDeclaration(node)) continue;
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

}
}
