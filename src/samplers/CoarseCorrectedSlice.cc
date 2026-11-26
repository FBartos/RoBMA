#include "CoarseCorrectedSlice.h"
#include "CoarseSurrogate.h"
#include "CorrectedSampler.h"

#include <BayesTools/sampler_declarations.h>

#include <sampler/Slicer.h>
#include <sampler/SingletonGraphView.h>
#include <sampler/MutableSampler.h>
#include <sampler/Sampler.h>
#include <graph/StochasticNode.h>
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
  // The kernel chosen at the first adaptOff(): the exact slice kernel replaces
  // the corrected update for the rest of this object's life.
  bool _exact;
  // Whether a surrogate evaluated during the current update used a normalizer
  // that is not the exact one (the update is then evidence).
  mutable bool _surrogate_not_exact;
  // Whether the next exact evaluation is the update's first (the current state).
  mutable bool _guard_next = false;
  // The adaptation flag, the evidence tally and the reporting window.
  CorrectedTally _tally;

  // The exact kernel's target: the full conditional, -Inf allowed (outside the
  // support), NaN and +Inf refused as in the corrected update. The first
  // evaluation of an update (the slicer's, at the current value) carries the
  // precision guard.
  double exactLogDensity() const
  {
    const double full = _guard_next ? currentFullConditional() : _graph->logFullConditional(_chain);
    if (std::isnan(full) || full == std::numeric_limits<double>::infinity()) {
      throwNodeError(_graph->node(), "Corrected slice encountered a non-finite full density");
    }
    return full;
  }

  // The full conditional of the chain's current state, as GraphView's
  // logFullConditional sums it (the prior, then the children in order), with
  // the precision guard of its likelihood terms (CorrectedSampler.h; after
  // adaptation only). A NaN sum is left to JAGS's own diagnosis.
  double currentFullConditional() const
  {
    _guard_next = false;
    LikelihoodMagnitude magnitude;
    double likelihood = 0.0;
    for (StochasticNode const *child : _graph->stochasticChildren()) {
      const double term = child->logDensity(_chain, PDF_LIKELIHOOD);
      magnitude.add(child, term);
      likelihood += term;
    }
    double prior = 0.0;
    prior += _graph->node()->logDensity(_chain, PDF_PRIOR);
    const double full = prior + likelihood;
    if (std::isnan(full)) return _graph->logFullConditional(_chain);
    check_likelihood_magnitude(magnitude, _tally.adapting(), corrected_kind_slice, "Corrected slice");
    return full;
  }
public:
  CoarseCorrectedSlice(SingletonGraphView const *graph, unsigned int chain,
                      SelNormCoarseSettings const &settings, double min_acceptance)
    : Slicer(1.0, 10), _graph(graph), _chain(chain), _settings(settings),
      _min_acceptance(min_acceptance), _surrogate_valid(false), _surrogate_value(0.0),
      _exact(false), _surrogate_not_exact(false), _tally(corrected_kind_slice)
  {
    for (StochasticNode *child : graph->stochasticChildren()) {
      const ChildKind kind = child_kind(child, settings);
      if (kind != ChildKind::other) _selected.emplace_back(child, kind);
      else _other.push_back(child);
    }
    if (_selected.empty()) throw std::logic_error("Corrected slice has no selection child.");
    graph->checkFinite(chain);
    register_corrected_sampler(this, corrected_kind_slice);
  }

  // A sampler destroyed before the report is read still counts for the worst
  // rate of its run.
  ~CoarseCorrectedSlice() override
  {
    unregister_corrected_sampler(this);
  }

  // Window counters of the report (read with the registry mutex held).
  std::uint64_t windowProposed() const override { return _tally.window_proposed(); }
  std::uint64_t windowAccepted() const override { return _tally.window_accepted(); }
  void clearWindow() override { _tally.clear_window(); }
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

  // The kernel is chosen once, at the first adaptOff(), from the evidence
  // proposals made while adapting (CorrectedTally::adapt_off): with the
  // fallback enabled, the exact kernel when at least 50 were made and their
  // acceptance rate is below min_acceptance. Both kernels leave the exact
  // target invariant, so the chain after adaptation is one homogeneous exact
  // kernel.
  void adaptOff() override
  {
    Slicer::adaptOff();
    if (_tally.adapt_off(_min_acceptance, 1.0)) {
      _exact = true;
      _surrogate_valid = false;
    }
  }

  void update(RNG *rng) override
  {
    _surrogate_valid = false;
    _surrogate_not_exact = false;
    corrected_kind_counters(corrected_kind_slice).updates.fetch_add(1, std::memory_order_relaxed);
    const double old_value = value();
    if (_exact) {
      // One ordinary slice update of the full conditional: no surrogate, no
      // correction, no counting.
      try {
        _guard_next = true;
        const bool stepped = updateStep(rng);
        _guard_next = false;
        if (!stepped) {
          throwNodeError(_graph->node(), state() == SLICER_POSINF ?
            "Corrected slice (exact kernel) is at a value with infinite density" :
            "Corrected slice (exact kernel) has a current value inconsistent with the data");
        }
      }
      catch (...) {
        _guard_next = false;
        setValue(old_value);
        throw;
      }
      return;
    }
    try {
      // Never reuse a full conditional from an earlier coordinate update.
      const double full_old = currentFullConditional();
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
      const double full_new = _graph->logFullConditional(_chain);
      if (std::isnan(full_new) || full_new == std::numeric_limits<double>::infinity()) {
        throwNodeError(_graph->node(), "Corrected slice proposed a non-finite full density");
      }
      if (full_new == -std::numeric_limits<double>::infinity()) {
        setValue(old_value);
        _tally.count(false, _surrogate_not_exact);
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
        _tally.count(false, _surrogate_not_exact);
      }
      else {
        _tally.count(true, _surrogate_not_exact);
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
      // Marker children take precedence over this factory.
      if (node->length() != 1 || node->isDiscreteValued() || !node->fullRank() ||
          BayesTools::sampler_declarations::has_sampler_marker_child(node)) continue;
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
