#include "CoarseCorrectedSlice.h"
#include "../distributions/DSELNORMMVSTEP.h"

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
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

std::atomic<std::uint64_t> coarse_proposed(0);
std::atomic<std::uint64_t> coarse_accepted(0);
std::atomic<std::uint64_t> coarse_correction_rejected(0);

struct SelectedChild {
  StochasticNode *node;
  DSELNORMMVSTEP const *distribution;
  std::vector<Node const *> parents;
  std::vector<unsigned int> lengths;
  mutable std::vector<double const *> parameters;
  bool controls_fixed = true;
  mutable bool controls_checked = false;
  mutable bool controls_valid = false;

  SelectedChild(StochasticNode *child, DSELNORMMVSTEP const *dist)
    : node(child), distribution(dist)
  {
    const unsigned int count = dist->npar();
    if (child->parents().size() < count) {
      throw std::logic_error("Selection child is missing distribution parameters.");
    }
    // StochasticNode appends truncation bounds after its distribution parents.
    parents.assign(child->parents().begin(), child->parents().begin() + count);
    lengths.resize(count);
    parameters.resize(count);
    for (unsigned int i = 0; i < count; ++i) {
      lengths[i] = parents[i]->length();
      // Only mean, covariance, weights, and kernel mode remain dynamic.
      if (i != 0 && i != 1 && i != 3 && i != 9 && !parents[i]->isFixed()) {
        controls_fixed = false;
      }
    }
  }

  double logDensity(unsigned int chain, SelNormCoarseSettings const &settings) const
  {
    for (unsigned int i = 0; i < parents.size(); ++i) {
      parameters[i] = parents[i]->value(chain);
    }
    bool valid;
    if (controls_fixed) {
      // The first surrogate call follows initialized-graph checkFinite().
      // This certificate belongs to one SelectedChild in one chain's sampler.
      if (!controls_checked) {
        controls_valid = distribution->checkControlParameterValue(parameters, lengths);
        controls_checked = true;
      }
      valid = controls_valid && distribution->checkDynamicParameterValue(parameters, lengths);
    } else {
      // Any nonfixed control retains the complete original node gate.
      valid = node->checkParentValues(chain);
    }
    if (!valid) return -std::numeric_limits<double>::infinity();
    return distribution->surrogateLogDensity(
      node->value(chain), node->length(), PDF_LIKELIHOOD, parameters, lengths,
      node->lowerLimit(chain), node->upperLimit(chain), settings
    );
  }
};

class CoarseCorrectedSlice : public Slicer {
  SingletonGraphView const *_graph;
  const unsigned int _chain;
  const SelNormCoarseSettings _settings;
  std::vector<SelectedChild> _selected;
  std::vector<StochasticNode *> _other;
  // Valid only between setValue calls inside one scalar update. Other sampled
  // coordinates may change between updates, so update() always invalidates it.
  mutable bool _surrogate_valid;
  mutable double _surrogate_value;
public:
  CoarseCorrectedSlice(SingletonGraphView const *graph, unsigned int chain,
                      SelNormCoarseSettings const &settings)
    : Slicer(1.0, 10), _graph(graph), _chain(chain), _settings(settings),
      _surrogate_valid(false), _surrogate_value(0.0)
  {
    for (StochasticNode *child : graph->stochasticChildren()) {
      DSELNORMMVSTEP const *dist =
        dynamic_cast<DSELNORMMVSTEP const *>(child->distribution());
      if (dist != nullptr) _selected.emplace_back(child, dist);
      else _other.push_back(child);
    }
    if (_selected.empty()) throw std::logic_error("Corrected slice has no selection child.");
    graph->checkFinite(chain);
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

  double logDensity() const
  {
    if (_surrogate_valid) return _surrogate_value;
    // The singleton GraphView decomposition preserves every original prior and
    // every non-selected child. Only selected children's log A is approximated.
    double result = _graph->logPrior(_chain);
    for (StochasticNode *child : _other) result += child->logDensity(_chain, PDF_LIKELIHOOD);
    for (SelectedChild const &child : _selected) result += child.logDensity(_chain, _settings);
    if (std::isnan(result) || result == std::numeric_limits<double>::infinity()) {
      throwNodeError(_graph->node(), "Corrected slice encountered a non-finite surrogate density");
    }
    _surrogate_value = result;
    _surrogate_valid = true;
    return result;
  }

  void update(RNG *rng)
  {
    _surrogate_valid = false;
    const double old_value = value();
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
      // adaptation is used only until JAGS invokes the inherited adaptOff().
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
      }
      else {
        coarse_accepted.fetch_add(1, std::memory_order_relaxed);
      }
    }
    catch (...) {
      setValue(old_value);
      throw;
    }
  }
};

bool hasSelectionChild(SingletonGraphView const &graph)
{
  for (StochasticNode *child : graph.stochasticChildren()) {
    if (dynamic_cast<DSELNORMMVSTEP const *>(child->distribution()) != nullptr) return true;
  }
  return false;
}

}

CoarseCorrectedSliceFactory::CoarseCorrectedSliceFactory(SelNormCoarseSettings const &settings)
{
  configure(false, settings);
}

void CoarseCorrectedSliceFactory::configure(bool enabled, SelNormCoarseSettings const &settings)
{
  if (!(settings.mean_step >= 0) || !std::isfinite(settings.mean_step) ||
      !(settings.diagonal_step >= 0) || !std::isfinite(settings.diagonal_step) ||
      !(settings.log_weight_step >= 0) || !std::isfinite(settings.log_weight_step) ||
      settings.max_rules < 3) {
    throw std::invalid_argument("Invalid corrected selection surrogate settings.");
  }
  const std::lock_guard<std::mutex> lock(_configuration_mutex);
  _configuration.enabled = enabled;
  _configuration.settings = settings;
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
          isObserved(node)) continue;
      std::unique_ptr<SingletonGraphView> view(new SingletonGraphView(node, graph));
      if (!hasSelectionChild(*view)) continue;
      std::vector<std::unique_ptr<MutableSampleMethod> > owned_methods;
      std::vector<MutableSampleMethod *> methods;
      for (unsigned int chain = 0; chain < node->nchain(); ++chain) {
        owned_methods.emplace_back(new CoarseCorrectedSlice(view.get(), chain, config.settings));
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

}
}
