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
// coarse rung: the product vector rule, a kernel mode not fixed at NORMAL, and,
// for a factor block, rule min(max_rules, rules) of its schedule affordable at
// the uncompacted rank (a compacted block or a forest only costs less). The
// decision reads fixed parents only, so it is the same in every chain. Any
// other cluster or factor child stays an exact term.
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
  if (!rule->isFixed() || rule->value(0)[0] != SELVECTOR_PRODUCT) return ChildKind::other;
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
  // and loadings), weights and kernel mode. Everything else is a control.
  static bool dynamic_parent(ChildKind kind, unsigned int i)
  {
    if (kind == ChildKind::mv) return i == 0 || i == 1 || i == 3 || i == 9;
    return i == 0 || i == 1 || i == 2 || i == 4 || i == 10;
  }

  SelectedChild(StochasticNode *child, ChildKind child_kind)
    : node(child), kind(child_kind), distribution(child->distribution())
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

  double surrogate(unsigned int chain, SelNormCoarseSettings const &settings) const
  {
    double const *x = node->value(chain);
    double const *lower = node->lowerLimit(chain);
    double const *upper = node->upperLimit(chain);
    switch (kind) {
    case ChildKind::mv:
      return static_cast<DSELNORMMVSTEP const *>(distribution)->surrogateLogDensity(
        x, node->length(), PDF_LIKELIHOOD, parameters, lengths, lower, upper, settings);
    case ChildKind::cluster:
      return static_cast<DSELNORMCLUSTERSTEP const *>(distribution)->surrogateLogDensity(
        x, node->length(), PDF_LIKELIHOOD, parameters, lengths, lower, upper, settings);
    default:
      return static_cast<DSELNORMFACTORSTEP const *>(distribution)->surrogateLogDensity(
        x, node->length(), PDF_LIKELIHOOD, parameters, lengths, lower, upper, settings);
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
        controls_valid = checkControl();
        controls_checked = true;
      }
      valid = controls_valid && checkDynamic();
    } else {
      // Any nonfixed control retains the complete original node gate.
      valid = node->checkParentValues(chain);
    }
    if (!valid) return -std::numeric_limits<double>::infinity();
    return surrogate(chain, settings);
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
      const ChildKind kind = child_kind(child, settings);
      if (kind != ChildKind::other) _selected.emplace_back(child, kind);
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
      if (!hasSelectionChild(*view, config.settings)) continue;
      std::vector<std::unique_ptr<MutableSampleMethod> > owned_methods;
      std::vector<MutableSampleMethod *> methods;
      for (unsigned int chain = 0; chain < node->nchain(); ++chain) {
        owned_methods.emplace_back(
          std::make_unique<CoarseCorrectedSlice>(view.get(), chain, config.settings));
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
