#include "../fp-contract-off.h"
#include "GatedCorrectedBlock.h"
#include "CoarseSurrogate.h"
#include "CorrectedBlockKernel.h"
#include "CorrectedSampler.h"
#include "GammaGroup.h"
#include "StepRowCache.h"

#include <BayesTools/factor_slice_kernel.h>
#include <BayesTools/sampler_declarations.h>

#include <distribution/Distribution.h>
#include <graph/DeterministicNode.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/MutableSampleMethod.h>
#include <sampler/MutableSampler.h>
#include <sampler/Sampler.h>

#include <cmath>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <JRmath.h>

namespace jags {
namespace RoBMA {
namespace {

const char *const gated_block_name = "RoBMA::GatedCorrectedBlock";
// redraws of a closed gate's prior draw before the previous values are kept
const unsigned int closed_draw_attempts = 100;

// Same kind of value: equal, or both NaN (two infinities of one sign are equal).
bool same_log_density(double before, double after)
{
  return before == after || (std::isnan(before) && std::isnan(after));
}

class GatedCorrectedBlock : public MutableSampleMethod {
  GraphView const *_graph;
  const unsigned int _chain;
  Node const *_active;
  CorrectedBlockKernel _kernel;
  std::vector<double> _lower, _upper;
  bool _updated;

  // Summed log density of the stochastic children of the block, the markers
  // included (the likelihood part of the full conditional).
  double children_log_density() const
  {
    double out = 0.0;
    for (StochasticNode const *child : _graph->stochasticChildren()) {
      out += child->logDensity(_chain, PDF_LIKELIHOOD);
    }
    return out;
  }

  bool descendants_finite() const
  {
    for (DeterministicNode const *child : _graph->deterministicChildren()) {
      double const *values = child->value(_chain);
      for (unsigned long i = 0; i < child->length(); ++i) {
        if (!std::isfinite(values[i])) return false;
      }
    }
    return true;
  }

  // Gate closed: by declaration every likelihood child is constant in the
  // members, so their full conditional is their own distribution given their
  // parents and a direct draw is an exact Gibbs update. A draw that puts a
  // member on (or outside) a bound of its support or makes a deterministic
  // descendant non-finite is outside the region where the model is defined
  // (and where the open update can start); it is redrawn, up to a fixed number
  // of times (each attempt an independent draw, so the kernel is a draw from
  // the restricted law or, with a probability that does not depend on the
  // state, no move), and otherwise the previous values are kept.
  // Any exception restores the block-entry values (and, through the graph
  // view, their deterministic descendants) before it propagates, as the open
  // update does.
  void closed_update(RNG *rng)
  {
    std::vector<StochasticNode *> const &nodes = _graph->nodes();
    const unsigned int d = static_cast<unsigned int>(nodes.size());
    std::vector<double> previous(d), drawn(d);
    _graph->getValue(previous, _chain);
    try {
      const bool adapting = _kernel.adapting();
      const double before = adapting ? children_log_density() : 0.0;
      bool valid = false;
      for (unsigned int attempt = 0; attempt < closed_draw_attempts && !valid; ++attempt) {
        for (unsigned int i = 0; i < d; ++i) {
          nodes[i]->randomSample(rng, _chain);
          drawn[i] = nodes[i]->value(_chain)[0];
        }
        _graph->setValue(drawn, _chain);
        valid = true;
        for (unsigned int i = 0; i < d && valid; ++i) {
          valid = std::isfinite(drawn[i]) && drawn[i] > _lower[i] && drawn[i] < _upper[i];
        }
        valid = valid && descendants_finite();
      }
      if (!valid) _graph->setValue(previous, _chain);
      // Diagnostic for declaration errors during adaptation only; the
      // declaration step in R establishes the condition.
      if (adapting && !same_log_density(before, children_log_density())) {
        throwNodeError(nodes[0], "gated node changes the likelihood while its gate is closed");
      }
    }
    catch (...) {
      _graph->setValue(previous, _chain);
      throw;
    }
  }

public:
  GatedCorrectedBlock(GraphView const *graph, unsigned int chain, Node const *active,
                      SelNormCoarseSettings const &settings, double min_acceptance,
                      CorrectedBlockOptions const &options)
    : _graph(graph), _chain(chain), _active(active),
      _kernel(graph, chain, settings, min_acceptance, options), _updated(false)
  {
    for (StochasticNode *node : graph->nodes()) {
      double lower, upper;
      node->support(&lower, &upper, 1, chain);
      _lower.push_back(lower);
      _upper.push_back(upper);
    }
  }

  void update(RNG *rng) override
  {
    _updated = true;
    corrected_kind_counters(corrected_kind_gated).updates.fetch_add(1, std::memory_order_relaxed);
    if (_active->value(_chain)[0] == 0.0) {
      closed_update(rng);
      return;
    }
    _kernel.update(rng);
  }

  bool isAdaptive() const override
  {
    return true;
  }

  void adaptOff() override
  {
    _kernel.adapt_off();
  }

  // Either gate state counts: a branch that was never open while adapting
  // keeps its prior-scale axes and does not hold up the adaptation.
  bool checkAdaptation() const override
  {
    return _updated;
  }
};

// Whether a member depends on another member (through any path).
bool members_depend(std::vector<StochasticNode *> const &members)
{
  for (StochasticNode const *member : members) {
    for (StochasticNode const *other : members) {
      if (member != other && BayesTools::sampler_declarations::depends_on(member, other)) return true;
    }
  }
  return false;
}

// A coarse-capable selection child or a singleton step-selection row exists,
// and every member reaches every such child only through the child's weights:
// no other parameter of such a child is a member or a deterministic descendant
// of the members.
bool weights_only(GraphView const &view, SelNormCoarseSettings const &settings)
{
  std::set<Node const *> reached(view.nodes().begin(), view.nodes().end());
  for (DeterministicNode const *node : view.deterministicChildren()) reached.insert(node);
  bool capable = false;
  const bool rows = !step_row_cache_test_flags().disable;
  for (StochasticNode const *child : view.stochasticChildren()) {
    const ChildKind kind = child_kind(child, settings);
    if (kind == ChildKind::other && rows && is_step_switch_row(child)) {
      // a singleton step-selection row of the row cache, reached only through
      // its weights (parameter 4)
      capable = true;
      std::vector<Node const *> const &parents = child->parents();
      for (unsigned int i = 0; i < parents.size(); ++i) {
        if (i != 4 && reached.count(parents[i]) > 0) return false;
      }
      continue;
    }
    if (kind == ChildKind::other) continue;
    capable = true;
    const unsigned int weights = SelectedChild::weights_parent(kind);
    std::vector<Node const *> const &parents = child->parents();
    for (unsigned int i = 0; i < parents.size(); ++i) {
      if (i != weights && reached.count(parents[i]) > 0) return false;
    }
  }
  return capable;
}

// 2.5 prior standard deviations of each direction of the block (the scale of
// the learned widths): log of a gamma member, logit of a beta member, the
// log-ratio coordinates of a gamma group; 1 for other members.
std::vector<double> prior_widths(std::vector<StochasticNode *> const &nodes,
                                 std::vector<GammaGroup> const &groups)
{
  const double multiplier = ::BayesTools::factor_slice::width_multiplier;
  std::vector<bool> grouped(nodes.size(), false);
  for (GammaGroup const &group : groups) {
    for (unsigned int index : group.members) grouped[index] = true;
  }
  std::vector<double> out;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (grouped[i]) continue;
    StochasticNode const *node = nodes[i];
    const std::string name = node->distribution()->name();
    std::vector<Node const *> const &parents = node->parents();
    double variance = NAN;
    if (name == "dgamma" && parents.size() >= 1) {
      variance = trigamma(parents[0]->value(0)[0]);
    } else if (name == "dbeta" && parents.size() >= 2) {
      variance = trigamma(parents[0]->value(0)[0]) + trigamma(parents[1]->value(0)[0]);
    }
    out.push_back(std::isfinite(variance) && variance > 0 ? multiplier * std::sqrt(variance) : 1.0);
  }
  for (GammaGroup const &group : groups) {
    const std::size_t m = group.shapes.size();
    for (std::size_t i = 0; i + 1 < m; ++i) {
      const double variance = trigamma(group.shapes[i]) + trigamma(group.shapes[m - 1]);
      out.push_back(std::isfinite(variance) && variance > 0 ? multiplier * std::sqrt(variance) : 1.0);
    }
  }
  return out;
}

}

std::vector<Sampler *> GatedCorrectedBlockFactory::makeSamplers(
    std::list<StochasticNode *> const &nodes, Graph const &graph) const
{
  std::vector<Sampler *> samplers;
  const CorrectedSamplerConfig config = corrected_sampler_configuration();
  if (!config.enabled) return samplers;
  const CorrectedTestFlags flags = corrected_test_flags();
  if (flags.disable_new || flags.disable_gated) return samplers;

  // The gated nodes of BayesTools' gated slice factory, grouped by their
  // <active> node in the order of first appearance.
  std::vector<std::pair<Node const *, std::vector<StochasticNode *> > > groups;
  for (StochasticNode *node : nodes) {
    if (node->length() != 1 || node->isDiscreteValued() || !node->fullRank() ||
        !isSupportFixed(node)) continue;
    Node const *active = BayesTools::sampler_declarations::gate_activity(node);
    if (active == nullptr) continue;
    double lower, upper;
    node->support(&lower, &upper, 1, 0);
    if (!(lower < upper)) continue;
    std::size_t g = 0;
    while (g < groups.size() && groups[g].first != active) ++g;
    if (g == groups.size()) groups.push_back(std::make_pair(active, std::vector<StochasticNode *>()));
    groups[g].second.push_back(node);
  }

  const std::string sampler_name(gated_block_name);
  try {
    for (auto const &group : groups) {
      if (members_depend(group.second)) continue;
      std::unique_ptr<GraphView> view(new GraphView(group.second, graph, true));
      if (!weights_only(*view, config.settings)) continue;
      CorrectedBlockOptions options;
      options.label = "Corrected slice (gated block)";
      options.kind = corrected_kind_gated;
      options.fixed_order = corrected_block_test_fixed_order();
      options.drop_sum_jacobian = flags.drop_sum_jacobian;
      options.groups = flags.disable_sum ? std::vector<GammaGroup>() : find_gamma_groups(view->nodes());
      options.initial_widths = prior_widths(view->nodes(), options.groups);
      options.row_cache = !step_row_cache_test_flags().disable;
      std::vector<std::unique_ptr<MutableSampleMethod> > owned_methods;
      std::vector<MutableSampleMethod *> methods;
      for (unsigned int chain = 0; chain < group.second[0]->nchain(); ++chain) {
        owned_methods.emplace_back(new GatedCorrectedBlock(view.get(), chain, group.first, config.settings,
                                                           config.min_acceptance, options));
        methods.push_back(owned_methods.back().get());
      }
      // The Sampler base adopts the graph, and owned_methods keeps the methods
      // until construction succeeds.
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

std::string GatedCorrectedBlockFactory::name() const
{
  return gated_block_name;
}

GatedCorrectedBlockFactory *gated_corrected_block_factory()
{
  static GatedCorrectedBlockFactory factory;
  return &factory;
}

}
}
