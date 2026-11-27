#include "../fp-contract-off.h"
#include "CoarseCorrectedBlock.h"
#include "CoarseSurrogate.h"
#include "CorrectedSampler.h"
#include "GammaGroup.h"

#include <BayesTools/sampler_declarations.h>

#include <graph/StochasticNode.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/MutableSampleMethod.h>
#include <sampler/MutableSampler.h>
#include <sampler/Sampler.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

const char *const block_sampler_name = "RoBMA::CoarseCorrectedBlock";

// One chain of a declared block: the corrected block kernel, adaptive in every
// update until the first adaptOff().
class CoarseCorrectedBlock : public MutableSampleMethod {
  CorrectedBlockKernel _kernel;
public:
  CoarseCorrectedBlock(GraphView const *graph, unsigned int chain, SelNormCoarseSettings const &settings,
                       double min_acceptance, CorrectedBlockOptions const &options)
    : _kernel(graph, chain, settings, min_acceptance, options) {}

  void update(RNG *rng) override
  {
    corrected_kind_counters(corrected_kind_block).updates.fetch_add(1, std::memory_order_relaxed);
    _kernel.update(rng);
  }

  bool isAdaptive() const override
  {
    return true;
  }

  // The kernel is chosen once, at the first adaptOff() (later calls change
  // nothing; JAGS calls adaptOff() repeatedly).
  void adaptOff() override
  {
    _kernel.adapt_off();
  }

  // As BayesTools' factor slice sampler: incomplete until one adaptive update.
  bool checkAdaptation() const override
  {
    return _kernel.axis_updates() > 0;
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
  const CorrectedTestFlags flags = corrected_test_flags();

  // The members of BayesTools' factor slice factory (BTFactorSlice.cc): scalar,
  // continuous, unobserved full-rank nodes with a fixed support with
  // lower < upper, without a gate marker, with one valid block marker.
  std::map<unsigned int, std::vector<StochasticNode *> > blocks;
  for (StochasticNode *node : nodes) {
    if (node->length() != 1 || node->isDiscreteValued() || !node->fullRank() ||
        !isSupportFixed(node)) continue;
    if (BayesTools::sampler_declarations::has_gate_marker(node)) continue;
    const unsigned int id = BayesTools::sampler_declarations::block_id(node);
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
      CorrectedBlockOptions options;
      options.label = "Corrected slice (block)";
      options.kind = corrected_kind_block;
      options.fixed_order = corrected_block_test_fixed_order();
      options.drop_sum_jacobian = flags.drop_sum_jacobian;
      // the view's member order (the GraphView sorts its nodes)
      options.groups = (flags.disable_new || flags.disable_sum) ? std::vector<GammaGroup>() : find_gamma_groups(view->nodes());
      std::vector<std::unique_ptr<MutableSampleMethod> > owned_methods;
      std::vector<MutableSampleMethod *> methods;
      for (unsigned int chain = 0; chain < block.second[0]->nchain(); ++chain) {
        owned_methods.emplace_back(new CoarseCorrectedBlock(view.get(), chain, config.settings,
                                                            config.min_acceptance, options));
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

}
}
