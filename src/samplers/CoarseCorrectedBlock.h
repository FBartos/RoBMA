#ifndef ROBMA_COARSE_CORRECTED_BLOCK_H
#define ROBMA_COARSE_CORRECTED_BLOCK_H

// The corrected block sampler (RoBMA::CoarseCorrectedBlock): a block of scalar
// parameters declared by BayesTools' dbt_block markers, with at least one
// coarse-capable selection child, updated by the corrected block kernel
// (CorrectedBlockKernel.h): one slice sweep on the surrogate density in a
// random direction order, corrected by one Metropolis step against the full
// joint conditional. A group of gamma members that enter the model only
// through their proportions (GammaGroup.h; the Dirichlet weights of a
// heterogeneity allocation) is held as its proportions, with its sum drawn
// from its gamma law.

#include "CorrectedBlockKernel.h"

#include <sampler/SamplerFactory.h>

namespace jags {
namespace RoBMA {

// Groups the still-unclaimed marked nodes by block id as BayesTools' factor
// slice factory does and claims a group when the shared corrected sampler
// configuration (CorrectedSampler.h) is enabled and the group has a
// coarse-capable selection child; every other group is left to BayesTools.
class CoarseCorrectedBlockFactory : public SamplerFactory {
public:
  std::vector<Sampler*> makeSamplers(std::list<StochasticNode*> const &nodes,
                                     Graph const &graph) const override;
  std::string name() const override;
};

// Borrowed process-lifetime factory. Module registration must not delete it.
CoarseCorrectedBlockFactory *coarse_corrected_block_factory();

}
}
#endif
