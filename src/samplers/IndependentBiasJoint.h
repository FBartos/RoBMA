#ifndef ROBMA_INDEPENDENT_BIAS_JOINT_H
#define ROBMA_INDEPENDENT_BIAS_JOINT_H

#include "JointSamplerTag.h"

#include <sampler/SamplerFactory.h>

namespace jags::RoBMA {

// Terminal deterministic tag. The parents are the two indicators,
// raw location, slab label, location sign, branch families, branch-specific
// raw slopes, and branch-major weights on the common selection partition.
// The moderated variant also takes the observed nodes in singleton order,
// the bound physical moderator coefficients and the row-major fixed design.
class IndependentBiasJointTag : public JointSamplerTag {
public:
  explicit IndependentBiasJointTag(bool moderated = false);
  bool checkParameterLength(std::vector<unsigned long> const &) const override;
};

class IndependentBiasJointFactory : public SamplerFactory {
public:
  std::vector<Sampler *> makeSamplers(std::list<StochasticNode *> const &nodes,
                                     Graph const &graph) const override;
  std::string name() const override;
};

} // namespace jags::RoBMA

#endif
