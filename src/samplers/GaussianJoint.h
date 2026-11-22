#ifndef ROBMA_GAUSSIAN_JOINT_H_
#define ROBMA_GAUSSIAN_JOINT_H_

#include "JointSamplerTag.h"

namespace jags {

class SamplerFactory;

namespace RoBMA {

// Terminal declaration. The six scalar stochastic parents are paired
// indicator/slab leaves in intercept, moderator, heterogeneity order. The
// remaining parents are the three active labels, moderator column, observed
// outcomes, sampling SEs, and the effect direction (+1 or -1).
class GaussianJointTag : public JointSamplerTag {
public:
  GaussianJointTag();
  bool checkParameterLength(const std::vector<unsigned long> &lengths) const override;
};

SamplerFactory *gaussian_joint_factory();

} // namespace RoBMA
} // namespace jags

#endif
